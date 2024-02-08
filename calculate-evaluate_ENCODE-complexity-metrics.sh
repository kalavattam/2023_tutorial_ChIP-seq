#!/bin/bash

#  calculate-evaluate_ENCODE-complexity-metrics.sh
#  KA

debug=true


#  Define functions ===========================================================
#TODO Function to print hep information for the script
function show_help() {
    cat << EOM
Usage: ${0} [OPTIONS]

This script runs analyses to calculate, evaluate, and return various ENCODE
library complexity metrics: encodeproject.org/data-standards/terms/#library

Options:
  -h, --help     Display this help message.
  -t, --threads  Number of threads to use (default: 1).
  -b, --bam      Path to the BAM file.
  -m, --mode     Mode to specify chromosome filter: SC (exclude
                 "SP_"), SP (only "SP_"), or both (all chromosomes) (default:
                 SC).
  -q, --mapq     MAPQ threshold for BAM filtering (default: 1).

Dependencies:
  samtools  Required for processing BAM files.
  sort      Involved in generating various alignment tallies.
  uniq      Involved in generating various alignment tallies.
  wc        Involved in generating various alignment tallies.
  awk       Used for filtering BAM file text streams and generating various
            alignment tallies.
  bc        Used for evaluating metric values.

Example:
  calculate-evaluate_ENCODE-complexity-metrics.sh \\
      --bam path/to/bam \\
      --threads 4 \\
      --mode SC \\
      --mapq 1
EOM
}


#  Function to check if PROGRAM; if not, exit with code 1
function check_program_in_path() {
    local program_name="$1"
    if ! command -v "${program_name}" &> /dev/null; then
        error_and_exit "${program_name} is not in PATH. Please install ${program_name} or add it to PATH to continue."
    fi
}


#  Function to exit with exit code 1, which stops the execution of all code,
#+ and return an error message
function error_and_exit() {
    echo "Error: ${1}" >&2
    exit 1
}


#  Function to calculate metric based on user-supplied conditions
function calculate_metric() {
    local threads="${1}"
    local mapq="${2}"
    local bam="${3}"
    local condition="${4}"
    local count_condition="${5}"

    samtools view -@ "${threads}" -f 2 -q "${mapq}" "${bam}" \
        | awk "${condition}" \
        | sort \
        | uniq -c \
        | awk "${count_condition}"
}


#  Function to evaluate NRF or PBC1/2 metrics
function evaluate_metric() {
    local metric="${1}"  # The metric to evaluate (NRF, PBC1, or PBC2)
    local metric_value="${2}"  # The value of the metric
    
    #  Below, use 'case' with 'bc' for floating-point comparisons to categorize
    #+ the metric value. 'bc -l' is a calculator that supports floating-point
    #+ arithmetic. Each 'case' option uses 'bc' to evaluate if the metric_value
    #+ meets a specific condition and echoes the corresponding evaluation
    #+ (e.g., "severe", "moderate", "mild", or "none")
    if [[ "${metric}" == "NRF" ]]; then
        # shellcheck disable=SC2194
        # shellcheck disable=SC2254
        case 1 in
            $(
                echo "${metric_value} < 0.5" | bc -l
            )) echo "concerning" ;;
            $(
                echo "${metric_value} >= 0.5 && ${metric_value} < 0.8" | bc -l
            )) echo "acceptable" ;;
            $(
                echo "${metric_value} >= 0.8 && ${metric_value} < 0.9" | bc -l
            )) echo "compliant" ;;
            $(
                echo "${metric_value} >= 0.9" | bc -l
            )) echo "ideal" ;;
        esac
    elif [[ "${metric}" == "PBC1" ]]; then
        #  For PBC1, different evaluation terms are returned

        # shellcheck disable=SC2194
        # shellcheck disable=SC2254
        case 1 in
            $(
                echo "${metric_value} < 0.5" | bc -l
            )) echo "severe" ;;
            $(
                echo "${metric_value} >= 0.5 && ${metric_value} < 0.8" | bc -l
            )) echo "moderate" ;;
            $(
                echo "${metric_value} >= 0.8 && ${metric_value} < 0.9" | bc -l
            )) echo "mild" ;;
            $(
                echo "${metric_value} >= 0.9" | bc -l
            )) echo "none" ;;
        esac
    elif [[ "${metric}" == "PBC2" ]]; then
        #  For PBC2, different thresholds are applied and different evaluation
        #+ terms are returned
        
        # shellcheck disable=SC2194
        # shellcheck disable=SC2254
        case 1 in
            $(
                echo "${metric_value} < 1" | bc -l
            )) echo "severe" ;;
            $(
                echo "${metric_value} >= 1 && ${metric_value} < 3" | bc -l
            )) echo "moderate" ;;
            $(
                echo "${metric_value} >= 3 && ${metric_value} < 10" | bc -l
            )) echo "mild" ;;
            $(
                echo "${metric_value} >= 10" | bc -l
            )) echo "none" ;;
        esac
    else
        error_and_exit "Invalid positional argument: ${metric}. Argument must be 'NRF', 'PBC1', or 'PBC2'."
    fi
}


#  Something ==================================================================
check_program_in_path "samtools"

if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
    show_help
    exit 0
fi

if ${debug}; then
    #  Values for debugging
    threads=1
    bam="03_bam/bowtie2/bam/in_G1_Hho1_6336.sort-coord.bam"
    mode="SC"
    mapq=1

    echo "
    threads=${threads}
    bam=${bam}
    mode=${mode}
    mapq=${mapq}
    "
else
    #  Default values
    threads=1
    bam=""
    mode="SC"
    mapq=1
fi

#  Parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case "${1}" in
        -t|--threads) threads="${2}"; shift 2 ;;
        -b|--bam) bam="${2}"; shift 2 ;;
        -m|--mode) mode="${2}"; shift 2 ;;
        -q|--mapq) mapq="${2}"; shift 2 ;;
        *) echo "Unknown parameter passed: ${1}"; exit 1 ;;
    esac
done

if [[ -z "${bam}" ]]; then
    error_and_exit "BAM file path is required. Use -b or --bam to specify it."
fi

if [[ ! -f "${bam}" ]]; then
    error_and_exit "Specified BAM file does not exist."
fi

#  Check mode and set awk command accordingly
# shellcheck disable=SC2016
case "${mode}" in
    SP) awk_cmd='{ if ($3 ~ /^SP_/) print $3, $4 }' ;;
    SC) awk_cmd='{ if ($3 !~ /^SP_/) print $3, $4 }' ;;
    both) awk_cmd='{ print $3, $4 }' ;;
    *) error_and_exit "Invalid mode: ${mode}. Mode must be 'SC', 'SP', or 'both'." ;;
esac

#  Determine M_1, the tally of genomic locations where exactly one read maps
#+ uniquely
# shellcheck disable=SC2016
M_1=$(
    calculate_metric \
        "${threads}" "${mapq}" "${bam}" \
        "${awk_cmd}" '$1 == 1 { count++ } END { print count }'
)

#  Determine M_2, the tally of genomic locations where exactly two reads map
#+ uniquely
# shellcheck disable=SC2016
M_2=$(
    calculate_metric \
        "${threads}" "${mapq}" "${bam}" \
        "${awk_cmd}" '$1 == 2 { count++ } END { print count }'
)

#  Determine M_DISTINCT, the count of the number of distinct genomic locations
#+ to which some read maps uniquely
# shellcheck disable=SC2016
M_distinct=$(
    calculate_metric \
        "${threads}" "${mapq}" "${bam}" \
        "${awk_cmd}" '{ print $0 }' \
            | wc -l
)

#  Calculate the ENCODE non-redundant fraction (NRF), the number of distinct
#+ alignments (after removing duplicates) divided by the total number of
#+ alignments
total="$(samtools view -@ "${threads}" -c -f 2 -q "${mapq}" "${bam}")"
dup="$(samtools view -@ "${threads}" -c -f 1024 -f 2 -q "${mapq}" "${bam}")"
nondup="$(( total - dup ))"

NRF="$(echo "scale=9; ${nondup} / ${total}" | bc)"

#  Calculate the ENCODE PCR bottleknecking coefficient #1 (PBC1):
#+ M_1 / M_DISTINCT
PBC1="$(echo "scale=9; ${M_1} / ${M_distinct}" | bc)"

#  Calculate the ENCODE PCR bottleknecking coefficient #2 (PBC2):
#+ M_1 / M_2
PBC2="$(echo "scale=9; ${M_1} / ${M_2}" | bc)"

#  Compare and record ENCODE QC evaluation based on the value of NRF
eval_NRF="$(evaluate_metric "NRF" "${NRF}")"

#  Compare and record ENCODE QC evaluation based on the value of PBC1
eval_PBC1="$(evaluate_metric "PBC1" "${PBC1}")"

#  Compare and record ENCODE QC evaluation based on the value of PBC2
eval_PBC2="$(evaluate_metric "PBC2" "${PBC2}")"

if ${debug}; then
    echo "
           bam=${bam}

         total=${total}
           dup=${dup}
        nondup=${nondup}

           M_1=${M_1}
           M_2=${M_2}
    M_distinct=${M_distinct}
           NRF=${NRF}
          PBC1=${PBC1}
          PBC2=${PBC2}

      eval_NRF=${eval_NRF}
     eval_PBC1=${eval_PBC1}
     eval_PBC2=${eval_PBC2}
    "
else
    echo "bam"$'\t'"total"$'\t'"dup"$'\t'"nondup"$'\t'"M_1"$'\t'"M_2"$'\t'"M_distinct"$'\t'"NRF"$'\t'"PBC1"$'\t'"PBC2"$'\t'"eval_NRF"$'\t'"eval_PBC1"$'\t'"eval_PBC2"
    echo "$(basename "${bam}")"$'\t'"${total}"$'\t'"${dup}"$'\t'"${nondup}"$'\t'"${M_1}"$'\t'"${M_2}"$'\t'"${M_distinct}"$'\t'"${NRF}"$'\t'"${PBC1}"$'\t'"${PBC2}"$'\t'"${eval_NRF}"$'\t'"${eval_PBC1}"$'\t'"${eval_PBC2}"
fi

#  SC mode
#
#         M_1_1=2902897
#         M_2_1=2189664
#  M_distinct_1=8747978
#         NRF_1=.685554002
#        PBC1_1=.331836339
#        PBC2_1=1.325727143
#
#        M_1_40=2725652
#        M_2_40=2020843
# M_distinct_40=8068906
#        NRF_40=.677960672
#       PBC1_40=.337796970
#       PBC2_40=1.348769795
