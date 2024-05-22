
`#draft_spike-in-processing.md`
<br />

<details>
<summary><b><font size="+2"><i>Table of contents</i></font></b></summary>
<br />
<!-- MarkdownTOC -->

1. [Examine and clock previous Bash code to tally/calculate alignments](#examine-and-clock-previous-bash-code-to-tallycalculate-alignments)
1. [Run integration testing for `tally_alignments.py`](#run-integration-testing-for-tally_alignmentspy)
1. [Run a simple single test of `calculate_scaling_factor.py`](#run-a-simple-single-test-of-calculate_scaling_factorpy)
1. [Run `calculate_scaling_factor.py` on all Hho1 and Hmo1 BAM files in an efficient manner](#run-calculate_scaling_factorpy-on-all-hho1-and-hmo1-bam-files-in-an-efficient-manner)
1. [Use deepTools `bamCompare` to generate BEDGRAPH files of spike-in-scaled coverage and siQ-ChIP-scaled coverage](#use-deeptools-bamcompare-to-generate-bedgraph-files-of-spike-in-scaled-coverage-and-siq-chip-scaled-coverage)

<!-- /MarkdownTOC -->
</details>
<br />

<a id="examine-and-clock-previous-bash-code-to-tallycalculate-alignments"></a>
## Examine and clock previous Bash code to tally/calculate alignments
<details>
<summary><i>Code: Examine and clock previous Bash code to tally/calculate alignments</i></summary>

```bash
#!/bin/bash

#  Define functions ===========================================================
calc_6f() { awk "BEGIN{ printf \"%.6f\n\", $* }"; }


calc_2f() { awk "BEGIN{ printf \"%.2f\n\", $* }"; }


tally_alignments() {
    local help
    help=$(
cat << EOM
Usage: tally_alignments

Tally alignments within BAM file.

Options:
  -h, --help     Display this help message
  -f, --file     Path to the BAM file (required) [chr]
  -t, --threads  Number of threads for parallel processing (default: 1) [int ≥
                 1]
  -y, --type     Types of alignments to tally: "all", "sc" (S. cerevisiae), or
                 "sp" (S. pombe, the "spike in") (default: "sc") [chr]
  -q, --mapq     Filter BAM to retain alignments with this MAPQ score or
                 greater (default: 1) [int ≥ 0]
  -m, --mito     Include mitochondrial alignments in tally (default: false)
                 [lgl]

Dependencies:
  Samtools: Required for tallying alignments

Example:
  tally_alignments
      --file "path/to/bam"
      --threads 4
      --type "sc"
      --mapq 1
      --mito true
EOM
    )

    if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
        echo "${help}"
        return 0
    fi

    #  Parse arguments
    local file=""
    local threads=1
    local type="sc"
    local mapq=1
    local mito="false"
    local valid_mito=false
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
            -f|--file) file="${2}"; shift 2 ;;
            -t|--threads) 
                threads="${2}"
                if ! [[ "${threads}" =~ ^[0-9]+$ ]] || [[ "${threads}" -lt 1 ]]; then
                    echo "Error: Argument 'threads' must be an int > 0." >&2
                    return 1
                fi
                shift 2
                ;;
            -y|--type) type="${2,,}"; shift 2 ;;
            -q|--mapq) 
                mapq="${2}"
                if ! [[ "${mapq}" =~ ^[0-9]+$ ]] || [[ "${mapq}" -lt 0 ]]; then
                    echo "Error: Argument 'mapq' must be an int ≥ 0." >&2
                    return 1
                fi
                shift 2
                ;;
            -m|--mito) mito="${2,,}"; shift 2
               # Validate mito input
                if [[ "${mito}" =~ ^(true|t|false|f)$ ]]; then
                    valid_mito=true
                else
                    echo "Error: Argument 'mito' must be true/t or false/f." >&2
                    return 1
                fi
                ;;
            *) echo "Error: Argument '${1}' is invalid." >&2
                return 1
                ;;
        esac
    done

    #  Validate that Samtools is in PATH
    if ! command -v samtools &> /dev/null; then
        echo "Error: Samtools is not in PATH." >&2
        return 1
    fi

    #  Validate mandatory argument 'file'
    if [[ -z "${file}" ]]; then
        echo "Error: Argument 'file' must be specified." >&2
        return 1
    fi

    #  Define function to determine inclusion of mitochondria in tally
    function include_mito() {
        [[ "${mito}" == "true" || "${mito}" == "t" ]]
    }

    #  Prepare mitochondrial inclusion based on type
    local SC_mito=""
    local SP_mito=""
    if include_mito; then
        SC_mito="Mito"
        SP_mito="SP_Mito"
    fi

    #  Execute Samtools based on type and include_mito
    case "${type}" in
        "all")
            if include_mito; then
                samtools view -c -@ "${threads}" -F 4 -q "${mapq}" "${file}"
            else
                samtools view -c -@ "${threads}" -F 4 -q "${mapq}" "${file}" \
                    I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI \
                    SP_II_TG SP_I SP_II SP_III SP_MTR
            fi
            ;;
        "sc")
            samtools view -c -@ "${threads}" -F 4 -q "${mapq}" "${file}" \
                I II III IV V VI VII VIII IX X \
                XI XII XIII XIV XV XVI ${SC_mito}
            ;;
        "sp")
            samtools view -c -@ "${threads}" -F 4 -q "${mapq}" "${file}" \
                SP_II_TG SP_I SP_II SP_III SP_MTR ${SP_mito}
            ;;
        *)
            echo "Error: Unrecognized type '${type}'" >&2
            return 1
            ;;
    esac
}


file="IP_Q_Hho1_6337.sort-coord.bam"
threads=4
type="sc"
mapq=1
mito="true"

time tally_alignments \
    --threads "${threads}" \
    --file "${file}" \
    --type "${type}" \
    --mapq "${mapq}" \
    --mito "${mito}"

#  threads=1
# real	0m14.637s
# user	0m18.712s
# sys	0m1.639s

#  threads=4
# real	0m5.686s
# user	0m18.599s
# sys	0m0.976s
```
</details>
<br />
<br />

<a id="run-integration-testing-for-tally_alignmentspy"></a>
## Run integration testing for `tally_alignments.py`
<details>
<summary><i>Code: Run integration testing for `tally_alignments.py`</i></summary>

```bash
#!/bin/bash

#  Function to run a test case
run_test() {
    local description="${1}"
    local file="${2}"
    local threads="${3}"
    local mapq="${4}"
    local include="${5}"
    local expected="${6}"

    if [[ -ne 6 ]]; then
        echo "
        Must provide the following six positional arguments:
        description="${1}"
        file="${2}"
        threads="${3}"
        mapq="${4}"
        include="${5}"
        expected="${6}"
        "

    echo "Running test: ${description}"
    if [[ -n "${include}" ]]; then
        output=$(
            python tally_alignments.py \
                --file "${file}" \
                --threads "${threads}" \
                --mapq "${mapq}" \
                --include "${include}"
            )
    else
        output=$(
            python tally_alignments.py \
                --file "${file}" \
                --threads "${threads}" \
                --mapq "${mapq}"
            )
    fi

    if [[ "${output}" -eq "${expected}" ]]; then
        echo "Test passed: Output (${output}) matches expected (${expected})"
    else
        echo "Test failed: Expected ${expected}, got ${output}"
    fi
}


file="/home/kalavatt/2023_tutorial_ChIP-seq/03_bam/bowtie2/bam/IP_Q_Hho1_6337.sort-coord.bam"

#  Test 1
description_1="Test with MAPQ of 20 with all main and spike-in chromosomes."
threads_1=8
mapq_1=20
include_1=""
expected_1=22649720

#  Test 2
description_2="Test with MAPQ of 1 and the chromosomes I, II, and III."
threads_2=8
mapq_2=1
include_2="I II III"
expected_2=2702762  # 28061262 for all, 0 for "I II III": #DONE Debug this

#  Execute tests
run_test \
    "${description_1}" \
    "${file}" \
    "${threads_1}" \
    "${mapq_1}" \
    "${include_1}" \
    "${expected_1}"

run_test \
    "${description_2}" \
    "${file}" \
    "${threads_2}" \
    "${mapq_2}" \
    "${include_2}" \
    "${expected_2}"


# python tally_alignments.py \
#     --file "${file}" \
#     --threads "${threads_2}" \
#     --mapq "${mapq_2}" #\
#     # --include "${include_5}"
#     # --include "${include_4}"
#     # --include "${include_3}"
#     # --include "${include_2}"

#  Expected (I II III): 2702762
#+ Expected (I): 471264
#+ Expected (all SC sans Mito): 27987372
#+ Expected (all SC): 27998574
#+ Expected (all SC and SP): 28061262

# include_3="I"
# include_4="I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI"
# include_5="${include_4} Mito"
```
</details>
<br />
<br />

<a id="run-a-simple-single-test-of-calculate_scaling_factorpy"></a>
## Run a simple single test of `calculate_scaling_factor.py`
<details>
<summary><i>Code: Run a simple single test of `calculate_scaling_factor.py`</i></summary>

```bash
#!/bin/bash

dir_proj=~/2023_tutorial_ChIP-seq
# bam_IP="${dir_proj}/03_bam/bowtie2/bam/IP_Q_Hho1_6337.sort-coord.bam"
# bam_in="${dir_proj}/03_bam/bowtie2/bam/in_Q_Hho1_6337.sort-coord.bam"

# threads=8
# mapq=1
# include_SC="I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI"
# include_SP="SP_II_TG SP_I SP_II SP_III SP_MTR"

# ml Pysam/0.22.0-GCC-13.2.0
# 
# IP_main=$(
#     python "${dir_proj}/tally_alignments.py" \
#         --file "${bam_IP}" \
#         --threads "${threads}" \
#         --mapq "${mapq}" \
#         --include "${include_SC}"
# )
# IP_spike_in=$(
#     python "${dir_proj}/tally_alignments.py" \
#         --file "${bam_IP}" \
#         --threads "${threads}" \
#         --mapq "${mapq}" \
#         --include "${include_SP}"
# )
# in_main=$(
#     python "${dir_proj}/tally_alignments.py" \
#         --file "${bam_in}" \
#         --threads "${threads}" \
#         --mapq "${mapq}" \
#         --include "${include_SC}"
# )
# in_spike_in=$(
#     python "${dir_proj}/tally_alignments.py" \
#         --file "${bam_in}" \
#         --threads "${threads}" \
#         --mapq "${mapq}" \
#         --include "${include_SP}"
# )

    IP_main=27987372
IP_spike_in=62488
    in_main=14005626
in_spike_in=422022

echo "
    IP_main=${IP_main}
IP_spike_in=${IP_spike_in}
    in_main=${in_main}
in_spike_in=${in_spike_in}
"

# touch calculate_scaling_factor.py
# vi calculate_scaling_factor.py

sf=$(
    python "${dir_proj}/calculate_scaling_factor.py" \
        -ip_m "${IP_main}" \
        -ip_s "${IP_spike_in}" \
        -in_m "${in_main}" \
        -in_s "${in_spike_in}"
)

echo "${sf}"  # 13.495782233484144
```
</details>
<br />
<br />

<a id="run-calculate_scaling_factorpy-on-all-hho1-and-hmo1-bam-files-in-an-efficient-manner"></a>
## Run `calculate_scaling_factor.py` on all Hho1 and Hmo1 BAM files in an efficient manner
<details>
<summary><i>Code: Run `calculate_scaling_factor.py` on all Hho1 and Hmo1 BAM files in an efficient manner</i></summary>

```bash
#!/bin/bash

#  Set location and file variables
dir_proj=~/2023_tutorial_ChIP-seq
dir_bam="${dir_proj}/03_bam/bowtie2/bam"

#  Set parameters
threads=4
mapq=1
include_SC="I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI"
include_SP="SP_II_TG SP_I SP_II SP_III SP_MTR"
results_file="${dir_proj}/scaling_factors_results.csv"

#  Define types (ChIP'd proteins) for iteration
types=("Hho1" "Hmo1")

#  Set the maximum number of processes to run in parallel
max_jobs=4  # 4 threads × 4 jobs = 16 threads in use at one time
job_count=0

#  Check variable assignments (optional)
check_variables=true
if ${check_variables}; then
    echo "
    dir_proj=${dir_proj}
    dir_bam=${dir_bam}
    
    threads=${threads}
    mapq=${mapq}
    include_SC=${include_SC}
    include_SP=${include_SP}
    results_file=${results_file}
    
    types=${types[*]}
    
    max_jobs=${max_jobs}
    job_count=${job_count}
    "
fi

#  Load the most up-to-date Pysam module  #TODO Make a Mamba environment
ml Pysam/0.22.0-GCC-13.2.0

for type in "${types[@]}"; do
    #  Initialize variable for results file
    results_file="${dir_proj}/scaling_factors_results_${type}.csv"
    
    #  Delete results file if present, then write its header
    [[ -f "${results_file}" ]] && rm "${results_file}"
    echo "sample,IP_main,IP_spike_in,in_main,in_spike_in,scaling_factor" \
        > "${results_file}"

    for IP_file in ${dir_bam}/IP_*_${type}*.sort-coord.bam; do
        sample_name="$(basename "${IP_file}" .sort-coord.bam)"
        in_file="${dir_bam}/in_${sample_name#IP_}.sort-coord.bam"

        if [[ ! -f "${in_file}" ]]; then
            echo "Input file for ${sample_name} not found, skipping..."
            continue
        fi

        echo "Processing ${sample_name#IP_} in background..."

        #  Start tally and scaling-factor calculation in the background
        (
            IP_main=$(
                python "${dir_proj}/tally_alignments.py" \
                    --file "${IP_file}" \
                    --threads "${threads}" \
                    --mapq "${mapq}" \
                    --include "${include_SC}"
            )
            IP_spike_in=$(
                python "${dir_proj}/tally_alignments.py" \
                    --file "${IP_file}" \
                    --threads "${threads}" \
                    --mapq "${mapq}" \
                    --include "${include_SP}"
            )
            in_main=$(
                python "${dir_proj}/tally_alignments.py" \
                    --file "${in_file}" \
                    --threads "${threads}" \
                    --mapq "${mapq}" \
                    --include "${include_SC}"
            )
            in_spike_in=$(
                python "${dir_proj}/tally_alignments.py" \
                    --file "${in_file}" \
                    --threads "${threads}" \
                    --mapq "${mapq}" \
                    --include "${include_SP}"
            )

            sf=$(
                python "${dir_proj}/calculate_scaling_factor.py" \
                    -ip_m "${IP_main}" \
                    -ip_s "${IP_spike_in}" \
                    -in_m "${in_main}" \
                    -in_s "${in_spike_in}"
            )

            #  Append tally and scaling factor information to results_file
            echo "${sample_name},${IP_main},${IP_spike_in},${in_main},${in_spike_in},${sf}" \
                >> "${results_file}"
        ) &

        #  Increment job count and wait if maximum parallel jobs are running
        (( job_count++ ))
        if (( job_count >= max_jobs )); then
            wait -n
            (( job_count-- ))
        fi
    done

    #  Wait for all background jobs to complete
    wait

    #  Normalize the scaling factors for a given ChIP'd protein by dividing
    #+ each one by the largest scaling factor in the group
    max_sf=$(cut -d ',' -f 6 "${results_file}" | sort -nr | head -n 1)
    cat "${results_file}" \
        | awk \
            -F ',' \
            -v max_sf="${max_sf}" \
            'BEGIN {
                OFS=","
            } NR == 1 {
                print $0,"normalized_sf"
            } NR > 1 {
                print $0,$6/max_sf
            }' \
                > "${results_file}.tmp"
    mv -f "${results_file}.tmp" "${results_file}"

    echo "Processing for ${type} complete. Results saved to ${results_file}."
done
```
</details>
<br />
<br />

<a id="use-deeptools-bamcompare-to-generate-bedgraph-files-of-spike-in-scaled-coverage-and-siq-chip-scaled-coverage"></a>
## Use deepTools `bamCompare` to generate BEDGRAPH files of spike-in-scaled coverage and siQ-ChIP-scaled coverage
<details>
<summary><i>Code: Use deepTools `bamCompare` to generate BEDGRAPH files of spike-in-scaled coverage and siQ-ChIP-scaled coverage</i></summary>

```bash
#!/bin/bash

#  Define functions ===========================================================
run_bamCompare() {
    local ip_file="${1}"
    local in_file="${2}"
    local sf="${3}"
    local bw="${4}"
    bamCompare \
        --numberOfProcessors "${threads}" \
        --bamfile1 "${ip_file}" \
        --bamfile2 "${in_file}" \
        --scaleFactors "${sf}:1" \
        --operation "${operation}" \
        --skipNAs \
        --skipZeroOverZero \
        --pseudocount "${pseudocount}" \
        --binSize "${bin_size}" \
        --outFileFormat "${format}" \
        --outFileName "${prefix_out}_${bw}"
}


#  Initialize variables and arrays ============================================
#  Set location and file variables
dir_proj=~/2023_tutorial_ChIP-seq
dir_bam="${dir_proj}/03_bam/bowtie2/bam"

#  Set parameters
threads=16
operation="ratio"
pseudocount=1
bin_size=30
format="bigwig"
prefix_out="test_${operation}"

#  Define types (ChIP'd proteins) for iteration
types=("Hho1" "Hmo1")
# types=("Hho1")

#  Set the maximum number of processes to run in parallel
max_jobs=1  # 16 threads × 1 jobs = 16 threads in use at one time
job_count=0


#  Do the main work ===========================================================
#  Set flags
check_functions=true
check_variables=true

#  Check functions (optional)
if ${check_functions}; then
    type run_bamCompare
fi

#  Check variable assignments (optional)
if ${check_variables}; then
    echo "
    dir_proj=${dir_proj}
    dir_bam=${dir_bam}
    
    threads=${threads}
    operation=${operation}
    pseudocount=${pseudocount}
    bin_size=${bin_size}
    format=${format}
    prefix_out=${prefix_out}
    
    types=(${types[*]})
    
    max_jobs=${max_jobs}
    job_count=${job_count}
    "
fi

#  Load the most up-to-date deepTools module  #TODO Make a Mamba environment
ml deepTools/3.5.4.post1-gfbf-2022b

for type in "${types[@]}"; do
    csv="${dir_proj}/scaling_factors_results_${type}.csv"

    for IP_file in ${dir_bam}/IP_*_${type}*.sort-coord.bam; do
        sample_name="$(basename "${IP_file}" .sort-coord.bam)"
        in_file="${dir_bam}/in_${sample_name#IP_}.sort-coord.bam"
        bw="${sample_name#IP_}.bw"

        if [[ ! -f "${in_file}" ]]; then
            echo "Input file for ${sample_name} not found, skipping..."
            continue
        fi

        normalized_sf=$(
            cat "${csv}" \
                | awk \
                    -F ',' \
                    -v sample_name="${sample_name}" \
                    '$1 == sample_name { print $7 }'
        )

        echo "Processing ${sample_name#IP_}..."

        #  Run bamCompare in the background
        run_bamCompare "${IP_file}" "${in_file}" "${normalized_sf}" "${bw}"
    done
done
```
</details>
<br />

<details>
<summary><i>Code: </i></summary>

```bash
#!/bin/bash

# f_in="IP_G1_Hho1_6336.bed.gz"
# f_ou="${f_in/.bed.gz/.nc.bedgraph}"
f_in="IP_G1_Hho1_6336.I.bed"
# f_ou="${f_in/.bed/.nc.bedgraph}"
f_ou="${f_in/.bed/.tc.bedgraph}"
b_sz=30
thrd=4

# python calculate_siQ-ChIP_normalized_coverage.py \
#     -i "${f_in}" \
#     -o "${f_ou}" \
#     -b "${b_sz}" \
#     -t "${thrd}"

python calculate_coverage.py \
    -i "${f_in}" \
    -o "${f_ou}" \
    -b "${b_sz}" \
    -t "${thrd}" #\
    #-n

#  Compare IP_G1_Hho1_6336.I.tc.bedgraph with the results from deepTools
#+ bamCoverage, IP_G1_Hho1_6336.I.bamCoverage-tc.bedgraph
bam="IP_G1_Hho1_6336.sort-coord.bam"
format="bedgraph"
b_sz=30
chr="I"
thrd=4
norm="None"
min_MAPQ=1

echo "
bamCoverage \\
    --bam \"${bam}\" \\
    --outFileName \"${bam/.bam/.${chr}.bamCoverage-tc.bedgraph}\" \\
    --outFileFormat \"bedgraph\" \\
    --binSize \"${b_sz}\" \\
    --region \"${chr}\" \\
    --numberOfProcessors \"${thrd}\" \\
    --normalizeUsing \"${norm}\" \\
    --minMappingQuality \"${min_MAPQ}\"
"

bamCoverage \
    --bam "${bam}" \
    --outFileName "${bam/.bam/.${chr}.bamCoverage-tc.bedgraph}" \
    --outFileFormat "bedgraph" \
    --binSize "${b_sz}" \
    --region "${chr}" \
    --numberOfProcessors "${thrd}" \
    --normalizeUsing "${norm}" \
    --minMappingQuality "${min_MAPQ}"

#  Continuing to suss out what is going on with `bamCoverage --normalizeUsing
#+ "None"`: Comparing `bamCoverage` to the custom script 
#+ `attempt-to-reproduce_bamCoverage-normalizeUsing-None.py`, which attempts to
#+ implement the default normalization performed by `bamCoverage` when argument
#+ `--normalizeUsing "None"` is specified (per the documentation), i.e., the
#+ evaluation of the following expression: `(total number of mapped reads *
#+ average fragment length) / effective genome size`
function attempt() {
    local infile="${1}"
    local outfile="${2}"
    local genome_size="${3}"
    local bin_size="${4}"
    local min_MAPQ="${5}"
    local region="${6}"

    ~/miniconda3/envs/deeptools_env/bin/python attempt-to-reproduce_bamCoverage-normalizeUsing-None.py \
        -i "${infile}" \
        -o "${outfile}" \
        -g "${genome_size}" \
        -b "${bin_size}" \
        -m "${min_MAPQ}" \
        -r "${region}"
}


bam="IP_G1_Hho1_6336.sort-coord.bam"
chr="I"
format="bedgraph"
genome_size=12157105
bin_size=30
min_MAPQ=1
region="I"
threads=4
normalization="None"

attempt \
    "${bam}" \
    "${bam%.bam}.${region}.attempt.${format}" \
    "${genome_size}" \
    "${bin_size}" \
    "${min_MAPQ}" \
    "${region}"

bamCoverage \
    --bam "${bam}" \
    --outFileName "${bam%.bam}.${region}.bamCoverage.${format}" \
    --outFileFormat "${format}" \
    --binSize "${bin_size}" \
    --region "${region}" \
    --numberOfProcessors "${threads}" \
    --normalizeUsing "${normalization}" \
    --effectiveGenomeSize "${genome_size}" \
    --minMappingQuality "${min_MAPQ}"

bamCoverage \
    --bam "${bam}" \
    --outFileName "${bam%.bam}.${region}.bamCoverage-exactScaling.${format}" \
    --outFileFormat "${format}" \
    --binSize "${bin_size}" \
    --region "${region}" \
    --numberOfProcessors "${threads}" \
    --normalizeUsing "${normalization}" \
    --effectiveGenomeSize "${genome_size}" \
    --minMappingQuality "${min_MAPQ}" \
    --exactScaling

bamCoverage \
    --bam "${bam}" \
    --outFileName "${bam%.bam}.${region}.bamCoverage-EGS-not-specified.${format}" \
    --outFileFormat "${format}" \
    --binSize "${bin_size}" \
    --region "${region}" \
    --numberOfProcessors "${threads}" \
    --normalizeUsing "${normalization}" \
    --minMappingQuality "${min_MAPQ}"


attempt \
    "${bam}" \
    "${bam%.bam}.all.attempt.${format}" \
    "${genome_size}" \
    "${bin_size}" \
    "${min_MAPQ}" \
    "all"

bamCoverage \
    --bam "${bam}" \
    --outFileName "${bam%.bam}.all.bamCoverage.${format}" \
    --outFileFormat "${format}" \
    --binSize "${bin_size}" \
    --numberOfProcessors "${threads}" \
    --normalizeUsing "${normalization}" \
    --effectiveGenomeSize "${genome_size}" \
    --minMappingQuality "${min_MAPQ}"
```
</details>
<br />