
`#draft_spike-in-processing.md`
<br />

<details>
<summary><b><font size="+2"><i>Table of contents</i></font></b></summary>
<br />
<!-- MarkdownTOC -->

1. [1. Tally/calculate alignments](#1-tallycalculate-alignments)
    1. [Code](#code)

<!-- /MarkdownTOC -->
</details>
<br />

<a id="1-tallycalculate-alignments"></a>
# 1. Tally/calculate alignments
<a id="code"></a>
## Code
<details>
<summary><i>Code: 1. Tally/calculate alignments</i></summary>

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

file="/home/kalavatt/2023_tutorial_ChIP-seq/03_bam/bowtie2/bam/IP_Q_Hho1_6337.sort-coord.bam"
threads=8
mapq=20
include="I II III"

tally_alignments.py \
    --file "${file}" \
    --threads "${threads}" \
    --mapq "${mapq}"

file="/home/kalavatt/2023_tutorial_ChIP-seq/03_bam/bowtie2/bam/IP_Q_Hho1_6337.sort-coord.bam"
threads=4
mapq=1
include="I II III"

tally_alignments.py \
    --file "${file}" \
    --threads "${threads}" \
    --mapq "${mapq}" \
    --include "${include}"

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

<details>
<summary><i>Code: Unit testing for tally_alignments.py</i></summary>

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

    echo "Running test: ${description}"
    if [[ -n "${include}" ]]; then
        output=$(
            python tally_alignments.py \
                --file "$file" \
                --threads "$threads" \
                --mapq "$mapq" \
                --include "$include"
            )
    else
        output=$(
            python tally_alignments.py \
                --file "$file" \
                --threads "$threads" \
                --mapq "$mapq"
            )
    fi

    if [[ "${output}" -eq "${expected}" ]]; then
        echo "Test passed: Output (${output}) matches expected (${expected})"
    else
        echo "Test failed: Expected ${expected}, got ${output}"
    fi
}

#  Test 1
description_1="Test with high MAPQ without specifying include chromosomes."
file_1="/home/kalavatt/2023_tutorial_ChIP-seq/03_bam/bowtie2/bam/IP_Q_Hho1_6337.sort-coord.bam"
threads_1=8
mapq_1=20
include_1=""
expected_1=500 # You need to set this based on your expected results

#  Test 2
description_2="Test with low MAPQ and specific include chromosomes."
file_2="/home/kalavatt/2023_tutorial_ChIP-seq/03_bam/bowtie2/bam/IP_Q_Hho1_6337.sort-coord.bam"
threads_2=4
mapq_2=1
include_2="I II III"
expected_2=150 # You need to set this based on your expected results

#  Execute tests
run_test \
    "${description_1}" \
    "${file_1}" \
    "${threads_1}" \
    "${mapq_1}" \
    "${include_1}" \
    "${expected_1}"

run_test \
    "${description_2}" \
    "${file_2}" \
    "${threads_2}" \
    "${mapq_2}" \
    "${include_2}" \
    "${expected_2}"
```
</details>
<br />
