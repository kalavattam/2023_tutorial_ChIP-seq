
`#work_align_tally_data-Brn1_genome-cat.md`
<br />
<br />

<details>
<summary><b><font size="+2"><i>Table of contents</i></font></b></summary>
<!-- MarkdownTOC -->

1. [Get situated](#get-situated)
    1. [Code](#code)
1. [Align the datasets](#align-the-datasets)
    1. [Run `bowtie2` alignment, etc.](#run-bowtie2-alignment-etc)
        1. [Get situated, set up variables and arrays](#get-situated-set-up-variables-and-arrays)
            1. [Code](#code-1)
        1. [Trim any adapter sequences present in the reads](#trim-any-adapter-sequences-present-in-the-reads)
            1. [Code](#code-2)
        1. [Align, sort, and index the sample datasets](#align-sort-and-index-the-sample-datasets)
            1. [`SmartMap`'s `Bowtie2` parameters](#smartmaps-bowtie2-parameters)
                1. [Code](#code-3)
            1. [Align \(etc.\) `atria`-trimmed `fastq`s](#align-etc-atria-trimmed-fastqs)
                1. [Code](#code-4)
            1. [Align \(etc.\) un-trimmed `fastq`s](#align-etc-un-trimmed-fastqs)
                1. [Code](#code-5)
    1. [Examine flags in bam outfiles](#examine-flags-in-bam-outfiles)
        1. [Initialize necessary functions](#initialize-necessary-functions)
            1. [Code](#code-6)
        1. [Initialize an array of bams](#initialize-an-array-of-bams)
            1. [Code](#code-7)
        1. [Check on flag information in bams](#check-on-flag-information-in-bams)
            1. [Code](#code-8)
1. [Tally/calculate alignments](#tallycalculate-alignments)
    1. [Tally/calculate alignments](#tallycalculate-alignments-1)
        1. [Initialize functions for doing floating point arithmetic, etc.](#initialize-functions-for-doing-floating-point-arithmetic-etc)
            1. [Code](#code-9)
        1. [Get situated, then initialize arrays, variables, etc.](#get-situated-then-initialize-arrays-variables-etc)
            1. [Code](#code-10)
        1. [Generate tab-separated table of alignment tallies/calculations](#generate-tab-separated-table-of-alignment-talliescalculations)
            1. [Code](#code-11)
        1. [Calculate CC/SS-styled scaling factors](#calculate-ccss-styled-scaling-factors)
            1. [Code](#code-12)
    1. [On calculating scaling factors](#on-calculating-scaling-factors)
        1. [Email from Christine \(edited by me\)](#email-from-christine-edited-by-me)
        1. [Notes on using Excel to calculate scaling factors](#notes-on-using-excel-to-calculate-scaling-factors)

<!-- /MarkdownTOC -->
</details>
<br />
<br />

<a id="get-situated"></a>
## Get situated
<a id="code"></a>
### Code
<details>
<summary><i>Code</i></summary>

```bash
#!/bin/bash

run_prev_sym_apch=FALSE
[[ "${run_prev_sym_apch}" == TRUE ]] &&
    {
        p_data="${HOME}/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802"
        p_rename="${p_data}/renamed"

        [[ ! -d "${p_rename}" ]] && mkdir -p "${p_rename}"

        unset data
        typeset -A data=(
            ["${p_data}/SRR7175375.fastq"]="${p_rename}/Brn1_Q_rep1_ChIP.fastq"
            ["${p_data}/SRR7175376.fastq"]="${p_rename}/Brn1_Q_rep1_input.fastq"
            ["${p_data}/SRR7175377.fastq"]="${p_rename}/Brn1_Q_rep2_ChIP.fastq"
            ["${p_data}/SRR7175378.fastq"]="${p_rename}/Brn1_Q_rep2_input.fastq"
            ["${p_data}/SRR7175379.fastq"]="${p_rename}/Brn1_Q_rep3_ChIP.fastq"
            ["${p_data}/SRR7175380.fastq"]="${p_rename}/Brn1_Q_rep3_input.fastq"
            ["${p_data}/SRR7175382.fastq"]="${p_rename}/Brn1_Q_all_input.fastq"
            ["${p_data}/SRR7175381.fastq"]="${p_rename}/Brn1_Q_all_ChIP.fastq"
            ["${p_data}/SRR7175367.fastq"]="${p_rename}/Brn1_log_rep1_ChIP.fastq"
            ["${p_data}/SRR7175368.fastq"]="${p_rename}/Brn1_log_rep1_input.fastq"
            ["${p_data}/SRR7175369.fastq"]="${p_rename}/Brn1_log_rep2_ChIP.fastq"
            ["${p_data}/SRR7175370.fastq"]="${p_rename}/Brn1_log_rep2_input.fastq"
            ["${p_data}/SRR7175371.fastq"]="${p_rename}/Brn1_log_rep3_ChIP.fastq"
            ["${p_data}/SRR7175372.fastq"]="${p_rename}/Brn1_log_rep3_input.fastq"
            ["${p_data}/SRR7175373.fastq"]="${p_rename}/Brn1_log_all_ChIP.fastq"
            ["${p_data}/SRR7175374.fastq"]="${p_rename}/Brn1_log_all_input.fastq"
        )

        #  Rather than rename the initial files, make "symbolic links" to the
        #+ files with new, better-detailed names
        for i in "${!data[@]}"; do
            echo """
              key  ${i}
            value  ${data[${i}]}
            """

            ln -s "${i}" "${data["${i}"]}"
        done
    }

#UPDATE
#  2023-0715-0716: fastqs are now compressed, and symlinked/renamed files are
#+ now stored in subdirectory sym/, e.g.,
p_data="${HOME}/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802"
p_sym="${p_data}/sym"
cd "${p_data}" || echo "cd'ing failed; check on this..."
ls -lhaFG "${p_sym}"

#NOTE
#  For more details on the symlinking/renaming process, see
#+ 2023_rDNA/results/2023-0228_work_fastqs_download/work_download-data.md
```
</details>
<br />
<br />

<a id="align-the-datasets"></a>
## Align the datasets
<a id="run-bowtie2-alignment-etc"></a>
### Run `bowtie2` alignment, etc.
<a id="get-situated-set-up-variables-and-arrays"></a>
#### Get situated, set up variables and arrays
<a id="code-1"></a>
##### Code
<details>
<summary><i>Code</i></summary>

```bash
#!/bin/bash

# tmux new -s a1  # tmux attach -t a1
# tmux new -s a2  # tmux attach -t a2
grabnode  # 8, 160, 1, N

module purge
ml SAMtools/1.16.1-GCC-11.2.0 Bowtie2/2.4.4-GCC-11.2.0

#  Set up work directory, directories for processed data
d_work="${HOME}/tsukiyamalab/Kris/2023_rDNA/results"  # cd "${d_work}"

[[ ! -d "${d_work}/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out" ]] && \
    mkdir -p "${d_work}/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out"

[[ ! -d "${d_work}/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out" ]] && \
    mkdir -p "${d_work}/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out"

#  Go to work directory
cd "${d_work}" || echo "cd'ing failed; check on this"

#  Initialize variables needed for alignment, etc.
  threads="${SLURM_CPUS_ON_NODE}"                                        # echo "${threads}"
  scratch="/fh/scratch/delete30/tsukiyama_t"                             # ., "${scratch}"
 d_genome="${HOME}/tsukiyamalab/Kris/genomes/combined_SC_SP"             # ls -lhaFG "${d_genome}"
 f_genome="${d_genome}/fasta/combined_SC_SP.fa"                          # ls -lhaFG "${f_genome}"
f_indices="${d_genome}/bowtie2/$(basename "${f_genome}" .fa)"            # ls -lhaFG "${f_indices}"*
  err_out="${d_work}/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out"  # ls -lhaFG "${err_out}"
   d_bams="$(dirname ${err_out})"                                        # ls -lhaFG "${d_bams}"

p_data="${HOME}/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym"        # ls -lhaFG "${p_data}"
unset fastqs
typeset -a fastqs=(
    "${p_data}/Ch_log_WT_Brn1_rep1.fastq.gz"
    "${p_data}/Ch_log_WT_Brn1_rep2.fastq.gz"
    "${p_data}/Ch_log_WT_Brn1_rep3.fastq.gz"
    "${p_data}/Ch_log_WT_Brn1_repM.fastq.gz"
    "${p_data}/Ch_Q_WT_Brn1_rep1.fastq.gz"
    "${p_data}/Ch_Q_WT_Brn1_rep2.fastq.gz"
    "${p_data}/Ch_Q_WT_Brn1_rep3.fastq.gz"
    "${p_data}/Ch_Q_WT_Brn1_repM.fastq.gz"
    "${p_data}/in_log_WT_Brn1_rep1.fastq.gz"
    "${p_data}/in_log_WT_Brn1_rep2.fastq.gz"
    "${p_data}/in_log_WT_Brn1_rep3.fastq.gz"
    "${p_data}/in_log_WT_Brn1_repM.fastq.gz"
    "${p_data}/in_Q_WT_Brn1_rep1.fastq.gz"
    "${p_data}/in_Q_WT_Brn1_rep2.fastq.gz"
    "${p_data}/in_Q_WT_Brn1_rep3.fastq.gz"
    "${p_data}/in_Q_WT_Brn1_repM.fastq.gz"
)

p_atria="${d_work}/2023-0406_tutorial_ChIP-seq_analysis/atria"
unset atria
typeset -a atria=(
    "${p_atria}/Ch_log_WT_Brn1_rep1.atria.fastq.gz"
    "${p_atria}/Ch_log_WT_Brn1_rep2.atria.fastq.gz"
    "${p_atria}/Ch_log_WT_Brn1_rep3.atria.fastq.gz"
    "${p_atria}/Ch_log_WT_Brn1_repM.atria.fastq.gz"
    "${p_atria}/Ch_Q_WT_Brn1_rep1.atria.fastq.gz"
    "${p_atria}/Ch_Q_WT_Brn1_rep2.atria.fastq.gz"
    "${p_atria}/Ch_Q_WT_Brn1_rep3.atria.fastq.gz"
    "${p_atria}/Ch_Q_WT_Brn1_repM.atria.fastq.gz"
    "${p_atria}/in_log_WT_Brn1_rep1.atria.fastq.gz"
    "${p_atria}/in_log_WT_Brn1_rep2.atria.fastq.gz"
    "${p_atria}/in_log_WT_Brn1_rep3.atria.fastq.gz"
    "${p_atria}/in_log_WT_Brn1_repM.atria.fastq.gz"
    "${p_atria}/in_Q_WT_Brn1_rep1.atria.fastq.gz"
    "${p_atria}/in_Q_WT_Brn1_rep2.atria.fastq.gz"
    "${p_atria}/in_Q_WT_Brn1_rep3.atria.fastq.gz"
    "${p_atria}/in_Q_WT_Brn1_repM.atria.fastq.gz"
)

run_check=TRUE
[[ "${run_check}" == TRUE ]] &&
    {
        echo '### echo "${threads}" ###'
        echo "${threads}"
        echo ""

        echo '### ls -lhaFG "${scratch} ###'
        ls -lhaFG "${scratch}"
        echo ""

        echo '### ls -lhaFG "${d_genome}" ###'
        ls -lhaFG "${d_genome}"
        echo ""

        echo '### ls -lhaFG "${f_genome}" ###'
        ls -lhaFG "${f_genome}"
        echo ""

        echo '### ls -lhaFG "${f_indices}"* ###'
        ls -lhaFG "${f_indices}"*
        echo ""

        echo '### ls -lhaFG "${err_out}" ###'
        ls -lhaFG "${err_out}"
        echo ""

        echo '### ls -lhaFG "${d_bams}" ###'
        ls -lhaFG "${d_bams}"
        echo ""

        echo '### ls -lhaFG "${p_data}" ###'
        ls -lhaFG "${p_data}"
        echo ""

        echo '### for i in "${fastqs[@]}"; do ls -lhaFG "${i}"; done ###'
        for i in "${fastqs[@]}"; do ls -lhaFG "${i}"; done
        echo ""

        echo '### for i in "${atria[@]}"; do echo "${i}"; done ###'
        for i in "${atria[@]}"; do echo "${i}"; done
    }
```
</details>
<br />

<a id="trim-any-adapter-sequences-present-in-the-reads"></a>
#### Trim any adapter sequences present in the reads
<a id="code-2"></a>
##### Code
<details>
<summary><i>Code</i></summary>

```bash
#!/bin/bash

#  Run print tests to check that the commands are correct/reasonable
print_test=TRUE
[[ "${print_test}" == TRUE ]] &&
    {
        for i in "${fastqs[@]}"; do
            echo """
            \"\${HOME}/tsukiyamalab/kalavatt/2023_rDNA/src/Atria/app-3.2.2/bin/atria\" \\
                -t \"${threads}\" \\
                -r \"${i}\" \\
                -o \"${p_atria}\" \\
                --length-range 35:500
            """
        done
    }

#  If things look good, then run the calls to atria (the program for trimming)
run=TRUE
[[ "${run}" == TRUE ]] &&
    {
        for i in ./2023-0406_tutorial_ChIP-seq_analysis/atria/*.fastq.gz; do
            if [[ ! -e "${i}" ]]; then
                echo "File does not exist; making file"
                #  Initialize environment for atria's dependencies
                [[ "${CONDA_DEFAULT_ENV}" != "atria_env" ]] && \
                    source activate atria_env
                
                for i in "${fastqs[@]}"; do
                    "${HOME}/tsukiyamalab/kalavatt/2023_rDNA/src/Atria/app-3.2.2/bin/atria" \
                        -t "${threads}" \
                        -r "${i}" \
                        -o "${p_atria}" \
                        --length-range 35:500
                done

                #  Move atria logs to atria/err_out/
                [[ $(find "${p_atria}" -type f \( -name "*.log" -o -name "*.json" \) | wc -l) -gt 0 ]] && \
                    mv ${p_atria}/*.{log,json} "${p_atria}/err_out"
            else
                echo "Atria-trimmed fastqs exist; skipping the running of Atria"
                break
            fi

        done
    }

#  Check that "${atria[@]}" array elements exist/are recognized by the system
run_check=TRUE
[[ "${run_check}" == TRUE ]] &&
    {
        echo '### for i in "${atria[@]}"; do ls -lhaFG "${i}"; done ###'
        for i in "${atria[@]}"; do ls -lhaFG "${i}"; done
    }
```
</details>
<br />

<a id="align-sort-and-index-the-sample-datasets"></a>
#### Align, sort, and index the sample datasets
<a id="smartmaps-bowtie2-parameters"></a>
##### `SmartMap`'s `Bowtie2` parameters
<a id="code-3"></a>
###### Code
<details>
<summary><i>Code</i></summary>

```bash
#!/bin/bash

#NOTE
#  SmartMap's Bowtie2 parameters used for the alignment of paired-end reads:
#+ ❯ bowtie2 \
#+ >     -p "${pcores}" \
#+ >     -x "${btindex}" \
#+ >     --end-to-end \
#+ >     --very-fast \
#+ >     --no-discordant \
#+ >     --no-mixed \
#+ >     -k "${kfilt}" \
#+ >     -I "${I}" \
#+ >     -X "${L}" \
#+ >     -1 "${R1fq}" \
#+ >     -2 "${R2fq}" \
#+ >         | awk '$2==99 || $2==163 || $2==355 || $2==419{ as = ""; ys = ""; for (i=12; i<=NF; i++) { if ($i ~/^AS/){ as = $i } if ($i ~/^YS/){ ys = $i; break } } print $3"\t"$4"\t"$4+$9-1"\t"$1"\t"as"\t"ys }' \
#+ >         | awk '$2>0 && $3>$2 && $5~/^AS:i:/ && $6~/^YS:i:/{ print }' \
#+ >         | sed -e "s/$s//g" \
#+ >             > "${outname}"
#+
#+ >      --end-to-end \  # single-end or paired-end data
#+ >      --very-fast \  # single-end or paired-end data
#+ >      --no-discordant \  # paired-end data only
#+ >      --no-mixed \  # paired-end data only

#NOTE
#  Previous Bowtie2 call for alignment of single-end data:
#+ ❯ bowtie2 \
#+ >     -p "${threads}" \
#+ >     -x "${f_indices}" \
#+ >     --very-sensitive-local \
#+ >     -U "${i}"

#NOTE
#  Updated Bowtie2 call for alignment of single-end data:
#+ ❯ bowtie2 \
#+ >     -p "${threads}" \
#+ >     -x "${f_indices}" \
#+ >     --end-to-end \
#+ >     --very-fast \
#+ >     -U "${i}"
```
</details>
<br />

<a id="align-etc-atria-trimmed-fastqs"></a>
##### Align (etc.) `atria`-trimmed `fastq`s
<a id="code-4"></a>
###### Code
<details>
<summary><i>Code</i></summary>

```bash
#!/bin/bash

#  Run print tests to check that the commands are correct/reasonable
print_test=TRUE
[[ "${print_test}" == TRUE ]] &&
    {
        for i in "${atria[@]}"; do
            # i="${atria[0]}"  # ., "${i}"
            in_fastq="${i}"  # ., "${in_fastq}"
            out_bam="${d_bams}/$(basename "${in_fastq}" .fastq.gz).bam"  # echo "${out_bam}"

            echo """
            #  ---------------------------------------------------------
            #  Align fastqs trimmed with atria; sort resulting bams
            {
                bowtie2 \\
                    -p ${threads} \\
                    -x ${f_indices} \\
                    --end-to-end \\
                    --very-fast \\
                    -U ${i} \\
                        | samtools sort \\
                            -@ ${threads} \\
                            -T ${scratch} \\
                            -O bam \\
                            -o ${out_bam}
            } \\
                 > >(tee -a ${err_out}/$(basename "${out_bam}" .bam).bowtie2_samtools-sort.stdout.txt) \\
                2> >(tee -a ${err_out}/$(basename "${out_bam}" .bam).bowtie2_samtools-sort.stderr.txt)

            #  Index the sorted bams
            if [[ -f ${out_bam} ]]; then
                samtools index \\
                    -@ ${threads} \\
                    ${out_bam} \\
                         > >(tee -a ${err_out}/$(basename ${out_bam} .bam).samtools-index.stdout.txt) \\
                        2> >(tee -a ${err_out}/$(basename ${out_bam} .bam).samtools-index.stderr.txt)
            fi
            """
        done
    }

run=TRUE
[[ "${run}" == TRUE ]] &&
    {
        for i in "${atria[@]}"; do
            # i="${atria[0]}"  # echo "${i}"
            in_fastq="${i}"  # ., "${in_fastq}"
            out_bam="${d_bams}/$(basename "${in_fastq}" .fastq.gz).bam"  # echo "${out_bam}"

            
            if [[ ! -e "${out_bam}" ]]; then
                echo "#### Aligning, sorting, indexing $(basename ${in_fastq}) ####"
                #  Align fastqs trimmed with atria; sort resulting bams
                {
                    bowtie2 \
                        -p "${threads}" \
                        --end-to-end \
                        --very-fast \
                        -x "${f_indices}" \
                        -U "${in_fastq}" \
                            | samtools sort \
                                -@ "${threads}" \
                                -T "${scratch}" \
                                -O "bam" \
                                -o "${out_bam}"
                } \
                     > "${err_out}/$(basename "${out_bam}" .bam).bowtie2_samtools-sort.stdout.txt" \
                    2> "${err_out}/$(basename "${out_bam}" .bam).bowtie2_samtools-sort.stderr.txt"

                #  Index the sorted bams
                if [[ -f "${out_bam}" ]]; then
                    samtools index \
                        -@ "${threads}" \
                        "${out_bam}" \
                             > >(tee -a "${err_out}/$(basename "${out_bam}" .bam).samtools-index.stdout.txt") \
                            2> >(tee -a "${err_out}/$(basename "${out_bam}" .bam).samtools-index.stderr.txt")
                fi

                if [[ -f "${out_bam}" && -f "${out_bam}.bai" ]]; then
                    echo "#DONE"
                else
                    echo "#PROBLEM"
                fi
                echo ""
            else
                echo "Bams exist; skipping the running of Bowtie 2, Samtools, etc."
                break
            fi
        done
    }
```
</details>
<br />

<a id="align-etc-un-trimmed-fastqs"></a>
##### Align (etc.) un-trimmed `fastq`s
<a id="code-5"></a>
###### Code
<details>
<summary><i>Code</i></summary>

```bash
#!/bin/bash

#  Run print tests to check that the commands are correct/reasonable
print_test=TRUE
[[ "${print_test}" == TRUE ]] &&
    {
        for i in "${fastqs[@]}"; do
            # i="${fastqs[0]}"  # ., "${i}"
            in_fastq="${i}"  # ., "${in_fastq}"
            out_bam="${d_bams}/$(basename "${in_fastq}" .fastq.gz).bam"  # echo "${out_bam}"

            echo """
            #  ---------------------------------------------------------
            #  Align untrimmed fastqs; sort resulting bams
            {
                bowtie2 \\
                    -p ${threads} \\
                    -x ${f_indices} \\
                    --very-sensitive-local \\
                    -U ${i} \\
                        | samtools sort \\
                            -@ ${threads} \\
                            -T ${scratch} \\
                            -O bam \\
                            -o ${out_bam}
            } \\
                 > >(tee -a ${err_out}/$(basename "${out_bam}" .bam).bowtie2_samtools-sort.stdout.txt) \\
                2> >(tee -a ${err_out}/$(basename "${out_bam}" .bam).bowtie2_samtools-sort.stderr.txt)

            #  Index the sorted bams
            if [[ -f ${out_bam} ]]; then
                samtools index \\
                    -@ ${threads} \\
                    ${out_bam} \\
                         > >(tee -a ${err_out}/$(basename ${out_bam} .bam).samtools-index.stdout.txt) \\
                        2> >(tee -a ${err_out}/$(basename ${out_bam} .bam).samtools-index.stderr.txt)
            fi
            """
        done
    }

run=TRUE
[[ "${run}" == TRUE ]] &&
    {
        for i in "${fastqs[@]}"; do
            # i="${fastqs[0]}"  # echo "${i}"
            in_fastq="${i}"  # ., "${in_fastq}"
            out_bam="${d_bams}/$(basename "${in_fastq}" .fastq.gz).bam"  # echo "${out_bam}"

            if [[ ! -e "${out_bam}" ]]; then
                [[ -f "${err_out}/$(basename "${out_bam}" .bam).bowtie2_samtools-sort.stdout.txt" ]] && \
                    rm "${err_out}/$(basename "${out_bam}" .bam).bowtie2_samtools-sort.stdout.txt"
                
                [[ -f "${err_out}/$(basename "${out_bam}" .bam).bowtie2_samtools-sort.stderr.txt" ]] && \
                    rm "${err_out}/$(basename "${out_bam}" .bam).bowtie2_samtools-sort.stderr.txt"

                echo "#### Aligning, sorting, indexing $(basename ${in_fastq}) ####"
                #  Align untrimmed fastqs; sort resulting bams
                {
                    bowtie2 \
                        -p "${threads}" \
                        --very-sensitive-local \
                        -x "${f_indices}" \
                        -U "${in_fastq}" \
                            | samtools sort \
                                -@ "${threads}" \
                                -T "${scratch}" \
                                -O "bam" \
                                -o "${out_bam}"
                } \
                     > "${err_out}/$(basename "${out_bam}" .bam).bowtie2_samtools-sort.stdout.txt" \
                    2> "${err_out}/$(basename "${out_bam}" .bam).bowtie2_samtools-sort.stderr.txt"

                #  Index the sorted bams
                if [[ -f "${out_bam}" ]]; then
                    samtools index \
                        -@ "${threads}" \
                        "${out_bam}" \
                             > >(tee -a "${err_out}/$(basename "${out_bam}" .bam).samtools-index.stdout.txt") \
                            2> >(tee -a "${err_out}/$(basename "${out_bam}" .bam).samtools-index.stderr.txt")
                fi

                if [[ -f "${out_bam}" && -f "${out_bam}.bai" ]]; then
                    echo "#DONE"
                else
                    echo "#PROBLEM"
                fi
                echo ""
            else
                echo "Bams exist; skipping the running of Bowtie 2, Samtools, etc."
                break
            fi
        done
    }
```
</details>
<br />

<a id="examine-flags-in-bam-outfiles"></a>
### Examine flags in bam outfiles
<a id="initialize-necessary-functions"></a>
#### Initialize necessary functions
<a id="code-6"></a>
##### Code
<details>
<summary><i>Code</i></summary>

```bash
#!/bin/bash

check_dependency() {
    what="""
    check_dependency()
    ------------------
    Check if a program is available in the current environment
    
    :param 1: program to check <chr>
    
    #DONE Check that params are not empty
    #TODO Check that params are appropriate formats/strings
    #TODO Print to stderr
    """

    warning="""
    WARNING: param 1 is empty; stopping the function
    """

    [[ -z "${1}" ]] &&
        {
            echo "${warning}"
            echo "${what}"

            return 1
        }

    command -v "${1}" &>/dev/null ||
        {
            echo "Warning: ${1} not found; install or load ${1}; stopping the function"
            return 1
        }
}


calculate_run_time() {
    what="""
    calculate_run_time()
    --------------------
    Calculate run time for chunk of code
    
    :param 1: start time in \$(date +%s) format
    :param 2: end time in \$(date +%s) format
    :param 3: message to be displayed when printing the run time <chr>
    
    #DONE Check that params are not empty
    #TODO Check that params are appropriate formats/strings
    #TODO Print to stderr
    """

    warning="""
    WARNING: param(s) 1, 2, and/or 3 is/are empty; stopping the function
    """

    [[ -z "${1}" || -z "${2}" || -z "${3}" ]] &&
        {
            echo "${warning}"
            echo "${what}"

            return 1
        }

    run_time="$(echo "${2}" - "${1}" | bc -l)"
    
    echo ""
    echo "${3}"
    printf 'Run time: %dh:%dm:%ds\n' \
        $(( run_time/3600 )) \
        $(( run_time%3600/60 )) \
        $(( run_time%60 ))
    echo ""
}


display_spinning_icon() {
    what="""
    display_spinning_icon()
    -----------------------
    Display \"spinning icon\" while a background process runs
    
    :param 1: PID of the last program the shell ran in the background <pos int>
    :param 2: message to be displayed next to the spinning icon <chr>

    #DONE Check that params are not empty
    #TODO Check that params are appropriate formats/strings
    #TODO Print to stderr
    """

    warning="""
    WARNING: param(s) 1 and/or 2 is/are empty; stopping the function
    """

    [[ -z "${1}" || -z "${2}" ]] &&
        {
            echo "${warning}"
            echo "${what}"

            return 1
        }

    spin="/|\\–"
    i=0
    while kill -0 "${1}" 2> /dev/null; do
        i=$(( (i + 1) % 4 ))
        printf "\r${spin:$i:1} %s" "${2}"
        sleep .15
    done
}


list_tally_flags() {
    what="""
    list_tally_flags()
    ------------------
    List and tally flags in a bam infile; function acts on a bam infile to
    perform piped commands (samtools view, cut, sort, uniq -c, sort -nr) that
    list and tally flags; function writes the results to a txt outfile, the
    name of which is derived from the infile
    
    :param 1: name of bam infile, including path (chr)

    #DONE Check that param string is not empty
    #DONE Check that param file exists
    #DONE Check that dependency is present
    #TODO Check that params are appropriate formats/strings
    #TODO Print to stderr
    """

    warning_param="""
    WARNING: param 1 is empty; stopping the function
    """

    warning_file="""
    WARNING: param 1 file not found; stopping the function
    """

    warning_depend="""
    WARNING: One or more dependencies not found; stopping the function
    """

    [[ -z "${1}" ]] &&
        {
            echo "${warning_param}"
            echo "${what}"
            return 1
        }

    [[ -f "${1}" ]] ||
        {
            echo "${warning_file}"
            echo "${what}"
            return 1
        }

    check_dependency samtools \
        && check_dependency sort \
        && check_dependency uniq \
        && check_dependency display_spinning_icon \
        && check_dependency calculate_run_time
    [[ $? -gt 0 ]] &&
        {
            echo "${warning_depend}"
            return 1
        }

    start="$(date +%s)"
    
    samtools view "${1}" \
        | cut -f 2 \
        | sort \
        | uniq -c \
        | sort -nr \
            > "${1/.bam/.flags.txt}" &
    display_spinning_icon $! \
    "Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on $(basename "${1}")... "
        
    end="$(date +%s)"
    echo ""
    calculate_run_time "${start}" "${end}"  \
    "List and tally flags in $(basename "${1}")."
}


list_tally_MAPQs() {
    what="""
    list_tally_MAPQs()
    ------------------
    List and tally MAPQ scores in a bam infile; function acts on a bam infile
    to perform piped commands (samtools view, cut, sort, uniq -c, sort -nr)
    that list and tally MAPQ scores; function writes the results to a txt
    outfile, the name of which is derived from the infile
    
    :param 1: name of bam infile, including path (chr)

    #DONE Check that param string is not empty
    #DONE Check that param file exists
    #DONE Check that dependency is present
    #TODO Check that params are appropriate formats/strings
    #TODO Print to stderr
    """

    warning_param="""
    WARNING: param 1 is empty; stopping the function
    """

    warning_file="""
    WARNING: param 1 file not found; stopping the function
    """

    warning_depend="""
    WARNING: One or more dependencies not found; stopping the function
    """

    [[ -z "${1}" ]] &&
        {
            echo "${warning_param}"
            echo "${what}"
            return 1
        }

    [[ -f "${1}" ]] ||
        {
            echo "${warning_file}"
            echo "${what}"
            return 1
        }

    check_dependency samtools \
        && check_dependency sort \
        && check_dependency uniq \
        && check_dependency display_spinning_icon \
        && check_dependency calculate_run_time
    [[ $? -gt 0 ]] &&
        {
            echo "${warning_depend}"
            return 1
        }

    start="$(date +%s)"
    
    samtools view "${1}" \
        | cut -f 5 \
        | sort \
        | uniq -c \
        | sort -nr \
            > "${1/.bam/.MAPQs.txt}" &
    display_spinning_icon $! \
    "Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on $(basename "${1}")... "
        
    end="$(date +%s)"
    echo ""
    calculate_run_time "${start}" "${end}"  \
    "List and tally MAPQ scores in $(basename "${1}")."
}
```
</details>
<br />

<a id="initialize-an-array-of-bams"></a>
#### Initialize an array of bams
<a id="code-7"></a>
##### Code
<details>
<summary><i>Code</i></summary>

```bash
#!/bin/bash

p_data="${HOME}/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams"
unset bams
typeset -a bams=(
    "${p_data}/Ch_log_WT_Brn1_rep1.atria.bam"
    "${p_data}/Ch_log_WT_Brn1_rep1.bam"
    "${p_data}/Ch_log_WT_Brn1_rep2.atria.bam"
    "${p_data}/Ch_log_WT_Brn1_rep2.bam"
    "${p_data}/Ch_log_WT_Brn1_rep3.atria.bam"
    "${p_data}/Ch_log_WT_Brn1_rep3.bam"
    "${p_data}/Ch_log_WT_Brn1_repM.atria.bam"
    "${p_data}/Ch_log_WT_Brn1_repM.bam"
    "${p_data}/Ch_Q_WT_Brn1_rep1.atria.bam"
    "${p_data}/Ch_Q_WT_Brn1_rep1.bam"
    "${p_data}/Ch_Q_WT_Brn1_rep2.atria.bam"
    "${p_data}/Ch_Q_WT_Brn1_rep2.bam"
    "${p_data}/Ch_Q_WT_Brn1_rep3.atria.bam"
    "${p_data}/Ch_Q_WT_Brn1_rep3.bam"
    "${p_data}/Ch_Q_WT_Brn1_repM.atria.bam"
    "${p_data}/Ch_Q_WT_Brn1_repM.bam"
    "${p_data}/in_log_WT_Brn1_rep1.atria.bam"
    "${p_data}/in_log_WT_Brn1_rep1.bam"
    "${p_data}/in_log_WT_Brn1_rep2.atria.bam"
    "${p_data}/in_log_WT_Brn1_rep2.bam"
    "${p_data}/in_log_WT_Brn1_rep3.atria.bam"
    "${p_data}/in_log_WT_Brn1_rep3.bam"
    "${p_data}/in_log_WT_Brn1_repM.atria.bam"
    "${p_data}/in_log_WT_Brn1_repM.bam"
    "${p_data}/in_Q_WT_Brn1_rep1.atria.bam"
    "${p_data}/in_Q_WT_Brn1_rep1.bam"
    "${p_data}/in_Q_WT_Brn1_rep2.atria.bam"
    "${p_data}/in_Q_WT_Brn1_rep2.bam"
    "${p_data}/in_Q_WT_Brn1_rep3.atria.bam"
    "${p_data}/in_Q_WT_Brn1_rep3.bam"
    "${p_data}/in_Q_WT_Brn1_repM.atria.bam"
    "${p_data}/in_Q_WT_Brn1_repM.bam"
)

run_check=TRUE
[[ "${run_check}" == TRUE ]] &&
    {
        for i in "${bams[@]}"; do ls -lhaFG "${i}"; done
    }
```
</details>
<br />

<a id="check-on-flag-information-in-bams"></a>
#### Check on flag information in bams
<a id="code-8"></a>
##### Code
<details>
<summary><i>Code</i></summary>

```bash
#!/bin/bash

for i in "${bams[@]}"; do
    # i="${bams[0]}"  # echo "${i}"
    echo "#### $(basename ${i}) ####"
    list_tally_flags "${i}"
    echo ""
done

for i in "${bams[@]}"; do
    if [[ -f "${i/.bam/.flags.txt}" ]]; then
        echo "#### $(basename ${i}) ####"
        cat "${i/.bam/.flags.txt}"
        echo ""
    fi
done

for i in "${bams[@]}"; do
    # i="${bams[0]}"  # echo "${i}"
    echo "#### $(basename ${i}) ####"
    list_tally_MAPQs "${i}"
    echo ""
done

for i in "${bams[@]}"; do
    if [[ -f "${i/.bam/.MAPQs.txt}" ]]; then
        echo "#### $(basename ${i}) ####"
        cat "${i/.bam/.MAPQs.txt}"
        echo ""
    fi
done
```
</details>
<br />
<br />

<a id="tallycalculate-alignments"></a>
## Tally/calculate alignments
<a id="tallycalculate-alignments-1"></a>
### Tally/calculate alignments
<a id="initialize-functions-for-doing-floating-point-arithmetic-etc"></a>
#### Initialize functions for doing floating point arithmetic, etc.
<a id="code-9"></a>
##### Code
<details>
<summary><i>Code</i></summary>

```bash
#!/bin/bash

calc_6f() { awk "BEGIN{ printf \"%.6f\n\", $* }"; }


calc_2f() { awk "BEGIN{ printf \"%.2f\n\", $* }"; }


tally_alignments() {
    local OPTIND opt type MAPQ file Mito
    while getopts "t:q:f:m:" opt; do
        case "${opt}" in
            t) type="${OPTARG}" ;;
            q) MAPQ="${OPTARG}" ;;
            f) file="${OPTARG}" ;;
            m) Mito="${OPTARG}" ;;
            *) return 1 ;;
        esac
    done
    
    [[ -z "${type}" ]] &&
        {
            echo "Warning: Argument \"type\" is empty; stopping the function"
            return 1
        }    
    [[ -z "${MAPQ}" ]] && MAPQ=0
    [[ -z "${file}" ]] &&
        {
            echo "Warning: Argument \"file\" is empty; stopping the function"
            return 1
        }
    [[ -z "${Mito}" ]] && Mito=FALSE

    debug=FALSE
    [[ "${debug}" == TRUE ]] &&
        {
            type="all"
            MAPQ=0
            file="${i}"
            Mito=FALSE
    
            echo "${type}"
            echo "${MAPQ}"
            echo "${file}"
            echo "${Mito}"
        }

    #  All alignments
    [[ "${type}" == "all" ]] &&
        {
            if [[ "${Mito}" == TRUE || "${Mito}" == T ]]; then
                tally="$(
                    samtools view -c \
                        -F 4 \
                        -q "${MAPQ}" \
                        "${file}"
                )"
            elif [[ "${Mito}" == FALSE || "${Mito}" == F ]]; then
                tally="$(
                    samtools view -c \
                        -F 4 \
                        -q "${MAPQ}" \
                        "${file}" \
                        I II III IV V VI VII VIII IX X \
                        XI XII XIII XIV XV XVI \
                        SP_II_TG SP_I SP_II SP_III SP_MTR
                )"
            fi

            echo "${tally}" && return 0
        }


    #  S. cerevisiae alignments only
    [[ "${type}" == "SC" || "${type}" == "sc" || "${type}" == "Sc" ]] && 
        {
            if [[ "${Mito}" == TRUE ]]; then
                tally="$(
                    samtools view -c \
                        -F 4 \
                        -q "${MAPQ}" \
                        "${file}" \
                        I II III IV V VI VII VIII IX X \
                        XI XII XIII XIV XV XVI Mito
                )"
            elif [[ "${Mito}" == FALSE || "${Mito}" == F ]]; then
                tally="$(
                    samtools view -c \
                        -F 4 \
                        -q "${MAPQ}" \
                        "${file}" \
                        I II III IV V VI VII VIII IX X \
                        XI XII XIII XIV XV XVI
                )"
            fi

            echo "${tally}" && return 0
        }

    #  S. pombe alignments only
    [[ "${type}" == "SP" || "${type}" == "sp" || "${type}" == "Sp" ]] &&
        {   
            if [[ "${Mito}" == TRUE ]]; then
                tally="$(
                    samtools view -c \
                        -F 4 \
                        -q "${MAPQ}" \
                        "${file}" \
                        SP_II_TG SP_I SP_II SP_III SP_MTR SP_Mito
                )"
            elif [[ "${Mito}" == FALSE || "${Mito}" == F ]]; then
                tally="$(
                    samtools view -c \
                        -F 4 \
                        -q "${MAPQ}" \
                        "${file}" \
                        SP_II_TG SP_I SP_II SP_III SP_MTR
                )"
            fi

            echo "${tally}" && return 0
        }
}


transpose_values() {
    matrix="$(
        cat "${1}" \
            | awk '
                {
                    #  Loop through each line and each field (column) in the
                    #+ input file, storing the value of each field in the array
                    #+ "a" indexed by the column number "i" and the row number
                    #+ "NR"
                    for (i = 1; i <= NF; i++) {
                        a[i, NR] = $i
                    }
                }
                END {
                    #  Loop through the stored values in the array "a" and
                    #+ print the transposed values
                    for (i = 1; i <= NF; i++) {
                        for (j = 1; j <= NR; j++) {
                            printf "%s", a[i, j]
                            if (j < NR) {
                                printf " "
                            }
                        }
                        print ""
                    }
                }
            ' 
    )"

    echo "${matrix}" && return 0
}


export -f calc_6f
export -f calc_2f
export -f tally_alignments
export -f transpose_values
```
</details>
<br />

<a id="get-situated-then-initialize-arrays-variables-etc"></a>
#### Get situated, then initialize arrays, variables, etc.
<a id="code-10"></a>
##### Code
<details>
<summary><i>Code</i></summary>

```bash
#!/bin/bash

#  Load necessary programs
module purge
ml SAMtools/1.16.1-GCC-11.2.0

cd "${HOME}/tsukiyamalab/Kris/2023_rDNA/results" ||
    echo "cd'ing failed; check on this..."

p_data="${HOME}/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams"
unset bams && typeset -a bams=(
    "${p_data}/in_log_WT_Brn1_rep1.atria.bam"
    "${p_data}/Ch_log_WT_Brn1_rep1.atria.bam"
    "${p_data}/in_log_WT_Brn1_rep2.atria.bam"
    "${p_data}/Ch_log_WT_Brn1_rep2.atria.bam"
    "${p_data}/in_log_WT_Brn1_rep3.atria.bam"
    "${p_data}/Ch_log_WT_Brn1_rep3.atria.bam"
    "${p_data}/in_log_WT_Brn1_repM.atria.bam"
    "${p_data}/Ch_log_WT_Brn1_repM.atria.bam"
    "${p_data}/in_Q_WT_Brn1_rep1.atria.bam"
    "${p_data}/Ch_Q_WT_Brn1_rep1.atria.bam"
    "${p_data}/in_Q_WT_Brn1_rep2.atria.bam"
    "${p_data}/Ch_Q_WT_Brn1_rep2.atria.bam"
    "${p_data}/in_Q_WT_Brn1_rep3.atria.bam"
    "${p_data}/Ch_Q_WT_Brn1_rep3.atria.bam"
    "${p_data}/in_Q_WT_Brn1_repM.atria.bam"
    "${p_data}/Ch_Q_WT_Brn1_repM.atria.bam"
)

run_check=FALSE
[[ "${run_check}" == TRUE ]] && \
    for i in "${bams[@]}"; do ls -lhaFG "${i}"; done

#  Initialize other necessary variables, etc.
    d_proj="${HOME}/tsukiyamalab/Kris/2023_rDNA"
     d_exp="results/2023-0406_tutorial_ChIP-seq_analysis"
 d_tallies="${d_proj}/${d_exp}"         # ls -lhaFG "${d_tallies}"
  d_ratios="${d_proj}/${d_exp}"         # ., "${d_ratios}"

     f_tmp="tmp.txt"
 f_tallies="tallies.tsv"
  f_header="header.txt"
f_ind_ctrl="ind_ctrl.txt"
f_ind_test="ind_test.txt"
    f_comb="comb_ctrl_test.txt"
  f_ratios="ratios.tsv"
 f_metrics="metrics.tsv"

     a_tmp="${d_tallies}/${f_tmp}"      # echo "${a_tmp}"     # ., "${a_tmp}"     # cat "${a_tmp}"
 a_tallies="${d_tallies}/${f_tallies}"  # ., "${a_tallies}"   # cat "${a_tallies}"
  a_header="${d_ratios}/${f_header}"    # echo "${a_header}"  # ., "${a_header}"  # cat "${a_header}"
a_ind_ctrl="${d_ratios}/${f_ind_ctrl}"  # ., "${a_ind_ctrl}"
a_ind_test="${d_ratios}/${f_ind_test}"  # ., "${a_ind_test}"
    a_comb="${d_ratios}/${f_comb}"      # ., "${a_comb}"
  a_ratios="${d_ratios}/${f_ratios}"    # ., "${a_ratios}"
 a_metrics="${d_ratios}/${f_metrics}"   # ., "${a_metrics}"

    n_head=1                                                  # echo "${n_head}"
n_rows_all=$(cat "${a_tallies}" | awk 'END { print NR }')     # echo "${n_rows_all}"
    n_rows=$(( n_rows_all - n_head ))                         # echo "${n_rows}"

unset     rows && typeset -a rows=( $(seq 2 ${n_rows_all}) )  # echo_test "${rows[@]}"
unset     cols && typeset -a cols=( $(seq 35 42) )            # echo_test "${cols[@]}"
unset rat_ctrl && typeset -a rat_ctrl
unset rat_test && typeset -a rat_test

#  If outfiles already exist, then delete them
[[ -f "${a_tmp}" ]]      && rm "${a_tmp}"
[[ -f "${a_tallies}" ]]  && rm "${a_tallies}"
[[ -f "${a_ind_ctrl}" ]] && rm "${a_ind_ctrl}"
[[ -f "${a_ind_test}" ]] && rm "${a_ind_test}"
[[ -f "${a_comb}" ]]     && rm "${a_comb}"
[[ -f "${a_header}" ]]   && rm "${a_header}"
[[ -f "${a_ratios}" ]]   && rm "${a_ratios}"
[[ -f "${a_metrics}" ]]  && rm "${a_metrics}"
```
</details>
<br />

<a id="generate-tab-separated-table-of-alignment-talliescalculations"></a>
#### Generate tab-separated table of alignment tallies/calculations
<a id="code-11"></a>
##### Code
<details>
<summary><i>Code</i></summary>

```bash
#!/bin/bash

#  How long does it take to run this chunk of code? Record the start time
start="$(date +%s)"

#  Touch a temporary file ("${a_tmp}") and add a file header
touch "${a_tmp}"
echo "sample unmapped mapped total_with_Mito total_SC_with_Mito total_SP_with_Mito total_sans_Mito total_SC_sans_Mito total_SP_sans_Mito total_SC_Mito total_SP_Mito total_1 total_2 total_3 total_23 total_30 total_40 total_42 SC_0 SC_1 SC_2 SC_3 SC_23 SC_30 SC_40 SC_42 SP_0 SP_1 SP_2 SP_3 SP_23 SP_30 SP_40 SP_42 CC_SP2SC_0 CC_SP2SC_1 CC_SP2SC_2 CC_SP2SC_3 CC_SP2SC_23 CC_SP2SC_30 CC_SP2SC_40 CC_SP2SC_42" \
    >> "${a_tmp}"  # cat "${a_tmp}"

#  Get the tallies, perform calculations
unset tallies && typeset -A tallies
for (( i=0; i < ${#bams[@]}; i++ )); do
    # i=1  # echo_test "${bams[@]}"
    file="${bams[${i}]}"                 # echo "${file}"    # ls -lhaFG "${file}"
    sample="$(basename "${file}" .bam)"  # echo "${sample}"
    
    #  Strip string '.atria' if present
    [[ "${sample}" =~ ".atria" ]] \
        && sample="${sample%.atria}"     # echo "${sample}"
    
    #  Print current sample
    echo "#### $(basename ${file}) ####"


    #  Define "tally types" -------------------------------
              unmapped="${sample}.unmapped"            # echo "${unmapped}"
                mapped="${sample}.mapped"              # echo "${mapped}"

       total_with_Mito="${sample}.total_with_Mito"     # echo "${total_with_Mito}"
    total_SC_with_Mito="${sample}.total_SC_with_Mito"  # echo "${total_SC_with_Mito}"
    total_SP_with_Mito="${sample}.total_SP_with_Mito"  # echo "${total_SP_with_Mito}"

       total_sans_Mito="${sample}.total_sans_Mito"     # echo "${total_sans_Mito}"
    total_SC_sans_Mito="${sample}.total_SC_sans_Mito"  # echo "${total_SC_sans_Mito}"
    total_SP_sans_Mito="${sample}.total_SP_sans_Mito"  # echo "${total_SP_sans_Mito}"
         total_SC_Mito="${sample}.SC_Mito"             # echo "${total_SC_Mito}"
         total_SP_Mito="${sample}.SP_Mito"             # echo "${total_SP_Mito}"

               total_1="${sample}.total_1"             # echo "${total_1}"
               total_2="${sample}.total_2"             # echo "${total_2}"
               total_3="${sample}.total_3"             # echo "${total_3}"
              total_23="${sample}.total_23"            # echo "${total_23}"
              total_30="${sample}.total_30"            # echo "${total_30}"
              total_40="${sample}.total_40"            # echo "${total_40}"
              total_42="${sample}.total_42"            # echo "${total_42}"

                 SC_0="${sample}.SC_0"                 # echo "${SC_0}"
                 SP_0="${sample}.SP_0"                 # echo "${SP_0}"
                 SC_1="${sample}.SC_1"                 # echo "${SC_1}"
                 SP_1="${sample}.SP_1"                 # echo "${SP_1}"
                 SC_2="${sample}.SC_2"                 # echo "${SC_2}"
                 SP_2="${sample}.SP_2"                 # echo "${SP_2}"
                 SC_3="${sample}.SC_3"                 # echo "${SC_3}"
                 SP_3="${sample}.SP_3"                 # echo "${SP_3}"
                SC_23="${sample}.SC_23"                # echo "${SC_23}"
                SP_23="${sample}.SP_23"                # echo "${SP_23}"
                SC_30="${sample}.SC_30"                # echo "${SC_30}"
                SP_30="${sample}.SP_30"                # echo "${SP_30}"
                SC_40="${sample}.SC_40"                # echo "${SC_40}"
                SP_40="${sample}.SP_40"                # echo "${SP_40}"
                SC_42="${sample}.SC_42"                # echo "${SC_42}"
                SP_42="${sample}.SP_42"                # echo "${SP_42}"

           CC_SP2SC_0="${sample}.CC_SP2SC_0"           # echo "${CC_SP2SC_0}"
           CC_SP2SC_1="${sample}.CC_SP2SC_1"           # echo "${CC_SP2SC_1}"
           CC_SP2SC_2="${sample}.CC_SP2SC_2"           # echo "${CC_SP2SC_2}"
           CC_SP2SC_3="${sample}.CC_SP2SC_3"           # echo "${CC_SP2SC_3}"
          CC_SP2SC_23="${sample}.CC_SP2SC_23"          # echo "${CC_SP2SC_23}"
          CC_SP2SC_30="${sample}.CC_SP2SC_30"          # echo "${CC_SP2SC_30}"
          CC_SP2SC_40="${sample}.CC_SP2SC_40"          # echo "${CC_SP2SC_40}"
          CC_SP2SC_42="${sample}.CC_SP2SC_42"          # echo "${CC_SP2SC_42}"

        #  CC_SC2SP_0="${sample}.CC_SC2SP_0"           # echo "${CC_SC2SP_0}"
        #  CC_SC2SP_1="${sample}.CC_SC2SP_1"           # echo "${CC_SC2SP_1}"
        #  CC_SC2SP_2="${sample}.CC_SC2SP_2"           # echo "${CC_SC2SP_2}"
        #  CC_SC2SP_3="${sample}.CC_SC2SP_3"           # echo "${CC_SC2SP_3}"
        # CC_SC2SP_23="${sample}.CC_SC2SP_23"          # echo "${CC_SC2SP_23}"
        # CC_SC2SP_30="${sample}.CC_SC2SP_30"          # echo "${CC_SC2SP_30}"
        # CC_SC2SP_40="${sample}.CC_SC2SP_40"          # echo "${CC_SC2SP_40}"
        # CC_SC2SP_42="${sample}.CC_SC2SP_42"          # echo "${CC_SC2SP_42}"


    #  Run print test -------------------------------------
    print_test=FALSE
    [[ "${print_test}" == TRUE ]] &&
        {
            echo """
            ${sample}
            ----------------------------------------
            ${unmapped}
            ${mapped}

            ${total_with_Mito}
            ${total_SC_with_Mito}
            ${total_SP_with_Mito}

            ${total_sans_Mito}
            ${total_SC_sans_Mito}
            ${total_SP_sans_Mito}
            ${total_SC_Mito}
            ${total_SP_Mito}

            ${total_1}
            ${total_2}
            ${total_3}
            ${total_23}
            ${total_30}
            ${total_40}
            ${total_42}

            ${SC_0}
            ${SP_0}
            ${SC_1}
            ${SP_1}
            ${SC_2}
            ${SP_2}
            ${SC_3}
            ${SP_3}
            ${SC_23}
            ${SP_23}
            ${SC_30}
            ${SP_30}
            ${SC_40}
            ${SP_40}
            ${SC_42}
            ${SC_42}

            ${CC_SP2SC_0}
            ${CC_SP2SC_1}
            ${CC_SP2SC_2}
            ${CC_SP2SC_3}
            ${CC_SP2SC_23}
            ${CC_SP2SC_30}
            ${CC_SP2SC_40}
            ${CC_SP2SC_42}
            """
        }


    #  Get basic tallies ----------------------------------
    #  Get total counts (including Mito and unmapped counts)
    tallies["${total_with_Mito}"]="$(
        samtools view -c "${file}"
    )"  # echo "${tallies[${total_with_Mito}]}"

    #  Get total unaligned counts
    tallies["${unmapped}"]="$(
        samtools view -c -f 4 "${file}"
    )"  # echo "${tallies[${unmapped}]}"

    #  Get total aligned counts
    tallies["${mapped}"]="$(
        tally_alignments -t "all" -q 0 -f "${file}" -m TRUE
    )"  # echo "${tallies[${mapped}]}"

    #+ Get total counts to S. cerevisiae (including Mito)
    tallies["${total_SC_with_Mito}"]="$(
        tally_alignments -t "SC" -q 0 -f "${file}" -m TRUE
    )"  # echo "${tallies[${total_SC_with_Mito}]}"

    #+ Get total counts to S. pombe (including Mito)
    tallies["${total_SP_with_Mito}"]="$(
        tally_alignments -t "SP" -q 0 -f "${file}" -m TRUE
    )"  # echo "${tallies[${total_SP_with_Mito}]}"


    #  Get main tallies/tallies for QC ------------------------
    #  Get total aligned counts sans Mito
    tallies["${total_sans_Mito}"]="$(
        tally_alignments -t "all" -q 0 -f "${file}" -m FALSE
    )"  # echo "${tallies[${total_sans_Mito}]}"

    #  Get total S. cerevisiae Mito counts
    tallies["${total_SC_sans_Mito}"]="$(
        tally_alignments -t "SC" -q 0 -f "${file}" -m FALSE
    )"  # echo "${tallies[${total_SC_sans_Mito}]}"

    #  Get total S. pombe Mito counts
    tallies["${total_SP_sans_Mito}"]="$(
        tally_alignments -t "SP" -q 0 -f "${file}" -m FALSE
    )"  # echo "${tallies[${total_SP_sans_Mito}]}"

    #  Get total S. cerevisiae Mito counts
    tallies["${total_SC_Mito}"]="$(
        samtools view -c "${file}" Mito
    )"  # echo "${tallies[${total_SC_Mito}]}"

    #  Get total S. pombe Mito counts
    tallies["${total_SP_Mito}"]="$(
        samtools view -c "${file}" SP_Mito
    )"  # echo "${tallies[${total_SP_Mito}]}"

    #  Get total aligned counts MAPQ >= 1
    tallies["${total_1}"]="$(
        tally_alignments -t "all" -q 1 -f "${file}" -m FALSE
    )"  # echo "${tallies[${total_1}]}"

    #  Get total aligned counts MAPQ >= 2
    tallies["${total_2}"]="$(
        tally_alignments -t "all" -q 2 -f "${file}" -m FALSE
    )"  # echo "${tallies[${total_2}]}"

    #  Get total aligned counts MAPQ >= 3
    tallies["${total_3}"]="$(
        tally_alignments -t "all" -q 3 -f "${file}" -m FALSE
    )"  # echo "${tallies[${total_3}]}"

    #  Get total aligned counts MAPQ >= 23
    tallies["${total_23}"]="$(
        tally_alignments -t "all" -q 23 -f "${file}" -m FALSE
    )"  # echo "${tallies[${total_23}]}"

    #  Get total aligned counts MAPQ >= 30
    tallies["${total_30}"]="$(
        tally_alignments -t "all" -q 30 -f "${file}" -m FALSE
    )"  # echo "${tallies[${total_30}]}"

    #  Get total aligned counts MAPQ >= 40
    tallies["${total_40}"]="$(
        tally_alignments -t "all" -q 40 -f "${file}" -m FALSE
    )"  # echo "${tallies[${total_40}]}"

    #  Get total aligned counts MAPQ >= 42
    tallies["${total_42}"]="$(
        tally_alignments -t "all" -q 42 -f "${file}" -m FALSE
    )"  # echo "${tallies[${total_42}]}"


    #  Get SC, SP tallies to calculate scaling factors ----
    #  Get total counts aligned MAPQ >= 0 -------
    #+ ...to S. cerevisiae
    tallies["${SC_0}"]="$(
        tally_alignments -t "SC" -q 0 -f "${file}" -m FALSE
    )"  # echo "${tallies[${SC_0}]}"

    #+ ...to S. pombe
    tallies["${SP_0}"]="$(
        tally_alignments -t "SP" -q 0 -f "${file}" -m FALSE
    )"  # echo "${tallies[${SP_0}]}"
    
    #  Get total counts aligned MAPQ >= 1 -------
    #+ ...to S. cerevisiae
    tallies["${SC_1}"]="$(
        tally_alignments -t "SC" -q 1 -f "${file}" -m FALSE
    )"  # echo "${tallies[${SC_1}]}"

    #+ ...to S. pombe
    tallies["${SP_1}"]="$(
        tally_alignments -t "SP" -q 1 -f "${file}" -m FALSE
    )"  # echo "${tallies[${SP_1}]}"

    #  Get total counts aligned MAPQ >= 2 -------
    #+ ...to S. cerevisiae
    tallies["${SC_2}"]="$(
        tally_alignments -t "SC" -q 2 -f "${file}" -m FALSE
    )"  # echo "${tallies[${SC_2}]}"

    #+ ...to S. pombe
    tallies["${SP_2}"]="$(
        tally_alignments -t "SP" -q 2 -f "${file}" -m FALSE
    )"  # echo "${tallies[${SP_2}]}"

    #  Get total counts aligned MAPQ >= 3 -------
    #+ ...to S. cerevisiae
    tallies["${SC_3}"]="$(
        tally_alignments -t "SC" -q 3 -f "${file}" -m FALSE
    )"  # echo "${tallies[${SC_3}]}"

    #+ ...to S. pombe
    tallies["${SP_3}"]="$(
        tally_alignments -t "SP" -q 3 -f "${file}" -m FALSE
    )"  # echo "${tallies[${SP_3}]}"

    #  Get total counts aligned MAPQ >= 23 ------
    #+ ...to S. cerevisiae
    tallies["${SC_23}"]="$(
        tally_alignments -t "SC" -q 23 -f "${file}" -m FALSE
    )"  # echo "${tallies[${SC_23}]}"

    #+ ...to S. pombe
    tallies["${SP_23}"]="$(
        tally_alignments -t "SP" -q 23 -f "${file}" -m FALSE
    )"  # echo "${tallies[${SP_23}]}"

    #  Get total counts aligned MAPQ >= 30 ------
    #+ ...to S. cerevisiae
    tallies["${SC_30}"]="$(
        tally_alignments -t "SC" -q 30 -f "${file}" -m FALSE
    )"  # echo "${tallies[${SC_30}]}"

    #+ ...to S. pombe
    tallies["${SP_30}"]="$(
        tally_alignments -t "SP" -q 30 -f "${file}" -m FALSE
    )"  # echo "${tallies[${SP_30}]}"

    #  Get total counts aligned MAPQ >= 40 ------
    #+ ...to S. cerevisiae
    tallies["${SC_40}"]="$(
        tally_alignments -t "SC" -q 40 -f "${file}" -m FALSE
    )"  # echo "${tallies[${SC_40}]}"

    #+ ...to S. pombe
    tallies["${SP_40}"]="$(
        tally_alignments -t "SP" -q 40 -f "${file}" -m FALSE
    )"  # echo "${tallies[${SP_40}]}"

    #  Get total counts aligned MAPQ >= 42 ------
    #+ ...to S. cerevisiae
    tallies["${SC_42}"]="$(
        tally_alignments -t "SC" -q 42 -f "${file}" -m FALSE
    )"  # echo "${tallies[${SC_42}]}"
    
    #+ ...to S. pombe
    tallies["${SP_42}"]="$(
        tally_alignments -t "SP" -q 42 -f "${file}" -m FALSE
    )"  # echo "${tallies[${SP_42}]}"


    #  Calculate SP:SC ratios -----------------------------
    #+ SP as numerator, MAPQ >= 0
    tallies["${CC_SP2SC_0}"]="$(
        calc_6f "${tallies["${SP_0}"]}"/"${tallies["${SC_0}"]}"
    )"  # echo "${tallies[${CC_SP2SC_0}]}"
    
    #+ SP as numerator, MAPQ >= 1
    tallies["${CC_SP2SC_1}"]="$(
        calc_6f "${tallies["${SP_1}"]}"/"${tallies["${SC_1}"]}"
    )"  # echo "${tallies[${CC_SP2SC_1}]}"
    
    #+ SP as numerator, MAPQ >= 2
    tallies["${CC_SP2SC_2}"]="$(
        calc_6f "${tallies["${SP_2}"]}"/"${tallies["${SC_2}"]}"
    )"  # echo "${tallies[${CC_SP2SC_2}]}"
    
    #+ SP as numerator, MAPQ >= 3
    tallies["${CC_SP2SC_3}"]="$(
        calc_6f "${tallies["${SP_3}"]}"/"${tallies["${SC_3}"]}"
    )"  # echo "${tallies[${CC_SP2SC_3}]}"
    
    #+ SP as numerator, MAPQ >= 23
    tallies["${CC_SP2SC_23}"]="$(
        calc_6f "${tallies["${SP_23}"]}"/"${tallies["${SC_23}"]}"
    )"  # echo "${tallies[${CC_SP2SC_23}]}"
    
    #+ SP as numerator, MAPQ >= 30
    tallies["${CC_SP2SC_30}"]="$(
        calc_6f "${tallies["${SP_30}"]}"/"${tallies["${SC_30}"]}"
    )"  # echo "${tallies[${CC_SP2SC_30}]}"
    
    #+ SP as numerator, MAPQ >= 40
    tallies["${CC_SP2SC_40}"]="$(
        calc_6f "${tallies["${SP_40}"]}"/"${tallies["${SC_40}"]}"
    )"  # echo "${tallies[${CC_SP2SC_40}]}"

    #+ SP as numerator, MAPQ >= 42
    tallies["${CC_SP2SC_42}"]="$(
        calc_6f "${tallies["${SP_42}"]}"/"${tallies["${SC_42}"]}"
    )"  # echo "${tallies[${CC_SP2SC_42}]}"
    
    # #+ SC as numerator, MAPQ >= 0
    # tallies["${CC_SC2SP_0}"]="$(
    #     calc_2f "${tallies["${SC_0}"]}"/"${tallies["${SP_0}"]}"
    # )"
    #
    # #+ SC as numerator, MAPQ >= 1
    # tallies["${CC_SC2SP_1}"]="$(
    #     calc_2f "${tallies["${SC_1}"]}"/"${tallies["${SP_1}"]}"
    # )"
    #
    # #+ SC as numerator, MAPQ >= 2
    # tallies["${CC_SC2SP_2}"]="$(
    #     calc_2f "${tallies["${SC_2}"]}"/"${tallies["${SP_2}"]}"
    # )"
    #
    # #+ SC as numerator, MAPQ >= 3
    # tallies["${CC_SC2SP_3}"]="$(
    #     calc_2f "${tallies["${SC_3}"]}"/"${tallies["${SP_3}"]}"
    # )"
    #
    # #+ SC as numerator, MAPQ >= 23
    # tallies["${CC_SC2SP_23}"]="$(
    #     calc_2f "${tallies["${SC_23}"]}"/"${tallies["${SP_23}"]}"
    # )"
    #
    # #+ SC as numerator, MAPQ >= 30
    # tallies["${CC_SC2SP_30}"]="$(
    #     calc_2f "${tallies["${SC_30}"]}"/"${tallies["${SP_30}"]}"
    # )"
    #
    # #+ SC as numerator, MAPQ >= 40
    # tallies["${CC_SC2SP_40}"]="$(
    #     calc_2f "${tallies["${SC_40}"]}"/"${tallies["${SP_40}"]}"
    # )"
    #
    # #+ SC as numerator, MAPQ >= 42
    # tallies["${CC_SC2SP_42}"]="$(
    #     calc_2f "${tallies["${SC_42}"]}"/"${tallies["${SP_42}"]}"
    # )"


    #  Run print test -------------------------------------
    print_test=TRUE
    [[ "${print_test}" == TRUE ]] &&
        {
            echo """
            ${sample}
            ------------------------------------------------------------
                      \"\${tallies[\${unmapped}]}\"  ${tallies[${unmapped}]}
                        \"\${tallies[\${mapped}]}\"  ${tallies[${mapped}]}

               \"\${tallies[\${total_with_Mito}]}\"  ${tallies[${total_with_Mito}]}
            \"\${tallies[\${total_SC_with_Mito}]}\"  ${tallies[${total_SC_with_Mito}]}
            \"\${tallies[\${total_SP_with_Mito}]}\"  ${tallies[${total_SP_with_Mito}]}

               \"\${tallies[\${total_sans_Mito}]}\"  ${tallies[${total_sans_Mito}]}
            \"\${tallies[\${total_SC_sans_Mito}]}\"  ${tallies[${total_SC_sans_Mito}]}
            \"\${tallies[\${total_SP_sans_Mito}]}\"  ${tallies[${total_SP_sans_Mito}]}

                 \"\${tallies[\${total_SC_Mito}]}\"  ${tallies[${total_SC_Mito}]}
                 \"\${tallies[\${total_SP_Mito}]}\"  ${tallies[${total_SP_Mito}]}

                       \"\${tallies[\${total_1}]}\"  ${tallies[${total_1}]}
                       \"\${tallies[\${total_2}]}\"  ${tallies[${total_2}]}
                       \"\${tallies[\${total_3}]}\"  ${tallies[${total_3}]}
                      \"\${tallies[\${total_23}]}\"  ${tallies[${total_23}]}
                      \"\${tallies[\${total_30}]}\"  ${tallies[${total_30}]}
                      \"\${tallies[\${total_40}]}\"  ${tallies[${total_40}]}
                      \"\${tallies[\${total_42}]}\"  ${tallies[${total_42}]}

                          \"\${tallies[\${SC_0}]}\"  ${tallies[${SC_0}]}
                          \"\${tallies[\${SC_1}]}\"  ${tallies[${SC_1}]}
                          \"\${tallies[\${SC_2}]}\"  ${tallies[${SC_2}]}
                          \"\${tallies[\${SC_3}]}\"  ${tallies[${SC_3}]}
                         \"\${tallies[\${SC_23}]}\"  ${tallies[${SC_23}]}
                         \"\${tallies[\${SC_30}]}\"  ${tallies[${SC_30}]}
                         \"\${tallies[\${SC_40}]}\"  ${tallies[${SC_40}]}
                         \"\${tallies[\${SC_42}]}\"  ${tallies[${SC_42}]}

                          \"\${tallies[\${SP_0}]}\"  ${tallies[${SP_0}]}
                          \"\${tallies[\${SP_1}]}\"  ${tallies[${SP_1}]}
                          \"\${tallies[\${SP_2}]}\"  ${tallies[${SP_2}]}
                          \"\${tallies[\${SP_3}]}\"  ${tallies[${SP_3}]}
                         \"\${tallies[\${SP_23}]}\"  ${tallies[${SP_23}]}
                         \"\${tallies[\${SP_30}]}\"  ${tallies[${SP_30}]}
                         \"\${tallies[\${SP_40}]}\"  ${tallies[${SP_40}]}
                         \"\${tallies[\${SP_42}]}\"  ${tallies[${SP_42}]}

                    \"\${tallies[\${CC_SP2SC_0}]}\"  ${tallies[${CC_SP2SC_0}]}
                    \"\${tallies[\${CC_SP2SC_1}]}\"  ${tallies[${CC_SP2SC_1}]}
                    \"\${tallies[\${CC_SP2SC_2}]}\"  ${tallies[${CC_SP2SC_2}]}
                    \"\${tallies[\${CC_SP2SC_3}]}\"  ${tallies[${CC_SP2SC_3}]}
                   \"\${tallies[\${CC_SP2SC_23}]}\"  ${tallies[${CC_SP2SC_23}]}
                   \"\${tallies[\${CC_SP2SC_30}]}\"  ${tallies[${CC_SP2SC_30}]}
                   \"\${tallies[\${CC_SP2SC_40}]}\"  ${tallies[${CC_SP2SC_40}]}
                   \"\${tallies[\${CC_SP2SC_42}]}\"  ${tallies[${CC_SP2SC_42}]}
            """
        }

    #  Run print test
    print_test=FALSE
    [[ "${print_test}" ==  TRUE ]] &&
        {
            echo "${sample} ${tallies[${unmapped}]} ${tallies[${mapped}]} ${tallies[${total_with_Mito}]} ${tallies[${total_SC_with_Mito}]} ${tallies[${total_SP_with_Mito}]} ${tallies[${total_sans_Mito}]} ${tallies[${total_SC_sans_Mito}]} ${tallies[${total_SP_sans_Mito}]} ${tallies[${total_SC_Mito}]} ${tallies[${total_SP_Mito}]} ${tallies[${total_1}]} ${tallies[${total_2}]} ${tallies[${total_3}]} ${tallies[${total_23}]} ${tallies[${total_30}]} ${tallies[${total_40}]} ${tallies[${total_42}]} ${tallies[${SC_0}]} ${tallies[${SC_1}]} ${tallies[${SC_2}]} ${tallies[${SC_3}]} ${tallies[${SC_23}]} ${tallies[${SC_30}]} ${tallies[${SC_40}]} ${tallies[${SC_42}]} ${tallies[${SP_0}]} ${tallies[${SP_1}]} ${tallies[${SP_2}]} ${tallies[${SP_3}]} ${tallies[${SP_23}]} ${tallies[${SP_30}]} ${tallies[${SP_40}]} ${tallies[${SP_42}]} ${tallies[${CC_SP2SC_0}]} ${tallies[${CC_SP2SC_1}]} ${tallies[${CC_SP2SC_2}]} ${tallies[${CC_SP2SC_3}]} ${tallies[${CC_SP2SC_23}]} ${tallies[${CC_SP2SC_30}]} ${tallies[${CC_SP2SC_40}]} ${tallies[${CC_SP2SC_42}]}"
            echo ""
        }

    #  Write the tallies (etc.) to the temporary file (line by line with each
    #+ iteration of the loop)
    echo "${sample} ${tallies[${unmapped}]} ${tallies[${mapped}]} ${tallies[${total_with_Mito}]} ${tallies[${total_SC_with_Mito}]} ${tallies[${total_SP_with_Mito}]} ${tallies[${total_sans_Mito}]} ${tallies[${total_SC_sans_Mito}]} ${tallies[${total_SP_sans_Mito}]} ${tallies[${total_SC_Mito}]} ${tallies[${total_SP_Mito}]} ${tallies[${total_1}]} ${tallies[${total_2}]} ${tallies[${total_3}]} ${tallies[${total_23}]} ${tallies[${total_30}]} ${tallies[${total_40}]} ${tallies[${total_42}]} ${tallies[${SC_0}]} ${tallies[${SC_1}]} ${tallies[${SC_2}]} ${tallies[${SC_3}]} ${tallies[${SC_23}]} ${tallies[${SC_30}]} ${tallies[${SC_40}]} ${tallies[${SC_42}]} ${tallies[${SP_0}]} ${tallies[${SP_1}]} ${tallies[${SP_2}]} ${tallies[${SP_3}]} ${tallies[${SP_23}]} ${tallies[${SP_30}]} ${tallies[${SP_40}]} ${tallies[${SP_42}]} ${tallies[${CC_SP2SC_0}]} ${tallies[${CC_SP2SC_1}]} ${tallies[${CC_SP2SC_2}]} ${tallies[${CC_SP2SC_3}]} ${tallies[${CC_SP2SC_23}]} ${tallies[${CC_SP2SC_30}]} ${tallies[${CC_SP2SC_40}]} ${tallies[${CC_SP2SC_42}]}" \
        >> "${a_tmp}"
done

#  Convert spaces to tabs
cat "${a_tmp}" | sed "s/\ /\t/g" > "${a_tallies}"
cat "${a_tallies}"

#  If "${a_tallies}" was written, then rm "${a_tmp}"
[[ -f "${a_tallies}" ]] && rm "${a_tmp}"

#  How long does it take to run this chunk of code? Record the end time, then
#+ calculate and print the total time
end="$(date +%s)"
calculate_run_time "${start}" "${end}"  \
"Generate tab-separated table of alignment tallies/calculations."
```
</details>
<br />

<a id="calculate-ccss-styled-scaling-factors"></a>
#### Calculate CC/SS-styled scaling factors
<a id="code-12"></a>
##### Code
<details>
<summary><i>Code</i></summary>

```bash
#!/bin/bash

#  How long does it take to run this chunk of code? Record the start time
start="$(date +%s)"

#  Calculate the scaling factors (ratios), writing them to arrays for control
#+ and test samples
for j in "${cols[@]}"; do
    for i in "${rows[@]}"; do
        #  If even row index, calculate control scaling factor
        if (( i % 2 == 0 )); then
            #  Using awk, extract the value at row i and column j from the
            #+ tallies file ("${a_tallies}"); the field separator (-F '\t') is
            #+ set to tab, and the variables i and j are passed to awk using
            #+ the -v flags; the 'NR == i { print $j }' condition prints the
            #+ value at the specified row and column
            ctrl="$(
                cat "${a_tallies}" \
                    | awk \
                        -F '\t' \
                        -v i="${i}" \
                        -v j="${j}" \
                        'NR == i { print $j }'
            )"

            #  Calculate scaling factor by dividing ctrl by ctrl (i.e., scaling
            #+ factor is 1)
            sf=$(echo "scale=0; ${ctrl} / ${ctrl}" | bc)
            [[ "${sf}" == "" ]] && sf="ND"

            #  Store scaling factor in rat_ctrl array
            rat_ctrl+=( "${sf}" )
        else
            #  Otherwise, if odd row index, then calculate test scaling factor
            #+ by incorporating information from odd row (test) and preceding
            #+ even row (ctrl)
            i_1=$(( i - 1 ))
            
            #  As above except row is i - 1
            ctrl="$(
                cat "${a_tallies}" \
                    | awk \
                        -F '\t' \
                        -v i="${i_1}" \
                        -v j="${j}" \
                        'NR == i { print $j }'
            )"

            #  As above
            test="$(
                cat "${a_tallies}" \
                    | awk \
                        -F '\t' \
                        -v i="${i}" \
                        -v j="${j}" \
                        'NR == i { print $j }'
            )"

            #  Calculate scaling factor by dividing ctrl by test
            sf=$(echo "scale=6; ${ctrl} / ${test}" | bc)
            [[ "${sf}" == "" ]] && sf="ND"

            #  Store scaling factor in rat_test array
            rat_test+=( "${sf}" )
        fi
    done
done

#  Run a print test for the arrays of ratios
print_test=FALSE
[[ ${print_test} == TRUE ]] &&
    {
        echo_test "${rat_ctrl[@]}" && echo ""
        echo_test "${rat_test[@]}"
    }

#  Print ctrl array elements, with 8 elements per column
[[ ! -f "${a_ind_ctrl}" ]] && touch "${a_ind_ctrl}"
for (( i = 0; i < ${#rat_ctrl[@]}; i++ )); do
    echo -n "${rat_ctrl[$i]} " >> "${a_ind_ctrl}"
    if (( (i + 1) % 8 == 0 )); then
        echo >> "${a_ind_ctrl}"
    fi
done
# cat "${a_ind_ctrl}"  # echo "${a_ind_ctrl}"

#  Transpose control values
transpose_values "${a_ind_ctrl}" > "${d_ratios}/tmp"
if [[ $? -eq 0 && -f "${d_ratios}/tmp" ]]; then
    mv -f "${d_ratios}/tmp" "${a_ind_ctrl}"
fi
# cat "${a_ind_ctrl}"  # echo "${a_ind_ctrl}"

#  Print test array elements, with 8 elements per column
[[ ! -f "${a_ind_test}" ]] && touch "${a_ind_test}"
for (( i = 0; i < ${#rat_test[@]}; i++ )); do
    echo -n "${rat_test[$i]} " >> "${a_ind_test}"
    if (( (i + 1) % 8 == 0 )); then
        echo >> "${a_ind_test}"
    fi
done
# cat "${a_ind_test}"  # echo "${a_ind_test}"

#  Transpose test values
transpose_values "${a_ind_test}" > "${d_ratios}/tmp"
if [[ $? -eq 0 && -f "${d_ratios}/tmp" ]]; then
    mv -f "${d_ratios}/tmp" "${a_ind_test}"
fi
# cat "${a_ind_test}"  # echo "${a_ind_test}"

#  Interleave the ctrl and test matrices
[[ -f "${comb_ctrl_test}" ]] && rm "${comb_ctrl_test}"
paste -d '\n' "${a_ind_ctrl}" "${a_ind_test}" > "${a_comb}"
# cat "${a_comb}"

#  Touch a temporary file with a header
touch "${a_header}"
echo "sf_0 sf_1 sf_2 sf_3 sf_23 sf_30 sf_40 sf_42" >> "${a_header}"
# cat "${a_header}"

#  Concatenate the header and interleaved-matrix files
cat "${a_header}" "${a_comb}" >> "${d_ratios}/tmp"
[[ $? -eq 0 && -f "${d_ratios}/tmp" ]] \
    && mv -f "${d_ratios}/tmp" "${a_comb}"
# cat "${a_header}"  # echo "${a_comb}"

#  Remove temporary files
if [[ $? -eq 0 && -f "${a_comb}" && ! -f "${d_ratios}/tmp" ]]; then
    rm "${a_ind_ctrl}" "${a_ind_test}" "${a_header}"
fi

#  Convert spaces to tabs, then column-bind the table of tallies with the table
#+ of ratios
cat "${a_comb}" | sed "s/\ /\t/g" > "${a_ratios}"
paste "${a_tallies}" "${a_ratios}" > "${a_metrics}"

#  If "${a_metrics}" was written, then rm other temporary files
[[ $? -eq 0 && -f "${a_metrics}" ]] \
    && rm "${a_comb}" "${a_ratios}"

#  Run print test
print_test=TRUE
[[ "${print_test}" == TRUE ]] &&
    {
        echo "----------------------------------------"
        cat "${a_tallies}"
        echo ""

        echo "----------------------------------------"
        cat "${a_metrics}"
    }

#  How long does it take to run this chunk of code? Record the end time, then
#+ calculate and print the total time
end="$(date +%s)"
calculate_run_time "${start}" "${end}"  \
"Calculate CC/SS-styled scaling factors."
```
</details>
<br />
<br />

<a id="on-calculating-scaling-factors"></a>
### On calculating scaling factors
<a id="email-from-christine-edited-by-me"></a>
#### Email from Christine (edited by me)
<details>
<summary><i>Email from Christine (edited by me)</i></summary>

Title: ChIP-seq analyses: Friendly reminder to check if, when using `deepTools`, you input the scaling factors as reciprocals  
From: Cucinotta, Christine E  
To: Alavattam, Kris  
Time: Tue 7/18/2023 1:38 PM

Hi Kris,

I take the `pombe:cerevisiae` ratio for all the samples. Then, I take the ratio of the sample we want to normalize to (e.g., WT or a time zero) and divide by the test sample.

I have trouble writing it out, so here’s an example:
```txt
 Q pombe:cerevisiae = 0.35  # Control ratio
5m pombe:cerevisiae = 0.26  # Test ratio

 Q scaling factor = 0.35/0.35 = 1     # Use this for the scaling factor in deepTools (control/control)  
5m Scaling factor = 0.35/0.26 = 1.35  # Use this number for the scaling factor in deepTools (control/test)
```
I hope this makes sense. Please let me know if not&mdash;happy to chat!

Thanks!

-Chrsitine
</details>
<br />

<a id="notes-on-using-excel-to-calculate-scaling-factors"></a>
#### Notes on using Excel to calculate scaling factors
<details>
<summary><i>Notes on using Excel to calculate scaling factors</i></summary>

Calculate the scaling factor by hand using the table of tallies/calculations
1. Calculate the quotient of a given sample's input SP-to-SC ratio and its ChIP SP-to-SC ratio&mdash;<i>this is the <b>scaling factor</i></b>
2. For example, for sample `Brn1_Q_rep1`, divide the `input` `SP-to-SC` value by the `ChIP` `SP-to-SC` value: 0.019802 ÷ 0.045415 = <b><i>0.436023</i></b>
3. Do this for all samples (e.g., sample-and-replicate-wise pairs of input and ChIP `bam`s)
</details>
<br />
