
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
        1. [Align, sort, and index the sample datasets](#align-sort-and-index-the-sample-datasets)
            1. [`SmartMap`'s `Bowtie2` parameters](#smartmaps-bowtie2-parameters)
                1. [Code](#code-2)
            1. [Align \(etc.\) untrimmed `fastq`s](#align-etc-untrimmed-fastqs)
                1. [Code](#code-3)
            1. [Align \(etc.\) un-trimmed `fastq`s](#align-etc-un-trimmed-fastqs)
                1. [Code](#code-4)
    1. [Examine flags in bam outfiles](#examine-flags-in-bam-outfiles)
        1. [Initialize necessary functions](#initialize-necessary-functions)
            1. [Code](#code-5)
        1. [Initialize an array of bams](#initialize-an-array-of-bams)
            1. [Code](#code-6)
        1. [Check on flag information in `bam`s](#check-on-flag-information-in-bams)
            1. [Code](#code-7)
1. [Tally/calculate alignments](#tallycalculate-alignments)
    1. [Tally/calculate alignments](#tallycalculate-alignments-1)
        1. [Initialize functions for doing floating point arithmetic, etc.](#initialize-functions-for-doing-floating-point-arithmetic-etc)
            1. [Code](#code-8)
        1. [Get situated, then initialize arrays, variables, etc.](#get-situated-then-initialize-arrays-variables-etc)
            1. [Code](#code-9)
        1. [Generate tab-separated table of alignment tallies/calculations](#generate-tab-separated-table-of-alignment-talliescalculations)
            1. [Code](#code-10)
        1. [Calculate CC/SS-styled scaling factors](#calculate-ccss-styled-scaling-factors)
            1. [Code](#code-11)
    1. [On calculating spike-in-derived scaling factors](#on-calculating-spike-in-derived-scaling-factors)
        1. [Email from Christine \(edited by me\)](#email-from-christine-edited-by-me)
        1. [Notes on using Excel to calculate spike-in-derived scaling factors](#notes-on-using-excel-to-calculate-spike-in-derived-scaling-factors)
        1. [Messages and notes associated with Biostars post on calculating spike-in-derived scaling factors](#messages-and-notes-associated-with-biostars-post-on-calculating-spike-in-derived-scaling-factors)
            1. [Initial post](#initial-post)
            1. [Notes associated with initial post](#notes-associated-with-initial-post)
            1. [Cleaned-up answer to initial post](#cleaned-up-answer-to-initial-post)
            1. [Question #1 in response to answer: How do things look if we take the fly-to-human ratio instead of the fly-to-all ratio?](#question-1-in-response-to-answer-how-do-things-look-if-we-take-the-fly-to-human-ratio-instead-of-the-fly-to-all-ratio)
            1. [Question #1: Response #1](#question-1-response-1)
            1. [Question #1: Response #2](#question-1-response-2)
            1. [Question #1: Response #3](#question-1-response-3)
            1. [Question #1: Response #4](#question-1-response-4)
            1. [Notes on Question #1: Response #4](#notes-on-question-1-response-4)
            1. [Question #1: Response #5](#question-1-response-5)
            1. [Question #2 in response to answer: On the use of alignment-ratio scaling factors versus `DESeq2::estimateSizeFactors()` scaling factors](#question-2-in-response-to-answer-on-the-use-of-alignment-ratio-scaling-factors-versus-deseq2estimatesizefactors-scaling-factors)
            1. [Question #2: Response #1](#question-2-response-1)
            1. [Question #3 in response to answer: On scaling input coverage and visualizing scaled coverage](#question-3-in-response-to-answer-on-scaling-input-coverage-and-visualizing-scaled-coverage)
            1. [Question #3: Response #1](#question-3-response-1)
        1. [How scaling factors are calculated by Egan et al., *PLOS One* 2016-1122](#how-scaling-factors-are-calculated-by-egan-et-al-plos-one-2016-1122)

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

run_previous_sym_approach=false
if ${run_previous_sym_approach}; then
    p_data="${HOME}/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802"
    p_rename="${p_data}/renamed"

    [[ ! -d "${p_rename}" ]] && mkdir -p "${p_rename}"

    unset data && typeset -A data=(
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
        echo "
          key  ${i}
        value  ${data[${i}]}
        "

        ln -s "${i}" "${data["${i}"]}"
    done
fi

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

#UPDATE
#  2024-0522: re-symlinking the data in new directory 2023_tutorial_ChIP-seq in
#+ the style of the rest of the data for the manuscript
run_new_sym_approach=true
if ${run_new_sym_approach}; then
    p_data="${HOME}/projects-etc/2023_rDNA/data/PRJNA471802/sym"
    p_sym="${HOME}/projects-etc/2023_tutorial_ChIP-seq/01_sym"

    unset data && typeset -A data=(
        ["${p_data}/Ch_log_WT_Brn1_rep1.fastq.gz"]="${p_sym}/IP_log_Brn1_rep1.fastq.gz"
        ["${p_data}/Ch_log_WT_Brn1_rep2.fastq.gz"]="${p_sym}/IP_log_Brn1_rep2.fastq.gz"
        ["${p_data}/Ch_log_WT_Brn1_rep3.fastq.gz"]="${p_sym}/IP_log_Brn1_rep3.fastq.gz"
        ["${p_data}/Ch_log_WT_Brn1_repM.fastq.gz"]="${p_sym}/IP_log_Brn1_repM.fastq.gz"
        ["${p_data}/Ch_Q_WT_Brn1_rep1.fastq.gz"]="${p_sym}/IP_Q_Brn1_rep1.fastq.gz"
        ["${p_data}/Ch_Q_WT_Brn1_rep2.fastq.gz"]="${p_sym}/IP_Q_Brn1_rep2.fastq.gz"
        ["${p_data}/Ch_Q_WT_Brn1_rep3.fastq.gz"]="${p_sym}/IP_Q_Brn1_rep3.fastq.gz"
        ["${p_data}/Ch_Q_WT_Brn1_repM.fastq.gz"]="${p_sym}/IP_Q_Brn1_repM.fastq.gz"
        ["${p_data}/in_log_WT_Brn1_rep1.fastq.gz"]="${p_sym}/in_log_Brn1_rep1.fastq.gz"
        ["${p_data}/in_log_WT_Brn1_rep2.fastq.gz"]="${p_sym}/in_log_Brn1_rep2.fastq.gz"
        ["${p_data}/in_log_WT_Brn1_rep3.fastq.gz"]="${p_sym}/in_log_Brn1_rep3.fastq.gz"
        ["${p_data}/in_log_WT_Brn1_repM.fastq.gz"]="${p_sym}/in_log_Brn1_repM.fastq.gz"
        ["${p_data}/in_Q_WT_Brn1_rep1.fastq.gz"]="${p_sym}/in_Q_Brn1_rep1.fastq.gz"
        ["${p_data}/in_Q_WT_Brn1_rep2.fastq.gz"]="${p_sym}/in_Q_Brn1_rep2.fastq.gz"
        ["${p_data}/in_Q_WT_Brn1_rep3.fastq.gz"]="${p_sym}/in_Q_Brn1_rep3.fastq.gz"
        ["${p_data}/in_Q_WT_Brn1_repM.fastq.gz"]="${p_sym}/in_Q_Brn1_repM.fastq.gz"
    )

    for i in "${!data[@]}"; do
        echo "
        Symlinking...
          key  ${i}
        value  ${data[${i}]}
        "

        ln -s "${i}" "${data["${i}"]}"
    done
fi
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
d_work="${HOME}/projects-etc/2023_tutorial_ChIP-seq"  # cd "${d_work}"
# ls -lhaFG "${d_work}"

if [[ ! -d "${d_work}/03_bam/bowtie2/err_out" ]]; then
    mkdir -p ${d_work}/03_bam/bowtie2/{bam,cvrg,err_out,macs3,ppqt,qc,siQ-ChIP}
fi

#  Go to work directory
cd "${d_work}" || echo "cd'ing failed; check on this"

#  Initialize variables needed for alignment, etc.
dir_base="${HOME}/tsukiyamalab"
dir_repo="Kris/2023_tutorial_ChIP-seq_analyses"
dir_untr="01_sym"
dir_trim="02_trim"
dir_bwt2="03_bam/bowtie2"
dir_genome=${HOME}/genomes/combined_SC_SP
dir_indx="${dir_genome}/bowtie2/combined_SC_SP"
file_fasta="${dir_genome}/fasta/combined_SC_SP.fa"
file_sizes="${dir_genome}/fasta/combined_SC_SP.chrom-info.tsv"

mapq=1

threads="${SLURM_CPUS_ON_NODE:-8}"                             # echo "${threads}"
time="8:00:00"                                                 # Job time for SLURM


p_data="${HOME}/projects-etc/2023_tutorial_ChIP-seq/01_sym"    # ls -lhaFG "${p_data}"
unset fastqs && typeset -a fastqs=(
    "${p_data}"  #INPROGRESS Moving all this to tutorial.md
)

run_check=true
if ${run_check}; then
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

    echo '### for i in "${fastqs[@]}"; do ls -lhaFG "${i}"*; done ###'
    for i in "${fastqs[@]}"; do ls -lhaFG "${i}"*; done
    echo ""
fi
```
</details>
<br />

<a id="align-sort-and-index-the-sample-datasets"></a>
#### Align, sort, and index the sample datasets
<a id="smartmaps-bowtie2-parameters"></a>
##### `SmartMap`'s `Bowtie2` parameters
<a id="code-2"></a>
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

#NOTE
#  Updated Bowtie2 call for alignment of single-end sequenced reads based on
#+ the call for paired-end sequenced reads in the following Bash script:
#+ `align-process-etc_fastqs_bowtie2.sh``
#+ ❯ bowtie2 \
#+ >     -p "${threads}" \
#+ >     -x "${f_indices}" \
#+ >     --very-sensitive-local 
#+ >     --no-unal \
#+ >     --phred33 \
#+ >     -U "${i}"
```
</details>
<br />

<a id="align-etc-untrimmed-fastqs"></a>
##### Align (etc.) untrimmed `fastq`s
<a id="code-3"></a>
###### Code
<details>
<summary><i>Code</i></summary>

```bash
#!/bin/bash

#  Run print tests to check that the commands are correct/reasonable
print_iteration=true
check_variables=true
for i in "${!fastqs[@]}"; do
    # i=0
    index="${i}"
    iter=$(( index + 1 ))
    file="${fastqs[${index}]}"            # echo "${file}"
    stem="$(basename ${file})"            # echo "${stem}"
    job_name="align-process-etc.${stem}"  # echo "${job_name}"
    
    fastq="${file}.fastq.gz"

    bam="${dir_bwt2}/bam/${stem}.bam"
    bam_coor="${bam/.bam/.sort-coord.bam}"
    bam_quer="${bam/.bam/.sort-qname.bam}"
    
    # bed_siQ="${dir_bwt2}/siQ-ChIP/${stem}.bed.gz"
    bed_etc="${dir_bwt2}/cvrg/${stem}"
    
    txt_met="${dir_bwt2}/qc/${stem}.picard-metrics.txt"
    txt_flg="${dir_bwt2}/qc/${stem}.samtools-flagstat.txt"
    txt_idx="${dir_bwt2}/qc/${stem}.samtools-idxstats.txt"
    txt_pre="${dir_bwt2}/qc/${stem}.preseq"

    #  Echo current iteration
    if ${print_iteration}; then
        echo "
        #  -------------------------------------
        ### ${iter} ###
        "
    fi

    #  Echo loop-dependent variables if check_variables is true
    if ${check_variables}; then
        echo "
        index=${index}
        iter=${iter}
        file=${file}
        stem=${stem}
        job_name=${job_name}
        
        fastq=${fastq}

        bam=${bam}
        bam_coor=${bam_coor}
        bam_quer=${bam_quer}
        # bed_siQ=${bed_siQ}
        bed_etc=${bed_etc}

        txt_met=${txt_met}
        txt_flg=${txt_flg}
        txt_idx=${txt_idx}
        txt_pre=${txt_pre}

        dir_base=${dir_base}
        dir_repo=${dir_repo}
        dir_work=${dir_work}
        dir_untr=${dir_untr}
        dir_trim=${dir_trim}
        dir_bwt2=${dir_bwt2}
        dir_indx=${dir_indx}
        file_fasta=${file_fasta}
        file_sizes=${file_sizes}

        mapq=${mapq}

        time=${time}
        threads=${threads}
        "
    fi
done
```
</details>
<br />

<a id="align-etc-un-trimmed-fastqs"></a>
##### Align (etc.) un-trimmed `fastq`s
<a id="code-4"></a>
###### Code
<details>
<summary><i>Code</i></summary>

```bash
#!/bin/bash

#  Run print tests to check that the commands are correct/reasonable
print_test=true
if ${print_test}; then
    for i in "${fastqs[@]}"; do
        # i="${fastqs[0]}"  # ., "${i}"
        in_fastq="${i}"  # ., "${in_fastq}"
        out_bam="${d_bams}/$(basename "${in_fastq}" .fastq.gz).bam"  # echo "${out_bam}"

        echo "
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
        "
    done
fi

run=true
if ${run}; then
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
fi
```
</details>
<br />

<a id="examine-flags-in-bam-outfiles"></a>
### Examine flags in bam outfiles
<a id="initialize-necessary-functions"></a>
#### Initialize necessary functions
<a id="code-5"></a>
##### Code
<details>
<summary><i>Code</i></summary>

```bash
#!/bin/bash

check_dependency() {
    what="
    check_dependency()
    ------------------
    Check if a program is available in the current environment
    
    :param 1: program to check <chr>
    
    #DONE Check that params are not empty
    #TODO Check that params are appropriate formats/strings
    #TODO Print to stderr
    "

    warning="
    WARNING: param 1 is empty; stopping the function
    "

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
    what="
    calculate_run_time()
    --------------------
    Calculate run time for chunk of code
    
    :param 1: start time in \$(date +%s) format
    :param 2: end time in \$(date +%s) format
    :param 3: message to be displayed when printing the run time <chr>
    
    #DONE Check that params are not empty
    #TODO Check that params are appropriate formats/strings
    #TODO Print to stderr
    "

    warning="
    WARNING: param(s) 1, 2, and/or 3 is/are empty; stopping the function
    "

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
    what="
    display_spinning_icon()
    -----------------------
    Display \"spinning icon\" while a background process runs
    
    :param 1: PID of the last program the shell ran in the background <pos int>
    :param 2: message to be displayed next to the spinning icon <chr>

    #DONE Check that params are not empty
    #TODO Check that params are appropriate formats/strings
    #TODO Print to stderr
    "

    warning="
    WARNING: param(s) 1 and/or 2 is/are empty; stopping the function
    "

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
    what="
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
    "

    warning_param="
    WARNING: param 1 is empty; stopping the function
    "

    warning_file="
    WARNING: param 1 file not found; stopping the function
    "

    warning_depend="
    WARNING: One or more dependencies not found; stopping the function
    "

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
    what="
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
    "

    warning_param="
    WARNING: param 1 is empty; stopping the function
    "

    warning_file="
    WARNING: param 1 file not found; stopping the function
    "

    warning_depend="
    WARNING: One or more dependencies not found; stopping the function
    "

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
<a id="code-6"></a>
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

run_check=true
if ${run_check}; then
    for i in "${bams[@]}"; do ls -lhaFG "${i}"; done
fi
```
</details>
<br />

<a id="check-on-flag-information-in-bams"></a>
#### Check on flag information in `bam`s
<a id="code-7"></a>
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
<a id="code-8"></a>
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
<a id="code-9"></a>
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

run_check=false
if ${run_check}; then
    for i in "${bams[@]}"; do ls -lhaFG "${i}"; done
fi

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
<a id="code-10"></a>
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
    print_test=false
    if ${print_test}; then
        echo "
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
        "
    fi


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
    print_test=true
    if ${print_test}; then
        echo "
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
        "
    fi

    #  Run print test
    print_test=false
    if ${print_test}; then
        echo "${sample} ${tallies[${unmapped}]} ${tallies[${mapped}]} ${tallies[${total_with_Mito}]} ${tallies[${total_SC_with_Mito}]} ${tallies[${total_SP_with_Mito}]} ${tallies[${total_sans_Mito}]} ${tallies[${total_SC_sans_Mito}]} ${tallies[${total_SP_sans_Mito}]} ${tallies[${total_SC_Mito}]} ${tallies[${total_SP_Mito}]} ${tallies[${total_1}]} ${tallies[${total_2}]} ${tallies[${total_3}]} ${tallies[${total_23}]} ${tallies[${total_30}]} ${tallies[${total_40}]} ${tallies[${total_42}]} ${tallies[${SC_0}]} ${tallies[${SC_1}]} ${tallies[${SC_2}]} ${tallies[${SC_3}]} ${tallies[${SC_23}]} ${tallies[${SC_30}]} ${tallies[${SC_40}]} ${tallies[${SC_42}]} ${tallies[${SP_0}]} ${tallies[${SP_1}]} ${tallies[${SP_2}]} ${tallies[${SP_3}]} ${tallies[${SP_23}]} ${tallies[${SP_30}]} ${tallies[${SP_40}]} ${tallies[${SP_42}]} ${tallies[${CC_SP2SC_0}]} ${tallies[${CC_SP2SC_1}]} ${tallies[${CC_SP2SC_2}]} ${tallies[${CC_SP2SC_3}]} ${tallies[${CC_SP2SC_23}]} ${tallies[${CC_SP2SC_30}]} ${tallies[${CC_SP2SC_40}]} ${tallies[${CC_SP2SC_42}]}"
        echo ""
    fi

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
<a id="code-11"></a>
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
print_test=false
if ${print_test}; then
    echo_test "${rat_ctrl[@]}" && echo ""
    echo_test "${rat_test[@]}"
fi

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

<a id="on-calculating-spike-in-derived-scaling-factors"></a>
### On calculating spike-in-derived scaling factors
<a id="email-from-christine-edited-by-me"></a>
#### Email from Christine (edited by me)
<details>
<summary><i>Email from Christine (edited by me)</i></summary>
<br />

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
5m scaling factor = 0.35/0.26 = 1.35  # Use this number for the scaling factor in deepTools (control/test)
```
I hope this makes sense. Please let me know if not&mdash;happy to chat!

Thanks!

-Christine
</details>
<br />

<a id="notes-on-using-excel-to-calculate-spike-in-derived-scaling-factors"></a>
#### Notes on using Excel to calculate spike-in-derived scaling factors
<details>
<summary><i>Notes on using Excel to calculate scaling factors</i></summary>
<br />

Calculate the scaling factor by hand using the table of tallies/calculations
1. Calculate the quotient of a given sample's input SP-to-SC ratio and its IP SP-to-SC ratio&mdash;<i>this is the <b>scaling factor</i></b>
2. For example, for sample `Brn1_Q_rep1`, divide the `input` `SP-to-SC` value by the `IP` `SP-to-SC` value: 0.019802 ÷ 0.045415 = <b><i>0.436023</i></b>
3. Do this for all samples (e.g., sample-and-replicate-wise pairs of input and IP `bam`s)
</details>
<br />

<a id="messages-and-notes-associated-with-biostars-post-on-calculating-spike-in-derived-scaling-factors"></a>
#### Messages and notes associated with [Biostars post](https://www.biostars.org/p/9572653/) on calculating spike-in-derived scaling factors
<a id="initial-post"></a>
##### [Initial post](https://www.biostars.org/p/9572653/#9572653)
<details>
<summary><i>Initial post</i></summary>
<br />

*[kalavattam](https://www.biostars.org/u/53721/)*:

1. For each sample, tally the numbers of genome-wide quality-checked (QC'd) alignments for input, immunoprecipitate (IP), spike-in input, and spike-in IP.
```R
## R ##
 main_input <- # Integer value of QC'd genome-wide alignments for input "main"
    main_IP <- # Integer value of QC'd genome-wide alignments for IP "main"

spike_input <- # Integer value of QC'd genome-wide alignments for input spike-in
   spike_IP <- # Integer value of QC'd genome-wide alignments for IP spike-in
```
*These alignment counts represent the raw data for each sample.*

2. Calculate the ratio of spike-in input to "main" input, and calculate the ratio spike-in IP to "main" IP.
```R
## R ##
ratio_input <- spike_input / main_input
   ratio_IP <- spike_IP / main_IP
```
*The assumption here is that the spike-in controls are proportional to the actual chromatin content and should ideally represent similar scaling factors.*

3. For the sample IP, calculate a scaling factor by dividing the input ratio by the IP ratio.
```R
## R ##
      sf_IP <- ratio_input / ratio_IP
```
*This scaling factor aims to correct for potential differences in IP efficiency and sequencing depth between the input and IP samples.*

4. For the sample input, the scaling factor is set to 1.
```R
## R ##
   sf_input <- 1 # (i.e., from rat_IP / rat_IP)
```
*Since the scaling factor for the input is calculated as 1, the input coverage will remain unchanged after scaling.*

5. Using deepTools `bamCoverage`, compute IP coverage by multiplying by the IP scaling factor.
```bash
## shell ##

main_IP_bam= # QC'd IP alignments to "main" genome (spike-in IP alignments have been excluded)
   bin_size= # For example, 1
      sf_IP= # Value calculated in step #3
 main_IP_bw= # Bigwig file of scaled coverage 

bamCoverage \
    -b "${main_IP_bam}" \
    --binSize "${bin_size}" \
    --scaleFactor ${scal_fctr} \
    -o "${main_IP_bw}"
```
*This should theoretically normalize the IP coverage based on the spike-in controls.*

I've done some searching, but so far I couldn't find any publications or resources that describe this particular normalization method. I'm curious to know if any of you have come across or used this approach in your ChIP-seq analyses. Are there any potential benefits or limitations associated with this method?

I will appreciate any insights and feedback on this topic, which will help me make informed decisions about normalization approaches for my ChIP-seq data. Thank you in advance for your time and assistance.
</details>
<br />

<a id="notes-associated-with-initial-post"></a>
##### Notes associated with [initial post](https://www.biostars.org/p/9572653/#9572653)
<details>
<summary><i>Notes associated with initial post</i></summary>
<br />

Strengths and assumptions for this normalization approach
- This method attempts to account for differences in sequencing depth and capture efficiency between IP and input samples using spike-in controls.
- It is based on the assumption that the spike-in controls accurately reflect the differences in IP efficiency and library preparation between samples.
- By scaling IP coverage based on input scaling factors, it theoretically reduces potential biases and enhances the comparability of different IP samples.

Weaknesses and considerations for this normalization approach
- The accuracy of the normalization heavily relies on the assumption that spike-in controls are consistent across all samples. If spike-ins vary in proportion or quality between samples, the normalization may introduce errors.
- This method does not account for biological variability or differences in IP efficiency that are not captured by spike-in controls.
    + Antibody variability: The efficiency of the IP step can vary between samples due to differences in antibody affinity, specificity, and the extent of cross-reactivity. Spike-in controls primarily address technical variations and sequencing depth but do not account for variations in antibody-antigen interactions, which can lead to differences in the IP efficiency for different samples.
    + Sample-specific conditions: Each sample might have unique characteristics that affect the efficiency of the IP process. These could include differences in chromatin accessibility, DNA fragmentation, and local chromatin structure, which may not be directly reflected by the spike-in controls.
    + Biological variability: Spike-in controls mainly address technical variation, and they might not capture inherent biological variability between samples. Biological differences in chromatin structure, epigenetic modifications, and transcription factor binding can lead to varying efficiencies in the IP process that are not addressed by the normalization based solely on spike-ins.
    + Non-specific binding: Non-specific binding of antibodies to unintended targets can also lead to differences in IP efficiency. Such non-specific interactions might not be consistently represented by spike-in controls.
- It assumes that the IP scaling factor applies uniformly across the entire genome, which may not be the case in regions with varying IP efficiencies.
- ~~The scaling factor can magnify noise and artifacts present in the spike-in controls, potentially introducing bias into the normalization process.~~
- ~~It might perform well when technical variations dominate the dataset, but biological variations may not be fully corrected.~~
</details>
<br />

<a id="cleaned-up-answer-to-initial-post"></a>
##### Cleaned-up [answer](https://www.biostars.org/p/9572653/#9572655) to [initial post](https://www.biostars.org/p/9572653/)
<details>
<summary><i>Cleaned-up answer to initial post</i></summary>
<br />

*[jared.andrews07](https://www.biostars.org/u/40195/)*:

This is largely how our group uses them to account for genome-wide shifts, as [done in this study](https://www.sciencedirect.com/science/article/pii/S1535610818305361?via%3Dihub#sec5). If we don't normalize the tracks this way, they look largely identical, as the binding profiles are largely similar, just very diminished genome-wide.

We also normalize the scaling factors to each other so that the highest scaling factor is set to 1. It feels better to scale down data than to scale up:

1. Calculate the percent *Drosophila* reads in the IP (from the read counts with duplicates removed)
2. Calculate the percent *Drosophila* reads in the input (from the read counts with duplicates removed)
3. Scaling factor = input fly % / IP fly %
4. For groups of samples that are being compared, all the scaling factors can be normalized by dividing by the highest or lowest scaling factor&mdash;this choice is made depending on how the downstream application will use the scaling factor. Ideally you should scale files down, not falsely inflate them. (For deepTools `bamCoverage`, you should divide all scaling factors by the highest scaling factor).

An example calculation:

Determine % of fly reads for each IP and input:

| Sample  | Human Reads | Fly Reads | Calculation           | % Fly Content |
|---------|-------------|-----------|-----------------------|---------------|
| IP_1    | 30,000,000  | 700,000   | 700,000/30,700,000    | 2.28%         |
| In_1    | 10,000,000  | 400,000   | 400,000/10,400,000    | 3.85%         |
| IP_2    | 40,000,000  | 3,000,000 | 3,000,000/43,000,000  | 6.98%         |
| In_2    | 10,000,000  | 450,000   | 450,000/10,450,000    | 4.31%         |
| IP_3    | 25,000,000  | 1,000,000 | 1,000,000/26,000,000  | 3.85%         |
| In_3    | 10,000,000  | 380,000   | 380,000/10,380,000    | 3.66%         |


Determine scaling factor for each IP (input fly % / IP fly % = scaling factor):

| Sample  | Calculation | Scaling Factor |
|---------|-------------|----------------|
| IP_1    | 3.85/2.28   | 1.69           |
| IP_2    | 4.31/6.98   | 0.62           |
| IP_3    | 3.66/3.85   | 0.95           |

Normalize scaling factors by dividing all scaling factors by highest scaling factor (this is used for `bamCoverage`&mdash;if the scaling factor is in the denominator instead, scale by the lowest scale factor):

| Sample  | Calculation | Normalized Scaling Factor |
|---------|-------------|---------------------------|
| IP_1    | 1.69/1.69   | 1.00                      |
| IP_2    | 0.62/1.69   | 0.367                     |
| IP_3    | 0.95/1.69   | 0.562                     |
</details>
<br />

<a id="question-1-in-response-to-answer-how-do-things-look-if-we-take-the-fly-to-human-ratio-instead-of-the-fly-to-all-ratio"></a>
##### [Question #1 in response to answer](https://www.biostars.org/p/9572653/): How do things look if we take the fly-to-human ratio instead of the fly-to-all ratio?
<details>
<summary><i>Question #1 in response to answer: How do things look if we take the fly-to-human ratio instead of the fly-to-all ratio?</i></summary>
<br />

*[kalavattam](https://www.biostars.org/u/53721/)*:

Thank you for this. In your example, you divide the fly counts by all counts (i.e., human plus fly counts) rather than human counts alone, which was what my colleagues described. Is there a reason you use the former denominator over the latter?

Because the spike-in alignment counts are low to begin with, the resulting scaling factors are similar if using (a) only human counts or (b) all counts for the denominator.

**(a) Using human counts in the denominator**

Determine % of fly reads for each IP and input:

| Sample  | Human Reads | Fly Reads | Calculation           | % Fly Content |
|---------|-------------|-----------|-----------------------|---------------|
| IP_1    | 30,000,000  | 700,000   | 700,000 / 30,000,000  | 2.33%         |
| Input_1 | 10,000,000  | 400,000   | 400,000 / 10,000,000  | 4.00%         |
| IP_2    | 40,000,000  | 3,000,000 | 3,000,000 / 40,000,000| 7.50%         |
| Input_2 | 10,000,000  | 450,000   | 450,000 / 10,000,000  | 4.50%         |
| IP_3    | 25,000,000  | 1,000,000 | 1,000,000 / 25,000,000| 4.00%         |
| Input_3 | 10,000,000  | 380,000   | 380,000 / 10,000,000  | 3.80%         |

Determine scaling factor for each IP (input fly % / IP fly % = scaling factor):

| Sample  | Calculation | Scaling Factor |
|---------|-------------|----------------|
| IP_1    | 4.00 / 2.33 | 1.72           |
| IP_2    | 4.50 / 7.50 | 0.600          |
| IP_3    | 3.80 / 4.00 | 0.950          |

Normalize scaling factors by dividing all scaling factors by highest scaling factor (this is used for `bamCoverage`&mdash;if the scaling factor is in the denominator instead, scale by the lowest scale factor):

| Sample  | Calculation | Normalized Scaling Factor |
|---------|-------------|---------------------------|
| IP_1    | 1.72 / 1.72 | 1.00                      |
| IP_2    | 0.60 / 1.72 | 0.349                     |
| IP_3    | 0.95 / 1.72 | 0.552                     |

**(b) Using all counts in the denominator**

Determine % of fly reads for each IP and input:

| Sample  | Human Reads | Fly Reads | Calculation          | % Fly Content |
|---------|-------------|-----------|----------------------|---------------|
| IP_1    | 30,000,000  | 700,000   | 700,000/30,700,000   | 2.28%         |
| Input_1 | 10,000,000  | 400,000   | 400,000/10,400,000   | 3.85%         |
| IP_2    | 40,000,000  | 3,000,000 | 3,000,000/43,000,000 | 6.98%         |
| Input_2 | 10,000,000  | 450,000   | 450,000/10,450,000   | 4.31%         |
| IP_3    | 25,000,000  | 1,000,000 | 1,000,000/26,000,000 | 3.85%         |
| Input_3 | 10,000,000  | 380,000   | 380,000/10,380,000   | 3.66%         |

Determine scaling factor for each IP (input fly % / IP fly % = scaling factor):

| Sample  | Calculation | Scaling Factor |
|---------|-------------|----------------|
| IP_1    | 3.85/2.28   | 1.69           |
| IP_2    | 4.31/6.98   | 0.62           |
| IP_3    | 3.66/3.85   | 0.95           |

Normalize scaling factors by dividing all scaling factors by highest scaling factor (this is used for `bamCoverage`&mdash;if the scaling factor is in the denominator instead, scale by the lowest scale factor):

| Sample  | Calculation | Normalized Scaling Factor |
|---------|-------------|---------------------------|
| IP_1    | 1.69/1.69   | 1.00                      |
| IP_2    | 0.62/1.69   | 0.367                     |
| IP_3    | 0.95/1.69   | 0.562                     |
</details>
<br />

<a id="question-1-response-1"></a>
##### [Question #1: Response #1](https://www.biostars.org/p/9572653/#9572681)
<details>
<summary><i>Question #1: Response #1</i></summary>
<br />

*[jared.andrews07](https://www.biostars.org/u/40195/)*:

> Thank you for this. In your example, you divide the fly counts by all counts (i.e., human plus fly counts) rather than human counts alone, which was what my colleagues described. Is there a reason you use the former denominator over the latter?

Well, the human and fly reads are still from one sample, so it makes sense to take them both into account to address total read depth differences (as one sample could have many more fly reads than another if the global shift is large). We can get pretty dramatic shifts where fly counts will be >5$\times$ higher in a sample with widespread depletion than a wildtype [sample] because the IP doesn't work as well due to much less protein. That represents a biological shift and is information that shouldn't be tossed out.

In this example, it doesn't make much of a difference in the final number, but in some experiments, I expect it'd be a more significant shift. In your example, you can already see a larger impact on samples with a greater proportion of fly reads&mdash;[e.g., a] change of 0.018 for `IP_2`, 0.01 for `IP_3`.

I can't think of any compelling reasons to not include the fly reads&mdash;it may be worth asking your colleagues why they do not.
</details>
<br />

<a id="question-1-response-2"></a>
##### [Question #1: Response #2](https://www.biostars.org/p/9572653/#9572685)
<details>
<summary><i>Question #1: Response #2</i></summary>
<br />

*[kalavattam](https://www.biostars.org/u/53721/)*:

Thank you, makes sense.

> In this example, it doesn't make much of a difference in the final number, but in some of experiments, I expect it'd be a more significant shift. In your example, you can already see a larger impact on samples with a greater proportion of fly reads&mdash;[e.g., a] change of 0.018 for `IP_2`, 0.01 for `IP_3`.
> 
> I can't think of any compelling reasons to not include the mouse reads&mdash;it may be worth asking your colleagues why they do not.

Good point, I will do so.
</details>
<br />

<a id="question-1-response-3"></a>
##### [Question #1: Response #3](https://www.biostars.org/p/9572653/#9574493)
<details>
<summary><i>Question #1: Response #3</i></summary>
<br />

*[kalavattam](https://www.biostars.org/u/53721/)*:

Thanks again, [jared.andrews07](https://www.biostars.org/u/40195/).

> In this example, it doesn't make much of a difference in the final number, but in some experiments, I expect it'd be a more significant shift. In your example, you can already see a larger impact on samples with a greater proportion of fly reads&mdash;[e.g., a] change of 0.018 for `IP_2`, 0.01 for `IP_3`.

In your previous work, did you happen to work with any ChIP-seq datasets where this shift is demonstrable?

One of the goals of my current work is&mdash;in addition to standard analyses&mdash;to provide my colleagues with a kind of tutorial for basic ChIP-seq analyses. I'd like to present and demonstrate the point that, with this approach to normalization, it's preferable to divide by the total counts (i.e, the "main" plus spike-in counts) rather than just the "main" counts to address read-depth differences. The shift is not apparent with the datasets I'm currently working with&mdash;i.e., dividing by either total or "main" yields similar scaling coefficients. I want to clearly show this because some colleagues have been dividing by only "main" counts and getting (apparently) reasonable results; I want to make the point that more extreme circumstances will yield inappropriate scaling coefficients.

So, if you happen to have encountered or know of any publicly available datasets useful for this purpose, would you mind to point me to their papers or accession numbers?
</details>
<br />

<a id="question-1-response-4"></a>
##### [Question #1: Response #4](https://www.biostars.org/p/9572653/#9574495)
<details>
<summary><i>Question #1: Response #4</i></summary>
<br />

*[jared.andrews07](https://www.biostars.org/u/40195/)*:

The study I linked to in my answer has some ChIP-seq data that will have some skew, as it contains data for a mutation that ablates H3K27me3 across the genome (with focal retention at certain loci, so not completely gone). So the samples with the mutation tend to have more fly reads than those without, as the IP just doesn't pull down as much.

If a rather extreme toy example is good enough, let's use the calculations in your comment above. I kicked the `IP_2` fly reads up to 9 million, and using all reads results in a normalized scaling factor of 0.144 while using only the human reads results in 0.116. Not a huge absolute difference, but noticeable compared to the difference between using all reads and only human reads seen for `IP_3` (0.028 for `IP_2` versus 0.01 for `IP_3`).
</details>
<br />

<a id="notes-on-question-1-response-4"></a>
##### Notes on [Question #1: Response #4](https://www.biostars.org/p/9572653/#9574495)
<details>
<summary><i>Notes on Question #1: Response #4</i></summary>
<br />

**(a) Using all counts in the denominator**

Determine % of fly reads for each IP and input:

| Sample  | Human Reads | Fly Reads | Calculation          | % Fly Content |
|---------|-------------|-----------|----------------------|---------------|
| IP_1    | 30,000,000  | 700,000   | 700,000/30,700,000   | 2.28%         |
| Input_1 | 10,000,000  | 400,000   | 400,000/10,400,000   | 3.85%         |
| IP_2    | 40,000,000  | <mark>9,000,000</mark> | <mark>9,000,000/49,000,000</mark> | <mark>18.37%</mark> |
| Input_2 | 10,000,000  | 450,000   | 450,000/10,450,000   | 4.31%         |
| IP_3    | 25,000,000  | 1,000,000 | 1,000,000/26,000,000 | 3.85%         |
| Input_3 | 10,000,000  | 380,000   | 380,000/10,380,000   | 3.66%         |

Determine scaling factor for each IP (input fly % / IP fly % = scaling factor):

| Sample  | Calculation | Scaling Factor |
|---------|-------------|----------------|
| IP_1    | 3.85/2.28   | 1.69           |
| IP_2    | 4.31/<mark>18.37</mark> | 0.23 |
| IP_3    | 3.66/3.85   | 0.95           |

Normalize scaling factors by dividing all scaling factors by highest scaling factor:

| Sample  | Calculation | Normalized Scaling Factor |
|---------|-------------|---------------------------|
| IP_1    | 1.69/1.69   | 1.00                      |
| IP_2    | <mark>0.23</mark>/1.69   | <mark>0.136</mark> |
| IP_3    | 0.95/1.69   | 0.562                     |

**(b) Using only human counts in the denominator**

Determine % of fly reads for each IP and input:

| Sample  | Human Reads | Fly Reads | Calculation          | % Fly Content |
|---------|-------------|-----------|----------------------|---------------|
| IP_1    | 30,000,000  | 700,000   | 700,000/30,000,000   | 2.33%         |
| Input_1 | 10,000,000  | 400,000   | 400,000/10,000,000   | 4.00%         |
| IP_2    | 40,000,000  | <mark>9,000,000</mark> | <mark>9,000,000/40,000,000</mark> | <mark>22.5%</mark> |
| Input_2 | 10,000,000  | 450,000   | 450,000/10,000,000   | 4.50%         |
| IP_3    | 25,000,000  | 1,000,000 | 1,000,000/25,000,000 | 4.00%         |
| Input_3 | 10,000,000  | 380,000   | 380,000/10,000,000   | 3.80%         |

Determine scaling factor for each IP (input fly % / IP fly % = scaling factor):

| Sample  | Calculation | Scaling Factor |
|---------|-------------|----------------|
| IP_1    | 4.00/2.33   | 1.72           |
| IP_2    | 4.50/<mark>22.5</mark> | 0.20 |
| IP_3    | 3.80/4.00   | 0.95           |

Normalize scaling factors by dividing all scaling factors by highest scaling factor:

| Sample  | Calculation | Normalized Scaling Factor |
|---------|-------------|---------------------------|
| IP_1    | 1.72/1.72   | 1.00                      |
| IP_2    | <mark>0.20</mark>/1.72   | <mark>0.116</mark> |
| IP_3    | 0.95/1.72   | 0.552                     |

**(c) Notes, comments**
> ...and using all reads results in a normalized scaling factor of 0.144 while using only the human reads results in 0.116.

Using all reads (main and spike-in), I'm not getting the 0.144 reported by Jared; instead, I am getting 0.136. However, using only main reads in the denominator of the scaling factor calculation, I am getting the 0.116 reported by Jared.
</details>
<br />

<a id="question-1-response-5"></a>
##### [Question #1: Response #5](https://www.biostars.org/p/9572653/#9574499)
<details>
<summary><i>Question #1: Response #5</i></summary>
<br />

*[kalavattam](https://www.biostars.org/u/53721/)*:

Sounds good, thank you so much!
</details>
<br />

<a id="question-2-in-response-to-answer-on-the-use-of-alignment-ratio-scaling-factors-versus-deseq2estimatesizefactors-scaling-factors"></a>
##### [Question #2 in response to answer](https://www.biostars.org/p/9572653/#9572677): On the use of alignment-ratio scaling factors versus `DESeq2::estimateSizeFactors()` scaling factors
<details>
<summary><i>Question #2 in response to answer: On the use of alignment-ratio scaling factors versus `DESeq2::estimateSizeFactors()` scaling factors</i></summary>
<br />

*[kalavattam](https://www.biostars.org/u/53721/)*:

Perhaps this question is better suited for a new post, but since it follows on this work, I will post it here for now and move it if requested.

If I wanted to use `DESeq2` for a differential binding analysis by generating a matrix of IP counts per *n*-bp bins, for example:

| chr   | start | stop  | IP_1_WT | IP_2_WT | IP_3_KO | IP_4_KO |
|-------|-------|-------|---------|---------|---------|---------|
| chr1  | 1     | 150   | 300     | 400     | 1200    | 1800    |
| chr1  | 151   | 300   | 250     | 433     | 1000    | 1100    |
| ...   | ...   | ...   | ...     | ...     | ...     | ...     |
| chrX  | 45001 | 45150 | 15000   | 18000   | 5500    | 4500    |
| chrX  | 45151 | 45300 | 12500   | 17500   | 5000    | 4500    |
| ...   | ...   | ...   | ...     | ...     | ...     | ...     |

Is it reasonable for me to supply the scaling factors calculated by the above-described method rather than use the values generated by `DESeq2::estimateSizeFactors()` (using either the human "main" counts or the fly spike-in counts)?
</details>
<br />

<a id="question-2-response-1"></a>
##### [Question #2: Response #1](https://www.biostars.org/p/9572653/#9572680)
<details>
<summary><i>Question #2: Response #1</i></summary>
<br />

*[jared.andrews07](https://www.biostars.org/u/40195/)*:

I'd probably use `DiffBind` and use the spike-ins in the [normalization function](https://rdrr.io/bioc/DiffBind/man/dba.normalize.html) (it uses `DESeq2` on the backend). It limits analysis to peak regions though, so depending on what you're trying to do, [you] may have to adjust.
</details>
<br />

<a id="question-3-in-response-to-answer-on-scaling-input-coverage-and-visualizing-scaled-coverage"></a>
##### [Question #3 in response to answer](https://www.biostars.org/p/9572653/#9572962): On scaling input coverage and visualizing scaled coverage
<details>
<summary><i>Question #3 in response to answer</i></summary>
<br />

*[kalavattam](https://www.biostars.org/u/53721/)*:

Hi [jared.andrews07](https://www.biostars.org/u/40195/) &ndash; thanks again for your advice and guidance.

A couple quick follow-up questions regarding these points:

> For groups of samples that are being compared, all the scaling factors can be normalized by dividing by the highest or lowest scaling factor...
> 
> ...
> 
> Normalized scaling factors by dividing all scaling factors by highest scaling factor (this used for bamCoverage - if the scaling factor is in the denominator instead, scale by the lowest scale factor):
> 
> IP 1 normalized scaling factor = 1.69/1.69 = 1
>
> IP 2 normalized scaling factor = 0.62/1.69 = 0.367
> 
> IP 3 normalized scaling factor = 0.95/1.69 = 0.562


\#1 &ndash; If we wanted to plot input coverage alongside IP coverage, would you also scale input down in this manner?

For example...
- Input 1 normalized scaling factor = 1/1.69 = 0.592
- Input 2 normalized scaling factor = 1/1.69 = 0.592
- Input 3 normalized scaling factor = 1/1.69 = 0.592


\#2 &ndash; And if we wanted to examine the $\log_2$ ratios of IP to input with, e.g., `bigwigCompare`, we'd use the scaled-down IP and input bigwigs? (However, to identify these kinds of positional changes in signal with confidence, I think a statistical approach (e.g., via `DiffBind`, `DESeq2`, `csaw`, etc.) is better suited than visualizing ratios of IP coverage to input coverage&mdash;but I'm asking because this is what my bench biologist colleagues are used to.)
</details>
<br />

<a id="question-3-response-1"></a>
##### [Question #3: Response #1](https://www.biostars.org/p/9572653/#9573065)
<details>
<summary><i>Question #3: Response #1</i></summary>
<br />

*[jared.andrews07](https://www.biostars.org/u/40195/)*:

1. I typically don't bother normalizing the input tracks, so kind of your call. I suppose I would scale them similarly.
2. I do the comparisons with the tools you mentioned to find actual differential loci, then use the normalized tracks if I need to actually show the profile/tracks in a figure. I don't do direct comparison[s] of the tracks&mdash;they are just a nice visual aid to go along with the statistical methods.
</details>
<br />


<a id="how-scaling-factors-are-calculated-by-egan-et-al-plos-one-2016-1122"></a>
#### How scaling factors are calculated by [Egan et al., *PLOS One* 2016-1122](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0166438)
<details>
<summary><i>How scaling factors are calculated by Egan et al., PLOS One 2016-1122</i></summary>
<br />

...
</details>
<br />
