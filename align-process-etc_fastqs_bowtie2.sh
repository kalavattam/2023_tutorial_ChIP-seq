#!/bin/bash

#  align-process-etc_fastqs_bowtie2.sh
#  KA

debug=true


#  Define functions ===========================================================
#  Function to return an error message and exit with exit code 1, which stops
#+ the execution of the script
function error_and_exit() {
    echo "Error: ${1}" >&2
    if ${debug}; then return 1; else exit 1; fi
}


#  Function to check if a necessary program 
function check_program_in_path() {
    local program_name="${1}"
    if ! command -v "${program_name}" &> /dev/null; then
        error_and_exit "${program_name} is not in PATH. Please install ${program_name} or add it to PATH to continue."
        exit 1
    fi
}


#  Function to print help information for the script
function show_help() {
    cat << EOM
Usage: ${0} [OPTIONS]

This script aligns paired- or single-end FASTQ files using Bowtie 2, processes
the output with SAMtools (coordinate sorting, queryname sorting, marking
duplicates, etc.), and outputs both a queryname-sorted BAM file and a
coordinate-sorted, duplicate-marked BAM file.

The script also (1) generates a BED file of fragments for use as input to siQ-
ChIP (for paired-end FASTQ input files only), (2) generates another BED file of
per-base coverage, and (3) uses the per-base coverage to calculate and output
bedGraph and bigWig files for reads per million (RPM)-scaled coverage.

Additionally, the script (1) calls picard CollectAlignmentSummaryMetrics to
assess the quality of read alignments as well as "the proportion of reads that
passed machine signal-to-noise threshold quality filters," (2) generates
Samtools flagstat and idxstats reports, and (3) calls preseq lc_extrap to
generate the expected yield for theoretical larger experiments; preseq
lc_extrap also gives bounds on the number of distinct reads in the library and
the associated confidence intervals.

Options:
   -h, --help      Display this help message.
   -t, --threads   Number of threads to use (default: 1).
   -i, --index     Path and stem for the Bowtie 2 index (required).
   -f, --fasta     PATH for genome FASTA file (required).
   -s, --sizes     Path for chromosome sizes file (required).
   -q, --mapq      MAPQ threshold for BAM filtering (default: 1).
   -m, --mode      Run Bowtie 2 in 'single' or 'paired' mode (default: paired).
  -f1, --fastq_1   Path to the first FASTQ file (required).
  -f2, --fastq_2   Path to the second FASTQ file (required).
   -b, --bam       Output path for the BAM file (required).
  -bc, --bam_coor  Output path for the coordinate-sored BAM file.
  -bq, --bam_quer  Output path for the queryname-sored BAM file.
  -bs, --bed_siQ   Output path for the BED file used as input to siQ-ChIP.
  -be, --bed_etc   Output path and stem for the BED file of per-base coverage,
                   related Mosdepth output files, and the RPM-scaled bedGraph
                   and bigWig files derived from the above BED file.
  -tm, --txt_met   Output path for the picard CollectAlignmentSummaryMetrics
                   TXT file.
  -tf, --txt_flg   Output path for the samtools flagstat TXT file.
  -ti, --txt_idx   Output path for the samtools idxstats TXT file.
  -tl, --txt_pre   Output path and stem for the preseq lc_extrap, bound_pop,
                   and c_curve TXT files.

Dependencies:
  - bedGraphToBigWig
  - bowtie2
  - mosdepth
  - picard
  - preseq
  - samtools

Note:
  - If not specified, then output paths for -bc, -bq, -bs, -be, -tm, -tf, -ti,
    and -ti are derived from the -b output path.
  - Arguments -f2 and -bs are not recognized when running in 'single' mode
    (-m 'single'). #TODO I may later implement an argument to allow for the
    input of a mean fragment length in order to generate a siQ-ChIP BED file
    from a single-end sequenced ChIP-seq alignments (or to override the use of
    TLEN with paired-end sequenced ChIP-seq alignments).

Examples:
  #  Running in 'paired' mode
  align-process-etc_fastqs_bowtie2.sh
      --threads 4
      --index path/to/indices/file_prefix
      --fasta path/to/fasta/file.fa
      --sizes path/to/sizes/file.tsv
      --mapq 1
      --mode paired
      --fastq_1 path/to/fastqs/file_1.fq.gz
      --fastq_2 path/to/fastqs/file_2.fq.gz
      --bam path/for/bams/file.bam
      --bam_coor path/for/bams/file.sort-coord.bam
      --bam_quer path/for/bams/file.sort-qname.bam
      --bed_siQ path/for/siq-ChIP/file.bed
      --bed_etc path/for/coverage/file_prefix
      --txt_met path/for/QC/file_metrics.txt
      --txt_flg path/for/QC/file_flagstat.txt
      --txt_idx path/for/QC/file_idxstats.txt
      --txt_pre path/for/QC/file_prefix

  #  Running in 'single' mode
  align-process-etc_fastqs_bowtie2.sh
      --threads 4
      --index path/to/indices/file_prefix
      --fasta path/to/fasta/file.fa
      --sizes path/to/sizes/file.tsv
      --mapq 1
      --mode single
      --fastq_1 path/to/fastqs/file.fastq.gz
      --bam path/for/bams/file.bam
      --bam_coor path/for/bams/file.sort-coord.bam
      --bam_quer path/for/bams/file.sort-qname.bam
      --bed_etc path/for/coverage/file_prefix
      --txt_met path/for/QC/file_metrics.txt
      --txt_flg path/for/QC/file_flagstat.txt
      --txt_idx path/for/QC/file_idxstats.txt
      --txt_pre path/for/QC/file_prefix
EOM
}


#  Something ==================================================================
check_program_in_path "bedGraphToBigWig"
check_program_in_path "bedSort"
check_program_in_path "bowtie2"
check_program_in_path "mosdepth"
check_program_in_path "picard"
check_program_in_path "samtools"

if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
    show_help
    if ${debug}; then :; else exit 0; fi
fi

if ${debug}; then
    #  Temporary absolute paths needed for testing
    dir_base="${HOME}/tsukiyamalab"
    dir_repo="Kris/2023_tutorial_ChIP-seq"
    dir_exp="${dir_base}/${dir_repo}"

    #  Values for debugging
    threads="${SLURM_CPUS_ON_NODE:-8}"
    index="${HOME}/genomes/combined_SC_SP/bowtie2/combined_SC_SP"
    fasta="${HOME}/genomes/combined_SC_SP/fasta/combined_SC_SP.fa"
    sizes="${HOME}/genomes/combined_SC_SP/fasta/combined_SC_SP.chrom-info.tsv"
    mapq=1

    # mode="paired"
    # fastq_1="${dir_exp}/02_trim/in_G1_Hho1_6336_R1.atria.fastq.gz"
    # fastq_2="${dir_exp}/02_trim/in_G1_Hho1_6336_R2.atria.fastq.gz"
    # bam="${dir_exp}/03_bam/bowtie2/bam/in_G1_Hho1_6336.bam"
    # bam_coor="${dir_exp}/03_bam/bowtie2/bam/in_G1_Hho1_6336.sort-coord.bam"
    # bam_quer="${dir_exp}/03_bam/bowtie2/bam/in_G1_Hho1_6336.sort-qname.bam"
    # bed_siQ="${dir_exp}/03_bam/bowtie2/siQ-ChIP/in_G1_Hho1_6336.bed.gz"
    # bed_etc="${dir_exp}/03_bam/bowtie2/cvrg/in_G1_Hho1_6336"
    # txt_met="${dir_exp}/03_bam/bowtie2/qc/in_G1_Hho1_6336.picard-metrics.txt"
    # txt_flg="${dir_exp}/03_bam/bowtie2/qc/in_G1_Hho1_6336.samtools-flagstat.txt"
    # txt_idx="${dir_exp}/03_bam/bowtie2/qc/in_G1_Hho1_6336.samtools-idxstats.txt"
    # txt_pre="${dir_exp}/03_bam/bowtie2/qc/in_G1_Hho1_6336.preseq"

    # mode="paired"
    # fastq_1="${dir_exp}/02_trim/IP_Q_Hmo1_7750_R1.atria.fastq.gz"
    # fastq_2="${dir_exp}/02_trim/IP_Q_Hmo1_7750_R2.atria.fastq.gz"
    # bam="${dir_exp}/03_bam/bowtie2/bam/IP_Q_Hmo1_7750.bam"
    # bam_coor="${dir_exp}/03_bam/bowtie2/bam/IP_Q_Hmo1_7750.sort-coord.bam"
    # bam_quer="${dir_exp}/03_bam/bowtie2/bam/IP_Q_Hmo1_7750.sort-qname.bam"
    # bed_siQ="${dir_exp}/03_bam/bowtie2/siQ-ChIP/IP_Q_Hmo1_7750.bed.gz"
    # bed_etc="${dir_exp}/03_bam/bowtie2/cvrg/IP_Q_Hmo1_7750"
    # txt_met="${dir_exp}/03_bam/bowtie2/qc/IP_Q_Hmo1_7750.picard-metrics.txt"
    # txt_flg="${dir_exp}/03_bam/bowtie2/qc/IP_Q_Hmo1_7750.samtools-flagstat.txt"
    # txt_idx="${dir_exp}/03_bam/bowtie2/qc/IP_Q_Hmo1_7750.samtools-idxstats.txt"
    # txt_pre="${dir_exp}/03_bam/bowtie2/qc/IP_Q_Hmo1_7750.preseq"

    mode="single"
    fastq_1="${dir_exp}/01_sym/IP_Q_Brn1_rep1.fastq.gz"
    if [[ -n "${fastq_2}" ]]; then unset fastq_2; fi
    bam="${dir_exp}/03_bam/bowtie2/bam/IP_Q_Brn1_rep1.bam"
    bam_coor="${dir_exp}/03_bam/bowtie2/bam/IP_Q_Brn1_rep1.sort-coord.bam"
    bam_quer="${dir_exp}/03_bam/bowtie2/bam/IP_Q_Brn1_rep1.sort-qname.bam"
    if [[ -n "${bed_siQ}" ]]; then unset bed_siQ; fi
    bed_etc="${dir_exp}/03_bam/bowtie2/cvrg/IP_Q_Brn1_rep1"
    txt_met="${dir_exp}/03_bam/bowtie2/qc/IP_Q_Brn1_rep1.picard-metrics.txt"
    txt_flg="${dir_exp}/03_bam/bowtie2/qc/IP_Q_Brn1_rep1.samtools-flagstat.txt"
    txt_idx="${dir_exp}/03_bam/bowtie2/qc/IP_Q_Brn1_rep1.samtools-idxstats.txt"
    txt_pre="${dir_exp}/03_bam/bowtie2/qc/IP_Q_Brn1_rep1.preseq"

    #  Optional: Check variable assignments 
    check_variables=true
    if ${check_variables}; then
        echo "
        ### Variable assignments for debugging ###

        threads=${threads}
        index=${index}
        fasta=${fasta}
        sizes=${sizes}
        mapq=${mapq}
        mode=${mode}
        fastq_1=${fastq_1}
        fastq_2=${fastq_2}
        bam=${bam}
        bam_coor=${bam_coor}
        bam_quer=${bam_quer}
        bed_siQ=${bed_siQ}
        bed_etc=${bed_etc}
        txt_met=${txt_met}
        txt_flg=${txt_flg}
        txt_idx=${txt_idx}
        txt_pre=${txt_pre}
        "
    fi
else
    #  Parse command line arguments
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
             -t|--threads)  threads="${2}";  shift 2 ;;
             -i|--index)    index="${2}";    shift 2 ;;
             -f|--fasta)    fasta="${2}";    shift 2 ;;
             -s|--sizes)    sizes="${2}";    shift 2 ;;
             -q|--mapq)     mapq="${2}";     shift 2 ;;
             -m|--mode)     mode="${2}";     shift 2 ;;
            -f1|--fastq_1)  fastq_1="${2}";  shift 2 ;;
            -f2|--fastq_2)  fastq_2="${2}";  shift 2 ;;
             -b|--bam)      bam="${2}";      shift 2 ;;
            -bc|--bam_coor) bam_coor="${2}"; shift 2 ;;
            -bq|--bam_quer) bam_quer="${2}"; shift 2 ;;
            -bs|--bed_siQ)  bed_siQ="${2}";  shift 2 ;;
            -be|--bed_etc)  bed_etc="${2}";  shift 2 ;;
            -tm|--txt_met)  txt_met="${2}";  shift 2 ;;
            -tf|--txt_flg)  txt_flg="${2}";  shift 2 ;;
            -ti|--txt_idx)  txt_idx="${2}";  shift 2 ;;
            -tl|--txt_pre)  txt_pre="${2}";  shift 2 ;;
            *) echo "Unknown parameter passed: ${1}"; exit 1 ;;
        esac
    done
fi

if [[ -z "${index}" ]]; then
    error_and_exit "File path and stem for Bowtie 2 indices are required. Use -i or --index to specify it."
fi

if [[ ! -d "$(dirname "${index}")" ]]; then
    error_and_exit "Specified directory for Bowtie 2 indices does not exist."
fi

if [[ -z "${sizes}" ]]; then
    error_and_exit "File path for chromosome sizes is required. Use -s or --sizes to specify it."
fi

if [[ ! -f "${sizes}" ]]; then
    error_and_exit "Specified chromosome sizes file does not exist."
fi

if [[ -z "${fasta}" ]]; then
    error_and_exit "File path for genome FASTA is required. Use -f or --fasta to specify it."
fi

if [[ ! -f "${fasta}" ]]; then
    error_and_exit "Specified genome FASTA file does not exist."
fi

if [[ -z "${mapq}" ]]; then
    mapq=1
fi

if [[ -z "${mode}" ]]; then
    mode="paired"
fi

case "${mode}" in
    single|paired) : ;;
    *) error_and_exit "Mode has been set to ${mode} but must be either 'single' or 'paired'." ;;
esac

if [[ -z "${fastq_1}" ]]; then
    error_and_exit "FASTQ_1 file path is required. Use -f1 or --fastq_1 to specify it."
fi

if [[ ! -f "${fastq_1}" ]]; then
    error_and_exit "Specified FASTQ_1 file does not exist."
fi

if [[ "${mode}" == "paired" ]]; then
    if [[ -z "${fastq_2}" ]]; then
        error_and_exit "FASTQ_2 file path is required. Use -f2 or --fastq_2 to specify it."
    fi

    if [[ ! -f "${fastq_2}" ]]; then
        error_and_exit "Specified FASTQ_2 file does not exist."
    fi
fi

if [[ -z "${bam}" ]]; then
    error_and_exit "BAM file path is required. Use -b or --bam to specify it."
fi

if [[ -z "${bam_coor}" ]]; then
    bam_coor="${bam/.bam/.sort-coord.bam}"
fi

if [[ -z "${bam_quer}" ]]; then
    bam_quer="${bam/.bam/.sort-qname.bam}"
fi

if [[ "${mode}" == "paired" ]]; then
    if [[ -z "${bed_siQ}" ]]; then
        bed_siQ="${bam/.bam/.siQ-ChIP.bed.gz}"
    fi
fi

if [[ -z "${bed_etc}" ]]; then
    bed_etc="${bam%.bam}"
fi

if [[ -z "${txt_met}" ]]; then
    txt_idx="${bam/.bam/.picard-metrics.txt}"
fi

if [[ -z "${txt_flg}" ]]; then
    txt_flg="${bam/.bam/.samtools-flagstat.txt}"
fi

if [[ -z "${txt_idx}" ]]; then
    txt_idx="${bam/.bam/.samtools-idxstats.txt}"
fi

if [[ -z "${txt_pre}" ]]; then
    txt_pre="${bam/.bam/.preseq}"
fi

#  Check if the BAM file exists; if not, perform alignment with Bowtie 2,
#+ converting the Bowtie 2 output to a BAM file with Samtools; in doing so,
#+ retain only properly paired reads (-f 2) that are greater than or equal to a
#+ user-supplied MAPQ score (-q "${mapq}") [by default, retain any-quality
#+ "maxi-reads" (-q 1)]
#+ 
#+ On how to interpret Bowtie 2 MAPQ scores:
#+ biofinysics.blogspot.com/2014/05/how-does-bowtie2-assign-mapq-scores.html
if [[
       ! -f "${bam}" \
    && ! -f "${bam_coor}" \
    && ! -f "${bam_quer}"
]]; then
    if [[ "${mode}" == "paired" ]]; then
        bowtie2 \
            -p "${threads}" \
            -x "${index}" \
            --very-sensitive-local \
            --no-unal \
            --no-mixed \
            --no-discordant \
            --no-overlap \
            --no-dovetail \
            --phred33 \
            -I 10 \
            -X 700 \
            -1 "${fastq_1}" \
            -2 "${fastq_2}" \
                | samtools view \
                    -@ "${threads}" \
                    -b \
                    -f 2 \
                    -q "${mapq}" \
                    -o "${bam}"
    elif [[ "${mode}" == "single" ]]; then
        bowtie2 \
            -p "${threads}" \
            -x "${index}" \
            --very-sensitive-local \
            --no-unal \
            --phred33 \
            -U "${fastq_1}" \
                | samtools view \
                    -@ "${threads}" \
                    -b \
                    -q "${mapq}" \
                    -o "${bam}"
    fi
else
    echo "Note: BAM file exists. Skipping alignment operations."
fi

#  Check if the BAM file exists to perform further operations, including sorting
#+ the BAM file by queryname
if [[
         -f "${bam}" \
    && ! -f "${bam_quer}"
]]; then
    samtools sort \
        -@ "${threads}" \
        -n \
        -o "${bam_quer}" \
        "${bam}"

    if [[ "${mode}" == "paired" ]]; then
        #  After sorting by queryname, fix the paired read-mate information, which
        #+ is required for subsequent operations 
        if [[ -f "${bam_quer}" ]]; then
            samtools fixmate \
                -@ "${threads}" \
                -m \
                "${bam_quer}" \
                "${bam_quer%.bam}.tmp.bam"

            #  Replace the original queryname-sorted BAM with queryname-sorted
            #+ mate-fixed BAM
            if [[ -f "${bam_quer%.bam}.tmp.bam" ]]; then
                mv -f \
                    "${bam_quer%.bam}.tmp.bam" \
                    "${bam_quer}"
            fi
        fi
    fi
else
    echo "Note: Queryname-sored BAM file exists. Skipping sort operations."
fi

#  For other downstream analyses, sort the queryname-sorted BAM by
#+ coordinates
if [[
         -f "${bam_quer}" \
    && ! -f "${bam_coor}"
]]; then
    samtools sort \
        -@ "${threads}" \
        -o "${bam_coor}" \
        "${bam_quer}"

    #  Index the coordinate-sorted BAM file
    if [[ ! -f "${bam_coor}.bai" ]]; then
        samtools index \
            -@ "${threads}" \
            "${bam_coor}"
    fi
else
    echo "Note: Coordinate-sored BAM file exists. Skipping sort operations."
fi

#  Mark duplicate alignments in the coordinate-sorted BAM file, and generate
#+ Picard, Samtools, and Preseq QC TXT files
#TODO Break this up; modularize it
if [[
         -f "${bam_coor}" \
    &&   -f "${bam_coor}.bai" \
    && ! -f "${txt_flg}" \
    && ! -f "${txt_idx}"
]]; then
    #  Mark duplicate alignments in the coordinate-sorted BAM file
    samtools markdup \
        -@ "${threads}" \
        "${bam_coor}" \
        "${bam_coor%.bam}.tmp.bam"

    #  Replace the original coordinate-sorted BAM with one in which
    #+ duplicates alignments are marked
    if [[ -f "${bam_coor%.bam}.tmp.bam" ]]; then
        mv -f \
            "${bam_coor%.bam}.tmp.bam" \
            "${bam_coor}"
    fi

    #+ Use picard CollectAlignmentSummaryMetrics to assess the quality of read
    #+ alignments as well as "the proportion of reads that passed machine
    #+ signal-to-noise threshold quality filters"
    # shellcheck disable=SC2181
    if [[ $? -eq 0 ]]; then
        #  Also, Picard CollectAlignmentSummaryMetrics requires that the
        #+ reference genome FASTA is indexed; so, generate an FAI file if
        #+ necessary
        if [[
                 -f "${bam_coor}" \
            &&   -f "${bam_coor}.bai" \
            && ! -f "${txt_met}"
        ]]; then
            if [[ "${fasta##*.}" == "gz" ]]; then
                decomp_fasta="${fasta%.gz}"
                if [[ ! -f "${decomp_fasta}.fai" ]]; then
                    if [[ ! -f "${decomp_fasta}" ]]; then
                        gunzip -c "${fasta}" > "${decomp_fasta}"
                    fi

                    samtools faidx -@ "${threads}" "${decomp_fasta}"
                fi

                picard CollectAlignmentSummaryMetrics \
                    --REFERENCE_SEQUENCE "${decomp_fasta}" \
                    --INPUT "${bam_coor}" \
                    --OUTPUT "${txt_met}"
            else
                if [[ ! -f "${fasta}.fai" ]]; then
                    samtools faidx -@ "${threads}" "${fasta}"
                fi

                picard CollectAlignmentSummaryMetrics \
                    --REFERENCE_SEQUENCE "${fasta}" \
                    --INPUT "${bam_coor}" \
                    --OUTPUT "${txt_met}"
            fi
        fi
    fi

    #  Generate Samtools flagstat and idxstats reports
    if [[
             -f "${bam_coor}" \
        &&   -f "${bam_coor}.bai" \
        && ! -f "${txt_flg}" \
        && ! -f "${txt_idx}"
    ]]; then
        samtools flagstat -@ "${threads}" "${bam_coor}" > "${txt_flg}"
        samtools idxstats "${bam_coor}" > "${txt_idx}"
    fi

    #  Use preseq lc_extrap to generate the expected yield for theoretical
    #+ larger experiments; preseq lc_extrap also gives bounds on the number of
    #+ distinct reads in the library and the associated confidence intervals
    if [[
             -f "${bam_coor}" \
        &&   -f "${bam_coor}.bai" \
        && ! -f "${txt_pre}"
    ]]; then
        preseq lc_extrap -v -P -B -r 24 \
            -o "${txt_pre}-lc-extrap.txt" \
            "${bam_coor}" \
                &> "${txt_pre%.txt}-lc-extrap-verbose.txt"

        preseq bound_pop -P -B -r 24 \
            -o "${txt_pre}-bound-pop.txt" \
            "${bam_coor}"

        preseq c_curve -P -B -r 24 \
            -o "${txt_pre}-c-curve.txt" \
            "${bam_coor}"
    fi
else
    echo "Note: Duplicate alignments are already marked in the coordinate-"
    echo "      sorted BAM file, and Picard and Samtools QC TXT files exist."
    echo "      Skipping operations."
fi

if [[ "${mode}" == "paired" ]]; then
    #  Generate a siQ-ChIP-input BED file from the queryname-sorted BAM
    if [[
             -f "${bam_quer}" \
        && ! -f "${bed_siQ}"
    ]]; then
        #  Extract fragment information to create the BED file, excluding
        #+ chromosomes starting with "SP_"
        samtools view "${bam_quer}" \
            | awk '{
                if (NR % 2 == 1) {
                    chr_1 = $3; 
                    start_1 = $4; 
                    len_1 = length($10);
                } else {
                    chr_2 = $3;
                    start_2 = $4; 
                    len_2 = length($10);
                    if (chr_1 == chr_2 && substr(chr_1, 1, 3) != "SP_") {
                        start = (start_1 < start_2) ? start_1 : start_2;
                        end = (start_1 < start_2) ? start_2 + len_2 - 1 : start_1 + len_1 - 1;
                        frag_length = end - start + 1; 
                        print chr_1, start, end, frag_length;
                    }
                }
            }' OFS='\t' \
            | sort -k1,1 -k2,2n \
            | gzip \
                > "${bed_siQ}"
    else
        echo "Note: BED file for siQ-ChIP exists. Skipping operations."
    fi
fi

#  Using the coordinate-sorted BAM, generate BED file of per-base coverage
#+ as well as RPM-scaled BG and BW files
if [[
         -f "${bam_coor}" \
    &&   -f "${bam_coor}.bai" \
    && ! -f "${bed_etc}.per-base.bed.gz"
]]; then
    #  Mosdepth throws errors if BAI is older than BAM
    if [[ "${bam_coor}.bai" -ot "${bam_coor}" ]]; then
        #  If BAI is older than BAM, then index the BAM again (current BAI will
        #+ be overwritten)
        samtools index -@ "${threads}" "${bam_coor}"
    fi

    mosdepth \
        "${bed_etc}" \
        "${bam_coor}"

    if [[
             -f "${bed_etc}.per-base.bed.gz" \
        && ! -f "${bed_etc}.rpm.bedgraph" \
        && ! -f "${bed_etc}.rpm.bedgraph.gz"
    ]]; then
        #TODO Error and exit for if total_reads is assigned 0
        total_reads="$(samtools view -c "${bam_coor}")"
        zcat "${bed_etc}.per-base.bed.gz" \
            | awk \
                -v total_reads="${total_reads}" \
                'BEGIN {
                    OFS="\t"
                } {
                    print $1, $2, $3, ($4 / total_reads) * 1000000
                }' \
                    > "${bed_etc}.rpm.bedgraph"
    fi

    if [[
             -f "${bed_etc}.rpm.bedgraph" \
        && ! -f "${bed_etc}.rpm.sort.bg"
    ]]; then
        bedSort \
            "${bed_etc}.rpm.bedgraph" \
            "${bed_etc}.rpm.sort.bg"
    fi

    if [[
           -f "${bed_etc}.rpm.bedgraph" \
        && -f "${bed_etc}.rpm.sort.bg"
    ]]; then
        mv -f \
            "${bed_etc}.rpm.sort.bg" \
            "${bed_etc}.rpm.bedgraph"
    fi

    if [[
             -f "${bed_etc}.rpm.bedgraph" \
        && ! -f "${bed_etc}.rpm.bigwig"
    ]]; then
        bedGraphToBigWig \
            "${bed_etc}.rpm.bedgraph" \
            "${sizes}" \
            "${bed_etc}.rpm.bigwig"
    fi

    if [[
             -f "${bed_etc}.rpm.bedgraph" \
        &&   -f "${bed_etc}.rpm.bigwig" \
        && ! -f "${bed_etc}.rpm.bedgraph.gz"
    ]]; then
        gzip "${bed_etc}.rpm.bedgraph"
    fi
else
    echo "Note: Coverage BED, BG, and BW files exist. Skipping operations."
fi

#  Remove the original BAM file (Bowtie 2 output piped to Samtools) if all
#+ other files have been successfully created
if [[
       -f "${bam}" \
    && -f "${bam_coor}" \
    && -f "${bam_quer}" \
    && -f "${bed_siQ}" \
    && -f "${bed_etc}.per-base.bed.gz" \
    && -f "${bed_etc}.rpm.bedgraph.gz" \
    && -f "${bed_etc}.rpm.bigwig"
]]; then
    rm "${bam}"
fi

#SCRATCH
# bamCoverage \
#     --bam "${bam_coor}" \
#     --numberOfProcessors "${threads}" \
#     --binSize "${bin_size}" \
#     --normalizeUsing BPM \
#     --outFileName "${bw}"
