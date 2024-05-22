#!/bin/bash

#  align-process-etc_fastqs_bowtie2_SE.sh
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

This script aligns single-end FASTQ files using Bowtie2, processes the output
with SAMtools (coordinate sorting, marking duplicates, etc.), and outputs both
a queryname-sorted BAM file and a coordinate-sorted, duplicate-marked BAM file.

This script aligns paired-end FASTQ files using Bowtie2, processes the output
with SAMtools (coordinate sorting, queryname sorting, marking duplicates,
etc.), and outputs both a queryname-sorted BAM file and a coordinate-sorted,
duplicate-marked BAM file.

The script also (1) generates a BED file of fragments for use as input to
siQ-ChIP, (2) generates another BED file of per-base coverage, and (3) uses the
per-base coverage to calculate and output BG and BW files for reads per million
(RPM)-scaled coverage.

Options:
   -h, --help      Display this help message.
   -t, --threads   Number of threads to use (default: 1).
   -i, --index     Path and stem for the Bowtie 2 index (required).
   -f, --fasta     PATH for genome FASTA file (required).
   -s, --sizes     Path for chromosome sizes file (required).
   -q, --mapq      MAPQ threshold for BAM filtering (default: 1).
  -fq, --fastq     Path to the FASTQ file (required).
   -b, --bam       Output path for the BAM file (required).
  -bc, --bam_coor  Output path for the coordinate-sored BAM file.
  -bq, --bam_quer  Output path for the queryname-sored BAM file.
  -be, --bed_etc   Output path and stem for the BED file of per-base coverage,
                   related Mosdepth output files, and the RPM-scaled BG and BW
                   files derived from the above BED file.
  -tm, --txt_met   Output path for the picard CollectAlignmentSummaryMetrics
                   TXT file.
  -tf, --txt_flg   Output path for the samtools flagstat TXT file.
  -ti, --txt_idx   Output path for the samtools idxstats TXT file.
  -tl, --txt_pre   Output path and stem for the preseq lc_extrap, bound_pop,
                   and c_curve TXT files.

#TODO Options:
  -bs, --bed_siQ   Output path for the BED file used as input to siQ-ChIP.

Dependencies:
  - bedGraphToBigWig
  - bowtie2
  - mosdepth
  - picard
  - preseq
  - samtools

Note:
  If not specified, then output paths for -bc, -bq, -bs, -be, -tm, -tf, -ti,
  and -ti are derived from the -b output path.

Example:
  align-process-etc_fastqs_bowtie2_SE.sh
      --threads 4
      --index path/to/indices/file_prefix
      --fasta path/to/fasta/file.fa
      --sizes path/to/sizes/file.tsv
      --mapq 1
      --fastq path/to/fastqs/file.fq.gz
      --bam path/for/bams/file.bam
      --bam_coor path/for/bams/file.sort-coord.bam
      --bam_quer path/for/bams/file.sort-qname.bam
      --bed_etc path/for/coverage/file_prefix
      --txt_met path/for/QC/file_metrics.txt
      --txt_flg path/for/QC/file_flagstat.txt
      --txt_idx path/for/QC/file_idxstats.txt
      --txt_pre path/for/QC/file_prefix
      
#TODO Option example:
      --bed_siQ path/for/siq-ChIP/file.bed
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
    threads=${SLURM_CPUS_ON_NODE:-1}
    index="${HOME}/genomes/combined_SC_SP/bowtie2/combined_SC_SP"
    fasta="${HOME}/genomes/combined_SC_SP/fasta/combined_SC_SP.fa"
    sizes="${HOME}/genomes/combined_SC_SP/fasta/combined_SC_SP.chrom-info.tsv"
    mapq=1

    fastq="${dir_exp}/01_sym/IP_log_Brn1_rep1.fastq.gz"
    bam="${dir_exp}/03_bam/bowtie2/bam/IP_log_Brn1_rep1.bam"
    bam_coor="${dir_exp}/03_bam/bowtie2/bam/IP_log_Brn1_rep1.sort-coord.bam"
    bam_quer="${dir_exp}/03_bam/bowtie2/bam/IP_log_Brn1_rep1.sort-qname.bam"
    # bed_siQ="${dir_exp}/03_bam/bowtie2/siQ-ChIP/IP_log_Brn1_rep1.bed.gz"
    bed_etc="${dir_exp}/03_bam/bowtie2/cvrg/IP_log_Brn1_rep1"
    txt_met="${dir_exp}/03_bam/bowtie2/qc/IP_log_Brn1_rep1.picard-metrics.txt"
    txt_flg="${dir_exp}/03_bam/bowtie2/qc/IP_log_Brn1_rep1.samtools-flagstat.txt"
    txt_idx="${dir_exp}/03_bam/bowtie2/qc/IP_log_Brn1_rep1.samtools-idxstats.txt"
    txt_pre="${dir_exp}/03_bam/bowtie2/qc/IP_log_Brn1_rep1.preseq"

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
            -fq|--fastq)    fastq="${2}";    shift 2 ;;
             -b|--bam)      bam="${2}";      shift 2 ;;
            -bc|--bam_coor) bam_coor="${2}"; shift 2 ;;
            -bq|--bam_quer) bam_quer="${2}"; shift 2 ;;
            -be|--bed_etc)  bed_etc="${2}";  shift 2 ;;
            -tm|--txt_met)  txt_met="${2}";  shift 2 ;;
            -tf|--txt_flg)  txt_flg="${2}";  shift 2 ;;
            -ti|--txt_idx)  txt_idx="${2}";  shift 2 ;;
            -tl|--txt_pre)  txt_pre="${2}";  shift 2 ;;
            *) echo "Unknown parameter passed: ${1}"; exit 1 ;;
        esac
        # -bs|--bed_siQ)  bed_siQ="${2}";  shift 2 ;;
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

if [[ -z "${fastq}" ]]; then
    error_and_exit "FASTQ file path is required. Use -fq or --fastq to specify it."
fi

if [[ ! -f "${fastq}" ]]; then
    error_and_exit "Specified FASTQ file does not exist."
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

# if [[ -z "${bed_siQ}" ]]; then
#     bed_siQ="${bam/.bam/.siQ-ChIP.bed.gz}"
# fi

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
#+ retain only aligned reads (-f 2) that are greater than or equal to a
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
    bowtie2 \
        -p "${threads}" \
        -x "${index}" \
        --very-sensitive-local \
        --no-unal \
        --phred33 \
        -U "${fastq}" \
            | samtools view \
                -@ "${threads}" \
                -b \
                -q "${mapq}" \
                -o "${bam}"
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

# #  Generate a siQ-ChIP-input BED file from the queryname-sorted BAM
# if [[
#          -f "${bam_quer}" \
#     && ! -f "${bed_siQ}"
# ]]; then
#     #  Extract fragment information to create the BED file, excluding
#     #+ chromosomes starting with "SP_"
#     samtools view "${bam_quer}" \
#         | awk '{
#             if (NR % 2 == 1) {
#                 chr_1 = $3; 
#                 start_1 = $4; 
#                 len_1 = length($10);
#             } else {
#                 chr_2 = $3;
#                 start_2 = $4; 
#                 len_2 = length($10);
#                 if (chr_1 == chr_2 && substr(chr_1, 1, 3) != "SP_") {
#                     start = (start_1 < start_2) ? start_1 : start_2;
#                     end = (start_1 < start_2) ? start_2 + len_2 - 1 : start_1 + len_1 - 1;
#                     frag_length = end - start + 1; 
#                     print chr_1, start, end, frag_length;
#                 }
#             }
#         }' OFS='\t' \
#         | sort -k1,1 -k2,2n \
#         | gzip \
#             > "${bed_siQ}"
# else
#     echo "Note: BED file for siQ-ChIP exists. Skipping operations."
# fi

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
        && ! -f "${bed_etc}.rpm.bg" \
        && ! -f "${bed_etc}.rpm.bg.gz"
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
                    > "${bed_etc}.rpm.bg"
    fi

    if [[
             -f "${bed_etc}.rpm.bg" \
        && ! -f "${bed_etc}.rpm.sort.bg"
    ]]; then
        bedSort \
            "${bed_etc}.rpm.bg" \
            "${bed_etc}.rpm.sort.bg"
    fi

    if [[
           -f "${bed_etc}.rpm.bg" \
        && -f "${bed_etc}.rpm.sort.bg"
    ]]; then
        mv -f \
            "${bed_etc}.rpm.sort.bg" \
            "${bed_etc}.rpm.bg"
    fi

    if [[
             -f "${bed_etc}.rpm.bg" \
        && ! -f "${bed_etc}.rpm.bw"
    ]]; then
        bedGraphToBigWig \
            "${bed_etc}.rpm.bg" \
            "${sizes}" \
            "${bed_etc}.rpm.bw"
    fi

    if [[
             -f "${bed_etc}.rpm.bg" \
        &&   -f "${bed_etc}.rpm.bw" \
        && ! -f "${bed_etc}.rpm.bg.gz"
    ]]; then
        gzip "${bed_etc}.rpm.bg"
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
    && -f "${bed_etc}.per-base.bed.gz" \
    && -f "${bed_etc}.rpm.bg.gz" \
    && -f "${bed_etc}.rpm.bw"
]]; then
    rm "${bam}"
fi
# && -f "${bed_siQ}" \

#SCRATCH
# bamCoverage \
#     --bam "${bam_coor}" \
#     --numberOfProcessors "${threads}" \
#     --binSize "${bin_size}" \
#     --normalizeUsing BPM \
#     --outFileName "${bw}"
