#!/bin/bash

#  align_process_etc.sh
#  KA

#  Set hardcoded general flags
interact=true
verbose=true
check_variables=true

#  Set hardcoded flags for quality checks
flag_fqc_fastq=true
flag_list_tally=true
flag_fqc_bam=true
flag_met=false
flag_flg=false
flag_idx=false
flag_lc_extrap=false
flag_bound_pop=false
flag_c_curve=false
flag_siq=true
flag_etc=true
flag_clean=true


#  Define functions ===========================================================
#  Function to return an error message to stderr and exit with exit code 1,
#+ which stops the execution of the script
function error_and_exit() {
    echo "Error: ${1}" >&2
    if ${interact}; then return 1; else exit 1; fi
}


#  Function to return a "note" message to stdout
function echo_note() { echo "Note: ${1}"; }


#  Function to check if a necessary program 
function check_program_in_path() {
    local program_name="${1}"
    if ! command -v "${program_name}" &> /dev/null; then
        error_and_exit "${program_name} is not in PATH. Please install ${program_name} or add it to PATH to continue."
    fi
}


#  Function to determine the extension of a FASTQ file
determine_fastq_extension() {
    local fastq="${1}"
    local extension=""

    if [[ -z "${fastq}" ]]; then
        echo "Error: No input file specified." >&2
        return 1
    fi

    if [[ ! -f "${fastq}" ]]; then
        echo "Error: Input file ${fastq} does not exist." >&2
        return 1
    fi

    if [[ "${fastq}" == *.fastq ]]; then
        extension=".fastq"
    elif [[ "${fastq}" == *.fastq.gz ]]; then
        extension=".fastq.gz"
    elif [[ "${fastq}" == *.fq ]]; then
        extension=".fq"
    elif [[ "${fastq}" == *.fq.gz ]]; then
        extension=".fq.gz"
    else
        error_and_exit "No valid FASTQ extension found for ${fastq}."
    fi

    echo "${extension}"
}


#  Function to list and tally flags, MAPQ values, or flags stratified by MAPQ
#+ values in a BAM file
list_tally_flags_mapq() {
    local infile=""
    local outfile=""
    local threads="${SLURM_CPUS_ON_NODE:-1}"
    local what="flags"
    local by_chr=false

    #  Parse keyword arguments
    if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
        echo "Usage: list_tally_flags --infile <path> [--outfile <path>]"
        echo "           [--threads <integer >= 1>]"
        echo "  -i, --infile   Specify the input BAM file (required)."
        echo "  -o, --outfile  Specify the output TSV file to write flag tallies (required)."
        echo "  -t, --threads  Specify the number of threads to use (default: ${threads})."
        echo "  -w, --what     Specify what to tally: 'flags', 'mapq', or 'both' (default: 'flags')."
        echo "  -b, --by_chr   Stratify the output by chromosome (optional)."
        return 1
    fi

    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
            -i|--infile)  infile="${2}";  shift 2 ;;
            -o|--outfile) outfile="${2}"; shift 2 ;;
            -t|--threads) threads="${2}"; shift 2 ;;
            -w|--what)    what="${2}";    shift 2 ;;
            -b|--by_chr)  by_chr=true; shift 1 ;;
            *) echo "Unknown parameter passed: ${1}"; return 1 ;;
        esac
    done

    #  Validate required parameter: infile
    if [[ -z "${infile}" ]]; then
        echo "Error: No input file specified." >&2
        return 1
    fi

    #  Check that input file exists
    if [[ ! -f "${infile}" ]]; then
        echo "Error: Input file does not exist." >&2
        return 1
    fi

    #  Check that input file is a BAM file
    if [[ ! "${infile}" =~ \.bam$ ]]; then
        echo "Error: Input file does not have a .bam extension." >&2
        return 1
    fi

    #  Check that threads is assigned to an integer >= 1
    if ! [[ "${threads}" =~ ^[0-9]+$ ]] || [[ "${threads}" -lt 1 ]]; then
        echo "Error: The number of threads must be an integer >= 1." >&2
        return 1
    fi

    #  Process according to the specified "what" option, checking and assigning
    #+ what if/as needed)
    case "${what}" in
        flags)
            if [[ -z "${outfile}" ]]; then
                outfile="${infile/.bam/.flags.txt}"
            fi

            if ${by_chr}; then
                if ! \
                    samtools view -@ "${threads}" "${infile}" \
                        | awk \
                            -v OFS='\t' \
                            '{
                                chr = $3; flag = $2;
                                count[chr,flag]++;
                            } END {
                                for (key in count) {
                                    split(key, a, SUBSEP);  # a[1] = chr, a[2] = flag
                                    print a[1], a[2], count[key];
                                }
                            }' \
                        | sort -k1,1 -k2,2n \
                            > "${outfile}";
                then
                    echo "Error: Failed to process BAM file with Samtools." >&2
                    return 1
                fi
            else
                if ! \
                    samtools view -@ "${threads}" "${infile}" \
                        | awk \
                            -v OFS='\t' \
                            '{
                                count[$2]++
                            } END {
                                for (flag in count) print flag, count[flag]
                            }' \
                        | sort -nr -k2 \
                            > "${outfile}";
                then
                    echo "Error: Failed to process BAM file with Samtools." >&2
                    return 1
                fi
            fi ;;
        mapq|mapqs|MAPQ|MAPQs)
            if [[ -z "${outfile}" ]]; then
                outfile="${infile/.bam/.MAPQs.txt}"
            fi

            if ${by_chr}; then
                if ! \
                    samtools view -@ "${threads}" "${infile}" \
                        | awk \
                            -v OFS='\t' \
                            '{
                                chr = $3; mapq = $5;
                                count[chr,mapq]++;
                            } END {
                                for (key in count) {
                                    split(key, a, SUBSEP);  # a[1] = chr, a[2] = mapq
                                    print a[1], a[2], count[key];
                                }
                            }' \
                        | sort -k1,1 -k2,2n \
                            > "${outfile}";
                then
                    echo "Error: Failed to process BAM file with Samtools." >&2
                    return 1
                fi
            else
                if ! \
                    samtools view -@ "${threads}" "${infile}" \
                        | awk \
                            -v OFS='\t' \
                            '{
                                count[$5]++
                            } END {
                                for (mapq in count) print mapq, count[mapq]
                            }' \
                        | sort -nr -k2 \
                            > "${outfile}";
                then
                    echo "Error: Failed to process BAM file with Samtools." >&2
                    return 1
                fi
            fi ;;
        both)
            if [[ -z "${outfile}" ]]; then
                outfile="${infile/.bam/.flags-by-MAPQs.txt}"
            fi

            if ${by_chr}; then
                if ! \
                    samtools view -@ "${threads}" "${infile}" \
                        | awk \
                            -v OFS='\t' \
                            '{
                                chr = $3; flag = $2; mapq = $5;
                                count[chr,flag,mapq]++;
                            } END {
                                for (key in count) {
                                    split(key, a, SUBSEP);  # a[1] = chr, a[2] = flag, a[3] = mapq
                                    print a[1], a[2], a[3], count[key];
                                }
                            }' \
                        | sort -k1,1 -k2,2n -k3,3n \
                            > "${outfile}";
                then
                    echo "Error: Failed to process BAM file with Samtools." >&2
                    return 1
                fi
            else
                if ! \
                    samtools view -@ "${threads}" "${infile}" \
                        | awk \
                            -v OFS='\t' \
                            '{
                                flag = $2; mapq = $5;
                                count[flag,mapq]++;
                            } END {
                                for (key in count)
                                    print key, count[key];
                            }' \
                        | awk \
                            -F'\t' \
                            '{
                                split($1, a, SUBSEP);  # a[1] = flag, a[2] = MAPQ
                                print a[1], a[2], $2;
                            }' \
                        | sort -k1,1n -k2,2n \
                            > "${outfile}"
                then
                    echo "Error: Failed to process BAM file with Samtools." >&2
                    return 1
                fi
            fi ;;
        *)
            echo "Error: Invalid option specified for what." >&2
            return 1
            ;;
    esac

    echo "list_tally_flags_mapq for ${what} tallied and written to ${outfile}."
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
   -m, --mode      Run Bowtie 2 in 'single' or 'paired' mode (default: paired).
   -q, --mapq      MAPQ threshold for BAM filtering (default: 1).
  -rf, --req_flag  Require flag 2 for BAM filtering (optional).
  -fl, --f_length  Fragment length in bp (integer > 0) for the generation of a
                   siQ-ChIP BED file; if specified and mode is 'paired', then
                   the value will override the calculation of fragment length
                   from aligned read mates (default for 'single' mode: 200;
                   the default for 'paired' mode is to use values stored in BAM
                   TLEN fields).
  -f1, --fastq_1   Path to the first FASTQ file (required).
  -f2, --fastq_2   Path to the second FASTQ file (required if mode is 'paired';
                   ignored if mode is 'single').
   -b, --bam       Output path for the BAM file (required).
  -bc, --bam_coor  Output path for the coordinate-sorted BAM file.
  -bq, --bam_quer  Output path for the queryname-sorted BAM file.
  -bs, --bed_siQ   Output path for the BED file used as input to siQ-ChIP.
  -be, --bed_etc   Output path and stem for the BED file of per-base coverage,
                   related Mosdepth output files, and the RPM-scaled bedGraph
                   and bigWig files derived from the above BED file.
  -dff, --d_fqc_f  Output directory for HTML file and compressed sub-directory
                   (ZIP) from running FastQC on FASTQ file(s) (required if
                   hardcoded flag_fqc_fastq=true).
  -dfb, --d_fqc_b  Output directory for HTML file and compressed sub-directory
                   (ZIP) from running FastQC on coordinate-sorted BAM file
                   (required if hardcoded flag_fqc_bam=true).
  -tl, --txt_list  Output path for the list_tally_flags_mapq TXT file.
  -tm, --txt_met   Output path for the picard CollectAlignmentSummaryMetrics
                   TXT file.
  -tf, --txt_flg   Output path for the samtools flagstat TXT file.
  -ti, --txt_idx   Output path for the samtools idxstats TXT file.
  -tp, --txt_pre   Output path and stem for the preseq lc_extrap, bound_pop,
                   and c_curve TXT files.

Dependencies:
  - bedGraphToBigWig
  - bowtie2
  - mosdepth
  - picard
  - preseq
  - samtools

Notes:
  - FASTQ infiles are expected to have one of the following extensions:
      + *.fastq.gz
      + *.fastq
      + *.fq.gz
      + *.fq
  - If not specified, then output paths for -bc, -bq, -bs, -be, -tl, -tm, -tf,
    -ti, and -tp are derived from the -b output path.
      + #TODO Get rid of this; if not specified, then simply not run; specified
              by assigning an outfile.
  - Argument -f2 is not recognized when running in 'single' mode (-m 'single').
  - #DONE I may later implement an argument to allow for the input of a
    mean fragment length in order to generate a siQ-ChIP BED file from
    single-end sequenced ChIP-seq alignments (or to override the use of TLEN
    with paired-end sequenced ChIP-seq alignments).
  - #TODO Test the above implementation.

Examples:
  #  Running in 'paired' mode
  align-process-etc_fastqs_bowtie2.sh
      --threads 4
      --index path/to/indices/file_prefix
      --fasta path/to/fasta/file.fa
      --sizes path/to/sizes/file.tsv
      --mode paired
      --mapq 1
      --req_flag
      --fastq_1 path/to/fastqs/file_1.fq.gz
      --fastq_2 path/to/fastqs/file_2.fq.gz
      --bam path/for/bams/file.bam
      --bam_coor path/for/bams/file.sort-coord.bam
      --bam_quer path/for/bams/file.sort-qname.bam
      --bed_siQ path/for/siq-ChIP/file.bed
      --bed_etc path/for/coverage/file_prefix
      --d_fqc_f path/for/QC/fastqc/fastq
      --d_fqc_b path/for/QC/fastqc/bam
      --txt_list path/for/QC/file_list-tally.txt
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
      --mode single
      --mapq 1
      --f_length 200
      --fastq_1 path/to/fastqs/file.fastq.gz
      --bam path/for/bams/file.bam
      --bam_coor path/for/bams/file.sort-coord.bam
      --bam_quer path/for/bams/file.sort-qname.bam
      --bed_siQ path/for/siq-ChIP/file.bed
      --bed_etc path/for/coverage/file_prefix
      --d_fqc_f path/for/QC/fastqc/fastq
      --d_fqc_b path/for/QC/fastqc/bam
      --txt_list path/for/QC/file_list-tally.txt
      --txt_met path/for/QC/file_metrics.txt
      --txt_flg path/for/QC/file_flagstat.txt
      --txt_idx path/for/QC/file_idxstats.txt
      --txt_pre path/for/QC/file_prefix
EOM
}


#  #TODO Name this section ====================================================
check_program_in_path "awk"
check_program_in_path "bedGraphToBigWig"
check_program_in_path "bedSort"
check_program_in_path "bowtie2"
check_program_in_path "fastqc"
check_program_in_path "mosdepth"
check_program_in_path "picard"
check_program_in_path "samtools"
check_program_in_path "sort"

if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
    show_help
    if ! ${interact}; then exit 0; fi
fi

if ${interact}; then
    #  Temporary absolute paths used for running tests
    dir_base="${HOME}/tsukiyamalab"
    dir_repo="Kris/2023_tutorial_ChIP-seq"
    dir_exp="${dir_base}/${dir_repo}"

    #  Values for running tests
    threads="${SLURM_CPUS_ON_NODE:-1}"
    index="${HOME}/genomes/combined_SC_SP/bowtie2/combined_SC_SP"
    fasta="${HOME}/genomes/combined_SC_SP/fasta/combined_SC_SP.fa"
    sizes="${HOME}/genomes/combined_SC_SP/fasta/combined_SC_SP.chrom-info.tsv"

    # mode="paired"
    # mapq=1
    # req_flag=true
    # if [[ -n "${f_length}" ]]; then unset f_length; fi
    # fastq_1="${dir_exp}/02_trim/in_G1_Hho1_6336_R1.atria.fastq.gz"
    # fastq_2="${dir_exp}/02_trim/in_G1_Hho1_6336_R2.atria.fastq.gz"
    # bam="${dir_exp}/03_bam/bowtie2/bam/in_G1_Hho1_6336.bam"
    # bam_coor="${dir_exp}/03_bam/bowtie2/bam/in_G1_Hho1_6336.sort-coord.bam"
    # bam_quer="${dir_exp}/03_bam/bowtie2/bam/in_G1_Hho1_6336.sort-qname.bam"
    # bed_siQ="${dir_exp}/03_bam/bowtie2/siQ-ChIP/in_G1_Hho1_6336.bed.gz"
    # d_fqc_f="${dir_exp}/03_bam/bowtie2/qc/fastqc/fastq"
    # d_fqc_b="${dir_exp}/03_bam/bowtie2/qc/fastqc/bam"
    # bed_etc="${dir_exp}/03_bam/bowtie2/cvrg/in_G1_Hho1_6336"
    # txt_list="${dir_exp}/03_bam/bowtie2/qc/in_G1_Hho1_6336.list-tally.txt"
    # txt_met="${dir_exp}/03_bam/bowtie2/qc/in_G1_Hho1_6336.picard-metrics.txt"
    # txt_flg="${dir_exp}/03_bam/bowtie2/qc/in_G1_Hho1_6336.samtools-flagstat.txt"
    # txt_idx="${dir_exp}/03_bam/bowtie2/qc/in_G1_Hho1_6336.samtools-idxstats.txt"
    # txt_pre="${dir_exp}/03_bam/bowtie2/qc/in_G1_Hho1_6336.preseq"

    mode="paired"
    mapq=1
    req_flag=true
    if [[ -n "${f_length}" ]]; then unset f_length; fi
    fastq_1="${dir_exp}/02_trim/IP_G1_Hho1_6336_R1.atria.fastq.gz"
    fastq_2="${dir_exp}/02_trim/IP_G1_Hho1_6336_R2.atria.fastq.gz"
    bam="${dir_exp}/03_bam/bowtie2/bam/IP_G1_Hho1_6336.bam"
    bam_coor="${dir_exp}/03_bam/bowtie2/bam/IP_G1_Hho1_6336.sort-coord.bam"
    bam_quer="${dir_exp}/03_bam/bowtie2/bam/IP_G1_Hho1_6336.sort-qname.bam"
    bed_siQ="${dir_exp}/03_bam/bowtie2/siQ-ChIP/IP_G1_Hho1_6336.bed.gz"
    bed_etc="${dir_exp}/03_bam/bowtie2/cvrg/IP_G1_Hho1_6336"
    d_fqc_f="${dir_exp}/03_bam/bowtie2/qc/fastqc/fastq"
    d_fqc_b="${dir_exp}/03_bam/bowtie2/qc/fastqc/bam"
    txt_list="${dir_exp}/03_bam/bowtie2/qc/IP_G1_Hho1_6336.list-tally.txt"
    txt_met="${dir_exp}/03_bam/bowtie2/qc/IP_G1_Hho1_6336.picard-metrics.txt"
    txt_flg="${dir_exp}/03_bam/bowtie2/qc/IP_G1_Hho1_6336.samtools-flagstat.txt"
    txt_idx="${dir_exp}/03_bam/bowtie2/qc/IP_G1_Hho1_6336.samtools-idxstats.txt"
    txt_pre="${dir_exp}/03_bam/bowtie2/qc/IP_G1_Hho1_6336.preseq"

    # mode="paired"
    # mapq=1
    # req_flag=true
    # if [[ -n "${f_length}" ]]; then unset f_length; fi
    # fastq_1="${dir_exp}/02_trim/IP_Q_Hmo1_7750_R1.atria.fastq.gz"
    # fastq_2="${dir_exp}/02_trim/IP_Q_Hmo1_7750_R2.atria.fastq.gz"
    # bam="${dir_exp}/03_bam/bowtie2/bam/IP_Q_Hmo1_7750.bam"
    # bam_coor="${dir_exp}/03_bam/bowtie2/bam/IP_Q_Hmo1_7750.sort-coord.bam"
    # bam_quer="${dir_exp}/03_bam/bowtie2/bam/IP_Q_Hmo1_7750.sort-qname.bam"
    # bed_siQ="${dir_exp}/03_bam/bowtie2/siQ-ChIP/IP_Q_Hmo1_7750.bed.gz"
    # bed_etc="${dir_exp}/03_bam/bowtie2/cvrg/IP_Q_Hmo1_7750"
    # d_fqc_f="${dir_exp}/03_bam/bowtie2/qc/fastqc/fastq"
    # d_fqc_b="${dir_exp}/03_bam/bowtie2/qc/fastqc/bam"
    # txt_list="${dir_exp}/03_bam/bowtie2/qc/IP_Q_Hmo1_7750.list-tally.txt"
    # txt_met="${dir_exp}/03_bam/bowtie2/qc/IP_Q_Hmo1_7750.picard-metrics.txt"
    # txt_flg="${dir_exp}/03_bam/bowtie2/qc/IP_Q_Hmo1_7750.samtools-flagstat.txt"
    # txt_idx="${dir_exp}/03_bam/bowtie2/qc/IP_Q_Hmo1_7750.samtools-idxstats.txt"
    # txt_pre="${dir_exp}/03_bam/bowtie2/qc/IP_Q_Hmo1_7750.preseq"

    # mode="single"
    # mapq=1
    # f_length=200
    # fastq_1="${dir_exp}/01_sym/IP_Q_Brn1_rep1.fastq.gz"
    # if [[ -n "${fastq_2}" ]]; then unset fastq_2; fi
    # bam="${dir_exp}/03_bam/bowtie2/bam/IP_Q_Brn1_rep1.bam"
    # bam_coor="${dir_exp}/03_bam/bowtie2/bam/IP_Q_Brn1_rep1.sort-coord.bam"
    # bam_quer="${dir_exp}/03_bam/bowtie2/bam/IP_Q_Brn1_rep1.sort-qname.bam"
    # bed_siQ="${dir_exp}/03_bam/bowtie2/siQ-ChIP/IP_Q_Brn1_rep1.bed.gz"
    # bed_etc="${dir_exp}/03_bam/bowtie2/cvrg/IP_Q_Brn1_rep1"
    # d_fqc_f="${dir_exp}/03_bam/bowtie2/qc/fastqc/fastq"
    # d_fqc_b="${dir_exp}/03_bam/bowtie2/qc/fastqc/bam"
    # txt_list="${dir_exp}/03_bam/bowtie2/qc/IP_Q_Brn1_rep1.list-tally.txt"
    # txt_met="${dir_exp}/03_bam/bowtie2/qc/IP_Q_Brn1_rep1.picard-metrics.txt"
    # txt_flg="${dir_exp}/03_bam/bowtie2/qc/IP_Q_Brn1_rep1.samtools-flagstat.txt"
    # txt_idx="${dir_exp}/03_bam/bowtie2/qc/IP_Q_Brn1_rep1.samtools-idxstats.txt"
    # txt_pre="${dir_exp}/03_bam/bowtie2/qc/IP_Q_Brn1_rep1.preseq"

    #  Optional: Check variable assignments 
    if ${check_variables}; then
        echo "
        ### Variable assignments for interactive testing ###

        threads=${threads}
        index=${index}
        fasta=${fasta}
        sizes=${sizes}
        mode=${mode}
        mapq=${mapq}
        req_flag=${req_flag}  # If mode=paired, then 'true' is filter for pairs, 'false' is do not.
        f_length=${f_length}  # If mode=paired, then empty or 0 means use TLEN.
        fastq_1=${fastq_1}
        fastq_2=${fastq_2}
        bam=${bam}
        bam_coor=${bam_coor}
        bam_quer=${bam_quer}
        bed_siQ=${bed_siQ}
        bed_etc=${bed_etc}
        d_fqc_f=${d_fqc_f}
        d_fqc_b=${d_fqc_b}
        txt_list=${txt_list}
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
             -rf|--req_flag) req_flag=true;   shift 1 ;;
              -m|--mode)     mode="${2}";     shift 2 ;;
             -f1|--fastq_1)  fastq_1="${2}";  shift 2 ;;
             -f2|--fastq_2)  fastq_2="${2}";  shift 2 ;;
              -b|--bam)      bam="${2}";      shift 2 ;;
             -bc|--bam_coor) bam_coor="${2}"; shift 2 ;;
             -bq|--bam_quer) bam_quer="${2}"; shift 2 ;;
             -bs|--bed_siQ)  bed_siQ="${2}";  shift 2 ;;
             -be|--bed_etc)  bed_etc="${2}";  shift 2 ;;
            -dff|--d_fqc_f)  d_fqc_f="${2}" ; shift 2 ;;
            -dfb|--d_fqc_b)  d_fqc_b="${2}" ; shift 2 ;;
             -tl|--txt_list) txt_list="${2}"; shift 2 ;;
             -tm|--txt_met)  txt_met="${2}";  shift 2 ;;
             -tf|--txt_flg)  txt_flg="${2}";  shift 2 ;;
             -ti|--txt_idx)  txt_idx="${2}";  shift 2 ;;
             -tp|--txt_pre)  txt_pre="${2}";  shift 2 ;;
            *) echo "Unknown parameter passed: ${1}"; exit 1 ;;
        esac
    done
fi

if [[ -z "${threads}" ]]; then
    threads="${SLURM_CPUS_ON_NODE:-1}"
    ${verbose} && echo_note "The number of threads has not been specified; defaulting to ${threads}."
fi

if ! [[ "${threads}" =~ ^[0-9]+$ ]] || [[ "${threads}" -lt 1 ]]; then
    error_and_exit "The number of threads must be an integer >= 1."
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
    ${verbose} && echo_note "MAPQ threshold has not been specified; defaulting to ${mapq}."
fi

if [[ -z "${req_flag}" || "${req_flag}" != "true" ]]; then
    req_flag=false
fi

if [[ -z "${mode}" ]]; then
    mode="paired"
    ${verbose} && echo_note "Mode has not been specified; defaulting to '${mode}' mode."
fi

case "${mode}" in
    single|paired) : ;;
    *) error_and_exit "Invalid mode '${mode}'. Must be either 'single' or 'paired'." ;;
esac

if [[ "${mode}" == "paired" && "${req_flag}" == "false" && ${verbose} ]]; then
    echo_note "Mode is 'paired' and BAM file will NOT be filtered for paired alignments (i.e., alignments for which flag bit 2 is set)."
else
    echo_note "Mode is 'paired' and BAM file will be filtered for paired alignments (i.e., alignments for which flag bit 2 is set)."
fi

if [[ -z "${f_length}" ]]; then
    if [[ "${mode}" == "single" ]]; then
        error_and_exit "Fragment length is required for 'single' mode. Use -fl or --f_length to specify it."
    else
        f_length=0
        ${verbose} && echo_note "Fragment length not specified for 'paired' mode; will determine from alignment (TLEN)."
    fi
elif ! [[ "${f_length}" =~ ^[0-9]+$ ]]; then
    error_and_exit "Fragment length must be a non-negative integer; received ${f_length}."
elif [[ "${f_length}" -eq 0 && "${mode}" == "single" ]]; then
    error_and_exit "Fragment length must be greater than 0 for 'single' mode."
elif [[ "${f_length}" -gt 0 && "${mode}" == "paired" && ${verbose} ]]; then
    echo_note "Fragment length specified as ${f_length} for 'paired' mode; will not determine from alignment (TLEN)."
elif [[ "${f_length}" -gt 0 && "${mode}" == "single" && ${verbose} ]]; then
    echo_note "Fragment length specified as ${f_length} for 'single' mode."
fi

if [[ -z "${fastq_1}" ]]; then
    error_and_exit "FASTQ_1 file path is required. Use -f1 or --fastq_1 to specify it."
fi

if [[ ! -f "${fastq_1}" ]]; then
    error_and_exit "Specified FASTQ_1 file does not exist."
fi

if [[ ! "${fastq_1}" =~ \.(fastq|fq)(\.gz)?$ ]]; then
    error_and_exit "FASTQ_1 does not have a valid extension, which must be one of the following: *.fastq, *.fastq.gz, *.fq, or *.fq.gz."
fi

if [[ "${mode}" == "paired" ]]; then
    if [[ -z "${fastq_2}" ]]; then
        error_and_exit "FASTQ_2 file path is required. Use -f2 or --fastq_2 to specify it."
    fi

    if [[ ! -f "${fastq_2}" ]]; then
        error_and_exit "Specified FASTQ_2 file does not exist."
    fi

    if [[ ! "${fastq_2}" =~ \.(fastq|fq)(\.gz)?$ ]]; then
        error_and_exit "FASTQ_2 does not have a valid extension, which must be one of the following: *.fastq, *.fastq.gz, *.fq, or *.fq.gz."
    fi
fi

if [[ "${mode}" == "single" && -n "${fastq_2}" && ${verbose} ]]; then
    echo_note "Mode is 'single' but FASTQ_2 file is assigned. FASTQ_2 will be ignored."
    echo "      fastq_2=\"${fastq_2}\""
fi

if [[ -z "${bam}" ]]; then
    error_and_exit "BAM file path is required. Use -b or --bam to specify it."
fi

if [[ -z "${bam_coor}" ]]; then
    bam_coor="${bam/.bam/.sort-coord.bam}"
    ${verbose} && echo_note "Coordinate-sorted BAM outfile has not been specified; defaulting to ${bam_coor}."
fi

if [[ -z "${bam_quer}" ]]; then
    bam_quer="${bam/.bam/.sort-qname.bam}"
    ${verbose} && echo_note "Queryname-sorted BAM outfile has not been specified; defaulting to ${bam_quer}."
fi

if [[ -z "${bed_siQ}" ]]; then
    bed_siQ="${bam/.bam/.siQ-ChIP.bed.gz}"
    ${verbose} && echo_note "siQ-ChIP BED outfile has not been specified; defaulting to ${bed_siQ}."
fi

if [[ "${bed_siQ}" != *.gz ]]; then
    bed_siQ="${bed_siQ}.gz"
    ${verbose} && echo_note "Appending suffix '.gz' to string assigned to bed_siQ: ${bed_siQ}" 
fi

if [[ -z "${bed_etc}" ]]; then
    bed_etc="${bam%.bam}"
    ${verbose} && echo_note "Stem for per-base coverage BED outfile, Mosdepth outfiles, etc. has not been specified; defaulting to ${bed_etc}."
fi

if [[ -z "${d_fqc_f}" && ${flag_fqc_fastq} ]]; then
    error_and_exit "FastQC directory is required. Use -dff or --d_fqc_f to specify it."
fi

if [[ ! -d "${d_fqc_f}" && ${flag_fqc_fastq} ]]; then
    error_and_exit "Specified directory for FastQC does not exist."
fi

if [[ -z "${d_fqc_b}" && ${flag_fqc_bam} ]]; then
    error_and_exit "FastQC directory is required. Use -dfb or --d_fqc_b to specify it."
fi

if [[ ! -d "${d_fqc_b}" && ${flag_fqc_bam} ]]; then
    error_and_exit "Specified directory for FastQC does not exist."
fi

if [[ -z "${txt_list}" ]]; then
    txt_idx="${bam/.bam/.list-tally.txt}"
    ${verbose} && echo_note "list_tally_flags_mapq outfile has not been specified; defaulting to ${txt_list}."
fi

if [[ -z "${txt_met}" ]]; then
    txt_idx="${bam/.bam/.picard-metrics.txt}"
    ${verbose} && echo_note "Picard CollectAlignmentSummaryMetrics outfile has not been specified; defaulting to ${txt_met}."
fi

if [[ -z "${txt_flg}" ]]; then
    txt_flg="${bam/.bam/.samtools-flagstat.txt}"
    ${verbose} && echo_note "Samtools flagstat outfile has not been specified; defaulting to ${txt_flg}."
fi

if [[ -z "${txt_idx}" ]]; then
    txt_idx="${bam/.bam/.samtools-idxstats.txt}"
    ${verbose} && echo_note "Samtools idxstats outfile has not been specified; defaulting to ${txt_idx}."
fi

if [[ -z "${txt_pre}" ]]; then
    txt_pre="${bam/.bam/.preseq}"
    ${verbose} && echo_note "Stem for preseq outfiles has not been specified; defaulting to ${txt_pre}."
fi


#  #TODO Name this section ====================================================
#  Step #0 ----------------------------
#+ Run FastQC on FASTQ file(s)
if ${flag_fqc_fastq}; then
    extension_1="$(determine_fastq_extension "${fastq_1}")" || {
        if ${interact}; then return 1; else exit 1; fi
    }
    file_fqc_1="${d_fqc_f}/$(basename "${fastq_1}" "${extension_1}")_fastqc.html"
    message_1="FastQC output already exists for FASTQ_1. Skipping FastQC operation (step #0)."
    
    if [[ "${mode}" == "single" ]]; then
        if [[ ! -f "${file_fqc_1}" ]]; then
            fastqc "${fastq_1}" -o "${d_fqc_f}"
        else
            echo_note "${message_1}"
        fi
    elif [[ "${mode}" == "paired" ]]; then
        extension_2="$(determine_fastq_extension "${fastq_2}")" || {
            if ${interact}; then return 1; else exit 1; fi
        }
        file_fqc_2="${d_fqc_f}/$(basename "${fastq_2}" "${extension_2}")_fastqc.html"
        message_2="FastQC output already exists for FASTQ_2. Skipping FastQC operation (step #0)."

        if [[ -f "${file_fqc_1}" || -f "${file_fqc_2}" ]]; then
            [[ -f "${file_fqc_1}" ]] && echo_note "${message_1}"
            [[ -f "${file_fqc_2}" ]] && echo_note "${message_2}"
        fi

        #  Execute FastQC if/as needed
        if [[ ! -f "${file_fqc_1}" && ! -f "${file_fqc_2}" ]]; then
            fastqc --threads 2 "${fastq_1}" "${fastq_2}" -o "${d_fqc_f}"
        elif [[ ! -f "${file_fqc_1}" && -f "${file_fqc_2}" ]]; then
            fastqc "${fastq_1}" -o "${d_fqc_f}"
        elif [[ -f "${file_fqc_1}" && ! -f "${file_fqc_2}" ]]; then
            fastqc "${fastq_2}" -o "${d_fqc_f}"
        fi
    fi
fi


#  Step #1 ----------------------------
#+ Check if the BAM files exist; if not, perform alignment with Bowtie 2,
#+ converting the Bowtie 2 output to a BAM file with Samtools; for paired-end
#+ sequenced data, retain only properly paired reads (-f 2) that are greater
#+ than or equal to a user-supplied MAPQ score (-q "${mapq}") [by default,
#+ retain any-quality "maxi-reads" (-q 1)]
#+ 
#+ On the definition of "maxi-reads" and how to interpret Bowtie 2 MAPQ scores:
#+ biofinysics.blogspot.com/2014/05/how-does-bowtie2-assign-mapq-scores.html
if [[
       ! -f "${bam}" \
    && ! -f "${bam_coor}" \
    && ! -f "${bam_quer}"
]]; then
    if [[ "${mode}" == "paired" ]]; then
        if ${req_flag}; then
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
        else
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
                        -q "${mapq}" \
                        -o "${bam}"
        fi
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
    echo_note "BAM file exists. Skipping alignment operations (step #1)."
fi


#  Step #2 ----------------------------
#+ Sort the BAM file by queryname after checking that the BAM infile exists and
#+ the queryname-sorted BAM outfile does not exist 
if [[ -f "${bam}" ]]; then
    if [[ ! -f "${bam_quer}" ]]; then
        if ${verbose}; then
            echo "Running step #2a: Sort the BAM file by querynames."
        fi

        samtools sort \
            -@ "${threads}" \
            -n \
            -o "${bam_quer}" \
            "${bam}"

        if [[ "${mode}" == "paired" ]]; then
            #  After sorting by queryname, fix the paired read mate
            #+ information; this information is required for subsequent
            #+ operations
            if ${verbose}; then
                echo "Running step #2b: Fix the aligned mate information in queryname-sorted BAM file."
            fi

            if [[ -f "${bam_quer}" ]]; then
                samtools fixmate \
                    -@ "${threads}" \
                    -c \
                    -m \
                    "${bam_quer}" \
                    "${bam_quer%.bam}.tmp.bam"

                #  Replace the original queryname-sorted BAM with
                #+ queryname-sorted mate-fixed BAM
                if [[ -f "${bam_quer%.bam}.tmp.bam" ]]; then
                    mv -f \
                        "${bam_quer%.bam}.tmp.bam" \
                        "${bam_quer}"
                fi
            fi
        elif [[ "${mode}" == "single" && ${verbose} ]]; then
            echo_note "Skipping step #2b (run samtools fixmate), which is not applicable for a BAM file derived from single-end sequencing data."
        fi
    else
        echo_note "Queryname-sorted BAM file exists. Skipping sort operations (step #2)."
    fi
else
    error_and_exit "BAM infile does not exist (step #2)."
fi


#  Step #3 ----------------------------
#+ For other downstream analyses, sort the queryname-sorted BAM by
#+ coordinates
if [[ -f "${bam_quer}" ]]; then
    if [[ ! -f "${bam_coor}" ]]; then
        if ${verbose}; then
            echo "Running step #3a: Sort the queryname-sorted BAM file by coordinates."
        fi

        samtools sort \
            -@ "${threads}" \
            -o "${bam_coor}" \
            "${bam_quer}"
    else
        echo_note "Coordinate-sorted BAM file exists. Skipping sort operation (step #3a)."
    fi
else
    error_and_exit "Queryname-sorted BAM infile does not exist (step #3a)."
fi

#  Index the coordinate-sorted BAM file
if [[ -f "${bam_coor}" ]]; then
    if [[ ! -f "${bam_coor}.bai" ]]; then
        if ${verbose}; then
            echo "Running step #3b: Index the coordinate-sorted BAM file."
        fi

        samtools index \
            -@ "${threads}" \
            "${bam_coor}"
    else
        echo_note "Index (BAI file) for coordinate-sorted BAM file exists. Skipping index operation (step #3b)."
    fi
else
    error_and_exit "Coordinate-sorted BAM infile does not exist (step #3b)."
fi


#  Step #4 ----------------------------
#+ Mark duplicate alignments in the coordinate-sorted BAM file
if [[ -f "${bam_coor}" ]]; then
    if [[ -f "${bam_coor}.bai" ]]; then
        if \
            samtools view -H "${bam_coor}" \
                | grep -q "@PG.*CL:samtools markdup";
        then
            echo_note "Duplicate alignments are already marked in coordinate-sorted BAM file. Skipping markdup operations (step #4)."
        else
            if ${verbose}; then
                echo "Running step #4: Mark duplicate alignments in the coordinate-sorted BAM file."
            fi

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
        fi
    else
        error_and_exit "Index (BAI file) for coordinate-sorted BAM file does not exist (step #4)."
    fi
else
    error_and_exit "Coordinate-sorted BAM file does not exist (step #4)."
fi


#  Step #5 ----------------------------
#+ Generate various quality-check files
if [[ -f "${bam_coor}" ]]; then
    if [[ -f "${bam_coor}.bai" ]]; then
        if \
            samtools view -H "${bam_coor}" \
               | grep -q "@PG.*CL:samtools markdup";
        then
            #  Step #5a -----
            file_fqc_bam="${d_fqc_b}/$(basename "${bam_coor}" .bam)_fastqc.html"
            if [[ ! -f "${file_fqc_bam}" && ${flag_fqc_bam} ]]; then
                ${verbose} \
                    && echo "Running step #5a: Run FastQC on coordinate-sorted BAM file." \
                    && echo "                  Output directory will be ${d_fqc_b}."
                fastqc "${bam_coor}" -o "${d_fqc_b}"
            else
                echo_note "FastQC output for coordinate-sorted BAM file already exists. Skipping FastQC operation (step #5a)."
            fi

            #  Step #5b -----
            if [[ ! -f "${txt_list}" && ${flag_list_tally} ]]; then
                ${verbose} \
                    && echo "Running step #5b: Run list_tally_flags_mapq on coordinate-sorted BAM file." \
                    && echo "                  Outfile will be ${txt_list}."
                list_tally_flags_mapq \
                    -t "${threads}" \
                    -i "${bam_coor}" \
                    -o "${txt_list}" \
                    -w "both" \
                    -b
            else
                echo_note "list_tally_flags_mapq TXT file already exists. Skipping list_tally_flags_mapq operation (step #5b)."
            fi

            #  Step #5c -----
            if [[ ! -f "${txt_met}" && ${flag_met} ]]; then
                #  Run picard CollectAlignmentSummaryMetrics to assess the
                #+ quality of read alignments as well as "the proportion of
                #+ reads that passed machine signal-to-noise threshold quality
                #+ filters"
                ${verbose} \
                    && echo "Running step #5c: Run picard CollectAlignmentSummaryMetrics on coordinate-sorted BAM file." \
                    && echo "                  Outfile will be ${txt_met}."

                if [[ "${fasta##*.}" == "gz" ]]; then
                    decomp_fasta="${fasta%.gz}"
                    if [[ ! -f "${decomp_fasta}.fai" ]]; then
                        #  Picard CollectAlignmentSummaryMetrics requires that
                        #+ the FASTA file is not compressed
                        ${verbose} && echo_note "Decompressing the reference genome FASTA file."

                        if [[ ! -f "${decomp_fasta}" ]]; then
                            gunzip -c "${fasta}" > "${decomp_fasta}"
                        fi

                        #  Also, Picard CollectAlignmentSummaryMetrics requires
                        #+ that the reference genome FASTA is indexed; so,
                        #+ generate an FAI file if necessary
                        ${verbose} && echo_note "Generating an index file (FAI) for the reference genome FASTA file."
                        samtools faidx -@ "${threads}" "${decomp_fasta}"
                    fi

                    picard CollectAlignmentSummaryMetrics \
                        --REFERENCE_SEQUENCE "${decomp_fasta}" \
                        --INPUT "${bam_coor}" \
                        --OUTPUT "${txt_met}"
                else
                    if [[ ! -f "${fasta}.fai" ]]; then
                        ${verbose} && echo_note "Generating an index file (FAI) for the reference genome FASTA file."

                        samtools faidx -@ "${threads}" "${fasta}"
                    fi

                    picard CollectAlignmentSummaryMetrics \
                        --REFERENCE_SEQUENCE "${fasta}" \
                        --INPUT "${bam_coor}" \
                        --OUTPUT "${txt_met}"
                fi
            else
                echo_note "picard CollectAlignmentSummaryMetrics TXT file already exists. Skipping picard CollectAlignmentSummaryMetrics operations (step #5c)."
            fi

            #  Step #5d -----
            if [[ ! -f "${txt_flg}" && ${flag_flg} ]]; then
                #  Generate Samtools flagstat report
                ${verbose} \
                    && echo "Running step #5d: Run samtools flagstat on coordinate-sorted BAM file." \
                    && echo "                  Outfile will be ${txt_flg}."

                samtools flagstat \
                    -@ "${threads}" \
                    "${bam_coor}" \
                        > "${txt_flg}"
            else
                echo_note "samtools flagstat TXT file already exists. Skipping samtools flagstat operation (step #5d)."
            fi

            #  Step #5e -----
            if [[ ! -f "${txt_idx}" && ${flag_idx} ]]; then
                #  Generate Samtools idxstats report
                ${verbose} \
                    && echo "Running step #5e: Run samtools idxstats on coordinate-sorted BAM file." \
                    && echo "                  Outfile will be ${txt_idx}."

                samtools idxstats "${bam_coor}" > "${txt_idx}"
            else
                echo_note "samtools idxstats TXT file already exists. Skipping samtools idxstats operation (step #5e)."
            fi

            #  Step #5f -----
            if [[ ! -f "${txt_pre}-lc-extrap.txt" && ${flag_lc_extrap} ]]; then
                #  Use preseq lc_extrap to generate the expected yield for theoretical
                #+ larger experiments; preseq lc_extrap also gives bounds on the number of
                #+ distinct reads in the library and the associated confidence intervals
                ${verbose} \
                    && echo "Running step #5f: Run preseq lc_extrap on coordinate-sorted BAM file." \
                    && echo "                  Outfile #1 will be ${txt_pre}-lc-extrap.txt." \
                    && echo "                  Outfile #2 will be ${txt_pre%.txt}-lc-extrap-verbose.txt."

                preseq lc_extrap -v -P -B -r 24 \
                    -o "${txt_pre}-lc-extrap.txt" \
                    "${bam_coor}" \
                        &> "${txt_pre%.txt}-lc-extrap-verbose.txt"
            else
                echo_note "preseq lc_extrap TXT files already exist. Skipping lc_extrap operation (step #5f)."
            fi

            #  Step #5g -----
            if [[ ! -f "${txt_pre}-bound-pop.txt" && ${flag_bound_pop} ]]; then
                #TODO Use preseq bound_pop to...
                ${verbose} \
                    && echo "Running step #5g: Run preseq bound_pop on coordinate-sorted BAM file." \
                    && echo "                  Outfile will be ${txt_pre}-bound-pop.txt."

                preseq bound_pop -P -B -r 24 \
                    -o "${txt_pre}-bound-pop.txt" \
                    "${bam_coor}"
            else
                echo_note "preseq bound_pop TXT files already exist. Skipping bound_pop operation (step #5g)."
            fi

            #  Step #5h -----
            if [[ ! -f "${txt_pre}-c-curve.txt" && ${flag_c_curve} ]]; then
                #TODO Use preseq c_curve to...
                ${verbose} \
                    && echo "Running step #5h: Run preseq c_curve on coordinate-sorted BAM file." \
                    && echo "                  Outfile will be ${txt_pre}-c-curve.txt."

                preseq c_curve -P -B -r 24 \
                    -o "${txt_pre}-c-curve.txt" \
                    "${bam_coor}"
            else
                echo_note "preseq c_curve TXT files already exist. Skipping c_curve operation (step #5h)."
            fi
        else
            error_and_exit "Duplicate alignments have not been marked in coordinate-sorted BAM file (step #5)."
        fi
    else
        error_and_exit "Index (BAI file) for coordinate-sorted BAM file does not exist (step #5)."
    fi
else
    error_and_exit "Coordinate-sorted BAM file does not exist (step #5)."
fi


#  Step #6 ----------------------------
if ${flag_siq} && [[ -f "${bam_quer}" ]]; then
    if [[ ! -f "${bed_siQ}" ]]; then
        if [[ "${mode}" == "paired" ]]; then
            if [[ "${f_length}" -eq 0 ]]; then
                samtools view -@ "${threads}" "${bam_quer}" \
                    | awk '{
                        if (NR % 2 == 1) {
                            chr_1 = $3; 
                            start_1 = $4; 
                            len_1 = length($10);
                        } else {
                            chr_2 = $3;
                            start_2 = $4; 
                            len_2 = length($10);

                            if (chr_1 == chr_2) {
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
            elif [[ "${f_length}" -gt 0 ]]; then
                samtools view -@ "${threads}" "${bam_quer}" \
                    | awk \
                        -v f_length="${f_length}" \
                        '{
                            #  Work with only first-in-pair alignments
                            if (and($2, 64)) {
                                chr = $3;
                                start = $4;
                                strand = and($2, 16);

                                if (strand == 16) {  # Reverse strand
                                    end = start;
                                    start = end - f_length + 1;
                                } else {  # Forward strand
                                    end = start + f_length - 1;
                                }

                                if (start < 1) start = 1;
                                frag_length = end - start + 1;
                                
                                print chr, start, end, frag_length;
                            }
                        }' \
                        OFS='\t' \
                    | sort -k1,1 -k2,2n \
                    | gzip \
                        > "${bed_siQ}"
            fi
        elif [[ "${mode}" == "single" ]]; then
            samtools view -@ "${threads}" "${bam_quer}" \
                | awk \
                    -v f_length="${f_length}" \
                    '{
                        chr = $3;
                        start = $4;
                        strand = and($2, 16);

                        if (strand == 16) {  # Reverse strand
                            end = start;
                            start = end - f_length + 1;
                        } else {  # Forward strand
                            end = start + f_length - 1;
                        }

                        if (start < 1) start = 1;
                        frag_length = end - start + 1;

                        print chr, start, end, frag_length;
                    }' \
                    OFS='\t' \
                | sort -k1,1 -k2,2n \
                | gzip \
                    > "${bed_siQ}"
        else
            error_and_exit "Invalid mode specified: ${mode}; expected 'single' or 'paired' (step #6)."
        fi
    else
        echo_note "BED file for siQ-ChIP already exists. Skipping operations (step #6)."
    fi
else
    error_and_exit "Queryname-sorted BAM file does not exist (step #6)."
fi


#  Step #7 ----------------------------
#+ Using the coordinate-sorted BAM file, generate a BED file of per-base
#+ coverage as well as RPM-scaled BEDGRAPH and BIGWIG files
if ${flag_etc} && [[ -f "${bam_coor}" ]]; then
    if [[ -f "${bam_coor}.bai" ]]; then
        #  Step #7a ---------
        if [[ ! -f "${bed_etc}.per-base.bed.gz" ]]; then
            #  Mosdepth throws errors if BAI is older than BAM
            if [[ "${bam_coor}.bai" -ot "${bam_coor}" ]]; then
                #  If BAI is older than BAM, then index the BAM again (current
                #+ BAI will be overwritten)
                samtools index \
                    -@ "${threads}" \
                    "${bam_coor}"
            fi

            mosdepth \
                "${bed_etc}" \
                "${bam_coor}"
        else
            echo_note "Per-base alignment coverage BED file exists. Skipping per-base alignment coverage operations (step #7a)."
        fi

        # Step #7b ---------
        if [[ -f "${bed_etc}.per-base.bed.gz" ]]; then
            if [[
                   ! -f "${bed_etc}.rpm.bedgraph" \
                || ! -f "${bed_etc}.rpm.bedgraph.gz"
            ]]; then
                # Step #7b1 -----
                total_reads="$(samtools view -c "${bam_coor}")"

                if [[ "${total_reads}" -eq 0 ]]; then
                    error_and_exit "In calculation for RPM alignment coverage, total_reads is ${total_reads} (step #7b1)."
                fi

                zcat "${bed_etc}.per-base.bed.gz" \
                    | awk -v total_reads="${total_reads}" '
                        BEGIN { OFS="\t" }
                        {
                            print $1, $2, $3, ($4 / total_reads) * 1000000
                        }' \
                            > "${bed_etc}.rpm.bedgraph"

                # Step #7b2 -----
                if [[ -f "${bed_etc}.rpm.bedgraph" ]]; then
                    if [[ ! -f "${bed_etc}.rpm.sort.bg" ]]; then
                        bedSort \
                            "${bed_etc}.rpm.bedgraph" \
                            "${bed_etc}.rpm.sort.bg"
                    else
                        echo_note "Sorted RPM alignment coverage BEDGRAPH file already exists. Skipping sorting operation (step #7b2)."
                    fi

                    mv -f \
                        "${bed_etc}.rpm.sort.bg" \
                        "${bed_etc}.rpm.bedgraph"
                else
                    error_and_exit "RPM alignment coverage BEDGRAPH file does not exist (step #7b2)."
                fi

                # Step #7b3 -----
                if [[ -f "${bed_etc}.rpm.bedgraph" ]]; then
                    if [[ ! -f "${bed_etc}.rpm.bigwig" ]]; then
                        bedGraphToBigWig "${bed_etc}.rpm.bedgraph" "${sizes}" "${bed_etc}.rpm.bigwig"
                    else
                        echo_note "RPM alignment coverage BIGWIG file already exists. Skipping sorting operation (step #7b3)."
                    fi
                else
                    error_and_exit "RPM alignment coverage BEDGRAPH file does not exist (step #7b3)."
                fi

                # Step #7b4 -----
                # Compress the BEDGRAPH file if both BEDGRAPH and BIGWIG files exist
                if [[
                       -f "${bed_etc}.rpm.bedgraph"
                    && -f "${bed_etc}.rpm.bigwig"
                ]]; then
                    if [[ ! -f "${bed_etc}.rpm.bedgraph.gz" ]]; then
                        gzip "${bed_etc}.rpm.bedgraph"
                    else
                        echo_note "Compressed RPM alignment coverage BEDGRAPH file already exists. Skipping compression (step #7b4)."
                    fi
                else
                    error_and_exit "Required RPM alignment coverage files do not exist (step #7b4)."
                fi
            else
                echo_note "RPM alignment coverage BEDGRAPH file exists. Skipping RPM alignment coverage operations (step #7b)."
            fi
        else
            error_and_exit "Per-base alignment coverage BED file does not exist (step #7b)."
        fi
    else
        error_and_exit "Index (BAI file) for coordinate-sorted BAM file does not exist (step #7)."
    fi
else
    error_and_exit "Coordinate-sorted BAM infile does not exist (step #7)."
fi


#  Step #8 ----------------------------
#+ Remove the original BAM file (Bowtie 2 output piped to Samtools) if all
#+ other files have been successfully created
if ${flag_clean}; then
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
fi
