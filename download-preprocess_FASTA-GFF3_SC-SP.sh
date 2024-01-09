#!/bin/bash

#  download-preprocess_FASTA-GFF3_SC-SP.sh
#  KA


#  Define functions ===========================================================
#  Function to exit with exit code 1, which stops the execution of all code,
#+ and return an error message
function error_and_exit() {
    echo "Error: ${1}" >&2
    exit 1
}


#  Function to return an error message with exit code 1, which stops the
#+ execution of the current code
function error_and_return() {
    echo "Error: ${1}" >&2
    return 1
}


#  Function to download a file
function download_file() {
    local url="${1}"
    local output_path="${2}"

    # if ! wget -q -O "${output_path}" "${url}"; then
    # if ! curl -s -o "${output_path}" "${url}"; then
    if ! curl -o "${output_path}" "${url}"; then
        error_and_exit "Failed to download ${url}; exiting."
    fi
}


#  Function to download and extract a tarball
function download_extract_tarball() {
    local url="${1}"
    local output_dir="${2}"

    # if ! wget -q -O - "${url}" | tar -xz -C "${output_dir}"; then
    # if ! curl -s "${url}" | tar -xz -C "${output_dir}"; then
    if ! curl "${url}" | tar -xz -C "${output_dir}"; then
        error_and_exit "Failed to download and extract tarball from ${url}; exiting."
    fi
}


#  Function to download and extract tarball of S. cerevisiae FASTA, GFF3, and
#+ SGD files
download_extract_SC() {
    local tarball_path="${dir_genomes}/${dir_SC}/${tarball_SC}"

    if [[ -d "${tarball_path%.tgz}" ]]; then
        error_and_return "Path ${tarball_path%.tgz} exists; skipping download and extraction of S. cerevisiae files."
    else
        download_extract_tarball \
            "${URL_SC}/${tarball_SC}" \
            "${dir_genomes}/${dir_SC}"
    fi
}


organize_SC_files() {
    local tarball_path="${dir_genomes}/${dir_SC}/${tarball_SC%.tgz}"

    if [[ -d "${tarball_path}" ]]; then
        #TODO Checks for the following
        mv "${tarball_path}/"*.sgd.gz "${dir_genomes}/${dir_SC}/sgd"
        mv "${tarball_path}/"*.gff.gz "${dir_genomes}/${dir_SC}/gff3"
        mv "${tarball_path}/"*.{fasta,fsa}.gz "${dir_genomes}/${dir_SC}/fasta"

        if find \
            "${tarball_path}" \
            -mindepth 1 \
            -maxdepth 1 \
            -print \
            -quit \
                | grep -q .
        then
            error_and_exit "Path ${tarball_path} is not empty; stopping directory removal and exiting."
        else
            rmdir "${tarball_path}"
        fi
    fi
}


download_SP_files() {
    local file_type="${1}"
    local array=("${!2}")
    local url_base="${3}"

    for file in "${array[@]}"; do
        local url="${url_base}/${file}"
        local output_file="${dir_genomes}/${dir_SP}/${file_type}/${file}"

        if [[ -f "${output_file}" ]]; then
            error_and_return "${output_file} already exists; skipping download of S. pombe files."
        else
            download_file "${url}" "${output_file}"
        fi
    done
}


#  Function to prepare S. cerevisiae FASTA for concatenation with S. pombe
#+ FASTA
process_SC_fasta() {
    local source_fasta="${dir_genomes}/${dir_SC}/fasta/S288C_reference_sequence_R64-3-1_20210421.fsa.gz"
    local processed_fasta="${dir_genomes}/${dir_SC}/fasta-processed/S288C_reference_sequence_R64-3-1_20210421.fa.gz"
    
    if [[ -f "${processed_fasta}" ]]; then
        error_and_return "${processed_fasta} already exists; skipping processing of S. cerevisiae FASTA file."
    elif [[ ! -f "${source_fasta}" ]]; then
        error_and_exit "${source_fasta} not found; exiting."
    else
        zcat "${source_fasta}" \
            | sed -r 's/^.*1133.*$/>I/g' \
            | sed -r 's/^.*1134.*$/>II/g' \
            | sed -r 's/^.*1135.*$/>III/g' \
            | sed -r 's/^.*1136.*$/>IV/g' \
            | sed -r 's/^.*1137.*$/>V/g' \
            | sed -r 's/^.*1138.*$/>VI/g' \
            | sed -r 's/^.*1139.*$/>VII/g' \
            | sed -r 's/^.*1140.*$/>VIII/g' \
            | sed -r 's/^.*1141.*$/>IX/g' \
            | sed -r 's/^.*1142.*$/>X/g' \
            | sed -r 's/^.*1143.*$/>XI/g' \
            | sed -r 's/^.*1144.*$/>XII/g' \
            | sed -r 's/^.*1145.*$/>XIII/g' \
            | sed -r 's/^.*1146.*$/>XIV/g' \
            | sed -r 's/^.*1147.*$/>XV/g' \
            | sed -r 's/^.*1148.*$/>XVI/g' \
            | sed -r 's/^.*1224.*$/>Mito/g' \
            | gzip \
                > "${processed_fasta}"
    fi
}


#  Function to prepare S. cerevisiae GFF3 for concatenation with S. pombe GFF3
process_SC_gff3() {
    local source_gff3="${dir_genomes}/${dir_SC}/gff3/saccharomyces_cerevisiae_R64-3-1_20210421.gff.gz"
    local processed_gff3="${dir_genomes}/${dir_SC}/gff3-processed/saccharomyces_cerevisiae_R64-3-1_20210421.gff3.gz"
    
    #  Remove "chr" prefixes and fasta sequences from SGD R64-3-1 S. cerevisiae
    #+ GFF3 file
    if [[ -f "${processed_gff3}" ]]; then
        error_and_return "${processed_gff3} already exists; skipping processing of S. cerevisiae GFF3 file."
    elif [[ ! -f "${source_gff3}" ]]; then
        error_and_exit "${source_gff3} not found; exiting."
    else
        zcat "${source_gff3}" \
            | sed 's/^chr//g' \
            | sed 's/^mt/Mito/g' \
            | sed '/^###/,$d' \
            | gzip \
                > "${processed_gff3}"
    fi
}


#  Function to prepare S. pombe FASTA for concatenation with S. cerevesiae
#+ FASTA
process_SP_fasta() {
    local source_fasta="${dir_genomes}/${dir_SP}/fasta/Schizosaccharomyces_pombe_all_chromosomes.fa.gz"
    local processed_fasta="${dir_genomes}/${dir_SP}/fasta-processed/Schizosaccharomyces_pombe_all_chromosomes.fa.gz"
    
    #  Append "SP_" prefix to S. pombe chromosome names
    if [[ -f "${processed_fasta}" ]]; then
        error_and_return "${processed_fasta} already exists; skipping processing of S. pombe FASTA file."
    elif [[ ! -f "${source_fasta}" ]]; then
        error_and_exit "${source_fasta} not found; exiting."
    else
        zcat "${source_fasta}" \
            | sed -r 's/^>chr_II_telomeric_gap\ .*$/>SP_II_TG/g' \
            | sed -r 's/^>I\ .*$/>SP_I/g' \
            | sed -r 's/^>II\ .*$/>SP_II/g' \
            | sed -r 's/^>III\ .*$/>SP_III/g' \
            | sed -r 's/^>mating_type_region\ .*$/>SP_MTR/g' \
            | sed -r 's/^>mitochondrial\ .*$/>SP_Mito/g' \
            | gzip \
                > "${processed_fasta}"
    fi
}


#  Function to prepare S. pombe GFF3 for concatenation with S. cerevesiae GFF3
process_SP_gff3() {
    local source_gff3="${dir_genomes}/${dir_SP}/gff3/Schizosaccharomyces_pombe_all_chromosomes.gff3.gz"
    local processed_gff3="${dir_genomes}/${dir_SP}/gff3-processed/Schizosaccharomyces_pombe_all_chromosomes.gff3.gz"

    #  Append "SP_" prefix to S. pombe chromosome names
    if [[ -f "${processed_gff3}" ]]; then
        error_and_return "${processed_gff3} already exists; skipping processing of S. pombe GFF3 file."
    elif [[ ! -f "${source_gff3}" ]]; then
        error_and_exit "${source_gff3} not found; exiting."
    else
        zcat "${source_gff3}" \
            | sed 's/^chr_II_telomeric_gap/SP_II_TG/g' \
            | sed 's/^I/SP_I/g' \
            | sed 's/^II/SP_II/g' \
            | sed 's/^III/SP_III/g' \
            | sed 's/^mating_type_region/SP_MTR/g' \
            | sed 's/^mitochondrial/SP_Mito/g' \
            | gzip \
                > "${processed_gff3}"
    fi
}


#  Initialize variables and arrays ============================================
#  Initialize variables
dir_genomes="${HOME}/genomes"
string_SC="Saccharomyces_cerevisiae"
string_SP="Schizosaccharomyces_pombe"
dir_SC="${string_SC}"
dir_SP="${string_SP}"

URL_SC="http://sgd-archive.yeastgenome.org/sequence/S288C_reference/genome_releases"
tarball_SC="S288C_reference_genome_R64-3-1_20210421.tgz"

URL_SP_fasta="https://www.pombase.org/data/genome_sequence_and_features/genome_sequence"
URL_SP_gff3="https://www.pombase.org/data/genome_sequence_and_features/gff3"

#  Initialize arrays
# shellcheck disable=SC2034
unset fasta_SP && typeset -a fasta_SP=(
    "${string_SP}_all_chromosomes.fa.gz"
    "${string_SP}_chr_II_telomeric_gap.fa.gz"
    "${string_SP}_chromosome_I.fa.gz"
    "${string_SP}_chromosome_II.fa.gz"
    "${string_SP}_chromosome_III.fa.gz"
    "${string_SP}_mating_type_region.fa.gz"
    "${string_SP}_mitochondrial_chromosome.fa.gz"
)

# shellcheck disable=SC2034
unset gff3_SP && typeset -a gff3_SP=(
    "${string_SP}_all_chromosomes.gff3.gz"
    "${string_SP}_chr_II_telomeric_gap.gff3.gz"
    "${string_SP}_chromosome_I.gff3.gz"
    "${string_SP}_chromosome_II.gff3.gz"
    "${string_SP}_chromosome_III.gff3.gz"
    "${string_SP}_mating_type_region.gff3.gz"
    "${string_SP}_mitochondrial_chromosome.gff3.gz"
)


#  Run main work ==============================================================
#  Handle S. cerevisiae files -------------------------------------------------
#  Create directories for storing essential FASTA and GFF3 files
mkdir -p \
    "${dir_genomes}/${dir_SC}/"{err_out,fasta,fasta-processed,gff3,gff3-processed,sgd}

#  Download and organize S. cerevisiae FASTA and GFF3 files
download_extract_SC  #TODO Print some output; remove curl -s option?
organize_SC_files

#  Prepare S. cerevisiae FASTA for concatenation with S. pombe FASTA
process_SC_fasta  #TODO Print some output

#  Prepare S. cerevisiae GFF3 for concatenation with S. pombe GFF3
process_SC_gff3  #TODO Print some output


#  Handle S. pombe files ------------------------------------------------------
#  Create directories for storing essential FASTA and GFF3 files
mkdir -p \
    "${dir_genomes}/${dir_SP}/"{err_out,fasta,fasta-processed,gff3,gff3-processed}

#  Download and organize S. pombe FASTA and GFF3 files
download_SP_files "fasta" "fasta_SP[@]" "${URL_SP_fasta}"  #TODO Print some output
download_SP_files "gff3" "gff3_SP[@]" "${URL_SP_gff3}"  #TODO Print some output

#  Prepare S. pombe FASTA for concatenation with S. cerevesiae FASTA
process_SP_fasta  #TODO Print some output

#  Prepare S. pombe GFF3 for concatenation with S. cerevesiae GFF3
process_SP_gff3  #TODO Print some output
