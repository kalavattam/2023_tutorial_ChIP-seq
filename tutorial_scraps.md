
<a id="3-align-trimmed-fastq-files"></a>
# 3. Align trimmed FASTQ files
<a id="a-generate-a-concatenated-annotated-assembly-of-the-s-cerevisiae-and-s-pombe-genomes"></a>
## a. Generate a concatenated annotated assembly of the *S. cerevisiae* and *S. pombe* genomes
<a id="code-3"></a>
### Code
<details>
<summary><i>Code: 3.a. Generate a concatenated annotated assembly of the *S. cerevisiae* and *S. pombe* genomes</i></summary>

```bash
#!/bin/bash

#  Define functions ===========================================================
#  Function to return an error message and exit code 1, which stops the
#+ interactive execution of code
function error_and_return() {
    echo "Error: ${1}" >&2
    return 1
}


#  Function to download a file
function download_file() {
    local url="${1}"
    local output_path="${2}"

    # if ! wget -q -O "${output_path}" "${url}"; then
    if ! curl -s -o "${output_path}" "${url}"; then
        error_and_return "Failed to download ${url}"
    fi
}
export -f download_file


#  Function to download and extract a tarball
function download_extract_tarball() {
    local url="${1}"
    local output_dir="${2}"

    # if ! wget -q -O - "${url}" | tar -xz -C "${output_dir}"; then
    if ! curl -s "${url}" | tar -xz -C "${output_dir}"; then
        error_and_return "Failed to download and extract tarball from ${url}"
    fi
}
export -f download_extract_tarball


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
unset fasta_SP && typeset -a fasta_SP=(
    ${string_SP}_all_chromosomes.fa.gz
    ${string_SP}_chr_II_telomeric_gap.fa.gz
    ${string_SP}_chromosome_I.fa.gz
    ${string_SP}_chromosome_II.fa.gz
    ${string_SP}_chromosome_III.fa.gz
    ${string_SP}_mating_type_region.fa.gz
    ${string_SP}_mitochondrial_chromosome.fa.gz
)

URL_SP_gff3="https://www.pombase.org/data/genome_sequence_and_features/gff3"
unset gff3_SP && typeset -a gff3_SP=(
    ${string_SP}_all_chromosomes.gff3.gz
    ${string_SP}_chr_II_telomeric_gap.gff3.gz
    ${string_SP}_chromosome_I.gff3.gz
    ${string_SP}_chromosome_II.gff3.gz
    ${string_SP}_chromosome_III.gff3.gz
    ${string_SP}_mating_type_region.gff3.gz
    ${string_SP}_mitochondrial_chromosome.gff3.gz
)

#  Time for jobs submitted to SLURM
time="4:00:00"  # Adjust as needed


#  Do the main work ===========================================================
#TODO Modularize the below code

#  Create directories for storing essential FASTA and GFF3 files --------------
if [[ ! -d "${dir_genomes}/${dir_SC}" ]]; then
    mkdir -p ${dir_genomes}/${dir_SC}/{err_out,fasta,fasta-processed,gff3,gff3-processed,sgd}
fi

if [[ ! -d "${dir_genomes}/${dir_SP}" ]]; then
    mkdir -p ${dir_genomes}/${dir_SP}/{err_out,fasta,fasta-processed,gff3,gff3-processed}
fi


#  Check variables and arrays (optional) --------------------------------------
check_variables=true   # Check variable assignments
check_arrays=true      # Check array assignments

if ${check_variables}; then
    echo "
    dir_genomes=\"${dir_genomes}\"
    string_SC=\"${string_SC}\"
    string_SP=\"${string_SP}\"
    dir_SC=\"${dir_SC}\"
    dir_SP=\"${dir_SP}\"

    URL_SC=\"${URL_SC}\"
    tarball_SC=\"${tarball_SC}\"

    URL_SP_fasta=\"${URL_SP_fasta}\"
    URL_SP_gff3=\"${URL_SP_gff3}\"
    "
fi

if ${check_arrays}; then
    for element in "${fasta_SP[@]}" "${gff3_SP[@]}"; do
        echo "${element}"
    done
fi


#  Download and store Saccharomyces cerevisiae FASTA and GFF3 files -----------
#  Set flags
check_operations=true  # Check operations to download genome files
run_operations=true    # Run operations to download genome files

#  Echo the download logic if check_operations is true
if ${check_operations}; then
    if [[ ! -d "${dir_genomes}/${dir_SC}/${tarball_SC%.tgz}" ]]; then
        job_name="download-extract.${tarball_SC%.tgz}"

        echo "
        sbatch << EOF
#!/bin/bash

#SBATCH --job-name=\"${job_name}\"
#SBATCH --nodes=1
#SBATCH --time=${time}
#SBATCH --error=\"${dir_genomes}/${dir_SC}/err_out/${job_name}.%A.stderr.txt\"
#SBATCH --output=\"${dir_genomes}/${dir_SC}/err_out/${job_name}.%A.stdout.txt\"

# Extract tarball command
download_extract_tarball \\
    \"${URL_SC}/${tarball_SC}\" \\
    \"${dir_genomes}/${dir_SC}\"
EOF
        "
    else
        echo "Warning: Path ${dir_genomes}/${dir_SC}/${tarball_SC%.tgz} exist; skipping download and extraction."
    fi
fi

#  Download and extract the tarball for Saccharomyces cerevisiae
if ${run_operations}; then
    if [[ ! -d "${dir_genomes}/${dir_SC}/${tarball_SC%.tgz}" ]]; then
        job_name="download-extract.${tarball_SC%.tgz}"

        sbatch << EOF
#!/bin/bash

#SBATCH --job-name="${job_name}"
#SBATCH --nodes=1
#SBATCH --time=${time}
#SBATCH --error="${dir_genomes}/${dir_SC}/err_out/${job_name}.%A.stderr.txt"
#SBATCH --output="${dir_genomes}/${dir_SC}/err_out/${job_name}.%A.stdout.txt"

# Extract tarball command
download_extract_tarball \
    "${URL_SC}/${tarball_SC}" \
    "${dir_genomes}/${dir_SC}"
EOF
    else
        echo "Warning: Path ${dir_genomes}/${dir_SC}/${tarball_SC%.tgz} exist; skipping download and extraction."
    fi
fi

#  Move the downloaded S. cerevisiae files to appropriate storage locations
if [[ -d "${dir_genomes}/${dir_SC}/${tarball_SC%.tgz}" ]]; then
    mv \
        "${dir_genomes}/${dir_SC}/${tarball_SC%.tgz}/"*.sgd.gz \
        "${dir_genomes}/${dir_SC}/sgd"

    mv \
        "${dir_genomes}/${dir_SC}/${tarball_SC%.tgz}/"*.gff.gz \
        "${dir_genomes}/${dir_SC}/gff3"

    mv \
        "${dir_genomes}/${dir_SC}/${tarball_SC%.tgz}/"*.{fasta,fsa}.gz \
        "${dir_genomes}/${dir_SC}/fasta"

    if [[ "$(
        find "${dir_genomes}/${dir_SC}/${tarball_SC%.tgz}" \
            -mindepth 1 \
            -maxdepth 1 \
            -print \
            -quit
    )" ]]; then
        echo "Warning: Path ${dir_genomes}/${dir_SC}/${tarball_SC%.tgz} is not empty; skipping directory removal."
    else
        rmdir "${dir_genomes}/${dir_SC}/${tarball_SC%.tgz}"
    fi
fi


#  Download and store Schizosaccharomyces pombe FASTA and GFF3 files ----------
#  Set flags
check_operations=true  # Check operations to download genome files
run_operations=true    # Run operations to download genome files

#  Loop through FASTA and GFF3 arrays for Schizosaccharomyces pombe
iter=0
for file_type in "fasta" "gff3"; do
    # file_type="fasta"

    eval array=( \"\${${file_type}_SP[@]}\" )
    url_base="URL_SP_${file_type}"

    for file in "${array[@]}"; do
        # file="${array[0]}"
        
        (( iter++ ))
        url="${!url_base}/${file}"
        output_file="${dir_genomes}/${dir_SP}/${file_type}/${file}"
        job_name="download.${dir_SP}.${file_type}.${iter}"

        #  Echo loop-dependent variables if check_variables is true
        if ${check_variables}; then
            echo "
            file=${file}
            url=${url}
            output_file=${output_file}
            "
        fi

        #  Echo the download logic if check_operations is true
        if ${check_operations}; then
            echo "
            ### ${iter} ###

            if \${run_operations} && [[ ! -f \"\${output_file}\" ]]; then
sbatch << EOF
#!/bin/bash

#SBATCH --job-name=\"${job_name}\"
#SBATCH --nodes=1
#SBATCH --time=${time}
#SBATCH --error=\"${dir_genomes}/${dir_SP}/err_out/${job_name}.stderr.txt\"
#SBATCH --output=\"${dir_genomes}/${dir_SP}/err_out/${job_name}.stdout.txt\"

# Download command
download_file \\
    \"${url}\" \\
    \"${output_file}\"
EOF
            fi
            "
        fi

        #  Download the file if run_operations is true
        if ${run_operations} && [[ ! -f "${output_file}" ]]; then
sbatch << EOF
#!/bin/bash

#SBATCH --job-name="${job_name}"
#SBATCH --nodes=1
#SBATCH --time=${time}
#SBATCH --error="${dir_genomes}/${dir_SP}/err_out/${job_name}.stderr.txt"
#SBATCH --output="${dir_genomes}/${dir_SP}/err_out/${job_name}.stdout.txt"

# Download command
download_file \
    "${url}" \
    "${output_file}"
EOF
        fi

        sleep 0.2
    done
done


#  Prepare S. cerevisiae fasta for concatenation with S. pombe fasta ----------
check_file=true
check_directory=true

#  Copy relevant S. cerevisiae fasta from fasta/ to fasta-processed/
cp \
    "${dir_genomes}/${dir_SC}/fasta/S288C_reference_sequence_R64-3-1_20210421.fsa.gz" \
    "${dir_genomes}/${dir_SC}/fasta-processed"

#  Check the chromosome names in the fasta
if ${check_file}; then
    zcat "${dir_genomes}/${dir_SC}/fasta-processed/"*.fsa.gz | grep "^>"
fi

#  Simplify the names of chromosomes in the S. cerevisiae fasta
if [[ -f "${dir_genomes}/${dir_SC}/fasta-processed/tmp.fa" ]]; then
    rm "${dir_genomes}/${dir_SC}/fasta-processed/tmp.fa"
fi

zcat "${dir_genomes}/${dir_SC}/fasta-processed/S288C_reference_sequence_R64-3-1_20210421.fsa.gz" \
    | sed 's/^>ref|NC_001133|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=I\]/>I/g;s/^>ref|NC_001134|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=II\]/>II/g;s/^>ref|NC_001135|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=III\]/>III/g;s/^>ref|NC_001136|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=IV\]/>IV/g;s/^>ref|NC_001137|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=V\]/>V/g;s/^>ref|NC_001138|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=VI\]/>VI/g;s/^>ref|NC_001139|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=VII\]/>VII/g;s/^>ref|NC_001140|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=VIII\]/>VIII/g;s/^>ref|NC_001141|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=IX\]/>IX/g;s/^>ref|NC_001142|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=X\]/>X/g;s/^>ref|NC_001143|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=XI\]/>XI/g;s/^>ref|NC_001144|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=XII\]/>XII/g;s/^>ref|NC_001145|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=XIII\]/>XIII/g;s/^>ref|NC_001146|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=XIV\]/>XIV/g;s/^>ref|NC_001147|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=XV\]/>XV/g;s/^>ref|NC_001148|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=XVI\]/>XVI/g;s/^>ref|NC_001224|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[location=mitochondrion\]\ \[top=circular\]/>Mito/g' \
        > "${dir_genomes}/${dir_SC}/fasta-processed/tmp.fa"

#  Check the new chromosome names in "tmp.fa"
if ${check_file}; then
    cat "${dir_genomes}/${dir_SC}/fasta-processed/tmp.fa" | grep "^>"
fi

#  Remove the intial infile
rm "${dir_genomes}/${dir_SC}/fasta-processed/S288C_reference_sequence_R64-3-1_20210421.fsa.gz"

#  Rename "tmp.fa" to "S288C_reference_sequence_R64-3-1_20210421.fa"
mv -f \
    "${dir_genomes}/${dir_SC}/fasta-processed/tmp.fa" \
    "${dir_genomes}/${dir_SC}/fasta-processed/S288C_reference_sequence_R64-3-1_20210421.fa"

#  Compress "S288C_reference_sequence_R64-3-1_20210421.fa"
gzip "${dir_genomes}/${dir_SC}/fasta-processed/S288C_reference_sequence_R64-3-1_20210421.fa"

#  Check the directory contents
if ${check_directory}; then
    ls -lhaFG "${dir_genomes}/${dir_SC}/fasta-processed"
fi

#  Check the chromosome names
if ${check_file}; then
    zcat "${dir_genomes}/${dir_SC}/fasta-processed/S288C_reference_sequence_R64-3-1_20210421.fa.gz" \
        | grep "^>"
fi


#  Prepare S. cerevisiae gff3 for concatenation with S. pombe gff3 ------------
check_file=true
check_directory=true

#  Copy relevant S. cerevisiae gff3 from gff3/ to gff3-processed/
cp \
    "${dir_genomes}/${dir_SC}/gff3/saccharomyces_cerevisiae_R64-3-1_20210421.gff.gz" \
    "${dir_genomes}/${dir_SC}/gff3-processed"

#  What do chromosome names look like?
if ${check_file}; then
    zcat "${dir_genomes}/${dir_SC}/gff3-processed/saccharomyces_cerevisiae_R64-3-1_20210421.gff.gz" \
        | cut -f 1 \
        | sort \
        | uniq \
        | head -100
fi
#NOTE There are fasta sequences in this file that need to be excluded

#  Decompress "saccharomyces_cerevisiae_R64-3-1_20210421.gff.gz"
gzip -cd \
    "${dir_genomes}/${dir_SC}/gff3-processed/saccharomyces_cerevisiae_R64-3-1_20210421.gff.gz" \
        > "${dir_genomes}/${dir_SC}/gff3-processed/saccharomyces_cerevisiae_R64-3-1_20210421.gff"

#  Remove everything after line containing ### (fasta chromosome sequences)
if [[ -f "${dir_genomes}/${dir_SC}/gff3-processed/tmp.gff3" ]]; then
    rm "${dir_genomes}/${dir_SC}/gff3-processed/tmp.gff3"
fi

sed -n '/###/q;p' \
    < "${dir_genomes}/${dir_SC}/gff3-processed/saccharomyces_cerevisiae_R64-3-1_20210421.gff" \
    > "${dir_genomes}/${dir_SC}/gff3-processed/tmp.gff3"

#  Check the chromosome names now that fasta sequences are gone
if ${check_file}; then
    cat "${dir_genomes}/${dir_SC}/gff3-processed/tmp.gff3" \
        | cut -f 1 \
        | sort \
        | uniq
fi

#  Use sed to rename chromosomes; save the results to "tmp_2.gff3"
cat "${dir_genomes}/${dir_SC}/gff3-processed/tmp.gff3" \
    | sed 's/^chr//g;s/^mt/Mito/g' \
        > "${dir_genomes}/${dir_SC}/gff3-processed/tmp_2.gff3"

#  Check the chromosome names in the updated file
if ${check_file}; then
    cat "${dir_genomes}/${dir_SC}/gff3-processed/tmp_2.gff3" \
        | cut -f 1 \
        | sort \
        | uniq
fi

#  Check on the file contents
if ${check_file}; then
    echo "### head ###"
    head -25 "${dir_genomes}/${dir_SC}/gff3-processed/tmp_2.gff3"
    echo ""

    echo "### tail ###"
    tail "${dir_genomes}/${dir_SC}/gff3-processed/tmp_2.gff3"
    echo ""
fi

#  Rename "tmp_2.gff3" to "saccharomyces_cerevisiae_R64-3-1_20210421.gff3"
mv \
    "${dir_genomes}/${dir_SC}/gff3-processed/tmp_2.gff3" \
    "${dir_genomes}/${dir_SC}/gff3-processed/saccharomyces_cerevisiae_R64-3-1_20210421.gff3"

#  Compress "saccharomyces_cerevisiae_R64-3-1_20210421.gff3"
gzip "${dir_genomes}/${dir_SC}/gff3-processed/saccharomyces_cerevisiae_R64-3-1_20210421.gff3"

#  Check the directory contents
if ${check_directory}; then
    ls -lhaFG "${dir_genomes}/${dir_SC}/gff3-processed"
fi

#  Remove unneeded files
rm \
    "${dir_genomes}/${dir_SC}/gff3-processed/saccharomyces_cerevisiae_R64-3-1_20210421.gff" \
    "${dir_genomes}/${dir_SC}/gff3-processed/saccharomyces_cerevisiae_R64-3-1_20210421.gff.gz" \
    "${dir_genomes}/${dir_SC}/gff3-processed/tmp.gff3"

#  Check the directory contents again
if ${check_directory}; then
    ls -lhaFG "${dir_genomes}/${dir_SC}/gff3-processed"
fi

#  Check chromosome names again
if ${check_file}; then
    zcat "${dir_genomes}/${dir_SC}/gff3-processed/saccharomyces_cerevisiae_R64-3-1_20210421.gff3.gz" \
        | cut -f 1 \
        | uniq
fi


#  Prepare S. pombe fasta for concatenation with S. cerevesiae fasta ----------
check_file=true
check_directory=true

#  Copy relevant S. pombe fasta from fasta/ to fasta-processed/
cp \
    "${dir_genomes}/${dir_SP}/fasta/Schizosaccharomyces_pombe_all_chromosomes.fa.gz" \
    "${dir_genomes}/${dir_SP}/fasta-processed"

#  What do chromosome names look like?
if ${check_file}; then
    zgrep "^>" "${dir_genomes}/${dir_SP}/fasta-processed/Schizosaccharomyces_pombe_all_chromosomes.fa.gz"
fi

#  Create a decompressed version of the fasta
gzip -cd \
    "${dir_genomes}/${dir_SP}/fasta-processed/Schizosaccharomyces_pombe_all_chromosomes.fa.gz" \
        > "${dir_genomes}/${dir_SP}/fasta-processed/Schizosaccharomyces_pombe_all_chromosomes.fa"

#  Use sed to rename chromosomes; save the results to "tmp.fa"
if [[ -f "${dir_genomes}/${dir_SP}/fasta-processed/tmp.fa" ]]; then
    rm "${dir_genomes}/${dir_SP}/fasta-processed/tmp.fa"
fi

cat "${dir_genomes}/${dir_SP}/fasta-processed/Schizosaccharomyces_pombe_all_chromosomes.fa" \
    | sed 's/^>chr_II_telomeric_gap\ Schizosaccharomyces_pombe/>SP_II_TG/g;s/^>I\ Schizosaccharomyces_pombe/>SP_I/g;s/^>II\ Schizosaccharomyces_pombe/>SP_II/g;s/^>III\ Schizosaccharomyces_pombe/>SP_III/g;s/^>mating_type_region\ Schizosaccharomyces_pombe/>SP_MTR/g;s/^>mitochondrial\ Schizosaccharomyces_pombe/>SP_Mito/g' \
        > "${dir_genomes}/${dir_SP}/fasta-processed/tmp.fa"

#  Check the new chromosome names
if ${check_file}; then
    cat "${dir_genomes}/${dir_SP}/fasta-processed/tmp.fa" | grep "^>"
fi

#  Overwrite the initial file with the contents of "tmp.fa"
mv -f \
    "${dir_genomes}/${dir_SP}/fasta-processed/tmp.fa" \
    "${dir_genomes}/${dir_SP}/fasta-processed/Schizosaccharomyces_pombe_all_chromosomes.fa"

#  Double check chromosome names
if ${check_file}; then
    cat "${dir_genomes}/${dir_SP}/fasta-processed/Schizosaccharomyces_pombe_all_chromosomes.fa" \
        | grep "^>"
fi

#  Remove the compressed initial file
rm "${dir_genomes}/${dir_SP}/fasta-processed/"*.gz

#  Compress the updated file (with new chromosome names)
gzip "${dir_genomes}/${dir_SP}/fasta-processed/"*.fa

#  Check the directory
if ${check_directory}; then
    ls -lhaFG "${dir_genomes}/${dir_SP}/fasta-processed"
fi

#  Check chromosome names again
if ${check_file}; then
    zcat "${dir_genomes}/${dir_SP}/fasta-processed/Schizosaccharomyces_pombe_all_chromosomes.fa.gz" \
        | grep "^>"
fi


#  Prepare S. pombe gff3 for concatenation with S. cerevesiae gff3 ------------
check_file=true
check_directory=true

#  Copy relevant S. pombe fasta from fasta/ to fasta-processed/
cp \
    "${dir_genomes}/${dir_SP}/gff3/Schizosaccharomyces_pombe_all_chromosomes.gff3.gz" \
    "${dir_genomes}/${dir_SP}/gff3-processed"

#  What do chromosome names look like?
if ${check_file}; then
    zcat "${dir_genomes}/${dir_SP}/gff3-processed/Schizosaccharomyces_pombe_all_chromosomes.gff3.gz" \
        | cut -f 1 \
        | sort \
        | uniq
fi

#  Use sed to rename chromosomes; save the results to "tmp.gff3"
zcat "${dir_genomes}/${dir_SP}/gff3-processed/Schizosaccharomyces_pombe_all_chromosomes.gff3.gz" \
    | sed 's/^chr_II_telomeric_gap/SP_II_TG/g;s/^I/SP_I/g;s/^II/SP_II/g;s/^III/SP_III/g;s/^mating_type_region/SP_MTR/g;s/^mitochondrial/SP_Mito/g' \
        > "${dir_genomes}/${dir_SP}/gff3-processed/tmp.gff3"

#  Check on the file contents
if ${check_file}; then
    echo "### head ###"
    head "${dir_genomes}/${dir_SP}/gff3-processed/tmp.gff3"
    echo ""

    echo "### tail ###"
    tail "${dir_genomes}/${dir_SP}/gff3-processed/tmp.gff3"
    echo ""

    echo "### tail (initial) ###"
    zcat "${dir_genomes}/${dir_SP}/gff3/Schizosaccharomyces_pombe_all_chromosomes.gff3.gz" \
        | tail
fi

#  Check on the chromosome names in "tmp.gff3"
if ${check_file}; then
    cat "${dir_genomes}/${dir_SP}/gff3-processed/tmp.gff3" \
        | cut -f 1 \
        | sort \
        | uniq
fi

#  Rename "tmp.gff3" to "Schizosaccharomyces_pombe_all_chromosomes.gff3"
mv -f \
    "${dir_genomes}/${dir_SP}/gff3-processed/tmp.gff3" \
    "${dir_genomes}/${dir_SP}/gff3-processed/Schizosaccharomyces_pombe_all_chromosomes.gff3"

#  Remove the initial gzipped file
rm "${dir_genomes}/${dir_SP}/gff3-processed/Schizosaccharomyces_pombe_all_chromosomes.gff3.gz"

#  Compress "Schizosaccharomyces_pombe_all_chromosomes.gff3"
gzip "${dir_genomes}/${dir_SP}/gff3-processed/Schizosaccharomyces_pombe_all_chromosomes.gff3"

#  Check the directory contents
if ${check_directory}; then
    ls -lhaFG "${dir_genomes}/${dir_SP}/gff3-processed"
fi

#  Check chromosome names again
if ${check_file}; then
    zcat "${dir_genomes}/${dir_SP}/gff3-processed/Schizosaccharomyces_pombe_all_chromosomes.gff3.gz" \
        | cut -f 1 \
        | sort \
        | uniq
fi
```
</details>
<br />
