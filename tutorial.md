
`#tutorial.md`
<br />

<details>
<summary><b><font size="+2"><i>Table of contents</i></font></b></summary>
<br />
<!-- MarkdownTOC -->

1. [1. Prepare FASTQ files of interest for processing and analyses](#1-prepare-fastq-files-of-interest-for-processing-and-analyses)
    1. [Code](#code)
1. [2. Adapter- and quality-trim the FASTQ files](#2-adapter--and-quality-trim-the-fastq-files)
    1. [a. Install Atria and dependencies](#a-install-atria-and-dependencies)
        1. [Code](#code-1)
    1. [b. Adapter- and quality-trim the FASTQ files using Atria](#b-adapter--and-quality-trim-the-fastq-files-using-atria)
        1. [Code](#code-2)
1. [3. Align trimmed FASTQ files](#3-align-trimmed-fastq-files)
    1. [a. Generate a concatenated annotated assembly of the *S. cerevisiae* and *S. pombe* genomes](#a-generate-a-concatenated-annotated-assembly-of-the-s-cerevisiae-and-s-pombe-genomes)
        1. [Code](#code-3)
    1. [b. Concatenate processed FASTA and GFF3 files in new directory `combined_SC_SP/`](#b-concatenate-processed-fasta-and-gff3-files-in-new-directory-combined_sc_sp)
        1. [Code](#code-4)
    1. [c. Create `bowtie2` and `bwa` indices for "`combined_SC_SP.fa.gz`"](#c-create-bowtie2-and-bwa-indices-for-combined_sc_spfagz)
        1. [Code](#code-5)
    1. [d. Use `bowtie2` to align the trimmed FASTQ files](#d-use-bowtie2-to-align-the-trimmed-fastq-files)
        1. [Code](#code-6)
    1. [e. Use `bwa` to align the trimmed FASTQ files](#e-use-bwa-to-align-the-trimmed-fastq-files)
        1. [Code](#code-7)
1. [X. Run `phantompeakqualtools` on aligned data](#x-run-phantompeakqualtools-on-aligned-data)
    1. [a. Install `phantompeakqualtools`](#a-install-phantompeakqualtools)
        1. [Code](#code-8)
        1. [Printed \(remote\)](#printed-remote)
        1. [Printed \(local\)](#printed-local)
    1. [b. Run `phantompeakqualtools`](#b-run-phantompeakqualtools)
        1. [Code](#code-9)
1. [4. Call peaks with MACS3](#4-call-peaks-with-macs3)
    1. [a. Install MACS3](#a-install-macs3)
        1. [Code](#code-10)
    1. [b. Run MACS3](#b-run-macs3)
        1. [Code](#code-11)
    1. [c. Run MACS3 with pooled replicates](#c-run-macs3-with-pooled-replicates)
        1. [Code](#code-12)
1. [5. Subset peaks](#5-subset-peaks)
    1. [a. Install environment for interactive R scripting](#a-install-environment-for-interactive-r-scripting)
        1. [Bash code](#bash-code)
        1. [R code](#r-code)
    1. [b. Perform set operations with peak intervals](#b-perform-set-operations-with-peak-intervals)
        1. [i. Get situated](#i-get-situated)
            1. [Bash code](#bash-code-1)
        1. [ii. Perform set operations with peak intervals, returning complete intervals from a given set](#ii-perform-set-operations-with-peak-intervals-returning-complete-intervals-from-a-given-set)
            1. [R code](#r-code-1)
        1. [iii. Perform set operations with peak intervals, returning partial intervals from a given set](#iii-perform-set-operations-with-peak-intervals-returning-partial-intervals-from-a-given-set)
            1. [R code](#r-code-2)
        1. [iv. Perform additional set operations with respect to three-way intersections](#iv-perform-additional-set-operations-with-respect-to-three-way-intersections)
            1. [R code](#r-code-3)
1. [6. Calculate sample scaling factors from *S. pombe* spike-ins](#6-calculate-sample-scaling-factors-from-s-pombe-spike-ins)
1. [7. Miscellaneous \(to be organized\)](#7-miscellaneous-to-be-organized)
    1. [x. Scratch](#x-scratch)
        1. [Code](#code-13)
        1. [Notes](#notes)
    1. [a. Determine the locations of low-complexity regions in *S. cerevisiae*](#a-determine-the-locations-of-low-complexity-regions-in-s-cerevisiae)
        1. [i. Install `sdust` via `minimap`](#i-install-sdust-via-minimap)
            1. [Code](#code-14)
            1. [Printed](#printed)
        1. [ii. Run `sdust` via `minimap`](#ii-run-sdust-via-minimap)
            1. [Code](#code-15)
    1. [b. Determine the effective genome size of *S. cerevisiae* \(50-mers\)](#b-determine-the-effective-genome-size-of-s-cerevisiae-50-mers)
        1. [i. Install `khmer`](#i-install-khmer)
            1. [Code](#code-16)
            1. [Printed](#printed-1)
        1. [ii. Run `khmer`](#ii-run-khmer)
            1. [Code](#code-17)
            1. [Printed](#printed-2)
    1. [b. Determine base statistics in *S. cerevisiae* FASTA files](#b-determine-base-statistics-in-s-cerevisiae-fasta-files)
        1. [i. Install `faCount`](#i-install-facount)
            1. [Code](#code-18)
        1. [ii. Run `faCount`](#ii-run-facount)
            1. [Code](#code-19)
            1. [Printed](#printed-3)

<!-- /MarkdownTOC -->
</details>
<br />

<a id="1-prepare-fastq-files-of-interest-for-processing-and-analyses"></a>
# 1. Prepare FASTQ files of interest for processing and analyses
<a id="code"></a>
## Code
<details>
<summary><i>Code: 1. Prepare FASTQ files of interest for processing and analyses</i></summary>

```bash
#!/bin/bash

#  Define functions ===========================================================
#  Function to return an error message and exit code 1, which stops the
#+ interactive execution of code
function error_and_return() {
    echo "Error: ${1}" >&2
    return 1
}


#  Initialize variables and arrays ============================================
#  Initialize variables for directories of interest
dir_base="${HOME}/tsukiyamalab"                          # Base directory where lab data is stored
dir_repo="Kris/2023_rDNA"                                # Repository directory for the specific project
dir_work="results/2023-0406_tutorial_ChIP-seq_analyses"  # Working directory for storing results of this tutorial
# dir_orig="Rina/ChIP-seq/230915_hho1_hmo1_rhirano"      # Directory where original ChIP-seq data is stored
dir_orig="Rachel/misc_data/experiments/ChIPs/HDAC_HAT_Q/230915_ChIPseq"  # Directory where original ChIP-seq data is stored
dir_sym="01_sym"                                         # Directory name where symbolic links will be stored

#  Create an associative array/hash map for renaming files through symbolic
#+ links. Original file stems are keys, and new file stems are values. This
#+ is the naming scheme for symlink-renamed files: assay_state_factor_strain
#+ 
#+ - Here, "IP" signifies ChIP-seq "immunoprecipitate" assay or experiment,
#+   i.e., the ChIP-seq data for the factor of interest
#+ - "in" denotes the ChIP-seq "input" assay, the... #TODO Write this
unset file_fastqs && typeset -A file_fastqs=(
    # ["6336_G1_in_S15"]="in_G1_Hho1_6336"
    # ["6336_G1_lP_S27"]="IP_G1_Hho1_6336"
    # ["6336_G2M_in_S17"]="in_G2M_Hho1_6336"
    # ["6336_G2M_lP_S29"]="IP_G2M_Hho1_6336"
    # ["6336_Q-in_S19"]="in_Q_Hho1_6336"
    # ["6336_Q_lP_S31"]="IP_Q_Hho1_6336"
    # ["6337_G1_in_S16"]="in_G1_Hho1_6337"
    # ["6337_GI_IP_S28"]="IP_G1_Hho1_6337"
    # ["6337_G2M_in_S18"]="in_G2M_Hho1_6337"
    # ["6337_G2M_lP_S30"]="IP_G2M_Hho1_6337"
    # ["6337_Q_in_S20"]="in_Q_Hho1_6337"
    # ["6337_Q_lP_S32"]="IP_Q_Hho1_6337"
    # ["7750_G1_in_S21"]="in_G1_Hmo1_7750"
    # ["7750_G1_lP_S33"]="IP_G1_Hmo1_7750"
    # ["7750_G2M_in_S23"]="in_G2M_Hmo1_7750"
    # ["7750_G2M_lP_S35"]="IP_G2M_Hmo1_7750"
    # ["7750_Q_ln_S25"]="in_Q_Hmo1_7750"
    # ["7750_Q_lP_S37"]="IP_Q_Hmo1_7750"
    # ["7751_G1_in_S22"]="in_G1_Hmo1_7751"
    # ["7751_G1_lP_S34"]="IP_G1_Hmo1_7751"
    # ["7751_G2M_in_S24"]="in_G2M_Hmo1_7751"
    # ["7751_G2M_lP_S36"]="IP_G2M_Hmo1_7751"
    # ["7751_Q_ln_S26"]="in_Q_Hmo1_7751"
    # ["7751_Q_lP_S38"]="IP_Q_Hmo1_7751"
    ["5781_input_S13"]="in_Q_untagged_5781"
    ["5781_IP_S14"]="IP_Q_untagged_5781"
    ["7041_input_S11"]="in_Q_Esa5_7041"
    ["7041_IP_S12"]="IP_Q_Esa5_7041"
    ["7568_input_S7"]="in_Q_Rpd3_7568"
    ["7568_IP_S8"]="IP_Q_Rpd3_7568"
    ["7569_input_S3"]="in_Q_Rpd3_7569"
    ["7569_IP_S4"]="IP_Q_Rpd3_7569"
    ["7691_input_S1"]="in_Q_Esa5_7691"
    ["7691_IP_S2"]="IP_Q_Esa5_7691"
    ["7692_input_S9"]="in_Q_Gcn5_7692"
    ["7692_IP_S10"]="IP_Q_Gcn5_7692"
    ["7709_input_S5"]="in_Q_Gcn5_7709"
    ["7709_IP_S6"]="IP_Q_Gcn5_7709"
)
# _R1_001.fastq.gz
# _R2_001.fastq.gz


#  Do the main work ===========================================================
#  Set flags to check variable and array assignments
check_variables=true
check_array=true

#  If check_variables is true, then echo the the variable assignments
if ${check_variables}; then
    echo "
    dir_base=${dir_base}
    dir_repo=${dir_repo}
    dir_work=${dir_work}
    dir_orig=${dir_orig}
    dir_sym=${dir_sym}
    "
fi

#  If check_array is true, then echo the the hash keys and values
if ${check_array}; then
    for i in "${!file_fastqs[@]}"; do
        key="${i}"
        value="${file_fastqs["${key}"]}"

        echo "
        key .......... ${key}
        value ........ ${value}
        "
    done
fi

#  Create the work directory if it does not exist
if [[ ! -d "${dir_base}/${dir_repo}/${dir_work}" ]]; then
    mkdir -p "${dir_base}/${dir_repo}/${dir_work}"
fi

#  Navigate to the work directory, and echo a message if navigation fails
cd "${dir_base}/${dir_repo}/${dir_work}" \
    || error_and_return "Failed to cd to ${dir_base}/${dir_repo}/${dir_work}."

#  If it doesn't exist, then create a directory to store symlinked FASTQ files
if [[ ! -d "${dir_sym}" ]]; then
    mkdir -p "${dir_sym}"
fi

#  Set flags for checking and running symlinking operations
check_operations=true
run_operations=true

#  Loop through each entry in the associative array to create symbolic links
for i in "${!file_fastqs[@]}"; do
    key="${i}"
    value="${file_fastqs["${key}"]}"
    
    #  If check_operations is true, then echo the operations that will be run
    if ${check_operations}; then
        echo "
        ### ${key} ###

        #  Check if the original file for read #1 exists before creating a symlink
        if [[ -f \"${dir_base}/${dir_orig}/${key}_R1_001.fastq.gz\" ]]; then
            ln -s \\
                \"${dir_base}/${dir_orig}/${key}_R1_001.fastq.gz\" \\
                \"${dir_sym}/${value}_R1.fastq.gz\"
        fi

        #  Check if the original file for read #2 exists before creating a symlink
        if [[ -f \"${dir_base}/${dir_orig}/${key}_R2_001.fastq.gz\" ]]; then
            ln -s \\
                \"${dir_base}/${dir_orig}/${key}_R2_001.fastq.gz\" \\
                \"${dir_sym}/${value}_R2.fastq.gz\"
        fi
        "
    fi

    #  If run_operations is true, then create the symbolic links
    if ${run_operations}; then
        #  Check if the original file for read #1 exists before creating a symlink
        if [[ -f "${dir_base}/${dir_orig}/${key}_R1_001.fastq.gz" ]]; then
            ln -s \
                "${dir_base}/${dir_orig}/${key}_R1_001.fastq.gz" \
                "${dir_sym}/${value}_R1.fastq.gz"
        fi

        #  Check if the original file for read #2 exists before creating a symlink
        if [[ -f "${dir_base}/${dir_orig}/${key}_R2_001.fastq.gz" ]]; then
            ln -s \
                "${dir_base}/${dir_orig}/${key}_R2_001.fastq.gz" \
                "${dir_sym}/${value}_R2.fastq.gz"
        fi
    fi
done

#  If check_operations is true, list the contents of the symlink storage
#+ directory
if ${check_operations}; then
    ls -lhaFG "${dir_sym}"
fi
```
</details>
<br />
<br />

<a id="2-adapter--and-quality-trim-the-fastq-files"></a>
# 2. Adapter- and quality-trim the FASTQ files
<a id="a-install-atria-and-dependencies"></a>
## a. Install Atria and dependencies
<a id="code-1"></a>
### Code
<details>
<summary><i>Code: 2.a. Install Atria and dependencies</i></summary>

```bash
#!/bin/bash

#  Define functions ===========================================================
#  Function to return an error message and exit code 1, which stops the
#+ interactive execution of code
function error_and_return() {
    echo "Error: ${1}" >&2
    return 1
}


#  Function to append PATH update to configuration file
function update_shell_config() {
    local config_file="${1}"
    local stem="${2}"
    local line_to_add="export PATH=\$PATH:\$HOME/${stem}/bin"

    #  Check if line already exists to avoid duplicates
    if ! grep -q "${line_to_add}" "${config_file}"; then
        echo "Appending PATH update to ${config_file}"
        echo "${line_to_add}" >> "${config_file}"
    else
        echo "PATH update already present in ${config_file}"
        return 1
    fi

    return 0
}


#  Function to check if Mamba is installed
function check_mamba_installed() {
    if ! : mamba &> /dev/null; then
        echo "Mamba is not installed on your system. Mamba is a package manager" 
        echo "that makes package installations faster and more reliable."
        echo ""
        echo "For installation instructions, please check the following link:"
        echo "https://github.com/mamba-org/mamba#installation"
        return 1
    fi
    
    return 0
}


#  Function to check if a specific Conda/Mamba environment is installed
function check_env_installed() {
    local env_name="${1}"

    if conda env list | grep -q "^${env_name} "; then
        return 0
    else
        echo "Environment \"${env_name}\" is not installed."
        return 1
    fi
}


#  Initialize variables and arrays ============================================
#  Initialize variables
URL_begin="https://julialang-s3.julialang.org/bin"
URL_mid=""
tarball=""

#  Detect operating system (OS) and system architecture
os="$(uname -s)"
arch="$(uname -m)"

#  Based on OS and architecture, set the URL and tarball name
case "${os}" in
    "Linux")
        URL_mid="linux/x64/1.8"
        tarball="julia-1.8.5-linux-x86_64.tar.gz"
        ;;
    "Darwin")
        if [[ "${arch}" = "x86_64" ]]; then
            URL_mid="mac/x64/1.9"
            tarball="julia-1.9.4-mac64.tar.gz"
        elif [[ "${arch}" = "arm64" ]]; then
            URL_mid="mac/aarch64/1.9"
            tarball="julia-1.9.4-macaarch64.tar.gz"
        else
            error_and_return "Unsupported architecture: ${arch}"
        fi
        ;;
    *)
        error_and_return "Unsupported operating system: ${os}"
        ;;
esac

#  Initialize variable for untarred directory in HOME
untarred="$(echo ${tarball} | awk -F '-' '{ print $1"-"$2 }')"

#  Initialize variable for name of conda/mamba environment for Atria
#+ dependencies
env_name="Atria_env"

#  Define the directory for Atria installation
atria_dir="${HOME}/Atria"


#  Do the main work ===========================================================
#  Install Julia --------------------------------------------------------------
#  Set flags for running echo tests, operations, etc.
check_variables=true  # Echo the variables assigned above
check_binary=true     # Check if the Julia binary is installed/in PATH
check_operation=true  # Check the operation to download and install Julia
run_operation=false   # Run the operation to download and install Julia
update_path=false     # Update PATH to include Julia binary

#  Check and echo variables
if ${check_variables}; then
    echo "
    URL_begin=${URL_begin}
    URL_mid=${URL_mid}
    tarball=${tarball}
    "
fi

#  Check if Julia binary is in PATH
if ${check_binary}; then
    if type julia &> /dev/null; then
        echo "Julia is in the PATH."
        echo "Available Julia binaries:"
        type -a julia
    else
        error_and_return "Julia is not in the PATH."
    fi
fi

#  Check the operation to download Julia
if ${check_operation}; then
    echo "
    curl \\
        -L \"${URL_begin}/${URL_mid}/${tarball}\" \\
        -o \"\${HOME}/${tarball}\"

    if [[ ! -d "\${HOME}/${untarred}" ]]; then
        tar zxf "\${HOME}/${tarball}"
        export PATH=\${PATH}:\${HOME}/${untarred}/bin
    fi

    if \${update_path}; then
        #  Determine which shell configuration file to update
        case \"${os}\" in
            \"Linux\")
                if [[ -f \"${HOME}/.bashrc\" ]]; then
                    shell_config=\"${HOME}/.bashrc\"
                elif [[ -f \"${HOME}/.bash_profile\" ]]; then
                    shell_config=\"${HOME}/.bash_profile\"
                fi
                ;;
            \"Darwin\")
                shell_config=\"${HOME}/.zshrc\"
                ;;
            *)
                error_and_return \"No known shell configuration file found.\"
                ;;
        esac

        # Call the function to update the configuration file
        update_shell_config \"\${shell_config}\" \"${untarred}\"

        echo \"To apply the update to PATH, please restart the terminal or\"
        echo \"source the configuration file.\"
    else
        echo \"For ${untarred} to remain in PATH, please export ${untarred} to\"
        echo \"PATH in the configuration file.\"
    fi
    "
fi

#  Run the operation to download Julia
if ${run_operation}; then
    curl \
        -L "${URL_begin}/${URL_mid}/${tarball}" \
        -o "${HOME}/${tarball}"

    if [[ ! -d "${HOME}/${untarred}" ]]; then
        tar zxf "${HOME}/${tarball}"
        export PATH=${PATH}:${HOME}/${untarred}/bin
    fi

    if ${update_path}; then
        #  Determine which shell configuration file to update
        case "${os}" in
            "Linux")
                if [[ -f "${HOME}/.bashrc" ]]; then
                    shell_config="${HOME}/.bashrc"
                elif [[ -f "${HOME}/.bash_profile" ]]; then
                    shell_config="${HOME}/.bash_profile"
                fi
                ;;
            "Darwin")
                shell_config="${HOME}/.zshrc"
                ;;
            *)
                error_and_return "No known shell configuration file found."
                ;;
        esac

        #  Call the function to update the configuration file
        update_shell_config "${shell_config}" "${untarred}"

        echo "To apply the update to PATH, please restart the terminal or"
        echo "source the configuration file."
    else
        echo "For ${untarred} to remain in PATH, please export ${untarred} to"
        echo "PATH in the configuration file."
    fi
fi


#  Install Atria dependencies -------------------------------------------------
#  Set flag(s)
create_mamba_env=true  # Install mamba environment if not detected
update_path=true       # Update PATH to include Atria binary  #TODO Write this

#  Check that Mamba is installed and in PATH
check_mamba_installed

#  Check that environment assigned to env_name is installed; if environment
#+ assigned to env_name is not installed, run the following
if [[ $(check_env_installed "${env_name}") -eq 0 ]]; then
    #  Handle the case when the environment is already installed
    echo "Activating environment ${env_name}"
    
    if ! mamba activate "${env_name}" &> /dev/null; then
        #  If `mamba activate` fails, try using `source activate`
        if ! conda activate "${env_name}" &> /dev/null; then
            if ! source activate "${env_name}" &> /dev/null; then
                #  If `source activate` also fails, return an error
                error_and_return "Failed to activate environment \"${env_name}\"."
            fi
        fi
    else
        echo "Environment \"${env_name}\" activated using mamba."
    fi
else
    #  Handle the case when the environment is not installed
    echo "Creating environment ${env_name}"
    
    if ${create_mamba_env}; then
        #  Switch `--yes` is set, which means no user input is required
        mamba create \
            --yes \
            --name "${env_name}" \
            --channel conda-forge \
                parallel \
                pbzip2 \
                pigz \
                r-argparse \
                r-ggsci \
                r-plotly \
                r-tidyverse
    fi
fi


#  Install Atria --------------------------------------------------------------
#  Set flags
install_atria=true  # Install Atria or not

if ${install_atria}; then
    #  Check if git and Julia are available
    if ! type git &> /dev/null; then
        error_and_return "git is not installed or not in the PATH."
    fi

    if ! type julia &> /dev/null; then
        error_and_return "Julia is not installed or not in the PATH."
    fi

    #  Clone the Atria repository if it doesn't already exist
    if [[ ! -d "${atria_dir}" ]]; then
        cd "$(dirname "${atria_dir}")" \
            || error_and_return "Failed to cd to $(dirname "${atria_dir}")."
        
        git clone "https://github.com/cihga39871/Atria.git" \
            || error_and_return "Failed to clone Atria repository."
    else
        echo "Atria directory already exists. Skipping git clone."
    fi

    #  Change to the Atria directory
    cd "${atria_dir}" \
        || error_and_return "Failed to change to Atria directory."

    #  Environment containing Atria dependencies must be activated prior to
    #+ installation of Atria
    if [[ "${CONDA_DEFAULT_ENV}" != "${env_name}" ]]; then
        if [[ "${CONDA_DEFAULT_ENV}" != "base" ]]; then
            mamba deactivate
        fi

        if ! mamba activate "${env_name}" &> /dev/null; then
            #  If `mamba activate` fails, try using `source activate`
            if ! conda activate "${env_name}" &> /dev/null; then
                #  If `conda activate` fails, try using `source activate`
                if ! source activate "${env_name}" &> /dev/null; then
                    #  If `source activate` also fails, return an error
                    error_and_return "Failed to activate environment \"${env_name}\"."
                fi
            fi
        fi
    fi

    #FIXME Installation issue on macOS: github.com/cihga39871/Atria/issues/14
    #  Run the Julia script to build Atria
    if ! julia build_atria.jl; then
        error_and_return "Failed to build Atria."
    fi

    echo "Atria installed successfully."

    #TODO
    #  Add the trimming program to PATH if not already present
    if ! grep -q "${trim_prog_dir}/bin" <<< "${PATH}"; then
        export PATH="${PATH}:${trim_prog_dir}/bin"
        local shell_config="${HOME}/.bashrc"  # Adjust based on the user's shell
        echo "export PATH=\"${PATH}:${trim_prog_dir}/bin\"" >> "${shell_config}"
        echo "Path updated in ${shell_config}. Please restart the terminal or source the configuration file."
    fi
fi
```
</details>
<br />

<a id="b-adapter--and-quality-trim-the-fastq-files-using-atria"></a>
## b. Adapter- and quality-trim the FASTQ files using Atria
<a id="code-2"></a>
### Code
<details>
<summary><i>Code: 2.b. Adapter- and quality-trim the FASTQ files using Atria</i></summary>

```bash
#!/bin/bash

#  Define function ============================================================
#  Function to return an error message and exit code 1, which stops the
#+ interactive execution of code
function error_and_return() {
    echo "Error: ${1}" >&2
    return 1
}


#  Initialize variables and arrays ============================================
dir_base="${HOME}/tsukiyamalab"                                    # Base directory for lab data
dir_repo="Kris/2023_rDNA"                                          # Repository directory
dir_work="results/2023-0406_tutorial_ChIP-seq_analyses"            # Work directory
dir_sym="01_sym"                                                   # Directory with symlinked FASTQs
dir_trim="02_trim"                                                 # Directory for trimmed FASTQs
env_atria="Atria_env"                                              # Conda environment for Atria
path_atria="${dir_base}/${dir_repo}/src/Atria/app-3.2.2/bin/atria" # Atria executable path
time="1:00:00"                                                     # Job time for SLURM
threads=4                                                          # Number of threads for SLURM jobs

#  Initialize an indexed array with FASTQ file stems
unset file_fastqs && typeset -a file_fastqs=(
    "in_G1_Hho1_6336"
    "IP_G1_Hho1_6336"
    "in_G2M_Hho1_6336"
    "IP_G2M_Hho1_6336"
    "in_Q_Hho1_6336"
    "IP_Q_Hho1_6336"
    "in_G1_Hho1_6337"
    "IP_G1_Hho1_6337"
    "in_G2M_Hho1_6337"
    "IP_G2M_Hho1_6337"
    "in_Q_Hho1_6337"
    "IP_Q_Hho1_6337"
    "in_G1_Hmo1_7750"
    "IP_G1_Hmo1_7750"
    "in_G2M_Hmo1_7750"
    "IP_G2M_Hmo1_7750"
    "in_Q_Hmo1_7750"
    "IP_Q_Hmo1_7750"
    "in_G1_Hmo1_7751"
    "IP_G1_Hmo1_7751"
    "in_G2M_Hmo1_7751"
    "IP_G2M_Hmo1_7751"
    "in_Q_Hmo1_7751"
    "IP_Q_Hmo1_7751"
    # "in_Q_untagged_5781"
    # "IP_Q_untagged_5781"
    # "in_Q_Esa5_7041"
    # "IP_Q_Esa5_7041"
    # "in_Q_Rpd3_7568"
    # "IP_Q_Rpd3_7568"
    # "in_Q_Rpd3_7569"
    # "IP_Q_Rpd3_7569"
    # "in_Q_Esa5_7691"
    # "IP_Q_Esa5_7691"
    # "in_Q_Gcn5_7692"
    # "IP_Q_Gcn5_7692"
    # "in_Q_Gcn5_7709"
    # "IP_Q_Gcn5_7709"
)


#  Do the main work ===========================================================
#  Set flags for checking variable and array assignments
check_variables=true
check_array=true

#  If check_variables is true, then echo the variable assignments
if ${check_variables}; then
    echo "
    dir_base=${dir_base}
    dir_repo=${dir_repo}
    dir_work=${dir_work}
    dir_sym=${dir_sym}
    dir_trim=${dir_trim}
    env_atria=${env_atria}
    path_atria=${path_atria}
    threads=${threads}
    "
fi

#  Echo array contents if check_array is true
if ${check_array}; then
    for i in "${file_fastqs[@]}"; do
        file="${i}"

        echo "
        read #1 ...... ${file}_R1.fastq.gz
        read #2 ...... ${file}_R2.fastq.gz
        "
    done
fi

#  If not already activated, the activate conda environment
if [[ "${CONDA_DEFAULT_ENV}" != "${env_atria}" ]]; then
    if [[ ${CONDA_DEFAULT_ENV} != "base" ]]; then
        conda deactivate
    fi

    source activate "${env_atria}"
fi

#  Navigate to the work directory
cd "${dir_base}/${dir_repo}/${dir_work}" \
    || error_and_return "Failed to cd to ${dir_base}/${dir_repo}/${dir_work}."

#  If it doesn't exist, then create a directory to store trimmed FASTQ files
if [[ ! -d "${dir_trim}" ]]; then
    mkdir -p "${dir_trim}/err_out"
fi

#  Set flags: checking variables, checking and submitting trimming jobs
check_variables=false
check_operations=true
run_operations=true

for i in "${!file_fastqs[@]}"; do
    index="${i}"
    iter=$(( index + 1 ))
    stem=${file_fastqs["${index}"]}
    job_name="${dir_trim}.${stem}"
    read_1="${dir_sym}/${stem}_R1.fastq.gz"
    read_2="${dir_sym}/${stem}_R2.fastq.gz"
    trim_1="${dir_trim}/${stem}_R1.atria.fastq.gz"
    trim_2="${dir_trim}/${stem}_R2.atria.fastq.gz"

    #  Echo loop-dependent variables if check_variables is true
    if ${check_variables}; then
        echo "
        index=${index}
        iter=${iter}
        stem=${stem}
        job_name=${job_name}
        read_1=${read_1}
        read_2=${read_2}
        trim_1=${trim_1}
        trim_2=${trim_2}
        "
    fi

    #  Echo the Atria trimming command if check_operations is true
    if ${check_operations}; then
        echo "
        #  -------------------------------------
        ### ${iter} ###

        if [[
                 -f \"${read_1}\" \\
            &&   -f \"${read_2}\" \\
            && ! -f \"${trim_1}\" \\
            && ! -f \"${trim_2}\"
        ]]; then
sbatch << EOF
#!/bin/bash

#SBATCH --job-name=\"${job_name}\"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=${threads}
#SBATCH --time=${time}
#SBATCH --error=\"${dir_trim}/err_out/${job_name}.%A.stderr.txt\"
#SBATCH --output=\"${dir_trim}/err_out/${job_name}.%A.stdout.txt\"

\"${path_atria}\" \\
    -t ${threads} \\
    -r \"${read_1}\" \\
    -R \"${read_2}\" \\
    -o \"${dir_trim}\" \\
    --length-range 35:500
        else
            echo \"
            Warning: Trimmed fastqs for $(basename ${read_1%_R1.fastq.gz}) exist; skipping trimming.
            \"
        fi
        "
    fi

    #  Submit the Atria trimming job if run_operations is true
    if ${run_operations}; then
        if [[
                 -f "${read_1}" \
            &&   -f "${read_2}" \
            && ! -f "${trim_1}" \
            && ! -f "${trim_2}"
        ]]; then
sbatch << EOF
#!/bin/bash

#SBATCH --job-name="${job_name}"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=${threads}
#SBATCH --time=${time}
#SBATCH --error="${dir_trim}/err_out/${job_name}.%A.stderr.txt"
#SBATCH --output="${dir_trim}/err_out/${job_name}.%A.stdout.txt"

"${path_atria}" \
    -t "${threads}" \
    -r "${read_1}" \
    -R "${read_2}" \
    -o "${dir_trim}" \
    --length-range 35:500
EOF
        else
            echo "
            Warning: Trimmed fastqs for $(basename ${read_1%_R1.fastq.gz}) exist; skipping trimming.
            "
        fi
    fi

    sleep 0.2  # Short pause to prevent rapid job-submission overload
done
```
</details>
<br />
<br />

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

#  Initialize variables =======================================================
dir_base="${HOME}/tsukiyamalab"                          # Base directory for lab data
dir_repo="Kris/2023_rDNA"                                # Repository directory
dir_work="results/2023-0406_tutorial_ChIP-seq_analyses"  # Work directory
threads=1
time="1:00:00"
job_name="download-preprocess_SC-SP"


#  Do the main work ===========================================================
cd "${dir_base}/${dir_repo}/${dir_work}" ||
    echo "cd'ing failed; check on this"

#  Submit to SLURM a script that downloads and preprocesses the FASTA and GFF3
#+ files for S. cerevisiae and S. pombe
sbatch \
    --job-name="${job_name}" \
    --nodes=1 \
    --cpus-per-task=${threads} \
    --time=${time} \
    --error="${job_name}.%A.stderr.txt" \
    --output="${job_name}.%A.stdout.txt" \
    download-preprocess_FASTA-GFF3_SC-SP.sh
```
</details>
<br />

<a id="b-concatenate-processed-fasta-and-gff3-files-in-new-directory-combined_sc_sp"></a>
## b. Concatenate processed FASTA and GFF3 files in new directory `combined_SC_SP/`
<a id="code-4"></a>
### Code
<details>
<summary><i>Code: 3.b. Concatenate processed FASTA and GFF3 files in new directory `combined_SC_SP/`</i></summary>

```bash
#!/bin/bash

#  Initialize variables =======================================================
dir_work="${HOME}/tsukiyamalab/kalavatt/genomes"
dir_SC="${dir_work}/Saccharomyces_cerevisiae"
dir_SP="${dir_work}/Schizosaccharomyces_pombe"
name_comb="combined_SC_SP"

dir_comb="${dir_work}/${name_comb}"
#TODO Four individual arguments for SC and SP FASTA and GFF3 files


#  Do the main work ===========================================================
#  Set flags for checking variable assignments and files
check_variables=true
check_file=true

#  If check_variables is true, then echo the variable assignments
if ${check_variables}; then
    echo "
    dir_work=${dir_work}
    dir_SC=${dir_SC}
    dir_SP=${dir_SP}
    name_comb=${name_comb}
    
    dir_comb=${dir_comb}
    "
fi

#  Get situated  #TODO Better description for this subsection
#TODO Four individual arguments for SC and SP FASTA and GFF3 files
if [[ ! -d "${dir_comb}" ]]; then
    mkdir -p "${dir_comb}/"{fasta,gff3}
fi

cp \
    "${dir_SC}/fasta-processed/S288C_reference_sequence_R64-3-1_20210421.fa.gz" \
    "${dir_comb}/fasta"
cp \
    "${dir_SP}/fasta-processed/Schizosaccharomyces_pombe_all_chromosomes.fa.gz" \
    "${dir_comb}/fasta"

cp \
    "${dir_SC}/gff3-processed/saccharomyces_cerevisiae_R64-3-1_20210421.gff3.gz" \
    "${dir_comb}/gff3"
cp \
    "${dir_SP}/gff3-processed/Schizosaccharomyces_pombe_all_chromosomes.gff3.gz" \
    "${dir_comb}/gff3"

#  Handle FASTA files
if [[ -f "${dir_comb}/fasta/${name_comb}.fa.gz" ]]; then
    rm "${dir_comb}/fasta/${name_comb}.fa.gz"
fi

cat \
    "${dir_comb}/fasta/S288C_reference_sequence_R64-3-1_20210421.fa.gz" \
    "${dir_comb}/fasta/Schizosaccharomyces_pombe_all_chromosomes.fa.gz" \
        > "${dir_comb}/fasta/${name_comb}.fa.gz"

if ${check_file}; then
    ls -lhaFG "${dir_comb}/fasta"
fi

if ${check_file}; then
    zcat "${dir_comb}/fasta/${name_comb}.fa.gz" | grep "^>"
fi

#  Handle GFF3 files
if [[ -f "${dir_comb}/gff3/${name_comb}.gff3.gz" ]]; then
    rm "${dir_comb}/gff3/${name_comb}.gff3.gz"
fi

cat \
    "${dir_comb}/gff3/saccharomyces_cerevisiae_R64-3-1_20210421.gff3.gz" \
    "${dir_comb}/gff3/Schizosaccharomyces_pombe_all_chromosomes.gff3.gz" \
        > "${dir_comb}/gff3/${name_comb}.gff3.gz"

#  Process the concatenated GFF3 file
#+ 
#+ 1. Initial pipe to awk: ARS_consensus_sequence records have strands labeled
#+    "0"; these need to be adjusted or else downstream programs will throw
#+    errors when encountering the "0" strands; thus, we change "0" to "."
#+ 
#+ 2. Subsequent pipe to sed: We need to replace or remove special characters
#+    that IGV can't handle; also, these characters can break the formation of
#+    data frames from gff3 via rtacklayer::import() or readr::read_tsv(), etc.
#+ 
#+ 3. Next pipe to awk: We want to create a concatenated S. cerevisiae/S. pombe
#+    gff3 file with "intelligible feature names": "Name=" value is "ID="
#+    value; this needs to be something that makes sense to a human; thus, we
#+    run an `awk` command that does the following:
#+    A. If field `$3` is `gene`, `blocked_reading_frame`, `ncRNA_gene`,
#+       `pseudogene`, `rRNA_gene`, `snRNA_gene`, `snoRNA_gene`, `tRNA_gene`,
#+       or `telomerase_RNA_gene`, then check for the presence of the `gene=`
#+       attribute; if present, then replace the `ID=` value with the `gene=`
#+       value
#+    B. If field `$3` is `ARS`, then check for the presence of the `Alias=`
#+       attribute; if present, then replace the `ID=` value with the `Alias=`
#+       value
zcat "${dir_comb}/gff3/${name_comb}.gff3.gz" \
    | awk -F "\t" '
        BEGIN { OFS = FS } {
            if ($7=="0") { $7="."; print $0 } else { print $0 }
        }
    ' \
    | sed -e 's/%20/-/g' \
          -e 's/%2C//g' \
          -e 's/%3B//g' \
          -e 's/%28//g' \
          -e 's/%29//g' \
          -e 's/%//g' \
    | awk '
        BEGIN { FS=OFS="\t" }
        $1 ~ /^#/ { print; next }  # Skip lines starting with #
        $3 ~ /^(blocked_reading_frame|gene|ncRNA_gene|pseudogene|rRNA_gene|snRNA_gene|snoRNA_gene|tRNA_gene|telomerase_RNA_gene)$/ {
            if (match($9, /gene=[^;]+/)) {
                #  Extract the gene name
                gene_name = substr($9, RSTART+5, RLENGTH-5)
                #  Replace Name value with gene value
                sub(/Name=[^;]+/, "Name=" gene_name, $9)
            }
        }
        $3 == "ARS" {
            if (match($9, /Alias=[^;]+/)) {
                #  Extract the alias name
                alias_name = substr($9, RSTART+6, RLENGTH-6)
                #  Replace Name value with alias value
                sub(/Name=[^;]+/, "Name=" alias_name, $9)
            }
        }
        { print }
    ' \
    | gzip \
        > "${dir_comb}/gff3/tmp.gff3.gz"

if [[ -f "${dir_comb}/gff3/tmp.gff3.gz" ]]; then
    mv -f \
        "${dir_comb}/gff3/tmp.gff3.gz" \
        "${dir_comb}/gff3/${name_comb}.gff3.gz"
fi

if ${check_file}; then
    ls -lhaFG "${dir_comb}/gff3"
fi

if ${check_file}; then
    zcat "${dir_comb}/gff3/${name_comb}.gff3.gz" \
        | cut -f 1 \
        | grep -v "#" \
        | sort \
        | uniq
fi
```
</details>
<br />
<br />

<a id="c-create-bowtie2-and-bwa-indices-for-combined_sc_spfagz"></a>
## c. Create `bowtie2` and `bwa` indices for "`combined_SC_SP.fa.gz`"
<a id="code-5"></a>
### Code
<details>
<summary><i>Code: 3.c. Create `bowtie2` and `bwa` indices for "`combined_SC_SP.fa.gz`"</i></summary>

```bash
#!/bin/bash

#  Define functions ===========================================================
#  Function to return an error message and exit code 1, which stops the
#+ interactive execution of code
function error_and_return() {
    echo "Error: ${1}" >&2
    return 1
}

#  Function to check if a SLURM module is loaded or not; if not, the module is
#+ loaded
function check_and_load_module() {
    local module_name="${1}"
    if ! module is-loaded "${module_name}" &> /dev/null; then
        echo "Loading module: ${module_name}"
        module load "${module_name}"
    else
        echo "Module already loaded: ${module_name}"
    fi
}


#  Initialize variables =======================================================
dir_work="${HOME}/tsukiyamalab/kalavatt/genomes/combined_SC_SP"
file_fa="${dir_work}/fasta/combined_SC_SP.fa.gz"
mod_samtools="SAMtools/1.16.1-GCC-11.2.0"
mod_bowtie2="Bowtie2/2.4.4-GCC-11.2.0"
mod_bwa="BWA/0.7.17-GCCcore-11.2.0"
time="8:00:00"

job_bowtie2="build-indices_bowtie2"
job_bwa="build-indices_bwa"


#  Do the main work ===========================================================
#  Set flags for checking variable assignments and files
check_variables=true
check_file=true

#  If check_variables is true, then echo the variable assignments
if ${check_variables}; then
    echo "
    dir_work=${dir_work}
    file_fa=${file_fa}
    mod_samtools=${mod_samtools}
    mod_bowtie2=${mod_bowtie2}
    mod_bwa=${mod_bwa}
    time=${time}

    job_bowtie2=${job_bowtie2}
    job_bwa=${job_bwa}
    "
fi


#  Index the FASTA file -----------------------------------
check_and_load_module "${mod_samtools}"

#  If necessary, unzip the FASTA file
if [[ "${file_fa}" == *.gz ]]; then
    gzip -cd "${file_fa}" > "${file_fa%.gz}"
    file_fa="${file_fa%.gz}"
fi

#  Check on the chromosomes in the FASTA file
if ${check_file}; then cat "${file_fa}" | grep "^>"; fi

#  If necessary, index the FASTA file
if [[ ! -f "${file_fa}.fai" ]]; then samtools faidx "${file_fa}"; fi

#  Create a CHROM-INFO file using the FASTA index
cut -f 1,2 "${file_fa}.fai" > "${file_fa/.fa/.chrom-info.tsv}"


#  Build Bowtie 2 indices ---------------------------------
check_and_load_module "${mod_bowtie2}"

#  Create a directory for storing the Bowtie 2 indices
if [[ ! -d "${dir_work}/bowtie2" ]]; then
    mkdir -p "${dir_work}/bowtie2/err_out"
fi

#  Using a HEREDOC, submit to SLURM a job for creating Bowtie 2 indices
sbatch << EOF
#!/bin/bash

#SBATCH --job-name="${job_bowtie2}"
#SBATCH --nodes=1
#SBATCH --time=${time}
#SBATCH --error="${dir_work}/bowtie2/err_out/${job_bowtie2}.%A.stderr.txt"
#SBATCH --output="${dir_work}/bowtie2/err_out/${job_bowtie2}.%A.stdout.txt"

bowtie2-build \
    "${file_fa}" \
    "${dir_work}/bowtie2/$(basename ${file_fa} .fa)"
EOF


#  Build BWA indices --------------------------------------
check_and_load_module "${mod_bwa}"

#  Create a directory for storing the BWA indices
if [[ ! -d "${dir_work}/bwa" ]]; then
    mkdir -p "${dir_work}/bwa/err_out"
fi

#  Copy the FASTA info the directory for storing BWA indices
if [[ ! -f "${dir_work}/bwa/$(basename ${file_fa})" ]]; then
    cp "${file_fa}" "${dir_work}/bwa/$(basename ${file_fa})"
fi

#  Relocate to the directory for storing BWA indices
if [[ "$(pwd)" != "${dir_work}/bwa" ]]; then
    cd "${dir_work}/bwa" \
        || error_and_return "cd'ing failed; check on this."
fi

#  Using a HEREDOC, submit to SLURM a job for creating BWA indices
sbatch << EOF
#!/bin/bash

#SBATCH --job-name="${job_bwa}"
#SBATCH --nodes=1
#SBATCH --time=${time}
#SBATCH --error="${dir_work}/bwa/err_out/${job_bwa}.%A.stderr.txt"
#SBATCH --output="${dir_work}/bwa/err_out/${job_bwa}.%A.stdout.txt"

bwa index "$(basename ${file_fa})"
EOF
```
</details>
<br />

<a id="d-use-bowtie2-to-align-the-trimmed-fastq-files"></a>
## d. Use `bowtie2` to align the trimmed FASTQ files
<a id="code-6"></a>
### Code
<details>
<summary><i>Code: 3.d. Use `bowtie2` to align the trimmed FASTQ files</i></summary>

```bash
#!/bin/bash

#  Define function ============================================================
#  Function to return an error message and exit code 1, which stops the
#+ interactive execution of code
function error_and_return() {
    echo "Error: ${1}" >&2
    return 1
}


#  Function to check if Mamba is installed
function check_mamba_installed() {
    if ! : mamba &> /dev/null; then
        echo "Mamba is not installed on your system. Mamba is a package manager" 
        echo "that makes package installations faster and more reliable."
        echo ""
        echo "For installation instructions, please check the following link:"
        echo "https://github.com/mamba-org/mamba#installation"
        return 1
    fi
    
    return 0
}


#  Function to deactivate a Conda/Mamba environment
function deactivate_env() {
    if [[ "${CONDA_DEFAULT_ENV}" != "base" ]]; then
        if ! mamba deactivate &> /dev/null; then
            if ! conda deactivate &> /dev/null; then
                if ! source deactivate &> /dev/null; then
                    error_and_return "Failed to deactivate environment."
                    return 1
                fi
            fi
        fi
    fi

    return 0
}


#  Function to check if a specific Conda/Mamba environment is installed
function check_env_installed() {
    local env_name="${1}"

    if conda env list | grep -q "^${env_name} "; then
        return 0
    else
        echo "Environment \"${env_name}\" is not installed."
        return 1
    fi
}


#  Function to activate a specific Conda/Mamba environment
function activate_env() {
    local env_name="${1}"

    if ! mamba activate "${env_name}" &> /dev/null; then
        if ! conda activate "${env_name}" &> /dev/null; then
            if ! source activate "${env_name}" &> /dev/null; then
                error_and_return "Failed to activate environment \"${env_name}\"."
                return 1
            fi
        fi
    fi
    
    echo "Environment \"${env_name}\" activated successfully."
    return 0
}


#  Initialize variables and arrays ============================================
dir_base="${HOME}/tsukiyamalab"                          # Base directory for lab data
dir_repo="Kris/2023_tutorial_ChIP-seq"                   # Repository directory
dir_untr="01_sym"                                        # Directory for initial, non-trimmed FASTQs
dir_trim="02_trim"                                       # Directory for trimmed FASTQs
dir_bwt2="03_bam/bowtie2"
dir_genome=${HOME}/genomes/combined_SC_SP
dir_indx="${dir_genome}/bowtie2/combined_SC_SP"
file_fasta="${dir_genome}/fasta/combined_SC_SP.fa"
file_sizes="${dir_genome}/fasta/combined_SC_SP.chrom-info.tsv"

mapq=1

time="8:00:00"                                           # Job time for SLURM
threads=8                                                # Number of threads for SLURM jobs

#  Initialize an indexed array with FASTQ file stems
unset file_fastqs && typeset -a file_fastqs=(
    # "${dir_trim}/in_G1_Hho1_6336"
    # "${dir_trim}/IP_G1_Hho1_6336"
    # "${dir_trim}/in_G2M_Hho1_6336"
    # "${dir_trim}/IP_G2M_Hho1_6336"
    # "${dir_trim}/in_Q_Hho1_6336"
    # "${dir_trim}/IP_Q_Hho1_6336"
    # "${dir_trim}/in_G1_Hho1_6337"
    # "${dir_trim}/IP_G1_Hho1_6337"
    # "${dir_trim}/in_G2M_Hho1_6337"
    # "${dir_trim}/IP_G2M_Hho1_6337"
    # "${dir_trim}/in_Q_Hho1_6337"
    # "${dir_trim}/IP_Q_Hho1_6337"
    # "${dir_trim}/in_G1_Hmo1_7750"
    # "${dir_trim}/IP_G1_Hmo1_7750"
    # "${dir_trim}/in_G2M_Hmo1_7750"
    # "${dir_trim}/IP_G2M_Hmo1_7750"
    # "${dir_trim}/in_Q_Hmo1_7750"
    # "${dir_trim}/IP_Q_Hmo1_7750"
    # "${dir_trim}/in_G1_Hmo1_7751"
    # "${dir_trim}/IP_G1_Hmo1_7751"
    # "${dir_trim}/in_G2M_Hmo1_7751"
    # "${dir_trim}/IP_G2M_Hmo1_7751"
    # "${dir_trim}/in_Q_Hmo1_7751"
    # "${dir_trim}/IP_Q_Hmo1_7751"
    # "${dir_untr}/in_Q_untagged_5781"
    # "${dir_untr}/IP_Q_untagged_5781"
    # "${dir_untr}/in_Q_Esa5_7041"
    # "${dir_untr}/IP_Q_Esa5_7041"
    # "${dir_untr}/in_Q_Esa5_7691"
    # "${dir_untr}/IP_Q_Esa5_7691"
    # "${dir_untr}/in_Q_Rpd3_7568"
    # "${dir_untr}/IP_Q_Rpd3_7568"
    # "${dir_untr}/in_Q_Rpd3_7569"
    # "${dir_untr}/IP_Q_Rpd3_7569"
    # "${dir_untr}/in_Q_Gcn5_7692"
    # "${dir_untr}/IP_Q_Gcn5_7692"
    # "${dir_untr}/in_Q_Gcn5_7709"
    # "${dir_untr}/IP_Q_Gcn5_7709"
    "${dir_untr}/IP_log_Brn1_rep1"
    "${dir_untr}/IP_log_Brn1_rep2"
    "${dir_untr}/IP_log_Brn1_rep3"
    "${dir_untr}/IP_log_Brn1_repM"
    "${dir_untr}/IP_Q_Brn1_rep1"
    "${dir_untr}/IP_Q_Brn1_rep2"
    "${dir_untr}/IP_Q_Brn1_rep3"
    "${dir_untr}/IP_Q_Brn1_repM"
    "${dir_untr}/in_log_Brn1_rep1"
    "${dir_untr}/in_log_Brn1_rep2"
    "${dir_untr}/in_log_Brn1_rep3"
    "${dir_untr}/in_log_Brn1_repM"
    "${dir_untr}/in_Q_Brn1_rep1"
    "${dir_untr}/in_Q_Brn1_rep2"
    "${dir_untr}/in_Q_Brn1_rep3"
    "${dir_untr}/in_Q_Brn1_repM"
)


#  Do the main work ===========================================================
#  Set flags for checking variable and array assignments, and to create the
#+ necessary mamba environment if not found
check_variables=true
check_array=true
create_mamba_env=false

#  If check_variables is true, then echo the variable assignments
if ${check_variables}; then
    echo "
    dir_base=${dir_base}
    dir_repo=${dir_repo}
    dir_untr=${dir_untr}
    dir_trim=${dir_trim}
    dir_bwt2=${dir_bwt2}
    dir_genome=${dir_genome}
    dir_indx=${dir_indx}
    file_fasta=${file_fasta}
    file_sizes=${file_sizes}
    
    mapq=${mapq}

    time=${time}
    threads=${threads}
    "
fi

#  Echo array contents if check_array is true
if ${check_array}; then
    for i in "${!file_fastqs[@]}"; do
        file="${file_fastqs[${i}]}"

        if [[ "$(dirname ${file})" == "${dir_trim}" ]]; then
            if [[ "$(basename ${file})" =~ "Brn1" ]]; then
                echo "read ...... ${file}.atria.fastq.gz"
            else
                echo "
                read #1 ...... ${file}_R1.atria.fastq.gz
                read #2 ...... ${file}_R2.atria.fastq.gz
                "
            fi
        elif [[ "$(dirname ${file})" == "${dir_untr}" ]]; then
            if [[ "$(basename ${file})" =~ "Brn1" ]]; then
                echo "read ...... ${file}.fastq.gz"
            else
                echo "
                read #1 ...... ${file}_R1.fastq.gz
                read #2 ...... ${file}_R2.fastq.gz
                "
            fi
        else
            error_and_return "Processing logic problem for ${file}."
        fi
    done
fi

#  Initialize conda/mamba environment containing necessary programs for
#+ alignment, quality checks, and post-processing
env_name="alignment-processing_env"

#  Check that Mamba is installed and in PATH
check_mamba_installed

#  If not in base environment, then deactivate current environment
deactivate_env

#  Check that environment assigned to env_name is installed; if environment
#+ assigned to env_name is not installed, run the following invocation of mamba
#+ to install it
if check_env_installed "${env_name}"; then
    #  Handle the case when the environment is already installed
    echo "Activating environment ${env_name}"
    
    activate_env "${env_name}"
    module load Java/17.0.6  #TEMP Until install up-to-date Java with Mamba
else
    #  Handle the case when the environment is not installed
    echo "Creating environment ${env_name}"
    
    if ${create_mamba_env}; then
        #  Switch `--yes` is not set, which means user input is required
        #NOTE Running this on FHCC Rhino; ergo, no CONDA_SUBDIR=osx-64
        mamba create \
            --name "${env_name}" \
            --channel bioconda \
                bamtools \
                bedtools \
                bowtie2 \
                bwa \
                fastqc \
                minimap \
                mosdepth \
                picard \
                preseq \
                samtools \
                ucsc-bedgraphtobigwig \
                ucsc-bedsort \
                ucsc-facount
        
        activate_env "${env_name}"
        module load Java/17.0.6  #TEMP Until install up-to-date Java with Mamba
    fi
fi

#  Navigate to the work directory
cd "${dir_base}/${dir_repo}" \
    || error_and_return "Failed to cd to ${dir_base}/${dir_repo}."

#  If it doesn't exist, create a directory to store Bowtie2-aligned BAM files
if [[ ! -d "${dir_bwt2}" ]]; then
    mkdir -p "${dir_bwt2}/"{bam,cvrg,err_out,qc,siQ-ChIP}
fi

#  Set flags: checking variables, checking and submitting Bowtie2 jobs
print_iteration=true
check_variables=true
check_operation=true
run_operation=true

for i in "${!file_fastqs[@]}"; do
    # i=0
    index="${i}"
    iter=$(( index + 1 ))
    file="${file_fastqs[${index}]}"                          # echo "${file}"
    stem="$(basename ${file})"                               # echo "${stem}"
    job_name="$(echo ${dir_bwt2} | sed 's:\/:_:g').${stem}"  # echo "${job_name}"
    script="align-process-etc_fastqs_bowtie2.sh"             #TODO Move out of loop
    
    #  Parse the files' source directory to determine which alignment script to
    #+ use and appropriate suffix for FASTQ file names
    if [[ "$(dirname ${file})" == "${dir_trim}" ]]; then
        if [[ "$(basename ${file})" =~ "Brn1" ]]; then
            mode="single"                                    # echo "${mode}"
            fastq_1="${file}.atria.fastq.gz"                 # echo "${fastq_1}"
            fastq_2="ignore"                                 # echo "${fastq_2}"
        else
            mode="paired"                                    # echo "${mode}"
            fastq_1="${file}_R1.atria.fastq.gz"              # echo "${fastq_1}"
            fastq_2="${file}_R2.atria.fastq.gz"              # echo "${fastq_2}"
        fi
    elif [[ "$(dirname ${file})" == "${dir_untr}" ]]; then
        if [[ "$(basename ${file})" =~ "Brn1" ]]; then
            mode="single"                                    # echo "${mode}"
            fastq_1="${file}.fastq.gz"                       # echo "${fastq_1}"
            fastq_2="ignore"                                 # echo "${fastq_2}"
        else
            mode="paired"                                    # echo "${mode}"
            fastq_1="${file}_R1.fastq.gz"                    # echo "${fastq_1}"
            fastq_2="${file}_R2.fastq.gz"                    # echo "${fastq_2}"
        fi
    else
        error_and_return "Processing logic problem for ${file}."
    fi
    
    bam="${dir_bwt2}/bam/${stem}.bam"
    bam_coor="${bam/.bam/.sort-coord.bam}"
    bam_quer="${bam/.bam/.sort-qname.bam}"
    
    if [[ "${mode}" == "paired" ]]; then
        bed_siQ="${dir_bwt2}/siQ-ChIP/${stem}.bed.gz"
    else
        bed_siQ="ignore"
    fi
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

        script=${script}
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

        dir_base=${dir_base}
        dir_repo=${dir_repo}
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

    if ${check_operation}; then
        echo "
sbatch << EOF
#!/bin/bash

#SBATCH --job-name=\"${job_name}\"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=${threads}
#SBATCH --time=${time}
#SBATCH --error=\"${dir_bwt2}/err_out/${job_name}.%A.stderr.txt\"
#SBATCH --output=\"${dir_bwt2}/err_out/${job_name}.%A.stdout.txt\"

bash ${script} \\
    --threads \"${threads}\" \\
    --index \"${dir_indx}\" \\
    --fasta \"${file_fasta}\" \\
    --sizes \"${file_sizes}\" \\
    --mapq \"${mapq}\" \\
    --mode \"${mode}\" \\
    --fastq_1 \"${fastq_1}\" \\
    --fastq_2 \"${fastq_2}\" \\
    --bam \"${bam}\" \\
    --bam_coor \"${bam_coor}\" \\
    --bam_quer \"${bam_quer}\" \\
    --bed_siQ \"${bed_siQ}\" \\
    --bed_etc \"${bed_etc}\" \\
    --txt_met \"${txt_met}\" \\
    --txt_flg \"${txt_flg}\" \\
    --txt_idx \"${txt_idx}\" \\
    --txt_pre \"${txt_pre}\"
EOF
        "
    fi

    if ${run_operation}; then
sbatch << EOF
#!/bin/bash

#SBATCH --job-name="${job_name}"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=${threads}
#SBATCH --time=${time}
#SBATCH --error="${dir_bwt2}/err_out/${job_name}.%A.stderr.txt"
#SBATCH --output="${dir_bwt2}/err_out/${job_name}.%A.stdout.txt"

bash ${script} \
    --threads "${threads}" \
    --index "${dir_indx}" \
    --fasta "${file_fasta}" \
    --sizes "${file_sizes}" \
    --mapq "${mapq}" \
    --mode "${mode}" \
    --fastq_1 "${fastq_1}" \
    --fastq_2 "${fastq_2}" \
    --bam "${bam}" \
    --bam_coor "${bam_coor}" \
    --bam_quer "${bam_quer}" \
    --bed_siQ "${bed_siQ}" \
    --bed_etc "${bed_etc}" \
    --txt_met "${txt_met}" \
    --txt_flg "${txt_flg}" \
    --txt_idx "${txt_idx}" \
    --txt_pre "${txt_pre}"
EOF
    fi

    sleep 0.2
done

#TODO Initial BAM outfiles (from just after alignment) wern't deleted upon completion; debug this
```
</details>
<br />

<a id="e-use-bwa-to-align-the-trimmed-fastq-files"></a>
## e. Use `bwa` to align the trimmed FASTQ files
`#TODO`

<a id="code-7"></a>
### Code
<details>
<summary><i>Code: 3.e. Use `bwa` to align the trimmed FASTQ files</i></summary>

```bash
#!/bin/bash

echo "
bwa mem \\
    -t ${threads} \\
    -k 19 -w 100 -d 100 -r 1.5 \\
    -M \\
    -I 10,700 \\
    \"${f_indices}\" \\
    \"${stem}_R1.atria.fastq.gz\" \\
    \"${stem}_R2.atria.fastq.gz\" \\
        | samtools sort \\
            -@ ${threads} \\
            -T \"${scratch}\" \\
            -O bam \\
            -o \"${out_bam}\"

#  Index the sorted bams
if [[ -f \"${out_bam}\" ]]; then
    samtools index \\
        -@ ${threads} \\
        \"${out_bam}\"
fi
"

bwa mem \
    -t ${threads} \
    -k 19 -w 100 -d 100 -r 1.5 \
    -M \
    -I 10,700 \
    "${f_indices}" \
    "${stem}_R1.atria.fastq.gz" \
    "${stem}_R2.atria.fastq.gz" \
        | samtools sort \
            -@ ${threads} \
            -O bam \
            -o "${out_bam}"


#  Index the sorted bams
if [[ -f "${out_bam}" ]]; then
    samtools index \
        -@ ${threads} \
        "${out_bam}"
fi
```
</details>
<br />
<br />

<a id="x-run-phantompeakqualtools-on-aligned-data"></a>
# X. Run `phantompeakqualtools` on aligned data
<a id="a-install-phantompeakqualtools"></a>
## a. Install `phantompeakqualtools`
<a id="code-8"></a>
### Code
<details>
<summary><i>Code: X.a. Install `phantompeakqualtools`</i></summary>

```bash
#!/bin/bash

#  Define functions ===========================================================
#  Function to return an error message and exit code 1, which stops the
#+ interactive execution of code
function error_and_return() {
    echo "Error: ${1}" >&2
    return 1
}


#  Function to check if Mamba is installed
function check_mamba_installed() {
    if ! : mamba &> /dev/null; then
        echo "Mamba is not installed on your system. Mamba is a package manager" 
        echo "that makes package installations faster and more reliable."
        echo ""
        echo "For installation instructions, please check the following link:"
        echo "https://github.com/mamba-org/mamba#installation"
        return 1
    fi
    
    return 0
}


#  Function to deactivate a Conda/Mamba environment
function deactivate_env() {
    if [[ "${CONDA_DEFAULT_ENV}" != "base" ]]; then
        if ! mamba deactivate &> /dev/null; then
            if ! conda deactivate &> /dev/null; then
                if ! source deactivate &> /dev/null; then
                    error_and_return "Failed to deactivate environment."
                    return 1
                fi
            fi
        fi
    fi

    return 0
}


#  Function to check if a specific Conda/Mamba environment is installed
function check_env_installed() {
    local env_name="${1}"

    if conda env list | grep -q "^${env_name} "; then
        return 0
    else
        echo "Environment \"${env_name}\" is not installed."
        return 1
    fi
}


#  Function to activate a specific Conda/Mamba environment
function activate_env() {
    local env_name="${1}"

    if ! mamba activate "${env_name}" &> /dev/null; then
        if ! conda activate "${env_name}" &> /dev/null; then
            if ! source activate "${env_name}" &> /dev/null; then
                error_and_return "Failed to activate environment \"${env_name}\"."
                return 1
            fi
        fi
    fi
    
    echo "Environment \"${env_name}\" activated successfully."
    return 0
}


#  Initialize variables and arrays ============================================
env_name="ppqt_env"


#  Do the main work ===========================================================
#  Set flag(s)
create_mamba_env=true  # Install mamba environment if not detected

#  Check that Mamba is installed and in PATH
check_mamba_installed

#  If not in base environment, then deactivate current environment
deactivate_env

#  Check that environment assigned to env_name is installed; if environment
#+ assigned to env_name is not installed, run the following
if check_env_installed "${env_name}"; then
    #  Handle the case when the environment is already installed
    echo "Activating environment ${env_name}"
    
    activate_env "${env_name}"
else
    #  Handle the case when the environment is not installed
    echo "Creating environment ${env_name}"
    
    if ${create_mamba_env}; then
        #  Switch `--yes` is set, which means no user input is required
        #NOTE Running this on FHCC Rhino; ergo, no CONDA_SUBDIR=osx-64
        mamba create \
            --yes \
            -n "${env_name}" \
            -c bioconda \
            -c conda-forge \
                phantompeakqualtools=1.2.2 \
                parallel
    fi
fi
```
</details>
<br />

<a id="printed-remote"></a>
### Printed (remote)
<details>
<summary><i>Printed: X.a. Install `phantompeakqualtools` (remote)</i></summary>

```txt
❯ if check_env_installed "${env_name}"; then
>     #  Handle the case when the environment is already installed
>     echo "Activating environment ${env_name}"
> 
>     activate_env "${env_name}"
> else
>     #  Handle the case when the environment is not installed
>     echo "Creating environment ${env_name}"
> 
>     if ${create_mamba_env}; then
>         #  Switch `--yes` is set, which means no user input is required
>         #NOTE Running this on FHCC Rhino; ergo, no CONDA_SUBDIR=osx-64
>         mamba create \
>             --yes \
>             -n "${env_name}" \
>             -c bioconda \
>                 phantompeakqualtools=1.2.2
>     fi
> fi
Environment "ppqt_env" is not installed.
Creating environment ppqt_env

                  __    __    __    __
                 /  \  /  \  /  \  /  \
                /    \/    \/    \/    \
███████████████/  /██/  /██/  /██/  /████████████████████████
              /  / \   / \   / \   / \  \____
             /  /   \_/   \_/   \_/   \    o \__,
            / _/                       \_____/  `
            |/
        ███╗   ███╗ █████╗ ███╗   ███╗██████╗  █████╗
        ████╗ ████║██╔══██╗████╗ ████║██╔══██╗██╔══██╗
        ██╔████╔██║███████║██╔████╔██║██████╔╝███████║
        ██║╚██╔╝██║██╔══██║██║╚██╔╝██║██╔══██╗██╔══██║
        ██║ ╚═╝ ██║██║  ██║██║ ╚═╝ ██║██████╔╝██║  ██║
        ╚═╝     ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝╚═════╝ ╚═╝  ╚═╝

        mamba (1.3.1) supported by @QuantStack

        GitHub:  https://github.com/mamba-org/mamba
        Twitter: https://twitter.com/QuantStack

█████████████████████████████████████████████████████████████


Looking for: ['phantompeakqualtools=1.2.2']

bioconda/linux-64                                           Using cache
bioconda/noarch                                             Using cache
conda-forge/linux-64                                        Using cache
conda-forge/noarch                                          Using cache
pkgs/main/linux-64                                            No change
pkgs/main/noarch                                              No change
pkgs/r/linux-64                                               No change
pkgs/r/noarch                                                 No change
Transaction

  Prefix: /home/kalavatt/miniconda3/envs/ppqt_env

  Updating specs:

   - phantompeakqualtools=1.2.2


  Package                               Version  Build                Channel                    Size
───────────────────────────────────────────────────────────────────────────────────────────────────────
  Install:
───────────────────────────────────────────────────────────────────────────────────────────────────────

  + _libgcc_mutex                           0.1  conda_forge          conda-forge/linux-64     Cached
  + _openmp_mutex                           4.5  2_gnu                conda-forge/linux-64     Cached
  + _r-mutex                              1.0.1  anacondar_1          conda-forge/noarch       Cached
  + argcomplete                           3.2.3  pyhd8ed1ab_0         conda-forge/noarch         40kB
  + binutils_impl_linux-64                 2.40  hf600244_0           conda-forge/linux-64     Cached
  + bioconductor-biocgenerics            0.48.1  r43hdfd78af_2        bioconda/noarch           658kB
  + bioconductor-biocparallel            1.36.0  r43hf17093f_1        bioconda/linux-64           2MB
  + bioconductor-biostrings              2.70.1  r43ha9d7317_1        bioconda/linux-64          15MB
  + bioconductor-data-packages         20231203  hdfd78af_0           bioconda/noarch          Cached
  + bioconductor-genomeinfodb            1.38.1  r43hdfd78af_1        bioconda/noarch             4MB
  + bioconductor-genomeinfodbdata        1.2.11  r43hdfd78af_1        bioconda/noarch          Cached
  + bioconductor-genomicranges           1.54.1  r43ha9d7317_1        bioconda/linux-64           2MB
  + bioconductor-iranges                 2.36.0  r43ha9d7317_1        bioconda/linux-64           3MB
  + bioconductor-rhtslib                  2.4.0  r43ha9d7317_1        bioconda/linux-64           2MB
  + bioconductor-rsamtools               2.18.0  r43hf17093f_1        bioconda/linux-64           4MB
  + bioconductor-s4vectors               0.40.2  r43ha9d7317_1        bioconda/linux-64           3MB
  + bioconductor-xvector                 0.42.0  r43ha9d7317_1        bioconda/linux-64         755kB
  + bioconductor-zlibbioc                1.48.0  r43ha9d7317_1        bioconda/linux-64          26kB
  + boost                                1.84.0  h8da182e_2           conda-forge/linux-64       16kB
  + bwidget                              1.9.14  ha770c72_1           conda-forge/linux-64     Cached
  + bzip2                                 1.0.8  hd590300_5           conda-forge/linux-64     Cached
  + c-ares                               1.28.1  hd590300_0           conda-forge/linux-64      169kB
  + ca-certificates                    2024.2.2  hbcca054_0           conda-forge/linux-64     Cached
  + cairo                                1.18.0  h3faef2a_0           conda-forge/linux-64     Cached
  + curl                                  8.7.1  hca28451_0           conda-forge/linux-64      165kB
  + expat                                 2.6.2  h59595ed_0           conda-forge/linux-64      138kB
  + font-ttf-dejavu-sans-mono              2.37  hab24e00_0           conda-forge/noarch       Cached
  + font-ttf-inconsolata                  3.000  h77eed37_0           conda-forge/noarch       Cached
  + font-ttf-source-code-pro              2.038  h77eed37_0           conda-forge/noarch       Cached
  + font-ttf-ubuntu                        0.83  h77eed37_1           conda-forge/noarch       Cached
  + fontconfig                           2.14.2  h14ed4e7_0           conda-forge/linux-64     Cached
  + fonts-conda-ecosystem                     1  0                    conda-forge/noarch       Cached
  + fonts-conda-forge                         1  0                    conda-forge/noarch       Cached
  + freetype                             2.12.1  h267a509_2           conda-forge/linux-64     Cached
  + fribidi                              1.0.10  h36c2ea0_0           conda-forge/linux-64     Cached
  + gawk                                  5.3.0  ha916aea_0           conda-forge/linux-64     Cached
  + gcc_impl_linux-64                    13.2.0  h338b0a0_5           conda-forge/linux-64     Cached
  + gettext                              0.21.1  h27087fc_0           conda-forge/linux-64     Cached
  + gfortran_impl_linux-64               13.2.0  h76e1118_5           conda-forge/linux-64     Cached
  + gmp                                   6.3.0  h59595ed_1           conda-forge/linux-64      570kB
  + graphite2                            1.3.13  h59595ed_1003        conda-forge/linux-64       97kB
  + gxx_impl_linux-64                    13.2.0  h338b0a0_5           conda-forge/linux-64     Cached
  + harfbuzz                              8.3.0  h3d44ed6_0           conda-forge/linux-64     Cached
  + icu                                    73.2  h59595ed_0           conda-forge/linux-64     Cached
  + jq                                      1.5  0                    bioconda/linux-64        Cached
  + kernel-headers_linux-64              2.6.32  he073ed8_17          conda-forge/noarch       Cached
  + keyutils                              1.6.1  h166bdaf_0           conda-forge/linux-64     Cached
  + krb5                                 1.21.2  h659d440_0           conda-forge/linux-64     Cached
  + ld_impl_linux-64                       2.40  h41732ed_0           conda-forge/linux-64     Cached
  + lerc                                  4.0.0  h27087fc_0           conda-forge/linux-64     Cached
  + libblas                               3.9.0  21_linux64_openblas  conda-forge/linux-64     Cached
  + libboost                             1.84.0  h8013b2b_2           conda-forge/linux-64        3MB
  + libboost-devel                       1.84.0  h00ab1b0_2           conda-forge/linux-64       39kB
  + libboost-headers                     1.84.0  ha770c72_2           conda-forge/linux-64       14MB
  + libboost-python                      1.84.0  py312hfb10629_2      conda-forge/linux-64      123kB
  + libboost-python-devel                1.84.0  py312h8da182e_2      conda-forge/linux-64       19kB
  + libcblas                              3.9.0  21_linux64_openblas  conda-forge/linux-64     Cached
  + libcurl                               8.7.1  hca28451_0           conda-forge/linux-64      398kB
  + libdeflate                             1.20  hd590300_0           conda-forge/linux-64       72kB
  + libedit                        3.1.20191231  he28a2e2_2           conda-forge/linux-64     Cached
  + libev                                  4.33  hd590300_2           conda-forge/linux-64     Cached
  + libexpat                              2.6.2  h59595ed_0           conda-forge/linux-64       74kB
  + libffi                                3.4.2  h7f98852_5           conda-forge/linux-64     Cached
  + libgcc                                7.2.0  h69d50b8_2           conda-forge/linux-64     Cached
  + libgcc-devel_linux-64                13.2.0  ha9c7c90_105         conda-forge/noarch       Cached
  + libgcc-ng                            13.2.0  h807b86a_5           conda-forge/linux-64     Cached
  + libgfortran-ng                       13.2.0  h69a702a_5           conda-forge/linux-64     Cached
  + libgfortran5                         13.2.0  ha4646dd_5           conda-forge/linux-64     Cached
  + libglib                              2.78.4  h783c2da_0           conda-forge/linux-64        3MB
  + libgomp                              13.2.0  h807b86a_5           conda-forge/linux-64     Cached
  + libiconv                               1.17  hd590300_2           conda-forge/linux-64     Cached
  + libjpeg-turbo                         3.0.0  hd590300_1           conda-forge/linux-64     Cached
  + liblapack                             3.9.0  21_linux64_openblas  conda-forge/linux-64     Cached
  + libnghttp2                           1.58.0  h47da74e_1           conda-forge/linux-64     Cached
  + libnsl                                2.0.1  hd590300_0           conda-forge/linux-64     Cached
  + libopenblas                          0.3.26  pthreads_h413a1c8_0  conda-forge/linux-64     Cached
  + libpng                               1.6.43  h2797004_0           conda-forge/linux-64     Cached
  + libsanitizer                         13.2.0  h7e041cc_5           conda-forge/linux-64     Cached
  + libsqlite                            3.45.2  h2797004_0           conda-forge/linux-64      857kB
  + libssh2                              1.11.0  h0841786_0           conda-forge/linux-64     Cached
  + libstdcxx-devel_linux-64             13.2.0  ha9c7c90_105         conda-forge/noarch       Cached
  + libstdcxx-ng                         13.2.0  h7e041cc_5           conda-forge/linux-64     Cached
  + libtiff                               4.6.0  h1dd3fc0_3           conda-forge/linux-64      283kB
  + libuuid                              2.38.1  h0b41bf4_0           conda-forge/linux-64     Cached
  + libwebp-base                          1.3.2  hd590300_0           conda-forge/linux-64     Cached
  + libxcb                                 1.15  h0b41bf4_0           conda-forge/linux-64     Cached
  + libxcrypt                            4.4.36  hd590300_1           conda-forge/linux-64     Cached
  + libzlib                              1.2.13  hd590300_5           conda-forge/linux-64     Cached
  + make                                    4.3  hd18ef5c_1           conda-forge/linux-64     Cached
  + mpfr                                  4.2.1  h9458935_0           conda-forge/linux-64     Cached
  + ncurses                        6.4.20240210  h59595ed_0           conda-forge/linux-64      896kB
  + numpy                                1.26.4  py312heda63a1_0      conda-forge/linux-64        7MB
  + openssl                               3.2.1  hd590300_1           conda-forge/linux-64        3MB
  + pango                                1.52.1  ha41ecd1_0           conda-forge/linux-64      444kB
  + pcre2                                 10.42  hcad00b1_0           conda-forge/linux-64     Cached
  + phantompeakqualtools                  1.2.2  hdfd78af_1           bioconda/noarch          Cached
  + pip                                    24.0  pyhd8ed1ab_0         conda-forge/noarch       Cached
  + pixman                               0.43.2  h59595ed_0           conda-forge/linux-64     Cached
  + pthread-stubs                           0.4  h36c2ea0_1001        conda-forge/linux-64     Cached
  + python                               3.12.2  hab00c5b_0_cpython   conda-forge/linux-64       32MB
  + python_abi                             3.12  4_cp312              conda-forge/linux-64     Cached
  + pyyaml                                6.0.1  py312h98912ed_1      conda-forge/linux-64      197kB
  + r-base                                4.3.3  hb8ee39d_0           conda-forge/linux-64       26MB
  + r-bh                               1.84.0_0  r43hc72bb7e_0        conda-forge/noarch       Cached
  + r-bitops                              1.0_7  r43h57805ef_2        conda-forge/linux-64     Cached
  + r-catools                            1.18.2  r43ha503ecb_2        conda-forge/linux-64     Cached
  + r-codetools                          0.2_20  r43hc72bb7e_0        conda-forge/noarch        108kB
  + r-cpp11                               0.4.7  r43hc72bb7e_0        conda-forge/noarch       Cached
  + r-crayon                              1.5.2  r43hc72bb7e_2        conda-forge/noarch       Cached
  + r-formatr                              1.14  r43hc72bb7e_1        conda-forge/noarch       Cached
  + r-futile.logger                       1.4.3  r43hc72bb7e_1005     conda-forge/noarch       Cached
  + r-futile.options                      1.0.1  r43hc72bb7e_1004     conda-forge/noarch       Cached
  + r-lambda.r                            1.2.4  r43hc72bb7e_3        conda-forge/noarch       Cached
  + r-rcpp                               1.0.12  r43h7df8631_0        conda-forge/linux-64     Cached
  + r-rcurl                           1.98_1.14  r43hf9611b0_0        conda-forge/linux-64     Cached
  + r-snow                                0.4_4  r43hc72bb7e_2        conda-forge/noarch       Cached
  + r-snowfall                         1.84_6.3  r43hc72bb7e_0        conda-forge/noarch       Cached
  + r-spp                                1.16.0  r43h21a89ab_9        bioconda/linux-64        Cached
  + readline                                8.2  h8228510_1           conda-forge/linux-64     Cached
  + samtools                                1.6  hc3601fc_10          bioconda/linux-64         513kB
  + sed                                     4.8  he412f7d_0           conda-forge/linux-64     Cached
  + setuptools                           69.2.0  pyhd8ed1ab_0         conda-forge/noarch        471kB
  + sysroot_linux-64                       2.12  he073ed8_17          conda-forge/noarch       Cached
  + tk                                   8.6.13  noxft_h4845f30_101   conda-forge/linux-64     Cached
  + tktable                                2.10  h0c5db8f_5           conda-forge/linux-64     Cached
  + toml                                 0.10.2  pyhd8ed1ab_0         conda-forge/noarch       Cached
  + tomlkit                              0.12.4  pyha770c72_0         conda-forge/noarch       Cached
  + tzdata                                2024a  h0c530f3_0           conda-forge/noarch       Cached
  + wheel                                0.43.0  pyhd8ed1ab_1         conda-forge/noarch         58kB
  + xmltodict                            0.13.0  pyhd8ed1ab_0         conda-forge/noarch       Cached
  + xorg-kbproto                          1.0.7  h7f98852_1002        conda-forge/linux-64     Cached
  + xorg-libice                           1.1.1  hd590300_0           conda-forge/linux-64     Cached
  + xorg-libsm                            1.2.4  h7391055_0           conda-forge/linux-64     Cached
  + xorg-libx11                           1.8.7  h8ee46fc_0           conda-forge/linux-64     Cached
  + xorg-libxau                          1.0.11  hd590300_0           conda-forge/linux-64     Cached
  + xorg-libxdmcp                         1.1.3  h7f98852_0           conda-forge/linux-64     Cached
  + xorg-libxext                          1.3.4  h0b41bf4_2           conda-forge/linux-64     Cached
  + xorg-libxrender                      0.9.11  hd590300_0           conda-forge/linux-64     Cached
  + xorg-libxt                            1.3.0  hd590300_1           conda-forge/linux-64     Cached
  + xorg-renderproto                     0.11.1  h7f98852_1002        conda-forge/linux-64     Cached
  + xorg-xextproto                        7.3.0  h0b41bf4_1003        conda-forge/linux-64     Cached
  + xorg-xproto                          7.0.31  h7f98852_1007        conda-forge/linux-64     Cached
  + xz                                    5.2.6  h166bdaf_0           conda-forge/linux-64     Cached
  + yaml                                  0.2.5  h7f98852_2           conda-forge/linux-64     Cached
  + yq                                    3.2.3  pyhd8ed1ab_0         conda-forge/noarch       Cached
  + zlib                                 1.2.13  hd590300_5           conda-forge/linux-64     Cached
  + zstd                                  1.5.5  hfc55251_0           conda-forge/linux-64     Cached

  Summary:

  Install: 147 packages

  Total download: 130MB

───────────────────────────────────────────────────────────────────────────────────────────────────────


libexpat                                            73.7kB @ 584.1kB/s  0.1s
gmp                                                569.9kB @   4.4MB/s  0.1s
c-ares                                             168.9kB @   1.2MB/s  0.2s
graphite2                                           96.9kB @ 574.5kB/s  0.2s
libboost-devel                                      38.5kB @ 215.5kB/s  0.1s
libboost                                             2.8MB @  12.9MB/s  0.1s
r-codetools                                        108.2kB @ 470.7kB/s  0.1s
samtools                                           512.9kB @   2.0MB/s  0.1s
bioconductor-zlibbioc                               25.6kB @  92.6kB/s  0.1s
libboost-python                                    123.3kB @ 332.5kB/s  0.1s
bioconductor-rhtslib                                 2.4MB @   5.2MB/s  0.3s
libboost-headers                                    13.7MB @  23.8MB/s  0.6s
bioconductor-genomicranges                           2.3MB @   3.7MB/s  0.3s
expat                                              137.6kB @ 211.7kB/s  0.1s
numpy                                                7.5MB @  10.5MB/s  0.5s
libtiff                                            282.7kB @ 385.8kB/s  0.1s
bioconductor-rsamtools                               4.2MB @   5.7MB/s  0.3s
setuptools                                         471.2kB @ 551.4kB/s  0.1s
argcomplete                                         40.2kB @  47.0kB/s  0.1s
libglib                                              2.7MB @   3.1MB/s  0.2s
boost                                               16.0kB @  15.3kB/s  0.2s
libsqlite                                          857.5kB @ 804.3kB/s  0.2s
bioconductor-iranges                                 2.6MB @   2.4MB/s  0.4s
r-base                                              25.7MB @  23.0MB/s  1.0s
bioconductor-genomeinfodb                            4.4MB @   3.9MB/s  0.3s
pango                                              444.2kB @ 361.3kB/s  0.2s
bioconductor-xvector                               755.5kB @ 606.7kB/s  0.1s
bioconductor-s4vectors                               2.6MB @   2.0MB/s  0.2s
wheel                                               58.0kB @  42.0kB/s  0.2s
pyyaml                                             196.6kB @ 141.7kB/s  0.1s
bioconductor-biocparallel                            1.7MB @   1.1MB/s  0.3s
ncurses                                            895.7kB @ 582.1kB/s  0.2s
openssl                                              2.9MB @   1.7MB/s  0.3s
libdeflate                                          71.5kB @  41.8kB/s  0.3s
libboost-python-devel                               19.5kB @  11.4kB/s  0.4s
libcurl                                            398.3kB @ 210.7kB/s  0.2s
curl                                               164.9kB @  87.1kB/s  0.2s
bioconductor-biocgenerics                          658.4kB @ 314.9kB/s  0.2s
bioconductor-biostrings                             14.6MB @   6.4MB/s  1.2s
python                                              32.3MB @  12.5MB/s  1.6s

Downloading and Extracting Packages

Preparing transaction: done
Verifying transaction: done
Executing transaction: done

To activate this environment, use

     $ mamba activate ppqt_env

To deactivate an active environment, use

     $ mamba deactivate
```
</details>
<br />

<a id="printed-local"></a>
### Printed (local)
<details>
<summary><i>Printed: X.a. Install `phantompeakqualtools` (local)</i></summary>

```txt
┌─[kalavattam][Kriss-MacBook-Pro][~]
└─▪  if check_env_installed "${env_name}"; then
└─▪     #  Handle the case when the environment is already installed
└─▪     echo "Activating environment ${env_name}"
└─▪
└─▪     activate_env "${env_name}"
└─▪ else
└─▪     #  Handle the case when the environment is not installed
└─▪     echo "Creating environment ${env_name}"
└─▪
└─▪     if ${create_mamba_env}; then
└─▪         #  Switch `--yes` is set, which means no user input is required
└─▪         #NOTE Running this on FHCC Rhino; ergo, no CONDA_SUBDIR=osx-64
└─▪         mamba create \
└─▪             --yes \
└─▪             -n "${env_name}" \
└─▪             -c bioconda \
└─▪                 phantompeakqualtools=1.2.2 \
└─▪                 parallel
└─▪     fi
└─▪ fi
Environment "ppqt_env" is not installed.
Creating environment ppqt_env

                  __    __    __    __
                 /  \  /  \  /  \  /  \
                /    \/    \/    \/    \
███████████████/  /██/  /██/  /██/  /████████████████████████
              /  / \   / \   / \   / \  \____
             /  /   \_/   \_/   \_/   \    o \__,
            / _/                       \_____/  `
            |/
        ███╗   ███╗ █████╗ ███╗   ███╗██████╗  █████╗
        ████╗ ████║██╔══██╗████╗ ████║██╔══██╗██╔══██╗
        ██╔████╔██║███████║██╔████╔██║██████╔╝███████║
        ██║╚██╔╝██║██╔══██║██║╚██╔╝██║██╔══██╗██╔══██║
        ██║ ╚═╝ ██║██║  ██║██║ ╚═╝ ██║██████╔╝██║  ██║
        ╚═╝     ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝╚═════╝ ╚═╝  ╚═╝

        mamba (0.25.0) supported by @QuantStack

        GitHub:  https://github.com/mamba-org/mamba
        Twitter: https://twitter.com/QuantStack

█████████████████████████████████████████████████████████████


Looking for: ['phantompeakqualtools=1.2.2', 'parallel']

r/osx-64                                                      No change
r/noarch                                                      No change
bioconda/osx-64                                      4.4MB @   2.4MB/s  1.9s
pkgs/r/osx-64                                                 No change
bioconda/noarch                                      5.2MB @   2.4MB/s  2.3s
pkgs/r/noarch                                                 No change
pkgs/main/noarch                                              No change
conda-forge/noarch                                  16.4MB @   4.3MB/s  4.1s
pkgs/main/osx-64                                     6.5MB @   1.6MB/s  2.4s
conda-forge/osx-64                                  35.0MB @   5.4MB/s  7.2s
Transaction

  Prefix: /Users/kalavattam/miniconda3/envs/ppqt_env

  Updating specs:

   - phantompeakqualtools=1.2.2
   - parallel


  Package                               Version  Build                 Channel                  Size
──────────────────────────────────────────────────────────────────────────────────────────────────────
  Install:
──────────────────────────────────────────────────────────────────────────────────────────────────────

  + _r-mutex                              1.0.1  anacondar_1           conda-forge/noarch     Cached
  + argcomplete                           3.2.3  pyhd8ed1ab_0          conda-forge/noarch     Cached
  + bioconductor-biocgenerics            0.48.1  r43hdfd78af_2         bioconda/noarch        Cached
  + bioconductor-biocparallel            1.36.0  r43hc0ef7c4_1         bioconda/osx-64        Cached
  + bioconductor-biostrings              2.70.1  r43h4c50009_1         bioconda/osx-64        Cached
  + bioconductor-data-packages         20231203  hdfd78af_0            bioconda/noarch        Cached
  + bioconductor-genomeinfodb            1.38.1  r43hdfd78af_1         bioconda/noarch        Cached
  + bioconductor-genomeinfodbdata        1.2.11  r43hdfd78af_1         bioconda/noarch        Cached
  + bioconductor-genomicranges           1.54.1  r43h4c50009_1         bioconda/osx-64        Cached
  + bioconductor-iranges                 2.36.0  r43h4c50009_1         bioconda/osx-64        Cached
  + bioconductor-rhtslib                  2.4.0  r43h4c50009_1         bioconda/osx-64        Cached
  + bioconductor-rsamtools               2.18.0  r43hc0ef7c4_1         bioconda/osx-64        Cached
  + bioconductor-s4vectors               0.40.2  r43h4c50009_1         bioconda/osx-64        Cached
  + bioconductor-xvector                 0.42.0  r43h4c50009_1         bioconda/osx-64        Cached
  + bioconductor-zlibbioc                1.48.0  r43h4c50009_1         bioconda/osx-64        Cached
  + boost                                1.84.0  hdce95a9_2            conda-forge/osx-64     Cached
  + bwidget                              1.9.14  h694c41f_1            conda-forge/osx-64     Cached
  + bzip2                                 1.0.8  h10d778d_5            conda-forge/osx-64     Cached
  + c-ares                               1.28.1  h10d778d_0            conda-forge/osx-64     Cached
  + ca-certificates                    2024.2.2  h8857fd0_0            conda-forge/osx-64     Cached
  + cairo                                1.18.0  h99e66fa_0            conda-forge/osx-64     Cached
  + cctools_osx-64                          986  h58a35ae_0            conda-forge/osx-64     Cached
  + clang                                18.1.2  default_h7151d67_1    conda-forge/osx-64     Cached
  + clang-18                             18.1.2  default_h7151d67_1    conda-forge/osx-64     Cached
  + clang_impl_osx-64                    18.1.2  h73f7f27_11           conda-forge/osx-64     Cached
  + clang_osx-64                         18.1.2  hb91bd55_11           conda-forge/osx-64     Cached
  + clangxx                              18.1.2  default_h7151d67_1    conda-forge/osx-64     Cached
  + clangxx_impl_osx-64                  18.1.2  hb14bd1d_11           conda-forge/osx-64     Cached
  + clangxx_osx-64                       18.1.2  hb91bd55_11           conda-forge/osx-64     Cached
  + compiler-rt                          18.1.2  ha38d28d_0            conda-forge/osx-64     Cached
  + compiler-rt_osx-64                   18.1.2  ha38d28d_0            conda-forge/noarch     Cached
  + curl                                  8.7.1  h726d00d_0            conda-forge/osx-64     Cached
  + expat                                 2.6.2  h73e2aa4_0            conda-forge/osx-64     Cached
  + font-ttf-dejavu-sans-mono              2.37  hab24e00_0            conda-forge/noarch     Cached
  + font-ttf-inconsolata                  3.000  h77eed37_0            conda-forge/noarch     Cached
  + font-ttf-source-code-pro              2.038  h77eed37_0            conda-forge/noarch     Cached
  + font-ttf-ubuntu                        0.83  h77eed37_1            conda-forge/noarch     Cached
  + fontconfig                           2.14.2  h5bb23bf_0            conda-forge/osx-64     Cached
  + fonts-conda-ecosystem                     1  0                     conda-forge/noarch     Cached
  + fonts-conda-forge                         1  0                     conda-forge/noarch     Cached
  + freetype                             2.12.1  h60636b9_2            conda-forge/osx-64     Cached
  + fribidi                              1.0.10  hbcb3906_0            conda-forge/osx-64     Cached
  + gawk                                  5.3.0  h2c496e9_0            conda-forge/osx-64     Cached
  + gettext                              0.21.1  h8a4c099_0            conda-forge/osx-64     Cached
  + gfortran_impl_osx-64                 12.3.0  hc328e78_3            conda-forge/osx-64     Cached
  + gfortran_osx-64                      12.3.0  h18f7dce_1            conda-forge/osx-64     Cached
  + gmp                                   6.3.0  h73e2aa4_1            conda-forge/osx-64     Cached
  + graphite2                            1.3.13  h73e2aa4_1003         conda-forge/osx-64     Cached
  + harfbuzz                              8.3.0  hf45c392_0            conda-forge/osx-64     Cached
  + htslib                               1.19.1  h365c357_2            bioconda/osx-64        Cached
  + icu                                    73.2  hf5e326d_0            conda-forge/osx-64     Cached
  + isl                                    0.26  imath32_h2e86a7b_101  conda-forge/osx-64     Cached
  + jq                                      1.5  0                     bioconda/osx-64        Cached
  + krb5                                 1.21.2  hb884880_0            conda-forge/osx-64     Cached
  + ld64_osx-64                             711  had5d0d3_0            conda-forge/osx-64     Cached
  + lerc                                  4.0.0  hb486fe8_0            conda-forge/osx-64     Cached
  + libblas                               3.9.0  21_osx64_openblas     conda-forge/osx-64     Cached
  + libboost                             1.84.0  h6ebd1c4_2            conda-forge/osx-64     Cached
  + libboost-devel                       1.84.0  hf2b3138_2            conda-forge/osx-64     Cached
  + libboost-headers                     1.84.0  h694c41f_2            conda-forge/osx-64     Cached
  + libboost-python                      1.84.0  py312h77b368e_2       conda-forge/osx-64     Cached
  + libboost-python-devel                1.84.0  py312hdce95a9_2       conda-forge/osx-64     Cached
  + libcblas                              3.9.0  21_osx64_openblas     conda-forge/osx-64     Cached
  + libclang-cpp18.1                     18.1.2  default_h7151d67_1    conda-forge/osx-64     Cached
  + libcurl                               8.7.1  h726d00d_0            conda-forge/osx-64     Cached
  + libcxx                               16.0.6  hd57cbcb_0            conda-forge/osx-64     Cached
  + libdeflate                             1.18  hac1461d_0            conda-forge/osx-64     Cached
  + libedit                        3.1.20191231  h0678c8f_2            conda-forge/osx-64     Cached
  + libev                                  4.33  h10d778d_2            conda-forge/osx-64     Cached
  + libexpat                              2.6.2  h73e2aa4_0            conda-forge/osx-64     Cached
  + libffi                                3.4.2  h0d85af4_5            conda-forge/osx-64     Cached
  + libgcc                                4.8.5  1                     conda-forge/osx-64     Cached
  + libgfortran                           5.0.0  13_2_0_h97931a8_3     conda-forge/osx-64     Cached
  + libgfortran-devel_osx-64             12.3.0  h0b6f5ec_3            conda-forge/noarch     Cached
  + libgfortran5                         13.2.0  h2873a65_3            conda-forge/osx-64     Cached
  + libglib                              2.78.1  h6d9ecee_0            conda-forge/osx-64     Cached
  + libiconv                               1.17  hd75f5a5_2            conda-forge/osx-64     Cached
  + libjpeg-turbo                       2.1.5.1  h0dc2134_1            conda-forge/osx-64     Cached
  + liblapack                             3.9.0  21_osx64_openblas     conda-forge/osx-64     Cached
  + libllvm18                            18.1.2  hbcf5fad_0            conda-forge/osx-64     Cached
  + libnghttp2                           1.58.0  h64cf6d3_1            conda-forge/osx-64     Cached
  + libopenblas                          0.3.26  openmp_hfef2a42_0     conda-forge/osx-64     Cached
  + libpng                               1.6.43  h92b6c6a_0            conda-forge/osx-64     Cached
  + libsqlite                            3.45.2  h92b6c6a_0            conda-forge/osx-64     Cached
  + libssh2                              1.11.0  hd019ec5_0            conda-forge/osx-64     Cached
  + libtiff                               4.6.0  hf955e92_0            conda-forge/osx-64     Cached
  + libwebp-base                          1.3.2  h0dc2134_0            conda-forge/osx-64     Cached
  + libxml2                              2.12.6  hc0ae0f7_1            conda-forge/osx-64     Cached
  + libzlib                              1.2.13  h8a1eda9_5            conda-forge/osx-64     Cached
  + llvm-openmp                          18.1.2  hb6ac08f_0            conda-forge/osx-64     Cached
  + llvm-tools                           18.1.2  hbcf5fad_0            conda-forge/osx-64     Cached
  + make                                    4.3  h22f3db7_1            conda-forge/osx-64     Cached
  + mpc                                   1.3.1  h81bd1dd_0            conda-forge/osx-64     Cached
  + mpfr                                  4.2.1  h0c69b56_0            conda-forge/osx-64     Cached
  + ncurses                        6.4.20240210  h73e2aa4_0            conda-forge/osx-64     Cached
  + numpy                                1.26.4  py312he3a82b2_0       conda-forge/osx-64     Cached
  + openssl                               3.2.1  hd75f5a5_1            conda-forge/osx-64     Cached
  + pango                               1.50.14  h19c1c8a_2            conda-forge/osx-64     Cached
  + parallel                           20170422  pl5.22.0_0            bioconda/osx-64           1MB
  + pcre2                                 10.40  h1c4e4bc_0            conda-forge/osx-64     Cached
  + perl                               5.22.0.1  0                     conda-forge/osx-64       15MB
  + phantompeakqualtools                  1.2.2  hdfd78af_1            bioconda/noarch        Cached
  + pip                                    24.0  pyhd8ed1ab_0          conda-forge/noarch     Cached
  + pixman                               0.43.4  h73e2aa4_0            conda-forge/osx-64     Cached
  + python                               3.12.2  h9f0c242_0_cpython    conda-forge/osx-64     Cached
  + python_abi                             3.12  4_cp312               conda-forge/osx-64     Cached
  + pyyaml                                6.0.1  py312h104f124_1       conda-forge/osx-64     Cached
  + r-base                                4.3.1  h61172b1_5            conda-forge/osx-64     Cached
  + r-bh                               1.84.0_0  r43hc72bb7e_0         conda-forge/noarch     Cached
  + r-bitops                              1.0_7  r43h6dc245f_2         conda-forge/osx-64     Cached
  + r-catools                            1.18.2  r43hac7d2d5_2         conda-forge/osx-64     Cached
  + r-codetools                          0.2_20  r43hc72bb7e_0         conda-forge/noarch     Cached
  + r-cpp11                               0.4.7  r43hc72bb7e_0         conda-forge/noarch     Cached
  + r-crayon                              1.5.2  r43hc72bb7e_2         conda-forge/noarch     Cached
  + r-formatr                              1.14  r43hc72bb7e_1         conda-forge/noarch     Cached
  + r-futile.logger                       1.4.3  r43hc72bb7e_1005      conda-forge/noarch     Cached
  + r-futile.options                      1.0.1  r43hc72bb7e_1004      conda-forge/noarch     Cached
  + r-lambda.r                            1.2.4  r43hc72bb7e_3         conda-forge/noarch     Cached
  + r-rcpp                               1.0.12  r43h29979af_0         conda-forge/osx-64     Cached
  + r-rcurl                           1.98_1.14  r43hbd64cb6_0         conda-forge/osx-64     Cached
  + r-snow                                0.4_4  r43hc72bb7e_2         conda-forge/noarch     Cached
  + r-snowfall                         1.84_6.3  r43hc72bb7e_0         conda-forge/noarch     Cached
  + r-spp                                1.16.0  r43h71b3455_9         bioconda/osx-64        Cached
  + readline                                8.2  h9e318b2_1            conda-forge/osx-64     Cached
  + samtools                             1.19.2  hd510865_1            bioconda/osx-64        Cached
  + setuptools                           69.2.0  pyhd8ed1ab_0          conda-forge/noarch     Cached
  + sigtool                               0.1.3  h88f4db0_0            conda-forge/osx-64     Cached
  + tapi                              1100.0.11  h9ce4665_0            conda-forge/osx-64     Cached
  + tk                                   8.6.13  h1abcd95_1            conda-forge/osx-64     Cached
  + tktable                                2.10  ha166976_5            conda-forge/osx-64     Cached
  + toml                                 0.10.2  pyhd8ed1ab_0          conda-forge/noarch     Cached
  + tomlkit                              0.12.4  pyha770c72_0          conda-forge/noarch     Cached
  + tzdata                                2024a  h0c530f3_0            conda-forge/noarch     Cached
  + wheel                                0.43.0  pyhd8ed1ab_1          conda-forge/noarch     Cached
  + xmltodict                            0.13.0  pyhd8ed1ab_0          conda-forge/noarch     Cached
  + xz                                    5.2.6  h775f41a_0            conda-forge/osx-64     Cached
  + yaml                                  0.2.5  h0d85af4_2            conda-forge/osx-64     Cached
  + yq                                    3.2.3  pyhd8ed1ab_0          conda-forge/noarch     Cached
  + zlib                                 1.2.13  h8a1eda9_5            conda-forge/osx-64     Cached
  + zstd                                  1.5.5  h829000d_0            conda-forge/osx-64     Cached

  Summary:

  Install: 140 packages

  Total download: 17MB

──────────────────────────────────────────────────────────────────────────────────────────────────────

parallel                                             1.2MB @ 901.0kB/s  1.4s
perl                                                15.3MB @   8.7MB/s  1.8s
Preparing transaction: done
Verifying transaction: done
Executing transaction: done

To activate this environment, use

     $ mamba activate ppqt_env

To deactivate an active environment, use

     $ mamba deactivate
```
</details>
<br />

<a id="b-run-phantompeakqualtools"></a>
## b. Run `phantompeakqualtools`
<a id="code-9"></a>
### Code
<details>
<summary><i>Code: X.b. Run `phantompeakqualtools`</i></summary>

```bash
#!/bin/bash

#  Define function ============================================================
#  Function to return an error message and exit code 1, which stops the
#+ interactive execution of code
function error_and_return() {
    echo "Error: ${1}" >&2
    return 1
}


#  Function to deactivate a Conda/Mamba environment
function deactivate_env() {
    if [[ "${CONDA_DEFAULT_ENV}" != "base" ]]; then
        if ! mamba deactivate &> /dev/null; then
            if ! conda deactivate &> /dev/null; then
                if ! source deactivate &> /dev/null; then
                    error_and_return "Failed to deactivate environment."
                    return 1
                fi
            fi
        fi
    fi

    return 0
}


#  Function to check if a specific Conda/Mamba environment is installed
function check_env_installed() {
    local env_name="${1}"

    if conda env list | grep -q "^${env_name} "; then
        return 0
    else
        echo "Environment \"${env_name}\" is not installed."
        return 1
    fi
}


#  Function to activate a specific Conda/Mamba environment
function activate_env() {
    local env_name="${1}"

    deactivate_env

    if ! mamba activate "${env_name}" &> /dev/null; then
        if ! conda activate "${env_name}" &> /dev/null; then
            if ! source activate "${env_name}" &> /dev/null; then
                error_and_return "Failed to activate environment \"${env_name}\"."
                return 1
            fi
        fi
    fi
    
    echo "Environment \"${env_name}\" activated successfully."
    return 0
}


#  Initialize variables and arrays ============================================
dir_base="${HOME}/tsukiyamalab"                          # Base directory for lab data
dir_repo="Kris/2023_tutorial_ChIP-seq"                   # Repository directory
dir_bams="03_bam/bowtie2/bam"                            # Directory for BAMs
dir_ppqt="03_bam/bowtie2/ppqt"                           # Directory for MACS3 outfiles

fdr=0.05
# gsize=12157105
# keep_dup="auto"
#
# time="4:00:00"                                           # Job time for SLURM (H:MM:SS)
# threads=1                                                # Number of threads for SLURM jobs

#  Initialize an indexed array of BAM file stems
unset file_bam_stems
typeset -a file_bam_stems=(
    "Q_untagged_5781"
    "Q_Esa5_7041"
    "Q_Esa5_7691"
    "Q_Rpd3_7568"
    "Q_Rpd3_7569"
    "Q_Gcn5_7692"
    "Q_Gcn5_7709"
    "G1_Hho1_6336"
    "G1_Hho1_6337"
    "G1_Hmo1_7750"
    "G1_Hmo1_7751"
    "G2M_Hho1_6336"
    "G2M_Hho1_6337"
    "G2M_Hmo1_7750"
    "G2M_Hmo1_7751"
    "Q_Hho1_6336"
    "Q_Hho1_6337"
    "Q_Hmo1_7750"
    "Q_Hmo1_7751"
)


#  Do the main work ===========================================================
#  Set flags for checking variable and array assignments
check_variables=true
check_array=true

#  If check_variables is true, then echo the variable assignments
if ${check_variables}; then
    echo "
    dir_base=${dir_base}
    dir_repo=${dir_repo}
    dir_bwt2=${dir_bams}
    dir_ppqt=${dir_ppqt}
    
    fdr=${fdr}
    # gsize=${gsize}
    # keep_dup=${keep_dup}
    #
    # time=${time}
    # threads=${threads}
    "
fi

#  Echo array contents if check_array is true
if ${check_array}; then
    for i in "${!file_bam_stems[@]}"; do
        file="${file_bam_stems[${i}]}"

        # echo "${file}"
        ls -lhaFG "${dir_base}/${dir_repo}/${dir_bams}/IP_${file}.sort-coord.bam"
        ls -lhaFG "${dir_base}/${dir_repo}/${dir_bams}/in_${file}.sort-coord.bam"
    done
fi

#  Initialize conda/mamba environment containing necessary programs for
#+ alignment, quality checks, and post-processing
env_name="ppqt_env"

check_env_installed "${env_name}"
activate_env "${env_name}"

#  Navigate to the work directory
cd "${dir_base}/${dir_repo}" \
    || error_and_return "Failed to cd to ${dir_base}/${dir_repo}."

#  If it doesn't exist, create a directory to store MACS3 outfiles
if [[ ! -d "${dir_ppqt}" ]]; then
    mkdir -p "${dir_ppqt}"
fi

#  Set flags: checking variables, checking and submitting Bowtie2 jobs
print_iteration=true
check_variables=true
check_operation=true
run_operation=true

for i in "${!file_bam_stems[@]}"; do
    # i=18
    index="${i}"
    iter=$(( index + 1 ))
    stem="${file_bam_stems[${index}]}"
    job_name="$(echo ${dir_ppqt} | sed 's:\/:_:g').${stem}"
    
    in="${dir_bams}/in_${stem}.sort-coord.bam"
    IP="${dir_bams}/IP_${stem}.sort-coord.bam"
    out="${dir_ppqt}/IP_${stem}"

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
        stem=${stem}
        job_name=${job_name}

        in=${in}
        IP=${IP}
        out=${out}
        "
    fi

    # echo "
    # run_spp.R \\
    #         -c=\"${IP}\" \\
    #         -savp \\
    #         -out=\"${out}\"
    # "
    #
    # run_spp.R \
    #     -c="${IP}" \
    #     -savp="${out}.pdf" \
    #     -out="${out}.txt"

    echo "
    run_spp.R \\
        -c=\"${IP}\" \\
        -i=\"${in}\" \\
        -fdr=\"${fdr}\" \\
        -odir=\"${dir_ppqt}\" \\
        -savn=\"${out}.narrowPeak\" \\
        -savr=\"${out}.regionPeak\" \\
        -savp=\"${out}.CCP.pdf\" \\
        -savd=\"${out}.Rdata\" \\
        -rf
    "

    run_spp.R \
        -c="${IP}" \
        -i="${in}" \
        -fdr="${fdr}" \
        -filtchr="SP|Mito" \
        -odir="${dir_ppqt}" \
        -savn="${out}.narrowPeak" \
        -savr="${out}.regionPeak" \
        -savp="${out}.CCP.pdf" \
        -savd="${out}.Rdata" \
        -rf

    # run_spp.R \
    #     -c=<ChIP_tagalign/BAM_file> \
    #     -i=<control_tagalign/BAM_file> \
    #     -npeak=<npeaks> \
    #     -odir=<peak_call_output_dir> \
    #     -savr \
    #     -savp \
    #     -savd \
    #     -rf


    # rm "${out}.txt"
    #     -i="${in}"

    if ${check_operation}; then
        echo "
sbatch << EOF
#!/bin/bash

#SBATCH --job-name=\"${job_name}\"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=${threads}
#SBATCH --time=${time}
#SBATCH --error=\"$(dirname ${dir_bams})/err_out/${job_name}.%A.stderr.txt\"
#SBATCH --output=\"$(dirname ${dir_bams})/err_out/${job_name}.%A.stdout.txt\"

macs3 callpeak \\
    --name \"${stem}\" \\
    --treatment \"${IP}\" \\
    --control \"${in}\" \\
    --format \"BAMPE\" \\
    --gsize \"${gsize}\" \\
    --keep-dup \"${keep_dup}\" \\
    --outdir \"${dir_macs}\" \\
    --bdg \\
    --SPMR \\
    --verbose 3

if [[ -f \"${dir_macs}/${stem}_summits.bed\" ]]; then
    find \"${dir_macs}\" \\
        -type f \\
        \( \\
               -name \"${stem}\"*\".bdg\" \\
            -o -name \"${stem}\"*\".narrowPeak\" \\
            -o -name \"${stem}\"*\".xls\" \\
            -o -name \"${stem}\"*\".bed\" \\
        \) \\
        -exec gzip {} \;
fi
EOF
        "
    fi

    if ${run_operation}; then
sbatch << EOF
#!/bin/bash

#SBATCH --job-name="${job_name}"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=${threads}
#SBATCH --time=${time}
#SBATCH --error="$(dirname ${dir_bams})/err_out/${job_name}.%A.stderr.txt"
#SBATCH --output="$(dirname ${dir_bams})/err_out/${job_name}.%A.stdout.txt"

macs3 callpeak \
    --name "${stem}" \
    --treatment "${IP}" \
    --control "${in}" \
    --format "BAMPE" \
    --gsize "${gsize}" \
    --keep-dup "${keep_dup}" \
    --outdir "${dir_macs}" \
    --bdg \
    --SPMR \
    --verbose 3

if [[ -f "${dir_macs}/${stem}_summits.bed" ]]; then
    find "${dir_macs}" \
        -type f \
        \( \
               -name "${stem}"*".bdg" \
            -o -name "${stem}"*".narrowPeak" \
            -o -name "${stem}"*".xls" \
            -o -name "${stem}"*".bed" \
        \) \
        -exec gzip {} \;
fi
EOF
    fi

    sleep 0.2
done
```
</details>
<br />

<a id="4-call-peaks-with-macs3"></a>
# 4. Call peaks with MACS3
<a id="a-install-macs3"></a>
## a. Install MACS3
<a id="code-10"></a>
### Code
<details>
<summary><i>Code: 4.a. Install MACS3</i></summary>

```bash
#!/bin/bash

#  Define functions ===========================================================
#  Function to return an error message and exit code 1, which stops the
#+ interactive execution of code
function error_and_return() {
    echo "Error: ${1}" >&2
    return 1
}


#  Function to check if Mamba is installed
function check_mamba_installed() {
    if ! : mamba &> /dev/null; then
        echo "Mamba is not installed on your system. Mamba is a package manager" 
        echo "that makes package installations faster and more reliable."
        echo ""
        echo "For installation instructions, please check the following link:"
        echo "https://github.com/mamba-org/mamba#installation"
        return 1
    fi
    
    return 0
}


#  Function to deactivate a Conda/Mamba environment
function deactivate_env() {
    if [[ "${CONDA_DEFAULT_ENV}" != "base" ]]; then
        if ! mamba deactivate &> /dev/null; then
            if ! conda deactivate &> /dev/null; then
                if ! source deactivate &> /dev/null; then
                    error_and_return "Failed to deactivate environment."
                    return 1
                fi
            fi
        fi
    fi

    return 0
}


#  Function to check if a specific Conda/Mamba environment is installed
function check_env_installed() {
    local env_name="${1}"

    if conda env list | grep -q "^${env_name} "; then
        return 0
    else
        echo "Environment \"${env_name}\" is not installed."
        return 1
    fi
}


#  Function to activate a specific Conda/Mamba environment
function activate_env() {
    local env_name="${1}"

    if ! mamba activate "${env_name}" &> /dev/null; then
        if ! conda activate "${env_name}" &> /dev/null; then
            if ! source activate "${env_name}" &> /dev/null; then
                error_and_return "Failed to activate environment \"${env_name}\"."
                return 1
            fi
        fi
    fi
    
    echo "Environment \"${env_name}\" activated successfully."
    return 0
}


#  Initialize variables and arrays ============================================
env_name="macs3_env"


#  Do the main work ===========================================================
#  Set flag(s)
create_mamba_env=true  # Install mamba environment if not detected

#  Check that Mamba is installed and in PATH
check_mamba_installed

#  If not in base environment, then deactivate current environment
deactivate_env

#  Check that environment assigned to env_name is installed; if environment
#+ assigned to env_name is not installed, run the following
if check_env_installed "${env_name}"; then
    #  Handle the case when the environment is already installed
    echo "Activating environment ${env_name}"
    
    activate_env "${env_name}"
else
    #  Handle the case when the environment is not installed
    echo "Creating environment ${env_name}"
    
    if ${create_mamba_env}; then
        #  Switch `--yes` is set, which means no user input is required
        #NOTE Running this on FHCC Rhino; ergo, no CONDA_SUBDIR=osx-64
        mamba create \
            --yes \
            -n "${env_name}" \
            -c conda-forge \
                python=3.10 \
                pip
        
        source activate "${env_name}"

        pip install macs3

        deactivate_env
    fi
fi
```
</details>
<br />

<a id="b-run-macs3"></a>
## b. Run MACS3
<a id="code-11"></a>
### Code
<details>
<summary><i>Code: Run MACS3</i></summary>

```bash
#!/bin/bash

#  Define function ============================================================
#  Function to return an error message and exit code 1, which stops the
#+ interactive execution of code
function error_and_return() {
    echo "Error: ${1}" >&2
    return 1
}


#  Function to deactivate a Conda/Mamba environment
function deactivate_env() {
    if [[ "${CONDA_DEFAULT_ENV}" != "base" ]]; then
        if ! mamba deactivate &> /dev/null; then
            if ! conda deactivate &> /dev/null; then
                if ! source deactivate &> /dev/null; then
                    error_and_return "Failed to deactivate environment."
                    return 1
                fi
            fi
        fi
    fi

    return 0
}


#  Function to check if a specific Conda/Mamba environment is installed
function check_env_installed() {
    local env_name="${1}"

    if conda env list | grep -q "^${env_name} "; then
        return 0
    else
        echo "Environment \"${env_name}\" is not installed."
        return 1
    fi
}


#  Function to activate a specific Conda/Mamba environment
function activate_env() {
    local env_name="${1}"

    deactivate_env

    if ! mamba activate "${env_name}" &> /dev/null; then
        if ! conda activate "${env_name}" &> /dev/null; then
            if ! source activate "${env_name}" &> /dev/null; then
                error_and_return "Failed to activate environment \"${env_name}\"."
                return 1
            fi
        fi
    fi
    
    echo "Environment \"${env_name}\" activated successfully."
    return 0
}


#  Initialize variables and arrays ============================================
dir_base="${HOME}/tsukiyamalab"                          # Base directory for lab data
dir_repo="Kris/2023_rDNA"                                # Repository directory
dir_work="results/2023-0406_tutorial_ChIP-seq_analyses"  # Work directory
dir_bams="03_bam/bowtie2/bam"                            # Directory for BAMs
dir_macs="03_bam/bowtie2/macs3"                          # Directory for MACS3 outfiles

gsize=12157105
keep_dup="auto"

time="4:00:00"                                           # Job time for SLURM (H:MM:SS)
threads=1                                                # Number of threads for SLURM jobs

#  Initialize an indexed array of BAM file stems
unset file_bam_stems && typeset -a file_bam_stems=(
    "Q_untagged_5781"
    "Q_Esa5_7041"
    "Q_Esa5_7691"
    "Q_Rpd3_7568"
    "Q_Rpd3_7569"
    "Q_Gcn5_7692"
    "Q_Gcn5_7709"
)


#  Do the main work ===========================================================
#  Set flags for checking variable and array assignments
check_variables=true
check_array=true

#  If check_variables is true, then echo the variable assignments
if ${check_variables}; then
    echo "
    dir_base=${dir_base}
    dir_repo=${dir_repo}
    dir_work=${dir_work}
    dir_bwt2=${dir_bams}
    dir_macs=${dir_macs}
    
    gsize=${gsize}
    keep_dup=${keep_dup}
    
    time=${time}
    threads=${threads}
    "
fi

#  Echo array contents if check_array is true
if ${check_array}; then
    for i in "${!file_bam_stems[@]}"; do
        file="${file_bam_stems[${i}]}"

        # echo "${file}"
        ls -lhaFG "${dir_base}/${dir_repo}/${dir_work}/${dir_bams}/IP_${file}.sort-coord.bam"
        ls -lhaFG "${dir_base}/${dir_repo}/${dir_work}/${dir_bams}/in_${file}.sort-coord.bam"
    done
fi

#  Initialize conda/mamba environment containing necessary programs for
#+ alignment, quality checks, and post-processing
env_name="macs3_env"

check_env_installed "${env_name}"
activate_env "${env_name}"

#  Navigate to the work directory
cd "${dir_base}/${dir_repo}/${dir_work}" \
    || error_and_return "Failed to cd to ${dir_base}/${dir_repo}/${dir_work}."

#  If it doesn't exist, create a directory to store MACS3 outfiles
if [[ ! -d "${dir_macs}" ]]; then
    mkdir -p "${dir_macs}"
fi

#  Set flags: checking variables, checking and submitting Bowtie2 jobs
print_iteration=true
check_variables=true
check_operation=true
run_operation=true

for i in "${!file_bam_stems[@]}"; do
    # i=0
    index="${i}"
    iter=$(( index + 1 ))
    stem="${file_bam_stems[${index}]}"
    job_name="$(echo ${dir_macs} | sed 's:\/:_:g').${stem}"
    
    in="${dir_bams}/in_${stem}.sort-coord.bam"
    IP="${dir_bams}/IP_${stem}.sort-coord.bam"

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
        stem=${stem}
        job_name=${job_name}

        in=${in}
        IP=${IP}
        "
    fi

    if ${check_operation}; then
        echo "
sbatch << EOF
#!/bin/bash

#SBATCH --job-name=\"${job_name}\"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=${threads}
#SBATCH --time=${time}
#SBATCH --error=\"$(dirname ${dir_bams})/err_out/${job_name}.%A.stderr.txt\"
#SBATCH --output=\"$(dirname ${dir_bams})/err_out/${job_name}.%A.stdout.txt\"

macs3 callpeak \\
    --name \"${stem}\" \\
    --treatment \"${IP}\" \\
    --control \"${in}\" \\
    --format \"BAMPE\" \\
    --gsize \"${gsize}\" \\
    --keep-dup \"${keep_dup}\" \\
    --outdir \"${dir_macs}\" \\
    --bdg \\
    --SPMR \\
    --verbose 3

if [[ -f \"${dir_macs}/${stem}_summits.bed\" ]]; then
    find \"${dir_macs}\" \\
        -type f \\
        \( \\
               -name \"${stem}\"*\".bdg\" \\
            -o -name \"${stem}\"*\".narrowPeak\" \\
            -o -name \"${stem}\"*\".xls\" \\
            -o -name \"${stem}\"*\".bed\" \\
        \) \\
        -exec gzip {} \;
fi
EOF
        "
    fi

    if ${run_operation}; then
sbatch << EOF
#!/bin/bash

#SBATCH --job-name="${job_name}"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=${threads}
#SBATCH --time=${time}
#SBATCH --error="$(dirname ${dir_bams})/err_out/${job_name}.%A.stderr.txt"
#SBATCH --output="$(dirname ${dir_bams})/err_out/${job_name}.%A.stdout.txt"

macs3 callpeak \
    --name "${stem}" \
    --treatment "${IP}" \
    --control "${in}" \
    --format "BAMPE" \
    --gsize "${gsize}" \
    --keep-dup "${keep_dup}" \
    --outdir "${dir_macs}" \
    --bdg \
    --SPMR \
    --verbose 3

if [[ -f "${dir_macs}/${stem}_summits.bed" ]]; then
    find "${dir_macs}" \
        -type f \
        \( \
               -name "${stem}"*".bdg" \
            -o -name "${stem}"*".narrowPeak" \
            -o -name "${stem}"*".xls" \
            -o -name "${stem}"*".bed" \
        \) \
        -exec gzip {} \;
fi
EOF
    fi

    sleep 0.2
done
```
</details>
<br />

<a id="c-run-macs3-with-pooled-replicates"></a>
## c. Run MACS3 with pooled replicates
<a id="code-12"></a>
### Code
<details>
<summary><i>Code: Run MACS3 with pooled replicates</i></summary>

```bash
#!/bin/bash

#  Define functions ===========================================================
#  Function to return an error message and exit code 1, which stops the
#+ interactive execution of code
function error_and_return() {
    echo "Error: ${1}" >&2
    return 1
}


#  Function to deactivate a Conda/Mamba environment
function deactivate_env() {
    if [[ "${CONDA_DEFAULT_ENV}" != "base" ]]; then
        if ! mamba deactivate &> /dev/null; then
            if ! conda deactivate &> /dev/null; then
                if ! source deactivate &> /dev/null; then
                    error_and_return "Failed to deactivate environment."
                    return 1
                fi
            fi
        fi
    fi

    return 0
}


#  Function to check if a specific Conda/Mamba environment is installed
function check_env_installed() {
    local env_name="${1}"

    if conda env list | grep -q "^${env_name} "; then
        return 0
    else
        echo "Environment \"${env_name}\" is not installed."
        return 1
    fi
}


#  Function to activate a specific Conda/Mamba environment
function activate_env() {
    local env_name="${1}"

    deactivate_env

    if ! mamba activate "${env_name}" &> /dev/null; then
        if ! conda activate "${env_name}" &> /dev/null; then
            if ! source activate "${env_name}" &> /dev/null; then
                error_and_return "Failed to activate environment \"${env_name}\"."
                return 1
            fi
        fi
    fi
    
    echo "Environment \"${env_name}\" activated successfully."
    return 0
}


#  Initialize variables and arrays ============================================
dir_base="${HOME}/tsukiyamalab"                          # Base directory for lab data
dir_repo="Kris/2023_rDNA"                                # Repository directory
dir_work="results/2023-0406_tutorial_ChIP-seq_analyses"  # Work directory
dir_bams="03_bam/bowtie2/bam"                            # Directory for BAMs
dir_macs="03_bam/bowtie2/macs3"                          # Directory for MACS3 outfiles

gsize=12157105
keep_dup="auto"

time="4:00:00"                                           # Job time for SLURM (H:MM:SS)
threads=1                                                # Number of threads for SLURM jobs

#  Initialize an indexed array of BAM file stems
unset file_bam_stems && typeset -A file_bam_stems=(
    ["Q_Esa5_7041"]="Q_Esa5_7691"
    ["Q_Rpd3_7568"]="Q_Rpd3_7569"
    ["Q_Gcn5_7692"]="Q_Gcn5_7709"
)


#  Do the main work ===========================================================
#  Set flags for checking variable and array assignments
check_variables=true
check_array=true

#  If check_variables is true, then echo the variable assignments
if ${check_variables}; then
    echo "
    dir_base=${dir_base}
    dir_repo=${dir_repo}
    dir_work=${dir_work}
    dir_bwt2=${dir_bams}
    dir_macs=${dir_macs}

    gsize=${gsize}
    keep_dup=${keep_dup}
    
    time=${time}
    threads=${threads}
    "
fi

#  Echo array contents if check_array is true
if ${check_array}; then
    for i in "${!file_bam_stems[@]}"; do
          key="${i}"
        value="${file_bam_stems[${key}]}"

        echo "[\"IP_${key}.sort-coord.bam\"]=\"IP_${value}.sort-coord.bam\""
        echo "[\"in_${key}.sort-coord.bam\"]=\"in_${value}.sort-coord.bam\""
        echo ""
        
        ls -lhaFG "${dir_base}/${dir_repo}/${dir_work}/${dir_bams}/IP_${key}.sort-coord.bam"
        ls -lhaFG "${dir_base}/${dir_repo}/${dir_work}/${dir_bams}/in_${key}.sort-coord.bam"
        ls -lhaFG "${dir_base}/${dir_repo}/${dir_work}/${dir_bams}/IP_${value}.sort-coord.bam"
        ls -lhaFG "${dir_base}/${dir_repo}/${dir_work}/${dir_bams}/in_${value}.sort-coord.bam"
        echo ""
        echo ""
    done
fi

#  Initialize conda/mamba environment containing necessary programs for
#+ alignment, quality checks, and post-processing
env_name="macs3_env"

check_env_installed "${env_name}"
activate_env "${env_name}"

#  Navigate to the work directory
cd "${dir_base}/${dir_repo}/${dir_work}" \
    || error_and_return "Failed to cd to ${dir_base}/${dir_repo}/${dir_work}."

#  If it doesn't exist, create a directory to store MACS3 outfiles
if [[ ! -d "${dir_macs}" ]]; then
    mkdir -p "${dir_macs}"
fi

#  Set flags: checking variables, checking and submitting Bowtie2 jobs
print_iteration=true
check_variables=true
check_operation=true
run_operation=false

for i in "${!file_bam_stems[@]}"; do
    # i="Q_Esa5_7041"
    iter=$(( index + 1 ))
    key="${i}"
    value="${file_bam_stems[${key}]}"
    stem="${key}.${value}"
    job_name="$(echo ${dir_macs} | sed 's:\/:_:g').${stem}"
    
    in_1="${dir_bams}/in_${key}.sort-coord.bam"
    in_2="${dir_bams}/in_${value}.sort-coord.bam"

    IP_1="${dir_bams}/IP_${key}.sort-coord.bam"
    IP_2="${dir_bams}/IP_${value}.sort-coord.bam"

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
        iter=${iter}
        key=${key}
        value=${value}
        stem=${stem}
        job_name=${job_name}

        in_1=${in_1}
        in_2=${in_2}
        IP_1=${IP_1}
        IP_2=${IP_2}
        "
    fi

    if ${check_operation}; then
        echo "
sbatch << EOF
#!/bin/bash

#SBATCH --job-name=\"${job_name}\"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=${threads}
#SBATCH --time=${time}
#SBATCH --error=\"$(dirname ${dir_bams})/err_out/${job_name}.%A.stderr.txt\"
#SBATCH --output=\"$(dirname ${dir_bams})/err_out/${job_name}.%A.stdout.txt\"

macs3 callpeak \\
    --name \"${stem}\" \\
    --treatment \"${IP_1}\" \"${IP_2}\" \\
    --control \"${in_1}\" \"${in_2}\" \\
    --format \"BAMPE\" \\
    --gsize \"${gsize}\" \\
    --keep-dup \"${keep_dup}\" \\
    --outdir \"${dir_macs}\" \\
    --bdg \\
    --SPMR \\
    --verbose 3

if [[ -f \"${dir_macs}/${stem}_summits.bed\" ]]; then
    find \"${dir_macs}\" \\
        -type f \\
        \( \\
               -name \"${stem}\"*\".bdg\" \\
            -o -name \"${stem}\"*\".narrowPeak\" \\
            -o -name \"${stem}\"*\".xls\" \\
            -o -name \"${stem}\"*\".bed\" \\
        \) \\
        -exec gzip {} \;
fi
EOF
        "
    fi

    if ${run_operation}; then
sbatch << EOF
#!/bin/bash

#SBATCH --job-name="${job_name}"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=${threads}
#SBATCH --time=${time}
#SBATCH --error="$(dirname ${dir_bams})/err_out/${job_name}.%A.stderr.txt"
#SBATCH --output="$(dirname ${dir_bams})/err_out/${job_name}.%A.stdout.txt"

macs3 callpeak \
    --name "${stem}" \
    --treatment "${IP_1}" "${IP_2}" \
    --control "${in_1}" "${in_2}" \
    --format "BAMPE" \
    --gsize "${gsize}" \
    --keep-dup "${keep_dup}" \
    --outdir "${dir_macs}" \
    --bdg \
    --SPMR \
    --verbose 3

if [[ -f "${dir_macs}/${stem}_summits.bed" ]]; then
    find "${dir_macs}" \
        -type f \
        \( \
               -name "${stem}"*".bdg" \
            -o -name "${stem}"*".narrowPeak" \
            -o -name "${stem}"*".xls" \
            -o -name "${stem}"*".bed" \
        \) \
        -exec gzip {} \;
fi
EOF
    fi

    sleep 0.2
done
```
</details>
<br />
<br />

<a id="5-subset-peaks"></a>
# 5. Subset peaks
<a id="a-install-environment-for-interactive-r-scripting"></a>
## a. Install environment for interactive R scripting
The purpose of the following two code chunks is to establish a computational environment, `R_env`, needed for interactive `R` scripting. This environment encompasses various programs and their dependencies necessary for our subsequent work with `narrowPeak` files generated by MACS3. The creation and setup of the `R_env` environment is expected to take between 10 to 20 minutes. After installation is complete, ensure the `R` interpreter is activated before proceeding to the second code chunk, in which we install the non-CRAN package `colorout`, which colorizes output in Unix-like terminal emulators, enhancing the readability of `R` operation outputs. Operations performed in this chunk:
- Mamba installation check
- Environment management
- Package installation
- Environment activation and custom package installation

<a id="bash-code"></a>
### Bash code
<details>
<summary><i>Bash code: 5.a. Install environment for interactive R scripting</i></summary>

```bash
#!/bin/bash

#  Define functions ===========================================================
#  Function to return an error message and exit code 1, which stops the
#+ interactive execution of code
function error_and_return() {
    echo "Error: ${1}" >&2
    return 1
}


#  Function to check if Mamba is installed
function check_mamba_installed() {
    if ! : mamba &> /dev/null; then
        echo "Mamba is not installed on your system. Mamba is a package manager" 
        echo "that makes package installations faster and more reliable."
        echo ""
        echo "For installation instructions, please check the following link:"
        echo "https://github.com/mamba-org/mamba#installation"
        return 1
    fi
    
    return 0
}


#  Function to deactivate a Conda/Mamba environment
function deactivate_env() {
    if [[ "${CONDA_DEFAULT_ENV}" != "base" ]]; then
        if ! mamba deactivate &> /dev/null; then
            if ! conda deactivate &> /dev/null; then
                if ! source deactivate &> /dev/null; then
                    error_and_return "Failed to deactivate environment."
                    return 1
                fi
            fi
        fi
    fi

    return 0
}


#  Function to check if a specific Conda/Mamba environment is installed
function check_env_installed() {
    local env_name="${1}"

    if conda env list | grep -q "^${env_name} "; then
        return 0
    else
        echo "Environment \"${env_name}\" is not installed."
        return 1
    fi
}


#  Function to activate a specific Conda/Mamba environment
function activate_env() {
    local env_name="${1}"

    if ! mamba activate "${env_name}" &> /dev/null; then
        if ! conda activate "${env_name}" &> /dev/null; then
            if ! source activate "${env_name}" &> /dev/null; then
                error_and_return "Failed to activate environment \"${env_name}\"."
                return 1
            fi
        fi
    fi
    
    echo "Environment \"${env_name}\" activated successfully."
    return 0
}


#  Initialize variables and arrays ============================================
env_name="R_env"


#  Do the main work ===========================================================
#  Set flag(s)
create_mamba_env=true  # Install mamba environment if not detected

#  Check that Mamba is installed and in PATH
check_mamba_installed

#  If not in base environment, then deactivate current environment
deactivate_env

#  Check that environment assigned to env_name is installed; if environment
#+ assigned to env_name is not installed, run the following
if check_env_installed "${env_name}"; then
    #  Handle the case when the environment is already installed
    echo "Activating environment ${env_name}"
    
    activate_env "${env_name}"
else
    #  Handle the case when the environment is not installed
    echo "Creating environment ${env_name}"
    
    #  In my experience, it takes 10-20 minutes to complete the following
    #+ installation
    if ${create_mamba_env}; then
        #  Switch `--yes` is set, which means no user input is required
        #NOTE Running this on FHCC Rhino; ergo, no CONDA_SUBDIR=osx-64
        mamba create \
            --yes \
            -n "${env_name}" \
            -c bioconda \
            -c conda-forge \
                bioconductor-annotationdbi \
                bioconductor-chipqc \
                bioconductor-chipseeker \
                bioconductor-clusterprofiler \
                bioconductor-deseq2 \
                bioconductor-diffbind \
                bioconductor-edger \
                bioconductor-enhancedvolcano \
                bioconductor-genomicfeatures \
                bioconductor-genomicranges \
                bioconductor-ihw \
                bioconductor-iranges \
                bioconductor-pcatools \
                bioconductor-sva \
                deeptools \
                phantompeakqualtools \
                r-argparse \
                r-dendextend \
                r-devtools \
                r-ggalt \
                r-ggpubr \
                r-ggrepel \
                r-pheatmap \
                r-readxl \
                r-rjson \
                r-tidyverse \
                r-upsetr \
                r-venneuler \
                r-writexl \
                r-xml2 \
                rename

        #  Activate the new environment and proceed to install a package
        #+ unavailable via CRAN and thus unavailable via bioconda, conda-forge,
        #+ etc.: colorout, which colorizes R output when running in a *nix
        #+ (e.g., Linux and OS X) terminal emulator
        source activate "${env_name}"

        #  To install colorout, first invoke the R interpreter 
        R
    fi
fi
```
</details>
<br />

<a id="r-code"></a>
### R code
<details>
<summary><i>R code: 5.a. Install environment for interactive R scripting</i></summary>

```r
#!/usr/bin/env Rscript

#  Check if the R package "colorout" is installed; if not, then install it from
#+ GitHub using the devtools package
if (!requireNamespace("colorout", quietly = TRUE)) {
    devtools::install_github("jalvesaq/colorout")
}

q()  # No need to save your workspace
```
</details>
<br />


<a id="b-perform-set-operations-with-peak-intervals"></a>
## b. Perform set operations with peak intervals
<a id="i-get-situated"></a>
### i. Get situated
<a id="bash-code-1"></a>
#### Bash code
The purpose of this `Bash` code chunk is to initialize an environment for interactive `R` scripting, `R_env`, that in turn allows us to access programs needed to process and perform set operations on `narrowPeak` files output by MACS3. Please modify the script to navigate (`cd`) into the directory containing the `narrowPeak` files, and ensure the `R` interpreter is activated before proceeding to the subsequent `R` code chunks. Operations performed in this chunk:
- Environment initialization
- Error handling
- Environment management
- Directory navigation
- Variable assignment checking
- Launching `R`

<details>
<summary><i>Bash code: 5.b.i. Get situated</i></summary>

```bash
#!/bin/bash

#  Define functions ===========================================================
#  Function to return an error message and exit code 1, which stops the
#+ interactive execution of code
function error_and_return() {
    echo "Error: ${1}" >&2
    return 1
}


#  Function to deactivate a Conda/Mamba environment
function deactivate_env() {
    if [[ "${CONDA_DEFAULT_ENV}" != "base" ]]; then
        if ! mamba deactivate &> /dev/null; then
            if ! conda deactivate &> /dev/null; then
                if ! source deactivate &> /dev/null; then
                    error_and_return "Failed to deactivate environment."
                    return 1
                fi
            fi
        fi
    fi

    return 0
}


#  Function to check if a specific Conda/Mamba environment is installed
function check_env_installed() {
    local env_name="${1}"

    if conda env list | grep -q "^${env_name} "; then
        return 0
    else
        echo "Environment \"${env_name}\" is not installed."
        return 1
    fi
}


#  Function to activate a specific Conda/Mamba environment
function activate_env() {
    local env_name="${1}"

    deactivate_env

    if ! mamba activate "${env_name}" &> /dev/null; then
        if ! conda activate "${env_name}" &> /dev/null; then
            if ! source activate "${env_name}" &> /dev/null; then
                error_and_return "Failed to activate environment \"${env_name}\"."
                return 1
            fi
        fi
    fi
    
    echo "Environment \"${env_name}\" activated successfully."
    return 0
}


#  Initialize variables and arrays ============================================
dir_base="${HOME}/tsukiyamalab"                          # Base directory for lab data
dir_repo="Kris/2023_rDNA"                                # Repository directory
dir_work="results/2023-0406_tutorial_ChIP-seq_analyses"  # Work directory
dir_macs="03_bam/bowtie2/macs3"                          # Directory for MACS3 outfiles

env_name="R_env"


#  Do the main work ===========================================================
#  Set flags for checking variable and array assignments
check_variables=true

#  If check_variables is true, then echo the variable assignments
if ${check_variables}; then
    echo "
    dir_base=${dir_base}
    dir_repo=${dir_repo}
    dir_work=${dir_work}
    dir_macs=${dir_macs}

    env_name=${env_name}
    "
fi

#  Initialize conda/mamba environment containing necessary programs for
#+ subsetting peaks
check_env_installed "${env_name}"
activate_env "${env_name}"

#  Navigate to the work directory
cd "${dir_base}/${dir_repo}/${dir_work}" \
    || error_and_return "Failed to cd to ${dir_base}/${dir_repo}/${dir_work}."

#  Invoke the R interpreter to begin the work to subset peaks
R
```
</details>
<br />

<a id="ii-perform-set-operations-with-peak-intervals-returning-complete-intervals-from-a-given-set"></a>
### ii. Perform set operations with peak intervals, returning complete intervals from a given set
The following chunk of `R` code is designed for interactive execution and focuses on performing set operations with genomic peak intervals, specifically those derived from MACS3 `narrowPeak` files. It makes use of the `GenomicRanges` library to manipulate and analyze genomic intervals, ultimately returning completed, "reduced" intervals from a given set. Please adjust the assignments to variables `dir_narrowPeak`, `peaks_Esa1`, `peaks_Gcn5`, and `peaks_Rpd3` as necessary. Prior to running the code in this chunk, it is necessary to first run the code in chunk [i](#i-get-situated). Key components of this chunk:
- Read MACS3 `narrowPeak` files
- Perform set operations on genomic intervals
- Perform pairwise evaluations of the specified peak sets
- Output complete intervals for given peak sets
- Write results to `BED` files

<a id="r-code-1"></a>
#### R code
<details>
<summary><i>R code: 5.b.ii. Perform set operations with peak intervals, returning complete intervals from a given set</i></summary>

```r
#!/usr/bin/env Rscript

#  Define functions ===========================================================
#  Function to read a MACS3 narrowPeak file and encode it as a GRanges object
read_narrowPeak_to_GRanges <- function(file_path) {
    df <- readr::read_tsv(
        file_path,
        col_names = c(  # See biostars.org/p/102710/ for more information
            "chr", "start", "stop", "name", "score", "strand", "signal_value",
            "p_value", "q_value", "peak"
        ),
        show_col_types = FALSE
    ) %>%
        as.data.frame() %>%
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    
    return(df)
}


#  Function to intersect the genomic ranges of two peak sets, then return a
#+ single peak set of complete reduced ranges based on those intersections
intersect_and_reduce_two_sets_to_one_complete_set <- function(
    peaks_1,
    peaks_2,
    prefix
) {
    #  Intersect the genomic ranges represented by peaks_1 and peaks_2 to
    #+ capture any two-way intersections; store the intersected ranges in
    #+ intxns
    intxns <- GenomicRanges::intersect(peaks_1, peaks_2)
    
    if (length(intxns) > 0) {
        #  Find overlaps between peaks_1 and intxns, storing the resulting
        #+ complete ranges in intxn_peaks_1
        overlaps_peaks_1 <- GenomicRanges::findOverlaps(peaks_1, intxns)
        intxns_peaks_1 <- peaks_1[queryHits(overlaps_peaks_1)]
        
        #  Find overlaps between peaks_2 and intxns, storing the resulting
        #+ complete ranges in intxns_peaks_2
        overlaps_peaks_2 <- GenomicRanges::findOverlaps(peaks_2, intxns)
        intxns_peaks_2 <- peaks_2[queryHits(overlaps_peaks_2)]
    } else {
        stop(paste0(
            "No intersections found between ", peaks_1, " and ", peaks_2, "."
        ))
    }
    
    #  Reduce the combined intersecting ranges (intxns_peaks_1 and
    #+ intxns_peaks_2) to remove any redundant intervals
    reduced <- GenomicRanges::reduce(c(intxns_peaks_1, intxns_peaks_2))
    mcols(reduced)$name <- paste(
        prefix, seq_len(length(reduced)), sep = "_"
    )
    
    return(list(
        intxns = intxns,
        intxns_peaks_1 = intxns_peaks_1,
        intxns_peaks_2 = intxns_peaks_2,
        reduced = reduced
    ))
}


#  Function to obtain peaks_1 (which was, for example, output by function
#+ intersect_and_reduce_two_sets_to_one_complete_set) that either intersect
#+ peaks_2 (get_intersect = TRUE) or are asymmetrically different than peaks_2
#+ (get_intersect = FALSE)
get_complete_reduced_ranges <- function(
    peaks_1,
    peaks_2,
    get_intersect = TRUE,
    prefix
) {
    #  Obtain complete ranges from processed_list that either intersect peaks
    #+ (queryHits) or are asymmetrically different than peaks (-queryHits)
    if (get_intersect) {
        ranges_selected <- peaks_1[queryHits(
            GenomicRanges::findOverlaps(peaks_1, peaks_2)
        )]
    } else {
        ranges_selected <- peaks_1[-queryHits(
            GenomicRanges::findOverlaps(peaks_1, peaks_2)
        )]
    }
    
    #  Reduce the overlapping or asymmetric ranges
    ranges_reduced <- GenomicRanges::reduce(ranges_selected)
    mcols(ranges_reduced)$name <- paste(
        prefix, seq_len(length(ranges_reduced)), sep = "_"
    )
    
    return(list(
        ranges_selected = ranges_selected,
        ranges_reduced = ranges_reduced
    ))
}


#  Function to write BED file from GRanges object
write_BED_from_GRanges <- function(GRanges, file_path) {
    GRanges %>%
        as.data.frame() %>%
        dplyr::select(seqnames, start, end, name) %>%
        readr::write_tsv(file_path, col_names = FALSE)
}


#  Function to process a combination of two peak sets
process_combination_of_two_peak_sets <- function(
    peaks_1,
    peaks_2,
    peaks_1_name,
    peaks_2_name
) {
    prefix_and <- paste0(
        "peak-subset_complete-reduced_", peaks_1_name, "-and-", peaks_2_name
    )
    prefix_not <- paste0(
        "peak-subset_complete-reduced_", peaks_1_name, "-not-", peaks_2_name
    )
    
    list_and <- get_complete_reduced_ranges(
        peaks_1,
        peaks_2,
        get_intersect = TRUE,
        prefix = prefix_and
    )
    list_not <- get_complete_reduced_ranges(
        peaks_1,
        peaks_2,
        get_intersect = FALSE,
        prefix = prefix_not
    )
    
    #  Return the lists for further processing
    return(list(
        list_and = list_and,
        list_not = list_not
    ))
}


#  Load required libraries ====================================================
library(GenomicRanges)
library(tidyverse)


#  Read in data ===============================================================
dir_narrowPeak <- "03_bam/bowtie2/macs3"

peaks_Esa1 <- file.path(
    dir_narrowPeak,
    "Q_Esa5_7041.Q_Esa5_7691_peaks.narrowPeak.gz"
)
peaks_Gcn5 <- file.path(
    dir_narrowPeak,
    "Q_Gcn5_7692.Q_Gcn5_7709_peaks.narrowPeak.gz"
)
peaks_Rpd3 <- file.path(
    dir_narrowPeak,
    "Q_Rpd3_7568.Q_Rpd3_7569_peaks.narrowPeak.gz"
)

#  Load the narrowPeak files as GRanges and data.frame objects
peaks_Esa1 <- read_narrowPeak_to_GRanges(peaks_Esa1)
peaks_Gcn5 <- read_narrowPeak_to_GRanges(peaks_Gcn5)
peaks_Rpd3 <- read_narrowPeak_to_GRanges(peaks_Rpd3)


#  Do the main work ===========================================================
#  Calculate and write out complete ranges... ---------------------------------
#+ (i.e., peak intervals) via set intersections and set asymmetric differences
#+ between one peak set and a second peak set

#  Perform pairwise evaluations ---------------------------
#  Initialize a list of peak data sets, where each entry contains genomic
#+ ranges for a specific protein stored as GRanges objects
peaks_list <- list(
    Esa1 = peaks_Esa1,
    Gcn5 = peaks_Gcn5,
    Rpd3 = peaks_Rpd3
)

#  Initialize a master list to hold all results from the below loops
results_list <- list()

#  Loop through each combination of peak sets
for (i in 1:(length(peaks_list) - 1)) {      # i goes to second-to-last element
    for (j in (i + 1):length(peaks_list)) {  # j starts from element next to i
        peaks_1_name <- names(peaks_list)[i]
        peaks_2_name <- names(peaks_list)[j]
        peaks_1 <- peaks_list[[peaks_1_name]]
        peaks_2 <- peaks_list[[peaks_2_name]]
        
        #  Run `process_combination_of_two_peak_sets` for the initial peak-set
        #+ pair
        results_key_initial <- paste(peaks_1_name, peaks_2_name, sep = "-")
        results_list[[results_key_initial]] <-
            process_combination_of_two_peak_sets(
                peaks_1,
                peaks_2,
                peaks_1_name,
                peaks_2_name
            )
        
        #  Run `process_combination_of_two_peak_sets` for the reversed peak-set
        #+ pair
        results_key_reversed <- paste(peaks_2_name, peaks_1_name, sep = "-")
        results_list[[results_key_reversed]] <-
            process_combination_of_two_peak_sets(
                peaks_2,
                peaks_1,
                peaks_2_name,
                peaks_1_name
            )
    }
}

#  After processing all combinations and fully populating results_list, loop 
#+ through results_list to write files
for (key in names(results_list)) {
    #  Extract the complete, reduced set intersection ("and") and complete,
    #+ reduced set asymmetric difference ("not") ranges for the current
    #+ combination
    ranges_and <- results_list[[key]]$intersect_list$reduced_ranges
    ranges_not <- results_list[[key]]$not_list$reduced_ranges
    
    #  Define file name prefix and suffix
    prefix <- paste0(
        "peak-subset_complete-reduced_", stringr::str_split_i(key, "-", 1)
    )
    suffix <- paste0(unlist(strsplit(key, "-"))[2])

    #  Define file paths
    file_path_and <- file.path(
        dir_narrowPeak,
        paste0(prefix, "-and-", suffix, ".bed.gz")
    )
    file_path_not <- file.path(
        dir_narrowPeak,
        paste0(prefix, "-not-", suffix, ".bed.gz")
    )
    
    #  Write the BED files
    write_BED_from_GRanges(ranges_and, file_path_and)
    write_BED_from_GRanges(ranges_not, file_path_not)
}

#  Evaluate combined HAT intervals with Rpd3 intervals ----
#  Combine Esa1-Gcn5 (HAT) genomic ranges (complete, reduced) and analyze their
#+ intersections and asymmetric differences with Rpd3 peaks, returning partial
#+ peak ranges
Esa1_Gcn5_list <- intersect_and_reduce_two_sets_to_one_complete_set(
    peaks_1 = peaks_Esa1,
    peaks_2 = peaks_Gcn5,
    prefix = "peak-subset_complete-reduced_Esa1-Gcn5"
)

Esa1_Gcn5_then_Rpd3_combo <- process_combination_of_two_peak_sets(
    Esa1_Gcn5_list$reduced,
    peaks_Rpd3,
    "Esa1-Gcn5",
    "Rpd3"
)
Rpd3_then_Esa1_Gcn5_combo <- process_combination_of_two_peak_sets(
    peaks_Rpd3,
    Esa1_Gcn5_list$reduced,
    "Rpd3",
    "Esa1-Gcn5"
)

#  Write BED files for the above results
write_BED_from_GRanges(
    Esa1_Gcn5_list$reduced,
    file.path(
        dir_narrowPeak,
        "peak-subset_complete-reduced_Esa1-Gcn5.bed.gz"
    )
)

write_BED_from_GRanges(
    Esa1_Gcn5_then_Rpd3_combo$intersect_list$reduced_ranges,
    file.path(
        dir_narrowPeak,
        "peak-subset_complete-reduced_Esa1-Gcn5-and-Rpd3.bed.gz"
    )
)
write_BED_from_GRanges(
    Esa1_Gcn5_then_Rpd3_combo$not_list$reduced_ranges,
    file.path(
        dir_narrowPeak,
        "peak-subset_complete-reduced_Esa1-Gcn5-not-Rpd3.bed.gz"
    )
)

write_BED_from_GRanges(
    Rpd3_then_Esa1_Gcn5_combo$intersect_list$reduced_ranges,
    file.path(
        dir_narrowPeak,
        "peak-subset_complete-reduced_Rpd3-and-Esa1-Gcn5.bed.gz"
    )
)
write_BED_from_GRanges(
    Rpd3_then_Esa1_Gcn5_combo$not_list$reduced_ranges,
    file.path(
        dir_narrowPeak,
        "peak-subset_complete-reduced_Rpd3-not-Esa1-Gcn5.bed.gz"
    )
)
```
</details>
<br />

<a id="iii-perform-set-operations-with-peak-intervals-returning-partial-intervals-from-a-given-set"></a>
### iii. Perform set operations with peak intervals, returning partial intervals from a given set
The `R` code in this chunk functions similarly to the code in chunk [ii](#ii-perform-set-operations-with-peak-intervals-returning-complete-intervals-from-a-given-set), but it returns partial genomic intervals from a given set. Again, please adjust the assignments to variables `dir_narrowPeak`, `peaks_Esa1`, `peaks_Gcn5`, and `peaks_Rpd3` as necessary. Prior to running the code in this chunk, it is necessary to first run the code in chunk [i](#i-get-situated), but it is not necessary to have run the code in chunk [ii](#ii-perform-set-operations-with-peak-intervals-returning-complete-intervals-from-a-given-set).

<a id="r-code-2"></a>
#### R code
<details>
<summary><i>R code: 5.b.iii. Perform set operations with peak intervals, returning partial intervals from a given set</i></summary>

```r
#!/usr/bin/env Rscript

#  Define functions ===========================================================
#  Function to read a MACS3 narrowPeak file and encode it as a GRanges object
read_narrowPeak_to_GRanges <- function(file_path) {
    df <- readr::read_tsv(
        file_path,
        col_names = c(  # See biostars.org/p/102710/ for more information
            "chr", "start", "stop", "name", "score", "strand", "signal_value",
            "p_value", "q_value", "peak"
        ),
        show_col_types = FALSE
    ) %>%
        as.data.frame() %>%
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    
    return(df)
}


#  Function to return the intersection between two peak sets; returns a GRanges
#+ of intersecting intervals and thus partial intervals relative to the initial
#+ sets
get_set_intersection <- function(
    peaks_1,
    peaks_2,
    peaks_1_name,
    peaks_2_name
) {
    peaks_1_2_intxn <- GenomicRanges::intersect(peaks_1, peaks_2)
    mcols(peaks_1_2_intxn)$name <- paste0(
        "peak-subset_partial_", peaks_1_name, "-and-", peaks_2_name, "_",
        seq_len(length(peaks_1_2_intxn))
    )

    return(peaks_1_2_intxn)
}


#  Function to return the A ← B and A → B set assymetric differences between
#+ two peak sets; returns a list of two GRanges objects: one for A ← B 
#+ asymmetric differences and one for A → B asymmetric differences; the
#+ asymmetric differences are partial intervals relative to the initial
#+ sets
get_set_asymmetric_difference <- function(
    peaks_1,
    peaks_2,
    peaks_1_name,
    peaks_2_name
) {
    peaks_1_2_asym_diff <- GenomicRanges::setdiff(peaks_1, peaks_2)
    mcols(peaks_1_2_asym_diff)$name <- paste0(
        "peak-subset_partial_", peaks_1_name, "-not-", peaks_2_name, "_",
        seq_len(length(peaks_1_2_asym_diff))
    )

    return(peaks_1_2_asym_diff)
}


#  Function to write BED file from GRanges object
write_BED_from_GRanges <- function(GRanges, file_path) {
    GRanges %>%
        as.data.frame() %>%
        dplyr::select(seqnames, start, end, name) %>%
        readr::write_tsv(file_path, col_names = FALSE)
}


#  Load required libraries ====================================================
library(GenomicRanges)
library(tidyverse)


#  Read in data ===============================================================
dir_narrowPeak <- "03_bam/bowtie2/macs3"

peaks_Esa1 <- file.path(
    dir_narrowPeak,
    "Q_Esa5_7041.Q_Esa5_7691_peaks.narrowPeak.gz"
)
peaks_Gcn5 <- file.path(
    dir_narrowPeak,
    "Q_Gcn5_7692.Q_Gcn5_7709_peaks.narrowPeak.gz"
)
peaks_Rpd3 <- file.path(
    dir_narrowPeak,
    "Q_Rpd3_7568.Q_Rpd3_7569_peaks.narrowPeak.gz"
)

#  Load the narrowPeak files as GRanges and data.frame objects
peaks_Esa1 <- read_narrowPeak_to_GRanges(peaks_Esa1)
peaks_Gcn5 <- read_narrowPeak_to_GRanges(peaks_Gcn5)
peaks_Rpd3 <- read_narrowPeak_to_GRanges(peaks_Rpd3)


#  Do the main work ===========================================================
#  Calculate and write out partial ranges... ----------------------------------
#+ (i.e., peak intervals) via set intersections and set asymmetric differences
#+ between one peak set and a second peak set

#  Perform pairwise evaluations ---------------------------
#  Initialize a list of peak data sets, where each entry contains genomic
#+ ranges for a specific protein stored as GRanges objects
peaks_list <- list(
    Esa1 = peaks_Esa1,
    Gcn5 = peaks_Gcn5,
    Rpd3 = peaks_Rpd3
)

#  Initialize a master list to hold all results from the below loops
results_list <- list()

#  Loop through each combination of peak sets
for (i in 1:(length(peaks_list) - 1)) {      # i goes to the second-to-last element
    for (j in (i + 1):length(peaks_list)) {  # j starts from element next to i
        peaks_1_name <- names(peaks_list)[i]
        peaks_2_name <- names(peaks_list)[j]
        peaks_1 <- peaks_list[[peaks_1_name]]
        peaks_2 <- peaks_list[[peaks_2_name]]
        
        #  Run `get_set_intersection` on the peak-set pair
        results_key_intxn <- paste0(
            "intersection_", peaks_1_name, "-", peaks_2_name
        )
        results_list[[results_key_intxn]] <-
            get_set_intersection(
                peaks_1,
                peaks_2,
                peaks_1_name,
                peaks_2_name
            )
        
        #  Run `get_set_asymmetric_difference` on the peak set pair
        results_key_asym_diff_1_2 <- paste0(
            "asym-diff_", peaks_1_name, "-", peaks_2_name
        )
        results_list[[results_key_asym_diff_1_2]] <-
            get_set_asymmetric_difference(
                peaks_1,
                peaks_2,
                peaks_1_name,
                peaks_2_name
            )

        #  Run `get_set_asymmetric_difference` on the reverse-ordered peak set
        #+ pair
        results_key_asym_diff_2_1 <- paste0(
            "asym-diff_", peaks_2_name, "-", peaks_1_name
        )
        results_list[[results_key_asym_diff_2_1]] <-
            get_set_asymmetric_difference(
                peaks_2,
                peaks_1,
                peaks_2_name,
                peaks_1_name
            )
    }
}  # str(results_list, max.level = 2)

#  After processing all combinations and fully populating results_list, loop 
#+ through results_list to write files
for (key in names(results_list)) {
    #  Extract the given GRanges object
    ranges <- results_list[[key]]

    #  Determine appropriate strings for file name
    if (isTRUE(stringr::str_detect(key, "intersection_"))) {
        name <- stringr::str_replace(key, "intersection_", "")
        bool <- "and"
    } else {
        name <- stringr::str_replace(key, "asym-diff_", "")
        bool <- "not"
    }
    peaks_1_name <- stringr::str_split_i(name, "-", 1)
    peaks_2_name <- stringr::str_split_i(name, "-", 2)
    
    #  Define file name prefix and suffix
    prefix <- paste0("peak-subset_partial")
    suffix <- paste(peaks_1_name, bool, peaks_2_name, sep = "-")

    #  Define file paths
    file_path <- file.path(
        dir_narrowPeak,
        paste0(prefix, "_", suffix, ".bed.gz")
    )
    
    #  Write the BED files
    write_BED_from_GRanges(ranges, file_path)
}

#  Evaluate combined HAT intervals with Rpd3 intervals ----
#  Combine Esa1-Gcn5 (HAT) genomic ranges (partial) and analyze their
#+ intersections and asymmetric differences with Rpd3 peaks, returning partial
#+ peak ranges
#+ 
#+ First, intersect the ranges from Esa1_peaks and Gcn5_peaks
Esa1_Gcn5 <- get_set_intersection(
    peaks_Esa1,
    peaks_Gcn5,
    "Esa1",
    "Gcn5"
)

#  Next, Return portions of the Esa1-Gcn5 intersecting ranges (partial) that
#+ (a) intersect with Rpd3 peak ranges and (b) do not intersect with Rpd3 peak
#+ ranges and vice versa
#+ 
#+ (a) Intersect the Esa1-Gcn5 intersections with peaks_Rpd3, capturing
#+ intersections
Esa1_Gcn5_and_Rpd3 <- get_set_intersection(
    Esa1_Gcn5,
    peaks_Rpd3,
    "Esa1-Gcn5",
    "Rpd3"
)

#  (b) Obtain the set asymmetric difference of Esa1_Gcn5 with respect to
#+     peaks_Rpd3 and vice versa
Esa1_Gcn5_not_Rpd3 <- get_set_asymmetric_difference(
    Esa1_Gcn5,
    peaks_Rpd3,
    "Esa1-Gcn5",
    "Rpd3"
)
Rpd3_not_Esa1_Gcn5 <- get_set_asymmetric_difference(
    peaks_Rpd3,
    Esa1_Gcn5,
    "Rpd3",
    "Esa1-Gcn5"
)

#  Write out BED files for the above set operations
write_BED_from_GRanges(
    Esa1_Gcn5_and_Rpd3,
    file.path(
        dir_narrowPeak,
        "peak-subset_partial_Esa1-Gcn5-and-Rpd3.bed.gz"
    )
)

write_BED_from_GRanges(
    Esa1_Gcn5_not_Rpd3,
    file.path(
        dir_narrowPeak,
        "peak-subset_partial_Esa1-Gcn5-not-Rpd3.bed.gz"
    )
)
write_BED_from_GRanges(
    Rpd3_not_Esa1_Gcn5,
    file.path(
        dir_narrowPeak,
        "peak-subset_partial_Rpd3-not-Esa1-Gcn5.bed.gz"
    )
)
```
</details>
<br />

<a id="iv-perform-additional-set-operations-with-respect-to-three-way-intersections"></a>
### iv. Perform additional set operations with respect to three-way intersections
The `R` code in this chunk functions similarly to the code in chunks [ii](#ii-perform-set-operations-with-peak-intervals-returning-complete-intervals-from-a-given-set) and [iii](#iii-perform-set-operations-with-peak-intervals-returning-partial-intervals-from-a-given-set), but it focuses not on output from pairwise comparisons but from three-way intersections. Again, please adjust the assignments to variables `dir_narrowPeak`, `peaks_Esa1`, `peaks_Gcn5`, and `peaks_Rpd3` as necessary. Prior to running the code in this chunk, it is necessary to first run the code in chunk [i](#i-get-situated), but it is not necessary to have run the code in chunks [ii](#ii-perform-set-operations-with-peak-intervals-returning-complete-intervals-from-a-given-set) and [iii](#iii-perform-set-operations-with-peak-intervals-returning-partial-intervals-from-a-given-set).


<a id="r-code-3"></a>
#### R code
<details>
<summary><i>R code: 5.b.iv. Perform additional set operations with respect to three-way intersections</i></summary>

```R
#!/usr/bin/env Rscript

#  Define functions ===========================================================
#  Function to read a MACS3 narrowPeak file and encode it as a GRanges object
read_narrowPeak_to_GRanges <- function(file_path) {
    df <- readr::read_tsv(
        file_path,
        col_names = c(  # See biostars.org/p/102710/ for more information
            "chr", "start", "stop", "name", "score", "strand", "signal_value",
            "p_value", "q_value", "peak"
        ),
        show_col_types = FALSE
    ) %>%
        as.data.frame() %>%
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    
    return(df)
}


#  Function to write BED file from GRanges object
write_BED_from_GRanges <- function(GRanges, file_path) {
    GRanges %>%
        as.data.frame() %>%
        dplyr::select(seqnames, start, end, name) %>%
        readr::write_tsv(file_path, col_names = FALSE)
}


#  Load required libraries ====================================================
library(GenomicRanges)
library(tidyverse)


#  Read in data ===============================================================
dir_narrowPeak <- "03_bam/bowtie2/macs3"

peaks_Esa1 <- file.path(
    dir_narrowPeak,
    "Q_Esa5_7041.Q_Esa5_7691_peaks.narrowPeak.gz"
)
peaks_Gcn5 <- file.path(
    dir_narrowPeak,
    "Q_Gcn5_7692.Q_Gcn5_7709_peaks.narrowPeak.gz"
)
peaks_Rpd3 <- file.path(
    dir_narrowPeak,
    "Q_Rpd3_7568.Q_Rpd3_7569_peaks.narrowPeak.gz"
)

#  Load the narrowPeak files as GRanges and data.frame objects
peaks_Esa1 <- read_narrowPeak_to_GRanges(peaks_Esa1)
peaks_Gcn5 <- read_narrowPeak_to_GRanges(peaks_Gcn5)
peaks_Rpd3 <- read_narrowPeak_to_GRanges(peaks_Rpd3)


#  Do the main work ===========================================================
#  (a) Obtain ONLY THE PORTIONS of peaks involved in a three-way intersection.
#+ (b) Identify the complete, specific intervals from each individual peak
#+ set that overlap with the three-way peak portion set, and then combine and
#+ reduce these overlapping intervals. And (c) obtain a combined, reduced peak
#+ set from all three peak sets without concern for what overlaps what and also
#+ without concern for retaining or discarding anything; i.e., reduction occurs
#+ only when overlaps are present, and non-overlapping intervals are retained
#+ in this peak set too.

#  (a) ------------------------------------------------------------------------
#  Capture any three-way intersections, i.e., only the portions of peak sets
#+ involved in a three-way intersection
intxn_three_way <- GenomicRanges::intersect(
    GenomicRanges::intersect(peaks_Esa1, peaks_Gcn5), peaks_Rpd3
)
mcols(intxn_three_way)$name <- paste0(
    "peak-subset_partial_three-way-intersection_",
    seq_len(length(intxn_three_way))
)

#  Write out BED file for the above set operation
write_BED_from_GRanges(
    intxn_three_way,
    file.path(
        dir_narrowPeak,
        "peak-subset_partial_three-way-intersection.bed.gz"
    )
)


#  (b) ------------------------------------------------------------------------
#  Return complete intervals for each of the three peak sets (Esa1, Gcn5, and
#+ Rpd3) based on the overlap of each of the three peaks sets with the GRanges
#+ object representing the three-way intersections
if (length(intxn_three_way) > 0) {
    intxn_three_way_Esa1 <- peaks_Esa1[
        queryHits(GenomicRanges::findOverlaps(
            peaks_Esa1, intxn_three_way
        ))
    ]
    
    intxn_three_way_Gcn5 <- peaks_Gcn5[
        queryHits(GenomicRanges::findOverlaps(
            peaks_Gcn5, intxn_three_way
        ))
    ]
    
    intxn_three_way_Rpd3 <- peaks_Rpd3[
        queryHits(GenomicRanges::findOverlaps(
            peaks_Rpd3, intxn_three_way
        ))
    ]
} else {
    stop("There are no three-way intersections.")
}

#  Generate a final GRanges object of reduced complete intervals derived
#+ from the three-way intersections of Esa1, Gcn5, and Rpd3
intxn_three_way_reduced <- GenomicRanges::reduce(c(
    intxn_three_way_Esa1,
    intxn_three_way_Gcn5,
    intxn_three_way_Rpd3
))
mcols(intxn_three_way_reduced)$name <- paste0(
    "peak-subset_complete-reduced_three-way-intersection_",
    seq_len(length(intxn_three_way_reduced))
)

#  Write out BED files for the above set operations
write_BED_from_GRanges(
    intxn_three_way_reduced,
    file.path(
        dir_narrowPeak,
        "peak-subset_complete-reduced_three-way-intersection.bed.gz"
    )
)

write_BED_from_GRanges(
    intxn_three_way_Esa1,
    file.path(
        dir_narrowPeak,
        "peak-subset_complete_Esa1-and-three-way-intersection.bed.gz"
    )
)
write_BED_from_GRanges(
    intxn_three_way_Gcn5,
    file.path(
        dir_narrowPeak,
        "peak-subset_complete_Gcn5-and-three-way-intersection.bed.gz"
    )
)
write_BED_from_GRanges(
    intxn_three_way_Rpd3,
    file.path(
        dir_narrowPeak,
        "peak-subset_complete_Rpd3-and-three-way-intersection.bed.gz"
    )
)


#  (c) ------------------------------------------------------------------------
#  Generate a GRanges object that is a combined, reduced peak set from all
#+ three peak sets; this means that reduction occurs only when overlaps are
#+ present, and non-overlapping intervals are retained in this peak set too.
peaks_all_reduced <- GenomicRanges::reduce(c(
    peaks_Esa1, peaks_Gcn5, peaks_Rpd3
))
mcols(peaks_all_reduced)$name <- paste0(
    "peak-subset_all-reduced_Esa1-Gcn5-Rpd3_",
    seq_len(length(peaks_all_reduced))
)

#  Write out BED file for the above set operation
write_BED_from_GRanges(
    peaks_all_reduced,
    file.path(
        dir_narrowPeak,
        "peak-subset_all-reduced_Esa1-Gcn5-Rpd3.bed.gz"
    )
)
```
</details>
<br />

<a id="6-calculate-sample-scaling-factors-from-s-pombe-spike-ins"></a>
# 6. Calculate sample scaling factors from *S. pombe* spike-ins
`#TODO`
<br />
<br />

<a id="7-miscellaneous-to-be-organized"></a>
# 7. Miscellaneous (to be organized)
<a id="x-scratch"></a>
## x. Scratch
<a id="code-13"></a>
### Code
<details>
<summary><i>Code: Scratch</i></summary>

```bash
#!/bin/bash

find . \
    -type f \
    -name 'subset-peaks_complete-reduced_*.bed.gz' \
    -print

find . \
    -type f \
    -name 'subset-peaks_complete-reduced_*.bed.gz' \
    -exec sh \
        -c 'mv "$0" "${0%/*}/bak.${0##*/}"' {} \;

find . \
    -type f \
    -name 'Ch_*.atria.*' \
    -exec bash \
        -c 'mv "$0" "${0/\/Ch_/\/IP_}"' {} \;
```
</details>
<br />

<a id="notes"></a>
### Notes
<details>
<summary><i>Notes: Scratch</i></summary>
<br />

Here's the breakdown of the command:
- `find . -type f -name 'subset-peaks_complete-reduced_*.bed.gz'`: This finds all files in the current directory (and subdirectories) that match the pattern 'subset-peaks_complete-reduced_*.bed.gz'.
- `-exec sh -c 'mv "$0" "${0%/*}/bak.${0##*/}"' {} \;`: For each file found, this executes a shell command to move (`mv`) the file. `$0` represents the found file (the `{}` placeholder is passed to the shell command as `$0`).
- `"${0%/*}/bak.${0##*/}"`: This constructs the new filename by appending "bak." to the basename of the file. `${0%/*}` extracts the directory part of the found file's path, and `${0##*/}` extracts the basename (the part after the last `/`). By combining these with `/bak.`, you prepend "bak." to the basename.
</details>
<br />

<a id="a-determine-the-locations-of-low-complexity-regions-in-s-cerevisiae"></a>
## a. Determine the locations of low-complexity regions in *S. cerevisiae*
<a id="i-install-sdust-via-minimap"></a>
### i. Install [`sdust`](https://pubmed.ncbi.nlm.nih.gov/16796549/) via [`minimap`](https://github.com/lh3/minimap/tree/master)
<a id="code-14"></a>
#### Code
<details>
<summary><i>Code: i. Install `sdust` via `minimap`</i></summary>

```bash
#!/bin/bash

#  Define functions ===========================================================
#  Function to return an error message and exit code 1, which stops the
#+ interactive execution of code
function error_and_return() {
    echo "Error: ${1}" >&2
    return 1
}


#  Function to check if Mamba is installed
function check_mamba_installed() {
    if ! : mamba &> /dev/null; then
        echo "Mamba is not installed on your system. Mamba is a package manager" 
        echo "that makes package installations faster and more reliable."
        echo ""
        echo "For installation instructions, please check the following link:"
        echo "https://github.com/mamba-org/mamba#installation"
        return 1
    fi
    
    return 0
}


function check_env_installed() {
    local env_name="${1}"

    if conda env list | grep -q "^${env_name} "; then
        return 0
    else
        return 1
    fi
}


#  Initialize variables and arrays ============================================
env_name="alignment-processing_env"


#  Do the main work ===========================================================
#  Set flag(s)
create_mamba_env=true  # Install mamba environment if not detected

#  Check that Mamba is installed and in PATH
check_mamba_installed

#  Check that environment assigned to env_name is installed; if environment
#+ assigned to env_name is not installed, run the following
if check_env_installed "${env_name}"; then
    #  Handle the case when the environment is already installed
    echo "Activating environment ${env_name}"
    
    #TODO Make the following a function
    if ! mamba activate "${env_name}" &> /dev/null; then
        #  If `mamba activate` fails, try using `source activate`
        if ! conda activate "${env_name}" &> /dev/null; then
            if ! source activate "${env_name}" &> /dev/null; then
                #  If `source activate` also fails, return an error
                error_and_return "Failed to activate environment \"${env_name}\"."
            fi
        fi
    else
        echo "Environment \"${env_name}\" activated using mamba."
    fi
else
    #  Handle the case when the environment is not installed
    echo "Creating environment ${env_name}"
    
    if ${create_mamba_env}; then
        #  Switch `--yes` is set, which means no user input is required
        #NOTE Running this on FHCC Rhino; ergo, no CONDA_SUBDIR=osx-64
        mamba create \
            --yes \
            --name "${env_name}" \
            --channel bioconda \
                bamtools \
                bedtools \
                bowtie2 \
                bwa \
                fastqc \
                minimap \
                mosdepth \
                picard \
                preseq \
                samtools \
                ucsc-bedgraphtobigwig \
                ucsc-facount
        
        # source activate "${env_name}"
    fi
fi

# mamba create \
#     --yes \
#     --name minimap_env \
#     --channel bioconda \
#         minimap

# source activate "alignment-processing_env"
# mamba install \
#     --yes \
#     --channel bioconda \
#         ucsc-facount
```
</details>
<br />

<a id="printed"></a>
#### Printed
<details>
<summary><i>Printed: i. Install `sdust` via `minimap`</i></summary>

```txt
❯ mamba create \
>     --yes \
>     --name minimap_env \
>     --channel bioconda \
>         minimap

                  __    __    __    __
                 /  \  /  \  /  \  /  \
                /    \/    \/    \/    \
███████████████/  /██/  /██/  /██/  /████████████████████████
              /  / \   / \   / \   / \  \____
             /  /   \_/   \_/   \_/   \    o \__,
            / _/                       \_____/  `
            |/
        ███╗   ███╗ █████╗ ███╗   ███╗██████╗  █████╗
        ████╗ ████║██╔══██╗████╗ ████║██╔══██╗██╔══██╗
        ██╔████╔██║███████║██╔████╔██║██████╔╝███████║
        ██║╚██╔╝██║██╔══██║██║╚██╔╝██║██╔══██╗██╔══██║
        ██║ ╚═╝ ██║██║  ██║██║ ╚═╝ ██║██████╔╝██║  ██║
        ╚═╝     ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝╚═════╝ ╚═╝  ╚═╝

        mamba (1.3.1) supported by @QuantStack

        GitHub:  https://github.com/mamba-org/mamba
        Twitter: https://twitter.com/QuantStack

█████████████████████████████████████████████████████████████


Looking for: ['minimap']

bioconda/linux-64                                           Using cache
bioconda/noarch                                             Using cache
conda-forge/noarch                                          Using cache
pkgs/main/noarch                                              No change
pkgs/main/linux-64                                            No change
pkgs/r/linux-64                                               No change
pkgs/r/noarch                                                 No change
conda-forge/linux-64                                37.7MB @   6.8MB/s  6.4s
Transaction

  Prefix: /home/kalavatt/miniconda3/envs/minimap_env

  Updating specs:

   - minimap


  Package          Version  Build        Channel                    Size
──────────────────────────────────────────────────────────────────────────
  Install:
──────────────────────────────────────────────────────────────────────────

  + _libgcc_mutex      0.1  conda_forge  conda-forge/linux-64     Cached
  + _openmp_mutex      4.5  2_gnu        conda-forge/linux-64     Cached
  + libgcc           7.2.0  h69d50b8_2   conda-forge/linux-64     Cached
  + libgcc-ng       13.2.0  h807b86a_4   conda-forge/linux-64      774kB
  + libgomp         13.2.0  h807b86a_4   conda-forge/linux-64      422kB
  + libstdcxx-ng    13.2.0  h7e041cc_4   conda-forge/linux-64        4MB
  + minimap            0.2  0            bioconda/linux-64         112kB

  Summary:

  Install: 7 packages

  Total download: 5MB

──────────────────────────────────────────────────────────────────────────


libgomp                                            422.1kB @   2.9MB/s  0.1s
libgcc-ng                                          773.8kB @   5.3MB/s  0.2s
libstdcxx-ng                                         3.8MB @  20.2MB/s  0.2s
minimap                                            112.5kB @ 316.6kB/s  0.4s

Downloading and Extracting Packages

Preparing transaction: done
Verifying transaction: done
Executing transaction: done

To activate this environment, use

     $ mamba activate minimap_env

To deactivate an active environment, use

     $ mamba deactivate


❯ mamba install \
>     --yes \
>     --channel bioconda \
>         ucsc-facount

                  __    __    __    __
                 /  \  /  \  /  \  /  \
                /    \/    \/    \/    \
███████████████/  /██/  /██/  /██/  /████████████████████████
              /  / \   / \   / \   / \  \____
             /  /   \_/   \_/   \_/   \    o \__,
            / _/                       \_____/  `
            |/
        ███╗   ███╗ █████╗ ███╗   ███╗██████╗  █████╗
        ████╗ ████║██╔══██╗████╗ ████║██╔══██╗██╔══██╗
        ██╔████╔██║███████║██╔████╔██║██████╔╝███████║
        ██║╚██╔╝██║██╔══██║██║╚██╔╝██║██╔══██╗██╔══██║
        ██║ ╚═╝ ██║██║  ██║██║ ╚═╝ ██║██████╔╝██║  ██║
        ╚═╝     ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝╚═════╝ ╚═╝  ╚═╝

        mamba (1.3.1) supported by @QuantStack

        GitHub:  https://github.com/mamba-org/mamba
        Twitter: https://twitter.com/QuantStack

█████████████████████████████████████████████████████████████


Looking for: ['ucsc-facount']

bioconda/noarch                                      5.1MB @   4.1MB/s  1.5s
bioconda/linux-64                                    5.3MB @   3.6MB/s  1.7s
pkgs/main/linux-64                                   6.7MB @   3.8MB/s  2.1s
pkgs/r/linux-64                                               No change
pkgs/r/noarch                                                 No change
pkgs/main/noarch                                   860.2kB @ 370.0kB/s  0.8s
conda-forge/noarch                                  15.9MB @   5.3MB/s  3.4s
conda-forge/linux-64                                38.4MB @   5.4MB/s  7.8s

Pinned packages:
  - python 3.10.*


Transaction

  Prefix: /home/kalavatt/miniconda3/envs/alignment-processing_env

  Updating specs:

   - ucsc-facount
   - ca-certificates
   - openssl


  Package                  Version  Build                Channel                    Size
──────────────────────────────────────────────────────────────────────────────────────────
  Install:
──────────────────────────────────────────────────────────────────────────────────────────

  + gsl                        2.7  he838d99_0           conda-forge/linux-64     Cached
  + libcblas                 3.9.0  21_linux64_openblas  conda-forge/linux-64       15kB
  + ucsc-facount               377  ha8a8165_3           bioconda/linux-64         137kB

  Change:
──────────────────────────────────────────────────────────────────────────────────────────

  - lcms2                     2.15  h7f713cb_2           conda-forge
  + lcms2                     2.15  haa2dc70_1           conda-forge/linux-64     Cached
  - libcups                  2.3.3  h4637d8d_4           conda-forge
  + libcups                  2.3.3  h36d4200_3           conda-forge/linux-64        5MB
  - mysql-connector-c       6.1.11  h659d440_1008        conda-forge
  + mysql-connector-c       6.1.11  h6eb9d5d_1007        conda-forge/linux-64     Cached
  - pango                  1.50.14  ha41ecd1_2           conda-forge
  + pango                  1.50.14  heaa33ce_1           conda-forge/linux-64     Cached

  Downgrade:
──────────────────────────────────────────────────────────────────────────────────────────

  - cairo                   1.18.0  h3faef2a_0           conda-forge
  + cairo                   1.16.0  hbbf8b49_1016        conda-forge/linux-64     Cached
  - curl                     8.5.0  hca28451_0           conda-forge
  + curl                    7.88.1  h37d81fd_2           pkgs/main/linux-64         82kB
  - harfbuzz                 8.3.0  h3d44ed6_0           conda-forge
  + harfbuzz                 7.3.0  hdb3a94d_0           conda-forge/linux-64     Cached
  - htslib                  1.19.1  h81da01d_1           bioconda
  + htslib                    1.17  h6bc39ce_1           bioconda/linux-64           2MB
  - icu                       73.2  h59595ed_0           conda-forge
  + icu                       72.1  hcb278e6_0           conda-forge/linux-64     Cached
  - krb5                    1.21.2  h659d440_0           conda-forge
  + krb5                    1.20.1  hf9c8cef_0           conda-forge/linux-64     Cached
  - libcurl                  8.5.0  hca28451_0           conda-forge
  + libcurl                 7.88.1  h91b91d3_2           pkgs/main/linux-64        393kB
  - libnghttp2              1.58.0  h47da74e_1           conda-forge
  + libnghttp2              1.52.0  ha637b67_1           pkgs/main/linux-64        687kB
  - libssh2                 1.11.0  h0841786_0           conda-forge
  + libssh2                 1.10.0  haa6b8db_3           conda-forge/linux-64     Cached
  - libtiff                  4.6.0  h8b53f26_0           conda-forge
  + libtiff                  4.5.1  h8b53f26_1           conda-forge/linux-64      417kB
  - libxml2                 2.12.5  h232c23b_0           conda-forge
  + libxml2                 2.11.5  h0d562d8_0           conda-forge/linux-64      705kB
  - mosdepth                 0.3.6  hd299d5a_0           bioconda
  + mosdepth                 0.3.3  hd299d5a_3           bioconda/linux-64        Cached
  - openjdk                 20.0.2  hfea2f88_1           conda-forge
  + openjdk                 11.0.1  h516909a_1016        conda-forge/linux-64      184MB
  - openssl                  3.2.1  hd590300_0           conda-forge
  + openssl                 1.1.1w  hd590300_0           conda-forge/linux-64        2MB
  - picard                   3.1.1  hdfd78af_0           bioconda
  + picard                   3.0.0  hdfd78af_0           bioconda/noarch            18MB
  - python                 3.10.13  hd12c33a_1_cpython   conda-forge
  + python                  3.10.8  h257c98d_0_cpython   conda-forge/linux-64     Cached
  - r-base                   4.3.1  h639d9d3_5           conda-forge
  + r-base                   4.2.3  h4a03800_2           conda-forge/linux-64       25MB
  - samtools                1.19.2  h50ea8bc_0           bioconda
  + samtools                  1.18  hd87286a_0           bioconda/linux-64         466kB
  - ucsc-bedgraphtobigwig      455  h2a80c09_0           bioconda
  + ucsc-bedgraphtobigwig      445  h954228d_0           bioconda/linux-64           2MB

  Summary:

  Install: 3 packages
  Change: 4 packages
  Downgrade: 19 packages

  Total download: 242MB

──────────────────────────────────────────────────────────────────────────────────────────


libcblas                                            14.6kB @ 115.9kB/s  0.1s
libtiff                                            416.5kB @   2.6MB/s  0.2s
openssl                                              2.0MB @  12.0MB/s  0.2s
libxml2                                            705.0kB @   4.3MB/s  0.2s
libcups                                              4.5MB @  20.5MB/s  0.2s
libnghttp2                                         687.2kB @   2.5MB/s  0.1s
samtools                                           466.2kB @   1.2MB/s  0.1s
libcurl                                            392.6kB @   1.0MB/s  0.2s
ucsc-facount                                       136.7kB @ 322.1kB/s  0.3s
curl                                                82.5kB @ 176.2kB/s  0.1s
htslib                                               2.5MB @   3.9MB/s  0.2s
r-base                                              25.1MB @  34.6MB/s  0.6s
ucsc-bedgraphtobigwig                                2.3MB @   2.4MB/s  0.6s
picard                                              18.3MB @  12.3MB/s  1.1s
openjdk                                            184.0MB @  90.6MB/s  2.0s

Downloading and Extracting Packages

Preparing transaction: done
Verifying transaction: done
Executing transaction: done
```
</details>
<br />

<a id="ii-run-sdust-via-minimap"></a>
### ii. Run [`sdust`](https://pubmed.ncbi.nlm.nih.gov/16796549/) via [`minimap`](https://github.com/lh3/minimap/tree/master)
<a id="code-15"></a>
#### Code
<details>
<summary><i>Code: ii. Run `sdust` via `minimap`</i></summary>

```bash
#!/bin/bash

grabnode  # 1, 20, 1, N


#  Define functions ===========================================================
#  ...

#  Initialize variables and arrays ============================================
d_fa="${HOME}/genomes/Saccharomyces_cerevisiae/fasta-processed"
f_fa="S288C_reference_sequence_R64-3-1_20210421.fa.gz"
a_fa="${d_fa}/${f_fa}"

f_bed="S288C_reference_sequence_R64-3-1_20210421.low_complexity.bed"
a_bed="${d_fa}/${f_bed}"

env_name="minimap_env"


#  Do the main work ===========================================================
#  Set flag(s)
check_variables=true
check_file_exists=true
check_command=true
run_command=true

#  Check variables
if ${check_variables}; then
    echo "
           d_fa=${d_fa}
           f_fa=${f_fa}
           a_fa=${a_fa}

          f_bed=${f_bed}
          a_bed=${a_bed}

       env_name=${env_name}
    "
fi

#  Check that input FASTQ file exists
if ${check_file_exists}; then ls -lhaFG "${a_fa}"; fi

#  If not already activated, the activate conda environment
if [[ "${CONDA_DEFAULT_ENV}" != "${env_name}" ]]; then
    if [[ ${CONDA_DEFAULT_ENV} != "base" ]]; then
        conda deactivate
    fi

    source activate "${env_name}"
fi

#  Create a BED file for S. cerevisiae regions of low complexity
if ${check_command}; then
    echo "
    if [[ ! -f \"${a_bed}\" ]]; then
        if [[ ! -f "${a_fa%.gz}" ]]; then
            gzip -cd \\
                \"${a_fa}\" \\
                    > \"${a_fa%.gz}\"
        fi
        
        if [[ -f \"${a_fa%.gz}\" ]]; then
            sdust \\
                \"${a_fa}\" \\
                    > \"${a_bed}\"
        fi
    fi
    "
fi

if ${run_command}; then
    if [[ ! -f "${a_bed}" ]]; then
        if [[ ! -f "${a_fa%.gz}" ]]; then
            gzip -cd "${a_fa}" > "${a_fa%.gz}"
        fi

        if [[ -f "${a_fa%.gz}" ]]; then
            sdust "${a_fa%.gz}" > "${a_bed}"
        fi
    fi
fi
```
</details>
<br />

<a id="b-determine-the-effective-genome-size-of-s-cerevisiae-50-mers"></a>
## b. Determine the effective genome size of *S. cerevisiae* (50-mers)
<a id="i-install-khmer"></a>
### i. Install [`khmer`](https://khmer.readthedocs.io/en/latest/)
<a id="code-16"></a>
#### Code
<details>
<summary><i>Code: i. Install `khmer`</i></summary>

```bash
#!/bin/bash

#  Define functions ===========================================================
#  Function to return an error message and exit code 1, which stops the
#+ interactive execution of code
function error_and_return() {
    echo "Error: ${1}" >&2
    return 1
}


#  Function to check if Mamba is installed
function check_mamba_installed() {
    if ! : mamba &> /dev/null; then
        echo "Mamba is not installed on your system. Mamba is a package manager" 
        echo "that makes package installations faster and more reliable."
        echo ""
        echo "For installation instructions, please check the following link:"
        echo "https://github.com/mamba-org/mamba#installation"
        return 1
    fi
    
    return 0
}


function check_env_installed() {
    local env_name="${1}"

    if conda env list | grep -q "^${env_name} "; then
        return 0
    else
        echo "Environment \"${env_name}\" is not installed."
        return 1
    fi
}


#  Initialize variables and arrays ============================================
env_name="khmer_env"


#  Do the main work ===========================================================
#  Install minimap ------------------------------------------------------------
#  Set flag(s)
create_mamba_env=true  # Install mamba environment if not detected

#  Check that Mamba is installed and in PATH
check_mamba_installed

#  Check that environment assigned to env_name is installed; if environment
#+ assigned to env_name is not installed, run the following
if [[ $(check_env_installed "${env_name}") -eq 0 ]]; then
    #  Handle the case when the environment is already installed
    echo "Activating environment ${env_name}"
    
    if ! mamba activate "${env_name}" &> /dev/null; then
        #  If `mamba activate` fails, try using `source activate`
        if ! conda activate "${env_name}" &> /dev/null; then
            if ! source activate "${env_name}" &> /dev/null; then
                #  If `source activate` also fails, return an error
                error_and_return "Failed to activate environment \"${env_name}\"."
            fi
        fi
    else
        echo "Environment \"${env_name}\" activated using mamba."
    fi
else
    #  Handle the case when the environment is not installed
    echo "Creating environment ${env_name}"
    
    if ${create_mamba_env}; then
        #  Switch `--yes` is set, which means no user input is required
        #NOTE Running this on FHCC Rhino; ergo, no CONDA_SUBDIR=osx-64
        mamba create \
            --yes \
            --name "${env_name}" \
            --channel bioconda \
                khmer
    fi
fi

#TODO Check this error message: -bash: [[: Environment "khmer_env" is not installed.: syntax error: invalid arithmetic operator (error token is ""khmer_env" is not installed.")
```
</details>
<br />

<a id="printed-1"></a>
#### Printed
<details>
<summary><i>Printed: ii. Install `khmer`</i></summary>

```txt
❯ mamba create \
>     --yes \
>     --name khmer_env \
>     --channel bioconda \
>         khmer

                  __    __    __    __
                 /  \  /  \  /  \  /  \
                /    \/    \/    \/    \
███████████████/  /██/  /██/  /██/  /████████████████████████
              /  / \   / \   / \   / \  \____
             /  /   \_/   \_/   \_/   \    o \__,
            / _/                       \_____/  `
            |/
        ███╗   ███╗ █████╗ ███╗   ███╗██████╗  █████╗
        ████╗ ████║██╔══██╗████╗ ████║██╔══██╗██╔══██╗
        ██╔████╔██║███████║██╔████╔██║██████╔╝███████║
        ██║╚██╔╝██║██╔══██║██║╚██╔╝██║██╔══██╗██╔══██║
        ██║ ╚═╝ ██║██║  ██║██║ ╚═╝ ██║██████╔╝██║  ██║
        ╚═╝     ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝╚═════╝ ╚═╝  ╚═╝

        mamba (1.3.1) supported by @QuantStack

        GitHub:  https://github.com/mamba-org/mamba
        Twitter: https://twitter.com/QuantStack

█████████████████████████████████████████████████████████████


Looking for: ['khmer']

bioconda/linux-64                                           Using cache
bioconda/noarch                                             Using cache
conda-forge/linux-64                                        Using cache
conda-forge/noarch                                          Using cache
pkgs/main/linux-64                                            No change
pkgs/main/noarch                                              No change
pkgs/r/linux-64                                               No change
pkgs/r/noarch                                                 No change
Transaction

  Prefix: /home/kalavatt/miniconda3/envs/khmer_env

  Updating specs:

   - khmer


  Package                Version  Build               Channel                    Size
───────────────────────────────────────────────────────────────────────────────────────
  Install:
───────────────────────────────────────────────────────────────────────────────────────

  + _libgcc_mutex            0.1  conda_forge         conda-forge/linux-64     Cached
  + _openmp_mutex            4.5  2_gnu               conda-forge/linux-64     Cached
  + bz2file                 0.98  py_0                conda-forge/noarch          9kB
  + bzip2                  1.0.8  hd590300_5          conda-forge/linux-64     Cached
  + ca-certificates   2023.11.17  hbcca054_0          conda-forge/linux-64     Cached
  + khmer                3.0.0a3  py38h94ffb2d_3      bioconda/linux-64          11MB
  + ld_impl_linux-64        2.40  h41732ed_0          conda-forge/linux-64     Cached
  + libffi                 3.4.2  h7f98852_5          conda-forge/linux-64     Cached
  + libgcc-ng             13.2.0  h807b86a_4          conda-forge/linux-64     Cached
  + libgomp               13.2.0  h807b86a_4          conda-forge/linux-64     Cached
  + libnsl                 2.0.1  hd590300_0          conda-forge/linux-64     Cached
  + libsqlite             3.44.2  h2797004_0          conda-forge/linux-64     Cached
  + libstdcxx-ng          13.2.0  h7e041cc_4          conda-forge/linux-64     Cached
  + libuuid               2.38.1  h0b41bf4_0          conda-forge/linux-64     Cached
  + libxcrypt             4.4.36  hd590300_1          conda-forge/linux-64      100kB
  + libzlib               1.2.13  hd590300_5          conda-forge/linux-64     Cached
  + ncurses                  6.4  h59595ed_2          conda-forge/linux-64     Cached
  + openssl                3.2.1  hd590300_0          conda-forge/linux-64        3MB
  + pip                   23.3.2  pyhd8ed1ab_0        conda-forge/noarch          1MB
  + python                3.8.18  hd12c33a_1_cpython  conda-forge/linux-64       24MB
  + python_abi               3.8  4_cp38              conda-forge/linux-64     Cached
  + readline                 8.2  h8228510_1          conda-forge/linux-64     Cached
  + screed                 1.0.4  py_0                bioconda/noarch            83kB
  + setuptools            69.0.3  pyhd8ed1ab_0        conda-forge/noarch        471kB
  + tk                    8.6.13  noxft_h4845f30_101  conda-forge/linux-64     Cached
  + wheel                 0.42.0  pyhd8ed1ab_0        conda-forge/noarch       Cached
  + xz                     5.2.6  h166bdaf_0          conda-forge/linux-64     Cached
  + zlib                  1.2.13  hd590300_5          conda-forge/linux-64     Cached

  Summary:

  Install: 28 packages

  Total download: 40MB

───────────────────────────────────────────────────────────────────────────────────────


pip                                                  1.4MB @  17.7MB/s  0.1s
libxcrypt                                          100.4kB @ 989.8kB/s  0.1s
setuptools                                         470.5kB @   3.4MB/s  0.1s
openssl                                              2.9MB @  16.7MB/s  0.2s
bz2file                                              8.8kB @  44.5kB/s  0.1s
screed                                              83.2kB @ 331.6kB/s  0.1s
python                                              24.0MB @  54.8MB/s  0.4s
khmer                                               10.8MB @  11.9MB/s  0.8s

Downloading and Extracting Packages

Preparing transaction: done
Verifying transaction: done
Executing transaction: done

To activate this environment, use

     $ mamba activate khmer_env

To deactivate an active environment, use

     $ mamba deactivate
```
</details>
<br />

<a id="ii-run-khmer"></a>
### ii. Run [`khmer`](https://khmer.readthedocs.io/en/latest/)
<a id="code-17"></a>
#### Code
<details>
<summary><i>Code: ii. Run `khmer`</i></summary>

```bash
#!/bin/bash

grabnode  # 1, 20, 1, N


#  Define functions ===========================================================
#  ...


#  Initialize variables and arrays ============================================
d_fa="${HOME}/genomes/Saccharomyces_cerevisiae/fasta-processed"
f_fa="S288C_reference_sequence_R64-3-1_20210421.fa.gz"
a_fa="${d_fa}/${f_fa}"

f_rep="S288C_reference_sequence_R64-3-1_20210421.unique-kmers_50.txt"
a_rep="${d_fa}/${f_rep}"

env_khmer="khmer_env"


#  Do the main work ===========================================================
#  Set flag(s)
check_variables=true
check_file_exists=true
check_command=true
run_command=true

#  Check variables
if ${check_variables}; then
    echo "
           d_fa=${d_fa}
           f_fa=${f_fa}
           a_fa=${a_fa}

          f_rep=${f_rep}
          a_rep=${a_rep}

    env_khmer=${env_khmer}
    "
fi

#  Check that input FASTQ file exists
if ${check_file_exists}; then ls -lhaFG "${a_fa}"; fi

#  If not already activated, the activate conda environment
if [[ "${CONDA_DEFAULT_ENV}" != "${env_khmer}" ]]; then
    if [[ ${CONDA_DEFAULT_ENV} != "base" ]]; then
        conda deactivate
    fi

    source activate "${env_khmer}"
fi

#  Estimate the number of unique 50-mers in S. cerevisiae
if ${check_command}; then
    echo "
        if [[ ! -f \"${a_rep}\" ]]; then
            if [[ ! -f "${a_fa%.gz}" ]]; then
                gzip -cd \\
                    \"${a_fa}\" \\
                        > \"${a_fa%.gz}\"
            fi
            
            if [[ -f \"${a_fa%.gz}\" ]]; then
                unique-kmers.py \\
                    -R \"${a_rep}\" \\
                    -k 50 \\
                    \"${a_fa}\"
            fi
        fi
    "
fi

if ${run_command}; then
    if [[ ! -f "${a_bed}" ]]; then
        if [[ ! -f "${a_fa%.gz}" ]]; then
            gzip -cd "${a_fa}" > "${a_fa%.gz}"
        fi

        if [[ -f "${a_fa%.gz}" ]]; then
            unique-kmers.py \
                -R "${a_rep}" \
                -k 50 \
                "${a_fa}"
        fi
    fi
fi
```
</details>
<br />

<a id="printed-2"></a>
#### Printed
<details>
<summary><i>Printed: ii. Run `khmer`</i></summary>

```txt
❯ if ${check_command}; then
>     echo "
>         if [[ ! -f \"${a_rep}\" ]]; then
>             if [[ ! -f "${a_fa%.gz}" ]]; then
>                 gzip -cd \\
>                     \"${a_fa}\" \\
>                         > \"${a_fa%.gz}\"
>             fi
> 
>             if [[ -f \"${a_fa%.gz}\" ]]; then
>                 unique-kmers.py \\
>                     -R \"${a_rep}\" \\
>                     -k 50 \\
>                     \"${a_fa}\"
>             fi
>         fi
>     "
> fi

    if [[ ! -f "/home/kalavatt/genomes/Saccharomyces_cerevisiae/fasta-processed/S288C_reference_sequence_R64-3-1_20210421.unique-kmers_50.txt" ]]; then
        if [[ ! -f /home/kalavatt/genomes/Saccharomyces_cerevisiae/fasta-processed/S288C_reference_sequence_R64-3-1_20210421.fa ]]; then
            gzip -cd \
                "/home/kalavatt/genomes/Saccharomyces_cerevisiae/fasta-processed/S288C_reference_sequence_R64-3-1_20210421.fa.gz" \
                    > "/home/kalavatt/genomes/Saccharomyces_cerevisiae/fasta-processed/S288C_reference_sequence_R64-3-1_20210421.fa"
        fi

        if [[ -f "/home/kalavatt/genomes/Saccharomyces_cerevisiae/fasta-processed/S288C_reference_sequence_R64-3-1_20210421.fa" ]]; then
            unique-kmers.py \
                -R "/home/kalavatt/genomes/Saccharomyces_cerevisiae/fasta-processed/S288C_reference_sequence_R64-3-1_20210421.unique-kmers_50.txt" \
                -k 50 \
                "/home/kalavatt/genomes/Saccharomyces_cerevisiae/fasta-processed/S288C_reference_sequence_R64-3-1_20210421.fa.gz"
        fi
    fi


❯ if ${run_command}; then
>     if [[ ! -f "${a_bed}" ]]; then
>         if [[ ! -f "${a_fa%.gz}" ]]; then
>             gzip -cd "${a_fa}" > "${a_fa%.gz}"
>         fi
> 
>         if [[ -f "${a_fa%.gz}" ]]; then
>             unique-kmers.py \
>                 -R "${a_rep}" \
>                 -k 50 \
>                 "${a_fa}"
>         fi
>     fi
> fi

|| This is the script unique-kmers.py in khmer.
|| You are running khmer version 3.0.0a3
|| You are also using screed version 1.0.4
||
|| If you use this script in a publication, please cite EACH of the following:
||
||   * MR Crusoe et al., 2015. https://doi.org/10.12688/f1000research.6924.1
||   * A. Döring et al. https://doi.org:80/10.1186/1471-2105-9-11
||   * Irber and Brown. https://doi.org/10.1101/056846
||
|| Please see http://khmer.readthedocs.io/en/latest/citations.html for details.

Estimated number of unique 50-mers in /home/kalavatt/genomes/Saccharomyces_cerevisiae/fasta-processed/S288C_reference_sequence_R64-3-1_20210421.fa.gz: 11624332
Total estimated number of unique 50-mers: 11624332
```
</details>
<br />

<a id="b-determine-base-statistics-in-s-cerevisiae-fasta-files"></a>
## b. Determine base statistics in *S. cerevisiae* FASTA files
<a id="i-install-facount"></a>
### i. Install [`faCount`](https://khmer.readthedocs.io/en/latest/)
<a id="code-18"></a>
#### Code
<details>
<summary><i>Code: i. Install `faCount`</i></summary>
<br />

*See mamba installation of "alignment-processing_env" [above](#i-install-sdust-via-minimap) (which includes the package ucsc-facount).*
</details>
<br />

<a id="ii-run-facount"></a>
### ii. Run `faCount`
<a id="code-19"></a>
#### Code
<details>
<summary><i>Code: Run `faCount`</i></summary>

```bash
#!/bin/bash

grabnode  # 1, 20, 1, N


#  Define functions ===========================================================
#  ...


#  Initialize variables and arrays ============================================
d_fa="${HOME}/genomes/Saccharomyces_cerevisiae/fasta-processed"
f_fa="S288C_reference_sequence_R64-3-1_20210421.fa.gz"
a_fa="${d_fa}/${f_fa}"

env_faCount="alignment-processing_env"


#  Do the main work ===========================================================
#  Set flag(s)
check_variables=true
check_file_exists=true
check_command=true
run_command=false

#  Check variables
if ${check_variables}; then
    echo "
           d_fa=${d_fa}
           f_fa=${f_fa}
           a_fa=${a_fa}

    env_faCount=${env_faCount}
    "
fi

#  Check that input FASTQ file exists
if ${check_file_exists}; then ls -lhaFG "${a_fa}"; fi

#  If not already activated, the activate conda environment
if [[ "${CONDA_DEFAULT_ENV}" != "${env_faCount}" ]]; then
    if [[ ${CONDA_DEFAULT_ENV} != "base" ]]; then
        conda deactivate
    fi

    source activate "${env_faCount}"
fi

#  Estimate the number of unique 50-mers in S. cerevisiae
if ${check_command}; then
    echo "faCount -summary \"${a_fa}\""
fi

faCount -summary "${a_fa}"
```
</details>
<br />

<a id="printed-3"></a>
#### Printed
<details>
<summary><i>Printed: Run `faCount`</i></summary>

```txt
❯ if ${check_command}; then
>     echo "faCount -summary \"${a_fa}\""
> fi
faCount -summary "/home/kalavatt/genomes/Saccharomyces_cerevisiae/fasta-processed/S288C_reference_sequence_R64-3-1_20210421.fa.gz"
 

❯ faCount -summary "${a_fa}"
#seq    len A   C   G   T   N   cpg
total   12157105    3766349 2320576 2317100 3753080 0   355299
prcnt   1.0     0.3098  0.1909  0.1906  0.3087  0.0000  0.0292
```
</details>
<br />

