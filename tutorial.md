
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
    1. [b. Concatenate processed fastas and gff3s in new directory `combined_SC_SP/`](#b-concatenate-processed-fastas-and-gff3s-in-new-directory-combined_sc_sp)
        1. [Code](#code-4)
    1. [c. Create `bowtie2` and `bwa` indices for "`combined_SC_SP.fa.gz`"](#c-create-bowtie2-and-bwa-indices-for-combined_sc_spfagz)
        1. [Code](#code-5)
        1. [Notes](#notes)
            1. [Breaking down the call to `mkdir -p`, which makes use of brace expansion](#breaking-down-the-call-to-mkdir--p-which-makes-use-of-brace-expansion)
    1. [d. Use `bowtie2` to align the trimmed FASTQ files](#d-use-bowtie2-to-align-the-trimmed-fastq-files)
        1. [Code](#code-6)
    1. [e. Use `bwa` to align the trimmed FASTQ files](#e-use-bwa-to-align-the-trimmed-fastq-files)
        1. [Code](#code-7)
1. [4. Call peaks with MACS3](#4-call-peaks-with-macs3)
    1. [a. Install MACS3](#a-install-macs3)
        1. [Code](#code-8)
    1. [b. Run MACS3](#b-run-macs3)
        1. [Code](#code-9)
1. [5. Miscellaneous \(to be organized\)](#5-miscellaneous-to-be-organized)
    1. [a. Determine the locations of low-complexity regions in *S. cerevisiae*](#a-determine-the-locations-of-low-complexity-regions-in-s-cerevisiae)
        1. [i. Install sdust via minimap](#i-install-sdust-via-minimap)
            1. [Code](#code-10)
            1. [Printed](#printed)
        1. [ii. Run sdust via minimap](#ii-run-sdust-via-minimap)
            1. [Code](#code-11)
    1. [b. Determine the effective genome size of *S. cerevisiae* \(50-mers\)](#b-determine-the-effective-genome-size-of-s-cerevisiae-50-mers)
        1. [i. Install khmer](#i-install-khmer)
            1. [Code](#code-12)
            1. [Printed](#printed-1)
        1. [ii. Run khmer](#ii-run-khmer)
            1. [Code](#code-13)
            1. [Printed](#printed-2)
    1. [b. Determine base statistics in *S. cerevisiae* FA files](#b-determine-base-statistics-in-s-cerevisiae-fa-files)
        1. [i. Install faCount](#i-install-facount)
            1. [Code](#code-14)
        1. [ii. Run faCount](#ii-run-facount)
            1. [Code](#code-15)
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

<a id="b-concatenate-processed-fastas-and-gff3s-in-new-directory-combined_sc_sp"></a>
## b. Concatenate processed fastas and gff3s in new directory `combined_SC_SP/`
<a id="code-4"></a>
### Code
<details>
<summary><i>Code: 3.b. Concatenate processed fastas and gff3s in new directory `combined_SC_SP/`</i></summary>

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
<br />

<a id="notes"></a>
### Notes
<details>
<summary><i>Notes: 3.a. Generate a concatenated annotated assembly of the *S. cerevisiae* and *S. pombe* genomes</i></summary>

<a id="breaking-down-the-call-to-mkdir--p-which-makes-use-of-brace-expansion"></a>
#### Breaking down the call to `mkdir -p`, which makes use of brace expansion
1. `[[ ! -d "${dir_genomes}" ]]`: This checks if the directory `"${dir_genomes}"` does not exist. `!` is used for negation, `-d` checks for the existence of a directory, and `${dir_genomes}` is a variable that should hold the path of the directory you're checking.
2. `mkdir -p "${HOME}/genomes/"{"${dir_SP}","${dir_SC}"}/{fasta,gff3}`: This command creates multiple directories in one go.
    + `mkdir -p`: The `mkdir` command is used to create directories, and the `-p` flag ensures that any necessary parent directories are also created (and also prevents an error if the directory already exists).
    + `"${dir_genomes}/"{"${dir_SP}","${dir_SC}"}/{fasta,gff3}`: This is an example of brace expansion. The `mkdir -p` command creates directories in the directory assigned to `dir_genomes`. For each of `${dir_SP}` and `${dir_SC}`, `mkdir -p` will create `fasta` and `gff3` subdirectories. For example, if `${dir_SP}` is "Schizosaccharomyces_pombe" and `${dir_SC}` is "Saccharomyces_cerevisiae", the command will create the following directories:
        - `${HOME}/genomes/Schizosaccharomyces_pombe/fasta`
        - `${HOME}/genomes/Schizosaccharomyces_pombe/gff3`
        - `${HOME}/genomes/Saccharomyces_cerevisiae/fasta`
        - `${HOME}/genomes/Saccharomyces_cerevisiae/gff3`

\#TODO To be continued after modularizing the above
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


# #  Function to check if a SLURM module is loaded or not; if not, the module is
# #+ loaded
# function check_and_load_module() {
#     local module_name="${1}"
#     if ! module is-loaded "${module_name}" &> /dev/null; then
#         echo "Loading module: ${module_name}"
#         module load "${module_name}"
#     else
#         echo "Module already loaded: ${module_name}"
#     fi
# }


#  Initialize variables and arrays ============================================
dir_base="${HOME}/tsukiyamalab"                          # Base directory for lab data
dir_repo="Kris/2023_rDNA"                                # Repository directory
dir_work="results/2023-0406_tutorial_ChIP-seq_analyses"  # Work directory
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

# mod_bowtie2="Bowtie2/2.4.4-GCC-11.2.0"
# mod_samtools="SAMtools/1.16.1-GCC-11.2.0"
# mod_bedtools="BEDTools/2.30.0-GCC-11.2.0"

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
    "${dir_untr}/in_Q_untagged_5781"
    "${dir_untr}/IP_Q_untagged_5781"
    "${dir_untr}/in_Q_Esa5_7041"
    "${dir_untr}/IP_Q_Esa5_7041"
    "${dir_untr}/in_Q_Rpd3_7568"
    "${dir_untr}/IP_Q_Rpd3_7568"
    "${dir_untr}/in_Q_Rpd3_7569"
    "${dir_untr}/IP_Q_Rpd3_7569"
    "${dir_untr}/in_Q_Esa5_7691"
    "${dir_untr}/IP_Q_Esa5_7691"
    "${dir_untr}/in_Q_Gcn5_7692"
    "${dir_untr}/IP_Q_Gcn5_7692"
    "${dir_untr}/in_Q_Gcn5_7709"
    "${dir_untr}/IP_Q_Gcn5_7709"
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

    # mod_bowtie2=${mod_bowtie2}
    # mod_samtools=${mod_samtools}
    # mod_bedtools=${mod_bedtools}
    "
fi

#  Echo array contents if check_array is true
if ${check_array}; then
    for i in "${!file_fastqs[@]}"; do
        file="${file_fastqs[${i}]}"

        if [[ "$(dirname ${file})" == "${dir_trim}" ]]; then
            echo "
            read #1 ...... ${file}_R1.atria.fastq.gz
            read #2 ...... ${file}_R2.atria.fastq.gz
            "
        elif [[ "$(dirname ${file})" == "${dir_untr}" ]]; then
            echo "
            read #1 ...... ${file}_R1.fastq.gz
            read #2 ...... ${file}_R2.fastq.gz
            "
        else
            error_and_return "Processing logic problem for ${file}."
        fi
    done
fi

# #  If the module for Bowtie2 is not loaded, then load it
# check_and_load_module "${mod_bowtie2}"
#
# #  If the corresponding Samtools module is not loaded, then load it
# check_and_load_module "${mod_samtools}"
#
# #  Load also BEDtools
# check_and_load_module "${mod_bedtools}"

#  Initialize conda/mamba environment containing necessary programs for
#+ alignment, quality checks, and post-processing
env_name="alignment-processing_env"

#TODO Make the following into a function
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

#  Navigate to the work directory
cd "${dir_base}/${dir_repo}/${dir_work}" \
    || error_and_return "Failed to cd to ${dir_base}/${dir_repo}/${dir_work}."

#  If it doesn't exist, create a directory to store Bowtie2-aligned BAM files
if [[ ! -d "${dir_bwt2}" ]]; then
    mkdir -p "${dir_bwt2}/"{err_out,bam,siQ-ChIP,cvrg,qc}
fi

#  Set flags: checking variables, checking and submitting Bowtie2 jobs
print_iteration=true
check_variables=true
check_operation=true
run_operation=true
# check_operations=false
# run_operations=false

for i in "${!file_fastqs[@]}"; do
    # i=0
    index="${i}"
    iter=$(( index + 1 ))
    file="${file_fastqs[${index}]}"
    stem="$(basename ${file})"
    job_name="$(echo ${dir_bwt2} | sed 's:\/:_:g').${stem}"
    
    #  Parse the files' source directory to determine appropriate suffix for
    #+ FASTQ file names
    if [[ "$(dirname ${file})" == "${dir_trim}" ]]; then
        fastq_1=${file}_R1.atria.fastq.gz
        fastq_2=${file}_R2.atria.fastq.gz
    elif [[ "$(dirname ${file})" == "${dir_untr}" ]]; then
        fastq_1=${file}_R1.fastq.gz
        fastq_2=${file}_R2.fastq.gz
    else
        error_and_return "Processing logic problem for ${file}."
    fi
    
    bam="${dir_bwt2}/bam/${stem}.bam"
    bam_coor="${bam/.bam/.sort-coord.bam}"
    bam_quer="${bam/.bam/.sort-qname.bam}"
    
    bed_siQ="${dir_bwt2}/siQ-ChIP/${stem}.bed.gz"
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

    if [[ ${check_operation} ]]; then
        echo "
sbatch << EOF
#!/bin/bash

#SBATCH --job-name=\"${job_name}\"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=${threads}
#SBATCH --time=${time}
#SBATCH --error=\"${dir_bwt2}/err_out/${job_name}.%A.stderr.txt\"
#SBATCH --output=\"${dir_bwt2}/err_out/${job_name}.%A.stdout.txt\"

bash align-process-etc_fastqs_bowtie2.sh \\
    --threads \"${threads}\" \\
    --index \"${dir_indx}\" \\
    --fasta \"${file_fasta}\" \\
    --sizes \"${file_sizes}\" \\
    --mapq \"${mapq}\" \\
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

    if [[ ${run_operation} ]]; then
sbatch << EOF
#!/bin/bash

#SBATCH --job-name="${job_name}"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=${threads}
#SBATCH --time=${time}
#SBATCH --error="${dir_bwt2}/err_out/${job_name}.%A.stderr.txt"
#SBATCH --output="${dir_bwt2}/err_out/${job_name}.%A.stdout.txt"

bash align-process-etc_fastqs_bowtie2.sh \
    --threads "${threads}" \
    --index "${dir_indx}" \
    --fasta "${file_fasta}" \
    --sizes "${file_sizes}" \
    --mapq "${mapq}" \
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

# for i in "${!file_fastqs[@]}"; do
#     # i=0
#     index="${i}"
#     iter=$(( index + 1 ))
#     file="${file_fastqs[${index}]}"
#     stem="$(basename ${file})"
#     job_name="$(echo ${dir_bwt2} | sed 's:\/:_:g').${stem}"
#    
#     #  Parse the files' source directory to determine appropriate suffix for
#     #+ FASTQ files
#     if [[ "$(dirname ${file})" == "${dir_trim}" ]]; then
#         fastq_1=${file}_R1.atria.fastq.gz
#         fastq_2=${file}_R2.atria.fastq.gz
#     elif [[ "$(dirname ${file})" == "${dir_untr}" ]]; then
#         fastq_1=${file}_R1.fastq.gz
#         fastq_2=${file}_R2.fastq.gz
#     else
#         error_and_return "Processing logic problem for ${file}."
#     fi
#    
#     bam="${dir_bwt2}/bam/${stem}.bam"
#     bam_coor="${bam/.bam/.sort-coord.bam}"
#     bam_quer="${bam/.bam/.sort-qname.bam}"
#    
#     bed_siQ="${dir_bwt2}/siQ-ChIP/${stem}.bed.gz"
#     bed_etc="${dir_bwt2}/cvrg/${stem}"
#    
#     txt_flg="${dir_bwt2}/qc/${stem}.samtools-flagstat.txt"
#     txt_idx="${dir_bwt2}/qc/${stem}.samtools-idxstats.txt"
#
#     #  Echo current iteration
#     if ${print_iteration}; then
#         echo "
#         #  -------------------------------------
#         ### ${iter} ###
#         "
#     fi
#    
#     #  Echo loop-dependent variables if check_variables is true
#     if ${check_variables}; then
#         echo "
#         index=${index}
#         iter=${iter}
#         file=${file}
#         stem=${stem}
#         job_name=${job_name}
#        
#         fastq_1=${fastq_1}
#         fastq_2=${fastq_2}
#
#         bam=${bam}
#         bam_coor=${bam_coor}
#         bam_quer=${bam_quer}
#         bed_siQ=${bed_siQ}
#         bed_etc=${bed_etc}
#
#         txt_flg=${txt_flg}
#         txt_idx=${txt_idx}
#         "
#     fi
#
#     #  Echo the Atria trimming command if check_operations is true
#     if ${check_operations}; then
#         echo "
#         if [[
#                  -f \"${fastq_1}\" \\
#             &&   -f \"${fastq_2}\" \\
#             && ! -f \"${bam_coor}\" \\
#             && ! -f \"${bam_quer}\" \\
#             && ! -f \"${bed}\"
#         ]]; then
# sbatch << EOF
# #!/bin/bash
#
# #SBATCH --job-name=\"${job_name}\"
# #SBATCH --nodes=1
# #SBATCH --cpus-per-task=${threads}
# #SBATCH --time=${time}
# #SBATCH --error=\"${dir_bwt2}/err_out/${job_name}.%A.stderr.txt\"
# #SBATCH --output=\"${dir_bwt2}/err_out/${job_name}.%A.stdout.txt\"
#
# #  Check if the BAM file exists; if not, perform alignment with Bowtie 2,
# #+ converting the Bowtie 2 output to a BAM file with Samtools; in doing so,
# #+ retain only properly paired reads (-f 2) that are high-quality multi-reads
# #+ or any-quality maxi-reads (-q 1)
# #+ 
# #+ On what multi- and maxi-reads are, and how to interpret Bowtie 2 MAPQ
# #+ scores:
# #+ biofinysics.blogspot.com/2014/05/how-does-bowtie2-assign-mapq-scores.html
# if [[ ! -f \"${bam}\" ]]; then
#     bowtie2 \\
#         -p ${threads} \\
#         -x \"${dir_indx}\" \\
#         --very-sensitive-local \\
#         --no-unal \\
#         --no-mixed \\
#         --no-discordant \\
#         --no-overlap \\
#         --no-dovetail \\
#         --phred33 \\
#         -I 10 \\
#         -X 700 \\
#         -1 \"${fastq_1}\" \\
#         -2 \"${fastq_2}\" \\
#             | samtools view \\
#                 -@ ${threads} \\
#                 -b \\
#                 -f 2 \\
#                 -q 1 \\
#                 -o \"${bam}\"
# fi
#
# #  Check if the BAM file exists to perform further operations
# if [[ -f \"${bam}\" ]]; then
#     #  Sort the BAM file by queryname if not already done
#     if [[ ! -f \"${bam_quer}\" ]]; then
#         samtools sort \\
#             -@ ${threads} \\
#             -n \\
#             -o \"${bam_quer}\" \\
#             \"${bam}\"
#     fi
#
#     #  Fix the paired read mate information, which is required after sorting by
#     #+ queryname for subsequent operations
#     if [[ -f \"${bam_quer}\" ]]; then
#         samtools fixmate \\
#             -@ ${threads} \\
#             -m \\
#             \"${bam_quer}\" \\
#             \"${bam_quer%.bam}.tmp.bam\"
#
#         #  Replace the original queryname-sorted BAM with queryname-sorted
#         #+ mate-fixed BAM
#         if [[ -f \"${bam_quer%.bam}.tmp.bam\" ]]; then
#             mv -f \\
#                 \"${bam_quer%.bam}.tmp.bam\" \\
#                 \"${bam_quer}\"
#         fi
#     fi
#
#     #  For downstream analyses, sort the queryname-sorted BAM by coordinates
#     if [[ ! -f \"${bam_coor}\" ]]; then
#         samtools sort \\
#             -@ ${threads} \\
#             -o \"${bam_coor}\" \\
#             \"${bam_quer}\"
#     fi
#
#     #  Index the coordinate-sorted BAM file
#     if [[ ! -f \"${bam_coor}.bai\" ]]; then
#         samtools index \\
#             -@ ${threads} \\
#             \"${bam_coor}\"
#     fi
#
#     #  Mark duplicate alignments in the coordinate-sorted BAM file
#     if [[
#            -f \"${bam_coor}\" \\
#         && -f \"${bam_coor}.bai\"
#     ]]; then
#         samtools markdup \\
#             -@ ${threads} \\
#             \"${bam_coor}\" \\
#             \"${bam_coor%.bam}.tmp.bam\"
#
#         #  Replace the original coordinate-sorted BAM with one in which
#         #+ duplicates alignments are marked
#         if [[ -f \"${bam_coor%.bam}.tmp.bam\" ]]; then
#             mv -f \\
#                 \"${bam_coor%.bam}.tmp.bam\" \\
#                 \"${bam_coor}\"
#         fi
#
#         #  If duplicate marking was successful, then generate flagstat and
#         #+ idxstats reports
#         if [[ $? -eq 0 ]]; then
#             samtools flagstat \\
#                 -@ ${threads} \\
#                 \"${bam_coor}\" \\
#                     > \"${txt_flg}\"
#
#             samtools idxstats \\
#                 \"${bam_coor}\" \\
#                     > \"${txt_idx}\"
#         fi
#     fi
#
#     #  Generate a BED file from the queryname-sorted BAM if it doesn't exist
#     if [[
#          ! -f \"${bed}\" \\
#         && -f \"${bam_quer}\"
#     ]]; then
#         #  Extract fragment information to create the BED file, excluding
#         #+ chromosomes starting with \"SP_\"
#         samtools view \"${bam_quer}\" \\
#             | awk '{
#                 if (NR % 2 == 1) {
#                     chr_1 = \$3; 
#                     start_1 = \$4; 
#                     len_1 = length(\$10);
#                 } else {
#                     chr_2 = \$3;
#                     start_2 = \$4; 
#                     len_2 = length(\$10);
#                     if (chr_1 == chr_2 && substr(chr_1, 1, 3) != \"SP_\") {
#                         start = (start_1 < start_2) ? start_1 : start_2;
#                         end = (start_1 < start_2) ? start_2 + len_2 - 1 : start_1 + len_1 - 1;
#                         frag_length = end - start + 1; 
#                         print chr_1, start, end, frag_length;
#                     }
#                 }
#             }' OFS='\t' \\
#             | sort -k1,1 -k2,2n \\
#             | gzip \\
#                 > \"${bed}\"
#     fi
#
#     #  Remove the original BAM file (output by Bowtie 2 piped to Samtools) if
#     #+ all other files have been successfully created
#     if [[ 
#            -f \"${bam_coor}\" \\
#         && -f \"${bam_quer}\" \\
#         && -f \"${bed}\"
#     ]]; then
#         rm \"${bam}\"
#     fi
# fi
# EOF
#         else
#             echo \"
#             Warning: ${stem} FASTQs do not appear to exist; skipping alignment and processing.
#             \"
#         fi
#         "
#     fi
#
#     #TODO Pick up here #TOMORROW: Write up/incluce code teaching how to
#     #+    generate Bowtie2 indices, and add variable above for the indices
#     #+    and the scratch directory, where sorting will take place
#
#     #  Submit the Atria trimming job if run_operations is true
#     if ${run_operations}; then
#         if [[
#                  -f "${fastq_1}" \
#             &&   -f "${fastq_2}" \
#             && ! -f "${bam_coor}" \
#             && ! -f "${bam_quer}" \
#             && ! -f "${bed}"
#         ]]; then
# sbatch << EOF
# #!/bin/bash
#
# #SBATCH --job-name="${job_name}"
# #SBATCH --nodes=1
# #SBATCH --cpus-per-task=${threads}
# #SBATCH --time=${time}
# #SBATCH --error="${dir_bwt2}/err_out/${job_name}.%A.stderr.txt"
# #SBATCH --output="${dir_bwt2}/err_out/${job_name}.%A.stdout.txt"
#
# #  Check if the BAM file exists; if not, perform alignment with Bowtie 2,
# #+ converting the Bowtie 2 output to a BAM file with Samtools; in doing so,
# #+ retain only properly paired reads (-f 2) that are high-quality multi-reads
# #+ or any-quality maxi-reads (-q 1)
# #+ 
# #+ On what multi- and maxi-reads are, and how to interpret Bowtie 2 MAPQ
# #+ scores:
# #+ biofinysics.blogspot.com/2014/05/how-does-bowtie2-assign-mapq-scores.html
# if [[ ! -f "${bam}" ]]; then
#     bowtie2 \
#         -p ${threads} \
#         -x "${dir_indx}" \
#         --very-sensitive-local \
#         --no-unal \
#         --no-mixed \
#         --no-discordant \
#         --no-overlap \
#         --no-dovetail \
#         --phred33 \
#         -I 10 \
#         -X 700 \
#         -1 "${fastq_1}" \
#         -2 "${fastq_2}" \
#             | samtools view \
#                 -@ ${threads} \
#                 -b \
#                 -f 2 \
#                 -q 1 \
#                 -o "${bam}"
# fi
#
# #  Check if the BAM file exists to perform further operations
# if [[ -f "${bam}" ]]; then
#     #  Sort the BAM file by queryname if not already done
#     if [[ ! -f "${bam_quer}" ]]; then
#         samtools sort \
#             -@ ${threads} \
#             -n \
#             -o "${bam_quer}" \
#             "${bam}"
#     fi
#
#     #  Fix the paired read mate information, which is required after sorting by
#     #+ queryname for subsequent operations
#     if [[ -f "${bam_quer}" ]]; then
#         samtools fixmate \
#             -@ ${threads} \
#             -m \
#             "${bam_quer}" \
#             "${bam_quer%.bam}.tmp.bam"
#
#         #  Replace the original queryname-sorted BAM with queryname-sorted
#         #+ mate-fixed BAM
#         if [[ -f "${bam_quer%.bam}.tmp.bam" ]]; then
#             mv -f \
#                 "${bam_quer%.bam}.tmp.bam" \
#                 "${bam_quer}"
#         fi
#     fi
#
#     #  For downstream analyses, sort the queryname-sorted BAM by coordinates
#     if [[ ! -f "${bam_coor}" ]]; then
#         samtools sort \
#             -@ ${threads} \
#             -o "${bam_coor}" \
#             "${bam_quer}"
#     fi
#
#     #  Index the coordinate-sorted BAM file
#     if [[ ! -f "${bam_coor}.bai" ]]; then
#         samtools index \
#             -@ ${threads} \
#             "${bam_coor}"
#     fi
#
#     #  Mark duplicate alignments in the coordinate-sorted BAM file
#     if [[
#            -f "${bam_coor}" \
#         && -f "${bam_coor}.bai"
#     ]]; then
#         samtools markdup \
#             -@ ${threads} \
#             "${bam_coor}" \
#             "${bam_coor%.bam}.tmp.bam"
#
#         #  Replace the original coordinate-sorted BAM with one in which
#         #+ duplicates alignments are marked
#         if [[ -f "${bam_coor%.bam}.tmp.bam" ]]; then
#             mv -f \
#                 "${bam_coor%.bam}.tmp.bam" \
#                 "${bam_coor}"
#         fi
#
#         #  If duplicate marking was successful, then generate flagstat and
#         #+ idxstats reports
#         if [[ $? -eq 0 ]]; then
#             samtools flagstat \
#                 -@ ${threads} \
#                 "${bam_coor}" \
#                     > "${txt_flg}"
#
#             samtools idxstats \
#                 "${bam_coor}" \
#                     > "${txt_idx}"
#         fi
#     fi
#
#     #  Generate a BED file from the queryname-sorted BAM if it doesn't exist
#     if [[
#          ! -f "${bed}" \
#         && -f "${bam_quer}"
#     ]]; then
#         #  Extract fragment information to create the BED file, excluding
#         #+ chromosomes starting with "SP_"
#         samtools view "${bam_quer}" \
#             | awk '{
#                 if (NR % 2 == 1) {
#                     chr_1 = \$3; 
#                     start_1 = \$4; 
#                     len_1 = length(\$10);
#                 } else {
#                     chr_2 = \$3;
#                     start_2 = \$4; 
#                     len_2 = length(\$10);
#                     if (chr_1 == chr_2 && substr(chr_1, 1, 3) != "SP_") {
#                         start = (start_1 < start_2) ? start_1 : start_2;
#                         end = (start_1 < start_2) ? start_2 + len_2 - 1 : start_1 + len_1 - 1;
#                         frag_length = end - start + 1; 
#                         print chr_1, start, end, frag_length;
#                     }
#                 }
#             }' OFS='\t' \
#             | sort -k1,1 -k2,2n \
#             | gzip \
#                 > "${bed}"
#     fi
#
#     #  Remove the original BAM file (output by Bowtie 2 piped to Samtools) if
#     #+ all other files have been successfully created
#     if [[ 
#            -f "${bam_coor}" \
#         && -f "${bam_quer}" \
#         && -f "${bed}"
#     ]]; then
#         rm "${bam}"
#     fi
# fi
# EOF
#         else
#             echo "
#             Warning: ${stem} FASTQs do not appear to exist; skipping alignment and processing.
#             "
#         fi
#     fi
#
#     sleep 0.2  # Short pause to prevent rapid job-submission overload
# done
```
</details>
<br />

<a id="e-use-bwa-to-align-the-trimmed-fastq-files"></a>
## e. Use `bwa` to align the trimmed FASTQ files
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

<a id="4-call-peaks-with-macs3"></a>
# 4. Call peaks with MACS3
<a id="a-install-macs3"></a>
## a. Install MACS3
<a id="code-8"></a>
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
            -n "macs3_env" \
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
<a id="code-9"></a>
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

    if [[ ${check_operation} ]]; then
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

    if [[ ${run_operation} ]]; then
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
<br />

<a id="5-miscellaneous-to-be-organized"></a>
# 5. Miscellaneous (to be organized)
<a id="a-determine-the-locations-of-low-complexity-regions-in-s-cerevisiae"></a>
## a. Determine the locations of low-complexity regions in *S. cerevisiae*
<a id="i-install-sdust-via-minimap"></a>
### i. Install [sdust](https://pubmed.ncbi.nlm.nih.gov/16796549/) via [minimap](https://github.com/lh3/minimap/tree/master)
<a id="code-10"></a>
#### Code
<details>
<summary><i>Code: i. Install sdust via minimap</i></summary>

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
<summary><i>Printed: i. Install sdust via minimap</i></summary>

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
### ii. Run [sdust](https://pubmed.ncbi.nlm.nih.gov/16796549/) via [minimap](https://github.com/lh3/minimap/tree/master)
<a id="code-11"></a>
#### Code
<details>
<summary><i>Code: ii. Run sdust via minimap</i></summary>

```bash
#!/bin/bash

grabnode  # 1, 20, 1, N


#  Define functions ===========================================================


#  Initialize variables and arrays ============================================
d_fa="${HOME}/genomes/Saccharomyces_cerevisiae/fasta-processed"
f_fa="S288C_reference_sequence_R64-3-1_20210421.fa.gz"
a_fa="${d_fa}/${f_fa}"

f_bed="S288C_reference_sequence_R64-3-1_20210421.low_complexity.bed"
a_bed="${d_fa}/${f_bed}"

env_minimap="minimap_env"


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

    env_minimap=${env_minimap}
    "
fi

#  Check that input FASTQ file exists
if ${check_file_exists}; then ls -lhaFG "${a_fa}"; fi

#  If not already activated, the activate conda environment
if [[ "${CONDA_DEFAULT_ENV}" != "${env_minimap}" ]]; then
    if [[ ${CONDA_DEFAULT_ENV} != "base" ]]; then
        conda deactivate
    fi

    source activate "${env_minimap}"
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
### i. Install [khmer](https://khmer.readthedocs.io/en/latest/)
<a id="code-12"></a>
#### Code
<details>
<summary><i>Code: i. Install khmer</i></summary>

```bash
#!/bin/bash

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
<summary><i>Printed: ii. Run khmer</i></summary>

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
### ii. Run [khmer](https://khmer.readthedocs.io/en/latest/)
<a id="code-13"></a>
#### Code
<details>
<summary><i>Code: ii. Install khmer</i></summary>

```bash
#!/bin/bash

grabnode  # 1, 20, 1, N


#  Define functions ===========================================================


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
<summary><i>Printed: ii. Install khmer</i></summary>

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

<a id="b-determine-base-statistics-in-s-cerevisiae-fa-files"></a>
## b. Determine base statistics in *S. cerevisiae* FA files
<a id="i-install-facount"></a>
### i. Install [faCount](https://khmer.readthedocs.io/en/latest/)
<a id="code-14"></a>
#### Code
<details>
<summary><i>Code: i. Install faCount</i></summary>
<br />

*See mamba installation of "alignment-processing_env" [above](#i-install-sdust-via-minimap) (which includes the package ucsc-facount).*
</details>
<br />

<a id="ii-run-facount"></a>
### ii. Run faCount
<a id="code-15"></a>
#### Code
<details>
<summary><i>Code: Run faCount</i></summary>

```bash
#!/bin/bash

grabnode  # 1, 20, 1, N


#  Define functions ===========================================================


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
<summary><i>Printed: Run faCount</i></summary>

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

