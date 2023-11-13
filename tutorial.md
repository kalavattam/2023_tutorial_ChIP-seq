
`#tutorial.md`
<br />

<details>
<summary><b><font size="+2"><i>Table of contents</i></font></b></summary>
<br />
<!-- MarkdownTOC -->

1. [1. Prepare FASTQ files of interest for processing and analyses](#1-prepare-fastq-files-of-interest-for-processing-and-analyses)
    1. [Code](#code)
    1. [Notes](#notes)
        1. [Comments \(`#`\)](#comments-)
        1. [Defining a function such as `error_and_return`](#defining-a-function-such-as-error_and_return)
        1. [The body of the function `error_and_return`](#the-body-of-the-function-error_and_return)
        1. [Variables](#variables)
        1. [The `${HOME}` variable](#the-%24home-variable)
        1. [Wrapping variables in curly braces `{}` and double quotes `""`](#wrapping-variables-in-curly-braces--and-double-quotes-)
        1. [Associative arrays \(hash maps\)](#associative-arrays-hash-maps)
        1. [Logical commands \(`true`, `false`\) assigned to "flag variables" \(flags\)](#logical-commands-true-false-assigned-to-flag-variables-flags)
        1. [`if` statements](#if-statements)
        1. [Conditional checks for directories \(`-d`\), files \(`-f`\), and logical negation \(`!`\)](#conditional-checks-for-directories--d-files--f-and-logical-negation-)
        1. [Calls to `ln`](#calls-to-ln)
        1. [Calls to `ls`](#calls-to-ls)
1. [2. Adapter- and quality-trim the FASTQ files](#2-adapter--and-quality-trim-the-fastq-files)
    1. [Install Atria and dependencies](#install-atria-and-dependencies)
        1. [Code](#code-1)
        1. [Notes](#notes-1)
    1. [Adapter- and quality-trim the FASTQ files using Atria](#adapter--and-quality-trim-the-fastq-files-using-atria)
        1. [Code](#code-2)
        1. [Notes](#notes-2)

<!-- /MarkdownTOC -->
</details>
<br />

<a id="1-prepare-fastq-files-of-interest-for-processing-and-analyses"></a>
## 1. Prepare FASTQ files of interest for processing and analyses
<a id="code"></a>
### Code
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
dir_orig="Rina/ChIP-seq/230915_hho1_hmo1_rhirano"        # Directory where original ChIP-seq data is stored
dir_sym="01_sym"                                         # Directory name where symbolic links will be stored

#  Create an associative array/hash map for renaming files through symbolic
#+ links. Original file stems are keys, and new file stems are values. This
#+ is naming scheme for symlink-renamed files: assay_state_factor_strain
unset file_fastqs && typeset -A file_fastqs=(
    ["6336_G1_in_S15"]="in_G1_Hho1_6336"
    ["6336_G1_lP_S27"]="IP_G1_Hho1_6336"
    ["6336_G2M_in_S17"]="in_G2M_Hho1_6336"
    ["6336_G2M_lP_S29"]="IP_G2M_Hho1_6336"
    ["6336_Q-in_S19"]="in_Q_Hho1_6336"
    ["6336_Q_lP_S31"]="IP_Q_Hho1_6336"
    ["6337_G1_in_S16"]="in_G1_Hho1_6337"
    ["6337_GI_IP_S28"]="IP_G1_Hho1_6337"
    ["6337_G2M_in_S18"]="in_G2M_Hho1_6337"
    ["6337_G2M_lP_S30"]="IP_G2M_Hho1_6337"
    ["6337_Q_in_S20"]="in_Q_Hho1_6337"
    ["6337_Q_lP_S32"]="IP_Q_Hho1_6337"
    ["7750_G1_in_S21"]="in_G1_Hmo1_7750"
    ["7750_G1_lP_S33"]="IP_G1_Hmo1_7750"
    ["7750_G2M_in_S23"]="in_G2M_Hmo1_7750"
    ["7750_G2M_lP_S35"]="IP_G2M_Hmo1_7750"
    ["7750_Q_ln_S25"]="in_Q_Hmo1_7750"
    ["7750_Q_lP_S37"]="IP_Q_Hmo1_7750"
    ["7751_G1_in_S22"]="in_G1_Hmo1_7751"
    ["7751_G1_lP_S34"]="IP_G1_Hmo1_7751"
    ["7751_G2M_in_S24"]="in_G2M_Hmo1_7751"
    ["7751_G2M_lP_S36"]="IP_G2M_Hmo1_7751"
    ["7751_Q_ln_S26"]="in_Q_Hmo1_7751"
    ["7751_Q_lP_S38"]="IP_Q_Hmo1_7751"
)


#  Do main work ===============================================================
#  Set flags to check variable and array assignments
check_variables=false
check_array=false

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
    || error_and_return "Failed to cd to ${dir_base}/${dir_repo}/${dir_work}"

#  If it doesn't exist, then create a directory to store symlinked FASTQ files
if [[ ! -d "${dir_sym}" ]]; then
    mkdir -p "${dir_sym}"
fi

#  Set flags for checking and running symlinking operations
check_operations=true
run_operations=false

#  Loop through each entry in the associative array to create symbolic links
for i in "${!file_fastqs[@]}"; do
    key="${i}"
    value="${file_fastqs["${key}"]}"
    
    #  If check_operations is true, then echo the operations that will be run
    if ${check_operations}; then
        echo "
        #  Check if the original file for read #1 exists before creating a symlink
        if [[ -f "${dir_base}/${dir_orig}/${key}_R1_001.fastq.gz" ]]; then
            ln -s \\
                \"${dir_base}/${dir_orig}/${key}_R1_001.fastq.gz\" \\
                \"${dir_sym}/${value}_R1.fastq.gz\"
        fi

        #  Check if the original file for read #2 exists before creating a symlink
        if [[ -f "${dir_base}/${dir_orig}/${key}_R2_001.fastq.gz" ]]; then
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

<a id="notes"></a>
### Notes
<details>
<summary><i>Notes: 1. Prepare FASTQ files of interest for processing and analyses</i></summary>
<br />

<a id="comments-"></a>
#### Comments (`#`)
The lines starting with `#` are comments. Comments are used to explain what the code does, but they are not executed. They're there to help anyone reading the code understand it better.

For example, `# Define functions ===========================================================` is a header comment, letting readers know that functions are defined below. The next comment explains what the function `error_and_return` does.

<a id="defining-a-function-such-as-error_and_return"></a>
#### Defining a function such as `error_and_return`
`function error_and_return()` starts the definition of a function named `error_and_return`. A function is a block of code that performs a specific task. Functions facilitate the reuse of the same code multiple times without having to write it out each time.

<a id="the-body-of-the-function-error_and_return"></a>
#### The body of the function `error_and_return`
Inside the function, we have the following code:
- `echo "Error: ${1}" >&2`
    + `echo` is a command used to display text.
    + `"Error: ${1}"` is the text to be displayed. Here, `${1}` is a placeholder for the first argument passed to the function. This means that whatever text is provided when calling `error_and_return` will be displayed after `"Error: "`.
    + `>&2` means that the error message is redirected to the standard error stream (stderr). This is typically used to output error messages.
- `return 1`
    + `return 1` exits the function and returns a value of 1. In Unix, Linux, and MacOS, returning a non-zero value generally indicates an error or abnormal condition. So, when `error_and_return` is called, it indicates that something went wrong.

<a id="variables"></a>
#### Variables
Variables (e.g., `dir_base`) are used to store and retrieve data (e.g., `dir_base="${HOME}/tsukiyamalab"`). They are like placeholders for values that can change over time. Variables make scripts flexible and reusable. For example, changing `dir_base` (to, say, `dir_base="${HOME}/Desktop"`) will update the base directory used throughout the script.

<a id="the-%24home-variable"></a>
#### The `${HOME}` variable
The `${HOME}` variable is an environmental variable that refers to the home directory of the current user. In Unix and Unix-like operating systems (e.g., MacOS), every user is assigned a unique directory where they can store personal files, configurations, and scripts. This directory is commonly referred to as the "home directory."
- In scripts and command lines, `${HOME}` is often used as a shorthand to access the user's home directory. For example, `cd ${HOME}` would change the current directory to the user's home directory.
- If the username is `john`, the home directory might be `/home/john` on Linux or `/Users/john` on macOS. In this case, `${HOME}` would be equivalent to `/home/john` or `/Users/john`, respectively.
- Using `${HOME}` makes scripts more portable and user-independent, as it automatically resolves to the home directory of the user running the script, without needing to hard-code the full path.

<a id="wrapping-variables-in-curly-braces--and-double-quotes-"></a>
#### Wrapping variables in curly braces `{}` and double quotes `""`
In shell scripting, it's considered good practice to wrap variables in curly braces (`{}`) and double quotes (`""`) for clarity and to prevent potential errors. This practice has several benefits:
- Curly braces help in clearly defining the boundary of a variable name. This is particularly useful when a variable is followed by text that could be misinterpreted as part of the variable name. For example, `${var}_suffix` clearly separates `var` from `_suffix`.
- Curly braces facilitate the concatenation of variables with strings. For instance, `${var}value` appends `value` to the contents of `var`.
- Enclosing variables in double quotes prevents word splitting and globbing. Word splitting can lead to unexpected behavior when variables contain spaces or newlines. For example, in `echo ${var}`, if `var` contains `file1 file2`, it will be split into two arguments. `echo "${var}"` would treat it as a single argument, preserving the intended behavior.
- When a variable is empty or undefined, using double quotes ensures that the script doesn't break or behave unpredictably. For example, `echo "${nonexistent}"` will safely print nothing, whereas `echo ${nonexistent}` might lead to unintended script behavior.
- Example:
    + With braces: `echo "${user}file"` makes it clear that `user` is the variable being referenced, and the string `"file"` is being appended to its contents.
    + Without braces: `echo "$userfile"` is misread as if there were a variable named `userfile`.

The use of curly braces and double quotes enhances the readability and reliability of shell scripts, making them less prone to errors.


<a id="associative-arrays-hash-maps"></a>
#### Associative arrays (hash maps)
Associative arrays (or hash maps) are collections of key-value pairs where each key is unique. They allow for more complex data structures, enabling you to map a unique key to a specific value. In the above chunk, keys are original file name stems and values are the new name stems for the symbolic links. For example, `file_fastqs["6336_G1_in_S15"]="in_G1_Hho1_6336"` maps an original file name stem, `"6336_G1_in_S15"`, to a new one, `"in_G1_Hho1_6336"`.

<a id="logical-commands-true-false-assigned-to-flag-variables-flags"></a>
#### Logical commands (`true`, `false`) assigned to "flag variables" (flags)
Flags are variables used to control the flow of the script. `check_variables`, `check_array`, `check_operations` and `run_operations` are flags. When set to `true`, they trigger specific operations like echoing commands or creating symbolic links. Conversely, setting them to `false` skips these operations.

<a id="if-statements"></a>
#### `if` statements
`if` statements are used for the conditional execution of code. They allow the script to make decisions and execute different blocks of code based on whether a condition is true or false. They follow the following logic: "`if` a given condition is `true`, execute a specific code block; `else` (optionally) execute a different code block or do nothing." For example, `if [[ -d "${directory}" ]]; then echo "Directory exists"; fi`.

<a id="conditional-checks-for-directories--d-files--f-and-logical-negation-"></a>
#### Conditional checks for directories (`-d`), files (`-f`), and logical negation (`!`)
- `[[ -d ... ]]`: Checks if a directory exists.
- `[[ -f ... ]]`: Checks if a file exists.
- `!`: Negates a condition, e.g., `[[ ! -d ... ]]` checks if a directory *does not* exist.

<a id="calls-to-ln"></a>
#### Calls to `ln`
The `ln` command creates symbolic links (symlinks), which are pointers to original files. The `-s` option creates a symlink. The command format is `ln -s original_file symlink_file`. Using `ln` like this maintains the integrity of raw data while simplifying access under new names.

<a id="calls-to-ls"></a>
#### Calls to `ls`
The `ls` command lists directory contents. Flags used include:
- `-l`: Long format with detailed information.
- `-h`: Human-readable file sizes.
- `-a`: Includes hidden files.
- `-F`: Appends a character indicating file type.
- `-G`: Colorizes output.

For more details on `ls` flags and shell commands in general, visit [ShellCheck](https://www.shellcheck.net/).
</details>
<br />
<br />

<a id="2-adapter--and-quality-trim-the-fastq-files"></a>
## 2. Adapter- and quality-trim the FASTQ files
<a id="install-atria-and-dependencies"></a>
### Install Atria and dependencies
<a id="code-1"></a>
#### Code
<details>
<summary><i>Code: Install Atria and dependencies</i></summary>

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
    if ! type -P mamba &>/dev/null; then
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
        echo "Environment '${env_name}' is not installed."
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
        URL_mid="linux/x64/1.9"
        tarball="julia-1.9.3-linux-x86_64.tar.gz"
        ;;
    "Darwin")
        if [[ "${arch}" = "x86_64" ]]; then
            URL_mid="mac/x64/1.9"
            tarball="julia-1.9.3-mac64.tar.gz"
        elif [[ "${arch}" = "arm64" ]]; then
            URL_mid="mac/aarch64/1.9"
            tarball="julia-1.9.3-macaarch64.tar.gz"
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


#  Install Julia ==============================================================
#  Set flags for running echo tests, operations, etc.
check_variables=true  # Echo the variables assigned above
check_binary=true     # Check if the Julia binary is installed/in PATH
check_operation=true  # Check the operation to download and install Julia
run_operation=false   # Run the operation to download and install Julia
update_path=true      # Update PATH to include Julia binary

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
    if type -P julia &>/dev/null; then
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
        local shell_config
        if [[ -f \"\${HOME}/.bashrc\" ]]; then
            shell_config=\"\${HOME}/.bashrc\"
        elif [[ -f \"\${HOME}/.bash_profile\" ]]; then
            shell_config=\"\${HOME}/.bash_profile\"
        else
            error_and_return \"No known shell configuration file found.\"
        fi

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
        local shell_config
        if [[ -f "${HOME}/.bashrc" ]]; then
            shell_config="${HOME}/.bashrc"
        elif [[ -f "${HOME}/.bash_profile" ]]; then
            shell_config="${HOME}/.bash_profile"
        else
            error_and_return "No known shell configuration file found."
        fi

        #  Call the function to update the configuration file
        update_shell_config "${shell_config}" "${untarred}"

        echo "To apply the update to PATH, please restart the terminal or"
        echo "source the configuration file."
    else
        echo "For ${untarred} to remain in PATH, please export ${untarred} to"
        echo "PATH in the configuration file."
    fi
fi


#  Install Atria dependencies =================================================
#  Set flag(s)
create_mamba_env=false  # Install mamba environment if not detected
update_path=true        # Update PATH to include Atria binary  #TODO

#  Check that Mamba is installed and in PATH
check_mamba_installed

#  Check that environment assigned to env_name is installed
check_env_installed "${env_name}"

#  If environment assigned to env_name is not installed, run the following
if [[ $? -eq 0 ]]; then
    #  Handle the case when the environment is already installed
    echo "Activating environment ${env_name}"
    
    if ! mamba activate "${env_name}" &> dev/null; then
        #  If `mamba activate` fails, try using `source activate`
        echo "Mamba activation failed. Trying with source activate..."

        if ! source activate "${env_name}" &>/dev/null; then
            #  If `source activate` also fails, return an error
            error_and_return "Failed to activate environment \"${env_name}\"."
        else
            echo "Environment \"${env_name}\" activated using source activate."
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


#  Install Atria ==============================================================
#  Set flags
install_atria=false  # Install Atria or not

if ${install_atria}; then
    #  Check if git and Julia are available
    if ! type -P git &>/dev/null; then
        error_and_return "Error: git is not installed or not in the PATH."
    fi
    if ! type -P julia &>/dev/null; then
        error_and_return "Error: Julia is not installed or not in the PATH."
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
        || error_and_return "Error: Failed to change to Atria directory."

    #  Environment containing Atria dependencies must be activated prior to
    #+ installation of Atria
    if [[ "${CONDA_DEFAULT_ENV}" != "${env_name}" ]]; then
        if [[ "${CONDA_DEFAULT_ENV}" != "base" ]]; then
            mamba deactivate
        fi

        if ! mamba activate "${env_name}" &> dev/null; then
            #  If `mamba activate` fails, try using `source activate`
            echo "Mamba activation failed. Trying with source activate..."
            if ! source activate "${env_name}" &>/dev/null; then
                #  If `source activate` also fails, return an error
                error_and_return "Failed to activate environment \"${env_name}\"."
            fi
        fi
    fi

    #  Run the Julia script to build Atria
    if ! julia build_atria.jl; then
        error_and_return "Failed to build Atria."
    fi

    echo "Atria installed successfully."

    #TODO
    #  Add the trimming program to PATH if not already present
    if ! grep -q "${trim_prog_dir}/bin" <<< "$PATH"; then
        export PATH="${PATH}:${trim_prog_dir}/bin"
        local shell_config="${HOME}/.bashrc"  # Adjust based on the user's shell
        echo "export PATH=\"${PATH}:${trim_prog_dir}/bin\"" >> "${shell_config}"
        echo "Path updated in ${shell_config}. Please restart the terminal or source the configuration file."
    fi
fi
```
</details>
<br />

<a id="notes-1"></a>
#### Notes
<details>
<summary><i>Notes: Install Atria and dependencies</i></summary>
<br />

`#TODO` Carefully explain all concepts in the above chunk.  
`#TODO` Carefully test the code in the above chunk.
</details>
<br />

<a id="adapter--and-quality-trim-the-fastq-files-using-atria"></a>
### Adapter- and quality-trim the FASTQ files using Atria
<a id="code-2"></a>
#### Code
<details>
<summary><i>Code: Adapter- and quality-trim the FASTQ files using Atria</i></summary>

```bash
#!/bin/bash

#  Initialize variables and arrays ============================================
dir_base="${HOME}/tsukiyamalab"                                    # Base directory for lab data
dir_repo="Kris/2023_rDNA"                                          # Repository directory
dir_work="results/2023-0406_tutorial_ChIP-seq_analyses"            # Work directory
dir_sym="01_sym"                                                   # Directory with symlinked FASTQs
dir_trim="02_trim"                                                 # Directory for trimmed FASTQs
env_atria="pairtools_env"                                          # Conda environment for Atria
path_atria="${dir_base}/${dir_repo}/src/Atria/app-3.2.2/bin/atria" # Atria executable path
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
)


#  Do main work ===============================================================
#  Set flags for checking variable and array assignments
check_variables=true
check_array=false

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
    || echo "cd'ing failed; check on this"

#  If it doesn't exist, then create a directory to store symlinked FASTQ files
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
#SBATCH --time=2:00:00
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
            Warning: Trimmed fastqs for $(basename ${read_1%_R1.fastq.gz})
                     exist; skipping trimming
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
#SBATCH --time=2:00:00
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
            Warning: Trimmed fastqs for $(basename ${read_1%_R1.fastq.gz})
                     exist; skipping trimming
            "
        fi
    fi

    sleep 0.2  # Short pause to prevent rapid job-submission overload
done
```
</details>
<br />

<a id="notes-2"></a>
#### Notes
<details>
<summary><i>Notes: Adapter- and quality-trim the FASTQ files using Atria</i></summary>


</details>
<br />
