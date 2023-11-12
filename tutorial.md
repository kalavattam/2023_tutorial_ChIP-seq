
`#tutorial.md`
<br />

<details>
<summary><b><font size="+2"><i>Table of contents</i></font></b></summary>
<br />
<!-- MarkdownTOC -->

1. [1. Prepare FASTQ files of interest for processing and analyses](#1-prepare-fastq-files-of-interest-for-processing-and-analyses)
    1. [Code](#code)
    1. [Notes](#notes)
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
    1. [Code](#code-1)

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
    || echo "cd'ing failed; check on this"

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
The `ln` command creates symbolic links (symlinks), which are pointers to original files. The `-s` option creates a symlink. The command format is `ln -s original_file symlink_file`. This maintains the integrity of raw data while simplifying access under new names.

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
<a id="code-1"></a>
### Code
<details>
<summary><i>Code: 2. Adapter- and quality-trim the FASTQ files</i></summary>

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

