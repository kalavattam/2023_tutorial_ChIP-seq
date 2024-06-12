
`#coverage-and-scaling_figures.md`
<br />

<details>
<summary><b><font size="+1"><i>Table of contents</i></font></b></summary>
<br />
<!-- MarkdownTOC -->

1. [Visualize raw coverage with bamCoverage: `flag-2_mapq-1`](#visualize-raw-coverage-with-bamcoverage-flag-2_mapq-1)
    1. [Bash code](#bash-code)

<!-- /MarkdownTOC -->
</details>
<br />

<a id="visualize-raw-coverage-with-bamcoverage-flag-2_mapq-1"></a>
## Visualize raw coverage with bamCoverage: `flag-2_mapq-1`
<a id="bash-code"></a>
### Bash code
<details>
<summary><i>Bash code: Raw coverage with deepTools</i></summary>

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
    if ! command -v mamba &> /dev/null; then
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
#  Set flags for filtering reads by flag and MAPQ, and for plotting coverage
#+ for read alignments or determined fragments
req_flag=false                # If true, filter BAM for flag bit 2 and MAPQ 1
req_frag=true                 # If true, draw coverage for fragments; if false,
                              # draw coverage from aligned reads

#  Set flags for checking variable assignmens, checking array assignments,
#+ checking jobs to submit, and submitting those jobs to SLURM
check_variables=true
check_array=false
print_jobs=true
submit_jobs=true

#  Assign job- and bamCoverage-related variables
env_name="deeptools_env"      # Virtual environment containing deepTools

time="0:30:00"                # Job time for SLURM
threads=4                     # Threads for deepTools/SLURM jobs

region="I"                    # Plot coverage for this region of the genome
format="bigwig"               # Format of coverage outfile
bin_size=1                    # Size of bins (bp) in outfile
normalization="None"          # "None" means don't normalize the outfile values

dir_proj="${HOME}/projects-etc/2023_tutorial_ChIP-seq"
if ${req_flag}; then
    dir_exp="flag-2_mapq-1"
else
    dir_exp="flag-NA_mapq-0"
fi
dir_bams="${dir_proj}/03_bam/bowtie2/${dir_exp}/bam"
dir_figs="${dir_proj}/04_figures/coverage-and-scaling/bowtie2/${dir_exp}"

unset file_bams && typeset -A file_bams=(
    ["IP_G1_Hho1_6336.sort-coord.bam"]="in_G1_Hho1_6336.sort-coord.bam"
    ["IP_G1_Hho1_6337.sort-coord.bam"]="in_G1_Hho1_6337.sort-coord.bam"
    ["IP_G1_Hmo1_7750.sort-coord.bam"]="in_G1_Hmo1_7750.sort-coord.bam"
    ["IP_G1_Hmo1_7751.sort-coord.bam"]="in_G1_Hmo1_7751.sort-coord.bam"
    ["IP_G2M_Hho1_6336.sort-coord.bam"]="in_G2M_Hho1_6336.sort-coord.bam"
    ["IP_G2M_Hho1_6337.sort-coord.bam"]="in_G2M_Hho1_6337.sort-coord.bam"
    ["IP_G2M_Hmo1_7750.sort-coord.bam"]="in_G2M_Hmo1_7750.sort-coord.bam"
    ["IP_G2M_Hmo1_7751.sort-coord.bam"]="in_G2M_Hmo1_7751.sort-coord.bam"
    ["IP_Q_Hho1_6336.sort-coord.bam"]="in_Q_Hho1_6336.sort-coord.bam"
    ["IP_Q_Hho1_6337.sort-coord.bam"]="in_Q_Hho1_6337.sort-coord.bam"
    ["IP_Q_Hmo1_7750.sort-coord.bam"]="in_Q_Hmo1_7750.sort-coord.bam"
    ["IP_Q_Hmo1_7751.sort-coord.bam"]="in_Q_Hmo1_7751.sort-coord.bam"
)


#  Do the main work ===========================================================
if ${check_variables}; then
    echo "
    req_flag=${req_flag}
    req_frag=${req_frag}
    
    env_name=${env_name}

    time=${time}
    threads=${threads}

    region=${region}
    format=${format}
    bin_size=${bin_size}
    normalization=${normalization}

    dir_proj=${dir_proj}
    dir_exp=${dir_exp}
    dir_bams=${dir_bams}
    dir_figs=${dir_figs}
    "
fi

if ${check_array}; then
    for key in "${!file_bams[@]}"; do
        echo "  key  ${key}"
        echo "value  ${file_bams[${key}]}"
        echo ""
    done
fi

run_chunk=false
if ${run_chunk}; then
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
    else
        error_and_return "Environment ${env_name} is not installed"

        if ${create_mamba_env}; then
            #  Handle the case when the environment is not installed
            echo "Creating environment ${env_name}"

            #  Switch `--yes` is not set, which means user input is required
            #NOTE Running this on FHCC Rhino; ergo, no CONDA_SUBDIR=osx-64
            # mamba create \
            #     --name "${env_name}" \
            #     --channel bioconda \
            #         bamtools \
            #         bedtools \
            #         bowtie2 \
            #         bwa \
            #         fastqc \
            #         minimap \
            #         mosdepth \
            #         picard \
            #         preseq \
            #         samtools \
            #         ucsc-bedgraphtobigwig \
            #         ucsc-bedsort \
            #         ucsc-facount
            #TODO Set up call to `mamba create` for deeptools_env
            
            activate_env "${env_name}"
        fi
    fi

    #  If not present, create the directories for outfiles
    if [[ ! -d "${dir_figs}" ]]; then
        mkdir -p ${dir_figs}/{reads,fragments}/err_out
    fi
fi

for key in "${!file_bams[@]}"; do
    IP="${dir_bams}/${key}"                   # ls -lhaFG "${IP}"
    input="${dir_bams}/${file_bams[${key}]}"  # ls -lhaFG "${input}"

    if ${req_frag}; then
        for file in "${IP}" "${input}"; do
            job_name="$(basename ${file} .bam)"

            if ${print_jobs}; then
                echo "
sbatch << EOF
#!/bin/bash

#SBATCH --job-name=\"${job_name}\"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=${threads}
#SBATCH --time=${time}
#SBATCH --error=\"${dir_figs}/fragments/err_out/${job_name}.%A.stderr.txt\"
#SBATCH --output=\"${dir_figs}/fragments/err_out/${job_name}.%A.stdout.txt\"

bamCoverage \\
    --numberOfProcessors \"${threads}\" \\
    --bam \"${file}\" \\
    --outFileName \"${dir_figs}/fragments/$(basename ${file} .bam).${region}.raw.${format}\" \\
    --outFileFormat \"${format}\" \\
    --binSize \"${bin_size}\" \\
    --region \"${region}\" \\
    --normalizeUsing \"${normalization}\" \\
    --samFlagInclude 64 \\
    --extendReads
EOF
                "
            fi
        
            if ${submit_jobs}; then
sbatch << EOF
#!/bin/bash

#SBATCH --job-name="${job_name}"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=${threads}
#SBATCH --time=${time}
#SBATCH --error="${dir_figs}/fragments/err_out/${job_name}.%A.stderr.txt"
#SBATCH --output="${dir_figs}/fragments/err_out/${job_name}.%A.stdout.txt"

bamCoverage \
    --numberOfProcessors "${threads}" \
    --bam "${file}" \
    --outFileName "${dir_figs}/fragments/$(basename ${file} .bam).${region}.raw.${format}" \
    --outFileFormat "${format}" \
    --binSize "${bin_size}" \
    --region "${region}" \
    --normalizeUsing "${normalization}" \
    --samFlagInclude 64 \
    --extendReads
EOF
            fi

            sleep 0.2
        done
    else
        for file in "${IP}" "${input}"; do
            job_name="$(basename ${file} .bam)"

            if ${print_jobs}; then
                echo "
sbatch << EOF
#!/bin/bash

#SBATCH --job-name=\"${job_name}\"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=${threads}
#SBATCH --time=${time}
#SBATCH --error=\"${dir_figs}/reads/err_out/${job_name}.%A.stderr.txt\"
#SBATCH --output=\"${dir_figs}/reads/err_out/${job_name}.%A.stdout.txt\"

bamCoverage \\
    --numberOfProcessors \"${threads}\" \\
    --bam \"${file}\" \\
    --outFileName \"${dir_figs}/reads/$(basename ${file} .bam).${region}.raw.${format}\" \\
    --outFileFormat \"${format}\" \\
    --binSize \"${bin_size}\" \\
    --region \"${region}\" \\
    --normalizeUsing \"${normalization}\"
EOF
                "
            fi
        
            if ${submit_jobs}; then
sbatch << EOF
#!/bin/bash

#SBATCH --job-name="${job_name}"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=${threads}
#SBATCH --time=${time}
#SBATCH --error="${dir_figs}/reads/err_out/${job_name}.%A.stderr.txt"
#SBATCH --output="${dir_figs}/reads/err_out/${job_name}.%A.stdout.txt"

bamCoverage \
    --numberOfProcessors "${threads}" \
    --bam "${file}" \
    --outFileName "${dir_figs}/reads/$(basename ${file} .bam).${region}.raw.${format}" \
    --outFileFormat "${format}" \
    --binSize "${bin_size}" \
    --region "${region}" \
    --normalizeUsing "${normalization}"
EOF
            fi

            sleep 0.2
        done
    fi
done
```

</details>
<br />