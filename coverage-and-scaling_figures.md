
`#coverage-and-scaling_figures.md`
<br />

<details>
<summary><b><font size="+1"><i>Table of contents</i></font></b></summary>
<br />
<!-- MarkdownTOC -->

1. [Raw coverage with deepTools](#raw-coverage-with-deeptools)
    1. [Bash code](#bash-code)

<!-- /MarkdownTOC -->
</details>
<br />

<a id="raw-coverage-with-deeptools"></a>
## Raw coverage with deepTools
<a id="bash-code"></a>
### Bash code
<details>
<summary><i>Bash code: Raw coverage with deepTools</i></summary>

```bash
#!/bin/bash

dir_proj="${HOME}/projects-etc/2023_tutorial_ChIP-seq"
dir_bams="${dir_proj}/03_bam/bowtie2/bam"
dir_figs="${dir_proj}/04_figures/coverage-and-scaling/bowtie2"

env="deeptools_env"

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

check_array=true
run_for_reads=false
run_for_fragments=true
print_jobs=true
submit_jobs=false

threads=4
time="0:30:00"

region="I"
format="bigwig"
bin_size=1
normalization="None"


cd "${dir_proj}" || echo "cd'ing failed; check on this"

if [[ ! -d "${dir_figs}" ]]; then
    mkdir -p ${dir_figs}/{reads,fragments}/err_out
fi

source activate "${env}"

if ${check_array}; then
    for key in "${!file_bams[@]}"; do
        echo "  key  ${key}"
        echo "value  ${file_bams[${key}]}"
        echo ""
    done
fi


for key in "${!file_bams[@]}"; do
    IP=${key}
    input=${file_bams[${key}]}

    if ${run_for_reads}; then
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
    --bam \"${dir_bams}/${file}\" \\
    --outFileName \"${dir_figs}/reads/${file/.bam/.${region}.raw.${format}}\" \\
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
    --bam "${dir_bams}/${file}" \
    --outFileName "${dir_figs}/reads/${file/.bam/.${region}.raw.${format}}" \
    --outFileFormat "${format}" \
    --binSize "${bin_size}" \
    --region "${region}" \
    --normalizeUsing "${normalization}"
EOF
            fi

            sleep 0.2
        done
    fi

    if 
done
```

</details>
<br />