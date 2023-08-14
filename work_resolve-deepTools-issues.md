
`#work_resolve-deepTools-issues.md`
<br />
<br />

<!-- MarkdownTOC -->

1. [Get situated](#get-situated)
    1. [Code](#code)
1. [Align the datasets](#align-the-datasets)
    1. [Submit job to make `bowtie2` indices for S288C R64-3-1](#submit-job-to-make-bowtie2-indices-for-s288c-r64-3-1)
        1. [Code](#code-1)
        1. [Printed](#printed)
    1. [Run `bowtie2` alignment](#run-bowtie2-alignment)
        1. [Code](#code-2)
        1. [Printed](#printed-1)
    1. [Examine flags in bam outfiles](#examine-flags-in-bam-outfiles)
        1. [Code](#code-3)
        1. [Printed](#printed-2)
    1. [Fix chromosome names in bams](#fix-chromosome-names-in-bams)
        1. [Get `jvarkit` `bamrenamechr` up and running](#get-jvarkit-bamrenamechr-up-and-running)
            1. [Code](#code-4)
            1. [Printed](#printed-3)
        1. [Run `jvarkit` `bamrenamechr` on problem bams](#run-jvarkit-bamrenamechr-on-problem-bams)
            1. [Code](#code-5)
            1. [Printed](#printed-4)
1. [Correct chromosome names in `S288C_reference_sequence_R64-3-1_20210421.fa`](#correct-chromosome-names-in-s288c_reference_sequence_r64-3-1_20210421fa)
    1. [Code](#code-6)
    1. [Printed](#printed-5)

<!-- /MarkdownTOC -->
<br />
<br />

<a id="get-situated"></a>
## Get situated
<a id="code"></a>
### Code
<details>
<summary><i>Code: Get situated</i></summary>

```bash
#!/bin/bash

p_data="${HOME}/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802"
p_rename="${p_data}/renamed"

if [[ ! -d "${p_rename}" ]]; then mkdir -p "${p_rename}"; fi

typeset -A data=(
    ["${p_data}/SRR7175375.fastq"]="${p_rename}/Brn1_Q_rep1_ChIP.fastq"
    ["${p_data}/SRR7175376.fastq"]="${p_rename}/Brn1_Q_rep1_input.fastq"
    ["${p_data}/SRR7175377.fastq"]="${p_rename}/Brn1_Q_rep2_ChIP.fastq"
    ["${p_data}/SRR7175378.fastq"]="${p_rename}/Brn1_Q_rep2_input.fastq"
    ["${p_data}/SRR7175379.fastq"]="${p_rename}/Brn1_Q_rep3_ChIP.fastq"
    ["${p_data}/SRR7175380.fastq"]="${p_rename}/Brn1_Q_rep3_input.fastq"
    ["${p_data}/SRR7175382.fastq"]="${p_rename}/Brn1_Q_all_input.fastq"
    ["${p_data}/SRR7175381.fastq"]="${p_rename}/Brn1_Q_all_ChIP.fastq"
    ["${p_data}/SRR7175367.fastq"]="${p_rename}/Brn1_log_rep1_ChIP.fastq"
    ["${p_data}/SRR7175368.fastq"]="${p_rename}/Brn1_log_rep1_input.fastq"
    ["${p_data}/SRR7175369.fastq"]="${p_rename}/Brn1_log_rep2_ChIP.fastq"
    ["${p_data}/SRR7175370.fastq"]="${p_rename}/Brn1_log_rep2_input.fastq"
    ["${p_data}/SRR7175371.fastq"]="${p_rename}/Brn1_log_rep3_ChIP.fastq"
    ["${p_data}/SRR7175372.fastq"]="${p_rename}/Brn1_log_rep3_input.fastq"
    ["${p_data}/SRR7175373.fastq"]="${p_rename}/Brn1_log_all_ChIP.fastq"
    ["${p_data}/SRR7175374.fastq"]="${p_rename}/Brn1_log_all_input.fastq"
)


#  Rather than rename the initial files, make "symbolic links" to the files
#+ with new, better-detailed names
for i in "${!data[@]}"; do
    echo """
      key: ${i}
    value: ${data[${i}]}
    """

    ln -s "${i}" "${data["${i}"]}"
done
```
</details>
<br />
<br />

<a id="align-the-datasets"></a>
## Align the datasets
<a id="submit-job-to-make-bowtie2-indices-for-s288c-r64-3-1"></a>
### Submit job to make `bowtie2` indices for S288C R64-3-1
<a id="code-1"></a>
#### Code
<details>
<summary><i>Code: Submit job to make bowtie2 indices for S288C R64-3-1</i></summary>

```bash
#!/bin/bash

source activate coverage_env
which bowtie2-build
bowtie2-build --help

job_name="bowtie2-build"  # echo "${job_name}"
threads=8  # echo "${threads}"
d_genome="${HOME}/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421"  # ., "${d_genome}"
f_genome="${d_genome}/fasta/S288C_reference_sequence_R64-3-1_20210421.fa"  # ., "${f_genome}"
o_indices="${d_genome}/bowtie2/$(basename "${f_genome}" .fa)"  # echo "${o_indices}"
sh_err_out="$(dirname "${d_genome}/bowtie2/$(basename "${f_genome}" .fa)")/sh_err_out"  # echo "${sh_err_out}"

if [[ ! -d "${d_genome}/bowtie2" ]]; then mkdir -p "${d_genome}/bowtie2/sh_err_out"; fi

echo """
sbatch \\
    --job-name=${job_name} \\
    --nodes=1 \\
    --cpus-per-task=${threads} \\
    --error=${sh_err_out}/${job_name}.%A.stderr.txt \\
    --output=${sh_err_out}/${job_name}.%A.stdout.txt \\
    bowtie2-build \\
        --threads ${threads} \\
        ${f_genome} \\
        ${o_indices} \\
             > >(tee -a ${sh_err_out}/${job_name}.stdout.txt) \\
            2> >(tee -a ${sh_err_out}/${job_name}.stderr.txt)
"""

if [[ -f "${sh_err_out}/${job_name}.sh" ]]; then
    rm "${sh_err_out}/${job_name}.sh"
fi
cat << script > "${sh_err_out}/${job_name}.sh"
#!/bin/bash

#SBATCH --job-name=${job_name}
#SBATCH --nodes=1
#SBATCH --cpus-per-task=${threads}
#SBATCH --error=${sh_err_out}/${job_name}.%A.stderr.txt
#SBATCH --output=${sh_err_out}/${job_name}.%A.stdout.txt

#  ${job_name}.sh
#  KA
#  $(date '+%Y-%m%d')

eval "\$(conda shell.bash hook)"
source activate coverage_env

bowtie2-build \\
    --threads ${threads} \\
    ${f_genome} \\
    ${o_indices}
script
# cat "${sh_err_out}/${job_name}.sh"

sbatch "${sh_err_out}/${job_name}.sh"
```
</details>
<br />

<a id="printed"></a>
#### Printed
<details>
<summary><i>Printed: Submit job to make bowtie2 indices for S288C R64-3-1</i></summary>

```txt
❯ source activate coverage_env


❯ which bowtie2-build
/home/kalavatt/miniconda3/envs/coverage_env/bin/bowtie2-build


❯ bowtie2-build --help
Bowtie 2 version 2.5.1 by Ben Langmead (langmea@cs.jhu.edu, www.cs.jhu.edu/~langmea)
Usage: bowtie2-build [options]* <reference_in> <bt2_index_base>
    reference_in            comma-separated list of files with ref sequences
    bt2_index_base          write bt2 data to files with this dir/basename
*** Bowtie 2 indexes will work with Bowtie v1.2.3 and later. ***
Options:
    -f                      reference files are Fasta (default)
    -c                      reference sequences given on cmd line (as
                            <reference_in>)
    --large-index           force generated index to be 'large', even if ref
                            has fewer than 4 billion nucleotides
    --debug                 use the debug binary; slower, assertions enabled
    --sanitized             use sanitized binary; slower, uses ASan and/or UBSan
    --verbose               log the issued command
    -a/--noauto             disable automatic -p/--bmax/--dcv memory-fitting
    -p/--packed             use packed strings internally; slower, less memory
    --bmax <int>            max bucket sz for blockwise suffix-array builder
    --bmaxdivn <int>        max bucket sz as divisor of ref len (default: 4)
    --dcv <int>             diff-cover period for blockwise (default: 1024)
    --nodc                  disable diff-cover (algorithm becomes quadratic)
    -r/--noref              don't build .3/.4 index files
    -3/--justref            just build .3/.4 index files
    -o/--offrate <int>      SA is sampled every 2^<int> BWT chars (default: 5)
    -t/--ftabchars <int>    # of chars consumed in initial lookup (default: 10)
    --threads <int>         # of threads
    --seed <int>            seed for random number generator
    -q/--quiet              verbose output (for debugging)
    --h/--help              print this message and quit
    --version               print version information and quit


❯ if [[ ! -d "${d_genome}/bowtie2" ]]; then mkdir -p "${d_genome}/bowtie2/sh_err_out"; fi
mkdir: created directory '/home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/bowtie2'
mkdir: created directory '/home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/bowtie2/sh_err_out'


❯ echo """
> sbatch \\
>     --job-name=${job_name} \\
>     --nodes=1 \\
>     --cpus-per-task=${threads} \\
>     --error=${sh_err_out}/${job_name}.%A.stderr.txt \\
>     --output=${sh_err_out}/${job_name}.%A.stdout.txt \\
>     bowtie2-build \\
>         --threads ${threads} \\
>         ${f_genome} \\
>         ${o_indices} \\
>              > >(tee -a ${sh_err_out}/${job_name}.stdout.txt) \\
>             2> >(tee -a ${sh_err_out}/${job_name}.stderr.txt)
> """
sbatch \
    --job-name=bowtie2-build \
    --nodes=1 \
    --cpus-per-task=8 \
    --error=/home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/bowtie2/sh_err_out/bowtie2-build.%A.stderr.txt \
    --output=/home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/bowtie2/sh_err_out/bowtie2-build.%A.stdout.txt \
    bowtie2-build \
        --threads 8 \
        /home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/fasta/S288C_reference_sequence_R64-3-1_20210421.fa \
        /home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/bowtie2/S288C_reference_sequence_R64-3-1_20210421 \
             > >(tee -a /home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/bowtie2/sh_err_out/bowtie2-build.stdout.txt) \
            2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/bowtie2/sh_err_out/bowtie2-build.stderr.txt)


❯ cat << script > "${sh_err_out}/${job_name}.sh"
> #!/bin/bash
> 
> #SBATCH --job-name=${job_name}
> #SBATCH --nodes=1
> #SBATCH --cpus-per-task=${threads}
> #SBATCH --error=${sh_err_out}/${job_name}.%A.stderr.txt
> #SBATCH --output=${sh_err_out}/${job_name}.%A.stdout.txt
> 
> #  ${job_name}.sh
> #  KA
> #  $(date '+%Y-%m%d')
> 
> eval "\$(conda shell.bash hook)"
> source activate coverage_env
> 
> bowtie2-build \\
>     --threads ${threads} \\
>     ${f_genome} \\
>     ${o_indices}
> script


❯ cat "${sh_err_out}/${job_name}.sh"
#!/bin/bash

#SBATCH --job-name=bowtie2-build
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --error=/home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/bowtie2/sh_err_out/bowtie2-build.%A.stderr.txt
#SBATCH --output=/home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/bowtie2/sh_err_out/bowtie2-build.%A.stdout.txt

#  bowtie2-build.sh
#  KA
#  2023-0406

eval "$(conda shell.bash hook)"
source activate coverage_env

bowtie2-build \
    --threads 8 \
    /home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/fasta/S288C_reference_sequence_R64-3-1_20210421.fa \
    /home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/bowtie2/S288C_reference_sequence_R64-3-1_20210421
```
</details>
<br />

<a id="run-bowtie2-alignment"></a>
### Run `bowtie2` alignment
<a id="code-2"></a>
#### Code
<details>
<summary><i>Code: Run bowtie2 alignment</i></summary>

```bash
#!/bin/bash

tmux new -s align
grabnode  # 16 cores, defaults

source activate coverage_env
which bowtie2
bowtie2 --help

d_work="${HOME}/tsukiyamalab/Kris/2023_rDNA/results"  # ., "${d_work}"
if [[ ! -d "${d_work}/2023-0406" ]]; then
    mkdir -p "${d_work}/2023-0406/bams/err_out"
fi

cd "${d_work}"
pwd

threads="${SLURM_CPUS_ON_NODE}"  # echo "${threads}"
d_genome="${HOME}/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421"  # ., "${d_genome}"
f_genome="${d_genome}/fasta/S288C_reference_sequence_R64-3-1_20210421.fa"  # ., "${f_genome}"
f_indices="${d_genome}/bowtie2/$(basename "${f_genome}" .fa)"  # ., "${f_indices}"*
err_out="${d_work}/2023-0406/bams/err_out"  # ., "${err_out}"
d_bams="$(dirname ${err_out})"

p_data="${HOME}/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed"
unset fastqs
typeset -a fastqs=(
    "${p_data}/Brn1_Q_rep1_ChIP.fastq"
    "${p_data}/Brn1_Q_rep1_input.fastq"
    "${p_data}/Brn1_Q_rep2_ChIP.fastq"
    "${p_data}/Brn1_Q_rep2_input.fastq"
    "${p_data}/Brn1_Q_rep3_ChIP.fastq"
    "${p_data}/Brn1_Q_rep3_input.fastq"
    "${p_data}/Brn1_Q_all_input.fastq"
    "${p_data}/Brn1_Q_all_ChIP.fastq"
    "${p_data}/Brn1_log_rep1_ChIP.fastq"
    "${p_data}/Brn1_log_rep1_input.fastq"
    "${p_data}/Brn1_log_rep2_ChIP.fastq"
    "${p_data}/Brn1_log_rep2_input.fastq"
    "${p_data}/Brn1_log_rep3_ChIP.fastq"
    "${p_data}/Brn1_log_rep3_input.fastq"
    "${p_data}/Brn1_log_all_ChIP.fastq"
    "${p_data}/Brn1_log_all_input.fastq"
)
for i in "${fastqs[@]}"; do ., "${i}"; done

for i in "${fastqs[@]}"; do
    # i="${fastqs[0]}"  # echo "${i}"
    
    # echo """
    # bowtie2 \\
    #     -p ${threads} \\
    #     --very-sensitive-local \\
    #     -x ${f_indices} \\
    #     -U ${i} \\
    #         | samtools sort \\
    #             -@ ${threads} \\
    #             -l 0 \\
    #             -O bam \\
    #         | samtools view \\
    #             -@ ${threads} \\
    #             -O bam \\
    #             -o ${i%.fastq}.bam  \\
    #                  > >(tee -a ${err_out}/bowtie2_samtools-sort_samtools-view.stdout.txt) \\
    #                 2> >(tee -a ${err_out}/bowtie2_samtools-sort_samtools-view.stderr.txt)
    # """

    in_fastq="${i}"  # ., "${in_fastq}"
    out_bam="${d_bams}/$(basename "${in_fastq}" .fastq).bam"  # echo "${out_bam}"
    echo "#### $(basename ${in_fastq}) ####"
    {
        bowtie2 \
            -p "${threads}" \
            --very-sensitive-local \
            -x "${f_indices}" \
            -U "${in_fastq}" \
                | samtools sort \
                    -@ "${threads}" \
                    -l 0 \
                    -O "bam" \
                | samtools view \
                    -@ "${threads}" \
                    -O "bam" \
                    -o "${out_bam}"
    }  \
         >> >(tee -a ${err_out}/$(basename "${out_bam}" .bam).bowtie2_samtools-sort_samtools-view.stdout.txt) \
        2>> >(tee -a ${err_out}/$(basename "${out_bam}" .bam).bowtie2_samtools-sort_samtools-view.stderr.txt)

    echo ""
    if [[ -f "${out_bam}" ]]; then
        samtools index \
            -@ "${threads}" \
            "${out_bam}" \
                 > >(tee -a ${err_out}/$(basename "${out_bam}" .bam).samtools-index.stdout.txt) \
                2> >(tee -a ${err_out}/$(basename "${out_bam}" .bam).samtools-index.stderr.txt)
    fi
    echo ""
done
```
</details>
<br />

<a id="printed-1"></a>
#### Printed
<details>
<summary><i>Printed: Run bowtie2 alignment</i></summary>

```txt
❯ which bowtie2
/home/kalavatt/miniconda3/envs/coverage_env/bin/bowtie2


❯ bowtie2 --help
Bowtie 2 version 2.5.1 by Ben Langmead (langmea@cs.jhu.edu, www.cs.jhu.edu/~langmea)
Usage:
  bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r> | --interleaved <i> | -b <bam>} [-S <sam>]
...


❯ d_work="${HOME}/tsukiyamalab/Kris/2023_rDNA/results"  # ., "${d_work}"
❯ if [[ ! -d "${d_work}/2023-0406" ]]; then
>     mkdir -p "${d_work}/2023-0406/bams/err_out"
> fi
mkdir: created directory '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406'
mkdir: created directory '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams'
mkdir: created directory '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/err_out'


❯ cd "${d_work}"
❯ pwd
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results


❯ ., "${d_genome}"
total 21M
drwxr-x---  5 kalavatt  548 Apr  6 08:04 ./
drwxrwx--- 12 kalavatt  482 Apr  6 07:53 ../
drwxrwx---  3 kalavatt  426 Apr  6 08:41 bowtie2/
drwxrwx---  2 kalavatt  329 Mar  7 14:25 bwa/
drwxrwx---  2 kalavatt   62 Mar  7 14:26 fasta/
-rw-r-----  1 kalavatt 3.6M Apr 27  2021 gene_association_R64-3-1_20210421.sgd.gz
-rw-r-----  1 kalavatt 1.1M Apr 21  2021 NotFeature_R64-3-1_20210421.fasta.gz
-rw-r-----  1 kalavatt 3.7M Apr 21  2021 orf_coding_all_R64-3-1_20210421.fasta.gz
-rw-r-----  1 kalavatt 2.6M Apr 21  2021 orf_trans_all_R64-3-1_20210421.fasta.gz
-rw-r-----  1 kalavatt 187K Apr 21  2021 other_features_genomic_R64-3-1_20210421.fasta.gz
-rw-r-----  1 kalavatt  42K Apr 27  2021 rna_coding_R64-3-1_20210421.fasta.gz
-rw-r-----  1 kalavatt 3.7M Apr 21  2021 S288C_reference_sequence_R64-3-1_20210421.fsa.gz
-rw-r-----  1 kalavatt 5.1M Apr 27  2021 saccharomyces_cerevisiae_R64-3-1_20210421.gff.gz


❯ ., "${f_genome}"
-rw-r----- 1 kalavatt 12M Apr 21  2021 /home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/fasta/S288C_reference_sequence_R64-3-1_20210421.fa


❯ ., "${f_indices}"*
-rw-rw---- 1 kalavatt 7.9M Apr  6 08:41 /home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/bowtie2/S288C_reference_sequence_R64-3-1_20210421.1.bt2
-rw-rw---- 1 kalavatt 2.9M Apr  6 08:41 /home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/bowtie2/S288C_reference_sequence_R64-3-1_20210421.2.bt2
-rw-rw---- 1 kalavatt  161 Apr  6 08:41 /home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/bowtie2/S288C_reference_sequence_R64-3-1_20210421.3.bt2
-rw-rw---- 1 kalavatt 2.9M Apr  6 08:41 /home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/bowtie2/S288C_reference_sequence_R64-3-1_20210421.4.bt2
-rw-rw---- 1 kalavatt 7.9M Apr  6 08:41 /home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/bowtie2/S288C_reference_sequence_R64-3-1_20210421.rev.1.bt2
-rw-rw---- 1 kalavatt 2.9M Apr  6 08:41 /home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/bowtie2/S288C_reference_sequence_R64-3-1_20210421.rev.2.bt2


❯ ., "${err_out}"
total 80K
drwxrws--- 2 kalavatt  0 Apr  6 08:58 ./
drwxrws--- 3 kalavatt 25 Apr  6 08:58 ../


❯ p_data="${HOME}/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed"


❯ unset fastqs


❯ typeset -a fastqs=(
>     "${p_data}/Brn1_Q_rep1_ChIP.fastq"
>     "${p_data}/Brn1_Q_rep1_input.fastq"
>     "${p_data}/Brn1_Q_rep2_ChIP.fastq"
>     "${p_data}/Brn1_Q_rep2_input.fastq"
>     "${p_data}/Brn1_Q_rep3_ChIP.fastq"
>     "${p_data}/Brn1_Q_rep3_input.fastq"
>     "${p_data}/Brn1_Q_all_input.fastq"
>     "${p_data}/Brn1_Q_all_ChIP.fastq"
>     "${p_data}/Brn1_log_rep1_ChIP.fastq"
>     "${p_data}/Brn1_log_rep1_input.fastq"
>     "${p_data}/Brn1_log_rep2_ChIP.fastq"
>     "${p_data}/Brn1_log_rep2_input.fastq"
>     "${p_data}/Brn1_log_rep3_ChIP.fastq"
>     "${p_data}/Brn1_log_rep3_input.fastq"
>     "${p_data}/Brn1_log_all_ChIP.fastq"
>     "${p_data}/Brn1_log_all_input.fastq"
> )


❯ for i in "${fastqs[@]}"; do ., "${i}"; done
lrwxrwxrwx 1 kalavatt 76 Apr  6 07:17 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_Q_rep1_ChIP.fastq -> /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/SRR7175375.fastq
lrwxrwxrwx 1 kalavatt 76 Apr  6 07:17 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_Q_rep1_input.fastq -> /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/SRR7175376.fastq
lrwxrwxrwx 1 kalavatt 76 Apr  6 07:17 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_Q_rep2_ChIP.fastq -> /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/SRR7175377.fastq
lrwxrwxrwx 1 kalavatt 76 Apr  6 07:17 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_Q_rep2_input.fastq -> /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/SRR7175378.fastq
lrwxrwxrwx 1 kalavatt 76 Apr  6 07:17 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_Q_rep3_ChIP.fastq -> /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/SRR7175379.fastq
lrwxrwxrwx 1 kalavatt 76 Apr  6 07:17 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_Q_rep3_input.fastq -> /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/SRR7175380.fastq
lrwxrwxrwx 1 kalavatt 76 Apr  6 07:17 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_Q_all_input.fastq -> /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/SRR7175382.fastq
lrwxrwxrwx 1 kalavatt 76 Apr  6 07:17 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_Q_all_ChIP.fastq -> /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/SRR7175381.fastq
lrwxrwxrwx 1 kalavatt 76 Apr  6 07:17 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_log_rep1_ChIP.fastq -> /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/SRR7175367.fastq
lrwxrwxrwx 1 kalavatt 76 Apr  6 07:17 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_log_rep1_input.fastq -> /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/SRR7175368.fastq
lrwxrwxrwx 1 kalavatt 76 Apr  6 07:17 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_log_rep2_ChIP.fastq -> /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/SRR7175369.fastq
lrwxrwxrwx 1 kalavatt 76 Apr  6 07:17 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_log_rep2_input.fastq -> /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/SRR7175370.fastq
lrwxrwxrwx 1 kalavatt 76 Apr  6 07:17 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_log_rep3_ChIP.fastq -> /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/SRR7175371.fastq
lrwxrwxrwx 1 kalavatt 76 Apr  6 07:17 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_log_rep3_input.fastq -> /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/SRR7175372.fastq
lrwxrwxrwx 1 kalavatt 76 Apr  6 07:17 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_log_all_ChIP.fastq -> /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/SRR7175373.fastq
lrwxrwxrwx 1 kalavatt 76 Apr  6 07:17 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_log_all_input.fastq -> /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/SRR7175374.fastq


❯ for i in "${fastqs[@]}"; do
>     echo """
>     bowtie2 \\
>         -p ${threads} \\
>         --very-sensitive-local \\
>         -x ${f_indices} \\
>         -U ${i} \\
>             | samtools sort \\
>                 -@ ${threads} \\
>                 -l 0 \\
>                 -O bam \\
>             | samtools view \\
>                 -@ ${threads} \\
>                 -O bam \\
>                 -o ${i%.fastq}.bam  \\
>                      > >(tee -a ${err_out}/bowtie2_samtools-sort_samtools-view.stdout.txt) \\
>                     2> >(tee -a ${err_out}/bowtie2_samtools-sort_samtools-view.stderr.txt)
>     """
> 
>     # bowtie2 \
>     #     -p "${threads}" \
>     #     --very-sensitive-local \
>     #     -x "${f_indices}" \
>     #     -U "${i}" \
>     #         | samtools sort \
>     #             -@ "${threads}" \
>     #             -l 0 \
>     #             -O "bam" \
>     #         | samtools view \
>     #             -@ "${threads}" \
>     #             -O "bam" \
>     #             -o "${i%.fastq}.bam"  \
>     #                  > >(tee -a ${err_out}/bowtie2_samtools-sort_samtools-view.stdout.txt) \
>     #                 2> >(tee -a ${err_out}/bowtie2_samtools-sort_samtools-view.stderr.txt)
> 
>     if [[ -f "${i%.fastq}.bam" ]]; then
>         samtools index \
>             -@ "${threads}" \
>             "${i%.fastq}.bam" \
>                  > >(tee -a ${err_out}/samtools-index.stdout.txt) \
>                 2> >(tee -a ${err_out}/samtools-index.stderr.txt)
>     fi
> done

    bowtie2 \
        -p 16 \
        --very-sensitive-local \
        -x /home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/bowtie2/S288C_reference_sequence_R64-3-1_20210421 \
        -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_Q_rep1_ChIP.fastq \
            | samtools sort \
                -@ 16 \
                -l 0 \
                -O bam \
            | samtools view \
                -@ 16 \
                -O bam \
                -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_Q_rep1_ChIP.bam  \
                     > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/err_out/bowtie2_samtools-sort_samtools-view.stdout.txt) \
                    2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/err_out/bowtie2_samtools-sort_samtools-view.stderr.txt)


    bowtie2 \
        -p 16 \
        --very-sensitive-local \
        -x /home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/bowtie2/S288C_reference_sequence_R64-3-1_20210421 \
        -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_Q_rep1_input.fastq \
            | samtools sort \
                -@ 16 \
                -l 0 \
                -O bam \
            | samtools view \
                -@ 16 \
                -O bam \
                -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_Q_rep1_input.bam  \
                     > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/err_out/bowtie2_samtools-sort_samtools-view.stdout.txt) \
                    2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/err_out/bowtie2_samtools-sort_samtools-view.stderr.txt)


    bowtie2 \
        -p 16 \
        --very-sensitive-local \
        -x /home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/bowtie2/S288C_reference_sequence_R64-3-1_20210421 \
        -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_Q_rep2_ChIP.fastq \
            | samtools sort \
                -@ 16 \
                -l 0 \
                -O bam \
            | samtools view \
                -@ 16 \
                -O bam \
                -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_Q_rep2_ChIP.bam  \
                     > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/err_out/bowtie2_samtools-sort_samtools-view.stdout.txt) \
                    2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/err_out/bowtie2_samtools-sort_samtools-view.stderr.txt)


    bowtie2 \
        -p 16 \
        --very-sensitive-local \
        -x /home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/bowtie2/S288C_reference_sequence_R64-3-1_20210421 \
        -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_Q_rep2_input.fastq \
            | samtools sort \
                -@ 16 \
                -l 0 \
                -O bam \
            | samtools view \
                -@ 16 \
                -O bam \
                -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_Q_rep2_input.bam  \
                     > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/err_out/bowtie2_samtools-sort_samtools-view.stdout.txt) \
                    2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/err_out/bowtie2_samtools-sort_samtools-view.stderr.txt)


    bowtie2 \
        -p 16 \
        --very-sensitive-local \
        -x /home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/bowtie2/S288C_reference_sequence_R64-3-1_20210421 \
        -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_Q_rep3_ChIP.fastq \
            | samtools sort \
                -@ 16 \
                -l 0 \
                -O bam \
            | samtools view \
                -@ 16 \
                -O bam \
                -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_Q_rep3_ChIP.bam  \
                     > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/err_out/bowtie2_samtools-sort_samtools-view.stdout.txt) \
                    2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/err_out/bowtie2_samtools-sort_samtools-view.stderr.txt)


    bowtie2 \
        -p 16 \
        --very-sensitive-local \
        -x /home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/bowtie2/S288C_reference_sequence_R64-3-1_20210421 \
        -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_Q_rep3_input.fastq \
            | samtools sort \
                -@ 16 \
                -l 0 \
                -O bam \
            | samtools view \
                -@ 16 \
                -O bam \
                -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_Q_rep3_input.bam  \
                     > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/err_out/bowtie2_samtools-sort_samtools-view.stdout.txt) \
                    2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/err_out/bowtie2_samtools-sort_samtools-view.stderr.txt)


    bowtie2 \
        -p 16 \
        --very-sensitive-local \
        -x /home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/bowtie2/S288C_reference_sequence_R64-3-1_20210421 \
        -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_Q_all_input.fastq \
            | samtools sort \
                -@ 16 \
                -l 0 \
                -O bam \
            | samtools view \
                -@ 16 \
                -O bam \
                -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_Q_all_input.bam  \
                     > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/err_out/bowtie2_samtools-sort_samtools-view.stdout.txt) \
                    2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/err_out/bowtie2_samtools-sort_samtools-view.stderr.txt)


    bowtie2 \
        -p 16 \
        --very-sensitive-local \
        -x /home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/bowtie2/S288C_reference_sequence_R64-3-1_20210421 \
        -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_Q_all_ChIP.fastq \
            | samtools sort \
                -@ 16 \
                -l 0 \
                -O bam \
            | samtools view \
                -@ 16 \
                -O bam \
                -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_Q_all_ChIP.bam  \
                     > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/err_out/bowtie2_samtools-sort_samtools-view.stdout.txt) \
                    2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/err_out/bowtie2_samtools-sort_samtools-view.stderr.txt)


    bowtie2 \
        -p 16 \
        --very-sensitive-local \
        -x /home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/bowtie2/S288C_reference_sequence_R64-3-1_20210421 \
        -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_log_rep1_ChIP.fastq \
            | samtools sort \
                -@ 16 \
                -l 0 \
                -O bam \
            | samtools view \
                -@ 16 \
                -O bam \
                -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_log_rep1_ChIP.bam  \
                     > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/err_out/bowtie2_samtools-sort_samtools-view.stdout.txt) \
                    2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/err_out/bowtie2_samtools-sort_samtools-view.stderr.txt)


    bowtie2 \
        -p 16 \
        --very-sensitive-local \
        -x /home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/bowtie2/S288C_reference_sequence_R64-3-1_20210421 \
        -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_log_rep1_input.fastq \
            | samtools sort \
                -@ 16 \
                -l 0 \
                -O bam \
            | samtools view \
                -@ 16 \
                -O bam \
                -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_log_rep1_input.bam  \
                     > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/err_out/bowtie2_samtools-sort_samtools-view.stdout.txt) \
                    2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/err_out/bowtie2_samtools-sort_samtools-view.stderr.txt)


    bowtie2 \
        -p 16 \
        --very-sensitive-local \
        -x /home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/bowtie2/S288C_reference_sequence_R64-3-1_20210421 \
        -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_log_rep2_ChIP.fastq \
            | samtools sort \
                -@ 16 \
                -l 0 \
                -O bam \
            | samtools view \
                -@ 16 \
                -O bam \
                -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_log_rep2_ChIP.bam  \
                     > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/err_out/bowtie2_samtools-sort_samtools-view.stdout.txt) \
                    2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/err_out/bowtie2_samtools-sort_samtools-view.stderr.txt)


    bowtie2 \
        -p 16 \
        --very-sensitive-local \
        -x /home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/bowtie2/S288C_reference_sequence_R64-3-1_20210421 \
        -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_log_rep2_input.fastq \
            | samtools sort \
                -@ 16 \
                -l 0 \
                -O bam \
            | samtools view \
                -@ 16 \
                -O bam \
                -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_log_rep2_input.bam  \
                     > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/err_out/bowtie2_samtools-sort_samtools-view.stdout.txt) \
                    2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/err_out/bowtie2_samtools-sort_samtools-view.stderr.txt)


    bowtie2 \
        -p 16 \
        --very-sensitive-local \
        -x /home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/bowtie2/S288C_reference_sequence_R64-3-1_20210421 \
        -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_log_rep3_ChIP.fastq \
            | samtools sort \
                -@ 16 \
                -l 0 \
                -O bam \
            | samtools view \
                -@ 16 \
                -O bam \
                -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_log_rep3_ChIP.bam  \
                     > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/err_out/bowtie2_samtools-sort_samtools-view.stdout.txt) \
                    2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/err_out/bowtie2_samtools-sort_samtools-view.stderr.txt)


    bowtie2 \
        -p 16 \
        --very-sensitive-local \
        -x /home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/bowtie2/S288C_reference_sequence_R64-3-1_20210421 \
        -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_log_rep3_input.fastq \
            | samtools sort \
                -@ 16 \
                -l 0 \
                -O bam \
            | samtools view \
                -@ 16 \
                -O bam \
                -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_log_rep3_input.bam  \
                     > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/err_out/bowtie2_samtools-sort_samtools-view.stdout.txt) \
                    2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/err_out/bowtie2_samtools-sort_samtools-view.stderr.txt)


    bowtie2 \
        -p 16 \
        --very-sensitive-local \
        -x /home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/bowtie2/S288C_reference_sequence_R64-3-1_20210421 \
        -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_log_all_ChIP.fastq \
            | samtools sort \
                -@ 16 \
                -l 0 \
                -O bam \
            | samtools view \
                -@ 16 \
                -O bam \
                -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_log_all_ChIP.bam  \
                     > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/err_out/bowtie2_samtools-sort_samtools-view.stdout.txt) \
                    2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/err_out/bowtie2_samtools-sort_samtools-view.stderr.txt)


    bowtie2 \
        -p 16 \
        --very-sensitive-local \
        -x /home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/bowtie2/S288C_reference_sequence_R64-3-1_20210421 \
        -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_log_all_input.fastq \
            | samtools sort \
                -@ 16 \
                -l 0 \
                -O bam \
            | samtools view \
                -@ 16 \
                -O bam \
                -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/renamed/Brn1_log_all_input.bam  \
                     > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/err_out/bowtie2_samtools-sort_samtools-view.stdout.txt) \
                    2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/err_out/bowtie2_samtools-sort_samtools-view.stderr.txt)


❯ for i in "${fastqs[@]}"; do
>     # i="${fastqs[0]}"  # echo "${i}"
> 
>     # echo """
>     # bowtie2 \\
>     #     -p ${threads} \\
>     #     --very-sensitive-local \\
>     #     -x ${f_indices} \\
>     #     -U ${i} \\
>     #         | samtools sort \\
>     #             -@ ${threads} \\
>     #             -l 0 \\
>     #             -O bam \\
>     #         | samtools view \\
>     #             -@ ${threads} \\
>     #             -O bam \\
>     #             -o ${i%.fastq}.bam  \\
>     #                  > >(tee -a ${err_out}/bowtie2_samtools-sort_samtools-view.stdout.txt) \\
>     #                 2> >(tee -a ${err_out}/bowtie2_samtools-sort_samtools-view.stderr.txt)
>     # """
> 
>     in_fastq="${i}"  # ., "${in_fastq}"
>     out_bam="${d_bams}/$(basename "${in_fastq}" .fastq).bam"  # echo "${out_bam}"
>     echo "#### $(basename ${in_fastq}) ####"
>     {
>         bowtie2 \
>               -p "${threads}" \
>               --very-sensitive-local \
>               -x "${f_indices}" \
>               -U "${in_fastq}" \
>                   | samtools sort \
>                       -@ "${threads}" \
>                       -l 0 \
>                       -O "bam" \
>                   | samtools view \
>                       -@ "${threads}" \
>                       -O "bam" \
>                       -o "${out_bam}"
>     }  \
>          >> >(tee -a ${err_out}/$(basename "${out_bam}" .bam).bowtie2_samtools-sort_samtools-view.stdout.txt) \
>         2>> >(tee -a ${err_out}/$(basename "${out_bam}" .bam).bowtie2_samtools-sort_samtools-view.stderr.txt)
> 
>     echo ""
>     if [[ -f "${out_bam}" ]]; then
>         samtools index \
>             -@ "${threads}" \
>             "${out_bam}" \
>                  > >(tee -a ${err_out}/$(basename "${out_bam}" .bam).samtools-index.stdout.txt) \
>                 2> >(tee -a ${err_out}/$(basename "${out_bam}" .bam).samtools-index.stderr.txt)
>     fi
>     echo ""
> done
#### Brn1_Q_rep1_ChIP.fastq ####
7627718 reads; of these:
  7627718 (100.00%) were unpaired; of these:
    2074925 (27.20%) aligned 0 times
    3727744 (48.87%) aligned exactly 1 time
    1825049 (23.93%) aligned >1 times
72.80% overall alignment rate
[bam_sort_core] merging from 0 files and 16 in-memory blocks...


#### Brn1_Q_rep1_input.fastq ####
11331936 reads; of these:
  11331936 (100.00%) were unpaired; of these:
    504957 (4.46%) aligned 0 times
    8094263 (71.43%) aligned exactly 1 time
    2732716 (24.12%) aligned >1 times
95.54% overall alignment rate
[bam_sort_core] merging from 0 files and 16 in-memory blocks...


#### Brn1_Q_rep2_ChIP.fastq ####
6881074 reads; of these:
  6881074 (100.00%) were unpaired; of these:
    1342180 (19.51%) aligned 0 times
    3424381 (49.77%) aligned exactly 1 time
    2114513 (30.73%) aligned >1 times
80.49% overall alignment rate
[bam_sort_core] merging from 0 files and 16 in-memory blocks...


#### Brn1_Q_rep2_input.fastq ####
13059890 reads; of these:
  13059890 (100.00%) were unpaired; of these:
    544018 (4.17%) aligned 0 times
    8941678 (68.47%) aligned exactly 1 time
    3574194 (27.37%) aligned >1 times
95.83% overall alignment rate
[bam_sort_core] merging from 0 files and 16 in-memory blocks...


#### Brn1_Q_rep3_ChIP.fastq ####
9660608 reads; of these:
  9660608 (100.00%) were unpaired; of these:
    780964 (8.08%) aligned 0 times
    5687748 (58.88%) aligned exactly 1 time
    3191896 (33.04%) aligned >1 times
91.92% overall alignment rate
[bam_sort_core] merging from 0 files and 16 in-memory blocks...


#### Brn1_Q_rep3_input.fastq ####
11160513 reads; of these:
  11160513 (100.00%) were unpaired; of these:
    367228 (3.29%) aligned 0 times
    8780918 (78.68%) aligned exactly 1 time
    2012367 (18.03%) aligned >1 times
96.71% overall alignment rate
[bam_sort_core] merging from 0 files and 16 in-memory blocks...


#### Brn1_Q_all_input.fastq ####
35552339 reads; of these:
  35552339 (100.00%) were unpaired; of these:
    1416168 (3.98%) aligned 0 times
    25816812 (72.62%) aligned exactly 1 time
    8319359 (23.40%) aligned >1 times
96.02% overall alignment rate
[bam_sort_core] merging from 0 files and 16 in-memory blocks...


#### Brn1_Q_all_ChIP.fastq ####
24169400 reads; of these:
  24169400 (100.00%) were unpaired; of these:
    4198054 (17.37%) aligned 0 times
    12839822 (53.12%) aligned exactly 1 time
    7131524 (29.51%) aligned >1 times
82.63% overall alignment rate
[bam_sort_core] merging from 0 files and 16 in-memory blocks...


#### Brn1_log_rep1_ChIP.fastq ####
6643463 reads; of these:
  6643463 (100.00%) were unpaired; of these:
    635238 (9.56%) aligned 0 times
    3058150 (46.03%) aligned exactly 1 time
    2950075 (44.41%) aligned >1 times
90.44% overall alignment rate
[bam_sort_core] merging from 0 files and 16 in-memory blocks...


#### Brn1_log_rep1_input.fastq ####
15256990 reads; of these:
  15256990 (100.00%) were unpaired; of these:
    782379 (5.13%) aligned 0 times
    11368165 (74.51%) aligned exactly 1 time
    3106446 (20.36%) aligned >1 times
94.87% overall alignment rate
[bam_sort_core] merging from 0 files and 16 in-memory blocks...


#### Brn1_log_rep2_ChIP.fastq ####
32074691 reads; of these:
  32074691 (100.00%) were unpaired; of these:
    3436279 (10.71%) aligned 0 times
    14504051 (45.22%) aligned exactly 1 time
    14134361 (44.07%) aligned >1 times
89.29% overall alignment rate
[bam_sort_core] merging from 0 files and 16 in-memory blocks...


#### Brn1_log_rep2_input.fastq ####
12144595 reads; of these:
  12144595 (100.00%) were unpaired; of these:
    593652 (4.89%) aligned 0 times
    8936412 (73.58%) aligned exactly 1 time
    2614531 (21.53%) aligned >1 times
95.11% overall alignment rate
[bam_sort_core] merging from 0 files and 16 in-memory blocks...


#### Brn1_log_rep3_ChIP.fastq ####
8783472 reads; of these:
  8783472 (100.00%) were unpaired; of these:
    430229 (4.90%) aligned 0 times
    5187717 (59.06%) aligned exactly 1 time
    3165526 (36.04%) aligned >1 times
95.10% overall alignment rate
[bam_sort_core] merging from 0 files and 16 in-memory blocks...


#### Brn1_log_rep3_input.fastq ####
8922323 reads; of these:
  8922323 (100.00%) were unpaired; of these:
    307004 (3.44%) aligned 0 times
    6865933 (76.95%) aligned exactly 1 time
    1749386 (19.61%) aligned >1 times
96.56% overall alignment rate
[bam_sort_core] merging from 0 files and 16 in-memory blocks...


#### Brn1_log_all_ChIP.fastq ####
47501626 reads; of these:
  47501626 (100.00%) were unpaired; of these:
    4501732 (9.48%) aligned 0 times
    22749937 (47.89%) aligned exactly 1 time
    20249957 (42.63%) aligned >1 times
90.52% overall alignment rate
[bam_sort_core] merging from 0 files and 16 in-memory blocks...


#### Brn1_log_all_input.fastq ####
36323908 reads; of these:
  36323908 (100.00%) were unpaired; of these:
    1683052 (4.63%) aligned 0 times
    27170475 (74.80%) aligned exactly 1 time
    7470381 (20.57%) aligned >1 times
95.37% overall alignment rate
[bam_sort_core] merging from 0 files and 16 in-memory blocks...
```
</details>
<br />

<a id="examine-flags-in-bam-outfiles"></a>
### Examine flags in bam outfiles
<a id="code-3"></a>
#### Code
<details>
<summary><i>Code: Examine flags in bam outfiles</i></summary>

```bash
#!/bin/bash

#  Still in coverage_env (need samtools)
calculate_run_time() {
    what="""
    calculate_run_time()
    --------------------
    Calculate run time for chunk of code
    
    :param 1: start time in \$(date +%s) format
    :param 2: end time in \$(date +%s) format
    :param 3: message to be displayed when printing the run time <chr>
    
    #TODO Check that params are not empty or inappropriate formats or strings
    """
    run_time="$(echo "${2}" - "${1}" | bc -l)"
    
    echo ""
    echo "${3}"
    printf 'Run time: %dh:%dm:%ds\n' \
        $(( run_time/3600 )) \
        $(( run_time%3600/60 )) \
        $(( run_time%60 ))
    echo ""
}


display_spinning_icon() {
    what="""
    display_spinning_icon()
    -----------------------
    Display \"spinning icon\" while a background process runs
    
    :param 1: PID of the last program the shell ran in the background (int)
    :param 2: message to be displayed next to the spinning icon (chr)

    #TODO Checks...
    """
    spin="/|\\–"
    i=0
    while kill -0 "${1}" 2> /dev/null; do
        i=$(( (i + 1) % 4 ))
        printf "\r${spin:$i:1} %s" "${2}"
        sleep .15
    done
}


list_tally_flags() {
    what="""
    list_tally_flags()
    ------------------
    List and tally flags in a bam infile; function acts on a bam infile to
    perform piped commands (samtools view, cut, sort, uniq -c, sort -nr) that
    list and tally flags; function writes the results to a txt outfile, the
    name of which is derived from the txt infile
    
    :param 1: name of bam infile, including path (chr)

    #TODO Checks...
    """
    start="$(date +%s)"
    
    samtools view "${1}" \
        | cut -d$'\t' -f 2 \
        | sort \
        | uniq -c \
        | sort -nr \
            > "${1/.bam/.flags.txt}" &
    display_spinning_icon $! \
    "Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on $(basename "${1}")... "
        
    end="$(date +%s)"
    echo ""
    calculate_run_time "${start}" "${end}"  \
    "List and tally flags in $(basename "${1}")."
}


p_data="${HOME}/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams"
unset bams
typeset -a bams=(
    "${p_data}/Brn1_Q_rep1_ChIP.bam"
    "${p_data}/Brn1_Q_rep1_input.bam"
    "${p_data}/Brn1_Q_rep2_ChIP.bam"
    "${p_data}/Brn1_Q_rep2_input.bam"
    "${p_data}/Brn1_Q_rep3_ChIP.bam"
    "${p_data}/Brn1_Q_rep3_input.bam"
    "${p_data}/Brn1_Q_all_input.bam"
    "${p_data}/Brn1_Q_all_ChIP.bam"
    "${p_data}/Brn1_log_rep1_ChIP.bam"
    "${p_data}/Brn1_log_rep1_input.bam"
    "${p_data}/Brn1_log_rep2_ChIP.bam"
    "${p_data}/Brn1_log_rep2_input.bam"
    "${p_data}/Brn1_log_rep3_ChIP.bam"
    "${p_data}/Brn1_log_rep3_input.bam"
    "${p_data}/Brn1_log_all_ChIP.bam"
    "${p_data}/Brn1_log_all_input.bam"
)
# for i in "${bams[@]}"; do echo "${i}"; done
# for i in "${bams[@]}"; do ., "${i}"; done

for i in "${bams[@]}"; do
    # i="${bams[0]}"  # echo "${i}"
    echo "#### $(basename ${i}) ####"
    list_tally_flags "${i}"
    echo ""
done

for i in "${bams[@]}"; do
    if [[ -f "${i/.bam/.flags.txt}" ]]; then
        echo "#### $(basename ${i}) ####"
        cat "${i/.bam/.flags.txt}"
        echo ""
    fi
done

samtools view Brn1_log_all_ChIP.bam | head
#NOTE Weird chromosome names; need to fix them
```
</details>
<br />

<a id="printed-2"></a>
#### Printed
<details>
<summary><i>Printed: Examine flags in bam outfiles</i></summary>

```txt
❯ typeset -a bams=(
>     "${p_data}/Brn1_Q_rep1_ChIP.bam"
>     "${p_data}/Brn1_Q_rep1_input.bam"
>     "${p_data}/Brn1_Q_rep2_ChIP.bam"
>     "${p_data}/Brn1_Q_rep2_input.bam"
>     "${p_data}/Brn1_Q_rep3_ChIP.bam"
>     "${p_data}/Brn1_Q_rep3_input.bam"
>     "${p_data}/Brn1_Q_all_input.bam"
>     "${p_data}/Brn1_Q_all_ChIP.bam"
>     "${p_data}/Brn1_log_rep1_ChIP.bam"
>     "${p_data}/Brn1_log_rep1_input.bam"
>     "${p_data}/Brn1_log_rep2_ChIP.bam"
>     "${p_data}/Brn1_log_rep2_input.bam"
>     "${p_data}/Brn1_log_rep3_ChIP.bam"
>     "${p_data}/Brn1_log_rep3_input.bam"
>     "${p_data}/Brn1_log_all_ChIP.bam"
>     "${p_data}/Brn1_log_all_input.bam"
> )


❯ for i in "${bams}"; do echo "${i}"; done
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep1_ChIP.bam


❯ for i in "${bams[@]}"; do echo "${i}"; done
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep1_ChIP.bam
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep1_input.bam
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep2_ChIP.bam
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep2_input.bam
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep3_ChIP.bam
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep3_input.bam
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_all_input.bam
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_all_ChIP.bam
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep1_ChIP.bam
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep1_input.bam
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep2_ChIP.bam
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep2_input.bam
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep3_ChIP.bam
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep3_input.bam
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_all_ChIP.bam
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_all_input.bam


❯ for i in "${bams[@]}"; do ., "${i}"; done
-rw-rw---- 1 kalavatt 175M Apr  6 09:40 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep1_ChIP.bam
-rw-rw---- 1 kalavatt 279M Apr  6 09:40 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep1_input.bam
-rw-rw---- 1 kalavatt 154M Apr  6 09:41 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep2_ChIP.bam
-rw-rw---- 1 kalavatt 318M Apr  6 09:42 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep2_input.bam
-rw-rw---- 1 kalavatt 226M Apr  6 09:43 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep3_ChIP.bam
-rw-rw---- 1 kalavatt 262M Apr  6 09:44 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep3_input.bam
-rw-rw---- 1 kalavatt 814M Apr  6 09:47 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_all_input.bam
-rw-rw---- 1 kalavatt 532M Apr  6 09:50 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_all_ChIP.bam
-rw-rw---- 1 kalavatt 137M Apr  6 09:50 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep1_ChIP.bam
-rw-rw---- 1 kalavatt 378M Apr  6 09:53 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep1_input.bam
-rw-rw---- 1 kalavatt 614M Apr  6 09:57 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep2_ChIP.bam
-rw-rw---- 1 kalavatt 303M Apr  6 09:58 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep2_input.bam
-rw-rw---- 1 kalavatt 192M Apr  6 09:59 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep3_ChIP.bam
-rw-rw---- 1 kalavatt 213M Apr  6 10:00 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep3_input.bam
-rw-rw---- 1 kalavatt 911M Apr  6 10:05 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_all_ChIP.bam
-rw-rw---- 1 kalavatt 848M Apr  6 10:10 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_all_input.bam


❯ for i in "${bams[@]}"; do
>     # i="${bams[0]}"  # echo "${i}"
>     echo "#### $(basename ${i}) ####"
>     list_tally_flags "${i}"
>     echo ""
> done
#### Brn1_Q_rep1_ChIP.bam ####
[1] 757
/ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Brn1_Q_rep1_ChIP.bam... [1]+  Done
samtools view "${1}" | cut -d'  ' -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"


List and tally flags in Brn1_Q_rep1_ChIP.bam.
Run time: 0h:0m:7s


#### Brn1_Q_rep1_input.bam ####
[1] 818
/ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Brn1_Q_rep1_input.bam... [1]+  Done
samtools view "${1}" | cut -d' ' -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"


List and tally flags in Brn1_Q_rep1_input.bam.
Run time: 0h:0m:10s


#### Brn1_Q_rep2_ChIP.bam ####
[1] 904
| Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Brn1_Q_rep2_ChIP.bam... [1]+  Done
samtools view "${1}" | cut -d'  ' -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"


List and tally flags in Brn1_Q_rep2_ChIP.bam.
Run time: 0h:0m:6s


#### Brn1_Q_rep2_input.bam ####
[1] 959
| Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Brn1_Q_rep2_input.bam... [1]+  Done
samtools view "${1}" | cut -d' ' -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"


List and tally flags in Brn1_Q_rep2_input.bam.
Run time: 0h:0m:11s


#### Brn1_Q_rep3_ChIP.bam ####
[1] 1118
| Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Brn1_Q_rep3_ChIP.bam... [1]+  Done
samtools view "${1}" | cut -d'  ' -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"


List and tally flags in Brn1_Q_rep3_ChIP.bam.
Run time: 0h:0m:8s


#### Brn1_Q_rep3_input.bam ####
[1] 1219
\ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Brn1_Q_rep3_input.bam... [1]+  Done
samtools view "${1}" | cut -d' ' -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"


List and tally flags in Brn1_Q_rep3_input.bam.
Run time: 0h:0m:10s


#### Brn1_Q_all_input.bam ####
[1] 1298
\ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Brn1_Q_all_input.bam... [1]+  Done
samtools view "${1}" | cut -d'  ' -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"


List and tally flags in Brn1_Q_all_input.bam.
Run time: 0h:0m:32s


#### Brn1_Q_all_ChIP.bam ####
[1] 1554
\ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Brn1_Q_all_ChIP.bam... [1]+  Done
samtools view "${1}" | cut -d'   ' -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"


List and tally flags in Brn1_Q_all_ChIP.bam.
Run time: 0h:0m:20s


#### Brn1_log_rep1_ChIP.bam ####
[1] 1709
| Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Brn1_log_rep1_ChIP.bam... [1]+  Done
samtools view "${1}" | cut -d'        ' -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"


List and tally flags in Brn1_log_rep1_ChIP.bam.
Run time: 0h:0m:6s


#### Brn1_log_rep1_input.bam ####
[1] 1760
/ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Brn1_log_rep1_input.bam... [1]+  Done
samtools view "${1}" | cut -d'       ' -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"


List and tally flags in Brn1_log_rep1_input.bam.
Run time: 0h:0m:13s


#### Brn1_log_rep2_ChIP.bam ####
[1] 1858
– Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Brn1_log_rep2_ChIP.bam... [1]+  Done
samtools view "${1}" | cut -d'        ' -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"


List and tally flags in Brn1_log_rep2_ChIP.bam.
Run time: 0h:0m:27s


#### Brn1_log_rep2_input.bam ####
[1] 2055
| Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Brn1_log_rep2_input.bam... [1]+  Done
samtools view "${1}" | cut -d'       ' -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"


List and tally flags in Brn1_log_rep2_input.bam.
Run time: 0h:0m:11s


#### Brn1_log_rep3_ChIP.bam ####
[1] 2202
| Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Brn1_log_rep3_ChIP.bam... [1]+  Done
samtools view "${1}" | cut -d'        ' -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"


List and tally flags in Brn1_log_rep3_ChIP.bam.
Run time: 0h:0m:8s


#### Brn1_log_rep3_input.bam ####
[1] 2270
– Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Brn1_log_rep3_input.bam... [1]+  Done
samtools view "${1}" | cut -d'       ' -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"


List and tally flags in Brn1_log_rep3_input.bam.
Run time: 0h:0m:8s


#### Brn1_log_all_ChIP.bam ####
[1] 2335
– Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Brn1_log_all_ChIP.bam... [1]+  Done
samtools view "${1}" | cut -d' ' -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"


List and tally flags in Brn1_log_all_ChIP.bam.
Run time: 0h:0m:41s


#### Brn1_log_all_input.bam ####
[1] 2644
| Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Brn1_log_all_input.bam... [1]+  Done
samtools view "${1}" | cut -d'        ' -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"


List and tally flags in Brn1_log_all_input.bam.
Run time: 0h:0m:33s


#### Brn1_Q_rep1_ChIP.bam ####
2782361 0
2770432 16
2074925 4

#### Brn1_Q_rep1_input.bam ####
5416610 0
5410369 16
 504957 4

#### Brn1_Q_rep2_ChIP.bam ####
2770893 0
2768001 16
1342180 4

#### Brn1_Q_rep2_input.bam ####
6262802 0
6253070 16
 544018 4

#### Brn1_Q_rep3_ChIP.bam ####
4441272 16
4438372 0
 780964 4

#### Brn1_Q_rep3_input.bam ####
5400653 0
5392632 16
 367228 4

#### Brn1_Q_all_input.bam ####
17080375 0
17055796 16
1416168 4

#### Brn1_Q_all_ChIP.bam ####
9991973 0
9979373 16
4198054 4

#### Brn1_log_rep1_ChIP.bam ####
3014203 16
2994022 0
 635238 4

#### Brn1_log_rep1_input.bam ####
7239390 0
7235221 16
 782379 4

#### Brn1_log_rep2_ChIP.bam ####
14322187 0
14316225 16
3436279 4

#### Brn1_log_rep2_input.bam ####
5783228 0
5767715 16
 593652 4

#### Brn1_log_rep3_ChIP.bam ####
4176724 0
4176519 16
 430229 4

#### Brn1_log_rep3_input.bam ####
4309260 16
4306059 0
 307004 4

#### Brn1_log_all_ChIP.bam ####
21506703 16
21493191 0
4501732 4

#### Brn1_log_all_input.bam ####
17328645 0
17312211 16
1683052 4


❯ samtools view Brn1_log_all_ChIP.bam | head
SRR7175373.380283   0   ref|NC_001133|  1   25  31M1I18M    *   0   0   CCACACCACACCCACACACCCACACACCACACCCACACACCACACCACAC  GGGGGIIIIIIIIIIIIIIIIIIIIIIIGIIIIIIIIIGGGGIGGGGGGI  AS:i:90 XS:i:62 XN:i:0  XM:i:0  XO:i:1  XG:i:1  NM:i:1  MD:Z:49 YT:Z:UU
SRR7175373.767277   0   ref|NC_001133|  1   14  1S31M1I17M  *   0   0   CCCACACCACACCCACACACCCACACACCACACCCACACACCACACCACA  AGGAAGGGGGG<GA.<GGGAAGGIGGGG....<GGGGGA.A.<AGGGGAA  AS:i:88 XS:i:82 XN:i:0  XM:i:0  XO:i:1  XG:i:1  NM:i:1  MD:Z:48 YT:Z:UU
SRR7175373.2091997  0   ref|NC_001133|  1   37  1S49M   *   0   0   CCCACACCACACCCACACACCCACACACCACACCACACACCACACCACAC  ...GA<AGAGGA<AAGGAGGGGA<AAAGGAGGIIGGAG<AGGIGAGGIGG  AS:i:98 XS:i:67 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:49 YT:Z:UU
SRR7175373.2270050  0   ref|NC_001133|  1   11  9S41M   *   0   0   ACACCACACCCACACCACACCCACACACCCACACACCACACCACACACCA  GGGGGIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIG  AS:i:82 XS:i:78 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:41 YT:Z:UU
SRR7175373.7025590  0   ref|NC_001133|  1   14  8S42M   *   0   0   CACCACACCCACACCACACCCACACACCCACACACCACACCACACACCAC  GAGGAGG<<AG.<.<<<<G<<<GGGGGGGGGGGGGGGGAGGAGAAAAAG#  AS:i:84 XS:i:78 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:42 YT:Z:UU
SRR7175373.8104795  0   ref|NC_001133|  1   21  4S46M   *   0   0   ACACCCACACCACACCCACACACCCACACACCACACCACACACCACACCA  GGGGGIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII  AS:i:92 XS:i:68 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:46 YT:Z:UU
SRR7175373.8877369  0   ref|NC_001133|  1   1   4S31M1I14M  *   0   0   ACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACC  AGA.AGGGGIAGAGGGGGGGGGGG.<.<GAGAGG.<.G<AAGGGG#####  AS:i:82 XS:i:82 XN:i:0  XM:i:0  XO:i:1  XG:i:1  NM:i:1  MD:Z:45 YT:Z:UU
SRR7175373.9602499  0   ref|NC_001133|  1   14  5S45M   *   0   0   CACACCCACACCACACCCACACACCCACACACCACACCACACACCACACC  GGAGAGGGIIIIGIIIIIIIIIIIIIIIIIIIIIIIIGIIIIIIIIIIII  AS:i:90 XS:i:82 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:45 YT:Z:UU
SRR7175373.9864638  0   ref|NC_001133|  1   14  7S43M   *   0   0   ACCACACCCACACCACACCCACACACCCACACACCACACCACACACCACA  GGAAGIGIIIIIGGGIIIIIGIGGIIGGIIIGGGGIIIGIIGGGGIGIGG  AS:i:86 XS:i:80 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:43 YT:Z:UU
SRR7175373.11488836 0   ref|NC_001133|  1   11  2S31M1I16M  *   0   0   ACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCAC  AAAAGAAAGGGGAGGGGIIII.AAAAGGGGGGIG.GGGI<GG<.<AAAGA  AS:i:86 XS:i:82 XN:i:0  XM:i:0  XO:i:1  XG:i:1  NM:i:1  MD:Z:47 YT:Z:UU
```
</details>
<br />

<a id="fix-chromosome-names-in-bams"></a>
### Fix chromosome names in bams
<a id="get-jvarkit-bamrenamechr-up-and-running"></a>
#### Get `jvarkit` `bamrenamechr` up and running
<a id="code-4"></a>
##### Code
<details>
<summary><i>Code: Get jvarkit bamrenamechr up and running</i></summary>

```bash
#!/bin/bash

#  Still in coverage_env
#+
#+ Use jvarkit from Pierre Lindenbaum: github.com/lindenb/jvarkit
#+ Pre-compiled jar here: uncloud.univ-nantes.fr/index.php/s/4sL77oWR2BFzSBH
#+
#+ In particular, jvarkit utility bamrenamechr

#NOTE Below chunk is not necessary
# module load picard/2.25.1-Java-11
#
# f_genome="${d_genome}/fasta/S288C_reference_sequence_R64-3-1_20210421.fa"  # ., "${f_genome}"
# f_dict="${f_genome/.fa/.dict}"
# java -jar "${EBROOTPICARD}/picard.jar" CreateSequenceDictionary \
#     --REFERENCE "${f_genome}" \
#     --OUTPUT "${f_dict}"
#
# ., "${f_dict}"
# cat "${f_dict}"
#
# module purge

module load Java/1.8.0_181

f_jar="${HOME}/tsukiyamalab/Kris/2023_rDNA/software/jvarkit/jvarkit.jar"  # ., ${f_jar}
java -jar "${f_jar}" bamrenamechr --help

pwd
if [[ ! -f test-jvarkit.bam ]]; then
    cp Brn1_log_rep1_ChIP.bam test-jvarkit.bam
    samtools index -@ "${SLURM_CPUS_ON_NODE}" test-jvarkit.bam
    # samtools idxstats test-jvarkit.bam \
    #     | grep -v -F '*' \
    #     | awk '{printf("%s\tchr%s\n",$1,$1);}' \
    #         >  test-jvarkit.txt
fi

#  Create a custom dictionary (tab-separated) mapping current chromosome names
#+ to desired chromosome names
cat test-jvarkit.txt

java -jar "${f_jar}" bamrenamechr --help
java -jar "${f_jar}" bamrenamechr \
    -f test-jvarkit.txt \
    test-jvarkit.bam \
        > test-jvarkit.renamed.bam

samtools view test-jvarkit.renamed.bam | head  # Looks good
samtools view test-jvarkit.renamed.bam | tail  # Looks good

if [[ -f test-jvarkit.renamed.bam ]]; then
    cp test-jvarkit.txt dictionary_chr-names.txt
    rm test-jvarkit.*
fi
```
</details>
<br />

<a id="printed-3"></a>
##### Printed
<details>
<summary><i>Printed: Get jvarkit bamrenamechr up and running</i></summary>

```txt
❯ module load picard/2.25.1-Java-11
To execute picard run: java -jar $EBROOTPICARD/picard.jar


❯ f_dict="${f_genome/.fa/.dict}"


❯ java -jar "${EBROOTPICARD}/picard.jar" CreateSequenceDictionary \
>     --REFERENCE "${f_genome}" \
>     --OUTPUT "${f_dict}"
11:30:23.947 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/app/software/picard/2.25.1-Java-11/picard.jar!/com/intel/gkl/native/libgkl_compression.so
[Thu Apr 06 11:30:23 PDT 2023] CreateSequenceDictionary --OUTPUT /home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/fasta/S288C_reference_sequence_R64-3-1_20210421.dict --REFERENCE /home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/fasta/S288C_reference_sequence_R64-3-1_20210421.fa --TRUNCATE_NAMES_AT_WHITESPACE true --NUM_SEQUENCES 2147483647 --VERBOSITY INFO --QUIET false --VALIDATION_STRINGENCY STRICT --COMPRESSION_LEVEL 5 --MAX_RECORDS_IN_RAM 500000 --CREATE_INDEX false --CREATE_MD5_FILE false --GA4GH_CLIENT_SECRETS client_secrets.json --help false --version false --showHidden false --USE_JDK_DEFLATER false --USE_JDK_INFLATER false
[Thu Apr 06 11:30:23 PDT 2023] Executing as kalavatt@gizmoj10 on Linux 4.15.0-192-generic amd64; OpenJDK 64-Bit Server VM 11.0.2+9; Deflater: Intel; Inflater: Intel; Provider GCS is not available; Picard version: Version:2.25.1
[Thu Apr 06 11:30:24 PDT 2023] picard.sam.CreateSequenceDictionary done. Elapsed time: 0.01 minutes.
Runtime.totalMemory()=2147483648


❯ ., "${f_dict}"
-rw-rw---- 1 kalavatt 3.5K Apr  6 11:30 /home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/fasta/S288C_reference_sequence_R64-3-1_20210421.dict


❯ cat "${f_dict}"
@HD     VN:1.6
@SQ     SN:ref|NC_001133|       LN:230218       M5:6681ac2f62509cfc220d78751b8dc524     UR:file:/home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/fasta/S288C_reference_sequence_R64-3-1_20210421.fa
@SQ     SN:ref|NC_001134|       LN:813184       M5:97a317c689cbdd7e92a5c159acd290d2     UR:file:/home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/fasta/S288C_reference_sequence_R64-3-1_20210421.fa
@SQ     SN:ref|NC_001135|       LN:316620       M5:54f4a74aa6392d9e19b82c38aa8ab345     UR:file:/home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/fasta/S288C_reference_sequence_R64-3-1_20210421.fa
@SQ     SN:ref|NC_001136|       LN:1531933      M5:74180788027e20df3de53dcb2367d9e3     UR:file:/home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/fasta/S288C_reference_sequence_R64-3-1_20210421.fa
@SQ     SN:ref|NC_001137|       LN:576874       M5:d2787193198c8d260f58f2097f9e1e39     UR:file:/home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/fasta/S288C_reference_sequence_R64-3-1_20210421.fa
@SQ     SN:ref|NC_001138|       LN:270161       M5:b7ebc601f9a7df2e1ec5863deeae88a3     UR:file:/home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/fasta/S288C_reference_sequence_R64-3-1_20210421.fa
@SQ     SN:ref|NC_001139|       LN:1090940      M5:a308c7ebf0b67c4926bc190dc4ba8ed8     UR:file:/home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/fasta/S288C_reference_sequence_R64-3-1_20210421.fa
@SQ     SN:ref|NC_001140|       LN:562643       M5:f66a4f8eef89fc3c3a393fe0210169f1     UR:file:/home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/fasta/S288C_reference_sequence_R64-3-1_20210421.fa
@SQ     SN:ref|NC_001141|       LN:439888       M5:4eae53ae7b2029b7e1075461c3eb9aac     UR:file:/home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/fasta/S288C_reference_sequence_R64-3-1_20210421.fa
@SQ     SN:ref|NC_001142|       LN:745751       M5:6757b8c7d9cca2c56401e2484cf5e2fb     UR:file:/home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/fasta/S288C_reference_sequence_R64-3-1_20210421.fa
@SQ     SN:ref|NC_001143|       LN:666816       M5:e72df2471be793f8aa06850348a896fa     UR:file:/home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/fasta/S288C_reference_sequence_R64-3-1_20210421.fa
@SQ     SN:ref|NC_001144|       LN:1078177      M5:77945d734ab92ad527d8920c9d78ac1c     UR:file:/home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/fasta/S288C_reference_sequence_R64-3-1_20210421.fa
@SQ     SN:ref|NC_001145|       LN:924431       M5:073f9ff1c599c1a6867de2c7e4355394     UR:file:/home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/fasta/S288C_reference_sequence_R64-3-1_20210421.fa
@SQ     SN:ref|NC_001146|       LN:784333       M5:188bca5110182a786cd42686ec6882c6     UR:file:/home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/fasta/S288C_reference_sequence_R64-3-1_20210421.fa
@SQ     SN:ref|NC_001147|       LN:1091291      M5:7e02090a38f05459102d1a9a83703534     UR:file:/home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/fasta/S288C_reference_sequence_R64-3-1_20210421.fa
@SQ     SN:ref|NC_001148|       LN:948066       M5:232475e9a61a5e07f9cb2df4a2dad757     UR:file:/home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/fasta/S288C_reference_sequence_R64-3-1_20210421.fa
@SQ     SN:ref|NC_001224|       LN:85779        M5:71c39cf065b8d574f636b654c274cf1b     UR:file:/home/kalavatt/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/fasta/S288C_reference_sequence_R64-3-1_20210421.fa


❯ module purge


❯ module load Java/1.8.0_181


❯ java -jar "${f_jar}" bamrenamechr --help
Usage: bamrenamechr [options] Files
  Options:
    --bamcompression
      Compression Level. 0: no compression. 9: max compression;
      Default: 5
    --dict
      Use this new dictionary A SAM Sequence dictionary source: it can be a
      *.dict file, a fasta file indexed with 'picard
      CreateSequenceDictionary', or any hts file containing a dictionary (VCF,
      BAM, CRAM, intervals...)
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -i, --ignore
      If the tool cannot convert a contig, skip the read
      Default: false
    -f, --mapping, -m
      load a custom name mapping. Format (chrom-source\tchrom-dest\n)+
    -o, --out
      Output file. Optional . Default: stdout
    -R, --reference
      For Reading CRAM. Indexed fasta Reference file. This file must be
      indexed with samtools faidx and with picard CreateSequenceDictionary
    --samoutputformat
      Sam output format.
      Default: SAM
      Possible Values: [BAM, SAM, CRAM]
    --version
      print version and exit


❯ pwd
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams


❯ if [[ ! -f test-jvarkit.bam ]]; then
>     cp Brn1_log_rep1_ChIP.bam test-jvarkit.bam
>     samtools index -@ ${SLURM_CPUS_ON_NODE} test-jvarkit.bam
>     samtools idxstats test-jvarkit.bam \
>         | grep -v -F '*' \
>         | awk '{printf("%s\tchr%s\n",$1,$1);}' \
>             >  test-jvarkit.txt
> fi
'Brn1_log_rep1_ChIP.bam' -> 'test-jvarkit.bam'


❯ #  Create a custom dictionary mapping current chromosome names to desired
❯ #+ chromosome names
❯ cat test-jvarkit.txt
ref|NC_001133|  I
ref|NC_001134|  II
ref|NC_001135|  III
ref|NC_001136|  IV
ref|NC_001137|  V
ref|NC_001138|  VI
ref|NC_001139|  VII
ref|NC_001140|  VIII
ref|NC_001141|  IX
ref|NC_001142|  X
ref|NC_001143|  XI
ref|NC_001144|  XII
ref|NC_001145|  XIII
ref|NC_001146|  XIV
ref|NC_001147|  XV
ref|NC_001148|  XVI
ref|NC_001224|  Mito


❯ java -jar "${f_jar}" bamrenamechr \
>     -f test-jvarkit.txt \
>     test-jvarkit.bam \
>         > test-jvarkit.renamed.bam
[WARN][ConvertBamChromosomes]num ignored read:0


❯ zcat ~/tsukiyamalab/Kris/genomes/S288C_reference_genome_R64-3-1_20210421/S288C_reference_sequence_R64-3-1_20210421.fsa.gz | grep "^>"
>ref|NC_001133| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=I]
>ref|NC_001134| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=II]
>ref|NC_001135| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=III]
>ref|NC_001136| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=IV]
>ref|NC_001137| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=V]
>ref|NC_001138| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=VI]
>ref|NC_001139| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=VII]
>ref|NC_001140| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=VIII]
>ref|NC_001141| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=IX]
>ref|NC_001142| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=X]
>ref|NC_001143| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XI]
>ref|NC_001144| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XII]
>ref|NC_001145| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XIII]
>ref|NC_001146| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XIV]
>ref|NC_001147| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XV]
>ref|NC_001148| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XVI]
>ref|NC_001224| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [location=mitochondrion] [top=circular]
```
</details>
<br />

<a id="run-jvarkit-bamrenamechr-on-problem-bams"></a>
#### Run `jvarkit` `bamrenamechr` on problem bams
<a id="code-5"></a>
##### Code
<details>
<summary><i>Code: Run jvarkit bamrenamechr on problem bams</i></summary>

```bash
#!/bin/bash

#  Still in coverage_env, and still have module Java/1.8.0_181 loaded
f_jar="${HOME}/tsukiyamalab/Kris/2023_rDNA/software/jvarkit/jvarkit.jar"  # ., ${f_jar}
p_data="${HOME}/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams"
unset bams
typeset -a bams=(
    "${p_data}/Brn1_Q_rep1_ChIP.bam"
    "${p_data}/Brn1_Q_rep1_input.bam"
    "${p_data}/Brn1_Q_rep2_ChIP.bam"
    "${p_data}/Brn1_Q_rep2_input.bam"
    "${p_data}/Brn1_Q_rep3_ChIP.bam"
    "${p_data}/Brn1_Q_rep3_input.bam"
    "${p_data}/Brn1_Q_all_input.bam"
    "${p_data}/Brn1_Q_all_ChIP.bam"
    "${p_data}/Brn1_log_rep1_ChIP.bam"
    "${p_data}/Brn1_log_rep1_input.bam"
    "${p_data}/Brn1_log_rep2_ChIP.bam"
    "${p_data}/Brn1_log_rep2_input.bam"
    "${p_data}/Brn1_log_rep3_ChIP.bam"
    "${p_data}/Brn1_log_rep3_input.bam"
    "${p_data}/Brn1_log_all_ChIP.bam"
    "${p_data}/Brn1_log_all_input.bam"
)
# for i in "${bams[@]}"; do echo "${i}"; done
# for i in "${bams[@]}"; do ., "${i}"; done

for i in "${bams[@]}"; do
    java -jar "${f_jar}" bamrenamechr \
        -f dictionary_chr-names.txt \
        "${i}" \
            > "${i%.bam}.renamed.bam"

    if [[ -f "${i%.bam}.renamed.bam" ]]; then
        if [[ -s "${i%.bam}.renamed.bam" ]]; then
            mv -f "${i%.bam}.renamed.bam" "${i}"
        fi
    fi
done

#NOTE 1/4 Did not set '--samoutputformat bam' when calling bamrenamechr, so the
#NOTE 2/4 outfiles are actually sams; thus, they need to have their extensions
#NOTE 3/4 changed, then I need to converted them to bams; finally, I need to
#NOTE 4/4 index the files
rm *.bai
rename 's/.bam/.sam/g' *.bam
for i in "${bams[@]}"; do
    samtools view -@ "${SLURM_CPUS_ON_NODE}" -S -b "${i/.bam/.sam}" > "${i}"
done
```
</details>
<br />

<a id="printed-4"></a>
##### Printed
<details>
<summary><i>Code: Run jvarkit bamrenamechr on problem bams</i></summary>

```txt
❯ for i in "${bams[@]}"; do
>     java -jar "${f_jar}" bamrenamechr \
>         -f dictionary_chr-names.txt \
>         "${i}" \
>             > "${i%.bam}.renamed.bam"
> 
>     if [[ -f "${i%.bam}.renamed.bam" ]]; then
>         if [[ -s "${i%.bam}.renamed.bam" ]]; then
>             mv -f "${i%.bam}.renamed.bam" "${i}"
>         fi
>     fi
> done
[WARN][ConvertBamChromosomes]num ignored read:0
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep1_ChIP.renamed.bam' -> '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep1_ChIP.bam'

[WARN][ConvertBamChromosomes]num ignored read:0
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep1_input.renamed.bam' -> '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep1_input.bam'

[WARN][ConvertBamChromosomes]num ignored read:0
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep2_ChIP.renamed.bam' -> '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep2_ChIP.bam'

[WARN][ConvertBamChromosomes]num ignored read:0
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep2_input.renamed.bam' -> '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep2_input.bam'

[WARN][ConvertBamChromosomes]num ignored read:0
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep3_ChIP.renamed.bam' -> '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep3_ChIP.bam'

[WARN][ConvertBamChromosomes]num ignored read:0
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep3_input.renamed.bam' -> '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep3_input.bam'

[WARN][ConvertBamChromosomes]num ignored read:0
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_all_input.renamed.bam' -> '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_all_input.bam'

[WARN][ConvertBamChromosomes]num ignored read:0
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_all_ChIP.renamed.bam' -> '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_all_ChIP.bam'

[WARN][ConvertBamChromosomes]num ignored read:0
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep1_ChIP.renamed.bam' -> '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep1_ChIP.bam'

[WARN][ConvertBamChromosomes]num ignored read:0
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep1_input.renamed.bam' -> '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep1_input.bam'

[WARN][ConvertBamChromosomes]num ignored read:0
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep2_ChIP.renamed.bam' -> '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep2_ChIP.bam'

[WARN][ConvertBamChromosomes]num ignored read:0
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep2_input.renamed.bam' -> '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep2_input.bam'

[WARN][ConvertBamChromosomes]num ignored read:0
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep3_ChIP.renamed.bam' -> '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep3_ChIP.bam'

[WARN][ConvertBamChromosomes]num ignored read:0
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep3_input.renamed.bam' -> '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep3_input.bam'

[WARN][ConvertBamChromosomes]num ignored read:0
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_all_ChIP.renamed.bam' -> '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_all_ChIP.bam'

[WARN][ConvertBamChromosomes]num ignored read:0
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_all_input.renamed.bam' -> '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_all_input.bam'


❯ rm *.bai
```
</details>
<br />

<a id="correct-chromosome-names-in-s288c_reference_sequence_r64-3-1_20210421fa"></a>
## Correct chromosome names in `S288C_reference_sequence_R64-3-1_20210421.fa`
<a id="code-6"></a>
### Code
<details>
<summary><i>Code: S288C_reference_sequence_R64-3-1_20210421.fa</i></summary>

```bash
#!/bin/bash

cd "${HOME}/genomes/S288C_reference_genome_R64-3-1_20210421/fasta"

cat S288C_reference_sequence_R64-3-1_20210421.fa \
    sed 

cp \
    S288C_reference_sequence_R64-3-1_20210421.fa \
    S288C_reference_sequence_R64-3-1_20210421.initial.fa

cat S288C_reference_sequence_R64-3-1_20210421.initial.fa \
    | sed "s:.*NC_001133.*:>I:g" \
    | sed "s:.*NC_001134.*:>II:g" \
    | sed "s:.*NC_001135.*:>III:g" \
    | sed "s:.*NC_001136.*:>IV:g" \
    | sed "s:.*NC_001137.*:>V:g" \
    | sed "s:.*NC_001138.*:>VI:g" \
    | sed "s:.*NC_001139.*:>VII:g" \
    | sed "s:.*NC_001140.*:>VIII:g" \
    | sed "s:.*NC_001141.*:>IX:g" \
    | sed "s:.*NC_001142.*:>X:g" \
    | sed "s:.*NC_001143.*:>XI:g" \
    | sed "s:.*NC_001144.*:>XII:g" \
    | sed "s:.*NC_001145.*:>XIII:g" \
    | sed "s:.*NC_001146.*:>XIV:g" \
    | sed "s:.*NC_001147.*:>XV:g" \
    | sed "s:.*NC_001148.*:>XVI:g" \
    | sed "s:.*NC_001224.*:>Mito:g" \
        > tmp

grep "^>" S288C_reference_sequence_R64-3-1_20210421.initial.fa
grep "^>" tmp

mv -f tmp S288C_reference_sequence_R64-3-1_20210421.fa
```
</details>
<br />

<a id="printed-5"></a>
### Printed
<details>
<summary><i>Printed: S288C_reference_sequence_R64-3-1_20210421.fa</i></summary>

```txt
❯ cp \
>     S288C_reference_sequence_R64-3-1_20210421.fa \
>     S288C_reference_sequence_R64-3-1_20210421.initial.fa
'S288C_reference_sequence_R64-3-1_20210421.fa' -> 'S288C_reference_sequence_R64-3-1_20210421.initial.fa'


❯ grep "^>" S288C_reference_sequence_R64-3-1_20210421.initial.fa
>ref|NC_001133| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=I]
>ref|NC_001134| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=II]
>ref|NC_001135| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=III]
>ref|NC_001136| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=IV]
>ref|NC_001137| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=V]
>ref|NC_001138| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=VI]
>ref|NC_001139| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=VII]
>ref|NC_001140| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=VIII]
>ref|NC_001141| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=IX]
>ref|NC_001142| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=X]
>ref|NC_001143| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XI]
>ref|NC_001144| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XII]
>ref|NC_001145| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XIII]
>ref|NC_001146| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XIV]
>ref|NC_001147| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XV]
>ref|NC_001148| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XVI]
>ref|NC_001224| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [location=mitochondrion] [top=circular]


❯ grep "^>" tmp
>I
>II
>III
>IV
>V
>VI
>VII
>VIII
>IX
>X
>XI
>XII
>XIII
>XIV
>XV
>XVI
>Mito


❯ mv -f tmp S288C_reference_sequence_R64-3-1_20210421.fa
renamed 'tmp' -> 'S288C_reference_sequence_R64-3-1_20210421.fa'


❯ for i in "${bams[@]}"; do
>     samtools view -@ "${SLURM_CPUS_ON_NODE}" -S -b "${i/.bam/.sam}" > "${i}"
> done
```
</details>
<br />

<details>
<summary><i>Code: </i></summary>

```bash
#!/bin/bash

threads="${SLURM_CPUS_ON_NODE}"
input="${p_data}/Brn1_Q_rep1_input.bam"
IP="${p_data}/Brn1_Q_rep1_ChIP.bam"

samtools index -@ "${threads}" "${IP}"
samtools index -@ "${threads}" "${input}"

bamCompare \
    -b1 "${IP}" \
    -b2 "${input}" \
    -o "${IP%%.IP.sort.bam}.input-normalized.bw" \
    --binSize 10 \
    --scaleFactorsMethod None \
    --normalizeUsing BPM \
    --numberOfProcessors "${threads}" \
        > >(tee -a "${IP%%.IP.sort.bam}.stdout.txt") \
        2> >(tee -a "${IP%%.IP.sort.bam}.stderr.txt")


```
</details>
<br />
