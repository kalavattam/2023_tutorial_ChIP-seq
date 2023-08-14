
`#work_align_tally_data-Brn1_genome-cat_printed.md`
<br />
<br />

<details>
<summary><b><font size="+2"><i>Table of contents</i></font></b></summary>
<!-- MarkdownTOC -->

1. [Get situated](#get-situated)
	1. [Printed](#printed)
1. [Align the datasets](#align-the-datasets)
	1. [Run `bowtie2` alignment, etc.](#run-bowtie2-alignment-etc)
		1. [Get situated, set up variables and arrays](#get-situated-set-up-variables-and-arrays)
			1. [Printed](#printed-1)
		1. [Trim any adapter sequences present in the reads](#trim-any-adapter-sequences-present-in-the-reads)
			1. [Printed](#printed-2)
		1. [Align, sort, and index the sample datasets](#align-sort-and-index-the-sample-datasets)
			1. [`SmartMap`'s `Bowtie2` parameters](#smartmaps-bowtie2-parameters)
			1. [Align \(etc.\) `atria`-trimmed `fastq`s](#align-etc-atria-trimmed-fastqs)
				1. [Printed](#printed-3)
			1. [Align \(etc.\) un-trimmed `fastq`s](#align-etc-un-trimmed-fastqs)
				1. [Printed](#printed-4)
	1. [Examine flags in bam outfiles](#examine-flags-in-bam-outfiles)
		1. [Initialize necessary functions](#initialize-necessary-functions)
			1. [Printed](#printed-5)
		1. [Initialize an array of bams](#initialize-an-array-of-bams)
			1. [Printed](#printed-6)
		1. [Check on flag information in bams](#check-on-flag-information-in-bams)
			1. [Printed](#printed-7)
1. [Tally/calculate alignments](#tallycalculate-alignments)
	1. [Tally/calculate alignments](#tallycalculate-alignments-1)
		1. [Initialize functions for doing floating point arithmetic, etc.](#initialize-functions-for-doing-floating-point-arithmetic-etc)
			1. [Printed](#printed-8)
		1. [Get situated, then initialize arrays, variables, etc.](#get-situated-then-initialize-arrays-variables-etc)
			1. [Printed](#printed-9)
		1. [Generate tab-separated table of alignment tallies/calculations](#generate-tab-separated-table-of-alignment-talliescalculations)
			1. [Printed](#printed-10)
		1. [Calculate CC/SS-styled scaling factors](#calculate-ccss-styled-scaling-factors)
	1. [On calculating scaling factors](#on-calculating-scaling-factors)
		1. [Email from Christine](#email-from-christine)
		1. [Notes on using Excel to calculate scaling factors](#notes-on-using-excel-to-calculate-scaling-factors)

<!-- /MarkdownTOC -->
</details>
<br />
<br />

<a id="get-situated"></a>
## Get situated
<a id="printed"></a>
### Printed
<details>
<summary><i>Printed</i></summary>

```txt
❯ #UPDATE
❯ #  2023-0715-0716: fastqs are now compressed, and symlinked/renamed files are
❯ #+ now stored in subdirectory sym/, e.g.,


❯ p_sym="${p_data}/sym"


❯ p_data="${HOME}/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802"


❯ cd "${p_data}" || echo "cd'ing failed; check on this..."


❯ ls -lhaFG "${p_sym}"
total 1.3M
drwxrws--- 2 kalavatt 2.0K Jul 15 09:18 ./
drwxrws--- 5 kalavatt 1.8K Jul 16 06:28 ../
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_log_WT_Brn1_rep1.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175367.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_log_WT_Brn1_rep2.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175369.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_log_WT_Brn1_rep3.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175371.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_log_WT_Brn1_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175373.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_log_WT_Rpb3_rep1.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175389.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_log_WT_Rpb3_rep2.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175391.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_log_WT_Rpb3_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175393.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_Q_msn2_Brn1_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR8200048.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_Q_SMC4-off_Brn1_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR8200046.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_Q_SMC4-off_Rpb3_rep1.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175383.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_Q_SMC4-off_Rpb3_rep2.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175385.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_Q_SMC4-off_Rpb3_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175387.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_Q_SMC4-off_Ser2-phos_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR8200052.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_Q_WT_Brn1_rep1.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175375.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_Q_WT_Brn1_rep2.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175377.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_Q_WT_Brn1_rep3.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175379.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_Q_WT_Brn1_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175381.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_Q_WT_Msn2_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR8200050.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_Q_WT_Rpb3_rep1.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175395.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_Q_WT_Rpb3_rep2.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175397.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_Q_WT_Rpb3_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175399.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_Q_WT_Ser2-phos_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR8200054.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_log_WT_Brn1_rep1.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175368.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_log_WT_Brn1_rep2.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175370.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_log_WT_Brn1_rep3.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175372.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_log_WT_Brn1_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175374.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_log_WT_Rpb3_rep1.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175390.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_log_WT_Rpb3_rep2.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175392.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_log_WT_Rpb3_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175394.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_Q_msn2_Brn1_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR8200049.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_Q_SMC4-off_Brn1_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR8200047.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_Q_SMC4-off_Rpb3_rep1.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175384.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_Q_SMC4-off_Rpb3_rep2.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175386.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_Q_SMC4-off_Rpb3_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175388.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_Q_SMC4-off_Ser2-phos_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR8200053.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_Q_WT_Brn1_rep1.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175376.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_Q_WT_Brn1_rep2.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175378.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_Q_WT_Brn1_rep3.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175380.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_Q_WT_Brn1_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175382.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_Q_WT_Msn2_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR8200051.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_Q_WT_Rpb3_rep1.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175396.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_Q_WT_Rpb3_rep2.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175398.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_Q_WT_Rpb3_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175400.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_Q_WT_Ser2-phos_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR8200055.fastq.gz


❯ #NOTE
❯ #  For more details on the symlinking/renaming process, see
❯ #+ 2023_rDNA/results/2023-0228_work_fastqs_download/work_download-data.md
```
</details>
<br />
<br />

<a id="align-the-datasets"></a>
## Align the datasets
<a id="run-bowtie2-alignment-etc"></a>
### Run `bowtie2` alignment, etc.
<a id="get-situated-set-up-variables-and-arrays"></a>
#### Get situated, set up variables and arrays
<a id="printed-1"></a>
##### Printed
<details>
<summary><i>Printed</i></summary>

```txt
❯ grabnode
--------------------------------------------------------------
     ** Please review options as the order has changed **

 - memory is requested first
 - CPU selections will be adjusted based on memory request
--------------------------------------------------------------

How much memory (GB) will you need?
Enter a value from 1 to 683 [20]: 20
How many cores (CPUs) would you like to grab on the node?
enter a value between 1 and 36 [36]: 8
Please enter the max number of days you would like to grab this node: [1-7] 1
Do you need a GPU ? [y/N]N

You have requested 8 CPUs on this node/server for 1 days or until you type exit.

WARNING: If you exit this shell before your jobs are finished, your jobs
on this node/server will be terminated. Please use sbatch for larger jobs.

Shared PI folders can be found in: /fh/fast, /fh/scratch and /fh/secure.

Requesting Queue: campus-new cores: 8 memory: 20 gpu: NONE
srun: job 24831181 queued and waiting for resources
srun: job 24831181 has been allocated resources


❯ ml SAMtools/1.16.1-GCC-11.2.0 Bowtie2/2.4.4-GCC-11.2.0


❯ which bowtie2
/app/software/Bowtie2/2.4.4-GCC-11.2.0/bin/bowtie2


❯ #  Set up work directory, directories for processed data


❯ d_work="${HOME}/tsukiyamalab/Kris/2023_rDNA/results"  # cd "${d_work}"


❯ [[ ! -d "${d_work}/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out" ]] && \
>    mkdir -p "${d_work}/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out"
mkdir: created directory '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria'
mkdir: created directory '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out'


❯ [[ ! -d "${d_work}/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out" ]] && \
>    mkdir -p "${d_work}/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out"
mkdir: created directory '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams'
mkdir: created directory '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out'


❯ #  Go to work directory


❯ cd "${d_work}" && pwd
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results


❯ #  Initialize variables needed for alignment, etc.


❯ threads="${SLURM_CPUS_ON_NODE}"  # echo "${threads}"


❯ d_genome="${HOME}/tsukiyamalab/Kris/genomes/combined_SC_SP"  # ., "${d_genome}"


❯ f_genome="${d_genome}/fasta/combined_SC_SP.fa"  # ., "${f_genome}"


❯ f_indices="${d_genome}/bowtie2/$(basename "${f_genome}" .fa)"  # ., "${f_indices}"*


❯ err_out="${d_work}/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out"  # ls -lhaFG "${err_out}"


❯ d_bams="$(dirname ${err_out})"  # ls -lhaFG "${d_bams}"


❯ p_data="${HOME}/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym"  # ., "${p_data}"


❯ typeset -a fastqs=(
>     "${p_data}/Ch_log_WT_Brn1_rep1.fastq.gz"
>     "${p_data}/Ch_log_WT_Brn1_rep2.fastq.gz"
>     "${p_data}/Ch_log_WT_Brn1_rep3.fastq.gz"
>     "${p_data}/Ch_log_WT_Brn1_repM.fastq.gz"
>     "${p_data}/Ch_Q_WT_Brn1_rep1.fastq.gz"
>     "${p_data}/Ch_Q_WT_Brn1_rep2.fastq.gz"
>     "${p_data}/Ch_Q_WT_Brn1_rep3.fastq.gz"
>     "${p_data}/Ch_Q_WT_Brn1_repM.fastq.gz"
>     "${p_data}/in_log_WT_Brn1_rep1.fastq.gz"
>     "${p_data}/in_log_WT_Brn1_rep2.fastq.gz"
>     "${p_data}/in_log_WT_Brn1_rep3.fastq.gz"
>     "${p_data}/in_log_WT_Brn1_repM.fastq.gz"
>     "${p_data}/in_Q_WT_Brn1_rep1.fastq.gz"
>     "${p_data}/in_Q_WT_Brn1_rep2.fastq.gz"
>     "${p_data}/in_Q_WT_Brn1_rep3.fastq.gz"
>     "${p_data}/in_Q_WT_Brn1_repM.fastq.gz"
> )


❯ run_check=TRUE


❯ [[ "${run_check}" == TRUE ]] &&
>     {
>         echo '### echo "${threads}" ###'
>         echo "${threads}"
>         echo ""
> 
>         echo '### ls -lhaFG "${d_genome}" ###'
>         ls -lhaFG "${d_genome}"
>         echo ""
> 
>         echo '### ls -lhaFG "${f_genome}" ###'
>         ls -lhaFG "${f_genome}"
>         echo ""
> 
>         echo '### ls -lhaFG "${f_indices}"* ###'
>         ls -lhaFG "${f_indices}"*
>         echo ""
> 
>         echo '### ls -lhaFG "${err_out}" ###'
>         ls -lhaFG "${err_out}"
>         echo ""
> 
>         echo '### ls -lhaFG "${d_bams}" ###'
>         ls -lhaFG "${d_bams}"
>         echo ""
> 
>         echo '### ls -lhaFG "${p_data}" ###'
>         ls -lhaFG "${p_data}"
>         echo ""
> 
>         echo '### for i in "${fastqs[@]}"; do ls -lhaFG "${i}"; done ###'
>         for i in "${fastqs[@]}"; do ls -lhaFG "${i}"; done
>
>         echo '### for i in "${atria[@]}"; do echo "${i}"; done ###'
>         for i in "${atria[@]}"; do echo "${i}"; done
>     }
### echo "${threads}" ###
8

### ls -lhaFG "${scratch} ###
total 512
drwxrws---   2 root   0 Jul 13 10:01 ./
drwxr-xr-x 250 root 248 Jul  5 07:00 ../

### ls -lhaFG "${d_genome}" ###
total 184K
drwxrwx---  5 kalavatt  70 May 29 11:47 ./
drwxrwx--- 14 kalavatt 537 Jul 11 15:19 ../
drwxrwx---  2 kalavatt 322 May 29 11:59 bowtie2/
drwxrwx---  2 kalavatt 242 May 29 11:48 fasta/
drwxrwx---  2 kalavatt 174 May 29 11:46 gff3/

### ls -lhaFG "${f_genome}" ###
-rw-rw---- 1 kalavatt 25M May 29 11:47 /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/fasta/combined_SC_SP.fa

### ls -lhaFG "${f_indices}"* ###
-rw-rw---- 1 kalavatt  12M May 29 11:59 /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP.1.bt2
-rw-rw---- 1 kalavatt 6.0M May 29 11:59 /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP.2.bt2
-rw-rw---- 1 kalavatt  269 May 29 11:59 /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP.3.bt2
-rw-rw---- 1 kalavatt 6.0M May 29 11:59 /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP.4.bt2
-rw-rw---- 1 kalavatt  12M May 29 12:00 /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP.rev.1.bt2
-rw-rw---- 1 kalavatt 6.0M May 29 12:00 /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP.rev.2.bt2
-rw-rw---- 1 kalavatt   23 May 29 11:59 /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP.stderr.txt
-rw-rw---- 1 kalavatt  13K May 29 12:00 /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP.stdout.txt

### ls -lhaFG "${err_out}" ###
total 72K
drwxrws--- 2 kalavatt  0 Jul 16 06:43 ./
drwxrws--- 3 kalavatt 25 Jul 16 06:43 ../

### ls -lhaFG "${d_bams}" ###
total 112K
drwxrws--- 3 kalavatt  25 Jul 16 06:43 ./
drwxrws--- 4 kalavatt 363 Jul 16 07:11 ../
drwxrws--- 2 kalavatt   0 Jul 16 06:43 err_out/

### ls -lhaFG "${p_data}" ###
total 1.3M
drwxrws--- 2 kalavatt 2.0K Jul 15 09:18 ./
drwxrws--- 5 kalavatt 1.8K Jul 16 06:28 ../
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_log_WT_Brn1_rep1.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175367.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_log_WT_Brn1_rep2.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175369.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_log_WT_Brn1_rep3.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175371.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_log_WT_Brn1_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175373.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_log_WT_Rpb3_rep1.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175389.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_log_WT_Rpb3_rep2.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175391.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_log_WT_Rpb3_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175393.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_Q_msn2_Brn1_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR8200048.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_Q_SMC4-off_Brn1_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR8200046.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_Q_SMC4-off_Rpb3_rep1.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175383.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_Q_SMC4-off_Rpb3_rep2.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175385.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_Q_SMC4-off_Rpb3_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175387.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_Q_SMC4-off_Ser2-phos_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR8200052.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_Q_WT_Brn1_rep1.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175375.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_Q_WT_Brn1_rep2.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175377.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_Q_WT_Brn1_rep3.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175379.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_Q_WT_Brn1_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175381.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_Q_WT_Msn2_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR8200050.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_Q_WT_Rpb3_rep1.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175395.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_Q_WT_Rpb3_rep2.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175397.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_Q_WT_Rpb3_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175399.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 Ch_Q_WT_Ser2-phos_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR8200054.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_log_WT_Brn1_rep1.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175368.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_log_WT_Brn1_rep2.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175370.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_log_WT_Brn1_rep3.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175372.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_log_WT_Brn1_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175374.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_log_WT_Rpb3_rep1.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175390.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_log_WT_Rpb3_rep2.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175392.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_log_WT_Rpb3_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175394.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_Q_msn2_Brn1_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR8200049.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_Q_SMC4-off_Brn1_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR8200047.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_Q_SMC4-off_Rpb3_rep1.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175384.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_Q_SMC4-off_Rpb3_rep2.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175386.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_Q_SMC4-off_Rpb3_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175388.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_Q_SMC4-off_Ser2-phos_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR8200053.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_Q_WT_Brn1_rep1.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175376.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_Q_WT_Brn1_rep2.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175378.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_Q_WT_Brn1_rep3.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175380.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_Q_WT_Brn1_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175382.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_Q_WT_Msn2_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR8200051.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_Q_WT_Rpb3_rep1.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175396.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_Q_WT_Rpb3_rep2.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175398.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_Q_WT_Rpb3_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175400.fastq.gz
lrwxrwxrwx 1 kalavatt   83 Jul 15 09:18 in_Q_WT_Ser2-phos_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR8200055.fastq.gz

### for i in "${fastqs[@]}"; do ls -lhaFG "${i}"; done ###
lrwxrwxrwx 1 kalavatt 83 Jul 15 09:18 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/Ch_log_WT_Brn1_rep1.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175367.fastq.gz
lrwxrwxrwx 1 kalavatt 83 Jul 15 09:18 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/Ch_log_WT_Brn1_rep2.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175369.fastq.gz
lrwxrwxrwx 1 kalavatt 83 Jul 15 09:18 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/Ch_log_WT_Brn1_rep3.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175371.fastq.gz
lrwxrwxrwx 1 kalavatt 83 Jul 15 09:18 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/Ch_log_WT_Brn1_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175373.fastq.gz
lrwxrwxrwx 1 kalavatt 83 Jul 15 09:18 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/Ch_Q_WT_Brn1_rep1.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175375.fastq.gz
lrwxrwxrwx 1 kalavatt 83 Jul 15 09:18 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/Ch_Q_WT_Brn1_rep2.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175377.fastq.gz
lrwxrwxrwx 1 kalavatt 83 Jul 15 09:18 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/Ch_Q_WT_Brn1_rep3.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175379.fastq.gz
lrwxrwxrwx 1 kalavatt 83 Jul 15 09:18 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/Ch_Q_WT_Brn1_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175381.fastq.gz
lrwxrwxrwx 1 kalavatt 83 Jul 15 09:18 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/in_log_WT_Brn1_rep1.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175368.fastq.gz
lrwxrwxrwx 1 kalavatt 83 Jul 15 09:18 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/in_log_WT_Brn1_rep2.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175370.fastq.gz
lrwxrwxrwx 1 kalavatt 83 Jul 15 09:18 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/in_log_WT_Brn1_rep3.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175372.fastq.gz
lrwxrwxrwx 1 kalavatt 83 Jul 15 09:18 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/in_log_WT_Brn1_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175374.fastq.gz
lrwxrwxrwx 1 kalavatt 83 Jul 15 09:18 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/in_Q_WT_Brn1_rep1.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175376.fastq.gz
lrwxrwxrwx 1 kalavatt 83 Jul 15 09:18 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/in_Q_WT_Brn1_rep2.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175378.fastq.gz
lrwxrwxrwx 1 kalavatt 83 Jul 15 09:18 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/in_Q_WT_Brn1_rep3.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175380.fastq.gz
lrwxrwxrwx 1 kalavatt 83 Jul 15 09:18 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/in_Q_WT_Brn1_repM.fastq.gz -> /home/kalavatt/tsukiyamalab/kalavatt/2023_rDNA/data/PRJNA471802/SRR7175382.fastq.gz

### for i in "${atria[@]}"; do echo "${i}"; done ###
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_log_WT_Brn1_rep1.atria.fastq.gz
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_log_WT_Brn1_rep2.atria.fastq.gz
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_log_WT_Brn1_rep3.atria.fastq.gz
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_log_WT_Brn1_repM.atria.fastq.gz
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_Q_WT_Brn1_rep1.atria.fastq.gz
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_Q_WT_Brn1_rep2.atria.fastq.gz
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_Q_WT_Brn1_rep3.atria.fastq.gz
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_Q_WT_Brn1_repM.atria.fastq.gz
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_log_WT_Brn1_rep1.atria.fastq.gz
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_log_WT_Brn1_rep2.atria.fastq.gz
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_log_WT_Brn1_rep3.atria.fastq.gz
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_log_WT_Brn1_repM.atria.fastq.gz
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_Q_WT_Brn1_rep1.atria.fastq.gz
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_Q_WT_Brn1_rep2.atria.fastq.gz
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_Q_WT_Brn1_rep3.atria.fastq.gz
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_Q_WT_Brn1_repM.atria.fastq.gz
```
</details>
<br />

<a id="trim-any-adapter-sequences-present-in-the-reads"></a>
#### Trim any adapter sequences present in the reads
<a id="printed-2"></a>
##### Printed
<details>
<summary><i>Printed</i></summary>

```txt
❯ #  Run print tests to check that the commands are correct/reasonable


❯ print_test=TRUE


❯ [[ "${print_test}" == TRUE ]] &&
>     {
>         for i in "${fastqs[@]}"; do
>             echo """
>             \"\${HOME}/tsukiyamalab/kalavatt/2023_rDNA/src/Atria/app-3.2.2/bin/atria\" \\
>                 -t \"${threads}\" \\
>                 -r \"${i}\" \\
>                 -o \"${p_atria}\" \\
>                 --length-range 35:500
>             """
>         done
>     }

            "${HOME}/tsukiyamalab/kalavatt/2023_rDNA/src/Atria/app-3.2.2/bin/atria" \
                -t "8" \
                -r "/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/Ch_log_WT_Brn1_rep1.fastq.gz" \
                -o "/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria" \
                --length-range 35:500


            "${HOME}/tsukiyamalab/kalavatt/2023_rDNA/src/Atria/app-3.2.2/bin/atria" \
                -t "8" \
                -r "/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/Ch_log_WT_Brn1_rep2.fastq.gz" \
                -o "/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria" \
                --length-range 35:500


            "${HOME}/tsukiyamalab/kalavatt/2023_rDNA/src/Atria/app-3.2.2/bin/atria" \
                -t "8" \
                -r "/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/Ch_log_WT_Brn1_rep3.fastq.gz" \
                -o "/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria" \
                --length-range 35:500


            "${HOME}/tsukiyamalab/kalavatt/2023_rDNA/src/Atria/app-3.2.2/bin/atria" \
                -t "8" \
                -r "/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/Ch_log_WT_Brn1_repM.fastq.gz" \
                -o "/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria" \
                --length-range 35:500


            "${HOME}/tsukiyamalab/kalavatt/2023_rDNA/src/Atria/app-3.2.2/bin/atria" \
                -t "8" \
                -r "/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/Ch_Q_WT_Brn1_rep1.fastq.gz" \
                -o "/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria" \
                --length-range 35:500


            "${HOME}/tsukiyamalab/kalavatt/2023_rDNA/src/Atria/app-3.2.2/bin/atria" \
                -t "8" \
                -r "/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/Ch_Q_WT_Brn1_rep2.fastq.gz" \
                -o "/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria" \
                --length-range 35:500


            "${HOME}/tsukiyamalab/kalavatt/2023_rDNA/src/Atria/app-3.2.2/bin/atria" \
                -t "8" \
                -r "/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/Ch_Q_WT_Brn1_rep3.fastq.gz" \
                -o "/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria" \
                --length-range 35:500


            "${HOME}/tsukiyamalab/kalavatt/2023_rDNA/src/Atria/app-3.2.2/bin/atria" \
                -t "8" \
                -r "/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/Ch_Q_WT_Brn1_repM.fastq.gz" \
                -o "/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria" \
                --length-range 35:500


            "${HOME}/tsukiyamalab/kalavatt/2023_rDNA/src/Atria/app-3.2.2/bin/atria" \
                -t "8" \
                -r "/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/in_log_WT_Brn1_rep1.fastq.gz" \
                -o "/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria" \
                --length-range 35:500


            "${HOME}/tsukiyamalab/kalavatt/2023_rDNA/src/Atria/app-3.2.2/bin/atria" \
                -t "8" \
                -r "/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/in_log_WT_Brn1_rep2.fastq.gz" \
                -o "/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria" \
                --length-range 35:500


            "${HOME}/tsukiyamalab/kalavatt/2023_rDNA/src/Atria/app-3.2.2/bin/atria" \
                -t "8" \
                -r "/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/in_log_WT_Brn1_rep3.fastq.gz" \
                -o "/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria" \
                --length-range 35:500


            "${HOME}/tsukiyamalab/kalavatt/2023_rDNA/src/Atria/app-3.2.2/bin/atria" \
                -t "8" \
                -r "/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/in_log_WT_Brn1_repM.fastq.gz" \
                -o "/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria" \
                --length-range 35:500


            "${HOME}/tsukiyamalab/kalavatt/2023_rDNA/src/Atria/app-3.2.2/bin/atria" \
                -t "8" \
                -r "/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/in_Q_WT_Brn1_rep1.fastq.gz" \
                -o "/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria" \
                --length-range 35:500


            "${HOME}/tsukiyamalab/kalavatt/2023_rDNA/src/Atria/app-3.2.2/bin/atria" \
                -t "8" \
                -r "/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/in_Q_WT_Brn1_rep2.fastq.gz" \
                -o "/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria" \
                --length-range 35:500


            "${HOME}/tsukiyamalab/kalavatt/2023_rDNA/src/Atria/app-3.2.2/bin/atria" \
                -t "8" \
                -r "/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/in_Q_WT_Brn1_rep3.fastq.gz" \
                -o "/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria" \
                --length-range 35:500


            "${HOME}/tsukiyamalab/kalavatt/2023_rDNA/src/Atria/app-3.2.2/bin/atria" \
                -t "8" \
                -r "/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/in_Q_WT_Brn1_repM.fastq.gz" \
                -o "/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria" \
                --length-range 35:500


❯ #  If things look good, then run the calls to atria (the program for trimming)


❯ run=TRUE


❯ [[ "${run}" == "${run}" ]] &&
>     {
>         #  Initialize environment for atria's dependencies
>         [[ "${CONDA_DEFAULT_ENV}" != "atria_env" ]] && source activate atria_env
> 
>         for i in "${fastqs[@]}"; do
>             "${HOME}/tsukiyamalab/kalavatt/2023_rDNA/src/Atria/app-3.2.2/bin/atria" \
>                 -t "${threads}" \
>                 -r "${i}" \
>                 -o "${p_atria}" \
>                 --length-range 35:500
>         done
>     }
┌ Info: ATRIA VERSIONS
│   atria = "v3.2.2"
└   julia = "v1.8.5"
┌ Info: ATRIA ARGUMENTS
└   command = `-t 8 -r /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/Ch_log_WT_Brn1_rep1.fastq.gz -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria --length-range 35:500`
┌ Info: ATRIA OUTPUT FILES
└   read = "/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_log_WT_Brn1_rep1.atria.fastq.gz"
┌ Info: ATRIA TRIMMERS AND FILTERS
│   tail_polyG_trimming = false
│   tail_polyT_trimming = false
│   tail_polyA_trimming = false
│   tail_polyC_trimming = false
│   adapter_trimming = true
│   consensus_calling = false
│   hard_clip_3_end = false
│   hard_clip_5_end = false
│   quality_trimming = true
│   tail_N_trimming = true
│   max_N_filtering = true
│   length_filtering = true
└   complexity_filtering = false
[ Info: Cycle 1: read 388237/388237 pairs; wrote 368390/368390; (copied 0/0)
[ Info: Cycle 2: read 385683/773920 pairs; wrote 365573/733963; (copied 0/0)
[ Info: Cycle 3: read 382096/1156016 pairs; wrote 360165/1094128; (copied 0/0)
[ Info: Cycle 4: read 377016/1533032 pairs; wrote 350980/1445108; (copied 0/0)
[ Info: Cycle 5: read 377016/1910048 pairs; wrote 355116/1800224; (copied 0/0)
... (see logs)
┌ Info: ATRIA COMPLETE
└   read = "/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_Q_WT_Brn1_repM.atria.fastq.gz"


❯ #  Move atria logs to atria/err_out/


❯ [[ $(find "${p_atria}" -type f \( -name "*.log" -o -name "*.json" \) | wc -l) -gt 0 ]] && \
>     mv ${p_atria}/*.{log,json} "${p_atria}/err_out"
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_log_WT_Brn1_rep1.atria.log' -> '/home/kalavatt/tsukiyamalab/Kris/2023_rD
NA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out/Ch_log_WT_Brn1_rep1.atria.log'
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_log_WT_Brn1_rep2.atria.log' -> '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out/Ch_log_WT_Brn1_rep2.atria.log'
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_log_WT_Brn1_rep3.atria.log' -> '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out/Ch_log_WT_Brn1_rep3.atria.log'
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_log_WT_Brn1_repM.atria.log' -> '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out/Ch_log_WT_Brn1_repM.atria.log'
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_Q_WT_Brn1_rep1.atria.log' -> '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out/Ch_Q_WT_Brn1_rep1.atria.log'
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_Q_WT_Brn1_rep2.atria.log' -> '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out/Ch_Q_WT_Brn1_rep2.atria.log'
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_Q_WT_Brn1_rep3.atria.log' -> '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out/Ch_Q_WT_Brn1_rep3.atria.log'
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_Q_WT_Brn1_repM.atria.log' -> '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out/Ch_Q_WT_Brn1_repM.atria.log'
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_log_WT_Brn1_rep1.atria.log' -> '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out/in_log_WT_Brn1_rep1.atria.log'
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_log_WT_Brn1_rep2.atria.log' -> '/home/kalavatt/tsukiyamalab/Kris/2023_rD
NA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out/in_log_WT_Brn1_rep2.atria.log'
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_log_WT_Brn1_rep3.atria.log' -> '/home/kalavatt/tsukiyamalab/Kris/2023_rD
NA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out/in_log_WT_Brn1_rep3.atria.log'
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_log_WT_Brn1_repM.atria.log' -> '/home/kalavatt/tsukiyamalab/Kris/2023_rD
NA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out/in_log_WT_Brn1_repM.atria.log'
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_Q_WT_Brn1_rep1.atria.log' -> '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA
/results/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out/in_Q_WT_Brn1_rep1.atria.log'
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_Q_WT_Brn1_rep2.atria.log' -> '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA
/results/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out/in_Q_WT_Brn1_rep2.atria.log'
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_Q_WT_Brn1_rep3.atria.log' -> '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA
/results/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out/in_Q_WT_Brn1_rep3.atria.log'
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_Q_WT_Brn1_repM.atria.log' -> '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA
/results/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out/in_Q_WT_Brn1_repM.atria.log'
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_log_WT_Brn1_rep1.atria.log.json' -> '/home/kalavatt/tsukiyamalab/Kris/20
23_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out/Ch_log_WT_Brn1_rep1.atria.log.json'
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_log_WT_Brn1_rep2.atria.log.json' -> '/home/kalavatt/tsukiyamalab/Kris/20
23_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out/Ch_log_WT_Brn1_rep2.atria.log.json'
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_log_WT_Brn1_rep3.atria.log.json' -> '/home/kalavatt/tsukiyamalab/Kris/20
23_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out/Ch_log_WT_Brn1_rep3.atria.log.json'
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_log_WT_Brn1_repM.atria.log.json' -> '/home/kalavatt/tsukiyamalab/Kris/20
23_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out/Ch_log_WT_Brn1_repM.atria.log.json'
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_Q_WT_Brn1_rep1.atria.log.json' -> '/home/kalavatt/tsukiyamalab/Kris/2023
_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out/Ch_Q_WT_Brn1_rep1.atria.log.json'
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_Q_WT_Brn1_rep2.atria.log.json' -> '/home/kalavatt/tsukiyamalab/Kris/2023
_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out/Ch_Q_WT_Brn1_rep2.atria.log.json'
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_Q_WT_Brn1_rep3.atria.log.json' -> '/home/kalavatt/tsukiyamalab/Kris/2023
_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out/Ch_Q_WT_Brn1_rep3.atria.log.json'
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_Q_WT_Brn1_repM.atria.log.json' -> '/home/kalavatt/tsukiyamalab/Kris/2023
_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out/Ch_Q_WT_Brn1_repM.atria.log.json'
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_log_WT_Brn1_rep1.atria.log.json' -> '/home/kalavatt/tsukiyamalab/Kris/20
23_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out/in_log_WT_Brn1_rep1.atria.log.json'
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_log_WT_Brn1_rep2.atria.log.json' -> '/home/kalavatt/tsukiyamalab/Kris/20
23_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out/in_log_WT_Brn1_rep2.atria.log.json'
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_log_WT_Brn1_rep3.atria.log.json' -> '/home/kalavatt/tsukiyamalab/Kris/20
23_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out/in_log_WT_Brn1_rep3.atria.log.json'
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_log_WT_Brn1_repM.atria.log.json' -> '/home/kalavatt/tsukiyamalab/Kris/20
23_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out/in_log_WT_Brn1_repM.atria.log.json'
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_Q_WT_Brn1_rep1.atria.log.json' -> '/home/kalavatt/tsukiyamalab/Kris/2023
_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out/in_Q_WT_Brn1_rep1.atria.log.json'
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_Q_WT_Brn1_rep2.atria.log.json' -> '/home/kalavatt/tsukiyamalab/Kris/2023
_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out/in_Q_WT_Brn1_rep2.atria.log.json'
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_Q_WT_Brn1_rep3.atria.log.json' -> '/home/kalavatt/tsukiyamalab/Kris/2023
_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out/in_Q_WT_Brn1_rep3.atria.log.json'
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_Q_WT_Brn1_repM.atria.log.json' -> '/home/kalavatt/tsukiyamalab/Kris/2023
_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/err_out/in_Q_WT_Brn1_repM.atria.log.json'


❯ #  Check that "${atria[@]}" array elements exist/are recognized by the system


❯ run_check=TRUE


❯ [[ "${run_check}" == TRUE ]] &&
>     {
>         echo '### for i in "${atria[@]}"; do ls -lhaFG "${i}"; done ###'
>         for i in "${atria[@]}"; do ls -lhaFG "${i}"; done
>     }
### for i in "${atria[@]}"; do ls -lhaFG "${i}"; done ###
-rw-rw---- 1 kalavatt 198M Jul 16 15:40 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_log_WT_Brn1_rep1.atria.fastq.gz
-rw-rw---- 1 kalavatt 959M Jul 16 15:41 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_log_WT_Brn1_rep2.atria.fastq.gz
-rw-rw---- 1 kalavatt 268M Jul 16 15:41 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_log_WT_Brn1_rep3.atria.fastq.gz
-rw-rw---- 1 kalavatt 1.4G Jul 16 15:42 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_log_WT_Brn1_repM.atria.fastq.gz
-rw-rw---- 1 kalavatt 240M Jul 16 15:43 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_Q_WT_Brn1_rep1.atria.fastq.gz
-rw-rw---- 1 kalavatt 215M Jul 16 15:43 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_Q_WT_Brn1_rep2.atria.fastq.gz
-rw-rw---- 1 kalavatt 303M Jul 16 15:43 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_Q_WT_Brn1_rep3.atria.fastq.gz
-rw-rw---- 1 kalavatt 760M Jul 16 15:44 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_Q_WT_Brn1_repM.atria.fastq.gz
-rw-rw---- 1 kalavatt 491M Jul 16 15:45 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_log_WT_Brn1_rep1.atria.fastq.gz
-rw-rw---- 1 kalavatt 391M Jul 16 15:45 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_log_WT_Brn1_rep2.atria.fastq.gz
-rw-rw---- 1 kalavatt 284M Jul 16 15:46 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_log_WT_Brn1_rep3.atria.fastq.gz
-rw-rw---- 1 kalavatt 1.2G Jul 16 15:47 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_log_WT_Brn1_repM.atria.fastq.gz
-rw-rw---- 1 kalavatt 364M Jul 16 15:47 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_Q_WT_Brn1_rep1.atria.fastq.gz
-rw-rw---- 1 kalavatt 421M Jul 16 15:48 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_Q_WT_Brn1_rep2.atria.fastq.gz
-rw-rw---- 1 kalavatt 358M Jul 16 15:48 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_Q_WT_Brn1_rep3.atria.fastq.gz
-rw-rw---- 1 kalavatt 1.2G Jul 16 15:49 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_Q_WT_Brn1_repM.atria.fastq.gz
```
</details>
<br />

<a id="align-sort-and-index-the-sample-datasets"></a>
#### Align, sort, and index the sample datasets
<a id="smartmaps-bowtie2-parameters"></a>
##### `SmartMap`'s `Bowtie2` parameters
<a id="align-etc-atria-trimmed-fastqs"></a>
##### Align (etc.) `atria`-trimmed `fastq`s
<a id="printed-3"></a>
###### Printed
<details>
<summary><i>Printed</i></summary>

```txt
❯ #  Run print tests to check that the commands are correct/reasonable


❯ print_test=TRUE


❯ [[ "${print_test}" == TRUE ]] &&
>     {
>         for i in "${atria[@]}"; do
>             # i="${atria[0]}"  # ., "${i}"
>             in_fastq="${i}"  # ., "${in_fastq}"
>             out_bam="${d_bams}/$(basename "${in_fastq}" .fastq.gz).bam"  # echo "${out_bam}"
> 
>             echo """
>             #  ---------------------------------------------------------
>             #  Align fastqs trimmed with atria; sort resulting bams
>             {
>                 bowtie2 \\
>                     -p ${threads} \\
>                     -x ${f_indices} \\
>                     --end-to-end \\
>                     --very-fast \\
>                     -U ${i} \\
>                         | samtools sort \\
>                             -@ ${threads} \\
>                             -T ${scratch} \\
>                             -O bam \\
>                             -o ${out_bam}
>             } \\
>                  > >(tee -a ${err_out}/$(basename "${out_bam}" .bam).bowtie2_samtools-sort.stdout.txt) \\
>                 2> >(tee -a ${err_out}/$(basename "${out_bam}" .bam).bowtie2_samtools-sort.stderr.txt)
> 
>             #  Index the sorted bams
>             if [[ -f ${out_bam} ]]; then
>                 samtools index \\
>                     -@ ${threads} \\
>                     ${out_bam} \\
>                          > >(tee -a ${err_out}/$(basename ${out_bam} .bam).samtools-index.stdout.txt) \\
>                         2> >(tee -a ${err_out}/$(basename ${out_bam} .bam).samtools-index.stderr.txt)
>             fi
>             """
>         done
            #  ---------------------------------------------------------
            #  Align fastqs trimmed with atria; sort resulting bams
            {
                bowtie2 \
                    -p 8 \
                    -x /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP \
                    --end-to-end \
                    --very-fast \
                    -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_log_WT_Brn1_rep1.atria.fastq.gz \
                        | samtools sort \
                            -@ 8 \
                            -T /fh/scratch/delete30/tsukiyama_t \
                            -O bam \
                            -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_log_WT_Brn1_rep1.atria.bam
            } \
                 > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_log_WT_Brn1_rep1.atria.bowtie2_samtools-sort.stdout.txt) \
                2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_log_WT_Brn1_rep1.atria.bowtie2_samtools-sort.stderr.txt)

            #  Index the sorted bams
            if [[ -f /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_log_WT_Brn1_rep1.atria.bam ]]; then
                samtools index \
                    -@ 8 \
                    /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_log_WT_Brn1_rep1.atria.bam \
                         > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_log_WT_Brn1_rep1.atria.samtools-index.stdout.txt) \
                        2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_log_WT_Brn1_rep1.atria.samtools-index.stderr.txt)
            fi


            #  ---------------------------------------------------------
            #  Align fastqs trimmed with atria; sort resulting bams
            {
                bowtie2 \
                    -p 8 \
                    -x /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP \
                    --end-to-end \
                    --very-fast \
                    -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_log_WT_Brn1_rep2.atria.fastq.gz \
                        | samtools sort \
                            -@ 8 \
                            -T /fh/scratch/delete30/tsukiyama_t \
                            -O bam \
                            -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_log_WT_Brn1_rep2.atria.bam
            } \
                 > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_log_WT_Brn1_rep2.atria.bowtie2_samtools-sort.stdout.txt) \
                2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_log_WT_Brn1_rep2.atria.bowtie2_samtools-sort.stderr.txt)

            #  Index the sorted bams
            if [[ -f /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_log_WT_Brn1_rep2.atria.bam ]]; then
                samtools index \
                    -@ 8 \
                    /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_log_WT_Brn1_rep2.atria.bam \
                         > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_log_WT_Brn1_rep2.atria.samtools-index.stdout.txt) \
                        2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_log_WT_Brn1_rep2.atria.samtools-index.stderr.txt)
            fi


            #  ---------------------------------------------------------
            #  Align fastqs trimmed with atria; sort resulting bams
            {
                bowtie2 \
                    -p 8 \
                    -x /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP \
                    --end-to-end \
                    --very-fast \
                    -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_log_WT_Brn1_rep3.atria.fastq.gz \
                        | samtools sort \
                            -@ 8 \
                            -T /fh/scratch/delete30/tsukiyama_t \
                            -O bam \
                            -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_log_WT_Brn1_rep3.atria.bam
            } \
                 > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_log_WT_Brn1_rep3.atria.bowtie2_samtools-sort.stdout.txt) \
                2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_log_WT_Brn1_rep3.atria.bowtie2_samtools-sort.stderr.txt)

            #  Index the sorted bams
            if [[ -f /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_log_WT_Brn1_rep3.atria.bam ]]; then
                samtools index \
                    -@ 8 \
                    /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_log_WT_Brn1_rep3.atria.bam \
                         > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_log_WT_Brn1_rep3.atria.samtools-index.stdout.txt) \
                        2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_log_WT_Brn1_rep3.atria.samtools-index.stderr.txt)
            fi


            #  ---------------------------------------------------------
            #  Align fastqs trimmed with atria; sort resulting bams
            {
                bowtie2 \
                    -p 8 \
                    -x /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP \
                    --end-to-end \
                    --very-fast \
                    -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_log_WT_Brn1_repM.atria.fastq.gz \
                        | samtools sort \
                            -@ 8 \
                            -T /fh/scratch/delete30/tsukiyama_t \
                            -O bam \
                            -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_log_WT_Brn1_repM.atria.bam
            } \
                 > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_log_WT_Brn1_repM.atria.bowtie2_samtools-sort.stdout.txt) \
                2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_log_WT_Brn1_repM.atria.bowtie2_samtools-sort.stderr.txt)

            #  Index the sorted bams
            if [[ -f /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_log_WT_Brn1_repM.atria.bam ]]; then
                samtools index \
                    -@ 8 \
                    /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_log_WT_Brn1_repM.atria.bam \
                         > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_log_WT_Brn1_repM.atria.samtools-index.stdout.txt) \
                        2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_log_WT_Brn1_repM.atria.samtools-index.stderr.txt)
            fi


            #  ---------------------------------------------------------
            #  Align fastqs trimmed with atria; sort resulting bams
            {
                bowtie2 \
                    -p 8 \
                    -x /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP \
                    --end-to-end \
                    --very-fast \
                    -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_Q_WT_Brn1_rep1.atria.fastq.gz \
                        | samtools sort \
                            -@ 8 \
                            -T /fh/scratch/delete30/tsukiyama_t \
                            -O bam \
                            -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_Q_WT_Brn1_rep1.atria.bam
            } \
                 > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_Q_WT_Brn1_rep1.atria.bowtie2_samtools-sort.stdout.txt) \
                2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_Q_WT_Brn1_rep1.atria.bowtie2_samtools-sort.stderr.txt)

            #  Index the sorted bams
            if [[ -f /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_Q_WT_Brn1_rep1.atria.bam ]]; then
                samtools index \
                    -@ 8 \
                    /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_Q_WT_Brn1_rep1.atria.bam \
                         > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_Q_WT_Brn1_rep1.atria.samtools-index.stdout.txt) \
                        2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_Q_WT_Brn1_rep1.atria.samtools-index.stderr.txt)
            fi


            #  ---------------------------------------------------------
            #  Align fastqs trimmed with atria; sort resulting bams
            {
                bowtie2 \
                    -p 8 \
                    -x /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP \
                    --end-to-end \
                    --very-fast \
                    -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_Q_WT_Brn1_rep2.atria.fastq.gz \
                        | samtools sort \
                            -@ 8 \
                            -T /fh/scratch/delete30/tsukiyama_t \
                            -O bam \
                            -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_Q_WT_Brn1_rep2.atria.bam
            } \
                 > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_Q_WT_Brn1_rep2.atria.bowtie2_samtools-sort.stdout.txt) \
                2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_Q_WT_Brn1_rep2.atria.bowtie2_samtools-sort.stderr.txt)

            #  Index the sorted bams
            if [[ -f /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_Q_WT_Brn1_rep2.atria.bam ]]; then
                samtools index \
                    -@ 8 \
                    /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_Q_WT_Brn1_rep2.atria.bam \
                         > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_Q_WT_Brn1_rep2.atria.samtools-index.stdout.txt) \
                        2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_Q_WT_Brn1_rep2.atria.samtools-index.stderr.txt)
            fi


            #  ---------------------------------------------------------
            #  Align fastqs trimmed with atria; sort resulting bams
            {
                bowtie2 \
                    -p 8 \
                    -x /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP \
                    --end-to-end \
                    --very-fast \
                    -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_Q_WT_Brn1_rep3.atria.fastq.gz \
                        | samtools sort \
                            -@ 8 \
                            -T /fh/scratch/delete30/tsukiyama_t \
                            -O bam \
                            -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_Q_WT_Brn1_rep3.atria.bam
            } \
                 > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_Q_WT_Brn1_rep3.atria.bowtie2_samtools-sort.stdout.txt) \
                2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_Q_WT_Brn1_rep3.atria.bowtie2_samtools-sort.stderr.txt)

            #  Index the sorted bams
            if [[ -f /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_Q_WT_Brn1_rep3.atria.bam ]]; then
                samtools index \
                    -@ 8 \
                    /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_Q_WT_Brn1_rep3.atria.bam \
                         > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_Q_WT_Brn1_rep3.atria.samtools-index.stdout.txt) \
                        2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_Q_WT_Brn1_rep3.atria.samtools-index.stderr.txt)
            fi


            #  ---------------------------------------------------------
            #  Align fastqs trimmed with atria; sort resulting bams
            {
                bowtie2 \
                    -p 8 \
                    -x /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP \
                    --end-to-end \
                    --very-fast \
                    -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/Ch_Q_WT_Brn1_repM.atria.fastq.gz \
                        | samtools sort \
                            -@ 8 \
                            -T /fh/scratch/delete30/tsukiyama_t \
                            -O bam \
                            -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_Q_WT_Brn1_repM.atria.bam
            } \
                 > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_Q_WT_Brn1_repM.atria.bowtie2_samtools-sort.stdout.txt) \
                2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_Q_WT_Brn1_repM.atria.bowtie2_samtools-sort.stderr.txt)

            #  Index the sorted bams
            if [[ -f /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_Q_WT_Brn1_repM.atria.bam ]]; then
                samtools index \
                    -@ 8 \
                    /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_Q_WT_Brn1_repM.atria.bam \
                         > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_Q_WT_Brn1_repM.atria.samtools-index.stdout.txt) \
                        2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_Q_WT_Brn1_repM.atria.samtools-index.stderr.txt)
            fi


            #  ---------------------------------------------------------
            #  Align fastqs trimmed with atria; sort resulting bams
            {
                bowtie2 \
                    -p 8 \
                    -x /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP \
                    --end-to-end \
                    --very-fast \
                    -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_log_WT_Brn1_rep1.atria.fastq.gz \
                        | samtools sort \
                            -@ 8 \
                            -T /fh/scratch/delete30/tsukiyama_t \
                            -O bam \
                            -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_log_WT_Brn1_rep1.atria.bam
            } \
                 > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_log_WT_Brn1_rep1.atria.bowtie2_samtools-sort.stdout.txt) \
                2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_log_WT_Brn1_rep1.atria.bowtie2_samtools-sort.stderr.txt)

            #  Index the sorted bams
            if [[ -f /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_log_WT_Brn1_rep1.atria.bam ]]; then
                samtools index \
                    -@ 8 \
                    /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_log_WT_Brn1_rep1.atria.bam \
                         > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_log_WT_Brn1_rep1.atria.samtools-index.stdout.txt) \
                        2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_log_WT_Brn1_rep1.atria.samtools-index.stderr.txt)
            fi


            #  ---------------------------------------------------------
            #  Align fastqs trimmed with atria; sort resulting bams
            {
                bowtie2 \
                    -p 8 \
                    -x /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP \
                    --end-to-end \
                    --very-fast \
                    -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_log_WT_Brn1_rep2.atria.fastq.gz \
                        | samtools sort \
                            -@ 8 \
                            -T /fh/scratch/delete30/tsukiyama_t \
                            -O bam \
                            -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_log_WT_Brn1_rep2.atria.bam
            } \
                 > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_log_WT_Brn1_rep2.atria.bowtie2_samtools-sort.stdout.txt) \
                2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_log_WT_Brn1_rep2.atria.bowtie2_samtools-sort.stderr.txt)

            #  Index the sorted bams
            if [[ -f /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_log_WT_Brn1_rep2.atria.bam ]]; then
                samtools index \
                    -@ 8 \
                    /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_log_WT_Brn1_rep2.atria.bam \
                         > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_log_WT_Brn1_rep2.atria.samtools-index.stdout.txt) \
                        2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_log_WT_Brn1_rep2.atria.samtools-index.stderr.txt)
            fi


            #  ---------------------------------------------------------
            #  Align fastqs trimmed with atria; sort resulting bams
            {
                bowtie2 \
                    -p 8 \
                    -x /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP \
                    --end-to-end \
                    --very-fast \
                    -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_log_WT_Brn1_rep3.atria.fastq.gz \
                        | samtools sort \
                            -@ 8 \
                            -T /fh/scratch/delete30/tsukiyama_t \
                            -O bam \
                            -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_log_WT_Brn1_rep3.atria.bam
            } \
                 > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_log_WT_Brn1_rep3.atria.bowtie2_samtools-sort.stdout.txt) \
                2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_log_WT_Brn1_rep3.atria.bowtie2_samtools-sort.stderr.txt)

            #  Index the sorted bams
            if [[ -f /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_log_WT_Brn1_rep3.atria.bam ]]; then
                samtools index \
                    -@ 8 \
                    /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_log_WT_Brn1_rep3.atria.bam \
                         > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_log_WT_Brn1_rep3.atria.samtools-index.stdout.txt) \
                        2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_log_WT_Brn1_rep3.atria.samtools-index.stderr.txt)
            fi


            #  ---------------------------------------------------------
            #  Align fastqs trimmed with atria; sort resulting bams
            {
                bowtie2 \
                    -p 8 \
                    -x /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP \
                    --end-to-end \
                    --very-fast \
                    -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_log_WT_Brn1_repM.atria.fastq.gz \
                        | samtools sort \
                            -@ 8 \
                            -T /fh/scratch/delete30/tsukiyama_t \
                            -O bam \
                            -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_log_WT_Brn1_repM.atria.bam
            } \
                 > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_log_WT_Brn1_repM.atria.bowtie2_samtools-sort.stdout.txt) \
                2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_log_WT_Brn1_repM.atria.bowtie2_samtools-sort.stderr.txt)

            #  Index the sorted bams
            if [[ -f /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_log_WT_Brn1_repM.atria.bam ]]; then
                samtools index \
                    -@ 8 \
                    /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_log_WT_Brn1_repM.atria.bam \
                         > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_log_WT_Brn1_repM.atria.samtools-index.stdout.txt) \
                        2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_log_WT_Brn1_repM.atria.samtools-index.stderr.txt)
            fi


            #  ---------------------------------------------------------
            #  Align fastqs trimmed with atria; sort resulting bams
            {
                bowtie2 \
                    -p 8 \
                    -x /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP \
                    --end-to-end \
                    --very-fast \
                    -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_Q_WT_Brn1_rep1.atria.fastq.gz \
                        | samtools sort \
                            -@ 8 \
                            -T /fh/scratch/delete30/tsukiyama_t \
                            -O bam \
                            -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_Q_WT_Brn1_rep1.atria.bam
            } \
                 > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_Q_WT_Brn1_rep1.atria.bowtie2_samtools-sort.stdout.txt) \
                2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_Q_WT_Brn1_rep1.atria.bowtie2_samtools-sort.stderr.txt)

            #  Index the sorted bams
            if [[ -f /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_Q_WT_Brn1_rep1.atria.bam ]]; then
                samtools index \
                    -@ 8 \
                    /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_Q_WT_Brn1_rep1.atria.bam \
                         > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_Q_WT_Brn1_rep1.atria.samtools-index.stdout.txt) \
                        2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_Q_WT_Brn1_rep1.atria.samtools-index.stderr.txt)
            fi


            #  ---------------------------------------------------------
            #  Align fastqs trimmed with atria; sort resulting bams
            {
                bowtie2 \
                    -p 8 \
                    -x /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP \
                    --end-to-end \
                    --very-fast \
                    -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_Q_WT_Brn1_rep2.atria.fastq.gz \
                        | samtools sort \
                            -@ 8 \
                            -T /fh/scratch/delete30/tsukiyama_t \
                            -O bam \
                            -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_Q_WT_Brn1_rep2.atria.bam
            } \
                 > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_Q_WT_Brn1_rep2.atria.bowtie2_samtools-sort.stdout.txt) \
                2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_Q_WT_Brn1_rep2.atria.bowtie2_samtools-sort.stderr.txt)

            #  Index the sorted bams
            if [[ -f /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_Q_WT_Brn1_rep2.atria.bam ]]; then
                samtools index \
                    -@ 8 \
                    /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_Q_WT_Brn1_rep2.atria.bam \
                         > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_Q_WT_Brn1_rep2.atria.samtools-index.stdout.txt) \
                        2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_Q_WT_Brn1_rep2.atria.samtools-index.stderr.txt)
            fi


            #  ---------------------------------------------------------
            #  Align fastqs trimmed with atria; sort resulting bams
            {
                bowtie2 \
                    -p 8 \
                    -x /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP \
                    --end-to-end \
                    --very-fast \
                    -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_Q_WT_Brn1_rep3.atria.fastq.gz \
                        | samtools sort \
                            -@ 8 \
                            -T /fh/scratch/delete30/tsukiyama_t \
                            -O bam \
                            -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_Q_WT_Brn1_rep3.atria.bam
            } \
                 > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_Q_WT_Brn1_rep3.atria.bowtie2_samtools-sort.stdout.txt) \
                2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_Q_WT_Brn1_rep3.atria.bowtie2_samtools-sort.stderr.txt)

            #  Index the sorted bams
            if [[ -f /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_Q_WT_Brn1_rep3.atria.bam ]]; then
                samtools index \
                    -@ 8 \
                    /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_Q_WT_Brn1_rep3.atria.bam \
                         > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_Q_WT_Brn1_rep3.atria.samtools-index.stdout.txt) \
                        2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_Q_WT_Brn1_rep3.atria.samtools-index.stderr.txt)
            fi


            #  ---------------------------------------------------------
            #  Align fastqs trimmed with atria; sort resulting bams
            {
                bowtie2 \
                    -p 8 \
                    -x /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP \
                    --end-to-end \
                    --very-fast \
                    -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/atria/in_Q_WT_Brn1_repM.atria.fastq.gz \
                        | samtools sort \
                            -@ 8 \
                            -T /fh/scratch/delete30/tsukiyama_t \
                            -O bam \
                            -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_Q_WT_Brn1_repM.atria.bam
            } \
                 > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_Q_WT_Brn1_repM.atria.bowtie2_samtools-sort.stdout.txt) \
                2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_Q_WT_Brn1_repM.atria.bowtie2_samtools-sort.stderr.txt)

            #  Index the sorted bams
            if [[ -f /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_Q_WT_Brn1_repM.atria.bam ]]; then
                samtools index \
                    -@ 8 \
                    /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_Q_WT_Brn1_repM.atria.bam \
                         > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_Q_WT_Brn1_repM.atria.samtools-index.stdout.txt) \
                        2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_Q_WT_Brn1_repM.atria.samtools-index.stderr.txt)
            fi


❯ run=TRUE


❯ [[ "${run}" == TRUE ]] &&
>     {
>         for i in "${atria[@]}"; do
>             # i="${atria[0]}"  # echo "${i}"
>             in_fastq="${i}"  # ., "${in_fastq}"
>             out_bam="${d_bams}/$(basename "${in_fastq}" .fastq.gz).bam"  # echo "${out_bam}"
> 
>             echo "#### Aligning, sorting, indexing $(basename ${in_fastq}) ####"
>             #  Align fastqs trimmed with atria; sort resulting bams
>             {
>                 bowtie2 \
>                     -p "${threads}" \
>                     --end-to-end \
>                     --very-fast \
>                     -x "${f_indices}" \
>                     -U "${in_fastq}" \
>                         | samtools sort \
>                             -@ "${threads}" \
>                             -T "${scratch}" \
>                             -O "bam" \
>                             -o "${out_bam}"
>             } \
>                  > "${err_out}/$(basename "${out_bam}" .bam).bowtie2_samtools-sort.stdout.txt" \
>                 2> "${err_out}/$(basename "${out_bam}" .bam).bowtie2_samtools-sort.stderr.txt"
> 
>             #  Index the sorted bams
>             if [[ -f "${out_bam}" ]]; then
>                 samtools index \
>                     -@ "${threads}" \
>                     "${out_bam}" \
>                          > >(tee -a "${err_out}/$(basename "${out_bam}" .bam).samtools-index.stdout.txt") \
>                         2> >(tee -a "${err_out}/$(basename "${out_bam}" .bam).samtools-index.stderr.txt")
>             fi
> 
>             if [[ -f "${out_bam}" && -f "${out_bam}.bai" ]]; then
>                 echo "#DONE"
>             else
>                 echo "#PROBLEM"
>             fi
>             echo ""
>         done
>     }
#### Aligning, sorting, indexing Ch_log_WT_Brn1_rep1.atria.fastq.gz ####
#DONE

#### Aligning, sorting, indexing Ch_log_WT_Brn1_rep2.atria.fastq.gz ####
#DONE

#### Aligning, sorting, indexing Ch_log_WT_Brn1_rep3.atria.fastq.gz ####
#DONE

#### Aligning, sorting, indexing Ch_log_WT_Brn1_repM.atria.fastq.gz ####
#DONE

#### Aligning, sorting, indexing Ch_Q_WT_Brn1_rep1.atria.fastq.gz ####
#DONE

#### Aligning, sorting, indexing Ch_Q_WT_Brn1_rep2.atria.fastq.gz ####
#DONE

#### Aligning, sorting, indexing Ch_Q_WT_Brn1_rep3.atria.fastq.gz ####
#DONE

#### Aligning, sorting, indexing Ch_Q_WT_Brn1_repM.atria.fastq.gz ####
#DONE

#### Aligning, sorting, indexing in_log_WT_Brn1_rep1.atria.fastq.gz ####
#DONE

#### Aligning, sorting, indexing in_log_WT_Brn1_rep2.atria.fastq.gz ####
#DONE

#### Aligning, sorting, indexing in_log_WT_Brn1_rep3.atria.fastq.gz ####
#DONE

#### Aligning, sorting, indexing in_log_WT_Brn1_repM.atria.fastq.gz ####
#DONE

#### Aligning, sorting, indexing in_Q_WT_Brn1_rep1.atria.fastq.gz ####
#DONE

#### Aligning, sorting, indexing in_Q_WT_Brn1_rep2.atria.fastq.gz ####
#DONE

#### Aligning, sorting, indexing in_Q_WT_Brn1_rep3.atria.fastq.gz ####
#DONE

#### Aligning, sorting, indexing in_Q_WT_Brn1_repM.atria.fastq.gz ####
#DONE
```
</details>
<br />

<a id="align-etc-un-trimmed-fastqs"></a>
##### Align (etc.) un-trimmed `fastq`s
<a id="printed-4"></a>
###### Printed
<details>
<summary><i>Printed</i></summary>

```txt
❯ #  Align un-trimmed fastqs --------------------------------


❯ #  Run print tests to check that the commands are correct/reasonable


❯ print_test=TRUE


❯ [[ "${print_test}" == TRUE ]] &&
>     {
>         for i in "${fastqs[@]}"; do
>             # i="${fastqs[0]}"  # ., "${i}"
>             in_fastq="${i}"  # ., "${in_fastq}"
>             out_bam="${d_bams}/$(basename "${in_fastq}" .fastq.gz).bam"  # echo "${out_bam}"
> 
>             echo """
>             #  ---------------------------------------------------------
>             #  Align untrimmed fastqs; sort resulting bams
>             {
>                 bowtie2 \\
>                     -p ${threads} \\
>                     -x ${f_indices} \\
>                     --very-sensitive-local \\
>                     -U ${i} \\
>                         | samtools sort \\
>                             -@ ${threads} \\
>                             -T ${scratch} \\
>                             -O bam \\
>                             -o ${out_bam}
>             } \\
>                  > >(tee -a ${err_out}/$(basename "${out_bam}" .bam).bowtie2_samtools-sort.stdout.txt) \\
>                 2> >(tee -a ${err_out}/$(basename "${out_bam}" .bam).bowtie2_samtools-sort.stderr.txt)
> 
>             #  Index the sorted bams
>             if [[ -f ${out_bam} ]]; then
>                 samtools index \\
>                     -@ ${threads} \\
>                     ${out_bam} \\
>                          > >(tee -a ${err_out}/$(basename ${out_bam} .bam).samtools-index.stdout.txt) \\
>                         2> >(tee -a ${err_out}/$(basename ${out_bam} .bam).samtools-index.stderr.txt)
>             fi
>             """
>         done
>     }

            #  ---------------------------------------------------------
            #  Align untrimmed fastqs; sort resulting bams
            {
                bowtie2 \
                    -p 8 \
                    -x /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP \
                    --very-sensitive-local \
                    -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/Ch_log_WT_Brn1_rep1.fastq.gz \
                        | samtools sort \
                            -@ 8 \
                            -T /fh/scratch/delete30/tsukiyama_t \
                            -O bam \
                            -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_log_WT_Brn1_rep1.bam
            } \
                 > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_log_WT_Brn1_rep1.bowtie2_samtools-sort.stdout.txt) \
                2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_log_WT_Brn1_rep1.bowtie2_samtools-sort.stderr.txt)

            #  Index the sorted bams
            if [[ -f /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_log_WT_Brn1_rep1.bam ]]; then
                samtools index \
                    -@ 8 \
                    /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_log_WT_Brn1_rep1.bam \
                         > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_log_WT_Brn1_rep1.samtools-index.stdout.txt) \
                        2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_log_WT_Brn1_rep1.samtools-index.stderr.txt)
            fi


            #  ---------------------------------------------------------
            #  Align untrimmed fastqs; sort resulting bams
            {
                bowtie2 \
                    -p 8 \
                    -x /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP \
                    --very-sensitive-local \
                    -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/Ch_log_WT_Brn1_rep2.fastq.gz \
                        | samtools sort \
                            -@ 8 \
                            -T /fh/scratch/delete30/tsukiyama_t \
                            -O bam \
                            -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_log_WT_Brn1_rep2.bam
            } \
                 > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_log_WT_Brn1_rep2.bowtie2_samtools-sort.stdout.txt) \
                2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_log_WT_Brn1_rep2.bowtie2_samtools-sort.stderr.txt)

            #  Index the sorted bams
            if [[ -f /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_log_WT_Brn1_rep2.bam ]]; then
                samtools index \
                    -@ 8 \
                    /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_log_WT_Brn1_rep2.bam \
                         > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_log_WT_Brn1_rep2.samtools-index.stdout.txt) \
                        2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_log_WT_Brn1_rep2.samtools-index.stderr.txt)
            fi


            #  ---------------------------------------------------------
            #  Align untrimmed fastqs; sort resulting bams
            {
                bowtie2 \
                    -p 8 \
                    -x /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP \
                    --very-sensitive-local \
                    -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/Ch_log_WT_Brn1_rep3.fastq.gz \
                        | samtools sort \
                            -@ 8 \
                            -T /fh/scratch/delete30/tsukiyama_t \
                            -O bam \
                            -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_log_WT_Brn1_rep3.bam
            } \
                 > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_log_WT_Brn1_rep3.bowtie2_samtools-sort.stdout.txt) \
                2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_log_WT_Brn1_rep3.bowtie2_samtools-sort.stderr.txt)

            #  Index the sorted bams
            if [[ -f /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_log_WT_Brn1_rep3.bam ]]; then
                samtools index \
                    -@ 8 \
                    /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_log_WT_Brn1_rep3.bam \
                         > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_log_WT_Brn1_rep3.samtools-index.stdout.txt) \
                        2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_log_WT_Brn1_rep3.samtools-index.stderr.txt)
            fi


            #  ---------------------------------------------------------
            #  Align untrimmed fastqs; sort resulting bams
            {
                bowtie2 \
                    -p 8 \
                    -x /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP \
                    --very-sensitive-local \
                    -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/Ch_log_WT_Brn1_repM.fastq.gz \
                        | samtools sort \
                            -@ 8 \
                            -T /fh/scratch/delete30/tsukiyama_t \
                            -O bam \
                            -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_log_WT_Brn1_repM.bam
            } \
                 > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_log_WT_Brn1_repM.bowtie2_samtools-sort.stdout.txt) \
                2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_log_WT_Brn1_repM.bowtie2_samtools-sort.stderr.txt)

            #  Index the sorted bams
            if [[ -f /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_log_WT_Brn1_repM.bam ]]; then
                samtools index \
                    -@ 8 \
                    /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_log_WT_Brn1_repM.bam \
                         > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_log_WT_Brn1_repM.samtools-index.stdout.txt) \
                        2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_log_WT_Brn1_repM.samtools-index.stderr.txt)
            fi


            #  ---------------------------------------------------------
            #  Align untrimmed fastqs; sort resulting bams
            {
                bowtie2 \
                    -p 8 \
                    -x /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP \
                    --very-sensitive-local \
                    -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/Ch_Q_WT_Brn1_rep1.fastq.gz \
                        | samtools sort \
                            -@ 8 \
                            -T /fh/scratch/delete30/tsukiyama_t \
                            -O bam \
                            -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_Q_WT_Brn1_rep1.bam
            } \
                 > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_Q_WT_Brn1_rep1.bowtie2_samtools-sort.stdout.txt) \
                2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_Q_WT_Brn1_rep1.bowtie2_samtools-sort.stderr.txt)

            #  Index the sorted bams
            if [[ -f /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_Q_WT_Brn1_rep1.bam ]]; then
                samtools index \
                    -@ 8 \
                    /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_Q_WT_Brn1_rep1.bam \
                         > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_Q_WT_Brn1_rep1.samtools-index.stdout.txt) \
                        2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_Q_WT_Brn1_rep1.samtools-index.stderr.txt)
            fi


            #  ---------------------------------------------------------
            #  Align untrimmed fastqs; sort resulting bams
            {
                bowtie2 \
                    -p 8 \
                    -x /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP \
                    --very-sensitive-local \
                    -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/Ch_Q_WT_Brn1_rep2.fastq.gz \
                        | samtools sort \
                            -@ 8 \
                            -T /fh/scratch/delete30/tsukiyama_t \
                            -O bam \
                            -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_Q_WT_Brn1_rep2.bam
            } \
                 > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_Q_WT_Brn1_rep2.bowtie2_samtools-sort.stdout.txt) \
                2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_Q_WT_Brn1_rep2.bowtie2_samtools-sort.stderr.txt)

            #  Index the sorted bams
            if [[ -f /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_Q_WT_Brn1_rep2.bam ]]; then
                samtools index \
                    -@ 8 \
                    /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_Q_WT_Brn1_rep2.bam \
                         > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_Q_WT_Brn1_rep2.samtools-index.stdout.txt) \
                        2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_Q_WT_Brn1_rep2.samtools-index.stderr.txt)
            fi


            #  ---------------------------------------------------------
            #  Align untrimmed fastqs; sort resulting bams
            {
                bowtie2 \
                    -p 8 \
                    -x /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP \
                    --very-sensitive-local \
                    -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/Ch_Q_WT_Brn1_rep3.fastq.gz \
                        | samtools sort \
                            -@ 8 \
                            -T /fh/scratch/delete30/tsukiyama_t \
                            -O bam \
                            -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_Q_WT_Brn1_rep3.bam
            } \
                 > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_Q_WT_Brn1_rep3.bowtie2_samtools-sort.stdout.txt) \
                2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_Q_WT_Brn1_rep3.bowtie2_samtools-sort.stderr.txt)

            #  Index the sorted bams
            if [[ -f /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_Q_WT_Brn1_rep3.bam ]]; then
                samtools index \
                    -@ 8 \
                    /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_Q_WT_Brn1_rep3.bam \
                         > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_Q_WT_Brn1_rep3.samtools-index.stdout.txt) \
                        2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_Q_WT_Brn1_rep3.samtools-index.stderr.txt)
            fi


            #  ---------------------------------------------------------
            #  Align untrimmed fastqs; sort resulting bams
            {
                bowtie2 \
                    -p 8 \
                    -x /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP \
                    --very-sensitive-local \
                    -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/Ch_Q_WT_Brn1_repM.fastq.gz \
                        | samtools sort \
                            -@ 8 \
                            -T /fh/scratch/delete30/tsukiyama_t \
                            -O bam \
                            -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_Q_WT_Brn1_repM.bam
            } \
                 > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_Q_WT_Brn1_repM.bowtie2_samtools-sort.stdout.txt) \
                2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_Q_WT_Brn1_repM.bowtie2_samtools-sort.stderr.txt)

            #  Index the sorted bams
            if [[ -f /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_Q_WT_Brn1_repM.bam ]]; then
                samtools index \
                    -@ 8 \
                    /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_Q_WT_Brn1_repM.bam \
                         > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_Q_WT_Brn1_repM.samtools-index.stdout.txt) \
                        2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/Ch_Q_WT_Brn1_repM.samtools-index.stderr.txt)
            fi


            #  ---------------------------------------------------------
            #  Align untrimmed fastqs; sort resulting bams
            {
                bowtie2 \
                    -p 8 \
                    -x /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP \
                    --very-sensitive-local \
                    -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/in_log_WT_Brn1_rep1.fastq.gz \
                        | samtools sort \
                            -@ 8 \
                            -T /fh/scratch/delete30/tsukiyama_t \
                            -O bam \
                            -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_log_WT_Brn1_rep1.bam
            } \
                 > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_log_WT_Brn1_rep1.bowtie2_samtools-sort.stdout.txt) \
                2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_log_WT_Brn1_rep1.bowtie2_samtools-sort.stderr.txt)

            #  Index the sorted bams
            if [[ -f /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_log_WT_Brn1_rep1.bam ]]; then
                samtools index \
                    -@ 8 \
                    /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_log_WT_Brn1_rep1.bam \
                         > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_log_WT_Brn1_rep1.samtools-index.stdout.txt) \
                        2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_log_WT_Brn1_rep1.samtools-index.stderr.txt)
            fi


            #  ---------------------------------------------------------
            #  Align untrimmed fastqs; sort resulting bams
            {
                bowtie2 \
                    -p 8 \
                    -x /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP \
                    --very-sensitive-local \
                    -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/in_log_WT_Brn1_rep2.fastq.gz \
                        | samtools sort \
                            -@ 8 \
                            -T /fh/scratch/delete30/tsukiyama_t \
                            -O bam \
                            -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_log_WT_Brn1_rep2.bam
            } \
                 > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_log_WT_Brn1_rep2.bowtie2_samtools-sort.stdout.txt) \
                2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_log_WT_Brn1_rep2.bowtie2_samtools-sort.stderr.txt)

            #  Index the sorted bams
            if [[ -f /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_log_WT_Brn1_rep2.bam ]]; then
                samtools index \
                    -@ 8 \
                    /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_log_WT_Brn1_rep2.bam \
                         > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_log_WT_Brn1_rep2.samtools-index.stdout.txt) \
                        2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_log_WT_Brn1_rep2.samtools-index.stderr.txt)
            fi


            #  ---------------------------------------------------------
            #  Align untrimmed fastqs; sort resulting bams
            {
                bowtie2 \
                    -p 8 \
                    -x /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP \
                    --very-sensitive-local \
                    -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/in_log_WT_Brn1_rep3.fastq.gz \
                        | samtools sort \
                            -@ 8 \
                            -T /fh/scratch/delete30/tsukiyama_t \
                            -O bam \
                            -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_log_WT_Brn1_rep3.bam
            } \
                 > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_log_WT_Brn1_rep3.bowtie2_samtools-sort.stdout.txt) \
                2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_log_WT_Brn1_rep3.bowtie2_samtools-sort.stderr.txt)

            #  Index the sorted bams
            if [[ -f /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_log_WT_Brn1_rep3.bam ]]; then
                samtools index \
                    -@ 8 \
                    /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_log_WT_Brn1_rep3.bam \
                         > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_log_WT_Brn1_rep3.samtools-index.stdout.txt) \
                        2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_log_WT_Brn1_rep3.samtools-index.stderr.txt)
            fi


            #  ---------------------------------------------------------
            #  Align untrimmed fastqs; sort resulting bams
            {
                bowtie2 \
                    -p 8 \
                    -x /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP \
                    --very-sensitive-local \
                    -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/in_log_WT_Brn1_repM.fastq.gz \
                        | samtools sort \
                            -@ 8 \
                            -T /fh/scratch/delete30/tsukiyama_t \
                            -O bam \
                            -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_log_WT_Brn1_repM.bam
            } \
                 > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_log_WT_Brn1_repM.bowtie2_samtools-sort.stdout.txt) \
                2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_log_WT_Brn1_repM.bowtie2_samtools-sort.stderr.txt)

            #  Index the sorted bams
            if [[ -f /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_log_WT_Brn1_repM.bam ]]; then
                samtools index \
                    -@ 8 \
                    /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_log_WT_Brn1_repM.bam \
                         > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_log_WT_Brn1_repM.samtools-index.stdout.txt) \
                        2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_log_WT_Brn1_repM.samtools-index.stderr.txt)
            fi


            #  ---------------------------------------------------------
            #  Align untrimmed fastqs; sort resulting bams
            {
                bowtie2 \
                    -p 8 \
                    -x /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP \
                    --very-sensitive-local \
                    -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/in_Q_WT_Brn1_rep1.fastq.gz \
                        | samtools sort \
                            -@ 8 \
                            -T /fh/scratch/delete30/tsukiyama_t \
                            -O bam \
                            -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_Q_WT_Brn1_rep1.bam
            } \
                 > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_Q_WT_Brn1_rep1.bowtie2_samtools-sort.stdout.txt) \
                2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_Q_WT_Brn1_rep1.bowtie2_samtools-sort.stderr.txt)

            #  Index the sorted bams
            if [[ -f /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_Q_WT_Brn1_rep1.bam ]]; then
                samtools index \
                    -@ 8 \
                    /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_Q_WT_Brn1_rep1.bam \
                         > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_Q_WT_Brn1_rep1.samtools-index.stdout.txt) \
                        2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_Q_WT_Brn1_rep1.samtools-index.stderr.txt)
            fi


            #  ---------------------------------------------------------
            #  Align untrimmed fastqs; sort resulting bams
            {
                bowtie2 \
                    -p 8 \
                    -x /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP \
                    --very-sensitive-local \
                    -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/in_Q_WT_Brn1_rep2.fastq.gz \
                        | samtools sort \
                            -@ 8 \
                            -T /fh/scratch/delete30/tsukiyama_t \
                            -O bam \
                            -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_Q_WT_Brn1_rep2.bam
            } \
                 > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_Q_WT_Brn1_rep2.bowtie2_samtools-sort.stdout.txt) \
                2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_Q_WT_Brn1_rep2.bowtie2_samtools-sort.stderr.txt)

            #  Index the sorted bams
            if [[ -f /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_Q_WT_Brn1_rep2.bam ]]; then
                samtools index \
                    -@ 8 \
                    /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_Q_WT_Brn1_rep2.bam \
                         > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_Q_WT_Brn1_rep2.samtools-index.stdout.txt) \
                        2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_Q_WT_Brn1_rep2.samtools-index.stderr.txt)
            fi


            #  ---------------------------------------------------------
            #  Align untrimmed fastqs; sort resulting bams
            {
                bowtie2 \
                    -p 8 \
                    -x /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP \
                    --very-sensitive-local \
                    -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/in_Q_WT_Brn1_rep3.fastq.gz \
                        | samtools sort \
                            -@ 8 \
                            -T /fh/scratch/delete30/tsukiyama_t \
                            -O bam \
                            -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_Q_WT_Brn1_rep3.bam
            } \
                 > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_Q_WT_Brn1_rep3.bowtie2_samtools-sort.stdout.txt) \
                2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_Q_WT_Brn1_rep3.bowtie2_samtools-sort.stderr.txt)

            #  Index the sorted bams
            if [[ -f /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_Q_WT_Brn1_rep3.bam ]]; then
                samtools index \
                    -@ 8 \
                    /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_Q_WT_Brn1_rep3.bam \
                         > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_Q_WT_Brn1_rep3.samtools-index.stdout.txt) \
                        2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_Q_WT_Brn1_rep3.samtools-index.stderr.txt)
            fi


            #  ---------------------------------------------------------
            #  Align untrimmed fastqs; sort resulting bams
            {
                bowtie2 \
                    -p 8 \
                    -x /home/kalavatt/tsukiyamalab/Kris/genomes/combined_SC_SP/bowtie2/combined_SC_SP \
                    --very-sensitive-local \
                    -U /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/data/PRJNA471802/sym/in_Q_WT_Brn1_repM.fastq.gz \
                        | samtools sort \
                            -@ 8 \
                            -T /fh/scratch/delete30/tsukiyama_t \
                            -O bam \
                            -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_Q_WT_Brn1_repM.bam
            } \
                 > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_Q_WT_Brn1_repM.bowtie2_samtools-sort.stdout.txt) \
                2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_Q_WT_Brn1_repM.bowtie2_samtools-sort.stderr.txt)

            #  Index the sorted bams
            if [[ -f /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_Q_WT_Brn1_repM.bam ]]; then
                samtools index \
                    -@ 8 \
                    /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_Q_WT_Brn1_repM.bam \
                         > >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_Q_WT_Brn1_repM.samtools-index.stdout.txt) \
                        2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/err_out/in_Q_WT_Brn1_repM.samtools-index.stderr.txt)
            fi


❯ run=TRUE


❯ [[ "${run}" == TRUE ]] &&
>     {
>         for i in "${fastqs[@]}"; do
>             # i="${fastqs[0]}"  # echo "${i}"
>             in_fastq="${i}"  # ., "${in_fastq}"
>             out_bam="${d_bams}/$(basename "${in_fastq}" .fastq.gz).bam"  # echo "${out_bam}"
> 
>             # [[ -f "${err_out}/$(basename "${out_bam}" .bam).bowtie2_samtools-sort.stdout.txt" ]] && \
>             #     rm "${err_out}/$(basename "${out_bam}" .bam).bowtie2_samtools-sort.stdout.txt"
>             #
>             # [[ -f "${err_out}/$(basename "${out_bam}" .bam).bowtie2_samtools-sort.stderr.txt" ]] && \
>             #     rm "${err_out}/$(basename "${out_bam}" .bam).bowtie2_samtools-sort.stderr.txt"
> 
>             echo "#### Aligning, sorting, indexing $(basename ${in_fastq}) ####"
>             #  Align untrimmed fastqs; sort resulting bams
>             {
>                 bowtie2 \
>                     -p "${threads}" \
>                     --very-sensitive-local \
>                     -x "${f_indices}" \
>                     -U "${in_fastq}" \
>                         | samtools sort \
>                             -@ "${threads}" \
>                             -T "${scratch}" \
>                             -O "bam" \
>                             -o "${out_bam}"
>             } \
>                  > "${err_out}/$(basename "${out_bam}" .bam).bowtie2_samtools-sort.stdout.txt" \
>                 2> "${err_out}/$(basename "${out_bam}" .bam).bowtie2_samtools-sort.stderr.txt"
> 
>             #  Index the sorted bams
>             if [[ -f "${out_bam}" ]]; then
>                 samtools index \
>                     -@ "${threads}" \
>                     "${out_bam}" \
>                          > >(tee -a "${err_out}/$(basename "${out_bam}" .bam).samtools-index.stdout.txt") \
>                         2> >(tee -a "${err_out}/$(basename "${out_bam}" .bam).samtools-index.stderr.txt")
>             fi
> 
>             if [[ -f "${out_bam}" && -f "${out_bam}.bai" ]]; then
>                 echo "#DONE"
>             else
>                 echo "#PROBLEM"
>             fi
>             echo ""
>         done
>     }
#### Aligning, sorting, indexing Ch_log_WT_Brn1_rep1.fastq.gz ####
#DONE

#### Aligning, sorting, indexing Ch_log_WT_Brn1_rep2.fastq.gz ####
#DONE

#### Aligning, sorting, indexing Ch_log_WT_Brn1_rep3.fastq.gz ####
#DONE

#### Aligning, sorting, indexing Ch_log_WT_Brn1_repM.fastq.gz ####
#DONE

#### Aligning, sorting, indexing Ch_Q_WT_Brn1_rep1.fastq.gz ####
#DONE

#### Aligning, sorting, indexing Ch_Q_WT_Brn1_rep2.fastq.gz ####
#DONE

#### Aligning, sorting, indexing Ch_Q_WT_Brn1_rep3.fastq.gz ####
#DONE

#### Aligning, sorting, indexing Ch_Q_WT_Brn1_repM.fastq.gz ####
#DONE

#### Aligning, sorting, indexing in_log_WT_Brn1_rep1.fastq.gz ####
#DONE

#### Aligning, sorting, indexing in_log_WT_Brn1_rep2.fastq.gz ####
#DONE

#### Aligning, sorting, indexing in_log_WT_Brn1_rep3.fastq.gz ####
#DONE

#### Aligning, sorting, indexing in_log_WT_Brn1_repM.fastq.gz ####
#DONE

#### Aligning, sorting, indexing in_Q_WT_Brn1_rep1.fastq.gz ####
#DONE

#### Aligning, sorting, indexing in_Q_WT_Brn1_rep2.fastq.gz ####
#DONE

#### Aligning, sorting, indexing in_Q_WT_Brn1_rep3.fastq.gz ####
#DONE

#### Aligning, sorting, indexing in_Q_WT_Brn1_repM.fastq.gz ####
#DONE
```
</details>
<br />

<a id="examine-flags-in-bam-outfiles"></a>
### Examine flags in bam outfiles
<a id="initialize-necessary-functions"></a>
#### Initialize necessary functions
<a id="printed-5"></a>
##### Printed
<details>
<summary><i>Printed</i></summary>

```txt
❯ check_dependency() {
>     what="""
>     check_dependency()
>     ------------------
>     Check if a program is available in the current environment
> 
>     :param 1: program to check <chr>
> 
>     #DONE Check that params are not empty
>     #TODO Check that params are appropriate formats/strings
>     #TODO Print to stderr
>     """
> 
>     warning="""
>     WARNING: param 1 is empty; stopping the function
>     """
> 
>     [[ -z "${1}" ]] &&
>         {
>             echo "${warning}"
>             echo "${what}"
> 
>             return 1
>         }
> 
>     command -v "${1}" &>/dev/null ||
>         {
>             echo "Warning: ${1} not found. Install or load ${1}. Stopping the function."
>             return 1
>         }
> }


❯ calculate_run_time() {
>     what="""
>     calculate_run_time()
>     --------------------
>     Calculate run time for chunk of code
> 
>     :param 1: start time in \$(date +%s) format
>     :param 2: end time in \$(date +%s) format
>     :param 3: message to be displayed when printing the run time <chr>
> 
>     #DONE Check that params are not empty
>     #TODO Check that params are appropriate formats/strings
>     #TODO Print to stderr
>     """
> 
>     warning="""
>     WARNING: param(s) 1, 2, and/or 3 is/are empty; stopping the function
>     """
> 
>     [[ -z "${1}" || -z "${2}" || -z "${3}" ]] &&
>         {
>             echo "${warning}"
>             echo "${what}"
> 
>             return 1
>         }
> 
>     run_time="$(echo "${2}" - "${1}" | bc -l)"
> 
>     echo ""
>     echo "${3}"
>     printf 'Run time: %dh:%dm:%ds\n' \
>         $(( run_time/3600 )) \
>         $(( run_time%3600/60 )) \
>         $(( run_time%60 ))
>     echo ""
> }


❯ display_spinning_icon() {
>     what="""
>     display_spinning_icon()
>     -----------------------
>     Display \"spinning icon\" while a background process runs
> 
>     :param 1: PID of the last program the shell ran in the background <pos int>
>     :param 2: message to be displayed next to the spinning icon <chr>
> 
>     #DONE Check that params are not empty
>     #TODO Check that params are appropriate formats/strings
>     #TODO Print to stderr
>     """
> 
>     warning="""
>     WARNING: param(s) 1 and/or 2 is/are empty; stopping the function
>     """
> 
>     [[ -z "${1}" || -z "${2}" ]] &&
>         {
>             echo "${warning}"
>             echo "${what}"
> 
>             return 1
>         }
> 
>     spin="/|\\–"
>     i=0
>     while kill -0 "${1}" 2> /dev/null; do
>         i=$(( (i + 1) % 4 ))
>         printf "\r${spin:$i:1} %s" "${2}"
>         sleep .15
>     done
> }


❯ list_tally_flags() {
>     what="""
>     list_tally_flags()
>     ------------------
>     List and tally flags in a bam infile; function acts on a bam infile to
>     perform piped commands (samtools view, cut, sort, uniq -c, sort -nr) that
>     list and tally flags; function writes the results to a txt outfile, the
>     name of which is derived from the infile
> 
>     :param 1: name of bam infile, including path (chr)
> 
>     #DONE Check that param string is not empty
>     #DONE Check that param file exists
>     #DONE Check that dependency is present
>     #TODO Check that params are appropriate formats/strings
>     #TODO Print to stderr
>     """
> 
>     warning_param="""
>     WARNING: param 1 is empty; stopping the function
>     """
> 
>     warning_file="""
>     WARNING: param 1 file not found; stopping the function
>     """
> 
>     warning_depend="""
>     WARNING: One or more dependencies not found; stopping the function
>     """
> 
>     [[ -z "${1}" ]] &&
>         {
>             echo "${warning_param}"
>             echo "${what}"
>             return 1
>         }
> 
>     [[ -f "${1}" ]] ||
>         {
>             echo "${warning_file}"
>             echo "${what}"
>             return 1
>         }
> 
>     check_dependency samtools \
>         && check_dependency sort \
>         && check_dependency uniq \
>         && check_dependency display_spinning_icon \
>         && check_dependency calculate_run_time
>     [[ $? -gt 0 ]] &&
>         {
>             echo "${warning_depend}"
>             return 1
>         }
> 
>     start="$(date +%s)"
> 
>     samtools view "${1}" \
>         | cut -f 2 \
>         | sort \
>         | uniq -c \
>         | sort -nr \
>             > "${1/.bam/.flags.txt}" &
>     display_spinning_icon $! \
>     "Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on $(basename "${1}")... "
> 
>     end="$(date +%s)"
>     echo ""
>     calculate_run_time "${start}" "${end}"  \
>     "List and tally flags in $(basename "${1}")."
> }


❯ list_tally_MAPQs() {
>     what="""
>     list_tally_MAPQs()
>     ------------------
>     List and tally MAPQ scores in a bam infile; function acts on a bam infile
>     to perform piped commands (samtools view, cut, sort, uniq -c, sort -nr)
>     that list and tally MAPQ scores; function writes the results to a txt
>     outfile, the name of which is derived from the infile
> 
>     :param 1: name of bam infile, including path (chr)
> 
>     #DONE Check that param string is not empty
>     #DONE Check that param file exists
>     #DONE Check that dependency is present
>     #TODO Check that params are appropriate formats/strings
>     #TODO Print to stderr
>     """
> 
>     warning_param="""
>     WARNING: param 1 is empty; stopping the function
>     """
> 
>     warning_file="""
>     WARNING: param 1 file not found; stopping the function
>     """
> 
>     warning_depend="""
>     WARNING: One or more dependencies not found; stopping the function
>     """
> 
>     [[ -z "${1}" ]] &&
>         {
>             echo "${warning_param}"
>             echo "${what}"
>             return 1
>         }
> 
>     [[ -f "${1}" ]] ||
>         {
>             echo "${warning_file}"
>             echo "${what}"
>             return 1
>         }
> 
>     check_dependency samtools \
>         && check_dependency sort \
>         && check_dependency uniq \
>         && check_dependency display_spinning_icon \
>         && check_dependency calculate_run_time
>     [[ $? -gt 0 ]] &&
>         {
>             echo "${warning_depend}"
>             return 1
>         }
> 
>     start="$(date +%s)"
> 
>     samtools view "${1}" \
>         | cut -f 5 \
>         | sort \
>         | uniq -c \
>         | sort -nr \
>             > "${1/.bam/.MAPQs.txt}" &
>     display_spinning_icon $! \
>     "Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on $(basename "${1}")... "
> 
>     end="$(date +%s)"
>     echo ""
>     calculate_run_time "${start}" "${end}"  \
>     "List and tally MAPQ scores in $(basename "${1}")."
> }
```
</details>
<br />

<a id="initialize-an-array-of-bams"></a>
#### Initialize an array of bams
<a id="printed-6"></a>
##### Printed
<details>
<summary><i>Printed</i></summary>

```txt
❯ p_data="${HOME}/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams"


❯ unset bams


❯ typeset -a bams=(
>     "${p_data}/Ch_log_WT_Brn1_rep1.atria.bam"
>     "${p_data}/Ch_log_WT_Brn1_rep1.bam"
>     "${p_data}/Ch_log_WT_Brn1_rep2.atria.bam"
>     "${p_data}/Ch_log_WT_Brn1_rep2.bam"
>     "${p_data}/Ch_log_WT_Brn1_rep3.atria.bam"
>     "${p_data}/Ch_log_WT_Brn1_rep3.bam"
>     "${p_data}/Ch_log_WT_Brn1_repM.atria.bam"
>     "${p_data}/Ch_log_WT_Brn1_repM.bam"
>     "${p_data}/Ch_Q_WT_Brn1_rep1.atria.bam"
>     "${p_data}/Ch_Q_WT_Brn1_rep1.bam"
>     "${p_data}/Ch_Q_WT_Brn1_rep2.atria.bam"
>     "${p_data}/Ch_Q_WT_Brn1_rep2.bam"
>     "${p_data}/Ch_Q_WT_Brn1_rep3.atria.bam"
>     "${p_data}/Ch_Q_WT_Brn1_rep3.bam"
>     "${p_data}/Ch_Q_WT_Brn1_repM.atria.bam"
>     "${p_data}/Ch_Q_WT_Brn1_repM.bam"
>     "${p_data}/in_log_WT_Brn1_rep1.atria.bam"
>     "${p_data}/in_log_WT_Brn1_rep1.bam"
>     "${p_data}/in_log_WT_Brn1_rep2.atria.bam"
>     "${p_data}/in_log_WT_Brn1_rep2.bam"
>     "${p_data}/in_log_WT_Brn1_rep3.atria.bam"
>     "${p_data}/in_log_WT_Brn1_rep3.bam"
>     "${p_data}/in_log_WT_Brn1_repM.atria.bam"
>     "${p_data}/in_log_WT_Brn1_repM.bam"
>     "${p_data}/in_Q_WT_Brn1_rep1.atria.bam"
>     "${p_data}/in_Q_WT_Brn1_rep1.bam"
>     "${p_data}/in_Q_WT_Brn1_rep2.atria.bam"
>     "${p_data}/in_Q_WT_Brn1_rep2.bam"
>     "${p_data}/in_Q_WT_Brn1_rep3.atria.bam"
>     "${p_data}/in_Q_WT_Brn1_rep3.bam"
>     "${p_data}/in_Q_WT_Brn1_repM.atria.bam"
>     "${p_data}/in_Q_WT_Brn1_repM.bam"
> )


❯ run_check=TRUE


❯ [[ "${run_check}" == TRUE ]] &&
>     {
>         for i in "${bams[@]}"; do ls -lhaFG "${i}"; done
>     }
-rw-rw---- 1 kalavatt 132M Jul 16 15:59 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_log_WT_Brn1_rep1.atria.bam
-rw-rw---- 1 kalavatt 142M Jul 16 14:29 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_log_WT_Brn1_rep1.bam
-rw-rw---- 1 kalavatt 585M Jul 16 20:31 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_log_WT_Brn1_rep2.atria.bam
-rw-rw---- 1 kalavatt 635M Jul 16 14:33 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_log_WT_Brn1_rep2.bam
-rw-rw---- 1 kalavatt 180M Jul 16 20:31 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_log_WT_Brn1_rep3.atria.bam
-rw-rw---- 1 kalavatt 196M Jul 16 14:34 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_log_WT_Brn1_rep3.bam
-rw-rw---- 1 kalavatt 863M Jul 16 21:17 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_log_WT_Brn1_repM.atria.bam
-rw-rw---- 1 kalavatt 940M Jul 16 14:39 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_log_WT_Brn1_repM.bam
-rw-rw---- 1 kalavatt 172M Jul 16 21:17 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_Q_WT_Brn1_rep1.atria.bam
-rw-rw---- 1 kalavatt 181M Jul 16 14:39 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_Q_WT_Brn1_rep1.bam
-rw-rw---- 1 kalavatt 150M Jul 16 21:18 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_Q_WT_Brn1_rep2.atria.bam
-rw-rw---- 1 kalavatt 159M Jul 16 14:40 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_Q_WT_Brn1_rep2.bam
-rw-rw---- 1 kalavatt 209M Jul 16 21:18 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_Q_WT_Brn1_rep3.atria.bam
-rw-rw---- 1 kalavatt 230M Jul 16 14:41 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_Q_WT_Brn1_rep3.bam
-rw-rw---- 1 kalavatt 506M Jul 16 21:20 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_Q_WT_Brn1_repM.atria.bam
-rw-rw---- 1 kalavatt 547M Jul 16 14:43 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_Q_WT_Brn1_repM.bam
-rw-rw---- 1 kalavatt 348M Jul 16 21:21 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_log_WT_Brn1_rep1.atria.bam
-rw-rw---- 1 kalavatt 385M Jul 16 14:45 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_log_WT_Brn1_rep1.bam
-rw-rw---- 1 kalavatt 279M Jul 16 21:22 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_log_WT_Brn1_rep2.atria.bam
-rw-rw---- 1 kalavatt 309M Jul 16 14:46 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_log_WT_Brn1_rep2.bam
-rw-rw---- 1 kalavatt 199M Jul 16 21:23 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_log_WT_Brn1_rep3.atria.bam
-rw-rw---- 1 kalavatt 216M Jul 16 14:47 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_log_WT_Brn1_rep3.bam
-rw-rw---- 1 kalavatt 781M Jul 16 22:20 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_log_WT_Brn1_repM.atria.bam
-rw-rw---- 1 kalavatt 865M Jul 16 14:52 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_log_WT_Brn1_repM.bam
-rw-rw---- 1 kalavatt 257M Jul 16 22:21 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_Q_WT_Brn1_rep1.atria.bam
-rw-rw---- 1 kalavatt 284M Jul 16 14:53 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_Q_WT_Brn1_rep1.bam
-rw-rw---- 1 kalavatt 294M Jul 16 22:22 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_Q_WT_Brn1_rep2.atria.bam
-rw-rw---- 1 kalavatt 325M Jul 16 14:54 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_Q_WT_Brn1_rep2.bam
-rw-rw---- 1 kalavatt 249M Jul 16 22:22 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_Q_WT_Brn1_rep3.atria.bam
-rw-rw---- 1 kalavatt 266M Jul 16 14:55 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_Q_WT_Brn1_rep3.bam
-rw-rw---- 1 kalavatt 756M Jul 16 22:50 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_Q_WT_Brn1_repM.atria.bam
-rw-rw---- 1 kalavatt 832M Jul 16 15:00 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_Q_WT_Brn1_repM.bam


❯ for i in "${bams[@]}"; do
>     # i="${bams[0]}"  # echo "${i}"
>     echo "#### $(basename ${i}) ####"
>     list_tally_MAPQs "${i}"
>     echo ""
> done
#### Ch_log_WT_Brn1_rep1.atria.bam ####
[1] 64378
| Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Ch_log_WT_Brn1_rep1.atria.bam... [1]+  Done

samtools view "${1}" | cut -f 5 | sort | uniq -c | sort -nr > "${1/.bam/.MAPQs.txt}"

List and tally MAPQ scores in Ch_log_WT_Brn1_rep1.atria.bam.
Run time: 0h:0m:7s


#### Ch_log_WT_Brn1_rep1.bam ####
[1] 64679
/ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Ch_log_WT_Brn1_rep1.bam... [1]+  Done

samtools view "${1}" | cut -f 5 | sort | uniq -c | sort -nr > "${1/.bam/.MAPQs.txt}"

List and tally MAPQ scores in Ch_log_WT_Brn1_rep1.bam.
Run time: 0h:0m:7s


#### Ch_log_WT_Brn1_rep2.atria.bam ####
[1] 64866
| Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Ch_log_WT_Brn1_rep2.atria.bam... [1]+  Done

samtools view "${1}" | cut -f 5 | sort | uniq -c | sort -nr > "${1/.bam/.MAPQs.txt}"

List and tally MAPQ scores in Ch_log_WT_Brn1_rep2.atria.bam.
Run time: 0h:0m:32s


#### Ch_log_WT_Brn1_rep2.bam ####
[1] 66137
– Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Ch_log_WT_Brn1_rep2.bam... [1]+  Done

samtools view "${1}" | cut -f 5 | sort | uniq -c | sort -nr > "${1/.bam/.MAPQs.txt}"

List and tally MAPQ scores in Ch_log_WT_Brn1_rep2.bam.
Run time: 0h:0m:33s


#### Ch_log_WT_Brn1_rep3.atria.bam ####
[1] 66684
/ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Ch_log_WT_Brn1_rep3.atria.bam... [1]+  Done

samtools view "${1}" | cut -f 5 | sort | uniq -c | sort -nr > "${1/.bam/.MAPQs.txt}"

List and tally MAPQ scores in Ch_log_WT_Brn1_rep3.atria.bam.
Run time: 0h:0m:10s


#### Ch_log_WT_Brn1_rep3.bam ####
[1] 67160
\ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Ch_log_WT_Brn1_rep3.bam... [1]+  Done

samtools view "${1}" | cut -f 5 | sort | uniq -c | sort -nr > "${1/.bam/.MAPQs.txt}"

List and tally MAPQ scores in Ch_log_WT_Brn1_rep3.bam.
Run time: 0h:0m:11s


#### Ch_log_WT_Brn1_repM.atria.bam ####
[1] 67858
\ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Ch_log_WT_Brn1_repM.atria.bam... [1]+  Done

samtools view "${1}" | cut -f 5 | sort | uniq -c | sort -nr > "${1/.bam/.MAPQs.txt}"

List and tally MAPQ scores in Ch_log_WT_Brn1_repM.atria.bam.
Run time: 0h:0m:49s


#### Ch_log_WT_Brn1_repM.bam ####
[1] 69105
– Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Ch_log_WT_Brn1_repM.bam... [1]+  Done

samtools view "${1}" | cut -f 5 | sort | uniq -c | sort -nr > "${1/.bam/.MAPQs.txt}"

List and tally MAPQ scores in Ch_log_WT_Brn1_repM.bam.
Run time: 0h:0m:54s


#### Ch_Q_WT_Brn1_rep1.atria.bam ####
[1] 71317
\ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Ch_Q_WT_Brn1_rep1.atria.bam... [1]+  Done

samtools view "${1}" | cut -f 5 | sort | uniq -c | sort -nr > "${1/.bam/.MAPQs.txt}"

List and tally MAPQ scores in Ch_Q_WT_Brn1_rep1.atria.bam.
Run time: 0h:0m:9s


#### Ch_Q_WT_Brn1_rep1.bam ####
[1] 71669
\ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Ch_Q_WT_Brn1_rep1.bam... [1]+  Done

samtools view "${1}" | cut -f 5 | sort | uniq -c | sort -nr > "${1/.bam/.MAPQs.txt}"

List and tally MAPQ scores in Ch_Q_WT_Brn1_rep1.bam.
Run time: 0h:0m:9s


#### Ch_Q_WT_Brn1_rep2.atria.bam ####
[1] 71867
– Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Ch_Q_WT_Brn1_rep2.atria.bam... [1]+  Done

samtools view "${1}" | cut -f 5 | sort | uniq -c | sort -nr > "${1/.bam/.MAPQs.txt}"

List and tally MAPQ scores in Ch_Q_WT_Brn1_rep2.atria.bam.
Run time: 0h:0m:9s


#### Ch_Q_WT_Brn1_rep2.bam ####
[1] 72049
| Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Ch_Q_WT_Brn1_rep2.bam... [1]+  Done

samtools view "${1}" | cut -f 5 | sort | uniq -c | sort -nr > "${1/.bam/.MAPQs.txt}"

List and tally MAPQ scores in Ch_Q_WT_Brn1_rep2.bam.
Run time: 0h:0m:8s


#### Ch_Q_WT_Brn1_rep3.atria.bam ####
[1] 72213
| Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Ch_Q_WT_Brn1_rep3.atria.bam... [1]+  Done

samtools view "${1}" | cut -f 5 | sort | uniq -c | sort -nr > "${1/.bam/.MAPQs.txt}"

List and tally MAPQ scores in Ch_Q_WT_Brn1_rep3.atria.bam.
Run time: 0h:0m:10s


#### Ch_Q_WT_Brn1_rep3.bam ####
[1] 72771
| Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Ch_Q_WT_Brn1_rep3.bam... [1]+  Done

samtools view "${1}" | cut -f 5 | sort | uniq -c | sort -nr > "${1/.bam/.MAPQs.txt}"

List and tally MAPQ scores in Ch_Q_WT_Brn1_rep3.bam.
Run time: 0h:0m:11s


#### Ch_Q_WT_Brn1_repM.atria.bam ####
[1] 73124
– Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Ch_Q_WT_Brn1_repM.atria.bam... [1]+  Done

samtools view "${1}" | cut -f 5 | sort | uniq -c | sort -nr > "${1/.bam/.MAPQs.txt}"

List and tally MAPQ scores in Ch_Q_WT_Brn1_repM.atria.bam.
Run time: 0h:0m:28s


#### Ch_Q_WT_Brn1_repM.bam ####
[1] 523
– Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Ch_Q_WT_Brn1_repM.bam... [1]+  Done

samtools view "${1}" | cut -f 5 | sort | uniq -c | sort -nr > "${1/.bam/.MAPQs.txt}"

List and tally MAPQ scores in Ch_Q_WT_Brn1_repM.bam.
Run time: 0h:0m:27s


#### in_log_WT_Brn1_rep1.atria.bam ####
[1] 1786
/ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on in_log_WT_Brn1_rep1.atria.bam... [1]+  Done

samtools view "${1}" | cut -f 5 | sort | uniq -c | sort -nr > "${1/.bam/.MAPQs.txt}"

List and tally MAPQ scores in in_log_WT_Brn1_rep1.atria.bam.
Run time: 0h:0m:16s


#### in_log_WT_Brn1_rep1.bam ####
[1] 2369
/ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on in_log_WT_Brn1_rep1.bam... [1]+  Done

samtools view "${1}" | cut -f 5 | sort | uniq -c | sort -nr > "${1/.bam/.MAPQs.txt}"

List and tally MAPQ scores in in_log_WT_Brn1_rep1.bam.
Run time: 0h:0m:18s


#### in_log_WT_Brn1_rep2.atria.bam ####
[1] 3124
– Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on in_log_WT_Brn1_rep2.atria.bam... [1]+  Done

samtools view "${1}" | cut -f 5 | sort | uniq -c | sort -nr > "${1/.bam/.MAPQs.txt}"

List and tally MAPQ scores in in_log_WT_Brn1_rep2.atria.bam.
Run time: 0h:0m:13s


#### in_log_WT_Brn1_rep2.bam ####
[1] 4012
/ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on in_log_WT_Brn1_rep2.bam... [1]+  Done

samtools view "${1}" | cut -f 5 | sort | uniq -c | sort -nr > "${1/.bam/.MAPQs.txt}"

List and tally MAPQ scores in in_log_WT_Brn1_rep2.bam.
Run time: 0h:0m:15s


#### in_log_WT_Brn1_rep3.atria.bam ####
[1] 4285
/ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on in_log_WT_Brn1_rep3.atria.bam... [1]+  Done

samtools view "${1}" | cut -f 5 | sort | uniq -c | sort -nr > "${1/.bam/.MAPQs.txt}"

List and tally MAPQ scores in in_log_WT_Brn1_rep3.atria.bam.
Run time: 0h:0m:9s


#### in_log_WT_Brn1_rep3.bam ####
[1] 4439
\ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on in_log_WT_Brn1_rep3.bam... [1]+  Done

samtools view "${1}" | cut -f 5 | sort | uniq -c | sort -nr > "${1/.bam/.MAPQs.txt}"

List and tally MAPQ scores in in_log_WT_Brn1_rep3.bam.
Run time: 0h:0m:11s


#### in_log_WT_Brn1_repM.atria.bam ####
[1] 4943
| Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on in_log_WT_Brn1_repM.atria.bam... [1]+  Done

samtools view "${1}" | cut -f 5 | sort | uniq -c | sort -nr > "${1/.bam/.MAPQs.txt}"

List and tally MAPQ scores in in_log_WT_Brn1_repM.atria.bam.
Run time: 0h:0m:39s


#### in_log_WT_Brn1_repM.bam ####
[1] 6697
/ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on in_log_WT_Brn1_repM.bam... [1]+  Done

samtools view "${1}" | cut -f 5 | sort | uniq -c | sort -nr > "${1/.bam/.MAPQs.txt}"

List and tally MAPQ scores in in_log_WT_Brn1_repM.bam.
Run time: 0h:0m:42s


#### in_Q_WT_Brn1_rep1.atria.bam ####
[1] 8825
\ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on in_Q_WT_Brn1_rep1.atria.bam... [1]+  Done

samtools view "${1}" | cut -f 5 | sort | uniq -c | sort -nr > "${1/.bam/.MAPQs.txt}"

List and tally MAPQ scores in in_Q_WT_Brn1_rep1.atria.bam.
Run time: 0h:0m:12s


#### in_Q_WT_Brn1_rep1.bam ####
[1] 9127
– Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on in_Q_WT_Brn1_rep1.bam... [1]+  Done

samtools view "${1}" | cut -f 5 | sort | uniq -c | sort -nr > "${1/.bam/.MAPQs.txt}"

List and tally MAPQ scores in in_Q_WT_Brn1_rep1.bam.
Run time: 0h:0m:15s


#### in_Q_WT_Brn1_rep2.atria.bam ####
[1] 10057
\ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on in_Q_WT_Brn1_rep2.atria.bam... [1]+  Done

samtools view "${1}" | cut -f 5 | sort | uniq -c | sort -nr > "${1/.bam/.MAPQs.txt}"

List and tally MAPQ scores in in_Q_WT_Brn1_rep2.atria.bam.
Run time: 0h:0m:15s


#### in_Q_WT_Brn1_rep2.bam ####
[1] 11168
| Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on in_Q_WT_Brn1_rep2.bam... [1]+  Done

samtools view "${1}" | cut -f 5 | sort | uniq -c | sort -nr > "${1/.bam/.MAPQs.txt}"

List and tally MAPQ scores in in_Q_WT_Brn1_rep2.bam.
Run time: 0h:0m:15s


#### in_Q_WT_Brn1_rep3.atria.bam ####
[1] 11680
| Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on in_Q_WT_Brn1_rep3.atria.bam... [1]+  Done

samtools view "${1}" | cut -f 5 | sort | uniq -c | sort -nr > "${1/.bam/.MAPQs.txt}"

List and tally MAPQ scores in in_Q_WT_Brn1_rep3.atria.bam.
Run time: 0h:0m:12s


#### in_Q_WT_Brn1_rep3.bam ####
[1] 12098
– Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on in_Q_WT_Brn1_rep3.bam... [1]+  Done

samtools view "${1}" | cut -f 5 | sort | uniq -c | sort -nr > "${1/.bam/.MAPQs.txt}"

List and tally MAPQ scores in in_Q_WT_Brn1_rep3.bam.
Run time: 0h:0m:14s


#### in_Q_WT_Brn1_repM.atria.bam ####
[1] 13193
– Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on in_Q_WT_Brn1_repM.atria.bam... [1]+  Done

samtools view "${1}" | cut -f 5 | sort | uniq -c | sort -nr > "${1/.bam/.MAPQs.txt}"

List and tally MAPQ scores in in_Q_WT_Brn1_repM.atria.bam.
Run time: 0h:0m:37s


#### in_Q_WT_Brn1_repM.bam ####
[1] 14570
– Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on in_Q_WT_Brn1_repM.bam... [1]+  Done

samtools view "${1}" | cut -f 5 | sort | uniq -c | sort -nr > "${1/.bam/.MAPQs.txt}"

List and tally MAPQ scores in in_Q_WT_Brn1_repM.bam.
Run time: 0h:0m:40s
```
</details>
<br />

<a id="check-on-flag-information-in-bams"></a>
#### Check on flag information in bams
<a id="printed-7"></a>
##### Printed
<details>
<summary><i>Printed</i></summary>

```txt
❯ for i in "${bams[@]}"; do
>     # i="${bams[0]}"  # echo "${i}"
>     echo "#### $(basename ${i}) ####"
>     list_tally_flags "${i}"
>     echo ""
> done
#### Ch_log_WT_Brn1_rep1.atria.bam ####
[1] 59460

| Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Ch_log_WT_Brn1_rep1.atria.bam... [1]+  Done
samtools view "${1}" | cut -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in Ch_log_WT_Brn1_rep1.atria.bam.
Run time: 0h:0m:10s


#### Ch_log_WT_Brn1_rep1.bam ####
[1] 59669
\ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Ch_log_WT_Brn1_rep1.bam... [1]+  Done
samtools view "${1}" | cut -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in Ch_log_WT_Brn1_rep1.bam.
Run time: 0h:0m:9s


#### Ch_log_WT_Brn1_rep2.atria.bam ####
[1] 60287
– Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Ch_log_WT_Brn1_rep2.atria.bam... [1]+  Done
samtools view "${1}" | cut -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in Ch_log_WT_Brn1_rep2.atria.bam.
Run time: 0h:0m:43s


#### Ch_log_WT_Brn1_rep2.bam ####
[1] 62542
\ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Ch_log_WT_Brn1_rep2.bam... [1]+  Done
samtools view "${1}" | cut -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in Ch_log_WT_Brn1_rep2.bam.
Run time: 0h:0m:48s


#### Ch_log_WT_Brn1_rep3.atria.bam ####
[1] 64344
– Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Ch_log_WT_Brn1_rep3.atria.bam... [1]+  Done
samtools view "${1}" | cut -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in Ch_log_WT_Brn1_rep3.atria.bam.
Run time: 0h:0m:13s


#### Ch_log_WT_Brn1_rep3.bam ####
[1] 64725
| Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Ch_log_WT_Brn1_rep3.bam... [1]+  Done
samtools view "${1}" | cut -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in Ch_log_WT_Brn1_rep3.bam.
Run time: 0h:0m:14s


#### Ch_log_WT_Brn1_repM.atria.bam ####
[1] 66124
– Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Ch_log_WT_Brn1_repM.atria.bam... [1]+  Done
samtools view "${1}" | cut -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in Ch_log_WT_Brn1_repM.atria.bam.
Run time: 0h:1m:17s


#### Ch_log_WT_Brn1_repM.bam ####
[1] 72459
| Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Ch_log_WT_Brn1_repM.bam... [1]+  Done
samtools view "${1}" | cut -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in Ch_log_WT_Brn1_repM.bam.
Run time: 0h:1m:11s


#### Ch_Q_WT_Brn1_rep1.atria.bam ####
[1] 3012
– Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Ch_Q_WT_Brn1_rep1.atria.bam... [1]+  Done
samtools view "${1}" | cut -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in Ch_Q_WT_Brn1_rep1.atria.bam.
Run time: 0h:0m:10s


#### Ch_Q_WT_Brn1_rep1.bam ####
[1] 3649
\ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Ch_Q_WT_Brn1_rep1.bam... [1]+  Done
samtools view "${1}" | cut -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in Ch_Q_WT_Brn1_rep1.bam.
Run time: 0h:0m:11s


#### Ch_Q_WT_Brn1_rep2.atria.bam ####
[1] 4117
| Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Ch_Q_WT_Brn1_rep2.atria.bam... [1]+  Done
samtools view "${1}" | cut -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in Ch_Q_WT_Brn1_rep2.atria.bam.
Run time: 0h:0m:10s


#### Ch_Q_WT_Brn1_rep2.bam ####
[1] 4288
– Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Ch_Q_WT_Brn1_rep2.bam... [1]+  Done
samtools view "${1}" | cut -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in Ch_Q_WT_Brn1_rep2.bam.
Run time: 0h:0m:10s


#### Ch_Q_WT_Brn1_rep3.atria.bam ####
[1] 4476
/ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Ch_Q_WT_Brn1_rep3.atria.bam... [1]+  Done
samtools view "${1}" | cut -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in Ch_Q_WT_Brn1_rep3.atria.bam.
Run time: 0h:0m:16s


#### Ch_Q_WT_Brn1_rep3.bam ####
[1] 5149
\ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Ch_Q_WT_Brn1_rep3.bam... [1]+  Done
samtools view "${1}" | cut -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in Ch_Q_WT_Brn1_rep3.bam.
Run time: 0h:0m:16s


#### Ch_Q_WT_Brn1_repM.atria.bam ####
[1] 6287
| Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Ch_Q_WT_Brn1_repM.atria.bam... [1]+  Done
samtools view "${1}" | cut -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in Ch_Q_WT_Brn1_repM.atria.bam.
Run time: 0h:0m:32s


#### Ch_Q_WT_Brn1_repM.bam ####
[1] 7080
| Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Ch_Q_WT_Brn1_repM.bam... [1]+  Done
samtools view "${1}" | cut -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in Ch_Q_WT_Brn1_repM.bam.
Run time: 0h:0m:36s


#### in_log_WT_Brn1_rep1.atria.bam ####
[1] 8664
/ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on in_log_WT_Brn1_rep1.atria.bam... [1]+  Done
samtools view "${1}" | cut -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in in_log_WT_Brn1_rep1.atria.bam.
Run time: 0h:0m:22s


#### in_log_WT_Brn1_rep1.bam ####
[1] 9160
\ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on in_log_WT_Brn1_rep1.bam... [1]+  Done
samtools view "${1}" | cut -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in in_log_WT_Brn1_rep1.bam.
Run time: 0h:0m:23s


#### in_log_WT_Brn1_rep2.atria.bam ####
[1] 10562
| Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on in_log_WT_Brn1_rep2.atria.bam... [1]+  Done
samtools view "${1}" | cut -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in in_log_WT_Brn1_rep2.atria.bam.
Run time: 0h:0m:19s


#### in_log_WT_Brn1_rep2.bam ####
[1] 11122
/ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on in_log_WT_Brn1_rep2.bam... [1]+  Done
samtools view "${1}" | cut -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in in_log_WT_Brn1_rep2.bam.
Run time: 0h:0m:20s


#### in_log_WT_Brn1_rep3.atria.bam ####
[1] 12039
| Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on in_log_WT_Brn1_rep3.atria.bam... [1]+  Done
samtools view "${1}" | cut -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in in_log_WT_Brn1_rep3.atria.bam.
Run time: 0h:0m:12s


#### in_log_WT_Brn1_rep3.bam ####
[1] 12451
– Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on in_log_WT_Brn1_rep3.bam... [1]+  Done
samtools view "${1}" | cut -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in in_log_WT_Brn1_rep3.bam.
Run time: 0h:0m:12s


#### in_log_WT_Brn1_repM.atria.bam ####
[1] 12688
| Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on in_log_WT_Brn1_repM.atria.bam... [1]+  Done
samtools view "${1}" | cut -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in in_log_WT_Brn1_repM.atria.bam.
Run time: 0h:0m:43s


#### in_log_WT_Brn1_repM.bam ####
[1] 18115
/ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on in_log_WT_Brn1_repM.bam... [1]+  Done
samtools view "${1}" | cut -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in in_log_WT_Brn1_repM.bam.
Run time: 0h:0m:45s


#### in_Q_WT_Brn1_rep1.atria.bam ####
[1] 19754
| Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on in_Q_WT_Brn1_rep1.atria.bam... [1]+  Done
samtools view "${1}" | cut -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in in_Q_WT_Brn1_rep1.atria.bam.
Run time: 0h:0m:13s


#### in_Q_WT_Brn1_rep1.bam ####
[1] 20117
\ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on in_Q_WT_Brn1_rep1.bam... [1]+  Done
samtools view "${1}" | cut -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in in_Q_WT_Brn1_rep1.bam.
Run time: 0h:0m:13s


#### in_Q_WT_Brn1_rep2.atria.bam ####
[1] 20516
| Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on in_Q_WT_Brn1_rep2.atria.bam... [1]+  Done
samtools view "${1}" | cut -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in in_Q_WT_Brn1_rep2.atria.bam.
Run time: 0h:0m:15s


#### in_Q_WT_Brn1_rep2.bam ####
[1] 21091
/ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on in_Q_WT_Brn1_rep2.bam... [1]+  Done
samtools view "${1}" | cut -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in in_Q_WT_Brn1_rep2.bam.
Run time: 0h:0m:16s


#### in_Q_WT_Brn1_rep3.atria.bam ####
[1] 21508
/ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on in_Q_WT_Brn1_rep3.atria.bam... [1]+  Done
samtools view "${1}" | cut -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in in_Q_WT_Brn1_rep3.atria.bam.
Run time: 0h:0m:13s


#### in_Q_WT_Brn1_rep3.bam ####
[1] 21734
| Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on in_Q_WT_Brn1_rep3.bam... [1]+  Done
samtools view "${1}" | cut -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in in_Q_WT_Brn1_rep3.bam.
Run time: 0h:0m:14s


#### in_Q_WT_Brn1_repM.atria.bam ####
[1] 22211
/ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on in_Q_WT_Brn1_repM.atria.bam... [1]+  Done
samtools view "${1}" | cut -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in in_Q_WT_Brn1_repM.atria.bam.
Run time: 0h:0m:37s


#### in_Q_WT_Brn1_repM.bam ####
[1] 23287
/ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on in_Q_WT_Brn1_repM.bam... [1]+  Done
samtools view "${1}" | cut -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in in_Q_WT_Brn1_repM.bam.
Run time: 0h:0m:44s


❯ for i in "${bams[@]}"; do
>     if [[ -f "${i/.bam/.flags.txt}" ]]; then
>         echo "#### $(basename ${i}) ####"
>         cat "${i/.bam/.flags.txt}"
>         echo ""
>     fi
> done
#### Ch_log_WT_Brn1_rep1.atria.bam ####
2880211 16
2872535 0
 499911 4

#### Ch_log_WT_Brn1_rep1.bam ####
3040301 16
3017984 0
 585178 4

#### Ch_log_WT_Brn1_rep2.atria.bam ####
13759369 0
13695345 16
2734535 4

#### Ch_log_WT_Brn1_rep2.bam ####
14450767 0
14447432 16
3176492 4

#### Ch_log_WT_Brn1_rep3.atria.bam ####
3933313 0
3918289 16
 316828 4

#### Ch_log_WT_Brn1_rep3.bam ####
4177590 16
4176208 0
 429674 4

#### Ch_log_WT_Brn1_repM.atria.bam ####
20568347 0
20490732 16
3551257 4

#### Ch_log_WT_Brn1_repM.bam ####
21661810 16
21648500 0
4191316 4

#### Ch_Q_WT_Brn1_rep1.atria.bam ####
2753517 0
2738599 16
1739181 4

#### Ch_Q_WT_Brn1_rep1.bam ####
2874061 0
2862832 16
1890825 4

#### Ch_Q_WT_Brn1_rep2.atria.bam ####
2713945 0
2707745 16
1066876 4

#### Ch_Q_WT_Brn1_rep2.bam ####
2846780 0
2843810 16
1190484 4

#### Ch_Q_WT_Brn1_rep3.atria.bam ####
4128527 0
4125918 16
 612389 4

#### Ch_Q_WT_Brn1_rep3.bam ####
4442444 16
4437645 0
 780519 4

#### Ch_Q_WT_Brn1_repM.atria.bam ####
9594496 0
9573809 16
3418392 4

#### Ch_Q_WT_Brn1_repM.bam ####
10156605 0
10150983 16
3861812 4

#### in_log_WT_Brn1_rep1.atria.bam ####
6779624 0
6769780 16
 333904 4

#### in_log_WT_Brn1_rep1.bam ####
7355182 0
7350937 16
 550871 4

#### in_log_WT_Brn1_rep2.atria.bam ####
5415076 0
5397393 16
 236372 4

#### in_log_WT_Brn1_rep2.bam ####
5876887 0
5862547 16
 405161 4

#### in_log_WT_Brn1_rep3.atria.bam ####
4033815 16
4033137 0
 195541 4

#### in_log_WT_Brn1_rep3.bam ####
4309499 16
4305931 0
 306893 4

#### in_log_WT_Brn1_repM.atria.bam ####
16227991 0
16200760 16
 765891 4

#### in_log_WT_Brn1_repM.bam ####
17538461 0
17522551 16
1262896 4

#### in_Q_WT_Brn1_rep1.atria.bam ####
5055045 0
5044226 16
 220686 4

#### in_Q_WT_Brn1_rep1.bam ####
5483653 0
5478220 16
 370063 4

#### in_Q_WT_Brn1_rep2.atria.bam ####
5865241 0
5848241 16
 217378 4

#### in_Q_WT_Brn1_rep2.bam ####
6351723 0
6342962 16
 365205 4

#### in_Q_WT_Brn1_rep3.atria.bam ####
5097671 0
5086953 16
 265277 4

#### in_Q_WT_Brn1_rep3.bam ####
5400769 0
5392613 16
 367131 4

#### in_Q_WT_Brn1_repM.atria.bam ####
16017543 0
15979878 16
 703297 4

#### in_Q_WT_Brn1_repM.bam ####
17235945 0
17214003 16
1102391 4


❯ for i in "${bams[@]}"; do
>     if [[ -f "${i/.bam/.MAPQs.txt}" ]]; then
>         echo "#### $(basename ${i}) ####"
>         cat "${i/.bam/.MAPQs.txt}"
>         echo ""
>     fi
> done
#### Ch_log_WT_Brn1_rep1.atria.bam ####
2945073 42
2623328 1
 518847 0
  63857 40
  22728 30
  10515 39
  10217 24
   9207 32
   8949 23
   6278 37
   6028 35
   4683 36
   3732 31
   3068 3
   3028 38
   2447 8
   2413 2
   1598 34
   1064 7
    688 15
    678 6
    617 16
    590 26
    535 4
    535 22
    511 11
    406 5
    355 14
    152 27
    149 12
    148 25
     98 18
     75 21
     39 17
     21 33

#### Ch_log_WT_Brn1_rep1.bam ####
3017395 44
2758358 1
 595062 0
  29139 42
  25740 31
  22044 37
  19791 39
  17535 11
  16904 41
  16764 38
  16123 34
  15319 40
  14456 32
  13845 33
  10208 35
   8198 36
   6686 14
   5281 25
   4871 2
   4842 28
   4571 21
   4002 22
   3774 17
   2851 18
   2481 24
   2293 9
   2041 12
   1764 16
   1125 19

#### Ch_log_WT_Brn1_rep2.atria.bam ####
14002100 42
12579782 1
2826048 0
 300446 40
 106109 30
  51017 39
  50294 24
  42217 23
  42138 32
  29795 35
  29525 37
  21726 36
  18074 31
  14773 3
  13545 38
  12122 8
  10770 2
   7806 34
   4660 7
   3175 15
   3174 6
   2835 22
   2738 26
   2721 16
   2483 11
   2295 4
   2041 5
   1655 14
    756 27
    708 12
    666 25
    411 18
    335 21
    193 17
    116 33

#### Ch_log_WT_Brn1_rep2.bam ####
14325384 44
13227785 1
3222823 0
 139813 42
 121118 31
 106859 37
  97425 39
  81467 41
  80225 11
  79924 38
  78030 40
  75369 34
  67565 32
  66338 33
  47971 35
  39465 36
  30622 14
  24762 25
  22870 2
  21809 21
  20084 22
  19636 28
  17792 17
  13323 18
  12300 24
  11040 9
   9112 12
   8303 16
   5477 19

#### Ch_log_WT_Brn1_rep3.atria.bam ####
4860491 42
2695978 1
 338356 0
 123359 40
  36526 30
  14795 24
  14210 39
  13180 23
  13128 32
   9932 37
   9491 35
   6539 36
   6056 31
   4909 38
   4134 3
   3325 8
   2853 34
   2295 2
   1526 7
    975 15
    926 26
    900 6
    803 16
    746 22
    656 11
    561 4
    521 5
    462 14
    205 27
    192 12
    144 25
    115 18
     71 21
     43 17
     27 33

#### Ch_log_WT_Brn1_rep3.bam ####
5062849 44
2874951 1
 440782 0
  52561 42
  41543 31
  30186 39
  29532 37
  29329 41
  27474 40
  24240 38
  24072 34
  20658 11
  20595 33
  19577 32
  14879 35
  12630 36
   8487 14
   6314 25
   5872 21
   5870 28
   5586 2
   5426 22
   4850 17
   3975 18
   3418 24
   2686 9
   2172 12
   1834 16
   1124 19

#### Ch_log_WT_Brn1_repM.atria.bam ####
21807635 42
17899082 1
3683398 0
 487693 40
 165362 30
  75660 39
  75273 24
  64525 32
  64267 23
  45954 37
  45164 35
  32967 36
  27827 31
  21885 3
  21454 38
  17921 8
  15498 2
  12208 34
   7195 7
   4832 6
   4827 15
   4194 16
   4172 26
   4114 22
   3759 11
   3350 4
   3019 5
   2385 14
   1158 27
   1076 12
    949 25
    628 18
    459 21
    283 17
    163 33

#### Ch_log_WT_Brn1_repM.bam ####
22405451 44
18861127 1
4258614 0
 221559 42
 188407 31
 158549 37
 147380 39
 127586 41
 121070 38
 120872 40
 118354 11
 115577 34
 101638 32
 100743 33
  72996 35
  60264 36
  46007 14
  36486 25
  33136 2
  32193 21
  30348 28
  29511 22
  26301 17
  20137 18
  18191 24
  16181 9
  13399 12
  11860 16
   7689 19

#### Ch_Q_WT_Brn1_rep1.atria.bam ####
3685023 42
1749695 0
1616895 1
  81918 40
  27594 30
   8835 23
   8748 24
   8408 32
   6889 39
   6003 35
   5562 37
   4992 31
   3155 36
   3038 38
   2804 3
   2140 8
   2102 34
   1562 2
   1089 7
    649 6
    641 15
    597 26
    518 16
    509 22
    468 11
    306 4
    305 5
    290 14
    143 12
    121 25
    113 27
     75 18
     46 21
     45 17
     19 33

#### Ch_Q_WT_Brn1_rep1.bam ####
3814348 44
1895799 0
1686784 1
  30200 31
  28293 42
  18435 41
  15980 40
  15572 39
  14899 37
  12912 33
  12847 34
  11867 38
  11467 11
   9762 32
   8001 35
   7321 36
   5408 14
   3719 22
   3453 2
   3273 28
   3057 21
   2890 17
   2830 25
   2423 18
   2412 24
   1581 9
    971 12
    768 16
    446 19

#### Ch_Q_WT_Brn1_rep2.atria.bam ####
3348796 42
1900257 1
1075903 0
  75405 40
  24088 30
   8905 24
   8449 23
   7436 32
   6110 39
   5285 35
   5012 37
   4464 31
   3018 36
   2709 38
   2484 3
   2068 8
   1844 34
   1167 2
    917 7
    602 6
    566 15
    487 26
    460 22
    424 16
    370 11
    299 5
    278 4
    226 14
    151 12
    104 25
     98 27
     79 18
     50 21
     42 17
     13 33

#### Ch_Q_WT_Brn1_rep2.bam ####
3481145 44
1988375 1
1196104 0
  29767 42
  26385 31
  18293 41
  14054 39
  13749 37
  13284 40
  11355 34
  11077 33
  10777 38
  10325 11
   9221 32
   7458 36
   7239 35
   4808 14
   4072 28
   3424 2
   3223 25
   3097 21
   3011 22
   2682 17
   2211 18
   2111 24
   1551 9
   1056 12
    716 16
    504 19

#### Ch_Q_WT_Brn1_rep3.atria.bam ####
5256157 42
2686308 1
 635770 0
 126865 40
  39140 30
  16327 24
  14030 23
  13605 39
  13571 32
   9939 35
   9590 37
   6740 36
   6715 31
   5449 38
   4720 3
   3880 2
   3697 8
   3159 34
   1929 7
   1300 6
   1217 15
   1077 16
   1004 26
    955 11
    921 22
    752 4
    669 14
    588 5
    209 27
    162 12
    149 25
     92 18
     77 21
     58 17
     13 33

#### Ch_Q_WT_Brn1_rep3.bam ####
5540008 44
2892269 1
 792488 0
  64790 42
  44453 31
  35829 41
  28370 40
  28286 39
  27757 37
  24769 11
  22852 34
  22547 38
  20278 33
  19853 32
  15468 36
  15274 35
  10938 14
   9039 28
   6068 22
   6055 17
   5760 25
   5655 21
   5596 2
   4386 18
   4165 24
   2770 9
   2048 12
   1708 16
   1129 19

#### Ch_Q_WT_Brn1_repM.atria.bam ####
12290058 42
6203473 1
3461395 0
 284215 40
  90822 30
  34010 24
  31325 23
  29512 32
  26544 39
  21233 35
  20208 37
  16118 31
  12895 36
  11100 38
   9957 3
   7897 8
   7117 34
   6588 2
   3889 7
   2540 6
   2449 15
   2005 16
   2002 26
   1922 22
   1827 11
   1360 4
   1200 14
   1173 5
    469 12
    404 27
    372 25
    246 18
    173 21
    149 17
     50 33

#### Ch_Q_WT_Brn1_repM.bam ####
12835510 44
6567349 1
3884286 0
 122810 42
 101090 31
  72506 41
  58024 39
  57626 40
  56437 37
  47052 34
  46719 11
  45154 38
  44204 33
  38777 32
  30546 35
  30248 36
  21081 14
  16408 28
  12767 22
  12595 2
  11873 25
  11732 21
  11702 17
   8928 18
   8699 24
   6034 9
   4076 12
   3133 16
   2034 19

#### in_log_WT_Brn1_rep1.atria.bam ####
10590105 42
2493180 1
 351181 0
 201739 40
  76537 30
  24065 24
  22332 32
  20440 23
  16185 35
  15150 39
  13824 37
  12460 31
   7766 38
   7398 36
   5698 3
   5310 34
   4982 8
   2728 7
   2044 2
   1841 6
   1303 15
   1148 16
   1136 26
    951 22
    852 11
    543 14
    514 5
    472 12
    444 4
    263 25
    251 27
    212 18
    112 17
     99 21
     43 33

#### in_log_WT_Brn1_rep1.bam ####
11318170 44
2708880 1
 557409 0
 133091 42
  83874 31
  65984 41
  46863 40
  40740 39
  34777 37
  33892 34
  32671 33
  30516 38
  27133 36
  26087 32
  22931 11
  22804 35
  12859 14
   9926 28
   7482 22
   6710 21
   6630 17
   5961 24
   5936 25
   5462 18
   4336 2
   2330 9
   1696 12
   1156 16
    684 19

#### in_log_WT_Brn1_rep2.atria.bam ####
8325296 42
2124441 1
 250385 0
 157054 40
  59273 30
  18644 24
  17432 32
  15871 23
  12553 35
  11497 39
  10716 37
   9740 31
   5963 38
   5746 36
   4357 3
   4134 34
   3830 8
   2133 7
   1726 2
   1472 6
    946 15
    884 26
    842 16
    775 22
    734 11
    433 14
    427 5
    400 12
    390 4
    195 25
    190 27
    157 18
     91 21
     81 17
     33 33

#### in_log_WT_Brn1_rep2.bam ####
8901995 44
2308538 1
 410638 0
 104567 42
  65476 31
  52057 41
  34581 40
  30608 39
  27097 37
  26411 34
  25505 33
  23185 38
  21630 36
  20223 32
  18124 11
  17068 35
  10132 14
   8505 28
   5638 22
   5291 17
   5274 21
   4707 24
   4672 25
   4338 18
   3655 2
   1921 9
   1318 12
    918 16
    523 19

#### in_log_WT_Brn1_rep3.atria.bam ####
6347819 42
1383922 1
 206733 0
 168296 40
  45490 30
  14609 23
  14291 24
  13519 32
  10166 35
   9992 39
   9363 37
   7849 31
   5174 38
   4857 36
   4134 3
   3405 34
   3176 8
   1769 7
   1388 2
   1092 6
    967 15
    797 16
    765 26
    615 11
    586 22
    394 14
    351 5
    305 4
    210 12
    132 27
    118 18
    100 25
     50 21
     44 17
     15 33

#### in_log_WT_Brn1_rep3.bam ####
6725442 44
1479637 1
 310826 0
  66145 42
  49972 31
  39259 41
  32714 40
  28157 39
  23644 37
  21493 34
  20390 33
  20344 38
  15347 36
  15105 32
  14601 35
  14134 11
   7805 14
   5954 28
   4924 22
   4256 21
   4092 17
   3980 25
   3641 24
   3523 18
   2999 2
   1567 9
   1128 12
    798 16
    446 19

#### in_log_WT_Brn1_repM.atria.bam ####
25263314 42
6001547 1
 808323 0
 527071 40
 181300 30
  57010 24
  53361 32
  50900 23
  38955 35
  36492 39
  33830 37
  30017 31
  18940 38
  18026 36
  14241 3
  12807 34
  12009 8
   6583 7
   5178 2
   4405 6
   3233 15
   2846 26
   2756 16
   2259 22
   2186 11
   1369 14
   1221 5
   1165 4
   1079 12
    591 27
    542 25
    503 18
    250 21
    237 17
     96 33

#### in_log_WT_Brn1_repM.bam ####
26945601 44
6497037 1
1278786 0
 303835 42
 199311 31
 157301 41
 114212 40
  99388 39
  85516 37
  81884 34
  78661 33
  74102 38
  64129 36
  61370 32
  55126 11
  54358 35
  30898 14
  24383 28
  18060 22
  16210 21
  15966 17
  14614 25
  14322 24
  13381 18
  10926 2
   5832 9
   4178 12
   2839 16
   1682 19

#### in_Q_WT_Brn1_rep1.atria.bam ####
7575923 42
2156379 1
 236725 0
 142143 40
  57353 30
  18051 24
  16972 32
  16731 39
  15637 23
  14602 35
  14485 37
   9703 31
   9020 38
   7813 36
   5128 3
   4982 34
   4371 8
   2100 2
   2071 7
   1414 6
   1361 26
   1156 15
   1139 16
   1004 22
    800 11
    593 14
    545 5
    438 4
    399 12
    272 27
    246 25
    178 18
     95 17
     90 21
     38 33

#### in_Q_WT_Brn1_rep1.bam ####
8021559 44
2340529 1
 375766 0
  96324 42
  62806 31
  55935 40
  47344 41
  47208 39
  36740 37
  33568 38
  32236 34
  26687 33
  24325 35
  22149 32
  20446 11
  20417 36
  10626 14
   8251 22
   7797 28
   6771 21
   6248 17
   6245 25
   5620 24
   5341 18
   4687 2
   2430 9
   1885 12
   1238 16
    758 19

#### in_Q_WT_Brn1_rep2.atria.bam ####
8374819 42
2948443 1
 234612 0
 157176 40
  62526 30
  19988 24
  18202 32
  16949 23
  15424 39
  14443 35
  13600 37
  10408 31
   8157 38
   7601 36
   5198 3
   4832 34
   4503 8
   2157 7
   2066 2
   1587 6
   1176 26
   1147 15
   1111 16
   1006 22
    859 11
    576 14
    548 5
    484 4
    390 12
    242 27
    227 25
    169 18
    113 21
     88 17
     33 33

#### in_Q_WT_Brn1_rep2.bam ####
8895664 44
3195672 1
 372196 0
 103244 42
  68319 31
  51209 41
  48861 40
  42396 39
  34991 37
  31310 34
  30963 38
  27848 33
  22707 32
  22471 35
  21323 11
  21205 36
  11271 14
   8946 28
   7653 22
   6954 21
   6629 17
   6475 25
   5422 18
   5335 24
   4585 2
   2426 9
   1864 12
   1226 16
    725 19

#### in_Q_WT_Brn1_rep3.atria.bam ####
8208677 42
1567566 1
 281494 0
 187892 40
  59173 30
  18455 23
  18257 24
  17781 32
  13576 35
  13428 39
  12352 37
  10305 31
   6996 38
   6323 36
   5270 3
   4637 34
   4012 8
   2475 2
   2364 7
   1365 6
   1175 15
   1064 26
   1064 16
    903 11
    802 22
    585 14
    488 5
    482 4
    257 12
    203 27
    170 25
    144 18
     85 21
     57 17
     24 33

#### in_Q_WT_Brn1_rep3.bam ####
8610317 44
1661701 1
 371212 0
  80249 42
  64518 31
  48035 41
  41716 40
  36891 39
  30420 37
  28838 34
  27089 33
  25983 38
  19947 32
  19888 11
  19405 35
  18496 36
  10393 14
   6682 28
   6220 22
   5365 17
   5316 21
   4626 24
   4625 25
   4349 18
   3494 2
   1924 9
   1321 12
    934 16
    559 19

#### in_Q_WT_Brn1_repM.atria.bam ####
24159430 42
6672387 1
 752871 0
 487229 40
 179054 30
  56317 24
  52937 32
  51040 23
  45838 39
  42612 35
  40441 37
  30411 31
  23917 38
  21786 36
  15667 3
  14416 34
  12876 8
   6551 2
   6541 7
   4397 6
   3624 26
   3444 15
   3354 16
   2830 22
   2550 11
   1725 14
   1557 5
   1412 4
   1025 12
    696 27
    640 25
    500 18
    301 21
    236 17
    106 33

#### in_Q_WT_Brn1_repM.bam ####
25527419 44
7197955 1
1119197 0
 279789 42
 195668 31
 146585 41
 146556 40
 126562 39
 101861 37
  92388 34
  90705 38
  81750 33
  66186 35
  64753 32
  61771 11
  60102 36
  32262 14
  23443 28
  22129 22
  19044 21
  18280 17
  17248 25
  15580 24
  15076 18
  12830 2
   6729 9
   5040 12
   3406 16
   2025 19
```
</details>
<br />
<br />

<a id="tallycalculate-alignments"></a>
## Tally/calculate alignments
<a id="tallycalculate-alignments-1"></a>
### Tally/calculate alignments
<a id="initialize-functions-for-doing-floating-point-arithmetic-etc"></a>
#### Initialize functions for doing floating point arithmetic, etc.
<a id="printed-8"></a>
##### Printed
<details>
<summary><i>Printed</i></summary>

```txt
❯ calc_6f() { awk "BEGIN{ printf \"%.6f\n\", $* }"; }


❯ calc_2f() { awk "BEGIN{ printf \"%.2f\n\", $* }"; }


❯ tally_alignments() {
>     local OPTIND opt type MAPQ file Mito
>     while getopts "t:q:f:m:" opt; do
>         case "${opt}" in
>             t) type="${OPTARG}" ;;
>             q) MAPQ="${OPTARG}" ;;
>             f) file="${OPTARG}" ;;
>             m) Mito="${OPTARG}" ;;
>             *) return 1 ;;
>         esac
>     done
> 
>     [[ -z "${type}" ]] &&
>         {
>             echo "Warning: Argument \"type\" is empty; stopping the function"
>             return 1
>         }
>     [[ -z "${MAPQ}" ]] && MAPQ=0
>     [[ -z "${file}" ]] &&
>         {
>             echo "Warning: Argument \"file\" is empty; stopping the function"
>             return 1
>         }
>     [[ -z "${Mito}" ]] && Mito=FALSE
> 
>     debug=FALSE
>     [[ "${debug}" == TRUE ]] &&
>         {
>             type="all"
>             MAPQ=0
>             file="${i}"
>             Mito=FALSE
> 
>             echo "${type}"
>             echo "${MAPQ}"
>             echo "${file}"
>             echo "${Mito}"
>         }
> 
>     #  All alignments
>     [[ "${type}" == "all" ]] &&
>         {
>             if [[ "${Mito}" == TRUE || "${Mito}" == T ]]; then
>                 tally="$(
>                     samtools view -c \
>                         -F 4 \
>                         -q "${MAPQ}" \
>                         "${file}"
>                 )"
>             elif [[ "${Mito}" == FALSE ]]; then
>                 tally="$(
>                     samtools view -c \
>                         -F 4 \
>                         -q "${MAPQ}" \
>                         "${file}" \
>                         I II III IV V VI VII VIII IX X \
>                         XI XII XIII XIV XV XVI \
>                         SP_II_TG SP_I SP_II SP_III SP_MTR
>                 )"
>             fi
> 
>             echo "${tally}" && return 0
>         }
> 
> 
>     #  S. cerevisiae alignments only
>     [[ "${type}" == "SC" || "${type}" == "sc" || "${type}" == "Sc" ]] &&
>         {
>             if [[ "${Mito}" == TRUE ]]; then
>                 tally="$(
>                     samtools view -c \
>                         -F 4 \
>                         -q "${MAPQ}" \
>                         "${file}" \
>                         I II III IV V VI VII VIII IX X \
>                         XI XII XIII XIV XV XVI Mito
>                 )"
>             elif [[ "${Mito}" == FALSE ]]; then
>                 tally="$(
>                     samtools view -c \
>                         -F 4 \
>                         -q "${MAPQ}" \
>                         "${file}" \
>                         I II III IV V VI VII VIII IX X \
>                         XI XII XIII XIV XV XVI
>                 )"
>             fi
> 
>             echo "${tally}" && return 0
>         }
> 
>     #  S. pombe alignments only
>     [[ "${type}" == "SP" || "${type}" == "sp" || "${type}" == "Sp" ]] &&
>         {
>             if [[ "${Mito}" == TRUE ]]; then
>                 tally="$(
>                     samtools view -c \
>                         -F 4 \
>                         -q "${MAPQ}" \
>                         "${file}" \
>                         SP_II_TG SP_I SP_II SP_III SP_MTR SP_Mito
>                 )"
>             elif [[ "${Mito}" == FALSE ]]; then
>                 tally="$(
>                     samtools view -c \
>                         -F 4 \
>                         -q "${MAPQ}" \
>                         "${file}" \
>                         SP_II_TG SP_I SP_II SP_III SP_MTR
>                 )"
>             fi
> 
>             echo "${tally}" && return 0
>         }
> }


❯ transpose_values() {
>     matrix="$(
>         awk '
>             {
>                 #  Loop through each line and each field (column) in the input
>                 #+ file, storing the value of each field in the array "a" indexed
>                 #+ by the column number "i" and the row number "NR"
>                 for (i = 1; i <= NF; i++) {
>                     a[i, NR] = $i
>                 }
>             }
>             END {
>                 #  Loop through the stored values in the array "a" and print the
>                 #+ transposed values
>                 for (i = 1; i <= NF; i++) {
>                     for (j = 1; j <= NR; j++) {
>                         printf "%s", a[i, j]
>                         if (j < NR) {
>                             printf " "
>                         }
>                     }
>                     print ""
>                 }
>             }
>         ' "${1}"
>     )"
> 
>     echo "${matrix}" && return 0
> }


❯ export -f calc_6f
❯ export -f calc_2f
❯ export -f tally_alignments
❯ export -f transpose_values
```
</details>
<br />

<a id="get-situated-then-initialize-arrays-variables-etc"></a>
#### Get situated, then initialize arrays, variables, etc.
<a id="printed-9"></a>
##### Printed
<details>
<summary><i>Printed</i></summary>

```txt
❯ module purge


❯ ml SAMtools/1.16.1-GCC-11.2.0 Bowtie2/2.4.4-GCC-11.2.0


❯ cd "${HOME}/tsukiyamalab/Kris/2023_rDNA/results" ||
>     echo "cd'ing failed; check on this..."


❯ p_data="${HOME}/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams"


❯ typeset -a bams=(
>     "${p_data}/in_log_WT_Brn1_rep1.atria.bam"
>     "${p_data}/Ch_log_WT_Brn1_rep1.atria.bam"
>     "${p_data}/in_log_WT_Brn1_rep2.atria.bam"
>     "${p_data}/Ch_log_WT_Brn1_rep2.atria.bam"
>     "${p_data}/in_log_WT_Brn1_rep3.atria.bam"
>     "${p_data}/Ch_log_WT_Brn1_rep3.atria.bam"
>     "${p_data}/in_log_WT_Brn1_repM.atria.bam"
>     "${p_data}/Ch_log_WT_Brn1_repM.atria.bam"
>     "${p_data}/in_Q_WT_Brn1_rep1.atria.bam"
>     "${p_data}/Ch_Q_WT_Brn1_rep1.atria.bam"
>     "${p_data}/in_Q_WT_Brn1_rep2.atria.bam"
>     "${p_data}/Ch_Q_WT_Brn1_rep2.atria.bam"
>     "${p_data}/in_Q_WT_Brn1_rep3.atria.bam"
>     "${p_data}/Ch_Q_WT_Brn1_rep3.atria.bam"
>     "${p_data}/in_Q_WT_Brn1_repM.atria.bam"
>     "${p_data}/Ch_Q_WT_Brn1_repM.atria.bam"
> )


❯ run_check=TRUE


❯ [[ "${run_check}" == TRUE ]] && \
>     for i in "${bams[@]}"; do ls -lhaFG "${i}"; done
-rw-rw---- 1 kalavatt 348M Jul 16 21:21 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_log_WT_Brn1_rep1.atria.bam
-rw-rw---- 1 kalavatt 132M Jul 16 15:59 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_log_WT_Brn1_rep1.atria.bam
-rw-rw---- 1 kalavatt 279M Jul 16 21:22 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_log_WT_Brn1_rep2.atria.bam
-rw-rw---- 1 kalavatt 585M Jul 16 20:31 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_log_WT_Brn1_rep2.atria.bam
-rw-rw---- 1 kalavatt 199M Jul 16 21:23 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_log_WT_Brn1_rep3.atria.bam
-rw-rw---- 1 kalavatt 180M Jul 16 20:31 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_log_WT_Brn1_rep3.atria.bam
-rw-rw---- 1 kalavatt 781M Jul 16 22:20 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_log_WT_Brn1_repM.atria.bam
-rw-rw---- 1 kalavatt 863M Jul 16 21:17 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_log_WT_Brn1_repM.atria.bam
-rw-rw---- 1 kalavatt 257M Jul 16 22:21 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_Q_WT_Brn1_rep1.atria.bam
-rw-rw---- 1 kalavatt 172M Jul 16 21:17 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_Q_WT_Brn1_rep1.atria.bam
-rw-rw---- 1 kalavatt 294M Jul 16 22:22 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_Q_WT_Brn1_rep2.atria.bam
-rw-rw---- 1 kalavatt 150M Jul 16 21:18 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_Q_WT_Brn1_rep2.atria.bam
-rw-rw---- 1 kalavatt 249M Jul 16 22:22 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_Q_WT_Brn1_rep3.atria.bam
-rw-rw---- 1 kalavatt 209M Jul 16 21:18 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_Q_WT_Brn1_rep3.atria.bam
-rw-rw---- 1 kalavatt 756M Jul 16 22:50 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/in_Q_WT_Brn1_repM.atria.bam
-rw-rw---- 1 kalavatt 506M Jul 16 21:20 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Ch_Q_WT_Brn1_repM.atria.bam
```
</details>
<br />

<a id="generate-tab-separated-table-of-alignment-talliescalculations"></a>
#### Generate tab-separated table of alignment tallies/calculations
<a id="printed-10"></a>
##### Printed
<details>
<summary><i>Printed</i></summary>

```txt
#### in_log_WT_Brn1_rep1.atria.bam ####

            in_log_WT_Brn1_rep1
            ------------------------------------------------------------
                      "${tallies[${unmapped}]}"  333904
                        "${tallies[${mapped}]}"  13549404

               "${tallies[${total_with_Mito}]}"  13883308
            "${tallies[${total_SC_with_Mito}]}"  13258095
            "${tallies[${total_SP_with_Mito}]}"  291309

               "${tallies[${total_sans_Mito}]}"  12972893
            "${tallies[${total_SC_sans_Mito}]}"  12706090
            "${tallies[${total_SP_sans_Mito}]}"  266803

                 "${tallies[${total_SC_Mito}]}"  552005
                 "${tallies[${total_SP_Mito}]}"  24506

                       "${tallies[${total_1}]}"  12959145
                       "${tallies[${total_2}]}"  10473556
                       "${tallies[${total_3}]}"  10472049
                      "${tallies[${total_23}]}"  10454533
                      "${tallies[${total_30}]}"  10413321
                      "${tallies[${total_40}]}"  10265586
                      "${tallies[${total_42}]}"  10074757

                          "${tallies[${SC_0}]}"  12706090
                          "${tallies[${SC_1}]}"  12692504
                          "${tallies[${SC_2}]}"  10300915
                          "${tallies[${SC_3}]}"  10299428
                         "${tallies[${SC_23}]}"  10282048
                         "${tallies[${SC_30}]}"  10240927
                         "${tallies[${SC_40}]}"  10095174
                         "${tallies[${SC_42}]}"  9905209

                          "${tallies[${SP_0}]}"  266803
                          "${tallies[${SP_1}]}"  266641
                          "${tallies[${SP_2}]}"  172641
                          "${tallies[${SP_3}]}"  172621
                         "${tallies[${SP_23}]}"  172485
                         "${tallies[${SP_30}]}"  172394
                         "${tallies[${SP_40}]}"  170412
                         "${tallies[${SP_42}]}"  169548

                    "${tallies[${CC_SP2SC_0}]}"  0.020998
                    "${tallies[${CC_SP2SC_1}]}"  0.021008
                    "${tallies[${CC_SP2SC_2}]}"  0.016760
                    "${tallies[${CC_SP2SC_3}]}"  0.016760
                   "${tallies[${CC_SP2SC_23}]}"  0.016775
                   "${tallies[${CC_SP2SC_30}]}"  0.016834
                   "${tallies[${CC_SP2SC_40}]}"  0.016881
                   "${tallies[${CC_SP2SC_42}]}"  0.017117

#### Ch_log_WT_Brn1_rep1.atria.bam ####

            Ch_log_WT_Brn1_rep1
            ------------------------------------------------------------
                      "${tallies[${unmapped}]}"  499911
                        "${tallies[${mapped}]}"  5752746

               "${tallies[${total_with_Mito}]}"  6252657
            "${tallies[${total_SC_with_Mito}]}"  5609097
            "${tallies[${total_SP_with_Mito}]}"  143649

               "${tallies[${total_sans_Mito}]}"  5740916
            "${tallies[${total_SC_sans_Mito}]}"  5597556
            "${tallies[${total_SP_sans_Mito}]}"  143360

                 "${tallies[${total_SC_Mito}]}"  11541
                 "${tallies[${total_SP_Mito}]}"  289

                       "${tallies[${total_1}]}"  5722062
                       "${tallies[${total_2}]}"  3098907
                       "${tallies[${total_3}]}"  3096507
                      "${tallies[${total_23}]}"  3085337
                      "${tallies[${total_30}]}"  3065407
                      "${tallies[${total_40}]}"  2998196
                      "${tallies[${total_42}]}"  2934584

                          "${tallies[${SC_0}]}"  5597556
                          "${tallies[${SC_1}]}"  5578873
                          "${tallies[${SC_2}]}"  3059048
                          "${tallies[${SC_3}]}"  3056672
                         "${tallies[${SC_23}]}"  3045575
                         "${tallies[${SC_30}]}"  3025669
                         "${tallies[${SC_40}]}"  2959163
                         "${tallies[${SC_42}]}"  2895786

                          "${tallies[${SP_0}]}"  143360
                          "${tallies[${SP_1}]}"  143189
                          "${tallies[${SP_2}]}"  39859
                          "${tallies[${SP_3}]}"  39835
                         "${tallies[${SP_23}]}"  39762
                         "${tallies[${SP_30}]}"  39738
                         "${tallies[${SP_40}]}"  39033
                         "${tallies[${SP_42}]}"  38798

                    "${tallies[${CC_SP2SC_0}]}"  0.025611
                    "${tallies[${CC_SP2SC_1}]}"  0.025666
                    "${tallies[${CC_SP2SC_2}]}"  0.013030
                    "${tallies[${CC_SP2SC_3}]}"  0.013032
                   "${tallies[${CC_SP2SC_23}]}"  0.013056
                   "${tallies[${CC_SP2SC_30}]}"  0.013134
                   "${tallies[${CC_SP2SC_40}]}"  0.013191
                   "${tallies[${CC_SP2SC_42}]}"  0.013398

#### in_log_WT_Brn1_rep2.atria.bam ####

            in_log_WT_Brn1_rep2
            ------------------------------------------------------------
                      "${tallies[${unmapped}]}"  236372
                        "${tallies[${mapped}]}"  10812469

               "${tallies[${total_with_Mito}]}"  11048841
            "${tallies[${total_SC_with_Mito}]}"  10571481
            "${tallies[${total_SP_with_Mito}]}"  240988

               "${tallies[${total_sans_Mito}]}"  10384551
            "${tallies[${total_SC_sans_Mito}]}"  10162156
            "${tallies[${total_SP_sans_Mito}]}"  222395

                 "${tallies[${total_SC_Mito}]}"  409325
                 "${tallies[${total_SP_Mito}]}"  18593

                       "${tallies[${total_1}]}"  10373037
                       "${tallies[${total_2}]}"  8254098
                       "${tallies[${total_3}]}"  8252807
                      "${tallies[${total_23}]}"  8238780
                      "${tallies[${total_30}]}"  8206403
                      "${tallies[${total_40}]}"  8090425
                      "${tallies[${total_42}]}"  7941054

                          "${tallies[${SC_0}]}"  10162156
                          "${tallies[${SC_1}]}"  10150752
                          "${tallies[${SC_2}]}"  8112719
                          "${tallies[${SC_3}]}"  8111449
                         "${tallies[${SC_23}]}"  8097546
                         "${tallies[${SC_30}]}"  8065242
                         "${tallies[${SC_40}]}"  7950863
                         "${tallies[${SC_42}]}"  7802230

                          "${tallies[${SP_0}]}"  222395
                          "${tallies[${SP_1}]}"  222285
                          "${tallies[${SP_2}]}"  141379
                          "${tallies[${SP_3}]}"  141358
                         "${tallies[${SP_23}]}"  141234
                         "${tallies[${SP_30}]}"  141161
                         "${tallies[${SP_40}]}"  139562
                         "${tallies[${SP_42}]}"  138824

                    "${tallies[${CC_SP2SC_0}]}"  0.021885
                    "${tallies[${CC_SP2SC_1}]}"  0.021898
                    "${tallies[${CC_SP2SC_2}]}"  0.017427
                    "${tallies[${CC_SP2SC_3}]}"  0.017427
                   "${tallies[${CC_SP2SC_23}]}"  0.017442
                   "${tallies[${CC_SP2SC_30}]}"  0.017502
                   "${tallies[${CC_SP2SC_40}]}"  0.017553
                   "${tallies[${CC_SP2SC_42}]}"  0.017793

#### Ch_log_WT_Brn1_rep2.atria.bam ####

            Ch_log_WT_Brn1_rep2
            ------------------------------------------------------------
                      "${tallies[${unmapped}]}"  2734535
                        "${tallies[${mapped}]}"  27454714

               "${tallies[${total_with_Mito}]}"  30189249
            "${tallies[${total_SC_with_Mito}]}"  26743282
            "${tallies[${total_SP_with_Mito}]}"  711432

               "${tallies[${total_sans_Mito}]}"  27353401
            "${tallies[${total_SC_sans_Mito}]}"  26644138
            "${tallies[${total_SP_sans_Mito}]}"  709263

                 "${tallies[${total_SC_Mito}]}"  99144
                 "${tallies[${total_SP_Mito}]}"  2169

                       "${tallies[${total_1}]}"  27262603
                       "${tallies[${total_2}]}"  14684393
                       "${tallies[${total_3}]}"  14673732
                      "${tallies[${total_23}]}"  14620998
                      "${tallies[${total_30}]}"  14525548
                      "${tallies[${total_40}]}"  14211401
                      "${tallies[${total_42}]}"  13913155

                          "${tallies[${SC_0}]}"  26644138
                          "${tallies[${SC_1}]}"  26554208
                          "${tallies[${SC_2}]}"  14478204
                          "${tallies[${SC_3}]}"  14467671
                         "${tallies[${SC_23}]}"  14415347
                         "${tallies[${SC_30}]}"  14320032
                         "${tallies[${SC_40}]}"  14009964
                         "${tallies[${SC_42}]}"  13713002

                          "${tallies[${SP_0}]}"  709263
                          "${tallies[${SP_1}]}"  708395
                          "${tallies[${SP_2}]}"  206189
                          "${tallies[${SP_3}]}"  206061
                         "${tallies[${SP_23}]}"  205651
                         "${tallies[${SP_30}]}"  205516
                         "${tallies[${SP_40}]}"  201437
                         "${tallies[${SP_42}]}"  200153

                    "${tallies[${CC_SP2SC_0}]}"  0.026620
                    "${tallies[${CC_SP2SC_1}]}"  0.026677
                    "${tallies[${CC_SP2SC_2}]}"  0.014241
                    "${tallies[${CC_SP2SC_3}]}"  0.014243
                   "${tallies[${CC_SP2SC_23}]}"  0.014266
                   "${tallies[${CC_SP2SC_30}]}"  0.014352
                   "${tallies[${CC_SP2SC_40}]}"  0.014378
                   "${tallies[${CC_SP2SC_42}]}"  0.014596

#### in_log_WT_Brn1_rep3.atria.bam ####

            in_log_WT_Brn1_rep3
            ------------------------------------------------------------
                      "${tallies[${unmapped}]}"  195541
                        "${tallies[${mapped}]}"  8066952

               "${tallies[${total_with_Mito}]}"  8262493
            "${tallies[${total_SC_with_Mito}]}"  8027800
            "${tallies[${total_SP_with_Mito}]}"  39152

               "${tallies[${total_sans_Mito}]}"  7587636
            "${tallies[${total_SC_sans_Mito}]}"  7548484
            "${tallies[${total_SP_sans_Mito}]}"  39152

                 "${tallies[${total_SC_Mito}]}"  479316
                 "${tallies[${total_SP_Mito}]}"  0

                       "${tallies[${total_1}]}"  7578953
                       "${tallies[${total_2}]}"  6200891
                       "${tallies[${total_3}]}"  6199862
                      "${tallies[${total_23}]}"  6188627
                      "${tallies[${total_30}]}"  6162439
                      "${tallies[${total_40}]}"  6074756
                      "${tallies[${total_42}]}"  5917888

                          "${tallies[${SC_0}]}"  7548484
                          "${tallies[${SC_1}]}"  7539870
                          "${tallies[${SC_2}]}"  6200843
                          "${tallies[${SC_3}]}"  6199830
                         "${tallies[${SC_23}]}"  6188614
                         "${tallies[${SC_30}]}"  6162429
                         "${tallies[${SC_40}]}"  6074756
                         "${tallies[${SC_42}]}"  5917888

                          "${tallies[${SP_0}]}"  39152
                          "${tallies[${SP_1}]}"  39083
                          "${tallies[${SP_2}]}"  48
                          "${tallies[${SP_3}]}"  32
                         "${tallies[${SP_23}]}"  13
                         "${tallies[${SP_30}]}"  10
                         "${tallies[${SP_40}]}"  0
                         "${tallies[${SP_42}]}"  0

                    "${tallies[${CC_SP2SC_0}]}"  0.005187
                    "${tallies[${CC_SP2SC_1}]}"  0.005184
                    "${tallies[${CC_SP2SC_2}]}"  0.000008
                    "${tallies[${CC_SP2SC_3}]}"  0.000005
                   "${tallies[${CC_SP2SC_23}]}"  0.000002
                   "${tallies[${CC_SP2SC_30}]}"  0.000002
                   "${tallies[${CC_SP2SC_40}]}"  0.000000
                   "${tallies[${CC_SP2SC_42}]}"  0.000000

#### Ch_log_WT_Brn1_rep3.atria.bam ####

            Ch_log_WT_Brn1_rep3
            ------------------------------------------------------------
                      "${tallies[${unmapped}]}"  316828
                        "${tallies[${mapped}]}"  7851602

               "${tallies[${total_with_Mito}]}"  8168430
            "${tallies[${total_SC_with_Mito}]}"  7768424
            "${tallies[${total_SP_with_Mito}]}"  83178

               "${tallies[${total_sans_Mito}]}"  7698275
            "${tallies[${total_SC_sans_Mito}]}"  7615097
            "${tallies[${total_SP_sans_Mito}]}"  83178

                 "${tallies[${total_SC_Mito}]}"  153327
                 "${tallies[${total_SP_Mito}]}"  0

                       "${tallies[${total_1}]}"  7677910
                       "${tallies[${total_2}]}"  4984373
                       "${tallies[${total_3}]}"  4982252
                      "${tallies[${total_23}]}"  4968750
                      "${tallies[${total_30}]}"  4941082
                      "${tallies[${total_40}]}"  4847398
                      "${tallies[${total_42}]}"  4728066

                          "${tallies[${SC_0}]}"  7615097
                          "${tallies[${SC_1}]}"  7594867
                          "${tallies[${SC_2}]}"  4984274
                          "${tallies[${SC_3}]}"  4982179
                         "${tallies[${SC_23}]}"  4968721
                         "${tallies[${SC_30}]}"  4941058
                         "${tallies[${SC_40}]}"  4847395
                         "${tallies[${SC_42}]}"  4728066

                          "${tallies[${SP_0}]}"  83178
                          "${tallies[${SP_1}]}"  83043
                          "${tallies[${SP_2}]}"  99
                          "${tallies[${SP_3}]}"  73
                         "${tallies[${SP_23}]}"  29
                         "${tallies[${SP_30}]}"  24
                         "${tallies[${SP_40}]}"  3
                         "${tallies[${SP_42}]}"  0

                    "${tallies[${CC_SP2SC_0}]}"  0.010923
                    "${tallies[${CC_SP2SC_1}]}"  0.010934
                    "${tallies[${CC_SP2SC_2}]}"  0.000020
                    "${tallies[${CC_SP2SC_3}]}"  0.000015
                   "${tallies[${CC_SP2SC_23}]}"  0.000006
                   "${tallies[${CC_SP2SC_30}]}"  0.000005
                   "${tallies[${CC_SP2SC_40}]}"  0.000001
                   "${tallies[${CC_SP2SC_42}]}"  0.000000

#### in_log_WT_Brn1_repM.atria.bam ####

            in_log_WT_Brn1_repM
            ------------------------------------------------------------
                      "${tallies[${unmapped}]}"  765891
                        "${tallies[${mapped}]}"  32428751

               "${tallies[${total_with_Mito}]}"  33194642
            "${tallies[${total_SC_with_Mito}]}"  31856761
            "${tallies[${total_SP_with_Mito}]}"  571990

               "${tallies[${total_sans_Mito}]}"  30945015
            "${tallies[${total_SC_sans_Mito}]}"  30416127
            "${tallies[${total_SP_sans_Mito}]}"  528888

                 "${tallies[${total_SC_Mito}]}"  1440634
                 "${tallies[${total_SP_Mito}]}"  43102

                       "${tallies[${total_1}]}"  30911132
                       "${tallies[${total_2}]}"  24928560
                       "${tallies[${total_3}]}"  24924677
                      "${tallies[${total_23}]}"  24881929
                      "${tallies[${total_30}]}"  24782146
                      "${tallies[${total_40}]}"  24430838
                      "${tallies[${total_42}]}"  23933786

                          "${tallies[${SC_0}]}"  30416127
                          "${tallies[${SC_1}]}"  30382591
                          "${tallies[${SC_2}]}"  24614485
                          "${tallies[${SC_3}]}"  24610669
                         "${tallies[${SC_23}]}"  24568200
                         "${tallies[${SC_30}]}"  24468582
                         "${tallies[${SC_40}]}"  24120863
                         "${tallies[${SC_42}]}"  23625412

                          "${tallies[${SP_0}]}"  528888
                          "${tallies[${SP_1}]}"  528541
                          "${tallies[${SP_2}]}"  314075
                          "${tallies[${SP_3}]}"  314008
                         "${tallies[${SP_23}]}"  313729
                         "${tallies[${SP_30}]}"  313564
                         "${tallies[${SP_40}]}"  309975
                         "${tallies[${SP_42}]}"  308374

                    "${tallies[${CC_SP2SC_0}]}"  0.017388
                    "${tallies[${CC_SP2SC_1}]}"  0.017396
                    "${tallies[${CC_SP2SC_2}]}"  0.012760
                    "${tallies[${CC_SP2SC_3}]}"  0.012759
                   "${tallies[${CC_SP2SC_23}]}"  0.012770
                   "${tallies[${CC_SP2SC_30}]}"  0.012815
                   "${tallies[${CC_SP2SC_40}]}"  0.012851
                   "${tallies[${CC_SP2SC_42}]}"  0.013053

#### Ch_log_WT_Brn1_repM.atria.bam ####

            Ch_log_WT_Brn1_repM
            ------------------------------------------------------------
                      "${tallies[${unmapped}]}"  3551257
                        "${tallies[${mapped}]}"  41059079

               "${tallies[${total_with_Mito}]}"  44610336
            "${tallies[${total_SC_with_Mito}]}"  40120529
            "${tallies[${total_SP_with_Mito}]}"  938550

               "${tallies[${total_sans_Mito}]}"  40792626
            "${tallies[${total_SC_sans_Mito}]}"  39856529
            "${tallies[${total_SP_sans_Mito}]}"  936097

                 "${tallies[${total_SC_Mito}]}"  264000
                 "${tallies[${total_SP_Mito}]}"  2453

                       "${tallies[${total_1}]}"  40662442
                       "${tallies[${total_2}]}"  22767532
                       "${tallies[${total_3}]}"  22752335
                      "${tallies[${total_23}]}"  22674868
                      "${tallies[${total_30}]}"  22531987
                      "${tallies[${total_40}]}"  22057002
                      "${tallies[${total_42}]}"  21575776

                          "${tallies[${SC_0}]}"  39856529
                          "${tallies[${SC_1}]}"  39727524
                          "${tallies[${SC_2}]}"  22521397
                          "${tallies[${SC_3}]}"  22506372
                         "${tallies[${SC_23}]}"  22429428
                         "${tallies[${SC_30}]}"  22286712
                         "${tallies[${SC_40}]}"  21816532
                         "${tallies[${SC_42}]}"  21336825

                          "${tallies[${SP_0}]}"  936097
                          "${tallies[${SP_1}]}"  934918
                          "${tallies[${SP_2}]}"  246135
                          "${tallies[${SP_3}]}"  245963
                         "${tallies[${SP_23}]}"  245440
                         "${tallies[${SP_30}]}"  245275
                         "${tallies[${SP_40}]}"  240470
                         "${tallies[${SP_42}]}"  238951

                    "${tallies[${CC_SP2SC_0}]}"  0.023487
                    "${tallies[${CC_SP2SC_1}]}"  0.023533
                    "${tallies[${CC_SP2SC_2}]}"  0.010929
                    "${tallies[${CC_SP2SC_3}]}"  0.010929
                   "${tallies[${CC_SP2SC_23}]}"  0.010943
                   "${tallies[${CC_SP2SC_30}]}"  0.011005
                   "${tallies[${CC_SP2SC_40}]}"  0.011022
                   "${tallies[${CC_SP2SC_42}]}"  0.011199

#### in_Q_WT_Brn1_rep1.atria.bam ####

            in_Q_WT_Brn1_rep1
            ------------------------------------------------------------
                      "${tallies[${unmapped}]}"  220686
                        "${tallies[${mapped}]}"  10099271

               "${tallies[${total_with_Mito}]}"  10319957
            "${tallies[${total_SC_with_Mito}]}"  9902653
            "${tallies[${total_SP_with_Mito}]}"  196618

               "${tallies[${total_sans_Mito}]}"  9062225
            "${tallies[${total_SC_sans_Mito}]}"  8879315
            "${tallies[${total_SP_sans_Mito}]}"  182910

                 "${tallies[${total_SC_Mito}]}"  1023338
                 "${tallies[${total_SP_Mito}]}"  13708

                       "${tallies[${total_1}]}"  9051915
                       "${tallies[${total_2}]}"  6908730
                       "${tallies[${total_3}]}"  6907543
                      "${tallies[${total_23}]}"  6895343
                      "${tallies[${total_30}]}"  6868153
                      "${tallies[${total_40}]}"  6768007
                      "${tallies[${total_42}]}"  6644434

                          "${tallies[${SC_0}]}"  8879315
                          "${tallies[${SC_1}]}"  8869108
                          "${tallies[${SC_2}]}"  6807940
                          "${tallies[${SC_3}]}"  6806770
                         "${tallies[${SC_23}]}"  6794672
                         "${tallies[${SC_30}]}"  6767541
                         "${tallies[${SC_40}]}"  6668609
                         "${tallies[${SC_42}]}"  6545581

                          "${tallies[${SP_0}]}"  182910
                          "${tallies[${SP_1}]}"  182807
                          "${tallies[${SP_2}]}"  100790
                          "${tallies[${SP_3}]}"  100773
                         "${tallies[${SP_23}]}"  100671
                         "${tallies[${SP_30}]}"  100612
                         "${tallies[${SP_40}]}"  99398
                         "${tallies[${SP_42}]}"  98853

                    "${tallies[${CC_SP2SC_0}]}"  0.020600
                    "${tallies[${CC_SP2SC_1}]}"  0.020612
                    "${tallies[${CC_SP2SC_2}]}"  0.014805
                    "${tallies[${CC_SP2SC_3}]}"  0.014805
                   "${tallies[${CC_SP2SC_23}]}"  0.014816
                   "${tallies[${CC_SP2SC_30}]}"  0.014867
                   "${tallies[${CC_SP2SC_40}]}"  0.014905
                   "${tallies[${CC_SP2SC_42}]}"  0.015102

#### Ch_Q_WT_Brn1_rep1.atria.bam ####

            Ch_Q_WT_Brn1_rep1
            ------------------------------------------------------------
                      "${tallies[${unmapped}]}"  1739181
                        "${tallies[${mapped}]}"  5492116

               "${tallies[${total_with_Mito}]}"  7231297
            "${tallies[${total_SC_with_Mito}]}"  5253236
            "${tallies[${total_SP_with_Mito}]}"  238880

               "${tallies[${total_sans_Mito}]}"  5369372
            "${tallies[${total_SC_sans_Mito}]}"  5137368
            "${tallies[${total_SP_sans_Mito}]}"  232004

                 "${tallies[${total_SC_Mito}]}"  115868
                 "${tallies[${total_SP_Mito}]}"  6876

                       "${tallies[${total_1}]}"  5359911
                       "${tallies[${total_2}]}"  3745001
                       "${tallies[${total_3}]}"  3743611
                      "${tallies[${total_23}]}"  3734792
                      "${tallies[${total_30}]}"  3717664
                      "${tallies[${total_40}]}"  3658456
                      "${tallies[${total_42}]}"  3579279

                          "${tallies[${SC_0}]}"  5137368
                          "${tallies[${SC_1}]}"  5128233
                          "${tallies[${SC_2}]}"  3597499
                          "${tallies[${SC_3}]}"  3596131
                         "${tallies[${SC_23}]}"  3587470
                         "${tallies[${SC_30}]}"  3570424
                         "${tallies[${SC_40}]}"  3513649
                         "${tallies[${SC_42}]}"  3435400

                          "${tallies[${SP_0}]}"  232004
                          "${tallies[${SP_1}]}"  231678
                          "${tallies[${SP_2}]}"  147502
                          "${tallies[${SP_3}]}"  147480
                         "${tallies[${SP_23}]}"  147322
                         "${tallies[${SP_30}]}"  147240
                         "${tallies[${SP_40}]}"  144807
                         "${tallies[${SP_42}]}"  143879

                    "${tallies[${CC_SP2SC_0}]}"  0.045160
                    "${tallies[${CC_SP2SC_1}]}"  0.045177
                    "${tallies[${CC_SP2SC_2}]}"  0.041001
                    "${tallies[${CC_SP2SC_3}]}"  0.041011
                   "${tallies[${CC_SP2SC_23}]}"  0.041066
                   "${tallies[${CC_SP2SC_30}]}"  0.041239
                   "${tallies[${CC_SP2SC_40}]}"  0.041213
                   "${tallies[${CC_SP2SC_42}]}"  0.041881

#### in_Q_WT_Brn1_rep2.atria.bam ####

            in_Q_WT_Brn1_rep2
            ------------------------------------------------------------
                      "${tallies[${unmapped}]}"  217378
                        "${tallies[${mapped}]}"  11713482

               "${tallies[${total_with_Mito}]}"  11930860
            "${tallies[${total_SC_with_Mito}]}"  11444841
            "${tallies[${total_SP_with_Mito}]}"  268641

               "${tallies[${total_sans_Mito}]}"  10908066
            "${tallies[${total_SC_sans_Mito}]}"  10657329
            "${tallies[${total_SP_sans_Mito}]}"  250737

                 "${tallies[${total_SC_Mito}]}"  787512
                 "${tallies[${total_SP_Mito}]}"  17904

                       "${tallies[${total_1}]}"  10895153
                       "${tallies[${total_2}]}"  7956792
                       "${tallies[${total_3}]}"  7955394
                      "${tallies[${total_23}]}"  7940826
                      "${tallies[${total_30}]}"  7908496
                      "${tallies[${total_40}]}"  7792753
                      "${tallies[${total_42}]}"  7649588

                          "${tallies[${SC_0}]}"  10657329
                          "${tallies[${SC_1}]}"  10644568
                          "${tallies[${SC_2}]}"  7822379
                          "${tallies[${SC_3}]}"  7821009
                         "${tallies[${SC_23}]}"  7806582
                         "${tallies[${SC_30}]}"  7774331
                         "${tallies[${SC_40}]}"  7660119
                         "${tallies[${SC_42}]}"  7517642

                          "${tallies[${SP_0}]}"  250737
                          "${tallies[${SP_1}]}"  250585
                          "${tallies[${SP_2}]}"  134413
                          "${tallies[${SP_3}]}"  134385
                         "${tallies[${SP_23}]}"  134244
                         "${tallies[${SP_30}]}"  134165
                         "${tallies[${SP_40}]}"  132634
                         "${tallies[${SP_42}]}"  131946

                    "${tallies[${CC_SP2SC_0}]}"  0.023527
                    "${tallies[${CC_SP2SC_1}]}"  0.023541
                    "${tallies[${CC_SP2SC_2}]}"  0.017183
                    "${tallies[${CC_SP2SC_3}]}"  0.017183
                   "${tallies[${CC_SP2SC_23}]}"  0.017196
                   "${tallies[${CC_SP2SC_30}]}"  0.017257
                   "${tallies[${CC_SP2SC_40}]}"  0.017315
                   "${tallies[${CC_SP2SC_42}]}"  0.017552

#### Ch_Q_WT_Brn1_rep2.atria.bam ####

            Ch_Q_WT_Brn1_rep2
            ------------------------------------------------------------
                      "${tallies[${unmapped}]}"  1066876
                        "${tallies[${mapped}]}"  5421690

               "${tallies[${total_with_Mito}]}"  6488566
            "${tallies[${total_SC_with_Mito}]}"  5202988
            "${tallies[${total_SP_with_Mito}]}"  218702

               "${tallies[${total_sans_Mito}]}"  5350248
            "${tallies[${total_SC_sans_Mito}]}"  5137697
            "${tallies[${total_SP_sans_Mito}]}"  212551

                 "${tallies[${total_SC_Mito}]}"  65291
                 "${tallies[${total_SP_Mito}]}"  6151

                       "${tallies[${total_1}]}"  5341736
                       "${tallies[${total_2}]}"  3442540
                       "${tallies[${total_3}]}"  3441459
                      "${tallies[${total_23}]}"  3433073
                      "${tallies[${total_30}]}"  3415828
                      "${tallies[${total_40}]}"  3360434
                      "${tallies[${total_42}]}"  3286560

                          "${tallies[${SC_0}]}"  5137697
                          "${tallies[${SC_1}]}"  5129454
                          "${tallies[${SC_2}]}"  3321989
                          "${tallies[${SC_3}]}"  3320927
                         "${tallies[${SC_23}]}"  3312680
                         "${tallies[${SC_30}]}"  3295510
                         "${tallies[${SC_40}]}"  3242161
                         "${tallies[${SC_42}]}"  3169083

                          "${tallies[${SP_0}]}"  212551
                          "${tallies[${SP_1}]}"  212282
                          "${tallies[${SP_2}]}"  120551
                          "${tallies[${SP_3}]}"  120532
                         "${tallies[${SP_23}]}"  120393
                         "${tallies[${SP_30}]}"  120318
                         "${tallies[${SP_40}]}"  118273
                         "${tallies[${SP_42}]}"  117477

                    "${tallies[${CC_SP2SC_0}]}"  0.041371
                    "${tallies[${CC_SP2SC_1}]}"  0.041385
                    "${tallies[${CC_SP2SC_2}]}"  0.036289
                    "${tallies[${CC_SP2SC_3}]}"  0.036295
                   "${tallies[${CC_SP2SC_23}]}"  0.036343
                   "${tallies[${CC_SP2SC_30}]}"  0.036510
                   "${tallies[${CC_SP2SC_40}]}"  0.036480
                   "${tallies[${CC_SP2SC_42}]}"  0.037070

#### in_Q_WT_Brn1_rep3.atria.bam ####

            in_Q_WT_Brn1_rep3
            ------------------------------------------------------------
                      "${tallies[${unmapped}]}"  265277
                        "${tallies[${mapped}]}"  10184624

               "${tallies[${total_with_Mito}]}"  10449901
            "${tallies[${total_SC_with_Mito}]}"  10143319
            "${tallies[${total_SP_with_Mito}]}"  41305

               "${tallies[${total_sans_Mito}]}"  9490476
            "${tallies[${total_SC_sans_Mito}]}"  9449171
            "${tallies[${total_SP_sans_Mito}]}"  41305

                 "${tallies[${total_SC_Mito}]}"  694148
                 "${tallies[${total_SP_Mito}]}"  0

                       "${tallies[${total_1}]}"  9477888
                       "${tallies[${total_2}]}"  7918342
                       "${tallies[${total_3}]}"  7916426
                      "${tallies[${total_23}]}"  7901664
                      "${tallies[${total_30}]}"  7868730
                      "${tallies[${total_40}]}"  7755327
                      "${tallies[${total_42}]}"  7581526

                          "${tallies[${SC_0}]}"  9449171
                          "${tallies[${SC_1}]}"  9436650
                          "${tallies[${SC_2}]}"  7918277
                          "${tallies[${SC_3}]}"  7916378
                         "${tallies[${SC_23}]}"  7901651
                         "${tallies[${SC_30}]}"  7868721
                         "${tallies[${SC_40}]}"  7755325
                         "${tallies[${SC_42}]}"  7581526

                          "${tallies[${SP_0}]}"  41305
                          "${tallies[${SP_1}]}"  41238
                          "${tallies[${SP_2}]}"  65
                          "${tallies[${SP_3}]}"  48
                         "${tallies[${SP_23}]}"  13
                         "${tallies[${SP_30}]}"  9
                         "${tallies[${SP_40}]}"  2
                         "${tallies[${SP_42}]}"  0

                    "${tallies[${CC_SP2SC_0}]}"  0.004371
                    "${tallies[${CC_SP2SC_1}]}"  0.004370
                    "${tallies[${CC_SP2SC_2}]}"  0.000008
                    "${tallies[${CC_SP2SC_3}]}"  0.000006
                   "${tallies[${CC_SP2SC_23}]}"  0.000002
                   "${tallies[${CC_SP2SC_30}]}"  0.000001
                   "${tallies[${CC_SP2SC_40}]}"  0.000000
                   "${tallies[${CC_SP2SC_42}]}"  0.000000

#### Ch_Q_WT_Brn1_rep3.atria.bam ####

            Ch_Q_WT_Brn1_rep3
            ------------------------------------------------------------
                      "${tallies[${unmapped}]}"  612389
                        "${tallies[${mapped}]}"  8254445

               "${tallies[${total_with_Mito}]}"  8866834
            "${tallies[${total_SC_with_Mito}]}"  8178828
            "${tallies[${total_SP_with_Mito}]}"  75617

               "${tallies[${total_sans_Mito}]}"  8110801
            "${tallies[${total_SC_sans_Mito}]}"  8035185
            "${tallies[${total_SP_sans_Mito}]}"  75616

                 "${tallies[${total_SC_Mito}]}"  143643
                 "${tallies[${total_SP_Mito}]}"  1

                       "${tallies[${total_1}]}"  8088510
                       "${tallies[${total_2}]}"  5404473
                       "${tallies[${total_3}]}"  5400771
                      "${tallies[${total_23}]}"  5384058
                      "${tallies[${total_30}]}"  5354027
                      "${tallies[${total_40}]}"  5256137
                      "${tallies[${total_42}]}"  5132737

                          "${tallies[${SC_0}]}"  8035185
                          "${tallies[${SC_1}]}"  8013026
                          "${tallies[${SC_2}]}"  5404377
                          "${tallies[${SC_3}]}"  5400693
                         "${tallies[${SC_23}]}"  5384039
                         "${tallies[${SC_30}]}"  5354012
                         "${tallies[${SC_40}]}"  5256132
                         "${tallies[${SC_42}]}"  5132735

                          "${tallies[${SP_0}]}"  75616
                          "${tallies[${SP_1}]}"  75484
                          "${tallies[${SP_2}]}"  96
                          "${tallies[${SP_3}]}"  78
                         "${tallies[${SP_23}]}"  19
                         "${tallies[${SP_30}]}"  15
                         "${tallies[${SP_40}]}"  5
                         "${tallies[${SP_42}]}"  2

                    "${tallies[${CC_SP2SC_0}]}"  0.009411
                    "${tallies[${CC_SP2SC_1}]}"  0.009420
                    "${tallies[${CC_SP2SC_2}]}"  0.000018
                    "${tallies[${CC_SP2SC_3}]}"  0.000014
                   "${tallies[${CC_SP2SC_23}]}"  0.000004
                   "${tallies[${CC_SP2SC_30}]}"  0.000003
                   "${tallies[${CC_SP2SC_40}]}"  0.000001
                   "${tallies[${CC_SP2SC_42}]}"  0.000000

#### in_Q_WT_Brn1_repM.atria.bam ####

            in_Q_WT_Brn1_repM
            ------------------------------------------------------------
                      "${tallies[${unmapped}]}"  703297
                        "${tallies[${mapped}]}"  31997421

               "${tallies[${total_with_Mito}]}"  32700718
            "${tallies[${total_SC_with_Mito}]}"  31490835
            "${tallies[${total_SP_with_Mito}]}"  506586

               "${tallies[${total_sans_Mito}]}"  29460746
            "${tallies[${total_SC_sans_Mito}]}"  28985770
            "${tallies[${total_SP_sans_Mito}]}"  474976

                 "${tallies[${total_SC_Mito}]}"  2505065
                 "${tallies[${total_SP_Mito}]}"  31610

                       "${tallies[${total_1}]}"  29424922
                       "${tallies[${total_2}]}"  22783840
                       "${tallies[${total_3}]}"  22779363
                      "${tallies[${total_23}]}"  22737818
                      "${tallies[${total_30}]}"  22645393
                      "${tallies[${total_40}]}"  22316141
                      "${tallies[${total_42}]}"  21875592

                          "${tallies[${SC_0}]}"  28985770
                          "${tallies[${SC_1}]}"  28950275
                          "${tallies[${SC_2}]}"  22548588
                          "${tallies[${SC_3}]}"  22544159
                         "${tallies[${SC_23}]}"  22502886
                         "${tallies[${SC_30}]}"  22410603
                         "${tallies[${SC_40}]}"  22084100
                         "${tallies[${SC_42}]}"  21644790

                          "${tallies[${SP_0}]}"  474976
                          "${tallies[${SP_1}]}"  474647
                          "${tallies[${SP_2}]}"  235252
                          "${tallies[${SP_3}]}"  235204
                         "${tallies[${SP_23}]}"  234932
                         "${tallies[${SP_30}]}"  234790
                         "${tallies[${SP_40}]}"  232041
                         "${tallies[${SP_42}]}"  230802

                    "${tallies[${CC_SP2SC_0}]}"  0.016387
                    "${tallies[${CC_SP2SC_1}]}"  0.016395
                    "${tallies[${CC_SP2SC_2}]}"  0.010433
                    "${tallies[${CC_SP2SC_3}]}"  0.010433
                   "${tallies[${CC_SP2SC_23}]}"  0.010440
                   "${tallies[${CC_SP2SC_30}]}"  0.010477
                   "${tallies[${CC_SP2SC_40}]}"  0.010507
                   "${tallies[${CC_SP2SC_42}]}"  0.010663

#### Ch_Q_WT_Brn1_repM.atria.bam ####

            Ch_Q_WT_Brn1_repM
            ------------------------------------------------------------
                      "${tallies[${unmapped}]}"  3418392
                        "${tallies[${mapped}]}"  19168305

               "${tallies[${total_with_Mito}]}"  22586697
            "${tallies[${total_SC_with_Mito}]}"  18636191
            "${tallies[${total_SP_with_Mito}]}"  532114

               "${tallies[${total_sans_Mito}]}"  18830521
            "${tallies[${total_SC_sans_Mito}]}"  18311441
            "${tallies[${total_SP_sans_Mito}]}"  519080

                 "${tallies[${total_SC_Mito}]}"  324750
                 "${tallies[${total_SP_Mito}]}"  13034

                       "${tallies[${total_1}]}"  18790189
                       "${tallies[${total_2}]}"  12591991
                       "${tallies[${total_3}]}"  12585830
                      "${tallies[${total_23}]}"  12551898
                      "${tallies[${total_30}]}"  12487541
                      "${tallies[${total_40}]}"  12275084
                      "${tallies[${total_42}]}"  11998618

                          "${tallies[${SC_0}]}"  18311441
                          "${tallies[${SC_1}]}"  18271807
                          "${tallies[${SC_2}]}"  12323823
                          "${tallies[${SC_3}]}"  12317741
                         "${tallies[${SC_23}]}"  12284157
                         "${tallies[${SC_30}]}"  12219965
                         "${tallies[${SC_40}]}"  12011998
                         "${tallies[${SC_42}]}"  11737260

                          "${tallies[${SP_0}]}"  519080
                          "${tallies[${SP_1}]}"  518382
                          "${tallies[${SP_2}]}"  268168
                          "${tallies[${SP_3}]}"  268089
                         "${tallies[${SP_23}]}"  267741
                         "${tallies[${SP_30}]}"  267576
                         "${tallies[${SP_40}]}"  263086
                         "${tallies[${SP_42}]}"  261358

                    "${tallies[${CC_SP2SC_0}]}"  0.028347
                    "${tallies[${CC_SP2SC_1}]}"  0.028371
                    "${tallies[${CC_SP2SC_2}]}"  0.021760
                    "${tallies[${CC_SP2SC_3}]}"  0.021764
                   "${tallies[${CC_SP2SC_23}]}"  0.021796
                   "${tallies[${CC_SP2SC_30}]}"  0.021897
                   "${tallies[${CC_SP2SC_40}]}"  0.021902
                   "${tallies[${CC_SP2SC_42}]}"  0.022267


❯ cat "${a_tallies}"
sample  unmapped        mapped  total_with_Mito total_SC_with_Mito      total_SP_with_Mito      total_sans_Mito total_SC_sans_Mito      total_SP_sans_Mito      total_SC_Mito   total_SP_Mito     total_1 total_2 total_3 total_23        total_30        total_40        total_42        SC_0    SC_1    SC_2    SC_3    SC_23   SC_30   SC_40   SC_42   SP_0    SP_1      SP_2    SP_3    SP_23   SP_30   SP_40   SP_42   CC_SP2SC_0      CC_SP2SC_1      CC_SP2SC_2      CC_SP2SC_3      CC_SP2SC_23     CC_SP2SC_30     CC_SP2SC_40     CC_SP2SC_42
in_log_WT_Brn1_rep1     333904  13549404        13883308        13258095        291309  12972893        12706090        266803  552005  24506   12959145        10473556        10472049  10454533        10413321        10265586        10074757        12706090        12692504        10300915        10299428        10282048        10240927        10095174 9905209  266803  266641  172641  172621  172485  172394  170412  169548  0.020998        0.021008        0.016760        0.016760        0.016775        0.016834        0.016881 0.017117
Ch_log_WT_Brn1_rep1     499911  5752746 6252657 5609097 143649  5740916 5597556 143360  11541   289     5722062 3098907 3096507 3085337 3065407 2998196 2934584 5597556 5578873 3059048   3056672 3045575 3025669 2959163 2895786 143360  143189  39859   39835   39762   39738   39033   38798   0.025611        0.025666        0.013030        0.013032        0.013056  0.013134        0.013191        0.013398
...


❯ calculate_run_time "${start}" "${end}"  \
>     "Generate tab-separated table of alignment tallies/calculations."

Generate tab-separated table of alignment tallies/calculations.
Run time: 0h:38m:39s
```
</details>
<br />

<a id="calculate-ccss-styled-scaling-factors"></a>
#### Calculate CC/SS-styled scaling factors
<br />
<br />

<a id="on-calculating-scaling-factors"></a>
### On calculating scaling factors
<a id="email-from-christine"></a>
#### Email from Christine
<details>
<summary><i>Email from Christine</i></summary>

Title: ChIP-seq analyses: Friendly reminder to check if, when using deepTools, you input the scaling factors as reciprocals  
From: Cucinotta, Christine E  
To: Alavattam, Kris  
Time: Tue 7/18/2023 1:38 PM

Hi Kris,

I take the Pombe:Sc ratio of all the samples. Then take the ratio of the sample we want to normalize to (e.g., WT or a time zero) and divide by the test sample.

I have trouble writing it out so here’s an example:  
Q Pombe:Sc = 0.35 (control)  
5m Pombe:Sc = 0.26 (test)

Scaling factor for Q = 0.35/0.35 = 1 (use this for the scaling factor in deeptools) (control/control)  
Scaling factor for 5m = 0.35/0.26 = 1.35 (use this number for the scaling factor in deeptools) (control/test)

I hope this makes sense. Please let me know if not – happy to chat!

Thanks!

-Chrsitine
</details>
<br />

<a id="notes-on-using-excel-to-calculate-scaling-factors"></a>
#### Notes on using Excel to calculate scaling factors
<details>
<summary><i>Notes on using Excel to calculate scaling factors</i></summary>

Calculate the scaling factor by hand using the table of tallies/calculations
1. Calculate the quotient of a given sample's input SP-to-SC ratio and its ChIP SP-to-SC ratio&mdash;<i>this is the <b>scaling factor</i></b>
2. For example, for sample `Brn1_Q_rep1`, divide the `input` `SP-to-SC` value by the `ChIP` `SP-to-SC` value: 0.019802 ÷ 0.045415 = <b><i>0.436023</i></b>
3. Do this for all samples (e.g., sample-and-replicate-wise pairs of input and ChIP `bam`s)
</details>
<br />
