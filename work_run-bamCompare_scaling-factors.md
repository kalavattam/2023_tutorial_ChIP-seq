
`#work_run-bamCompare_scaling-factors-CC-style.md`

## Run `bamCompare` using scaling factors calculated
*...as performed by Christine Cuccinota*
### Code
<details>
<summary><i>Code: Run `bamCompare` using scaling factors calculated</i></summary>

```bash
#!/bin/bash

#  Get situated, set up variables and arrays ----------------------------------
tmux new -s bws
grabnode  # 16 cores, defaults

ml deepTools/3.5.1-foss-2021b
# which bamCoverage

d_work="${HOME}/tsukiyamalab/Kris/2023_rDNA/results"  # ., "${d_work}"
[[ ! -d "${d_work}/2023-0406/bws" ]] &&
	{
    	mkdir -p "${d_work}/2023-0406/bws/err_out"
	}

cd "${d_work}" || echo "cd'ing failed; check on this"


#  Set variables and create a "dictionary" for ChIP samples -------------------
#  ...set up to associate ChIP samples with input samples and size factors
p_data="${HOME}/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams"  # echo "${p_data}"
p_out="${d_work}/2023-0406/bws"  # echo "${p_out}"
p_eo="${d_work}/2023-0406/bws/err_out"  # echo "${p_eo}"
threads="${SLURM_CPUS_ON_NODE}"  # echo "${threads}"

unset dictionary
typeset -A dictionary
dictionary["${p_data}/Brn1_Q_rep1_ChIP"]="${p_data}/Brn1_Q_rep1_input:0.43602334"
dictionary["${p_data}/Brn1_Q_rep2_ChIP"]="${p_data}/Brn1_Q_rep2_input:0.558518908"
dictionary["${p_data}/Brn1_Q_rep3_ChIP"]="${p_data}/Brn1_Q_rep3_input:0.433102714"
dictionary["${p_data}/Brn1_Q_all_ChIP"]="${p_data}/Brn1_Q_all_input:0.568815116"
dictionary["${p_data}/Brn1_log_rep1_ChIP"]="${p_data}/Brn1_log_rep1_input:0.850833236"
dictionary["${p_data}/Brn1_log_rep2_ChIP"]="${p_data}/Brn1_log_rep2_input:0.851440329"
dictionary["${p_data}/Brn1_log_rep3_ChIP"]="${p_data}/Brn1_log_rep3_input:0.451079005"
dictionary["${p_data}/Brn1_log_all_ChIP"]="${p_data}/Brn1_log_all_input:0.764050676"

#  Check the dictionary of key-value pairs
for i in "${!dictionary[@]}"; do
	echo "  key (ChIP) .............. ${i}"
	echo "value (input) ............. ${dictionary["${i}"]%:*}"
	echo "value (scaling factor) .... ${dictionary["${i}"]#*:}"
	echo ""
done

#  Run echo tests for bamCompare with ChIP and input bams, and scaling factors
for i in "${!dictionary[@]}"; do
	echo "#### $(basename ${i}) ####"
	echo "  key (ChIP) .............. ${i}"
	echo "value (input) ............. ${dictionary["${i}"]%:*}"
	echo "value (scaling factor) .... ${dictionary["${i}"]#*:}"

	ChIP="${i}.bam"  # echo "${ChIP}"
	input="${dictionary["${i}"]%:*}.bam"  # echo "${input}"
	SF="${dictionary["${i}"]#*:}"  # echo "${SF}"
	outbw="${p_out}/$(basename "${ChIP}" .bam).bamCompare-to-input.SF-coef-${SF}.bw"  # echo "${outbw}"
	stderr="${p_eo}/$(basename "${ChIP}" .bam).bamCompare-to-input.SF-coef-${SF}.stderr.txt"  # echo "${stderr}"
	stdout="${p_eo}/$(basename "${ChIP}" .bam).bamCompare-to-input.SF-coef-${SF}.stdout.txt"  # echo "${stdout}"

	echo """
	bamCompare \\
		-p "${threads}" \\
		-b1 "${ChIP}" \\
		-b2 "${input}" \\
		-o "${outbw}" \\
		--scaleFactorsMethod None \\
		--scaleFactors ${SF}:1 \\
		--normalizeUsing None \\
		--minMappingQuality 0 \\
		--binSize 1 \\
		--smoothLength 4 \\
		--extendReads 150 \\
		--centerReads \\
			1> >(tee -a ${stderr}) \\
			2> >(tee -a ${stdout})
	
	"""
done

#  Run the actual calls to bamCompare
for i in "${!dictionary[@]}"; do
	echo "#### $(basename ${i}) ####"
	echo "  key (ChIP) .............. ${i}"
	echo "value (input) ............. ${dictionary["${i}"]%:*}"
	echo "value (scaling factor) .... ${dictionary["${i}"]#*:}"

	ChIP="${i}.bam"  # echo "${ChIP}"
	input="${dictionary["${i}"]%:*}.bam"  # echo "${input}"
	SF="${dictionary["${i}"]#*:}"  # echo "${SF}"
	outbw="${p_out}/$(basename "${ChIP}" .bam).bamCompare-to-input.SF-coef-${SF}.bw"  # echo "${outbw}"
	stderr="${p_eo}/$(basename "${ChIP}" .bam).bamCompare-to-input.SF-coef-${SF}.stderr.txt"  # echo "${stderr}"
	stdout="${p_eo}/$(basename "${ChIP}" .bam).bamCompare-to-input.SF-coef-${SF}.stdout.txt"  # echo "${stdout}"

	echo """
	bamCompare \\
		-p "${threads}" \\
		-b1 "${ChIP}" \\
		-b2 "${input}" \\
		-o "${outbw}" \\
		--scaleFactorsMethod None \\
		--scaleFactors ${SF}:1 \\
		--normalizeUsing None \\
		--minMappingQuality 0 \\
		--binSize 1 \\
		--smoothLength 4 \\
		--extendReads 150 \\
		--centerReads \\
			1> >(tee -a ${stderr}) \\
			2> >(tee -a ${stdout})
	
	"""

	bamCompare \
		-p "${threads}" \
		-b1 "${ChIP}" \
		-b2 "${input}" \
		-o "${outbw}" \
		--scaleFactorsMethod None \
		--scaleFactors ${SF}:1 \
		--normalizeUsing None \
		--minMappingQuality 0 \
		--binSize 1 \
		--smoothLength 4 \
		--extendReads 150 \
		--centerReads \
			1> >(tee -a ${stderr}) \
			2> >(tee -a ${stdout})
done

# bamCoverage \
# 	-b ${SAMPLE}/${SAMPLE}.nodup.norm.bam \
# 	-o /scratch/${PBS_JOBID}/${SAMPLE}/${SAMPLE}.bw \
# 	-p 1 \
# 	--normalizeTo1x ${genomesize} \
# 	--scaleFactor ${SPIKE} \
# 	--ignoreForNormalization chrX chrY \
# 	-bl ${ENCODEBED} \
# 	--smoothLength 40 \
# 	--binSize 1 \
# 	--centerReads \
# 	--extendReads 150
```
</details>
<br />

### Printed
<details>
<summary><i>Printed: Run `bamCompare` using scaling factors calculated</i></summary>

```txt
❯ d_work="${HOME}/tsukiyamalab/Kris/2023_rDNA/results"  # cd "${d_work}"


❯ ., "${d_work}"
total 208K
drwxrws--- 4 kalavatt  110 May 30 14:03 ./
drwxrws--- 7 kalavatt  174 Apr 27 13:59 ../
drwxrws--- 2 kalavatt  412 Mar  3 00:01 2023-0228/
drwxrws--- 3 kalavatt   78 May 30 14:03 2023-0406/
-rw-rw---- 1 kalavatt 4.0K May 30 14:03 ._.DS_Store
-rw-rw---- 1 kalavatt 6.1K May 30 14:03 .DS_Store


❯ if [[ ! -d "${d_work}/2023-0406/bws" ]]; then
>    mkdir -p "${d_work}/2023-0406/bws/err_out"
>fi
mkdir: created directory '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws'
mkdir: created directory '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/err_out'


❯ cd "${d_work}" && pwd
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results


❯ ml deepTools/3.5.1-foss-2021b


❯ which bamCoverage
/app/software/deepTools/3.5.1-foss-2021b/bin/bamCoverage


❯ grabnode
--------------------------------------------------------------
     ** Please review options as the order has changed **

 - memory is requested first
 - CPU selections will be adjusted based on memory request
--------------------------------------------------------------

How much memory (GB) will you need?
Enter a value from 1 to 683 [20]:
How many cores (CPUs) would you like to grab on the node?
enter a value between 1 and 36 [36]: 16
Please enter the max number of days you would like to grab this node: [1-7] 1
Do you need a GPU ? [y/N]N

You have requested 16 CPUs on this node/server for 1 days or until you type exit.

Warning: If you exit this shell before your jobs are finished, your jobs
on this node/server will be terminated. Please use sbatch for larger jobs.

Shared PI folders can be found in: /fh/fast, /fh/scratch and /fh/secure.

Requesting Queue: campus-new cores: 16 memory: 20 gpu: NONE
srun: job 21838445 queued and waiting for resources
srun: job 21838445 has been allocated resources


❯ ml deepTools/3.5.1-foss-2021b


❯ which bamCoverage
/app/software/deepTools/3.5.1-foss-2021b/bin/bamCoverage


❯ d_work="${HOME}/tsukiyamalab/Kris/2023_rDNA/results"  # ., "${d_work}"


❯ if [[ ! -d "${d_work}/2023-0406/bws" ]]; then
>    mkdir -p "${d_work}/2023-0406/bws/err_out"
>fi


❯ cd "${d_work}" && pwd
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results


❯ p_data="${HOME}/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams"


❯ p_out="${d_work}/2023-0406/bws"


❯ p_eo="${d_work}/2023-0406/bws/err_out"


❯ unset dictionary


❯ typeset -A dictionary


❯ dictionary["${p_data}/Brn1_Q_rep1_ChIP"]="${p_data}/Brn1_Q_rep1_input:0.43602334"


❯ dictionary["${p_data}/Brn1_Q_rep2_ChIP"]="${p_data}/Brn1_Q_rep2_input:0.558518908"


❯ dictionary["${p_data}/Brn1_Q_rep3_ChIP"]="${p_data}/Brn1_Q_rep3_input:0.433102714"


❯ dictionary["${p_data}/Brn1_Q_all_ChIP"]="${p_data}/Brn1_Q_all_input:0.568815116"


❯ dictionary["${p_data}/Brn1_log_rep1_ChIP"]="${p_data}/Brn1_log_rep1_input:0.850833236"


❯ dictionary["${p_data}/Brn1_log_rep2_ChIP"]="${p_data}/Brn1_log_rep2_input:0.851440329"


❯ dictionary["${p_data}/Brn1_log_rep3_ChIP"]="${p_data}/Brn1_log_rep3_input:0.451079005"


❯ dictionary["${p_data}/Brn1_log_all_ChIP"]="${p_data}/Brn1_log_all_input:0.764050676"


❯ for i in "${!dictionary[@]}"; do
>    echo "  key (ChIP) .............. ${i}"
>    echo "value (input) ............. ${dictionary["${i}"]%:*}"
>    echo "value (scaling factor) .... ${dictionary["${i}"]#*:}"
>    echo ""
>done
  key (ChIP) .............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep2_ChIP
value (input) ............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep2_input
value (scaling factor) .... 0.851440329

  key (ChIP) .............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep2_ChIP
value (input) ............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep2_input
value (scaling factor) .... 0.558518908

  key (ChIP) .............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep3_ChIP
value (input) ............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep3_input
value (scaling factor) .... 0.433102714

  key (ChIP) .............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_all_ChIP
value (input) ............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_all_input
value (scaling factor) .... 0.764050676

  key (ChIP) .............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep1_ChIP
value (input) ............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep1_input
value (scaling factor) .... 0.850833236

  key (ChIP) .............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep1_ChIP
value (input) ............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep1_input
value (scaling factor) .... 0.43602334

  key (ChIP) .............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep3_ChIP
value (input) ............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep3_input
value (scaling factor) .... 0.451079005

  key (ChIP) .............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_all_ChIP
value (input) ............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_all_input
value (scaling factor) .... 0.568815116


❯ for i in "${!dictionary[@]}"; do
>    echo "#### $(basename ${i}) ####"
>    echo "  key (ChIP) .............. ${i}"
>    echo "value (input) ............. ${dictionary["${i}"]%:*}"
>    echo "value (scaling factor) .... ${dictionary["${i}"]#*:}"
>
>    ChIP="${i}.bam"  # echo "${ChIP}"
>    input="${dictionary["${i}"]%:*}.bam"  # echo "${input}"
>    SF="${dictionary["${i}"]#*:}"  # echo "${SF}"
>    outbw="${p_out}/$(basename "${ChIP}" .bam).bamCompare-to-input.SF-coef-${SF}.bw"  # echo "${outbw}"
>    stderr="${p_eo}/$(basename "${ChIP}" .bam).bamCompare-to-input.SF-coef-${SF}.stderr.txt"  # echo "${stderr}"
>    stdout="${p_eo}/$(basename "${ChIP}" .bam).bamCompare-to-input.SF-coef-${SF}.stdout.txt"  # echo "${stdout}"
>
>    echo """
>    bamCompare \\
>        -p "${threads}" \\
>        -b1 "${ChIP}" \\
>        -b2 "${input}" \\
>        -o "${outbw}" \\
>        --scaleFactorsMethod None \\
>        --scaleFactors ${SF}:1 \\
>        --normalizeUsing None \\
>        --minMappingQuality 0 \\
>        --binSize 1 \\
>        --smoothLength 4 \\
>        --extendReads 150 \\
>        --centerReads \\
>            1> >(tee -a ${stderr}) \\
>            2> >(tee -a ${stdout})
>
>    """
>done
#### Brn1_log_rep2_ChIP ####
  key (ChIP) .............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep2_ChIP
value (input) ............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep2_input
value (scaling factor) .... 0.851440329

    bamCompare \
        -p 16 \
        -b1 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep2_ChIP.bam \
        -b2 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep2_input.bam \
        -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/Brn1_log_rep2_ChIP.bamCompare-to-input.SF-coef-0.851440329.bw \
        --scaleFactorsMethod None \
        --scaleFactors 0.851440329:1 \
        --normalizeUsing None \
        --minMappingQuality 0 \
        --binSize 1 \
        --smoothLength 4 \
        --extendReads 150 \
        --centerReads \
            1> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/err_out/Brn1_log_rep2_ChIP.bamCompare-to-input.SF-coef-0.851440329.stderr.txt) \
            2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/err_out/Brn1_log_rep2_ChIP.bamCompare-to-input.SF-coef-0.851440329.stdout.txt)


#### Brn1_Q_rep2_ChIP ####
  key (ChIP) .............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep2_ChIP
value (input) ............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep2_input
value (scaling factor) .... 0.558518908

    bamCompare \
        -p 16 \
        -b1 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep2_ChIP.bam \
        -b2 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep2_input.bam \
        -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/Brn1_Q_rep2_ChIP.bamCompare-to-input.SF-coef-0.558518908.bw \
        --scaleFactorsMethod None \
        --scaleFactors 0.558518908:1 \
        --normalizeUsing None \
        --minMappingQuality 0 \
        --binSize 1 \
        --smoothLength 4 \
        --extendReads 150 \
        --centerReads \
            1> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/err_out/Brn1_Q_rep2_ChIP.bamCompare-to-input.SF-coef-0.558518908.stderr.txt) \
            2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/err_out/Brn1_Q_rep2_ChIP.bamCompare-to-input.SF-coef-0.558518908.stdout.txt)


#### Brn1_Q_rep3_ChIP ####
  key (ChIP) .............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep3_ChIP
value (input) ............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep3_input
value (scaling factor) .... 0.433102714

    bamCompare \
        -p 16 \
        -b1 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep3_ChIP.bam \
        -b2 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep3_input.bam \
        -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/Brn1_Q_rep3_ChIP.bamCompare-to-input.SF-coef-0.433102714.bw \
        --scaleFactorsMethod None \
        --scaleFactors 0.433102714:1 \
        --normalizeUsing None \
        --minMappingQuality 0 \
        --binSize 1 \
        --smoothLength 4 \
        --extendReads 150 \
        --centerReads \
            1> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/err_out/Brn1_Q_rep3_ChIP.bamCompare-to-input.SF-coef-0.433102714.stderr.txt) \
            2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/err_out/Brn1_Q_rep3_ChIP.bamCompare-to-input.SF-coef-0.433102714.stdout.txt)


#### Brn1_log_all_ChIP ####
  key (ChIP) .............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_all_ChIP
value (input) ............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_all_input
value (scaling factor) .... 0.764050676

    bamCompare \
        -p 16 \
        -b1 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_all_ChIP.bam \
        -b2 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_all_input.bam \
        -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/Brn1_log_all_ChIP.bamCompare-to-input.SF-coef-0.764050676.bw \
        --scaleFactorsMethod None \
        --scaleFactors 0.764050676:1 \
        --normalizeUsing None \
        --minMappingQuality 0 \
        --binSize 1 \
        --smoothLength 4 \
        --extendReads 150 \
        --centerReads \
            1> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/err_out/Brn1_log_all_ChIP.bamCompare-to-input.SF-coef-0.764050676.stderr.txt) \
            2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/err_out/Brn1_log_all_ChIP.bamCompare-to-input.SF-coef-0.764050676.stdout.txt)


#### Brn1_log_rep1_ChIP ####
  key (ChIP) .............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep1_ChIP
value (input) ............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep1_input
value (scaling factor) .... 0.850833236

    bamCompare \
        -p 16 \
        -b1 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep1_ChIP.bam \
        -b2 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep1_input.bam \
        -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/Brn1_log_rep1_ChIP.bamCompare-to-input.SF-coef-0.850833236.bw \
        --scaleFactorsMethod None \
        --scaleFactors 0.850833236:1 \
        --normalizeUsing None \
        --minMappingQuality 0 \
        --binSize 1 \
        --smoothLength 4 \
        --extendReads 150 \
        --centerReads \
            1> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/err_out/Brn1_log_rep1_ChIP.bamCompare-to-input.SF-coef-0.850833236.stderr.txt) \
            2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/err_out/Brn1_log_rep1_ChIP.bamCompare-to-input.SF-coef-0.850833236.stdout.txt)


#### Brn1_Q_rep1_ChIP ####
  key (ChIP) .............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep1_ChIP
value (input) ............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep1_input
value (scaling factor) .... 0.43602334

    bamCompare \
        -p 16 \
        -b1 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep1_ChIP.bam \
        -b2 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep1_input.bam \
        -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/Brn1_Q_rep1_ChIP.bamCompare-to-input.SF-coef-0.43602334.bw \
        --scaleFactorsMethod None \
        --scaleFactors 0.43602334:1 \
        --normalizeUsing None \
        --minMappingQuality 0 \
        --binSize 1 \
        --smoothLength 4 \
        --extendReads 150 \
        --centerReads \
            1> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/err_out/Brn1_Q_rep1_ChIP.bamCompare-to-input.SF-coef-0.43602334.stderr.txt) \
            2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/err_out/Brn1_Q_rep1_ChIP.bamCompare-to-input.SF-coef-0.43602334.stdout.txt)


#### Brn1_log_rep3_ChIP ####
  key (ChIP) .............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep3_ChIP
value (input) ............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep3_input
value (scaling factor) .... 0.451079005

    bamCompare \
        -p 16 \
        -b1 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep3_ChIP.bam \
        -b2 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep3_input.bam \
        -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/Brn1_log_rep3_ChIP.bamCompare-to-input.SF-coef-0.451079005.bw \
        --scaleFactorsMethod None \
        --scaleFactors 0.451079005:1 \
        --normalizeUsing None \
        --minMappingQuality 0 \
        --binSize 1 \
        --smoothLength 4 \
        --extendReads 150 \
        --centerReads \
            1> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/err_out/Brn1_log_rep3_ChIP.bamCompare-to-input.SF-coef-0.451079005.stderr.txt) \
            2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/err_out/Brn1_log_rep3_ChIP.bamCompare-to-input.SF-coef-0.451079005.stdout.txt)


#### Brn1_Q_all_ChIP ####
  key (ChIP) .............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_all_ChIP
value (input) ............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_all_input
value (scaling factor) .... 0.568815116

    bamCompare \
        -p 16 \
        -b1 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_all_ChIP.bam \
        -b2 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_all_input.bam \
        -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/Brn1_Q_all_ChIP.bamCompare-to-input.SF-coef-0.568815116.bw \
        --scaleFactorsMethod None \
        --scaleFactors 0.568815116:1 \
        --normalizeUsing None \
        --minMappingQuality 0 \
        --binSize 1 \
        --smoothLength 4 \
        --extendReads 150 \
        --centerReads \
            1> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/err_out/Brn1_Q_all_ChIP.bamCompare-to-input.SF-coef-0.568815116.stderr.txt) \
            2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/err_out/Brn1_Q_all_ChIP.bamCompare-to-input.SF-coef-0.568815116.stdout.txt)


❯ for i in "${!dictionary[@]}"; do
>    echo "#### $(basename ${i}) ####"
>    echo "  key (ChIP) .............. ${i}"
>    echo "value (input) ............. ${dictionary["${i}"]%:*}"
>    echo "value (scaling factor) .... ${dictionary["${i}"]#*:}"
>
>    ChIP="${i}.bam"  # echo "${ChIP}"
>    input="${dictionary["${i}"]%:*}.bam"  # echo "${input}"
>    SF="${dictionary["${i}"]#*:}"  # echo "${SF}"
>    outbw="${p_out}/$(basename "${ChIP}" .bam).bamCompare-to-input.SF-coef-${SF}.bw"  # echo "${outbw}"
>    stderr="${p_eo}/$(basename "${ChIP}" .bam).bamCompare-to-input.SF-coef-${SF}.stderr.txt"  # echo "${stderr}"
>    stdout="${p_eo}/$(basename "${ChIP}" .bam).bamCompare-to-input.SF-coef-${SF}.stdout.txt"  # echo "${stdout}"
>
>    echo """
>    bamCompare \\
>        -p "${threads}" \\
>        -b1 "${ChIP}" \\
>        -b2 "${input}" \\
>        -o "${outbw}" \\
>        --scaleFactorsMethod None \\
>        --scaleFactors ${SF}:1 \\
>        --normalizeUsing None \\
>        --minMappingQuality 0 \\
>        --binSize 1 \\
>        --smoothLength 4 \\
>        --extendReads 150 \\
>        --centerReads \\
>            1> >(tee -a ${stderr}) \\
>            2> >(tee -a ${stdout})
>
>    """
>
>    bamCompare \
>        -p "${threads}" \
>        -b1 "${ChIP}" \
>        -b2 "${input}" \
>        -o "${outbw}" \
>        --scaleFactorsMethod None \
>        --scaleFactors ${SF}:1 \
>        --normalizeUsing None \
>        --minMappingQuality 0 \
>        --binSize 1 \
>        --smoothLength 4 \
>        --extendReads 150 \
>        --centerReads \
>            1> >(tee -a ${stderr}) \
>            2> >(tee -a ${stdout})
>done
#### Brn1_log_rep2_ChIP ####
  key (ChIP) .............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep2_ChIP
value (input) ............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep2_input
value (scaling factor) .... 0.851440329

    bamCompare \
        -p 16 \
        -b1 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep2_ChIP.bam \
        -b2 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep2_input.bam \
        -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/Brn1_log_rep2_ChIP.bamCompare-to-input.SF-coef-0.851440329.bw \
        --scaleFactorsMethod None \
        --scaleFactors 0.851440329:1 \
        --normalizeUsing None \
        --minMappingQuality 0 \
        --binSize 1 \
        --smoothLength 4 \
        --extendReads 150 \
        --centerReads \
            1> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/err_out/Brn1_log_rep2_ChIP.bamCompare-to-input.SF-coef-0.851440329.stderr.txt) \
            2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/err_out/Brn1_log_rep2_ChIP.bamCompare-to-input.SF-coef-0.851440329.stdout.txt)


bamFilesList: ['/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep2_ChIP.bam', '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep2_input.bam']
binLength: 1
numberOfSamples: 0
blackListFileName: None
skipZeroOverZero: False
bed_and_bin: False
genomeChunkSize: None
defaultFragmentLength: 150
numberOfProcessors: 16
verbose: False
region: None
bedFile: None
minMappingQuality: 0
ignoreDuplicates: False
chrsToSkip: []
stepSize: 1
center_read: True
samFlag_include: None
samFlag_exclude: None
minFragmentLength: 0
maxFragmentLength: 0
zerosToNans: False
smoothLength: 4
save_data: False
out_file_for_raw_data: None
maxPairedFragmentLength: 600


#### Brn1_Q_rep2_ChIP ####
  key (ChIP) .............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep2_ChIP
value (input) ............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep2_input
value (scaling factor) .... 0.558518908

    bamCompare \
        -p 16 \
        -b1 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep2_ChIP.bam \
        -b2 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep2_input.bam \
        -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/Brn1_Q_rep2_ChIP.bamCompare-to-input.SF-coef-0.558518908.bw \
        --scaleFactorsMethod None \
        --scaleFactors 0.558518908:1 \
        --normalizeUsing None \
        --minMappingQuality 0 \
        --binSize 1 \
        --smoothLength 4 \
        --extendReads 150 \
        --centerReads \
            1> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/err_out/Brn1_Q_rep2_ChIP.bamCompare-to-input.SF-coef-0.558518908.stderr.txt) \
            2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/err_out/Brn1_Q_rep2_ChIP.bamCompare-to-input.SF-coef-0.558518908.stdout.txt)


bamFilesList: ['/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep2_ChIP.bam', '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep2_input.bam']
binLength: 1
numberOfSamples: 0
blackListFileName: None
skipZeroOverZero: False
bed_and_bin: False
genomeChunkSize: None
defaultFragmentLength: 150
numberOfProcessors: 16
verbose: False
region: None
bedFile: None
minMappingQuality: 0
ignoreDuplicates: False
chrsToSkip: []
stepSize: 1
center_read: True
samFlag_include: None
samFlag_exclude: None
minFragmentLength: 0
maxFragmentLength: 0
zerosToNans: False
smoothLength: 4
save_data: False
out_file_for_raw_data: None
maxPairedFragmentLength: 600


#### Brn1_Q_rep3_ChIP ####
  key (ChIP) .............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep3_ChIP
value (input) ............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep3_input
value (scaling factor) .... 0.433102714

    bamCompare \
        -p 16 \
        -b1 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep3_ChIP.bam \
        -b2 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep3_input.bam \
        -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/Brn1_Q_rep3_ChIP.bamCompare-to-input.SF-coef-0.433102714.bw \
        --scaleFactorsMethod None \
        --scaleFactors 0.433102714:1 \
        --normalizeUsing None \
        --minMappingQuality 0 \
        --binSize 1 \
        --smoothLength 4 \
        --extendReads 150 \
        --centerReads \
            1> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/err_out/Brn1_Q_rep3_ChIP.bamCompare-to-input.SF-coef-0.433102714.stderr.txt) \
            2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/err_out/Brn1_Q_rep3_ChIP.bamCompare-to-input.SF-coef-0.433102714.stdout.txt)


bamFilesList: ['/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep3_ChIP.bam', '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep3_input.bam']
binLength: 1
numberOfSamples: 0
blackListFileName: None
skipZeroOverZero: False
bed_and_bin: False
genomeChunkSize: None
defaultFragmentLength: 150
numberOfProcessors: 16
verbose: False
region: None
bedFile: None
minMappingQuality: 0
ignoreDuplicates: False
chrsToSkip: []
stepSize: 1
center_read: True
samFlag_include: None
samFlag_exclude: None
minFragmentLength: 0
maxFragmentLength: 0
zerosToNans: False
smoothLength: 4
save_data: False
out_file_for_raw_data: None
maxPairedFragmentLength: 600


#### Brn1_log_all_ChIP ####
  key (ChIP) .............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_all_ChIP
value (input) ............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_all_input
value (scaling factor) .... 0.764050676

    bamCompare \
        -p 16 \
        -b1 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_all_ChIP.bam \
        -b2 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_all_input.bam \
        -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/Brn1_log_all_ChIP.bamCompare-to-input.SF-coef-0.764050676.bw \
        --scaleFactorsMethod None \
        --scaleFactors 0.764050676:1 \
        --normalizeUsing None \
        --minMappingQuality 0 \
        --binSize 1 \
        --smoothLength 4 \
        --extendReads 150 \
        --centerReads \
            1> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/err_out/Brn1_log_all_ChIP.bamCompare-to-input.SF-coef-0.764050676.stderr.txt) \
            2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/err_out/Brn1_log_all_ChIP.bamCompare-to-input.SF-coef-0.764050676.stdout.txt)


bamFilesList: ['/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_all_ChIP.bam', '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_all_input.bam']
binLength: 1
numberOfSamples: 0
blackListFileName: None
skipZeroOverZero: False
bed_and_bin: False
genomeChunkSize: None
defaultFragmentLength: 150
numberOfProcessors: 16
verbose: False
region: None
bedFile: None
minMappingQuality: 0
ignoreDuplicates: False
chrsToSkip: []
stepSize: 1
center_read: True
samFlag_include: None
samFlag_exclude: None
minFragmentLength: 0
maxFragmentLength: 0
zerosToNans: False
smoothLength: 4
save_data: False
out_file_for_raw_data: None
maxPairedFragmentLength: 600


#### Brn1_log_rep1_ChIP ####
  key (ChIP) .............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep1_ChIP
value (input) ............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep1_input
value (scaling factor) .... 0.850833236

    bamCompare \
        -p 16 \
        -b1 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep1_ChIP.bam \
        -b2 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep1_input.bam \
        -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/Brn1_log_rep1_ChIP.bamCompare-to-input.SF-coef-0.850833236.bw \
        --scaleFactorsMethod None \
        --scaleFactors 0.850833236:1 \
        --normalizeUsing None \
        --minMappingQuality 0 \
        --binSize 1 \
        --smoothLength 4 \
        --extendReads 150 \
        --centerReads \
            1> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/err_out/Brn1_log_rep1_ChIP.bamCompare-to-input.SF-coef-0.850833236.stderr.txt) \
            2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/err_out/Brn1_log_rep1_ChIP.bamCompare-to-input.SF-coef-0.850833236.stdout.txt)


bamFilesList: ['/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep1_ChIP.bam', '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep1_input.bam']
binLength: 1
numberOfSamples: 0
blackListFileName: None
skipZeroOverZero: False
bed_and_bin: False
genomeChunkSize: None
defaultFragmentLength: 150
numberOfProcessors: 16
verbose: False
region: None
bedFile: None
minMappingQuality: 0
ignoreDuplicates: False
chrsToSkip: []
stepSize: 1
center_read: True
samFlag_include: None
samFlag_exclude: None
minFragmentLength: 0
maxFragmentLength: 0
zerosToNans: False
smoothLength: 4
save_data: False
out_file_for_raw_data: None
maxPairedFragmentLength: 600


#### Brn1_Q_rep1_ChIP ####
  key (ChIP) .............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep1_ChIP
value (input) ............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep1_input
value (scaling factor) .... 0.43602334

    bamCompare \
        -p 16 \
        -b1 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep1_ChIP.bam \
        -b2 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep1_input.bam \
        -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/Brn1_Q_rep1_ChIP.bamCompare-to-input.SF-coef-0.43602334.bw \
        --scaleFactorsMethod None \
        --scaleFactors 0.43602334:1 \
        --normalizeUsing None \
        --minMappingQuality 0 \
        --binSize 1 \
        --smoothLength 4 \
        --extendReads 150 \
        --centerReads \
            1> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/err_out/Brn1_Q_rep1_ChIP.bamCompare-to-input.SF-coef-0.43602334.stderr.txt) \
            2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/err_out/Brn1_Q_rep1_ChIP.bamCompare-to-input.SF-coef-0.43602334.stdout.txt)


bamFilesList: ['/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep1_ChIP.bam', '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_rep1_input.bam']
binLength: 1
numberOfSamples: 0
blackListFileName: None
skipZeroOverZero: False
bed_and_bin: False
genomeChunkSize: None
defaultFragmentLength: 150
numberOfProcessors: 16
verbose: False
region: None
bedFile: None
minMappingQuality: 0
ignoreDuplicates: False
chrsToSkip: []
stepSize: 1
center_read: True
samFlag_include: None
samFlag_exclude: None
minFragmentLength: 0
maxFragmentLength: 0
zerosToNans: False
smoothLength: 4
save_data: False
out_file_for_raw_data: None
maxPairedFragmentLength: 600


#### Brn1_log_rep3_ChIP ####
  key (ChIP) .............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep3_ChIP
value (input) ............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep3_input
value (scaling factor) .... 0.451079005

    bamCompare \
        -p 16 \
        -b1 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep3_ChIP.bam \
        -b2 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep3_input.bam \
        -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/Brn1_log_rep3_ChIP.bamCompare-to-input.SF-coef-0.451079005.bw \
        --scaleFactorsMethod None \
        --scaleFactors 0.451079005:1 \
        --normalizeUsing None \
        --minMappingQuality 0 \
        --binSize 1 \
        --smoothLength 4 \
        --extendReads 150 \
        --centerReads \
            1> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/err_out/Brn1_log_rep3_ChIP.bamCompare-to-input.SF-coef-0.451079005.stderr.txt) \
            2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/err_out/Brn1_log_rep3_ChIP.bamCompare-to-input.SF-coef-0.451079005.stdout.txt)


bamFilesList: ['/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep3_ChIP.bam', '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_log_rep3_input.bam']
binLength: 1
numberOfSamples: 0
blackListFileName: None
skipZeroOverZero: False
bed_and_bin: False
genomeChunkSize: None
defaultFragmentLength: 150
numberOfProcessors: 16
verbose: False
region: None
bedFile: None
minMappingQuality: 0
ignoreDuplicates: False
chrsToSkip: []
stepSize: 1
center_read: True
samFlag_include: None
samFlag_exclude: None
minFragmentLength: 0
maxFragmentLength: 0
zerosToNans: False
smoothLength: 4
save_data: False
out_file_for_raw_data: None
maxPairedFragmentLength: 600


#### Brn1_Q_all_ChIP ####
  key (ChIP) .............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_all_ChIP
value (input) ............. /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_all_input
value (scaling factor) .... 0.568815116

    bamCompare \
        -p 16 \
        -b1 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_all_ChIP.bam \
        -b2 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_all_input.bam \
        -o /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/Brn1_Q_all_ChIP.bamCompare-to-input.SF-coef-0.568815116.bw \
        --scaleFactorsMethod None \
        --scaleFactors 0.568815116:1 \
        --normalizeUsing None \
        --minMappingQuality 0 \
        --binSize 1 \
        --smoothLength 4 \
        --extendReads 150 \
        --centerReads \
            1> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/err_out/Brn1_Q_all_ChIP.bamCompare-to-input.SF-coef-0.568815116.stderr.txt) \
            2> >(tee -a /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bws/err_out/Brn1_Q_all_ChIP.bamCompare-to-input.SF-coef-0.568815116.stdout.txt)


bamFilesList: ['/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_all_ChIP.bam', '/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406/bams/Brn1_Q_all_input.bam']
binLength: 1
numberOfSamples: 0
blackListFileName: None
skipZeroOverZero: False
bed_and_bin: False
genomeChunkSize: None
defaultFragmentLength: 150
numberOfProcessors: 16
verbose: False
region: None
bedFile: None
minMappingQuality: 0
ignoreDuplicates: False
chrsToSkip: []
stepSize: 1
center_read: True
samFlag_include: None
samFlag_exclude: None
minFragmentLength: 0
maxFragmentLength: 0
zerosToNans: False
smoothLength: 4
save_data: False
out_file_for_raw_data: None
maxPairedFragmentLength: 600
```
</details>
<br />
