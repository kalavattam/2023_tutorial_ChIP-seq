<details>
<summary><i>Printed: </i></summary>

```txt
❯ if [[
>        ! -f "${bam}" \
>     && ! -f "${bam_coor}" \
>     && ! -f "${bam_quer}"
> ]]; then
>     if [[ "${mode}" == "paired" ]]; then
>         bowtie2 \
>             -p "${threads}" \
>             -x "${index}" \
>             --very-sensitive-local \
>             --no-unal \
>             --no-mixed \
>             --no-discordant \
>             --no-overlap \
>             --no-dovetail \
>             --phred33 \
>             -I 10 \
>             -X 700 \
>             -1 "${fastq_1}" \
>             -2 "${fastq_2}" \
>                 | samtools view \
>                     -@ "${threads}" \
>                     -b \
>                     -f 2 \
>                     -q "${mapq}" \
>                     -o "${bam}"
>     elif [[ "${mode}" == "single" ]]; then
>         bowtie2 \
>             -p "${threads}" \
>             -x "${index}" \
>             --very-sensitive-local \
>             --no-unal \
>             --phred33 \
>             -U "${fastq_1}" \
>                 | samtools view \
>                     -@ "${threads}" \
>                     -b \
>                     -q "${mapq}" \
>                     -o "${bam}"
>     fi
> else
>     echo_note "BAM file exists. Skipping alignment operations (step #1)."
> fi
15419397 reads; of these:
  15419397 (100.00%) were paired; of these:
    901946 (5.85%) aligned concordantly 0 times
    9532537 (61.82%) aligned concordantly exactly 1 time
    4984914 (32.33%) aligned concordantly >1 times
94.15% overall alignment rate


❯ if [[ -f "${bam}" ]]; then
>     if [[ ! -f "${bam_quer}" ]]; then>         samtools sort \
>             -@ "${threads}" \
>             -n \
>             -o "${bam_quer}" \
>             "${bam}"
> 
>         if [[ "${mode}" == "paired" ]]; then
>             #  After sorting by queryname, fix the paired read-mate information, which
>             #+ is required for subsequent operations
>             if [[ -f "${bam_quer}" ]]; then
>                 samtools fixmate \
>                     -@ "${threads}" \>                     -m \
>                     "${bam_quer}" \
>                     "${bam_quer%.bam}.tmp.bam"
> 
>                 #  Replace the original queryname-sorted BAM with queryname-sorted
>                 #+ mate-fixed BAM
>                 if [[ -f "${bam_quer%.bam}.tmp.bam" ]]; then
>                     mv -f \
>                         "${bam_quer%.bam}.tmp.bam" \
>                         "${bam_quer}">                 fi
>             fi
>         fi
>     else
>         echo_note "Queryname-sored BAM file exists. Skipping sort operations (step #2)."
>     fi
> else
>     error_and_exit "BAM infile does not exist (step #2)."
> fi
[bam_sort_core] merging from 1 files and 8 in-memory blocks...
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_tutorial_ChIP-seq/03_bam/bowtie2/bam/IP_G1_Hho1_6336.sort-qname.tmp.bam' -> '/home/kalavatt/tsukiyamalab/Kris/2023_tutorial_ChIP-seq/03_bam/bowtie2/bam/IP_G1_Hho1_6336.sort-qname.bam'


❯ if [[ -f "${bam_quer}" ]]; then
>     if [[ ! -f "${bam_coor}" ]]; then>         samtools sort \
>             -@ "${threads}" \
>             -o "${bam_coor}" \
>             "${bam_quer}"
> 
>         #  Index the coordinate-sorted BAM file
>         if [[ ! -f "${bam_coor}.bai" ]]; then
>             samtools index \
>                 -@ "${threads}" \
>                 "${bam_coor}"
>         fi
>     else
>         echo_note "Coordinate-sored BAM file exists. Skipping sort operations (step #3)."
>     fi
> else
>     error_and_exit "Queryname-sorted BAM infile does not exist (step #3)."
> fi
[bam_sort_core] merging from 1 files and 8 in-memory blocks...


❯ if [[ -f "${bam_coor}" && -f "${bam_coor}.bai" ]]; then
>         #  Mark duplicate alignments in the coordinate-sorted BAM file
>         samtools markdup \
>             -@ "${threads}" \
>             "${bam_coor}" \
>             "${bam_coor%.bam}.tmp.bam"
> 
>         #  Replace the original coordinate-sorted BAM with one in which
>         #+ duplicates alignments are marked
>         if [[ -f "${bam_coor%.bam}.tmp.bam" ]]; then
>             mv -f \
>                 "${bam_coor%.bam}.tmp.bam" \
>                 "${bam_coor}"
>         fi
> fi
renamed '/home/kalavatt/tsukiyamalab/Kris/2023_tutorial_ChIP-seq/03_bam/bowtie2/bam/IP_G1_Hho1_6336.sort-coord.tmp.bam' -> '/home/kalavatt/tsukiyamalab/Kris/2023_tutorial_ChIP-seq/03_bam/bowtie2/bam/IP_G1_Hho1_6336.sort-coord.bam'


❯ if [[ -f "${bam_coor}" ]]; then
>     if [[ -f "${bam_coor}.bai" ]]; then
>         if \
>             samtools view -H "${bam_coor}" \
>                | grep -q "@PG.*CL:samtools markdup";
>         then
>             if [[ ! -f "${txt_met}" ]]; then
>                 #  Run picard CollectAlignmentSummaryMetrics to assess the
>                 #+ quality of read alignments as well as "the proportion of
>                 #+ reads that passed machine signal-to-noise threshold quality
>                 #+ filters"
>                 if ${verbose}; then
>                     echo "Running step #5: Generate Picard CollectAlignmentSummaryMetrics TXT file."
>                 fi
> 
>                 if [[ "${fasta##*.}" == "gz" ]]; then
>                     decomp_fasta="${fasta%.gz}"
>                     if [[ ! -f "${decomp_fasta}.fai" ]]; then
>                         #  Picard CollectAlignmentSummaryMetrics requires that
>                         #+ the FASTA file is not compressed
>                         if ${verbose}; then
>                             echo_note "Decompressing the reference genome FASTA file."
>                         fi
> 
>                         if [[ ! -f "${decomp_fasta}" ]]; then
>                             gunzip -c "${fasta}" > "${decomp_fasta}"
>                         fi
> 
>                         #  Also, Picard CollectAlignmentSummaryMetrics requires
>                         #+ that the reference genome FASTA is indexed; so,
>                         #+ generate an FAI file if necessary
>                         if ${verbose}; then
>                             echo_note "Generating an index file (FAI) for the reference genome FASTA file."
>                         fi
>                         samtools faidx -@ "${threads}" "${decomp_fasta}"
>                     fi
> 
>                     picard CollectAlignmentSummaryMetrics \
>                         --REFERENCE_SEQUENCE "${decomp_fasta}" \
>                         --INPUT "${bam_coor}" \
>                         --OUTPUT "${txt_met}"
>                 else
>                     if [[ ! -f "${fasta}.fai" ]]; then
>                         if ${verbose}; then
>                             echo_note "Generating an index file (FAI) for the reference genome FASTA file."
>                         fi
> 
>                         samtools faidx -@ "${threads}" "${fasta}"
>                     fi
> 
>                     picard CollectAlignmentSummaryMetrics \
>                         --REFERENCE_SEQUENCE "${fasta}" \
>                         --INPUT "${bam_coor}" \
>                         --OUTPUT "${txt_met}"
>                 fi
>             else
>                 echo_note "Picard CollectAlignmentSummaryMetrics TXT file already exists. Skipping CollectAlignmentSummaryMetrics operations (step #5)."
>             fi
>         else
>             error_and_exit "Duplicate alignments have not been marked in coordinate-sorted BAM file (step #5)."
>         fi
>     else
>         error_and_exit "Index (BAI file) for coordinate-sorted BAM file does not exist (step #5)."
>     fi
> else
>     error_and_exit "Coordinate-sorted BAM file does not exist (step #5)."
> fi
Running step #5: Generate Picard CollectAlignmentSummaryMetrics TXT file.
Picked up JAVA_TOOL_OPTIONS: -Djava.io.tmpdir=/loc/scratch/43305178
Error: LinkageError occurred while loading main class picard.cmdline.PicardCommandLine
	java.lang.UnsupportedClassVersionError: picard/cmdline/PicardCommandLine has been compiled by a more recent version of the Java Runtime (class file version 61.0), this version of the Java Runtime only recognizes class file versions up to 55.0
```
</details>
<br />