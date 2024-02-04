
`#notes.md`
<br />

<details>
<summary><b><font size="+2"><i>Table of contents</i></font></b></summary>
<br />
<!-- MarkdownTOC -->

1. [Notes on teaching materials at the Harvard Chan Bioinformatics Core: Intro-to-ChIPseq-flipped](#notes-on-teaching-materials-at-the-harvard-chan-bioinformatics-core-intro-to-chipseq-flipped)
	1. [lessons/05_filtering_BAM_files](#lessons05_filtering_bam_files)
		1. [Filtering reads](#filtering-reads)
		1. [Multimapping reads \(multimappers\)](#multimapping-reads-multimappers)
			1. [Conversation with ChatGPT](#conversation-with-chatgpt)
				1. [Me](#me)
				1. [ChatGPT](#chatgpt)
			1. [How multimapping reads \(multimappers\) are handled in CUT&RUN and ATAC-seq analyses](#how-multimapping-reads-multimappers-are-handled-in-cutrun-and-atac-seq-analyses)
			1. [Notes to self](#notes-to-self)
		1. [Duplicate reads](#duplicate-reads)
			1. [Conversation with ChatGPT](#conversation-with-chatgpt-1)
				1. [Me](#me-1)
				1. [ChatGPT](#chatgpt-1)
				1. [Me](#me-2)
				1. [ChatGPT](#chatgpt-2)
					1. [Explanation](#explanation)
					1. [Note](#note)
				1. [Me](#me-3)
					1. [ChatGPT](#chatgpt-3)
					1. [Explanation](#explanation-1)
					1. [Note](#note-1)
			1. [Some additional notes on duplicates](#some-additional-notes-on-duplicates)
			1. [On how duplicates are treated in CUT&RUN and ATAC-seq analyses](#on-how-duplicates-are-treated-in-cutrun-and-atac-seq-analyses)
		1. [Note to self: Takeaways for multimappers and duplicates](#note-to-self-takeaways-for-multimappers-and-duplicates)
		1. [Filtering workflow](#filtering-workflow)
			1. [Note to self: Installation of necessary programs](#note-to-self-installation-of-necessary-programs)
			1. [1. Sort BAM files by genomic coordinates](#1-sort-bam-files-by-genomic-coordinates)
			1. [2. Filter the reads to keep only uniquely mapping reads](#2-filter-the-reads-to-keep-only-uniquely-mapping-reads)
				1. [Additional filtering considerations for CUT&RUN](#additional-filtering-considerations-for-cutrun)
				1. [Additional filtering considerations for ATAC-seq data](#additional-filtering-considerations-for-atac-seq-data)
			1. [Filtering out blacklisted regions](#filtering-out-blacklisted-regions)
1. [Notes on teaching materials at teaching materials at the Harvard Chan Bioinformatics Core: Intro-to-ChIPseq](#notes-on-teaching-materials-at-teaching-materials-at-the-harvard-chan-bioinformatics-core-intro-to-chipseq)
	1. [lessons/07_handling-replicates-idr](#lessons07_handling-replicates-idr)

<!-- /MarkdownTOC -->
</details>
<br />

<a id="notes-on-teaching-materials-at-the-harvard-chan-bioinformatics-core-intro-to-chipseq-flipped"></a>
# Notes on [teaching materials at the Harvard Chan Bioinformatics Core](https://github.com/hbctraining): [Intro-to-ChIPseq-flipped](https://github.com/hbctraining/Intro-to-ChIPseq-flipped/tree/main)
<a id="lessons05_filtering_bam_files"></a>
## [lessons](https://github.com/hbctraining/Intro-to-ChIPseq-flipped/tree/main/lessons)/[05_filtering_BAM_files](https://github.com/hbctraining/Intro-to-ChIPseq-flipped/blob/main/lessons/05_filtering_BAM_files.md)
<a id="filtering-reads"></a>
### [Filtering reads](https://github.com/hbctraining/Intro-to-ChIPseq-flipped/blob/main/lessons/05_filtering_BAM_files.md#filtering-reads)
A key issue when working with a ChIP-seq data is to move forward with only the uniquely mapping reads. Allowing for multi-mapped reads increases the number of usable reads and the sensitivity of peak detection; however, the number of false positives may also increase [1](https://www.ncbi.nlm.nih.gov/pubmed/21779159/). To increase our confidence in peak calling and improve data reproducibility, we need to filter out both multi-mapping reads and duplicate reads.

<a id="multimapping-reads-multimappers"></a>
### [Multimapping reads (multimappers)](https://github.com/hbctraining/Intro-to-ChIPseq-flipped/blob/main/lessons/05_filtering_BAM_files.md#multi-mapping-reads)
Multi-mapping reads are reads that are mapping to multiple loci on the reference genome.

<a id="conversation-with-chatgpt"></a>
#### Conversation with ChatGPT
<a id="me"></a>
##### Me
How would you define "multimapping reads" if you were writing a text book?

<a id="chatgpt"></a>
##### ChatGPT
In a textbook, "multimapping reads" would be defined as follows: Multimapping reads refer to sequences in next-generation sequencing data that <i>align to multiple locations in a reference genome with similar levels of confidence</i>. This occurs due to the presence of repetitive elements, paralogous genes, or other sequences with high similarity within the genome. Such reads pose challenges in accurate alignment and quantification, as <i>it is difficult to ascertain the true origin of these sequences</i>. They are often handled separately in data analysis, with some methods distributing them proportionally across all matched locations or excluding them to increase confidence in mapping results.

<a id="how-multimapping-reads-multimappers-are-handled-in-cutrun-and-atac-seq-analyses"></a>
#### How multimapping reads (multimappers) are handled in CUT&RUN and ATAC-seq analyses
- ATAC-seq: Use of uniquely mapping reads is critical for ATAC-seq analysis. A unique mapping rate over 80% is typical for a successful experiment. Multimappers are always removed.
- CUT&RUN: There are no mentions of a multimapper removal step in CUT&RUN analysis approaches. Perhaps because Bowtie2 defaults to search for multiple alignments, yet only reports the best one (i.e., not using the `-k` option).

<a id="notes-to-self"></a>
#### Notes to self
Thus far, we're running Bowtie 2 as follows:
```bash
bowtie2 \
    -p ${threads} \
    -x "${dir_indx}" \
    --very-sensitive-local \
    --no-unal \
    --no-mixed \
    --no-discordant \
    --no-overlap \
    --no-dovetail \
    --phred33 \
    -I 10 \
    -X 700 \
    -1 "${trim_1}" \
    -2 "${trim_2}" \
        | samtools view \
            -@ ${threads} \
            -b \
            -o "${bam}"
```

We're not running in `-k` mode, which means Bowtie 2 searches for multiple alignments but only selects the single "best" (presumably it breaks ties when necessary&mdash;but I need to check on this). Thus, essentially, the multimappers "have become" unimappers and there are no secondary alignments in our data.


<a id="duplicate-reads"></a>
### [Duplicate reads](https://github.com/hbctraining/Intro-to-ChIPseq-flipped/blob/main/lessons/05_filtering_BAM_files.md#duplicate-reads)
Duplicate reads are reads that map at the exact same location, with the same coordinates and the same strand. These duplicates can arise from experimental artifacts, but can also contribute to genuine ChIP-signal.

- The bad kind of duplicates: If initial starting material is low, this can lead to overamplification of this material before sequencing. Any biases in PCR will compound this problem and can lead to artificially enriched regions.
- The good kind of duplicates: You can expect some biological duplicates with ChIP-seq since you are only sequencing a small part of the genome. This number can increase if your depth of coverage is excessive or if your protein only binds to few sites. If there are a good proportion of biological dupicates, removal can lead to an underestimation of the ChIP signal.
To get an idea on what to expect in terms of duplication rate, we encourage you to take a look at the [ENCODE quality metrics for complexity](https://www.encodeproject.org/data-standards/terms/#library). Different metrics are described and there is also a table which describes how to classify a sample as good, moderate or bad, based on these values.

<a id="conversation-with-chatgpt-1"></a>
#### Conversation with ChatGPT
<a id="me-1"></a>
##### Me
On [this page](https://www.encodeproject.org/data-standards/terms/#library), I see the ENCODE quality metrics for complexity. Here, different metrics are described and there is also a table which describes how to classify a sample as good, moderate or bad, based on these values.

![ENCODE quality metrics for library complexity](./pngs/ENCODE-library-complexity.png)

Do you know how I can run the quality metrics on my data? It seems I need the following values:
- M<sub>1</sub>: number of genomic locations where exactly one read maps uniquely
- M<sub>DISTINCT</sub>: number of distinct genomic locations to which some read maps uniquely
- M<sub>2</sub>: number of genomic locations where two reads map uniquely
- Number of distinct uniquely mapping reads (i.e., after removing duplicates)
- Total number of reads

<a id="chatgpt-1"></a>
##### ChatGPT
To calculate the ENCODE quality metrics for library complexity from your ChIP-seq or ATAC-seq data, you'll need to determine several values, such as M<sub>1</sub>, M<sub>2</sub>, and M<sub>DISTINCT</sub>, the number of distinct uniquely mapping reads, and the total number of reads. These metrics are used to compute PBC1 (PCR Bottlenecking Coefficient 1), PBC2 (PCR Bottlenecking Coefficient 2), and NRF (Non-Redundant Fraction).

Here are general steps to calculate these metrics:

1. Alignment: Align your reads to the reference genome using an aligner like BWA or Bowtie2.
2. Filtering for unique reads: Use tools like samtools to filter uniquely mapped reads (often reads with a mapping quality above a certain threshold, e.g., q30).
3. Calculate M<sub>1</sub>, M<sub>2</sub>, and M<sub>DISTINCT</sub>:
    - M<sub>1</sub>: Count genomic locations where exactly one read maps uniquely. This can be done using samtools and custom scripts or awk commands.
    - M<sub>2</sub>: Count genomic locations where exactly two reads map uniquely.
    - M<sub>DISTINCT</sub>: Count the number of distinct genomic locations to which some read maps uniquely.
4. Calculate NRF: This is the number of distinct uniquely mapping reads (after removing duplicates) divided by the total number of reads. Tools like samtools can help identify and filter out duplicates.
5. Compute PBC1 and PBC2:
    - PBC1: Calculate as M<sub>1</sub> divided by M<sub>DISTINCT</sub>.
    - PBC2: Calculate as M<sub>1</sub> divided by M<sub>2</sub>.

Tools like samtools, bedtools, and awk can be used in combination to calculate these metrics. Specific commands depend on the exact format of your data and the tools available in your environment.

For more detailed information and guidelines, you can refer to the [ENCODE Project's Standards](https://www.encodeproject.org/data-standards/terms/#library).

Please note that the exact commands and scripts to calculate these metrics can be quite specific and may require custom adaptation for your dataset.

<a id="me-2"></a>
##### Me
What would/could be code to "count genomic locations where exactly one read maps uniquely?"

<a id="chatgpt-2"></a>
##### ChatGPT
To count the number of genomic locations where exactly one read maps uniquely (M<sub>1</sub>), you can use a combination of samtools and awk. The process involves filtering uniquely mapped reads and then counting the locations with only one read. Here's an example code snippet:

```bash
#  1. Filter uniquely mapped reads (assuming SAM flag 0x4 for unmapped and 0x10
#+    for not primary alignment)
#+ 2. Extract chromosome and position
#+ 3. Sort and count unique occurrences
#+ 4. Filter locations with exactly one read

samtools view -F 0x4 -F 0x10 -q 30 yourfile.bam \
	| awk '{ print $3, $4 }' \
	| sort \
	| uniq -c \
	| awk '$1 == 1 { count++ } END { print count }'
```

<a id="explanation"></a>
###### Explanation
- `samtools view -F 0x4 -F 0x10 -q 30`: Filters out unmapped reads, secondary alignments, and reads with a mapping quality score less than 30.
- `awk '{ print $3, $4 }'`: Extracts the chromosome and start position of each read.
- `sort | uniq -c`: Sorts the coordinates and counts the number of occurrences of each unique coordinate.
- `awk '$1 == 1 { count++ } END { print count }'`: Counts the number of genomic locations where the count is exactly one.
<a id="note"></a>
###### Note
This example assumes a typical format for BAM files and a certain threshold for mapping quality (`q30`). You may need to adjust parameters based on your specific data and requirements.

<a id="me-3"></a>
##### Me
And what would/could code be to determine M<sub>DISTINCT</sub>, the count of the number of distinct genomic locations to which some read maps uniquely?

<a id="chatgpt-3"></a>
###### ChatGPT
To determine M<sub>DISTINCT</sub>, which is the count of the number of distinct genomic locations to which some read maps uniquely, you can again use samtools to filter uniquely mapped reads and then awk and sort to find unique genomic locations. Here's an example command sequence:
```bash
#  1. Filter for uniquely mapped reads
#+ 2. Extract chromosome and position
#+ 3. Sort and get unique coordinates
#+ 4. Count unique locations

samtools view -F 0x4 -F 0x10 -q 30 yourfile.bam \
	| awk '{ print $3, $4 }' \
	| sort \
	| uniq \
	| wc -l
```

<a id="explanation-1"></a>
###### Explanation
- `samtools view -F 0x4 -F 0x10 -q 30`: This filters the BAM file for uniquely mapped reads, excluding unmapped reads, secondary alignments, and reads below a certain mapping quality (here, 30).
- `awk '{ print $3, $4 }'`: This extracts the chromosome ($3) and start position ($4) of each read.
- `sort | uniq`: This sorts the coordinates and filters them to get only unique coordinates.
- `wc -l`: This counts the number of lines, which corresponds to the number of unique genomic locations.

<a id="note-1"></a>
###### Note
As with M<sub>1</sub>, ensure that the command is adjusted according to your specific data and the alignment method used.

<a id="some-additional-notes-on-duplicates"></a>
#### Some additional notes on duplicates
Most peak calling algorithms also implement methods to deal with duplicate reads. While they are commonly removed prior to peak calling, another option is to leave them now and deal with them later. Skip the duplicate filtering at this step if:
- You are planning on performing a differential binding analysis.
- <mark>You are expecting binding in repetitive regions (also, use paired-end sequencing).</mark> (<i>Note to self: with our analyses of the</i> S. cerevisiae <i>rDNA locus, we are.</i>)
- You have included UMIs into your experimental setup.

<a id="on-how-duplicates-are-treated-in-cutrun-and-atac-seq-analyses"></a>
#### On how duplicates are treated in CUT&RUN and ATAC-seq analyses
- ATAC-seq: Duplicate removal is always performed for ATAC-seq data.
- CUT&RUN: Duplicate removal is an optional step in many CUT&RUN analysis approaches. The default is usually to keep duplicates, because CUT&RUN increases the likelihood of biological duplicates. More specifically, nuclease cleavage of chromatin by its stereotypical nature is influenced by conformation of chromatin and/or nuclease bias, increasing the likelihood of identical reads that are originated from different cells. Therefore, we should remove duplicate with caution. <i>Assess the library complexity first, and then check if there are unreasonbly high amount of duplications. If not, and your experiment does not over-amplify, you might not want to remove the duplicate.</i>

<a id="note-to-self-takeaways-for-multimappers-and-duplicates"></a>
### Note to self: Takeaways for multimappers and duplicates
- Exclude or otherwise "handle" multimappers if present
- Retain duplicates if library complexity is satisfactory

<a id="filtering-workflow"></a>
### Filtering workflow
The older version of Bowtie 2 had an argument that allowed us to easily perform filtering during the alignment process. But the latest Bowtie 2 does not have this option. As a result, the filtering will be done with the use of a tool called sambamba. Sambamba is an open source tool that provides methods for working with SAM/BAM files, similar to samtools, except with faster processing times and in some cases added functionality.

This lesson will consist of two steps:
1. Sort BAM files by genomic coordinates (using samtools).
2. Filter the reads to keep only uniquely mapping reads (using sambamba). This will also remove any unmapped reads.

<a id="note-to-self-installation-of-necessary-programs"></a>
#### Note to self: Installation of necessary programs
- [sambamba via conda (bioconda)](https://anaconda.org/bioconda/sambamba)
- [samtools via conda (bioconda)](https://anaconda.org/bioconda/samtools)

<a id="1-sort-bam-files-by-genomic-coordinates"></a>
#### 1. Sort BAM files by genomic coordinates
Before we can do the filtering, we need to sort our BAM alignment files by genomic coordinates (instead of by name). To perform the sorting, we could use Samtools, a tool we previously used when coverting our SAM file to a BAM file.

The command we use this time is `samtools sort` with the parameter `-o`, indicating the path to the output file. The example code looks like below (please do not run):
```bash
# DO NOT RUN
$ cd ~/chipseq_workshop/results/bowtie2/
$ samtools sort wt_sample2_chip.bam -o wt_sample2_chip_sorted.bam
```

NOTE: To run the above code, you will need the BAM file generated from the alignment lesson. If you do not have this file, you could copy over the BAM file to your directory:
```bash
$ cp /n/groups/hbctraining/harwell-datasets/workshop_material/results/bowtie2/wt_sample2_chip.bam ~/chipseq_workshop/results/bowtie2/wt_sample2_chip.bam
```
The samtools sort code above takes about 7 min to finish. Instead of running it in the class, we have generated the output BAM file. Please copy over the BAM file to your directory:
```bash
$ cd ~/chipseq_workshop/results/bowtie2/
$ cp /n/groups/hbctraining/harwell-datasets/workshop_material/results/bowtie2/wt_sample2_chip_sorted.bam ~/chipseq_workshop/results/bowtie2/wt_sample2_chip_sorted.bam
```

We could take a glimpse of the sorted BAM file using `samtools view`:
```bash
$ samtools view wt_sample2_chip_sorted.bam | less
```

<a id="2-filter-the-reads-to-keep-only-uniquely-mapping-reads"></a>
#### 2. Filter the reads to keep only uniquely mapping reads
Next, we can filter the sorted BAM files to keep only uniquely mapping reads. We will use the sambamba view command with the following parameters:
- `-t`: number of threads(cores)
- `-h`: print SAM header before reads
- `-f`: format of output file (default is SAM)
- `-F`: set custom filter&mdash;we will be using the filter to remove duplicates, multimappers and unmapped reads.

```bash
$ sambamba view -h -t 2 -f bam \
-F "[XS] == null and not unmapped and not duplicate" \
wt_sample2_chip_sorted.bam > wt_sample2_chip_final.bam
```

We filter out unmapped reads by specifying in the filter `not unmapped`, and we filter out duplicates with `not duplicate`. Also, among the reads that are aligned, we filter out multimappers by specifying `[XS] == null`. `XS` is a tag generated by Bowtie 2 that gives an alignment score for the second-best alignment, and it is only present if the read is aligned and more than one alignment is found for the read.

<a id="additional-filtering-considerations-for-cutrun"></a>
##### Additional filtering considerations for CUT&RUN
Once the CUT&RUN sequence reads have been aligned to the genome, the resulting BAM files can be filtered by fragment size. Fragments can be divided into ≤120-bp and >120-bp fractions. For transcription factors or proteins with an expected punctate binding profile, you can use the ≤120-bp fraction, which is likely to contain binding sites. The range can be increased depending on the protein of interest, and alternatively BAM files without filtering can also be used. Example code for filtering BAM files by fragment size:
```bash
sambamba view \
	--format bam \
	--nthreads 6 \
    -F "((template_length > 0 and template_length < 120) or (template_length < 0 and template_length > -120))" "${file}" \
     | samtools view -b \
     	> "bams_sizeSelect/${s}-sizeSelect.bam"
```

<a id="additional-filtering-considerations-for-atac-seq-data"></a>
##### Additional filtering considerations for ATAC-seq data
There are two additional filtering steps that need to be performed for ATAC-seq data analysis:
1. Filtering mitochondrial reads.
    - The mitochondrial genome, which is more accessible due to the lack of chromatin packaging will result in extremely high read coverage.
    - These reads should be discarded.
    - Since there are no ATAC-seq peaks of interest in the mitochondrial genome, these reads are discarded.
    - The Omni-ATAC method uses detergents to remove mitochondria from the samples prior to sequencing and is another option to deal with this issue.
2. Filtering BAM files based on fragment size.
    - Typically, a successful ATAC-seq experiment should generate a fragment size distribution plot with decreasing and periodical peaks corresponding to the nucleosome-free regions (NFR) (<100 bp) and mono-, di-, and tri-nucleosomes (~200, 400, and 600 bp, respectively)
    - Fragments from the NFR are expected to be enriched around the transcription start site (TSS).
    - Fragments from nucleosome-bound regions are expected to be depleted at TSS with a slight enrichment of flanking regions around TSS.
    - Use example sambamba code from above (CUT&RUN filtering) to filter out fragments by size.
    	+ Create BAM files for NFR, mono-nuc, di-nuc, tri-nuc.
    	+ Typically, the NFR BAM is used for peak calling.
    - Shift the reads in the BAM file. To account for the 9-bp duplication created by DNA repair of the nick by Tn5 transposase, reads should be shifted +4 bp and −5 bp for, respectively, the positive and negative strands.

<a id="filtering-out-blacklisted-regions"></a>
#### Filtering out blacklisted regions
Although we do not perform this step, it is common practice to apply an additional level of filtering to our BAM files. That is, we remove alignments that occur within defined blacklisted regions. <mark>We will filter out blacklist regions post-peak calling</mark>.

Blacklisted regions represent artifact regions that tend to show artificially high signal&dash;i.e., excessive unstructured anomalous reads mapping. These regions are often found at specific types of repeats such as centromeres, telomeres, and satellite repeats, and typically appear uniquely mappable so the simple mappability filters applied above do not remove them. The ENCODE and modENCODE consortia have compiled blacklists for various species and genome versions including human, mouse, worm and fly. These blacklisted regions (coordinate files) can be filtered out from our alignment files before proceeding to peak calling.

If we wanted to filter blacklist regions at this point in our workflow, we would use the following code:
```bash
# DO NOT RUN
$ bedtools intersect -v \
	-abam wt_sample2_chip_final.bam \
	-b mm10-blacklist.v2.bed \
		> wt_sample2_chip_final_blacklist_filtered.bam
```

bedtools is a suite of tools that we will discuss in more detail in a later lesson when blacklist filtering is applied.

<a id="notes-on-teaching-materials-at-teaching-materials-at-the-harvard-chan-bioinformatics-core-intro-to-chipseq"></a>
# Notes on teaching materials at [teaching materials at the Harvard Chan Bioinformatics Core](https://github.com/hbctraining): [Intro-to-ChIPseq](https://github.com/hbctraining/Intro-to-ChIPseq)
<a id="lessons07_handling-replicates-idr"></a>
## [lessons](https://github.com/hbctraining/Intro-to-ChIPseq/blob/master/schedule/2-day.md)/[07_handling-replicates-idr](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/07_handling-replicates-idr.html)
