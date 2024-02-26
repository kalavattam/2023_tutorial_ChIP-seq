
`#notes.md`
<br />

<details>
<summary><b><font size="+2"><i>Table of contents</i></font></b></summary>
<br />
<!-- MarkdownTOC -->

1. [Open tabs: what they are, what they contain, why they're important or useful, etc.](#open-tabs-what-they-are-what-they-contain-why-theyre-important-or-useful-etc)
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
1. [Open tabs related to the filtering of duplicates](#open-tabs-related-to-the-filtering-of-duplicates)
1. [Notes on MACS3 documentation](#notes-on-macs3-documentation)
    1. [Advanced step-by-step peak calling using MACS3 commands](#advanced-step-by-step-peak-calling-using-macs3-commands)
        1. [Step 1: Filter duplicates](#step-1-filter-duplicates)
            1. [Documentation for `macs3 filterdup`](#documentation-for-macs3-filterdup)
            1. [My proposed invocation of `macs3 filterdup`](#my-proposed-invocation-of-macs3-filterdup)
                1. [Code](#code)
                1. [Printed](#printed)
            1. [Conversation with GPT4](#conversation-with-gpt4)
                1. [Me](#me-4)
                1. [GPT4](#gpt4)
                    1. [What is the binomial test?](#what-is-the-binomial-test)
                    1. [Why is it applicable in this circumstance?](#why-is-it-applicable-in-this-circumstance)
                    1. [Why not remove all duplicates wholesale?](#why-not-remove-all-duplicates-wholesale)
                    1. [Conclusion](#conclusion)
                1. [Me](#me-5)
                1. [GPT4](#gpt4-1)
                    1. [Scenario](#scenario)
                    1. [Binomial test in R](#binomial-test-in-r)
                    1. [Interpretation](#interpretation)
                    1. [Note](#note-2)
                1. [Me](#me-6)
                1. [GPT4](#gpt4-2)
                    1. [Null hypothesis \(H0\)](#null-hypothesis-h0)
                    1. [Alternative Hypothesis \(H1\)](#alternative-hypothesis-h1)
                1. [Additional conversation](#additional-conversation)
        1. [Step 2: Decide the fragment length `d`](#step-2-decide-the-fragment-length-d)
            1. [Documentation for `macs predictd`](#documentation-for-macs-predictd)
            1. [My proposed invocation](#my-proposed-invocation)
            1. [Conversation with GPT4](#conversation-with-gpt4-1)
                1. [Me](#me-7)
                1. [GPT4](#gpt4-3)
        1. [Step 3: Extend ChIP sample to get ChIP coverage track](#step-3-extend-chip-sample-to-get-chip-coverage-track)
            1. [Documentation for `macs3 pileup`](#documentation-for-macs3-pileup)
            1. [On trying to understand what `--format BAMPE` does in `macs3 pileup`](#on-trying-to-understand-what---format-bampe-does-in-macs3-pileup)
            1. [My proposed invocation of `macs3 pileup`](#my-proposed-invocation-of-macs3-pileup)
            1. [Step 4: Build a local bias track from the control sample](#step-4-build-a-local-bias-track-from-the-control-sample)
                1. [The `d` background](#the-d-background)
                1. [The slocal background](#the-slocal-background)
                1. [The llocal background](#the-llocal-background)
                1. [The genome background](#the-genome-background)
                1. [Combine and generate the maximum background noise](#combine-and-generate-the-maximum-background-noise)
                    1. [On how to do this with paired-end data](#on-how-to-do-this-with-paired-end-data)
                    1. [d2 was defined here...](#d2-was-defined-here)
                    1. [Note to self](#note-to-self)

<!-- /MarkdownTOC -->
</details>
<br />

<a id="open-tabs-what-they-are-what-they-contain-why-theyre-important-or-useful-etc"></a>
# Open tabs: what they are, what they contain, why they're important or useful, etc.
- [nf-core ChIP-seq (2.0.0)](https://nf-co.re/chipseq/2.0.0): The nf-core ChIP-seq pipeline; I'm drawing on this for the processing of alignments (e.g., marking duplicates, etc.), the estimation of library complexity, the creation of RPM-normalized BW files, etc. It is essentially an example gameplan for the many steps we want to include and document in the in-progress tutorial.
- [nf-core MNase-seq (1.0.0)](https://nf-co.re/mnaseseq/1.0.0): This pipeline serves as a supplemental reference for the nf-core ChIP-seq pipeline, which is described below.
- [nf-core ChIP-seq pipeline Github repository search results for "macs"](https://github.com/search?q=repo%3Anf-core%2Fchipseq+macs&type=code): I sought to understand how the authors of the nf-core ChIP-seq pipeline called macs2. From [here](https://github.com/nf-core/chipseq/blob/51eba00b32885c4d0bec60db3cb0a45eb61e34c5/modules/nf-core/modules/macs2/callpeak/main.nf), we see the answer:
```python
#!/usr/bin/env python3

"""
macs2 \\
    callpeak \\
    ${args_list.join(' ')} \\
    --gsize $macs2_gsize \\
    --format $format \\
    --name $prefix \\
    --treatment $ipbam \\
    $control
"""
```
- [nf-core ChIP-seq pipeline GitHub repository search results for "bedgraphtobigwig"](https://github.com/search?q=repo%3Anf-core%2Fchipseq+bedgraphtobigwig&type=code): I was examining these results to understand how the authors of the nf-core ChIP-seq pipeline went about generating RPM-normalized coverage tracks in the form of BW files.
- [nf-core ChIP-seq pipeline GitHub repository search results for "deepTools"](https://github.com/search?q=repo%3Anf-core%2Fchipseq+deepTools&type=code): I was searching for the term "deepTools" to understand how the authors of nf-core ChIP-seq pipeline used the program in their pipeline; indeed, they seem to use it in different ways at different steps. I think I initially was looking into this for the sake of running quality-check analyses, such as [plotFingerprint](https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html), but see now how this can be useful for describing and running other analyses&mdash;such as [computeMatrix](https://deeptools.readthedocs.io/en/develop/content/tools/computeMatrix.html), [plotHeatmap](https://deeptools.readthedocs.io/en/develop/content/tools/plotHeatmap.html), and [plotProfile](https://deeptools.readthedocs.io/en/develop/content/tools/plotProfile.html)&mdash;downstream in the tutorial.
- [ENCODE (old) quality metrics](https://genome.ucsc.edu/ENCODE/qualityMetrics.html): This is an outdated ENCODE material on the quality metrics used in the first phase of the project (2003-2012). It includes spreadsheets of quality metrics for datasets in the first phase and definitions of the quality metrics. The updated ENCODE quality metrics can be found [here](https://www.encodeproject.org/data-standards/terms/).
- [hbctraining (new) 06_peak_calling_macs](https://github.com/hbctraining/Intro-to-ChIPseq-flipped/blob/main/lessons/06_peak_calling_macs.md): The Harvard Chan Bioinformatics Core lesson on calling peaks with MACS2. This resource goes over (briefly) how MACS2 works, including associated theory. There is a discussion of important parameters and how those parameters should be changed in case paired-end reads are used, in case we are analyzing CUT&RUN data, and in case we are analyzing ATAC-seq data. This is an important resource for our work. Here is an example of how MACS2 is called [in the tutorial](https://github.com/hbctraining/Intro-to-ChIPseq-flipped/blob/main/lessons/06_peak_calling_macs.md#macs2-parameters):
```bash
#!/bin/bash

#  For single-end ChIP-seq data
macs2 callpeak \
    -t /n/groups/hbctraining/harwell-datasets/workshop_material/results/bowtie2/wt_sample1_chip_final.bam \
    -c /n/groups/hbctraining/harwell-datasets/workshop_material/results/bowtie2/wt_sample1_input_final.bam \
    -f BAM \
    -g mm \
    -n wt_sample1 \
    --outdir macs2 \
        2> macs2/wt_sample1_macs2.log
```
- [hbctraining (old) two-day schedule](https://hbctraining.github.io/Intro-to-ChIPseq/schedule/2-day.html): This is a table of contents for the various modules (lessons) in the Harvard Chan Bioinformatics Core Intro to ChIP-seq course (an older version of the course referenced above). This is a good supplemental reference for the newer course (again, referenced above).
- [hbctraining (old) 05_peak_calling_macs](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html): The older Harvard Chan Bioinformatics Core Intro to ChIP-seq course on calling peaks with MACS2. This is a good supplemental reference for the newer course material on calling peaks with MACS2. Found here are some details not found in the older material. Here is an example call to MACS2 from the tutorial:
```bash
#!/bin/bash

#  For single-end ChIP-seq data
macs2 callpeak \
    -t bowtie2/H1hesc_Nanog_Rep1_aln.bam \
	-c bowtie2/H1hesc_Input_Rep1_aln.bam \
 	-f BAM -g 1.3e+8 \
	-n Nanog-rep1 \
	--outdir macs2 \
	    2> macs2/Nanog-rep1-macs2.log
```
- [hbctraining (old) 03_align_and_filtering](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/03_align_and_filtering.html): The older Harvard Chan Bioinformatics Core Intro to ChIP-seq course on aligning reads with Bowtie 2, understanding the information in SAM/BAM files, and filtering aligned reads to keep only uniquely mapped ones.
- [hbctraining (old)  06_combine_chipQC_and_metrics](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/06_combine_chipQC_and_metrics.html): The older Harvard Chan Bioinformatics Core Intro to ChIP-seq course on generating a report containing QC metrics using [ChIPQC](https://bioconductor.org/packages/release/bioc/html/ChIPQC.html). The lesson material describes specific QC metrics, how to interpret values associated with them, and points out sources of low quality data.
- [hbctraining (old)  07_handling-replicates-idr](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/07_handling-replicates-idr.html): The older Harvard Chan Bioinformatics Core Intro to ChIP-seq course on combining replicates to only get the highly reproducible peaks using the IDR method. This briefly goes over different strategies for handling replicates in ChIP-seq before diving into the IDR method, including (a) the evaluation of peak consistency between true replicates, (b) the evaluation peak consistency between pooled pseudo-replicates, and (c) the evaluation self-consistency for each individual replicate. Note that running the IDR method entails running MACS2/3 in a less stringent manner, as the "IDR algorithm requires sampling of both signal and noise distributions to separate the peaks into two groups". Code for (a) is provided with descriptions; ~~only descriptions~~ code and descriptions are provided for (b), and only descriptions are provided for (c); however, the descriptions give sufficient detail such that I could draft code to implement the work without much trouble.
- [Biostar Handbook ChIP-seq alignment](https://www.biostarhandbook.com/chip-seq-alignments.html): Part of the Biostars Handbook, this is a detailed, systematic lesson on aligning ChIP-seq data. It's pretty basic, but provides context and details.
- [Biostar Handbook ChIP-seq peak-calling](https://www.biostarhandbook.com/chip-peak-calling.html): Also part of the Biostars Handbook, the information on peak-calling here is detailed but basic; what is useful is the discussion and code related to the use of [sdust](https://github.com/lh3/sdust) (for determining low-complexity regions) and [bedtools](https://bedtools.readthedocs.io/en/latest/) (for keeping high-scoring values that do not fall on low-complexity regions and selecting genomic regions upstream of genes). Here's an example call to MACS2 from this resource:
```bash
#!/bin/bash

macs2 callpeak \
    -t "${ETH1}" "${ETH2}" \
    -c "${GLU1}" "${GLU2}" \
    --gsize 1E7 \
    --name ethanol \
    --outdir ethanol/ 
```
- [Biostar Handbook Ming Tang's Guide to ChIP-seq Analysis](https://www.biostarhandbook.com/ming-tangs-guide-to-chip-seq-analysis.html): Also part of the Biostars Handbook, I'm including this link here for quick access to the ChIP-seq analysis tutorial by Ming "Tommy" Tang. However, this resource only discusses the outdated MACS program, not MACS2/3.
- [crazyhottommy's ChIP-seq analysis notes: Peak calling with MACS2](https://github.com/crazyhottommy/ChIP-seq-analysis/blob/master/part1.3_MACS2_peak_calling_details.md): From Ming "Tommy" Tang's ChIP-seq notes, some details on/about setting MACS2 parameters&mdash;particularly distinct parameters for calling narrow and broad/gapped peaks. There is [a similar, subsequent page covering this material in crazyhottommy's ChIP-seq analysis notes too](https://github.com/crazyhottommy/ChIP-seq-analysis/blob/master/part1_peak_calling.md):
```bash
#!/bin/bash

#  Old ENCODE approach to generating narrow peaks and preliminary signal tracks
macs2 callpeak \
    -t %s -c %s %(experiment.name, control.name) \
	-f BED -n %s/%s %(peaks_dirname, prefix) \
	-g %s %(genomesize) \
	-p 1e-2 \
	--nomodel \
	--shift 0 \
	--extsize %s %(fraglen) \
	--keep-dup all \
	-B \
	--SPMR

#  Old ENCODE approach to generating broad and gapped peaks
macs2 callpeak \
    -t %s -c %s %(experiment.name, control.name) \
    -f BED -n %s/%s %(peaks_dirname, prefix) \
    -g %s %(genomesize) \
    -p 1e-2 \
    --broad \
    --nomodel \
    --shift 0 \
    --extsize %s %(fraglen) \
    --keep-dup all

#  Use --nomodel and provide the --extsize of either 147 bp or the fragment
#+ length predicted by strand cross-correlation analysis
#+
#+ For narrow peaks
macs2 callpeak \
    -t IP.bam \
    -c Input.bam \
    -n test \
    -p 0.01 \
    --nomodel \
    --extsize fragment_length \
    --keep-dup all \
    -g hs

#+ For borad regions
macs2 callpeak \
    -t IP.bam \
    -c Input.bam \
    -n test \
    --broad \
    -p 0.01 \
    --nomodel \
    --extsize fragment_length \
    --keep-dup all \
    -g hs

#  Call peaks with macs2 using --broad, building model
macs2 callpeak \
    -t ../data/wgEncodeSydhHistonePanc1H3k27acUcdAlnRep1.bam \
    -c ../data/wgEncodeSydhHistonePanc1InputUcdAlnRep1.bam \
    --broad \
    -g hs \
    --broad-cutoff 0.1 \
    -n panc1H3k27acRep1 \
    --outdir panc1H3k27acRep1_with_model_broad

#  Call peaks with macs2 using --broad, bypass the model
macs2 callpeak \
    -t ../data/wgEncodeSydhHistonePanc1H3k27acUcdAlnRep1.bam \
    -c ../data/wgEncodeSydhHistonePanc1InputUcdAlnRep1.bam \
    --broad \
    -g hs \
    --broad-cutoff 0.1 \
    -n panc1H3k27acRep1 \
    --outdir panc1H3k27acRep1_without_model_broad \
    --nomodel \
    --extsize 146
#  --nomodel and --extsize 146 tell MACS2 to use 146 bp as the fragment size
#+ for the piling up of alignments
```
- [Description of ChiLin quality metrics](http://cistrome.org/chilin/_downloads/instructions.pdf): Taken from [ChiLin ChIP-seq data quality and analysis pipeline](http://cistrome.org/chilin/), this PDF file provides a well-organized, thorough overview of useful quality metrics for the assessment of ChIP-seq datasets. (Currently, we calculate many of these metrics; we may seek to obtain more.)
- [ENCODE4 histone ChIP-seq data standards and processing pipeline (2020-07)](https://www.encodeproject.org/chip-seq/histone-encode4/#standards): This page provides ENCODE assay and pipeline overviews for both replicated and unreplicated experiments, as well as reference files, pipeline guidelines (e.g., minimal lengths of reads, etc.), and quality-metric standards (e.g., the preferred values for NRF, PBC1, and PBC2 are as follows: NRF > 0.9, PBC1 > 0.9, and PBC2 > 10; assuming the use of human or mouse model organisms, the minimum ENCODE standard for each replicate in a ChIP-seq experiment targeting broad marks and investigated as a broad histone mark is 20 million usable fragments, meaning the recommended value sequencing depth is >45 million paired-end reads, although >35 million is acceptable; etc.).
- [GitHub repository for the ENCODE ChIP-seq pipeline (2.0)](https://github.com/ENCODE-DCC/chip-seq-pipeline2). An example call to MACS2 from [this script for ChIP-seq](https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/26eeda81a0540dc793fc69b0c390d232ca7ca50a/src/encode_task_macs2_signal_track_chip.py#L90) and [this script for ATAC-seq](https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/master/src/encode_task_macs2_atac.py):
```python
#!/usr/bin/env python3

#  ChIP-seq
run_shell_cmd(
    ' macs2 callpeak '
    '-t {ta} {ctl_param} -f BED -n {prefix} -g {gensz} -p {pval_thresh} '
    '--nomodel --shift {shiftsize} --extsize {extsize} --keep-dup all -B --SPMR'.format(
        ta=ta,
        ctl_param='-c {ctl_ta}'.format(ctl_ta=ctl_ta) if ctl_ta else '',
        prefix=prefix,
        gensz=gensz,
        pval_thresh=pval_thresh,
        shiftsize=0,
        extsize=fraglen,
    )
)

#  ATAC-seq
run_shell_cmd(
        'macs2 callpeak '
        '-t {ta} -f BED -n {prefix} -g {gensz} -p {pval_thresh} '
        '--shift {shiftsize} --extsize {extsize} '
        '--nomodel -B --SPMR --keep-dup all --call-summits'.format(
            ta=ta,
            prefix=prefix,
            gensz=gensz,
            pval_thresh=pval_thresh,
            shiftsize=shiftsize,
            extsize=smooth_win,
        )
    )
```
- [F1000 article: "From reads to regions: a Bioconductor workflow to detect differential binding in ChIP-seq data"](https://f1000research.com/articles/4-1080): This article is a tutorial for ChIP-seq differential binding analyses by authors Aaron Lun and Gordon Smyth. "Most conventional approaches to ChIP-seq data analysis involve the detection of the absolute presence (or absence) of a binding site. However, an alternative strategy is to identify changes in the binding intensity between two biological conditions, i.e., differential binding (DB). This may yield more relevant results than conventional analyses, as changes in binding can be associated with the biological difference being investigated. The aim of this article is to facilitate the implementation of DB analyses, by comprehensively describing a computational workflow for the detection of DB regions from ChIP-seq data. ... \[D\]etection of DB regions will be conducted using the counts for sliding windows from the csaw package, with statistical modelling performed using methods in the edgeR package."
- [NBIS Epigenomics Workshop 2023: ChIP-seq](https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/chipseq_tutorials.html): Very detailed material on ChIP-seq, including...
	+ [data processing, quality checks, and narrow peak detection](https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/chipseqProc/lab-chipseq-processing.html)
	+ [broad peak detection](https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/chipseqBroadPeaks/lab-broadpeaks.html)
	+ [peak annotation](https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/atac-chip-downstream/lab-PeakAnnot.html)
	+ [the detection of differential occupancy in ChIP-seq data using DiffBind](https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/diffBind/lab-diffBinding.html)
	
```bash
#!/bin/bash

#  Narrow peak calling
macs2 callpeak \
    -t ENCFF000PED.preproc.bam \
    -c ENCFF000PET.preproc.bam \
    -f BAM \
    -g 4.9e8 \
    -n hela_1_REST.chr12.macs2 \
    -q 0.01 \
        &> macs.log

#  Broad peak calling
macs3 callpeak \
    --broad \
    -t ENCFF395DAJ.chr12.MAPQ30.blcklst.rh.sorted.bam \
    -c ENCFF956GLJ.chr12.MAPQ30.blcklst.rh.sorted.bam \
    -f BAMPE \
    -g 04.9e8 \
    --broad-cutoff 0.1 \
    -n neuroGM23338_macs3_rep1
```

- [Documentation for the Bioconductor package ChIPSeqSpike](https://bioconductor.org/packages/3.12/bioc/vignettes/ChIPSeqSpike/inst/doc/ChIPSeqSpike.pdf): "ChIPSeqSpike: ChIP-seq data scaling with spike-in control"
	+ [NBIS Epigenomic Workshop 2023 entry on ChIPSeqSpike](https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/quantitativeChip/lab-exospike.html)
	+ I might be wrong, but I think modern versions of this package are no longer maintained as R and Bioconductor are updated.
- [Documentation for the Bioconductor package ChIPQC](https://bioconductor.org/packages/3.18/bioc/vignettes/ChIPQC/inst/doc/ChIPQC.pdf): Definitely draw on this.
- [PDF presentation on ChIP-seq peak calling by Stéphanie Le Gras at IGBMC (France)](http://genomeast.igbmc.fr/wiki/lib/exe/fetch.php?media=training:1709ngsbioinfo_12a_chipseq_peakcalling.pdf): Includes useful information on the use of spike-ins.
- The example [usage page]((https://deeptools.readthedocs.io/en/develop/content/example_usage.html)) in the deepTools documentation.
- [Biostars post: "Purpose of bamCoverage RPKM normalization method"](https://www.biostars.org/p/9474318/): I never understood why RPKM normalization is used by some in ChIP-seq analyses (let alone in RNA-seq analyses). This post includes a comment from ATPoint echoing my concerns and advancing that, in this context, use of the CPM normalization seems more appropriate than use of the RPKM normalization.
- [Biostars post: "How to plot coverage and depth statistics of a bam file"](https://www.biostars.org/p/104063/): Here, lots of advice on how to go about calculating coverage and depth statistics for BAM files. Because of information presented here and elsewhere, I have moved forward with [mosdepth](https://github.com/brentp/mosdepth) for calculating coverage and depth statistics.
- [MACS2/3 issue: "How to incorporate ChIP-seq spike in normalisation into peak calling?"](https://github.com/macs3-project/MACS/issues/356): Not very scrutable advice on including ChIP-seq spike-in scaling factors in MACS2/3 peak calling. I want to better understand, implement, and test this approach, perhaps including it in the final tutorial (or at least the linked notebook).
<br />
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
Multi-mapping reads are reads that map to multiple loci on the reference genome.

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

We're not running in `-k` mode, which means Bowtie 2 searches for multiple alignments but only selects the single "best" (presumably it breaks ties when necessary&mdash;but I need to check on this). Thus, essentially, the multimappers "have become" unimappers and there are no secondary alignments in our data. Still, we can identify multimappers in SAM/BAM files by the presence of the XS tag, "a tag generated by Bowtie2 that gives an alignment score for the second-best alignment;" it is "only present if the read is aligned and more than one alignment is found for the read."

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
On [this page](https://www.encodeproject.org/data-standards/terms/#library), I see the ENCODE quality metrics for complexity. Here, different metrics are described and there is also a table which describes how to classify a sample as good, moderate, or bad, based on these values.

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
<br />
<br />

<a id="notes-on-teaching-materials-at-teaching-materials-at-the-harvard-chan-bioinformatics-core-intro-to-chipseq"></a>
# Notes on teaching materials at [teaching materials at the Harvard Chan Bioinformatics Core](https://github.com/hbctraining): [Intro-to-ChIPseq](https://github.com/hbctraining/Intro-to-ChIPseq)
<a id="lessons07_handling-replicates-idr"></a>
## [lessons](https://github.com/hbctraining/Intro-to-ChIPseq/blob/master/schedule/2-day.md)/[07_handling-replicates-idr](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/07_handling-replicates-idr.html)
<br />
<br />

<a id="open-tabs-related-to-the-filtering-of-duplicates"></a>
# Open tabs related to the filtering of duplicates
*`#TODO` Take notes on this stuff, similar to what you did above*
- [The Duplicates Dilemma of ChIP-seq](https://www.biostars.org/p/9543379/)
- [Tian et al. (Yan, Slager), Identification of factors associated with duplicate rate in ChIP-seq data, PLoS ONE 2019-0403](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0214723)
- [If we use MACS2, do we need to remove duplicate sequences with `samtools rmdup`?](https://www.biostars.org/p/300640/)
- [Did you remove ChIP-seq duplicates?](https://www.biostars.org/p/318974/)
- [Duplicated reads in ChIP-seq](https://www.biostars.org/p/225609/)
- [Handing duplicates in DiffBind](https://support.bioconductor.org/p/118792/)
- [Leleu, Lefebrve, Rougemont, Processing and analyzing ChIP-seq data: from short reads to regulatory interactions, Brief Funct Genomics, 2010-12](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3080774/)
- [To filter or not to filter duplicate reads – ChIP-seq](https://ethanomics.wordpress.com/2012/01/06/to-filter-or-not-to-filter-duplicate-reads-chip-seq/)
- [Duplicate Reads in ChIP-seq](https://www.seqanswers.com/forum/applications-forums/sample-prep-library-generation/39917-duplicate-reads-in-chip-seq)
- [Removing duplicates is it really necessary?](https://www.seqanswers.com/forum/bioinformatics/bioinformatics-aa/5774-removing-duplicates-is-it-really-necessary)
- [csaw - workflow to incorporate input/control samples?](https://support.bioconductor.org/p/82099/)
- [csaw: when not to dedup](https://support.bioconductor.org/p/82288/)
- [csaw: removing PCR duplicates using dedup flag](https://support.bioconductor.org/p/87874/)
- [csaw: when not to dedup](https://support.bioconductor.org/p/82288/)

*To be mentioned elsewhere; this comment from Tao Liu seems quite relevant to the above searches* ([source](https://github.com/macs3-project/MACS/issues/145), specifically [here](https://github.com/macs3-project/MACS/issues/145#issuecomment-662532470)):
> `--keep-dup all` is recommended if the duplicated alignments have been removed in preprocessing steps with other tools. The idea is that we shall not perform the same task with two different tools in our pipeline. For example, if you rely on picard to remove duplicates, then you just leave the task to picard and let macs not to remove duplicates at all. In this way, you will have better control of your pipeline and avoid unexpected errors.
<br />
<br />

<a id="notes-on-macs3-documentation"></a>
# Notes on [MACS3 documentation](https://github.com/macs3-project/MACS)
*I've edited the initial writing of Tao Liu for clarity.*

<a id="advanced-step-by-step-peak-calling-using-macs3-commands"></a>
## [Advanced step-by-step peak calling using MACS3 commands](https://github.com/macs3-project/MACS/blob/master/docs/Advanced_Step-by-step_Peak_Calling.md)
Over the years, I have gotten many emails from users asking if they can analyze their "<i>x</i>-seq" (not ChIP-seq) data with MACS, or if they can turn on or off some of the features in `callpeak`. In most cases, I would simply reply that they may have to find a dedicated tool for the various data in question, because the `callpeak` module is specifically designed and tuned for ChIP-seq data. However, MACS3 in fact contains a suite of subcommands, and one can control every single step with a self-implemented pipeline that combines them, analyzing the data of interest in a highly customized way. In this tutorial, I show how `callpeak`, the main MACS3 function, can be decomposed into a pipeline containing various MACS3 subcommands, including `filterdup`, `predictd`, `pileup`, `bdgcmp`, `bdgopt`, and `bdgpeakcall` (or `bdgbroadcall` in the case of broad marks). To analyze special data in a special way, some of these steps may need to be skipped, or some of their parameters may need to be tweaked.

Now, let's suppose we are dealing with the two test files, CTCF_ChIP_200K.bed.gz and CTCF_Control_200K.bed.gz, that are found in the MACS3 GitHub repository.

<i>Note: This tutorial is written with single-end datasets in mind. For use with for paired-end data, the various calls to the command line will be need to be modified.</i>

<a id="step-1-filter-duplicates"></a>
### Step 1: Filter duplicates
In the first step of ChIP-seq analysis by `callpeak`, ChIP and control data need to be read and the redundant reads at each genomic loci have to be removed. I won't go over the rationale, but just tell you how this can be done by filterdup subcommand. By default, the maximum number of allowed duplicated reads is 1, or `--keep-dup=1` for callpeak. To simulate this behavior, do the following:
```bash
#!/bin/bash

macs3 filterdup \
    -i CTCF_ChIP_200K.bed.gz \
    --keep-dup=1 \
    -o CTCF_ChIP_200K_filterdup.bed

macs3 filterdup \
    -i CTCF_Control_200K.bed.gz \
    --keep-dup=1 \
    -o CTCF_Control_200K_filterdup.bed
```

You can set different number for `--keep-dup` or let MACS3 automatically decide the maximum allowed duplicated reads for each genomic loci for ChIP and control separately. Check `macs3 filterdup -h` for detail, and remember if you go with `auto` way, the genome size should be set accordingly. Note, in the output, MACS3 will give you the final number of reads kept after filtering, you'd better write the numbers down since we need them when we have to scale the ChIP and control signals to the same depth. In this case, the number is 199,583 for ChIP and 199,867 for control, and the ratio between them is 199,583/199,867 = 0.99858.

<a id="documentation-for-macs3-filterdup"></a>
#### Documentation for `macs3 filterdup`
<details>
<summary><i>Text: Documentation for `macs3 filterdup`</i></summary>

```txt
#!/bin/bash

❯ macs3 filterdup -h
usage: macs3 filterdup [-h] -i IFILE [IFILE ...] [-f {AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE,BEDPE}] [-g GSIZE]
                       [-s TSIZE] [-p PVALUE] [--keep-dup KEEPDUPLICATES] [--buffer-size BUFFER_SIZE] [--verbose VERBOSE] [--outdir OUTDIR]
                       [-o OUTPUTFILE] [-d]

options:
  -h, --help            show this help message and exit
  -i IFILE [IFILE ...], --ifile IFILE [IFILE ...]
                        Alignment file. If multiple files are given as '-t A B C', then they will all be read and combined. REQUIRED.
  -f {AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE,BEDPE}, --format {AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE,BEDPE}
                        Format of tag file, "AUTO", "BED" or "ELAND" or "ELANDMULTI" or "ELANDEXPORT" or "SAM" or "BAM" or "BOWTIE" or
                        "BAMPE" or "BEDPE". The default AUTO option will let 'macs3 filterdup' decide which format the file is. Please check
                        the definition in README file if you choose ELAND/ELANDMULTI/ELANDEXPORT/SAM/BAM/BOWTIE or BAMPE/BEDPE. DEFAULT:
                        "AUTO"
  -g GSIZE, --gsize GSIZE
                        Effective genome size. It can be 1.0e+9 or 1000000000, or shortcuts:'hs' for human (2,913,022,398), 'mm' for mouse
                        (2,652,783,500), 'ce' for C. elegans (100,286,401) and 'dm' for fruitfly (142,573,017), Default:hs. The effective
                        genome size numbers for the above four species are collected from Deeptools
                        https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html Please refer to deeptools to
                        define the best genome size you plan to use.
  -s TSIZE, --tsize TSIZE
                        Tag size. This will override the auto detected tag size. DEFAULT: Not set
  -p PVALUE, --pvalue PVALUE
                        Pvalue cutoff for binomial distribution test. DEFAULT:1e-5
  --keep-dup KEEPDUPLICATES
                        It controls the 'macs3 filterdup' behavior towards duplicate tags/pairs at the exact same location -- the same
                        coordination and the same strand. The 'auto' option makes 'macs3 filterdup' calculate the maximum tags at the exact
                        same location based on binomal distribution using given -p as pvalue cutoff; and the 'all' option keeps every tags
                        (useful if you only want to convert formats). If an integer is given, at most this number of tags will be kept at the
                        same location. Note, MACS3 callpeak function uses KEEPDUPLICATES=1 as default. Note, if you've used samtools or
                        picard to flag reads as 'PCR/Optical duplicate' in bit 1024, MACS3 will still read them although the reads may be
                        decided by MACS3 as duplicate later. Default: auto
  --buffer-size BUFFER_SIZE
                        Buffer size for incrementally increasing internal array size to store reads alignment information. In most cases, you
                        don't have to change this parameter. However, if there are large number of chromosomes/contigs/scaffolds in your
                        alignment, it's recommended to specify a smaller buffer size in order to decrease memory usage (but it will take
                        longer time to read alignment files). Minimum memory requested for reading an alignment file is about # of CHROMOSOME
                        * BUFFER_SIZE * 8 Bytes. DEFAULT: 100000
  --verbose VERBOSE     Set verbose level. 0: only show critical message, 1: show additional warning message, 2: show process information, 3:
                        show debug messages. If you want to know where are the duplicate reads, use 3. DEFAULT:2
  --outdir OUTDIR       If specified all output files will be written to that directory. Default: the current working directory
  -o OUTPUTFILE, --ofile OUTPUTFILE
                        Output BED file name. If not specified, will write to standard output. Note, if the input format is BAMPE or BEDPE,
                        the output will be in BEDPE format. DEFAULT: stdout
  -d, --dry-run         When set, filterdup will only output numbers instead of writing output files, including maximum allowable duplicates,
                        total number of reads before filtering, total number of reads after filtering, and redundant rate. Default: not set
```
</details>
<br />

<a id="my-proposed-invocation-of-macs3-filterdup"></a>
#### My proposed invocation of `macs3 filterdup`
Based on the above information, this is the command I propose to call:
<a id="code"></a>
##### Code
<details>
<summary><i>Code: My proposed invocation of `macs3 filterdup`</i></summary>

```bash
#!/bin/bash

file_bam="path/to/file.bam"
gsize=11624332
dir_out="path/for/outfiles"
file_out="file.bed"

macs3 filterdup \
    --ifile "${file_bam}" \
    --format "BAMPE" \
    --gsize "${gsize}" \
    --pvalue "1e-5" \
    --keep-dup "auto" \
    --buffer-size 100000 \
    --outdir "${dir_out}" \
    --ofile "${file_out}" \
    --verbose 3

#  Also, note that the use of --dry-run could make for an interesting and
#+ informative quality control metric

#  Rough-draft test code for IP, input
file_IP_bam="IP_Q_Esa5_7041.sort-coord.bam"
file_IP_bed="IP_Q_Esa5_7041.sort-coord.bed"
file_in_bam="in_${file_IP_bam##IP_}"
file_in_bed="${file_in/.bam/.bed}"
gsize=11624332
dir_out="macs3_test"

source activate macs3_env

cd "${HOME}/projects-etc/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/03_bam/bowtie2/bam"

if [[ ! -d "${dir_out}" ]]; then mkdir "${dir_out}"; fi

macs3 filterdup \
    --ifile "${file_IP_bam}" \
    --format "BAMPE" \
    --gsize "${gsize}" \
    --pvalue "1e-5" \
    --keep-dup "auto" \
    --buffer-size 100000 \
    --outdir "${dir_out}" \
    --ofile "${file_IP_bed}" \
    --verbose 3 \
    --dry-run

#  It seems that running the above still results in an empty BEDPE file; need
#+ to report this as a bug

macs3 filterdup \
    --ifile "${file_IP_bam}" \
    --format "BAMPE" \
    --gsize "${gsize}" \
    --pvalue "1e-5" \
    --keep-dup "auto" \
    --buffer-size 100000 \
    --outdir "${dir_out}" \
    --ofile "${file_IP_bed}" \
    --verbose 3

#  It seems that the binomial test is an evaluation of the max number of
#+ duplicates allowed?

macs3 filterdup \
    --ifile "${file_in_bam}" \
    --format "BAMPE" \
    --gsize "${gsize}" \
    --pvalue "1e-5" \
    --keep-dup "auto" \
    --buffer-size 100000 \
    --outdir "${dir_out}" \
    --ofile "${file_in_bed}" \
    --verbose 3
```
</details>
<br />

<a id="printed"></a>
##### Printed
<details>
<summary><i>Printed: My proposed invocation of `macs3 filterdup`</i></summary>

```txt
❯ file_IP_bam="IP_Q_Esa5_7041.sort-coord.bam"


❯ file_IP_bed="IP_Q_Esa5_7041.sort-coord.bed"


❯ file_in_bam="in_${file_IP_bam##IP_}"


❯ file_in_bed="${file_in/.bam/.bed}"


❯ gsize=11624332


❯ dir_out="macs3_test"


❯ cd "${HOME}/projects-etc/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/03_bam/bowtie2/bam"


❯ source activate macs3_env


❯ if [[ ! -d macs3_test ]]; then mkdir macs3_test; fi
mkdir: created directory 'macs3_test'


❯ macs3 filterdup \
>     --ifile "${file_IP_bam}" \
>     --format "BAMPE" \
>     --gsize "${gsize}" \
>     --pvalue "1e-5" \
>     --keep-dup "auto" \
>     --buffer-size 100000 \
>     --outdir "${dir_out}" \
>     --ofile "${file_IP_bed}" \
>     --verbose 3 \
>     --dry-run
INFO  @ 23 Feb 2024 08:56:54: [61 MB] # read input file in Paired-end mode.
INFO  @ 23 Feb 2024 08:56:54: [61 MB] # read treatment fragments...
INFO  @ 23 Feb 2024 08:56:56: [74 MB]  1000000 fragments parsed
INFO  @ 23 Feb 2024 08:56:58: [84 MB]  2000000 fragments parsed
INFO  @ 23 Feb 2024 08:57:01: [91 MB]  3000000 fragments parsed
INFO  @ 23 Feb 2024 08:57:03: [99 MB]  4000000 fragments parsed
INFO  @ 23 Feb 2024 08:57:05: [105 MB]  5000000 fragments parsed
INFO  @ 23 Feb 2024 08:57:08: [110 MB]  6000000 fragments parsed
INFO  @ 23 Feb 2024 08:57:10: [118 MB]  7000000 fragments parsed
INFO  @ 23 Feb 2024 08:57:12: [126 MB]  8000000 fragments parsed
INFO  @ 23 Feb 2024 08:57:15: [135 MB]  9000000 fragments parsed
INFO  @ 23 Feb 2024 08:57:17: [143 MB]  10000000 fragments parsed
INFO  @ 23 Feb 2024 08:57:18: [146 MB] 10385714 fragments have been read.
INFO  @ 23 Feb 2024 08:57:24: [146 MB] # total fragments/pairs in alignment file: 10385714
INFO  @ 23 Feb 2024 08:57:24: [146 MB] calculate max duplicate tags in single position based on binomal distribution...
INFO  @ 23 Feb 2024 08:57:24: [146 MB]  max_dup_tags based on binomal = 7
INFO  @ 23 Feb 2024 08:57:24: [146 MB] filter out redundant tags at the same location and the same strand by allowing at most 7 tag(s)
INFO  @ 23 Feb 2024 08:57:35: [162 MB]  tags after filtering in alignment file: 10382147
INFO  @ 23 Feb 2024 08:57:35: [162 MB]  Redundant rate of alignment file: 0.00
INFO  @ 23 Feb 2024 08:57:35: [162 MB] Dry-run is finished!


❯ macs3 filterdup \
>        --ifile "${file_IP_bam}" \
>        --format "BAMPE" \
>        --gsize "${gsize}" \
>        --pvalue "1e-5" \
>        --keep-dup "auto" \
>        --buffer-size 100000 \
>        --outdir "${dir_out}" \
>        --ofile "${file_IP_bed}" \
>        --verbose 3
INFO  @ 23 Feb 2024 09:01:11: [61 MB] # read input file in Paired-end mode.
INFO  @ 23 Feb 2024 09:01:11: [61 MB] # read treatment fragments...
INFO  @ 23 Feb 2024 09:01:13: [74 MB]  1000000 fragments parsed
INFO  @ 23 Feb 2024 09:01:15: [84 MB]  2000000 fragments parsed
INFO  @ 23 Feb 2024 09:01:18: [91 MB]  3000000 fragments parsed
INFO  @ 23 Feb 2024 09:01:20: [99 MB]  4000000 fragments parsed
INFO  @ 23 Feb 2024 09:01:22: [105 MB]  5000000 fragments parsed
INFO  @ 23 Feb 2024 09:01:24: [110 MB]  6000000 fragments parsed
INFO  @ 23 Feb 2024 09:01:27: [118 MB]  7000000 fragments parsed
INFO  @ 23 Feb 2024 09:01:29: [126 MB]  8000000 fragments parsed
INFO  @ 23 Feb 2024 09:01:31: [135 MB]  9000000 fragments parsed
INFO  @ 23 Feb 2024 09:01:34: [143 MB]  10000000 fragments parsed
INFO  @ 23 Feb 2024 09:01:34: [145 MB] 10385714 fragments have been read.
INFO  @ 23 Feb 2024 09:01:40: [145 MB] # total fragments/pairs in alignment file: 10385714
INFO  @ 23 Feb 2024 09:01:40: [145 MB] calculate max duplicate tags in single position based on binomal distribution...
INFO  @ 23 Feb 2024 09:01:41: [145 MB]  max_dup_tags based on binomal = 7
INFO  @ 23 Feb 2024 09:01:41: [145 MB] filter out redundant tags at the same location and the same strand by allowing at most 7 tag(s)
INFO  @ 23 Feb 2024 09:01:51: [155 MB]  tags after filtering in alignment file: 10382147
INFO  @ 23 Feb 2024 09:01:51: [155 MB]  Redundant rate of alignment file: 0.00
INFO  @ 23 Feb 2024 09:01:51: [155 MB] Write to BED file
INFO  @ 23 Feb 2024 09:02:07: [156 MB] finished! Check IP_Q_Esa5_7041.sort-coord.bed.


❯ macs3 filterdup \
>     --ifile "${file_in_bam}" \
>     --format "BAMPE" \
>     --gsize "${gsize}" \
>     --pvalue "1e-5" \
>     --keep-dup "auto" \
>     --buffer-size 100000 \
>     --outdir "${dir_out}" \
>     --ofile "${file_in_bed}" \
>     --verbose 3
INFO  @ 23 Feb 2024 13:08:37: [61 MB] # read input file in Paired-end mode.
INFO  @ 23 Feb 2024 13:08:37: [61 MB] # read treatment fragments...
INFO  @ 23 Feb 2024 13:08:43: [75 MB]  1000000 fragments parsed
INFO  @ 23 Feb 2024 13:08:49: [78 MB]  2000000 fragments parsed
INFO  @ 23 Feb 2024 13:08:57: [88 MB]  3000000 fragments parsed
INFO  @ 23 Feb 2024 13:09:04: [92 MB]  4000000 fragments parsed
INFO  @ 23 Feb 2024 13:09:13: [101 MB]  5000000 fragments parsed
INFO  @ 23 Feb 2024 13:09:20: [109 MB]  6000000 fragments parsed
INFO  @ 23 Feb 2024 13:09:26: [117 MB]  7000000 fragments parsed
INFO  @ 23 Feb 2024 13:09:34: [124 MB]  8000000 fragments parsed
INFO  @ 23 Feb 2024 13:09:41: [133 MB]  9000000 fragments parsed
INFO  @ 23 Feb 2024 13:09:49: [140 MB]  10000000 fragments parsed
INFO  @ 23 Feb 2024 13:09:57: [149 MB]  11000000 fragments parsed
INFO  @ 23 Feb 2024 13:10:04: [157 MB]  12000000 fragments parsed
INFO  @ 23 Feb 2024 13:10:10: [165 MB]  13000000 fragments parsed
INFO  @ 23 Feb 2024 13:10:15: [172 MB] 13741759 fragments have been read.
INFO  @ 23 Feb 2024 13:10:39: [172 MB] # total fragments/pairs in alignment file: 13741759
INFO  @ 23 Feb 2024 13:10:39: [172 MB] calculate max duplicate tags in single position based on binomal distribution...
INFO  @ 23 Feb 2024 13:10:39: [172 MB]  max_dup_tags based on binomal = 8
INFO  @ 23 Feb 2024 13:10:39: [172 MB] filter out redundant tags at the same location and the same strand by allowing at most 8 tag(s)
INFO  @ 23 Feb 2024 13:11:10: [197 MB]  tags after filtering in alignment file: 13736308
INFO  @ 23 Feb 2024 13:11:10: [197 MB]  Redundant rate of alignment file: 0.00
INFO  @ 23 Feb 2024 13:11:10: [197 MB] Write to BED file
INFO  @ 23 Feb 2024 13:12:04: [198 MB] finished! Check in_Q_Esa5_7041.sort-coord.bed.
```
</details>
<br />

<a id="conversation-with-gpt4"></a>
#### Conversation with GPT4
<a id="me-4"></a>
##### Me
I see in this in the documentation for `macs3 filterdup`:
```txt
It controls the 'macs3 filterdup' behavior towards duplicate tags/pairs at the exact same location -- the same
                        coordination and the same strand. The 'auto' option makes 'macs3 filterdup' calculate the maximum tags at the exact
                        same location based on binomal distribution using given -p as pvalue cutoff; and the 'all' option keeps every tags
                        (useful if you only want to convert formats). If an integer is given, at most this number of tags will be kept at the
                        same location. Note, MACS3 callpeak function uses KEEPDUPLICATES=1 as default. Note, if you've used samtools or
                        picard to flag reads as 'PCR/Optical duplicate' in bit 1024, MACS3 will still read them although the reads may be
                        decided by MACS3 as duplicate later. Default: auto
```

Can you explain to me what is going on with this binomial test? What is this test, why is it applicable in this circumstance, and why are not all duplicates just being removed wholesale, i.e., why is this test being run instead of the wholesale removal of identified duplicates (fragments with same 5' positions)?

<a id="gpt4"></a>
##### GPT4
The `macs3 filterdup` command's use of a binomial test to manage duplicate tags at the exact same location is a sophisticated method for distinguishing between true biological signal and technical duplicates arising from PCR amplification or sequencing. Here's an explanation of the concepts involved:

<a id="what-is-the-binomial-test"></a>
###### What is the binomial test?
The binomial test is a statistical test that determines the probability of observing a given number of "successes" in a fixed number of trials, with a known probability of success on each trial. In the context of ChIP-seq data analysis, a "success" could be considered the occurrence of a tag (read or read pair) at a particular genomic location, and the fixed number of trials could be the total number of reads mapping to that location.

<a id="why-is-it-applicable-in-this-circumstance"></a>
###### Why is it applicable in this circumstance?
In ChIP-seq experiments, duplicates can arise due to the overamplification of certain fragments during PCR, leading to an artificially high number of reads mapping to the same genomic location. However, not all duplicates are artifacts. Some duplicates represent true biological signal, especially in regions of high protein-DNA binding affinity. The binomial test helps distinguish between these scenarios by assessing whether the observed number of duplicates at a specific location is significantly higher than expected by chance, given the overall read depth and the expected distribution of reads across the genome.

<a id="why-not-remove-all-duplicates-wholesale"></a>
###### Why not remove all duplicates wholesale?
Simply removing all duplicates could potentially eliminate valuable biological signals, especially in high-coverage areas where true binding events may generate multiple reads mapping to the same location. By using a statistical test like the binomial distribution, MACS3 can retain a certain number of duplicates that are statistically justified, thus preserving the true signal while reducing the impact of technical artifacts.

This approach reflects a balance between completely removing all duplicates (which might discard true signals) and keeping all duplicates (which might include many artifacts). By setting a threshold based on the binomial distribution, macs3 filterdup aims to optimize the balance between sensitivity (detecting true binding events) and specificity (minimizing the impact of PCR duplicates).

<a id="conclusion"></a>
###### Conclusion
The use of a binomial test in `macs3 filterdup` is a nuanced approach to managing duplicates in ChIP-seq data. It allows researchers to fine-tune their analysis to distinguish between true biological signals and technical artifacts, ensuring that the resulting peaks represent genuine protein-DNA interactions. This statistical method adds a layer of sophistication to ChIP-seq analysis, helping to improve the reliability and interpretability of the results.

<a id="me-5"></a>
##### Me
Can you give me an example of what such a binomial test could like (in R code)? Try to set up and contextualize the test based on what was discussed above (e.g., "The binomial test helps distinguish between these scenarios by assessing whether the observed number of duplicates at a specific location is significantly higher than expected by chance, given the overall read depth and the expected distribution of reads across the genome.") This will help me to better understand what's happening under the hood (so to speak).

<a id="gpt4-1"></a>
##### GPT4
To illustrate how a binomial test could be applied in the context of filtering duplicates in ChIP-seq data, let's set up a hypothetical scenario and then go through the R code to perform the test.

<a id="scenario"></a>
###### Scenario
Imagine you're analyzing ChIP-seq data and have come across a genomic location where 100 reads are mapped. From prior knowledge or based on the overall data, you estimate that, on average, 5% of the genomic locations should have reads due to true biological signals (i.e., the probability of a read being at a specific location due to a true signal is 0.05). You want to test if the observed number of reads (duplicates) at this location significantly exceeds what you would expect by chance, using a significance level of 0.05.

<a id="binomial-test-in-r"></a>
###### Binomial test in R
The binomial test will assess whether the observed number of "successes" (in this case, reads at a specific location, which we're treating as duplicates) significantly deviates from the expected number under a binomial distribution. Here, "success" is the event of a read mapping to the specific location due to a true biological signal, and the number of trials is the total number of reads mapped to that location.

```r
# Parameters
observed_duplicates <- 100  # Observed no. reads at the location
total_reads <- 100          # Total reads mapped to this location
prob_success <- 0.05        # Prob. read at this location b/c/o/true signal

# Perform the binomial test
binom_test_result <- binom.test(
    observed_duplicates, total_reads, prob_success, alternative = "greater"
)

# Print the result
print(binom_test_result)
```

<a id="interpretation"></a>
###### Interpretation
The output of `binom.test()` will give you several pieces of information, including the p-value. The p-value tells you the probability of observing as many or more duplicates as you did if the true rate of duplicates due to biological signals was indeed 5%.
- If the p-value is less than your significance level (e.g., 0.05), you reject the null hypothesis and conclude that the observed number of duplicates is significantly higher than expected by chance alone. This suggests that not all duplicates at this location are due to technical artifacts; some may represent true biological signals.
- If the p-value is greater than your significance level, you fail to reject the null hypothesis, suggesting that the observed number of duplicates could be explained by the expected distribution of reads across the genome, and they might primarily be technical duplicates.

<a id="note-2"></a>
###### Note
This example simplifies the actual complexities involved in ChIP-seq data analysis and assumes all reads at a location are duplicates, which may not always be the case. The real scenario might involve distinguishing between different types of duplicates and considering other biological and technical factors.

<a id="me-6"></a>
##### Me
Can you explicitly clarify what null and alternative hypotheses are for the above-mentioned example binomial test?

<a id="gpt4-2"></a>
##### GPT4
In the context of the binomial test example provided for analyzing ChIP-seq data duplicates, the null and alternative hypotheses can be defined as follows:

<a id="null-hypothesis-h0"></a>
###### Null hypothesis (H<sub>0</sub>)
The null hypothesis states that the observed number of reads (duplicates) at a specific genomic location is consistent with the expected number based on a binomial distribution, with the success probability equal to the probability of a read mapping to that location due to a true biological signal. In simpler terms, it suggests that the observed duplicates are no more frequent than would be expected by chance alone, given the overall read depth and the probability of a read being at that location due to a true biological signal. For the given example, this means assuming the true rate of reads mapping to the location due to biological signals is 5%.

<b>H<sub>0</sub></b>: The number of observed duplicates at this location is not significantly higher than expected by chance, based on a binomial distribution with a success probability of 0.05.

<a id="alternative-hypothesis-h1"></a>
###### Alternative Hypothesis (H<sub>1</sub>)
The alternative hypothesis posits that the observed number of reads (duplicates) at the specific location is significantly greater than what would be expected based on a binomial distribution with the given success probability. This would imply that the observed number of duplicates cannot be solely attributed to chance and may indicate an enrichment of reads beyond what is expected from the background read distribution, suggesting potential biological significance.

<b>H<sub>1</sub></b>: The number of observed duplicates at this location is significantly higher than expected by chance, indicating potential biological significance beyond mere technical duplication.

The binomial test then evaluates these hypotheses by calculating the probability (p-value) of observing the given number of duplicates (or more) under the null hypothesis. If this probability is lower than the predetermined significance level (e.g., α = 0.05), the null hypothesis is rejected in favor of the alternative, suggesting that the observed duplication is unlikely to be solely due to chance.

<a id="additional-conversation"></a>
##### Additional conversation
<a id="me-7"></a>
<b><i>Me</i></b><br />
Can you teach me what a binomial distribution is? Please contextualize your lesson with the R code example (etc.) above.

<a id="gpt4-3"></a>
<b><i>GPT4</i></b><br />
The binomial distribution is used in statistics to describe the probability of obtaining a specific number of successes in a fixed number of independent trials, each with the same probability of success. In the context of your question about `macs3 filterdup`, a binomial test could assess whether the number of duplicates at a location significantly exceeds what would be expected by chance, given the total read depth and uniform distribution assumption across the genome. This method allows for distinguishing genuine biological enrichment from technical duplicates.

<a id="me-7"></a>
<b><i>Me</i></b><br />
So, in this particular example, the number of duplicates at a given location is the "successes", and the "fixed number of independent trials" is some value obtained from calculating "what would be expected by chance, given the total read depth and uniform distribution assumption across the genome."

If this is correct, can you give me some example R code demonstrating the above?

<a id="gpt4-4"></a>
<b><i>GPT4</i></b><br />
In R, you can perform a binomial test to evaluate if the number of duplicates at a specific genomic location significantly exceeds what's expected by chance. This involves using the `binom.test` function, where the number of successes is the observed number of duplicates, and the number of trials is based on the expected distribution, considering the total read depth and genome uniformity. The success probability in each trial (p) could be the expected fraction of reads at any location if reads were evenly distributed. Here's a simplified example:
```r
#  Observed number of duplicates at a location
observed_duplicates <- 10

#  Total number of reads (trials)
total_reads <- 100000

#  Assuming uniform distribution, the probability of any read mapping to this
#+ specific location
p_success <- 1 / genome_size

#  Perform binomial test
binom.test(observed_duplicates, total_reads, p_success)
```

This test helps determine if the observed duplication at a location is statistically significant, considering the genome-wide read distribution.

<a id="me-8"></a>
<b><i>Me</i></b><br />
In example R code, how might the expected distribution be calculated?

<a id="gpt4-5"></a>
<b><i>GPT4</i></b><br />
The expected distribution of reads across the genome can be calculated by assuming a uniform distribution of reads. This means if you have a total number of reads obtained from sequencing and a known genome size, you can calculate the expected number of reads per genomic location (base pair, for instance) by dividing the total number of reads by the genome size. This gives you an average read density or the probability of a read mapping to any specific location in the genome, assuming reads are evenly distributed.

<a id="me-9"></a>
<b><i>Me</i></b><br />
Let's say the total estimated number of unique 50-mers in my genome  of interest, *S. cerevisiae*, is 11,624,332, and 80,000,000 paired-end ChIP-seq reads were sequenced. What would example R code to calculate the expected distribution look like?

Furthermore, are we only concerned with the mean of the expected distribution? We don't also try to calculate other summary statistics?

<a id="gpt4-6"></a>
<b><i>GPT4</i></b><br />
Given your scenario, we can use R to calculate the expected number of reads per unique 50-mer (genomic location) under the assumption of a uniform distribution. The mean of the expected distribution is the primary statistic of interest here because it represents the average expected coverage across the genome, which is crucial for assessing whether the observed duplicates at a specific location significantly exceed this expectation.

Here's how you might calculate it:
```r
total_unique_50mers <- 11624332
total_reads <- 80000000  # Assuming each read represents a 50-mer

#  Calculate the expected number of reads per 50-mer
expected_reads_per_50mer <- total_reads / total_unique_50mers

print(expected_reads_per_50mer)
```

This code calculates the average (mean) number of reads per unique 50-mer across the genome. In this context, other summary statistics like variance or standard deviation are less relevant because we're primarily interested in comparing the observed number of reads at specific locations to this mean expectation to assess duplication.

<a id="me-10"></a>
<b><i>Me</i></b><br />
Please connect the example to calculate `expected_reads_per_50mer` with the above code examples you described,
```r
# Parameters
observed_duplicates <- 100  # Observed no. reads at the location
total_reads <- 100          # Total reads mapped to this location
prob_success <- 0.05        # Prob. read at this location b/c/o/true signal

# Perform the binomial test
binom_test_result <- binom.test(
    observed_duplicates, total_reads, prob_success, alternative = "greater"
)

# Print the result
print(binom_test_result)
```
and
```r
#  Observed number of duplicates at a location
observed_duplicates <- 10

#  Total number of reads (trials)
total_reads <- 100000

#  Assuming uniform distribution, the probability of any read mapping to this
#+ specific location
p_success <- 1 / genome_size

#  Perform binomial test
binom.test(observed_duplicates, total_reads, p_success)
```

I think that doing so will help me to get a bigger picture of what is going on with duplicate-handling by MACS3.

<a id="gpt4-7"></a>
<b><i>GPT4</i></b><br />
The calculation of `expected_reads_per_50mer` helps us understand the baseline or average number of reads we might expect to map to any given unique 50-mer in the genome, assuming a uniform distribution of reads. This average serves as a reference point when we perform a binomial test to assess the significance of observing a specific number of duplicates at a particular genomic location.

In the context of duplicate handling by MACS3, this understanding aids in distinguishing between true biological signal and potential technical artifacts. The binomial test evaluates whether the observed number of duplicates at a location significantly exceeds what we would expect based on the average distribution (i.e., `expected_reads_per_50mer`). If the probability of observing this many or more duplicates by chance (given the total reads and the expected probability of any read mapping to this specific location) is low, it suggests that the observed duplicates may not be solely due to random sequencing depth but could be indicative of true enrichment or technical biases in the library preparation.

Thus, by setting `prob_success` in the binomial test to the inverse of the total unique 50-mers (or using the calculated `expected_reads_per_50mer` to inform our expectation), we directly connect the expected uniform distribution across the genome with the statistical assessment of observed duplicates, providing a more nuanced approach to identifying significant enrichments or technical artifacts in ChIP-seq data.

<a id="me-11"></a>
<b><i>Me</i></b><br />
Please show it all to me with example R code.

<a id="gpt4-8"></a>
<b><i>GPT4</i></b><br />
To connect the calculation of expected_reads_per_50mer with the binomial test examples provided, consider a scenario where the total number of unique 50-mers in the *S. cerevisiae* genome is 11,624,332, and 80,000,000 paired-end ChIP-seq reads were sequenced. Here's how the overall process might look in R code:
```r
#  Parameters
genome_size <- 11624332  # Total estimated unique 50-mers in S. cerevisiae
total_reads <- 80000000  # Total paired-end ChIP-seq reads (50-bp) sequenced

#  Calculate the expected reads per unique 50-mer assuming uniform distribution
expected_reads_per_50mer <- total_reads / genome_size

#  For a specific location, assuming observed duplicates and total reads mapped
observed_duplicates <- 100
total_reads_at_location <- 1000

#  Assuming uniform distribution, calculate the probability of any read mapping
#+ to this specific location
prob_success <- 1 / genome_size

#  Perform the binomial test to assess if the observed number of duplicates at
#+ the location is significantly higher than expected by chance
binom_test_result <- binom.test(
    observed_duplicates, total_reads_at_location, prob_success,
    alternative = "greater"
)

#  Print the results
print(binom_test_result)
```

This example integrates the concept of `expected_reads_per_50mer` into the assessment of whether the observed number of duplicates at a specific location is unusually high given the total read depth and the assumed uniform distribution of reads across the genome. The binomial test's null hypothesis is that the observed number of duplicates comes from a distribution with success probability equal to `prob_success`, and the alternative hypothesis is that the true probability of success is greater, indicating potential enrichment or technical artifacts.

<a id="me-12"></a>
<b><i>Me</i></b><br />
Actually, it does not look like `expected_reads_per_50mer` is integrated into the binomial test. It is initialized/calculated, but not referenced in subsequent code in any way.

<a id="gpt4-9"></a>
<b><i>GPT4</i></b><br />
You're correct, and I appreciate your patience. The `expected_reads_per_50mer` was calculated to give us an idea of the expected number of reads per unique 50-mer in the genome, which serves as a basis for understanding how reads might be distributed across the genome if they were uniformly distributed. However, it wasn't directly used in the binomial test examples provided.

The correct approach to integrate `expected_reads_per_50mer` into the binomial test would involve using it to calculate the probability of success (`p_success`) in the binomial distribution model. Essentially, `p_success` represents the probability that a given read will map to a specific 50-mer location if reads are uniformly distributed across the genome. Here's an example of how you might calculate this and use it in a binomial test:
```r
#  Total estimated unique 50-mers in the yeast genome
unique_50mers <- 11624332

#  Total reads sequenced
total_reads_sequenced <- 80000000

#  Calculate the expected reads per 50-mer
expected_reads_per_50mer <- total_reads_sequenced / unique_50mers

#  Assuming uniform distribution, calculate the probability of a read mapping
#+ to a specific 50-mer
p_success <- 1 / unique_50mers

#  Observed number of duplicates at a specific 50-mer location
observed_duplicates <- 100

#  Perform the binomial test using the observed duplicates and the calculated
#+ probability of success
binom_test_result <- binom.test(
    observed_duplicates, expected_reads_per_50mer, p_success,
    alternative = "greater"
)

#  Print the result
print(binom_test_result)
```

In this revised example, `expected_reads_per_50mer` is used to provide the number of trials in the binomial test (`binom.test` function in R), which is more aligned with the actual setup of the problem. This approach better reflects the expected behavior of the sequencing experiment given the genome complexity and read distribution, making the binomial test's application and interpretation more accurate.

<a id="me-13"></a>
<b><i>Me</i></b><br />
Thank you. This makes a bit more sense, particularly
```r
binom_test_result <- binom.test(
    observed_duplicates, expected_reads_per_50mer, ...
)
```

However, I think it's reasonable for the user to supply a fixed alpha threshold to `p_success` (e.g., 0.00001) rather than calculate it as you have (which results in the value 8.60264487e-8, which in turn strikes me as perhaps too stringent).

<a id="gpt4-10"></a>
<b><i>GPT4</i></b><br />
Setting a fixed alpha threshold for `p_success` in the binomial test instead of calculating it based on the uniform distribution assumption is indeed a practical approach, especially in contexts where the uniform distribution assumption may not hold true due to various factors like sequencing biases or non-uniform genomic features. Using a fixed alpha threshold allows for a more flexible interpretation of duplication significance, accommodating the variability and complexity of genomic data. This approach can help in differentiating between significant enrichment of duplicates due to biological reasons versus technical artifacts.

<a id="me-14"></a>
<b><i>Me</i></b><br />
Speaking of which, why is a uniform distribution assumed given that this is an analysis of ChIP-seq data, which has the expectation that only a subset of the whole genome—the protein binding sites—will be enriched by the experiment?

Also, is a uniform distribution indeed assumed by `macs3 filterdup`?

<a id="gpt4-11"></a>
<b><i>GPT4</i></b><br />
The assumption of a uniform distribution for ChIP-seq data analysis, including in tools like `macs3 filterdup`, is a simplification that facilitates computational analysis. While ChIP-seq data indeed targets specific genomic regions (protein-binding sites), assuming a uniform distribution allows for a baseline expectation against which significant deviations (enrichments) can be detected. This assumption helps in identifying regions with unexpectedly high read counts, which may indicate biological significance. However, it's crucial to consider that in practice, ChIP-seq data will not be uniformly distributed due to biological variances and technical biases.

<a id="me-15"></a>
<b><i>Me</i></b><br />
In looking at the following printed output from calling `macs3 filterdups`,
```txt
❯ macs3 filterdup \
>        --ifile "${file_IP_bam}" \
>        --format "BAMPE" \
>        --gsize "${gsize}" \
>        --pvalue "1e-5" \
>        --keep-dup "auto" \
>        --buffer-size 100000 \
>        --outdir "${dir_out}" \
>        --ofile "${file_IP_bed}" \
>        --verbose 3
INFO  @ 23 Feb 2024 09:01:11: [61 MB] # read input file in Paired-end mode.
INFO  @ 23 Feb 2024 09:01:11: [61 MB] # read treatment fragments...
INFO  @ 23 Feb 2024 09:01:13: [74 MB]  1000000 fragments parsed
INFO  @ 23 Feb 2024 09:01:15: [84 MB]  2000000 fragments parsed
INFO  @ 23 Feb 2024 09:01:18: [91 MB]  3000000 fragments parsed
INFO  @ 23 Feb 2024 09:01:20: [99 MB]  4000000 fragments parsed
INFO  @ 23 Feb 2024 09:01:22: [105 MB]  5000000 fragments parsed
INFO  @ 23 Feb 2024 09:01:24: [110 MB]  6000000 fragments parsed
INFO  @ 23 Feb 2024 09:01:27: [118 MB]  7000000 fragments parsed
INFO  @ 23 Feb 2024 09:01:29: [126 MB]  8000000 fragments parsed
INFO  @ 23 Feb 2024 09:01:31: [135 MB]  9000000 fragments parsed
INFO  @ 23 Feb 2024 09:01:34: [143 MB]  10000000 fragments parsed
INFO  @ 23 Feb 2024 09:01:34: [145 MB] 10385714 fragments have been read.
INFO  @ 23 Feb 2024 09:01:40: [145 MB] # total fragments/pairs in alignment file: 10385714
INFO  @ 23 Feb 2024 09:01:40: [145 MB] calculate max duplicate tags in single position based on binomal distribution...
INFO  @ 23 Feb 2024 09:01:41: [145 MB]  max_dup_tags based on binomal = 7
INFO  @ 23 Feb 2024 09:01:41: [145 MB] filter out redundant tags at the same location and the same strand by allowing at most 7 tag(s)
INFO  @ 23 Feb 2024 09:01:51: [155 MB]  tags after filtering in alignment file: 10382147
INFO  @ 23 Feb 2024 09:01:51: [155 MB]  Redundant rate of alignment file: 0.00
INFO  @ 23 Feb 2024 09:01:51: [155 MB] Write to BED file
INFO  @ 23 Feb 2024 09:02:07: [156 MB] finished! Check IP_Q_Esa5_7041.sort-coord.bed.
```
it seems that the binomial test is an evaluation of the max number of duplicates allowed? Or is it more like what we discussed and the program has just done some kind of after-the-fact test to see what the maximum number of duplicate tags at a location is after running the (multiple) binomial tests?

<a id="gpt4-12"></a>
<b><i>GPT4</i></b><br />
The `macs3 filterdup` command uses a binomial distribution to decide the maximum number of duplicate tags (reads or fragments) at the same location that can be considered as not being due to PCR duplication artifacts. The binomial distribution models the probability of observing a given number of duplicates purely by chance, given the total number of reads and the assumption of a uniform distribution of reads across the genome. The `auto` option lets MACS3 determine this maximum number dynamically based on a specified p-value (probability threshold), which helps differentiate between genuine biological signals and technical duplicates. This approach allows for a more nuanced handling of duplicates, recognizing that not all duplicates are artifacts, especially in high-throughput sequencing data where some level of duplication is expected due to the nature of the biological samples and the sequencing process.

<a id="me-16"></a>
<b><i>Me</i></b><br />
I guess my assumption was that the test is being run on a locus-by-locus-by-locus-etc. basis to decide whether given duplicates at a given locus can remain purely on the basis of the given locus?

<a id="gpt4-3"></a>
<b><i>GPT4</i></b><br />
Yes, your understanding is correct. The binomial test in `macs3 filterdup` is applied locus by locus to decide whether the observed number of duplicates at each specific location exceeds what would be expected by chance, based on the total read depth and a uniform distribution assumption across the genome. This method allows MACS3 to discern between likely technical duplicates (PCR artifacts) and genuine biological duplicates (true signal enrichment), enabling a more refined approach to handling duplicates in ChIP-seq data analysis.

<a id="me-17"></a>
<b><i>Me</i></b><br />
Then where does the cutoff of 7 duplicates enter into the picture?

<a id="gpt4-4"></a>
<b><i>GPT4</i></b><br />
The cutoff of 7 duplicates determined by `macs3 filterdup` using a binomial distribution test represents the maximum number of duplicate reads allowed at a single genomic location before being considered redundant. This threshold is calculated based on the total number of reads, the genome size, and a specified p-value, aiming to balance between retaining genuine biological signal and removing technical duplicates. The process involves statistically estimating the likelihood of observing a given number of duplicates purely by chance, thus setting a data-driven threshold for filtering out excessive duplicates.

<a id="me-18"></a>
<b><i>Me</i></b><br />
So it's not on a per-locus basis&mdash;it's on an all-loci basis then?

<a id="gpt4-5"></a>
<b><i>GPT4</i></b><br />
Yes, the threshold of 7 duplicates established by macs3 filterdup is applied across all loci in the dataset, rather than being calculated for each locus individually. This global threshold is statistically derived to maximize the retention of genuine biological signals while minimizing the impact of technical duplicates. By setting a uniform threshold based on the entire dataset, macs3 filterdup ensures a consistent approach to filtering duplicates across the genome, rather than varying the threshold from one genomic location to another.

<a id="me-19"></a>
<b><i>Me</i></b><br />
How?

<a id="gpt4-3"></a>
<b><i>GPT4</i></b><br />
The threshold of 7 duplicates is determined by a binomial test that evaluates the entire dataset's characteristics, such as the total number of reads and the genome size, to statistically estimate the maximum number of duplicates that could be present at any location by chance. This approach aims to identify a balance where genuine biological signals are preserved while minimizing the impact of technical artifacts. By applying a uniform threshold across all loci, `macs3 filterdup` ensures consistency in filtering, which is crucial for maintaining the integrity of the dataset and the reliability of downstream analyses.

<a id="step-2-decide-the-fragment-length-d"></a>
### Step 2: Decide the fragment length `d`
This is an important step for MACS3 to analyze ChIP-seq and also for other types of data since the location of sequenced read may only tell you the end of a DNA fragment that you are interested in (such as TFBS or DNA hypersensitive regions), and you have to estimate how long this DNA fragment is in order to recover the actual enrichment. You can also regard this as a data smoothing technic. You know that `macs3 callpeak` will output something like, it can identify certain number of pairs of peaks and it can predict the fragment length, or `d` in MACS3 terminology, using cross-correlation. In fact, you can also do this using `predictd` module. Normally, we only need to do this for ChIP data:
```bash
#!/bin/bash

macs3 predictd \
    -i CTCF_ChIP_200K_filterdup.bed \
    -g hs \
    -m 5 50
```

Here the `-g` (the genome size) need to be set according to your sample, and the `mfold` parameters have to be set reasonably. To simulate the default behavior of macs3 callpeak, set `-m 5 50`. Of course, you can tweak it. The output from `predictd` will tell you the fragment length `d`, and in this example, it is 254. Write the number down on your notebook since we will need it in the next step. Of course, if you do not want to extend the reads or you have a better estimation on fragment length, you can simply skip this step.

<a id="documentation-for-macs-predictd"></a>
#### Documentation for `macs predictd`
<details>
<summary><i>Text: Documentation for `macs predictd`</i></summary>

```txt
❯ macs3 predictd -h
usage: macs3 predictd [-h] -i IFILE [IFILE ...] [-f {AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE,BEDPE}] [-g GSIZE] [-s TSIZE] [--bw BW] [--d-min D_MIN] [-m MFOLD MFOLD] [--outdir OUTDIR] [--rfile RFILE] [--buffer-size BUFFER_SIZE] [--verbose VERBOSE]

options:
  -h, --help            show this help message and exit
  -i IFILE [IFILE ...], --ifile IFILE [IFILE ...]
                        ChIP-seq alignment file. If multiple files are given as '-t A B C', then they will all be read and combined. REQUIRED.
  -f {AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE,BEDPE}, --format {AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE,BEDPE}
                        Format of tag file, "AUTO", "BED" or "ELAND" or "ELANDMULTI" or "ELANDEXPORT" or "SAM" or "BAM" or "BOWTIE" or "BAMPE" or "BEDPE". The default AUTO option will let MACS decide which format the file is. However, if you want to decide the average insertion size/fragment
                        size from PE data such as BEDPE or BAMPE, please specify the format as BAMPE or BEDPE since MACS3 won't automatically recognize three two formats with -f AUTO. Please be aware that in PE mode, -g, -s, --bw, --d-min, -m, and --rfile have NO effect. DEFAULT: "AUTO"
  -g GSIZE, --gsize GSIZE
                        Effective genome size. It can be 1.0e+9 or 1000000000, or shortcuts:'hs' for human (2,913,022,398), 'mm' for mouse (2,652,783,500), 'ce' for C. elegans (100,286,401) and 'dm' for fruitfly (142,573,017), Default:hs. The effective genome size numbers for the above four
                        species are collected from Deeptools https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html Please refer to deeptools to define the best genome size you plan to use.
  -s TSIZE, --tsize TSIZE
                        Tag size. This will override the auto detected tag size. DEFAULT: Not set
  --bw BW               Band width for picking regions to compute fragment size. This value is only used while building the shifting model. DEFAULT: 300
  --d-min D_MIN         Minimum fragment size in basepair. Any predicted fragment size less than this will be excluded. DEFAULT: 20
  -m MFOLD MFOLD, --mfold MFOLD MFOLD
                        Select the regions within MFOLD range of high-confidence enrichment ratio against background to build model. Fold-enrichment in regions must be lower than upper limit, and higher than the lower limit. Use as "-m 10 30". DEFAULT:5 50
  --outdir OUTDIR       If specified all output files will be written to that directory. Default: the current working directory
  --rfile RFILE         PREFIX of filename of R script for drawing X-correlation figure. DEFAULT:'predictd_model.R' and R file will be predicted_model.R
  --buffer-size BUFFER_SIZE
                        Buffer size for incrementally increasing internal array size to store reads alignment information. In most cases, you don't have to change this parameter. However, if there are large number of chromosomes/contigs/scaffolds in your alignment, it's recommended to
                        specify a smaller buffer size in order to decrease memory usage (but it will take longer time to read alignment files). Minimum memory requested for reading an alignment file is about # of CHROMOSOME * BUFFER_SIZE * 8 Bytes. DEFAULT: 100000
  --verbose VERBOSE     Set verbose level of runtime message. 0: only show critical message, 1: show additional warning message, 2: show process information, 3: show debug messages. DEFAULT:2
```
</details>
<br />

<details>
<summary><i>Code: </i></summary>

```bash
#!/bin/bash

#  If you want to pile up control data as single-end like callpeak, then the
#+ ${ctrl_prefix}_filterdup.bed file has to be broken down into single end BED
#+ file...

macs3 predictd \
    --ifile "${infile}" \
    --format "BED" \
    --gsize "${gsize}" \
    --outdir "${dir_out}" \
    --verbose 3

```
</details>
<br />

<a id="my-proposed-invocation"></a>
#### My proposed invocation
<details>
<summary><i>Code: </i></summary>

```bash
#!/bin/bash

#  Rough-draft test code for IP, input
file_IP_bam="IP_Q_Esa5_7041.sort-coord.bam"
file_IP_bed="IP_Q_Esa5_7041.sort-coord.bed"
file_in_bam="in_${file_IP_bam##IP_}"
file_in_bed="${file_in_bam/.bam/.bed}"
file_in_bed_se="${file_in_bed/.bed/.se.bed}"
gsize=11624332
kmer=50
dir_out="macs3_test"

source activate macs3_env

cd "${HOME}/projects-etc/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/03_bam/bowtie2/bam"

if [[ ! -d "${dir_out}" ]]; then mkdir "${dir_out}"; fi

if [[ ! -f "${dir_out}/${file_IP_bed}" ]]; then
    macs3 filterdup \
        --ifile "${file_IP_bam}" \
        --format "BAMPE" \
        --gsize "${gsize}" \
        --pvalue "1e-5" \
        --keep-dup "auto" \
        --buffer-size 100000 \
        --outdir "${dir_out}" \
        --ofile "${file_IP_bed}" \
        --verbose 3
fi

if [[ ! -f "${dir_out}/${file_in_bed}" ]]; then
    macs3 filterdup \
        --ifile "${file_in_bam}" \
        --format "BAMPE" \
        --gsize "${gsize}" \
        --pvalue "1e-5" \
        --keep-dup "auto" \
        --buffer-size 100000 \
        --outdir "${dir_out}" \
        --ofile "${file_in_bed}" \
        --verbose 3
fi

if [[ ! -f "${dir_out}/${file_in_bed_se}" ]]; then
    cat "${dir_out}/${file_in_bed}" \
        | awk \
            -v OFS="\t" \
            -v kmer="${kmer}" \
            '{ print $1,$2,$2+kmer,".",".","+\n"$1,$3-kmer,$3,".",".","-" }' \
                > "${dir_out}/${file_in_bed_se}"
fi

check_bed_se=true
if ${check_bed_se}; then
    less "${dir_out}/${file_in_bed_se}"
    # q
fi

macs3 predictd \
    --ifile "${dir_out}/${file_in_bed_se}" \
    --format "BED" \
    --gsize "${gsize}" \
    --outdir "${dir_out}" \
    --verbose 3

macs3 callpeak \
    --treatment "${file_IP_bam}" \
    --control "${file_in_bam}" \
    --format "BAMPE" \
    --gsize "${gsize}" \
    --keep-dup "auto" \
    --outdir "${dir_out}" \
    --bdg \
    --SPMR \
    --verbose 3
```
</details>
<br />

<details>
<summary><i>Printed: </i></summary>

https://github.com/macs3-project/MACS/issues/353

```txt
❯ macs3 predictd \
>     --ifile "${dir_out}/${file_in_bed_se}" \
>     --format "BED" \
>     --gsize "${gsize}" \
>     --outdir "${dir_out}" \
>     --verbose 3
INFO  @ 24 Feb 2024 11:05:37: [61 MB] # read alignment files...
INFO  @ 24 Feb 2024 11:05:37: [61 MB] # read treatment tags...
INFO  @ 24 Feb 2024 11:05:38: [119 MB]  1000000 reads parsed
INFO  @ 24 Feb 2024 11:05:39: [126 MB]  2000000 reads parsed
INFO  @ 24 Feb 2024 11:05:39: [136 MB]  3000000 reads parsed
INFO  @ 24 Feb 2024 11:05:40: [136 MB]  4000000 reads parsed
INFO  @ 24 Feb 2024 11:05:40: [138 MB]  5000000 reads parsed
INFO  @ 24 Feb 2024 11:05:41: [148 MB]  6000000 reads parsed
INFO  @ 24 Feb 2024 11:05:41: [158 MB]  7000000 reads parsed
INFO  @ 24 Feb 2024 11:05:42: [160 MB]  8000000 reads parsed
INFO  @ 24 Feb 2024 11:05:43: [167 MB]  9000000 reads parsed
INFO  @ 24 Feb 2024 11:05:43: [169 MB]  10000000 reads parsed
INFO  @ 24 Feb 2024 11:05:44: [169 MB]  11000000 reads parsed
INFO  @ 24 Feb 2024 11:05:44: [178 MB]  12000000 reads parsed
INFO  @ 24 Feb 2024 11:05:45: [178 MB]  13000000 reads parsed
INFO  @ 24 Feb 2024 11:05:46: [183 MB]  14000000 reads parsed
INFO  @ 24 Feb 2024 11:05:46: [191 MB]  15000000 reads parsed
INFO  @ 24 Feb 2024 11:05:47: [191 MB]  16000000 reads parsed
INFO  @ 24 Feb 2024 11:05:47: [206 MB]  17000000 reads parsed
INFO  @ 24 Feb 2024 11:05:48: [212 MB]  18000000 reads parsed
INFO  @ 24 Feb 2024 11:05:48: [212 MB]  19000000 reads parsed
INFO  @ 24 Feb 2024 11:05:49: [214 MB]  20000000 reads parsed
INFO  @ 24 Feb 2024 11:05:50: [214 MB]  21000000 reads parsed
INFO  @ 24 Feb 2024 11:05:50: [218 MB]  22000000 reads parsed
INFO  @ 24 Feb 2024 11:05:51: [219 MB]  23000000 reads parsed
INFO  @ 24 Feb 2024 11:05:51: [227 MB]  24000000 reads parsed
INFO  @ 24 Feb 2024 11:05:52: [232 MB]  25000000 reads parsed
INFO  @ 24 Feb 2024 11:05:52: [233 MB]  26000000 reads parsed
INFO  @ 24 Feb 2024 11:05:53: [237 MB]  27000000 reads parsed
INFO  @ 24 Feb 2024 11:05:53: [237 MB] tag size is determined as 50 bps
INFO  @ 24 Feb 2024 11:05:53: [237 MB] # tag size = 50
INFO  @ 24 Feb 2024 11:05:53: [237 MB] # total tags in alignment file: 27472616
INFO  @ 24 Feb 2024 11:05:53: [237 MB] # Build Peak Model...
DEBUG @ 24 Feb 2024 11:05:53: [237 MB] #2 min_tags: 3545; max_tags:35451;
INFO  @ 24 Feb 2024 11:05:53: [237 MB] #2 looking for paired plus/minus strand peaks...
DEBUG @ 24 Feb 2024 11:05:53: [237 MB] Chromosome: b'I'
DEBUG @ 24 Feb 2024 11:05:53: [237 MB] Number of unique tags on + strand: 185204
DEBUG @ 24 Feb 2024 11:05:53: [237 MB] Number of peaks in + strand: 0
DEBUG @ 24 Feb 2024 11:05:53: [237 MB] Number of unique tags on - strand: 185204
DEBUG @ 24 Feb 2024 11:05:53: [237 MB] Number of peaks in - strand: 0
DEBUG @ 24 Feb 2024 11:05:53: [237 MB] Chrom b'I' is discarded!
DEBUG @ 24 Feb 2024 11:05:53: [237 MB] Chromosome: b'II'
DEBUG @ 24 Feb 2024 11:05:53: [237 MB] Number of unique tags on + strand: 671328
DEBUG @ 24 Feb 2024 11:05:53: [237 MB] Number of peaks in + strand: 0
DEBUG @ 24 Feb 2024 11:05:53: [237 MB] Number of unique tags on - strand: 671328
DEBUG @ 24 Feb 2024 11:05:53: [237 MB] Number of peaks in - strand: 0
DEBUG @ 24 Feb 2024 11:05:53: [237 MB] Chrom b'II' is discarded!
DEBUG @ 24 Feb 2024 11:05:53: [237 MB] Chromosome: b'III'
DEBUG @ 24 Feb 2024 11:05:53: [237 MB] Number of unique tags on + strand: 261050
DEBUG @ 24 Feb 2024 11:05:53: [237 MB] Number of peaks in + strand: 1
DEBUG @ 24 Feb 2024 11:05:53: [237 MB] plus peaks: first - (299, 3918.0) ... last - (299, 3918.0)
DEBUG @ 24 Feb 2024 11:05:53: [237 MB] Number of unique tags on - strand: 261050
DEBUG @ 24 Feb 2024 11:05:53: [237 MB] Number of peaks in - strand: 1
DEBUG @ 24 Feb 2024 11:05:53: [237 MB] minus peaks: first - (408, 3831.0) ... last - (408, 3831.0)
DEBUG @ 24 Feb 2024 11:05:53: [237 MB] ip_max: 1; im_max: 1
DEBUG @ 24 Feb 2024 11:05:53: [237 MB] Paired centers: first - 353 ... second - 353
DEBUG @ 24 Feb 2024 11:05:53: [237 MB] Number of paired peaks in this chromosome: 1
DEBUG @ 24 Feb 2024 11:05:53: [237 MB] Chromosome: b'IV'
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of unique tags on + strand: 1280955
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of peaks in + strand: 0
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of unique tags on - strand: 1280955
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of peaks in - strand: 0
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Chrom b'IV' is discarded!
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Chromosome: b'IX'
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of unique tags on + strand: 370116
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of peaks in + strand: 0
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of unique tags on - strand: 370116
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of peaks in - strand: 0
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Chrom b'IX' is discarded!
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Chromosome: b'Mito'
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of unique tags on + strand: 1689248
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of peaks in + strand: 14
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] plus peaks: first - (2559, 15966.0) ... last - (79479, 25502.0)
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of unique tags on - strand: 1689248
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of peaks in - strand: 14
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] minus peaks: first - (2865, 16193.0) ... last - (79751, 25583.0)
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] ip_max: 14; im_max: 14
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Paired centers: first - 2712 ... second - 79615
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of paired peaks in this chromosome: 14
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Chromosome: b'SP_I'
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of unique tags on + strand: 23541
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of peaks in + strand: 0
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of unique tags on - strand: 23541
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of peaks in - strand: 0
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Chrom b'SP_I' is discarded!
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Chromosome: b'SP_II'
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of unique tags on + strand: 19214
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of peaks in + strand: 0
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of unique tags on - strand: 19214
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of peaks in - strand: 0
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Chrom b'SP_II' is discarded!
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Chromosome: b'SP_III'
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of unique tags on + strand: 24548
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of peaks in + strand: 0
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of unique tags on - strand: 24548
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of peaks in - strand: 0
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Chrom b'SP_III' is discarded!
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Chromosome: b'SP_II_TG'
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of unique tags on + strand: 64
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of peaks in + strand: 0
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of unique tags on - strand: 64
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of peaks in - strand: 0
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Chrom b'SP_II_TG' is discarded!
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Chromosome: b'SP_MTR'
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of unique tags on + strand: 123
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of peaks in + strand: 0
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of unique tags on - strand: 123
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of peaks in - strand: 0
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Chrom b'SP_MTR' is discarded!
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Chromosome: b'SP_Mito'
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of unique tags on + strand: 7808
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of peaks in + strand: 0
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of unique tags on - strand: 7808
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of peaks in - strand: 0
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Chrom b'SP_Mito' is discarded!
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Chromosome: b'V'
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of unique tags on + strand: 483332
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of peaks in + strand: 0
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of unique tags on - strand: 483332
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of peaks in - strand: 0
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Chrom b'V' is discarded!
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Chromosome: b'VI'
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of unique tags on + strand: 226515
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of peaks in + strand: 0
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of unique tags on - strand: 226515
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of peaks in - strand: 0
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Chrom b'VI' is discarded!
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Chromosome: b'VII'
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of unique tags on + strand: 901719
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of peaks in + strand: 0
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of unique tags on - strand: 901719
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of peaks in - strand: 0
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Chrom b'VII' is discarded!
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Chromosome: b'VIII'
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of unique tags on + strand: 490159
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of peaks in + strand: 2
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] plus peaks: first - (214099, 4754.0) ... last - (215362, 3679.0)
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of unique tags on - strand: 490159
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of peaks in - strand: 2
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] minus peaks: first - (214247, 4689.0) ... last - (215697, 3662.0)
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] ip_max: 2; im_max: 2
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Paired centers: first - 214173 ... second - 215529
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of paired peaks in this chromosome: 2
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Chromosome: b'X'
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of unique tags on + strand: 621946
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of peaks in + strand: 0
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of unique tags on - strand: 621946
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Number of peaks in - strand: 0
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Chrom b'X' is discarded!
DEBUG @ 24 Feb 2024 11:05:54: [237 MB] Chromosome: b'XI'
DEBUG @ 24 Feb 2024 11:05:55: [237 MB] Number of unique tags on + strand: 551808
DEBUG @ 24 Feb 2024 11:05:55: [237 MB] Number of peaks in + strand: 0
DEBUG @ 24 Feb 2024 11:05:55: [237 MB] Number of unique tags on - strand: 551808
DEBUG @ 24 Feb 2024 11:05:55: [237 MB] Number of peaks in - strand: 0
DEBUG @ 24 Feb 2024 11:05:55: [237 MB] Chrom b'XI' is discarded!
DEBUG @ 24 Feb 2024 11:05:55: [237 MB] Chromosome: b'XII'
DEBUG @ 24 Feb 2024 11:05:55: [237 MB] Number of unique tags on + strand: 2824287
DEBUG @ 24 Feb 2024 11:05:55: [237 MB] Number of peaks in + strand: 1
DEBUG @ 24 Feb 2024 11:05:55: [237 MB] plus peaks: first - (490051, 23123.0) ... last - (490051, 23123.0)
DEBUG @ 24 Feb 2024 11:05:55: [237 MB] Number of unique tags on - strand: 2824287
DEBUG @ 24 Feb 2024 11:05:55: [237 MB] Number of peaks in - strand: 1
DEBUG @ 24 Feb 2024 11:05:55: [237 MB] minus peaks: first - (490245, 23158.0) ... last - (490245, 23158.0)
DEBUG @ 24 Feb 2024 11:05:55: [237 MB] ip_max: 1; im_max: 1
DEBUG @ 24 Feb 2024 11:05:55: [237 MB] Paired centers: first - 490148 ... second - 490148
DEBUG @ 24 Feb 2024 11:05:55: [237 MB] Number of paired peaks in this chromosome: 1
DEBUG @ 24 Feb 2024 11:05:55: [237 MB] Chromosome: b'XIII'
DEBUG @ 24 Feb 2024 11:05:55: [237 MB] Number of unique tags on + strand: 763883
DEBUG @ 24 Feb 2024 11:05:55: [237 MB] Number of peaks in + strand: 0
DEBUG @ 24 Feb 2024 11:05:55: [237 MB] Number of unique tags on - strand: 763883
DEBUG @ 24 Feb 2024 11:05:55: [237 MB] Number of peaks in - strand: 0
DEBUG @ 24 Feb 2024 11:05:55: [237 MB] Chrom b'XIII' is discarded!
DEBUG @ 24 Feb 2024 11:05:55: [237 MB] Chromosome: b'XIV'
DEBUG @ 24 Feb 2024 11:05:55: [237 MB] Number of unique tags on + strand: 639487
DEBUG @ 24 Feb 2024 11:05:55: [237 MB] Number of peaks in + strand: 0
DEBUG @ 24 Feb 2024 11:05:55: [237 MB] Number of unique tags on - strand: 639487
DEBUG @ 24 Feb 2024 11:05:55: [237 MB] Number of peaks in - strand: 0
DEBUG @ 24 Feb 2024 11:05:55: [237 MB] Chrom b'XIV' is discarded!
DEBUG @ 24 Feb 2024 11:05:55: [237 MB] Chromosome: b'XV'
DEBUG @ 24 Feb 2024 11:05:55: [237 MB] Number of unique tags on + strand: 910215
DEBUG @ 24 Feb 2024 11:05:55: [237 MB] Number of peaks in + strand: 0
DEBUG @ 24 Feb 2024 11:05:55: [237 MB] Number of unique tags on - strand: 910215
DEBUG @ 24 Feb 2024 11:05:55: [237 MB] Number of peaks in - strand: 0
DEBUG @ 24 Feb 2024 11:05:55: [237 MB] Chrom b'XV' is discarded!
DEBUG @ 24 Feb 2024 11:05:55: [237 MB] Chromosome: b'XVI'
DEBUG @ 24 Feb 2024 11:05:55: [237 MB] Number of unique tags on + strand: 789758
DEBUG @ 24 Feb 2024 11:05:55: [237 MB] Number of peaks in + strand: 0
DEBUG @ 24 Feb 2024 11:05:56: [237 MB] Number of unique tags on - strand: 789758
DEBUG @ 24 Feb 2024 11:05:56: [237 MB] Number of peaks in - strand: 0
DEBUG @ 24 Feb 2024 11:05:56: [237 MB] Chrom b'XVI' is discarded!
INFO  @ 24 Feb 2024 11:05:56: [237 MB] #2 Total number of paired peaks: 18
WARNING @ 24 Feb 2024 11:05:56: [237 MB] #2 MACS3 needs at least 100 paired peaks at + and - strand to build the model, but can only find 18! Please make your MFOLD range broader and try again. If MACS3 still can't build the model, we suggest to use --nomodel and --extsize 147 or other fixed number instead.
WARNING @ 24 Feb 2024 11:05:56: [237 MB] #2 Process for pairing-model is terminated!
WARNING @ 24 Feb 2024 11:05:56: [237 MB] # Can't find enough pairs of symmetric peaks to build model!


❯ macs3 callpeak \
>     --treatment "${file_IP_bam}" \
>     --control "${file_in_bam}" \
>     --format "BAMPE" \
>     --gsize "${gsize}" \
>     --keep-dup "auto" \
>     --outdir "${dir_out}" \
>     --bdg \
>     --SPMR \
>     --verbose 3
INFO  @ 24 Feb 2024 11:40:15: [61 MB]
# Command line: callpeak --treatment IP_Q_Esa5_7041.sort-coord.bam --control in_Q_Esa5_7041.sort-coord.bam --format BAMPE --gsize 11624332 --keep-dup auto --outdir macs3_test --bdg --SPMR --verbose 3
# ARGUMENTS LIST:
# name = NA
# format = BAMPE
# ChIP-seq file = ['IP_Q_Esa5_7041.sort-coord.bam']
# control file = ['in_Q_Esa5_7041.sort-coord.bam']
# effective genome size = 1.16e+07
# band width = 300
# model fold = [5, 50]
# qvalue cutoff = 5.00e-02
# The maximum gap between significant sites is assigned as the read length/tag size.
# The minimum length of peaks is assigned as the predicted fragment length "d".
# Larger dataset will be scaled towards smaller dataset.
# Range for calculating regional lambda is: 1000 bps and 10000 bps
# Broad region calling is off
# Paired-End mode is on
# MACS will save fragment pileup signal per million reads

INFO  @ 24 Feb 2024 11:40:15: [61 MB] #1 read fragment files...
INFO  @ 24 Feb 2024 11:40:15: [61 MB] #1 read treatment fragments...
INFO  @ 24 Feb 2024 11:40:17: [77 MB]  1000000 fragments parsed
INFO  @ 24 Feb 2024 11:40:19: [81 MB]  2000000 fragments parsed
INFO  @ 24 Feb 2024 11:40:22: [92 MB]  3000000 fragments parsed
INFO  @ 24 Feb 2024 11:40:24: [100 MB]  4000000 fragments parsed
INFO  @ 24 Feb 2024 11:40:26: [110 MB]  5000000 fragments parsed
INFO  @ 24 Feb 2024 11:40:28: [116 MB]  6000000 fragments parsed
INFO  @ 24 Feb 2024 11:40:30: [125 MB]  7000000 fragments parsed
INFO  @ 24 Feb 2024 11:40:32: [133 MB]  8000000 fragments parsed
INFO  @ 24 Feb 2024 11:40:35: [141 MB]  9000000 fragments parsed
INFO  @ 24 Feb 2024 11:40:37: [146 MB]  10000000 fragments parsed
INFO  @ 24 Feb 2024 11:40:38: [154 MB] 10385714 fragments have been read.
INFO  @ 24 Feb 2024 11:40:43: [154 MB] #1.2 read input fragments...
INFO  @ 24 Feb 2024 11:40:46: [155 MB]  1000000 fragments parsed
INFO  @ 24 Feb 2024 11:40:48: [163 MB]  2000000 fragments parsed
INFO  @ 24 Feb 2024 11:40:50: [171 MB]  3000000 fragments parsed
INFO  @ 24 Feb 2024 11:40:52: [179 MB]  4000000 fragments parsed
INFO  @ 24 Feb 2024 11:40:55: [188 MB]  5000000 fragments parsed
INFO  @ 24 Feb 2024 11:40:57: [196 MB]  6000000 fragments parsed
INFO  @ 24 Feb 2024 11:40:59: [204 MB]  7000000 fragments parsed
INFO  @ 24 Feb 2024 11:41:01: [212 MB]  8000000 fragments parsed
INFO  @ 24 Feb 2024 11:41:03: [220 MB]  9000000 fragments parsed
INFO  @ 24 Feb 2024 11:41:05: [228 MB]  10000000 fragments parsed
INFO  @ 24 Feb 2024 11:41:08: [236 MB]  11000000 fragments parsed
INFO  @ 24 Feb 2024 11:41:10: [244 MB]  12000000 fragments parsed
INFO  @ 24 Feb 2024 11:41:12: [252 MB]  13000000 fragments parsed
INFO  @ 24 Feb 2024 11:41:14: [260 MB] 13741759 fragments have been read.
INFO  @ 24 Feb 2024 11:41:21: [260 MB] #1 mean fragment size is determined as 199.3 bp from treatment
INFO  @ 24 Feb 2024 11:41:21: [260 MB] #1 note: mean fragment size in control is 188.1 bp -- value ignored
INFO  @ 24 Feb 2024 11:41:21: [260 MB] #1 fragment size = 199.3
INFO  @ 24 Feb 2024 11:41:21: [260 MB] #1  total fragments in treatment: 10385714
INFO  @ 24 Feb 2024 11:41:21: [260 MB] #1 calculate max duplicate fragments in single position based on binomial distribution...
INFO  @ 24 Feb 2024 11:41:21: [260 MB] #1  max_dup_tags based on binomial = 7
INFO  @ 24 Feb 2024 11:41:21: [260 MB] #1 filter out redundant fragments by allowing at most 7 identical fragment(s)
INFO  @ 24 Feb 2024 11:41:31: [276 MB] #1  fragments after filtering in treatment: 10382147
INFO  @ 24 Feb 2024 11:41:31: [276 MB] #1  Redundant rate of treatment: 0.00
INFO  @ 24 Feb 2024 11:41:31: [276 MB] #1  total fragments in control: 13741759
INFO  @ 24 Feb 2024 11:41:31: [276 MB] #1  for control, calculate max duplicate fragments in single position based on binomial distribution...
INFO  @ 24 Feb 2024 11:41:31: [276 MB] #1  max_dup_tags based on binomial = 8
INFO  @ 24 Feb 2024 11:41:31: [276 MB] #1 filter out redundant fragments by allowing at most 7 identical fragment(s)
INFO  @ 24 Feb 2024 11:41:44: [298 MB] #1  fragments after filtering in control: 13731787
INFO  @ 24 Feb 2024 11:41:44: [298 MB] #1  Redundant rate of control: 0.00
INFO  @ 24 Feb 2024 11:41:44: [298 MB] #1 finished!
INFO  @ 24 Feb 2024 11:41:44: [298 MB] #2 Build Peak Model...
INFO  @ 24 Feb 2024 11:41:44: [298 MB] #2 Skipped...
INFO  @ 24 Feb 2024 11:41:44: [298 MB] #3 Call peaks...
INFO  @ 24 Feb 2024 11:41:44: [298 MB] #3 Pre-compute pvalue-qvalue table...
INFO  @ 24 Feb 2024 11:42:53: [408 MB] #3 In the peak calling step, the following will be performed simultaneously:
INFO  @ 24 Feb 2024 11:42:53: [408 MB] #3   Write bedGraph files for treatment pileup (after scaling if necessary)... NA_treat_pileup.bdg
INFO  @ 24 Feb 2024 11:42:53: [408 MB] #3   Write bedGraph files for control lambda (after scaling if necessary)... NA_control_lambda.bdg
INFO  @ 24 Feb 2024 11:42:53: [408 MB] #3   --SPMR is requested, so pileup will be normalized by sequencing depth in million reads.
INFO  @ 24 Feb 2024 11:42:53: [408 MB] #3 Call peaks for each chromosome...
INFO  @ 24 Feb 2024 11:42:59: [408 MB] #4 Write output xls file... macs3_test/NA_peaks.xls
INFO  @ 24 Feb 2024 11:42:59: [408 MB] #4 Write peak in narrowPeak format file... macs3_test/NA_peaks.narrowPeak
INFO  @ 24 Feb 2024 11:42:59: [408 MB] #4 Write summits bed file... macs3_test/NA_summits.bed
INFO  @ 24 Feb 2024 11:42:59: [408 MB] Done!
```
</details>
<br />

<a id="conversation-with-gpt4-1"></a>
#### Conversation with GPT4
<a id="me-7"></a>
##### Me
My understanding is that running `macs3 predictd` is not necessary for paired-end ChIP-seq immunoprecipitate data. This is because&mdash;if I am not mistaken&mdash;the fragment length is inferred from the alignment of concordant (or otherwise proper) read pairs; `macs3 predictd` is applicable for single-end ChIP-seq data presumably sequenced to a sufficient depth that a cross-strand correlation analysis can be run to infer an average fragment length for the aligned reads. Is this all correct? Do you have anything to change or add?

<a id="gpt4-3"></a>
##### GPT4
Your understanding is correct. For paired-end ChIP-seq data, the fragment length (or insert size) is typically inferred directly from the aligned read pairs, which makes the `macs3 predictd` step unnecessary. This command is more relevant for single-end data, where fragment lengths can't be directly determined from the reads themselves. In those cases, `macs3 predictd` can help estimate the distribution of fragment lengths by analyzing the distances between peaks on opposite DNA strands, a method that relies on the assumption of a uniform distribution of fragment lengths across the genome. For paired-end data, the direct measurement of the distance between reads provides a more accurate and straightforward method of determining fragment lengths.

<a id="step-3-extend-chip-sample-to-get-chip-coverage-track"></a>
### Step 3: Extend ChIP sample to get ChIP coverage track
Now that you have estimated the fragment length, we next can use the `MACS3 pileup` subcommand to generate a pileup track in BEDGRAPH format for the ChIP sample. Since we are dealing with ChIP-seq data in this tutorial, we need to extend reads in 5' to 3' direction, which is the default behavior of `pileup` function. If you are dealing with some DNase-seq data or you think the cutting site, which is detected by short read sequencing, is in the middle of the fragment you are interested in, you can use the `-B` option to extend the read in both directions. Here is the command to simulate `callpeak` behavior:

```bash
#!/bin/bash

macs3 pileup \
    -i CTCF_ChIP_200K_filterdup.bed \
    -o CTCF_ChIP_200K_filterdup.pileup.bdg \
    --extsize 254
```

The file CTCF_ChIP_200K_filterdup.pileup.bdg now contains the fragment pileup signals for ChIP sample.

<a id="documentation-for-macs3-pileup"></a>
#### Documentation for `macs3 pileup`
<details>
<summary><i>Text: Documentation for `macs3 pileup`</i></summary>

```txt
❯ macs3 pileup -h
usage: macs3 pileup [-h] -i IFILE [IFILE ...] -o OUTPUTFILE [--outdir OUTDIR]
                    [-f {AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE,BEDPE}] [-B] [--extsize EXTSIZE]
                    [--buffer-size BUFFER_SIZE] [--verbose VERBOSE]

options:
  -h, --help            show this help message and exit
  -i IFILE [IFILE ...], --ifile IFILE [IFILE ...]
                        Alignment file. If multiple files are given as '-t A B C', then they will all be read and combined. REQUIRED.
  -o OUTPUTFILE, --ofile OUTPUTFILE
                        Output bedGraph file name. If not specified, will write to standard output. REQUIRED.
  --outdir OUTDIR       If specified all output files will be written to that directory. Default: the current working directory
  -f {AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE,BEDPE}, --format {AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE,BEDPE}
                        Format of tag file, "AUTO", "BED", "ELAND", "ELANDMULTI", "ELANDEXPORT", "SAM", "BAM", "BOWTIE", "BAMPE", or "BEDPE".
                        The default AUTO option will let 'macs3 pileup' decide which format the file is. DEFAULT: "AUTO", MACS3 will pick a
                        format from "AUTO", "BED", "ELAND", "ELANDMULTI", "ELANDEXPORT", "SAM", "BAM" and "BOWTIE". If the format is BAMPE or
                        BEDPE, please specify it explicitly. Please note that when the format is BAMPE or BEDPE, the -B and --extsize options
                        would be ignored.
  -B, --both-direction  By default, any read will be extended towards downstream direction by extension size. So it's [0,size-1] (1-based
                        index system) for plus strand read and [-size+1,0] for minus strand read where position 0 is 5' end of the aligned
                        read. Default behavior can simulate MACS3 way of piling up ChIP sample reads where extension size is set as fragment
                        size/d. If this option is set as on, aligned reads will be extended in both upstream and downstream directions by
                        extension size. It means [-size,size] where 0 is the 5' end of a aligned read. It can partially simulate MACS3 way of
                        piling up control reads. However MACS3 local bias is calculated by maximizing the expected pileup over a ChIP
                        fragment size/d estimated from 10kb, 1kb, d and whole genome background. This option will be ignored when the format
                        is set as BAMPE or BEDPE. DEFAULT: False
  --extsize EXTSIZE     The extension size in bps. Each alignment read will become a EXTSIZE of fragment, then be piled up. Check description
                        for -B for detail. It's twice the `shiftsize` in old MACSv1 language. This option will be ignored when the format is
                        set as BAMPE or BEDPE. DEFAULT: 200
  --buffer-size BUFFER_SIZE
                        Buffer size for incrementally increasing internal array size to store reads alignment information. In most cases, you
                        don't have to change this parameter. However, if there are large number of chromosomes/contigs/scaffolds in your
                        alignment, it's recommended to specify a smaller buffer size in order to decrease memory usage (but it will take
                        longer time to read alignment files). Minimum memory requested for reading an alignment file is about # of CHROMOSOME
                        * BUFFER_SIZE * 8 Bytes. DEFAULT: 100000
  --verbose VERBOSE     Set verbose level. 0: only show critical message, 1: show additional warning message, 2: show process information, 3:
                        show debug messages. If you want to know where are the duplicate reads, use 3. DEFAULT:2
```
</details>
<br />

<a id="on-trying-to-understand-what---format-bampe-does-in-macs3-pileup"></a>
#### On trying to understand what `--format BAMPE` does in `macs3 pileup`
The best I can tell, invoking `macs3 pileup` with `--format BAMPE` facilitates the piling up of whole fragments in general ([source](https://github.com/macs3-project/MACS/issues/145), specifically [here](https://github.com/macs3-project/MACS/issues/145#issuecomment-299528856)):
> With `-f BAMPE` on, MACS2 read the left mate and the insertion length information from BAM file, and discard right mate. With `-f BAM`, MACS2 only keeps the left mate, so it's definitely not a feature you want for your paired-end data.

But what is the coded evidence for this?
- [Pileup.pyx (pileup_and_write_pe)](https://github.com/macs3-project/MACS/blob/3cfb5dc3874e0f2ca3aa287ac5bc764ca0193df3/MACS3/Signal/Pileup.pyx#L165)
- [pileup_cmd.py (pileup_and_write_pe)](https://github.com/macs3-project/MACS/blob/master/MACS3/Commands/pileup_cmd.py#L52C9-L52C28)

<a id="my-proposed-invocation-of-macs3-pileup"></a>
#### My proposed invocation of `macs3 pileup`
<details>
<summary><i>Code: My proposed invocation of `macs3 pileup`</i></summary>

```bash
#!/bin/bash

file_bed="macs3_test/IP_Q_Esa5_7041.sort-coord.bed"
file_bg="${file_bed/.bed/.pileup.bg}"

macs3 pileup \
    --format "BEDPE" \
    --ifile "${file_bed}" \
    --ofile "${file_bg}" \
    --buffer-size 100000 \
    --verbose 3
```
</details>
<br />

<details>
<summary><i>Printed: My proposed invocation of `macs3 pileup`</i></summary>

```txt
❯ file_bed="macs3_test/IP_Q_Esa5_7041.sort-coord.bed"


❯ file_bg="${file_bed/.bed/.pileup.bg}"


❯ macs3 pileup \
>     --format "BEDPE" \
>     --ifile "${file_bed}" \
>     --ofile "${file_bg}" \
>     --buffer-size 100000 \
>     --verbose 3
INFO  @ 23 Feb 2024 12:09:40: [61 MB] # read alignment files...
INFO  @ 23 Feb 2024 12:09:40: [61 MB] # read input file in Paired-end mode.
INFO  @ 23 Feb 2024 12:09:40: [61 MB] # read treatment fragments...
INFO  @ 23 Feb 2024 12:09:42: [139 MB]  1000000 fragments parsed
INFO  @ 23 Feb 2024 12:09:44: [154 MB]  2000000 fragments parsed
INFO  @ 23 Feb 2024 12:09:46: [161 MB]  3000000 fragments parsed
INFO  @ 23 Feb 2024 12:09:48: [181 MB]  4000000 fragments parsed
INFO  @ 23 Feb 2024 12:09:50: [185 MB]  5000000 fragments parsed
INFO  @ 23 Feb 2024 12:09:52: [189 MB]  6000000 fragments parsed
INFO  @ 23 Feb 2024 12:09:54: [209 MB]  7000000 fragments parsed
INFO  @ 23 Feb 2024 12:09:56: [209 MB]  8000000 fragments parsed
INFO  @ 23 Feb 2024 12:09:58: [209 MB]  9000000 fragments parsed
INFO  @ 23 Feb 2024 12:10:00: [222 MB]  10000000 fragments parsed
INFO  @ 23 Feb 2024 12:10:19: [222 MB] # total fragments/pairs in alignment file: 10382147
INFO  @ 23 Feb 2024 12:10:19: [222 MB] # Pileup paired-end alignment file.
INFO  @ 23 Feb 2024 12:10:28: [222 MB] # Done! Check macs3_test/IP_Q_Esa5_7041.sort-coord.pileup.bg
```
</details>
<br />

<a id="step-4-build-a-local-bias-track-from-the-control-sample"></a>
#### Step 4: Build a local bias track from the control sample
By default, MACS3 `callpeak` function computes the local bias by taking the maximum bias from surrounding 1 kb (set by `--slocal`), 10 kb (set by `--llocal`), the size of fragment length `d` (predicted as what you got from `predictd`), and the whole genome background. Here I show you how each of the bias is calculated and how they can be combined using the subcommands.

<a id="the-d-background"></a>
##### The `d` background
Basically, to create the background noise track, you need to extend the control read to both sides (-B option) using pileup function. The idea is that the cutting site from control sample contains the noise representing a region surrounding it. To do this, take half of d you got from `predictd`, 127 (1/2 of 254) for our example, then:
```r
macs3 pileup -i CTCF_Control_200K_filterdup.bed -B --extsize 127 -o d_bg.bdg
```

The file d_bg.bdg contains the `d` background from control.

<a id="the-slocal-background"></a>
##### The slocal background
Next, you can create a background noise track of slocal local window, or 1 kb window by default. Simply imagine that each cutting site (sequenced read) represent a 1 kb (default, you can tweak it) surrounding noise. So:
```r
macs3 pileup -i CTCF_Control_200K_filterdup.bed -B --extsize 500 -o 1k_bg.bdg
```

Note, here 500 is the 1/2 of 1k. Because the ChIP signal track was built by extending reads into `d` size fragments, we have to normalize the 1 kb noise by multiplying the values by `d/slocal`, which is 254/1000 = 0.254 in our example. To do so, use the `bdgopt` subcommand:
```r
macs3 bdgopt -i 1k_bg.bdg -m multiply -p 0.254 -o 1k_bg_norm.bdg
```

The file1k_bg_norm.bdg contains the slocal background from control. Note, we don't have to do this for `d` background because the multiplier is simply 1.

<a id="the-llocal-background"></a>
##### The llocal background
The background noise from larger region can be generated in the same way as slocal backgound. The only difference is the size for extension. MACS3 callpeak by default asks for 10 kb (you can change this value) surrounding window, so:
```r
macs3 pileup -i CTCF_Control_200K_filterdup.bed -B --extsize 5000 -o 10k_bg.bdg
```

The extsize has to be set as 1/2 of llocal. Then, the multiplier now is `d/llocal`, or 0.0254 in our example.
```r
macs3 bdgopt -i 10k_bg.bdg -m multiply -p 0.0254 -o 10k_bg_norm.bdg
```

The file 10k_bg_norm.bdg now contains the slocal background from control.

<a id="the-genome-background"></a>
##### The genome background
The whole genome background can be calculated as (the_number_of_control_reads \* fragment_length) / genome_size. You don't need to run subcommands to build a genome background track since it's just a single value.

<a id="combine-and-generate-the-maximum-background-noise"></a>
##### Combine and generate the maximum background noise
Now all the above background noises have to be combined and the maximum bias for each genomic location need be computed. This is the default behavior of MACS3 `callpeak`, but you can have your own pipeline to include some of them or even make more noise (such as 5 k or 50 k background) then include more tracks. Here is the way to combine and get the maximum bias.

Take the maximum between slocal (1 k) and llocal (10 k) background:
```r
macs3 bdgcmp -m max -t 1k_bg_norm.bdg -c 10k_bg_norm.bdg -o 1k_10k_bg_norm.bdg
```

Then, take the maximum then by comparing with d background:
```r
macs3 bdgcmp -m max -t 1k_10k_bg_norm.bdg -c d_bg.bdg -o d_1k_10k_bg_norm.bdg
```

Finally, combine with the genome wide background using bdgopt subcommand
```r
macs3 bdgopt -i d_1k_10k_bg_norm.bdg -m max -p .0188023 -o local_bias_raw.bdg
```

Now the file local_bias_raw.bdg is a BEDGRAPH file containing the raw local bias from control data.

<a id="on-how-to-do-this-with-paired-end-data"></a>
###### [On how to do this with paired-end data](https://github.com/macs3-project/MACS/issues/356)
@olechnwin Sorry for the confusion! I am discussing the part for ## build local background track from control when you pile up the control data.

```bash
# d background
macs2 pileup -f BEDPE -B -i ${ctrl_prefix}_filterdup.bed 
    -o ${ctrl_prefix}_${prefix}_filterdup_d_bg.pileup.bdg --extsize $d2
# slocal background
macs2 pileup -f BEDPE -B --extsize 500 -i ${ctrl_prefix}_filterdup.bed 
     -o ${ctrl_prefix}_${prefix}_filterdup_1k_bg.pileup.bdg
# large local background default 10k (--extsize 5000)
macs2 pileup -f BEDPE -i ${ctrl_prefix}_filterdup.bed -B --extsize 5000 
    -o ${ctrl_prefix}_${prefix}_filterdup_10k_bg.bdg
```

But let me state clearly. I am not suggesting that you have to replicate the same process in callpeak and this is the beauty of your pipeline since you can control each step. As I mentioned, the default behavior pile up control data as single-end is based on some assumption. And if you have a good reason to pile up control data as paired-end, feel free to do so. My guess is that the difference will be small.

If you want to pile up control data as single-end like callpeak, then the ${ctrl_prefix}_filterdup.bed file has to be broken down into single end BED file like:
```bash
awk -v OFS="\t" -v kmer=$read_len '{print $1,$2,$2+kmer,".",".","+\n"$1,$3-kmer,$3,".",".","-"}' \
  ${ctrl_prefix}_filterdup.bed > ${ctrl_prefix}_filterdup.as_single_end.bed
```

Then:
```bash
# d background
macs2 pileup -f BED -B -i ${ctrl_prefix}_filterdup.as_single_end.bed \
  -o ${ctrl_prefix}_${prefix}_filterdup_d_bg.pileup.bdg --extsize $d2
# slocal background
macs2 pileup -f BED -B -i ${ctrl_prefix}_filterdup.as_single_end.bed \
  -o ${ctrl_prefix}_${prefix}_filterdup_1k_bg.pileup.bdg --extsize 500
# large local background default 10k (--extsize 5000)
macs2 pileup -f BED -B -i ${ctrl_prefix}_filterdup.as_single_end.bed \
  -o ${ctrl_prefix}_${prefix}_filterdup_10k_bg.pileup.bdg --extsize 5000
# divide by 2
awk -v OFS="\t" '{print $1,$2,$3/2}' ${ctrl_prefix}_${prefix}_filterdup_d_bg.pileup.bdg \
  > ${ctrl_prefix}_${prefix}_filterdup_d_bg.pileup.halved.bdg
mv ${ctrl_prefix}_${prefix}_filterdup_d_bg.pileup.halved.bdg  ${ctrl_prefix}_${prefix}_filterdup_d_bg.pileup.bdg

awk -v OFS="\t" '{print $1,$2,$3/2}' ${ctrl_prefix}_${prefix}_filterdup_1k_bg.pileup.bdg \
  > ${ctrl_prefix}_${prefix}_filterdup_1k_bg.pileup.halved.bdg; 
mv ${ctrl_prefix}_${prefix}_filterdup_1k_bg.pileup.halved.bdg  ${ctrl_prefix}_${prefix}_filterdup_1k_bg.pileup.bdg

awk -v OFS="\t" '{print $1,$2,$3/2}' ${ctrl_prefix}_${prefix}_filterdup_10k_bg.pileup.bdg \
  > ${ctrl_prefix}_${prefix}_filterdup_10k_bg.pileup.halved.bdg; 
mv ${ctrl_prefix}_${prefix}_filterdup_10k_bg.pileup.halved.bdg  ${ctrl_prefix}_${prefix}_filterdup_10k_bg.pileup.bdg
```
After that, just keep using the same code in your pipeline. (btw, I didn't test my code, so there may be a bug...)


<a id="d2-was-defined-here"></a>
###### d2 was defined here...
@ taoliu,
Thank you so much for looking into this. Sorry for the delay. I got distracted...
See below for my pipeline.
If I understand you correctly, you were suggesting to do the subcommands instead of callpeak? Is there a way to get pileup and fold enrichment when doing the subcommand?

My pipeline, some commands were broken up into multiple lines for easy viewing:
```bash
# specify read length
read_len=150

# specify peaks cut off
log_qval_cutoff=1.301
log_pval_cutoff=2

# convert sam to bam for ChIP sample
samtools view -S -b $file | samtools sort - -o ${prefix}.sorted.bam
# remove duplicates
macs2 filterdup -f BAMPE -i ${prefix}.sorted.bam --keep-dup=1 -o ${prefix}_filterdup.bed

# convert sam to bam for control sample
samtools view -S -b $ctrl_file | samtools sort - -o ${ctrl_prefix}.sorted.bam
# remove duplicates
macs2 filterdup -f BAMPE -i ${ctrl_prefix}.sorted.bam --keep-dup=1 
    -o ${ctrl_prefix}_filterdup.bed

## generate ChIP coverage tracks
# get d (fragment size) from each peak file then calculate p normalization factor
# for background
d=$(cat $peaks_dir/${prefix}_peaks.xls | awk '/# d/ {print $4}')
d2=$(echo $d/2 | bc)
macs2 pileup -f BEDPE -B -i ${prefix}_filterdup.bed 
     -o ${prefix}_filterdup.pileup.bdg --extsize $d

## build local background track from control
# d background
macs2 pileup -f BEDPE -B -i ${ctrl_prefix}_filterdup.bed 
    -o ${ctrl_prefix}_${prefix}_filterdup_d_bg.pileup.bdg --extsize $d2
# slocal background
macs2 pileup -f BEDPE -B --extsize 500 -i ${ctrl_prefix}_filterdup.bed 
     -o ${ctrl_prefix}_${prefix}_filterdup_1k_bg.pileup.bdg
# normalize the small local (slocal) background
p=$(echo "scale=4;"$d/1000 | bc)
macs2 bdgopt -i ${ctrl_prefix}_${prefix}_filterdup_1k_bg.pileup.bdg -m multiply -p $p 
    -o ${ctrl_prefix}_${prefix}_filterdup_1k_bg_norm.bdg
# large local background default 10k (--extsize 5000)
macs2 pileup -f BEDPE -i ${ctrl_prefix}_filterdup.bed -B --extsize 5000 
    -o ${ctrl_prefix}_${prefix}_filterdup_10k_bg.bdg
# normalize the large local (llocal) background
p=$(echo "scale=4;"$d/10000 | bc)
macs2 bdgopt -i ${ctrl_prefix}_${prefix}_filterdup_10k_bg.bdg -m multiply -p $p 
    -o ${ctrl_prefix}_${prefix}_filterdup_10k_bg_norm.bdg
# combine and generate maximum background noise
macs2 bdgcmp -m max -t ${ctrl_prefix}_${prefix}_filterdup_1k_bg_norm.bdg 
    -c ${ctrl_prefix}_${prefix}_filterdup_10k_bg_norm.bdg 
    -o ${ctrl_prefix}_${prefix}_filterdup_1k_10k_bg_norm.bdg
# Then, take the maximum then by comparing with d background:
macs2 bdgcmp -m max -t ${ctrl_prefix}_${prefix}_filterdup_1k_10k_bg_norm.bdg 
    -c ${ctrl_prefix}_${prefix}_filterdup_d_bg.pileup.bdg 
    -o ${ctrl_prefix}_${prefix}_filterdup_d_1k_10k_bg_norm.bdg

## finally combine with genome wide background
#p=the_number_of_control_reads*fragment_length/genome_size
ctrl_num_read=$(wc -l ${ctrl_prefix}_filterdup.bed | awk '{print $1}') 
genome_size=2700000000
p=$(echo "scale=4;"$ctrl_num_read*$d/$genome_size | bc)
echo $ctrl_num_read $genome_size $p
macs2 bdgopt -i ${ctrl_prefix}_${prefix}_filterdup_d_1k_10k_bg_norm.bdg -m max -p $p
    -o ${ctrl_prefix}_${prefix}_filterdup_local_bias_raw.bdg        

# scaling chip and control based on spike in
treat_num_read=$(wc -l ${prefix}_filterdup.bed | awk '{print $1}')
treat_spike_num_read=$(wc -l ${prefix}.ecoli_filterdup.bed | awk '{print $1}')
ctrl_spike_num_read=$(wc -l ${ctrl_prefix}.ecoli_filterdup.bed | awk '{print $1}')
min_spike_num_read=$(wc -l *ecoli_filterdup.bed | sort -nk 1 | head -n 1 | awk '{print $1}')
scale_factor_spike_treat=$(echo "scale=4;"$min_spike_num_read/$treat_spike_num_read | bc)
scale_factor_spike_ctrl=$(echo "scale=4;"$min_spike_num_read/$ctrl_spike_num_read | bc)
p=$(echo "scale=4;"$scale_factor_spike_treat/$scale_factor_spike_ctrl | bc)
echo spike in scale factor=$p
if (( $(echo $p '>' 1 |bc -l) )); then
    p=$(echo "scale=10;"1/$p | bc)
    echo scale down control...
    # scale down control
    $run macs2 bdgopt -i ${ctrl_prefix}_${prefix}_filterdup_local_bias_raw.bdg -m multiply -p $p \
            -o ${ctrl_prefix}_${prefix}_filterdup_local_lambda.bdg
    $run macs2 bdgopt -i ${prefix}_filterdup.pileup.bdg -m multiply -p 1 -o ${prefix}_filterdup_scale.pileup.bdg
    echo scaled spike-in control= $(echo $p*$ctrl_spike_num_read | bc)
    echo scaled spike-in treat= $(echo $treat_spike_num_read | bc)
else
    echo scale down treatment
    # scale down treatment
    $run macs2 bdgopt -i ${prefix}_filterdup.pileup.bdg -m multiply -p $p -o ${prefix}_filterdup_scale.pileup.bdg
    $run macs2 bdgopt -i ${ctrl_prefix}_${prefix}_filterdup_local_bias_raw.bdg -m multiply -p 1 \
            -o ${ctrl_prefix}_${prefix}_filterdup_local_lambda.bdg
    echo scaled spike-in control= $(echo $p*$treat_spike_num_read | bc)
    echo scaled spike-in treat= $(echo $ctrl_spike_num_read | bc)
fi

# compare chip and control and calculate q and pvalue
macs2 bdgcmp -t ${prefix}_filterdup_scale.pileup.bdg 
       -c  ${ctrl_prefix}_${prefix}_filterdup_local_lambda.bdg -m qpois 
        -o ${prefix}_qvalue.bdg
macs2 bdgcmp -t ${prefix}_filterdup_scale.pileup.bdg 
       -c  ${ctrl_prefix}_${prefix}_filterdup_local_lambda.bdg -m ppois 
       -o ${prefix}_pvalue.bdg
#Call peaks on score track using a cutoff
macs2 bdgpeakcall -i ${prefix}_qvalue.bdg -c $log_qval_cutoff -l $d -g $read_len 
        -o ${prefix}_qval_peaks.bed
macs2 bdgpeakcall -i ${prefix}_pvalue.bdg -c $log_pval_cutoff -l $d -g $read_len 
        -o ${prefix}_pval_peaks.bed
```

edited: to fix mixed up scaling.
p.s: with `-f BEDPE`, you don't really need to specify `--extsize` anymore.

<a id="note-to-self"></a>
###### Note to self
This is weird as hell. The guy calculates d2 from already having run `callpeaks`. But the whole point of doing this step-by-step work is to avoid calling `callpeaks`.


Up to this point, it seems like the steps for the control ("in") data should be:
```bash
#!/bin/bash

file_in_bam="in_${file_IP_bam##IP_}"
file_in_bed="${file_in/.bam/.bed}"
gsize=11624332
dir_out="macs3_test"

#  Step 1: Filter duplicates
macs3 filterdup \
    --ifile "${file_in_bam}" \
    --format "BAMPE" \
    --gsize "${gsize}" \
    --pvalue "1e-5" \
    --keep-dup "auto" \
    --buffer-size 100000 \
    --outdir "${dir_out}" \
    --ofile "${file_in_bed}" \
    --verbose 3


```

