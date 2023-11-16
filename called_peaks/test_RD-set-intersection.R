#!/usr/bin/env Rscript

#  test_RD-set-intersection.R
#  KA

#  This script is designed to be called in directory called_peaks with RStudio 
if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
    dir_script <- dirname(rstudioapi::getSourceEditorContext()$path)
}

setwd(dir_script)  # Set the wrk directory to directory called_peaks
# getwd()

#  Load required libraries
library(GenomicRanges)
library(IRanges)
library(tidyverse)
#NOTE We load venneuler below

#  Load in narrowPeak files that are in directory called_peaks
Esa1_peaks_file <- "Esa1_merged_IP_peaks.narrowPeak"
Gcn5_peaks_file <- "Gcn5_merged_IP_peaks.narrowPeak"
Rpd3_peaks_file <- "Rpd3_merged_IP_peaks.narrowPeak"


#  Read in narrowPeak files ===================================================
#  Define function to read narrowPeak files and convert them to GRanges
col_names <- c(  # See biostars.org/p/102710/ for more information
    "chr", "start", "stop", "name", "score", "strand", "signal_value",
    "p_value", "q_value", "peak"
)

read_narrowPeak <- function(file_path) {
    readr::read_tsv(
        file_path, col_names = col_names, show_col_types = FALSE
    ) %>%
        as.data.frame() %>%
        makeGRangesFromDataFrame(keep.extra.columns = TRUE)
}

#  Load the narrowPeak files as GRanges objects
Esa1_peaks <- read_narrowPeak(Esa1_peaks_file)
Gcn5_peaks <- read_narrowPeak(Gcn5_peaks_file)
Rpd3_peaks <- read_narrowPeak(Rpd3_peaks_file)


#  Example #1: Standard overlap analysis ======================================
#  Example: finding overlaps between Esa1 and Rpd3
overlaps_Esa1_Rpd3 <- GenomicRanges::findOverlaps(Esa1_peaks, Rpd3_peaks)

#  To get the actual overlapping ranges, use the queryHits function
Esa1_overlap_Rpd3_ranges <- Esa1_peaks[queryHits(overlaps_Esa1_Rpd3)]
#  Esa1_overlap_Rpd3_ranges will contain those regions in Esa1_peaks that
#+ overlap with Rpd3_peaks.
#+ 
#+ Here's how it works:
#+ 1. findOverlaps(Esa1_peaks, Rpd3_peaks) returns an object that contains the
#+    indices of the overlapping ranges between Esa1_peaks and Rpd3_peaks.
#+    Specifically, it gives two sets of indices:
#+    a. Query hits (the indices of Esa1_peaks that have overlaps)
#+    b. Subject hits (the indices of Rpd3_peaks that overlap with Esa1_peaks)
#+ 
#+ 2. queryHits(overlaps_Esa1_Rpd3) extracts just the indices of Esa1_peaks
#+    that have overlaps with Rpd3_peaks.
#+ 
#+ 3. Esa1_peaks[queryHits(overlaps_Esa1_Rpd3)] then uses these indices to
#+    subset Esa1_peaks, giving you the actual GRanges objects (genomic ranges)
#+    from Esa1_peaks that overlap with those in Rpd3_peaks.
#+ 
#+ So, Esa1_overlap_Rpd3_ranges will contain the subset of Esa1_peaks that have
#+ some degree of overlap with Rpd3_peaks. This allows you to focus on those
#+ specific regions of the genome where both Esa1 and Rpd3 are present, which
#+ is useful for the type of comparative analysis described in your research
#+ questions.

#  You can also count the number of overlaps
overlaps_Esa1_Rpd3_tally <- length(overlaps_Esa1_Rpd3)  # 1179


#  Example #2: An overlap-with-slack analysis =================================
#  To allow for some "slack" or flexibility in defining overlaps, you can use
#+ the maxgap and minoverlap arguments in the findOverlaps function. These
#+ parameters define the conditions under which two ranges are considered to be
#+ overlapping. Here's a brief overview of each:
#+     1. maxgap: This parameter defines the maximum gap allowed between two
#+        ranges for them to be considered as overlapping. For instance, if
#+        maxgap is set to 10, two ranges will be considered overlapping even
#+        if they are up to 10 base pairs apart. This can be useful when you
#+        want to include interactions or effects that might occur between
#+        genomic regions that are close but not exactly adjacent. By default,
#+        maxgap is -1.
#+ 
#+     2. minoverlap: This parameter specifies the minimum number of base pairs
#+        that must overlap for two ranges to be considered as overlapping. For
#+        instance, if minoverlap is set to 5, then two ranges need to overlap
#+        by at least 5 base pairs to be counted as an overlap. By default,
#+        this value is 0. (The reason it is set to 0 is because maxgap is -1
#+        by default.)

#  Define the amount of slack for overlaps
min_overlap <- 500  # Minimum overlap of 500 base pairs

#  Find overlaps with slack
overlaps_Esa1_Rpd3_min_500 <- findOverlaps(
    Esa1_peaks, Rpd3_peaks, minoverlap = min_overlap
)
overlaps_Esa1_Rpd3_min_500_tally <- length(overlaps_Esa1_Rpd3_min_500)  # 479
#  When using the default type parameter (which is "any"), you cannot
#+ simultaneously set non-default values for both maxgap and minoverlap. The
#+ function requires that at least one of these parameters be left at its
#+ default setting.

#  Here's a breakdown of what you can do:
#+ 1. Use default for one parameter: You can choose to set either maxgap or
#+    minoverlap to a non-default value, but not both. Decide which parameter
#+    is more critical for your analysis and set the other one to its default
#+    value.
#+    - If minoverlap is more important (i.e., you want to ensure a minimum
#+      number of overlapping bases), set minoverlap to your desired value and
#+      leave maxgap as its default.
#+    - If maxgap is more important (i.e., you're interested in regions that
#+      are close, even if they don't overlap by a minimum number of bases),
#+      set maxgap to your desired value and leave minoverlap as its default.
#+ 
#+ 2. Change the overlap type: The constraint is specific to the default
#+    overlap type ("any"). You can work with other type values, like
#+    "within", "start", "end", etc., which allow for setting both maxgap and
#+    minoverlap to non-default values. However, be aware that changing the
#+    type changes the definition of what constitutes an overlap, which might
#+    or might not be suitable for your analysis.


#  Example #3: A more specific overlap-with-slack analysis ====================
type <- "start"
min_overlap <- 500  # Minimum overlap of 500 base pairs
max_gap <- 10

overlaps_Esa1_Rpd3_start_min_500_max_10 <- findOverlaps(
    Esa1_peaks, Rpd3_peaks,
    type = type, minoverlap = min_overlap, maxgap = max_gap
)

#  Here, findOverlaps is configured to identify overlaps based on the following
#+ criteria:
#+ 1. type (type = "start"): This setting specifies that the function should
#+    consider an overlap if the start positions of the ranges in Esa1_peaks
#+    and Rpd3_peaks are within a specified distance of each other. The "start"
#+    type focuses on the starting points of the ranges.
#+ 
#+ 2. Minimum overlap (minoverlap = 500): This criterion requires that for an
#+    overlap to be considered valid, the overlapping regions must be at least
#+    500 base pairs long. This means that even if the start positions are
#+    within the specified distance, the function will only report an overlap
#+    if the actual overlapping segment is at least 500 bp.
#+ 
#+ 3. Maximum gap (maxgap = 10): This parameter allows for a gap of up to 10
#+    base pairs between the start positions of the ranges. In other words,
#+    even if the start positions of the ranges in Esa1_peaks and Rpd3_peaks
#+    are not exactly aligned, they can be considered overlapping if they are
#+    within 10 base pairs of each other.
#+ 
#+ The resulting "Hits" object (which is assigned to variable
#+ overlaps_Esa1_Rpd3_start_min_500_max_10) represents the pairs of ranges
#+ (from Esa1_peaks and Rpd3_peaks) that meet these criteria. The "queryHits"
#+ and "subjectHits" indicate the indices of the ranges in Esa1_peaks (the
#+ "query") and Rpd3_peaks (the "subject"), respectively, that form these
#+ overlapping pairs. (See below for more information on "query" and "subject"
#+ assignments below.)
#+ 
#+ Using the queryLength and subjectLength functions, you can get information
#+ on the total number of ranges in each of the two sets (Esa1_peaks and
#+ Rpd3_peaks).
#+ 
#+ Note that the "query" is the first first argument passed to the findOverlaps
#+ function. In our case, that's Esa1_peaks. That is, GRanges object Esa1_peaks
#+ is the set of genomic ranges that we are testing for overlaps against
#+ another set (Rpd3_peaks).
#+ 
#+ On the other hand, "subject" is the second argument in the findOverlaps
#+ function. In our scenario, GRanges object Rpd3_peaks is the subject. It's
#+ the set of genomic ranges against which the query is being compared to find
#+ overlaps.


#  Example #4: Let's isolate intersecting--and not--peaks (Example #1) ========
Esa1_overlap_Rpd3_ranges <- Esa1_peaks[queryHits(overlaps_Esa1_Rpd3)]
Rpd3_overlap_Esa1_ranges <- Rpd3_peaks[subjectHits(overlaps_Esa1_Rpd3)]
#  These are the same.

peaks_in_Esa1_not_Rpd3 <- GenomicRanges::setdiff(Esa1_peaks, Rpd3_peaks)
peaks_in_Rpd3_not_Esa1 <- GenomicRanges::setdiff(Rpd3_peaks, Esa1_peaks)

#  Like above, it's possible to allow for slack in identifying non-overlapping
#+ regions between two sets of genomic ranges. The mechanism is a bit different
#+ though. We have to adjust the Esa1_peaks and Rpd3_peaks before performing
#+ the GenomicRanges::setdiff operation. For example, we can to expand the
#+ Rpd3_peaks by a certain number of base pairs (the "slack") on both sides,
#+ and then use setdiff to find Esa1_peaks that do not overlap with these
#+ expanded Rpd3_peaks.
#+ 
#+ Here's how we can implement this:
#+ 1. Define the slack: Decide how many base pairs we want to add to each side
#+    of the Rpd3_peaks.
#+ 
#+ 2. Expand Rpd3_peaks: Use the GenomicRanges::resize function, specifying the
#+    amount of slack and ensuring that the resizing is centered.
#+ 
#+ 3. Perform setdiff: Use the expanded Rpd3_peaks in the setdiff function.
#+ 
#+ For example:
#  Define slack
slack <- 100

#  Expand Rpd3 peaks by slack
expanded_Rpd3_peaks <- GenomicRanges::resize(
    Rpd3_peaks, width = width(Rpd3_peaks) + 2 * slack, fix = "center"
)

#  Find Esa1 peaks that do not overlap the expanded ("slackened") Rpd3 peaks
peaks_in_Esa1_not_expanded_Rpd3 <- GenomicRanges::setdiff(
    Esa1_peaks, expanded_Rpd3_peaks
)


#  Example #5: Contextualizing ranges calculated above ========================
#  Load in gff3 file for S. cerevisiae/S. pombe -------------------------------
col_names <- c(
    "chr", "source", "type", "start", "end", "frame", "strand", "score",
    "attributes"
)
exclude_types <- c(
    "CDS", "mRNA", "ncRNA", "rRNA", "snoRNA", "snRNA", "telomerase_RNA",
    "transposable_element", "tRNA"
)

f_gtf <- "combined_SC_SP.clean.intelligible.gff3.gz"
t_gtf <- readr::read_tsv(
    f_gtf,
    comment = "#",
    col_names = col_names,
    show_col_types = FALSE
) %>%
    dplyr::select(-c(source, frame, score)) %>%
    dplyr::mutate(to = strsplit(attributes, ";")) %>%
    tidyr::unnest(to) %>%
    dplyr::filter(stringr::str_detect(to, "ID=")) %>%
    dplyr::group_by(attributes) %>%
    dplyr::mutate(row = dplyr::row_number()) %>%
    tidyr::spread(row, to) %>%
    dplyr::mutate(features = stringr::str_remove_all(`1`, "ID=")) %>%
    dplyr::ungroup() %>%
    dplyr::select(-c(attributes, `1`)) %>%
    dplyr::relocate(c(chr, start, end, strand, features, type))

check_tibble <- FALSE
if(base::isTRUE(check_tibble)) {
    check <- t_gtf %>%
        dplyr::group_by(type) %>%
        dplyr::summarize(no = dplyr::n())
}

t_gtf <- t_gtf %>%
    dplyr::filter(!type %in% exclude_types) %>%
    dplyr::filter(!grepl("_ARS", features)) %>%
    dplyr::filter(!grepl("_blocked_reading_frame", features)) %>%
    dplyr::filter(!grepl("_centromere", features)) %>%
    dplyr::filter(!grepl("_mating_type_region", features)) %>%
    dplyr::filter(!grepl("_silent_mating_type_cassette_array", features)) %>%
    dplyr::filter(!grepl("_telomere", features)) %>%
    dplyr::arrange(chr, start, end, strand)

t_gtf$type <- gsub("_gene", "", t_gtf$type)

check_tibble <- FALSE
if(base::isTRUE(check_tibble)) {
    check_after <- t_gtf %>%
        dplyr::group_by(type) %>%
        dplyr::summarize(no = dplyr::n())
}

rm(col_names, exclude_types)


#  Convert tibble for S. cerevisiae/S. pombe features to GRanges object -------
G_gtf <- t_gtf %>%
    as.data.frame() %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)


#  Use findOverlaps to associate feature info with ChIP-seq ranges ------------
#  Define function to concatenate gene names for overlapping ranges
concatenate_genes <- function(indices, gtf) {
    gene_names <- mcols(gtf[indices])$features
    return(paste(gene_names, collapse = "; "))
}

#  Find overlaps between Esa1_overlap_Rpd3_ranges and G_gtf
overlaps <- findOverlaps(Esa1_overlap_Rpd3_ranges, G_gtf)

#  Initialize a vector to store gene names
gene_names <- vector("character", length(Esa1_overlap_Rpd3_ranges))

#  Iterate over each range in Esa1_overlap_Rpd3_ranges
for (i in seq_along(Esa1_overlap_Rpd3_ranges)) {
    # Find indices of overlapping genes in G_gtf
    overlapping_indices <- subjectHits(overlaps[queryHits(overlaps) == i])

    if (length(overlapping_indices) > 0) {
        #  If one or more genes overlap a range, then associate those genes
        #+ with the appropriate index and, if more than one gene, concatenate
        #+ gene names
        gene_names[i] <- concatenate_genes(overlapping_indices, G_gtf)
    } else {
        #  Otherwise, assign placeholder to signify the range does not overlap
        #+ any genes
        gene_names[i] <- "No overlap"
    }
}

#  Assign gene names to Esa1_overlap_Rpd3_ranges
Esa1_overlap_Rpd3_ranges$gene_names <- gene_names

#  View the updated Esa1_overlap_Rpd3_ranges with gene information
Esa1_overlap_Rpd3_ranges


#  Example #6: Draw a Venn diagram with base R graphics =======================
#  (But this way to do it sucks.)

overlapping_peaks_Rpd3_Esa1 <- Esa1_peaks[queryHits(overlaps_Esa1_Rpd3)]
# overlapping_peaks_Rpd3_Esa1 <- Rpd3_peaks[subjectHits(overlaps_Esa1_Rpd3)]
# Same as above.

peaks_in_Esa1_not_Rpd3 <- GenomicRanges::setdiff(Esa1_peaks, Rpd3_peaks)
peaks_in_Rpd3_not_Esa1 <- GenomicRanges::setdiff(Rpd3_peaks, Esa1_peaks)

#  Prepare the numbers for the Venn diagram
count_Esa1_only <- length(peaks_in_Esa1_not_Rpd3)
count_Rpd3_only <- length(peaks_in_Rpd3_not_Esa1)
count_overlap <- length(overlapping_peaks_Rpd3_Esa1)

#  Basic plot setup
plot(1, type = "n", xlab = "", ylab = "", xlim = c(0, 3), ylim = c(0, 2))
rect(0, 0, 3, 2, col = "white")

#  Draw the circles
symbols(
    c(1, 2), c(1, 1),
    circles = c(1, 1),
    add = TRUE, inches = FALSE,
    fg = "black"
)

#  Add text labels
text(1, 1.5, "Esa1", cex = 0.8)
text(2, 1.5, "Rpd3", cex = 0.8)
text(1, 0.5, count_Esa1_only, cex = 1)
text(2, 0.5, count_Rpd3_only, cex = 1)
text(1.5, 1, count_overlap, cex = 1)


#  Example #7: Draw a Venn diagram with package venneuler =====================
#  (This way does not suck. This is the program Alison used to make the nice
#+ plots in the NAR paper.)
library(venneuler)

#  Counts for each set and intersection
venn_data <- c(
    "Esa1" = length(peaks_in_Esa1_not_Rpd3),
    "Rpd3" = length(peaks_in_Rpd3_not_Esa1),
    "Esa1&Rpd3" = length(overlapping_peaks_Rpd3_Esa1)
)

#  Generate the Venn diagram
venn_diagram <- venneuler(venn_data)

#  Plot the diagram: This looks much nicer and much more sensible
plot(venn_diagram)
