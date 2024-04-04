#!/usr/bin/env Rscript

#  calculate_coverage_log2-frac-fIPx-finx.R
#  KA

#  Load necessary libraries ===================================================
library(GenomicRanges)
library(tidyverse)


#  Define functions ===========================================================
#  Function to read in bedGraph as dataframe
read_bedGraph <- function(bedGraph) {
    df <- read.table(
        bedGraph,
        header = FALSE,
        skip = 1,
        col.names = c("chr", "start", "end", "coverage")
    ) %>%
        tidyr::unite(
            combined,
            c(chr, start, end),
            sep = "_",
            remove = FALSE
        ) %>%
        dplyr::relocate(combined, .after = coverage)
    
    return(df)
}


#  Function to organize bedGraph dataframe as GRanges
generate_GRanges <- function(df) {
    gr <- GenomicRanges::GRanges(
        seqnames = df$chr,
        ranges = IRanges::IRanges(start = df$start, end = df$end),
        coverage = df$coverage
    )
}


#  Do the main work ===========================================================
#  Initialize variable to set appropriate working directory
d_repo <- "/Users/kalavattam/repos"
d_proj <- "2023_tutorial_ChIP-seq"
d_exp <- "siQ-ChIP_2024-0329/seqL"

#  Set appropriate working directory
setwd(file.path(d_repo, d_proj, d_exp))
# getwd()  # list.files()

#  Load data
df_IP <- read_bedGraph("NormCovIP_G1_7750_Hmo1_seqL.bedGraph")
# gr_IP <- generate_GRanges(df_IP)

df_in <- read_bedGraph("NormCovIN_G1_7750_Hmo1_seqL.bedGraph")
# gr_in <- generate_GRanges(df_in)

test <- dplyr::full_join(df_IP, df_in, by = "combined")

#  Combine and reduce the intervals from both datasets
all_intervals <- GenomicRanges::reduce(c(gr_IP, gr_in))

#  Initialize a default coverage of 0 for all intervals
mcols(all_intervals)$cvg_IP <- 0
mcols(all_intervals)$cvg_in <- 0

#  Fill in actual coverage values from the original datasets
mcols(all_intervals)$cvg_IP[match(gr_IP, all_intervals)] <-
    mcols(gr_IP)$coverage
mcols(all_intervals)$cvg_in[match(gr_in, all_intervals)] <-
    mcols(gr_in)$coverage
