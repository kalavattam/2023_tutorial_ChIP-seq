#!/usr/bin/env Rscript

#  add_cluster-number_matrix.R
#  KA

#  Draws on work in test_RD-data-viz.md; can use this code as a template to
#+ associate cluster group numbers with rows and columns in the matrices

#  Load necessary library
library(tidyverse)

#  Set up variables
p_proj <- "/Users/kalavatt/projects-etc" 
p_repo <- "2023_rDNA"
p_exp <- "results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz"

f_bed <- "all_merged_Steinmetz_kmeans-7.bed"
f_mat <- "all_merged_Steinmetz_kmeans-7.mat.gz"
col_mat <- c(
    "chr", "start", "stop", "name", "score", "strand",
    paste("bin", c(1:600), sep = "_")  # Apparently, 600 values/row (feature)
)

#  Go to work directory
file.path(p_proj, p_repo, p_exp) %>% setwd()

#  Load bed and matrix files as tibble data objects
t_bed <- readr::read_tsv(f_bed, show_col_types = FALSE)
t_mat <- readr::read_tsv(
    f_mat, col_names = col_mat, skip = 1, show_col_types = FALSE
)

#  Associate cluster group numbers with matrix values
t_full <- dplyr::right_join(
    t_bed[, colnames(t_bed) %in% c("name", "deepTools_group")],
    t_mat,
    by = "name"
) %>%
    dplyr::relocate(c("name", "deepTools_group"), .after = "stop")
