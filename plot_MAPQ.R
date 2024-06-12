#!/usr/bin/env Rscript

#  plot_MAPQ.R
#  KA

#  Run in "interactive mode" or not
interact <- TRUE

#  Load required libraries
library(argparse)
library(scales)
library(tidyverse)

#  Set options
options(scipen = 999)

#  Parse arguments in "interactive mode" or "command line mode"
if (interact) {
    dir_repo <- "/Users/kalavattam/repos/2023_tutorial_ChIP-seq"
    dir_exp <- file.path(dir_repo, "alignment-tallies_bwa")
    infile <- "IP_G1_Hho1_6336.list-tally.txt"
    plot_type <- "flag"  # "MAPQ"  # "MAPQ-by-flag"
    outfile <- file.path(
        dir_exp,
        paste(tools::file_path_sans_ext(infile), plot_type, "pdf", sep = ".")
    )
    height <- 6
    width <- 10
    per_chr <- TRUE
    chr_list <- c("I", "XII", "SP_III")  # chr_list <- NULL
    zoom_pct <- 100  # zoom_pct <- 1.5
    
    parse_args <- list(
        infile = file.path(dir_exp, infile),
        plot_type = plot_type,
        outfile = outfile,
        height = height,
        width = width,
        per_chr = per_chr,
        chr_list = chr_list,
        zoom_pct = zoom_pct
    )
} else {
    #  Function to parse arguments from the command line
    parse_args <- function() {
        parser <- ArgumentParser(
            description = "Plot ... from a tab-separated text file of....")
        parser$add_argument(
            "-i", "--infile",
            type = "character",
            help = "Path to the input tab-separated text file."
        )
        parser$add_argument(
            "-o", "--outfile",
            type = "character",
            help = paste(
                "Path to save the output plot. Extension can be any handled",
                "by ggplot2."
            )
        )
        parser$add_argument(
            "-p", "--plot_type",
            type = "character",
            default = "MAPQ",
            choices = c("MAPQ", "MAPQ-by-flag", "flags"),
            help = "Type of plot: 'MAPQ', 'MAPQ-by-flag', or 'flags'."
        )
        parser$add_argument(
            "-h", "--height",
            type = "integer",
            default = 6,
            help = "Height of the plot in inches (default: 6)."
        )
        parser$add_argument(
            "-w", "--width",
            type = "integer",
            default = 10,
            help = "Width of the plot in inches (default: 10)."
        )
        parser$add_argument(
            "-c", "--per_chr",
            action = "store_true",
            help = "Generate per-chromosome plots."
        )
        parser$add_argument(
            "-l", "--chr_list",
            type = "character",
            nargs = "*",
            help = paste(
                "List of chromosomes to include in per-chromosome plots. If",
                "not specified, all chromosomes will be included."
            )
        )
        parser$add_argument(
            "-z", "--zoom_pct",
            type = "numeric",
            default = 100,
            help = "Zoom percentage for the y-axis (default: 100)."
        )
        args <- parser$parse_args()
        return(args)
    }
}


#  Function to apply zoom to the y-axis
apply_zoom <- function(plot, zoom_pct) {
    if (!is.null(zoom_pct) && zoom_pct < 100) {
        max_y <- max(plot$data$tally, na.rm = TRUE)
        plot <- plot + coord_cartesian(ylim = c(0, max_y * (zoom_pct / 100)))
    }
    
    return(plot)
}


#  Main function to process data and create plot
main <- function() {
    if (interact) {
        args <- parse_args
    } else {
        args <- parse_args()
    }
    
    #  Read the input file
    data <- readr::read_tsv(
        args$infile,
        col_names = c("chr", "flag", "MAPQ", "tally"),
        show_col_types = FALSE
    )
    
    #  Generate the plot based on the selected plot type
    if (args$plot_type == "MAPQ") {
        plot <- ggplot2::ggplot(
            data, aes(x = as.factor(MAPQ), y = tally)
        ) +
            geom_bar(stat = "identity") +
            labs(title = "MAPQ distribution", x = "MAPQ", y = "tally") +
            theme_minimal() +
            scale_x_discrete(
                breaks = function(x) x[as.numeric(x) %% 5 == 0]
            ) +
            scale_y_continuous(labels = scales::comma)
    } else if (args$plot_type == "MAPQ-by-flag") {
        plot <- ggplot2::ggplot(
            data, aes(x = as.factor(MAPQ), y = tally, fill = as.factor(flag))
        ) +
            geom_bar(stat = "identity", position = "dodge") +
            labs(
                title = "MAPQ distribution by flag",
                x = "MAPQ",
                y = "tally",
                fill = "flag"
            ) +
            theme_minimal() +
            scale_x_discrete(
                breaks = function(x) x[as.numeric(x) %% 5 == 0]
            ) +
            scale_y_continuous(labels = scales::comma)
    } else if (args$plot_type == "flag") {
        plot <- ggplot2::ggplot(
            data, aes(x = as.factor(flag), y = tally)
        ) +
            geom_bar(stat = "identity") +
            labs(title = "Flag distribution", x = "flag", y = "tally") +
            theme_minimal() +
            scale_y_continuous(labels = scales::comma)
    }
    
    #  Apply zoom if necessary
    if (zoom_pct < 100) {
        plot <- apply_zoom(plot, args$zoom_pct)
    }
    
    #  Save the plot to a file, changing outfile name if zoom is applied
    ggplot2::ggsave(
        if (zoom_pct < 100) {
            gsub(".pdf$", paste0(".", zoom_pct, ".pdf"), args$outfile)
        } else {
            args$outfile
        },
        plot,
        height = args$height,
        width = args$width
    )
    
    #  Generate per-chromosome plots if requested
    if (args$per_chr) {
        if (is.null(args$chr_list) || length(args$chr_list) == 0) {
            chromosomes <- unique(data$chr)
        } else {
            chromosomes <- args$chr_list
        }
        
        for (chromosome in chromosomes) {
            chr_data <- data %>% filter(chr == chromosome)
            
            if (args$plot_type == "MAPQ") {
                chr_plot <- ggplot2::ggplot(
                    chr_data, aes(x = as.factor(MAPQ), y = tally)
                ) +
                    geom_bar(stat = "identity") +
                    labs(
                        title = paste(
                            "MAPQ distribution for chromosome", chromosome
                        ),
                        x = "MAPQ",
                        y = "tally"
                    ) +
                    theme_minimal() +
                    scale_x_discrete(
                        breaks = function(x) x[as.numeric(x) %% 5 == 0]
                    ) +
                    scale_y_continuous(labels = scales::comma)
            } else if (args$plot_type == "MAPQ-by-flag") {
                chr_plot <- ggplot2::ggplot(
                    chr_data, aes(
                        x = as.factor(MAPQ), y = tally, fill = as.factor(flag)
                    )
                ) +
                    geom_bar(stat = "identity", position = "dodge") +
                    labs(
                        title = paste(
                            "MAPQ distribution by flag for chromosome",
                            chromosome
                        ),
                        x = "MAPQ",
                        y = "tally",
                        fill = "flag"
                    ) +
                    theme_minimal() +
                    scale_x_discrete(
                        breaks = function(x) x[as.numeric(x) %% 5 == 0]
                    ) +
                    scale_y_continuous(labels = scales::comma)
            } else if (args$plot_type == "flags") {
                chr_plot <- ggplot2::ggplot(
                    chr_data, aes(x = as.factor(flag), y = tally)
                ) +
                    geom_bar(stat = "identity") +
                    labs(
                        title = paste(
                            "Flag distribution for chromosome", chromosome
                        ),
                        x = "Flag",
                        y = "tally"
                    ) +
                    theme_minimal() +
                    scale_y_continuous(labels = scales::comma)
            }
            
            #  Apply zoom if necessary
            if (args$zoom_pct < 100) {
                chr_plot <- apply_zoom(chr_plot, args$zoom_pct)
            }
            
            #  Determine outfile name
            if (args$zoom_pct < 100) {
                chr_outfile <- gsub(
                    ".pdf$",
                    paste0(".", chromosome, ".", args$zoom_pct, ".pdf"),
                    args$outfile
                )
            } else {
                chr_outfile <- gsub(
                    ".pdf$", paste0(".", chromosome, ".pdf"), args$outfile
                )
            }
            
            #  Save the plot to a file
            ggplot2::ggsave(
                chr_outfile,
                chr_plot,
                height = args$height,
                width = args$width
            )
        }
    }
}

#  Run the main function
main()
