#!/usr/bin/env Rscript

#  generate_TxDb-from-gff3.R
#  KA

#  Initialize functions -------------------------------------------------------
quietly_load_library <- function(library_name) {
    #  Function to load a library quietly
    suppressPackageStartupMessages(library(
        library_name, character.only = TRUE
    ))
}


get_gff3_connection <- function(gff3_file) {
    #  Function to check if the file is gzipped and return the appropriate file
    #+ connection or path
    if (grepl("\\.gz$", gff3_file)) {
        return(gzfile(gff3_file))
    } else {
        return(gff3_file)
    }
}


parse_args <- function() {
    # parse_args function
    #
    # This function defines and parses the command-line arguments using the
    # argparse library.
    #
    # It sets up the expected input arguments for the GFF3 file path and the
    # output database file path. It returns a list of arguments that can be
    # used in the script.
    #
    # Returns:
    #   A list containing the values of the parsed command-line arguments.
    parser <- ArgumentParser(
        description = "
            Generate TxDb from GFF3 file and save as a database. A TxDb object
            is an SQLite-based annotation data package that stores
            transcript-centric annotations for a genome. It includes
            information about genes, transcripts, exons, cds, and much more.
        "
    )

    #  Define arguments
    parser$add_argument(
        "-g", "--gff3",
        type = "character",
        required = TRUE,
        help = "Path to the GFF3 file"
    )
    parser$add_argument(
        "-d", "--dbfile",
        type = "character",
        required = TRUE,
        help = "
            Full path and name for the database output file (with .db
            extension)
        "
    )

    #  Parse the arguments
    args <- parser$parse_args()
    return(args)
}


main <- function() {
    #  Get command line arguments
    arguments <- parse_args()

    gff3_file <- arguments$gff3
    db_filepath <- arguments$dbfile
    
    #  Check if the GFF3 file is gzipped and get the appropriate file
    #+ connection or path
    gff3_connection <- get_gff3_connection(gff3_file)

    #  Generate TxDb from gff3 file
    TxDb <- GenomicFeatures::makeTxDbFromGFF(
        gff3_file,
        format = "gff3"
    )

    #  Save the TxDb to the specified database file
    AnnotationDbi::saveDb(TxDb, db_filepath)

    #  Print success message
    cat("Database saved to:", db_filepath, "\n")
}


#  Load required libraries ----------------------------------------------------
quietly_load_library("argparse")
quietly_load_library("GenomicFeatures")
quietly_load_library("AnnotationDbi")
quietly_load_library("stringr")


#  Run the main function if the script is not being sourced -------------------
if (!interactive() && sys.nframe() == 0) {
    main()
}
