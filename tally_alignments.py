#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright 2024 by Kris Alavattam
# Email: kalavattam@gmail.com
# Usage: python tally_alignments.py \
#            -f <bam_file> \
#            -t <number_of_threads> \
#            -q <minimum_MAPQ> \
#            -i <chromosomes_to_include>
# Description: This script tallies alignments in a BAM file based on specified
# mapping quality (MAPQ values). Users can optionally specify chromosomes to
# include; if none are specified, all chromosomes are included. The script
# supports parallel processing using multiple threads.
#
# Distributed under terms of the MIT license.

import argparse
import pysam
import sys


def is_included_chromosome(read, included_chromosomes=None):
    """
    Check if a read should be included based on chromosome name.

    Args:
        read (pysam.AlignedSegment): A single read from a BAM file.
        included_chromosomes (set, optional): Set of chromosomes to include.
                                              If None, includes all
                                              chromosomes.

    Returns:
        bool: True if the read's chromosome is not in the excluded list, False
              otherwise.
    """
    if included_chromosomes is None:
        return True  # Include all chromosomes if no specific list is provided
    else:
        return read.reference_name in included_chromosomes


def split_chromosomes(chrom_str):
    """
    Split a string into a list of chromosome names, handling space-separated,
    comma-separated, and comma-and-space-separated strings.

    Args:
        chrom_str (str): A string containing chromosome names such as
                         'I II III', 'I,II,III', or 'I, II, III'.

    Returns:
        list: A list of strings, where each string is a chromosome name
              extracted from the input string.
    """
    if ',' in chrom_str:
        return chrom_str.split(',')
    elif ', ' in chrom_str:
        return chrom_str.split(', ')
    else:
        return chrom_str.split()


def tally_alignments(file, threads, mapq, included_chromosomes=None):
    """
    Tally alignments based on MAPQ score and specified chromosome inclusions.
    
    Args:
        file (str): Path to the BAM file.
        threads (int): Number of threads for processing.
        mapq (int): Minimum MAPQ score to consider.
        included_chromosomes (set, optional): Specific chromosomes to include.
    
    Returns:
        int: Number of alignments meeting the criteria.
    """
    try:
        bamfile = pysam.AlignmentFile(file, "rb", threads=threads)
        count = 0
        for read in bamfile.fetch():
            if (
                read.mapping_quality >= mapq and
                is_included_chromosome(read, included_chromosomes)
            ):
                count += 1
        bamfile.close()
        return count
    except Exception as e:
        print(f"Error processing the BAM file: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    """
    Execute the primary control flow of the script, handling command line
    parsing and facilitating the tallying of alignments based on user-provided
    parameters. This function sets up logging, validates input, and calls a
    function to tally alignments if all parameters are correctly specified.

    Command-line Arguments:
        -f/--file (str): Path to the BAM file. [required]
        -t/--threads (int): Number of processing threads. [default: 1]
        -q/--mapq (int): Minimum MAPQ score to consider. [default: 1]
        -i/--include (list): Chromosomes to include in the tally, separated by
                             spaces, commas, or commas-and-one-space.
                             [default: None]

    Returns:
        None, but outputs the number of alignments to standard output.

    Exits:
        1 on error with a message to standard error,
        0 on successful completion.
    """
    #  Define arguments
    parser = argparse.ArgumentParser(
        description='Tally alignments within BAM file.'
    )
    parser.add_argument(
        '-f',
        '--file',
        required=True,
        type=str,
        help='Path to the BAM file [required]'
    )
    parser.add_argument(
        '-t',
        '--threads',
        type=int,
        default=1,
        help='Number of threads for parallel processing (default: 1) [int ≥ 1]'
    )
    parser.add_argument(
        '-q',
        '--mapq',
        type=int,
        default=1,
        help=(
            'Filter BAM to retain alignments with this MAPQ score or greater '
            '(default: 1) [int ≥ 0]'
        )
    )
    parser.add_argument(
        '-i',
        '--include',
        type=split_chromosomes,
        default=None,
        help='List of chromosomes to include in tallying (default: all)'
    )
    
    #  Display help and exit if no arguments were provided
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    #  Parse arguments
    args = parser.parse_args()
    
    #  Validate arguments 'threads' and 'mapq'
    if args.threads < 1:
        parser.error('Argument "threads" must be an int > 0.')
    if args.mapq < 0:
        parser.error('Argument "mapq" must be an int ≥ 0.')

    #  Run tally_alignments, returning results
    result = tally_alignments(args.file, args.threads, args.mapq, args.include)
    print(f"{result}")


if __name__ == "__main__":
    main()
