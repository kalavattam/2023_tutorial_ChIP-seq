#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright 2024 by Kris Alavattam
# Email: kalavattam@gmail.com
# Usage: python tally_alignments.py \
#            -f <bam_file> \
#            -t <no_threads> \
#            -q <minimum_MAPQ> \
#            -i <chromosomes_to_include>
# Description: This script tallies alignments in a BAM file based on specified
#              mapping quality (MAPQ values). Users can optionally specify
#              chromosomes to include; if none are specified, all chromosomes
#              are included. The script supports parallel processing using
#              multiple threads.

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
    #  Parse arguments
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
        nargs='*',
        default=None,
        help='List of chromosomes to include in tallying (default: all)'
    )
    
    args = parser.parse_args()
    
    #  Validate arguments 'threads' and 'mapq'
    if args.threads < 1:
        parser.error('Argument "threads" must be an int > 0.')
    if args.mapq < 0:
        parser.error('Argument "mapq" must be an int ≥ 0.')

    included_chromosomes = set(args.include) if args.include else None
    result = tally_alignments(
        args.file, args.threads, args.mapq, included_chromosomes
    )
    print(f"{result}")


if __name__ == "__main__":
    main()
