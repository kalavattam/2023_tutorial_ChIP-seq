#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright 2024 by Kris Alavattam
# Email: kalavattam@gmail.com
# Usage: python calculate_scaling_factor.py \
#            --ip_main <no_ip_main_alignments> \
#            --ip_spike_in <no_ip_spike_in_alignments> \
#            --input_main <no_input_main_alignments> \
#            --input_spike_in <no_input_spike_in_alignments>
# Description: This script calculates a spike-in-derived scaling factor for
# ChIP-seq datasets. It requires counts of "main" and spike-in alignments for
# corresponding IP and input samples. The script calculates the ratio of
# spike-in to main alignments for IP and input samples, and then computes a
# final scaling factor by dividing the input sample ratio by the IP sample
# ratio.
#
# Distributed under terms of the MIT license.

import argparse
import sys


def calculate_scaling_factor(ip_main, ip_spike_in, input_main, input_spike_in):
    """
    Calculate the spike-in-derived scaling factor based on alignment counts.

    This function computes two ratios: the ratio of spike-in to "main"
    alignments for the IP sample and the same ratio for the input sample. It
    then returns the scaling factor derived by dividing the input sample ratio
    by the IP sample ratio.

    Args:
        ip_main (int): The number of "main" alignments in the IP sample.
        ip_spike_in (int): The number of spike-in alignments in the IP sample.
        input_main (int): The number of "main" alignments in the input sample.
        input_spike_in (int): The number of spike-in alignments in the input
                              sample.

    Returns:
        float: The scaling factor calculated as the ratio of input sample
               ratio to IP sample ratio.
    
    Raises:
        ValueError: If any of the counts are negative.
        TypeError: If any of the counts are not integers.
        ZeroDivisionError: If any of the main alignment counts are zero.
    """
    #  Validate that all inputs are integers
    for count in [ip_main, ip_spike_in, input_main, input_spike_in]:
        if not isinstance(count, int):
            raise TypeError(f'Expected integer, got {type(count).__name__}')

    #  Validate that all inputs are non-negative
    if any(
        count < 0 for count in [
            ip_main, ip_spike_in, input_main, input_spike_in
        ]
    ):
        raise ValueError('Alignment counts cannot be negative.')

    #  Check for division by zero
    if ip_main == 0 or input_main == 0:
        raise ZeroDivisionError(
            'Main alignment counts cannot be zero to avoid division by zero.'
        )

    ratio_ip = ip_spike_in / ip_main
    ratio_input = input_spike_in / input_main
    
    #  Calculate the scaling factor
    return ratio_input / ratio_ip


def main():
    """
    Execute the primary control flow of the script.

    This function facilitates the calculation of a spike-in-derived scaling
    factor for ChIP-seq samples. It parses command-line arguments for the
    number of alignments from both IP and input samples, calculates ratios of
    spike-in to main alignments for both sample types, and then computes a
    final scaling factor by dividing the input sample ratio by the IP sample
    ratio.

    Command-line Arguments:
        -ip_m/--ip_main (int): The number of "main" alignments for the ChIP-seq
                               IP sample.
        -ip_s/--ip_spike_in (int): The number of spike-in alignments for the
                                   ChIP-seq IP sample.
        -in_m/--input_main (int): The number of "main" alignments for the
                                  corresponding ChIP-seq input sample.
        -in_s/--input_spike_in (int): The number of spike-in alignments for the
                                      corresponding ChIP-seq input sample.

    The function will terminate early and print an error message if any of the
    main alignment counts are zero, preventing division by zero errors.

    Returns:
        None, but outputs the final scaling factor to standard output, which is
        the ratio of the input sample spike-in to main alignments ratio divided
        by the IP sample spike-in to main alignments ratio.
    """
    #  Define arguments
    parser = argparse.ArgumentParser(description=(
        'Calculate spike-in-derived scaling factor for a ChIP-seq sample with '
        'IP and input data.'
    ))
    parser.add_argument(
        '-ip_m',
        '--ip_main',
        type=int,
        help='Number of "main" alignments for the ChIP-seq IP sample.'
    )
    parser.add_argument(
        '-ip_s',
        '--ip_spike_in',
        type=int,
        help='Number of spike-in alignments for the ChIP-seq IP sample.'
    )
    parser.add_argument(
        '-in_m',
        '--input_main',
        type=int,
        help=(
            'Number of "main" alignments for the corresponding ChIP-seq input '
            'sample.'
        )
    )
    parser.add_argument(
        '-in_s',
        '--input_spike_in',
        type=int,
        help=(
            'Number of spike-in alignments for the corresponding ChIP-seq '
            'input sample.'
        )
    )

    #  Display help and exit if no arguments were provided
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    #  Calculate the scaling factor, handling exceptions as necessary
    try:
        result = calculate_scaling_factor(
            args.ip_main,
            args.ip_spike_in,
            args.input_main,
            args.input_spike_in
        )
        print(f'{result}')
    except ValueError as e:
        print(f'Error: {e}', file=sys.stderr)
        sys.exit(1)
    except TypeError as e:
        print(f'Error: {e}', file=sys.stderr)
        sys.exit(1)
    except ZeroDivisionError as e:
        print(f'Error: {e}', file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
