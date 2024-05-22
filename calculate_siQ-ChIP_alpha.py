#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright 2024 by Kris Alavattam
# Email: kalavattam@gmail.com
# Usage: python calculate_siQ-ChIP_alpha.py \
#            --mass_ip <mass_of_IP_sample> \
#            --mass_in <mass_of_input_sample> \
#            --volume_ip <volume_of_IP_sample> \
#            --volume_in <volume_of_input_sample> \
#            --depth_ip <sequencing_depth_of_IP_sample> \
#            --depth_in <sequencing_depth_of_input_sample> \
#            --length_ip <mean_fragment_length_of_IP_sample> \
#            --length_in <mean_fragment_length_of_input_sample>
# Description: This script calculates the siQ-ChIP alpha scaling factor for
# ChIP-seq datasets. It uses the provided mass, volume, sequencing depth, and
# fragment length of IP and input samples to compute the alpha scaling factor
# based on a mathematical formula per the details described here:
# - pubmed.ncbi.nlm.nih.gov/32994221/
# - pubmed.ncbi.nlm.nih.gov/37160995/ (particularly this)
#
# Distributed under terms of the MIT license.

import argparse
import sys


def calculate_alpha(
    mass_ip, mass_in,
    volume_ip, volume_in,
    depth_ip, depth_in,
    length_ip, length_in
):
    """
    Calculate a siQ-ChIP 'alpha' scaling factor using the provided values.
    
    Args:
        mass_ip (float): Mass of IP sample.
        mass_in (float): Mass of input sample.
        volume_ip (float): Volume of IP sample.
        volume_in (float): Volume of input sample.
        depth_ip (float): Sequencing depth of IP sample.
        depth_in (float): Sequencing depth of input sample.
        length_ip (float): Mean fragment length of IP sample.
        length_in (float): Mean fragment length of input sample.

    Returns:
        float: The calculated alpha scaling factor.
    """
    alpha = (
        (mass_ip / mass_in) *
        (volume_in / volume_ip) *
        (depth_in / depth_ip) *
        (length_in / length_ip)
    )
    return alpha


def main():
    """
    Execute the main functionality of the script.

    This function orchestrates the calculation of the siQ-ChIP alpha scaling
    factor for ChIP-seq datasets by parsing command-line arguments for the
    required experimental values of IP and input samples. It then calculates
    the alpha scaling factor based on these inputs and prints the result.

    Command-line Arguments:
        -mi/--mass_ip (float): Required. Specifies the mass of the IP sample.
        -mn/--mass_in (float): Required. Specifies the mass of the input
                               sample.
        -vi/--volume_ip (float): Required. Specifies the volume of the IP
                                 sample.
        -vn/--volume_in (float): Required. Specifies the volume of the input
                                 sample.
        -di/--depth_ip (float): Required. Specifies the sequencing depth of the
                                IP sample.
        -dn/--depth_in (float): Required. Specifies the sequencing depth of the
                                input sample.
        -li/--length_ip (float): Required. Specifies the mean fragment length
                                 of the IP sample.
        -ln/--length_in (float): Required. Specifies the mean fragment length
                                 of the input sample.

    Returns:
        None. Outputs the calculated alpha scaling factor to stdout.
    """
    parser = argparse.ArgumentParser(
        description='Calculate a siQ-ChIP alpha scaling factor.'
    )
    parser.add_argument(
        '-mi',
        '--mass_ip',
        type=float,
        required=True,
        help='Mass of IP sample'
    )
    parser.add_argument(
        '-mn',
        '--mass_in',
        type=float,
        required=True,
        help='Mass of input sample'
    )
    parser.add_argument(
        '-vi',
        '--volume_ip',
        type=float,
        required=True,
        help='Volume of IP sample'
    )
    parser.add_argument(
        '-vn',
        '--volume_in',
        type=float,
        required=True,
        help='Volume of input sample'
    )
    parser.add_argument(
        '-di',
        '--depth_ip',
        type=float,
        required=True,
        help='Sequencing depth of IP sample'
    )
    parser.add_argument(
        '-dn',
        '--depth_in',
        type=float,
        required=True,
        help='Sequencing depth of input sample'
    )
    parser.add_argument(
        '-li',
        '--length_ip',
        type=float,
        required=True,
        help='Mean fragment length of IP sample'
    )
    parser.add_argument(
        '-ln',
        '--length_in',
        type=float,
        required=True,
        help='Mean fragment length of input sample'
    )

    #  Display help and exit if no arguments were provided
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    alpha = calculate_alpha(
        args.mass_ip, args.mass_in, args.volume_ip, args.volume_in,
        args.depth_ip, args.depth_in, args.length_ip, args.length_in
    )
    print(f"{alpha}")


if __name__ == "__main__":
    main()
