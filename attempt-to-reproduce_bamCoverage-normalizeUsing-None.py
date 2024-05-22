import argparse
import pysam
import sys
from collections import defaultdict


def calculate_coverage(
    bam_file_path, effective_genome_size, bin_size, min_mapq, region
):
    """
    Calculate the normalized and binned coverage from a BAM file using the
    following formula: `(total number of mapped reads * average fragment
    length) / effective genome size`.
    """
    bamfile = pysam.AlignmentFile(bam_file_path, "rb")
    coverage = defaultdict(float)
    total_reads = 0
    total_fragment_length = 0
    num_fragments = 0

    #  Define the region to fetch
    if region.lower() == 'all':
        region = None  # `None` will fetch the whole file

    for read in bamfile.fetch(region=region):
        if (
            not read.is_unmapped and
            not read.is_secondary and
            not read.is_supplementary and
            read.mapping_quality >= min_mapq
        ):
            total_reads += 1
            if read.is_proper_pair and read.is_read1:
                fragment_length = read.template_length
                if fragment_length > 0:
                    total_fragment_length += fragment_length
                    num_fragments += 1
                    start = read.reference_start
                    end = read.reference_end
                    chrom = bamfile.get_reference_name(read.reference_id)
                    for pos in range(start, end):
                        bin_start = (pos // bin_size) * bin_size
                        coverage[(chrom, bin_start)] += 1

    average_fragment_length = (
        total_fragment_length / num_fragments if num_fragments else 0
    )
    scaling_factor = (
        (total_reads * average_fragment_length) / effective_genome_size
    )

    for key in coverage:
        coverage[key] /= scaling_factor

    return coverage


def write_bedgraph(coverage, output_file, bin_size):
    """
    Write the normalized and binned coverage data to a BEDGRAPH file.
    """
    with open(output_file, 'w') as outfile:
        for (chrom, bin_start), value in sorted(
            coverage.items(), key=lambda x: (x[0], x[1])
        ):
            outfile.write(
                f"{chrom}\t{bin_start}\t{bin_start + bin_size}\t{value:.6f}\n"
            )


def is_bam_indexed(bam_file_path):
    """
    Check if the BAM file is indexed by trying to fetch from it. This function
    assumes the BAM file is coordinate-sorted.
    """
    try:
        bamfile = pysam.AlignmentFile(bam_file_path, "rb")
        try:
            bamfile.fetch()
            return True
        except ValueError:
            return False
    finally:
        bamfile.close()


def main():
    #  Parse arguments
    parser = argparse.ArgumentParser(description=(
        "Calculate normalized and binned coverage from a BAM file using the "
        "formula for when deepTools `bamCoverage` is called with "
        "`--normalizeUsing 'None'`."
    ))
    parser.add_argument(
        '-i',
        '--input_bam',
        required=True,
        help="Input BAM file path (coordinate-sorted and indexed)."
    )
    parser.add_argument(
        '-o',
        '--output_bedgraph',
        required=True,
        help="Output BEDGRAPH file path."
    )
    parser.add_argument(
        '-g',
        '--genome_size',
        type=int,
        required=True,
        help="Effective genome size."
    )
    parser.add_argument(
        '-b',
        '--bin_size',
        type=int,
        default=30,
        help="Bin size for coverage calculation in base pairs."
    )
    parser.add_argument(
        '-m',
        '--min_mapq',
        type=int,
        default=1,
        help="Minimum mapping quality for reads to be considered."
    )
    parser.add_argument(
        '-r',
        '--region',
        default='all',
        help=(
            "Region of the BAM to analyze (format: 'chr:start-end' or 'chr' "
            "or 'all')."
        )
    )

    #  Display help and exit if no arguments were provided
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    #  Check if BAM file is indexed
    if not is_bam_indexed(args.input_bam):
        print(
            "Error: BAM file must be indexed. Please create an index (e.g., "
            "using `samtools index`)."
        )
        sys.exit(1)

    coverage = calculate_coverage(
        args.input_bam, args.genome_size, args.bin_size, args.min_mapq,
        args.region
    )
    write_bedgraph(
        coverage, args.output_bedgraph, args.bin_size
    )


if __name__ == "__main__":
    main()
