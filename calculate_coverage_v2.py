import argparse
import pysam
import sys

from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor


def parse_bam_file(bam_file_path):
    """
    Read a BAM file and return a dictionary of lists containing (chromosome,
    start, end). Entries are grouped by chromosome. Only considers reads with
    FLAGs 99 and 163.
    """
    bam_entries = defaultdict(list)
    bam_entries = defaultdict(list)
    with pysam.AlignmentFile(bam_file_path, 'rb') as bam_file:
        for read in bam_file.fetch():
            if read.flag in {99, 163}:
                chrom = bam_file.get_reference_name(read.reference_id)
                start = read.reference_start
                end = start + read.template_length
                frag_length = read.template_length
                bam_entries[chrom].append((start, end, frag_length))

    return bam_entries


def calculate_coverage_for_chromosome(
    chromosome, bam_entries, total_reads, bin_size, normalize
):
    """
    Calculate the binned coverage for a set of BED entries, either 'normalized'
    or 'traditional'.
    """
    if bin_size <= 0:
        raise ValueError("bin_size must be greater than 0")

    coverage = defaultdict(float)
    for start, end, frag_length in bam_entries:
        contribution = 1 / frag_length if normalize else 1
        for pos in range(start, end + 1):
            coverage[(chromosome, pos)] += contribution

    if normalize:
        for key in coverage:
            coverage[key] /= total_reads

    binned_coverage = defaultdict(float)
    for (chrom, pos), value in coverage.items():
        bin_start = (pos // bin_size) * bin_size
        binned_coverage[(chrom, bin_start)] += value / bin_size

    return binned_coverage


def calculate_coverage_task(data):
    return calculate_coverage_for_chromosome(*data)


def write_bedgraph(coverage, output_file, bin_size):
    """
    Write the normalized and binned coverage data to a BEDGRAPH file.
    """
    with open(output_file, 'w') as outfile:
        for (chrom, bin_start), value in sorted(
            coverage.items(), key=lambda x: (x[0][0], x[0][1])
        ):
            output_line = (
                f"{chrom}\t{bin_start}\t{bin_start + bin_size}\t{value:.24f}\n"
            )
            outfile.write(output_line)


def main():
    #  Parse arguments
    parser = argparse.ArgumentParser(
        description=(
            'Calculate binned coverage, normalized or not, from BAM file.'
        )
    )
    parser.add_argument(
        '-i',
        '--input_file',
        required=True,
        help='Path to the input BAM file.'
    )
    parser.add_argument(
        '-o',
        '--output_file',
        required=True,
        help='Path to the output BEDGRAPH file.'
    )
    parser.add_argument(
        '-b',
        '--bin_size',
        type=int, default=30,
        help='Bin size for coverage calculation in base pairs.'
    )
    parser.add_argument(
        '-t',
        '--threads',
        type=int,
        default=1,
        help='Number of threads for parallel processing.'
    )
    parser.add_argument(
        '-n',
        '--normalize',
        action='store_true',
        help=(
            'Normalize coverage by fragment length and total reads, ' +
            'generating "normalized coverage" as described in Dickson et ' +
            'al., Sci Rep 2023. If not specified, then "typical/raw ' +
            'coverage" is calculated.'
        )
    )

    #  Display help and exit if no arguments are provided
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    print(f"Input file: {args.input_file}")
    print(f"Output file: {args.output_file}")
    print(f"Bin size: {args.bin_size}")
    print(f"Threads: {args.threads}")
    print(f"Normalize: {args.normalize}")

    if args.bin_size <= 0:
        raise ValueError("bin_size must be greater than 0")

    #  Parse and process BAM file
    bam_entries = parse_bam_file(args.input_file)
    total_reads = len(bam_entries)
    combined_coverage = defaultdict(float)

    #  Prepare data for parallel processing
    task_data = [
        (chrom, entries, total_reads, args.bin_size, args.normalize)
        for chrom, entries in bam_entries.items()
    ]

    #  Execute parallel processing
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        results = executor.map(calculate_coverage_task, task_data)

    #  Combine results from all processes
    for result in results:
        for key, value in result.items():
            combined_coverage[key] += value

    #  Write output to BEDGRAPH file
    write_bedgraph(combined_coverage, args.output_file, args.bin_size)


if __name__ == "__main__":
    main()
