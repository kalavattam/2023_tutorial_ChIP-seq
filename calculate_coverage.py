import argparse
import gzip
import sys

from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor


def parse_bed_file(bed_file_path):
    """
    Read a BED file and return a list of tuples containing (chromosome, start,
    end, fragment_length). Automatically detects if the file is gzip-
    compressed.
    """
    bed_entries = []
    open_func = gzip.open if bed_file_path.endswith('.gz') else open

    with open_func(bed_file_path, 'rt') as file:
        for line in file:
            fields = line.strip().split()
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            frag_length = int(fields[3])
            bed_entries.append((chrom, start, end, frag_length))

    return bed_entries


def calculate_coverage_for_chromosome(args):
    """
    Calculate the binned coverage for a set of BED entries, either 'normalized'
    or 'traditional'.
    """
    bed_entries, total_reads, bin_size, normalize = args
    coverage = defaultdict(float)
    for chrom, start, end, frag_length in bed_entries:
        if normalize:
            contribution = 1 / frag_length
        else:
            contribution = 1
        for pos in range(start, end + 1):
            coverage[(chrom, pos)] += contribution

    if normalize:
        #  Normalize the coverage by the number of observations
        for key in coverage:
            coverage[key] /= total_reads

    #  Bin the coverage
    binned_coverage = defaultdict(float)
    for (chrom, pos), value in coverage.items():
        bin_start = (pos // bin_size) * bin_size
        binned_coverage[(chrom, bin_start)] += value / bin_size

    return binned_coverage


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
        description='Calculate binned normalized coverage from BED file.'
    )
    parser.add_argument(
        '-i',
        '--input_file',
        required=True,
        help='Path to the input BED file.'
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
        type=int,
        default=30,
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
            'al., Sci Rep 2023.'
        )
    )

    #  Display help and exit if no arguments are provided
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    #  Parse and process BED file
    bed_entries = parse_bed_file(args.input_file)
    total_reads = len(bed_entries)
    entries_by_chromosome = defaultdict(list)
    for entry in bed_entries:
        entries_by_chromosome[entry[0]].append(entry)

    #  Prepare data for parallel processing
    task_data = [
        (entries, total_reads, args.bin_size, args.normalize)
        for entries in entries_by_chromosome.values()
    ]

    #  Execute parallel processing
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        results = executor.map(calculate_coverage_for_chromosome, task_data)

    #  Combine results from all processes
    combined_coverage = defaultdict(float)
    for result in results:
        for key, value in result.items():
            combined_coverage[key] += value

    #  Write output
    write_bedgraph(combined_coverage, args.output_file, args.bin_size)


if __name__ == "__main__":
    main()
