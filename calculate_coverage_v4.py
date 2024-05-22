import argparse
import numpy as np
import pysam
import sys

from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor


def get_fragment_from_read(
    read, extend_reads, center_reads, default_fragment_length,
    max_paired_fragment_length
):
    """
    Get read start and end position of a read in the style of deepTools
    `get_fragment_from_read`.
    """
    if not extend_reads and default_fragment_length == 'read length':
        return read.get_blocks()
    
    if not extend_reads:
        return [(read.reference_start, read.reference_end)]
    
    if (
        read.is_paired and read.is_proper_pair
    ) and (
        abs(read.template_length) < max_paired_fragment_length
    ):
        if read.is_reverse:
            fragment_start = read.next_reference_start
            fragment_end = read.reference_end
        else:
            fragment_start = read.reference_start
            fragment_end = read.reference_start + abs(read.template_length)
    else:
        if read.is_reverse:
            fragment_start = read.reference_end - default_fragment_length
            fragment_end = read.reference_end
        else:
            fragment_start = read.reference_start
            fragment_end = read.reference_start + default_fragment_length

    if center_reads:
        fragment_center = (
            fragment_end - (fragment_end - fragment_start) / 2
        )
        fragment_start = int(
            fragment_center - read.infer_query_length(always=False) / 2
        )
        fragment_end = (
            fragment_start + read.infer_query_length(always=False)
        )

    assert fragment_start < fragment_end, f"Fragment start greater than fragment end for read {read.query_name}"
    
    return [(fragment_start, fragment_end)]


def parse_bam_file(
    bam_file_path, sam_flag_include, min_mapq, extend_reads, center_reads,
    default_fragment_length, max_paired_fragment_length
):
    bam_entries = defaultdict(list)
    with pysam.AlignmentFile(bam_file_path, 'rb') as bam_file:
        for read in bam_file.fetch():
            if (
                read.flag & sam_flag_include
            ) and (
                read.mapping_quality >= min_mapq
            ):
                chrom = bam_file.get_reference_name(read.reference_id)
                fragments = get_fragment_from_read(
                    read,
                    extend_reads,
                    center_reads,
                    default_fragment_length,
                    max_paired_fragment_length
                )
                for fragment_start, fragment_end in fragments:
                    frag_length = fragment_end - fragment_start
                    bam_entries[chrom].append(
                        (fragment_start, fragment_end, frag_length)
                    )

    return bam_entries


def calculate_coverage_for_chromosome(
    chromosome, bam_entries, total_reads, bin_size, normalize
):
    if bin_size <= 0:
        raise ValueError("bin_size must be greater than 0")

    max_end = max(end for _, end, _ in bam_entries)
    n_reg_bins = (max_end // bin_size) + 1
    coverages = np.zeros(n_reg_bins, dtype=float)

    last_e_idx = None  # Used for managing overlaps
    prev_end = 0  # Track the end of the last fragment processed

    for start, end, frag_length in bam_entries:
        contribution = 1 / frag_length if normalize else 1

        fragment_start = max(start, 0)
        fragment_end = min(end, max_end)

        s_idx = (fragment_start // bin_size)
        e_idx = (fragment_end // bin_size) + 1

        print(
            f"Processing fragment from {start} to {end}, indices {s_idx} " +
            f"to {e_idx}, contribution {contribution}"
        )

        if fragment_start > prev_end:
            last_e_idx = None  # Reset last_e_idx as there's no overlap

        if last_e_idx is not None and s_idx < last_e_idx:
            s_idx = last_e_idx
            print(f"Adjusting start index due to overlap: new s_idx {s_idx}")

        if s_idx >= e_idx:
            print(
                f"Skipping fragment with adjusted start index {s_idx} >= " +
                f"end index {e_idx}"
            )
            continue

        coverages[s_idx:e_idx] += contribution
        print(
            f"Incremented coverage from indices {s_idx} to {e_idx} by " +
            f"{contribution}"
        )
        
        last_e_idx = e_idx
        prev_end = fragment_end  # Update end of last processed fragment

    binned_coverage = {
        (chromosome, i * bin_size): val
        for i, val in enumerate(coverages)
        if val > 0
    }

    if normalize:
        for key in binned_coverage:
            binned_coverage[key] /= total_reads

    return binned_coverage


# def calculate_coverage_for_chromosome(
#     chromosome, bam_entries, total_reads, bin_size, normalize
# ):
#     if bin_size <= 0:
#         raise ValueError("bin_size must be greater than 0")
#
#     #  Initialize the coverages array based on the maximum end position
#     #  provided in bam_entries
#     max_end = max(end for _, end, _ in bam_entries)
#     n_reg_bins = (max_end // bin_size) + 1
#     coverages = np.zeros(n_reg_bins, dtype=float)
#
#     last_e_idx = None  # Used for managing overlaps
#
#     for start, end, frag_length in bam_entries:
#         contribution = 1 / frag_length if normalize else 1
#
#         #  Adjust fragment start and end to fit within bins
#         fragment_start = max(start, 0)  # Assuming 0 is the chr/region start
#         fragment_end = min(end, max_end)
#
#         s_idx = (fragment_start // bin_size)
#         e_idx = (fragment_end // bin_size) + 1
#
#         #  Ensure no double counting for overlapping fragments
#         if last_e_idx is not None and s_idx < last_e_idx:
#             s_idx = last_e_idx
#         if s_idx >= e_idx:
#             continue
#
#         #  Increment coverage
#         coverages[s_idx:e_idx] += contribution
#         last_e_idx = e_idx
#
#     #  Convert coverages to a dictionary for compatibility with your existing
#     #  format
#     binned_coverage = {
#         (chromosome, i * bin_size):
#         val for i, val in enumerate(coverages) if val > 0
#     }
#
#     #  Normalize if necessary
#     if normalize:
#         for key in binned_coverage:
#             binned_coverage[key] /= total_reads
#
#     return binned_coverage


# def calculate_coverage_for_chromosome(
#     chromosome, bam_entries, total_reads, bin_size, normalize
# ):
#     if bin_size <= 0:
#         raise ValueError("bin_size must be greater than 0")
#
#     coverage = defaultdict(float)
#     for start, end, frag_length in bam_entries:
#         contribution = 1 / frag_length if normalize else 1
#         for pos in range(start, end + 1):
#             coverage[(chromosome, pos)] += contribution
#
#     if normalize:
#         for key in coverage:
#             coverage[key] /= total_reads
#
#     binned_coverage = defaultdict(float)
#     for (chrom, pos), value in coverage.items():
#         bin_start = (pos // bin_size) * bin_size
#         binned_coverage[(chrom, bin_start)] += value / bin_size
#
#     return binned_coverage


def calculate_coverage_task(data):
    return calculate_coverage_for_chromosome(*data)


def write_bedgraph(coverage, output_file, bin_size):
    with open(output_file, 'w') as outfile:
        for (chrom, bin_start), value in sorted(
            coverage.items(), key=lambda x: (x[0][0], x[0][1])
        ):
            output_line = (
                f"{chrom}\t{bin_start}\t{bin_start + bin_size}\t{value:.24f}\n"
            )
            outfile.write(output_line)


def main():
    parser = argparse.ArgumentParser(
        description=(
            'Calculate binned coverage, normalized or not, from BAM file.'
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '-i', '--input_file',
        required=True,
        help='Path to the input BAM file.'
    )
    parser.add_argument(
        '-o', '--output_file',
        required=True,
        help='Path to the output BEDGRAPH file.'
    )
    parser.add_argument(
        '-b', '--bin_size',
        type=int,
        default=30,
        help='Bin size for coverage calculation in base pairs.'
    )
    parser.add_argument(
        '-t', '--threads',
        type=int,
        default=1,
        help='Number of threads for parallel processing.'
    )
    parser.add_argument(
        '-n', '--normalize',
        action='store_true',
        help=(
            'Normalize coverage by fragment length and total reads, ' +
            'generating "normalized coverage" as described in Dickson et ' +
            'al., Sci Rep 2023. If not specified, then "typical/raw ' +
            'coverage" is calculated.'
        )
    )
    parser.add_argument(
        '-s', '--SAM_flag_include',
        type=int,
        default=0,
        help='SAM flag to include reads.'
    )
    parser.add_argument(
        '-m', '--min_MAPQ',
        type=int,
        default=0,
        help='Minimum mapping quality (MAPQ).'
    )
    parser.add_argument(
        '-e', '--extend_reads',
        action='store_true',
        help='Extend reads to fragment length (in the manner of deepTools).'
    )
    parser.add_argument(
        '-c', '--center_reads',
        action='store_true',
        help=(
            'Center reads on fragment length (in the manner of ' +
            'deepTools).'
        )
    )
    parser.add_argument(
        '-d', '--default_fragment_length',
        type=int,
        default=200,
        help='Default fragment length for extending reads.'
    )
    parser.add_argument(
        '-l', '--max_paired_fragment_length',
        type=int,
        default=1000,
        help='Maximum fragment length for considering "proper pairs."'
    )
    
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    
    args = parser.parse_args()

    print("#  Parameters --------------------------")
    print(f"                Input file: {args.input_file}")
    print(f"               Output file: {args.output_file}")
    print(f"                  Bin size: {args.bin_size}")
    print(f"                   Threads: {args.threads}")
    print(f"                 Normalize: {args.normalize}")
    print(f"          SAM flag include: {args.SAM_flag_include}")
    print(f"              Minimum MAPQ: {args.min_MAPQ}")
    print(f"              Extend reads: {args.extend_reads}")
    print(f"              Center reads: {args.center_reads}")
    print(f"   Default fragment length: {args.default_fragment_length}")
    print(f"Max paired fragment length: {args.max_paired_fragment_length}")

    if args.bin_size <= 0:
        raise ValueError("bin_size must be greater than 0")

    bam_entries = parse_bam_file(
        args.input_file,
        args.SAM_flag_include,
        args.min_MAPQ,
        args.extend_reads,
        args.center_reads,
        args.default_fragment_length,
        args.max_paired_fragment_length
    )

    total_reads = sum(len(entries) for entries in bam_entries.values())
    combined_coverage = defaultdict(float)

    task_data = [
        (chrom, entries, total_reads, args.bin_size, args.normalize)
        for chrom, entries in bam_entries.items()
    ]

    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        results = executor.map(calculate_coverage_task, task_data)

    for result in results:
        for key, value in result.items():
            combined_coverage[key] += value

    write_bedgraph(combined_coverage, args.output_file, args.bin_size)


if __name__ == "__main__":
    main()
