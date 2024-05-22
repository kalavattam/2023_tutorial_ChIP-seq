
`#TBD.md`
<br />

<details>
<summary><b><font size="+1"><i>Table of contents</i></font></b></summary>
<br />
<!-- MarkdownTOC -->

1. [Testing methods to generate "typical coverage"](#testing-methods-to-generate-typical-coverage)
    1. [Bash code: Get situated, downsample a BAM file for tests](#bash-code-get-situated-downsample-a-bam-file-for-tests)
    1. [List and tally FLAGs in BAM files of interest](#list-and-tally-flags-in-bam-files-of-interest)
    1. [Coverage tests, 2024-0503-0506](#coverage-tests-2024-0503-0506)
1. [Assessments/comparisons of BEDGRAPH files generated from custom code and `bamCoverage` \(`coverage-tests_2024-0503-0506`\)](#assessmentscomparisons-of-bedgraph-files-generated-from-custom-code-and-bamcoverage-coverage-tests_2024-0503-0506)
    1. [Images](#images)
    1. [Breakdown of rows in images](#breakdown-of-rows-in-images)
    1. [Observations](#observations)
        1. [On the `--normalizeUsing` documentation issue \(i.e., regarding `None` and `RPGC`\)](#on-the---normalizeusing-documentation-issue-ie-regarding-none-and-rpgc)
1. [Snippets: How different normalizations are calculated in the deepTools codebase](#snippets-how-different-normalizations-are-calculated-in-the-deeptools-codebase)
1. [Reproducing what deepTools `bamCoverage` does via simplified versions of its methods](#reproducing-what-deeptools-bamcoverage-does-via-simplified-versions-of-its-methods)
    1. [Drafting a simplified version of `count_reads_in_region`](#drafting-a-simplified-version-of-count_reads_in_region)
        1. [The simplified version of `count_reads_in_region`](#the-simplified-version-of-count_reads_in_region)
        1. [Breaking down simplified version of `count_reads_in_region`'s parameter `step_size`](#breaking-down-simplified-version-of-count_reads_in_regions-parameter-step_size)
    1. [Drafting simplified versions of `is_proper_pair` and `get_fragment_from_read`](#drafting-simplified-versions-of-is_proper_pair-and-get_fragment_from_read)
        1. [The simplified version of `is_proper_pair`](#the-simplified-version-of-is_proper_pair)
        1. [The simplified version of `get_fragment_from_read`](#the-simplified-version-of-get_fragment_from_read)
        1. [Setting up tests for `is_proper_pair` and `get_fragment_from_read`](#setting-up-tests-for-is_proper_pair-and-get_fragment_from_read)
        1. [Additional breakdown of `get_fragment_from_read`](#additional-breakdown-of-get_fragment_from_read)
    1. [Drafting a simplified version of `get_read_and_fragment_length` \(renamed to `determine_fragment_and_read_lengths`\)](#drafting-a-simplified-version-of-get_read_and_fragment_length-renamed-to-determine_fragment_and_read_lengths)
        1. [The original function](#the-original-function)
        1. [The simplified version](#the-simplified-version)
        1. [On including multiprocessing in the simplified version](#on-including-multiprocessing-in-the-simplified-version)
        1. [On running integration tests for `process_reads` and `determine_fragment_and_read_lengths`](#on-running-integration-tests-for-process_reads-and-determine_fragment_and_read_lengths)
    1. [Refactoring `get_fragment_from_read` such that each mate in a pair of reads is extended to match the fragment length, etc.](#refactoring-get_fragment_from_read-such-that-each-mate-in-a-pair-of-reads-is-extended-to-match-the-fragment-length-etc)
1. [Testing the simplified functions](#testing-the-simplified-functions)
1. [Test `calculate_coverage_v2.py`](#test-calculate_coverage_v2py)

<!-- /MarkdownTOC -->
</details>
<br />

<a id="testing-methods-to-generate-typical-coverage"></a>
## Testing methods to generate "typical coverage"
<a id="bash-code-get-situated-downsample-a-bam-file-for-tests"></a>
### Bash code: Get situated, downsample a BAM file for tests
<details>
<summary><i>Bash code: Get situated, downsample a BAM file for tests</i></summary>

```bash
#!/bin/bash

#  Go to work directory
cd ~/repos/2023_tutorial_ChIP-seq || echo "cd'ing failed; check on this"

#  Source appropriate environment
conda activate deeptools_env

#  Downsample a BAM file for interactive testing
bam_init="IP_G1_Hho1_6336.sort-coord.bam"
bam_down="${bam_init%.bam}.ds.bam"

threads=4
samp_frac=0.1
samp_seed=24

if [[ ! -f "${bam_down}" ]]; then
    samtools view \
        -@ "${threads}" \
        --require-flags 2 \
        --subsample "${samp_frac}" \
        --subsample-seed "${samp_seed}" \
        "${bam_init}" \
        -o "${bam_down}"

    samtools sort \
        -@ "${threads}" \
        -n \
        -o "${bam_down%.bam}.qname.bam" \
        "${bam_down}"

    if [[ -f "${bam_down%.bam}.qname.bam" ]]; then
        mv -f \
            "${bam_down%.bam}.qname.bam" \
            "${bam_down}"
    fi

    samtools fixmate \
        -@ "${threads}" \
        -m \
        "${bam_down}" \
        "${bam_down%.bam}.fixmt.bam"

    if [[ -f "${bam_down%.bam}.fixmt.bam" ]]; then
        mv -f \
            "${bam_down%.bam}.fixmt.bam" \
            "${bam_down}"
    fi

    samtools sort \
        -@ "${threads}" \
        -o "${bam_down%.bam}.coord.bam" \
        "${bam_down}"

    if [[ -f "${bam_down%.bam}.coord.bam" ]]; then
        mv -f \
            "${bam_down%.bam}.coord.bam" \
            "${bam_down}"
    fi
    
    samtools index \
        -@ "${threads}" \
        "${bam_down}"
fi
```
</details>
<br />

<a id="list-and-tally-flags-in-bam-files-of-interest"></a>
### List and tally FLAGs in BAM files of interest
<details>
<summary><i>Bash code: List and tally FLAGs in BAM files of interest</i></summary>

```bash
#!/bin/bash

#  Define functions
calculate_run_time() {
    what="""
    calculate_run_time()
    --------------------
    Calculate run time for chunk of code
    
    :param 1: start time in \$(date +%s) format
    :param 2: end time in \$(date +%s) format
    :param 3: message to be displayed when printing the run time <chr>
    """
    run_time="$(echo "${2}" - "${1}" | bc -l)"
    
    echo ""
    echo "${3}"
    printf 'Run time: %dh:%dm:%ds\n' \
        $(( run_time/3600 )) \
        $(( run_time%3600/60 )) \
        $(( run_time%60 ))
    echo ""
}


display_spinning_icon() {
    what="""
    display_spinning_icon()
    -----------------------
    Display \"spinning icon\" while a background process runs
    
    :param 1: PID of the last program the shell ran in the background (int)
    :param 2: message to be displayed next to the spinning icon (chr)
    """
    spin="/|\\–"
    i=0
    while kill -0 "${1}" 2> /dev/null; do
        i=$(( (i + 1) % 4 ))
        printf "\r${spin:$i:1} %s" "${2}"
        sleep .15
    done
}


list_tally_flags() {
    what="""
    list_tally_flags()
    ------------------
    List and tally flags in a bam infile; function acts on a bam infile to
    perform piped commands (samtools view, cut, sort, uniq -c, sort -nr) that
    list and tally flags; function writes the results to a txt outfile, the
    name of which is derived from the txt infile
    
    :param 1: name of bam infile, including path (chr)
    """
    start="$(date +%s)"
    
    samtools view "${1}" \
        | cut -d$'\t' -f 2 \
        | sort \
        | uniq -c \
        | sort -nr \
            > "${1/.bam/.flags.txt}" &
    display_spinning_icon $! \
        "Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on $(basename "${1}")... "
        
    end="$(date +%s)"
    echo ""
    calculate_run_time "${start}" "${end}"  \
    "List and tally flags in $(basename "${1}")."
}


#  Define variables
bam_full="IP_G1_Hho1_6336.sort-coord.bam"
bam_down="IP_G1_Hho1_6336.sort-coord.ds.bam"


#  Do main work
list_tally_flags "${bam_full}"
list_tally_flags "${bam_down}"

cat "${bam_full%.bam}.flags.txt"
cat "${bam_down%.bam}.flags.txt"
```

</details>
<br />

<details>
<summary><i>Printed: List and tally FLAGs in BAM files of interest</i></summary>

```txt
┌─[kalavattam][Kriss-MacBook-Pro][±][main S:5 U:5 ?:29 ✗][~/repos/2023_tutorial_ChIP-seq]
└─▪  calculate_run_time() {
└─▪     what="""
└─▪     calculate_run_time()
└─▪     --------------------
└─▪     Calculate run time for chunk of code
└─▪
└─▪     :param 1: start time in \$(date +%s) format
└─▪     :param 2: end time in \$(date +%s) format
└─▪     :param 3: message to be displayed when printing the run time <chr>
└─▪     """
└─▪     run_time="$(echo "${2}" - "${1}" | bc -l)"
└─▪
└─▪     echo ""
└─▪     echo "${3}"
└─▪     printf 'Run time: %dh:%dm:%ds\n' \
└─▪         $(( run_time/3600 )) \
└─▪         $(( run_time%3600/60 )) \
└─▪         $(( run_time%60 ))
└─▪     echo ""
└─▪ }


┌─[kalavattam][Kriss-MacBook-Pro][±][main S:5 U:5 ?:29 ✗][~/repos/2023_tutorial_ChIP-seq]
└─▪  display_spinning_icon() {
└─▪     what="""
└─▪     display_spinning_icon()
└─▪     -----------------------
└─▪     Display \"spinning icon\" while a background process runs
└─▪
└─▪     :param 1: PID of the last program the shell ran in the background (int)
└─▪     :param 2: message to be displayed next to the spinning icon (chr)
└─▪     """
└─▪     spin="/|\\–"
└─▪     i=0
└─▪     while kill -0 "${1}" 2> /dev/null; do
└─▪         i=$(( (i + 1) % 4 ))
└─▪         printf "\r${spin:$i:1} %s" "${2}"
└─▪         sleep .15
└─▪     done
└─▪ }


┌─[kalavattam][Kriss-MacBook-Pro][±][main S:5 U:5 ?:29 ✗][~/repos/2023_tutorial_ChIP-seq]
└─▪  list_tally_flags() {
└─▪     what="""
└─▪     list_tally_flags()
└─▪     ------------------
└─▪     List and tally flags in a bam infile; function acts on a bam infile to
└─▪     perform piped commands (samtools view, cut, sort, uniq -c, sort -nr) that
└─▪     list and tally flags; function writes the results to a txt outfile, the
└─▪     name of which is derived from the txt infile
└─▪
└─▪     :param 1: name of bam infile, including path (chr)
└─▪     """
└─▪     start="$(date +%s)"
└─▪
└─▪     samtools view "${1}" \
└─▪         | cut -d$'\t' -f 2 \
└─▪         | sort \
└─▪         | uniq -c \
└─▪         | sort -nr \
└─▪             > "${1/.bam/.flags.txt}" &
└─▪     display_spinning_icon $! \
└─▪         "Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on $(basename "${1}")... "
└─▪
└─▪     end="$(date +%s)"
└─▪     echo ""
└─▪     calculate_run_time "${start}" "${end}"  \
└─▪     "List and tally flags in $(basename "${1}")."
└─▪ }


┌─[kalavattam][Kriss-MacBook-Pro][±][main S:5 U:5 ?:29 ✗][~/repos/2023_tutorial_ChIP-seq]
└─▪  bam_full="IP_G1_Hho1_6336.sort-coord.bam"


┌─[kalavattam][Kriss-MacBook-Pro][±][main S:5 U:5 ?:29 ✗][~/repos/2023_tutorial_ChIP-seq]
└─▪  bam_down="IP_G1_Hho1_6336.sort-coord.ds.bam"


┌─[kalavattam][Kriss-MacBook-Pro][±][main S:5 U:5 ?:29 ✗][~/repos/2023_tutorial_ChIP-seq]
└─▪  list_tally_flags "${bam_full}"
[1] 31348
– Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on IP_G1_Hho1_6336.sort-coord.bam... [1]+  Done
samtools view "${1}" | cut -d'    ' -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in IP_G1_Hho1_6336.sort-coord.bam.
Run time: 0h:2m:57s


┌─[kalavattam][Kriss-MacBook-Pro][±][main S:5 U:5 ?:31 ✗][~/repos/2023_tutorial_ChIP-seq]
└─▪  cat "${bam_full%.bam}.flags.txt"
7259170 99
7259170 147
7252868 83
7252868 163


┌─[kalavattam][Kriss-MacBook-Pro][±][main S:5 U:5 ?:31 ✗][~/repos/2023_tutorial_ChIP-seq]
└─▪  cat "${bam_down%.bam}.flags.txt"
725512 99
725512 147
723811 83
723811 163
```
</details>
<br />

<a id="coverage-tests-2024-0503-0506"></a>
### Coverage tests, 2024-0503-0506
<details>
<summary><i>Bash code: Coverage tests, 2024-0503-0506</i></summary>

```bash
#!/bin/bash

#  Define variables shared by each test call to `bamCoverage`
bam="IP_G1_Hho1_6336.sort-coord.bam"
outdir="coverage-tests_2024-0503"
chr="I"
format="bedgraph"
genome_size=12157105
bin_size=30
min_MAPQ=1
region="I"
threads=4
normalization="None"
flag=2


#  #1
#+ `bamCoverage` default behavior with specification of...
#+     - binSize: 30
#+     - region: I
#+     - normalizeUsing: None
#+     - effectiveGenomeSize: 12157105
#+     - miMappingQuality: 1
bamCoverage \
    --bam "${bam}" \
    --outFileName "${outdir}/${bam%.bam}.${region}.bamCoverage.${flag}.${format}" \
    --outFileFormat "${format}" \
    --binSize "${bin_size}" \
    --region "${region}" \
    --numberOfProcessors "${threads}" \
    --normalizeUsing "${normalization}" \
    --effectiveGenomeSize "${genome_size}" \
    --minMappingQuality "${min_MAPQ}"


#  #2
#+ `bamCoverage` default behavior with specification of FLAG 64: Use of only
#+ the first mate in proper pairs
flag=64

bamCoverage \
    --bam "${bam}" \
    --outFileName "${outdir}/${bam%.bam}.${region}.bamCoverage.${flag}.${format}" \
    --outFileFormat "${format}" \
    --binSize "${bin_size}" \
    --region "${region}" \
    --numberOfProcessors "${threads}" \
    --normalizeUsing "${normalization}" \
    --effectiveGenomeSize "${genome_size}" \
    --minMappingQuality "${min_MAPQ}" \
    --samFlagInclude "${flag}"


#  #3
#+ `bamCoverage` default behavior with specification of FLAG 128: Use of only
#+ the second mate in proper pairs
flag=128

bamCoverage \
    --bam "${bam}" \
    --outFileName "${outdir}/${bam%.bam}.${region}.bamCoverage.${flag}.${format}" \
    --outFileFormat "${format}" \
    --binSize "${bin_size}" \
    --region "${region}" \
    --numberOfProcessors "${threads}" \
    --normalizeUsing "${normalization}" \
    --effectiveGenomeSize "${genome_size}" \
    --minMappingQuality "${min_MAPQ}" \
    --samFlagInclude "${flag}"


#  #4
#+ `bamCoverage` with invocation of argument `--extendReads` and specification
#+ of...
#+     - binSize: 30
#+     - region: I
#+     - normalizeUsing: None
#+     - effectiveGenomeSize: 12157105
#+     - miMappingQuality: 1
flag=2

bamCoverage \
    --bam "${bam}" \
    --outFileName "${outdir}/${bam%.bam}.${region}.bamCoverage.${flag}.extendReads.${format}" \
    --outFileFormat "${format}" \
    --binSize "${bin_size}" \
    --region "${region}" \
    --numberOfProcessors "${threads}" \
    --normalizeUsing "${normalization}" \
    --effectiveGenomeSize "${genome_size}" \
    --minMappingQuality "${min_MAPQ}" \
    --extendReads


#  #5
#+ `bamCoverage` with invocation of argument `--extendReads` and specification
#+ of FLAG 64: Use of only the first mate in proper pairs
flag=64

bamCoverage \
    --bam "${bam}" \
    --outFileName "${outdir}/${bam%.bam}.${region}.bamCoverage.${flag}.extendReads.${format}" \
    --outFileFormat "${format}" \
    --binSize "${bin_size}" \
    --region "${region}" \
    --numberOfProcessors "${threads}" \
    --normalizeUsing "${normalization}" \
    --effectiveGenomeSize "${genome_size}" \
    --minMappingQuality "${min_MAPQ}" \
    --samFlagInclude "${flag}" \
    --extendReads


#  #6
#+ `bamCoverage` with invocation of argument `--extendReads` and specification
#+ of FLAG 128: Use of only the second mate in proper pairs
flag=128

bamCoverage \
    --bam "${bam}" \
    --outFileName "${outdir}/${bam%.bam}.${region}.bamCoverage.${flag}.extendReads.${format}" \
    --outFileFormat "${format}" \
    --binSize "${bin_size}" \
    --region "${region}" \
    --numberOfProcessors "${threads}" \
    --normalizeUsing "${normalization}" \
    --effectiveGenomeSize "${genome_size}" \
    --minMappingQuality "${min_MAPQ}" \
    --samFlagInclude "${flag}" \
    --extendReads
```
</details>
<br />

<details>
<summary><i>Printed: Coverage tests, 2024-0503</i></summary>

```txt
❯ bam="IP_G1_Hho1_6336.sort-coord.bam"
❯ outdir="coverage-tests_2024-0503"
❯ chr="I"
❯ format="bedgraph"
❯ genome_size=12157105
❯ bin_size=30
❯ min_MAPQ=1
❯ region="I"
❯ threads=4
❯ normalization="None"
❯ flag=2


❯ bamCoverage \
>     --bam "${bam}" \
>     --outFileName "${outdir}/${bam%.bam}.${region}.bamCoverage.${flag}.${format}" \
>     --outFileFormat "${format}" \
>     --binSize "${bin_size}" \
>     --region "${region}" \
>     --numberOfProcessors "${threads}" \
>     --normalizeUsing "${normalization}" \
>     --effectiveGenomeSize "${genome_size}" \
>     --minMappingQuality "${min_MAPQ}"
bamFilesList: ['IP_G1_Hho1_6336.sort-coord.bam']
binLength: 30
numberOfSamples: None
blackListFileName: None
skipZeroOverZero: False
bed_and_bin: False
genomeChunkSize: None
defaultFragmentLength: read length
numberOfProcessors: 4
verbose: False
region: I:30
bedFile: None
minMappingQuality: 1
ignoreDuplicates: False
chrsToSkip: []
stepSize: 30
center_read: False
samFlag_include: None
samFlag_exclude: None
minFragmentLength: 0
maxFragmentLength: 0
zerosToNans: False
smoothLength: None
save_data: False
out_file_for_raw_data: None
maxPairedFragmentLength: 1000


❯ flag=64


❯ bamCoverage \
>     --bam "${bam}" \
>     --outFileName "${outdir}/${bam%.bam}.${region}.bamCoverage.${flag}.${format}" \
>     --outFileFormat "${format}" \
>     --binSize "${bin_size}" \
>     --region "${region}" \
>     --numberOfProcessors "${threads}" \
>     --normalizeUsing "${normalization}" \
>     --effectiveGenomeSize "${genome_size}" \
>     --minMappingQuality "${min_MAPQ}" \
>     --samFlagInclude "${flag}"
bamFilesList: ['IP_G1_Hho1_6336.sort-coord.bam']
binLength: 30
numberOfSamples: None
blackListFileName: None
skipZeroOverZero: False
bed_and_bin: False
genomeChunkSize: None
defaultFragmentLength: read length
numberOfProcessors: 4
verbose: False
region: I:30
bedFile: None
minMappingQuality: 1
ignoreDuplicates: False
chrsToSkip: []
stepSize: 30
center_read: False
samFlag_include: 64
samFlag_exclude: None
minFragmentLength: 0
maxFragmentLength: 0
zerosToNans: False
smoothLength: None
save_data: False
out_file_for_raw_data: None
maxPairedFragmentLength: 1000


❯ flag=128


❯ bamCoverage \
>     --bam "${bam}" \
>     --outFileName "${outdir}/${bam%.bam}.${region}.bamCoverage.${flag}.${format}" \
>     --outFileFormat "${format}" \
>     --binSize "${bin_size}" \
>     --region "${region}" \
>     --numberOfProcessors "${threads}" \
>     --normalizeUsing "${normalization}" \
>     --effectiveGenomeSize "${genome_size}" \
>     --minMappingQuality "${min_MAPQ}" \
>     --samFlagInclude "${flag}"
bamFilesList: ['IP_G1_Hho1_6336.sort-coord.bam']
binLength: 30
numberOfSamples: None
blackListFileName: None
skipZeroOverZero: False
bed_and_bin: False
genomeChunkSize: None
defaultFragmentLength: read length
numberOfProcessors: 4
verbose: False
region: I:30
bedFile: None
minMappingQuality: 1
ignoreDuplicates: False
chrsToSkip: []
stepSize: 30
center_read: False
samFlag_include: 128
samFlag_exclude: None
minFragmentLength: 0
maxFragmentLength: 0
zerosToNans: False
smoothLength: None
save_data: False
out_file_for_raw_data: None
maxPairedFragmentLength: 1000
```
</details>
<br />

<details>
<summary><i>Printed: Coverage tests continued (`--extendReads`), 2024-0506</i></summary>

```txt
❯ bam="IP_G1_Hho1_6336.sort-coord.bam"
❯ outdir="coverage-tests_2024-0503"
❯ chr="I"
❯ format="bedgraph"
❯ genome_size=12157105
❯ bin_size=30
❯ min_MAPQ=1
❯ region="I"
❯ threads=4
❯ normalization="None"
❯ flag=2


❯ bamCoverage \
>     --bam "${bam}" \
>     --outFileName "${outdir}/${bam%.bam}.${region}.bamCoverage.${flag}.extendReads.${format}" \
>     --outFileFormat "${format}" \
>     --binSize "${bin_size}" \
>     --region "${region}" \
>     --numberOfProcessors "${threads}" \
>     --normalizeUsing "${normalization}" \
>     --effectiveGenomeSize "${genome_size}" \
>     --minMappingQuality "${min_MAPQ}" \
>     --extendReads
bamFilesList: ['IP_G1_Hho1_6336.sort-coord.bam']
binLength: 30
numberOfSamples: None
blackListFileName: None
skipZeroOverZero: False
bed_and_bin: False
genomeChunkSize: None
defaultFragmentLength: 210
numberOfProcessors: 4
verbose: False
region: I:30
bedFile: None
minMappingQuality: 1
ignoreDuplicates: False
chrsToSkip: []
stepSize: 30
center_read: False
samFlag_include: None
samFlag_exclude: None
minFragmentLength: 0
maxFragmentLength: 0
zerosToNans: False
smoothLength: None
save_data: False
out_file_for_raw_data: None
maxPairedFragmentLength: 840


❯ bamCoverage \
>     --bam "${bam}" \
>     --outFileName "${outdir}/${bam%.bam}.${region}.bamCoverage.${flag}.extendReads.${format}" \
>     --outFileFormat "${format}" \
>     --binSize "${bin_size}" \
>     --region "${region}" \
>     --numberOfProcessors "${threads}" \
>     --normalizeUsing "${normalization}" \
>     --effectiveGenomeSize "${genome_size}" \
>     --minMappingQuality "${min_MAPQ}" \
>     --samFlagInclude "${flag}" \
>     --extendReads
bamFilesList: ['IP_G1_Hho1_6336.sort-coord.bam']
binLength: 30
numberOfSamples: None
blackListFileName: None
skipZeroOverZero: False
bed_and_bin: False
genomeChunkSize: None
defaultFragmentLength: 210
numberOfProcessors: 4
verbose: False
region: I:30
bedFile: None
minMappingQuality: 1
ignoreDuplicates: False
chrsToSkip: []
stepSize: 30
center_read: False
samFlag_include: 64
samFlag_exclude: None
minFragmentLength: 0
maxFragmentLength: 0
zerosToNans: False
smoothLength: None
save_data: False
out_file_for_raw_data: None
maxPairedFragmentLength: 840


❯ bamCoverage \
>     --bam "${bam}" \
>     --outFileName "${outdir}/${bam%.bam}.${region}.bamCoverage.${flag}.extendReads.${format}" \
>     --outFileFormat "${format}" \
>     --binSize "${bin_size}" \
>     --region "${region}" \
>     --numberOfProcessors "${threads}" \
>     --normalizeUsing "${normalization}" \
>     --effectiveGenomeSize "${genome_size}" \
>     --minMappingQuality "${min_MAPQ}" \
>     --samFlagInclude "${flag}" \
>     --extendReads
bamFilesList: ['IP_G1_Hho1_6336.sort-coord.bam']
binLength: 30
numberOfSamples: None
blackListFileName: None
skipZeroOverZero: False
bed_and_bin: False
genomeChunkSize: None
defaultFragmentLength: 210
numberOfProcessors: 4
verbose: False
region: I:30
bedFile: None
minMappingQuality: 1
ignoreDuplicates: False
chrsToSkip: []
stepSize: 30
center_read: False
samFlag_include: 128
samFlag_exclude: None
minFragmentLength: 0
maxFragmentLength: 0
zerosToNans: False
smoothLength: None
save_data: False
out_file_for_raw_data: None
maxPairedFragmentLength: 840
```
</details>
<br />

<a id="assessmentscomparisons-of-bedgraph-files-generated-from-custom-code-and-bamcoverage-coverage-tests_2024-0503-0506"></a>
## Assessments/comparisons of BEDGRAPH files generated from custom code and `bamCoverage` (`coverage-tests_2024-0503-0506`)
<a id="images"></a>
### Images
<details>
<summary><i>IGV snapshot of BEDGRAPH files, I:1-12,000</i></summary>

<img src="./coverage-tests_2024-0503-0506/igv_snapshot_1.png" width="2000"/>
</details>
<br />

<details>
<summary><i>IGV snapshot of BEDGRAPH files, I:64,000-165,870</i></summary>

<img src="./coverage-tests_2024-0503-0506/igv_snapshot_2.png" width="2000"/>
</details>
<br />

<a id="breakdown-of-rows-in-images"></a>
### Breakdown of rows in images
<details>
<summary><i>Text: Breakdown of rows in images</i></summary>
<br />

- Rows 1&ndash;6: BEDGRAPH files from calls to deepTools `bamCoverage` without specification of normalization.
    + Rows 1&ndash;3: *Sans* read extension.
        - Row 1: Both FLAGs `64` and `128` (first and second in pair, respectively)
        - Row 2: Only FLAG `64` (first in pair)
        - Row 3: Only FLAG `128` (second in pair)
    + Rows 4&ndash;6: *With* read extension (`--extendReads`).
        - Row 4: Both FLAGs `64` and `128` (first and second in pair, respectively)
        - Row 5: Only FLAG `64` (first in pair)
        - Row 6: Only FLAG `128` (second in pair)
- Rows 7 and 8: BEDGRAPH files from calls to custom script `calculate_coverage.py`, which calculates "traditional coverage" and "normalized coverage" per Dickson et al., *Sci Rep* 2023.
    + Row 7: Traditional coverage.
    + Row 8: Normalized coverage.
</details>
<br />

<a id="observations"></a>
### Observations
<details>
<summary><i>Text: Observations</i></summary>
<br />

- The BEDGRAPH track of custom-calculated (`calculate_coverage.py`) tradtional coverage (row 7) is quite similar to the calls to deepTools `bamCoverage --extendReads` that make use of reads with only FLAG `64` (first in pair) or `128` (second in pair) (rows 5 and 6). However, why are the y-axis **~70 units lower** for the custom-calculated traditional coverage (row 7) versus the coverage from calls deepTools `bamCoverage --extendReads` (rows 5 and 6)?
- A call to `bamCoverage` without specifying a normalization (`--normalizeUsing None`) *does not* make use of effective genome size (`--effectiveGenomeSize`), despite what the documentation says:
> --normalizeUsing
> Possible choices: RPKM, CPM, BPM, RPGC, None
> 
> Use one of the entered methods to normalize the number of reads per bin. By default, no normalization is performed. RPKM = Reads Per Kilobase per Million mapped reads; CPM = Counts Per Million mapped reads, same as CPM in RNA-seq; BPM = Bins Per Million mapped reads, same as TPM in RNA-seq; RPGC = reads per genomic content (1x normalization); Mapped reads are considered after blacklist filtering (if applied). RPKM (per bin) = number of reads per bin / (number of mapped reads (in millions) * bin length (kb)). CPM (per bin) = number of reads per bin / number of mapped reads (in millions). BPM (per bin) = number of reads per bin / sum of all reads per bin (in millions). RPGC (per bin) = number of reads per bin / scaling factor for 1x average coverage. None = the default and equivalent to not setting this option at all. This scaling factor, in turn, is determined from the sequencing depth: (total number of mapped reads * fragment length) / effective genome size. The scaling factor used is the inverse of the sequencing depth computed for the sample to match the 1x coverage. This option requires --effectiveGenomeSize. Each read is considered independently, if you want to only count one mate from a pair in paired-end data, then use the --samFlagInclude/--samFlagExclude options. (Default: None)
</details>
<br />

<a id="on-the---normalizeusing-documentation-issue-ie-regarding-none-and-rpgc"></a>
#### On the `--normalizeUsing` documentation issue (i.e., regarding `None` and `RPGC`)
<details>
<summary><i>Me</i></summary>
<br />

Please consider the following portion of documentation for deepTools `bamCoverage` and let me know what you think deepTools `bamCoverage` is doing when `--normalizeUsing None` is invoked by a user:
> --normalizeUsing
> Possible choices: RPKM, CPM, BPM, RPGC, None
> 
> Use one of the entered methods to normalize the number of reads per bin. By default, no normalization is performed. RPKM = Reads Per Kilobase per Million mapped reads; CPM = Counts Per Million mapped reads, same as CPM in RNA-seq; BPM = Bins Per Million mapped reads, same as TPM in RNA-seq; RPGC = reads per genomic content (1x normalization); Mapped reads are considered after blacklist filtering (if applied). RPKM (per bin) = number of reads per bin / (number of mapped reads (in millions) * bin length (kb)). CPM (per bin) = number of reads per bin / number of mapped reads (in millions). BPM (per bin) = number of reads per bin / sum of all reads per bin (in millions). RPGC (per bin) = number of reads per bin / scaling factor for 1x average coverage. None = the default and equivalent to not setting this option at all. This scaling factor, in turn, is determined from the sequencing depth: (total number of mapped reads * fragment length) / effective genome size. The scaling factor used is the inverse of the sequencing depth computed for the sample to match the 1x coverage. This option requires --effectiveGenomeSize. Each read is considered independently, if you want to only count one mate from a pair in paired-end data, then use the --samFlagInclude/--samFlagExclude options. (Default: None)
</details>
<br />

<details>
<summary><i>GPT4</i></summary>
<br />

The `bamCoverage` documentation suggests that when the `--normalizeUsing None` option is invoked, no normalization is performed on the read counts per bin by default. This means that the output values you see are simply the raw counts of reads overlapping each genomic bin specified by the user.

However, the documentation seems to have a confusing statement regarding how "None" is described:

> "None = the default and equivalent to not setting this option at all. This scaling factor, in turn, is determined from the sequencing depth: (total number of mapped reads * fragment length) / effective genome size. The scaling factor used is the inverse of the sequencing depth computed for the sample to match the 1x coverage. This option requires --effectiveGenomeSize."

It appears there is a contradiction or an error in the documentation here. Under "None," no normalization or scaling should occur, meaning the raw read counts are reported without any adjustments for genome size, read length, or total read number. The discussion of a scaling factor, which seems to relate to normalization approaches like RPGC where coverage is adjusted to reflect a per-genome coverage normalization, doesn't align with the "None" setting.

From a practical perspective, using `--normalizeUsing None` in `bamCoverage` should result in output where each bin's value is simply the count of reads overlapping that bin, without any further transformation or scaling. This is typically used when you want to assess the raw coverage data, perhaps for visualization purposes where absolute counts are more informative than normalized data, or when downstream analyses require raw counts.

If you're seeing or expecting any other behavior (like unintended normalization), it might be due to a misunderstanding from the documentation or possibly a software bug. It might be helpful to directly test the output of `bamCoverage` with and without the `--normalizeUsing None` option on a small sample data set to see if there are any differences in the output. If there is unexpected behavior, it might warrant further investigation or a query to the deepTools development team.
</details>
<br />

<details>
<summary><i>Me</i></summary>
<br />

It seems to me that, in the documentation for deepTools `bamCoverage` argument `--normalizeUsing`, the following three sentences are out of place:
> This scaling factor, in turn, is determined from the sequencing depth: (total number of mapped reads * fragment length) / effective genome size. The scaling factor used is the inverse of the sequencing depth computed for the sample to match the 1x coverage. This option requires --effectiveGenomeSize.

Because these sentences follow this sentence on option "None",
> None = the default and equivalent to not setting this option at all.

it seems like the three sentences are describing what happens when the user specifies `--normalizeUsing "None"`. However, it is clear in [the source code](https://github.com/deeptools/deepTools/blob/master/deeptools/getScaleFactor.py#L251-L254) that these three sentences are describing what happens when the user specifies `--normalizeUsing "RPGC"`:
```python
current_coverage = \
            float(bam_mapped * fragment_length) / args.effectiveGenomeSize
# the scaling sets the coverage to match 1x
scale_factor *= 1.0 / current_coverage
```

I think it would be clearer if the documentation for deepTools `bamCoverage` argument `--normalizeUsing` read like this (note that the bolded text is to emphasize the new location of the sentence describing what `--normalizeUsing None` does; apart from that, nothing else in the documentation has been changed):
> RPGC (per bin) = number of reads per bin / scaling factor for 1x average coverage. None = the default and equivalent to not setting this option at all. This scaling factor, in turn, is determined from the sequencing depth: (total number of mapped reads * fragment length) / effective genome size. The scaling factor used is the inverse of the sequencing depth computed for the sample to match the 1x coverage. This option requires --effectiveGenomeSize. **None = the default and equivalent to not setting this option at all.** Each read is considered independently, if you want to only count one mate from a pair in paired-end data, then use the --samFlagInclude/--samFlagExclude options. (Default: None)

I'd like to write all this up for a GitHub Issue to be submitted to the maintainers of deepTools. If you agree with what I have written, then please help me rephrase all of this clearly and concisely for the GitHub Issue.
</details>
<br />

<details>
<summary><i>GPT4</i></summary>
<br />

You've made a very clear observation regarding the potential misplacement of the documentation concerning the --normalizeUsing None and RPGC options in the deepTools bamCoverage tool. Rearranging the documentation as you've suggested would indeed clarify the intended functionality and prevent confusion about the behavior of the None normalization setting.

Here's a concise rephrasing for a GitHub issue submission that captures your findings and suggestion:

**Title: Documentation Clarification Needed for --normalizeUsing Option in bamCoverage**

***Body:***

Hello deepTools team,

I've noticed a potential inconsistency in the `bamCoverage` documentation regarding the `--normalizeUsing` option. Specifically, the descriptions for the `None` and `RPGC` settings appear to be intertwined, which could lead to confusion about the behavior of these options.

***Current Documentation:***

> --normalizeUsing
> Possible choices: RPKM, CPM, BPM, RPGC, None
> 
> Use one of the entered methods to normalize the number of reads per bin. By default, no normalization is performed. RPKM = Reads Per Kilobase per Million mapped reads; CPM = Counts Per Million mapped reads, same as CPM in RNA-seq; BPM = Bins Per Million mapped reads, same as TPM in RNA-seq; RPGC = reads per genomic content (1x normalization); Mapped reads are considered after blacklist filtering (if applied). RPKM (per bin) = number of reads per bin / (number of mapped reads (in millions) * bin length (kb)). CPM (per bin) = number of reads per bin / number of mapped reads (in millions). BPM (per bin) = number of reads per bin / sum of all reads per bin (in millions). RPGC (per bin) = number of reads per bin / scaling factor for 1x average coverage. None = the default and equivalent to not setting this option at all. This scaling factor, in turn, is determined from the sequencing depth: (total number of mapped reads * fragment length) / effective genome size. The scaling factor used is the inverse of the sequencing depth computed for the sample to match the 1x coverage. This option requires --effectiveGenomeSize. Each read is considered independently, if you want to only count one mate from a pair in paired-end data, then use the --samFlagInclude/--samFlagExclude options. (Default: None)

***Issue:***
The description suggests that "None" is equivalent to not setting the normalization option, yet it immediately follows with a detailed explanation of calculating a scaling factor, which seems more appropriate under the `RPGC` option. This might incorrectly imply that `None` involves some form of normalization.

***Suggested Change:***
It would be clearer to move the detailed explanation about the scaling factor under the RPGC option directly, and reaffirm that `None` truly means no normalization. Perhaps something like this:

> RPGC (per bin) = number of reads per bin / scaling factor for 1x average coverage. This scaling factor is determined from the sequencing depth: (total number of mapped reads * fragment length) / effective genome size. The scaling factor used is the inverse of the sequencing depth computed for the sample to match the 1x coverage. This option requires --effectiveGenomeSize.
None = the default and equivalent to not setting this option at all.
Each read is considered independently, if you want to only count one mate from a pair in paired-end data, then use the --samFlagInclude/--samFlagExclude options. (Default: None)

This rearrangement would help clarify the documentation and ensure users correctly understand the behavior of each normalization option.

Thank you for considering this clarification to enhance the utility and user-friendliness of deepTools.

This draft clearly states the issue and provides a constructive suggestion for improving the documentation. You can adjust the technical details to match your exact understanding or any additional insights you wish to include.
</details>
<br />

<details>
<summary><i>Me</i></summary>
<br />

Submission of edited version of the above text to the deepTools team via [a GitHub Issue](https://github.com/deeptools/deepTools/issues/1311).
</details>
</details>
<br />
<br />

<a id="snippets-how-different-normalizations-are-calculated-in-the-deeptools-codebase"></a>
## Snippets: How different normalizations are calculated in the deepTools codebase
<details>
<summary><i>Text and code</i></summary>
<br />

- Calculation of coverage if `--normalizeUsing "RPGC"`:
    + [`getScaleFactor.py#L216-L254`](https://github.com/deeptools/deepTools/blob/master/deeptools/getScaleFactor.py#L216-L251)
    + [`getScaleFactor.py#L251-L254`](https://github.com/deeptools/deepTools/blob/master/deeptools/getScaleFactor.py#L251-L254)
    ```python
current_coverage = \
        float(bam_mapped * fragment_length) / args.effectiveGenomeSize
# the scaling sets the coverage to match 1x
scale_factor *= 1.0 / current_coverage 
    ```
- Calculation of coverage if `--normalizeUsing "RPKM"`
    + [`getScaleFactor.py#L263-L268`](https://github.com/deeptools/deepTools/blob/master/deeptools/getScaleFactor.py#L263-L268)
    + [`getScaleFactor.py#L259-L268`](https://github.com/deeptools/deepTools/blob/master/deeptools/getScaleFactor.py#L259-L268)
    ```python
# the RPKM is the # reads per tile / \
#    ( total reads (in millions) * tile length in Kb)
million_reads_mapped = float(bam_mapped) / 1e6
tile_len_in_kb = float(args.binSize) / 1000

scale_factor *= 1.0 / (million_reads_mapped * tile_len_in_kb)
    ```
- Calculation of coverage if `--normalizeUsing "CPM"`
    + [`getScaleFactor.py#L277-L279`](https://github.com/deeptools/deepTools/blob/master/deeptools/getScaleFactor.py#L277-L279)
    + [`getScaleFactor.py#L273-L282`](https://github.com/deeptools/deepTools/blob/master/deeptools/getScaleFactor.py#L273-L282)
    ```python
# the CPM (norm is based on post-filtering total counts of reads in BAM "bam_mapped")
million_reads_mapped = float(bam_mapped) / 1e6
scale_factor *= 1.0 / (million_reads_mapped)
    ```
- Calculation of coverage if `--normalizeUsing "BPM"`
    + [`getScaleFactor.py#L284-L294`](https://github.com/deeptools/deepTools/blob/master/deeptools/getScaleFactor.py#L284-L294)
    + [`getScaleFactor.py#L287-L292`](https://github.com/deeptools/deepTools/blob/master/deeptools/getScaleFactor.py#L287-L292)
    ```python
# the BPM (norm is based on post-filtering total counts of reads in BAM "bam_mapped")
# sampled_bins_sum = getSampledSum(args.bam)
tile_len_in_kb = float(args.binSize) / 1000
tpm_scaleFactor = (bam_mapped / tile_len_in_kb) / 1e6

scale_factor *= 1 / (tpm_scaleFactor * tile_len_in_kb)
    ```
</details>
<br />
<br />

<a id="reproducing-what-deeptools-bamcoverage-does-via-simplified-versions-of-its-methods"></a>
## Reproducing what deepTools `bamCoverage` does via simplified versions of its methods
<a id="drafting-a-simplified-version-of-count_reads_in_region"></a>
### Drafting a simplified version of [`count_reads_in_region`](https://github.com/deeptools/deepTools/blob/master/deeptools/countReadsPerBin.py#L400)
<details>
<summary><i>Me</i></summary>
<br />

I'd like to write my own version of the following function that makes use of `pysam`:

<details>
<summary><i>Python code</i></summary>

```python
def count_reads_in_region(self, chrom, start, end, bed_regions_list=None):
    """Counts the reads in each bam file at each 'stepSize' position
    within the interval (start, end) for a window or bin of size binLength.

    The stepSize controls the distance between bins. For example,
    a step size of 20 and a bin size of 20 will create bins next to
    each other. If the step size is smaller than the bin size the
    bins will overlap.

    If a list of bedRegions is given, then the number of reads
    that overlaps with each region is counted.

    Parameters
    ----------
    chrom : str
        Chrom name
    start : int
        start coordinate
    end : int
        end coordinate
    bed_regions_list: list
        List of list of tuples of the form (start, end)
        corresponding to bed regions to be processed.
        If not bed file was passed to the object constructor
        then this list is empty.

    Returns
    -------
    numpy array
        The result is a numpy array that as rows each bin
        and as columns each bam file.


    Examples
    --------
    Initialize some useful values

    >>> test = Tester()
    >>> c = CountReadsPerBin([test.bamFile1, test.bamFile2], 25, 0, stepSize=50)

    The transpose is used to get better looking numbers. The first line
    corresponds to the number of reads per bin in the first bamfile.

    >>> _array, __ = c.count_reads_in_region(test.chrom, 0, 200)
    >>> _array
    array([[0., 0.],
           [0., 1.],
           [1., 1.],
           [1., 2.]])

    """

    if start > end:
        raise NameError("start %d bigger that end %d" % (start, end))

    if self.stepSize is None and bed_regions_list is None:
        raise ValueError("stepSize is not set!")
    # array to keep the read counts for the regions
    subnum_reads_per_bin = []

    start_time = time.time()

    bam_handles = []
    for fname in self.bamFilesList:
        try:
            bam_handles.append(bamHandler.openBam(fname))
        except SystemExit:
            sys.exit(sys.exc_info()[1])
        except:
            bam_handles.append(pyBigWig.open(fname))

    blackList = None
    if self.blackListFileName is not None:
        blackList = GTF(self.blackListFileName)

    # A list of lists of tuples
    transcriptsToConsider = []
    if bed_regions_list is not None:
        if self.bed_and_bin:
            transcriptsToConsider.append([(x[1][0][0], x[1][0][1], self.binLength) for x in bed_regions_list])
        else:
            transcriptsToConsider = [x[1] for x in bed_regions_list]
    else:
        if self.stepSize == self.binLength:
            transcriptsToConsider.append([(start, end, self.binLength)])
        else:
            for i in range(start, end, self.stepSize):
                if i + self.binLength > end:
                    break
                if blackList is not None and blackList.findOverlaps(chrom, i, i + self.binLength):
                    continue
                transcriptsToConsider.append([(i, i + self.binLength)])

    if self.save_data:
        _file = open(deeptools.utilities.getTempFileName(suffix='.bed'), 'w+t')
        _file_name = _file.name
    else:
        _file_name = ''

    for bam in bam_handles:
        for trans in transcriptsToConsider:
            tcov = self.get_coverage_of_region(bam, chrom, trans)
            if bed_regions_list is not None and not self.bed_and_bin:
                subnum_reads_per_bin.append(np.sum(tcov))
            else:
                subnum_reads_per_bin.extend(tcov)

    subnum_reads_per_bin = np.concatenate([subnum_reads_per_bin]).reshape(-1, len(self.bamFilesList), order='F')

    if self.save_data:
        idx = 0
        for i, trans in enumerate(transcriptsToConsider):
            if len(trans[0]) != 3:
                starts = ",".join([str(x[0]) for x in trans])
                ends = ",".join([str(x[1]) for x in trans])
                _file.write("\t".join([chrom, starts, ends]) + "\t")
                _file.write("\t".join(["{}".format(x) for x in subnum_reads_per_bin[i, :]]) + "\n")
            else:
                for exon in trans:
                    for startPos in range(exon[0], exon[1], exon[2]):
                        if idx >= subnum_reads_per_bin.shape[0]:
                            # At the end of chromosomes (or due to blacklisted regions), there are bins smaller than the bin size
                            # Counts there are added to the bin before them, but range() will still try to include them.
                            break
                        _file.write("{0}\t{1}\t{2}\t".format(chrom, startPos, min(startPos + exon[2], exon[1])))
                        _file.write("\t".join(["{}".format(x) for x in subnum_reads_per_bin[idx, :]]) + "\n")
                        idx += 1
        _file.close()

    if self.verbose:
        endTime = time.time()
        rows = subnum_reads_per_bin.shape[0]
        print("%s countReadsInRegions_worker: processing %d "
              "(%.1f per sec) @ %s:%s-%s" %
              (multiprocessing.current_process().name,
               rows, rows / (endTime - start_time), chrom, start, end))

    return subnum_reads_per_bin, _file_name
```
</details>
<br />

We don't need to implement the BED file thing, nor do we need to implement the handling of a blacklist. What do you think? Is this doable?
</details>
<br />

<details>
<summary><i>GPT4</i></summary>
<br />

Certainly, you can simplify the `count_reads_in_region` function by stripping out the functionality related to BED files and blacklist handling, focusing purely on counting reads within specified genomic intervals using `pysam`. Below is a streamlined version of your function using `pysam` to count reads directly:

<details>
<summary><i>Python code</i></summary>

```python
import pysam
import numpy as np

def count_reads_in_region(bam_file_path, chrom, start, end, step_size, bin_length):
    """
    Counts reads in each specified region of a BAM file.

    Parameters:
    - bam_file_path (str): Path to the BAM file.
    - chrom (str): Chromosome name.
    - start (int): Start coordinate of the region.
    - end (int): End coordinate of the region.
    - step_size (int): Distance between consecutive bins.
    - bin_length (int): Length of each bin.

    Returns:
    - numpy array: Array where each element is the count of reads in each bin.
    """

    # Open the BAM file using pysam.
    bam = pysam.AlignmentFile(bam_file_path, "rb")

    # Calculate the number of bins based on step size and region length.
    num_bins = (end - start) // step_size
    read_counts = np.zeros(num_bins, dtype=int)

    # Iterate through each bin and count reads.
    for i in range(num_bins):
        bin_start = start + i * step_size
        bin_end = bin_start + bin_length
        if bin_end > end:
            bin_end = end

        # Fetch reads within the current bin.
        for read in bam.fetch(chrom, bin_start, bin_end):
            if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:
                read_counts[i] += 1

    return read_counts

# Example usage:
bam_path = 'example.bam'
chromosome = 'chr1'
region_start = 100000
region_end = 101000
bin_size = 50
step = 50

# Get read counts
counts = count_reads_in_region(bam_path, chromosome, region_start, region_end, step, bin_size)
print("Size of array:", counts.size)
print("Read counts per bin:", counts)
```
</details>
<br />

**Explanation:**
1. *Function Definition:* The function is defined with parameters to specify the BAM file, chromosome, start and end of the region, step size between bins, and bin length.
2. *BAM File Handling:* Opens the BAM file using `pysam`.
3. *Bins Calculation:* Calculates the number of bins based on the `step_size` and the length of the region. Initializes an array `read_counts` to store read counts for each bin.
4. *Read Counting:* Iterates over each bin, using `pysam.fetch()` to retrieve reads within the bin coordinates. It increments the count in the respective bin if the read is properly mapped (ignoring unmapped, secondary, or supplementary alignments).
5. *Return:* Returns the numpy array containing the read counts for each bin.

This function efficiently calculates the number of reads per bin over a specified genomic region, using step size and bin length parameters to determine how the region is divided. Adjust the `step` and `bin_length` as needed based on the desired resolution for counting reads.
</details>
<br />

<a id="the-simplified-version-of-count_reads_in_region"></a>
#### The simplified version of `count_reads_in_region`
<details>
<summary><i>Python code: The simplified version of `count_reads_in_region`</i></summary>

```python
import pysam
import numpy as np

def count_reads_in_region(bam_file_path, chrom, start, end, step_size, bin_length):
    """
    Counts reads in each specified region of a BAM file.

    Parameters:
    - bam_file_path (str): Path to the BAM file.
    - chrom (str): Chromosome name.
    - start (int): Start coordinate of the region.
    - end (int): End coordinate of the region.
    - step_size (int): Distance between consecutive bins.
    - bin_length (int): Length of each bin.

    Returns:
    - numpy array: Array where each element is the count of reads in each bin.
    """

    # Open the BAM file using pysam.
    bam = pysam.AlignmentFile(bam_file_path, "rb")

    # Calculate the number of bins based on step size and region length.
    num_bins = (end - start) // step_size
    read_counts = np.zeros(num_bins, dtype=int)

    # Iterate through each bin and count reads.
    for i in range(num_bins):
        bin_start = start + i * step_size
        bin_end = bin_start + bin_length
        if bin_end > end:
            bin_end = end

        # Fetch reads within the current bin.
        for read in bam.fetch(chrom, bin_start, bin_end):
            if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:
                read_counts[i] += 1

    return read_counts
```
</details>
<br />

<a id="breaking-down-simplified-version-of-count_reads_in_regions-parameter-step_size"></a>
#### Breaking down simplified version of `count_reads_in_region`'s parameter `step_size`
<details>
<summary><i>Me</i></summary>
<br />

Let's take a step back and consider the code in this function and how the function operates:

<details>
<summary><i>Python code</i></summary>

```python
import pysam
import numpy as np

def count_reads_in_region(bam_file_path, chrom, start, end, step_size, bin_length):
    """
    Counts reads in each specified region of a BAM file.

    Parameters:
    - bam_file_path (str): Path to the BAM file.
    - chrom (str): Chromosome name.
    - start (int): Start coordinate of the region.
    - end (int): End coordinate of the region.
    - step_size (int): Distance between consecutive bins.
    - bin_length (int): Length of each bin.

    Returns:
    - numpy array: Array where each element is the count of reads in each bin.
    """

    # Open the BAM file using pysam.
    bam = pysam.AlignmentFile(bam_file_path, "rb")

    # Calculate the number of bins based on step size and region length.
    num_bins = (end - start) // step_size
    read_counts = np.zeros(num_bins, dtype=int)

    # Iterate through each bin and count reads.
    for i in range(num_bins):
        bin_start = start + i * step_size
        bin_end = bin_start + bin_length
        if bin_end > end:
            bin_end = end

        # Fetch reads within the current bin.
        for read in bam.fetch(chrom, bin_start, bin_end):
            if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:
                read_counts[i] += 1

    return read_counts
```
</details>
<br />

What would happen if `step_size` was a smaller number than `bin_length`? And what if `step_size` was a larger number than `bin_length`? And what if both `step_size` and `bin_length` had the same size?
</details>
<br />

<details>
<summary><i>GPT4</i></summary>
<br />

The `count_reads_in_region` function you've written utilizes `step_size` and `bin_length` to control how reads are counted across a specified region in a BAM file. These parameters influence the placement and overlap of bins where reads are counted. Let's break down what happens under different scenarios regarding the values of step_size and `bin_length`:

**1. Step Size Smaller Than Bin Length**  
When `step_size` is smaller than `bin_length`, the bins in which reads are counted will overlap. This scenario allows for a more detailed and continuous coverage profile across the region, as each bin's start is less than the length of the previous bin's length. This overlap can be useful for detecting patterns and ensuring that features aren't missed due to the placement of bin boundaries.

*Example:*
- `bin_length` = 100
- `step_size` = 50
- This setup results in 50 base pairs of overlap between consecutive bins.

**2. Step Size Larger Than Bin Length**  
If `step_size` is larger than `bin_length`, there will be gaps between each bin where no reads are counted. This might be useful in scenarios where you only need periodic snapshots of the coverage but could miss important features in the regions between bins.

*Example:*
- `bin_length` = 50
- `step_size` = 100
- This setup results in 50 base pairs between bins where reads are not counted.

**3. Step Size Equal to Bin Length**  
When `step_size` equals `bin_length`, each bin is adjacent to the next without overlap or gaps. This is a common setup for generating histograms or coverage maps where each segment of the genome is counted once and only once.

*Example:*
- `bin_length` = 50
- `step_size` = 50
- This setup results in bins that cover the region continuously and completely without any overlap or gaps.

**Implementation Notes**  
- Efficiency Considerations: The function iterates through the BAM file multiple times (once for each bin). This could be inefficient for large regions or high coverage data. It might be beneficial to optimize this by counting reads in a single pass and allocating counts to bins as reads are encountered.
- Edge Cases: The function adjusts the `bin_end` to not exceed the `end` of the region. This adjustment is crucial for the last bin, especially when the region's length isn't an exact multiple of `step_size`.
</details>
<br />

<a id="drafting-simplified-versions-of-is_proper_pair-and-get_fragment_from_read"></a>
### Drafting simplified versions of `is_proper_pair` and `get_fragment_from_read`
<details>
<summary><i>Me</i></summary>
<br />

With your help, I'd like to write simplified versions of the following two functions. Care to assist me with this? Keep in mind that the only major packages we should be drawing from are `pysam` and `numpy`.

<details>
<summary><i>Python code</i></summary>

```python
def is_proper_pair(read, maxPairedFragmentLength):
    """
    Checks if a read is proper pair meaning that both mates are facing each other and are in
    the same chromosome and are not to far away. The sam flag for proper pair can not
    always be trusted. Note that if the fragment size is > maxPairedFragmentLength (~2kb
    usually) that False will be returned.
    :return: bool

    >>> import pysam
    >>> import os
    >>> from deeptools.countReadsPerBin import CountReadsPerBin as cr
    >>> root = os.path.dirname(os.path.abspath(__file__)) + "/test/test_data/"
    >>> bam = pysam.AlignmentFile("{}/test_proper_pair_filtering.bam".format(root))
    >>> iter = bam.fetch()
    >>> read = next(iter)
    >>> cr.is_proper_pair(read, 1000) # "keep" read
    True
    >>> cr.is_proper_pair(read, 200) # "keep" read, but maxPairedFragmentLength is too short
    False
    >>> read = next(iter)
    >>> cr.is_proper_pair(read, 1000) # "improper pair"
    False
    >>> read = next(iter)
    >>> cr.is_proper_pair(read, 1000) # "mismatch chr"
    False
    >>> read = next(iter)
    >>> cr.is_proper_pair(read, 1000) # "same orientation1"
    False
    >>> read = next(iter)
    >>> cr.is_proper_pair(read, 1000) # "same orientation2"
    False
    >>> read = next(iter)
    >>> cr.is_proper_pair(read, 1000) # "rev first"
    False
    >>> read = next(iter)
    >>> cr.is_proper_pair(read, 1000) # "rev first OK"
    True
    >>> read = next(iter)
    >>> cr.is_proper_pair(read, 1000) # "for first"
    False
    >>> read = next(iter)
    >>> cr.is_proper_pair(read, 1000) # "for first"
    True
    """
    if not read.is_proper_pair:
        return False
    if read.reference_id != read.next_reference_id:
        return False
    if abs(read.template_length) > maxPairedFragmentLength:
        return False
    # check that the mates face each other (inward)
    if read.is_reverse is read.mate_is_reverse:
        return False
    if read.is_reverse:
        if read.reference_start >= read.next_reference_start:
            return True
    else:
        if read.reference_start <= read.next_reference_start:
            return True
    return False


def get_fragment_from_read(self, read):
    """Get read start and end position of a read.
    If given, the reads are extended as follows:
    If reads are paired end, each read mate is extended to match
    the fragment length, otherwise, a default fragment length
    is used. If reads are split (give by the CIGAR string) then
    the multiple positions of the read are returned.
    When reads are extended the cigar information is
    skipped.

    Parameters
    ----------
    read: pysam object.

    The following values are defined (for forward reads)::


             |--          -- read.tlen --              --|
             |-- read.alen --|
        -----|===============>------------<==============|----
             |               |            |
        read.reference_start
                    read.reference_end  read.pnext

          and for reverse reads


             |--             -- read.tlen --           --|
                                         |-- read.alen --|
        -----|===============>-----------<===============|----
             |                           |               |
          read.pnext           read.reference_start  read.reference_end

    this is a sketch of a pair-end reads

    The function returns the fragment start and end, either
    using the paired end information (if available) or
    extending the read in the appropriate direction if this
    is single-end.

    Parameters
    ----------
    read : pysam read object


    Returns
    -------
    list of tuples
        [(fragment start, fragment end)]


    >>> test = Tester()
    >>> c = CountReadsPerBin([], 1, 1, 200, extendReads=True)
    >>> c.defaultFragmentLength=100
    >>> c.get_fragment_from_read(test.getRead("paired-forward"))
    [(5000000, 5000100)]
    >>> c.get_fragment_from_read(test.getRead("paired-reverse"))
    [(5000000, 5000100)]
    >>> c.defaultFragmentLength = 200
    >>> c.get_fragment_from_read(test.getRead("single-forward"))
    [(5001491, 5001691)]
    >>> c.get_fragment_from_read(test.getRead("single-reverse"))
    [(5001536, 5001736)]
    >>> c.defaultFragmentLength = 'read length'
    >>> c.get_fragment_from_read(test.getRead("single-forward"))
    [(5001491, 5001527)]
    >>> c.defaultFragmentLength = 'read length'
    >>> c.extendReads = False
    >>> c.get_fragment_from_read(test.getRead("paired-forward"))
    [(5000000, 5000036)]

    Tests for read centering.

    >>> c = CountReadsPerBin([], 1, 1, 200, extendReads=True, center_read=True)
    >>> c.defaultFragmentLength = 100
    >>> assert c.get_fragment_from_read(test.getRead("paired-forward")) == [(5000032, 5000068)]
    >>> c.defaultFragmentLength = 200
    >>> assert c.get_fragment_from_read(test.getRead("single-reverse")) == [(5001618, 5001654)]
    """
    # if no extension is needed, use pysam get_blocks
    # to identify start and end reference positions.
    # get_blocks return a list of start and end positions
    # based on the CIGAR if skipped regions are found.
    # E.g for a cigar of 40M260N22M
    # get blocks return two elements for the first 40 matches
    # and the for the last 22 matches.
    if self.defaultFragmentLength == 'read length':
        return read.get_blocks()

    else:
        if self.is_proper_pair(read, self.maxPairedFragmentLength):
            if read.is_reverse:
                fragmentStart = read.next_reference_start
                fragmentEnd = read.reference_end
            else:
                fragmentStart = read.reference_start
                # the end of the fragment is defined as
                # the start of the forward read plus the insert length
                fragmentEnd = read.reference_start + abs(read.template_length)

        # Extend using the default fragment length
        else:
            if read.is_reverse:
                fragmentStart = read.reference_end - self.defaultFragmentLength
                fragmentEnd = read.reference_end
            else:
                fragmentStart = read.reference_start
                fragmentEnd = read.reference_start + self.defaultFragmentLength

    if self.center_read:
        fragmentCenter = fragmentEnd - (fragmentEnd - fragmentStart) / 2
        fragmentStart = int(fragmentCenter - read.infer_query_length(always=False) / 2)
        fragmentEnd = fragmentStart + read.infer_query_length(always=False)

    assert fragmentStart < fragmentEnd, "fragment start greater than fragment" \
                                        "end for read {}".format(read.query_name)
    return [(fragmentStart, fragmentEnd)]
```
</details>
</details>
<br />

<details>
<summary><i>GPT4</i></summary>
<br />

Let's work together to create simplified versions of the two functions you've listed using pysam and numpy. We'll keep the functions as intuitive and minimalistic as possible while ensuring they perform their intended tasks efficiently.

**Function 1: Check if a Read is a Proper Pair**
The first function, `is_proper_pair`, checks if a read from a BAM file is properly paired according to several criteria, including mate orientation, chromosome alignment, and fragment length. Here's a simplified version:

<details>
<summary><i>Python code</i></summary>

```python
import pysam

def is_proper_pair(read, max_fragment_length):
    """
    Checks if the BAM read is a proper pair.

    Args:
    read (pysam.AlignedSegment): A single read from a BAM file.
    max_fragment_length (int): Maximum allowed fragment length.

    Returns:
    bool: True if the read is a proper pair, False otherwise.
    """
    # Check basic flags and chromosome alignment
    if not read.is_proper_pair or read.reference_id != read.next_reference_id:
        return False
    
    # Check fragment size
    if abs(read.template_length) > max_fragment_length:
        return False

    # Check orientation for paired-end reads
    if read.is_reverse == read.mate_is_reverse:
        return False

    # Additional check for correct mate orientation
    if read.is_reverse:
        return read.reference_start >= read.next_reference_start
    else:
        return read.reference_start <= read.next_reference_start
```
</details>
<br />

**Function 2: Get Fragment from Read**
The second function, `get_fragment_from_read`, determines the start and end positions of a read fragment, extending it if necessary based on certain conditions. Here's a streamlined version:

<details>
<summary><i>Python code</i></summary>

```python
def get_fragment_from_read(read, default_fragment_length, extend_reads=True, center_read=False):
    """
    Retrieves or calculates the start and end positions of the read or its fragment.

    Args:
    read (pysam.AlignedSegment): A single read from a BAM file.
    default_fragment_length (int): Default length to which single-end reads or unmatched pairs should be extended.
    extend_reads (bool): Whether to extend reads to the fragment length.
    center_read (bool): Whether to center the read around the midpoint of its fragment.

    Returns:
    tuple: The start and end positions of the read or its fragment.
    """
    if not extend_reads:
        return read.get_blocks()

    # If read is properly paired and aligned
    if read.is_proper_pair and read.reference_id == read.next_reference_id:
        if read.is_reverse:
            start = read.next_reference_start
            end = read.reference_end
        else:
            start = read.reference_start
            end = read.reference_start + abs(read.template_length)
    else:
        # Handle single reads or improperly paired reads
        if read.is_reverse:
            start = max(0, read.reference_end - default_fragment_length)
            end = read.reference_end
        else:
            start = read.reference_start
            end = read.reference_start + default_fragment_length

    # Adjust for centering the read, if requested
    if center_read:
        center = start + (end - start) // 2
        half_length = default_fragment_length // 2
        start, end = center - half_length, center + half_length

    return start, end
```
</details>
<br />

**Integration and Testing**
To integrate and test these functions:
1. Integration: You can incorporate these functions into a script that processes a BAM file. For example, you might iterate over reads in a BAM file and apply these functions to each read to analyze read pairing and calculate fragment sizes.
2. Testing: For testing, you could use a small subset of your BAM data or synthetic BAM files created to represent specific scenarios (proper pairs, improper pairs, single ends, etc.). This would allow you to verify that the functions behave as expected under controlled conditions.

Let me know if you need further assistance with implementing these functions or setting up tests!
</details>
<br />

<a id="the-simplified-version-of-is_proper_pair"></a>
#### The simplified version of `is_proper_pair`
<details>
<summary><i>Python code: The simplified version of `is_proper_pair`</i></summary>

```python
import pysam

def is_proper_pair(read, max_fragment_length):
    """
    Checks if the BAM read is a proper pair.

    Args:
    read (pysam.AlignedSegment): A single read from a BAM file.
    max_fragment_length (int): Maximum allowed fragment length.

    Returns:
    bool: True if the read is a proper pair, False otherwise.
    """
    # Check basic flags and chromosome alignment
    if not read.is_proper_pair or read.reference_id != read.next_reference_id:
        return False
    
    # Check fragment size
    if abs(read.template_length) > max_fragment_length:
        return False

    # Check orientation for paired-end reads
    if read.is_reverse == read.mate_is_reverse:
        return False

    # Additional check for correct mate orientation
    if read.is_reverse:
        return read.reference_start >= read.next_reference_start
    else:
        return read.reference_start <= read.next_reference_start
```
</details>
<br />

<a id="the-simplified-version-of-get_fragment_from_read"></a>
#### The simplified version of `get_fragment_from_read`
<details>
<summary><i>Python code: The simplified version of `get_fragment_from_read`</i></summary>

```python
def get_fragment_from_read(read, default_fragment_length, extend_reads=True, center_read=False):
    """
    Retrieves or calculates the start and end positions of the read or its fragment.

    Args:
    read (pysam.AlignedSegment): A single read from a BAM file.
    default_fragment_length (int): Default length to which single-end reads or unmatched pairs should be extended.
    extend_reads (bool): Whether to extend reads to the fragment length.
    center_read (bool): Whether to center the read around the midpoint of its fragment.

    Returns:
    tuple: The start and end positions of the read or its fragment.
    """
    if not extend_reads:
        return read.get_blocks()

    # If read is properly paired and aligned
    if read.is_proper_pair and read.reference_id == read.next_reference_id:
        if read.is_reverse:
            start = read.next_reference_start
            end = read.reference_end
        else:
            start = read.reference_start
            end = read.reference_start + abs(read.template_length)
    else:
        # Handle single reads or improperly paired reads
        if read.is_reverse:
            start = max(0, read.reference_end - default_fragment_length)
            end = read.reference_end
        else:
            start = read.reference_start
            end = read.reference_start + default_fragment_length

    # Adjust for centering the read, if requested
    if center_read:
        center = start + (end - start) // 2
        half_length = default_fragment_length // 2
        start, end = center - half_length, center + half_length

    return start, end
```
</details>
<br />

<a id="setting-up-tests-for-is_proper_pair-and-get_fragment_from_read"></a>
#### Setting up tests for `is_proper_pair` and `get_fragment_from_read`
<details>
<summary><i>Me</i></summary>
<br />

Yes, please assist with setting up tests. For example, here's how I set up a very simple test for `count_reads_in_region`:

<details>
<summary><i>Python code</i></summary>

```python
#  Example usage
bam_path = 'IP_G1_Hho1_6336.sort-coord.ds.bam'
chromosome = 'I'
region_start = 1
region_end = 230218
bin_size = 20
step = 20

#  Get read counts
counts = count_reads_in_region(
    bam_path, chromosome, region_start, region_end, step, bin_size
)
print("Number of elements in array:", counts.size)
print("Read counts per bin:", counts)
```
</details>
<br />

Something as simple as this for the two functions would be fine to start with (as far as I am concerned).
</details>
<br />

<details>
<summary><i>GPT4</i></summary>
<br />

To set up simple tests for the functions `is_proper_pair` and `get_fragment_from_read`, you can follow a similar approach to your example for `count_reads_in_region`. Let's create a test setup that uses a small BAM file and checks a few reads to validate the functionality of each function.

**Setup for Testing**
First, ensure that you have a small BAM file that includes a variety of reads (properly paired, improperly paired, single-end, etc.) to thoroughly test different scenarios. For testing, you can extract a small portion of your BAM file with samtools or create one manually.

**Sample Test Code**
Here’s how you might set up the tests for the two functions:

<details>
<summary><i>Python code</i></summary>

```python
import pysam

# Load the functions you've defined or imported.
from your_module import is_proper_pair, get_fragment_from_read

# Define the path to your BAM file
bam_path = 'small_sample.bam'  # Ensure this BAM file exists and is indexed

# Open the BAM file
bam_file = pysam.AlignmentFile(bam_path, "rb")

# Fetch a few reads to test
reads = list(bam_file.fetch('I', 100, 200))  # Adjust coordinates based on your BAM file content

# Test is_proper_pair
print("Testing is_proper_pair function:")
for read in reads:
    print(f"Read {read.query_name} proper pair status: {is_proper_pair(read, 2000)}")

# Test get_fragment_from_read
print("\nTesting get_fragment_from_read function:")
for read in reads:
    fragment = get_fragment_from_read(read, 500, extend_reads=True)
    print(f"Read {read.query_name} fragment: Start {fragment[0]}, End {fragment[1]}")

# Close the BAM file
bam_file.close()
```
</details>
<br />

**What This Test Does:**
1. Load Functions: Import the functions you want to test from wherever they are defined.
2. Open BAM File: Open a BAM file that includes relevant test data.
3. Fetch Reads: Fetch a set of reads from a specific region. Adjust the `'I'`, `100`, `200` parameters to match a region in your BAM file that will provide a good variety of reads.
4. Test `is_proper_pair`: Iterate through the fetched reads and print out whether each read is considered a proper pair according to the specified maximum fragment length.
5. Test `get_fragment_from_read`: Similarly, print out the calculated fragment start and end for each read, considering the default fragment length and whether to extend the reads.

**Adjustments:**
- Regions and BAM Path: Ensure that the chromosome names and coordinates match what is available in your BAM file.
- Function Parameters: Adjust the default fragment length and maximum fragment length according to what makes sense for your data.

By running this script, you'll be able to see immediate outputs that indicate whether each function is behaving as expected. This method is a good starting point for more complex automated testing or integration into a larger pipeline. If you need to set up more formal unit tests or require further customization, let me know!
</details>
<br />

<a id="additional-breakdown-of-get_fragment_from_read"></a>
#### Additional breakdown of `get_fragment_from_read`
<details>
<summary><i>Me</i></summary>
<br />

We wrote the following function:
<details>
<summary><i>Python code</i></summary>

```python
def get_fragment_from_read(read, default_fragment_length, extend_reads=True, center_read=False):
    """
    Retrieves or calculates the start and end positions of the read or its fragment.

    Args:
    read (pysam.AlignedSegment): A single read from a BAM file.
    default_fragment_length (int): Default length to which single-end reads or unmatched pairs should be extended.
    extend_reads (bool): Whether to extend reads to the fragment length.
    center_read (bool): Whether to center the read around the midpoint of its fragment.

    Returns:
    tuple: The start and end positions of the read or its fragment.
    """
    if not extend_reads:
        return read.get_blocks()

    # If read is properly paired and aligned
    if read.is_proper_pair and read.reference_id == read.next_reference_id:
        if read.is_reverse:
            start = read.next_reference_start
            end = read.reference_end
        else:
            start = read.reference_start
            end = read.reference_start + abs(read.template_length)
    else:
        # Handle single reads or improperly paired reads
        if read.is_reverse:
            start = max(0, read.reference_end - default_fragment_length)
            end = read.reference_end
        else:
            start = read.reference_start
            end = read.reference_start + default_fragment_length

    # Adjust for centering the read, if requested
    if center_read:
        center = start + (end - start) // 2
        half_length = default_fragment_length // 2
        start, end = center - half_length, center + half_length

    return start, end
```
</details>
<br />

Please explain to me what is happening in the function line by line.
</details>
<br />

<details>
<summary><i>GPT4</i></summary>
<br />

Here's a line-by-line explanation of the function `get_fragment_from_read`:
1. *Function Definition:* The function `get_fragment_from_read` is defined with four parameters: `read`, `default_fragment_length`, `extend_reads`, and `center_read`. It takes a read from a BAM file and potentially extends it based on the given parameters.
2. *Docstring:* The docstring describes the function, its parameters, and what it returns. It clearly explains the purpose and usage of each parameter.
3. *Check if Extending Reads:* The first conditional (`if not extend_reads:`) checks if the `extend_reads` parameter is `False`. If so, the function simply returns the blocks of the read (start and end positions of the aligned segments from the CIGAR string) without extending them.
4. *Properly Paired Read Handling:* The main if condition (`if read.is_proper_pair and read.reference_id == read.next_reference_id:`) checks if the read is part of a properly paired alignment and if both reads are mapped to the same chromosome.
    + *Reverse Strand Logic:* If the read is on the reverse strand (`if read.is_reverse:`), the function sets `start` to the start of the mate (next read) and `end` to the end of the current read.
    + *Forward Strand Logic:* Otherwise, for forward strand reads, it sets `start` to the start of the current read and `end` to the start plus the template length (distance between paired reads).
5. *Single or Improperly Paired Read Handling:* If the read is not properly paired or not paired at all:
    + *Reverse Strand Single Read:* For reads on the reverse strand, it calculates `start` by subtracting the `default_fragment_length` from the read's end, ensuring it doesn’t go below zero, and sets end to the read's end.
    + *Forward Strand Single Read:* For forward reads, it sets `start` to the read's start and `end` to the start plus the `default_fragment_length`.
6. *Center Read Adjustment:* If `center_read` is `True`, the function recalculates the `start` and `end` to center the read around the middle of the fragment defined by `default_fragment_length`. It calculates the center, then adjusts `start` and `end` to be equally distant from this center, making the fragment centered around the midpoint of the original or extended alignment.
7. *Return:* The function returns a tuple `(start, end)`, representing the start and end positions of the read or its fragment after any adjustments for extension and centering.

This function is flexible in handling both paired and single-end reads and allows for extensions and adjustments based on the specific needs of the analysis, such as extending reads to cover more genomic territory or centering reads for more uniform coverage analysis.
</details>
<br />

<details>
<summary><i>Me</i></summary>
<br />

What is the `//` operation in the following snippet from `get_fragment_from_read`?
<details>
<summary><i>Python code</i></summary>

```python
    # Adjust for centering the read, if requested
    if center_read:
        center = start + (end - start) // 2
        half_length = default_fragment_length // 2
        start, end = center - half_length, center + half_length
```
</details>
</details>
<br />

<details>
<summary><i>GPT4</i></summary>
<br />

The `//` operator in Python performs floor division. Unlike the standard division operator `/`, which returns a floating point result, floor division returns the largest possible integer that is less than or equal to the division result. This means it discards any fractional part and returns an integer result.

In the context of the code snippet you provided:
<details>
<summary><i>Python code</i></summary>

```python
    # Adjust for centering the read, if requested
    if center_read:
        center = start + (end - start) // 2
        half_length = default_fragment_length // 2
        start, end = center - half_length, center + half_length
```
</details>
<br />

- `center = start + (end - start) // 2`: This calculates the midpoint between start and end. `(end - start) // 2` computes half the distance between start and end but rounds down to the nearest integer, ensuring that the center point remains a whole number.
- `half_length = default_fragment_length // 2`: This calculates half of the `default_fragment_length`, again rounding down to ensure it remains an integer.

This use of floor division ensures that all calculations remain in integer space, which is important when dealing with genomic coordinates or array indices that must be whole numbers.
</details>
<br />
<br />

<a id="drafting-a-simplified-version-of-get_read_and_fragment_length-renamed-to-determine_fragment_and_read_lengths"></a>
### Drafting a simplified version of `get_read_and_fragment_length` (renamed to `determine_fragment_and_read_lengths`)
<a id="the-original-function"></a>
#### The original function
<details>
<summary><i>Me</i></summary>
<br />

Let's write a simplified version of the following function. The function is taken from [this script](https://github.com/deeptools/deepTools/blob/master/deeptools/getFragmentAndReadSize.py) in the deepTools source code.

<details>
<summary><i>Python code</i></summary>

```python
def get_read_and_fragment_length(bamFile, return_lengths=False, blackListFileName=None,
                                 binSize=50000, distanceBetweenBins=1000000,
                                 numberOfProcessors=None, verbose=False):
    """
    Estimates the fragment length and read length through sampling

    Parameters
    ----------
    bamFile : str
        BAM file name
    return_lengths : bool
    numberOfProcessors : int
    verbose : bool
    binSize : int
    distanceBetweenBins : int

    Returns
    -------
    d : dict
        tuple of two dictionaries, one for the fragment length and the other
for the read length. The dictionaries summarise the mean, median etc. values

    """

    bam_handle = bamHandler.openBam(bamFile)
    chrom_sizes = list(zip(bam_handle.references, bam_handle.lengths))

    distanceBetweenBins *= 2
    fl = []

    # Fix issue #522, allow distanceBetweenBins == 0
    if distanceBetweenBins == 0:
        imap_res = mapReduce.mapReduce((bam_handle.filename, distanceBetweenBins),
                                       getFragmentLength_wrapper,
                                       chrom_sizes,
                                       genomeChunkLength=binSize,
                                       blackListFileName=blackListFileName,
                                       numberOfProcessors=numberOfProcessors,
                                       verbose=verbose)
        fl = np.concatenate(imap_res)

    # Try to ensure we have at least 1000 regions from which to compute statistics, halving the intra-bin distance as needed
    while len(fl) < 1000 and distanceBetweenBins > 1:
        distanceBetweenBins /= 2
        stepsize = binSize + distanceBetweenBins
        imap_res = mapReduce.mapReduce((bam_handle.filename, distanceBetweenBins),
                                       getFragmentLength_wrapper,
                                       chrom_sizes,
                                       genomeChunkLength=stepsize,
                                       blackListFileName=blackListFileName,
                                       numberOfProcessors=numberOfProcessors,
                                       verbose=verbose)

        fl = np.concatenate(imap_res)

    if len(fl):
        fragment_length = fl[:, 0]
        read_length = fl[:, 1]
        if fragment_length.mean() > 0:
            fragment_len_dict = {'sample_size': len(fragment_length),
                                 'min': fragment_length.min(),
                                 'qtile25': np.percentile(fragment_length, 25),
                                 'mean': np.mean(fragment_length),
                                 'median': np.median(fragment_length),
                                 'qtile75': np.percentile(fragment_length, 75),
                                 'max': fragment_length.max(),
                                 'std': np.std(fragment_length),
                                 'mad': np.median(np.abs(fragment_length - np.median(fragment_length))),
                                 'qtile10': np.percentile(fragment_length, 10),
                                 'qtile20': np.percentile(fragment_length, 20),
                                 'qtile30': np.percentile(fragment_length, 30),
                                 'qtile40': np.percentile(fragment_length, 40),
                                 'qtile60': np.percentile(fragment_length, 60),
                                 'qtile70': np.percentile(fragment_length, 70),
                                 'qtile80': np.percentile(fragment_length, 80),
                                 'qtile90': np.percentile(fragment_length, 90),
                                 'qtile99': np.percentile(fragment_length, 99)}
        else:
            fragment_len_dict = None

        if return_lengths and fragment_len_dict is not None:
            fragment_len_dict['lengths'] = fragment_length

        read_len_dict = {'sample_size': len(read_length),
                         'min': read_length.min(),
                         'qtile25': np.percentile(read_length, 25),
                         'mean': np.mean(read_length),
                         'median': np.median(read_length),
                         'qtile75': np.percentile(read_length, 75),
                         'max': read_length.max(),
                         'std': np.std(read_length),
                         'mad': np.median(np.abs(read_length - np.median(read_length))),
                         'qtile10': np.percentile(read_length, 10),
                         'qtile20': np.percentile(read_length, 20),
                         'qtile30': np.percentile(read_length, 30),
                         'qtile40': np.percentile(read_length, 40),
                         'qtile60': np.percentile(read_length, 60),
                         'qtile70': np.percentile(read_length, 70),
                         'qtile80': np.percentile(read_length, 80),
                         'qtile90': np.percentile(read_length, 90),
                         'qtile99': np.percentile(read_length, 99)}
        if return_lengths:
            read_len_dict['lengths'] = read_length
    else:
        fragment_len_dict = None
        read_len_dict = None

    return fragment_len_dict, read_len_dict
```
</details>
</details>
<br />

<a id="the-simplified-version"></a>
#### The simplified version
<details>
<summary><i>GPT4</i></summary>
<br />

To simplify the `get_read_and_fragment_length` function from deepTools, we can focus on the core task of estimating fragment and read lengths using pysam, without the advanced sampling and multiprocessing used in the original function. Here's a streamlined version that uses direct calculation for a given number of reads from a BAM file:
<details>
<summary><i>Python code</i></summary>

```python
import pysam
import numpy as np

def determine_fragment_and_read_lengths(bam_file_path, sample_size=1000):
    """
    Estimate the mean fragment and read length by sampling reads from a BAM file.

    Args:
    bam_file_path (str): Path to the BAM file.
    sample_size (int): Number of reads to sample for estimating lengths.

    Returns:
    dict: Dictionary containing mean and median fragment and read lengths.
    """
    bam = pysam.AlignmentFile(bam_file_path, "rb")
    fragment_lengths = []
    read_lengths = []

    # Sample a fixed number of reads
    for read in bam.fetch():
        if read.is_proper_pair:
            # Calculate fragment length
            if read.is_reverse:
                fragment_length = read.reference_start - read.next_reference_start
            else:
                fragment_length = read.next_reference_start - read.reference_start
            
            fragment_lengths.append(abs(fragment_length))
            read_lengths.append(read.infer_query_length())

        if len(fragment_lengths) >= sample_size:
            break

    # Compute statistics
    fragment_data = {
        'mean_fragment_length': np.mean(fragment_lengths),
        'median_fragment_length': np.median(fragment_lengths),
        'mean_read_length': np.mean(read_lengths),
        'median_read_length': np.median(read_lengths)
    }

    return fragment_data

# Example usage:
bam_path = 'path_to_your.bam'
result = determine_fragment_and_read_lengths(bam_path)
print(result)
```
</details>
<br />

**Key Changes:**
1. Direct Fetch and Calculation: The function directly fetches reads from the BAM file and calculates the lengths without the need for multiprocessing or complex binning strategies.
2. Fragment Length Calculation: It computes the absolute fragment length based on paired-end information available directly from the reads.
3. Statistics Calculation: It calculates mean and median lengths for both fragments and reads based only on the sampled data.

This simplified approach ensures you can still gather useful metrics about your sequencing data without the complexity of the original function, making it more accessible for direct usage or integration into other scripts.
</details>
<br />

<a id="on-including-multiprocessing-in-the-simplified-version"></a>
#### On including multiprocessing in the simplified version
<details>
<summary><i>Me</i></summary>
<br />

If possible, I would like to include some form of multiprocessing.
</details>
<br />

<details>
<summary><i>GPT4</i></summary>
<br />

To include multiprocessing in the function `determine_fragment_and_read_lengths`, we can utilize Python's multiprocessing library to parallelize the reading and processing of the BAM file. This approach can significantly speed up the processing for large BAM files by utilizing multiple CPU cores to compute fragment and read lengths concurrently.

Here's an updated version of the function that incorporates multiprocessing:
<details>
<summary><i>Python code</i></summary>

```python
import pysam
import numpy as np
from multiprocessing import Pool, cpu_count

def process_reads(bam_file_path, regions):
    """
    Process a subset of regions in the BAM file to compute read and fragment lengths.

    Args:
    bam_file_path (str): Path to the BAM file.
    regions (list): List of regions (tuples of (contig, start, end)) to process.

    Returns:
    list: Fragment lengths and read lengths from the specified regions.
    """
    bam = pysam.AlignmentFile(bam_file_path, "rb")
    fragment_lengths = []
    read_lengths = []

    for contig, start, end in regions:
        for read in bam.fetch(contig, start, end):
            if read.is_proper_pair and not read.is_secondary and not read.is_supplementary:
                # Calculate fragment length
                fragment_length = abs(read.template_length)
                fragment_lengths.append(fragment_length)
                read_lengths.append(read.infer_query_length())

    return fragment_lengths, read_lengths

def determine_fragment_and_read_lengths(bam_file_path, num_processes=None):
    """
    Estimate the mean fragment and read length by sampling reads from a BAM file with multiprocessing.

    Args:
    bam_file_path (str): Path to the BAM file.
    num_processes (int): Number of processes to use.

    Returns:
    dict: Dictionary containing mean and median fragment and read lengths.
    """
    if num_processes is None:
        num_processes = cpu_count()

    # Open the BAM file to determine regions to process
    bam = pysam.AlignmentFile(bam_file_path, "rb")
    total_length = sum(bam.lengths)
    chunk_size = total_length // num_processes
    regions = []

    current_pos = 0
    for contig, length in zip(bam.references, bam.lengths):
        while current_pos < length - chunk_size:
            regions.append((contig, current_pos, current_pos + chunk_size))
            current_pos += chunk_size
        current_pos = 0

    # Process each chunk using multiple processes
    with Pool(num_processes) as pool:
        results = pool.starmap(process_reads, [(bam_file_path, region) for region in regions])

    # Combine results from all processes
    fragment_lengths = []
    read_lengths = []
    for fl, rl in results:
        fragment_lengths.extend(fl)
        read_lengths.extend(rl)

    # Compute statistics
    fragment_data = {
        'mean_fragment_length': np.mean(fragment_lengths),
        'median_fragment_length': np.median(fragment_lengths),
        'mean_read_length': np.mean(read_lengths),
        'median_read_length': np.median(read_lengths)
    }

    return fragment_data

# Example usage:
bam_path = 'path_to_your.bam'
result = determine_fragment_and_read_lengths(bam_path)
print(result)
```
</details>
<br />

**Modifications:**
1. *Chunking the BAM File:* The BAM file is divided into chunks based on the total length and the number of processes. Each chunk corresponds to a specific region of the genome.
2. *Multiprocessing:* The `multiprocessing.Pool` is used to process each chunk concurrently. The function `process_reads` is designed to be called by each process, handling a subset of regions and computing the read and fragment lengths.
3. *Combining Results:* After all processes complete, the results are aggregated to compute the overall statistics.

This approach allows efficient use of multiple CPU cores, thus reducing the computation time for large datasets.
</details>
<br />

<a id="on-running-integration-tests-for-process_reads-and-determine_fragment_and_read_lengths"></a>
#### On running integration tests for `process_reads` and `determine_fragment_and_read_lengths`
<details>
<summary><i>Me</i></summary>
<br />

Let's set up simple integration tests for `process_reads` and `determine_fragment_and_read_lengths`.
</details>
<br />

<details>
<summary><i>GPT4</i></summary>
<br />

To set up simple integration tests for `process_reads` and `determine_fragment_and_read_lengths`, you'll need a testing framework like `pytest` or `unittest` in Python. Here, I'll outline a basic approach using pytest for both functions. Ensure you have a test BAM file and corresponding indexed file (`.bam` and `.bai`) available for testing.

**1. Setting Up the Environment**
First, install pytest if it isn't already installed:
<details>
<summary><i>Bash code</i></summary>

```bash
pip install pytest
```
</details>
<br />

**2. Test File Setup**
Create a Python test file, for example `test_read_processing.py`.

**3. Writing Tests for `determine_fragment_and_read_lengths`**
Assuming `determine_fragment_and_read_lengths` function estimates the average read and fragment lengths from a BAM file, your test might look something like this:
<details>
<summary><i>Python code</i></summary>

```python
import pytest
from your_module import determine_fragment_and_read_lengths

def test_get_read_and_fragment_length():
    bam_file_path = "path/to/test.bam"  # Ensure this is a small BAM file for quick tests
    # Expected values should be known from manual inspection or previous reliable outputs
    expected_fragment_length = 150  # example value
    expected_read_length = 100  # example value

    fragment_length, read_length = determine_fragment_and_read_lengths(bam_file_path)
    
    assert fragment_length == expected_fragment_length, f"Expected fragment length {expected_fragment_length}, got {fragment_length}"
    assert read_length == expected_read_length, f"Expected read length {expected_read_length}, got {read_length}"
```
</details>
<br />

**4. Writing Tests for `process_reads`**
If `process_reads` is expected to process reads in some manner (e.g., filtering, mapping modifications), a test could look like this:
<details>
<summary><i>Python code</i></summary>

```python
from your_module import process_reads

def test_process_reads():
    bam_file_path = "path/to/test.bam"
    processed_data = process_reads(bam_file_path)
    
    # Verify some conditions that should be met by the processed data
    # This could be counts of reads, specific flags set/unset, or other data characteristics
    assert len(processed_data) > 0, "No data processed."
    # Add more specific assertions based on what process_reads is expected to do
```
</details>
<br />

**5. Running the Tests**
Execute the tests by running `pytest` in your terminal:
<details>
<summary><i>Bash code</i></summary>

```bash
pytest test_read_processing.py
```
</details>
<br />

**Considerations**
- *Test Data:* Make sure your test BAM file is representative of the data your functions are meant to handle but small enough to keep tests fast.
- *Expected Results:* Your tests should have known expected results. These could be hard-coded into your tests as shown, or you could use a fixture to load expected data from a file.
- *Environment:* Ensure that the environment where tests are run has all the necessary dependencies installed and configured (e.g., `pysam`, `numpy`).

This setup provides a straightforward way to integrate tests into your development process, improving code reliability and facilitating future modifications.
</details>
<br />

<a id="refactoring-get_fragment_from_read-such-that-each-mate-in-a-pair-of-reads-is-extended-to-match-the-fragment-length-etc"></a>
### Refactoring `get_fragment_from_read` such that each mate in a pair of reads is extended to match the fragment length, etc.
<details>
<summary><i>Me</i></summary>
<br />

We need to adapt and otherwise refactor `get_fragment_from_read` such that if reads are from paired-end sequencing, then each mate in a pair of reads is extended to match the fragment length for that pair of reads. Otherwise, a default fragment length should be used, and that value should have been determined from running `determine_fragment_and_read_lengths`—that is, the median fragment length returned by that function . If reads are split (as determined by interpreting the CIGAR string), then the multiple positions of the read should be returned. However, when reads are extended, then the assessment of CIGAR information should be skipped.

Current working version of `get_fragment_from_read`:
<details>
<summary><i>Python code</i></summary>

```python
def get_fragment_from_read(
    read, default_fragment_length, extend_reads=True, center_read=False
):
    """
    Retrieves or calculates the start and end positions of the read or its
    fragment.

    Args:
    read (pysam.AlignedSegment): A single read from a BAM file.
    default_fragment_length (int): Default length to which single-end reads or
        unmatched pairs should be extended.
    extend_reads (bool): Whether to extend reads to the fragment length.
    center_read (bool): Whether to center the read around the midpoint of its
        fragment.

    Returns:
    tuple: The start and end positions of the read or its fragment.
    """
    if not extend_reads:
        return read.get_blocks()

    # If read is properly paired and aligned
    if (
        read.is_proper_pair
        and read.reference_id == read.next_reference_id
    ):
        if read.is_reverse:
            start = read.next_reference_start
            end = read.reference_end
        else:
            start = read.reference_start
            end = read.reference_start + abs(read.template_length)
    else:
        # Handle single reads or improperly paired reads
        if read.is_reverse:
            start = max(0, read.reference_end - default_fragment_length)
            end = read.reference_end
        else:
            start = read.reference_start
            end = read.reference_start + default_fragment_length

    # Adjust for centering the read, if requested
    if center_read:
        center = start + (end - start) // 2
        half_length = default_fragment_length // 2
        start, end = center - half_length, center + half_length

    return start, end
```
</details>
</details>
<br />

<details>
<summary><i>GPT4</i></summary>
<br />

To adapt the `get_fragment_from_read` function based on your requirements, the following modifications and enhancements will ensure that it properly extends paired-end reads to match the fragment length, skips CIGAR information when extending, and utilizes the median fragment length from a previous analysis as the default length. The function will also handle split reads when not extending by returning all aligned segments based on the CIGAR string.

Here's how you can refactor `get_fragment_from_read`:

<details>
<summary><i>Python code</i></summary>

```python
def get_fragment_from_read(read, fragment_length_dict, extend_reads=True, center_read=False):
    """
    Retrieves or calculates the start and end positions of the read or its
    fragment based on whether the reads are extended to match the fragment length.

    Args:
    read (pysam.AlignedSegment): A single read from a BAM file.
    fragment_length_dict (dict): Dictionary containing 'median_fragment_length'
        which is used as the default length for extending single-end reads or unmatched pairs.
    extend_reads (bool): Whether to extend reads to the fragment length.
    center_read (bool): Whether to center the read around the midpoint of its fragment.

    Returns:
    list of tuples: Each tuple is a start and end position of a fragment.
    """
    if not extend_reads:
        # Return all aligned blocks based on CIGAR string
        return read.get_blocks()

    # Use median fragment length from previous analysis if no specific length provided
    default_fragment_length = fragment_length_dict.get('median_fragment_length', 200)  # Default if not provided

    fragments = []
    if read.is_proper_pair and read.reference_id == read.next_reference_id:
        # For paired-end reads where mates are properly paired, extend to match the fragment length
        if read.is_reverse:
            fragment_start = read.next_reference_start
            fragment_end = read.reference_end
        else:
            fragment_start = read.reference_start
            fragment_end = read.reference_start + abs(read.template_length)
        fragments.append((fragment_start, fragment_end))
    else:
        # Handle single or improperly paired reads
        if read.is_reverse:
            fragment_start = max(0, read.reference_end - default_fragment_length)
            fragment_end = read.reference_end
        else:
            fragment_start = read.reference_start
            fragment_end = read.reference_start + default_fragment_length
        fragments.append((fragment_start, fragment_end))

    if center_read:
        # If centering is requested, adjust fragments to be centered around the read's midpoint
        centered_fragments = []
        for start, end in fragments:
            center = start + (end - start) // 2
            half_length = default_fragment_length // 2
            centered_fragments.append((center - half_length, center + half_length))
        return centered_fragments

    return fragments

# Usage example
read_data = pysam.AlignedSegment()  # This would be a read object from pysam
fragment_length_info = {'median_fragment_length': 300}  # Example median fragment length obtained from analysis
fragments = get_fragment_from_read(read_data, fragment_length_info, extend_reads=True, center_read=False)
print(fragments)
```
</details>
<br />

**Explanation:**
1. *Handling CIGAR Information:* When `extend_reads` is `False`, the function returns all alignment blocks for the read as specified by the CIGAR string, allowing detailed analysis of split reads.
2. *Extension Logic:* For properly paired reads, each mate's alignment is extended to the full fragment length. For single-end or improperly paired reads, a default fragment length is used.
3. *Centering Adjustments:* If `center_read` is enabled, the function recalculates the start and end positions to center the fragment around the midpoint, using the provided or default fragment length.

This version of the function ensures flexibility in handling different types of read alignments while providing options for extending and centering reads based on the sequencing data's characteristics.
</details>
<br />

<a id="testing-the-simplified-functions"></a>
## Testing the simplified functions
<details>
<summary><i>Python code: Testing the simplified functions</i></summary>

```python
#  Import libraries, define functions =========================================
import pysam
import numpy as np
from multiprocessing import Pool, cpu_count


def count_reads_in_region(
    bam_file_path, chrom, start, end, step_size, bin_length
):
    """
    Counts reads in each specified region of a BAM file.

    Args:
    bam_file_path (str): Path to the BAM file.
    chrom (str): Chromosome name.
    start (int): Start coordinate of the region.
    end (int): End coordinate of the region.
    step_size (int): Distance between consecutive bins.
    bin_length (int): Length of each bin.

    Returns:
    read_counts (numpy array): Array where each element is the count of reads
        in each bin.
    """

    #  Open the BAM file using pysam.
    bam = pysam.AlignmentFile(bam_file_path, "rb")

    #  Calculate the number of bins based on step size and region length.
    num_bins = (end - start) // step_size
    read_counts = np.zeros(num_bins, dtype=int)

    #  Iterate through each bin and count reads.
    for i in range(num_bins):
        bin_start = start + i * step_size
        bin_end = bin_start + bin_length
        if bin_end > end:
            bin_end = end

        #  Fetch reads within the current bin.
        for read in bam.fetch(chrom, bin_start, bin_end):
            if (
                not read.is_unmapped
                and not read.is_secondary
                and not read.is_supplementary
            ):
                read_counts[i] += 1

    return read_counts


def is_proper_pair(read, max_fragment_length):
    """
    Checks if the BAM read is a proper pair.

    Args:
    read (pysam.AlignedSegment): A single read from a BAM file.
    max_fragment_length (int): Maximum allowed fragment length.

    Returns:
    bool: True if the read is a proper pair, False otherwise.
    """
    #  Check basic flags and chromosome alignment
    if (
        not read.is_proper_pair
        or read.reference_id != read.next_reference_id
    ):
        return False
    
    #  Check fragment size
    if abs(read.template_length) > max_fragment_length:
        return False

    #  Check orientation for paired-end reads
    if read.is_reverse == read.mate_is_reverse:
        return False

    #  Additional check for correct mate orientation
    if read.is_reverse:
        return read.reference_start >= read.next_reference_start
    else:
        return read.reference_start <= read.next_reference_start


def process_reads(bam_file_path, regions):
    """
    Process a subset of regions in the BAM file to compute read and fragment
    lengths.

    Args:
    bam_file_path (str): Path to the BAM file.
    regions (list): List of regions (tuples of (contig, start, end)) to
        process.

    Returns:
    fragment_lengths, read_lengths (lists): Fragment lengths and read lengths
        from the specified regions.
    """
    bam = pysam.AlignmentFile(bam_file_path, "rb")
    fragment_lengths = []
    read_lengths = []

    for contig, start, end in regions:
        for read in bam.fetch(contig, start, end):
            if (
                read.is_proper_pair
                and not read.is_secondary
                and not read.is_supplementary
            ):
                #  Calculate fragment length
                fragment_length = abs(read.template_length)
                fragment_lengths.append(fragment_length)
                read_lengths.append(read.infer_query_length())

    return fragment_lengths, read_lengths


def determine_fragment_and_read_lengths(bam_file_path, num_processes=None):
    """
    Estimate the mean fragment and read length by sampling reads from a BAM
    file with multiprocessing.

    Args:
    bam_file_path (str): Path to the BAM file.
    num_processes (int): Number of processes to use.

    Returns:
    fragment_read_data (dict): Dictionary containing mean and median fragment
        and read lengths.
    """
    if num_processes is None:
        num_processes = cpu_count()

    # Open the BAM file to determine regions to process
    bam = pysam.AlignmentFile(bam_file_path, "rb")
    total_length = sum(bam.lengths)
    chunk_size = total_length // num_processes
    regions = []

    # current_pos = 0
    # for contig, length in zip(bam.references, bam.lengths):
    #     while current_pos < length - chunk_size:
    #         regions.append((contig, current_pos, current_pos + chunk_size))
    #         current_pos += chunk_size
    #     current_pos = 0
    current_pos = 0
    for contig, length in zip(bam.references, bam.lengths):
        while current_pos + chunk_size < length:
            regions.append((contig, current_pos, current_pos + chunk_size))
            current_pos += chunk_size
        # Handle the last chunk which might be smaller than chunk_size
        if current_pos < length:
            regions.append((contig, current_pos, length))
        current_pos = 0  # reset for the next contig

    # Process each chunk using multiple processes
    with Pool(num_processes) as pool:
        results = pool.starmap(
            process_reads, [(bam_file_path, region) for region in regions]
        )

    # Combine results from all processes
    fragment_lengths = []
    read_lengths = []
    for fl, rl in results:
        fragment_lengths.extend(fl)
        read_lengths.extend(rl)

    # Compute statistics
    fragment_read_data = {
        'mean_fragment_length': np.mean(fragment_lengths),
        'median_fragment_length': np.median(fragment_lengths),
        'mean_read_length': np.mean(read_lengths),
        'median_read_length': np.median(read_lengths)
    }

    return fragment_read_data


def get_fragment_from_read(
    read, fragment_length_dict, extend_reads=True, center_read=False
):
    """
    Retrieves or calculates the start and end positions of the read or its
    fragment based on whether the reads are extended to match the fragment
    length.

    Args:
    read (pysam.AlignedSegment): A single read from a BAM file.
    fragment_length_dict (dict): Dictionary containing 'median_fragment_length'
        which is used as the default length for extending single-end reads or
        unmatched pairs.
    extend_reads (bool): Whether to extend reads to the fragment length.
    center_read (bool): Whether to center the read around the midpoint of its
        fragment.

    Returns:
    fragments (list of tuples): Each tuple is a start and end position of a
        fragment.
    """
    if not extend_reads:
        # Return all aligned blocks based on CIGAR string
        return read.get_blocks()

    #  Use median fragment length from previous analysis if no specific length
    #  provided
    default_fragment_length = fragment_length_dict.get(
        'median_fragment_length', 200
    )  # Default if not provided

    fragments = []
    if read.is_proper_pair and read.reference_id == read.next_reference_id:
        #  For paired-end reads where mates are properly paired, extend to
        #  match the fragment length
        if read.is_reverse:
            fragment_start = read.next_reference_start
            fragment_end = read.reference_end
        else:
            fragment_start = read.reference_start
            fragment_end = read.reference_start + abs(read.template_length)
        fragments.append((fragment_start, fragment_end))
    else:
        #  Handle single or improperly paired reads
        if read.is_reverse:
            fragment_start = max(
                0, read.reference_end - default_fragment_length
            )
            fragment_end = read.reference_end
        else:
            fragment_start = read.reference_start
            fragment_end = read.reference_start + default_fragment_length
        fragments.append((fragment_start, fragment_end))

    if center_read:
        #  If centering is requested, adjust fragments to be centered around
        #  the read's midpoint
        centered_fragments = []
        for start, end in fragments:
            center = start + (end - start) // 2
            half_length = default_fragment_length // 2
            centered_fragments.append((
                center - half_length, center + half_length
            ))
        return centered_fragments

    return fragments


#  Define variables ===========================================================
bam_path = 'IP_G1_Hho1_6336.sort-coord.ds.bam'
chromosome = 'I'
# region_start = 1
# region_end = 230218
region_start = 100000
region_end = 102000  # 200000
bin_size = 20
step = 20
max_fragment_length = 2000
default_fragment_length = 500


#  Run tests ==================================================================
#  Test count_reads_in_region -------------------------------------------------
counts = count_reads_in_region(
    bam_path, chromosome, region_start, region_end, step, bin_size
)
print("Testing count_reads_in_region function:")
print("    - Number of elements in array:", counts.size)
print("    - Read counts per bin:", counts)


#  Test is_proper_pair --------------------------------------------------------
#  Open the BAM file
bam_file = pysam.AlignmentFile(bam_path, "rb")

#  Fetch reads to test
reads = list(bam_file.fetch(chromosome, region_start, region_end))

print("Testing is_proper_pair function:")
for read in reads:
    print(
        f"    - Read {read.query_name} proper pair status: "
        f"{is_proper_pair(read, max_fragment_length)}"
    )

#  Test process_reads ---------------------------------------------------------
#  List of regions on chromosome I of S. cerevisiae
regions = [
    ('I', 25000, 26000),
    ('I', 50000, 51000),
    ('I', 75000, 76000),
    ('I', 100000, 101000),
    ('I', 125000, 126000),
    ('I', 150000, 151000),
    ('I', 175000, 176000),
    ('I', 200000, 201000),
    ('I', 225000, 226000)
]

print("Testing process_reads function:")
for region in regions:
    fragment_lengths, read_lengths = process_reads(bam_path, [region])
    print(f"- Region: {region}")
    print(f"- Fragment Lengths: {fragment_lengths}")
    print(f"- Read Lengths: {read_lengths}\n")


#  determine_fragment_and_read_lengths ----------------------------------------
# if __name__ == '__main__':
determine_fragment_and_read_lengths(bam_path, num_processes=1)
# determine_fragment_and_read_lengths(bam_path, num_processes=4)


import concurrent.futures
from concurrent.futures import ProcessPoolExecutor


def determine_fragment_and_read_lengths(bam_file_path, num_processes=None):
    if num_processes is None:
        num_processes = cpu_count()

    bam = pysam.AlignmentFile(bam_file_path, "rb")
    total_length = sum(bam.lengths)
    chunk_size = total_length // num_processes
    regions = []

    for contig, length in zip(bam.references, bam.lengths):
        current_pos = 0
        while current_pos + chunk_size < length:
            regions.append((contig, current_pos, current_pos + chunk_size))
            current_pos += chunk_size
        if current_pos < length:
            regions.append((contig, current_pos, length))

    results = []
    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        future_to_process = {executor.submit(process_reads, bam_file_path, region): region for region in regions}
        for future in concurrent.futures.as_completed(future_to_process):
            fl, rl = future.result()
            results.append((fl, rl))

    # Combine results from all processes
    fragment_lengths = []
    read_lengths = []
    for fl, rl in results:
        fragment_lengths.extend(fl)
        read_lengths.extend(rl)

    # Compute statistics
    fragment_read_data = {
        'mean_fragment_length': np.mean(fragment_lengths),
        'median_fragment_length': np.median(fragment_lengths),
        'mean_read_length': np.mean(read_lengths),
        'median_read_length': np.median(read_lengths)
    }

    return fragment_read_data


determine_fragment_and_read_lengths(bam_path, num_processes=1)


# #  Test get_fragment_from_read
# print("Testing get_fragment_from_read function:")
# for read in reads:
#     fragment = get_fragment_from_read(
#         read, default_fragment_length, extend_reads=True
#     )
#     print(
#         f"Read {read.query_name} fragment: "
#         f"Start {fragment[0]}, End {fragment[1]}"
#     )

# #  Usage example
# read_data = pysam.AlignedSegment()  # This would be a read object from pysam
# fragment_length_info = {'median_fragment_length': 300}  # Example median fragment length obtained from analysis
# fragments = get_fragment_from_read(read_data, fragment_length_info, extend_reads=True, center_read=False)
# print(fragments)

#  Close the BAM file
bam_file.close()


# #  Example usage:
# bam_path = 'path_to_your.bam'
# result = determine_fragment_and_read_lengths(bam_path)
# print(result)
```
</details>
<br />


<a id="test-calculate_coverage_v2py"></a>
## Test `calculate_coverage_v2.py`
<details>
<summary><i>Code: Test `calculate_coverage_v2.py`</i></summary>

```bash
#!/bin/bash

#  Run tests
dir_repo="/Users/kalavattam/repos/2023_tutorial_ChIP-seq"
bam_file_path="${dir_repo}/IP_G1_Hho1_6336.sort-coord.bam"
# output_file="${dir_repo}/IP_G1_Hho1_6336.sort-coord.test-v2.bedgraph"
# output_file="${dir_repo}/IP_G1_Hho1_6336.sort-coord.test-v3b.extend-F.bedgraph"
# output_file="${dir_repo}/IP_G1_Hho1_6336.sort-coord.test-v3c.extend-T.bedgraph"
output_file="${dir_repo}/IP_G1_Hho1_6336.sort-coord.test-v4.extend-T.bedgraph"
bin_size=30
# threads=4
threads=6
normalize=False
SAM_flag_include=64
min_MAPQ=1

# python calculate_coverage_v2.py \
#     -i "${bam_file_path}" \
#     -o "${output_file}" \
#     -b "${bin_size}" \
#     -t "${threads}"

# python calculate_coverage_v3.py \
#     -i "${bam_file_path}" \
#     -o "${output_file}" \
#     -b "${bin_size}" \
#     -t "${threads}" \
#     -s "${SAM_flag_include}" \
#     -m "${min_MAPQ}"

# python calculate_coverage_v3.py \
#     -i "${bam_file_path}" \
#     -o "${output_file}" \
#     -b "${bin_size}" \
#     -t "${threads}" \
#     -s "${SAM_flag_include}" \
#     -m "${min_MAPQ}" \
#     -e

# python calculate_coverage_v3.py \
#     -i "${bam_file_path}" \
#     -o "${output_file}" \
#     -b "${bin_size}" \
#     -t "${threads}" \
#     -s "${SAM_flag_include}" \
#     -m "${min_MAPQ}" \
#     -e \
#     -d 210 \
#     -l 840

python calculate_coverage_v4.py \
    -i "${bam_file_path}" \
    -o "${output_file}" \
    -b "${bin_size}" \
    -t "${threads}" \
    -s "${SAM_flag_include}" \
    -m "${min_MAPQ}" \
    -e \
    -d 210 \
    -l 840

dir_repo="/Users/kalavattam/repos/2023_tutorial_ChIP-seq"
bam_file_path="${dir_repo}/IP_G1_Hho1_6336.sort-coord.bam"
output_file="${bam_file_path%.bam}.test-bamCoverage-again.skipNAs.bedgraph"
chr="I"
format="bedgraph"
bin_size=30
min_MAPQ=1
threads=6
normalization="None"
flag=64

bamCoverage \
    --bam "${bam_file_path}" \
    --outFileName "${output_file}" \
    --outFileFormat "${format}" \
    --binSize "${bin_size}" \
    --numberOfProcessors "${threads}" \
    --normalizeUsing "${normalization}" \
    --minMappingQuality "${min_MAPQ}" \
    --samFlagInclude "${flag}" \
    --extendReads \
    --skipNonCoveredRegions
```
</details>
<br />