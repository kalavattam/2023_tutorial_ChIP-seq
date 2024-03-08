
`#scraps.md`
<br />
<br />

<a id="previous-work"></a>
#### Previous work
<a id="printed-8"></a>
##### Printed
<details>
<summary><i>Printed</i></summary>

```txt
❯ display_spinning_icon() {
>     what="""
>     display_spinning_icon()
>     -----------------------
>     Display \"spinning icon\" while a background process runs
> 
>     :param 1: PID of the last program the shell ran in the background (int)
>     :param 2: message to be displayed next to the spinning icon (chr)
> 
>     #TODO Checks...
>     """
>     spin="/|\\–"
>     i=0
>     while kill -0 "${1}" 2> /dev/null; do
>         i=$(( (i + 1) % 4 ))
>         printf "\r${spin:$i:1} %s" "${2}"
>         sleep .15
>     done
> }


❯ list_tally_flags() {
>     what="""
>     list_tally_flags()
>     ------------------
>     List and tally flags in a bam infile; function acts on a bam infile to
>     perform piped commands (samtools view, cut, sort, uniq -c, sort -nr) that
>     list and tally flags; function writes the results to a txt outfile, the
>     name of which is derived from the txt infile
> 
>     :param 1: name of bam infile, including path (chr)
> 
>     #TODO Checks...
>     """
>     start="$(date +%s)"
> 
>     samtools view "${1}" \
>         | cut -f 2 \
>         | sort \
>         | uniq -c \
>         | sort -nr \
>             > "${1/.bam/.flags.txt}" &
>     display_spinning_icon $! \
>     "Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on $(basename "${1}")... "
> 
>     end="$(date +%s)"
>     echo ""
>     calculate_run_time "${start}" "${end}"  \
>     "List and tally flags in $(basename "${1}")."
> }


❯ p_data="${HOME}/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams"


❯ unset bams


❯ typeset -a bams=(
>     "${p_data}/Brn1_Q_rep1_ChIP.bam"
>     "${p_data}/Brn1_Q_rep1_input.bam"
>     "${p_data}/Brn1_Q_rep2_ChIP.bam"
>     "${p_data}/Brn1_Q_rep2_input.bam"
>     "${p_data}/Brn1_Q_rep3_ChIP.bam"
>     "${p_data}/Brn1_Q_rep3_input.bam"
>     "${p_data}/Brn1_Q_all_input.bam"
>     "${p_data}/Brn1_Q_all_ChIP.bam"
>     "${p_data}/Brn1_log_rep1_ChIP.bam"
>     "${p_data}/Brn1_log_rep1_input.bam"
>     "${p_data}/Brn1_log_rep2_ChIP.bam"
>     "${p_data}/Brn1_log_rep2_input.bam"
>     "${p_data}/Brn1_log_rep3_ChIP.bam"
>     "${p_data}/Brn1_log_rep3_input.bam"
>     "${p_data}/Brn1_log_all_ChIP.bam"
>     "${p_data}/Brn1_log_all_input.bam"
> )


❯ for i in "${bams[@]}"; do echo "${i}"; done
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_Q_rep1_ChIP.bam
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_Q_rep1_input.bam
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_Q_rep2_ChIP.bam
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_Q_rep2_input.bam
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_Q_rep3_ChIP.bam
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_Q_rep3_input.bam
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_Q_all_input.bam
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_Q_all_ChIP.bam
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_log_rep1_ChIP.bam
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_log_rep1_input.bam
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_log_rep2_ChIP.bam
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_log_rep2_input.bam
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_log_rep3_ChIP.bam
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_log_rep3_input.bam
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_log_all_ChIP.bam
/home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_log_all_input.bam


❯ for i in "${bams[@]}"; do ., "${i}"; done
-rw-rw---- 1 kalavatt 1.2G May 29 12:54 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_Q_rep1_ChIP.bam
-rw-rw---- 1 kalavatt 1.8G May 29 12:55 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_Q_rep1_input.bam
-rw-rw---- 1 kalavatt 1.1G May 29 12:55 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_Q_rep2_ChIP.bam
-rw-rw---- 1 kalavatt 2.1G May 29 12:56 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_Q_rep2_input.bam
-rw-rw---- 1 kalavatt 1.6G May 29 12:57 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_Q_rep3_ChIP.bam
-rw-rw---- 1 kalavatt 1.8G May 29 12:57 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_Q_rep3_input.bam
-rw-rw---- 1 kalavatt 5.7G May 29 12:59 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_Q_all_input.bam
-rw-rw---- 1 kalavatt 3.8G May 29 13:01 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_Q_all_ChIP.bam
-rw-rw---- 1 kalavatt 1.1G May 29 13:01 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_log_rep1_ChIP.bam
-rw-rw---- 1 kalavatt 2.5G May 29 13:02 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_log_rep1_input.bam
-rw-rw---- 1 kalavatt 5.1G May 29 13:04 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_log_rep2_ChIP.bam
-rw-rw---- 1 kalavatt 2.0G May 29 13:04 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_log_rep2_input.bam
-rw-rw---- 1 kalavatt 1.4G May 29 13:05 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_log_rep3_ChIP.bam
-rw-rw---- 1 kalavatt 1.5G May 29 13:06 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_log_rep3_input.bam
-rw-rw---- 1 kalavatt 7.6G May 29 13:08 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_log_all_ChIP.bam
-rw-rw---- 1 kalavatt 5.8G May 29 13:11 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_log_all_input.bam


❯ for i in "${bams[@]}"; do
>    # i="${bams[0]}"  # echo "${i}"
>    echo "#### $(basename ${i}) ####"
>    list_tally_flags "${i}"
>    echo ""
>done
#### Brn1_Q_rep1_ChIP.bam ####
[1] 9115
| Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Brn1_Q_rep1_ChIP.bam...
– Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Brn1_Q_rep1_ChIP.bam... [1]+  Done
samtools view "${1}" | cut -d'  ' -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in Brn1_Q_rep1_ChIP.bam.
Run time: 0h:0m:6s


#### Brn1_Q_rep1_input.bam ####
[1] 9204
/ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Brn1_Q_rep1_input.bam... [1]+  Done
samtools view "${1}" | cut -d' ' -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in Brn1_Q_rep1_input.bam.
Run time: 0h:0m:10s


#### Brn1_Q_rep2_ChIP.bam ####
[1] 9299
/ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Brn1_Q_rep2_ChIP.bam... [1]+  Done
samtools view "${1}" | cut -d'  ' -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in Brn1_Q_rep2_ChIP.bam.
Run time: 0h:0m:5s


#### Brn1_Q_rep2_input.bam ####
[1] 9363
– Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Brn1_Q_rep2_input.bam... [1]+  Done
samtools view "${1}" | cut -d' ' -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in Brn1_Q_rep2_input.bam.
Run time: 0h:0m:11s


#### Brn1_Q_rep3_ChIP.bam ####
[1] 9479
– Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Brn1_Q_rep3_ChIP.bam... [1]+  Done
samtools view "${1}" | cut -d'  ' -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in Brn1_Q_rep3_ChIP.bam.
Run time: 0h:0m:8s


#### Brn1_Q_rep3_input.bam ####
[1] 9564
– Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Brn1_Q_rep3_input.bam... [1]+  Done
samtools view "${1}" | cut -d' ' -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in Brn1_Q_rep3_input.bam.
Run time: 0h:0m:10s


#### Brn1_Q_all_input.bam ####
[1] 9663
\ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Brn1_Q_all_input.bam... [1]+  Done
samtools view "${1}" | cut -d'  ' -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in Brn1_Q_all_input.bam.
Run time: 0h:0m:31s


#### Brn1_Q_all_ChIP.bam ####
[1] 10031
– Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Brn1_Q_all_ChIP.bam... [1]+  Done
samtools view "${1}" | cut -d'   ' -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in Brn1_Q_all_ChIP.bam.
Run time: 0h:0m:19s


#### Brn1_log_rep1_ChIP.bam ####
[1] 10209
– Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Brn1_log_rep1_ChIP.bam... [1]+  Done
samtools view "${1}" | cut -d'        ' -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in Brn1_log_rep1_ChIP.bam.
Run time: 0h:0m:6s


#### Brn1_log_rep1_input.bam ####
[1] 10271
/ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Brn1_log_rep1_input.bam... [1]+  Done
samtools view "${1}" | cut -d'       ' -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in Brn1_log_rep1_input.bam.
Run time: 0h:0m:13s


#### Brn1_log_rep2_ChIP.bam ####
[1] 10400
| Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Brn1_log_rep2_ChIP.bam... [1]+  Done
samtools view "${1}" | cut -d'        ' -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in Brn1_log_rep2_ChIP.bam.
Run time: 0h:0m:26s


#### Brn1_log_rep2_input.bam ####
[1] 10642
/ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Brn1_log_rep2_input.bam... [1]+  Done
samtools view "${1}" | cut -d'       ' -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in Brn1_log_rep2_input.bam.
Run time: 0h:0m:11s


#### Brn1_log_rep3_ChIP.bam ####
[1] 10794
\ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Brn1_log_rep3_ChIP.bam... [1]+  Done
samtools view "${1}" | cut -d'        ' -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in Brn1_log_rep3_ChIP.bam.
Run time: 0h:0m:7s


#### Brn1_log_rep3_input.bam ####
[1] 10868
– Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Brn1_log_rep3_input.bam... [1]+  Done
samtools view "${1}" | cut -d'       ' -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in Brn1_log_rep3_input.bam.
Run time: 0h:0m:8s


#### Brn1_log_all_ChIP.bam ####
[1] 10947
/ Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Brn1_log_all_ChIP.bam... [1]+  Done
samtools view "${1}" | cut -d' ' -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in Brn1_log_all_ChIP.bam.
Run time: 0h:0m:41s


#### Brn1_log_all_input.bam ####
[1] 11316
| Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on Brn1_log_all_input.bam... [1]+  Done
samtools view "${1}" | cut -d'        ' -f 2 | sort | uniq -c | sort -nr > "${1/.bam/.flags.txt}"

List and tally flags in Brn1_log_all_input.bam.
Run time: 0h:0m:31s


❯ for i in "${bams[@]}"; do
>    if [[ -f "${i/.bam/.flags.txt}" ]]; then
>        echo "#### $(basename ${i}) ####"
>        cat "${i/.bam/.flags.txt}"
>        echo ""
>    fi
>done
#### Brn1_Q_rep1_ChIP.bam ####
2874061 0
2862832 16
1890825 4

#### Brn1_Q_rep1_input.bam ####
5483653 0
5478220 16
 370063 4

#### Brn1_Q_rep2_ChIP.bam ####
2846780 0
2843810 16
1190484 4

#### Brn1_Q_rep2_input.bam ####
6351723 0
6342962 16
 365205 4

#### Brn1_Q_rep3_ChIP.bam ####
4442444 16
4437645 0
 780519 4

#### Brn1_Q_rep3_input.bam ####
5400769 0
5392613 16
 367131 4

#### Brn1_Q_all_input.bam ####
17235945 0
17214003 16
1102391 4

#### Brn1_Q_all_ChIP.bam ####
10156605 0
10150983 16
3861812 4

#### Brn1_log_rep1_ChIP.bam ####
3040301 16
3017984 0
 585178 4

#### Brn1_log_rep1_input.bam ####
7355182 0
7350937 16
 550871 4

#### Brn1_log_rep2_ChIP.bam ####
14450767 0
14447432 16
3176492 4

#### Brn1_log_rep2_input.bam ####
5876887 0
5862547 16
 405161 4

#### Brn1_log_rep3_ChIP.bam ####
4177590 16
4176208 0
 429674 4

#### Brn1_log_rep3_input.bam ####
4309499 16
4305931 0
 306893 4

#### Brn1_log_all_ChIP.bam ####
21661810 16
21648500 0
4191316 4

#### Brn1_log_all_input.bam ####
17538461 0
17522551 16
1262896 4


❯ samtools view "${p_data}/Brn1_log_all_ChIP.bam" | head
SRR7175373.380283       0       I       1       25      31M1I18M        *       0       0       CCACACCACACCCACACACCCACACACCACACCCACACACCACACCACAC      GGGGGIIIIIIIIIIIIIIIIIIIIIIIGIIIIIIIIIGGGGIGGGGGGI AS:i:90  XS:i:62 XN:i:0  XM:i:0  XO:i:1  XG:i:1  NM:i:1  MD:Z:49 YT:Z:UU
SRR7175373.767277       0       I       1       14      1S31M1I17M      *       0       0       CCCACACCACACCCACACACCCACACACCACACCCACACACCACACCACA      AGGAAGGGGGG<GA.<GGGAAGGIGGGG....<GGGGGA.A.<AGGGGAA AS:i:88  XS:i:82 XN:i:0  XM:i:0  XO:i:1  XG:i:1  NM:i:1  MD:Z:48 YT:Z:UU
SRR7175373.2091997      0       I       1       37      1S49M   *       0       0       CCCACACCACACCCACACACCCACACACCACACCACACACCACACCACAC      ...GA<AGAGGA<AAGGAGGGGA<AAAGGAGGIIGGAG<AGGIGAGGIGG      AS:i:98     XS:i:67 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:49 YT:Z:UU
SRR7175373.2270050      0       I       1       11      9S41M   *       0       0       ACACCACACCCACACCACACCCACACACCCACACACCACACCACACACCA      GGGGGIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIG      AS:i:82     XS:i:78 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:41 YT:Z:UU
SRR7175373.7025590      0       I       1       14      8S42M   *       0       0       CACCACACCCACACCACACCCACACACCCACACACCACACCACACACCAC      GAGGAGG<<AG.<.<<<<G<<<GGGGGGGGGGGGGGGGAGGAGAAAAAG#      AS:i:84     XS:i:78 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:42 YT:Z:UU
SRR7175373.8104795      0       I       1       21      4S46M   *       0       0       ACACCCACACCACACCCACACACCCACACACCACACCACACACCACACCA      GGGGGIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII      AS:i:92     XS:i:68 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:46 YT:Z:UU
SRR7175373.8877369      0       I       1       1       4S31M1I14M      *       0       0       ACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACC      AGA.AGGGGIAGAGGGGGGGGGGG.<.<GAGAGG.<.G<AAGGGG##### AS:i:82  XS:i:82 XN:i:0  XM:i:0  XO:i:1  XG:i:1  NM:i:1  MD:Z:45 YT:Z:UU
SRR7175373.9602499      0       I       1       14      5S45M   *       0       0       CACACCCACACCACACCCACACACCCACACACCACACCACACACCACACC      GGAGAGGGIIIIGIIIIIIIIIIIIIIIIIIIIIIIIGIIIIIIIIIIII      AS:i:90     XS:i:82 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:45 YT:Z:UU
SRR7175373.9864638      0       I       1       14      7S43M   *       0       0       ACCACACCCACACCACACCCACACACCCACACACCACACCACACACCACA      GGAAGIGIIIIIGGGIIIIIGIGGIIGGIIIGGGGIIIGIIGGGGIGIGG      AS:i:86     XS:i:80 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:43 YT:Z:UU
SRR7175373.11488836     0       I       1       11      2S31M1I16M      *       0       0       ACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCAC      AAAAGAAAGGGGAGGGGIIII.AAAAGGGGGGIG.GGGI<GG<.<AAAGA AS:i:86  XS:i:82 XN:i:0  XM:i:0  XO:i:1  XG:i:1  NM:i:1  MD:Z:47 YT:Z:UU


❯ #  Initialize functions for doing floating point arithmetic -------------------


❯ calc_6f() { awk "BEGIN{ printf \"%.6f\n\", $* }"; }


❯ calc_2f() { awk "BEGIN{ printf \"%.2f\n\", $* }"; }


❯ ml SAMtools/1.16.1-GCC-11.2.0 Bowtie2/2.4.4-GCC-11.2.0


❯ #  Get situated, initialize array of sample datasets --------------------------


❯ ml SAMtools/1.16.1-GCC-11.2.0 Bowtie2/2.4.4-GCC-11.2.0


❯ cd "${HOME}/tsukiyamalab/Kris/2023_rDNA/results" ||
>    echo "cd'ing failed; check on this..."


❯ p_data="${HOME}/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams"


❯ unset bams


❯ typeset -a bams=(
>    "${p_data}/Brn1_Q_rep1_ChIP.bam"
>    "${p_data}/Brn1_Q_rep1_input.bam"
>    "${p_data}/Brn1_Q_rep2_ChIP.bam"
>    "${p_data}/Brn1_Q_rep2_input.bam"
>    "${p_data}/Brn1_Q_rep3_ChIP.bam"
>    "${p_data}/Brn1_Q_rep3_input.bam"
>    "${p_data}/Brn1_Q_all_input.bam"
>    "${p_data}/Brn1_Q_all_ChIP.bam"
>    "${p_data}/Brn1_log_rep1_ChIP.bam"
>    "${p_data}/Brn1_log_rep1_input.bam"
>    "${p_data}/Brn1_log_rep2_ChIP.bam"
>    "${p_data}/Brn1_log_rep2_input.bam"
>    "${p_data}/Brn1_log_rep3_ChIP.bam"
>    "${p_data}/Brn1_log_rep3_input.bam"
>    "${p_data}/Brn1_log_all_ChIP.bam"
>    "${p_data}/Brn1_log_all_input.bam"
>)


❯ for i in "${bams[@]}"; do ., "${i}"; done
-rw-rw---- 1 kalavatt 1.2G May 29 12:54 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_Q_rep1_ChIP.bam
-rw-rw---- 1 kalavatt 1.8G May 29 12:55 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_Q_rep1_input.bam
-rw-rw---- 1 kalavatt 1.1G May 29 12:55 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_Q_rep2_ChIP.bam
-rw-rw---- 1 kalavatt 2.1G May 29 12:56 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_Q_rep2_input.bam
-rw-rw---- 1 kalavatt 1.6G May 29 12:57 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_Q_rep3_ChIP.bam
-rw-rw---- 1 kalavatt 1.8G May 29 12:57 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_Q_rep3_input.bam
-rw-rw---- 1 kalavatt 5.7G May 29 12:59 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_Q_all_input.bam
-rw-rw---- 1 kalavatt 3.8G May 29 13:01 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_Q_all_ChIP.bam
-rw-rw---- 1 kalavatt 1.1G May 29 13:01 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_log_rep1_ChIP.bam
-rw-rw---- 1 kalavatt 2.5G May 29 13:02 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_log_rep1_input.bam
-rw-rw---- 1 kalavatt 5.1G May 29 13:04 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_log_rep2_ChIP.bam
-rw-rw---- 1 kalavatt 2.0G May 29 13:04 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_log_rep2_input.bam
-rw-rw---- 1 kalavatt 1.4G May 29 13:05 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_log_rep3_ChIP.bam
-rw-rw---- 1 kalavatt 1.5G May 29 13:06 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_log_rep3_input.bam
-rw-rw---- 1 kalavatt 7.6G May 29 13:08 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_log_all_ChIP.bam
-rw-rw---- 1 kalavatt 5.8G May 29 13:11 /home/kalavatt/tsukiyamalab/Kris/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analysis/bams/Brn1_log_all_input.bam


❯ unset tallies


❯ typeset -A tallies


❯ for i in "${bams[@]}"; do
>    # i="${bams[0]}"  # echo "${i}"
>    echo "#### $(basename ${i}) ####"
>
>    sample="$(basename "${i}" .bam)"  # echo "${sample}"
>    total="${sample}_total"  # echo "${total}"
>    mapped="${sample}_mapped"  # echo "${mapped}"
>    unmapped="${sample}_unmapped"  # echo "${unmapped}"
>    SC="${sample}_SC"  # echo "${SC}"
>    SP="${sample}_SP"  # echo "${SP}"
>    SP2SC="${sample}_SP_to_SC"  # echo "${SP2SC}"
>    SC2SP="${sample}_SC_to_SP"  # echo "${SC2SP}"
>
>    tallies["${total}"]="$(samtools view -c "${i}")"  # echo "${tallies["${total}"]}"
>    tallies["${mapped}"]="$(samtools view -F 4 -c "${i}")"  # echo "${tallies["${mapped}"]}"
>    tallies["${unmapped}"]="$(samtools view -f 4 -c "${i}")"  # echo "${tallies["${unmapped}"]}"
>    tallies["${SC}"]="$(samtools view -F 4 -c "${i}" I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI Mito)"  # echo "${tallies["${SC}"]}"
>    tallies["${SP}"]="$(samtools view -F 4 -c "${i}" SP_II_TG SP_I SP_II SP_III SP_MTR SP_Mito)"  # echo "${tallies["${SP}"]}"
>    tallies["${SP2SC}"]="$(calc_6f "${tallies["${SP}"]}"/"${tallies["${SC}"]}")"  # echo "${tallies["${SP2SC}"]}"
>    tallies["${SC2SP}"]="$(calc_2f "${tallies["${SC}"]}"/"${tallies["${SP}"]}")"  # echo "${tallies["${SC2SP}"]}"
>
>    echo "${sample} ${tallies["${total}"]} ${tallies["${mapped}"]} ${tallies["${unmapped}"]} ${tallies["${SC}"]} ${tallies["${SP}"]} ${tallies["${SP2SC}"]} ${tallies["${SC2SP}"]}" \
>        >> "2023-0406_tutorial_ChIP-seq_analysis/bams/tmp.txt"
>
>    echo "${sample} ${tallies["${total}"]} ${tallies["${mapped}"]} ${tallies["${unmapped}"]} ${tallies["${SC}"]} ${tallies["${SP}"]} ${tallies["${SP2SC}"]} ${tallies["${SC2SP}"]}"
>    echo ""
>done
#### Brn1_Q_rep1_ChIP.bam ####
Brn1_Q_rep1_ChIP 7627718 5736893 1890825 5487671 249222 0.045415 22.02

#### Brn1_Q_rep1_input.bam ####
Brn1_Q_rep1_input 11331936 10961873 370063 10749022 212851 0.019802 50.50

#### Brn1_Q_rep2_ChIP.bam ####
Brn1_Q_rep2_ChIP 6881074 5690590 1190484 5462233 228357 0.041807 23.92

#### Brn1_Q_rep2_input.bam ####
Brn1_Q_rep2_input 13059890 12694685 365205 12405031 289654 0.023350 42.83

#### Brn1_Q_rep3_ChIP.bam ####
Brn1_Q_rep3_ChIP 9660608 8880089 780519 8797438 82651 0.009395 106.44

#### Brn1_Q_rep3_input.bam ####
Brn1_Q_rep3_input 11160513 10793382 367131 10749642 43740 0.004069 245.76

#### Brn1_Q_all_input.bam ####
Brn1_Q_all_input 35552339 34449948 1102391 33903904 546044 0.016106 62.09

#### Brn1_Q_all_ChIP.bam ####
Brn1_Q_all_ChIP 24169400 20307588 3861812 19748404 559184 0.028315 35.32

#### Brn1_log_rep1_ChIP.bam ####
Brn1_log_rep1_ChIP 6643463 6058285 585178 5906242 152043 0.025743 38.85

#### Brn1_log_rep1_input.bam ####
Brn1_log_rep1_input 15256990 14706119 550871 14390918 315201 0.021903 45.66

#### Brn1_log_rep2_ChIP.bam ####
Brn1_log_rep2_ChIP 32074691 28898199 3176492 28145866 752333 0.026730 37.41

#### Brn1_log_rep2_input.bam ####
Brn1_log_rep2_input 12144595 11739434 405161 11478196 261238 0.022759 43.94

#### Brn1_log_rep3_ChIP.bam ####
Brn1_log_rep3_ChIP 8783472 8353798 429674 8263429 90369 0.010936 91.44

#### Brn1_log_rep3_input.bam ####
Brn1_log_rep3_input 8922323 8615430 306893 8573141 42289 0.004933 202.73

#### Brn1_log_all_ChIP.bam ####
Brn1_log_all_ChIP 47501626 43310310 4191316 42314977 995333 0.023522 42.51

#### Brn1_log_all_input.bam ####
Brn1_log_all_input 36323908 35061012 1262896 34442035 618977 0.017972 55.64


❯ cat "2023-0406_tutorial_ChIP-seq_analysis/bams/tmp.txt" | sed 's/\ /\t/g' > "2023-0406_tutorial_ChIP-seq_analysis/bams/tallies.tsv"


❯ cat "2023-0406_tutorial_ChIP-seq_analysis/bams/tallies.tsv"
sample  total   mapped  unmapped    SC  SP  SP-to-SC    SC-to-SP
Brn1_Q_rep1_ChIP    7627718 5736893 1890825 5487671 249222      0.045415    22.02
Brn1_Q_rep1_input   11331936    10961873    370063  10749022    212851      0.019802    50.50
Brn1_Q_rep2_ChIP    6881074 5690590 1190484 5462233 228357      0.041807    23.92
Brn1_Q_rep2_input   13059890    12694685    365205  12405031    289654      0.023350    42.83
Brn1_Q_rep3_ChIP    9660608 8880089 780519  8797438 82651       0.009395    106.44
Brn1_Q_rep3_input   11160513    10793382    367131  10749642    43740       0.004069    245.76
Brn1_Q_all_input    35552339    34449948    1102391 33903904    546044      0.016106    62.09
Brn1_Q_all_ChIP 24169400    20307588    3861812 19748404    559184      0.028315    35.32
Brn1_log_rep1_ChIP  6643463 6058285 585178  5906242 152043      0.025743    38.85
Brn1_log_rep1_input 15256990    14706119    550871  14390918    315201      0.021903    45.66
Brn1_log_rep2_ChIP  32074691    28898199    3176492 28145866    752333      0.026730    37.41
Brn1_log_rep2_input 12144595    11739434    405161  11478196    261238      0.022759    43.94
Brn1_log_rep3_ChIP  8783472 8353798 429674  8263429 90369       0.010936    91.44
Brn1_log_rep3_input 8922323 8615430 306893  8573141 42289       0.004933    202.73
Brn1_log_all_ChIP   47501626    43310310    4191316 42314977    995333      0.023522    42.51
Brn1_log_all_input  36323908    35061012    1262896 34442035    618977      0.017972    55.64


❯ if [[ -f "2023-0406_tutorial_ChIP-seq_analysis/bams/tallies.tsv" ]]; then
>    rm "2023-0406_tutorial_ChIP-seq_analysis/bams/tmp.txt"
>fi
```
</details>
<br />

#### Installation of `R_env` on the FHCC cluster
<details>
<summary><i>Printed: </i></summary>

```txt
Environment "R_env" is not installed.
Creating environment R_env

                  __    __    __    __
                 /  \  /  \  /  \  /  \
                /    \/    \/    \/    \
███████████████/  /██/  /██/  /██/  /████████████████████████
              /  / \   / \   / \   / \  \____
             /  /   \_/   \_/   \_/   \    o \__,
            / _/                       \_____/  `
            |/
        ███╗   ███╗ █████╗ ███╗   ███╗██████╗  █████╗
        ████╗ ████║██╔══██╗████╗ ████║██╔══██╗██╔══██╗
        ██╔████╔██║███████║██╔████╔██║██████╔╝███████║
        ██║╚██╔╝██║██╔══██║██║╚██╔╝██║██╔══██╗██╔══██║
        ██║ ╚═╝ ██║██║  ██║██║ ╚═╝ ██║██████╔╝██║  ██║
        ╚═╝     ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝╚═════╝ ╚═╝  ╚═╝

        mamba (1.3.1) supported by @QuantStack

        GitHub:  https://github.com/mamba-org/mamba
        Twitter: https://twitter.com/QuantStack

█████████████████████████████████████████████████████████████


Looking for: ['bioconductor-annotationdbi', 'bioconductor-chipqc', 'bioconductor-chipseeker', 'bioconductor-clusterprofiler', 'bioconductor-deseq2', 'bioconductor-diffbind', 'bioconductor-edger', 'bioconductor-enhancedvolcano', 'bioconductor-genomicfeatures', 'bioconductor-genomicranges', 'bioconductor-ihw', 'bioconductor-iranges', 'bioconductor-pcatools', 'bioconductor-sva', 'deeptools', 'phantompeakqualtools', 'r-argparse', 'r-dendextend', 'r-devtools', 'r-ggalt', 'r-ggpubr', 'r-ggrepel', 'r-pheatmap', 'r-readxl', 'r-rjson', 'r-tidyverse', 'r-upsetr', 'r-venneuler', 'r-writexl', 'r-xml2', 'rename']

bioconda/linux-64                                           Using cache
bioconda/noarch                                             Using cache
conda-forge/linux-64                                        Using cache
conda-forge/noarch                                          Using cache
pkgs/r/noarch                                                 No change
pkgs/r/linux-64                                               No change
pkgs/main/noarch                                              No change
pkgs/main/linux-64                                            No change
Transaction

  Prefix: /home/kalavatt/miniconda3/envs/R_env

  Updating specs:

   - bioconductor-annotationdbi
   - bioconductor-chipqc
   - bioconductor-chipseeker
   - bioconductor-clusterprofiler
   - bioconductor-deseq2
   - bioconductor-diffbind
   - bioconductor-edger
   - bioconductor-enhancedvolcano
   - bioconductor-genomicfeatures
   - bioconductor-genomicranges
   - bioconductor-ihw
   - bioconductor-iranges
   - bioconductor-pcatools
   - bioconductor-sva
   - deeptools
   - phantompeakqualtools
   - r-argparse
   - r-dendextend
   - r-devtools
   - r-ggalt
   - r-ggpubr
   - r-ggrepel
   - r-pheatmap
   - r-readxl
   - r-rjson
   - r-tidyverse
   - r-upsetr
   - r-venneuler
   - r-writexl
   - r-xml2
   - rename


  Package                                                  Version  Build                Channel                    Size
──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
  Install:
──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

  + _libgcc_mutex                                              0.1  conda_forge          conda-forge/linux-64     Cached
  + _openmp_mutex                                              4.5  2_gnu                conda-forge/linux-64     Cached
  + _r-mutex                                                 1.0.1  anacondar_1          conda-forge/noarch       Cached
  + alsa-lib                                                1.2.11  hd590300_1           conda-forge/linux-64      555kB
  + argcomplete                                              3.2.2  pyhd8ed1ab_0         conda-forge/noarch         40kB
  + aws-c-auth                                              0.7.16  h79b3bcb_6           conda-forge/linux-64      103kB
  + aws-c-cal                                               0.6.10  hb29e0c7_1           conda-forge/linux-64       55kB
  + aws-c-common                                            0.9.13  hd590300_0           conda-forge/linux-64      226kB
  + aws-c-compression                                       0.2.18  hecc5fa9_1           conda-forge/linux-64       19kB
  + aws-c-event-stream                                       0.4.2  hf9b2f7b_4           conda-forge/linux-64       54kB
  + aws-c-http                                               0.8.1  h5d7533a_5           conda-forge/linux-64      195kB
  + aws-c-io                                                0.14.5  h50678d4_1           conda-forge/linux-64      157kB
  + aws-c-mqtt                                              0.10.2  hf479d2b_4           conda-forge/linux-64      164kB
  + aws-c-s3                                                 0.5.2  h4ad9680_0           conda-forge/linux-64      105kB
  + aws-c-sdkutils                                          0.1.15  hecc5fa9_1           conda-forge/linux-64       55kB
  + aws-checksums                                           0.1.18  hecc5fa9_1           conda-forge/linux-64       50kB
  + aws-crt-cpp                                             0.26.2  h19f5d62_7           conda-forge/linux-64      334kB
  + aws-sdk-cpp                                           1.11.267  h5606698_1           conda-forge/linux-64        4MB
  + binutils_impl_linux-64                                    2.40  hf600244_0           conda-forge/linux-64     Cached
  + bioconductor-annotate                                   1.78.0  r43hdfd78af_0        bioconda/noarch             2MB
  + bioconductor-annotationdbi                              1.62.2  r43hdfd78af_0        bioconda/noarch             5MB
  + bioconductor-apeglm                                     1.22.1  r43hf17093f_0        bioconda/linux-64           1MB
  + bioconductor-beachmat                                   2.16.0  r43hf17093f_0        bioconda/linux-64           1MB
  + bioconductor-biobase                                    2.60.0  r43ha9d7317_0        bioconda/linux-64           3MB
  + bioconductor-biocfilecache                               2.8.0  r43hdfd78af_0        bioconda/noarch           591kB
  + bioconductor-biocgenerics                               0.46.0  r43hdfd78af_0        bioconda/noarch           655kB
  + bioconductor-biocio                                     1.10.0  r43hdfd78af_0        bioconda/noarch           469kB
  + bioconductor-biocparallel                               1.34.2  r43hf17093f_0        bioconda/linux-64           2MB
  + bioconductor-biocsingular                               1.16.0  r43hf17093f_0        bioconda/linux-64         991kB
  + bioconductor-biomart                                    2.56.1  r43hdfd78af_0        bioconda/noarch           927kB
  + bioconductor-biostrings                                 2.68.1  r43ha9d7317_0        bioconda/linux-64          14MB
  + bioconductor-bsgenome                                   1.68.0  r43hdfd78af_0        bioconda/noarch             7MB
  + bioconductor-chipqc                                     1.36.0  r43hdfd78af_0        bioconda/noarch             2MB
  + bioconductor-chipseeker                                 1.36.0  r43hdfd78af_0        bioconda/noarch             7MB
  + bioconductor-chipseq                                    1.50.0  r43ha9d7317_0        bioconda/linux-64           3MB
  + bioconductor-clusterprofiler                             4.8.1  r43hdfd78af_0        bioconda/noarch           856kB
  + bioconductor-data-packages                            20231203  hdfd78af_0           bioconda/noarch           416kB
  + bioconductor-delayedarray                               0.26.6  r43ha9d7317_0        bioconda/linux-64           2MB
  + bioconductor-delayedmatrixstats                         1.22.1  r43hdfd78af_0        bioconda/noarch           802kB
  + bioconductor-deseq2                                     1.40.2  r43hf17093f_0        bioconda/linux-64           3MB
  + bioconductor-diffbind                                   3.10.0  r43hf17093f_0        bioconda/linux-64           7MB
  + bioconductor-dose                                       3.26.1  r43hdfd78af_0        bioconda/noarch             7MB
  + bioconductor-edger                                      3.42.4  r43hf17093f_0        bioconda/linux-64           3MB
  + bioconductor-enhancedvolcano                            1.20.0  r43hdfd78af_0        bioconda/noarch             5MB
  + bioconductor-enrichplot                                 1.20.0  r43hdfd78af_0        bioconda/noarch           368kB
  + bioconductor-fgsea                                      1.26.0  r43hf17093f_0        bioconda/linux-64           5MB
  + bioconductor-genefilter                                 1.82.1  r43ha1e849b_0        bioconda/linux-64           1MB
  + bioconductor-genomeinfodb                               1.36.1  r43hdfd78af_0        bioconda/noarch             4MB
  + bioconductor-genomeinfodbdata                           1.2.11  r43hdfd78af_1        bioconda/noarch             9kB
  + bioconductor-genomicalignments                          1.36.0  r43ha9d7317_0        bioconda/linux-64           2MB
  + bioconductor-genomicfeatures                            1.52.1  r43hdfd78af_0        bioconda/noarch             2MB
  + bioconductor-genomicranges                              1.52.0  r43ha9d7317_0        bioconda/linux-64           2MB
  + bioconductor-ggtree                                      3.8.0  r43hdfd78af_0        bioconda/noarch           905kB
  + bioconductor-go.db                                      3.17.0  r43hdfd78af_0        bioconda/noarch             9kB
  + bioconductor-gosemsim                                   2.26.0  r43hf17093f_0        bioconda/linux-64         951kB
  + bioconductor-greylistchip                               1.32.0  r43hdfd78af_0        bioconda/noarch           874kB
  + bioconductor-hdo.db                                     0.99.1  r43hdfd78af_1        bioconda/noarch             9kB
  + bioconductor-ihw                                        1.28.0  r43hdfd78af_0        bioconda/noarch             4MB
  + bioconductor-iranges                                    2.34.1  r43ha9d7317_0        bioconda/linux-64           3MB
  + bioconductor-keggrest                                   1.40.0  r43hdfd78af_0        bioconda/noarch           202kB
  + bioconductor-limma                                      3.56.2  r43ha9d7317_0        bioconda/linux-64           3MB
  + bioconductor-lpsymphony                                 1.28.1  r43hf17093f_0        bioconda/linux-64          11MB
  + bioconductor-matrixgenerics                             1.12.2  r43hdfd78af_0        bioconda/noarch           471kB
  + bioconductor-pcatools                                   2.12.0  r43hf17093f_0        bioconda/linux-64           5MB
  + bioconductor-qvalue                                     2.32.0  r43hdfd78af_0        bioconda/noarch             3MB
  + bioconductor-rhtslib                                     2.2.0  r43ha9d7317_0        bioconda/linux-64           2MB
  + bioconductor-rsamtools                                  2.16.0  r43hf17093f_0        bioconda/linux-64           4MB
  + bioconductor-rtracklayer                                1.60.0  r43ha9d7317_0        bioconda/linux-64           6MB
  + bioconductor-s4arrays                                    1.0.4  r43ha9d7317_0        bioconda/linux-64         804kB
  + bioconductor-s4vectors                                  0.38.1  r43ha9d7317_0        bioconda/linux-64           3MB
  + bioconductor-scaledmatrix                                1.8.1  r43hdfd78af_0        bioconda/noarch           659kB
  + bioconductor-shortread                                  1.58.0  r43hf17093f_0        bioconda/linux-64           6MB
  + bioconductor-sparsematrixstats                          1.12.2  r43hf17093f_0        bioconda/linux-64           1MB
  + bioconductor-summarizedexperiment                       1.30.2  r43hdfd78af_0        bioconda/noarch             2MB
  + bioconductor-sva                                        3.48.0  r43ha9d7317_0        bioconda/linux-64         501kB
  + bioconductor-systempiper                                 2.6.3  r43hdfd78af_0        bioconda/noarch             7MB
  + bioconductor-treeio                                     1.24.1  r43hdfd78af_0        bioconda/noarch           874kB
  + bioconductor-txdb.celegans.ucsc.ce6.ensgene              3.2.2  r43hdfd78af_15       bioconda/noarch            10kB
  + bioconductor-txdb.dmelanogaster.ucsc.dm3.ensgene         3.2.2  r43hdfd78af_16       bioconda/noarch            10kB
  + bioconductor-txdb.hsapiens.ucsc.hg18.knowngene           3.2.2  r43hdfd78af_15       bioconda/noarch            10kB
  + bioconductor-txdb.hsapiens.ucsc.hg19.knowngene           3.2.2  r43hdfd78af_15       bioconda/noarch            10kB
  + bioconductor-txdb.mmusculus.ucsc.mm10.knowngene         3.10.0  r43hdfd78af_8        bioconda/noarch            10kB
  + bioconductor-txdb.mmusculus.ucsc.mm9.knowngene           3.2.2  r43hdfd78af_15       bioconda/noarch            10kB
  + bioconductor-txdb.rnorvegicus.ucsc.rn4.ensgene           3.2.2  r43hdfd78af_15       bioconda/noarch            10kB
  + bioconductor-xvector                                    0.40.0  r43ha9d7317_0        bioconda/linux-64         757kB
  + bioconductor-zlibbioc                                   1.46.0  r43ha9d7317_0        bioconda/linux-64          26kB
  + boost                                                   1.84.0  h17c5347_1           conda-forge/linux-64       13kB
  + brotli                                                   1.1.0  hd590300_1           conda-forge/linux-64     Cached
  + brotli-bin                                               1.1.0  hd590300_1           conda-forge/linux-64     Cached
  + bwidget                                                 1.9.14  ha770c72_1           conda-forge/linux-64     Cached
  + bzip2                                                    1.0.8  hd590300_5           conda-forge/linux-64     Cached
  + c-ares                                                  1.27.0  hd590300_0           conda-forge/linux-64      164kB
  + ca-certificates                                       2024.2.2  hbcca054_0           conda-forge/linux-64     Cached
  + cairo                                                   1.18.0  h3faef2a_0           conda-forge/linux-64     Cached
  + certifi                                               2024.2.2  pyhd8ed1ab_0         conda-forge/noarch        161kB
  + contourpy                                                1.2.0  py310hd41b1e2_0      conda-forge/linux-64     Cached
  + curl                                                     8.5.0  hca28451_0           conda-forge/linux-64     Cached
  + cycler                                                  0.12.1  pyhd8ed1ab_0         conda-forge/noarch       Cached
  + deeptools                                                3.5.4  pyhdfd78af_1         bioconda/noarch           151kB
  + deeptoolsintervals                                       0.1.9  py310h8472f5a_7      bioconda/linux-64          78kB
  + expat                                                    2.5.0  hcb278e6_1           conda-forge/linux-64     Cached
  + font-ttf-dejavu-sans-mono                                 2.37  hab24e00_0           conda-forge/noarch       Cached
  + font-ttf-inconsolata                                     3.000  h77eed37_0           conda-forge/noarch       Cached
  + font-ttf-source-code-pro                                 2.038  h77eed37_0           conda-forge/noarch       Cached
  + font-ttf-ubuntu                                           0.83  h77eed37_1           conda-forge/noarch       Cached
  + fontconfig                                              2.14.2  h14ed4e7_0           conda-forge/linux-64     Cached
  + fonts-conda-ecosystem                                        1  0                    conda-forge/noarch       Cached
  + fonts-conda-forge                                            1  0                    conda-forge/noarch       Cached
  + fonttools                                               4.49.0  py310h2372a71_0      conda-forge/linux-64        2MB
  + freetype                                                2.12.1  h267a509_2           conda-forge/linux-64     Cached
  + fribidi                                                 1.0.10  h36c2ea0_0           conda-forge/linux-64     Cached
  + gawk                                                     5.3.0  ha916aea_0           conda-forge/linux-64        1MB
  + gcc_impl_linux-64                                       13.2.0  h338b0a0_5           conda-forge/linux-64     Cached
  + gettext                                                 0.21.1  h27087fc_0           conda-forge/linux-64     Cached
  + gflags                                                   2.2.2  he1b5a44_1004        conda-forge/linux-64     Cached
  + gfortran_impl_linux-64                                  13.2.0  h76e1118_5           conda-forge/linux-64     Cached
  + giflib                                                   5.2.1  h0b41bf4_3           conda-forge/linux-64     Cached
  + glog                                                     0.7.0  hed5481d_0           conda-forge/linux-64      144kB
  + glpk                                                       5.0  h445213a_0           conda-forge/linux-64     Cached
  + gmp                                                      6.3.0  h59595ed_0           conda-forge/linux-64     Cached
  + graphite2                                               1.3.13  h58526e2_1001        conda-forge/linux-64     Cached
  + gxx_impl_linux-64                                       13.2.0  h338b0a0_5           conda-forge/linux-64     Cached
  + harfbuzz                                                 8.3.0  h3d44ed6_0           conda-forge/linux-64     Cached
  + htslib                                                  1.19.1  h81da01d_2           bioconda/linux-64           3MB
  + icu                                                       73.2  h59595ed_0           conda-forge/linux-64     Cached
  + importlib-metadata                                       7.0.1  pyha770c72_0         conda-forge/noarch         26kB
  + jq                                                         1.5  0                    bioconda/linux-64        Cached
  + kernel-headers_linux-64                                 2.6.32  he073ed8_17          conda-forge/noarch        711kB
  + keyutils                                                 1.6.1  h166bdaf_0           conda-forge/linux-64     Cached
  + kiwisolver                                               1.4.5  py310hd41b1e2_1      conda-forge/linux-64     Cached
  + krb5                                                    1.21.2  h659d440_0           conda-forge/linux-64     Cached
  + lcms2                                                     2.15  h7f713cb_2           conda-forge/linux-64     Cached
  + ld_impl_linux-64                                          2.40  h41732ed_0           conda-forge/linux-64     Cached
  + lerc                                                     4.0.0  h27087fc_0           conda-forge/linux-64     Cached
  + libabseil                                           20240116.1  cxx17_h59595ed_2     conda-forge/linux-64        1MB
  + libarrow                                                14.0.2  h5001e6d_10_cpu      conda-forge/linux-64        8MB
  + libarrow-acero                                          14.0.2  h59595ed_10_cpu      conda-forge/linux-64      577kB
  + libarrow-dataset                                        14.0.2  h59595ed_10_cpu      conda-forge/linux-64      581kB
  + libarrow-substrait                                      14.0.2  h469e5c9_10_cpu      conda-forge/linux-64      519kB
  + libblas                                                  3.9.0  21_linux64_openblas  conda-forge/linux-64     Cached
  + libboost                                                1.84.0  h8013b2b_1           conda-forge/linux-64        3MB
  + libboost-devel                                          1.84.0  h00ab1b0_1           conda-forge/linux-64       35kB
  + libboost-headers                                        1.84.0  ha770c72_1           conda-forge/linux-64       14MB
  + libboost-python                                         1.84.0  py310hcb52e73_1      conda-forge/linux-64      119kB
  + libboost-python-devel                                   1.84.0  py310h17c5347_1      conda-forge/linux-64       17kB
  + libbrotlicommon                                          1.1.0  hd590300_1           conda-forge/linux-64     Cached
  + libbrotlidec                                             1.1.0  hd590300_1           conda-forge/linux-64     Cached
  + libbrotlienc                                             1.1.0  hd590300_1           conda-forge/linux-64     Cached
  + libcblas                                                 3.9.0  21_linux64_openblas  conda-forge/linux-64     Cached
  + libcrc32c                                                1.1.2  h9c3ff4c_0           conda-forge/linux-64     Cached
  + libcups                                                  2.3.3  h4637d8d_4           conda-forge/linux-64     Cached
  + libcurl                                                  8.5.0  hca28451_0           conda-forge/linux-64     Cached
  + libdeflate                                                1.18  h0b41bf4_0           conda-forge/linux-64     Cached
  + libedit                                           3.1.20191231  he28a2e2_2           conda-forge/linux-64     Cached
  + libev                                                     4.33  hd590300_2           conda-forge/linux-64     Cached
  + libevent                                                2.1.12  hf998b51_1           conda-forge/linux-64     Cached
  + libexpat                                                 2.5.0  hcb278e6_1           conda-forge/linux-64     Cached
  + libffi                                                   3.4.2  h7f98852_5           conda-forge/linux-64     Cached
  + libgcc                                                   7.2.0  h69d50b8_2           conda-forge/linux-64     Cached
  + libgcc-devel_linux-64                                   13.2.0  ha9c7c90_105         conda-forge/noarch       Cached
  + libgcc-ng                                               13.2.0  h807b86a_5           conda-forge/linux-64     Cached
  + libgfortran-ng                                          13.2.0  h69a702a_5           conda-forge/linux-64     Cached
  + libgfortran5                                            13.2.0  ha4646dd_5           conda-forge/linux-64     Cached
  + libgit2                                                  1.7.1  hca3a8ce_0           conda-forge/linux-64      857kB
  + libglib                                                 2.78.1  hebfc3b9_0           conda-forge/linux-64     Cached
  + libgomp                                                 13.2.0  h807b86a_5           conda-forge/linux-64     Cached
  + libgoogle-cloud                                         2.21.0  h72bcb37_2           conda-forge/linux-64        1MB
  + libgoogle-cloud-storage                                 2.21.0  hc7a4891_2           conda-forge/linux-64      750kB
  + libgrpc                                                 1.61.1  h42401df_1           conda-forge/linux-64        7MB
  + libiconv                                                  1.17  hd590300_2           conda-forge/linux-64     Cached
  + libjpeg-turbo                                          2.1.5.1  hd590300_1           conda-forge/linux-64     Cached
  + liblapack                                                3.9.0  21_linux64_openblas  conda-forge/linux-64     Cached
  + libnghttp2                                              1.58.0  h47da74e_1           conda-forge/linux-64     Cached
  + libnsl                                                   2.0.1  hd590300_0           conda-forge/linux-64     Cached
  + libopenblas                                             0.3.26  pthreads_h413a1c8_0  conda-forge/linux-64     Cached
  + libparquet                                              14.0.2  h352af49_10_cpu      conda-forge/linux-64        1MB
  + libpng                                                  1.6.43  h2797004_0           conda-forge/linux-64      288kB
  + libprotobuf                                             4.25.2  h08a7969_1           conda-forge/linux-64        3MB
  + libre2-11                                           2023.09.01  h5a48ba9_2           conda-forge/linux-64      233kB
  + libsanitizer                                            13.2.0  h7e041cc_5           conda-forge/linux-64     Cached
  + libsqlite                                               3.45.1  h2797004_0           conda-forge/linux-64     Cached
  + libssh2                                                 1.11.0  h0841786_0           conda-forge/linux-64     Cached
  + libstdcxx-devel_linux-64                                13.2.0  ha9c7c90_105         conda-forge/noarch       Cached
  + libstdcxx-ng                                            13.2.0  h7e041cc_5           conda-forge/linux-64     Cached
  + libthrift                                               0.19.0  hb90f79a_1           conda-forge/linux-64     Cached
  + libtiff                                                  4.6.0  h8b53f26_0           conda-forge/linux-64     Cached
  + libutf8proc                                              2.8.0  h166bdaf_0           conda-forge/linux-64     Cached
  + libuuid                                                 2.38.1  h0b41bf4_0           conda-forge/linux-64     Cached
  + libwebp-base                                             1.3.2  hd590300_0           conda-forge/linux-64     Cached
  + libxcb                                                    1.15  h0b41bf4_0           conda-forge/linux-64     Cached
  + libxcrypt                                               4.4.36  hd590300_1           conda-forge/linux-64     Cached
  + libxml2                                                 2.12.5  h232c23b_0           conda-forge/linux-64     Cached
  + libzlib                                                 1.2.13  hd590300_5           conda-forge/linux-64     Cached
  + lz4-c                                                    1.9.4  hcb278e6_0           conda-forge/linux-64     Cached
  + make                                                       4.3  hd18ef5c_1           conda-forge/linux-64     Cached
  + matplotlib-base                                          3.8.3  py310h62c0568_0      conda-forge/linux-64        7MB
  + mpfr                                                     4.2.1  h9458935_0           conda-forge/linux-64      642kB
  + munkres                                                  1.0.7  py_1                 bioconda/noarch          Cached
  + ncurses                                                    6.4  h59595ed_2           conda-forge/linux-64     Cached
  + nlopt                                                    2.7.1  py310h46d55b5_4      conda-forge/linux-64      398kB
  + numpy                                                   1.26.4  py310hb13e2d6_0      conda-forge/linux-64        7MB
  + openjdk                                                 20.0.2  hfea2f88_1           conda-forge/linux-64     Cached
  + openjpeg                                                 2.5.2  h488ebb8_0           conda-forge/linux-64      342kB
  + openssl                                                  3.2.1  hd590300_0           conda-forge/linux-64     Cached
  + orc                                                      1.9.2  h00e871a_2           conda-forge/linux-64        1MB
  + packaging                                                 23.2  pyhd8ed1ab_0         conda-forge/noarch       Cached
  + pandoc                                                3.1.12.2  ha770c72_0           conda-forge/linux-64       21MB
  + pango                                                  1.50.14  ha41ecd1_2           conda-forge/linux-64     Cached
  + pcre2                                                    10.40  hc3806b6_0           conda-forge/linux-64     Cached
  + perl                                                    5.32.1  7_hd590300_perl5     conda-forge/linux-64     Cached
  + phantompeakqualtools                                     1.2.2  hdfd78af_1           bioconda/noarch            14kB
  + pillow                                                  10.0.1  py310h29da1c1_1      conda-forge/linux-64     Cached
  + pip                                                       24.0  pyhd8ed1ab_0         conda-forge/noarch       Cached
  + pixman                                                  0.43.2  h59595ed_0           conda-forge/linux-64     Cached
  + plotly                                                  5.19.0  pyhd8ed1ab_0         conda-forge/noarch          6MB
  + proj                                                     9.3.1  h1d62c97_0           conda-forge/linux-64        3MB
  + pthread-stubs                                              0.4  h36c2ea0_1001        conda-forge/linux-64     Cached
  + py2bit                                                   0.3.0  py310h4b81fae_8      bioconda/linux-64          26kB
  + pybigwig                                                0.3.22  py310h79000e5_1      bioconda/linux-64        Cached
  + pyparsing                                                3.1.2  pyhd8ed1ab_0         conda-forge/noarch         89kB
  + pysam                                                   0.22.0  py310h41dec4a_1      bioconda/linux-64           5MB
  + python                                                 3.10.13  hd12c33a_1_cpython   conda-forge/linux-64     Cached
  + python-dateutil                                          2.9.0  pyhd8ed1ab_0         conda-forge/noarch        223kB
  + python_abi                                                3.10  4_cp310              conda-forge/linux-64     Cached
  + pyyaml                                                   6.0.1  py310h2372a71_1      conda-forge/linux-64     Cached
  + r-abind                                                  1.4_5  r43hc72bb7e_1005     conda-forge/noarch         78kB
  + r-amap                                                  0.8_19  r43hcf54a89_1        conda-forge/linux-64      191kB
  + r-ape                                                    5.7_1  r43h08d816e_1        conda-forge/linux-64        3MB
  + r-aplot                                                  0.2.2  r43hc72bb7e_0        conda-forge/noarch         76kB
  + r-argparse                                               2.2.2  r43hc72bb7e_1        conda-forge/noarch       Cached
  + r-arrow                                                 14.0.1  r43h59595ed_0        conda-forge/linux-64        4MB
  + r-ash                                                   1.0_15  r43h61816a4_1008     conda-forge/linux-64       47kB
  + r-ashr                                                  2.2_63  r43ha503ecb_0        conda-forge/linux-64        1MB
  + r-askpass                                                1.2.0  r43h57805ef_0        conda-forge/linux-64     Cached
  + r-assertthat                                             0.2.1  r43hc72bb7e_4        conda-forge/noarch       Cached
  + r-backports                                              1.4.1  r43h57805ef_2        conda-forge/linux-64     Cached
  + r-base                                                   4.3.1  h639d9d3_5           conda-forge/linux-64     Cached
  + r-base64enc                                              0.1_3  r43h57805ef_1006     conda-forge/linux-64     Cached
  + r-bbmle                                               1.0.25.1  r43hc72bb7e_0        conda-forge/noarch        831kB
  + r-bdsmatrix                                              1.3_7  r43h57805ef_0        conda-forge/linux-64      317kB
  + r-bh                                                  1.84.0_0  r43hc72bb7e_0        conda-forge/noarch         11MB
  + r-bit                                                    4.0.5  r43h57805ef_1        conda-forge/linux-64     Cached
  + r-bit64                                                  4.0.5  r43h57805ef_2        conda-forge/linux-64     Cached
  + r-bitops                                                 1.0_7  r43h57805ef_2        conda-forge/linux-64       44kB
  + r-blob                                                   1.2.4  r43hc72bb7e_1        conda-forge/noarch       Cached
  + r-boot                                                  1.3_30  r43hc72bb7e_0        conda-forge/noarch        628kB
  + r-brew                                                  1.0_10  r43hc72bb7e_0        conda-forge/noarch         68kB
  + r-brio                                                   1.1.4  r43h57805ef_0        conda-forge/linux-64       43kB
  + r-broom                                                  1.0.5  r43hc72bb7e_1        conda-forge/noarch       Cached
  + r-bslib                                                  0.6.1  r43hc72bb7e_0        conda-forge/noarch       Cached
  + r-cachem                                                 1.0.8  r43h57805ef_1        conda-forge/linux-64     Cached
  + r-callr                                                  3.7.5  r43hc72bb7e_0        conda-forge/noarch        421kB
  + r-car                                                    3.1_2  r43hc72bb7e_1        conda-forge/noarch          2MB
  + r-cardata                                                3.0_5  r43hc72bb7e_2        conda-forge/noarch          2MB
  + r-caret                                                 6.0_94  r43h57805ef_1        conda-forge/linux-64        4MB
  + r-catools                                               1.18.2  r43ha503ecb_2        conda-forge/linux-64      221kB
  + r-cellranger                                             1.1.0  r43hc72bb7e_1006     conda-forge/noarch       Cached
  + r-class                                                 7.3_22  r43h57805ef_1        conda-forge/linux-64      107kB
  + r-cli                                                    3.6.2  r43ha503ecb_0        conda-forge/linux-64        1MB
  + r-clipr                                                  0.8.0  r43hc72bb7e_2        conda-forge/noarch       Cached
  + r-clock                                                  0.7.0  r43ha503ecb_1        conda-forge/linux-64        2MB
  + r-coda                                                0.19_4.1  r43hc72bb7e_0        conda-forge/noarch        337kB
  + r-codetools                                             0.2_19  r43hc72bb7e_1        conda-forge/noarch        108kB
  + r-colorspace                                             2.1_0  r43h57805ef_1        conda-forge/linux-64     Cached
  + r-commonmark                                             1.9.1  r43h57805ef_0        conda-forge/linux-64      138kB
  + r-conflicted                                             1.2.0  r43h785f33e_1        conda-forge/noarch       Cached
  + r-conquer                                                1.3.3  r43h08d816e_2        conda-forge/linux-64      537kB
  + r-corrplot                                                0.92  r43hc72bb7e_2        conda-forge/noarch          4MB
  + r-cowplot                                                1.1.3  r43hc72bb7e_0        conda-forge/noarch          1MB
  + r-cpp11                                                  0.4.7  r43hc72bb7e_0        conda-forge/noarch       Cached
  + r-crayon                                                 1.5.2  r43hc72bb7e_2        conda-forge/noarch       Cached
  + r-credentials                                            2.0.1  r43hc72bb7e_0        conda-forge/noarch        226kB
  + r-crosstalk                                              1.2.1  r43hc72bb7e_0        conda-forge/noarch       Cached
  + r-curl                                                   5.1.0  r43hf9611b0_0        conda-forge/linux-64     Cached
  + r-data.table                                            1.15.2  r43h029312a_0        conda-forge/linux-64        2MB
  + r-dbi                                                    1.2.2  r43hc72bb7e_0        conda-forge/noarch        851kB
  + r-dbplyr                                                 2.4.0  r43hc72bb7e_0        conda-forge/noarch       Cached
  + r-deldir                                                 2.0_4  r43h61816a4_0        conda-forge/linux-64      285kB
  + r-dendextend                                            1.17.1  r43hc72bb7e_1        conda-forge/noarch          4MB
  + r-desc                                                   1.4.3  r43hc72bb7e_0        conda-forge/noarch        334kB
  + r-devtools                                               2.4.5  r43hc72bb7e_2        conda-forge/noarch        425kB
  + r-diagram                                                1.6.5  r43ha770c72_2        conda-forge/noarch        671kB
  + r-diffobj                                                0.3.5  r43h57805ef_2        conda-forge/linux-64      988kB
  + r-digest                                                0.6.34  r43ha503ecb_0        conda-forge/linux-64      201kB
  + r-downlit                                                0.4.3  r43hc72bb7e_0        conda-forge/noarch        120kB
  + r-downloader                                               0.4  r43hc72bb7e_1005     conda-forge/noarch         36kB
  + r-dplyr                                                  1.1.4  r43ha503ecb_0        conda-forge/linux-64     Cached
  + r-dqrng                                                  0.3.2  r43ha503ecb_0        conda-forge/linux-64      164kB
  + r-dtplyr                                                 1.3.1  r43hc72bb7e_1        conda-forge/noarch       Cached
  + r-e1071                                                 1.7_14  r43ha503ecb_0        conda-forge/linux-64      580kB
  + r-ellipsis                                               0.3.2  r43h57805ef_2        conda-forge/linux-64     Cached
  + r-emdbook                                               1.3.13  r43hc72bb7e_0        conda-forge/noarch        222kB
  + r-etrunct                                                  0.1  r43hc72bb7e_1005     conda-forge/noarch         19kB
  + r-evaluate                                                0.23  r43hc72bb7e_0        conda-forge/noarch       Cached
  + r-extrafont                                               0.19  r43ha770c72_1        conda-forge/noarch         68kB
  + r-extrafontdb                                              1.0  r43hc72bb7e_1005     conda-forge/noarch         22kB
  + r-fansi                                                  1.0.6  r43h57805ef_0        conda-forge/linux-64     Cached
  + r-farver                                                 2.1.1  r43ha503ecb_2        conda-forge/linux-64     Cached
  + r-fastmap                                                1.1.1  r43ha503ecb_1        conda-forge/linux-64     Cached
  + r-fastmatch                                              1.1_4  r43h57805ef_0        conda-forge/linux-64       48kB
  + r-fdrtool                                               1.2.17  r43h57805ef_2        conda-forge/linux-64      156kB
  + r-filelock                                               1.0.3  r43h57805ef_0        conda-forge/linux-64       33kB
  + r-findpython                                             1.0.8  r43hc72bb7e_1        conda-forge/noarch       Cached
  + r-fontawesome                                            0.5.2  r43hc72bb7e_0        conda-forge/noarch       Cached
  + r-forcats                                                1.0.0  r43hc72bb7e_1        conda-forge/noarch       Cached
  + r-foreach                                                1.5.2  r43hc72bb7e_2        conda-forge/noarch        139kB
  + r-foreign                                               0.8_86  r43h57805ef_0        conda-forge/linux-64      267kB
  + r-formatr                                                 1.14  r43hc72bb7e_1        conda-forge/noarch        165kB
  + r-fs                                                     1.6.3  r43ha503ecb_0        conda-forge/linux-64     Cached
  + r-futile.logger                                          1.4.3  r43hc72bb7e_1005     conda-forge/noarch        105kB
  + r-futile.options                                         1.0.1  r43hc72bb7e_1004     conda-forge/noarch         29kB
  + r-future                                                1.33.1  r43hc72bb7e_0        conda-forge/noarch        624kB
  + r-future.apply                                          1.11.1  r43hc72bb7e_0        conda-forge/noarch        166kB
  + r-gargle                                                 1.5.2  r43h785f33e_0        conda-forge/noarch       Cached
  + r-generics                                               0.1.3  r43hc72bb7e_2        conda-forge/noarch       Cached
  + r-gert                                                   2.0.1  r43hc25a090_0        conda-forge/linux-64      258kB
  + r-ggalt                                                  0.4.0  r43ha770c72_4        conda-forge/noarch          2MB
  + r-ggforce                                                0.4.2  r43ha503ecb_0        conda-forge/linux-64        2MB
  + r-ggfun                                                  0.1.4  r43hc72bb7e_0        conda-forge/noarch        203kB
  + r-ggnewscale                                            0.4.10  r43hc72bb7e_0        conda-forge/noarch        354kB
  + r-ggplot2                                                3.5.0  r43hc72bb7e_0        conda-forge/noarch          5MB
  + r-ggplotify                                              0.1.2  r43hc72bb7e_0        conda-forge/noarch        147kB
  + r-ggpubr                                                 0.6.0  r43hc72bb7e_1        conda-forge/noarch          2MB
  + r-ggraph                                                 2.1.0  r43ha503ecb_2        conda-forge/linux-64        4MB
  + r-ggrepel                                                0.9.5  r43ha503ecb_0        conda-forge/linux-64      271kB
  + r-ggsci                                                  3.0.1  r43hc72bb7e_0        conda-forge/noarch          2MB
  + r-ggsignif                                               0.6.4  r43hc72bb7e_1        conda-forge/noarch        572kB
  + r-gh                                                     1.4.0  r43hc72bb7e_1        conda-forge/noarch        108kB
  + r-gitcreds                                               0.1.2  r43hc72bb7e_2        conda-forge/noarch         95kB
  + r-globals                                               0.16.2  r43hc72bb7e_1        conda-forge/noarch        121kB
  + r-glue                                                   1.7.0  r43h57805ef_0        conda-forge/linux-64      155kB
  + r-googledrive                                            2.1.1  r43hc72bb7e_1        conda-forge/noarch       Cached
  + r-googlesheets4                                          1.1.1  r43h785f33e_1        conda-forge/noarch       Cached
  + r-gower                                                  1.0.1  r43h57805ef_1        conda-forge/linux-64      225kB
  + r-gplots                                               3.1.3.1  r43hc72bb7e_0        conda-forge/noarch        604kB
  + r-graphlayouts                                           1.1.0  r43ha503ecb_0        conda-forge/linux-64        4MB
  + r-gridextra                                                2.3  r43hc72bb7e_1005     conda-forge/noarch          1MB
  + r-gridgraphics                                           0.5_1  r43hc72bb7e_2        conda-forge/noarch        262kB
  + r-gson                                                   0.1.0  r43hc72bb7e_1        conda-forge/noarch        227kB
  + r-gtable                                                 0.3.4  r43hc72bb7e_0        conda-forge/noarch       Cached
  + r-gtools                                                 3.9.5  r43h57805ef_0        conda-forge/linux-64      367kB
  + r-hardhat                                                1.3.1  r43hc72bb7e_0        conda-forge/noarch        810kB
  + r-haven                                                  2.5.4  r43ha503ecb_0        conda-forge/linux-64     Cached
  + r-hexbin                                                1.28.3  r43h61816a4_1        conda-forge/linux-64     Cached
  + r-highr                                                   0.10  r43hc72bb7e_1        conda-forge/noarch       Cached
  + r-hms                                                    1.1.3  r43hc72bb7e_1        conda-forge/noarch       Cached
  + r-htmltools                                              0.5.7  r43ha503ecb_0        conda-forge/linux-64     Cached
  + r-htmlwidgets                                            1.6.4  r43hc72bb7e_1        conda-forge/noarch        425kB
  + r-httpuv                                                1.6.14  r43ha503ecb_0        conda-forge/linux-64      798kB
  + r-httr                                                   1.4.7  r43hc72bb7e_0        conda-forge/noarch       Cached
  + r-httr2                                                  1.0.0  r43hc72bb7e_0        conda-forge/noarch        539kB
  + r-hwriter                                              1.3.2.1  r43hc72bb7e_2        conda-forge/noarch        123kB
  + r-ids                                                    1.0.1  r43hc72bb7e_3        conda-forge/noarch       Cached
  + r-igraph                                                 2.0.2  r43hbec7d4a_0        conda-forge/linux-64        5MB
  + r-ini                                                    0.3.1  r43hc72bb7e_1005     conda-forge/noarch         33kB
  + r-interp                                                 1.1_6  r43hcf54a89_0        conda-forge/linux-64        1MB
  + r-invgamma                                                 1.1  r43hc72bb7e_3        conda-forge/noarch         35kB
  + r-ipred                                                 0.9_14  r43h57805ef_1        conda-forge/linux-64      392kB
  + r-irlba                                                2.3.5.1  r43h316c678_1        conda-forge/linux-64      308kB
  + r-isoband                                                0.2.7  r43ha503ecb_2        conda-forge/linux-64     Cached
  + r-iterators                                             1.0.14  r43hc72bb7e_2        conda-forge/noarch        349kB
  + r-jpeg                                                  0.1_10  r43h0de940f_3        conda-forge/linux-64       54kB
  + r-jquerylib                                              0.1.4  r43hc72bb7e_2        conda-forge/noarch       Cached
  + r-jsonlite                                               1.8.8  r43h57805ef_0        conda-forge/linux-64     Cached
  + r-kernsmooth                                           2.23_22  r43h13b3f57_0        conda-forge/linux-64      100kB
  + r-knitr                                                   1.45  r43hc72bb7e_1        conda-forge/noarch          1MB
  + r-labeling                                               0.4.3  r43hc72bb7e_0        conda-forge/noarch       Cached
  + r-lambda.r                                               1.2.4  r43hc72bb7e_3        conda-forge/noarch        119kB
  + r-later                                                  1.3.2  r43ha503ecb_0        conda-forge/linux-64     Cached
  + r-lattice                                               0.22_5  r43h57805ef_0        conda-forge/linux-64     Cached
  + r-latticeextra                                          0.6_30  r43hc72bb7e_2        conda-forge/noarch          2MB
  + r-lava                                                   1.7.3  r43hc72bb7e_0        conda-forge/noarch          3MB
  + r-lazyeval                                               0.2.2  r43h57805ef_4        conda-forge/linux-64     Cached
  + r-lifecycle                                              1.0.4  r43hc72bb7e_0        conda-forge/noarch       Cached
  + r-listenv                                                0.9.1  r43hc72bb7e_0        conda-forge/noarch        121kB
  + r-lme4                                                1.1_35.1  r43ha503ecb_0        conda-forge/linux-64        4MB
  + r-locfit                                               1.5_9.9  r43h57805ef_0        conda-forge/linux-64      548kB
  + r-lubridate                                              1.9.3  r43h57805ef_0        conda-forge/linux-64     Cached
  + r-magrittr                                               2.0.3  r43h57805ef_2        conda-forge/linux-64     Cached
  + r-maps                                                   3.4.2  r43h57805ef_0        conda-forge/linux-64        2MB
  + r-maptools                                               1.1_8  r43h57805ef_0        conda-forge/linux-64        2MB
  + r-mass                                                  7.3_60  r43h57805ef_1        conda-forge/linux-64     Cached
  + r-matrix                                                 1.6_5  r43h316c678_0        conda-forge/linux-64        4MB
  + r-matrixmodels                                           0.5_3  r43hc72bb7e_0        conda-forge/noarch        390kB
  + r-matrixstats                                            1.2.0  r43h57805ef_0        conda-forge/linux-64      466kB
  + r-memoise                                                2.0.1  r43hc72bb7e_2        conda-forge/noarch       Cached
  + r-mgcv                                                   1.9_1  r43h316c678_0        conda-forge/linux-64        3MB
  + r-mime                                                    0.12  r43h57805ef_2        conda-forge/linux-64     Cached
  + r-miniui                                               0.1.1.1  r43hc72bb7e_1004     conda-forge/noarch         54kB
  + r-minqa                                                  1.2.6  r43hcf54a89_0        conda-forge/linux-64      143kB
  + r-mixsqp                                                0.3_54  r43h08d816e_0        conda-forge/linux-64      202kB
  + r-modelmetrics                                         1.2.2.2  r43ha503ecb_3        conda-forge/linux-64      159kB
  + r-modelr                                                0.1.11  r43hc72bb7e_1        conda-forge/noarch       Cached
  + r-munsell                                                0.5.0  r43hc72bb7e_1006     conda-forge/noarch       Cached
  + r-mvtnorm                                                1.2_4  r43hd9ac46e_0        conda-forge/linux-64      740kB
  + r-nlme                                                 3.1_164  r43h61816a4_0        conda-forge/linux-64     Cached
  + r-nloptr                                                 2.0.3  r43hcf54a89_2        conda-forge/linux-64      412kB
  + r-nnet                                                  7.3_19  r43h57805ef_1        conda-forge/linux-64      131kB
  + r-nozzle.r1                                            1.1_1.1  r43ha770c72_2        conda-forge/noarch        370kB
  + r-numderiv                                          2016.8_1.1  r43hc72bb7e_5        conda-forge/noarch        127kB
  + r-openssl                                                2.1.1  r43hb353fa6_0        conda-forge/linux-64     Cached
  + r-parallelly                                            1.37.1  r43h57805ef_0        conda-forge/linux-64      370kB
  + r-patchwork                                              1.2.0  r43hc72bb7e_0        conda-forge/noarch          3MB
  + r-pbkrtest                                               0.5.2  r43hc72bb7e_1        conda-forge/noarch        195kB
  + r-pheatmap                                              1.0.12  r43hc72bb7e_4        conda-forge/noarch         91kB
  + r-pillar                                                 1.9.0  r43hc72bb7e_1        conda-forge/noarch       Cached
  + r-pkgbuild                                               1.4.2  r43hc72bb7e_0        conda-forge/noarch        205kB
  + r-pkgconfig                                              2.0.3  r43hc72bb7e_3        conda-forge/noarch       Cached
  + r-pkgdown                                                2.0.7  r43hc72bb7e_1        conda-forge/noarch        722kB
  + r-pkgload                                                1.3.4  r43hc72bb7e_0        conda-forge/noarch        196kB
  + r-plogr                                                  0.2.0  r43hc72bb7e_1005     conda-forge/noarch         22kB
  + r-plotly                                                4.10.4  r43hc72bb7e_0        conda-forge/noarch          3MB
  + r-plotrix                                                3.8_4  r43hc72bb7e_0        conda-forge/noarch          1MB
  + r-plyr                                                   1.8.9  r43ha503ecb_0        conda-forge/linux-64      824kB
  + r-png                                                    0.1_8  r43h81d01c5_1        conda-forge/linux-64       60kB
  + r-polyclip                                              1.10_6  r43ha503ecb_0        conda-forge/linux-64      124kB
  + r-polynom                                                1.4_1  r43hc72bb7e_2        conda-forge/noarch        396kB
  + r-praise                                                 1.0.0  r43hc72bb7e_1007     conda-forge/noarch         25kB
  + r-prettyunits                                            1.2.0  r43hc72bb7e_0        conda-forge/noarch       Cached
  + r-proc                                                  1.18.5  r43ha503ecb_0        conda-forge/linux-64      839kB
  + r-processx                                               3.8.3  r43h57805ef_0        conda-forge/linux-64      323kB
  + r-prodlim                                           2023.08.28  r43ha503ecb_0        conda-forge/linux-64      433kB
  + r-profvis                                                0.3.8  r43h57805ef_3        conda-forge/linux-64      206kB
  + r-progress                                               1.2.3  r43hc72bb7e_0        conda-forge/noarch       Cached
  + r-progressr                                             0.14.0  r43hc72bb7e_0        conda-forge/noarch        350kB
  + r-proj4                                                 1.0_14  r43h595fb24_0        conda-forge/linux-64       40kB
  + r-promises                                               1.2.1  r43ha503ecb_0        conda-forge/linux-64     Cached
  + r-proxy                                                 0.4_27  r43h57805ef_2        conda-forge/linux-64      185kB
  + r-ps                                                     1.7.6  r43h57805ef_0        conda-forge/linux-64      313kB
  + r-purrr                                                  1.0.2  r43h57805ef_0        conda-forge/linux-64     Cached
  + r-quantreg                                                5.97  r43hd9ac46e_0        conda-forge/linux-64        2MB
  + r-r.methodss3                                            1.8.2  r43hc72bb7e_2        conda-forge/noarch         97kB
  + r-r.oo                                                  1.26.0  r43hc72bb7e_0        conda-forge/noarch        975kB
  + r-r.utils                                               2.12.3  r43hc72bb7e_0        conda-forge/noarch          1MB
  + r-r6                                                     2.5.1  r43hc72bb7e_2        conda-forge/noarch       Cached
  + r-ragg                                                   1.2.5  r43hd759249_3        conda-forge/linux-64      438kB
  + r-rappdirs                                               0.3.3  r43h57805ef_2        conda-forge/linux-64     Cached
  + r-rcmdcheck                                              1.4.0  r43h785f33e_2        conda-forge/noarch        177kB
  + r-rcolorbrewer                                           1.1_3  r43h785f33e_2        conda-forge/noarch       Cached
  + r-rcpp                                                  1.0.12  r43h7df8631_0        conda-forge/linux-64        2MB
  + r-rcpparmadillo                                     0.12.8.1.0  r43h08d816e_0        conda-forge/linux-64      917kB
  + r-rcppeigen                                          0.3.4.0.0  r43h08d816e_0        conda-forge/linux-64        1MB
  + r-rcppnumerical                                          0.6_0  r43ha503ecb_0        conda-forge/linux-64      225kB
  + r-rcurl                                              1.98_1.14  r43hf9611b0_0        conda-forge/linux-64      826kB
  + r-readr                                                  2.1.5  r43ha503ecb_0        conda-forge/linux-64      788kB
  + r-readxl                                                 1.4.3  r43ha5c9fba_0        conda-forge/linux-64     Cached
  + r-recipes                                               1.0.10  r43hc72bb7e_0        conda-forge/noarch          2MB
  + r-rematch                                                2.0.0  r43hc72bb7e_0        conda-forge/noarch       Cached
  + r-rematch2                                               2.1.2  r43hc72bb7e_3        conda-forge/noarch       Cached
  + r-remotes                                              2.4.2.1  r43hc72bb7e_0        conda-forge/noarch        402kB
  + r-reprex                                                 2.1.0  r43hc72bb7e_0        conda-forge/noarch        499kB
  + r-reshape2                                               1.4.4  r43ha503ecb_3        conda-forge/linux-64      122kB
  + r-restfulr                                              0.0.15  r43h56115f1_3        bioconda/linux-64         441kB
  + r-rio                                                    1.0.1  r43hc72bb7e_0        conda-forge/noarch        551kB
  + r-rjava                                                 1.0_11  r43h57805ef_0        conda-forge/linux-64      782kB
  + r-rjson                                                 0.2.21  r43ha503ecb_3        conda-forge/linux-64      155kB
  + r-rlang                                                  1.1.3  r43ha503ecb_0        conda-forge/linux-64        2MB
  + r-rmarkdown                                               2.26  r43hc72bb7e_0        conda-forge/noarch          2MB
  + r-roxygen2                                               7.3.1  r43ha503ecb_0        conda-forge/linux-64      685kB
  + r-rpart                                                 4.1.23  r43h57805ef_0        conda-forge/linux-64      700kB
  + r-rprojroot                                              2.0.4  r43hc72bb7e_0        conda-forge/noarch        109kB
  + r-rsqlite                                                2.3.4  r43ha503ecb_0        conda-forge/linux-64        1MB
  + r-rstatix                                                0.7.2  r43hc72bb7e_1        conda-forge/noarch        616kB
  + r-rstudioapi                                            0.15.0  r43hc72bb7e_0        conda-forge/noarch       Cached
  + r-rsvd                                                   1.0.5  r43hc72bb7e_2        conda-forge/noarch          4MB
  + r-rttf2pt1                                              1.3.12  r43h57805ef_1        conda-forge/linux-64      109kB
  + r-rversions                                              2.1.2  r43hc72bb7e_2        conda-forge/noarch         73kB
  + r-rvest                                                  1.0.4  r43hc72bb7e_0        conda-forge/noarch        298kB
  + r-sass                                                   0.4.8  r43ha503ecb_0        conda-forge/linux-64     Cached
  + r-scales                                                 1.3.0  r43hc72bb7e_0        conda-forge/noarch       Cached
  + r-scatterpie                                             0.2.1  r43hc72bb7e_1        conda-forge/noarch        141kB
  + r-selectr                                                0.4_2  r43hc72bb7e_3        conda-forge/noarch       Cached
  + r-sessioninfo                                            1.2.2  r43hc72bb7e_2        conda-forge/noarch        198kB
  + r-shadowtext                                             0.1.3  r43hc72bb7e_0        conda-forge/noarch        221kB
  + r-shape                                                1.4.6.1  r43ha770c72_0        conda-forge/noarch        761kB
  + r-shiny                                                  1.8.0  r43h785f33e_0        conda-forge/noarch          4MB
  + r-sitmo                                                  2.0.2  r43ha503ecb_2        conda-forge/linux-64      142kB
  + r-slam                                                  0.1_50  r43h1df0287_3        conda-forge/linux-64      197kB
  + r-snow                                                   0.4_4  r43hc72bb7e_2        conda-forge/noarch        115kB
  + r-snowfall                                            1.84_6.3  r43hc72bb7e_0        conda-forge/noarch        261kB
  + r-sourcetools                                          0.1.7_1  r43ha503ecb_1        conda-forge/linux-64       54kB
  + r-sp                                                     2.1_3  r43h57805ef_0        conda-forge/linux-64        2MB
  + r-sparsem                                                 1.81  r43h61816a4_2        conda-forge/linux-64      937kB
  + r-spp                                                   1.16.0  r43h21a89ab_9        bioconda/linux-64         359kB
  + r-squarem                                               2021.1  r43hc72bb7e_2        conda-forge/noarch        195kB
  + r-statmod                                                1.5.0  r43hd8f1df9_1        conda-forge/linux-64      313kB
  + r-stringi                                                1.8.3  r43h9facbd6_0        conda-forge/linux-64      898kB
  + r-stringr                                                1.5.1  r43h785f33e_0        conda-forge/noarch       Cached
  + r-survival                                               3.5_8  r43h57805ef_0        conda-forge/linux-64        6MB
  + r-sys                                                    3.4.2  r43h57805ef_1        conda-forge/linux-64     Cached
  + r-systemfonts                                            1.0.5  r43haf97adc_0        conda-forge/linux-64     Cached
  + r-testthat                                               3.2.1  r43ha503ecb_0        conda-forge/linux-64        2MB
  + r-textshaping                                            0.3.7  r43hd87b9d6_0        conda-forge/linux-64     Cached
  + r-tibble                                                 3.2.1  r43h57805ef_2        conda-forge/linux-64     Cached
  + r-tidygraph                                              1.3.0  r43ha503ecb_0        conda-forge/linux-64      557kB
  + r-tidyr                                                  1.3.1  r43ha503ecb_0        conda-forge/linux-64        1MB
  + r-tidyselect                                             1.2.0  r43hc72bb7e_1        conda-forge/linux-64     Cached
  + r-tidytree                                               0.4.6  r43hc72bb7e_0        conda-forge/noarch        348kB
  + r-tidyverse                                              2.0.0  r43h785f33e_1        conda-forge/noarch       Cached
  + r-timechange                                             0.3.0  r43ha503ecb_0        conda-forge/linux-64      191kB
  + r-timedate                                            4032.109  r43hc72bb7e_0        conda-forge/noarch          1MB
  + r-tinytex                                                 0.49  r43hc72bb7e_1        conda-forge/noarch        147kB
  + r-truncnorm                                              1.0_9  r43h57805ef_1        conda-forge/linux-64       41kB
  + r-tweenr                                                 2.0.3  r43ha503ecb_0        conda-forge/linux-64      442kB
  + r-tzdb                                                   0.4.0  r43ha503ecb_1        conda-forge/linux-64     Cached
  + r-upsetr                                                 1.4.0  r43hc72bb7e_4        conda-forge/noarch          3MB
  + r-urlchecker                                             1.0.1  r43hc72bb7e_2        conda-forge/noarch         52kB
  + r-usethis                                                2.2.3  r43hc72bb7e_0        conda-forge/noarch        858kB
  + r-utf8                                                   1.2.4  r43h57805ef_0        conda-forge/linux-64     Cached
  + r-uuid                                                   1.2_0  r43h57805ef_0        conda-forge/linux-64       55kB
  + r-vctrs                                                  0.6.5  r43ha503ecb_0        conda-forge/linux-64     Cached
  + r-venneuler                                              1.1_4  r43hc72bb7e_0        conda-forge/noarch         62kB
  + r-viridis                                                0.6.5  r43hc72bb7e_0        conda-forge/noarch          3MB
  + r-viridislite                                            0.4.2  r43hc72bb7e_1        conda-forge/noarch       Cached
  + r-vroom                                                  1.6.5  r43ha503ecb_0        conda-forge/linux-64     Cached
  + r-waldo                                                  0.5.2  r43hc72bb7e_0        conda-forge/noarch        112kB
  + r-whisker                                                0.4.1  r43hc72bb7e_1        conda-forge/noarch         82kB
  + r-withr                                                  3.0.0  r43hc72bb7e_0        conda-forge/noarch        249kB
  + r-writexl                                                1.5.0  r43h57805ef_0        conda-forge/linux-64      140kB
  + r-xfun                                                    0.42  r43ha503ecb_0        conda-forge/linux-64      467kB
  + r-xml                                              3.99_0.16.1  r43hc6530ce_0        conda-forge/linux-64        2MB
  + r-xml2                                                   1.3.6  r43hbfba7a4_1        conda-forge/linux-64      346kB
  + r-xopen                                                  1.0.0  r43hc72bb7e_1005     conda-forge/noarch         30kB
  + r-xtable                                                 1.8_4  r43hc72bb7e_5        conda-forge/noarch        697kB
  + r-yaml                                                   2.3.8  r43h57805ef_0        conda-forge/linux-64      118kB
  + r-yulab.utils                                            0.1.4  r43hc72bb7e_0        conda-forge/noarch         87kB
  + r-zip                                                    2.3.1  r43h57805ef_0        conda-forge/linux-64      134kB
  + re2                                                 2023.09.01  h7f4b329_2           conda-forge/linux-64       27kB
  + readline                                                   8.2  h8228510_1           conda-forge/linux-64     Cached
  + rename                                                   1.601  hdfd78af_1           bioconda/noarch          Cached
  + s2n                                                      1.4.5  h06160fa_0           conda-forge/linux-64      338kB
  + samtools                                                1.19.2  h50ea8bc_1           bioconda/linux-64         463kB
  + scipy                                                   1.12.0  py310hb13e2d6_2      conda-forge/linux-64       16MB
  + sed                                                        4.8  he412f7d_0           conda-forge/linux-64     Cached
  + setuptools                                              69.1.1  pyhd8ed1ab_0         conda-forge/noarch        470kB
  + six                                                     1.16.0  pyh6c4a22f_0         conda-forge/noarch       Cached
  + snappy                                                  1.1.10  h9fff704_0           conda-forge/linux-64     Cached
  + sqlite                                                  3.45.1  h2c6b66d_0           conda-forge/linux-64      848kB
  + sysroot_linux-64                                          2.12  he073ed8_17          conda-forge/noarch         15MB
  + tenacity                                                 8.2.3  pyhd8ed1ab_0         conda-forge/noarch         23kB
  + tk                                                      8.6.13  noxft_h4845f30_101   conda-forge/linux-64     Cached
  + tktable                                                   2.10  h0c5db8f_5           conda-forge/linux-64     Cached
  + toml                                                    0.10.2  pyhd8ed1ab_0         conda-forge/noarch       Cached
  + tomlkit                                                 0.12.4  pyha770c72_0         conda-forge/noarch         37kB
  + tzdata                                                   2024a  h0c530f3_0           conda-forge/noarch       Cached
  + unicodedata2                                            15.1.0  py310h2372a71_0      conda-forge/linux-64     Cached
  + wheel                                                   0.42.0  pyhd8ed1ab_0         conda-forge/noarch       Cached
  + xmltodict                                               0.13.0  pyhd8ed1ab_0         conda-forge/noarch       Cached
  + xorg-fixesproto                                            5.0  h7f98852_1002        conda-forge/linux-64     Cached
  + xorg-inputproto                                          2.3.2  h7f98852_1002        conda-forge/linux-64     Cached
  + xorg-kbproto                                             1.0.7  h7f98852_1002        conda-forge/linux-64     Cached
  + xorg-libice                                              1.1.1  hd590300_0           conda-forge/linux-64     Cached
  + xorg-libsm                                               1.2.4  h7391055_0           conda-forge/linux-64     Cached
  + xorg-libx11                                              1.8.7  h8ee46fc_0           conda-forge/linux-64     Cached
  + xorg-libxau                                             1.0.11  hd590300_0           conda-forge/linux-64     Cached
  + xorg-libxdmcp                                            1.1.3  h7f98852_0           conda-forge/linux-64     Cached
  + xorg-libxext                                             1.3.4  h0b41bf4_2           conda-forge/linux-64     Cached
  + xorg-libxfixes                                           5.0.3  h7f98852_1004        conda-forge/linux-64     Cached
  + xorg-libxi                                              1.7.10  h7f98852_0           conda-forge/linux-64     Cached
  + xorg-libxrender                                         0.9.11  hd590300_0           conda-forge/linux-64     Cached
  + xorg-libxt                                               1.3.0  hd590300_1           conda-forge/linux-64     Cached
  + xorg-libxtst                                             1.2.3  h7f98852_1002        conda-forge/linux-64     Cached
  + xorg-recordproto                                        1.14.2  h7f98852_1002        conda-forge/linux-64     Cached
  + xorg-renderproto                                        0.11.1  h7f98852_1002        conda-forge/linux-64     Cached
  + xorg-xextproto                                           7.3.0  h0b41bf4_1003        conda-forge/linux-64     Cached
  + xorg-xproto                                             7.0.31  h7f98852_1007        conda-forge/linux-64     Cached
  + xz                                                       5.2.6  h166bdaf_0           conda-forge/linux-64     Cached
  + yaml                                                     0.2.5  h7f98852_2           conda-forge/linux-64     Cached
  + yq                                                       3.2.3  pyhd8ed1ab_0         conda-forge/noarch         23kB
  + zipp                                                    3.17.0  pyhd8ed1ab_0         conda-forge/noarch       Cached
  + zlib                                                    1.2.13  hd590300_5           conda-forge/linux-64     Cached
  + zstd                                                     1.5.5  hfc55251_0           conda-forge/linux-64     Cached

  Summary:

  Install: 572 packages

  Total download: 491MB

──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────


c-ares                                             163.6kB @ 848.9kB/s  0.2s
aws-c-common                                       225.6kB @ 895.3kB/s  0.3s
libre2-11                                          232.6kB @ 913.5kB/s  0.1s
s2n                                                337.9kB @   1.1MB/s  0.1s
aws-c-cal                                           55.1kB @ 172.6kB/s  0.1s
libabseil                                            1.3MB @   3.4MB/s  0.4s
aws-c-io                                           157.2kB @ 414.7kB/s  0.1s
aws-c-http                                         195.0kB @ 448.4kB/s  0.1s
libboost-headers                                    13.7MB @  27.2MB/s  0.5s
aws-crt-cpp                                        333.8kB @ 598.8kB/s  0.1s
gawk                                                 1.2MB @   2.1MB/s  0.3s
libgrpc                                              7.0MB @  11.6MB/s  0.3s
proj                                                 3.0MB @   4.8MB/s  0.2s
libarrow-acero                                     577.2kB @ 816.0kB/s  0.2s
libarrow-dataset                                   580.9kB @ 818.4kB/s  0.2s
r-rpart                                            700.1kB @ 869.9kB/s  0.1s
pandoc                                              21.0MB @  25.6MB/s  0.8s
htslib                                               3.0MB @   3.7MB/s  0.2s
r-gower                                            225.0kB @ 252.2kB/s  0.2s
r-foreign                                          266.9kB @ 298.6kB/s  0.1s
r-mvtnorm                                          740.2kB @ 825.4kB/s  0.1s
r-bdsmatrix                                        316.7kB @ 330.1kB/s  0.2s
r-gtools                                           367.1kB @ 378.3kB/s  0.1s
r-filelock                                          33.0kB @  31.3kB/s  0.2s
r-xml                                                1.8MB @   1.7MB/s  0.2s
r-commonmark                                       138.2kB @ 130.9kB/s  0.1s
r-ps                                               312.8kB @ 295.6kB/s  0.1s
r-proj4                                             40.2kB @  34.8kB/s  0.1s
r-uuid                                              55.2kB @  47.7kB/s  0.1s
r-kernsmooth                                       100.3kB @  81.0kB/s  0.2s
r-ash                                               46.9kB @  37.7kB/s  0.2s
sysroot_linux-64                                    15.1MB @  11.3MB/s  0.7s
r-rlang                                              1.5MB @   1.1MB/s  0.2s
r-rcpp                                               2.0MB @   1.5MB/s  0.2s
r-class                                            106.8kB @  79.6kB/s  0.1s
r-minqa                                            143.0kB @  91.1kB/s  0.4s
r-rcpparmadillo                                    916.6kB @ 558.7kB/s  0.4s
r-plyr                                             823.8kB @ 490.7kB/s  0.4s
r-maptools                                           1.8MB @   1.0MB/s  0.5s
r-proc                                             839.5kB @ 464.8kB/s  0.2s
r-interp                                             1.5MB @ 819.8kB/s  0.2s
r-rcppnumerical                                    225.3kB @ 122.8kB/s  0.2s
r-matrix                                             3.9MB @   2.0MB/s  0.7s
r-r.methodss3                                       96.9kB @  49.6kB/s  0.2s
r-timedate                                           1.3MB @ 622.9kB/s  0.2s
r-invgamma                                          35.1kB @  17.4kB/s  0.2s
r-squarem                                          194.8kB @  96.8kB/s  0.3s
r-ini                                               33.0kB @  15.6kB/s  0.2s
r-whisker                                           81.5kB @  37.5kB/s  0.2s
r-codetools                                        108.4kB @  49.7kB/s  0.2s
r-hwriter                                          123.2kB @  56.5kB/s  0.2s
r-rprojroot                                        109.3kB @  48.0kB/s  0.1s
r-remotes                                          401.9kB @ 176.2kB/s  0.2s
r-extrafontdb                                       21.6kB @   9.5kB/s  0.2s
r-polynom                                          396.2kB @ 163.3kB/s  0.3s
r-urlchecker                                        51.6kB @  21.3kB/s  0.2s
r-rversions                                         72.7kB @  29.9kB/s  0.2s
r-credentials                                      225.8kB @  92.8kB/s  0.3s
r-sessioninfo                                      198.0kB @  77.0kB/s  0.2s
r-desc                                             333.6kB @ 128.3kB/s  0.2s
r-latticeextra                                       2.2MB @ 833.0kB/s  0.2s
r-corrplot                                           3.8MB @   1.4MB/s  0.7s
r-futile.logger                                    105.3kB @  38.5kB/s  0.2s
r-r.utils                                            1.4MB @ 510.7kB/s  0.2s
deeptoolsintervals                                  77.9kB @  27.6kB/s  0.1s
r-extrafont                                         68.4kB @  24.2kB/s  0.5s
bioconductor-limma                                   2.9MB @   1.0MB/s  0.3s
plotly                                               6.1MB @   2.1MB/s  0.4s
bioconductor-biocparallel                            1.7MB @ 573.0kB/s  0.3s
r-timechange                                       191.2kB @  63.4kB/s  0.1s
r-dqrng                                            164.2kB @  53.8kB/s  0.2s
r-prodlim                                          432.6kB @ 138.4kB/s  0.1s
r-httpuv                                           797.7kB @ 254.8kB/s  0.1s
r-ipred                                            392.1kB @ 121.9kB/s  0.1s
libboost-python                                    118.6kB @  35.9kB/s  0.2s
numpy                                                7.0MB @   2.1MB/s  0.4s
r-lme4                                               4.2MB @   1.2MB/s  0.2s
r-restfulr                                         441.0kB @ 126.0kB/s  0.2s
bioconductor-biocgenerics                          655.4kB @ 186.9kB/s  0.8s
r-httr2                                            539.0kB @ 150.8kB/s  0.1s
r-pkgbuild                                         204.5kB @  55.5kB/s  0.1s
r-gh                                               107.9kB @  29.3kB/s  0.1s
bioconductor-beachmat                                1.4MB @ 355.9kB/s  0.5s
r-pkgload                                          196.1kB @  50.8kB/s  0.2s
r-usethis                                          858.0kB @ 221.7kB/s  0.2s
bioconductor-data-packages                         416.2kB @ 103.1kB/s  0.2s
bioconductor-delayedarray                            2.4MB @ 581.4kB/s  0.7s
bioconductor-delayedmatrixstats                    802.2kB @ 196.5kB/s  0.3s
scipy                                               16.5MB @   4.0MB/s  1.1s
r-profvis                                          205.6kB @  49.5kB/s  0.2s
bioconductor-scaledmatrix                          658.5kB @ 154.8kB/s  0.5s
r-pkgdown                                          722.4kB @ 168.7kB/s  0.1s
r-rvest                                            298.2kB @  67.0kB/s  0.3s
r-waldo                                            112.0kB @  25.2kB/s  0.2s
bioconductor-genomicranges                           2.3MB @ 510.7kB/s  0.5s
r-patchwork                                          3.3MB @ 716.9kB/s  0.3s
r-viridis                                            3.0MB @ 652.9kB/s  0.2s
r-aplot                                             75.7kB @  15.8kB/s  0.3s
r-plotly                                             2.9MB @ 583.9kB/s  0.5s
bioconductor-summarizedexperiment                    1.9MB @ 383.0kB/s  0.5s
bioconductor-qvalue                                  2.8MB @ 566.8kB/s  0.4s
r-ggrepel                                          271.5kB @  53.3kB/s  0.2s
bioconductor-biostrings                             14.3MB @   2.7MB/s  1.1s
r-conquer                                          536.7kB @ 103.0kB/s  0.1s
r-caret                                              3.6MB @ 682.3kB/s  0.2s
bioconductor-go.db                                   8.9kB @   1.7kB/s  0.4s
bioconductor-ggtree                                904.8kB @ 168.9kB/s  0.6s
r-scatterpie                                       141.5kB @  26.0kB/s  0.2s
bioconductor-genomicalignments                       2.4MB @ 446.6kB/s  0.3s
bioconductor-apeglm                                  1.3MB @ 228.9kB/s  0.4s
bioconductor-genomicfeatures                         2.2MB @ 385.5kB/s  0.3s
bioconductor-txdb.hsapiens.ucsc.hg19.knowngene       9.5kB @   1.6kB/s  0.1s
r-car                                                1.7MB @ 292.1kB/s  0.5s
bioconductor-txdb.celegans.ucsc.ce6.ensgene          9.5kB @   1.6kB/s  0.3s
bioconductor-rtracklayer                             5.6MB @ 924.0kB/s  0.8s
bioconductor-txdb.mmusculus.ucsc.mm10.knowngene      9.5kB @   1.6kB/s  0.3s
aws-checksums                                       50.1kB @   8.2kB/s  0.2s
libpng                                             288.2kB @  46.2kB/s  0.2s
re2                                                 26.6kB @   4.3kB/s  0.2s
bioconductor-dose                                    7.1MB @   1.1MB/s  0.9s
libgit2                                            857.2kB @ 134.6kB/s  0.1s
aws-c-event-stream                                  54.0kB @   8.5kB/s  0.2s
aws-c-mqtt                                         164.4kB @  25.4kB/s  0.2s
aws-c-auth                                         103.4kB @  15.9kB/s  0.2s
libgoogle-cloud-storage                            749.9kB @ 114.0kB/s  0.3s
samtools                                           463.4kB @  69.0kB/s  0.2s
r-parallelly                                       370.1kB @  55.1kB/s  0.2s
r-proxy                                            184.8kB @  27.3kB/s  0.2s
bioconductor-chipseeker                              7.3MB @   1.1MB/s  1.0s
r-statmod                                          312.7kB @  45.3kB/s  0.2s
r-jpeg                                              53.7kB @   7.8kB/s  0.2s
bioconductor-diffbind                                7.5MB @   1.1MB/s  1.0s
r-data.table                                         2.0MB @ 289.1kB/s  0.2s
r-writexl                                          139.9kB @  19.7kB/s  0.1s
r-catools                                          221.5kB @  31.1kB/s  0.2s
r-sparsem                                          937.4kB @ 130.7kB/s  0.5s
r-maps                                               2.4MB @ 329.5kB/s  0.4s
r-locfit                                           547.8kB @  75.2kB/s  0.2s
r-nnet                                             130.7kB @  17.9kB/s  0.2s
r-rttf2pt1                                         109.0kB @  14.7kB/s  0.6s
r-rcppeigen                                          1.5MB @ 200.6kB/s  0.2s
r-ape                                                2.9MB @ 386.5kB/s  0.3s
r-iterators                                        348.9kB @  46.0kB/s  0.1s
r-futile.options                                    28.5kB @   3.7kB/s  0.2s
r-irlba                                            308.2kB @  39.8kB/s  0.5s
r-survival                                           6.2MB @ 796.6kB/s  0.5s
r-gridgraphics                                     262.3kB @  33.8kB/s  0.3s
r-snow                                             114.8kB @  14.5kB/s  0.3s
r-dbi                                              851.2kB @ 106.7kB/s  0.2s
r-withr                                            249.2kB @  31.2kB/s  0.2s
r-coda                                             336.8kB @  42.2kB/s  0.3s
r-nozzle.r1                                        369.9kB @  45.6kB/s  0.6s
r-matrixmodels                                     389.6kB @  48.0kB/s  0.2s
r-diagram                                          671.2kB @  82.5kB/s  0.2s
r-gplots                                           604.2kB @  73.7kB/s  0.3s
r-yulab.utils                                       87.3kB @  10.5kB/s  0.2s
pyparsing                                           89.5kB @  10.8kB/s  0.2s
certifi                                            160.6kB @  19.4kB/s  0.2s
r-rsvd                                               3.6MB @ 429.1kB/s  0.5s
python-dateutil                                    222.7kB @  26.2kB/s  0.2s
importlib-metadata                                  26.4kB @   3.1kB/s  0.2s
r-future.apply                                     166.1kB @  19.5kB/s  0.2s
r-lava                                               2.5MB @ 288.5kB/s  0.4s
r-ashr                                               1.1MB @ 128.7kB/s  0.2s
r-diffobj                                          988.2kB @ 113.2kB/s  0.3s
r-processx                                         323.3kB @  37.0kB/s  0.3s
r-nloptr                                           411.6kB @  46.1kB/s  0.2s
libboost-python-devel                               16.6kB @   1.9kB/s  0.2s
r-callr                                            420.6kB @  45.9kB/s  0.2s
r-downlit                                          119.9kB @  13.1kB/s  0.2s
r-arrow                                              3.6MB @ 393.4kB/s  0.5s
r-reprex                                           499.5kB @  54.0kB/s  0.1s
bioconductor-ihw                                     4.0MB @ 434.8kB/s  0.7s
matplotlib-base                                      7.0MB @ 759.0kB/s  0.6s
r-reshape2                                         121.7kB @  13.0kB/s  0.4s
r-roxygen2                                         684.9kB @  72.2kB/s  0.4s
r-gson                                             227.1kB @  23.6kB/s  0.3s
bioconductor-biocio                                469.4kB @  48.4kB/s  0.7s
r-ggnewscale                                       354.2kB @  36.1kB/s  0.2s
bioconductor-rsamtools                               4.2MB @ 425.9kB/s  0.6s
r-ggplotify                                        146.9kB @  14.9kB/s  0.3s
r-shadowtext                                       221.4kB @  22.4kB/s  0.3s
r-dendextend                                         3.7MB @ 370.9kB/s  0.3s
bioconductor-biocfilecache                         590.6kB @  57.4kB/s  0.4s
bioconductor-annotate                                1.8MB @ 177.2kB/s  0.4s
r-upsetr                                             3.2MB @ 309.5kB/s  0.7s
bioconductor-genefilter                              1.4MB @ 132.6kB/s  0.4s
bioconductor-shortread                               5.6MB @ 527.7kB/s  0.4s
bioconductor-sva                                   501.0kB @  46.6kB/s  0.6s
bioconductor-txdb.dmelanogaster.ucsc.dm3.ensgene     9.5kB @ 878.0 B/s  0.3s
r-ggpubr                                             2.1MB @ 190.9kB/s  0.5s
bioconductor-annotationdbi                           5.2MB @ 472.0kB/s  1.1s
aws-c-sdkutils                                      55.4kB @   5.0kB/s  0.1s
sqlite                                             848.2kB @  76.5kB/s  0.1s
libboost-devel                                      35.2kB @   3.2kB/s  0.1s
bioconductor-txdb.mmusculus.ucsc.mm9.knowngene       9.5kB @ 853.0 B/s  0.4s
bioconductor-txdb.rnorvegicus.ucsc.rn4.ensgene       9.5kB @ 851.0 B/s  0.4s
openjpeg                                           341.6kB @  30.3kB/s  0.3s
bioconductor-systempiper                             6.6MB @ 576.2kB/s  0.9s
kernel-headers_linux-64                            710.6kB @  61.9kB/s  0.1s
r-sourcetools                                       54.2kB @   4.7kB/s  0.1s
r-truncnorm                                         41.0kB @   3.6kB/s  0.2s
aws-sdk-cpp                                          3.6MB @ 311.1kB/s  0.5s
r-png                                               60.0kB @   5.1kB/s  0.2s
r-matrixstats                                      465.8kB @  39.9kB/s  0.3s
r-bitops                                            44.1kB @   3.8kB/s  0.3s
r-digest                                           201.2kB @  17.2kB/s  0.3s
r-cli                                                1.3MB @ 106.2kB/s  0.3s
r-rjava                                            782.2kB @  65.8kB/s  0.3s
r-rjson                                            155.2kB @  13.1kB/s  0.3s
r-rcurl                                            826.1kB @  69.1kB/s  0.3s
r-sp                                                 1.7MB @ 142.2kB/s  0.4s
r-listenv                                          121.3kB @  10.1kB/s  0.1s
r-formatr                                          165.0kB @  13.7kB/s  0.1s
r-plogr                                             22.1kB @   1.8kB/s  0.2s
r-etrunct                                           19.4kB @   1.6kB/s  0.3s
r-praise                                            25.3kB @   2.1kB/s  0.2s
r-progressr                                        349.9kB @  28.5kB/s  0.2s
r-downloader                                        36.3kB @   3.0kB/s  0.4s
r-tinytex                                          146.8kB @  12.0kB/s  0.4s
tenacity                                            22.8kB @   1.8kB/s  0.3s
tomlkit                                             37.2kB @   3.0kB/s  0.4s
argcomplete                                         40.1kB @   3.2kB/s  0.4s
r-emdbook                                          222.3kB @  17.7kB/s  0.4s
r-future                                           623.5kB @  49.6kB/s  0.4s
r-gert                                             257.6kB @  20.2kB/s  0.3s
r-tweenr                                           441.5kB @  34.3kB/s  0.3s
r-clock                                              1.7MB @ 132.4kB/s  0.4s
fonttools                                            2.3MB @ 176.5kB/s  0.4s
r-xopen                                             30.3kB @   2.3kB/s  0.2s
yq                                                  23.4kB @   1.8kB/s  0.2s
r-rsqlite                                            1.2MB @  95.0kB/s  0.3s
r-igraph                                             5.0MB @ 380.3kB/s  0.9s
r-hardhat                                          809.6kB @  60.3kB/s  0.3s
r-spp                                              358.8kB @  26.7kB/s  0.5s
bioconductor-genomeinfodb                            4.4MB @ 324.3kB/s  0.8s
r-ggplot2                                            4.6MB @ 339.1kB/s  0.7s
r-ggsignif                                         572.1kB @  41.7kB/s  0.4s
r-ggalt                                              2.2MB @ 162.4kB/s  0.6s
bioconductor-biomart                               927.3kB @  66.9kB/s  0.6s
r-quantreg                                           1.6MB @ 112.7kB/s  0.4s
bioconductor-deseq2                                  2.8MB @ 201.8kB/s  0.4s
r-devtools                                         425.0kB @  30.0kB/s  0.2s
r-rstatix                                          615.9kB @  43.4kB/s  0.2s
mpfr                                               641.5kB @  44.8kB/s  0.4s
libgoogle-cloud                                      1.2MB @  83.0kB/s  0.3s
aws-c-s3                                           105.2kB @   7.2kB/s  0.4s
libboost                                             2.7MB @ 181.6kB/s  0.5s
r-deldir                                           285.1kB @  19.2kB/s  0.3s
r-amap                                             190.6kB @  12.8kB/s  0.4s
bioconductor-fgsea                                   4.7MB @ 317.1kB/s  1.3s
r-brio                                              43.0kB @   2.9kB/s  0.3s
r-xfun                                             466.8kB @  30.9kB/s  0.3s
r-yaml                                             118.3kB @   7.8kB/s  0.3s
bioconductor-enhancedvolcano                         5.1MB @ 337.9kB/s  1.1s
r-e1071                                            579.6kB @  37.9kB/s  0.2s
r-mixsqp                                           202.1kB @  13.2kB/s  0.2s
setuptools                                         469.6kB @  30.6kB/s  0.2s
r-globals                                          121.0kB @   7.8kB/s  0.3s
r-foreach                                          138.6kB @   8.9kB/s  0.3s
r-mgcv                                               3.3MB @ 209.6kB/s  0.6s
r-plotrix                                            1.1MB @  72.4kB/s  0.6s
r-snowfall                                         261.5kB @  16.6kB/s  0.4s
nlopt                                              398.0kB @  25.1kB/s  0.2s
bioconductor-zlibbioc                               25.5kB @   1.6kB/s  0.6s
r-ragg                                             438.3kB @  26.8kB/s  0.6s
r-rmarkdown                                          2.1MB @ 126.9kB/s  0.5s
r-rcmdcheck                                        176.7kB @  10.8kB/s  0.4s
r-graphlayouts                                       3.6MB @ 219.6kB/s  0.8s
bioconductor-genomeinfodbdata                        8.5kB @ 516.0 B/s  0.2s
r-readr                                            788.0kB @  47.6kB/s  0.3s
r-tidyr                                              1.1MB @  68.1kB/s  0.3s
r-ggfun                                            203.5kB @  12.2kB/s  0.3s
r-cowplot                                            1.3MB @  77.7kB/s  0.3s
r-pbkrtest                                         195.4kB @  11.6kB/s  0.4s
r-ggsci                                              2.1MB @ 124.4kB/s  0.5s
bioconductor-lpsymphony                             10.9MB @ 635.0kB/s  1.6s
bioconductor-gosemsim                              950.9kB @  55.0kB/s  0.6s
libprotobuf                                          2.8MB @ 162.2kB/s  0.4s
libparquet                                           1.2MB @  66.4kB/s  0.3s
bioconductor-chipseq                                 2.6MB @ 146.6kB/s  0.7s
r-polyclip                                         124.5kB @   7.0kB/s  0.3s
r-fastmatch                                         48.2kB @   2.7kB/s  0.3s
r-slam                                             197.0kB @  11.1kB/s  0.3s
bioconductor-pcatools                                5.0MB @ 280.0kB/s  1.1s
r-modelmetrics                                     159.1kB @   8.9kB/s  0.3s
r-sitmo                                            142.3kB @   7.9kB/s  0.3s
r-shape                                            761.4kB @  42.2kB/s  0.3s
r-numderiv                                         127.4kB @   7.1kB/s  0.4s
r-bbmle                                            831.2kB @  45.4kB/s  0.3s
r-knitr                                              1.3MB @  71.6kB/s  0.3s
r-gridextra                                          1.0MB @  57.3kB/s  0.4s
r-pheatmap                                          91.0kB @   5.0kB/s  0.4s
bioconductor-sparsematrixstats                       1.3MB @  67.7kB/s  0.4s
bioconductor-s4arrays                              804.3kB @  43.1kB/s  0.4s
pysam                                                4.6MB @ 245.8kB/s  0.7s
bioconductor-biobase                                 2.5MB @ 135.8kB/s  0.5s
r-recipes                                            1.6MB @  82.3kB/s  0.2s
bioconductor-keggrest                              202.0kB @  10.7kB/s  0.3s
r-rio                                              550.8kB @  29.0kB/s  0.3s
glog                                               143.6kB @   7.5kB/s  0.2s
bioconductor-iranges                                 2.6MB @ 134.8kB/s  0.8s
r-stringi                                          897.9kB @  46.4kB/s  0.3s
r-abind                                             78.3kB @   4.0kB/s  0.2s
r-fdrtool                                          156.3kB @   8.1kB/s  0.4s
bioconductor-bsgenome                                7.3MB @ 377.4kB/s  0.7s
bioconductor-chipqc                                  2.3MB @ 115.8kB/s  0.6s
r-boot                                             628.3kB @  32.1kB/s  0.3s
r-brew                                              67.5kB @   3.4kB/s  0.3s
r-r.oo                                             974.9kB @  49.6kB/s  0.3s
r-shiny                                              3.5MB @ 176.7kB/s  0.3s
boost                                               13.2kB @ 662.0 B/s  0.3s
bioconductor-s4vectors                               2.6MB @ 128.4kB/s  0.5s
bioconductor-xvector                               756.5kB @  37.4kB/s  0.6s
aws-c-compression                                   19.2kB @ 947.0 B/s  0.1s
bioconductor-matrixgenerics                        471.5kB @  23.3kB/s  0.6s
r-ggraph                                             3.8MB @ 189.1kB/s  0.4s
libarrow-substrait                                 519.5kB @  25.4kB/s  0.2s
r-zip                                              133.7kB @   6.5kB/s  0.3s
bioconductor-greylistchip                          874.4kB @  42.7kB/s  0.6s
r-xtable                                           696.9kB @  33.8kB/s  0.2s
r-tidygraph                                        557.0kB @  27.0kB/s  0.2s
r-gitcreds                                          94.8kB @   4.6kB/s  0.4s
r-venneuler                                         61.6kB @   3.0kB/s  0.4s
bioconductor-hdo.db                                  9.0kB @ 431.0 B/s  0.4s
r-ggforce                                            1.8MB @  86.6kB/s  0.5s
bioconductor-txdb.hsapiens.ucsc.hg18.knowngene       9.5kB @ 450.0 B/s  0.4s
alsa-lib                                           554.7kB @  26.0kB/s  0.4s
r-miniui                                            53.8kB @   2.5kB/s  0.2s
r-htmlwidgets                                      424.9kB @  19.7kB/s  0.2s
bioconductor-clusterprofiler                       856.3kB @  39.7kB/s  0.8s
r-bh                                                11.4MB @ 526.6kB/s  1.2s
orc                                                  1.0MB @  46.9kB/s  0.3s
r-lambda.r                                         119.4kB @   5.5kB/s  0.3s
bioconductor-rhtslib                                 2.5MB @ 113.0kB/s  0.5s
r-tidytree                                         348.4kB @  15.8kB/s  0.1s
bioconductor-enrichplot                            368.4kB @  16.8kB/s  0.4s
bioconductor-treeio                                874.5kB @  39.5kB/s  0.3s
r-glue                                             155.5kB @   7.0kB/s  0.2s
r-xml2                                             346.2kB @  15.6kB/s  0.2s
libarrow                                             8.1MB @ 362.9kB/s  0.8s
phantompeakqualtools                                14.3kB @ 639.0 B/s  0.5s
r-testthat                                           1.7MB @  74.4kB/s  0.3s
py2bit                                              25.6kB @   1.1kB/s  0.4s
r-cardata                                            1.8MB @  80.9kB/s  0.5s
bioconductor-edger                                   2.8MB @ 124.0kB/s  0.5s
deeptools                                          151.3kB @   6.7kB/s  0.3s
bioconductor-biocsingular                          990.6kB @  43.7kB/s  0.4s

Downloading and Extracting Packages

Preparing transaction: done
Verifying transaction: done
Executing transaction: done

To activate this environment, use

     $ mamba activate R_env

To deactivate an active environment, use

     $ mamba deactivate


R version 4.3.1 (2023-06-16) -- "Beagle Scouts"
Copyright (C) 2023 The R Foundation for Statistical Computing
Platform: x86_64-conda-linux-gnu (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

  Natural language support but running in an English locale

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> if (!requireNamespace("colorout", quietly = TRUE)) {
+     devtools::install_github("jalvesaq/colorout")
+ }
Downloading GitHub repo jalvesaq/colorout@HEAD
── R CMD build ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
✔  checking for file ‘/tmp/RtmpczZKwR/remotes125c74cbd5e2/jalvesaq-colorout-c6113a2/DESCRIPTION’ (550ms)
─  preparing ‘colorout’:
✔  checking DESCRIPTION meta-information
─  cleaning src
─  checking for LF line-endings in source and make files and shell scripts
─  checking for empty or unneeded directories
─  building ‘colorout_1.3-0.2.tar.gz’
   Warning: invalid uid value replaced by that for user 'nobody'
   Warning: invalid gid value replaced by that for user 'nobody'

* installing *source* package ‘colorout’ ...
** using staged installation
** libs
using C compiler: ‘x86_64-conda-linux-gnu-cc (conda-forge gcc 13.2.0-5) 13.2.0’
x86_64-conda-linux-gnu-cc -I"/home/kalavatt/miniconda3/envs/R_env/lib/R/include" -DNDEBUG   -DNDEBUG -D_FORTIFY_SOURCE=2 -O2 -isystem /home/kalavatt/miniconda3/envs/R_env/include -I/home/kalavatt/miniconda3/envs/R_env/include -Wl,-rpath-link,/home/kalavatt/miniconda3/envs/R_env/lib    -fpic  -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /home/kalavatt/miniconda3/envs/R_env/include -fdebug-prefix-map=/home/conda/feedstock_root/build_artifacts/r-base-split_1695319594214/work=/usr/local/src/conda/r-base-4.3.1 -fdebug-prefix-map=/home/kalavatt/miniconda3/envs/R_env=/usr/local/src/conda-prefix  -c colorout.c -o colorout.o
x86_64-conda-linux-gnu-cc -shared -L/home/kalavatt/miniconda3/envs/R_env/lib/R/lib -Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,--disable-new-dtags -Wl,--gc-sections -Wl,--allow-shlib-undefined -Wl,-rpath,/home/kalavatt/miniconda3/envs/R_env/lib -Wl,-rpath-link,/home/kalavatt/miniconda3/envs/R_env/lib -L/home/kalavatt/miniconda3/envs/R_env/lib -o colorout.so colorout.o -L/home/kalavatt/miniconda3/envs/R_env/lib/R/lib -lR
installing to /home/kalavatt/miniconda3/envs/R_env/lib/R/library/00LOCK-colorout/00new/colorout/libs
** R
** inst
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
*** copying figures
** building package indices
** testing if installed package can be loaded from temporary location
** checking absolute paths in shared objects and dynamic libraries
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (colorout)
```
</details>
<br />