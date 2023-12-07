
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
