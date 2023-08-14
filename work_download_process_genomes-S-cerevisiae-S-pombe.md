
`#work_download-process_genomes-S-cerevisiae-S-pombe.md`
<br />
<br />

<!-- MarkdownTOC -->

1. [Get situated](#get-situated)
    1. [Code](#code)
1. [Download *S. pombe* fastas, gff3s](#download-s-pombe-fastas-gff3s)
    1. [Code](#code-1)
    1. [Printed](#printed)
    1. [Notes](#notes)
1. [Download *S. cerevisiae* fastas, gff3](#download-s-cerevisiae-fastas-gff3)
    1. [Code](#code-2)
    1. [Printed](#printed-1)
1. [Prepare *S. pombe* fasta, gff3 for concatenation with *S. cerevisiae*](#prepare-s-pombe-fasta-gff3-for-concatenation-with-s-cerevisiae)
    1. [Code](#code-3)
    1. [Printed](#printed-2)
1. [Prepare *S. cerevisiae* fasta, gff3 for concatenation with *S. pombe*](#prepare-s-cerevisiae-fasta-gff3-for-concatenation-with-s-pombe)
    1. [Code](#code-4)
    1. [Printed](#printed-3)
1. [Concatenate processed fastas and gff3s in new directory `combined_SC_SP/`](#concatenate-processed-fastas-and-gff3s-in-new-directory-combined_sc_sp)
    1. [Code](#code-5)
1. [Create `bowtie2` indices for "`combined_SC_SP.fa.gz`"](#create-bowtie2-indices-for-combined_sc_spfagz)
    1. [Code](#code-6)
    1. [Printed](#printed-4)
1. [Copy files to Rina and Rachel](#copy-files-to-rina-and-rachel)
    1. [Get situated, copy files](#get-situated-copy-files)
        1. [Code](#code-7)
        1. [Printed](#printed-5)

<!-- /MarkdownTOC -->
<br />
<br />

<a id="get-situated"></a>
## Get situated
<a id="code"></a>
### Code
<details>
<summary><i>Code: Get situated</i></summary>

```bash
#!/bin/bash

cd "${HOME}/tsukiyamalab/kalavatt/genomes" ||
    echo "cd'ing failed; check on this..."

d_pombe="Schizosaccharomyces_pombe/"
if [[ ! -d "${d_pombe}" ]]; then mkdir "${d_pombe}"; fi

d_cerevisiae="Saccharomyces_cerevisiae/"
if [[ ! -d "${d_cerevisiae}" ]]; then mkdir "${d_cerevisiae}"; fi

alias .,="ls -lhaFG"
alias .,s="ls -lhaFG ./*"
```
</details>
<br />
<br />

<a id="download-s-pombe-fastas-gff3s"></a>
## Download *S. pombe* fastas, gff3s
<a id="code-1"></a>
### Code
<details>
<summary><i>Code: Download S. pombe fastas, gff3s</i></summary>

```bash
#!/bin/bash

cd "${HOME}/genomes/" ||
    echo "cd'ing failed; check on this..."

cd "${d_pombe}" ||
    echo "cd'ing failed; check on this..."


#  Download fastas ------------------------------------------------------------
if [[ ! -d "fasta/" ]]; then mkdir "fasta/"; fi
cd "fasta/"

u_fa="https://www.pombase.org/data/genome_sequence_and_features/genome_sequence"
f_fa=(
    Schizosaccharomyces_pombe_all_chromosomes.fa.gz
    Schizosaccharomyces_pombe_chr_II_telomeric_gap.fa.gz
    Schizosaccharomyces_pombe_chromosome_I.fa.gz
    Schizosaccharomyces_pombe_chromosome_II.fa.gz
    Schizosaccharomyces_pombe_chromosome_III.fa.gz
    Schizosaccharomyces_pombe_mating_type_region.fa.gz
    Schizosaccharomyces_pombe_mitochondrial_chromosome.fa.gz
)

#  Download the files
for i in "${f_fa[@]}"; do curl "${u_fa}/${i}" > "${i}"; done

#  Check the directory
.,

#  How do the chromosome names look?
zcat Schizosaccharomyces_pombe_all_chromosomes.fa.gz | grep "^>"


#  Download gff3s -------------------------------------------------------------
if [[ ! -d "gff3" ]]; then mkdir "gff3"; fi
cd "gff3/"

u_gff3="https://www.pombase.org/data/genome_sequence_and_features/gff3/"
f_gff3=(
    Schizosaccharomyces_pombe_all_chromosomes.gff3.gz
    Schizosaccharomyces_pombe_chr_II_telomeric_gap.gff3.gz
    Schizosaccharomyces_pombe_chromosome_I.gff3.gz
    Schizosaccharomyces_pombe_chromosome_II.gff3.gz
    Schizosaccharomyces_pombe_chromosome_III.gff3.gz
    Schizosaccharomyces_pombe_mating_type_region.gff3.gz
    Schizosaccharomyces_pombe_mitochondrial_chromosome.gff3.gz
)

#  Download the files
for i in "${f_gff3[@]}"; do curl "${u_gff3}/${i}" > "${i}"; done

#  Check the directory
.,

#  How do the chromosome names look?
zcat Schizosaccharomyces_pombe_all_chromosomes.gff3.gz \
    | cut -f 1 \
    | sort \
    | uniq
```
</details>
<br />

<a id="printed"></a>
### Printed
<details>
<summary><i>Printed: Download S. pombe fastas, gff3s</i></summary>

```txt
❯ if [[ ! -d "fasta/" ]]; then mkdir "fasta/"; fi
mkdir: created directory 'fasta/'


❯ cd "fasta/"
/home/kalavatt/tsukiyamalab/kalavatt/genomes/Schizosaccharomyces_pombe/fasta


❯ u_fa="https://www.pombase.org/data/genome_sequence_and_features/genome_sequence"


❯ f_fa=(
>    Schizosaccharomyces_pombe_all_chromosomes.fa.gz
>    Schizosaccharomyces_pombe_chr_II_telomeric_gap.fa.gz
>    Schizosaccharomyces_pombe_chromosome_I.fa.gz
>    Schizosaccharomyces_pombe_chromosome_II.fa.gz
>    Schizosaccharomyces_pombe_chromosome_III.fa.gz
>    Schizosaccharomyces_pombe_mating_type_region.fa.gz
>    Schizosaccharomyces_pombe_mitochondrial_chromosome.fa.gz
>)


❯ for i in "${f_fa[@]}"; do curl "${u_fa}/${i}" > "${i}"; done
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100 3685k  100 3685k    0     0  1848k      0  0:00:01  0:00:01 --:--:-- 1847k
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100  6568  100  6568    0     0  13295      0 --:--:-- --:--:-- --:--:-- 13295
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100 1633k  100 1633k    0     0  1031k      0  0:00:01  0:00:01 --:--:-- 1030k
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100 1325k  100 1325k    0     0   923k      0  0:00:01  0:00:01 --:--:--  923k
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100  710k  100  710k    0     0   555k      0  0:00:01  0:00:01 --:--:--  555k
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100  6563  100  6563    0     0  13285      0 --:--:-- --:--:-- --:--:-- 13285
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100  6229  100  6229    0     0  12609      0 --:--:-- --:--:-- --:--:-- 12660


❯ .,
total 9.3M
drwxrwx--- 2 kalavatt  466 May 29 08:48 ./
drwxrwx--- 3 kalavatt   23 May 29 08:47 ../
-rw-rw---- 1 kalavatt 3.6M May 29 08:48 Schizosaccharomyces_pombe_all_chromosomes.fa.gz
-rw-rw---- 1 kalavatt 6.5K May 29 08:48 Schizosaccharomyces_pombe_chr_II_telomeric_gap.fa.gz
-rw-rw---- 1 kalavatt 1.6M May 29 08:48 Schizosaccharomyces_pombe_chromosome_I.fa.gz
-rw-rw---- 1 kalavatt 1.3M May 29 08:48 Schizosaccharomyces_pombe_chromosome_II.fa.gz
-rw-rw---- 1 kalavatt 711K May 29 08:48 Schizosaccharomyces_pombe_chromosome_III.fa.gz
-rw-rw---- 1 kalavatt 6.5K May 29 08:48 Schizosaccharomyces_pombe_mating_type_region.fa.gz
-rw-rw---- 1 kalavatt 6.1K May 29 08:48 Schizosaccharomyces_pombe_mitochondrial_chromosome.fa.gz


❯ zcat Schizosaccharomyces_pombe_all_chromosomes.fa.gz | grep "^>"
>chr_II_telomeric_gap Schizosaccharomyces_pombe
>I Schizosaccharomyces_pombe
>II Schizosaccharomyces_pombe
>III Schizosaccharomyces_pombe
>mating_type_region Schizosaccharomyces_pombe
>mitochondrial Schizosaccharomyces_pombe


❯ if [[ ! -d "gff3" ]]; then mkdir "gff3"; fi
mkdir: created directory 'gff3'


❯ cd "gff3/"
/home/kalavatt/tsukiyamalab/kalavatt/genomes/Schizosaccharomyces_pombe/gff3


❯ u_gff3="https://www.pombase.org/data/genome_sequence_and_features/gff3/"


❯ f_gff3=(
>    Schizosaccharomyces_pombe_all_chromosomes.gff3.gz
>    Schizosaccharomyces_pombe_chr_II_telomeric_gap.gff3.gz
>    Schizosaccharomyces_pombe_chromosome_I.gff3.gz
>    Schizosaccharomyces_pombe_chromosome_II.gff3.gz
>    Schizosaccharomyces_pombe_chromosome_III.gff3.gz
>    Schizosaccharomyces_pombe_mating_type_region.gff3.gz
>    Schizosaccharomyces_pombe_mitochondrial_chromosome.gff3.gz
>)


❯ for i in "${f_gff3[@]}"; do curl "${u_gff3}/${i}" > "${i}"; done
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100  600k  100  600k    0     0   393k      0  0:00:01  0:00:01 --:--:--  393k
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100   297  100   297    0     0    602      0 --:--:-- --:--:-- --:--:--   602
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100  266k  100  266k    0     0   238k      0  0:00:01  0:00:01 --:--:--  237k
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100  220k  100  220k    0     0   197k      0  0:00:01  0:00:01 --:--:--  197k
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100  113k  100  113k    0     0   118k      0 --:--:-- --:--:-- --:--:--  118k
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100   399  100   399    0     0    810      0 --:--:-- --:--:-- --:--:--   810
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100  1478  100  1478    0     0   2991      0 --:--:-- --:--:-- --:--:--  2991


❯ .,
total 2.6M
drwxrwx--- 2 kalavatt  480 May 29 09:02 ./
drwxrwx--- 4 kalavatt   45 May 29 09:01 ../
-rw-rw---- 1 kalavatt 601K May 29 09:02 Schizosaccharomyces_pombe_all_chromosomes.gff3.gz
-rw-rw---- 1 kalavatt  297 May 29 09:02 Schizosaccharomyces_pombe_chr_II_telomeric_gap.gff3.gz
-rw-rw---- 1 kalavatt 267K May 29 09:02 Schizosaccharomyces_pombe_chromosome_I.gff3.gz
-rw-rw---- 1 kalavatt 221K May 29 09:02 Schizosaccharomyces_pombe_chromosome_II.gff3.gz
-rw-rw---- 1 kalavatt 114K May 29 09:02 Schizosaccharomyces_pombe_chromosome_III.gff3.gz
-rw-rw---- 1 kalavatt  399 May 29 09:02 Schizosaccharomyces_pombe_mating_type_region.gff3.gz
-rw-rw---- 1 kalavatt 1.5K May 29 09:02 Schizosaccharomyces_pombe_mitochondrial_chromosome.gff3.gz


❯ zcat Schizosaccharomyces_pombe_all_chromosomes.gff3.gz \
>    | cut -f 1 \
>    | sort \
>    | uniq
chr_II_telomeric_gap
##gff-version 3
I
II
III
mating_type_region
mitochondrial
```
</details>
<br />

<a id="notes"></a>
### Notes
<details>
<summary><i>Notes: Download S. pombe fastas, gff3s</i></summary>

Fasta file information as of 2023-0529, the date of downloading (no `README` in directory):
- Schizosaccharomyces_pombe_all_chromosomes.fa.gz 2023-05-11 02:56 3.6M
- Schizosaccharomyces_pombe_chr_II_telomeric_gap.fa.gz 2023-05-11 02:56 6.4K
- Schizosaccharomyces_pombe_chromosome_I.fa.gz 2023-05-11 02:56 1.6M
- Schizosaccharomyces_pombe_chromosome_II.fa.gz 2023-05-11 02:56 1.3M
- Schizosaccharomyces_pombe_chromosome_III.fa.gz 2023-05-11 02:56 711K
- Schizosaccharomyces_pombe_mating_type_region.fa.gz 2023-05-11 02:56 6.4K
- Schizosaccharomyces_pombe_mitochondrial_chromosome.fa.gz 2023-05-11 02:56 6.1K

Gff3 file information as of 2023-0529, the date of downloading (no `README` in directory):
- Schizosaccharomyces_pombe_all_chromosomes.gff3.gz 2023-05-28 02:12 601K
- Schizosaccharomyces_pombe_chr_II_telomeric_gap.gff3.gz 2023-05-28 02:12 297
- Schizosaccharomyces_pombe_chromosome_I.gff3.gz 2023-05-28 02:12 267K
- Schizosaccharomyces_pombe_chromosome_II.gff3.gz 2023-05-28 02:12 220K
- Schizosaccharomyces_pombe_chromosome_III.gff3.gz 2023-05-28 02:12 114K
- Schizosaccharomyces_pombe_mating_type_region.gff3.gz 2023-05-28 02:12 401
- Schizosaccharomyces_pombe_mitochondrial_chromosome.gff3.gz 2023-05-28 02:12 1.4K
</details>
<br />
<br />

<a id="download-s-cerevisiae-fastas-gff3"></a>
## Download *S. cerevisiae* fastas, gff3
<a id="code-2"></a>
### Code
<details>
<summary><i>Code: Download S. cerevisiae fastas, gff3</i></summary>

```bash
#!/bin/bash

cd "${HOME}/genomes/" ||
    echo "cd'ing failed; check on this..."

cd "${d_cerevisiae}" ||
    echo "cd'ing failed; check on this..."


#  Download fastas, gff3 ------------------------------------------------------
u_tgz="http://sgd-archive.yeastgenome.org/sequence/S288C_reference/genome_releases"
f_tgz="S288C_reference_genome_R64-3-1_20210421.tgz"

#  Download the tarfile of fastas and gff3
curl "${u_tgz}/${f_tgz}" > "${f_tgz}"

#  Decompress the directory
tar -xvzf "${f_tgz}"

#  Rename the directory to "fasta/"
mv S288C_reference_genome_R64-3-1_20210421/ fasta/

#  Create a "gff3/" directory and move the gff3 file from "fasta/" to "gff3/"
mkdir gff3/ && \
    mv fasta/saccharomyces_cerevisiae_R64-3-1_20210421.gff.gz gff3/

#  Check how everything looks
.,s
.,
```
</details>
<br />

<a id="printed-1"></a>
### Printed
<details>
<summary><i>Printed: Download S. cerevisiae fastas, gff3</i></summary>

```txt
❯ cd "${d_cerevisiae}" ||
>    echo "cd'ing failed; check on this..."
/home/kalavatt/genomes/Saccharomyces_cerevisiae


❯ u_tgz="http://sgd-archive.yeastgenome.org/sequence/S288C_reference/genome_releases"


❯ f_tgz="S288C_reference_genome_R64-3-1_20210421.tgz"


❯ curl "${u_tgz}/${f_tgz}" > "${f_tgz}"
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100 19.6M  100 19.6M    0     0  28.2M      0 --:--:-- --:--:-- --:--:-- 28.2M


❯ tar -xvzf "${f_tgz}"
S288C_reference_genome_R64-3-1_20210421/
S288C_reference_genome_R64-3-1_20210421/orf_coding_all_R64-3-1_20210421.fasta.gz
S288C_reference_genome_R64-3-1_20210421/other_features_genomic_R64-3-1_20210421.fasta.gz
S288C_reference_genome_R64-3-1_20210421/orf_trans_all_R64-3-1_20210421.fasta.gz
S288C_reference_genome_R64-3-1_20210421/S288C_reference_sequence_R64-3-1_20210421.fsa.gz
S288C_reference_genome_R64-3-1_20210421/NotFeature_R64-3-1_20210421.fasta.gz
S288C_reference_genome_R64-3-1_20210421/gene_association_R64-3-1_20210421.sgd.gz
S288C_reference_genome_R64-3-1_20210421/saccharomyces_cerevisiae_R64-3-1_20210421.gff.gz
S288C_reference_genome_R64-3-1_20210421/rna_coding_R64-3-1_20210421.fasta.gz


❯ mv S288C_reference_genome_R64-3-1_20210421/ fasta/
renamed 'S288C_reference_genome_R64-3-1_20210421/' -> 'fasta/'


❯ mkdir gff3/ && \
>    mv fasta/saccharomyces_cerevisiae_R64-3-1_20210421.gff.gz gff3/
mkdir: created directory 'gff3/'
renamed 'fasta/saccharomyces_cerevisiae_R64-3-1_20210421.gff.gz' -> 'gff3/saccharomyces_cerevisiae_R64-3-1_20210421.gff.gz'


❯ .,s
-rw-rw---- 1 kalavatt 20M May 29 08:53 ./S288C_reference_genome_R64-3-1_20210421.tgz

./fasta:
total 19M
drwxr-x--- 2 kalavatt  413 May 29 08:54 ./
drwxrwx--- 4 kalavatt  106 May 29 08:54 ../
-rw-r----- 1 kalavatt 3.6M Apr 27  2021 gene_association_R64-3-1_20210421.sgd.gz
-rw-r----- 1 kalavatt 1.1M Apr 21  2021 NotFeature_R64-3-1_20210421.fasta.gz
-rw-r----- 1 kalavatt 3.7M Apr 21  2021 orf_coding_all_R64-3-1_20210421.fasta.gz
-rw-r----- 1 kalavatt 2.6M Apr 21  2021 orf_trans_all_R64-3-1_20210421.fasta.gz
-rw-r----- 1 kalavatt 187K Apr 21  2021 other_features_genomic_R64-3-1_20210421.fasta.gz
-rw-r----- 1 kalavatt  42K Apr 27  2021 rna_coding_R64-3-1_20210421.fasta.gz
-rw-r----- 1 kalavatt 3.7M Apr 21  2021 S288C_reference_sequence_R64-3-1_20210421.fsa.gz

./gff3:
total 6.0M
drwxrwx--- 2 kalavatt   66 May 29 08:54 ./
drwxrwx--- 4 kalavatt  106 May 29 08:54 ../
-rw-r----- 1 kalavatt 5.1M Apr 27  2021 saccharomyces_cerevisiae_R64-3-1_20210421.gff.gz


❯ .,
total 23M
drwxrwx---  4 kalavatt 106 May 29 08:54 ./
drwxrwx--- 16 kalavatt 578 May 29 08:46 ../
drwxr-x---  2 kalavatt 413 May 29 08:54 fasta/
drwxrwx---  2 kalavatt  66 May 29 08:54 gff3/
-rw-rw----  1 kalavatt 20M May 29 08:53 S288C_reference_genome_R64-3-1_20210421.tgz
```
</details>
<br />
<br />

<a id="prepare-s-pombe-fasta-gff3-for-concatenation-with-s-cerevisiae"></a>
## Prepare *S. pombe* fasta, gff3 for concatenation with *S. cerevisiae*
<a id="code-3"></a>
### Code
<details>
<summary><i>Code: Prepare S. pombe fasta, gff3 for concatenation with S. cerevisiae</i></summary>

```bash
#!/bin/bash

cd "${HOME}/genomes/" ||
    echo "cd'ing failed; check on this..."

cd "${d_pombe}" ||
    echo "cd'ing failed; check on this..."


#  Process fasta --------------------------------------------------------------
#  Create a directory for the processed fasta
if [[ ! -d "fasta-processed" ]]; then mkdir "fasta-processed"; fi

#  Copy over file to process, then cd into directory
cp \
    "fasta/Schizosaccharomyces_pombe_all_chromosomes.fa.gz" \
    "fasta-processed/"

cd "fasta-processed" ||
    echo "cd'ing failed; check on this..."

#  What do chromosome names look like?
zgrep "^>" "Schizosaccharomyces_pombe_all_chromosomes.fa.gz"

#  Create a decompressed version of the fasta
gzip -cd "Schizosaccharomyces_pombe_all_chromosomes.fa.gz" \
    > "Schizosaccharomyces_pombe_all_chromosomes.fa"

#  Use sed to rename chromosomes; save the results to "tmp.fa"
if [[ -f "tmp.fa" ]]; then rm "tmp.fa"; fi
sed 's/^>chr_II_telomeric_gap\ Schizosaccharomyces_pombe/>SP_II_TG/g;s/^>I\ Schizosaccharomyces_pombe/>SP_I/g;s/^>II\ Schizosaccharomyces_pombe/>SP_II/g;s/^>III\ Schizosaccharomyces_pombe/>SP_III/g;s/^>mating_type_region\ Schizosaccharomyces_pombe/>SP_MTR/g;s/^>mitochondrial\ Schizosaccharomyces_pombe/>SP_Mito/g' "Schizosaccharomyces_pombe_all_chromosomes.fa" \
    > "tmp.fa"

#  Check the new chromosome names
cat "tmp.fa" | grep "^>"

#  Overwrite the initial file with the contents of "tmp.fa"
mv -f "tmp.fa" "Schizosaccharomyces_pombe_all_chromosomes.fa"

#  Double check chromosome names
cat "Schizosaccharomyces_pombe_all_chromosomes.fa" | grep "^>"

#  Remove the compressed initial file
rm *.gz

#  Compress the updated file (with new chromosome names)
gzip *.fa

#  Check the directory
.,

#  Check chromosome names again
zcat "Schizosaccharomyces_pombe_all_chromosomes.fa.gz" | grep "^>"


#  Process gff3 ---------------------------------------------------------------
cd .. && pwd

#  Create a directory for the processed gff3
if [[ ! -d "gff3-processed/" ]]; then mkdir "gff3-processed/"; fi

#  Copy over file to process, then cd into directory
cp "gff3/Schizosaccharomyces_pombe_all_chromosomes.gff3.gz" "gff3-processed/"

cd "gff3-processed/" ||
    echo "cd'ing failed; check on this..."

#  What do chromosome names look like?
zcat "Schizosaccharomyces_pombe_all_chromosomes.gff3.gz" \
    | cut -f 1 \
    | sort \
    | uniq

#  Use sed to rename chromosomes; save the results to "tmp.gff3"
zcat "Schizosaccharomyces_pombe_all_chromosomes.gff3.gz" \
    | sed 's/^chr_II_telomeric_gap/SP_II_TG/g;s/^I/SP_I/g;s/^II/SP_II/g;s/^III/SP_III/g;s/^mating_type_region/SP_MTR/g;s/^mitochondrial/SP_Mito/g' \
        > "tmp.gff3"

#  Check on the file contents
head "tmp.gff3"
tail "tmp.gff3"
zcat "../gff3/Schizosaccharomyces_pombe_all_chromosomes.gff3.gz" | tail

#  Check on the chromosome names in "tmp.gff3"
cat "tmp.gff3" \
    | cut -f 1 \
    | sort \
    | uniq

#  Rename "tmp.gff3" to "Schizosaccharomyces_pombe_all_chromosomes.gff3"
mv "tmp.gff3" "Schizosaccharomyces_pombe_all_chromosomes.gff3"

#  Remove the initial gzipped file
rm "Schizosaccharomyces_pombe_all_chromosomes.gff3.gz"

#  Compress "Schizosaccharomyces_pombe_all_chromosomes.gff3"
gzip "Schizosaccharomyces_pombe_all_chromosomes.gff3"

#  Check the directory contents
.,

#  Check chromosome names again
zcat "Schizosaccharomyces_pombe_all_chromosomes.gff3.gz" \
    | cut -f 1 \
    | sort \
    | uniq
```
</details>
<br />

<a id="printed-2"></a>
### Printed
<details>
<summary><i>Printed: Prepare S. pombe fasta, gff3 for concatenation with S. cerevisiae</i></summary>

```txt
❯ cd "${d_pombe}" ||
>    echo "cd'ing failed; check on this..."
/home/kalavatt/genomes/Schizosaccharomyces_pombe


❯ if [[ ! -d "fasta-processed" ]]; then mkdir "fasta-processed"; fi
mkdir: created directory 'fasta-processed'


❯ cp \
>    "fasta/Schizosaccharomyces_pombe_all_chromosomes.fa.gz" \
>    "fasta-processed/"
'fasta/Schizosaccharomyces_pombe_all_chromosomes.fa.gz' -> 'fasta-processed/Schizosaccharomyces_pombe_all_chromosomes.fa.gz'


❯ cd "fasta-processed" ||
>    echo "cd'ing failed; check on this..."
/home/kalavatt/genomes/Schizosaccharomyces_pombe/fasta-processed


❯ zgrep "^>" "Schizosaccharomyces_pombe_all_chromosomes.fa.gz"
>chr_II_telomeric_gap Schizosaccharomyces_pombe
>I Schizosaccharomyces_pombe
>II Schizosaccharomyces_pombe
>III Schizosaccharomyces_pombe
>mating_type_region Schizosaccharomyces_pombe
>mitochondrial Schizosaccharomyces_pombe


❯ gzip -cd "Schizosaccharomyces_pombe_all_chromosomes.fa.gz" \
>    > "Schizosaccharomyces_pombe_all_chromosomes.fa"


❯ if [[ -f "tmp.fa" ]]; then rm "tmp.fa"; fi


❯ sed 's/^>chr_II_telomeric_gap\ Schizosaccharomyces_pombe/>SP_II_TG/g;s/^>I\ Schizosaccharomyces_pombe/>SP_I/g;s/^>II\ Schizosaccharomyces_pombe/>SP_II/g;s/^>III\ Schizosaccharomyces_pombe/>SP_III/g;s/^>mating_type_region\ Schizosaccharomyces_pombe/>SP_MTR/g;s/^>mitochondrial\ Schizosaccharomyces_pombe/SP_Mito/g' "Schizosaccharomyces_pombe_all_chromosomes.fa" \
>    > "tmp.fa"


❯ cat "tmp.fa" | grep "^>"
>SP_II_TG
>SP_I
>SP_II
>SP_III
>SP_MTR
>SP_Mito


❯ mv -f "tmp.fa" "Schizosaccharomyces_pombe_all_chromosomes.fa"
renamed 'tmp.fa' -> 'Schizosaccharomyces_pombe_all_chromosomes.fa'


❯ cat "Schizosaccharomyces_pombe_all_chromosomes.fa" | grep "^>"
>SP_II_TG
>SP_I
>SP_II
>SP_III
>SP_MTR
>SP_Mito


❯ rm *.gz


❯ gzip *.fa


❯ .,
total 4.6M
drwxrwx--- 2 kalavatt   65 May 29 09:25 ./
drwxrwx--- 5 kalavatt   78 May 29 09:12 ../
-rw-rw---- 1 kalavatt 3.8M May 29 09:21 Schizosaccharomyces_pombe_all_chromosomes.fa.gz


❯ zcat "Schizosaccharomyces_pombe_all_chromosomes.fa.gz" | grep "^>"
>SP_II_TG
>SP_I
>SP_II
>SP_III
>SP_MTR
>SP_Mito


❯ cd .. && pwd
/home/kalavatt/genomes/Schizosaccharomyces_pombe


❯ if [[ ! -d "gff3-processed/" ]]; then mkdir "gff3-processed/"; fi
mkdir: created directory 'gff3-processed/'


❯ cp "gff3/Schizosaccharomyces_pombe_all_chromosomes.gff3.gz" "gff3-processed/"
'gff3/Schizosaccharomyces_pombe_all_chromosomes.gff3.gz' -> 'gff3-processed/Schizosaccharomyces_pombe_all_chromosomes.gff3.gz'


❯ cd "gff3-processed/" ||
>    echo "cd'ing failed; check on this..."
/home/kalavatt/genomes/Schizosaccharomyces_pombe/gff3-processed


❯ zcat "Schizosaccharomyces_pombe_all_chromosomes.gff3.gz" \
>    | cut -f 1 \
>    | sort \
>    | uniq
chr_II_telomeric_gap
##gff-version 3
I
II
III
mating_type_region
mitochondrial


❯ zcat "Schizosaccharomyces_pombe_all_chromosomes.gff3.gz" \
>    | sed 's/^chr_II_telomeric_gap/SP_II_TG/g;s/^I/SP_I/g;s/^II/SP_II/g;s/^III/SP_III/g;s/^mating_type_region/SP_MTR/g;s/^mitochondrial/SP_Mito/g' \
>        > "tmp.gff3"


❯ head "tmp.gff3"
##gff-version 3
SP_I    PomBase gene    1798347 1798835 .   +   .   ID=SPAC1002.01;Name=mrx11
SP_I    PomBase mRNA    1798347 1798835 .   +   .   ID=SPAC1002.01.1;Parent=SPAC1002.01
SP_I    PomBase CDS 1798347 1798835 .   +   0   ID=SPAC1002.01.1:exon:1;Parent=SPAC1002.01.1
SP_I    PomBase gene    1799014 1800053 .   +   .   ID=SPAC1002.02;Name=pom34
SP_I    PomBase mRNA    1799014 1800053 .   +   .   ID=SPAC1002.02.1;Parent=SPAC1002.02
SP_I    PomBase five_prime_UTR  1799014 1799127 .   +   .   ID=SPAC1002.02.1:five_prime_UTR:1;Parent=SPAC1002.02.1
SP_I    PomBase CDS 1799128 1799817 .   +   0   ID=SPAC1002.02.1:exon:1;Parent=SPAC1002.02.1
SP_I    PomBase three_prime_UTR 1799818 1800053 .   +   .   ID=SPAC1002.02.1:three_prime_UTR:1;Parent=SPAC1002.02.1
SP_I    PomBase gene    1799915 1803070 .   -   .   ID=SPAC1002.03c;Name=gls2


❯ tail "tmp.gff3"
SP_I    PomBase TR_box  4514798 4514815 .   -   .   ID=CU329670_TR_box_4514798..4514815
SP_III  PomBase long_terminal_repeat    2319921 2320269 .   -   .   ID=SPLTRC.71
SP_II   PomBase long_terminal_repeat    2339945 2340297 .   -   .   ID=SPLTRB.34
SP_I    PomBase long_terminal_repeat    4525577 4525926 .   +   .   ID=SPLTRA.71
SP_III  PomBase gene_group  2450422 2452883 .   +   .   ID=CU329672_gene_group_2450422..2452883
SP_III  PomBase long_terminal_repeat    2108288 2108631 .   -   .   ID=SPLTRC.58
SP_III  PomBase long_terminal_repeat    782301  782649  .   -   .   ID=SPLTRC.27
SP_III  PomBase region  814799  814896  .   -   .   ID=CU329672_region_814799..814896
SP_III  PomBase dh_repeat   1087567 1091508 .   +   .   ID=SPRPTCENC.9
SP_I    PomBase long_terminal_repeat    32863   33057   .   -   .   ID=SPLTRA.6


❯ zcat "../gff3/Schizosaccharomyces_pombe_all_chromosomes.gff3.gz" | tail
I   PomBase TR_box  4514798 4514815 .   -   .   ID=CU329670_TR_box_4514798..4514815
III PomBase long_terminal_repeat    2319921 2320269 .   -   .   ID=SPLTRC.71
II  PomBase long_terminal_repeat    2339945 2340297 .   -   .   ID=SPLTRB.34
I   PomBase long_terminal_repeat    4525577 4525926 .   +   .   ID=SPLTRA.71
III PomBase gene_group  2450422 2452883 .   +   .   ID=CU329672_gene_group_2450422..2452883
III PomBase long_terminal_repeat    2108288 2108631 .   -   .   ID=SPLTRC.58
III PomBase long_terminal_repeat    782301  782649  .   -   .   ID=SPLTRC.27
III PomBase region  814799  814896  .   -   .   ID=CU329672_region_814799..814896
III PomBase dh_repeat   1087567 1091508 .   +   .   ID=SPRPTCENC.9
I   PomBase long_terminal_repeat    32863   33057   .   -   .   ID=SPLTRA.6


❯ cat "tmp.gff3" \
>    | cut -f 1 \
>    | sort \
>    | uniq
##gff-version 3
SP_I
SP_II
SP_III
SP_II_TG
SP_Mito
SP_MTR


❯ mv "tmp.gff3" "Schizosaccharomyces_pombe_all_chromosomes.gff3"
renamed 'tmp.gff3' -> 'Schizosaccharomyces_pombe_all_chromosomes.gff3'


❯ rm "Schizosaccharomyces_pombe_all_chromosomes.gff3.gz"


❯ gzip "Schizosaccharomyces_pombe_all_chromosomes.gff3"


❯ .,
total 984K
drwxrwx--- 2 kalavatt   67 May 29 09:39 ./
drwxrwx--- 6 kalavatt  110 May 29 09:31 ../
-rw-rw---- 1 kalavatt 648K May 29 09:36 Schizosaccharomyces_pombe_all_chromosomes.gff3.gz


❯ zcat "Schizosaccharomyces_pombe_all_chromosomes.gff3.gz" \
>    | cut -f 1 \
>    | sort \
>    | uniq
##gff-version 3
SP_I
SP_II
SP_III
SP_II_TG
SP_Mito
SP_MTR
```
</details>
<br />
<br />

<a id="prepare-s-cerevisiae-fasta-gff3-for-concatenation-with-s-pombe"></a>
## Prepare *S. cerevisiae* fasta, gff3 for concatenation with *S. pombe*
<a id="code-4"></a>
### Code
<details>
<summary><i>Code: Prepare S. cerevisiae fasta, gff3 for concatenation with S. pombe</i></summary>

```bash
#!/bin/bash

cd "${HOME}/genomes/" ||
    echo "cd'ing failed; check on this..."

cd "${d_cerevisiae}" ||
    echo "cd'ing failed; check on this..."


#  Process fasta --------------------------------------------------------------
if [[ ! -d "fasta-processed" ]]; then mkdir "fasta-processed"; fi

#  Copy over file to process, then cd into directory
cp \
    "fasta/S288C_reference_sequence_R64-3-1_20210421.fsa.gz" \
    "fasta-processed/"

cd "fasta-processed" ||
    echo "cd'ing failed; check on this..."

#  Check the chromosome names in the fasta
zcat "S288C_reference_sequence_R64-3-1_20210421.fsa.gz" | grep "^>"

#  Simplify the names of chromosomes in the S. cerevisiae fasta
if [[ -f "tmp.fa" ]]; then rm "tmp.fa"; fi
zcat "S288C_reference_sequence_R64-3-1_20210421.fsa.gz" \
    | sed 's/^>ref|NC_001133|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=I\]/>I/g;s/^>ref|NC_001134|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=II\]/>II/g;s/^>ref|NC_001135|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=III\]/>III/g;s/^>ref|NC_001136|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=IV\]/>IV/g;s/^>ref|NC_001137|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=V\]/>V/g;s/^>ref|NC_001138|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=VI\]/>VI/g;s/^>ref|NC_001139|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=VII\]/>VII/g;s/^>ref|NC_001140|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=VIII\]/>VIII/g;s/^>ref|NC_001141|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=IX\]/>IX/g;s/^>ref|NC_001142|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=X\]/>X/g;s/^>ref|NC_001143|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=XI\]/>XI/g;s/^>ref|NC_001144|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=XII\]/>XII/g;s/^>ref|NC_001145|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=XIII\]/>XIII/g;s/^>ref|NC_001146|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=XIV\]/>XIV/g;s/^>ref|NC_001147|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=XV\]/>XV/g;s/^>ref|NC_001148|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=XVI\]/>XVI/g;s/^>ref|NC_001224|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[location=mitochondrion\]\ \[top=circular\]/>Mito/g' \
        > "tmp.fa"

#  Check the new chromosome names in "tmp.fa"
cat "tmp.fa" | grep "^>"

#  Remove the intial infile
rm "S288C_reference_sequence_R64-3-1_20210421.fsa.gz"

#  Rename "tmp.fa" to "S288C_reference_sequence_R64-3-1_20210421.fa"
mv -f "tmp.fa" "S288C_reference_sequence_R64-3-1_20210421.fa"

#  Compress "S288C_reference_sequence_R64-3-1_20210421.fa"
gzip "S288C_reference_sequence_R64-3-1_20210421.fa"

#  Check the directory contents
.,

#  Check the chromosome names
zcat "S288C_reference_sequence_R64-3-1_20210421.fa.gz" | grep "^>"


#  Process gff3 ---------------------------------------------------------------
cd .. && pwd

#  Create a directory for the processed gff3
if [[ ! -d "gff3-processed/" ]]; then mkdir "gff3-processed/"; fi

#  Copy over file to process, then cd into directory
cp "gff3/saccharomyces_cerevisiae_R64-3-1_20210421.gff.gz" "gff3-processed/"

cd "gff3-processed/" ||
    echo "cd'ing failed; check on this..."

#  What do chromosome names look like?
zcat "saccharomyces_cerevisiae_R64-3-1_20210421.gff.gz" \
    | cut -f 1 \
    | sort \
    | uniq
#NOTE There are fasta sequences in this file that need to be excluded

#  Decompress "saccharomyces_cerevisiae_R64-3-1_20210421.gff.gz"
gzip -cd "saccharomyces_cerevisiae_R64-3-1_20210421.gff.gz" \
    > "saccharomyces_cerevisiae_R64-3-1_20210421.gff"

#  Remove everything after line containing ### (fasta chromosome sequences)
if [[ -f "tmp.gff3" ]]; then rm "tmp.gff3"; fi
sed -n '/###/q;p' < saccharomyces_cerevisiae_R64-3-1_20210421.gff > tmp.gff3

#  Check the chromosome names now that fasta sequences are gone
cat "tmp.gff3" \
    | cut -f 1 \
    | sort \
    | uniq

#  Use sed to rename chromosomes; save the results to "tmp.gff3"
cat "tmp.gff3" \
    | sed 's/^chr//g;s/^mt/Mito/g' \
        > "tmp.2.gff3"

#  Check the chromosome names in the updated file
cat "tmp.2.gff3" \
    | cut -f 1 \
    | sort \
    | uniq

#  Check on the file contents
head "tmp.2.gff3"
tail "tmp.2.gff3"

#  Rename "tmp.2.gff3" to "saccharomyces_cerevisiae_R64-3-1_20210421.gff3"
mv "tmp.2.gff3" "saccharomyces_cerevisiae_R64-3-1_20210421.gff3"

#  Compress "saccharomyces_cerevisiae_R64-3-1_20210421.gff3"
gzip "saccharomyces_cerevisiae_R64-3-1_20210421.gff3"

#  Check the directory contents
.,

#  Remove unneeded files
rm \
    "saccharomyces_cerevisiae_R64-3-1_20210421.gff" \
    "saccharomyces_cerevisiae_R64-3-1_20210421.gff.gz" \
    "tmp.gff3"

#  Check the directory contents again
.,

#  Check chromosome names again
zcat "saccharomyces_cerevisiae_R64-3-1_20210421.gff3.gz" \
    | cut -f 1 \
    | sort \
    | uniq
```
</details>
<br />

<a id="printed-3"></a>
### Printed
<details>
<summary><i>Printed: Prepare S. cerevisiae fasta, gff3 for concatenation with S. pombe</i></summary>

```txt
❯ cd "${HOME}/genomes/" ||
>    echo "cd'ing failed; check on this..."


❯ cd "${d_cerevisiae}" ||
>    echo "cd'ing failed; check on this..."
/home/kalavatt/genomes/Saccharomyces_cerevisiae


❯ if [[ ! -d "fasta-processed" ]]; then mkdir "fasta-processed"; fi
mkdir: created directory 'fasta-processed'


❯ cp \
>    "fasta/S288C_reference_sequence_R64-3-1_20210421.fsa.gz" \
>    "fasta-processed/"
'fasta/S288C_reference_sequence_R64-3-1_20210421.fsa.gz' -> 'fasta-processed/S288C_reference_sequence_R64-3-1_20210421.fsa.gz'


❯ cd "fasta-processed" ||
>    echo "cd'ing failed; check on this..."
/home/kalavatt/genomes/Saccharomyces_cerevisiae/fasta-processed


❯ zcat "S288C_reference_sequence_R64-3-1_20210421.fsa.gz" | grep "^>"
>ref|NC_001133| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=I]
>ref|NC_001134| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=II]
>ref|NC_001135| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=III]
>ref|NC_001136| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=IV]
>ref|NC_001137| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=V]
>ref|NC_001138| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=VI]
>ref|NC_001139| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=VII]
>ref|NC_001140| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=VIII]
>ref|NC_001141| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=IX]
>ref|NC_001142| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=X]
>ref|NC_001143| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XI]
>ref|NC_001144| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XII]
>ref|NC_001145| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XIII]
>ref|NC_001146| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XIV]
>ref|NC_001147| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XV]
>ref|NC_001148| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XVI]
>ref|NC_001224| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [location=mitochondrion] [top=circular]


❯ zcat "S288C_reference_sequence_R64-3-1_20210421.fsa.gz" \
>    | sed 's/^>ref|NC_001133|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=I\]/>I/g;s/^>ref|NC_001134|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=II\]/>II/g;s/^>ref|NC_001135|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=III\]/>III/g;s/^>ref|NC_001136|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=IV\]/>IV/g;s/^>ref|NC_001137|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=V\]/>V/g;s/^>ref|NC_001138|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=VI\]/>VI/g;s/^>ref|NC_001139|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=VII\]/>VII/g;s/^>ref|NC_001140|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=VIII\]/>VIII/g;s/^>ref|NC_001141|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=IX\]/>IX/g;s/^>ref|NC_001142|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=X\]/>X/g;s/^>ref|NC_001143|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=XI\]/>XI/g;s/^>ref|NC_001144|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=XII\]/>XII/g;s/^>ref|NC_001145|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=XIII\]/>XIII/g;s/^>ref|NC_001146|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=XIV\]/>XIV/g;s/^>ref|NC_001147|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=XV\]/>XV/g;s/^>ref|NC_001148|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=XVI\]/>XVI/g;s/^>ref|NC_001224|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[location=mitochondrion\]\ \[top=circular\]/>Mito/g' \
>        > "tmp.fa"


❯ cat "tmp.fa" | grep "^>"
>I
>II
>III
>IV
>V
>VI
>VII
>VIII
>IX
>X
>XI
>XII
>XIII
>XIV
>XV
>XVI
>Mito


❯ rm "S288C_reference_sequence_R64-3-1_20210421.fsa.gz"


❯ mv -f "tmp.fa" "S288C_reference_sequence_R64-3-1_20210421.fa"
renamed 'tmp.fa' -> 'S288C_reference_sequence_R64-3-1_20210421.fa'


❯ gzip "S288C_reference_sequence_R64-3-1_20210421.fa"


❯ .,
total 4.5M
drwxrwx--- 2 kalavatt   65 May 29 10:07 ./
drwxrwx--- 5 kalavatt  139 May 29 09:53 ../
-rw-rw---- 1 kalavatt 3.7M May 29 09:58 S288C_reference_sequence_R64-3-1_20210421.fa.gz


❯ zcat "S288C_reference_sequence_R64-3-1_20210421.fa.gz" | grep "^>"
>I
>II
>III
>IV
>V
>VI
>VII
>VIII
>IX
>X
>XI
>XII
>XIII
>XIV
>XV
>XVI
>Mito


❯ cd .. && pwd
/home/kalavatt/genomes/Saccharomyces_cerevisiae


❯ if [[ ! -d "gff3-processed/" ]]; then mkdir "gff3-processed/"; fi
mkdir: created directory 'gff3-processed/'


❯ cp "gff3/saccharomyces_cerevisiae_R64-3-1_20210421.gff.gz" "gff3-processed/"
'gff3/saccharomyces_cerevisiae_R64-3-1_20210421.gff.gz' -> 'gff3-processed/saccharomyces_cerevisiae_R64-3-1_20210421.gff.gz'


❯ cd "gff3-processed/" ||
>    echo "cd'ing failed; check on this..."
/home/kalavatt/genomes/Saccharomyces_cerevisiae/gff3-processed


❯ gzip -cd "saccharomyces_cerevisiae_R64-3-1_20210421.gff.gz" \
>    > "saccharomyces_cerevisiae_R64-3-1_20210421.gff"


❯ if [[ -f "tmp.gff3" ]]; then rm "tmp.gff3"; fi


❯ sed -n '/###/q;p' < saccharomyces_cerevisiae_R64-3-1_20210421.gff > tmp.gff3


❯ cat "tmp.gff3" \
>    | cut -f 1 \
>    | sort \
>    | uniq
#
#!assembly R64-3-1
chrI
chrII
chrIII
chrIV
chrIX
chrmt
chrV
chrVI
chrVII
chrVIII
chrX
chrXI
chrXII
chrXIII
chrXIV
chrXV
chrXVI
# Created by Saccharomyces Genome Database (http://www.yeastgenome.org/)
#!data-source SGD
#!date-produced 2021-04-27 10:49:32
# Features from the 16 nuclear chromosomes labeled chrI to chrXVI,
##gff-version 3
# https://downloads.yeastgenome.org/latest/saccharomyces_cerevisiae.gff.gz
# Please send comments and suggestions to sgd-helpdesk@lists.stanford.edu
# plus the mitochondrial genome labeled chrmt.
#!refseq-version GCF_000146045.2
# Saccharomyces cerevisiae S288C genome (version=R64-3-1)
# SGD is funded as a National Human Genome Research Institute Biomedical Informatics Resource from
# the U. S. National Institutes of Health to Stanford University.
# Weekly updates of this file are available for download from:


❯ cat "tmp.gff3" \
>    | sed 's/^chr//g;s/^mt/Mito/g' \
>        > "tmp.2.gff3"


❯ cat "tmp.2.gff3" \
>    | cut -f 1 \
>    | sort \
>    | uniq
#
#!assembly R64-3-1
# Created by Saccharomyces Genome Database (http://www.yeastgenome.org/)
#!data-source SGD
#!date-produced 2021-04-27 10:49:32
# Features from the 16 nuclear chromosomes labeled chrI to chrXVI,
##gff-version 3
# https://downloads.yeastgenome.org/latest/saccharomyces_cerevisiae.gff.gz
I
II
III
IV
IX
Mito
# Please send comments and suggestions to sgd-helpdesk@lists.stanford.edu
# plus the mitochondrial genome labeled chrmt.
#!refseq-version GCF_000146045.2
# Saccharomyces cerevisiae S288C genome (version=R64-3-1)
# SGD is funded as a National Human Genome Research Institute Biomedical Informatics Resource from
# the U. S. National Institutes of Health to Stanford University.
V
VI
VII
VIII
# Weekly updates of this file are available for download from:
X
XI
XII
XIII
XIV
XV
XVI


❯ head "tmp.2.gff3"
##gff-version 3
#!date-produced 2021-04-27 10:49:32
#!data-source SGD
#!assembly R64-3-1
#!refseq-version GCF_000146045.2
#
# Saccharomyces cerevisiae S288C genome (version=R64-3-1)
#
# Features from the 16 nuclear chromosomes labeled chrI to chrXVI,
# plus the mitochondrial genome labeled chrmt.


❯ tail "tmp.2.gff3"
Mito    SGD origin_of_replication   82329   82600   .   +   .   ID=ORI5;Name=ORI5;gene=ORI5;Alias=ORI5;Note=Mitochondrial%20origin%20of%20replication;display=Mitochondrial%20origin%20of%20replication;dbxref=SGD:S000029671;curie=SGD:S000029671
Mito    SGD tRNA_gene   85035   85112   .   +   .   ID=YNCQ0026W;Name=YNCQ0026W;Alias=tM%28CAU%29Q2,tRNA-fMet;Ontology_term=GO:0005739,GO:0030533,GO:0070125,SO:0000704;Note=Mitochondrial%20formylated%20methionine%20tRNA%20%28tRNA-fMet%29%3B%20predicted%20by%20tRNAscan-SE%20analysis;display=Mitochondrial%20formylated%20methionine%20tRNA%20%28tRNA-fMet%29;dbxref=SGD:S000007326;curie=SGD:S000007326
Mito    SGD noncoding_exon  85035   85112   .   +   .   Parent=YNCQ0026W_tRNA;Name=YNCQ0026W_noncoding_exon
Mito    SGD tRNA    85035   85112   .   +   .   ID=YNCQ0026W_tRNA;Name=YNCQ0026W_tRNA;Parent=YNCQ0026W
Mito    SGD ncRNA_gene  85295   85777   .   +   .   ID=YNCQ0027W;Name=YNCQ0027W;gene=RPM1;Alias=RPM1,Q0285;Ontology_term=GO:0001682,GO:0004526,GO:0005739,GO:0008033,GO:0030678,SO:0000704;Note=RNA%20component%20of%20mitochondrial%20RNase%20P%3B%20mitochondrial%20RNase%20P%20also%20contains%20the%20protein%20subunit%20Rpm2p%3B%20RNase%20P%20removes%205'%20extensions%20from%20mitochondrial%20tRNA%20precursors%3B%20RPM1%20is%20conserved%20in%20bacteria%2C%20fungi%2C%20and%20protozoa;display=RNA%20component%20of%20mitochondrial%20RNase%20P;dbxref=SGD:S000029023;curie=SGD:S000029023
Mito    SGD noncoding_exon  85295   85777   .   +   .   Parent=YNCQ0027W_ncRNA;Name=YNCQ0027W_noncoding_exon
Mito    SGD ncRNA   85295   85777   .   +   .   ID=YNCQ0027W_ncRNA;Name=YNCQ0027W_ncRNA;Parent=YNCQ0027W
Mito    SGD gene    85554   85709   .   +   .   ID=Q0297;Name=Q0297;Alias=ORF12;Ontology_term=GO:0003674,GO:0005575,GO:0008150,SO:0000704;Note=Dubious%20open%20reading%20frame%3B%20unlikely%20to%20encode%20a%20functional%20protein%2C%20based%20on%20available%20experimental%20and%20comparative%20sequence%20data%3B%20partially%20overlaps%20the%20verified%20gene%20RPM1;display=Dubious%20open%20reading%20frame;dbxref=SGD:S000007284;orf_classification=Dubious;curie=SGD:S000007284
Mito    SGD CDS 85554   85709   .   +   0   Parent=Q0297_mRNA;Name=Q0297_CDS;orf_classification=Dubious;protein_id=UniProtKB:Q9ZZV8
Mito    SGD mRNA    85554   85709   .   +   .   ID=Q0297_mRNA;Name=Q0297_mRNA;Parent=Q0297


❯ mv "tmp.2.gff3" "saccharomyces_cerevisiae_R64-3-1_20210421.gff3"
renamed 'tmp.2.gff3' -> 'saccharomyces_cerevisiae_R64-3-1_20210421.gff3'


❯ gzip "saccharomyces_cerevisiae_R64-3-1_20210421.gff3"


❯ mv "tmp.2.gff3" "saccharomyces_cerevisiae_R64-3-1_20210421.gff3"
renamed 'tmp.2.gff3' -> 'saccharomyces_cerevisiae_R64-3-1_20210421.gff3'


❯ .,
total 34M
drwxrwx--- 2 kalavatt  222 May 29 10:41 ./
drwxrwx--- 6 kalavatt  171 May 29 10:15 ../
-rw-rw---- 1 kalavatt  20M May 29 10:21 saccharomyces_cerevisiae_R64-3-1_20210421.gff
-rw-rw---- 1 kalavatt 1.6M May 29 10:37 saccharomyces_cerevisiae_R64-3-1_20210421.gff3.gz
-rw-r----- 1 kalavatt 5.1M May 29 10:16 saccharomyces_cerevisiae_R64-3-1_20210421.gff.gz
-rw-rw---- 1 kalavatt 7.4M May 29 10:34 tmp.gff3
 

❯ rm \
>    "saccharomyces_cerevisiae_R64-3-1_20210421.gff" \
>    "saccharomyces_cerevisiae_R64-3-1_20210421.gff.gz" \
>    "tmp.gff3"


❯ .,
total 1.9M
drwxrwx--- 2 kalavatt   67 May 29 10:42 ./
drwxrwx--- 6 kalavatt  171 May 29 10:15 ../
-rw-rw---- 1 kalavatt 1.6M May 29 10:37 saccharomyces_cerevisiae_R64-3-1_20210421.gff3.gz


❯ zcat "saccharomyces_cerevisiae_R64-3-1_20210421.gff3.gz" \
>    | cut -f 1 \
>    | sort \
>    | uniq
#
#!assembly R64-3-1
# Created by Saccharomyces Genome Database (http://www.yeastgenome.org/)
#!data-source SGD
#!date-produced 2021-04-27 10:49:32
# Features from the 16 nuclear chromosomes labeled chrI to chrXVI,
##gff-version 3
# https://downloads.yeastgenome.org/latest/saccharomyces_cerevisiae.gff.gz
I
II
III
IV
IX
Mito
# Please send comments and suggestions to sgd-helpdesk@lists.stanford.edu
# plus the mitochondrial genome labeled chrmt.
#!refseq-version GCF_000146045.2
# Saccharomyces cerevisiae S288C genome (version=R64-3-1)
# SGD is funded as a National Human Genome Research Institute Biomedical Informatics Resource from
# the U. S. National Institutes of Health to Stanford University.
V
VI
VII
VIII
# Weekly updates of this file are available for download from:
X
XI
XII
XIII
XIV
XV
XVI
```
</details>
<br />
<br />

<a id="concatenate-processed-fastas-and-gff3s-in-new-directory-combined_sc_sp"></a>
## Concatenate processed fastas and gff3s in new directory `combined_SC_SP/`
<a id="code-5"></a>
### Code
<details>
<summary><i>Code: Concatenate processed fastas and gff3s in new directory combined_SC_SP/</i></summary>

```bash
#!/bin/bash

cd "${HOME}/tsukiyamalab/kalavatt/genomes/" ||
    echo "cd'ing failed; check on this..."

if [[ ! -d "combined_SC_SP" ]]; then mkdir -p combined_SC_SP/{fasta,gff3}; fi

cp \
    "Saccharomyces_cerevisiae/gff3-processed/saccharomyces_cerevisiae_R64-3-1_20210421.gff3.gz" \
    "combined_SC_SP/gff3/"
cp \
    "Schizosaccharomyces_pombe/gff3-processed/Schizosaccharomyces_pombe_all_chromosomes.gff3.gz" \
    "combined_SC_SP/gff3/"

cp \
    "Saccharomyces_cerevisiae/fasta-processed/S288C_reference_sequence_R64-3-1_20210421.fa.gz" \
    "combined_SC_SP/fasta/"
cp \
    "Schizosaccharomyces_pombe/fasta-processed/Schizosaccharomyces_pombe_all_chromosomes.fa.gz" \
    "combined_SC_SP/fasta/"

cd "combined_SC_SP/" ||
    echo "cd'ing failed; check on this..."

cd "gff3/" ||
    echo "cd'ing failed; check on this..."

cat \
    "saccharomyces_cerevisiae_R64-3-1_20210421.gff3.gz" \
    "Schizosaccharomyces_pombe_all_chromosomes.gff3.gz" \
        > "combined_SC_SP.gff3.gz"

.,

zcat "combined_SC_SP.gff3.gz" \
    | cut -f 1 \
    | sort \
    | uniq

cd "../fasta" ||
    echo "cd'ing failed; check on this..."

cat \
    "S288C_reference_sequence_R64-3-1_20210421.fa.gz" \
    "Schizosaccharomyces_pombe_all_chromosomes.fa.gz" \
        > "combined_SC_SP.fa.gz"

.,

zcat "combined_SC_SP.fa.gz" | grep "^>"
```
</details>
<br />

<details>
<summary><i>Printed: Concatenate processed fastas and gff3s in new directory combined_SC_SP/</i></summary>

```txt
❯ cd "${HOME}/tsukiyamalab/kalavatt/genomes/combined_SC_SP" ||
>    echo "cd'ing failed; check on this..."


❯ if [[ ! -d "combined_SC_SP" ]]; then mkdir -p combined_SC_SP/{fasta,gff3}; fi
mkdir: created directory 'combined_SC_SP'
mkdir: created directory 'combined_SC_SP/fasta'
mkdir: created directory 'combined_SC_SP/gff3'


❯ cp \
>    "Saccharomyces_cerevisiae/gff3-processed/saccharomyces_cerevisiae_R64-3-1_20210421.gff3.gz" \
>    "combined_SC_SP/gff3/"
'Saccharomyces_cerevisiae/gff3-processed/saccharomyces_cerevisiae_R64-3-1_20210421.gff3.gz' -> 'combined_SC_SP/gff3/saccharomyces_cerevisiae_R64-3-1_20210421.gff3.gz'


❯ cp \
>    "Schizosaccharomyces_pombe/gff3-processed/Schizosaccharomyces_pombe_all_chromosomes.gff3.gz" \
>    "combined_SC_SP/gff3/"
'Schizosaccharomyces_pombe/gff3-processed/Schizosaccharomyces_pombe_all_chromosomes.gff3.gz' -> 'combined_SC_SP/gff3/Schizosaccharomyces_pombe_all_chromosomes.gff3.gz'


❯ cp \
>    "Saccharomyces_cerevisiae/fasta-processed/S288C_reference_sequence_R64-3-1_20210421.fa.gz" \
>    "combined_SC_SP/fasta/"
'Saccharomyces_cerevisiae/fasta-processed/S288C_reference_sequence_R64-3-1_20210421.fa.gz' -> 'combined_SC_SP/fasta/S288C_reference_sequence_R64-3-1_20210421.fa.gz'


❯ cp \
>    "Schizosaccharomyces_pombe/fasta-processed/Schizosaccharomyces_pombe_all_chromosomes.fa.gz" \
>    "combined_SC_SP/fasta/"
'Schizosaccharomyces_pombe/fasta-processed/Schizosaccharomyces_pombe_all_chromosomes.fa.gz' -> 'combined_SC_SP/fasta/Schizosaccharomyces_pombe_all_chromosomes.fa.gz'


❯ cd "combined_SC_SP/" ||
>    echo "cd'ing failed; check on this..."
/home/kalavatt/tsukiyamalab/kalavatt/genomes/combined_SC_SP


❯ cd "gff3/" ||
>    echo "cd'ing failed; check on this..."
/home/kalavatt/tsukiyamalab/kalavatt/genomes/combined_SC_SP/gff3


❯ cat \
>    "saccharomyces_cerevisiae_R64-3-1_20210421.gff3.gz" \
>    "Schizosaccharomyces_pombe_all_chromosomes.gff3.gz" \
>        > "combined_SC_SP.gff3.gz"


❯ .,
total 5.9M
drwxrwx--- 2 kalavatt  174 May 29 10:56 ./
drwxrwx--- 4 kalavatt   45 May 29 10:52 ../
-rw-rw---- 1 kalavatt 2.2M May 29 10:56 combined_SC_SP.gff3.gz
-rw-rw---- 1 kalavatt 1.6M May 29 10:52 saccharomyces_cerevisiae_R64-3-1_20210421.gff3.gz
-rw-rw---- 1 kalavatt 648K May 29 10:52 Schizosaccharomyces_pombe_all_chromosomes.gff3.gz


❯ zcat "combined_SC_SP.gff3.gz" \
>    | cut -f 1 \
>    | sort \
>    | uniq
#
#!assembly R64-3-1
# Created by Saccharomyces Genome Database (http://www.yeastgenome.org/)
#!data-source SGD
#!date-produced 2021-04-27 10:49:32
# Features from the 16 nuclear chromosomes labeled chrI to chrXVI,
##gff-version 3
# https://downloads.yeastgenome.org/latest/saccharomyces_cerevisiae.gff.gz
I
II
III
IV
IX
Mito
# Please send comments and suggestions to sgd-helpdesk@lists.stanford.edu
# plus the mitochondrial genome labeled chrmt.
#!refseq-version GCF_000146045.2
# Saccharomyces cerevisiae S288C genome (version=R64-3-1)
# SGD is funded as a National Human Genome Research Institute Biomedical Informatics Resource from
SP_I
SP_II
SP_III
SP_II_TG
SP_Mito
SP_MTR
# the U. S. National Institutes of Health to Stanford University.
V
VI
VII
VIII
# Weekly updates of this file are available for download from:
X
XI
XII
XIII
XIV
XV
XVI


❯ cd "../fasta" ||
>    echo "cd'ing failed; check on this..."


❯ cat \
>    "S288C_reference_sequence_R64-3-1_20210421.fa.gz" \
>    "Schizosaccharomyces_pombe_all_chromosomes.fa.gz" \
>        > "combined_SC_SP.fa.gz"


❯ .,
total 18M
drwxrwx--- 2 kalavatt  168 May 29 11:01 ./
drwxrwx--- 4 kalavatt   45 May 29 10:52 ../
-rw-rw---- 1 kalavatt 7.4M May 29 11:01 combined_SC_SP.fa.gz
-rw-rw---- 1 kalavatt 3.7M May 29 10:54 S288C_reference_sequence_R64-3-1_20210421.fa.gz
-rw-rw---- 1 kalavatt 3.8M May 29 10:54 Schizosaccharomyces_pombe_all_chromosomes.fa.gz


❯ zcat "combined_SC_SP.fa.gz" | grep "^>"
>I
>II
>III
>IV
>V
>VI
>VII
>VIII
>IX
>X
>XI
>XII
>XIII
>XIV
>XV
>XVI
>Mito
>SP_II_TG
>SP_I
>SP_II
>SP_III
>SP_MTR
>SP_Mito
```
</details>
<br />
<br />

<a id="create-bowtie2-indices-for-combined_sc_spfagz"></a>
## Create `bowtie2` indices for "`combined_SC_SP.fa.gz`"
<a id="code-6"></a>
### Code
<details>
<summary><i>Code: Create bowtie2 indices for "combined_SC_SP.fa.gz"</i></summary>

```bash
#!/bin/bash

cd "${HOME}/genomes/combined_SC_SP" ||
    echo "cd'ing failed; check on this..."

if [[ ! -d "bowtie2/" ]]; then mkdir "bowtie2/"; fi

#  Index the fasta file
cd "fasta/"

gzip -cd "combined_SC_SP.fa.gz" > "combined_SC_SP.fa"

ml SAMtools/1.16.1-GCC-11.2.0 Bowtie2/2.4.4-GCC-11.2.0

cat "combined_SC_SP.fa" | grep "^>"

samtools faidx "combined_SC_SP.fa"

#  Create a "chrom-info" file
cut -f 1,2 "combined_SC_SP.fa.fai" > "combined_SC_SP.chrom-info.tsv"

#  Build the indices
cd .. && pwd

bowtie2-build fasta/combined_SC_SP.fa bowtie2/combined_SC_SP \
    1> >(tee -a bowtie2/combined_SC_SP.stdout.txt) \
    2> >(tee -a bowtie2/combined_SC_SP.stderr.txt)
```
</details>
<br />

<a id="printed-4"></a>
### Printed
<details>
<summary><i>Printed: Create bowtie2 indices for "combined_SC_SP.fa.gz"</i></summary>

```txt
❯ cd "${HOME}/genomes/combined_SC_SP" ||
>    echo "cd'ing failed; check on this..."


❯ if [[ ! -d "bowtie2/" ]]; then mkdir "bowtie2/"; fi
mkdir: created directory 'bowtie2/'


❯ #  Index the fasta file


❯ cd "fasta/"
/home/kalavatt/genomes/combined_SC_SP/fasta


❯ gzip -cd "combined_SC_SP.fa.gz" > "combined_SC_SP.fa"


❯ ml SAMtools/1.16.1-GCC-11.2.0 Bowtie2/2.4.4-GCC-11.2.0


❯ cat "combined_SC_SP.fa" | grep "^>"
>I
>II
>III
>IV
>V
>VI
>VII
>VIII
>IX
>X
>XI
>XII
>XIII
>XIV
>XV
>XVI
>Mito
>SP_II_TG
>SP_I
>SP_II
>SP_III
>SP_MTR
>SP_Mito


❯ samtools faidx "combined_SC_SP.fa"


❯ #  Create a "chrom-info" file


❯ cut -f 1,2 "combined_SC_SP.fa.fai" > "combined_SC_SP.chrom-info.tsv"


❯ #  Build the indices


❯ cd .. && pwd
/home/kalavatt/genomes/combined_SC_SP


❯ bowtie2-build fasta/combined_SC_SP.fa bowtie2/combined_SC_SP \
>    1> >(tee -a bowtie2/combined_SC_SP.stdout.txt) \
>    2> >(tee -a bowtie2/combined_SC_SP.stderr.txt)
Settings:
  Output files: "bowtie2/combined_SC_SP.*.bt2"
  Line rate: 6 (line is 64 bytes)
  Lines per side: 1 (side is 64 bytes)
  Offset rate: 4 (one in 16)
  FTable chars: 10
  Strings: unpacked
  Max bucket size: default
  Max bucket size, sqrt multiplier: default
  Max bucket size, len divisor: 4
  Difference-cover sample period: 1024
  Endianness: little
  Actual local endianness: little
  Sanity checking: disabled
  Assertions: disabled
  Random seed: 0
  Sizeofs: void*:8, int:4, long:8, size_t:8
Input files DNA, FASTA:
  fasta/combined_SC_SP.fa
Reading reference sizes
Building a SMALL index
  Time reading reference sizes: 00:00:00
Calculating joined length
Writing header
Reserving space for joined string
Joining reference sequences
  Time to join reference sequences: 00:00:00
bmax according to bmaxDivN setting: 6197021
Using parameters --bmax 4647766 --dcv 1024
  Doing ahead-of-time memory usage test
  Passed!  Constructing with these parameters: --bmax 4647766 --dcv 1024
Constructing suffix-array element generator
Building DifferenceCoverSample
  Building sPrime
  Building sPrimeOrder
  V-Sorting samples
  V-Sorting samples time: 00:00:00
  Allocating rank array
  Ranking v-sort output
  Ranking v-sort output time: 00:00:00
  Invoking Larsson-Sadakane on ranks
  Invoking Larsson-Sadakane on ranks time: 00:00:00
  Sanity-checking and returning
Building samples
Reserving space for 12 sample suffixes
Generating random suffixes
QSorting 12 sample offsets, eliminating duplicates
QSorting sample offsets, eliminating duplicates time: 00:00:00
Multikey QSorting 12 samples
  (Using difference cover)
  Multikey QSorting samples time: 00:00:00
Calculating bucket sizes
Splitting and merging
  Splitting and merging time: 00:00:00
Avg bucket size: 2.75423e+06 (target: 4647765)
Converting suffix-array elements to index image
Allocating ftab, absorbFtab
Entering Ebwt loop
Getting block 1 of 9
  Reserving size (4647766) for bucket 1
  Calculating Z arrays for bucket 1
  Entering block accumulator loop for bucket 1:
  bucket 1: 10%
  bucket 1: 20%
  bucket 1: 30%
  bucket 1: 40%
  bucket 1: 50%
  bucket 1: 60%
  bucket 1: 70%
  bucket 1: 80%
  bucket 1: 90%
  bucket 1: 100%
  Sorting block of length 3260001 for bucket 1
  (Using difference cover)
  Sorting block time: 00:00:01
Returning block of 3260002 for bucket 1
Getting block 2 of 9
  Reserving size (4647766) for bucket 2
  Calculating Z arrays for bucket 2
  Entering block accumulator loop for bucket 2:
  bucket 2: 10%
  bucket 2: 20%
  bucket 2: 30%
  bucket 2: 40%
  bucket 2: 50%
  bucket 2: 60%
  bucket 2: 70%
  bucket 2: 80%
  bucket 2: 90%
  bucket 2: 100%
  Sorting block of length 2867770 for bucket 2
  (Using difference cover)
  Sorting block time: 00:00:00
Returning block of 2867771 for bucket 2
Getting block 3 of 9
  Reserving size (4647766) for bucket 3
  Calculating Z arrays for bucket 3
  Entering block accumulator loop for bucket 3:
  bucket 3: 10%
  bucket 3: 20%
  bucket 3: 30%
  bucket 3: 40%
  bucket 3: 50%
  bucket 3: 60%
  bucket 3: 70%
  bucket 3: 80%
  bucket 3: 90%
  bucket 3: 100%
  Sorting block of length 2560185 for bucket 3
  (Using difference cover)
  Sorting block time: 00:00:00
Returning block of 2560186 for bucket 3
Getting block 4 of 9
  Reserving size (4647766) for bucket 4
  Calculating Z arrays for bucket 4
  Entering block accumulator loop for bucket 4:
  bucket 4: 10%
  bucket 4: 20%
  bucket 4: 30%
  bucket 4: 40%
  bucket 4: 50%
  bucket 4: 60%
  bucket 4: 70%
  bucket 4: 80%
  bucket 4: 90%
  bucket 4: 100%
  Sorting block of length 2972413 for bucket 4
  (Using difference cover)
  Sorting block time: 00:00:00
Returning block of 2972414 for bucket 4
Getting block 5 of 9
  Reserving size (4647766) for bucket 5
  Calculating Z arrays for bucket 5
  Entering block accumulator loop for bucket 5:
  bucket 5: 10%
  bucket 5: 20%
  bucket 5: 30%
  bucket 5: 40%
  bucket 5: 50%
  bucket 5: 60%
  bucket 5: 70%
  bucket 5: 80%
  bucket 5: 90%
  bucket 5: 100%
  Sorting block of length 2671929 for bucket 5
  (Using difference cover)
  Sorting block time: 00:00:01
Returning block of 2671930 for bucket 5
Getting block 6 of 9
  Reserving size (4647766) for bucket 6
  Calculating Z arrays for bucket 6
  Entering block accumulator loop for bucket 6:
  bucket 6: 10%
  bucket 6: 20%
  bucket 6: 30%
  bucket 6: 40%
  bucket 6: 50%
  bucket 6: 60%
  bucket 6: 70%
  bucket 6: 80%
  bucket 6: 90%
  bucket 6: 100%
  Sorting block of length 4280047 for bucket 6
  (Using difference cover)
  Sorting block time: 00:00:01
Returning block of 4280048 for bucket 6
Getting block 7 of 9
  Reserving size (4647766) for bucket 7
  Calculating Z arrays for bucket 7
  Entering block accumulator loop for bucket 7:
  bucket 7: 10%
  bucket 7: 20%
  bucket 7: 30%
  bucket 7: 40%
  bucket 7: 50%
  bucket 7: 60%
  bucket 7: 70%
  bucket 7: 80%
  bucket 7: 90%
  bucket 7: 100%
  Sorting block of length 736213 for bucket 7
  (Using difference cover)
  Sorting block time: 00:00:00
Returning block of 736214 for bucket 7
Getting block 8 of 9
  Reserving size (4647766) for bucket 8
  Calculating Z arrays for bucket 8
  Entering block accumulator loop for bucket 8:
  bucket 8: 10%
  bucket 8: 20%
  bucket 8: 30%
  bucket 8: 40%
  bucket 8: 50%
  bucket 8: 60%
  bucket 8: 70%
  bucket 8: 80%
  bucket 8: 90%
  bucket 8: 100%
  Sorting block of length 4277836 for bucket 8
  (Using difference cover)
  Sorting block time: 00:00:01
Returning block of 4277837 for bucket 8
Getting block 9 of 9
  Reserving size (4647766) for bucket 9
  Calculating Z arrays for bucket 9
  Entering block accumulator loop for bucket 9:
  bucket 9: 10%
  bucket 9: 20%
  bucket 9: 30%
  bucket 9: 40%
  bucket 9: 50%
  bucket 9: 60%
  bucket 9: 70%
  bucket 9: 80%
  bucket 9: 90%
  bucket 9: 100%
  Sorting block of length 1161682 for bucket 9
  (Using difference cover)
  Sorting block time: 00:00:00
Returning block of 1161683 for bucket 9
Exited Ebwt loop
fchr[A]: 0
fchr[C]: 7802696
fchr[G]: 12398259
fchr[T]: 16992544
fchr[$]: 24788084
Exiting Ebwt::buildToDisk()
Returning from initFromVector
Wrote 12457710 bytes to primary EBWT file: bowtie2/combined_SC_SP.1.bt2
Wrote 6197028 bytes to secondary EBWT file: bowtie2/combined_SC_SP.2.bt2
Re-opening _in1 and _in2 as input streams
Returning from Ebwt constructor
Headers:
    len: 24788084
    bwtLen: 24788085
    sz: 6197021
    bwtSz: 6197022
    lineRate: 6
    offRate: 4
    offMask: 0xfffffff0
    ftabChars: 10
    eftabLen: 20
    eftabSz: 80
    ftabLen: 1048577
    ftabSz: 4194308
    offsLen: 1549256
    offsSz: 6197024
    lineSz: 64
    sideSz: 64
    sideBwtSz: 48
    sideBwtLen: 192
    numSides: 129105
    numLines: 129105
    ebwtTotLen: 8262720
    ebwtTotSz: 8262720
    color: 0
    reverse: 0
Total time for call to driver() for forward index: 00:00:11
Reading reference sizes
  Time reading reference sizes: 00:00:00
Calculating joined length
Writing header
Reserving space for joined string
Joining reference sequences
  Time to join reference sequences: 00:00:01
  Time to reverse reference sequence: 00:00:00
bmax according to bmaxDivN setting: 6197021
Using parameters --bmax 4647766 --dcv 1024
  Doing ahead-of-time memory usage test
  Passed!  Constructing with these parameters: --bmax 4647766 --dcv 1024
Constructing suffix-array element generator
Building DifferenceCoverSample
  Building sPrime
  Building sPrimeOrder
  V-Sorting samples
  V-Sorting samples time: 00:00:00
  Allocating rank array
  Ranking v-sort output
  Ranking v-sort output time: 00:00:00
  Invoking Larsson-Sadakane on ranks
  Invoking Larsson-Sadakane on ranks time: 00:00:00
  Sanity-checking and returning
Building samples
Reserving space for 12 sample suffixes
Generating random suffixes
QSorting 12 sample offsets, eliminating duplicates
QSorting sample offsets, eliminating duplicates time: 00:00:00
Multikey QSorting 12 samples
  (Using difference cover)
  Multikey QSorting samples time: 00:00:00
Calculating bucket sizes
Splitting and merging
  Splitting and merging time: 00:00:00
Split 1, merged 7; iterating...
Splitting and merging
  Splitting and merging time: 00:00:00
Avg bucket size: 3.54115e+06 (target: 4647765)
Converting suffix-array elements to index image
Allocating ftab, absorbFtab
Entering Ebwt loop
Getting block 1 of 7
  Reserving size (4647766) for bucket 1
  Calculating Z arrays for bucket 1
  Entering block accumulator loop for bucket 1:
  bucket 1: 10%
  bucket 1: 20%
  bucket 1: 30%
  bucket 1: 40%
  bucket 1: 50%
  bucket 1: 60%
  bucket 1: 70%
  bucket 1: 80%
  bucket 1: 90%
  bucket 1: 100%
  Sorting block of length 4364643 for bucket 1
  (Using difference cover)
  Sorting block time: 00:00:00
Returning block of 4364644 for bucket 1
Getting block 2 of 7
  Reserving size (4647766) for bucket 2
  Calculating Z arrays for bucket 2
  Entering block accumulator loop for bucket 2:
  bucket 2: 10%
  bucket 2: 20%
  bucket 2: 30%
  bucket 2: 40%
  bucket 2: 50%
  bucket 2: 60%
  bucket 2: 70%
  bucket 2: 80%
  bucket 2: 90%
  bucket 2: 100%
  Sorting block of length 3463592 for bucket 2
  (Using difference cover)
  Sorting block time: 00:00:01
Returning block of 3463593 for bucket 2
Getting block 3 of 7
  Reserving size (4647766) for bucket 3
  Calculating Z arrays for bucket 3
  Entering block accumulator loop for bucket 3:
  bucket 3: 10%
  bucket 3: 20%
  bucket 3: 30%
  bucket 3: 40%
  bucket 3: 50%
  bucket 3: 60%
  bucket 3: 70%
  bucket 3: 80%
  bucket 3: 90%
  bucket 3: 100%
  Sorting block of length 3603783 for bucket 3
  (Using difference cover)
  Sorting block time: 00:00:00
Returning block of 3603784 for bucket 3
Getting block 4 of 7
  Reserving size (4647766) for bucket 4
  Calculating Z arrays for bucket 4
  Entering block accumulator loop for bucket 4:
  bucket 4: 10%
  bucket 4: 20%
  bucket 4: 30%
  bucket 4: 40%
  bucket 4: 50%
  bucket 4: 60%
  bucket 4: 70%
  bucket 4: 80%
  bucket 4: 90%
  bucket 4: 100%
  Sorting block of length 4546286 for bucket 4
  (Using difference cover)
  Sorting block time: 00:00:01
Returning block of 4546287 for bucket 4
Getting block 5 of 7
  Reserving size (4647766) for bucket 5
  Calculating Z arrays for bucket 5
  Entering block accumulator loop for bucket 5:
  bucket 5: 10%
  bucket 5: 20%
  bucket 5: 30%
  bucket 5: 40%
  bucket 5: 50%
  bucket 5: 60%
  bucket 5: 70%
  bucket 5: 80%
  bucket 5: 90%
  bucket 5: 100%
  Sorting block of length 1146309 for bucket 5
  (Using difference cover)
  Sorting block time: 00:00:00
Returning block of 1146310 for bucket 5
Getting block 6 of 7
  Reserving size (4647766) for bucket 6
  Calculating Z arrays for bucket 6
  Entering block accumulator loop for bucket 6:
  bucket 6: 10%
  bucket 6: 20%
  bucket 6: 30%
  bucket 6: 40%
  bucket 6: 50%
  bucket 6: 60%
  bucket 6: 70%
  bucket 6: 80%
  bucket 6: 90%
  bucket 6: 100%
  Sorting block of length 3789355 for bucket 6
  (Using difference cover)
  Sorting block time: 00:00:01
Returning block of 3789356 for bucket 6
Getting block 7 of 7
  Reserving size (4647766) for bucket 7
  Calculating Z arrays for bucket 7
  Entering block accumulator loop for bucket 7:
  bucket 7: 10%
  bucket 7: 20%
  bucket 7: 30%
  bucket 7: 40%
  bucket 7: 50%
  bucket 7: 60%
  bucket 7: 70%
  bucket 7: 80%
  bucket 7: 90%
  bucket 7: 100%
  Sorting block of length 3874110 for bucket 7
  (Using difference cover)
  Sorting block time: 00:00:01
Returning block of 3874111 for bucket 7
Exited Ebwt loop
fchr[A]: 0
fchr[C]: 7802696
fchr[G]: 12398259
fchr[T]: 16992544
fchr[$]: 24788084
Exiting Ebwt::buildToDisk()
Returning from initFromVector
Wrote 12457710 bytes to primary EBWT file: bowtie2/combined_SC_SP.rev.1.bt2
Wrote 6197028 bytes to secondary EBWT file: bowtie2/combined_SC_SP.rev.2.bt2
Re-opening _in1 and _in2 as input streams
Returning from Ebwt constructor
Headers:
    len: 24788084
    bwtLen: 24788085
    sz: 6197021
    bwtSz: 6197022
    lineRate: 6
    offRate: 4
    offMask: 0xfffffff0
    ftabChars: 10
    eftabLen: 20
    eftabSz: 80
    ftabLen: 1048577
    ftabSz: 4194308
    offsLen: 1549256
    offsSz: 6197024
    lineSz: 64
    sideSz: 64
    sideBwtSz: 48
    sideBwtLen: 192
    numSides: 129105
    numLines: 129105
    ebwtTotLen: 8262720
    ebwtTotSz: 8262720
    color: 0
    reverse: 1
Total time for backward call to driver() for mirror index: 00:00:12
```
</details>
<br />
<br />

<a id="copy-files-to-rina-and-rachel"></a>
## Copy files to Rina and Rachel
<a id="get-situated-copy-files"></a>
### Get situated, copy files
<a id="code-7"></a>
#### Code
<details>
<summary><i>Code: Get situated, copy files</i></summary>

```bash
#!/bin/bash

cd "${HOME}/tsukiyamalab/kalavatt/genomes" ||
    echo "cd'ing failed; check on this..."

dir_SC_SP="combined_SC_SP/"
dir_Rachel="${HOME}/tsukiyamalab/Rachel"
dir_Rina="${HOME}/tsukiyamalab/Rina"

cp -r "${dir_SC_SP}" "${dir_Rachel}"
cp -r "${dir_SC_SP}" "${dir_Rina}"
```
</details>
<br />

<a id="printed-5"></a>
#### Printed
<details>
<summary><i>Printed: Get situated, copy files</i></summary>

```txt
❯ cd "${HOME}/tsukiyamalab/kalavatt/genomes" ||
>    echo "cd'ing failed; check on this..."


❯ dir_SC_SP="combined_SC_SP/"


❯ dir_Rachel="${HOME}/tsukiyamalab/Rachel"


❯ dir_Rina="${HOME}/tsukiyamalab/Rina"


❯ cp -r "${dir_SC_SP}" "${dir_Rachel}"
'combined_SC_SP/' -> '/home/kalavatt/tsukiyamalab/Rachel/combined_SC_SP'
'combined_SC_SP/gff3' -> '/home/kalavatt/tsukiyamalab/Rachel/combined_SC_SP/gff3'
'combined_SC_SP/gff3/saccharomyces_cerevisiae_R64-3-1_20210421.gff3.gz' -> '/home/kalavatt/tsukiyamalab/Rachel/combined_SC_SP/gff3/saccharomyces_cerevisiae_R64-3-1_20210421.gff3.gz'
'combined_SC_SP/gff3/combined_SC_SP.gff3.gz' -> '/home/kalavatt/tsukiyamalab/Rachel/combined_SC_SP/gff3/combined_SC_SP.gff3.gz'
'combined_SC_SP/gff3/Schizosaccharomyces_pombe_all_chromosomes.gff3.gz' -> '/home/kalavatt/tsukiyamalab/Rachel/combined_SC_SP/gff3/Schizosaccharomyces_pombe_all_chromosomes.gff3.gz'
'combined_SC_SP/bowtie2' -> '/home/kalavatt/tsukiyamalab/Rachel/combined_SC_SP/bowtie2'
'combined_SC_SP/bowtie2/combined_SC_SP.stdout.txt' -> '/home/kalavatt/tsukiyamalab/Rachel/combined_SC_SP/bowtie2/combined_SC_SP.stdout.txt'
'combined_SC_SP/bowtie2/combined_SC_SP.stderr.txt' -> '/home/kalavatt/tsukiyamalab/Rachel/combined_SC_SP/bowtie2/combined_SC_SP.stderr.txt'
'combined_SC_SP/bowtie2/combined_SC_SP.3.bt2' -> '/home/kalavatt/tsukiyamalab/Rachel/combined_SC_SP/bowtie2/combined_SC_SP.3.bt2'
'combined_SC_SP/bowtie2/combined_SC_SP.4.bt2' -> '/home/kalavatt/tsukiyamalab/Rachel/combined_SC_SP/bowtie2/combined_SC_SP.4.bt2'
'combined_SC_SP/bowtie2/combined_SC_SP.1.bt2' -> '/home/kalavatt/tsukiyamalab/Rachel/combined_SC_SP/bowtie2/combined_SC_SP.1.bt2'
'combined_SC_SP/bowtie2/combined_SC_SP.2.bt2' -> '/home/kalavatt/tsukiyamalab/Rachel/combined_SC_SP/bowtie2/combined_SC_SP.2.bt2'
'combined_SC_SP/bowtie2/combined_SC_SP.rev.1.bt2' -> '/home/kalavatt/tsukiyamalab/Rachel/combined_SC_SP/bowtie2/combined_SC_SP.rev.1.bt2'
'combined_SC_SP/bowtie2/combined_SC_SP.rev.2.bt2' -> '/home/kalavatt/tsukiyamalab/Rachel/combined_SC_SP/bowtie2/combined_SC_SP.rev.2.bt2'
'combined_SC_SP/fasta' -> '/home/kalavatt/tsukiyamalab/Rachel/combined_SC_SP/fasta'
'combined_SC_SP/fasta/combined_SC_SP.fa' -> '/home/kalavatt/tsukiyamalab/Rachel/combined_SC_SP/fasta/combined_SC_SP.fa'
'combined_SC_SP/fasta/S288C_reference_sequence_R64-3-1_20210421.fa.gz' -> '/home/kalavatt/tsukiyamalab/Rachel/combined_SC_SP/fasta/S288C_reference_sequence_R64-3-1_20210421.fa.gz'
'combined_SC_SP/fasta/Schizosaccharomyces_pombe_all_chromosomes.fa.gz' -> '/home/kalavatt/tsukiyamalab/Rachel/combined_SC_SP/fasta/Schizosaccharomyces_pombe_all_chromosomes.fa.gz'
'combined_SC_SP/fasta/combined_SC_SP.fa.gz' -> '/home/kalavatt/tsukiyamalab/Rachel/combined_SC_SP/fasta/combined_SC_SP.fa.gz'
'combined_SC_SP/fasta/combined_SC_SP.fa.fai' -> '/home/kalavatt/tsukiyamalab/Rachel/combined_SC_SP/fasta/combined_SC_SP.fa.fai'


❯ cp -r "${dir_SC_SP}" "${dir_Rina}"
'combined_SC_SP/' -> '/home/kalavatt/tsukiyamalab/Rina/combined_SC_SP'
'combined_SC_SP/gff3' -> '/home/kalavatt/tsukiyamalab/Rina/combined_SC_SP/gff3'
'combined_SC_SP/gff3/saccharomyces_cerevisiae_R64-3-1_20210421.gff3.gz' -> '/home/kalavatt/tsukiyamalab/Rina/combined_SC_SP/gff3/saccharomyces_cerevisiae_R64-3-1_20210421.gff3.gz'
'combined_SC_SP/gff3/combined_SC_SP.gff3.gz' -> '/home/kalavatt/tsukiyamalab/Rina/combined_SC_SP/gff3/combined_SC_SP.gff3.gz'
'combined_SC_SP/gff3/Schizosaccharomyces_pombe_all_chromosomes.gff3.gz' -> '/home/kalavatt/tsukiyamalab/Rina/combined_SC_SP/gff3/Schizosaccharomyces_pombe_all_chromosomes.gff3.gz'
'combined_SC_SP/bowtie2' -> '/home/kalavatt/tsukiyamalab/Rina/combined_SC_SP/bowtie2'
'combined_SC_SP/bowtie2/combined_SC_SP.stdout.txt' -> '/home/kalavatt/tsukiyamalab/Rina/combined_SC_SP/bowtie2/combined_SC_SP.stdout.txt'
'combined_SC_SP/bowtie2/combined_SC_SP.stderr.txt' -> '/home/kalavatt/tsukiyamalab/Rina/combined_SC_SP/bowtie2/combined_SC_SP.stderr.txt'
'combined_SC_SP/bowtie2/combined_SC_SP.3.bt2' -> '/home/kalavatt/tsukiyamalab/Rina/combined_SC_SP/bowtie2/combined_SC_SP.3.bt2'
'combined_SC_SP/bowtie2/combined_SC_SP.4.bt2' -> '/home/kalavatt/tsukiyamalab/Rina/combined_SC_SP/bowtie2/combined_SC_SP.4.bt2'
'combined_SC_SP/bowtie2/combined_SC_SP.1.bt2' -> '/home/kalavatt/tsukiyamalab/Rina/combined_SC_SP/bowtie2/combined_SC_SP.1.bt2'
'combined_SC_SP/bowtie2/combined_SC_SP.2.bt2' -> '/home/kalavatt/tsukiyamalab/Rina/combined_SC_SP/bowtie2/combined_SC_SP.2.bt2'
'combined_SC_SP/bowtie2/combined_SC_SP.rev.1.bt2' -> '/home/kalavatt/tsukiyamalab/Rina/combined_SC_SP/bowtie2/combined_SC_SP.rev.1.bt2'
'combined_SC_SP/bowtie2/combined_SC_SP.rev.2.bt2' -> '/home/kalavatt/tsukiyamalab/Rina/combined_SC_SP/bowtie2/combined_SC_SP.rev.2.bt2'
'combined_SC_SP/fasta' -> '/home/kalavatt/tsukiyamalab/Rina/combined_SC_SP/fasta'
'combined_SC_SP/fasta/combined_SC_SP.fa' -> '/home/kalavatt/tsukiyamalab/Rina/combined_SC_SP/fasta/combined_SC_SP.fa'
'combined_SC_SP/fasta/S288C_reference_sequence_R64-3-1_20210421.fa.gz' -> '/home/kalavatt/tsukiyamalab/Rina/combined_SC_SP/fasta/S288C_reference_sequence_R64-3-1_20210421.fa.gz'
'combined_SC_SP/fasta/Schizosaccharomyces_pombe_all_chromosomes.fa.gz' -> '/home/kalavatt/tsukiyamalab/Rina/combined_SC_SP/fasta/Schizosaccharomyces_pombe_all_chromosomes.fa.gz'
'combined_SC_SP/fasta/combined_SC_SP.fa.gz' -> '/home/kalavatt/tsukiyamalab/Rina/combined_SC_SP/fasta/combined_SC_SP.fa.gz'
'combined_SC_SP/fasta/combined_SC_SP.fa.fai' -> '/home/kalavatt/tsukiyamalab/Rina/combined_SC_SP/fasta/combined_SC_SP.fa.fai'
```
</details>
<br />
