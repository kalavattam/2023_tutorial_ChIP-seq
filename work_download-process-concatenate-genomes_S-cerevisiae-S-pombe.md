
`#work_download-process_genomes-S-cerevisiae-S-pombe.md`
<br />
<br />

<details>
<summary><font size="+2"><b><i>Table of Contents</i></b></font></summary>
<!-- MarkdownTOC -->

1. [Get situated](#get-situated)
    1. [Code](#code)
1. [Download *S. pombe* fastas, gff3s](#download-s-pombe-fastas-gff3s)
    1. [Code](#code-1)
    1. [Notes](#notes)
1. [Download *S. cerevisiae* fastas, gff3](#download-s-cerevisiae-fastas-gff3)
    1. [Code](#code-2)
1. [Prepare *S. pombe* fasta, gff3 for concatenation with *S. cerevisiae*](#prepare-s-pombe-fasta-gff3-for-concatenation-with-s-cerevisiae)
    1. [Code](#code-3)
1. [Prepare *S. cerevisiae* fasta, gff3 for concatenation with *S. pombe*](#prepare-s-cerevisiae-fasta-gff3-for-concatenation-with-s-pombe)
    1. [Code](#code-4)
1. [Concatenate processed fastas and gff3s in new directory `combined_SC_SP/`](#concatenate-processed-fastas-and-gff3s-in-new-directory-combined_sc_sp)
    1. [Code](#code-5)
1. [Create `bowtie2` indices for "`combined_SC_SP.fa.gz`"](#create-bowtie2-indices-for-combined_sc_spfagz)
    1. [Code](#code-6)
1. [Copy files to Rina and Rachel](#copy-files-to-rina-and-rachel)
    1. [Code](#code-7)

<!-- /MarkdownTOC -->
</details>
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
<summary><i>Code: Download *S. pombe* fastas, gff3s</i></summary>

```bash
#!/bin/bash

cd "${HOME}/genomes/" || echo "cd'ing failed; check on this..."
cd "${d_pombe}" || echo "cd'ing failed; check on this..."


#  Download fastas ------------------------------------------------------------
if [[ ! -d "fasta/" ]]; then mkdir "fasta/"; fi

cd "fasta/" || echo "cd'ing failed; check on this..."

u_fa="https://www.pombase.org/data/genome_sequence_and_features/genome_sequence"
unset f_fa && typeset -a f_fa=(
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
run_check=true
if ${run_check}; then .,; fi

#  How do the chromosome names look?
run_check=true
if ${run_check}; then
    zcat Schizosaccharomyces_pombe_all_chromosomes.fa.gz | grep "^>"
fi


#  Download gff3s -------------------------------------------------------------
if [[ ! -d "gff3" ]]; then mkdir "gff3"; fi

cd "gff3/" || echo "cd'ing failed; check on this..."

u_gff3="https://www.pombase.org/data/genome_sequence_and_features/gff3/"
unset f_gff3 && typeset -a f_gff3=(
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
run_check=true
if ${run_check}; then .,; fi

#  How do the chromosome names look?
run_check=true
if ${run_check}; then
    zcat Schizosaccharomyces_pombe_all_chromosomes.gff3.gz \
        | cut -f 1 \
        | sort \
        | uniq
fi
```
</details>
<br />

<a id="notes"></a>
### Notes
<details>
<summary><i>Notes: Download *S. pombe* fastas, gff3s</i></summary>

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
<summary><i>Code: Download *S. cerevisiae* fastas, gff3</i></summary>

```bash
#!/bin/bash

cd "${HOME}/genomes/" || echo "cd'ing failed; check on this..."
cd "${d_cerevisiae}" || echo "cd'ing failed; check on this..."


#  Download fastas, gff3 ------------------------------------------------------
u_tgz="http://sgd-archive.yeastgenome.org/sequence/S288C_reference/genome_releases"
f_tgz="S288C_reference_genome_R64-3-1_20210421.tgz"

#  Download the tarfile of fastas and gff3
curl "${u_tgz}/${f_tgz}" > "${f_tgz}"

#  Decompress the directory
tar -xvzf "${f_tgz}"

#  Rename the directory to "fasta/"
mv "S288C_reference_genome_R64-3-1_20210421/" "fasta/"

#  Create a "gff3/" directory and move the gff3 file from "fasta/" to "gff3/"
mkdir gff3/ && \
    mv "fasta/saccharomyces_cerevisiae_R64-3-1_20210421.gff.gz" "gff3/"

#  Check how everything looks
run_check=true
if ${run_check}; then
    .,s
    .,
fi
```
</details>
<br />
<br />


<a id="prepare-s-pombe-fasta-gff3-for-concatenation-with-s-cerevisiae"></a>
## Prepare *S. pombe* fasta, gff3 for concatenation with *S. cerevisiae*
<a id="code-3"></a>
### Code
<details>
<summary><i>Code: Prepare *S. pombe* fasta, gff3 for concatenation with *S. cerevisiae*</i></summary>

```bash
#!/bin/bash

cd "${HOME}/genomes/" || echo "cd'ing failed; check on this..."
cd "${d_pombe}" || echo "cd'ing failed; check on this..."


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
run_check=true
if ${run_check}; then
    zgrep "^>" "Schizosaccharomyces_pombe_all_chromosomes.fa.gz"
fi

#  Create a decompressed version of the fasta
gzip -cd "Schizosaccharomyces_pombe_all_chromosomes.fa.gz" \
    > "Schizosaccharomyces_pombe_all_chromosomes.fa"

#  Use sed to rename chromosomes; save the results to "tmp.fa"
if [[ -f "tmp.fa" ]]; then rm "tmp.fa"; fi
sed 's/^>chr_II_telomeric_gap\ Schizosaccharomyces_pombe/>SP_II_TG/g;s/^>I\ Schizosaccharomyces_pombe/>SP_I/g;s/^>II\ Schizosaccharomyces_pombe/>SP_II/g;s/^>III\ Schizosaccharomyces_pombe/>SP_III/g;s/^>mating_type_region\ Schizosaccharomyces_pombe/>SP_MTR/g;s/^>mitochondrial\ Schizosaccharomyces_pombe/>SP_Mito/g' "Schizosaccharomyces_pombe_all_chromosomes.fa" \
    > "tmp.fa"

#  Check the new chromosome names
run_check=true
if ${run_check}; then cat "tmp.fa" | grep "^>"; fi

#  Overwrite the initial file with the contents of "tmp.fa"
mv -f "tmp.fa" "Schizosaccharomyces_pombe_all_chromosomes.fa"

#  Double check chromosome names
cat "Schizosaccharomyces_pombe_all_chromosomes.fa" | grep "^>"

#  Remove the compressed initial file
rm *.gz

#  Compress the updated file (with new chromosome names)
gzip *.fa

#  Check the directory
run_check=true
if ${run_check}; then .,; fi

#  Check chromosome names again
run_check=true
if ${run_check}; then
    zcat "Schizosaccharomyces_pombe_all_chromosomes.fa.gz" | grep "^>"
fi


#  Process gff3 ---------------------------------------------------------------
cd .. && pwd

#  Create a directory for the processed gff3
if [[ ! -d "gff3-processed/" ]]; then mkdir "gff3-processed/"; fi

#  Copy over file to process, then cd into directory
cp "gff3/Schizosaccharomyces_pombe_all_chromosomes.gff3.gz" "gff3-processed/"

cd "gff3-processed/" || echo "cd'ing failed; check on this..."

#  What do chromosome names look like?
run_check=true
if ${run_check}; then
    zcat "Schizosaccharomyces_pombe_all_chromosomes.gff3.gz" \
        | cut -f 1 \
        | sort \
        | uniq
fi

#  Use sed to rename chromosomes; save the results to "tmp.gff3"
zcat "Schizosaccharomyces_pombe_all_chromosomes.gff3.gz" \
    | sed 's/^chr_II_telomeric_gap/SP_II_TG/g;s/^I/SP_I/g;s/^II/SP_II/g;s/^III/SP_III/g;s/^mating_type_region/SP_MTR/g;s/^mitochondrial/SP_Mito/g' \
        > "tmp.gff3"

#  Check on the file contents
run_check=true
if ${run_check}; then
    head "tmp.gff3"
    tail "tmp.gff3"
    zcat "../gff3/Schizosaccharomyces_pombe_all_chromosomes.gff3.gz" | tail
fi

#  Check on the chromosome names in "tmp.gff3"
run_check=true
if ${run_check}; then cat "tmp.gff3" | cut -f 1 | sort | uniq; fi

#  Rename "tmp.gff3" to "Schizosaccharomyces_pombe_all_chromosomes.gff3"
mv "tmp.gff3" "Schizosaccharomyces_pombe_all_chromosomes.gff3"

#  Remove the initial gzipped file
rm "Schizosaccharomyces_pombe_all_chromosomes.gff3.gz"

#  Compress "Schizosaccharomyces_pombe_all_chromosomes.gff3"
gzip "Schizosaccharomyces_pombe_all_chromosomes.gff3"

#  Check the directory contents
run_check=true
if ${run_check}; then .,; fi

#  Check chromosome names again
run_check=true
if ${run_check}; then
    zcat "Schizosaccharomyces_pombe_all_chromosomes.gff3.gz" \
        | cut -f 1 \
        | sort \
        | uniq
fi
    ```
</details>
<br />
<br />

<a id="prepare-s-cerevisiae-fasta-gff3-for-concatenation-with-s-pombe"></a>
## Prepare *S. cerevisiae* fasta, gff3 for concatenation with *S. pombe*
<a id="code-4"></a>
### Code
<details>
<summary><i>Code: Prepare *S. cerevisiae* fasta, gff3 for concatenation with *S. pombe*</i></summary>

```bash
#!/bin/bash

cd "${HOME}/genomes/" || echo "cd'ing failed; check on this..."
cd "${d_cerevisiae}" || echo "cd'ing failed; check on this..."


#  Process fasta --------------------------------------------------------------
if [[ ! -d "fasta-processed" ]]; then mkdir "fasta-processed"; fi

#  Copy over file to process, then cd into directory
cp \
    "fasta/S288C_reference_sequence_R64-3-1_20210421.fsa.gz" \
    "fasta-processed/"

cd "fasta-processed" || echo "cd'ing failed; check on this..."

#  Check the chromosome names in the fasta
zcat "S288C_reference_sequence_R64-3-1_20210421.fsa.gz" | grep "^>"

#  Simplify the names of chromosomes in the S. cerevisiae fasta
if [[ -f "tmp.fa" ]]; then rm "tmp.fa"; fi
zcat "S288C_reference_sequence_R64-3-1_20210421.fsa.gz" \
    | sed 's/^>ref|NC_001133|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=I\]/>I/g;s/^>ref|NC_001134|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=II\]/>II/g;s/^>ref|NC_001135|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=III\]/>III/g;s/^>ref|NC_001136|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=IV\]/>IV/g;s/^>ref|NC_001137|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=V\]/>V/g;s/^>ref|NC_001138|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=VI\]/>VI/g;s/^>ref|NC_001139|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=VII\]/>VII/g;s/^>ref|NC_001140|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=VIII\]/>VIII/g;s/^>ref|NC_001141|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=IX\]/>IX/g;s/^>ref|NC_001142|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=X\]/>X/g;s/^>ref|NC_001143|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=XI\]/>XI/g;s/^>ref|NC_001144|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=XII\]/>XII/g;s/^>ref|NC_001145|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=XIII\]/>XIII/g;s/^>ref|NC_001146|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=XIV\]/>XIV/g;s/^>ref|NC_001147|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=XV\]/>XV/g;s/^>ref|NC_001148|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[chromosome=XVI\]/>XVI/g;s/^>ref|NC_001224|\ \[org=Saccharomyces\ cerevisiae\]\ \[strain=S288C\]\ \[moltype=genomic\]\ \[location=mitochondrion\]\ \[top=circular\]/>Mito/g' \
        > "tmp.fa"

#  Check the new chromosome names in "tmp.fa"
run_check=true
if ${run_check}; then cat "tmp.fa" | grep "^>"; fi

#  Remove the intial infile
rm "S288C_reference_sequence_R64-3-1_20210421.fsa.gz"

#  Rename "tmp.fa" to "S288C_reference_sequence_R64-3-1_20210421.fa"
mv -f "tmp.fa" "S288C_reference_sequence_R64-3-1_20210421.fa"

#  Compress "S288C_reference_sequence_R64-3-1_20210421.fa"
gzip "S288C_reference_sequence_R64-3-1_20210421.fa"

#  Check the directory contents
run_check=true
if ${run_check}; then .,; fi

#  Check the chromosome names
run_check=true
if ${run_check}; then
    zcat "S288C_reference_sequence_R64-3-1_20210421.fa.gz" | grep "^>"
fi


#  Process gff3 ---------------------------------------------------------------
cd .. && pwd

#  Create a directory for the processed gff3
if [[ ! -d "gff3-processed/" ]]; then mkdir "gff3-processed/"; fi

#  Copy over file to process, then cd into directory
cp "gff3/saccharomyces_cerevisiae_R64-3-1_20210421.gff.gz" "gff3-processed/"

cd "gff3-processed/" || echo "cd'ing failed; check on this..."

#  What do chromosome names look like?
run_check=true
if ${run_check}; then
    zcat "saccharomyces_cerevisiae_R64-3-1_20210421.gff.gz" \
        | cut -f 1 \
        | sort \
        | uniq
fi
#NOTE There are fasta sequences in this file that need to be excluded

#  Decompress "saccharomyces_cerevisiae_R64-3-1_20210421.gff.gz"
gzip -cd "saccharomyces_cerevisiae_R64-3-1_20210421.gff.gz" \
    > "saccharomyces_cerevisiae_R64-3-1_20210421.gff"

#  Remove everything after line containing ### (fasta chromosome sequences)
if [[ -f "tmp.gff3" ]]; then rm "tmp.gff3"; fi
sed -n '/###/q;p' < saccharomyces_cerevisiae_R64-3-1_20210421.gff > tmp.gff3

#  Check the chromosome names now that fasta sequences are gone
run_check=true
if ${run_check}; then cat "tmp.gff3" | cut -f 1 | sort | uniq; fi

#  Use sed to rename chromosomes; save the results to "tmp.gff3"
cat "tmp.gff3" \
    | sed 's/^chr//g;s/^mt/Mito/g' \
        > "tmp.2.gff3"

#  Check the chromosome names in the updated file
run_check=true
if ${run_check}; then cat "tmp.2.gff3" | cut -f 1 | sort | uniq; fi

#  Check on the file contents
head "tmp.2.gff3"
tail "tmp.2.gff3"

#  Rename "tmp.2.gff3" to "saccharomyces_cerevisiae_R64-3-1_20210421.gff3"
mv "tmp.2.gff3" "saccharomyces_cerevisiae_R64-3-1_20210421.gff3"

#  Compress "saccharomyces_cerevisiae_R64-3-1_20210421.gff3"
gzip "saccharomyces_cerevisiae_R64-3-1_20210421.gff3"

#  Check the directory contents
run_check=true
if ${run_check}; then .,; fi

#  Remove unneeded files
rm \
    "saccharomyces_cerevisiae_R64-3-1_20210421.gff" \
    "saccharomyces_cerevisiae_R64-3-1_20210421.gff.gz" \
    "tmp.gff3"

#  Check the directory contents again
run_check=true
if ${run_check}; then .,; fi

#  Check chromosome names again
run_check=true
if ${run_check}; then
    zcat "saccharomyces_cerevisiae_R64-3-1_20210421.gff3.gz" \
        | cut -f 1 \
        | sort \
        | uniq
fi
```
</details>
<br />
<br />

<a id="concatenate-processed-fastas-and-gff3s-in-new-directory-combined_sc_sp"></a>
## Concatenate processed fastas and gff3s in new directory `combined_SC_SP/`
<a id="code-5"></a>
### Code
<details>
<summary><i>Code: Concatenate processed fastas and gff3s in new directory `combined_SC_SP/`</i></summary>

```bash
#!/bin/bash

cd "${HOME}/tsukiyamalab/kalavatt/genomes/" \
    || echo "cd'ing failed; check on this..."


#  Get situated ---------------------------------------------------------------
[[ ! -d "combined_SC_SP" ]] && mkdir -p combined_SC_SP/{fasta,gff3} || true

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


#  gff3 -----------------------------------------------------------------------
cd "combined_SC_SP/gff3/" || echo "cd'ing failed; check on this..."

[[ -f "combined_SC_SP.gff3.gz" ]] && rm "combined_SC_SP.gff3.gz" || true

cat \
    "saccharomyces_cerevisiae_R64-3-1_20210421.gff3.gz" \
    "Schizosaccharomyces_pombe_all_chromosomes.gff3.gz" \
        > "combined_SC_SP.gff3.gz"

#  ARS_consensus_sequence records have strands labeled "0"; these need to be
#+ adjusted or else downstream programs will throw errors when encountering the
#+ "0" strands
zcat "combined_SC_SP.gff3.gz" \
    | awk -F "\t" '
        BEGIN { OFS = FS } {
            if ($7=="0") { $7="."; print $0 } else { print $0 }
        }
    ' \
    | gzip \
        > "tmp.gff3.gz"

mv -f "tmp.gff3.gz" "combined_SC_SP.gff3.gz"

#  Replace or remove special characters that IGV can't handle
zcat "combined_SC_SP.gff3.gz" \
    | sed -e 's/%20/-/g' \
          -e 's/%2C//g' \
          -e 's/%3B//g' \
          -e 's/%28//g' \
          -e 's/%29//g' \
          -e 's/%//g' \
    | gzip \
        > "combined_SC_SP.clean.gff3.gz"

#  Create S. cerevisiae/S. pombe concatenated gff3 file with "intelligible
#+ feature names": "Name=" value is "ID=" value; this needs to be something
#+ that makes sense to a human
#+ 
#+ Need to write an `awk` command to
#+ 1. If field `$3` is `gene`, `blocked_reading_frame`, `ncRNA_gene`,
#+ `pseudogene`, `rRNA_gene`, `snRNA_gene`, `snoRNA_gene`, `tRNA_gene`,
#+ `telomerase_RNA_gene` then check for the presence of the `gene=` attribute;
#+ if present, then replace the `ID=` value with the `gene=` value
#+ 2. If field `$3` is `ARS`, then check for the presence of the `Alias=`
#+ attribute; if present, then replace the `ID=` value with the `Alias=` value
zcat "combined_SC_SP.clean.gff3.gz" \
    | awk '
        BEGIN { FS=OFS="\t" }
        $1 ~ /^#/ { print; next }  # Skip lines starting with #
        $3 ~ /^(blocked_reading_frame|gene|ncRNA_gene|pseudogene|rRNA_gene|snRNA_gene|snoRNA_gene|tRNA_gene|telomerase_RNA_gene)$/ {
            if (match($9, /gene=[^;]+/)) {
                #  Extract the gene name
                gene_name = substr($9, RSTART+5, RLENGTH-5)
                #  Replace Name value with gene value
                sub(/Name=[^;]+/, "Name=" gene_name, $9)
            }
        }
        $3 == "ARS" {
            if (match($9, /Alias=[^;]+/)) {
                #  Extract the alias name
                alias_name = substr($9, RSTART+6, RLENGTH-6)
                #  Replace Name value with alias value
                sub(/Name=[^;]+/, "Name=" alias_name, $9)
            }
        }
        { print }
    ' \
    | gzip \
        > "combined_SC_SP.clean.intelligible.gff3.gz"

run_check=true
if ${run_check}; then .,; fi

run_check=false
if ${run_check}; then
    zcat "combined_SC_SP.clean.intelligible.gff3.gz" | cut -f 1 | sort | uniq
fi


#  fasta ----------------------------------------------------------------------
cd "../fasta" || echo "cd'ing failed; check on this..."

cat \
    "S288C_reference_sequence_R64-3-1_20210421.fa.gz" \
    "Schizosaccharomyces_pombe_all_chromosomes.fa.gz" \
        > "combined_SC_SP.fa.gz"

run_check=true
if ${run_check}; then ls -lhaFG; fi

run_check=true
if ${run_check}; then zcat "combined_SC_SP.fa.gz" | grep "^>"; fi
```
</details>
<br />
<br />

<a id="create-bowtie2-indices-for-combined_sc_spfagz"></a>
## Create `bowtie2` indices for "`combined_SC_SP.fa.gz`"
<a id="code-6"></a>
### Code
<details>
<summary><i>Code: Create `bowtie2` indices for "`combined_SC_SP.fa.gz`"</i></summary>

```bash
#!/bin/bash

cd "${HOME}/genomes/combined_SC_SP" || echo "cd'ing failed; check on this..."

if [[ ! -d "bowtie2/" ]]; then mkdir "bowtie2/"; fi

#  Index the fasta file
cd "fasta/"

gzip -cd "combined_SC_SP.fa.gz" > "combined_SC_SP.fa"

module load SAMtools/1.16.1-GCC-11.2.0 Bowtie2/2.4.4-GCC-11.2.0

run_check=true
if ${run_check}; then cat "combined_SC_SP.fa" | grep "^>"; fi

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
<br />

<a id="copy-files-to-rina-and-rachel"></a>
## Copy files to Rina and Rachel
<a id="code-7"></a>
### Code
<details>
<summary><i>Code: Copy files to Rina and Rachel</i></summary>

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
