
`#set-up_siQ-ChIP_assessment.md`
<br />
<br />

## Clone the forked repository for siQ-ChIP
### Code
<details>
<summary><i>Code: Clone the forked repository for siQ-ChIP</i></summary>

```bash
#!/bin/bash

cd "${HOME}/repos/2023_rDNA/src" \
    || echo "cd'ing failed; check on this"

#  Clone forked siQ-ChIP repository
if [[ ! -d siQ-ChIP ]]; then
    gh repo clone kalavattam/siQ-ChIP
fi

cd siQ-ChIP || echo "cd'ing failed; check on this"

git checkout -b cerevisiae

git push origin cerevisiae
```
</details>
<br />

### Printed
<details>
<summary><i>Printed: Clone the forked repository for siQ-ChIP</i></summary>

```txt
❯ cd "${HOME}/repos/2023_rDNA/src" \
>     || echo "cd'ing failed; check on this"


❯ if [[ ! -d siQ-ChIP ]]; then
then>     gh repo clone kalavattam/siQ-ChIP
then> fi
Cloning into 'siQ-ChIP'...
remote: Enumerating objects: 155, done.
remote: Counting objects: 100% (114/114), done.
remote: Compressing objects: 100% (106/106), done.
remote: Total 155 (delta 56), reused 27 (delta 7), pack-reused 41
Receiving objects: 100% (155/155), 33.79 MiB | 8.76 MiB/s, done.
Resolving deltas: 100% (73/73), done.
Updating upstream
From https://github.com/BradleyDickson/siQ-ChIP
 * [new branch]      master     -> upstream/master


A new release of gh is available: 2.11.3 → 2.39.2
To upgrade, run: brew update && brew upgrade gh
https://github.com/cli/cli/releases/tag/v2.39.2


❯ cd siQ-ChIP && git checkout -b cerevisiae
Switched to a new branch 'cerevisiae'


❯ git push origin cerevisiae
Total 0 (delta 0), reused 0 (delta 0), pack-reused 0
remote:
remote: Create a pull request for 'cerevisiae' on GitHub by visiting:
remote:      https://github.com/kalavattam/siQ-ChIP/pull/new/cerevisiae
remote:
To https://github.com/kalavattam/siQ-ChIP.git
 * [new branch]      cerevisiae -> cerevisiae
```
</details>
<br />

## [Dependencies]((https://github.com/kalavattam/siQ-ChIP#dependencies-and-assumptions))
### Text
<details>
<summary><i>Text: Dependencies</i></summary>

- [bc](https://anaconda.org/conda-forge/bc)
- [gfortran](https://anaconda.org/conda-forge/gfortran)
- [gnuplot](https://anaconda.org/conda-forge/gnuplot)
</details>
<br />

## Lines of code that need to be changed
### Text
<details>
<summary><i>Text: Lines of code that need to be changed</i></summary>

- [`WGfrechet.sh`](https://github.com/kalavattam/siQ-ChIP/blob/master/WGfrechet.sh)
    + [`L12`](https://github.com/kalavattam/siQ-ChIP/blob/master/WGfrechet.sh#L12)
    + [`L15-L16`](https://github.com/kalavattam/siQ-ChIP/blob/master/WGfrechet.sh#L15-L16)
    + [`L22-L34`](https://github.com/kalavattam/siQ-ChIP/blob/master/WGfrechet.sh#L22-L34)
    + `#MAYBE` [`L42-L43`](https://github.com/kalavattam/siQ-ChIP/blob/master/WGfrechet.sh#L42-L43)
- [`frechet.f90`](https://github.com/kalavattam/siQ-ChIP/blob/master/frechet.f90): It's not clear to me that anything needs to be changed here&mdash;but I am not sure
- [`Bins.f90`](https://github.com/kalavattam/siQ-ChIP/blob/master/Bins.f90): Again, it's not clear to me that anything needs to be changed here&mdash;but I am not sure
- [`getalpha.f90`](https://github.com/kalavattam/siQ-ChIP/blob/master/getalpha.f90): Again, it's not clear to me that anything needs to be changed here&mdash;but I am not sure
- [`binReads.f90`](https://github.com/kalavattam/siQ-ChIP/blob/master/binReads.f90): It's not clear to me that anything needs to be changed here&mdash;but I am not sure
</details>
<br />

## Working through the [`README`](https://github.com/BradleyDickson/siQ-ChIP#readme)
### Text
<details>
<summary><i>Text: Working through the `README`</i></summary>

#### Philosophy of use
...

To give you a quick sense of what siQ will require, you will need the following starting ingredients. We will cover all details for each of these points below.
1. bed files of aligned sequencing data
2. parameter files for siQ scaling (IP mass, input mass, average fragment length)
3. build EXPlayout file for your experiment (to define which ip, input, and parameters go together)
4. Link any annotations you want:  `ln -s PATH/your_annotations.bed ./Annotations.bed` ---> You gotta use `Annotations.bed `for the name you link to! If you have none, use `touch Annotations.bed`
5. Execute `./getsiq.sh > errorfile` or whatever is appropriate for your HPC

Each of these steps (save for generating your aligned bed files) is discussed below.

#### To perform siQ-ChIP
At this point you have determined your antibody:chromatin isotherm and managed to demonstrate clear observation of signal (captured DNA mass). Or maybe you just did ChIP without this isotherm, and that's ok. You can still apply the analysis here, but you inherit some caveats. Your samples have all been sequenced and you have mapped your data to your target genome. You will need the bed files from your alignment and you will need to prepare the following parameter files for all of your samples. (Bed files are to be sorted as usual with `sort -k1,1 -k2,2n`)

<b>Bed file format</b>: A bed file containing all QC'd paired-end reads for an IP and an INPUT sample. An example of the first few lines of one of these files are as follows where the chr, start, stop and length of reads is listed:

```txt
chr1    100000041       100000328       287
chr1    100000189       100000324       135
chr1    10000021        10000169        148
chr1    100000389       100000596       207
chr1    100000748       100001095       347
chr1    100000917       100001015       98
chr1    100000964       100001113       149
chr1    10000122        10000449        327
chr1    100001232       100001602       370
```

#### Parameter files
Each ChIP reaction has its own parameter file that contains the information required to compute the siQ-ChIP quantitative scale. Each parameter file <b><i>must</i></b> have the following information in the following format (example given below):
```txt
input sample volume
total volume before removal of input
input DNA mass (ng)
IP DNA mass (ng)
IP average fragment length (from Bioanalyzer)
input average fragment length (from Bioanalyzer)
```

An example file would look like this:
```txt
50
250
135
10
400
382
```

You may take input sample, then split the remaining chromatin for separate IPs. That's fine. Just be sure to enter the IP volume plus the input volume for the total volume in this parameter file.

At this point, you have a parameter file for each of your samples and you have a bed file for each sample (IP and input). Next, you need to build a "layout file" to tell the siQ-ChIP scripts which files go together and which samples should be quantitatively compared.

#### The EXPlayout file
siQ-ChIP enables the comparison of two or more ChIP-seq experiments. So we assume you have two IP datasets, two input datasets, and two sets of measurements required to evaluate the quantitative scale for each of these IP-cases. It is fine if you have one input that was shared for different IP.

The siQ-ChIP track for experiment 1 is built by combining IP1.bed input1.bed params1.in Likewise, the second experiment is processed using IP2.bed input2.bed params2.in. This is to say that the IP, input, and measurements (params) will be integrated to produce one track (at quantified scale) for each experiment.

To build all the siQ-ChIP tracks and to compare annotated fragment distributions and tracks, we only need to build the following EXPlayout file and make sure our params files are defined (see below). <b><i>No dashes in file names.</i></b>
```txt
#getTracks: IP.bed input.bed params output_name
IP1.bed input1.bed params1.in exp1siq
IP2.bed input2.bed params2.in exp2siq
#getResponse: CNTR.bed EXP.bed output_name
exp1siq.bed exp2siq.bed exp1TOexp2
#getFracts: data any order, last is output_name
IP1.bed IP2.bed input1.bed input2.bed MyExperiment
#END
```

If you only have one IP and one input, the you may make an EXPlayout like this:
```txt
#getTracks: IP.bed input.bed params output_name
IP1.bed input1.bed params1.in exp1siq
#getResponse: CNTR.bed EXP.bed output_name
#getFracts: data any order, last is output_name
#END
```

This will only build the siQ-scaled track that you have data for.

Likewise, you could compare two tracks that you've already built with siQ, as follows:
```txt
#getTracks: IP.bed input.bed params output_name
#getResponse: CNTR.bed EXP.bed output_name
exp1siq.bed exp2siq.bed exp1TOexp2
#getFracts: data any order, last is output_name
#END
```


</details>
<br />
