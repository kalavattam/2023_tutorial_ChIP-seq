
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
- [`frechet.f90`](https://github.com/kalavattam/siQ-ChIP/blob/master/frechet.f90): It's not clear to me that anything needs to be changed here&mdash;but I am not sure.
- [`Bins.f90`](https://github.com/kalavattam/siQ-ChIP/blob/master/Bins.f90): Again, it's not clear to me that anything needs to be changed here&mdash;but I am not sure.
- [`getalpha.f90`](https://github.com/kalavattam/siQ-ChIP/blob/master/getalpha.f90): Again, it's not clear to me that anything needs to be changed here&mdash;but I am not sure.
- [`binReads.f90`](https://github.com/kalavattam/siQ-ChIP/blob/master/binReads.f90): It's not clear to me that anything needs to be changed here&mdash;but I am not sure.
</details>
<br />

## Working through the [`README`](https://github.com/BradleyDickson/siQ-ChIP#readme)
### Text
<details>
<summary><i>Text: Working through the `README`</i></summary>

#### Philosophy of use
...

To give you a quick sense of what siQ will require, you will need the following starting ingredients. We will cover all details for each of these points below.
1. <mark>bed files of aligned sequencing data</mark>
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

This is useful if you thought of a comparison to make after you built your siQ scaled tracks, or if you acquired new data at a later time and don't need to rebuild all siQ-tracks.

The key here is that you need these section headers in EXPlayout, but the sections can be empty. Finally, you have to use the name EXPlayout. If you have to redo something or add a new comparison, save your EXPlayout to a meaningful name. Then edit it. This might get improved at a later time.

The getFracts section outputs datafiles named MyExperiment.
</details>
<br />

## Create siQ-ChIP environment (draft)
### Code
<details>
<summary><i>Code: Create siQ-ChIP environment (draft)</i></summary>

```bash
#!/bin/bash

unset envs && typeset -a envs
while IFS=$'\n' read -r line; do
    env_name=$(echo "${line}" | awk '{ print $1 }')
    envs+=( "${env_name}" )
done < <(conda info -e | grep -v "^#")
# echo_test "${envs[@]}"

EOI="siQ-ChIP_env"
found=0

for env in "${envs[@]}"; do
    if [[ "${env}" == *"${EOI}"* ]]; then
        echo "Env found amidst conda envs: ${env}."
        found=1
        break
    fi
done

if [[ ${found} -eq 0 ]]; then
    echo "Env not found amidst conda envs. Installing ${EOI}."
    mamba create \
        -n "${EOI}" \
        -c conda-forge \
        -c bioconda \
            bc \
            bedtools \
            gfortran \
            gnuplot
fi
```
</details>
<br />

### Printed
<details>
<summary><i>Printed: Create siQ-ChIP environment (draft)</i></summary>

```txt
❯ unset envs && typeset -a envs
❯ while IFS=$'\n' read -r line; do
>     env_name=$(echo "${line}" | awk '{ print $1 }')
>     envs+=( "${env_name}" )
> done < <(conda info -e | grep -v "^#")


❯ EOI="siQ-ChIP_env"
❯ found=0


❯ for env in "${envs[@]}"; do
>     if [[ "${env}" == *"${EOI}"* ]]; then
>         echo "Env found amidst conda envs: ${env}."
>         found=1
>         break
>     fi
> done
 

❯ if [[ ${found} -eq 0 ]]; then
>     echo "Env not found amidst conda envs. Installing ${EOI}."
>     mamba create \
>         -n "${EOI}" \
>         -c conda-forge \
>         -c bioconda \
>             bc \
>             bedtools \
>             gfortran \
>             gnuplot
> fi
Env not found amidst conda envs. Installing .

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


Looking for: ['bc', 'bedtools', 'gfortran', 'gnuplot']

bioconda/noarch                                      4.9MB @   3.4MB/s  1.6s
bioconda/linux-64                                    5.2MB @   2.2MB/s  2.5s
pkgs/main/linux-64                                   6.4MB @   2.7MB/s  2.8s
pkgs/r/linux-64                                               No change
pkgs/main/noarch                                              No change
pkgs/r/noarch                                                 No change
conda-forge/noarch                                  14.9MB @   4.5MB/s  3.9s
conda-forge/linux-64                                36.5MB @   5.6MB/s  7.2s
Transaction

  Prefix: /home/kalavatt/miniconda3/envs/siQ-ChIP_env

  Updating specs:

   - bc
   - bedtools
   - gfortran
   - gnuplot


  Package                           Version  Build               Channel                    Size
──────────────────────────────────────────────────────────────────────────────────────────────────
  Install:
──────────────────────────────────────────────────────────────────────────────────────────────────

  + _libgcc_mutex                       0.1  conda_forge         conda-forge/linux-64     Cached
  + _openmp_mutex                       4.5  2_gnu               conda-forge/linux-64     Cached
  + alsa-lib                         1.2.10  hd590300_0          conda-forge/linux-64     Cached
  + atk-1.0                          2.38.0  hd4edc92_1          conda-forge/linux-64     Cached
  + attr                              2.5.1  h166bdaf_1          conda-forge/linux-64     Cached
  + bc                               1.07.1  h7f98852_0          conda-forge/linux-64      103kB
  + bedtools                         2.31.1  hf5e1c6e_0          bioconda/linux-64        Cached
  + binutils_impl_linux-64             2.40  hf600244_0          conda-forge/linux-64     Cached
  + bzip2                             1.0.8  hd590300_5          conda-forge/linux-64     Cached
  + ca-certificates              2023.11.17  hbcca054_0          conda-forge/linux-64     Cached
  + cairo                            1.18.0  h3faef2a_0          conda-forge/linux-64     Cached
  + chrpath                            0.16  h7f98852_1002       conda-forge/linux-64       30kB
  + dbus                             1.13.6  h5008d03_3          conda-forge/linux-64     Cached
  + expat                             2.5.0  hcb278e6_1          conda-forge/linux-64     Cached
  + font-ttf-dejavu-sans-mono          2.37  hab24e00_0          conda-forge/noarch       Cached
  + font-ttf-inconsolata              3.000  h77eed37_0          conda-forge/noarch       Cached
  + font-ttf-source-code-pro          2.038  h77eed37_0          conda-forge/noarch       Cached
  + font-ttf-ubuntu                    0.83  h77eed37_1          conda-forge/noarch       Cached
  + fontconfig                       2.14.2  h14ed4e7_0          conda-forge/linux-64     Cached
  + fonts-conda-ecosystem                 1  0                   conda-forge/noarch       Cached
  + fonts-conda-forge                     1  0                   conda-forge/noarch       Cached
  + freetype                         2.12.1  h267a509_2          conda-forge/linux-64     Cached
  + fribidi                          1.0.10  h36c2ea0_0          conda-forge/linux-64     Cached
  + gcc                              13.2.0  h574f8da_2          conda-forge/linux-64       27kB
  + gcc_impl_linux-64                13.2.0  h338b0a0_3          conda-forge/linux-64     Cached
  + gdk-pixbuf                      2.42.10  h829c605_4          conda-forge/linux-64      572kB
  + gettext                          0.21.1  h27087fc_0          conda-forge/linux-64     Cached
  + gfortran                         13.2.0  h0584b13_2          conda-forge/linux-64       27kB
  + gfortran_impl_linux-64           13.2.0  h76e1118_3          conda-forge/linux-64     Cached
  + giflib                            5.2.1  h0b41bf4_3          conda-forge/linux-64     Cached
  + glib                             2.78.2  hfc55251_0          conda-forge/linux-64      488kB
  + glib-tools                       2.78.2  hfc55251_0          conda-forge/linux-64      112kB
  + gnuplot                           5.4.8  h142138f_0          conda-forge/linux-64        1MB
  + graphite2                        1.3.13  h58526e2_1001       conda-forge/linux-64     Cached
  + gst-plugins-base                 1.22.7  h8e1006c_0          conda-forge/linux-64     Cached
  + gstreamer                        1.22.7  h98fc4e7_0          conda-forge/linux-64     Cached
  + gtk2                            2.24.33  h90689f9_2          conda-forge/linux-64     Cached
  + harfbuzz                          8.3.0  h3d44ed6_0          conda-forge/linux-64     Cached
  + icu                                73.2  h59595ed_0          conda-forge/linux-64     Cached
  + kernel-headers_linux-64          2.6.32  he073ed8_16         conda-forge/noarch       Cached
  + keyutils                          1.6.1  h166bdaf_0          conda-forge/linux-64     Cached
  + krb5                             1.21.2  h659d440_0          conda-forge/linux-64     Cached
  + lame                              3.100  h166bdaf_1003       conda-forge/linux-64     Cached
  + ld_impl_linux-64                   2.40  h41732ed_0          conda-forge/linux-64     Cached
  + lerc                              4.0.0  h27087fc_0          conda-forge/linux-64     Cached
  + libcap                             2.69  h0f662aa_0          conda-forge/linux-64     Cached
  + libclang                         15.0.7  default_hb11cfb5_4  conda-forge/linux-64     Cached
  + libclang13                       15.0.7  default_ha2b6cf4_4  conda-forge/linux-64     Cached
  + libcups                           2.3.3  h4637d8d_4          conda-forge/linux-64     Cached
  + libdeflate                         1.19  hd590300_0          conda-forge/linux-64     Cached
  + libedit                    3.1.20191231  he28a2e2_2          conda-forge/linux-64     Cached
  + libevent                         2.1.12  hf998b51_1          conda-forge/linux-64     Cached
  + libexpat                          2.5.0  hcb278e6_1          conda-forge/linux-64     Cached
  + libffi                            3.4.2  h7f98852_5          conda-forge/linux-64     Cached
  + libflac                           1.4.3  h59595ed_0          conda-forge/linux-64     Cached
  + libgcc-devel_linux-64            13.2.0  ha9c7c90_103        conda-forge/noarch       Cached
  + libgcc-ng                        13.2.0  h807b86a_3          conda-forge/linux-64     Cached
  + libgcrypt                        1.10.3  hd590300_0          conda-forge/linux-64     Cached
  + libgd                             2.3.3  h119a65a_9          conda-forge/linux-64      224kB
  + libgfortran5                     13.2.0  ha4646dd_3          conda-forge/linux-64     Cached
  + libglib                          2.78.2  h783c2da_0          conda-forge/linux-64        3MB
  + libgomp                          13.2.0  h807b86a_3          conda-forge/linux-64     Cached
  + libgpg-error                       1.47  h71f35ed_0          conda-forge/linux-64     Cached
  + libiconv                           1.17  h166bdaf_0          conda-forge/linux-64     Cached
  + libjpeg-turbo                     3.0.0  hd590300_1          conda-forge/linux-64     Cached
  + libllvm15                        15.0.7  h5cf9203_3          conda-forge/linux-64     Cached
  + libnsl                            2.0.1  hd590300_0          conda-forge/linux-64     Cached
  + libogg                            1.3.4  h7f98852_1          conda-forge/linux-64     Cached
  + libopus                           1.3.1  h7f98852_1          conda-forge/linux-64     Cached
  + libpng                           1.6.39  h753d276_0          conda-forge/linux-64     Cached
  + libpq                              16.1  hfc447b1_0          conda-forge/linux-64        3MB
  + libsanitizer                     13.2.0  h7e041cc_3          conda-forge/linux-64     Cached
  + libsndfile                        1.2.2  hc60ed4a_1          conda-forge/linux-64     Cached
  + libsqlite                        3.44.2  h2797004_0          conda-forge/linux-64     Cached
  + libstdcxx-ng                     13.2.0  h7e041cc_3          conda-forge/linux-64     Cached
  + libsystemd0                         255  h3516f8a_0          conda-forge/linux-64      404kB
  + libtiff                           4.6.0  ha9c0a0a_2          conda-forge/linux-64     Cached
  + libuuid                          2.38.1  h0b41bf4_0          conda-forge/linux-64     Cached
  + libvorbis                         1.3.7  h9c3ff4c_0          conda-forge/linux-64     Cached
  + libwebp                           1.3.2  h658648e_1          conda-forge/linux-64       85kB
  + libwebp-base                      1.3.2  hd590300_0          conda-forge/linux-64     Cached
  + libxcb                             1.15  h0b41bf4_0          conda-forge/linux-64     Cached
  + libxkbcommon                      1.6.0  h5d7e998_0          conda-forge/linux-64     Cached
  + libxml2                          2.11.6  h232c23b_0          conda-forge/linux-64     Cached
  + libzlib                          1.2.13  hd590300_5          conda-forge/linux-64     Cached
  + lz4-c                             1.9.4  hcb278e6_0          conda-forge/linux-64     Cached
  + mpg123                           1.32.3  h59595ed_0          conda-forge/linux-64     Cached
  + mysql-common                     8.0.33  hf1915f5_6          conda-forge/linux-64     Cached
  + mysql-libs                       8.0.33  hca2cd23_6          conda-forge/linux-64     Cached
  + ncurses                             6.4  h59595ed_2          conda-forge/linux-64     Cached
  + nspr                               4.35  h27087fc_0          conda-forge/linux-64     Cached
  + nss                                3.95  h1d7d5a4_0          conda-forge/linux-64     Cached
  + openssl                           3.1.4  hd590300_0          conda-forge/linux-64     Cached
  + pango                           1.50.14  ha41ecd1_2          conda-forge/linux-64     Cached
  + pcre2                             10.42  hcad00b1_0          conda-forge/linux-64        1MB
  + pip                              23.3.1  pyhd8ed1ab_0        conda-forge/noarch       Cached
  + pixman                           0.42.2  h59595ed_0          conda-forge/linux-64     Cached
  + pthread-stubs                       0.4  h36c2ea0_1001       conda-forge/linux-64     Cached
  + pulseaudio-client                  16.1  hb77b528_5          conda-forge/linux-64     Cached
  + python                           3.12.0  hab00c5b_0_cpython  conda-forge/linux-64     Cached
  + qt-main                          5.15.8  h82b777d_17         conda-forge/linux-64       61MB
  + readline                            8.2  h8228510_1          conda-forge/linux-64     Cached
  + setuptools                       68.2.2  pyhd8ed1ab_0        conda-forge/noarch       Cached
  + sysroot_linux-64                   2.12  he073ed8_16         conda-forge/noarch       Cached
  + tk                               8.6.13  noxft_h4845f30_101  conda-forge/linux-64     Cached
  + tzdata                            2023c  h71feb2d_0          conda-forge/noarch       Cached
  + wheel                            0.42.0  pyhd8ed1ab_0        conda-forge/noarch       Cached
  + xcb-util                          0.4.0  hd590300_1          conda-forge/linux-64     Cached
  + xcb-util-image                    0.4.0  h8ee46fc_1          conda-forge/linux-64     Cached
  + xcb-util-keysyms                  0.4.0  h8ee46fc_1          conda-forge/linux-64     Cached
  + xcb-util-renderutil               0.3.9  hd590300_1          conda-forge/linux-64     Cached
  + xcb-util-wm                       0.4.1  h8ee46fc_1          conda-forge/linux-64     Cached
  + xkeyboard-config                   2.40  hd590300_0          conda-forge/linux-64     Cached
  + xorg-kbproto                      1.0.7  h7f98852_1002       conda-forge/linux-64     Cached
  + xorg-libice                       1.1.1  hd590300_0          conda-forge/linux-64     Cached
  + xorg-libsm                        1.2.4  h7391055_0          conda-forge/linux-64     Cached
  + xorg-libx11                       1.8.7  h8ee46fc_0          conda-forge/linux-64     Cached
  + xorg-libxau                      1.0.11  hd590300_0          conda-forge/linux-64     Cached
  + xorg-libxdmcp                     1.1.3  h7f98852_0          conda-forge/linux-64     Cached
  + xorg-libxext                      1.3.4  h0b41bf4_2          conda-forge/linux-64     Cached
  + xorg-libxrender                  0.9.11  hd590300_0          conda-forge/linux-64     Cached
  + xorg-libxt                        1.3.0  hd590300_1          conda-forge/linux-64     Cached
  + xorg-renderproto                 0.11.1  h7f98852_1002       conda-forge/linux-64     Cached
  + xorg-xextproto                    7.3.0  h0b41bf4_1003       conda-forge/linux-64     Cached
  + xorg-xf86vidmodeproto             2.3.1  h7f98852_1002       conda-forge/linux-64     Cached
  + xorg-xproto                      7.0.31  h7f98852_1007       conda-forge/linux-64     Cached
  + xz                                5.2.6  h166bdaf_0          conda-forge/linux-64     Cached
  + zlib                             1.2.13  hd590300_5          conda-forge/linux-64     Cached
  + zstd                              1.5.5  hfc55251_0          conda-forge/linux-64     Cached

  Summary:

  Install: 129 packages

  Total download: 71MB

──────────────────────────────────────────────────────────────────────────────────────────────────


Confirm changes: [Y/n] Y
pcre2                                                1.0MB @  11.6MB/s  0.1s
libglib                                              2.7MB @  27.4MB/s  0.1s
chrpath                                             29.6kB @ 245.7kB/s  0.1s
gcc                                                 27.0kB @ 206.7kB/s  0.1s
libpq                                                2.5MB @  15.3MB/s  0.2s
libwebp                                             84.9kB @ 482.4kB/s  0.1s
libgd                                              224.4kB @   1.1MB/s  0.1s
gfortran                                            26.5kB @ 127.6kB/s  0.1s
glib-tools                                         112.0kB @ 385.9kB/s  0.1s
bc                                                 102.7kB @ 342.5kB/s  0.3s
libsystemd0                                        404.3kB @   1.1MB/s  0.1s
gdk-pixbuf                                         572.0kB @   1.5MB/s  0.1s
glib                                               488.2kB @   1.2MB/s  0.2s
gnuplot                                              1.2MB @   2.1MB/s  0.3s
qt-main                                             61.1MB @  60.4MB/s  1.1s

Downloading and Extracting Packages

Preparing transaction: done
Verifying transaction: done
Executing transaction: done

To activate this environment, use

     $ mamba activate siQ-ChIP_env

To deactivate an active environment, use

     $ mamba deactivate
```
</details>
<br />
