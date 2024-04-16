

<details>
<summary><i>Printed: Installation of phantompeakqualtools environment, etc. (MacBook)</i></summary>

```txt
❯ mamba search phantompeakqualtools

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

        mamba (0.25.0) supported by @QuantStack

        GitHub:  https://github.com/mamba-org/mamba
        Twitter: https://twitter.com/QuantStack

█████████████████████████████████████████████████████████████

Loading channels: done
# Name                       Version           Build  Channel
phantompeakqualtools             1.2               0  bioconda
phantompeakqualtools             1.2               1  bioconda
phantompeakqualtools           1.2.1               0  bioconda
phantompeakqualtools         1.2.1.1               0  bioconda
phantompeakqualtools           1.2.2               0  bioconda
phantompeakqualtools           1.2.2      hdfd78af_1  bioconda

    ~/repos/2023_tutorial_ChIP-seq    main !1 ?10 ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── 1m 57s   base   15:33:15 
❯ mamba create -n ppqt_env -c bioconda phantompeakqualtools=1.2.2

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

        mamba (0.25.0) supported by @QuantStack

        GitHub:  https://github.com/mamba-org/mamba
        Twitter: https://twitter.com/QuantStack

█████████████████████████████████████████████████████████████


Looking for: ['phantompeakqualtools=1.2.2']

bioconda/osx-64                                             Using cache
bioconda/noarch                                             Using cache
conda-forge/osx-64                                          Using cache
conda-forge/noarch                                          Using cache
r/osx-64                                                    Using cache
r/noarch                                                    Using cache
pkgs/main/osx-64                                              No change
pkgs/r/osx-64                                                 No change
pkgs/main/noarch                                              No change
pkgs/r/noarch                                                 No change
Transaction

  Prefix: /Users/kalavattam/miniconda3/envs/ppqt_env

  Updating specs:

   - phantompeakqualtools=1.2.2


  Package                               Version  Build                 Channel                  Size
──────────────────────────────────────────────────────────────────────────────────────────────────────
  Install:
──────────────────────────────────────────────────────────────────────────────────────────────────────

  + _r-mutex                              1.0.1  anacondar_1           conda-forge/noarch     Cached
  + argcomplete                           3.2.3  pyhd8ed1ab_0          conda-forge/noarch     Cached
  + bioconductor-biocgenerics            0.48.1  r43hdfd78af_2         bioconda/noarch        Cached
  + bioconductor-biocparallel            1.36.0  r43hc0ef7c4_1         bioconda/osx-64        Cached
  + bioconductor-biostrings              2.70.1  r43h4c50009_1         bioconda/osx-64        Cached
  + bioconductor-data-packages         20231203  hdfd78af_0            bioconda/noarch        Cached
  + bioconductor-genomeinfodb            1.38.1  r43hdfd78af_1         bioconda/noarch        Cached
  + bioconductor-genomeinfodbdata        1.2.11  r43hdfd78af_1         bioconda/noarch        Cached
  + bioconductor-genomicranges           1.54.1  r43h4c50009_1         bioconda/osx-64        Cached
  + bioconductor-iranges                 2.36.0  r43h4c50009_1         bioconda/osx-64        Cached
  + bioconductor-rhtslib                  2.4.0  r43h4c50009_1         bioconda/osx-64        Cached
  + bioconductor-rsamtools               2.18.0  r43hc0ef7c4_1         bioconda/osx-64        Cached
  + bioconductor-s4vectors               0.40.2  r43h4c50009_1         bioconda/osx-64        Cached
  + bioconductor-xvector                 0.42.0  r43h4c50009_1         bioconda/osx-64        Cached
  + bioconductor-zlibbioc                1.48.0  r43h4c50009_1         bioconda/osx-64        Cached
  + boost                                1.84.0  hdce95a9_2            conda-forge/osx-64       16kB
  + bwidget                              1.9.14  h694c41f_1            conda-forge/osx-64     Cached
  + bzip2                                 1.0.8  h10d778d_5            conda-forge/osx-64     Cached
  + c-ares                               1.28.1  h10d778d_0            conda-forge/osx-64     Cached
  + ca-certificates                    2024.2.2  h8857fd0_0            conda-forge/osx-64     Cached
  + cairo                                1.18.0  h99e66fa_0            conda-forge/osx-64     Cached
  + cctools_osx-64                          986  h58a35ae_0            conda-forge/osx-64     Cached
  + clang                                18.1.2  default_h7151d67_1    conda-forge/osx-64     Cached
  + clang-18                             18.1.2  default_h7151d67_1    conda-forge/osx-64     Cached
  + clang_impl_osx-64                    18.1.2  h73f7f27_11           conda-forge/osx-64     Cached
  + clang_osx-64                         18.1.2  hb91bd55_11           conda-forge/osx-64     Cached
  + clangxx                              18.1.2  default_h7151d67_1    conda-forge/osx-64     Cached
  + clangxx_impl_osx-64                  18.1.2  hb14bd1d_11           conda-forge/osx-64     Cached
  + clangxx_osx-64                       18.1.2  hb91bd55_11           conda-forge/osx-64     Cached
  + compiler-rt                          18.1.2  ha38d28d_0            conda-forge/osx-64     Cached
  + compiler-rt_osx-64                   18.1.2  ha38d28d_0            conda-forge/noarch     Cached
  + curl                                  8.7.1  h726d00d_0            conda-forge/osx-64     Cached
  + expat                                 2.6.2  h73e2aa4_0            conda-forge/osx-64     Cached
  + font-ttf-dejavu-sans-mono              2.37  hab24e00_0            conda-forge/noarch     Cached
  + font-ttf-inconsolata                  3.000  h77eed37_0            conda-forge/noarch     Cached
  + font-ttf-source-code-pro              2.038  h77eed37_0            conda-forge/noarch     Cached
  + font-ttf-ubuntu                        0.83  h77eed37_1            conda-forge/noarch     Cached
  + fontconfig                           2.14.2  h5bb23bf_0            conda-forge/osx-64     Cached
  + fonts-conda-ecosystem                     1  0                     conda-forge/noarch     Cached
  + fonts-conda-forge                         1  0                     conda-forge/noarch     Cached
  + freetype                             2.12.1  h60636b9_2            conda-forge/osx-64     Cached
  + fribidi                              1.0.10  hbcb3906_0            conda-forge/osx-64     Cached
  + gawk                                  5.3.0  h2c496e9_0            conda-forge/osx-64     Cached
  + gettext                              0.21.1  h8a4c099_0            conda-forge/osx-64     Cached
  + gfortran_impl_osx-64                 12.3.0  hc328e78_3            conda-forge/osx-64     Cached
  + gfortran_osx-64                      12.3.0  h18f7dce_1            conda-forge/osx-64     Cached
  + gmp                                   6.3.0  h73e2aa4_1            conda-forge/osx-64     Cached
  + graphite2                            1.3.13  h73e2aa4_1003         conda-forge/osx-64     Cached
  + harfbuzz                              8.3.0  hf45c392_0            conda-forge/osx-64     Cached
  + htslib                               1.19.1  h365c357_2            bioconda/osx-64        Cached
  + icu                                    73.2  hf5e326d_0            conda-forge/osx-64     Cached
  + isl                                    0.26  imath32_h2e86a7b_101  conda-forge/osx-64     Cached
  + jq                                      1.5  0                     bioconda/osx-64        Cached
  + krb5                                 1.21.2  hb884880_0            conda-forge/osx-64     Cached
  + ld64_osx-64                             711  had5d0d3_0            conda-forge/osx-64     Cached
  + lerc                                  4.0.0  hb486fe8_0            conda-forge/osx-64     Cached
  + libblas                               3.9.0  21_osx64_openblas     conda-forge/osx-64     Cached
  + libboost                             1.84.0  h6ebd1c4_2            conda-forge/osx-64        2MB
  + libboost-devel                       1.84.0  hf2b3138_2            conda-forge/osx-64       39kB
  + libboost-headers                     1.84.0  h694c41f_2            conda-forge/osx-64       14MB
  + libboost-python                      1.84.0  py312h77b368e_2       conda-forge/osx-64      107kB
  + libboost-python-devel                1.84.0  py312hdce95a9_2       conda-forge/osx-64       20kB
  + libcblas                              3.9.0  21_osx64_openblas     conda-forge/osx-64     Cached
  + libclang-cpp18.1                     18.1.2  default_h7151d67_1    conda-forge/osx-64     Cached
  + libcurl                               8.7.1  h726d00d_0            conda-forge/osx-64     Cached
  + libcxx                               16.0.6  hd57cbcb_0            conda-forge/osx-64     Cached
  + libdeflate                             1.18  hac1461d_0            conda-forge/osx-64     Cached
  + libedit                        3.1.20191231  h0678c8f_2            conda-forge/osx-64     Cached
  + libev                                  4.33  h10d778d_2            conda-forge/osx-64     Cached
  + libexpat                              2.6.2  h73e2aa4_0            conda-forge/osx-64     Cached
  + libffi                                3.4.2  h0d85af4_5            conda-forge/osx-64     Cached
  + libgcc                                4.8.5  1                     conda-forge/osx-64     Cached
  + libgfortran                           5.0.0  13_2_0_h97931a8_3     conda-forge/osx-64     Cached
  + libgfortran-devel_osx-64             12.3.0  h0b6f5ec_3            conda-forge/noarch     Cached
  + libgfortran5                         13.2.0  h2873a65_3            conda-forge/osx-64     Cached
  + libglib                              2.78.1  h6d9ecee_0            conda-forge/osx-64     Cached
  + libiconv                               1.17  hd75f5a5_2            conda-forge/osx-64     Cached
  + libjpeg-turbo                       2.1.5.1  h0dc2134_1            conda-forge/osx-64     Cached
  + liblapack                             3.9.0  21_osx64_openblas     conda-forge/osx-64     Cached
  + libllvm18                            18.1.2  hbcf5fad_0            conda-forge/osx-64     Cached
  + libnghttp2                           1.58.0  h64cf6d3_1            conda-forge/osx-64     Cached
  + libopenblas                          0.3.26  openmp_hfef2a42_0     conda-forge/osx-64     Cached
  + libpng                               1.6.43  h92b6c6a_0            conda-forge/osx-64     Cached
  + libsqlite                            3.45.2  h92b6c6a_0            conda-forge/osx-64     Cached
  + libssh2                              1.11.0  hd019ec5_0            conda-forge/osx-64     Cached
  + libtiff                               4.6.0  hf955e92_0            conda-forge/osx-64     Cached
  + libwebp-base                          1.3.2  h0dc2134_0            conda-forge/osx-64     Cached
  + libxml2                              2.12.6  hc0ae0f7_1            conda-forge/osx-64     Cached
  + libzlib                              1.2.13  h8a1eda9_5            conda-forge/osx-64     Cached
  + llvm-openmp                          18.1.2  hb6ac08f_0            conda-forge/osx-64     Cached
  + llvm-tools                           18.1.2  hbcf5fad_0            conda-forge/osx-64     Cached
  + make                                    4.3  h22f3db7_1            conda-forge/osx-64     Cached
  + mpc                                   1.3.1  h81bd1dd_0            conda-forge/osx-64     Cached
  + mpfr                                  4.2.1  h0c69b56_0            conda-forge/osx-64     Cached
  + ncurses                        6.4.20240210  h73e2aa4_0            conda-forge/osx-64     Cached
  + numpy                                1.26.4  py312he3a82b2_0       conda-forge/osx-64        7MB
  + openssl                               3.2.1  hd75f5a5_1            conda-forge/osx-64     Cached
  + pango                               1.50.14  h19c1c8a_2            conda-forge/osx-64     Cached
  + pcre2                                 10.40  h1c4e4bc_0            conda-forge/osx-64     Cached
  + phantompeakqualtools                  1.2.2  hdfd78af_1            bioconda/noarch        Cached
  + pip                                    24.0  pyhd8ed1ab_0          conda-forge/noarch     Cached
  + pixman                               0.43.4  h73e2aa4_0            conda-forge/osx-64     Cached
  + python                               3.12.2  h9f0c242_0_cpython    conda-forge/osx-64     Cached
  + python_abi                             3.12  4_cp312               conda-forge/osx-64     Cached
  + pyyaml                                6.0.1  py312h104f124_1       conda-forge/osx-64     Cached
  + r-base                                4.3.1  h61172b1_5            conda-forge/osx-64     Cached
  + r-bh                               1.84.0_0  r43hc72bb7e_0         conda-forge/noarch     Cached
  + r-bitops                              1.0_7  r43h6dc245f_2         conda-forge/osx-64     Cached
  + r-catools                            1.18.2  r43hac7d2d5_2         conda-forge/osx-64     Cached
  + r-codetools                          0.2_20  r43hc72bb7e_0         conda-forge/noarch     Cached
  + r-cpp11                               0.4.7  r43hc72bb7e_0         conda-forge/noarch     Cached
  + r-crayon                              1.5.2  r43hc72bb7e_2         conda-forge/noarch     Cached
  + r-formatr                              1.14  r43hc72bb7e_1         conda-forge/noarch     Cached
  + r-futile.logger                       1.4.3  r43hc72bb7e_1005      conda-forge/noarch     Cached
  + r-futile.options                      1.0.1  r43hc72bb7e_1004      conda-forge/noarch     Cached
  + r-lambda.r                            1.2.4  r43hc72bb7e_3         conda-forge/noarch     Cached
  + r-rcpp                               1.0.12  r43h29979af_0         conda-forge/osx-64     Cached
  + r-rcurl                           1.98_1.14  r43hbd64cb6_0         conda-forge/osx-64     Cached
  + r-snow                                0.4_4  r43hc72bb7e_2         conda-forge/noarch     Cached
  + r-snowfall                         1.84_6.3  r43hc72bb7e_0         conda-forge/noarch      261kB
  + r-spp                                1.16.0  r43h71b3455_9         bioconda/osx-64        Cached
  + readline                                8.2  h9e318b2_1            conda-forge/osx-64     Cached
  + samtools                             1.19.2  hd510865_1            bioconda/osx-64        Cached
  + setuptools                           69.2.0  pyhd8ed1ab_0          conda-forge/noarch     Cached
  + sigtool                               0.1.3  h88f4db0_0            conda-forge/osx-64     Cached
  + tapi                              1100.0.11  h9ce4665_0            conda-forge/osx-64     Cached
  + tk                                   8.6.13  h1abcd95_1            conda-forge/osx-64     Cached
  + tktable                                2.10  ha166976_5            conda-forge/osx-64     Cached
  + toml                                 0.10.2  pyhd8ed1ab_0          conda-forge/noarch     Cached
  + tomlkit                              0.12.4  pyha770c72_0          conda-forge/noarch     Cached
  + tzdata                                2024a  h0c530f3_0            conda-forge/noarch     Cached
  + wheel                                0.43.0  pyhd8ed1ab_1          conda-forge/noarch     Cached
  + xmltodict                            0.13.0  pyhd8ed1ab_0          conda-forge/noarch     Cached
  + xz                                    5.2.6  h775f41a_0            conda-forge/osx-64     Cached
  + yaml                                  0.2.5  h0d85af4_2            conda-forge/osx-64     Cached
  + yq                                    3.2.3  pyhd8ed1ab_0          conda-forge/noarch     Cached
  + zlib                                 1.2.13  h8a1eda9_5            conda-forge/osx-64     Cached
  + zstd                                  1.5.5  h829000d_0            conda-forge/osx-64     Cached

  Summary:

  Install: 138 packages

  Total download: 23MB

──────────────────────────────────────────────────────────────────────────────────────────────────────

Confirm changes: [Y/n] Y
libboost-devel                                      39.4kB @  79.3kB/s  0.5s
boost                                               16.2kB @  25.2kB/s  0.1s
libboost-python                                    107.2kB @ 132.7kB/s  0.2s
libboost-python-devel                               19.6kB @  20.5kB/s  0.1s
libboost                                             2.1MB @   2.1MB/s  1.0s
r-snowfall                                         261.5kB @ 259.6kB/s  1.0s
numpy                                                7.0MB @   4.5MB/s  1.6s
libboost-headers                                    13.8MB @   6.7MB/s  2.1s
Preparing transaction: done
Verifying transaction: done
Executing transaction: done

To activate this environment, use

     $ mamba activate ppqt_env

To deactivate an active environment, use

     $ mamba deactivate


❯ zaliases


❯ alias R-test="(
    conda activate ppqt_env
    open -na /Applications/RStudio.app
)"


❯ R-test


❯ conda activate ppqt_env


❯ run_spp.R --help
Loading required package: Rcpp
Warning message:
package ‘Rcpp’ was built under R version 4.3.2
Error in parse.arguments(args) : Illegal argument --help
Execution halted


❯ run_spp.R -h
Loading required package: Rcpp
Warning message:
package ‘Rcpp’ was built under R version 4.3.2
Error in parse.arguments(args) : Illegal argument -h
Execution halted


❯ run_spp.R
Loading required package: Rcpp
Warning message:
package ‘Rcpp’ was built under R version 4.3.2
Usage: Rscript run_spp.R <options>
MANDATORY ARGUMENTS
-c=<ChIP_alignFile>, full path and name (or URL) of tagAlign/BAM file (can be gzipped) (FILE EXTENSION MUST BE tagAlign.gz, tagAlign, bam or bam.gz)
MANDATORY ARGUMENTS FOR PEAK CALLING
-i=<Input_alignFile>, full path and name (or URL) of tagAlign/BAM file (can be gzipped) (FILE EXTENSION MUST BE tagAlign.gz, tagAlign, bam or bam.gz)
OPTIONAL ARGUMENTS
-s=<min>:<step>:<max> , strand shifts at which cross-correlation is evaluated, default=-500:5:1500
-speak=<strPeak>, user-defined cross-correlation peak strandshift
-x=<min>:<max>, strand shifts to exclude (This is mainly to avoid region around phantom peak) default=10:(readlen+10)
-p=<nodes> , number of parallel processing nodes, default=0
-fdr=<falseDisoveryRate> , false discovery rate threshold for peak calling
-npeak=<numPeaks>, threshold on number of peaks to call
-tmpdir=<tempdir> , Temporary directory (if not specified R function tempdir() is used)
-filtchr=<chrnamePattern> , Pattern to use to remove tags that map to specific chromosomes e.g. _ will remove all tags that map to chromosomes with _ in their name
OUTPUT ARGUMENTS
-odir=<outputDirectory> name of output directory (If not set same as ChIP file directory is used)
-savn=<narrowpeakfilename> OR -savn NarrowPeak file name (fixed width peaks)
-savr=<regionpeakfilename> OR -savr RegionPeak file name (variable width peaks with regions of enrichment)
-savd=<rdatafile> OR -savd, save Rdata file
-savp=<plotdatafile> OR -savp, save cross-correlation plot
-out=<resultfile>, append peakshift/phantomPeak results to a file
     format:Filename<tab>numReads<tab>estFragLen<tab>corr_estFragLen<tab>PhantomPeak<tab>corr_phantomPeak<tab>argmin_corr<tab>min_corr<tab>Normalized SCC (NSC)<tab>Relative SCC (RSC)<tab>QualityTag)
-rf, if plot or rdata or narrowPeak file exists replace it. If not used then the run is aborted if the plot or Rdata or narrowPeak file exists
-clean, if present will remove the original chip and control files after reading them in. CAUTION: Use only if the script calling run_spp.R is creating temporary files
```
</details>
<br />
