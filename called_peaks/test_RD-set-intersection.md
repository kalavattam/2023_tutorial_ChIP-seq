
`#test_RD-set-intersection.md`
<br />

<details>
<summary><b><font size="+2"><i>Table of contents</i></font></b></summary>
<!-- MarkdownTOC -->

1. [Copy data from RD to KA](#copy-data-from-rd-to-ka)
    1. [Code](#code)
1. [Set up environment for running R script, `test_RD-set-intersection.R`](#set-up-environment-for-running-r-script-test_rd-set-intersectionr)
    1. [Configure the package management program, Mamba/Conda; then, install packages](#configure-the-package-management-program-mambaconda-then-install-packages)
        1. [Code](#code-1)
        1. [Printed](#printed)

<!-- /MarkdownTOC -->
</details>
<br />
<br />

<a id="copy-data-from-rd-to-ka"></a>
## Copy data from RD to KA
<a id="code"></a>
### Code
<details>
<summary><i>Code: Copy data from RD to KA</i></summary>

```bash
#!/bin/bash

grabnode  # 1, 20, 1, N


#  Initialize necessary variables =============================================
    p_RD="${HOME}/tsukiyamalab/Rachel"
p_RD_exp="misc_data/experiments/ChIPs/HDAC_HAT_Q/results"

    p_KA="${HOME}/tsukiyamalab/Kris"
p_KA_exp="2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses"

#  Using GSE95556_Sc.cerevisiae.feature.anno_Steinmetz_2013.gtf, Christine's
#+ BED file
dir_MACS2="called_peaks"  # ls -lhaFG "${bed}"

check_variables=true
if ${check_variables}; then
    echo """
          p_RD=${p_RD}
      p_RD_exp=${p_RD_exp}

          p_KA=${p_KA}
      p_KA_exp=${p_KA_exp}

    dir_MACS2=${dir_MACS2}
    """
fi


#  Run operations ===========================================================
#  Go to RD's work directory
cd "${p_RD}/${p_RD_exp}" || echo "cd'ing failed; check on this"

if [[ ! -d "${p_KA}/${p_KA_exp}/${dir_MACS2}" ]]; then
    cp -r "${dir_MACS2}" "${p_KA}/${p_KA_exp}"
fi

cd "${p_KA}/${p_KA_exp}" || echo "cd'ing failed; check on this"

ls -1 "${dir_MACS2}"
```
</details>
<br />
<br />

<a id="set-up-environment-for-running-r-script-test_rd-set-intersectionr"></a>
## Set up environment for running R script, [`test_RD-set-intersection.R`](./test_RD-set-intersection.R)
<a id="configure-the-package-management-program-mambaconda-then-install-packages"></a>
### Configure the package management program, Mamba/Conda; then, install packages
<a id="code-1"></a>
#### Code
<details>
<summary><i>Code: Configure the package management program, Mamba/Conda; then, install packages</i></summary>

```bash
#!/bin/bash

check_settings=false  # Flag variable to check conda/mamba settings
if ${check_settings}; then
    #  Display current channel configuration
    conda config --show channels

    #  Describe channel_priority setting
    conda config --describe channel_priority
fi

change_channel_priority=true  # Flag to change channel_priority
if ${change_channel_priority}; then
    #  We can do it manually:
    # open ~/.condarc  # Change channel_priority from to disabled; save file.
    
    #  Or we can do it from the command line (which is easier/more
    #+ reproducible):
    conda config --set channel_priority disabled

    #  When channel_priority is set to "disabled", this helps Mamba/Conda
    #+ resolve issues regarding package versions between the "conda-forge" and
    #+ "bioconda" channels.
    #+ 
    #+ This needs to be set only once. The settings are saved thereafter.
    #+ 
    #+ More information on conda channels:
    #+ conda.io/projects/conda/en/latest/user-guide/tasks/manage-channels.html
fi

set_channel_priority_order=false  # Flag to configure channel order
if ${set_channel_priority_order}; then    
    conda config --add channels defaults
    conda config --add channels r
    conda config --add channels bioconda
    conda config --add channels conda-forge
fi

#  After running the above, the channel-order priority will be, from highest to
#+ lowest, (1) "conda-forge", (2) "bioconda", (3) "r", and then (4) "defaults".
#+ 
#+ It is important that the channels are in this order even though we have set
#+ channel_priority to "disabled"; the channel-order priority directs
#+ Conda/Mamba to check and pull certain packages from some before checking
#+ others.
#+ 
#+ This needs to be set only once. The settings are saved thereafter.

#  Once the channel_priority is disabled and the channel order is established,
#+ let's make an environment for running test_RD-set-intersection.R.

#  First, let's make sure we're in Conda/Mamba's default environment, "base":
if [[ "${CONDA_DEFAULT_ENV}" != "base" ]]; then
    conda deactivate
fi

#  Now, let's make the environment:
create_environment=true      # Flag to install environment
env_name="test_overlap_env"  # Change to preferred name for new environment
if ${create_environment}; then
    CONDA_SUBDIR=osx-64 \
    mamba create \
        -n "${env_name}" \
        -c conda-forge \
        -c bioconda \
            bioconductor-genomicranges \
            bioconductor-iranges \
            r-tidyverse \
            r-venneuler
fi

#  Note that bioconda packages are not compiled for the ARM architecture used
#+ by Apple machines from 2020 onward, so we set the Conda/Mamba environmental
#+ variable CONDA_SUBDIR to osx-64, which essentially tells Conda/Mamba to
#+ install versions of the packages compiled for the previous Intel 64-bit
#+ architecture used by Apple machines prior to 2020. (Modern Macs do this
#+ through very efficient emulation with the program Rosetta. It's all
#+ happening "under the hood", so you do not need to worry about this.)
#+ 
#+ I think you (Rachel) are currently on an older Apple machine, so it's not
#+ strictly necessary to declare CONDA_SUBDIR=osx-64 ; however, it's still good
#+ practice to include CONDA_SUBDIR=osx-64 in that case too.
#+ 
#+ That said, do *not* include CONDA_SUBDIR=osx-64 when using Conda/Mamba to
#+ install package environments on the FHCC cluster, which uses a Linux
#+ operating system with an Intel 64-bit architecture.

#  Now, let's create a shell (Bash, Zsh, etc.) alias to open RStudio in our new
#+ environment, ensuring that RStudio's interactive R environment can access
#+ the programs we just installed:
alias R-overlap="(
    conda activate ${env_name}
    open -na /Applications/RStudio.app
)"

#  Remember that aliases defined in a terminal session are not persistent
#+ across sessions. To make this alias permanent, you should add it to your
#+ ~/.bashrc, ~/.bash_profile, ~/.zshrc, etc. file (depending on your shell,
#+ operating system, and general setup for these things).
#+ 
#+ When adding the alias to ~/.bashrc, ~/.bash_profile, ~/.zshrc, etc.,
#+ explicitly specify the environment name, i.e., do not use a variable as
#+ above, e.g.,
#+     alias R-overlap="(
#+         mamba activate test_overlap_env
#+         open -na /Applications/RStudio.app
#+     )"
#+ 
#+ In order for this change to take effect, you need to either restart the
#+ terminal or source ~/.bashrc, ~/.bash_profile, ~/.zshrc on the command line,
#+ e.g.,
#+     source ~/.bashrc
#+     source ~/.bash_profile
#+     source ~/.zshrc
#+ 
#+ Also, the open -na /Applications/RStudio.app command is specific to macOS
#+ and will work as expected on that platform. It opens a new instance of
#+ RStudio. You'll need to use a different command if on Linux or Windows, but
#+ no need to worry about that now.

#  Now, let's open RStudio with an interactive R environment connected to
#+ "test_overlap_env" (or whatever you initialized variable env_name as). You
#+ can open RStudio by typing R-overlap in the command line:
open_RStudio_with_env=true  # Flag to do the above
if ${open_RStudio_with_env}; then
    R-overlap
fi
```
</details>
<br />

<a id="printed"></a>
#### Printed
*It's a good idea to keep a record of the installation.*

<details>
<summary><i>Printed: Configure the package management program, Mamba/Conda; then, install packages</i></summary>

```txt
❯ if ${create_environment}; then
then>     env_name="test_overlap_env"
then>     CONDA_SUBDIR=osx-64 \
then>     mamba create \
then>         -n "${env_name}" \
then>         -c conda-forge \
then>         -c bioconda \
then>             bioconductor-genomicranges \
then>             bioconductor-iranges \
then>             r-tidyverse \
then>             r-venneuler
then> fi

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

        mamba (1.4.1) supported by @QuantStack

        GitHub:  https://github.com/mamba-org/mamba
        Twitter: https://twitter.com/QuantStack

█████████████████████████████████████████████████████████████


Looking for: ['bioconductor-genomicranges', 'bioconductor-iranges', 'r-tidyverse', 'r-venneuler']

r/osx-64                                                      No change
r/noarch                                                      No change
pkgs/r/osx-64                                                 No change
pkgs/r/noarch                                                 No change
bioconda/osx-64                                      4.3MB @   4.0MB/s  1.2s
bioconda/noarch                                      4.7MB @   3.9MB/s  1.3s
pkgs/main/noarch                                              No change
pkgs/main/osx-64                                     6.0MB @   3.6MB/s  1.5s
conda-forge/noarch                                  14.6MB @   6.7MB/s  2.4s
conda-forge/osx-64                                  32.2MB @   6.8MB/s  5.3s
Transaction

  Prefix: /Users/kalavatt/mambaforge/envs/test_overlap_env

  Updating specs:

   - bioconductor-genomicranges
   - bioconductor-iranges
   - r-tidyverse
   - r-venneuler


  Package                               Version  Build               Channel                  Size
────────────────────────────────────────────────────────────────────────────────────────────────────
  Install:
────────────────────────────────────────────────────────────────────────────────────────────────────

  + _r-mutex                              1.0.1  anacondar_1         conda-forge/noarch     Cached
  + argcomplete                           3.1.6  pyhd8ed1ab_0        conda-forge/noarch       40kB
  + bioconductor-biocgenerics            0.46.0  r43hdfd78af_0       bioconda/noarch        Cached
  + bioconductor-data-packages         20230718  hdfd78af_1          bioconda/noarch        Cached
  + bioconductor-genomeinfodb            1.36.1  r43hdfd78af_0       bioconda/noarch        Cached
  + bioconductor-genomeinfodbdata        1.2.10  r43hdfd78af_0       bioconda/noarch        Cached
  + bioconductor-genomicranges           1.52.0  r43h4c50009_0       bioconda/osx-64        Cached
  + bioconductor-iranges                 2.34.1  r43h4c50009_0       bioconda/osx-64        Cached
  + bioconductor-s4vectors               0.38.1  r43h4c50009_0       bioconda/osx-64        Cached
  + bioconductor-xvector                 0.40.0  r43h4c50009_0       bioconda/osx-64        Cached
  + bioconductor-zlibbioc                1.46.0  r43h4c50009_0       bioconda/osx-64        Cached
  + bwidget                              1.9.14  h694c41f_1          conda-forge/osx-64     Cached
  + bzip2                                 1.0.8  h10d778d_5          conda-forge/osx-64     Cached
  + c-ares                               1.21.0  h10d778d_0          conda-forge/osx-64     Cached
  + ca-certificates                  2023.08.22  hecd8cb5_0          pkgs/main/osx-64       Cached
  + cairo                                1.18.0  h99e66fa_0          conda-forge/osx-64     Cached
  + cctools_osx-64                      973.0.1  ha1c5b94_15         conda-forge/osx-64     Cached
  + clang                                16.0.6  hc177806_2          conda-forge/osx-64       21kB
  + clang-16                             16.0.6  default_h762fdd7_2  conda-forge/osx-64      695kB
  + clang_impl_osx-64                    16.0.6  h8787910_6          conda-forge/osx-64     Cached
  + clang_osx-64                         16.0.6  hb91bd55_6          conda-forge/osx-64     Cached
  + clangxx                              16.0.6  default_h762fdd7_2  conda-forge/osx-64       21kB
  + clangxx_impl_osx-64                  16.0.6  h1b7723c_6          conda-forge/osx-64     Cached
  + clangxx_osx-64                       16.0.6  h2fff7d5_6          conda-forge/osx-64     Cached
  + compiler-rt                          16.0.6  he1888fc_1          conda-forge/osx-64     Cached
  + compiler-rt_osx-64                   16.0.6  he1888fc_1          conda-forge/noarch     Cached
  + curl                                  8.4.0  h726d00d_0          conda-forge/osx-64     Cached
  + expat                                 2.5.0  hf0c8a7f_1          conda-forge/osx-64     Cached
  + font-ttf-dejavu-sans-mono              2.37  hab24e00_0          conda-forge/noarch     Cached
  + font-ttf-inconsolata                  3.000  h77eed37_0          conda-forge/noarch     Cached
  + font-ttf-source-code-pro              2.038  h77eed37_0          conda-forge/noarch     Cached
  + font-ttf-ubuntu                        0.83  hab24e00_0          conda-forge/noarch     Cached
  + fontconfig                           2.14.2  h5bb23bf_0          conda-forge/osx-64     Cached
  + fonts-conda-ecosystem                     1  0                   conda-forge/noarch     Cached
  + fonts-conda-forge                         1  0                   conda-forge/noarch     Cached
  + freetype                             2.12.1  h60636b9_2          conda-forge/osx-64     Cached
  + fribidi                              1.0.10  hbcb3906_0          conda-forge/osx-64     Cached
  + gettext                              0.21.1  h8a4c099_0          conda-forge/osx-64     Cached
  + gfortran_impl_osx-64                 12.3.0  h54fd467_1          conda-forge/osx-64     Cached
  + gfortran_osx-64                      12.3.0  h18f7dce_1          conda-forge/osx-64     Cached
  + gmp                                   6.3.0  h93d8f39_0          conda-forge/osx-64      521kB
  + graphite2                            1.3.14  he9d5cce_1          pkgs/main/osx-64       Cached
  + harfbuzz                              8.3.0  hf45c392_0          conda-forge/osx-64        1MB
  + icu                                    73.2  hf5e326d_0          conda-forge/osx-64     Cached
  + isl                                    0.25  hb486fe8_0          conda-forge/osx-64     Cached
  + jq                                      1.6  hc929b4f_1000       conda-forge/osx-64     Cached
  + krb5                                 1.21.2  hb884880_0          conda-forge/osx-64     Cached
  + ld64_osx-64                             609  ha20a434_15         conda-forge/osx-64     Cached
  + lerc                                  4.0.0  hb486fe8_0          conda-forge/osx-64     Cached
  + libblas                               3.9.0  19_osx64_openblas   conda-forge/osx-64     Cached
  + libclang-cpp16                       16.0.6  default_h762fdd7_2  conda-forge/osx-64       13MB
  + libcurl                               8.4.0  h726d00d_0          conda-forge/osx-64     Cached
  + libcxx                               16.0.6  hd57cbcb_0          conda-forge/osx-64     Cached
  + libdeflate                             1.19  ha4e1b8e_0          conda-forge/osx-64     Cached
  + libedit                        3.1.20221030  h6c40b1e_0          pkgs/main/osx-64       Cached
  + libev                                  4.33  haf1e3a3_1          conda-forge/osx-64     Cached
  + libexpat                              2.5.0  hf0c8a7f_1          conda-forge/osx-64     Cached
  + libffi                                3.4.4  hecd8cb5_0          pkgs/main/osx-64       Cached
  + libgfortran                           5.0.0  13_2_0_h97931a8_1   conda-forge/osx-64     Cached
  + libgfortran-devel_osx-64             12.3.0  h0b6f5ec_1          conda-forge/noarch     Cached
  + libgfortran5                         13.2.0  h2873a65_1          conda-forge/osx-64     Cached
  + libglib                              2.78.1  h6d9ecee_0          conda-forge/osx-64     Cached
  + libiconv                               1.17  hac89ed1_0          conda-forge/osx-64     Cached
  + libjpeg-turbo                         3.0.0  h0dc2134_1          conda-forge/osx-64     Cached
  + liblapack                             3.9.0  19_osx64_openblas   conda-forge/osx-64     Cached
  + libllvm16                            16.0.6  he4b1e75_2          conda-forge/osx-64     Cached
  + libnghttp2                           1.58.0  h64cf6d3_0          conda-forge/osx-64     Cached
  + libopenblas                          0.3.24  openmp_h48a4ad5_0   conda-forge/osx-64     Cached
  + libpng                               1.6.39  ha978bb4_0          conda-forge/osx-64     Cached
  + libsqlite                            3.44.0  h92b6c6a_0          conda-forge/osx-64     Cached
  + libssh2                              1.11.0  hd019ec5_0          conda-forge/osx-64     Cached
  + libtiff                               4.6.0  h684deea_2          conda-forge/osx-64     Cached
  + libwebp-base                          1.3.2  h0dc2134_0          conda-forge/osx-64     Cached
  + libxml2                              2.11.5  h3346baf_1          conda-forge/osx-64     Cached
  + libzlib                              1.2.13  h8a1eda9_5          conda-forge/osx-64     Cached
  + llvm-openmp                          17.0.5  hb6ac08f_0          conda-forge/osx-64      300kB
  + llvm-tools                           16.0.6  he4b1e75_2          conda-forge/osx-64     Cached
  + make                                    4.3  h22f3db7_1          conda-forge/osx-64     Cached
  + mpc                                   1.3.1  h81bd1dd_0          conda-forge/osx-64     Cached
  + mpfr                                  4.2.1  h0c69b56_0          conda-forge/osx-64     Cached
  + ncurses                                 6.4  h93d8f39_2          conda-forge/osx-64     Cached
  + oniguruma                             6.9.9  h10d778d_0          conda-forge/osx-64     Cached
  + openjdk                              21.0.1  hf4d7fad_0          conda-forge/osx-64     Cached
  + openssl                               3.1.4  hd75f5a5_0          conda-forge/osx-64     Cached
  + pandoc                                3.1.3  h9d075a6_0          conda-forge/osx-64     Cached
  + pango                               1.50.14  h19c1c8a_2          conda-forge/osx-64     Cached
  + pcre2                                 10.40  h1c4e4bc_0          conda-forge/osx-64     Cached
  + pip                                  23.3.1  pyhd8ed1ab_0        conda-forge/noarch     Cached
  + pixman                               0.42.2  he965462_0          conda-forge/osx-64     Cached
  + python                               3.12.0  h30d4d87_0_cpython  conda-forge/osx-64     Cached
  + python_abi                             3.12  4_cp312             conda-forge/osx-64     Cached
  + pyyaml                                6.0.1  py312h104f124_1     conda-forge/osx-64     Cached
  + r-askpass                             1.2.0  r43h6dc245f_0       conda-forge/osx-64     Cached
  + r-assertthat                          0.2.1  r43hc72bb7e_4       conda-forge/noarch     Cached
  + r-backports                           1.4.1  r43h6dc245f_2       conda-forge/osx-64     Cached
  + r-base                                4.3.2  hf996cbc_0          conda-forge/osx-64       25MB
  + r-base64enc                           0.1_3  r43h6dc245f_1006    conda-forge/osx-64     Cached
  + r-bit                                 4.0.5  r43h6dc245f_1       conda-forge/osx-64     Cached
  + r-bit64                               4.0.5  r43h6dc245f_2       conda-forge/osx-64     Cached
  + r-bitops                              1.0_7  r43h6dc245f_2       conda-forge/osx-64     Cached
  + r-blob                                1.2.4  r43hc72bb7e_1       conda-forge/noarch     Cached
  + r-broom                               1.0.5  r43hc72bb7e_1       conda-forge/noarch     Cached
  + r-bslib                               0.5.1  r43hc72bb7e_0       conda-forge/noarch     Cached
  + r-cachem                              1.0.8  r43h6dc245f_1       conda-forge/osx-64     Cached
  + r-callr                               3.7.3  r43hc72bb7e_1       conda-forge/noarch     Cached
  + r-cellranger                          1.1.0  r43hc72bb7e_1006    conda-forge/noarch     Cached
  + r-cli                                 3.6.1  r43hac7d2d5_1       conda-forge/osx-64     Cached
  + r-clipr                               0.8.0  r43hc72bb7e_2       conda-forge/noarch     Cached
  + r-colorspace                          2.1_0  r43h6dc245f_1       conda-forge/osx-64     Cached
  + r-conflicted                          1.2.0  r43h785f33e_1       conda-forge/noarch     Cached
  + r-cpp11                               0.4.6  r43hc72bb7e_0       conda-forge/noarch     Cached
  + r-crayon                              1.5.2  r43hc72bb7e_2       conda-forge/noarch     Cached
  + r-curl                                5.1.0  r43h0100ac3_0       conda-forge/osx-64     Cached
  + r-data.table                         1.14.8  r43h7eccc33_2       conda-forge/osx-64     Cached
  + r-dbi                                 1.1.3  r43hc72bb7e_2       conda-forge/noarch     Cached
  + r-dbplyr                              2.4.0  r43hc72bb7e_0       conda-forge/noarch     Cached
  + r-digest                             0.6.33  r43hac7d2d5_0       conda-forge/osx-64     Cached
  + r-dplyr                               1.1.3  r43hac7d2d5_0       conda-forge/osx-64     Cached
  + r-dtplyr                              1.3.1  r43hc72bb7e_1       conda-forge/noarch     Cached
  + r-ellipsis                            0.3.2  r43h6dc245f_2       conda-forge/osx-64     Cached
  + r-evaluate                             0.23  r43hc72bb7e_0       conda-forge/noarch     Cached
  + r-fansi                               1.0.5  r43hb2c329c_0       conda-forge/osx-64     Cached
  + r-farver                              2.1.1  r43hac7d2d5_2       conda-forge/osx-64     Cached
  + r-fastmap                             1.1.1  r43hac7d2d5_1       conda-forge/osx-64     Cached
  + r-fontawesome                         0.5.2  r43hc72bb7e_0       conda-forge/noarch     Cached
  + r-forcats                             1.0.0  r43hc72bb7e_1       conda-forge/noarch     Cached
  + r-fs                                  1.6.3  r43hac7d2d5_0       conda-forge/osx-64     Cached
  + r-gargle                              1.5.2  r43h785f33e_0       conda-forge/noarch     Cached
  + r-generics                            0.1.3  r43hc72bb7e_2       conda-forge/noarch     Cached
  + r-ggplot2                             3.4.4  r43hc72bb7e_0       conda-forge/noarch     Cached
  + r-glue                                1.6.2  r43h6dc245f_2       conda-forge/osx-64     Cached
  + r-googledrive                         2.1.1  r43hc72bb7e_1       conda-forge/noarch     Cached
  + r-googlesheets4                       1.1.1  r43h785f33e_1       conda-forge/noarch     Cached
  + r-gtable                              0.3.4  r43hc72bb7e_0       conda-forge/noarch     Cached
  + r-haven                               2.5.3  r43hac7d2d5_0       conda-forge/osx-64     Cached
  + r-highr                                0.10  r43hc72bb7e_1       conda-forge/noarch     Cached
  + r-hms                                 1.1.3  r43hc72bb7e_1       conda-forge/noarch     Cached
  + r-htmltools                           0.5.7  r43h64b2c41_0       conda-forge/osx-64     Cached
  + r-httr                                1.4.7  r43hc72bb7e_0       conda-forge/noarch     Cached
  + r-ids                                 1.0.1  r43hc72bb7e_3       conda-forge/noarch     Cached
  + r-isoband                             0.2.7  r43hac7d2d5_2       conda-forge/osx-64     Cached
  + r-jquerylib                           0.1.4  r43hc72bb7e_2       conda-forge/noarch     Cached
  + r-jsonlite                            1.8.7  r43h6dc245f_0       conda-forge/osx-64     Cached
  + r-knitr                                1.45  r43hc72bb7e_0       conda-forge/noarch     Cached
  + r-labeling                            0.4.3  r43hc72bb7e_0       conda-forge/noarch     Cached
  + r-lattice                            0.22_5  r43hb2c329c_0       conda-forge/osx-64     Cached
  + r-lifecycle                           1.0.4  r43hc72bb7e_0       conda-forge/noarch     Cached
  + r-lubridate                           1.9.3  r43h6dc245f_0       conda-forge/osx-64     Cached
  + r-magrittr                            2.0.3  r43h6dc245f_2       conda-forge/osx-64     Cached
  + r-mass                               7.3_60  r43h6dc245f_1       conda-forge/osx-64     Cached
  + r-matrix                              1.6_3  r43h39411cf_0       conda-forge/osx-64        4MB
  + r-memoise                             2.0.1  r43hc72bb7e_2       conda-forge/noarch     Cached
  + r-mgcv                                1.9_0  r43h9c380c6_0       conda-forge/osx-64     Cached
  + r-mime                                 0.12  r43h6dc245f_2       conda-forge/osx-64     Cached
  + r-modelr                             0.1.11  r43hc72bb7e_1       conda-forge/noarch     Cached
  + r-munsell                             0.5.0  r43hc72bb7e_1006    conda-forge/noarch     Cached
  + r-nlme                              3.1_163  r43hfe07776_0       conda-forge/osx-64     Cached
  + r-openssl                             2.1.1  r43hc61a7e2_0       conda-forge/osx-64     Cached
  + r-pillar                              1.9.0  r43hc72bb7e_1       conda-forge/noarch     Cached
  + r-pkgconfig                           2.0.3  r43hc72bb7e_3       conda-forge/noarch     Cached
  + r-prettyunits                         1.2.0  r43hc72bb7e_0       conda-forge/noarch     Cached
  + r-processx                            3.8.2  r43h6dc245f_0       conda-forge/osx-64     Cached
  + r-progress                            1.2.2  r43hc72bb7e_4       conda-forge/noarch     Cached
  + r-ps                                  1.7.5  r43h6dc245f_1       conda-forge/osx-64     Cached
  + r-purrr                               1.0.2  r43h6dc245f_0       conda-forge/osx-64     Cached
  + r-r6                                  2.5.1  r43hc72bb7e_2       conda-forge/noarch     Cached
  + r-ragg                                1.2.6  r43h1914a8a_0       conda-forge/osx-64      397kB
  + r-rappdirs                            0.3.3  r43h6dc245f_2       conda-forge/osx-64     Cached
  + r-rcolorbrewer                        1.1_3  r43h785f33e_2       conda-forge/noarch     Cached
  + r-rcpp                               1.0.11  r43hac7d2d5_0       conda-forge/osx-64     Cached
  + r-rcurl                           1.98_1.13  r43hbd64cb6_0       conda-forge/osx-64     Cached
  + r-readr                               2.1.4  r43hac7d2d5_1       conda-forge/osx-64     Cached
  + r-readxl                              1.4.3  r43h88814b1_0       conda-forge/osx-64     Cached
  + r-rematch                             2.0.0  r43hc72bb7e_0       conda-forge/noarch     Cached
  + r-rematch2                            2.1.2  r43hc72bb7e_3       conda-forge/noarch     Cached
  + r-reprex                              2.0.2  r43hc72bb7e_2       conda-forge/noarch     Cached
  + r-rjava                               1.0_6  r43hc02fff5_7       conda-forge/osx-64     Cached
  + r-rlang                               1.1.2  r43h64b2c41_0       conda-forge/osx-64     Cached
  + r-rmarkdown                            2.25  r43hc72bb7e_0       conda-forge/noarch     Cached
  + r-rstudioapi                         0.15.0  r43hc72bb7e_0       conda-forge/noarch     Cached
  + r-rvest                               1.0.3  r43hc72bb7e_2       conda-forge/noarch     Cached
  + r-sass                                0.4.7  r43hac7d2d5_0       conda-forge/osx-64     Cached
  + r-scales                              1.2.1  r43hc72bb7e_2       conda-forge/noarch     Cached
  + r-selectr                             0.4_2  r43hc72bb7e_3       conda-forge/noarch     Cached
  + r-stringi                             1.8.1  r43hb55eabf_0       conda-forge/osx-64      854kB
  + r-stringr                             1.5.1  r43h785f33e_0       conda-forge/noarch      296kB
  + r-sys                                 3.4.2  r43h6dc245f_1       conda-forge/osx-64     Cached
  + r-systemfonts                         1.0.5  r43hac8360d_0       conda-forge/osx-64     Cached
  + r-textshaping                         0.3.7  r43hfe21e11_0       conda-forge/osx-64     Cached
  + r-tibble                              3.2.1  r43h6dc245f_2       conda-forge/osx-64     Cached
  + r-tidyr                               1.3.0  r43hac7d2d5_1       conda-forge/osx-64     Cached
  + r-tidyselect                          1.2.0  r43hbe3e9c8_1       conda-forge/osx-64     Cached
  + r-tidyverse                           2.0.0  r43h785f33e_1       conda-forge/noarch     Cached
  + r-timechange                          0.2.0  r43hac7d2d5_1       conda-forge/osx-64     Cached
  + r-tinytex                              0.48  r43hc72bb7e_1       conda-forge/noarch     Cached
  + r-tzdb                                0.4.0  r43hac7d2d5_1       conda-forge/osx-64     Cached
  + r-utf8                                1.2.4  r43hb2c329c_0       conda-forge/osx-64     Cached
  + r-uuid                                1.1_1  r43h6dc245f_0       conda-forge/osx-64     Cached
  + r-vctrs                               0.6.4  r43h64b2c41_0       conda-forge/osx-64     Cached
  + r-venneuler                           1.1_3  r43hc72bb7e_2       conda-forge/noarch     Cached
  + r-viridislite                         0.4.2  r43hc72bb7e_1       conda-forge/noarch     Cached
  + r-vroom                               1.6.4  r43hac7d2d5_0       conda-forge/osx-64     Cached
  + r-withr                               2.5.2  r43hc72bb7e_0       conda-forge/noarch     Cached
  + r-xfun                                 0.41  r43h64b2c41_0       conda-forge/osx-64     Cached
  + r-xml2                                1.3.5  r43h2e0d1c5_0       conda-forge/osx-64     Cached
  + r-yaml                                2.3.7  r43h6dc245f_1       conda-forge/osx-64     Cached
  + readline                                8.2  h9e318b2_1          conda-forge/osx-64     Cached
  + setuptools                           68.2.2  pyhd8ed1ab_0        conda-forge/noarch     Cached
  + sigtool                               0.1.3  h88f4db0_0          conda-forge/osx-64     Cached
  + tapi                              1100.0.11  h9ce4665_0          conda-forge/osx-64     Cached
  + tk                                   8.6.13  h1abcd95_1          conda-forge/osx-64     Cached
  + tktable                                2.10  ha166976_5          conda-forge/osx-64     Cached
  + toml                                 0.10.2  pyhd8ed1ab_0        conda-forge/noarch     Cached
  + tomlkit                              0.12.3  pyha770c72_0        conda-forge/noarch       37kB
  + tzdata                                2023c  h71feb2d_0          conda-forge/noarch     Cached
  + wheel                                0.41.3  pyhd8ed1ab_0        conda-forge/noarch     Cached
  + xmltodict                            0.13.0  pyhd8ed1ab_0        conda-forge/noarch     Cached
  + xz                                    5.4.2  h6c40b1e_0          pkgs/main/osx-64       Cached
  + yaml                                  0.2.5  h0d85af4_2          conda-forge/osx-64     Cached
  + yq                                    3.2.3  pyhd8ed1ab_0        conda-forge/noarch     Cached
  + zlib                                 1.2.13  h8a1eda9_5          conda-forge/osx-64     Cached
  + zstd                                  1.5.5  h829000d_0          conda-forge/osx-64     Cached

  Summary:

  Install: 222 packages

  Total download: 47MB

────────────────────────────────────────────────────────────────────────────────────────────────────


Confirm changes: [Y/n] Y
gmp                                                520.9kB @   6.5MB/s  0.1s
llvm-openmp                                        299.8kB @   3.1MB/s  0.1s
clang-16                                           695.1kB @   6.1MB/s  0.1s
harfbuzz                                             1.3MB @  11.3MB/s  0.1s
r-stringr                                          296.0kB @   2.1MB/s  0.0s
tomlkit                                             37.1kB @ 251.6kB/s  0.0s
clangxx                                             21.3kB @ 111.4kB/s  0.1s
libclang-cpp16                                      12.8MB @  50.5MB/s  0.3s
r-matrix                                             4.0MB @  15.3MB/s  0.1s
r-stringi                                          853.9kB @   3.3MB/s  0.2s
r-ragg                                             397.0kB @   1.4MB/s  0.1s
argcomplete                                         39.6kB @ 133.5kB/s  0.0s
clang                                               21.1kB @  70.0kB/s  0.0s
r-base                                              25.3MB @  40.8MB/s  0.4s

Downloading and Extracting Packages

Preparing transaction: done
Verifying transaction: done
Executing transaction: /
To install TinyTeX with `tinytex::install_tinytex()` the system must have a functional Perl
installation with a `File::Find` module. Most end-user systems will already satisfy this
requirement; however, some minimal contexts (e.g., containers) may not. Perl is available
via Conda Forge as the package `perl`. See https://github.com/rstudio/tinytex/issues/419


done

To activate this environment, use

     $ mamba activate test_overlap_env

To deactivate an active environment, use

     $ mamba deactivate
```
</details>
<br />
