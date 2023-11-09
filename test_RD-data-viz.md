
`#test_RD-data-viz.md`
<br />

<details>
<summary><b><font size="+2"><i>Table of contents</i></font></b></summary>
<!-- MarkdownTOC -->

1. [Name of section](#name-of-section)
    1. [Code](#code)
    1. [Printed](#printed)

<!-- /MarkdownTOC -->
</details>
<br />
<br />

<a id="name-of-section"></a>
## Name of section
<a id="code"></a>
### Code
<details>
<summary><i>Code: Name of section</i></summary>

```bash
#!/bin/bash

grabnode  # 4, 80, 1, N


#  Initialize necessary variables =============================================
   p_RD="${HOME}/tsukiyamalab/Rachel"             # echo "${p_data}"
  p_exp="misc_data/experiments/ChIPs/HDAC_HAT_Q"  # echo "${p_data}"
 d_base="${p_RD}/${p_exp}"                        # ls -lhaFG "${d_base}"
 p_data="230915_ChIPseq"                          # echo "${p_data}"
 d_work="${d_base}/results"                       # ls -lhaFG "${d_work}"

   p_eo="${d_work}/err_out"                       # ls -lhaFG "${p_eo}"
  d_bws="${d_work}/bws"                           # ls -lhaFG "${d_bws}"
 d_bams="${d_work}/bams"                          # ls -lhaFG "${d_bams}"
d_plots="${d_work}/plots"                         # ls -lhaFG "${d_plots}"

threads="${SLURM_CPUS_ON_NODE}"                   # echo "${threads}"

#  Using GSE95556_Sc.cerevisiae.feature.anno_Steinmetz_2013.gtf, Christine's
#+ BED file
bed="GSE95556_Sc.cerevisiae.feature.anno_Steinmetz_2013.gtf"  # ls -lhaFG "${bed}"

check_variables=true
if ${check_variables}; then
    echo """
       p_RD=${p_RD}
      p_exp=${p_exp}
     d_base=${d_base}
     p_data=${p_data}
     d_work=${d_work}

       p_eo=${p_eo}
      d_bws=${d_bws}
     d_bams=${d_bams}
    d_plots=${d_plots}

    threads=${threads}
    """
fi



#  Run operations ===========================================================
#  Go to RD's work directory
cd "${p_RD}/${p_exp}/${d_work}" || echo "cd'ing failed; check on this"

test_KA=true
if ${test_KA}; then
    #  For KA's tests, change locations to store outfiles and stderr/stdout
    #+ files
       p_KA="${HOME}/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses"
    d_plots="${p_KA}/test_RD-data-viz"
       p_eo="${d_plots}/err_out"

    check_variables=true
    if ${check_variables}; then
        echo """
           p_RD=${p_RD}
           p_KA=${p_KA}
          p_exp=${p_exp}
         d_base=${d_base}
         p_data=${p_data}
         d_work=${d_work}

           p_eo=${p_eo}
          d_bws=${d_bws}
         d_bams=${d_bams}
        d_plots=${d_plots}

        threads=${threads}
        """
    fi

    if [[ ! -d "${d_plots}" ]]; then mkdir -p "${d_plots}/err_out"; fi
else
    #  Otherwise, create directory for outfiles in RD's work directory
    if [[ ! -d ${plots} ]]; then mkdir plots; fi
fi

#  Source deepTools if haven't already
if [[ $(command -v deeptools > /dev/null) ]]; then
    module purge
    ml deepTools/3.5.1-foss-2021b
fi


#  Make the prerequisite files for plots --------------------------------------
print_test=true
if ${print_test}; then
    echo """
    computeMatrix reference-point \\
        --referencePoint TSS \\
        -b 1000 -a 1000 \\
        -R \"${d_base}/${bed}\" \\
        -S \"${d_bws}/5781_IP.bamCompare-to-input_no_norm.bw\" \\
           \"${d_bws}/7692_IP.bamCompare-to-input_no_norm.bw\" \\
           \"${d_bws}/7709_IP.bamCompare-to-input_no_norm.bw\" \\
        --skipZeros \\
        -o \"${d_plots}/Gcn5_reps_Steinmetz.gz\" \\
        -p \"${threads}\" \\
        --outFileSortedRegions \"${d_plots}/Gcn5_reps_Steinmetz.bed\"
    """
fi

run=true
if ${run}; then
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 1000 -a 1000 \
        -R "${d_base}/${bed}" \
        -S "${d_bws}/5781_IP.bamCompare-to-input_no_norm.bw" \
           "${d_bws}/7692_IP.bamCompare-to-input_no_norm.bw" \
           "${d_bws}/7709_IP.bamCompare-to-input_no_norm.bw" \
        --skipZeros \
        -o "${d_plots}/Gcn5_reps_Steinmetz.gz" \
        -p "${threads}" \
        --outFileSortedRegions "${d_plots}/Gcn5_reps_Steinmetz.bed"
fi
#  ^ gives the following error: TypeError: cannot pickle '_io.TextIOWrapper'
#+ object
#+ 
#+ Internet search suggests that this is due to the version of python and the
#+ multiple input files fixed the variable name d_bw to d_bws, which fixed both
#+ this and the issue mentioned below
#+ 
#+ KA, 2023-1109: All good.

computeMatrix reference-point \
    --referencePoint TSS \
    -b 1000 -a 1000 \
    -R "${p_RD}/${p_exp}/+1Nuc_all_sacCer3.bed" \
    -S "${d_bws}/7692_IP.bamCompare-to-input_no_norm.bw" \
    --skipZeros \
    -o "${d_plots}/Gcn5_reps_+1Nuc.gz" \
    -p "${threads}" \
    --outFileSortedRegions "${d_plots}/Gcn5_7692_+1Nuc.bed"

ls -lhaFG "${d_bws}/"*"_no_norm.bw"
computeMatrix reference-point \
    --referencePoint TSS \
    -b 1000 -a 1000 \
    -R "GSE95556_Sc.cerevisiae.feature.anno_Steinmetz_2013.gtf" \
    -S "${d_bws}/"*"_no_norm.bw" \
    --skipZeros \
    -o "${d_plots}/all_bws_Steinmetz.gz" \
    -p "${threads}" \
    --outFileSortedRegions "${d_plots}/all_bws_Steinmetz.bed"
#  ^ both the above pieces of code give the following error:
#+ xmlrpc.client.ProtocolError: <ProtocolError for deepblue.mpi-inf.mpg.de/xmlrpc: 503 Service Unavailable>
#+ 
#+ Quick internet search suggests that this is due to an incorrect path to the
#+ reference file (BED or GTF), but I know those files are there...
#+ Upon closer reading of the error code, it could be an issue with the bw
#+  files
#+ 
#+ -> turns out I mistyped the d_bws variable as d_bw in my code... so it
#+ couldn't find the files
#+ 
#+ KA, 2023-1109: Both all good.

#  For just Esa1 replicates with untagged control
computeMatrix reference-point \
    --referencePoint TSS \
    -b 1000 -a 1000 \
    -R "${d_base}/GSE95556_Sc.cerevisiae.feature.anno_Steinmetz_2013.gtf" \
    -S "${d_bws}/5781_IP.bamCompare-to-input_no_norm.bw" \
       "${d_bws}/7041_IP.bamCompare-to-input_no_norm.bw" \
       "${d_bws}/7691_IP.bamCompare-to-input_no_norm.bw" \
    --skipZeros \
    -o "${d_plots}/Esa1_reps_Steinmetz.gz" \
    -p "${threads}" \
    --outFileSortedRegions "${d_plots}/Esa1_reps_Steinmetz.bed"

#  For just Rpd3 replicates with untagged control
computeMatrix reference-point --referencePoint TSS \
    -b 1000 -a 1000 \
    -R "${d_base}/GSE95556_Sc.cerevisiae.feature.anno_Steinmetz_2013.gtf" \
    -S "${d_bws}/5781_IP.bamCompare-to-input_no_norm.bw" \
        "${d_bws}/7568_IP.bamCompare-to-input_no_norm.bw" \
        "${d_bws}/7569_IP.bamCompare-to-input_no_norm.bw" \
    --skipZeros \
    -o "${d_plots}/Rpd3_reps_Steinmetz.gz" \
    -p "${threads}" \
    --outFileSortedRegions "${d_plots}/Rpd3_reps_Steinmetz.bed"


#  Make a heat map with all files ---------------------------------------------
plotHeatmap \
    -m "${d_plots}/all_bws_Steinmetz.gz" \
    -out "${d_plots}/all_bws_Steinmetz.png" \
    --samplesLabel 5781_untagged 7041_Esa1 7568_Rpd3 7569_Rpd3 7691_Esa1 7692_Gcn5 7709_Gcn5 \
    --plotTitle untagged_Gcn5_Esa1_Rpd3_Q \
    --dpi 600
#  KA, 2023-1109: Does not work. Try again.

#  KA, 2023-1109: Try again w/o --samplesLabel.
plotHeatmap \
    -m "${d_plots}/all_bws_Steinmetz.gz" \
    -out "${d_plots}/all_bws_Steinmetz.png" \
    --plotTitle untagged_Gcn5_Esa1_Rpd3_Q \
    --dpi 600
#  KA, 2023-1109: Now it works.

plotHeatmap \
    -m "${d_plots}/Gcn5_reps_Steinmetz.gz" \
    -out "${d_plots}/Gcn5_reps_Steinmetz.png" \
    --samplesLabel 5781_untagged 7692_Gcn5 7709_Gcn5 \
    --plotTitle untagged_Gcn5_Q \
    --dpi 600

plotHeatmap \
    -m "${d_plots}/Esa1_reps_Steinmetz.gz" \
    -out "${d_plots}/Esa1_reps_Steinmetz.png" \
    --samplesLabel 5781_untagged 7041_Esa1 7691_Esa1 \
    --plotTitle untagged_Esa1_Q \
    --dpi 600

plotHeatmap \
    -m "${d_plots}/Rpd3_reps_Steinmetz.gz" \
    -out "${d_plots}/Rpd3_reps_Steinmetz.png" \
    --samplesLabel 5781_untagged 7568_Rpd3 7569_Rpd3 \
    --plotTitle untagged_Rpd3_Q \
    --dpi 600


#  ----------------------------------------------------------------------------
### 24 October 2023 - now working with the merged bio rep files ###

#  All merged replicates
computeMatrix reference-point \
    --referencePoint TSS \
    -b 1000 -a 1000 \
    -R "${d_base}/GSE95556_Sc.cerevisiae.feature.anno_Steinmetz_2013.gtf" \
    -S "${d_bws}/Esa1_merged_IP.bamCompare-to-input_no_norm.bw" \
        "${d_bws}/Gcn5_merged_IP.bamCompare-to-input_no_norm.bw" \
        "${d_bws}/Rpd3_merged_IP.bamCompare-to-input_no_norm.bw" \
    --skipZeros \
    -o "${d_plots}/all_merged_Steinmetz.gz" \
    -p "${threads}" \
    --outFileSortedRegions "${d_plots}/all_merged_Steinmetz.bed"

plotHeatmap -m "${d_plots}/all_merged_Steinmetz.gz" \
    -out "${d_plots}/all_merged_Steinmetz.png" \
    --samplesLabel Esa1 Gcn5 Rpd3 \
    --plotTitle HATs_HDAC_merged_reps_Q \
    --dpi 600

#  Each Gcn5 replicate and the merged file
computeMatrix reference-point \
    --referencePoint TSS \
    -b 1000 -a 1000 \
    -R "${d_base}/GSE95556_Sc.cerevisiae.feature.anno_Steinmetz_2013.gtf" \
    -S "${d_bws}/7692_IP.bamCompare-to-input_no_norm.bw" \
        "${d_bws}/7709_IP.bamCompare-to-input_no_norm.bw" \
        "${d_bws}/Gcn5_merged_IP.bamCompare-to-input_no_norm.bw" \
    --skipZeros \
    -o "${d_plots}/Gcn5_reps_merge_Steinmetz.gz" \
    -p "${threads}" \
    --outFileSortedRegions "${d_plots}/Gcn5_reps_merge_Steinmetz.bed"

plotHeatmap \
    -m "${d_plots}/Gcn5_reps_merge_Steinmetz.gz" \
    -out "${d_plots}/Gcn5_reps_merge_Steinmetz.png" \
    --samplesLabel 7692 7709 Gcn5_merge \
    --plotTitle Gcn5_bio_reps_merge \
    --dpi 600

#  Each Esa1 replicate and the merged file
computeMatrix reference-point \
    --referencePoint TSS \
    -b 1000 -a 1000 \
    -R "${d_base}/GSE95556_Sc.cerevisiae.feature.anno_Steinmetz_2013.gtf" \
    -S "${d_bws}/7041_IP.bamCompare-to-input_no_norm.bw" \
        "${d_bws}/7691_IP.bamCompare-to-input_no_norm.bw" \
        "${d_bws}/Esa1_merged_IP.bamCompare-to-input_no_norm.bw" \
    --skipZeros \
    -o "${d_plots}/Esa1_reps_merge_Steinmetz.gz" \
    -p "${threads}" \
    --outFileSortedRegions "${d_plots}/Esa1_reps_merge_Steinmetz.bed"

plotHeatmap \
    -m "${d_plots}/Esa1_reps_merge_Steinmetz.gz" \
    -out "${d_plots}/Esa1_reps_merge_Steinmetz.png" \
    --samplesLabel 7041 7691 Esa1_merge \
    --plotTitle Esa1_bio_reps_merge \
    --dpi 600

#  Each Rpd3 replicate and the merged file
computeMatrix reference-point --referencePoint TSS \
    -b 1000 -a 1000 \
    -R "${d_base}/GSE95556_Sc.cerevisiae.feature.anno_Steinmetz_2013.gtf" \
    -S "${d_bws}/7568_IP.bamCompare-to-input_no_norm.bw" \
        "${d_bws}/7569_IP.bamCompare-to-input_no_norm.bw" \
        "${d_bws}/Rpd3_merged_IP.bamCompare-to-input_no_norm.bw" \
    --skipZeros \
    -o "${d_plots}/Rpd3_reps_merge_Steinmetz.gz" \
    -p "${threads}" \
    --outFileSortedRegions "${d_plots}/Rpd3_reps_merge_Steinmetz.bed"

plotHeatmap \
    -m "${d_plots}"/Rpd3_reps_merge_Steinmetz.gz \
    -out "${d_plots}"/Rpd3_reps_merge_Steinmetz.png \
    --samplesLabel 7568 7569 Rpd3_merge \
    --plotTitle Rpd3_bio_reps_merge \
    --dpi 600

#  All merged replicates along with untagged control
computeMatrix reference-point \
    --referencePoint TSS \
    -b 1000 -a 1000 \
    -R "${d_base}/GSE95556_Sc.cerevisiae.feature.anno_Steinmetz_2013.gtf" \
    -S "${d_bws}/5781_IP.bamCompare-to-input_no_norm.bw" \
        "${d_bws}/Esa1_merged_IP.bamCompare-to-input_no_norm.bw" \
        "${d_bws}/Gcn5_merged_IP.bamCompare-to-input_no_norm.bw" \
        "${d_bws}/Rpd3_merged_IP.bamCompare-to-input_no_norm.bw" \
    --skipZeros \
    -o "${d_plots}/all_merged_untagged_Steinmetz.gz" \
    -p "${threads}" \
    --outFileSortedRegions "${d_plots}/all_merged_untagged_Steinmetz.bed"

plotHeatmap \
    -m "${d_plots}"/all_merged_untagged_Steinmetz.gz \
    -out "${d_plots}"/all_merged_untagged_Steinmetz.png \
    --samplesLabel untagged Esa1 Gcn5 Rpd3 \
    --plotTitle HATs_HDAC_merged_untagged_reps_Q \
    --dpi 600

#  n.b. all Steinmetz computeMatrix files give the following output:
#+ Skipping CUT442, due to being absent in the computeMatrix output.
#+ Skipping XUT_14R-344, due to being absent in the computeMatrix output.
#+ Skipping SUT760, due to being absent in the computeMatrix output.
#+ Skipping XUT_14F-324, due to being absent in the computeMatrix output.
#+ Skipping YGL263W, due to being absent in the computeMatrix output.


#  ----------------------------------------------------------------------------
### 31 October 2023 - now trying k-means clustering of heat maps ###

#  KA, 2023-1109: Minor refactoring: Addition of for loops for k
print_test=true
run_command=true

iter=0
unset stems && typeset -a stems=(
    "all_merged_untagged_Steinmetz"  # Merged factors w/untagged control
    "all_merged_Steinmetz"           # Merged factors w/o untagged control
)
for k in 2 3 4 5 6 7; do
    for stem in "${stems[@]}"; do
        (( iter++ ))

        if [[ "${stem}" == "all_merged_untagged_Steinmetz" ]]; then
            samples=(untagged Esa1 Gcn5 Rpd3)
        elif [[ "${stem}" == "all_merged_Steinmetz" ]]; then
            samples=(Esa1 Gcn5 Rpd3)
        fi

        if ${print_test}; then
            echo """
            ### ${iter} ###
            plotHeatmap \\
                -m \"${d_plots}/${stem}.gz\" \\
                -out \"${d_plots}/${stem}_kmeans-${k}.png\" \\
                --outFileSortedRegions \"${d_plots}/${stem}_kmeans-${k}.bed\" \\
                --outFileNameMatrix \"${d_plots}/${stem}_kmeans-${k}.mat.gz\" \\
                --samplesLabel ${samples[*]} \\
                --plotTitle \"HATs_HDAC_${stem##all_}_reps_Q\" \\
                --dpi 600 \\
                --kmeans ${k}
            """
        fi

        if ${run_command}; then
            plotHeatmap \
                -m "${d_plots}/${stem}.gz" \
                -out "${d_plots}/${stem}_kmeans-${k}.png" \
                --outFileSortedRegions "${d_plots}/${stem}_kmeans-${k}.bed" \
                --outFileNameMatrix "${d_plots}/${stem}_kmeans-${k}.mat.gz" \
                --samplesLabel ${samples[*]} \
                --plotTitle "HATs_HDAC_${stem##all_}_reps_Q" \
                --dpi 600 \
                --kmeans ${k}
        fi
    done
done


#  ----------------------------------------------------------------------------
### 1 November 2023 - making heatmaps from normalized bws ###

#  All unmerged files with MAPQ1 normalization
computeMatrix reference-point \
    --referencePoint TSS \
    -b 1000 -a 1000 \
    -R GSE95556_Sc.cerevisiae.feature.anno_Steinmetz_2013.gtf \
    -S "${d_bws}/5781_IP.bamCompare-to-input.SF1-coef-0.112792.bw" \
       "${d_bws}/7041_IP.bamCompare-to-input.SF1-coef-0.788518.bw" \
       "${d_bws}/7691_IP.bamCompare-to-input.SF1-coef-0.96677.bw" \
       "${d_bws}/7568_IP.bamCompare-to-input.SF1-coef-0.649159.bw" \
       "${d_bws}/7569_IP.bamCompare-to-input.SF1-coef-1.208895.bw" \
       "${d_bws}/7692_IP.bamCompare-to-input.SF1-coef-0.639539.bw" \
       "${d_bws}/7709_IP.bamCompare-to-input.SF1-coef-0.784.bw" \
    --skipZeros \
    -o "${d_plots}/all_bws_SF1_Steinmetz.gz" \
    -p "${threads}" \
    --outFileSortedRegions "${d_plots}/all_bws_SF1_Steinmetz.bed"

plotHeatmap \
    -m "${d_plots}/all_bws_SF1_Steinmetz.gz" \
    -out "${d_plots}/all_bws_SF1_Steinmetz.png" \
    --outFileSortedRegions "${d_plots}/all_bws_SF1_Steinmetz.bed" \
    --outFileNameMatrix "${d_plots}/all_bws_SF1_Steinmetz.mat.gz" \
    --samplesLabel 5781_untagged 7041_Esa1 7691_Esa1 7568_Rpd3 7569_Rpd3 7692_Gcn5 7709_Gcn5 \
    --plotTitle untagged_Esa1_Rpd3_Gcn5_SF1_Q \
    --dpi 600

#  All merged files with MAPQ1 normalization
computeMatrix reference-point \
    --referencePoint TSS \
    -b 1000 -a 1000 \
    -R "${d_base}/GSE95556_Sc.cerevisiae.feature.anno_Steinmetz_2013.gtf" \
    -S "${d_bws}/Esa1_merged_IP.bamCompare-to-input.SF1-coef-0.880503.bw" \
       "${d_bws}/Rpd3_merged_IP.bamCompare-to-input.SF1-coef-0.969325.bw" \
       "${d_bws}/Gcn5_merged_IP.bamCompare-to-input.SF1-coef-0.695127.bw" \
    --skipZeros \
    -o "${d_plots}/all_merged_SF1_Steinmetz.gz" \
    -p "${threads}" \
    --outFileSortedRegions "${d_plots}/all_merged_SF1_Steinmetz.bed"

plotHeatmap \
    -m "${d_plots}/all_merged_SF1_Steinmetz.gz" \
    -out "${d_plots}/all_merged_SF1_Steinmetz.png" \
    --outFileSortedRegions "${d_plots}/all_merged_SF1_Steinmetz.bed" \
    --outFileNameMatrix "${d_plots}/all_merged_SF1_Steinmetz.mat.gz" \
    --samplesLabel Esa1 Rpd3 Gcn5 \
    --plotTitle HATs_HDAC_merged_SF1_reps_Q \
    --dpi 600


#  All unmerged files with MAPQ2 normalization
computeMatrix reference-point \
    --referencePoint TSS \
    -b 1000 -a 1000 \
    -R GSE95556_Sc.cerevisiae.feature.anno_Steinmetz_2013.gtf \
    -S "${d_bws}/5781_IP.bamCompare-to-input.SF2-coef-0.109063.bw" \
       "${d_bws}/7041_IP.bamCompare-to-input.SF2-coef-0.837016.bw" \
       "${d_bws}/7691_IP.bamCompare-to-input.SF2-coef-0.949774.bw" \
       "${d_bws}/7568_IP.bamCompare-to-input.SF2-coef-0.649534.bw" \
       "${d_bws}/7569_IP.bamCompare-to-input.SF2-coef-1.310394.bw" \
       "${d_bws}/7692_IP.bamCompare-to-input.SF2-coef-0.616505.bw" \
       "${d_bws}/7709_IP.bamCompare-to-input.SF2-coef-0.831745.bw" \
    --skipZeros \
    -o "${d_plots}/all_bws_SF2_Steinmetz.gz" \
    -p "${threads}" \
    --outFileSortedRegions "${d_plots}/all_bws_SF2_Steinmetz.bed"

plotHeatmap \
    -m "${d_plots}/all_bws_SF2_Steinmetz.gz" \
    -out "${d_plots}/all_bws_SF2_Steinmetz.png" \
    --outFileSortedRegions "${d_plots}/all_bws_SF2_Steinmetz.bed" \
    --outFileNameMatrix "${d_plots}/all_bws_SF2_Steinmetz.mat.gz" \
    --samplesLabel 5781_untagged 7041_Esa1 7691_Esa1 7568_Rpd3 7569_Rpd3 7692_Gcn5 7709_Gcn5 \
    --plotTitle untagged_Esa1_Rpd3_Gcn5_SF2_Q \
    --dpi 600

#  All merged files with MAPQ2 normalization
computeMatrix reference-point \
    --referencePoint TSS \
    -b 1000 -a 1000 \
    -R "${d_base}/GSE95556_Sc.cerevisiae.feature.anno_Steinmetz_2013.gtf" \
    -S "${d_bws}/Esa1_merged_IP.bamCompare-to-input.SF2-coef-0.898049.bw" \
       "${d_bws}/Rpd3_merged_IP.bamCompare-to-input.SF2-coef-1.021932.bw" \
       "${d_bws}/Gcn5_merged_IP.bamCompare-to-input.SF2-coef-0.699061.bw" \
    --skipZeros \
    -o "${d_plots}/all_merged_SF2_Steinmetz.gz" \
    -p "${threads}" \
    --outFileSortedRegions "${d_plots}/all_merged_SF2_Steinmetz.bed"

plotHeatmap \
    -m "${d_plots}/all_merged_SF2_Steinmetz.gz" \
    -out "${d_plots}/all_merged_SF2_Steinmetz.png" \
    --outFileSortedRegions "${d_plots}/all_merged_SF2_Steinmetz.bed" \
    --outFileNameMatrix "${d_plots}/all_merged_SF2_Steinmetz.mat.gz" \
    --samplesLabel Esa1 Rpd3 Gcn5 \
    --plotTitle HATs_HDAC_merged_SF2_reps_Q \
    --dpi 600

#  All unmerged files with MAPQ3 normalization
computeMatrix reference-point \
    --referencePoint TSS \
    -b 1000 -a 1000 \
    -R GSE95556_Sc.cerevisiae.feature.anno_Steinmetz_2013.gtf \
    -S "${d_bws}/5781_IP.bamCompare-to-input.SF3-coef-0.109045.bw" \
       "${d_bws}/7041_IP.bamCompare-to-input.SF3-coef-0.837038.bw" \
       "${d_bws}/7691_IP.bamCompare-to-input.SF3-coef-0.949745.bw" \
       "${d_bws}/7568_IP.bamCompare-to-input.SF3-coef-0.649477.bw" \
       "${d_bws}/7569_IP.bamCompare-to-input.SF3-coef-1.310305.bw" \
       "${d_bws}/7692_IP.bamCompare-to-input.SF3-coef-0.616394.bw" \
       "${d_bws}/7709_IP.bamCompare-to-input.SF3-coef-0.831602.bw" \
    --skipZeros \
    -o "${d_plots}/all_bws_SF3_Steinmetz.gz" \
    -p "${threads}" \
    --outFileSortedRegions "${d_plots}/all_bws_SF3_Steinmetz.bed"

plotHeatmap \
    -m "${d_plots}/all_bws_SF3_Steinmetz.gz" \
    -out "${d_plots}/all_bws_SF3_Steinmetz.png" \
    --outFileSortedRegions "${d_plots}/all_bws_SF3_Steinmetz.bed" \
    --outFileNameMatrix "${d_plots}/all_bws_SF3_Steinmetz.mat.gz" \
    --samplesLabel 5781_untagged 7041_Esa1 7691_Esa1 7568_Rpd3 7569_Rpd3 7692_Gcn5 7709_Gcn5 \
    --plotTitle untagged_Esa1_Rpd3_Gcn5_SF3_Q \
    --dpi 600

#  All merged files with MAPQ3 normalization
computeMatrix reference-point \
    --referencePoint TSS \
    -b 1000 -a 1000 \
    -R "${d_base}/GSE95556_Sc.cerevisiae.feature.anno_Steinmetz_2013.gtf" \
    -S "${d_bws}/Esa1_merged_IP.bamCompare-to-input.SF3-coef-0.898045.bw" \
       "${d_bws}/Rpd3_merged_IP.bamCompare-to-input.SF3-coef-1.021858.bw" \
       "${d_bws}/Gcn5_merged_IP.bamCompare-to-input.SF3-coef-0.698936.bw" \
    --skipZeros \
    -o "${d_plots}/all_merged_SF3_Steinmetz.gz" \
    -p "${threads}" \
    --outFileSortedRegions "${d_plots}/all_merged_SF3_Steinmetz.bed"

plotHeatmap \
    -m "${d_plots}/all_merged_SF3_Steinmetz.gz" \
    -out "${d_plots}/all_merged_SF3_Steinmetz.png" \
    --outFileSortedRegions "${d_plots}/all_merged_SF3_Steinmetz.bed" \
    --outFileNameMatrix "${d_plots}/all_merged_SF3_Steinmetz.mat.gz" \
    --samplesLabel Esa1 Rpd3 Gcn5 \
    --plotTitle HATs_HDAC_merged_SF3_reps_Q \
    --dpi 600


#  ----------------------------------------------------------------------------
### 3 November 2023 - making k-means-clustered heatmaps with normalized merged files ###
print_test=true
run_command=true
iter=0
for MAPQ in 1 2 3; do
    for k in 2 3 4 5 6 7; do
        (( iter++ ))
        if ${print_test}; then
            echo """
            ### ${iter} ###
            plotHeatmap \\
                -m \"${d_plots}/all_merged_SF${MAPQ}_Steinmetz.gz\" \\
                -out \"${d_plots}/all_merged_SF${MAPQ}_Steinmetz_kmeans-${k}.png\" \\
                --outFileSortedRegions \"${d_plots}/all_merged_SF${MAPQ}_Steinmetz_kmeans-${k}.bed\" \\
                --outFileNameMatrix \"${d_plots}/all_merged_SF${MAPQ}_Steinmetz_kmeans-${k}.mat.gz\" \\
                --samplesLabel Esa1 Rpd3 Gcn5 \\
                --plotTitle \"HATs_HDAC_merged_SF${MAPQ}_reps_Q\" \\
                --dpi 600 \\
                --kmeans ${k}
            """
        fi

        if ${run_command}; then
            plotHeatmap \
                -m "${d_plots}/all_merged_SF${MAPQ}_Steinmetz.gz" \
                -out "${d_plots}/all_merged_SF${MAPQ}_Steinmetz_kmeans-${k}.png" \
                --outFileSortedRegions "${d_plots}/all_merged_SF${MAPQ}_Steinmetz_kmeans-${k}.bed" \
                --outFileNameMatrix "${d_plots}/all_merged_SF${MAPQ}_Steinmetz_kmeans-${k}.mat.gz" \
                --samplesLabel Esa1 Rpd3 Gcn5 \
                --plotTitle "HATs_HDAC_merged_SF${MAPQ}_reps_Q" \
                --dpi 600 \
                --kmeans ${k}
        fi
    done
done


#  ----------------------------------------------------------------------------
### 6 November 2023 - trying plotCorrelation with merged files ###

#  Compare all MAPQ2 files to each other
multiBigwigSummary bins \
    -b "${d_bws}/5781_IP.bamCompare-to-input.SF2-coef-0.109063.bw" \
       "${d_bws}/7041_IP.bamCompare-to-input.SF2-coef-0.837016.bw" \
       "${d_bws}/7691_IP.bamCompare-to-input.SF2-coef-0.949774.bw" \
       "${d_bws}/7568_IP.bamCompare-to-input.SF2-coef-0.649534.bw" \
       "${d_bws}/7569_IP.bamCompare-to-input.SF2-coef-1.310394.bw" \
       "${d_bws}/7692_IP.bamCompare-to-input.SF2-coef-0.616505.bw" \
       "${d_bws}/7709_IP.bamCompare-to-input.SF2-coef-0.831745.bw" \
    --labels untagged 7041_Esa1 7691_Esa1 7568_Rpd3 7569_Rpd3 7692_Gcn5 7709_Gcn5 \
    -o "${d_plots}/all_MAPQ2file_compare.npz" \
    -p "${threads}"

#  Compare all merged (and untagged control) MAPQ2 files to each other
multiBigwigSummary bins \
    -b "${d_bws}/5781_IP.bamCompare-to-input.SF2-coef-0.109063.bw" \
       "${d_bws}/Esa1_merged_IP.bamCompare-to-input.SF2-coef-0.898049.bw" \
       "${d_bws}/Rpd3_merged_IP.bamCompare-to-input.SF2-coef-1.021932.bw" \
       "${d_bws}/Gcn5_merged_IP.bamCompare-to-input.SF2-coef-0.699061.bw" \
    --labels untagged Esa1 Rpd3 Gcn5 \
    -o "${d_plots}/merged_untagged_MAPQ2_compare.npz" \
    -p "${threads}"

#  Compare all merged MAPQ2 files to each other (no untagged control)
multiBigwigSummary bins \
    -b "${d_bws}/Esa1_merged_IP.bamCompare-to-input.SF2-coef-0.898049.bw" \
       "${d_bws}/Rpd3_merged_IP.bamCompare-to-input.SF2-coef-1.021932.bw" \
       "${d_bws}/Gcn5_merged_IP.bamCompare-to-input.SF2-coef-0.699061.bw" \
    --labels Esa1 Rpd3 Gcn5 \
    -o "${d_plots}/merged_MAPQ2_compare.npz" \
    -p "${threads}"

#  Draw plots
print_test=true
run_command=true
iter=0

unset stems && typeset -a stems=(
    "all_MAPQ2file_compare"
    "merged_untagged_MAPQ2_compare"
    "merged_MAPQ2_compare"
)
for cor in "pearson" "spearman"; do
    for what in "scatterplot" "heatmap"; do
        for stem in "${stems[@]}"; do
            (( iter++ ))
            
            if ${print_test}; then
                echo """
                ### ${iter} ###
                plotCorrelation \\
                    --corData \"${d_plots}/${stem}.npz\" \\
                    --corMethod \"${cor}\" \\
                    --whatToPlot \"${what}\" \\
                    --plotFile \"${d_plots}/${stem}_${cor}_${what}.png\"
                """
            fi

            if ${run_command}; then
                plotCorrelation \
                    --corData "${d_plots}/${stem}.npz" \
                    --corMethod "${cor}" \
                    --whatToPlot "${what}" \
                    --plotFile "${d_plots}/${stem}_${cor}_${what}.png"
            fi
        done
    done
done


#  Compress the directory in which all the outfiles are stored ================
if [[ -d "${d_plots}" && ! -f "${d_plots}.tar.gz" ]]; then
    tar czf "${d_plots}.tar.gz" "${d_plots}"
fi

if [[ -d "${d_plots}" && -f "${d_plots}.tar.gz" ]]; then
    rm -rf "${d_plots}"
fi
```
</details>
<br />

<a id="printed"></a>
### Printed
<details>
<summary><i>Printed: Name of section (selected operations)</i></summary>

```txt
❯ if ${check_variables}; then
>     echo """
>        p_RD=${p_RD}
>       p_exp=${p_exp}
>      d_base=${d_base}
>      p_data=${p_data}
>      d_work=${d_work}
> 
>        p_eo=${p_eo}
>       d_bws=${d_bws}
>      d_bams=${d_bams}
>     d_plots=${d_plots}
> 
>     threads=${threads}
>     """
> fi

       p_RD=/home/kalavatt/tsukiyamalab/Rachel
      p_exp=misc_data/experiments/ChIPs/HDAC_HAT_Q
     d_base=/home/kalavatt/tsukiyamalab/Rachel/misc_data/experiments/ChIPs/HDAC_HAT_Q
     p_data=230915_ChIPseq
     d_work=/home/kalavatt/tsukiyamalab/Rachel/misc_data/experiments/ChIPs/HDAC_HAT_Q/results

       p_eo=/home/kalavatt/tsukiyamalab/Rachel/misc_data/experiments/ChIPs/HDAC_HAT_Q/results/err_out
      d_bws=/home/kalavatt/tsukiyamalab/Rachel/misc_data/experiments/ChIPs/HDAC_HAT_Q/results/bws
     d_bams=/home/kalavatt/tsukiyamalab/Rachel/misc_data/experiments/ChIPs/HDAC_HAT_Q/results/bams
    d_plots=/home/kalavatt/tsukiyamalab/Rachel/misc_data/experiments/ChIPs/HDAC_HAT_Q/results/plots

    threads=4


❯ if ${test_KA}; then
>     #  For KA's tests, change locations to store outfiles and stderr/stdout
>     #+ files
>        p_KA="${HOME}/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses"
>     d_plots="${p_KA}/test_RD-data-viz"
>        p_eo="${d_plots}/err_out"
> 
>     check_variables=true
>     if ${check_variables}; then
>         echo """
>            p_RD=${p_RD}
>            p_KA=${p_KA}
>           p_exp=${p_exp}
>          d_base=${d_base}
>          p_data=${p_data}
>          d_work=${d_work}
> 
>            p_eo=${p_eo}
>           d_bws=${d_bws}
>          d_bams=${d_bams}
>         d_plots=${d_plots}
> 
>         threads=${threads}
>         """
>     fi
> 
>     if [[ ! -d "${d_plots}" ]]; then mkdir -p "${d_plots}/err_out"; fi
> else
>     #  Otherwise, create directory for outfiles in RD's work directory
>     if [[ ! -d ${plots} ]]; then mkdir plots; fi
> fi

           p_RD=/home/kalavatt/tsukiyamalab/Rachel
           p_KA=/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses
          p_exp=misc_data/experiments/ChIPs/HDAC_HAT_Q
         d_base=/home/kalavatt/tsukiyamalab/Rachel/misc_data/experiments/ChIPs/HDAC_HAT_Q
         p_data=230915_ChIPseq
         d_work=/home/kalavatt/tsukiyamalab/Rachel/misc_data/experiments/ChIPs/HDAC_HAT_Q/results

           p_eo=/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/err_out
          d_bws=/home/kalavatt/tsukiyamalab/Rachel/misc_data/experiments/ChIPs/HDAC_HAT_Q/results/bws
         d_bams=/home/kalavatt/tsukiyamalab/Rachel/misc_data/experiments/ChIPs/HDAC_HAT_Q/results/bams
        d_plots=/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz

        threads=4


❯ #  Make the prerequisite files for plots --------------------------------------
❯ print_test=true
❯ if ${print_test}; then
>     echo """
>     computeMatrix reference-point \\
>         --referencePoint TSS \\
>         -b 1000 -a 1000 \\
>         -R \"${d_base}/${bed}\" \\
>         -S \"${d_bws}/5781_IP.bamCompare-to-input_no_norm.bw\" \\
>            \"${d_bws}/7692_IP.bamCompare-to-input_no_norm.bw\" \\
>            \"${d_bws}/7709_IP.bamCompare-to-input_no_norm.bw\" \\
>         --skipZeros \\
>         -o \"${d_plots}/Gcn5_reps_Steinmetz.gz\" \\
>         -p \"${threads}\" \\
>         --outFileSortedRegions \"${d_plots}/Gcn5_reps_Steinmetz.bed\"
>     """
> fi

    computeMatrix reference-point \
        --referencePoint TSS \
        -b 1000 -a 1000 \
        -R "/home/kalavatt/tsukiyamalab/Rachel/misc_data/experiments/ChIPs/HDAC_HAT_Q/GSE95556_Sc.cerevisiae.feature.anno_Steinmetz_2013.gtf" \
        -S "/home/kalavatt/tsukiyamalab/Rachel/misc_data/experiments/ChIPs/HDAC_HAT_Q/results/bws/5781_IP.bamCompare-to-input_no_norm.bw" \
           "/home/kalavatt/tsukiyamalab/Rachel/misc_data/experiments/ChIPs/HDAC_HAT_Q/results/bws/7692_IP.bamCompare-to-input_no_norm.bw" \
           "/home/kalavatt/tsukiyamalab/Rachel/misc_data/experiments/ChIPs/HDAC_HAT_Q/results/bws/7709_IP.bamCompare-to-input_no_norm.bw" \
        --skipZeros \
        -o "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/Gcn5_reps_Steinmetz.gz" \
        -p "4" \
        --outFileSortedRegions "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/Gcn5_reps_Steinmetz.bed"


❯ run=true
❯ if ${run}; then
>     computeMatrix reference-point \
>         --referencePoint TSS \
>         -b 1000 -a 1000 \
>         -R "${d_base}/${bed}" \
>         -S "${d_bws}/5781_IP.bamCompare-to-input_no_norm.bw" \
>            "${d_bws}/7692_IP.bamCompare-to-input_no_norm.bw" \
>            "${d_bws}/7709_IP.bamCompare-to-input_no_norm.bw" \
>         --skipZeros \
>         -o "${d_plots}/Gcn5_reps_Steinmetz.gz" \
>         -p "${threads}" \
>         --outFileSortedRegions "${d_plots}/Gcn5_reps_Steinmetz.bed"
> fi
Skipping CUT442, due to being absent in the computeMatrix output.
Skipping XUT_14R-344, due to being absent in the computeMatrix output.
Skipping SUT760, due to being absent in the computeMatrix output.
Skipping XUT_14F-324, due to being absent in the computeMatrix output.
Skipping YGL263W, due to being absent in the computeMatrix output.


❯ computeMatrix reference-point \
>     --referencePoint TSS \
>     -b 1000 -a 1000 \
>     -R "${p_RD}/${p_exp}/+1Nuc_all_sacCer3.bed" \
>     -S "${d_bws}/7692_IP.bamCompare-to-input_no_norm.bw" \
>     --skipZeros \
>     -o "${d_plots}/Gcn5_reps_+1Nuc.gz" \
>     -p "${threads}" \
>     --outFileSortedRegions "${d_plots}/Gcn5_7692_+1Nuc.bed"


❯ ls -lhaFG "${d_bws}/"*"_no_norm.bw"
-rw-rw---- 1 rdell  93M Oct  3 12:39 /home/kalavatt/tsukiyamalab/Rachel/misc_data/experiments/ChIPs/HDAC_HAT_Q/results/bws/5781_IP.bamCompare-to-input_no_norm.bw
-rw-rw---- 1 rdell 101M Oct  3 12:45 /home/kalavatt/tsukiyamalab/Rachel/misc_data/experiments/ChIPs/HDAC_HAT_Q/results/bws/7041_IP.bamCompare-to-input_no_norm.bw
-rw-rw---- 1 rdell 100M Oct  3 12:48 /home/kalavatt/tsukiyamalab/Rachel/misc_data/experiments/ChIPs/HDAC_HAT_Q/results/bws/7568_IP.bamCompare-to-input_no_norm.bw
-rw-rw---- 1 rdell 104M Oct  3 12:53 /home/kalavatt/tsukiyamalab/Rachel/misc_data/experiments/ChIPs/HDAC_HAT_Q/results/bws/7569_IP.bamCompare-to-input_no_norm.bw
-rw-rw---- 1 rdell 102M Oct  3 12:37 /home/kalavatt/tsukiyamalab/Rachel/misc_data/experiments/ChIPs/HDAC_HAT_Q/results/bws/7691_IP.bamCompare-to-input_no_norm.bw
-rw-rw---- 1 rdell 102M Oct  3 12:42 /home/kalavatt/tsukiyamalab/Rachel/misc_data/experiments/ChIPs/HDAC_HAT_Q/results/bws/7692_IP.bamCompare-to-input_no_norm.bw
-rw-rw---- 1 rdell  99M Oct  3 12:50 /home/kalavatt/tsukiyamalab/Rachel/misc_data/experiments/ChIPs/HDAC_HAT_Q/results/bws/7709_IP.bamCompare-to-input_no_norm.bw
-rw-rw---- 1 rdell 108M Oct 24 10:50 /home/kalavatt/tsukiyamalab/Rachel/misc_data/experiments/ChIPs/HDAC_HAT_Q/results/bws/Esa1_merged_IP.bamCompare-to-input_no_norm.bw
-rw-rw---- 1 rdell 107M Oct 24 10:42 /home/kalavatt/tsukiyamalab/Rachel/misc_data/experiments/ChIPs/HDAC_HAT_Q/results/bws/Gcn5_merged_IP.bamCompare-to-input_no_norm.bw
-rw-rw---- 1 rdell 109M Oct 24 10:46 /home/kalavatt/tsukiyamalab/Rachel/misc_data/experiments/ChIPs/HDAC_HAT_Q/results/bws/Rpd3_merged_IP.bamCompare-to-input_no_norm.bw


❯ computeMatrix reference-point \
>     --referencePoint TSS \
>     -b 1000 -a 1000 \
>     -R "GSE95556_Sc.cerevisiae.feature.anno_Steinmetz_2013.gtf" \
>     -S "${d_bws}/"*"_no_norm.bw" \
>     --skipZeros \
>     -o "${d_plots}/all_bws_Steinmetz.gz" \
>     -p "${threads}" \
>     --outFileSortedRegions "${d_plots}/all_bws_Steinmetz.bed"
Skipping CUT442, due to being absent in the computeMatrix output.
Skipping XUT_14R-344, due to being absent in the computeMatrix output.
Skipping SUT760, due to being absent in the computeMatrix output.
Skipping XUT_14F-324, due to being absent in the computeMatrix output.
Skipping YGL263W, due to being absent in the computeMatrix output.


❯ #  For just Esa1 replicates with untagged control
❯ computeMatrix reference-point \
>     --referencePoint TSS \
>     -b 1000 -a 1000 \
>     -R "${d_base}/GSE95556_Sc.cerevisiae.feature.anno_Steinmetz_2013.gtf" \
>     -S "${d_bws}/5781_IP.bamCompare-to-input_no_norm.bw" \
>        "${d_bws}/7041_IP.bamCompare-to-input_no_norm.bw" \
>        "${d_bws}/7691_IP.bamCompare-to-input_no_norm.bw" \
>     --skipZeros \
>     -o "${d_plots}/Esa1_reps_Steinmetz.gz" \
>     -p "${threads}" \
>     --outFileSortedRegions "${d_plots}/Esa1_reps_Steinmetz.bed"
Skipping CUT442, due to being absent in the computeMatrix output.
Skipping XUT_14R-344, due to being absent in the computeMatrix output.
Skipping SUT760, due to being absent in the computeMatrix output.
Skipping XUT_14F-324, due to being absent in the computeMatrix output.
Skipping YGL263W, due to being absent in the computeMatrix output.


❯ #  For just Rpd3 replicates with untagged control
❯ computeMatrix reference-point --referencePoint TSS \
>     -b 1000 -a 1000 \
>     -R "${d_base}/GSE95556_Sc.cerevisiae.feature.anno_Steinmetz_2013.gtf" \
>     -S "${d_bws}/5781_IP.bamCompare-to-input_no_norm.bw" \
>         "${d_bws}/7568_IP.bamCompare-to-input_no_norm.bw" \
>         "${d_bws}/7569_IP.bamCompare-to-input_no_norm.bw" \
>     --skipZeros \
>     -o "${d_plots}/Rpd3_reps_Steinmetz.gz" \
>     -p "${threads}" \
>     --outFileSortedRegions "${d_plots}/Rpd3_reps_Steinmetz.bed"
Skipping CUT442, due to being absent in the computeMatrix output.
Skipping XUT_14R-344, due to being absent in the computeMatrix output.
Skipping SUT760, due to being absent in the computeMatrix output.
Skipping XUT_14F-324, due to being absent in the computeMatrix output.
Skipping YGL263W, due to being absent in the computeMatrix output.


❯ #  Make a heat map with all files ---------------------------------------------
❯ plotHeatmap \
>     -m "${d_plots}"/all_bws_Steinmetz.gz \
>     -out "${d_plots}"/all_bws_Steinmetz.png \
>     --samplesLabel 5781_untagged 7041_Esa1 7568_Rpd3 7569_Rpd3 7691_Esa1 7692_Gcn5 7709_Gcn5 \
>     --plotTitle untagged_Gcn5_Esa1_Rpd3_Q \
>     --dpi 600
Traceback (most recent call last):
  File "/app/software/deepTools/3.5.1-foss-2021b/bin/plotHeatmap", line 12, in <module>
    main(args)
  File "/app/software/deepTools/3.5.1-foss-2021b/lib/python3.9/site-packages/deeptools/plotHeatmap.py", line 840, in main
    hm.matrix.set_sample_labels(args.samplesLabel)
  File "/app/software/deepTools/3.5.1-foss-2021b/lib/python3.9/site-packages/deeptools/heatmapper.py", line 1164, in set_sample_labels
    raise ValueError("length new labels != length original labels")
ValueError: length new labels != length original labels
❯ #  KA, 2023-1109: Does not work. Try again.

❯ #  KA, 2023-1109: Try again w/o --samplesLabel.
❯ plotHeatmap \
>     -m "${d_plots}/all_bws_Steinmetz.gz" \
>     -out "${d_plots}/all_bws_Steinmetz.png" \
>     --plotTitle untagged_Gcn5_Esa1_Rpd3_Q \
>     --dpi 600


❯ plotHeatmap \
>     -m "${d_plots}/Gcn5_reps_Steinmetz.gz" \
>     -out "${d_plots}/Gcn5_reps_Steinmetz.png" \
>     --samplesLabel 5781_untagged 7692_Gcn5 7709_Gcn5 \
>     --plotTitle untagged_Gcn5_Q \
>     --dpi 600


❯ plotHeatmap \
>     -m "${d_plots}/Esa1_reps_Steinmetz.gz" \
>     -out "${d_plots}/Esa1_reps_Steinmetz.png" \
>     --samplesLabel 5781_untagged 7041_Esa1 7691_Esa1 \
>     --plotTitle untagged_Esa1_Q \
>     --dpi 600


❯ plotHeatmap \
>     -m "${d_plots}/Rpd3_reps_Steinmetz.gz" \
>     -out "${d_plots}/Rpd3_reps_Steinmetz.png" \
>     --samplesLabel 5781_untagged 7568_Rpd3 7569_Rpd3 \
>     --plotTitle untagged_Rpd3_Q \
>     --dpi 600


❯ #  ----------------------------------------------------------------------------
❯ ### 24 October 2023 - now working with the merged bio rep files ###


❯ #  All merged replicates
❯ computeMatrix reference-point \
>     --referencePoint TSS \
>     -b 1000 -a 1000 \
>     -R "${d_base}/GSE95556_Sc.cerevisiae.feature.anno_Steinmetz_2013.gtf" \
>     -S "${d_bws}/Esa1_merged_IP.bamCompare-to-input_no_norm.bw" \
>         "${d_bws}/Gcn5_merged_IP.bamCompare-to-input_no_norm.bw" \
>         "${d_bws}/Rpd3_merged_IP.bamCompare-to-input_no_norm.bw" \
>     --skipZeros \
>     -o "${d_plots}/all_merged_Steinmetz.gz" \
>     -p "${threads}" \
>     --outFileSortedRegions "${d_plots}/all_merged_Steinmetz.bed"
Skipping CUT442, due to being absent in the computeMatrix output.
Skipping XUT_14R-344, due to being absent in the computeMatrix output.
Skipping SUT760, due to being absent in the computeMatrix output.
Skipping XUT_14F-324, due to being absent in the computeMatrix output.
Skipping YGL263W, due to being absent in the computeMatrix output.


❯ plotHeatmap -m "${d_plots}/all_merged_Steinmetz.gz" \
>     -out "${d_plots}/all_merged_Steinmetz.png" \
>     --samplesLabel Esa1 Gcn5 Rpd3 \
>     --plotTitle HATs_HDAC_merged_reps_Q \
>     --dpi 600


❯ #  Each Gcn5 replicate and the merged file
❯ computeMatrix reference-point \
>     --referencePoint TSS \
>     -b 1000 -a 1000 \
>     -R "${d_base}/GSE95556_Sc.cerevisiae.feature.anno_Steinmetz_2013.gtf" \
>     -S "${d_bws}/7692_IP.bamCompare-to-input_no_norm.bw" \
>         "${d_bws}/7709_IP.bamCompare-to-input_no_norm.bw" \
>         "${d_bws}/Gcn5_merged_IP.bamCompare-to-input_no_norm.bw" \
>     --skipZeros \
>     -o "${d_plots}/Gcn5_reps_merge_Steinmetz.gz" \
>     -p "${threads}" \
>     --outFileSortedRegions "${d_plots}/Gcn5_reps_merge_Steinmetz.bed"
Skipping CUT442, due to being absent in the computeMatrix output.
Skipping XUT_14R-344, due to being absent in the computeMatrix output.
Skipping SUT760, due to being absent in the computeMatrix output.
Skipping XUT_14F-324, due to being absent in the computeMatrix output.
Skipping YGL263W, due to being absent in the computeMatrix output.


❯ plotHeatmap \
>     -m "${d_plots}"/Gcn5_reps_merge_Steinmetz.gz \
>     -out "${d_plots}"/Gcn5_reps_merge_Steinmetz.png \
>     --samplesLabel 7692 7709 Gcn5_merge \
>     --plotTitle Gcn5_bio_reps_merge \
>     --dpi 600


❯ #  Each Esa1 replicate and the merged file
❯ computeMatrix reference-point \
>     --referencePoint TSS \
>     -b 1000 -a 1000 \
>     -R "${d_base}/GSE95556_Sc.cerevisiae.feature.anno_Steinmetz_2013.gtf" \
>     -S "${d_bws}/7041_IP.bamCompare-to-input_no_norm.bw" \
>         "${d_bws}/7691_IP.bamCompare-to-input_no_norm.bw" \
>         "${d_bws}/Esa1_merged_IP.bamCompare-to-input_no_norm.bw" \
>     --skipZeros \
>     -o "${d_plots}/Esa1_reps_merge_Steinmetz.gz" \
>     -p "${threads}" \
>     --outFileSortedRegions "${d_plots}/Esa1_reps_merge_Steinmetz.bed"
Skipping CUT442, due to being absent in the computeMatrix output.
Skipping XUT_14R-344, due to being absent in the computeMatrix output.
Skipping SUT760, due to being absent in the computeMatrix output.
Skipping XUT_14F-324, due to being absent in the computeMatrix output.
Skipping YGL263W, due to being absent in the computeMatrix output.


❯ plotHeatmap \
>     -m "${d_plots}/Esa1_reps_merge_Steinmetz.gz" \
>     -out "${d_plots}/Esa1_reps_merge_Steinmetz.png" \
>     --samplesLabel 7041 7691 Esa1_merge \
>     --plotTitle Esa1_bio_reps_merge \
>     --dpi 600


❯ #  Each Rpd3 replicate and the merged file
❯ computeMatrix reference-point --referencePoint TSS \
>     -b 1000 -a 1000 \
>     -R "${d_base}/GSE95556_Sc.cerevisiae.feature.anno_Steinmetz_2013.gtf" \
>     -S "${d_bws}/7568_IP.bamCompare-to-input_no_norm.bw" \
>         "${d_bws}/7569_IP.bamCompare-to-input_no_norm.bw" \
>         "${d_bws}/Rpd3_merged_IP.bamCompare-to-input_no_norm.bw" \
>     --skipZeros \
>     -o "${d_plots}/Rpd3_reps_merge_Steinmetz.gz" \
>     -p "${threads}" \
>     --outFileSortedRegions "${d_plots}/Rpd3_reps_merge_Steinmetz.bed"
Skipping CUT442, due to being absent in the computeMatrix output.
Skipping XUT_14R-344, due to being absent in the computeMatrix output.
Skipping SUT760, due to being absent in the computeMatrix output.
Skipping XUT_14F-324, due to being absent in the computeMatrix output.
Skipping YGL263W, due to being absent in the computeMatrix output.


❯ plotHeatmap \
>     -m "${d_plots}"/Rpd3_reps_merge_Steinmetz.gz \
>     -out "${d_plots}"/Rpd3_reps_merge_Steinmetz.png \
>     --samplesLabel 7568 7569 Rpd3_merge \
>     --plotTitle Rpd3_bio_reps_merge \
>     --dpi 600


❯ #  All merged replicates along with untagged control
❯ computeMatrix reference-point \
>     --referencePoint TSS \
>     -b 1000 -a 1000 \
>     -R "${d_base}/GSE95556_Sc.cerevisiae.feature.anno_Steinmetz_2013.gtf" \
>     -S "${d_bws}/5781_IP.bamCompare-to-input_no_norm.bw" \
>         "${d_bws}/Esa1_merged_IP.bamCompare-to-input_no_norm.bw" \
>         "${d_bws}/Gcn5_merged_IP.bamCompare-to-input_no_norm.bw" \
>         "${d_bws}/Rpd3_merged_IP.bamCompare-to-input_no_norm.bw" \
>     --skipZeros \
>     -o "${d_plots}/all_merged_untagged_Steinmetz.gz" \
>     -p "${threads}" \
>     --outFileSortedRegions "${d_plots}/all_merged_untagged_Steinmetz.bed"
Skipping CUT442, due to being absent in the computeMatrix output.
Skipping XUT_14R-344, due to being absent in the computeMatrix output.
Skipping SUT760, due to being absent in the computeMatrix output.
Skipping XUT_14F-324, due to being absent in the computeMatrix output.
Skipping YGL263W, due to being absent in the computeMatrix output.


❯ plotHeatmap \
>     -m "${d_plots}"/all_merged_untagged_Steinmetz.gz \
>     -out "${d_plots}"/all_merged_untagged_Steinmetz.png \
>     --samplesLabel untagged Esa1 Gcn5 Rpd3 \
>     --plotTitle HATs_HDAC_merged_untagged_reps_Q \
>     --dpi 600


❯ #  ----------------------------------------------------------------------------
❯ ### 31 October 2023 - now trying k-means clustering of heat maps ###

❯ #  KA, 2023-1109: Minor refactoring: Addition of for loops for k

❯ print_test=true
❯ run_command=true
 
❯ iter=0
❯ unset stems && typeset -a stems=(
>     "all_merged_untagged_Steinmetz"  # Merged factors w/untagged control
>     "all_merged_Steinmetz"           # Merged factors w/o untagged control
> )


❯ for k in 2 3 4 5 6 7; do
>     for stem in "${stems[@]}"; do
>         (( iter++ ))
>         if ${print_test}; then
>             echo """
>             ### ${iter} ###
>             plotHeatmap \\
>                 -m \"${d_plots}/${stem}.gz\" \\
>                 -out \"${d_plots}/${stem}_kmeans-${k}.png\" \\
>                 --outFileSortedRegions \"${d_plots}/${stem}_kmeans-${k}.bed\" \\
>                 --outFileNameMatrix \"${d_plots}/${stem}_kmeans-${k}.mat.gz\" \\
>                 --samplesLabel untagged Esa1 Gcn5 Rpd3 \\
>                 --plotTitle \"HATs_HDAC_${stem##all_}_reps_Q\" \\
>                 --dpi 600 \\
>                 --kmeans ${k}
>             """
>         fi
> 
>         if ${run_command}; then
>             plotHeatmap \
>                 -m "${d_plots}/${stem}.gz" \
>                 -out "${d_plots}/${stem}_kmeans-${k}.png" \
>                 --outFileSortedRegions "${d_plots}/${stem}_kmeans-${k}.bed" \
>                 --outFileNameMatrix "${d_plots}/${stem}_kmeans-${k}.mat.gz" \
>                 --samplesLabel untagged Esa1 Gcn5 Rpd3 \
>                 --plotTitle "HATs_HDAC_${stem##all_}_reps_Q" \
>                 --dpi 600 \
>                 --kmeans ${k}
>         fi
>     done
> done

            ### 1 ###
            plotHeatmap \
                -m "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_untagged_Steinmetz.gz" \
                -out "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_untagged_Steinmetz_kmeans-2.png" \
                --outFileSortedRegions "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_untagged_Steinmetz_kmeans-2.bed" \
                --outFileNameMatrix "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_untagged_Steinmetz_kmeans-2.mat.gz" \
                --samplesLabel untagged Esa1 Gcn5 Rpd3 \
                --plotTitle "HATs_HDAC_merged_untagged_Steinmetz_reps_Q" \
                --dpi 600 \
                --kmeans 2


            ### 2 ###
            plotHeatmap \
                -m "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_Steinmetz.gz" \
                -out "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_Steinmetz_kmeans-2.png" \
                --outFileSortedRegions "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_Steinmetz_kmeans-2.bed" \
                --outFileNameMatrix "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_Steinmetz_kmeans-2.mat.gz" \
                --samplesLabel Esa1 Gcn5 Rpd3 \
                --plotTitle "HATs_HDAC_merged_Steinmetz_reps_Q" \
                --dpi 600 \
                --kmeans 2


            ### 3 ###
            plotHeatmap \
                -m "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_untagged_Steinmetz.gz" \
                -out "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_untagged_Steinmetz_kmeans-3.png" \
                --outFileSortedRegions "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_untagged_Steinmetz_kmeans-3.bed" \
                --outFileNameMatrix "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_untagged_Steinmetz_kmeans-3.mat.gz" \
                --samplesLabel untagged Esa1 Gcn5 Rpd3 \
                --plotTitle "HATs_HDAC_merged_untagged_Steinmetz_reps_Q" \
                --dpi 600 \
                --kmeans 3


            ### 4 ###
            plotHeatmap \
                -m "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_Steinmetz.gz" \
                -out "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_Steinmetz_kmeans-3.png" \
                --outFileSortedRegions "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_Steinmetz_kmeans-3.bed" \
                --outFileNameMatrix "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_Steinmetz_kmeans-3.mat.gz" \
                --samplesLabel Esa1 Gcn5 Rpd3 \
                --plotTitle "HATs_HDAC_merged_Steinmetz_reps_Q" \
                --dpi 600 \
                --kmeans 3


            ### 5 ###
            plotHeatmap \
                -m "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_untagged_Steinmetz.gz" \
                -out "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_untagged_Steinmetz_kmeans-4.png" \
                --outFileSortedRegions "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_untagged_Steinmetz_kmeans-4.bed" \
                --outFileNameMatrix "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_untagged_Steinmetz_kmeans-4.mat.gz" \
                --samplesLabel untagged Esa1 Gcn5 Rpd3 \
                --plotTitle "HATs_HDAC_merged_untagged_Steinmetz_reps_Q" \
                --dpi 600 \
                --kmeans 4


            ### 6 ###
            plotHeatmap \
                -m "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_Steinmetz.gz" \
                -out "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_Steinmetz_kmeans-4.png" \
                --outFileSortedRegions "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_Steinmetz_kmeans-4.bed" \
                --outFileNameMatrix "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_Steinmetz_kmeans-4.mat.gz" \
                --samplesLabel Esa1 Gcn5 Rpd3 \
                --plotTitle "HATs_HDAC_merged_Steinmetz_reps_Q" \
                --dpi 600 \
                --kmeans 4


            ### 7 ###
            plotHeatmap \
                -m "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_untagged_Steinmetz.gz" \
                -out "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_untagged_Steinmetz_kmeans-5.png" \
                --outFileSortedRegions "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_untagged_Steinmetz_kmeans-5.bed" \
                --outFileNameMatrix "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_untagged_Steinmetz_kmeans-5.mat.gz" \
                --samplesLabel untagged Esa1 Gcn5 Rpd3 \
                --plotTitle "HATs_HDAC_merged_untagged_Steinmetz_reps_Q" \
                --dpi 600 \
                --kmeans 5


            ### 8 ###
            plotHeatmap \
                -m "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_Steinmetz.gz" \
                -out "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_Steinmetz_kmeans-5.png" \
                --outFileSortedRegions "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_Steinmetz_kmeans-5.bed" \
                --outFileNameMatrix "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_Steinmetz_kmeans-5.mat.gz" \
                --samplesLabel Esa1 Gcn5 Rpd3 \
                --plotTitle "HATs_HDAC_merged_Steinmetz_reps_Q" \
                --dpi 600 \
                --kmeans 5


            ### 9 ###
            plotHeatmap \
                -m "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_untagged_Steinmetz.gz" \
                -out "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_untagged_Steinmetz_kmeans-6.png" \
                --outFileSortedRegions "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_untagged_Steinmetz_kmeans-6.bed" \
                --outFileNameMatrix "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_untagged_Steinmetz_kmeans-6.mat.gz" \
                --samplesLabel untagged Esa1 Gcn5 Rpd3 \
                --plotTitle "HATs_HDAC_merged_untagged_Steinmetz_reps_Q" \
                --dpi 600 \
                --kmeans 6


            ### 10 ###
            plotHeatmap \
                -m "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_Steinmetz.gz" \
                -out "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_Steinmetz_kmeans-6.png" \
                --outFileSortedRegions "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_Steinmetz_kmeans-6.bed" \
                --outFileNameMatrix "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_Steinmetz_kmeans-6.mat.gz" \
                --samplesLabel Esa1 Gcn5 Rpd3 \
                --plotTitle "HATs_HDAC_merged_Steinmetz_reps_Q" \
                --dpi 600 \
                --kmeans 6


            ### 11 ###
            plotHeatmap \
                -m "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_untagged_Steinmetz.gz" \
                -out "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_untagged_Steinmetz_kmeans-7.png" \
                --outFileSortedRegions "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_untagged_Steinmetz_kmeans-7.bed" \
                --outFileNameMatrix "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_untagged_Steinmetz_kmeans-7.mat.gz" \
                --samplesLabel untagged Esa1 Gcn5 Rpd3 \
                --plotTitle "HATs_HDAC_merged_untagged_Steinmetz_reps_Q" \
                --dpi 600 \
                --kmeans 7


            ### 12 ###
            plotHeatmap \
                -m "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_Steinmetz.gz" \
                -out "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_Steinmetz_kmeans-7.png" \
                --outFileSortedRegions "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_Steinmetz_kmeans-7.bed" \
                --outFileNameMatrix "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_Steinmetz_kmeans-7.mat.gz" \
                --samplesLabel Esa1 Gcn5 Rpd3 \
                --plotTitle "HATs_HDAC_merged_Steinmetz_reps_Q" \
                --dpi 600 \
                --kmeans 7


❯ #  ----------------------------------------------------------------------------
❯ ### 1 November 2023 - making heatmaps from normalized bws ###

❯ #  All unmerged files with MAPQ1 normalization
❯ computeMatrix reference-point \
>     --referencePoint TSS \
>     -b 1000 -a 1000 \
>     -R GSE95556_Sc.cerevisiae.feature.anno_Steinmetz_2013.gtf \
>     -S "${d_bws}/5781_IP.bamCompare-to-input.SF1-coef-0.112792.bw" \
>        "${d_bws}/7041_IP.bamCompare-to-input.SF1-coef-0.788518.bw" \
>        "${d_bws}/7691_IP.bamCompare-to-input.SF1-coef-0.96677.bw" \
>        "${d_bws}/7568_IP.bamCompare-to-input.SF1-coef-0.649159.bw" \
>        "${d_bws}/7569_IP.bamCompare-to-input.SF1-coef-1.208895.bw" \
>        "${d_bws}/7692_IP.bamCompare-to-input.SF1-coef-0.639539.bw" \
>        "${d_bws}/7709_IP.bamCompare-to-input.SF1-coef-0.784.bw" \
>     --skipZeros \
>     -o "${d_plots}/all_bws_SF1_Steinmetz.gz" \
>     -p "${threads}" \
>     --outFileSortedRegions "${d_plots}/all_bws_SF1_Steinmetz.bed"
Skipping CUT442, due to being absent in the computeMatrix output.
Skipping XUT_14R-344, due to being absent in the computeMatrix output.
Skipping SUT760, due to being absent in the computeMatrix output.
Skipping XUT_14F-324, due to being absent in the computeMatrix output.
Skipping YGL263W, due to being absent in the computeMatrix output.


❯ plotHeatmap \
>     -m "${d_plots}/all_bws_SF1_Steinmetz.gz" \
>     -out "${d_plots}/all_bws_SF1_Steinmetz.png" \
>     --outFileSortedRegions "${d_plots}/all_bws_SF1_Steinmetz.bed" \
>     --outFileNameMatrix "${d_plots}/all_bws_SF1_Steinmetz.mat.gz" \
>     --samplesLabel 5781_untagged 7041_Esa1 7691_Esa1 7568_Rpd3 7569_Rpd3 7692_Gcn5 7709_Gcn5 \
>     --plotTitle untagged_Esa1_Rpd3_Gcn5_SF1_Q \
>     --dpi 600


❯ #  All merged files with MAPQ1 normalization
❯ computeMatrix reference-point \
>     --referencePoint TSS \
>     -b 1000 -a 1000 \
>     -R "${d_base}/GSE95556_Sc.cerevisiae.feature.anno_Steinmetz_2013.gtf" \
>     -S "${d_bws}/Esa1_merged_IP.bamCompare-to-input.SF1-coef-0.880503.bw" \
>        "${d_bws}/Rpd3_merged_IP.bamCompare-to-input.SF1-coef-0.969325.bw" \
>        "${d_bws}/Gcn5_merged_IP.bamCompare-to-input.SF1-coef-0.695127.bw" \
>     --skipZeros \
>     -o "${d_plots}/all_merged_SF1_Steinmetz.gz" \
>     -p "${threads}" \
>     --outFileSortedRegions "${d_plots}/all_merged_SF1_Steinmetz.bed"
Skipping CUT442, due to being absent in the computeMatrix output.
Skipping XUT_14R-344, due to being absent in the computeMatrix output.
Skipping SUT760, due to being absent in the computeMatrix output.
Skipping XUT_14F-324, due to being absent in the computeMatrix output.
Skipping YGL263W, due to being absent in the computeMatrix output.


❯ plotHeatmap \
>     -m "${d_plots}/all_merged_SF1_Steinmetz.gz" \
>     -out "${d_plots}/all_merged_SF1_Steinmetz.png" \
>     --outFileSortedRegions "${d_plots}/all_merged_SF1_Steinmetz.bed" \
>     --outFileNameMatrix "${d_plots}/all_merged_SF1_Steinmetz.mat.gz" \
>     --samplesLabel Esa1 Rpd3 Gcn5 \
>     --plotTitle HATs_HDAC_merged_SF1_reps_Q \
>     --dpi 600


❯ #  All unmerged files with MAPQ2 normalization
❯ computeMatrix reference-point \
>     --referencePoint TSS \
>     -b 1000 -a 1000 \
>     -R GSE95556_Sc.cerevisiae.feature.anno_Steinmetz_2013.gtf \
>     -S "${d_bws}/5781_IP.bamCompare-to-input.SF2-coef-0.109063.bw" \
>        "${d_bws}/7041_IP.bamCompare-to-input.SF2-coef-0.837016.bw" \
>        "${d_bws}/7691_IP.bamCompare-to-input.SF2-coef-0.949774.bw" \
>        "${d_bws}/7568_IP.bamCompare-to-input.SF2-coef-0.649534.bw" \
>        "${d_bws}/7569_IP.bamCompare-to-input.SF2-coef-1.310394.bw" \
>        "${d_bws}/7692_IP.bamCompare-to-input.SF2-coef-0.616505.bw" \
>        "${d_bws}/7709_IP.bamCompare-to-input.SF2-coef-0.831745.bw" \
>     --skipZeros \
>     -o "${d_plots}/all_bws_SF2_Steinmetz.gz" \
>     -p "${threads}" \
>     --outFileSortedRegions "${d_plots}/all_bws_SF2_Steinmetz.bed"
Skipping CUT442, due to being absent in the computeMatrix output.
Skipping XUT_14R-344, due to being absent in the computeMatrix output.
Skipping SUT760, due to being absent in the computeMatrix output.
Skipping XUT_14F-324, due to being absent in the computeMatrix output.
Skipping YGL263W, due to being absent in the computeMatrix output.


❯ plotHeatmap \
>     -m "${d_plots}/all_bws_SF2_Steinmetz.gz" \
>     -out "${d_plots}/all_bws_SF2_Steinmetz.png" \
>     --outFileSortedRegions "${d_plots}/all_bws_SF2_Steinmetz.bed" \
>     --outFileNameMatrix "${d_plots}/all_bws_SF2_Steinmetz.mat.gz" \
>     --samplesLabel 5781_untagged 7041_Esa1 7691_Esa1 7568_Rpd3 7569_Rpd3 7692_Gcn5 7709_Gcn5 \
>     --plotTitle untagged_Esa1_Rpd3_Gcn5_SF2_Q \
>     --dpi 600


❯ #  All merged files with MAPQ2 normalization
❯ computeMatrix reference-point \
>     --referencePoint TSS \
>     -b 1000 -a 1000 \
>     -R "${d_base}/GSE95556_Sc.cerevisiae.feature.anno_Steinmetz_2013.gtf" \
>     -S "${d_bws}/Esa1_merged_IP.bamCompare-to-input.SF2-coef-0.898049.bw" \
>        "${d_bws}/Rpd3_merged_IP.bamCompare-to-input.SF2-coef-1.021932.bw" \
>        "${d_bws}/Gcn5_merged_IP.bamCompare-to-input.SF2-coef-0.699061.bw" \
>     --skipZeros \
>     -o "${d_plots}/all_merged_SF2_Steinmetz.gz" \
>     -p "${threads}" \
>     --outFileSortedRegions "${d_plots}/all_merged_SF2_Steinmetz.bed"
Skipping CUT442, due to being absent in the computeMatrix output.
Skipping XUT_14R-344, due to being absent in the computeMatrix output.
Skipping SUT760, due to being absent in the computeMatrix output.
Skipping XUT_14F-324, due to being absent in the computeMatrix output.
Skipping YGL263W, due to being absent in the computeMatrix output.


❯ plotHeatmap \
>     -m "${d_plots}/all_merged_SF2_Steinmetz.gz" \
>     -out "${d_plots}/all_merged_SF2_Steinmetz.png" \
>     --outFileSortedRegions "${d_plots}/all_merged_SF2_Steinmetz.bed" \
>     --outFileNameMatrix "${d_plots}/all_merged_SF2_Steinmetz.mat.gz" \
>     --samplesLabel Esa1 Rpd3 Gcn5 \
>     --plotTitle HATs_HDAC_merged_SF2_reps_Q \
>     --dpi 600


❯ #  All unmerged files with MAPQ3 normalization
❯ computeMatrix reference-point \
>     --referencePoint TSS \
>     -b 1000 -a 1000 \
>     -R GSE95556_Sc.cerevisiae.feature.anno_Steinmetz_2013.gtf \
>     -S "${d_bws}/5781_IP.bamCompare-to-input.SF3-coef-0.109045.bw" \
>        "${d_bws}/7041_IP.bamCompare-to-input.SF3-coef-0.837038.bw" \
>        "${d_bws}/7691_IP.bamCompare-to-input.SF3-coef-0.949745.bw" \
>        "${d_bws}/7568_IP.bamCompare-to-input.SF3-coef-0.649477.bw" \
>        "${d_bws}/7569_IP.bamCompare-to-input.SF3-coef-1.310305.bw" \
>        "${d_bws}/7692_IP.bamCompare-to-input.SF3-coef-0.616394.bw" \
>        "${d_bws}/7709_IP.bamCompare-to-input.SF3-coef-0.831602.bw" \
>     --skipZeros \
>     -o "${d_plots}/all_bws_SF3_Steinmetz.gz" \
>     -p "${threads}" \
>     --outFileSortedRegions "${d_plots}/all_bws_SF3_Steinmetz.bed"
Skipping CUT442, due to being absent in the computeMatrix output.
Skipping XUT_14R-344, due to being absent in the computeMatrix output.
Skipping SUT760, due to being absent in the computeMatrix output.
Skipping XUT_14F-324, due to being absent in the computeMatrix output.
Skipping YGL263W, due to being absent in the computeMatrix output.
 

❯ plotHeatmap \
>     -m "${d_plots}/all_bws_SF3_Steinmetz.gz" \
>     -out "${d_plots}/all_bws_SF3_Steinmetz.png" \
>     --outFileSortedRegions "${d_plots}/all_bws_SF3_Steinmetz.bed" \
>     --outFileNameMatrix "${d_plots}/all_bws_SF3_Steinmetz.mat.gz" \
>     --samplesLabel 5781_untagged 7041_Esa1 7691_Esa1 7568_Rpd3 7569_Rpd3 7692_Gcn5 7709_Gcn5 \
>     --plotTitle untagged_Esa1_Rpd3_Gcn5_SF3_Q \
>     --dpi 600
 

❯ #  All merged files with MAPQ3 normalization
❯ computeMatrix reference-point \
>     --referencePoint TSS \
>     -b 1000 -a 1000 \
>     -R "${d_base}/GSE95556_Sc.cerevisiae.feature.anno_Steinmetz_2013.gtf" \
>     -S "${d_bws}/Esa1_merged_IP.bamCompare-to-input.SF3-coef-0.898045.bw" \
>        "${d_bws}/Rpd3_merged_IP.bamCompare-to-input.SF3-coef-1.021858.bw" \
>        "${d_bws}/Gcn5_merged_IP.bamCompare-to-input.SF3-coef-0.698936.bw" \
>     --skipZeros \
>     -o "${d_plots}/all_merged_SF3_Steinmetz.gz" \
>     -p "${threads}" \
>     --outFileSortedRegions "${d_plots}/all_merged_SF3_Steinmetz.bed"
Skipping CUT442, due to being absent in the computeMatrix output.
Skipping XUT_14R-344, due to being absent in the computeMatrix output.
Skipping SUT760, due to being absent in the computeMatrix output.
Skipping XUT_14F-324, due to being absent in the computeMatrix output.
Skipping YGL263W, due to being absent in the computeMatrix output.


❯ plotHeatmap \
>     -m "${d_plots}/all_merged_SF3_Steinmetz.gz" \
>     -out "${d_plots}/all_merged_SF3_Steinmetz.png" \
>     --outFileSortedRegions "${d_plots}/all_merged_SF3_Steinmetz.bed" \
>     --outFileNameMatrix "${d_plots}/all_merged_SF3_Steinmetz.mat.gz" \
>     --samplesLabel Esa1 Rpd3 Gcn5 \
>     --plotTitle HATs_HDAC_merged_SF3_reps_Q \
>     --dpi 600


❯ #  ----------------------------------------------------------------------------
❯ ### 3 November 2023 - making k-means-clustered heatmaps with normalized merged files ###


❯ for MAPQ in 1 2 3; do
>     for k in 2 3 4 5 6 7; do
>         (( iter++ ))
>         if ${print_test}; then
>             echo """
>             ### ${iter} ###
>             plotHeatmap \\
>                 -m \"${d_plots}/all_merged_SF${MAPQ}_Steinmetz.gz\" \\
>                 -out \"${d_plots}/all_merged_SF${MAPQ}_Steinmetz_kmeans-${k}.png\" \\
>                 --outFileSortedRegions \"${d_plots}/all_merged_SF${MAPQ}_Steinmetz_kmeans-${k}.bed\" \\
>                 --outFileNameMatrix \"${d_plots}/all_merged_SF${MAPQ}_Steinmetz_kmeans-${k}.mat.gz\" \\
>                 --samplesLabel Esa1 Rpd3 Gcn5 \\
>                 --plotTitle \"HATs_HDAC_merged_SF${MAPQ}_reps_Q\" \\
>                 --dpi 600 \\
>                 --kmeans ${k}
>             """
>         fi
> 
>         if ${run_command}; then
>             plotHeatmap \
>                 -m "${d_plots}/all_merged_SF${MAPQ}_Steinmetz.gz" \
>                 -out "${d_plots}/all_merged_SF${MAPQ}_Steinmetz_kmeans-${k}.png" \
>                 --outFileSortedRegions "${d_plots}/all_merged_SF${MAPQ}_Steinmetz_kmeans-${k}.bed" \
>                 --outFileNameMatrix "${d_plots}/all_merged_SF${MAPQ}_Steinmetz_kmeans-${k}.mat.gz" \
>                 --samplesLabel Esa1 Rpd3 Gcn5 \
>                 --plotTitle "HATs_HDAC_merged_SF${MAPQ}_reps_Q" \
>                 --dpi 600 \
>                 --kmeans ${k}
>         fi
>     done
> done

            ### 1 ###
            plotHeatmap \
                -m "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF1_Steinmetz.gz" \
                -out "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF1_Steinmetz_kmeans-2.png" \
                --outFileSortedRegions "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF1_Steinmetz_kmeans-2.bed" \
                --outFileNameMatrix "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF1_Steinmetz_kmeans-2.mat.gz" \
                --samplesLabel Esa1 Rpd3 Gcn5 \
                --plotTitle "HATs_HDAC_merged_SF1_reps_Q" \
                --dpi 600 \
                --kmeans 2


            ### 2 ###
            plotHeatmap \
                -m "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF1_Steinmetz.gz" \
                -out "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF1_Steinmetz_kmeans-3.png" \
                --outFileSortedRegions "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF1_Steinmetz_kmeans-3.bed" \
                --outFileNameMatrix "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF1_Steinmetz_kmeans-3.mat.gz" \
                --samplesLabel Esa1 Rpd3 Gcn5 \
                --plotTitle "HATs_HDAC_merged_SF1_reps_Q" \
                --dpi 600 \
                --kmeans 3


            ### 3 ###
            plotHeatmap \
                -m "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF1_Steinmetz.gz" \
                -out "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF1_Steinmetz_kmeans-4.png" \
                --outFileSortedRegions "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF1_Steinmetz_kmeans-4.bed" \
                --outFileNameMatrix "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF1_Steinmetz_kmeans-4.mat.gz" \
                --samplesLabel Esa1 Rpd3 Gcn5 \
                --plotTitle "HATs_HDAC_merged_SF1_reps_Q" \
                --dpi 600 \
                --kmeans 4


            ### 4 ###
            plotHeatmap \
                -m "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF1_Steinmetz.gz" \
                -out "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF1_Steinmetz_kmeans-5.png" \
                --outFileSortedRegions "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF1_Steinmetz_kmeans-5.bed" \
                --outFileNameMatrix "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF1_Steinmetz_kmeans-5.mat.gz" \
                --samplesLabel Esa1 Rpd3 Gcn5 \
                --plotTitle "HATs_HDAC_merged_SF1_reps_Q" \
                --dpi 600 \
                --kmeans 5


            ### 5 ###
            plotHeatmap \
                -m "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF1_Steinmetz.gz" \
                -out "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF1_Steinmetz_kmeans-6.png" \
                --outFileSortedRegions "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF1_Steinmetz_kmeans-6.bed" \
                --outFileNameMatrix "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF1_Steinmetz_kmeans-6.mat.gz" \
                --samplesLabel Esa1 Rpd3 Gcn5 \
                --plotTitle "HATs_HDAC_merged_SF1_reps_Q" \
                --dpi 600 \
                --kmeans 6


            ### 6 ###
            plotHeatmap \
                -m "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF1_Steinmetz.gz" \
                -out "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF1_Steinmetz_kmeans-7.png" \
                --outFileSortedRegions "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF1_Steinmetz_kmeans-7.bed" \
                --outFileNameMatrix "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF1_Steinmetz_kmeans-7.mat.gz" \
                --samplesLabel Esa1 Rpd3 Gcn5 \
                --plotTitle "HATs_HDAC_merged_SF1_reps_Q" \
                --dpi 600 \
                --kmeans 7


            ### 7 ###
            plotHeatmap \
                -m "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF2_Steinmetz.gz" \
                -out "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF2_Steinmetz_kmeans-2.png" \
                --outFileSortedRegions "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF2_Steinmetz_kmeans-2.bed" \
                --outFileNameMatrix "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF2_Steinmetz_kmeans-2.mat.gz" \
                --samplesLabel Esa1 Rpd3 Gcn5 \
                --plotTitle "HATs_HDAC_merged_SF2_reps_Q" \
                --dpi 600 \
                --kmeans 2


            ### 8 ###
            plotHeatmap \
                -m "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF2_Steinmetz.gz" \
                -out "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF2_Steinmetz_kmeans-3.png" \
                --outFileSortedRegions "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF2_Steinmetz_kmeans-3.bed" \
                --outFileNameMatrix "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF2_Steinmetz_kmeans-3.mat.gz" \
                --samplesLabel Esa1 Rpd3 Gcn5 \
                --plotTitle "HATs_HDAC_merged_SF2_reps_Q" \
                --dpi 600 \
                --kmeans 3


            ### 9 ###
            plotHeatmap \
                -m "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF2_Steinmetz.gz" \
                -out "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF2_Steinmetz_kmeans-4.png" \
                --outFileSortedRegions "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF2_Steinmetz_kmeans-4.bed" \
                --outFileNameMatrix "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF2_Steinmetz_kmeans-4.mat.gz" \
                --samplesLabel Esa1 Rpd3 Gcn5 \
                --plotTitle "HATs_HDAC_merged_SF2_reps_Q" \
                --dpi 600 \
                --kmeans 4


            ### 10 ###
            plotHeatmap \
                -m "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF2_Steinmetz.gz" \
                -out "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF2_Steinmetz_kmeans-5.png" \
                --outFileSortedRegions "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF2_Steinmetz_kmeans-5.bed" \
                --outFileNameMatrix "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF2_Steinmetz_kmeans-5.mat.gz" \
                --samplesLabel Esa1 Rpd3 Gcn5 \
                --plotTitle "HATs_HDAC_merged_SF2_reps_Q" \
                --dpi 600 \
                --kmeans 5


            ### 11 ###
            plotHeatmap \
                -m "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF2_Steinmetz.gz" \
                -out "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF2_Steinmetz_kmeans-6.png" \
                --outFileSortedRegions "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF2_Steinmetz_kmeans-6.bed" \
                --outFileNameMatrix "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF2_Steinmetz_kmeans-6.mat.gz" \
                --samplesLabel Esa1 Rpd3 Gcn5 \
                --plotTitle "HATs_HDAC_merged_SF2_reps_Q" \
                --dpi 600 \
                --kmeans 6


            ### 12 ###
            plotHeatmap \
                -m "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF2_Steinmetz.gz" \
                -out "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF2_Steinmetz_kmeans-7.png" \
                --outFileSortedRegions "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF2_Steinmetz_kmeans-7.bed" \
                --outFileNameMatrix "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF2_Steinmetz_kmeans-7.mat.gz" \
                --samplesLabel Esa1 Rpd3 Gcn5 \
                --plotTitle "HATs_HDAC_merged_SF2_reps_Q" \
                --dpi 600 \
                --kmeans 7


            ### 13 ###
            plotHeatmap \
                -m "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF3_Steinmetz.gz" \
                -out "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF3_Steinmetz_kmeans-2.png" \
                --outFileSortedRegions "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF3_Steinmetz_kmeans-2.bed" \
                --outFileNameMatrix "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF3_Steinmetz_kmeans-2.mat.gz" \
                --samplesLabel Esa1 Rpd3 Gcn5 \
                --plotTitle "HATs_HDAC_merged_SF3_reps_Q" \
                --dpi 600 \
                --kmeans 2


            ### 14 ###
            plotHeatmap \
                -m "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF3_Steinmetz.gz" \
                -out "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF3_Steinmetz_kmeans-3.png" \
                --outFileSortedRegions "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF3_Steinmetz_kmeans-3.bed" \
                --outFileNameMatrix "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF3_Steinmetz_kmeans-3.mat.gz" \
                --samplesLabel Esa1 Rpd3 Gcn5 \
                --plotTitle "HATs_HDAC_merged_SF3_reps_Q" \
                --dpi 600 \
                --kmeans 3


            ### 15 ###
            plotHeatmap \
                -m "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF3_Steinmetz.gz" \
                -out "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF3_Steinmetz_kmeans-4.png" \
                --outFileSortedRegions "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF3_Steinmetz_kmeans-4.bed" \
                --outFileNameMatrix "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF3_Steinmetz_kmeans-4.mat.gz" \
                --samplesLabel Esa1 Rpd3 Gcn5 \
                --plotTitle "HATs_HDAC_merged_SF3_reps_Q" \
                --dpi 600 \
                --kmeans 4


            ### 16 ###
            plotHeatmap \
                -m "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF3_Steinmetz.gz" \
                -out "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF3_Steinmetz_kmeans-5.png" \
                --outFileSortedRegions "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF3_Steinmetz_kmeans-5.bed" \
                --outFileNameMatrix "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF3_Steinmetz_kmeans-5.mat.gz" \
                --samplesLabel Esa1 Rpd3 Gcn5 \
                --plotTitle "HATs_HDAC_merged_SF3_reps_Q" \
                --dpi 600 \
                --kmeans 5


            ### 17 ###
            plotHeatmap \
                -m "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF3_Steinmetz.gz" \
                -out "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF3_Steinmetz_kmeans-6.png" \
                --outFileSortedRegions "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF3_Steinmetz_kmeans-6.bed" \
                --outFileNameMatrix "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF3_Steinmetz_kmeans-6.mat.gz" \
                --samplesLabel Esa1 Rpd3 Gcn5 \
                --plotTitle "HATs_HDAC_merged_SF3_reps_Q" \
                --dpi 600 \
                --kmeans 6


            ### 18 ###
            plotHeatmap \
                -m "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF3_Steinmetz.gz" \
                -out "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF3_Steinmetz_kmeans-7.png" \
                --outFileSortedRegions "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF3_Steinmetz_kmeans-7.bed" \
                --outFileNameMatrix "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_merged_SF3_Steinmetz_kmeans-7.mat.gz" \
                --samplesLabel Esa1 Rpd3 Gcn5 \
                --plotTitle "HATs_HDAC_merged_SF3_reps_Q" \
                --dpi 600 \
                --kmeans 7


❯ #  ----------------------------------------------------------------------------
❯ ### 6 November 2023 - trying plotCorrelation with merged files ###

❯ #  Compare all MAPQ2 files to each other
❯ multiBigwigSummary bins \
>     -b "${d_bws}/5781_IP.bamCompare-to-input.SF2-coef-0.109063.bw" \
>        "${d_bws}/7041_IP.bamCompare-to-input.SF2-coef-0.837016.bw" \
>        "${d_bws}/7691_IP.bamCompare-to-input.SF2-coef-0.949774.bw" \
>        "${d_bws}/7568_IP.bamCompare-to-input.SF2-coef-0.649534.bw" \
>        "${d_bws}/7569_IP.bamCompare-to-input.SF2-coef-1.310394.bw" \
>        "${d_bws}/7692_IP.bamCompare-to-input.SF2-coef-0.616505.bw" \
>        "${d_bws}/7709_IP.bamCompare-to-input.SF2-coef-0.831745.bw" \
>     --labels untagged 7041_Esa1 7691_Esa1 7568_Rpd3 7569_Rpd3 7692_Gcn5 7709_Gcn5 \
>     -o "${d_plots}/all_MAPQ2file_compare.npz" \
>     -p "${threads}"
Number of bins found: 2490


❯ #  Compare all merged (and untagged control) MAPQ2 files to each other
❯ multiBigwigSummary bins \
>     -b "${d_bws}/5781_IP.bamCompare-to-input.SF2-coef-0.109063.bw" \
>        "${d_bws}/Esa1_merged_IP.bamCompare-to-input.SF2-coef-0.898049.bw" \
>        "${d_bws}/Rpd3_merged_IP.bamCompare-to-input.SF2-coef-1.021932.bw" \
>        "${d_bws}/Gcn5_merged_IP.bamCompare-to-input.SF2-coef-0.699061.bw" \
>     --labels untagged Esa1 Rpd3 Gcn5 \
>     -o "${d_plots}/merged_untagged_MAPQ2_compare.npz" \
>     -p "${threads}"
Number of bins found: 2490


❯ #  Compare all merged MAPQ2 files to each other (no untagged control)
❯ multiBigwigSummary bins \
>     -b "${d_bws}/Esa1_merged_IP.bamCompare-to-input.SF2-coef-0.898049.bw" \
>        "${d_bws}/Rpd3_merged_IP.bamCompare-to-input.SF2-coef-1.021932.bw" \
>        "${d_bws}/Gcn5_merged_IP.bamCompare-to-input.SF2-coef-0.699061.bw" \
>     --labels Esa1 Rpd3 Gcn5 \
>     -o "${d_plots}/merged_MAPQ2_compare.npz" \
>     -p "${threads}"
Number of bins found: 2490


❯ #  Draw plots
❯ print_test=true

❯ run_command=true
❯ iter=0
❯ unset stems && typeset -a stems=(
>     "all_MAPQ2file_compare"
>     "merged_untagged_MAPQ2_compare"
>     "merged_MAPQ2_compare"
> )


❯ for cor in "pearson" "spearman"; do
>     for what in "scatterplot" "heatmap"; do
>         for stem in "${stems[@]}"; do
>             (( iter++ ))
> 
>             if ${print_test}; then
>                 echo """
>                 ### ${iter} ###
>                 plotCorrelation \\
>                     --corData \"${d_plots}/${stem}.npz\" \\
>                     --corMethod \"${cor}\" \\
>                     --whatToPlot \"${what}\" \\
>                     --plotFile \"${d_plots}/${stem}_${cor}_${what}.png\"
>                 """
>             fi
> 
>             if ${run_command}; then
>                 plotCorrelation \
>                     --corData "${d_plots}/${stem}.npz" \
>                     --corMethod "${cor}" \
>                     --whatToPlot "${what}" \
>                     --plotFile "${d_plots}/${stem}_${cor}_${what}.png"
>             fi
>         done
>     done
> done

                ### 1 ###
                plotCorrelation \
                    --corData "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_MAPQ2file_compare.npz" \
                    --corMethod "pearson" \
                    --whatToPlot "scatterplot" \
                    --plotFile "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_MAPQ2file_compare_pearson_scatterplot.png"


                ### 2 ###
                plotCorrelation \
                    --corData "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/merged_untagged_MAPQ2_compare.npz" \
                    --corMethod "pearson" \
                    --whatToPlot "scatterplot" \
                    --plotFile "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/merged_untagged_MAPQ2_compare_pearson_scatterplot.png"


                ### 3 ###
                plotCorrelation \
                    --corData "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/merged_MAPQ2_compare.npz" \
                    --corMethod "pearson" \
                    --whatToPlot "scatterplot" \
                    --plotFile "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/merged_MAPQ2_compare_pearson_scatterplot.png"


                ### 4 ###
                plotCorrelation \
                    --corData "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_MAPQ2file_compare.npz" \
                    --corMethod "pearson" \
                    --whatToPlot "heatmap" \
                    --plotFile "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_MAPQ2file_compare_pearson_heatmap.png"


                ### 5 ###
                plotCorrelation \
                    --corData "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/merged_untagged_MAPQ2_compare.npz" \
                    --corMethod "pearson" \
                    --whatToPlot "heatmap" \
                    --plotFile "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/merged_untagged_MAPQ2_compare_pearson_heatmap.png"


                ### 6 ###
                plotCorrelation \
                    --corData "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/merged_MAPQ2_compare.npz" \
                    --corMethod "pearson" \
                    --whatToPlot "heatmap" \
                    --plotFile "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/merged_MAPQ2_compare_pearson_heatmap.png"


                ### 7 ###
                plotCorrelation \
                    --corData "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_MAPQ2file_compare.npz" \
                    --corMethod "spearman" \
                    --whatToPlot "scatterplot" \
                    --plotFile "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_MAPQ2file_compare_spearman_scatterplot.png"


                ### 8 ###
                plotCorrelation \
                    --corData "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/merged_untagged_MAPQ2_compare.npz" \
                    --corMethod "spearman" \
                    --whatToPlot "scatterplot" \
                    --plotFile "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/merged_untagged_MAPQ2_compare_spearman_scatterplot.png"


                ### 9 ###
                plotCorrelation \
                    --corData "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/merged_MAPQ2_compare.npz" \
                    --corMethod "spearman" \
                    --whatToPlot "scatterplot" \
                    --plotFile "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/merged_MAPQ2_compare_spearman_scatterplot.png"


                ### 10 ###
                plotCorrelation \
                    --corData "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_MAPQ2file_compare.npz" \
                    --corMethod "spearman" \
                    --whatToPlot "heatmap" \
                    --plotFile "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/all_MAPQ2file_compare_spearman_heatmap.png"


                ### 11 ###
                plotCorrelation \
                    --corData "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/merged_untagged_MAPQ2_compare.npz" \
                    --corMethod "spearman" \
                    --whatToPlot "heatmap" \
                    --plotFile "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/merged_untagged_MAPQ2_compare_spearman_heatmap.png"


                ### 12 ###
                plotCorrelation \
                    --corData "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/merged_MAPQ2_compare.npz" \
                    --corMethod "spearman" \
                    --whatToPlot "heatmap" \
                    --plotFile "/home/kalavatt/2023_rDNA/results/2023-0406_tutorial_ChIP-seq_analyses/test_RD-data-viz/merged_MAPQ2_compare_spearman_heatmap.png"


❯ #  Compress the directory in which all the outfiles are stored ================
❯ if [[ -d "${d_plots}" && ! -f "${d_plots}.tar.gz" ]]; then
>     tar czf "${d_plots}.tar.gz" "${d_plots}"
> fi
tar: Removing leading `/' from member names


❯ if [[ -d "${d_plots}" && -f "${d_plots}.tar.gz" ]]; then
>     rm -rf "${d_plots}"
> fi
```
</details>
<br />
