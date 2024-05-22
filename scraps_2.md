
`scraps_2.md`
<br />

$$\text{RPKM}_{i} = \frac{R_{i}}{\frac{L_{i}}{10^3} \times \frac{\hat{R}}{10^6}}$$
$$\text{RPKM}_{i} = \frac{R_{i} \times 10^9}{L_{i} \times \hat{R}}$$

where $i$ is a given bin, $R_{i}$ is the number of aligned reads per bin $i$, $L_{i}$ is the bin length in bases, and $\hat{R}$ is the total number of aligned reads.

$$\text{CPM}_{i} = \frac{R_{i}}{\frac{\hat{R}}{10^6}}$$
$$\text{CPM}_{i} = \frac{R_{i}}{\hat{R}} \times 10^6$$

where $i$ is a given bin, $R_{i}$ is the number of aligned reads per bin $i$, and $\hat{R}$ is the total number of aligned reads.

$$\frac{\frac{\hat{R}}{\frac{L_{i}}{10^3}}}{10^6}$$

$$\text{BPM}_{i} = \frac{R_{i}}{\frac{\sum_{j} R_{j}}{10^6}}$$

where $i$ is a given bin, $R_{i}$ is the number of aligned reads per bin $i$, and $\sum_{j} R_{j}$ is the total number of aligned reads across all bins.

$$\text{BPM}_{i} = \frac{R_{i} \times 10^3}{L_{i}} \times \frac{1}{\sum_{j=1}^n \frac{R_{j} \times 10^3}{L_{j}}} \times 10^6$$

$$\frac{1}{\sum_{j=1}^n \frac{R_{j} \times 10^3}{L_{j}}}$$

*My question to GPT4o about the above.*

<details>
<summary><i>Me</i></summary>
<br />

Does this look correct to you?

$$\text{BPM}_{i} = \frac{R_{i} \times 10^3}{L_{i}} \times \frac{1}{\sum_{j=1}^n \frac{R_{j} \times 10^3}{L_{j}}} \times 10^6$$
</details>
<br />

*This answer from GPT4o is wrong.*

<details>
<summary><i>GPT4o</i></summary>
<br />

The formula you provided is close, but it can be simplified. Let's break it down to ensure clarity and correctness.

The expression $\frac{R_{i} \times 10^3}{L_{i}}$ scales the reads per bin by the bin length, and you want to normalize this by the total across all bins. Here is the corrected and simplified formula for BPM (Bins Per Million):
$$\text{BPM}_{i} = \frac{R_{i}}{\sum_{j=1}^n R_{j}}$$

However, it seems you are mixing bin length normalization (which is typically for RPKM) with the BPM concept, which is simpler.

But if you need to incorporate bin length, use RPKM instead, which is:
</details>
<br />

<details>
<summary><i>Me</i></summary>
<br />

BPM is essentially just TPM except gene lengths (which differ from gene to gene) are replaced with bin lengths (which do not differ (except for the ones at the very ends of chromosomes)).
</details>
<br />

<details>
<summary><i>GPT4o</i></summary>
<br />

Understood. In that case, the expression you initially wrote can be slightly refined for clarity and correctness. Here’s the refined formula:

$$\text{BPM}_{i} = \frac{\frac{R_{i} \times 10^3}{L_{i}}}{\sum_{j=1}^n \frac{R_{j} \times 10^3}{L_{j}}} \times 10^6$$

We can simplify the formula by canceling out the $10^3$ factor in the numerator and the denominator:

$$\text{BPM}_{i} = \frac{\frac{R_{i}}{L_{i}}}{\sum_{j=1}^n \frac{R_{j}}{L_{j}}} \times 10^6$$
Here is the complete and simplified notation:

$$\text{BPM}_{i} = \frac{\frac{R_{i}}{L_{i}}}{\sum_{j=1}^n \frac{R_{j}}{L_{j}}} \times 10^6$$

where
- $i$ is a given bin,
- $R_{i}$ is the number of aligned reads per bin 
- $L_{i}$ is the bin length in bases,
- $\sum_{j=1}^n \frac{R_j}{L_j}$ is the sum of normalized read counts for all bins.

This formula reflects the relationship between the read counts, bin lengths, and the total normalized read count, scaled to per million.
</details>
<br />


$$\text{RPGC}_{i} = \frac{R_{i}}{\frac{\hat{R} \times F}{G}}$$, which can be rewritten as $$\text{RPGC}_{i} = \frac{R_{i}}{\hat{R} \times F} \times G$$



In the context of processed ChIP-seq signal data, "coverage" refers to the density of reads, or fragments inferred from reads (see Note \_), aligned to regions of a reference genome. It is useful to frame the discussion of coverage in terms of two separate-but-related elements:
1. A histogram intersecting a genomic range:
    a. Coverage is, in part, a representation of the distribution of reads or fragments across bins, which are consecutive windows of a defined size, starting at one position of a genome and ending at a second, further position.
    b. This is a representation with two dimensions: genomic coordinates on the x-axis and the number of reads or fragments on the y-axis.
2. A scaling factor for normalization:
    a. The calculation of coverage involves a value used to adjust the raw read or fragment counts that compose a histogram in a way that accounts for sequencing depth, library size, or other technical variations. This adjustment can be multiplicative or divisional, which is to say that the histogram y-axis values can be multiplied or divided by the scaling factor.
    b. The primary purpose of applying a scaling factor is to normalize the coverage, making it possible to accurately compare results from within a sample or between different samples (see Section \_ and Note \_).

Please explain $\cap$ versus $\cup$.

$x \cap r \neq \emptyset$
