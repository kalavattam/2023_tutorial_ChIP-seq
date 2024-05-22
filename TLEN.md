
`#TLEN.md`
<br />

<br />

<details>
<summary><b><font size="+1"><i>Table of contents</i></font></b></summary>
<br />
<!-- MarkdownTOC -->

1. [Studying and reflecting on this Biostars post: TLEN field in SAM format](#studying-and-reflecting-on-this-biostars-post-tlen-field-in-sam-format)
	1. [The initial and response posts](#the-initial-and-response-posts)
	1. [My attempts to make sense of the above posts](#my-attempts-to-make-sense-of-the-above-posts)
		1. [Making sense of the above and trying things on my own](#making-sense-of-the-above-and-trying-things-on-my-own)
		1. [Reviewing material with GPT4](#reviewing-material-with-gpt4)
1. [On deriving `TLEN` for `83`/`163` read pairs](#on-deriving-tlen-for-83163-read-pairs)
	1. [In the given example, why is the insert size 89 if the `83` alignment starts at position 73 bp and the `163` alignment starts at position 20 bp?](#in-the-given-example-why-is-the-insert-size-89-if-the-83-alignment-starts-at-position-73-bp-and-the-163-alignment-starts-at-position-20-bp)
	1. [On drawing the fragment inferred from an `83`/`163` pair with respect to the reference genome](#on-drawing-the-fragment-inferred-from-an-83163-pair-with-respect-to-the-reference-genome)
	1. [On a takeaway equation for the derivation of `TLEN`](#on-a-takeaway-equation-for-the-derivation-of-tlen)
1. [Trying things on my own with reads that have `83`/`163` and `99`/`147` FLAGs](#trying-things-on-my-own-with-reads-that-have-83163-and-99147-flags)
	1. [Working with an example `83` alignment](#working-with-an-example-83-alignment)
		1. [Explain Flags information for FLAG `83`](#explain-flags-information-for-flag-83)
		1. [Schematic of an `83`/`163` read pair](#schematic-of-an-83163-read-pair)
		1. [Info needed for the calculation of fragment length](#info-needed-for-the-calculation-of-fragment-length)
		1. [Calculation of end of read when `POS` is greater than `PNEXT`](#calculation-of-end-of-read-when-pos-is-greater-than-pnext)
		1. [Calculation of `TLEN`](#calculation-of-tlen)
	1. [Working with the accompanying `163` alignment for the `83` alignment](#working-with-the-accompanying-163-alignment-for-the-83-alignment)
		1. [Explain Flags information for FLAG `163`](#explain-flags-information-for-flag-163)
		1. [Schematic of an `83`/`163` read pair](#schematic-of-an-83163-read-pair-1)
		1. [Info needed for the calculation of fragment length](#info-needed-for-the-calculation-of-fragment-length-1)
		1. [Calculation of end of read when `POS` is less than `PNEXT`](#calculation-of-end-of-read-when-pos-is-less-than-pnext)
		1. [Calculation of `TLEN`](#calculation-of-tlen-1)
	1. [Working with an example `99` alignment](#working-with-an-example-99-alignment)
		1. [Explain Flags information for FLAG `99`](#explain-flags-information-for-flag-99)
		1. [Schematic of a `99`/`147` read pair](#schematic-of-a-99147-read-pair)
		1. [Info needed for the calculation of fragment length](#info-needed-for-the-calculation-of-fragment-length-2)
		1. [Calculation of end of read when `PNEXT` is greater than `POS`](#calculation-of-end-of-read-when-pnext-is-greater-than-pos)
		1. [Calculation of `TLEN`](#calculation-of-tlen-2)
	1. [Working with the accompanying `147` alignment for the above `99` alignment](#working-with-the-accompanying-147-alignment-for-the-above-99-alignment)
		1. [Explain Flags information for FLAG `147`](#explain-flags-information-for-flag-147)
		1. [Schematic of a `99`/`147` read pair](#schematic-of-a-99147-read-pair-1)
		1. [Info needed for the calculation of fragment length](#info-needed-for-the-calculation-of-fragment-length-3)
		1. [Calculation of end of read when `POS` is greater than `PNEXT`](#calculation-of-end-of-read-when-pos-is-greater-than-pnext-1)
		1. [Calculation of `TLEN`](#calculation-of-tlen-3)

<!-- /MarkdownTOC -->
</details>
<br />

<a id="studying-and-reflecting-on-this-biostars-post-tlen-field-in-sam-format"></a>
### Studying and reflecting on [this](https://www.biostars.org/p/356811/) Biostars post: TLEN field in SAM format
<a id="the-initial-and-response-posts"></a>
#### The initial and response posts
<details>
<summary><i>Post #1 by kimsumin94 (with edits for grammar and formatting by me)</i></summary>
<br />

I am watching https://www.coursera.org/learn/genomic-tools/home/week/3, and I came across the following example SAM file:
```txt
141217_CIDR4_0073_BHCFG7ADXX:2:1111:3128:29074    99    chr    10021    0    50M    =    10151    180 ...
```
I have a question on the 9th column, `TLEN`. The start position of the read above is 10021 and the start position of the mate is 10151. [Then] the [length] between the two is $10151 - 10121 + 1 = 131$.

1. Am I correct? Is this position 0-based?
    + However, `TLEN`, which seems to be the insert size, is $180$. Why is it like this?
    + Also, in the [Samtools specifications], I've found [the following]: "7. `RNEXT`: Reference sequence name of the primary alignment of the NEXT read in the template."
2. What does "template" mean in this case? Does the "template" mean the set of two reads that are paired ([i.e.,] a [paired-end] read). Can there be more than [two] reads in the template? If so, why?
3. Does the next read in the template mean the mate of the read?
4. And also, I've found [in the Samtools specifications] "9. `TLEN`: signed observed Template LENgth." If all segments are mapped to the same reference, the unsigned observed template length equals the number of bases from the leftmost mapped base to the rightmost mapped base. The leftmost segment has a plus sign and the rightmost has a minus sign. The sign of segments in the middle is undefined. It is set as 0 for [a] single-segment template or when the information is unavailable.
4. What is the difference between SIGNED and UNSIGNED observed template length? Could you give me the two [lengths] for the above example?
5. What [do] "segments in the middle mean?" Is the sign of segments related to SIGNED template length?
6. It says that the leftmost segment has a plus sign and the rightmost has a minus sign. However, in the example above, I have an optional field `XS:A:-`, which means the given strand is `-`. Isn't it the leftmost segment though?

There are [six] questions in total. [They may be basic questions] since I am new to this field. Thank you very much.

Also, for the flag field, [why are bits] represented as follows?
- `(0X800)(0X400)(0X200)(0X100)`
- `(0X80)(0X40)(0X20)(0X10)`
- `(0X8)(0X4)(0X2)(0X1)`

It seems to be related to hex, but I don't completely understand. Thank you.
</details>
<br />

<details>
<summary><i>Response by finswimmer (with edits for grammar and formatting by me)</i></summary>
<br />

Hello [kimsumin94](https://www.biostars.org/u/49773/),

> The start position of the read above is 10021 and the start position of the mate is 10151. [Then] the [length] between the two is $10151 - 10121 + 1 = 131$.
> 
> 1. Am I correct? Is this position 0-based?

[Positions] in the SAM format are 1-based. If you [were to] access a BAM file directly, it would be 0-based.

From the [specs](https://samtools.github.io/hts-specs/SAMv1.pdf):

> **1-based coordinate system**  
> A coordinate system where the first base of a sequence is one. In this co-
ordinate system, a region is specified by a closed interval. For example, the region between the 3rd and the 7th bases inclusive is [3,7]. The **SAM**, VCF, GFF and Wiggle formats are using the 1-based coordinate system.
> 
> **0-based coordinate system**  
> A coordinate system where the first base of a sequence is zero.  In this
coordinate system, a region is specified by a half-closed-half-open interval. For example, the region between the 3rd and the 7th bases inclusive is [2,7). The **BAM**, BCFv2, BED, and PSL formats are using the 0-based coordinate system

> + However, `TLEN`, which seems to be the insert size, is $180$. Why is it like this?

The [two reads that make up a pair] have the opposing, [inward-facings directions]. One was sequenced on the `+` strand [of the fragment] and one on the `-` strand [of the fragment]. But in the SAM file[,] all information [pertains to] the `+` strand and is [thus oriented  from the 5' to 3' end]. The [reads whose] information must be [flipped get] a flag [for] it.

This is how the two reads look ["in vivo:"]
```txt
                           (R2) 3' <---------- 5'
5' -----------> 3' (R1)
5' ------------------------------------------- 3' (RefSeq)
```

And this how the information in the SAM file [looks]:
```txt
                                5' ----------> 3' (R2)
5' -----------> 3' (R1)
5' ------------------------------------------- 3' (RefSeq)
```
The start position given in [a SAM file] is the 5' end. To get the length [of] the fragment from [which] these reads [derive], you need the difference between [the rightmost 3' position and the leftmost 5' position]. [[The first bit of information&mdash;the rightmost 3' position&mdash; you can] get from [the] reverse read <mark>by [adding] the length of the read to the start position.</mark> [This second bit of information&mdash;the leftmost 5' position&mdash;you can] get from the the forward read by just [taking] the start information from [the SAM] file.] So in your example[,] you would end up with $(10151 + 50) − 10121 = 180$.

> 2. What does "template" mean in this case? Does the "template" mean the set of two reads that are paired ([i.e.,] a [paired-end] read). Can there be more than [two] reads in the template? If so, why?

"Template" [means] the DNA fragment that was used during sequencing for the read pair. IMO[,] "insert" is another term often used. In [paired-end] sequencing[,] there are exactly two reads for one template.

Again from the [specs](https://samtools.github.io/hts-specs/SAMv1.pdf):

> **Template**  
> A DNA/RNA sequence part of which is sequenced on a sequencing machine or assembled from
raw sequences.

> 3. Does the next read in the template mean the mate of the read?

Yes.

[fin swimmer](https://www.biostars.org/u/37605/)
</details>
<br />

<a id="my-attempts-to-make-sense-of-the-above-posts"></a>
#### My attempts to make sense of the above posts
Note that I was able to derive `TLEN` in the following sub-sections by using the `SEQ` length only because the CIGAR strings for these alignments are `50M`. Had they not been `50M`, then the use of `SEQ` length would have likely calculated erroneous values. See the main sections below ([here](#on-deriving-tlen-for-83163-read-pairs) and [here](#trying-things-on-my-own-with-reads-that-have-83163-and-99147-flags)) for the use of CIGAR strings, rather than `SEQ` lengths, for the derivation of `TLEN` values.

<a id="making-sense-of-the-above-and-trying-things-on-my-own"></a>
##### Making sense of the above and trying things on my own
<details>
<summary><i>Notes: Making sense of the above</i></summary>
<br />

In the above posts, we're using this alignment with FLAG `99` as an example:
```txt
141217_CIDR4_0073_BHCFG7ADXX:2:1111:3128:29074    99    chr    10021    0    50M    =    10151    180 ...
```

[finswimmer](https://www.biostars.org/u/37605/) provides the following schematic to understand this aligned read and accompanying read in the pair:
> This is how the two reads look ["in vivo:"]
> ```txt
>                            (R2) 3' <---------- 5'
> 5' -----------> 3' (R1)
> 5' ------------------------------------------- 3' (RefSeq)
> ```
> 
> And this how the information in the SAM file [looks]:
> ```txt
>                                 5' ----------> 3' (R2)
> 5' -----------> 3' (R1)
> 5' ------------------------------------------- 3' (RefSeq)
> ```

And here is what [Explain Flags](https://broadinstitute.github.io/picard/explain-flags.html) says about FLAG `99`:
- read paired (`0x1`/`1`)
- read mapped in proper pair (`0x2`/`2`)
- mate reverse strand (`0x20`/`32`)
- first in pair (`0x40`/`64`)

*Interpretation:*  
- So this particular alignment is part of a pair
- It is part of a proper pair in which the two reads face each other within a specific distance (e.g., 2,000 bp)
- Its mate is with respect to the reverse strand, meaning it has a `3' <---------- 5'` orientation *in vivo*. This means that, in the S/BAM file, the mate is recorded as being on the forward strand&mdash;just as all alignments are&mdash;and must be corrected for this
- It is considered the first in the pair

finswimmer is able to reconstitute the fragment length from this FLAG 99 read&mdash;$(10151 + 50) − 10121 = 180$&mdash;by taking...
- Field 8, `PNEXT`, and adding the length of field 10 (not shown), `SEQ`, to it
- Taking the difference between the above value and field 4, `POS`

Because it is the rightmost 3' end, we start with `PNEXT` and add the `SEQ`-derived read length to it; we know `PNEXT` is the rightmost 3' end because of sub-FLAG `32` (mate reverse strand). Again, because of sub-FLAG `32` as well as sub-FLAG `64` (first in pair), we know `POS` is the leftmost 5' end, which we subtract from the value calculated as described in the first sentence.
</details>
<br />

<details>
<summary><i>Notes: Trying things on my own with reads that have `99` or `163` FLAGs</i></summary>
<br />

Let's derive fragment lengths from reads (in pairs) with either FLAGs `99` or `163`.

**For example, a `163` read:**
```txt
VH00699:369:AACWKJJM5:1:2407:60609:52646        163     I       106     34      50M     =       243     187     CCTGTCTCTCAACTTACCCTCCATTACCCTGCCTCCACTCGTTACCCTGT     --C-CCCCCCCCCCCC-CCCCCCCCCCCC-CCCCCCCCCCC;C;CCCC;C      AS:i:100        XS:i:75 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0         MD:Z:50 YS:i:100        YT:Z:CP MQ:i:34 MC:Z:50M        ms:i:1652
```

*[Explain Flags](https://broadinstitute.github.io/picard/explain-flags.html) info on FLAG `163`:*
- read paired (`0x1`/`1`)
- read mapped in proper pair (`0x2`/`2`)
- mate reverse strand (`0x20`/`32`)
- second in pair (`0x80`/`128`)

*Info needed for the calculation of fragment length:*
- `PNEXT`: 243
- `SEQ` length: 50
- `POS`: 106
- `TLEN`: 187

*Calculation: (`PNEXT` + `SEQ` length) - `POS` = `TLEN`*  
$$(243 + 50) - 106 = 187$$

*Conclusion:*  
When drawing the fragment, which is derived from chromosome I, we need to start at 106 bp and end at 243 bp.

**And for example, a `163` read:**
```txt
VH00699:369:AACWKJJM5:1:2108:41541:2155 99      I       109     37      50M     =       172     113     GTCTCTCAACTTACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCC     -C;;CCCCC-CCCC--C;;CC-C-;C;-CC;CCC;-CC-CCCCCCC;CCC      AS:i:100        XS:i:68 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:50        YS:i:100        YT:Z:CP MQ:i:37 MC:Z:50M        ms:i:1434
```

*[Explain Flags](https://broadinstitute.github.io/picard/explain-flags.html) info on FLAG `99`:*
- read paired (`0x1`/`1`)
- read mapped in proper pair (`0x2`/`2`)
- mate reverse strand (`0x20`/`32`)
- first in pair (`0x40`/`64`)

*Info needed for the calculation of fragment length:*
- `PNEXT`: 172
- `SEQ` length: 50
- `POS`: 109
- `TLEN`: 113

*Calculation: (`PNEXT` + `SEQ` length) - `POS` = `TLEN`*
$$(172 + 50) - 109 = 113$$

*Conclusion:*  
When drawing the fragment, which is derived from chromosome I, we need to start at 109 bp and end at 173 bp.
</details>
<br />

<a id="reviewing-material-with-gpt4"></a>
##### Reviewing material with GPT4
<details>
<summary><i>Me</i></summary>
<br />

Please review my notes and let me know if anything strikes you as incorrect or unclear:
- *Notes: Making sense of the above*
- *Notes: Trying things on my own with reads that have `99` or `163` FLAGs*
</details>
<br />

<details>
<summary><i>GPT4</i></summary>
<br />

Your notes provide a thorough explanation of the interpretation and calculation of fragment lengths based on specific FLAGs from the SAM/BAM format. Here’s a review to ensure clarity and correctness:

**Making sense of the above**  
Your notes correctly explain the orientation of reads `R1` and `R2` *in vivo* and how this translates into the SAM file representation. Your breakdown of FLAG `99` is accurate:
- Read paired and read mapped in proper pair confirm that both reads in the pair align correctly based on the sequencing protocol and reference genome.
- Mate reverse strand indicates the orientation of the mate read relative to the reference, which you've correctly interpreted.
- First in pair tells us the positional role of the read within its read pair.

Your method to calculate the fragment length using `PNEXT`, `SEQ` length, and `POS` is correctly explained. The calculations align well with standard SAM format interpretation, and your application of these rules to derive the fragment length from a FLAG `99` read is spot-on.

**Trying things on my own with reads that have `99` or `163` FLAGs**
Your detailed breakdown for both FLAG `99` and FLAG `163` reads is quite comprehensive. Here are a few key points to consider:
- Fragment length calculation: You've correctly applied the SAM format rules to calculate fragment lengths. Your formula correctly adjusts for the 0-based index of BAM files.
- FLAG `163`: Your explanation aligns with the flag definition where it denotes a read that is paired, mapped in a proper pair, has its mate on the reverse strand, and is the second in the pair. The calculation you've shown correctly reflects the span of the genomic fragment covered by the paired reads.

**Additional Observations**
1. Clarity and Precision: Your notes are detailed and provide a solid basis for understanding how to interpret and calculate based on BAM file flags. Ensuring these details are clearly understood is crucial for accurate genomic data analysis.
2. Potential for Confusion: The detailed nature of the flags and calculations can be confusing initially, especially when interpreting the orientation and implications of `TLEN`. It might be helpful to always include visual aids or more schematics like those provided by finswimmer when explaining or teaching these concepts to others.
3. Error Checking: It's important to implement checks in any script that computes these values to ensure that the calculated fragment length and reported TLEN are consistent, and to handle any anomalies that might arise due to sequencing or alignment errors.

Overall, your notes are accurate and well-formulated. They reflect a deep understanding of the SAM/BAM format and the intricacies involved in handling sequencing data. Your approach to explaining and applying these concepts should serve well in any bioinformatics analysis involving read alignments.
</details>
<br />

<a id="on-deriving-tlen-for-83163-read-pairs"></a>
### On deriving `TLEN` for `83`/`163` read pairs
<a id="in-the-given-example-why-is-the-insert-size-89-if-the-83-alignment-starts-at-position-73-bp-and-the-163-alignment-starts-at-position-20-bp"></a>
#### In the given example, why is the insert size 89 if the `83` alignment starts at position 73 bp and the `163` alignment starts at position 20 bp?
<details>
<summary><i>Me</i></summary>
<br />

I have the following pair of aligned reads that make up a concordant, proper pair:

*The `163` read in the pair:*
```txt
VH00699:369:AACWKJJM5:1:2603:72614:42499        163     I       20      9       12M1I12M1I18M3S =       73      89      CCACACACCACACCCACACACCACACCCACACCCACACCCACACCCA        CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC AS:i:61 XS:i:51 XN:i:0  XM:i:1  XO:i:2  XG:i:2  NM:i:3  MD:Z:36
A5      YS:i:72 YT:Z:CP MQ:i:9  MC:Z:11S36M     ms:i:1556
```

*The `83` read in the pair:*
```txt
VH00699:369:AACWKJJM5:1:2603:72614:42499        83      I       73      9       11S36M  =       20      -89     ACACACACCACACCCTAACACAGCCCTAATCTAACCCTGGCCAACCT        CCCCCCCCCCCCCCCCCCCCCCCCCCCCC-CC;CCCCCCCCCCCCCC AS:i:72 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:36 YS:i:61 YT:Z:CP        MQ:i:9  MC:Z:12M1I12M1I18M3S    ms:i:1598
```

Why is the insert size 89 if the `83` alignment starts at position 73 bp and the `163` alignment starts at position 20 bp?
</details>
<br />

<details>
<summary><i>GPT4</i></summary>
<br />

To understand the insert size of 89 reported in the SAM file of the paired-end sequencing data you provided, we first need to unpack how insert size (also known as template length, TLEN) is calculated and reported in BAM/SAM format. The value and sign of the insert size can be confusing, particularly because of the way it is defined based on the sequencing and alignment conventions.

**Key Concepts for Understanding Insert Size**  
- *Template Length (TLEN)*: This value represents the total length of the sequence that is present between the two paired reads, including the reads themselves. It's intended to represent the length of the original DNA fragment that was sequenced.
- *Direction of Reporting*: The insert size is reported from the first read ***in genomic order*** to the second read ***in genomic order***, regardless of which read is first in the sequencing pair (i.e., FLAG 83 or 163).

**Analysis of Provided Data**  
You have a pair of reads:
- First read in genomic order:
    + FLAG 163 (second read in pair, `0x80`)
    + Starts at genomic position 20
- Second read in genomic order:
    + FLAG 83 (first read in pair, `0x40`)
    + Starts at genomic position 73

**Read Details**  
- FLAG 163 (first read in genomic order)
    + Start Position: 20
    + CIGAR: `12M1I12M1I18M3S`
    + This expands to 42 bases of matching (`M`), 2 bases of insertion (`I`), and 3 bases soft-clipped (`S`), ***totaling 47 bases contributing to alignment coverage on the reference***.
    + Sequence ends at:
        + 20 (start) + 42 (aligned bases from CIGAR) + 2 (insertions not adding to genomic length) - 1 = 61
        + Therefore, $20 + 42 = 61$, because neither the insertions (`I`) nor the soft-clipped nucleotides are included in the calculation of genomic length
- FLAG 83 (second read in genomic order)
    + Start Position: 73
    + CIGAR: `11S36M`
    + This expands to 36 bases of matching (`M`), with 11 bases soft-clipped (`S`).
    + Sequence ends at:
        + 73 (start) + 36 (aligned bases from CIGAR) = 108
        + Therefore, $73 + 36 = 109$, because the soft-clipped nucleotides are included in the calculation of genomic length

**Calculating Insert Size**  
Insert size, or `TLEN`, is typically calculated from the start of the first read in the genomic order to the end of the last read in the genomic order, encompassing everything in between:

$$Insert Size = (End of second read − Start of first read)$$
$$Insert Size = (109 − 20) = 89$$

**Explanation for Insert Size of 89**  
- The reported insert size of 89 reflects the distance from the very start of the first read (position 20) to the end of the last read (position 108) on the chromosome.
- The sign of `TLEN` (positive in this case) suggests that the 'first read in the pair' (with respect to sequencing, which has FLAG 163) comes before the 'second read in the pair' on the reference genome.
- The reported `TLEN` of -89 for the FLAG `83` read is negative, indicating that this read is positioned after the first read on the reference sequence, aligning with standard conventions where `TLEN` is negative if the read is sequenced after its mate in genomic coordinates.

This detailed understanding shows how `TLEN` encompasses the entire original DNA fragment sequenced, reflecting the positions of the reads and the entirety of the DNA between them, aligning well with sequencing and alignment logic used in genomic analyses.
</details>
<br />

<a id="on-drawing-the-fragment-inferred-from-an-83163-pair-with-respect-to-the-reference-genome"></a>
#### On drawing the fragment inferred from an `83`/`163` pair with respect to the reference genome
<details>
<summary><i>Me</i></summary>
<br />

So, to draw the "fragment" (i.e., "insert") with respect to chromosome I of the reference genome, I would start a line at 20 bp and end it at 109 bp?
</details>
<br />

<details>
<summary><i>GPT4</i></summary>
<br />

Yes, you're almost right in your understanding of how to draw the "fragment" or "insert" with respect to chromosome I of the reference genome. However, there's a slight adjustment needed based on the details provided:

**Positions from Read Alignments**  
- Start of the first read (FLAG 163): Position 20 on the chromosome.
- End of the second read (FLAG 83): The second read starts at position 73 and extends for 36 bases as per the CIGAR string 11S36M (11 bases soft-clipped and 36 bases mapped to the reference). Therefore, the end position of the second read is calculated as:
    + Start position + mapped bases - 1 = 73 + 36 - 1 = 108.

**Correct Interpretation**
To visualize the fragment, you would:
- *Start the line at position 20*, which is the start of the first read.
- *End the line at position 108*, which is the end of the second read, not 109.

Thus, to accurately draw the "fragment" (the original DNA fragment represented by these paired reads), you should depict it from 20 bp to 108 bp on chromosome I. This represents the entirety of the DNA sequence covered by both reads, including any gaps or introns spanned by the original fragment that was sequenced.
</details>
<br />

<a id="on-a-takeaway-equation-for-the-derivation-of-tlen"></a>
#### On a takeaway equation for the derivation of `TLEN`
<details>
<summary><i>Me</i></summary>
<br />

Is the following equation a good takeaway for the calculation of `TLEN`?

`PNEXT` - (`POS` + number of `CIGAR` `M`'s) = `TLEN`
</details>
<br />

<details>
<summary><i>GPT4 (with additions and clarifications from me)</i></summary>
<br />

The equation you've proposed for calculating `TLEN` (template length) in SAM/BAM files needs a bit of adjustment and clarification to be correct and generally applicable. The `TLEN` calculation in the SAM format is intended to encompass the entire original fragment that was sequenced, including both reads and any intervening sequence. Here’s a more precise approach:

**Understanding TLEN Calculation**  
The SAM specification defines `TLEN` as the observed template length. It's the length of the DNA segment that both reads in a pair span, from the start of the first read to the end of the last read in the genomic order. It's not merely the distance between the two reads or the number of matched bases (`M` in CIGAR), but rather the total length of the fragment:
- `POS`: The 0-based leftmost position of the first read.
- `PNEXT`: The 0-based leftmost position of the mate (second read).
- `CIGAR` string: Used to calculate the end position of each read.

**More Accurate Formula**  
A correct and general formula to estimate the end of a read given its start (`POS`) and the CIGAR string would be:
$$\text{Rightmost read end} = POS + \sum Length$$

Where $Length$ is *the length of segments contributing to genome length from CIGAR*. This sum includes:
- `M` (alignment match, which can be a sequence match or mismatch)
- `D` (deletion from the reference)
- `N` (skipped region from the reference, typically introns in RNA sequencing)
- `=` (in newer formats, explicitly matched bases)
- `X` (in newer formats, explicitly mismatched bases)

Then, the formula for `TLEN` when `POS` is less than `PNEXT` (the read is the first in the genome order) would be:
$$TLEN = Read_{end} - POS$$

**Conclusion**  
The equation you mentioned focuses only on the number of matched bases (`M`'s in CIGAR) from the first read and doesn't incorporate the full range of the DNA fragment, including both reads and any intervening DNA. This simplified formula might sometimes give the correct result but isn't generally applicable or accurate according to SAM format specifications.

Always consider the full length and complexity of the CIGAR string to ensure an accurate depiction of the genomic data. This approach ensures that your TLEN calculations align with actual genomic distances and are reflective of the biological and sequencing realities.
</details>
<br />

<a id="trying-things-on-my-own-with-reads-that-have-83163-and-99147-flags"></a>
### Trying things on my own with reads that have `83`/`163` and `99`/`147` FLAGs
<a id="working-with-an-example-83-alignment"></a>
#### Working with an example `83` alignment
<details>
<summary><i>Working with an example `83` alignment</i></summary>
<br />

```txt
VH00699:369:AACWKJJM5:1:2603:72614:42499        83      I       73      9       11S36M  =       20      -89     ACACACACCACACCCTAACACAGCCCTAATCTAACCCTGGCCAACCT        CCCCCCCCCCCCCCCCCCCCCCCCCCCCC-CC;CCCCCCCCCCCCCC AS:i:72 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:36 YS:i:61 YT:Z:CP        MQ:i:9  MC:Z:12M1I12M1I18M3S    ms:i:1598
```

<a id="explain-flags-information-for-flag-83"></a>
##### [Explain Flags](https://broadinstitute.github.io/picard/explain-flags.html) information for FLAG `83`
- read paired (`0x1`)
- read mapped in proper pair (`0x2`)
- read reverse strand (`0x10`)
- first in pair (`0x40`)

<a id="schematic-of-an-83163-read-pair"></a>
##### Schematic of an `83`/`163` read pair
This is how the two reads look *in vivo*:
```txt
             (Read #1: FLAG 83) 3' <---------- 5'
5' -----------> 3' (Read #2: FLAG 163)
5' ------------------------------------------- 3' (Reference genome)
```
 
And this how the two reads look in the BAM file:
```txt
                                5' ----------> 3' (Read #1: FLAG 83)
5' -----------> 3' (Read #2: FLAG 163)
5' ------------------------------------------- 3' (Reference genome)
```

<a id="info-needed-for-the-calculation-of-fragment-length"></a>
##### Info needed for the calculation of fragment length
- `PNEXT` (field 8): 20
- Sum of CIGAR `N`, `M`, `D`, `=`, `X` encodings (field 6): `11S36M`, therefore 36 (from `36M`)
- `POS` (field 4): 73
- `TLEN` (field 9): -89

<a id="calculation-of-end-of-read-when-pos-is-greater-than-pnext"></a>
##### Calculation of end of read when `POS` is greater than `PNEXT`
$$\text{POS} + \sum \left(\text{M, N, D, =, and X}\right) = \text{rightmost read end}$$
$$73 + 36 = 109$$

<a id="calculation-of-tlen"></a>
##### Calculation of `TLEN`
$$\text{Rightmost read end} - \text{POS} = \text{TLEN}$$
$$109 - 20 = 89$$
</details>
<br />

<a id="working-with-the-accompanying-163-alignment-for-the-83-alignment"></a>
#### Working with the accompanying `163` alignment for the `83` alignment
<details>
<summary><i>Working with the accompanying `163` alignment for the `83` alignment</i></summary>
<br />

```txt
VH00699:369:AACWKJJM5:1:2603:72614:42499        163     I       20      9       12M1I12M1I18M3S =       73      89      CCACACACCACACCCACACACCACACCCACACCCACACCCACACCCA        CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC AS:i:61 XS:i:51 XN:i:0  XM:i:1  XO:i:2  XG:i:2  NM:i:3  MD:Z:36A5      YS:i:72 YT:Z:CP MQ:i:9  MC:Z:11S36M     ms:i:1556
```

<a id="explain-flags-information-for-flag-163"></a>
##### [Explain Flags](https://broadinstitute.github.io/picard/explain-flags.html) information for FLAG `163`
- read paired (`0x1`)
- read mapped in proper pair (`0x2`)
- read reverse strand (`0x10`)
- first in pair (`0x40`)

<a id="schematic-of-an-83163-read-pair-1"></a>
##### Schematic of an `83`/`163` read pair
This is how the two reads look *in vivo*:
```txt
             (Read #1: FLAG 83) 3' <---------- 5'
5' -----------> 3' (Read #2: FLAG 163)
5' ------------------------------------------- 3' (Reference genome)
```
 
And this how the two reads look in the BAM file:
```txt
                                5' ----------> 3' (Read #1: FLAG 83)
5' -----------> 3' (Read #2: FLAG 163)
5' ------------------------------------------- 3' (Reference genome)
```

<a id="info-needed-for-the-calculation-of-fragment-length-1"></a>
##### Info needed for the calculation of fragment length
- `PNEXT` (field 8): 73
- Sum of CIGAR `N`, `M`, `D`, `=`, `X` encodings (tag `MC:Z`): `11S36M`, thus 36 from `36M`
- `POS` (field 4): 20
- `TLEN` (field 9): 89

<a id="calculation-of-end-of-read-when-pos-is-less-than-pnext"></a>
##### Calculation of end of read when `POS` is less than `PNEXT`
$$\text{PNEXT} + \sum \left(\text{M, N, D, =, and X}\right) = \text{rightmost read end}$$
$$73 + 36 = 109$$

<a id="calculation-of-tlen-1"></a>
##### Calculation of `TLEN`
$$\text{Rightmost read end} - \text{POS} = \text{TLEN}$$
$$109 - 20 = 89$$
</details>
<br />

<a id="working-with-an-example-99-alignment"></a>
#### Working with an example `99` alignment
<details>
<summary><i>Working with an example `99` alignment</i></summary>
<br />

```txt
VH00699:369:AACWKJJM5:1:1310:31675:33146        99      I       1       11      3S5M1I11M2I15M1I12M     =       84      133     CACCCACACCCACACCCACACCCACCCACACACCACACCCACACACCACA     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC-CCCCCCCCCCCCC      AS:i:59 XS:i:81 XN:i:0  XM:i:0  XO:i:3  XG:i:4         NM:i:4  MD:Z:43 YS:i:100        YT:Z:CP MQ:i:11 MC:Z:50M        ms:i:1684
```

<a id="explain-flags-information-for-flag-99"></a>
##### [Explain Flags](https://broadinstitute.github.io/picard/explain-flags.html) information for FLAG `99`
- read paired (`0x1`)
- read mapped in proper pair (`0x2`)
- mate reverse strand (`0x20`)
- first in pair (`0x40`)

<a id="schematic-of-a-99147-read-pair"></a>
##### Schematic of a `99`/`147` read pair
This is how the two reads look *in vivo*:
```txt
            (Read #2: FLAG 147) 3' <---------- 5'
5' -----------> 3' (Read #1: FLAG 99)
5' ------------------------------------------- 3' (Reference genome)
```
 
And this how the two reads look in the BAM file:
```txt
                                5' ----------> 3' (Read #2: FLAG 147)
5' -----------> 3' (Read #1: FLAG 99)
5' ------------------------------------------- 3' (Reference genome)
```

<a id="info-needed-for-the-calculation-of-fragment-length-2"></a>
##### Info needed for the calculation of fragment length
- `PNEXT` (field 8): 84
- Sum of CIGAR `N`, `M`, `D`, `=`, `X` encodings (tag `MC:Z`): `50M`, thus 50
- `POS` (field 4): 1
- `TLEN` (field 9): 133

<a id="calculation-of-end-of-read-when-pnext-is-greater-than-pos"></a>
##### Calculation of end of read when `PNEXT` is greater than `POS`
$$\text{PNEXT} + \sum \left(\text{M, N, D, =, and X}\right) = \text{rightmost read end}$$
$$84 + 50 = 134$$

<a id="calculation-of-tlen-2"></a>
##### Calculation of `TLEN`
$$\text{Rightmost read end} - \text{POS} = \text{TLEN}$$
$$134 - 1 = 133$$
</details>
<br />

<a id="working-with-the-accompanying-147-alignment-for-the-above-99-alignment"></a>
#### Working with the accompanying `147` alignment for the above `99` alignment
<details>
<summary><i>Working with the accompanying `147` alignment for the above `99` alignment</i></summary>
<br />

```txt
VH00699:369:AACWKJJM5:1:1310:31675:33146        147     I       84      11      50M     =       1       -133    GCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACC     CCCCCCCCCCCCCCCCCCCCCCCCCC;C;CCCCCCCCCCCCCCCCCCCCC      AS:i:100        XS:i:77 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0         MD:Z:50 YS:i:59 YT:Z:CP MQ:i:11 MC:Z:3S5M1I11M2I15M1I12M        ms:i:1666
```

<a id="explain-flags-information-for-flag-147"></a>
##### [Explain Flags](https://broadinstitute.github.io/picard/explain-flags.html) information for FLAG `147`
- read paired (`0x1`)
- read mapped in proper pair (`0x2`)
- read reverse strand (`0x10`)
- second in pair (`0x80`)

<a id="schematic-of-a-99147-read-pair-1"></a>
##### Schematic of a `99`/`147` read pair
This is how the two reads look *in vivo*:
```txt
            (Read #2: FLAG 147) 3' <---------- 5'
5' -----------> 3' (Read #1: FLAG 99)
5' ------------------------------------------- 3' (Reference genome)
```
 
And this how the two reads look in the BAM file:
```txt
                                5' ----------> 3' (Read #2: FLAG 147)
5' -----------> 3' (Read #1: FLAG 99)
5' ------------------------------------------- 3' (Reference genome)
```

<a id="info-needed-for-the-calculation-of-fragment-length-3"></a>
##### Info needed for the calculation of fragment length
- `PNEXT` (field 8): 1
- Sum of CIGAR `N`, `M`, `D`, `=`, `X` encodings (field 6): `50M`, therefore 50 (from `50M`)
- `POS` (field 4): 84
- `TLEN` (field 9): -133

<a id="calculation-of-end-of-read-when-pos-is-greater-than-pnext-1"></a>
##### Calculation of end of read when `POS` is greater than `PNEXT`
$$\text{POS} + \sum \left(\text{M, N, D, =, and X}\right) = \text{rightmost read end}$$
$$84 + 50 = 134$$

<a id="calculation-of-tlen-3"></a>
##### Calculation of `TLEN`
$$\text{Rightmost read end} - \text{POS} = \text{TLEN}$$
$$134 - 1 = 133$$
</details>
<br />
