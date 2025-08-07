# Project 2: Genomic Analysis - Variant Calling in Arabidopsis thaliana

## Overview
This project focuses on developing variant calling pipelines for the plant Arabidopsis thaliana to determine genetic variation. The analysis includes read alignment, variant detection, and statistical analysis of genomic data.

## Table of Contents
- [Docker Setup](#docker-setup)
- [Index Creation](#index-creation)
- [Read Alignment](#read-alignment)
- [Variant Calling](#variant-calling)
- [Analysis Results](#analysis-results)

## Docker Setup

```bash
# Access Docker container with genomic tools
docker run -it --rm -v /Users/labuser3/Downloads/genomic_data_science:/workspace gencommand_image:v1.0.0 /bin/bash
```

## Index Creation

### Question 1: How many index files did the operation create?

Before Bowtie2 can align reads to a reference genome, it needs the reference genome to be in a special, efficient format — that's what the index is.

```bash
bowtie2-build wu_0_A.fasta wu_0
```

**Answer:** 6 index files were created.

## Read Analysis

### Question 2: How many sequences are in the raw sequence files?

```bash
wc -l wu_0_A_wgs.fastq | awk '{print $1/4}'
```

**Answer:** 147,354 sequences

## Read Alignment

### Full Match Alignment

Now that we have the indexed files, we can perform alignment. Note: specify the directory for index files without any extension, include thread number, and convert to SAM format.

```bash
bowtie2 -p 4 -x wu_0 wu_0_A_wgs.fastq -S exome.bt2.sam
```

**Alignment Statistics:**
- Total reads: 147,354
- Unpaired reads: 147,354 (100.00%)
- Aligned 0 times: 9,635 (6.54%)
- Aligned exactly 1 time: 93,780 (63.64%)
- Aligned >1 times: 43,939 (29.82%)
- **Overall alignment rate: 93.46%**

### Local Match Alignment

```bash
bowtie2 -p 4 --local -x wu_0 wu_0_A_wgs.fastq -S exome.local.bt2.sam
```

**Alignment Statistics:**
- Total reads: 147,354
- Unpaired reads: 147,354 (100.00%)
- Aligned 0 times: 6,310 (4.28%)
- Aligned exactly 1 time: 84,939 (57.64%)
- Aligned >1 times: 56,105 (38.07%)
- **Overall alignment rate: 95.72%**

### Question 3: How many reads were mapped in the full match scenario?

**Answer:** 137,719 reads

*Note: This includes reads that aligned exactly once (93,780) plus reads that aligned multiple times (43,939).*

### Question 4: How many reads were mapped in the local match scenario?

**Answer:** 141,044 reads

*Note: This includes reads that aligned exactly once (84,939) plus reads that aligned multiple times (56,105).*

## Variant Calling

### Question 5: How many entries were reported for Chr3?

**Step 1:** Convert SAM to BAM and sort
```bash
samtools view -bT wu_0.v7.fas exome.bt2.sam > exome.bt2.bam
samtools sort exome.bt2.bam exome.bt2.sorted
```

**Step 2:** Generate VCF file using mpileup
```bash
samtools mpileup -f wu_0.v7.fas -uv exome.bt2.sorted.bam > exome.full.mpileup.vcf
```

**Step 3:** Count Chr3 entries
```bash
cut -f1 exome.full.mpileup.vcf | grep "Chr3" | wc -l
```

**Answer:** 360,296 entries

> **⚠️ Important:** It is critical to sort the BAM file before analyzing it with samtools.

### Question 6: How many entries have 'A' as the corresponding genome letter?

**Correct approach:**
```bash
cat exome.full.mpileup.vcf | grep -v "^#" | cut -f4 | grep -c -P "^A$"
```

**Answer:** 1,150,985 entries

*Note: The pattern `^A$` ensures we match exactly "A" (nothing else), excluding header lines.*

### Question 7: How many entries have exactly 20 supporting reads (read depth)?

Read depth is indicated by the 'DP=' field in column 8:

```bash
cut -f8 exome.full.mpileup.vcf | grep -c "DP=20"
```

**Answer:** 1,816 entries

### Question 8: How many alignments contained insertions and/or deletions in full match alignment?

```bash
samtools view exome.bt2.sam | awk '{if($6 ~ /[ID]/) count++} END {print count}'
```

**Answer:** 2,782 alignments

### Question 9: How many alignments contained insertions and/or deletions in local alignment?

```bash
samtools view exome.local.bt2.sam | awk '{if($6 ~ /[ID]/) count++} END {print count}'
```

**Answer:** 2,613 alignments

### Question 10: How many entries represent indels?

```bash
bcftools view -v indels exome.full.mpileup.vcf | grep -v "^#" | wc -l
```

**Answer:** 1,972 entries

### Question 11: How many entries are reported for position 175672 on Chr1?

**Step 1:** Filter for Chr1 entries
```bash
cut -f1,2 exome.full.mpileup.vcf | grep -v "^#" | grep "Chr1" > chr.txt
```

**Step 2:** Count entries at specific position
```bash
grep "175672" chr.txt | wc -l
```

**Answer:** 24 entries

## Advanced Variant Calling

### Question 12: How many variants are called on Chr3?

**Step 1:** Generate BCF file
```bash
samtools mpileup -f wu_0.v7.fas -g exome.bt2.sorted.bam > out.full.mpileup.bcf
```

**Step 2:** Call variants
```bash
bcftools call -m -v -O v out.full.mpileup.bcf > out.final.vcf
```

**Step 3:** Count Chr3 variants
```bash
cat out.final.vcf | grep -v "^#" | cut -f1 | sort | uniq -c | grep "Chr3"
```

**Answer:** 398 variants

### Question 13: How many variants represent an A→T SNP?

```bash
cat out.final.vcf | grep -v "^#" | cut -f4,5 | grep -P "^A\tT$" | wc -l
```

**Answer:** 392 variants

### Question 14: How many entries are indels?

```bash
bcftools view -v indels out.final.vcf | grep -v "^#" | wc -l
```

**Answer:** 320 entries

### Question 15: How many entries have precisely 20 supporting reads (read depth)?

```bash
grep -v "^#" out.final.vcf | grep -w "DP=20" | wc -l
```

**Answer:** 2 entries

### Question 16: What type of variant is called at position 11937923 on Chr3?

```bash
grep -P "^Chr3\t11937923\t" out.final.vcf
```

**Answer:** SNP

## Summary

This analysis successfully:
- Processed 147,354 raw sequencing reads
- Achieved 93.46% alignment rate with full match and 95.72% with local alignment
- Identified 360,296 potential variant sites on Chr3
- Called 398 final variants on Chr3, including 392 A→T SNPs
- Detected 320 indels in the final variant set

## Tools Used
- **Bowtie2**: Read alignment
- **SAMtools**: File format conversion and mpileup
- **BCFtools**: Variant calling and filtering
- **Docker**: Containerized environment

## Files Generated
- `exome.bt2.sam`: Full match alignment
- `exome.local.bt2.sam`: Local alignment
- `exome.bt2.sorted.bam`: Sorted BAM file
- `exome.full.mpileup.vcf`: Initial variant calls
- `out.final.vcf`: Final variant calls
