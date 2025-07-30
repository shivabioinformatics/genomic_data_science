# Genomic Analysis Commands
## Johns Hopkins University - Data Science Course Project

This document contains the specific commands used to analyze the Arabidopsis thaliana Wu_0_A strain genomic data.

## Data Files
- `athal_wu_0_A.bam` - Aligned sequencing reads
- `athal_wu_0_A_annot.gtf` - Gene annotations  
- `athal_wu_0_A.bam.bai` - BAM index file

## SAMtools Commands

### 1. View BAM Header
```bash
samtools view -H athal_wu_0_A.bam | head -20
```
**Purpose**: Examine the BAM file header to understand reference genome and metadata.

### 2. Count Total Reads
```bash
samtools view -c athal_wu_0_A.bam
```
**Purpose**: Get the total number of reads in the BAM file.

### 3. Count Mapped Reads
```bash
samtools view -c -F 4 athal_wu_0_A.bam
```
**Purpose**: Count reads that are mapped (flag 4 is unmapped, so -F 4 excludes unmapped reads).

### 4. Count Unmapped Reads
```bash
samtools view -c -f 4 athal_wu_0_A.bam
```
**Purpose**: Count reads that are unmapped (flag 4 indicates unmapped reads).

### 5. View Sample Alignments
```bash
samtools view athal_wu_0_A.bam | head -10
```
**Purpose**: Examine the first 10 alignment records to understand the data format.

### 6. Get Quality Statistics
```bash
samtools flagstat athal_wu_0_A.bam
```
**Purpose**: Get comprehensive alignment statistics including mapping rates and quality metrics.

### 7. Chromosome/Contig Statistics
```bash
samtools idxstats athal_wu_0_A.bam
```
**Purpose**: Get read counts for each chromosome or contig in the reference genome.

### 8. Coverage Analysis
```bash
samtools depth athal_wu_0_A.bam | head -20
```
**Purpose**: Calculate coverage at each position in the genome.

## BEDtools Commands

### 9. Convert BAM to BED Format
```bash
bedtools bamtobed -i athal_wu_0_A.bam > reads.bed
```
**Purpose**: Convert BAM file to BED format for interval-based analysis.

### 10. Intersect Reads with Genes
```bash
bedtools intersect -a reads.bed -b athal_wu_0_A_annot.gtf > reads_in_genes.bed
```
**Purpose**: Find reads that overlap with annotated genes.

### 11. Count Overlapping Features
```bash
bedtools intersect -c -a reads.bed -b athal_wu_0_A_annot.gtf
```
**Purpose**: Count how many gene features each read overlaps with.

### 12. Coverage of Genes
```bash
bedtools coverage -a athal_wu_0_A_annot.gtf -b reads.bed
```
**Purpose**: Calculate coverage statistics for each gene annotation.

## Unix Command Line Tools

### 13. Count Lines in Files
```bash
wc -l athal_wu_0_A_annot.gtf
```
**Purpose**: Count the number of annotation lines in the GTF file.

### 14. Extract Specific Columns
```bash
cut -f1,3,4,5 athal_wu_0_A_annot.gtf | head -10
```
**Purpose**: Extract specific columns (chromosome, feature type, start, end) from GTF file.

### 15. Filter by Feature Type
```bash
grep "gene" athal_wu_0_A_annot.gtf | wc -l
```
**Purpose**: Count the number of gene features in the annotation file.

### 16. Sort and Count Unique Values
```bash
cut -f3 athal_wu_0_A_annot.gtf | sort | uniq -c
```
**Purpose**: Count the frequency of each feature type in the GTF file.

## Analysis Workflow

1. **Initial Exploration**: Use SAMtools to understand the BAM file structure
2. **Quality Assessment**: Check mapping rates and alignment quality
3. **Data Conversion**: Convert BAM to BED format for interval analysis
4. **Feature Analysis**: Use BEDtools to analyze read-gene relationships
5. **Statistics**: Generate summary statistics for the analysis

## Expected Output Examples

### Mapping Statistics
```
1000000 + 0 in total (QC-passed reads + QC-failed reads)
950000 + 0 mapped (95.00% : N/A)
50000 + 0 unmapped (5.00% : N/A)
```

### Chromosome Distribution
```
Chr1    100000  0       100000
Chr2    80000   0       80000
Chr3    75000   0       75000
```

### Gene Coverage
```
Chr1    gene    1000    2000    gene1   100     50      0.5
Chr1    gene    3000    4000    gene2   200     150     0.75
```

## Notes

- All commands assume you're in the directory containing the data files
- The BAM file must be indexed (BAI file present) for some SAMtools commands
- BEDtools commands require BED format input files
- Results may vary depending on the specific questions being addressed 