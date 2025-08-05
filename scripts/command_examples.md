# Genomic Analysis Commands

This document contains the specific commands used to analyze the Arabidopsis thaliana Wu_0_A strain genomic data.

## Data Files
- `athal_wu_0_A.bam` - Aligned sequencing reads
- `athal_wu_0_A_annot.gtf` - Gene annotations  
- `athal_wu_0_A.bam.bai` - BAM index file


### 1. Samtools 

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

Practice Questions



 First, a couple of introductory comments. A BAM file contains alignments for a set of input reads. Each read can have 0 (none), 1 or multiple alignments on the genome. These questions explore the relationships between reads and alignments.
 The number of alignments is the number of entries, excluding the header, contained in the BAM file, or equivalently in its SAM conversion.

#### 1. How many alignments does the set contain?  

```bash
samtools flagstats athal_wu_0_A.bam
```
or 
```bash 
samtools view athal_wu_0_A.bam | wc -l
```
**output**: 
flagstats athal_wu_0_A.bam
221372 + 0 in total (QC-passed reads + QC-failed reads) # numbers of alignmnet 
221372 + 0 primary
0 + 0 secondary
0 + 0 supplementary
20419 + 0 duplicates
20419 + 0 primary duplicates
218469 + 0 mapped (98.69% : N/A)
218469 + 0 primary mapped (98.69% : N/A)
155851 + 0 paired in sequencing # how many of them are properly paired
77895 + 0 read1
77956 + 0 read2
142663 + 0 properly paired (91.54% : N/A)
150045 + 0 with itself and mate mapped
2903 + 0 singletons (1.86% : N/A)
4938 + 0 with mate mapped to a different chr
2120 + 0 with mate mapped to a different chr (mapQ>=5) 

#### 2. How many alignments show the read’s mate unmapped?
 
```bash
samtools view athal_wu_0_A.bam | cut -f7 | grep -c "*"
```
output: 65521

An alignment with an unmapped mate is marked with a ‘*’ in column 7. Note that the question asks for alignments, not reads, so we simply count the number of lines in the SAM file with a ‘*’ in column 7:
= (equals sign):

Means the mate is mapped to the same chromosome as the current read
Example: If current read is on Chr3, mate is also on Chr3

* (asterisk):

Means the mate is unmapped (has no valid alignment)
The mate read couldn't be aligned to the reference genome

#### 3. How many alignments contain a deletion (D)? 

In the column 6 of the alignment (bam file), we have CIGAR operations which is stands for ***(Compact Idiosyncratic Gapped Alignment Report)***.
```bash
samtools view athal_wu_0_A.bam | cut -f6 | grep -c 'D'
```
output: 2451

Common operations:

M = Match/mismatch (alignment match)
I = Insertion (bases in read, not in reference)
D = Deletion (bases in reference, not in read)
S = Soft clipping (bases in read, not aligned)
H = Hard clipping (bases not in read)
N = Skipped region (for RNA-seq, represents introns)

Examples:

51M = 51 consecutive matching/mismatching bases
36M = 36 consecutive matching/mismatching bases
50M10I40M = 50 matches, 10 inserted bases, 40 matches
75M25S = 75 matches, 25 bases soft-clipped at end
20S80M = 20 bases soft-clipped at start, 80 matches

All the reads show simple CIGAR strings like 51M or 36M, meaning they are perfect alignments with no insertions, deletions, or clipping - just straight matches to the reference genome.

#### 4. How many alignments show the read’s mate mapped to the same chromosome?
Alignments with the read’s mate mapped to the same chromosome are marked with a ‘=’ in column 7:

```bash
samtools view athal_wu_0_A.bam | cut –f7 | grep –c ‘=’
```
output: 150913

#### 5. How many alignments are spliced?
A spliced alignment will be marked with an ‘N’ (intron gap) in the CIGAR field:


```bash
samtools view athal_wu_0_A.bam | cut -f6 | grep -c 'N'
```
output: 0 
#### 6. How many alignments does the set contain?
We first need to construct the reduced set, i.e. to extract from the original set only those alignments in the specified region. For this, we need to sort and index the file:
This will create the file ‘athal_wu_0_A.sorted.bam’. We then index this file:
This will create the index file ‘athal_wu_0_A.sorted.bam.bai’ in the current directory. Lastly, we extract alignments in the specified range:
The option ‘-b’ will generate output in BAM format. The resulting BAM file will be sorted, so it can be indexed directly if needed.Common pitfalls: make sure to specify the correct reference sequence (‘Chr3’, not ‘chr3’) and exclude ‘,’ when representing the query coordinates. Also, make sure to use the sorted and index BAM file. To determine the number of alignments in the new (region) file, we can use the same commands as for Q1, e.g.:

```bash
samtools view -b athal_wu_0_A.bam Chr3:11777000-11794000 > region_subset.bam
#index the subset
samtools index region_subset.bam

#view the alignment
samtools view -c region_subset.bam Chr3:11777000-11794000
```
output: 7081

#### 7. How many alignments are spliced?
```bash
samtools view region_subset.bam | grep -c 'N'
```
output: 1983

#### 8. How many alignments contain a deletion (D)?
```bash
samtools view region_subset.bam | cut -f6 | grep -c 'D'
```
output: 1019
#### 9. How many alignments show the read’s mate mapped to the same chromosome?

```bash
samtools view region_subset.bam | cut -f7 | grep -c "="
```
output: 4670

#### 10 How many alignments are spliced?

```bash
samtools view region_subset.bam | cut -f6 | grep -c 'N'
```
output: 0 

##### 11: How many sequences are in the genome file? 
This information can be found in the header of the BAM file. Starting with the original BAM file, we use samtools to display the header information and count the number of lines describing the sequences in the reference genome:
```bash
samtools view –H athal_wu_0_A.bam | grep –c “SN:”
```
#### 12: What is the length of the first sequence in the genome file?
The length information is stored alongside the sequence identifier in the header (pattern ‘LN:seq_length’):

```bash
samtools view –H athal_wu_0_A.bam | grep “SN:” | more
```




#### 14: What is the read identifier (name) for the first alignment?
GAII05_0002:1:113:7822:3886#0

#### 15: What is the start position of this read’s mate on the genome? Give this as ‘chrom:pos’ if the read was mapped, or ‘*” if unmapped.
output: chr3 11700332

#### 16: Using BEDtools, examine how many of the alignments at point 2 in the specified range overlap exons at the locus of interest. Use the BEDtools ‘-wo’ option to only report non-zero overlaps. The list of exons is given in the included ‘athal_wu_0_A_annot.gtf’ GTF file. 

```bash 
% cut –f22 overlaps.bed | sort –nrk1 > lengths
```
What the -wo option does:

-w: Reports the number of base pairs of overlap
-o: Only reports overlaps (non-zero)
Combined -wo: Reports only entries with non-zero overlap AND shows the overlap length

output: 3101 
or 
```bash
bedtools intersect –abam athal_wu_0_A.region.bam –b athal_wu_0_A_annot.gtf –bed -wo > overlaps.bed
```
We start by running BEDtools on the alignment set restricted to the specified region (Chr3:11777000-11794000) and the GTF annotation file listed above. To allow the input to be read directly from the BAM file, we use the option ‘-abam’; in this case we will need to also specify ‘-bed’ for the BAM alignment information to be shown in BED column format in the output:
This will create a file with the following format: Columns 1-12 : alignment information, converted to BED format Columns 13-21 : annotation (exon) information, from the GTF file Column 22 : length of the overlap.

Alternatively, we could first convert the BAM file to BED format using ‘bedtools bamtobed’ then use the resulting file in the ‘bedtools intersect’ command. To answer the question, the number of overlaps reported is precisely the number of lines in the file (because only entries in the first file that have overlaps in file B are reported, according to the option ‘-wo’):

#### 17: How many of these are 10 bases or longer?
The size of the overlap is listed in column 22 of the ‘overlaps.bed’ file. To determine those longer than 10 bases, we extract the column, sort numerically in decreasing order, and simply determine by visual inspection of the file the number of such records. For instance, in ‘vim’ we search for the first line listing ‘9’ (‘:/9’), then determine its line number (Ctrl+g). Alternatively, one can use grep with option ‘-n’ to list the lines and corresponding line numbers:

```bash
cut –f22 overlaps.bed | sort –nrk1 > lengths
```
or
```bash
cut –f22 overlaps.bed | sort –nrk1 | grep –n “^9” | head -1
```
For the latter, the last ”10” line will be immediately above the first “9”, so subtract 1 from the answer.

#### 18: How many alignments overlap the annotations?
Columns 1-12 define the alignments:

```bash
cut –f1-12 overlaps.bed | sort –u | wc -l
```
Potential pitfalls: Multiple reads may map at the same coordinates, so the information in columns 1-3 is insufficient. The minimum information needed to define the alignments is contained in columns 1-5, which include the read ID and the flag, specifying whether this is read 1 or read 2 in a pair with the same read ID).

#### 19: Conversely, how many exons have reads mapped to them?
```bash
cut –f13-21 overlaps.bed | sort –u | wc -l
```

#### 20: If you were to convert the transcript annotations in the file “athal_wu_0_A_annot.gtf” into BED format, how many BED records would be generated?

```bash 
cut –f9 athal_wu_0_A.annot.gtf | cut –d ‘ ‘ –f1,2 | sort –u | wc -l
```
output: 4