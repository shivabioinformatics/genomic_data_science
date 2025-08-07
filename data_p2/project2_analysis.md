# Project 2 Genomic Analysis Scripts
## Developing Variant Calling in n the plant Arabidopsis thaliana, one of the strain to determine genetic variaiton. 

## Docker Setup
```bash
# Access Docker container with genomic tools
docker run -it --rm -v /Users/labuser3/Downloads/genomic_data_science:/workspace gencommand_image:v1.0.0 /bin/bash
```

**Question 1:** How many index files did the operation create?

Before Bowtie2 can align reads to a reference genome, it needs the reference genome to be in a special, efficient format — that's what the index is.

```bash
bowtie2-build wu_0_A.fasta wu_0
```
operation: 6 

**Question 2:** We then examine how many sequences we have in the raw sequences files before we proceed to any analysis. 

```bash
 wc -l wu_0_A_wgs.fastq | awk '{print $1/4}'
```
output: 147354 


Now that I have the indexed files, I dont need to specifiy the extension, need the directory for the insex files without any extension (which was my error), include the threads number and convert it to sam file. You will get a summary statistics when it is done. 

```bash
bowtie2 -p 4 -x wu_0  wu_0_A_wgs.fastq -S exome.bt2.sam
```
-U for single-end seuencing (Unpaired)

output: 
147354 reads; of these: # number of reads
  147354 (100.00%) were unpaired; of these:
    9635 (6.54%) aligned 0 times
    93780 (63.64%) aligned exactly 1 time
    43939 (29.82%) aligned >1 times
93.46% overall alignment rate

for local match 
```bash
bowtie2 -p 4 --local  -x wu_0  wu_0_A_wgs.fastq -S exome.local.bt2.sam
```
147354 reads; of these:
  147354 (100.00%) were unpaired; of these:
    6310 (4.28%) aligned 0 times
    84939 (57.64%) aligned exactly 1 time
    56105 (38.07%) aligned >1 times
95.72% overall alignment rate

**Question 3:** How many reads were mapped in the scenario in Question 7?
137719

This question explores the difference between reads and alignments. As a reminder, a read may have 0 (unmapped), 1 (unique), or multiple alignments in a SAM file. The information on the number of reads in each category is contained in the 'mapped.full.stats' and 'mapped.local.stats' summary files above: simply add the numbers of reads reported to map exactly once and those reported to match multiple times.

**Question 4:** How many reads were mapped in the scenario in Question 8?
141044

This question explores the difference between reads and alignments. As a reminder, a read may have 0 (unmapped), 1 (unique), or multiple alignments in a SAM file. The information on the number of reads in each category is contained in the 'mapped.full.stats' and 'mapped.local.stats' summary files above: simply add the numbers of reads reported to map exactly once and those reported to match multiple times.

Note that by default bowtie2 reports only one match per read, so in this case the number of mapped reads at Q10 and the number of alignments at Q8 will be the same.

Note that by default bowtie2 reports only one match per read, so in this case the number of mapped reads at Q9 and the number of alignments at Q7 will be the same.

**Question 5:** How many entries were reported for Chr3?

```bash
samtools view -bT wu_0.v7.fas exome.bt2.sam > exome.bt2.bam
#then sorting it 
samtools sort exome.bt2.bam exome.bt2.sorted
```
This will create the BAM file 'exom.bt2.sorted.bam', which will be used to determine sites of variation.Step 2: Determine candidate sites using 'samtools mpileup', providing the reference fasta genome (option '-f') and using the option '-uv' to report the output in uncompressed VCF format:

```bash
samtools mpileup -f wu_0.v7.fas -uv exome.bt2.sorted.bam > exome.full.mpileup.vcf
```
Step 3: Count the number of entries in the VCF file located on Chr3. The chromosome information is listed in column 1, once we filter out the header lines (marked with "#"):

```bash
cut -f1 exome.full.mpileup.vcf | grep "Chr3" | wc -l
```
output:360296
***Potential pitfalls: It is critical to sort the BAM file before analyzing it with samtools.***

**Question 6:** How many entries have 'A' as the corresponding genome letter?

I used this command and listed these options: 
A
A
A
A
A
A
A
A
A
GA
A
A

So probably not a good command to use in this case

```bash 
cut -f4 exome.full.mpileup.vcf | grep "A" 
```
What it does:

Extracts REF column (field 4) from all lines (including headers)
Shows lines containing "A" anywhere in the field
Alternative would be : 

```bash
cat exome.full.mpileup.vcf | grep -v "^#" | cut -f4 | grep -c -P "^A$"
```
What it does:
where "^A$" tells 'grep' to look for patterns that consist exclusively of 'A' (i.e., between the start '^' and end '$' of the line).
Excludes header lines (grep -v "^#")
Extracts REF column (field 4)
Counts lines where REF is exactly "A" (nothing else)
Returns a count number

output: 1150985

**Question 7:** How many entries have exactly 20 supporting reads (read depth)?
Read depth is indicated by the 'DP=' field in column 8:

```bash
cut -f8 exome.full.mpileup.vcf | grep -c "DP=20"
```
or 
```bash
cat exome.full.mpileup.vcf | grep -v "^#" | grep -c "DP=20;"
```

output:1816

**Question 8:** How many entries represent indels?

```bash

```

**Question 8:** How many alignments contained insertions and/or deletions, in the scenario in full match alignmetn (no loca)?

```bash
samtools view exome.bt2.sam | awk '{if($6 ~ /[ID]/) count++} END {print count}'
```
output: 2782

**Question 9:** How many alignments contained insertions and/or deletions, in the scenario for local alignment?

```bash
samtools view exome.local.bt2.sam | awk '{if($6 ~ /[ID]/) count++} END {print count}'
```
output: 2613

***Question10*** How many entries represent indels?
```bash
bcftools view -v indels exome.full.mpileup.vcf | grep -v "^#" | wc -l
```
output: 1972

***Question 11*** How many entries are reported for position 175672 on Chr1?

First I prefer to filter with only Chr1 
```bash
cut -f1,2 exome.full.mpileup.vcf | grep -v "^#" | grep "Chr1"  > chr.txt


grep "175672" chr.txt | wc -l
```
output: 24

***Question 12*** How many variants are called on Chr3?

Step 1: First re-run ‘SAMtools mpileup’ with the BCF output option ‘-g’:

My mistake from before was I did not use sorted.bam as a input file
```bash
samtools mpileup -f wu_0.v7.fas -g exome.bt2.sorted.bam > out.full.mpileup.bcf
```
then call variants using ‘BCFtools call’ with the multi-allelic caller (option ‘-m’), showing only variant sites (‘-v’) and presenting the output in uncompressed VCF format (‘-O v’), as instructed:

```bash
bcftools call –m –v –O v out.full.mpileup.bcf > out.final.vcf
```
Step 2: To answer the question, we count all reported variants that show ‘Chr3’ in column 1:
```bash
cat out.final.vcf | grep -v "^#" | cut -f1 | sort | uniq -c | grep "Chr3"

```

output: 398


***Question 13*** How many variants represent an A->T SNP? If useful, you can use ‘grep –P’ to allow tabular spaces in the search term.

```bash
cat out.final.vcf | grep -v "^#" | cut -f4,5 | grep -P "^A\tT$" | wc -l

```
output: 392


***Question 14*** How many entries are indels? 

```bash
bcftools view -v indels out.final.vcf | grep -v "^#" | wc -l

```
output: 320


***Question 15** 
