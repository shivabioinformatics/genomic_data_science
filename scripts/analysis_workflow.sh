#!/bin/bash

# Genomic Data Analysis Workflow
# Johns Hopkins University - Data Science Course Project
# Analysis of Arabidopsis thaliana Wu_0_A strain

echo "=== Genomic Data Analysis Workflow ==="
echo "Analyzing Arabidopsis thaliana Wu_0_A strain data"
echo ""

# Set data directory
DATA_DIR="../data/gencommand_proj2_data"
BAM_FILE="$DATA_DIR/athal_wu_0_A.bam"
GTF_FILE="$DATA_DIR/athal_wu_0_A_annot.gtf"

echo "Data files:"
echo "BAM file: $BAM_FILE"
echo "GTF file: $GTF_FILE"
echo ""

# 1. Basic BAM file information
echo "=== 1. Basic BAM File Information ==="
echo "Command: samtools view -H $BAM_FILE | head -20"
echo "Purpose: View BAM file header information"
echo ""

# 2. Count total reads
echo "=== 2. Count Total Reads ==="
echo "Command: samtools view -c $BAM_FILE"
echo "Purpose: Count total number of reads in BAM file"
echo ""

# 3. Count mapped reads
echo "=== 3. Count Mapped Reads ==="
echo "Command: samtools view -c -F 4 $BAM_FILE"
echo "Purpose: Count reads that are mapped (not unmapped)"
echo ""

# 4. Count unmapped reads
echo "=== 4. Count Unmapped Reads ==="
echo "Command: samtools view -c -f 4 $BAM_FILE"
echo "Purpose: Count reads that are unmapped"
echo ""

# 5. View first few alignments
echo "=== 5. View Sample Alignments ==="
echo "Command: samtools view $BAM_FILE | head -10"
echo "Purpose: View first 10 alignment records"
echo ""

# 6. Convert BAM to BED format
echo "=== 6. Convert BAM to BED Format ==="
echo "Command: bedtools bamtobed -i $BAM_FILE > reads.bed"
echo "Purpose: Convert BAM file to BED format for interval analysis"
echo ""

# 7. Intersect reads with genes
echo "=== 7. Intersect Reads with Genes ==="
echo "Command: bedtools intersect -a reads.bed -b $GTF_FILE > reads_in_genes.bed"
echo "Purpose: Find reads that overlap with annotated genes"
echo ""

# 8. Count reads per chromosome
echo "=== 8. Count Reads Per Chromosome ==="
echo "Command: samtools idxstats $BAM_FILE"
echo "Purpose: Get read counts for each chromosome/contig"
echo ""

# 9. Coverage analysis
echo "=== 9. Coverage Analysis ==="
echo "Command: samtools depth $BAM_FILE | head -20"
echo "Purpose: Calculate coverage at each position"
echo ""

# 10. Quality statistics
echo "=== 10. Quality Statistics ==="
echo "Command: samtools flagstat $BAM_FILE"
echo "Purpose: Get comprehensive alignment statistics"
echo ""

echo "=== Analysis Complete ==="
echo "Note: This script demonstrates the workflow. Actual commands may vary based on specific questions." 