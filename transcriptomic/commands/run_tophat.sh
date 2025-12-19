#!/bin/bash

# TopHat RNA-seq analysis script
# Analysis of Day8 and Day16 samples

# Create tophat directory if it doesn't exist
mkdir -p tophat

# Change to tophat directory as working directory
cd tophat

# Set paths to data files
DATA_DIR="../data"
REFERENCE_GENOME="$DATA_DIR/athal_chr.fa"
GTF_FILE="$DATA_DIR/athal_genes.gtf"
SAMPLE1="$DATA_DIR/Day8.fastq"
SAMPLE2="$DATA_DIR/Day16.fastq"

echo "Starting TopHat analysis..."
echo "Reference genome: $REFERENCE_GENOME"
echo "GTF file: $GTF_FILE"
echo "Samples: Day8.fastq, Day16.fastq"

# Run TopHat for Day8 sample
echo "Processing Day8 sample..."
tophat -G "$GTF_FILE" -o Day8_tophat_out "$REFERENCE_GENOME" "$SAMPLE1"

# Run TopHat for Day16 sample
echo "Processing Day16 sample..."
tophat -G "$GTF_FILE" -o Day16_tophat_out "$REFERENCE_GENOME" "$SAMPLE2"

echo "TopHat analysis completed!"
echo "Results are in:"
echo "  - Day8_tophat_out/"
echo "  - Day16_tophat_out/"