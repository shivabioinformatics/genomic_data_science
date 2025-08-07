# Genomic Analysis Project: Variant Calling in Arabidopsis thaliana

## Overview
This project demonstrates a complete variant calling pipeline for genomic analysis using Arabidopsis thaliana as a model organism. The analysis includes read alignment, variant detection, and statistical analysis of genomic data.

## Project Description
This repository contains the analysis scripts and documentation for performing variant calling on genomic sequencing data. The workflow processes raw sequencing reads, aligns them to a reference genome, and identifies genetic variants including SNPs and indels.

## Key Features
- **Read Alignment**: Using Bowtie2 for both full and local alignment strategies
- **Variant Calling**: Comprehensive pipeline using SAMtools and BCFtools
- **Quality Control**: Multiple filtering steps and statistical analysis
- **Documentation**: Detailed step-by-step analysis with results

## Files Included
- `project2_analysis.md`: Complete analysis workflow with 16 questions and answers
- `history_command`: Command history from the analysis session
- Various genomic data files (excluded from git due to size)

## Analysis Results Summary
- **Total Reads Processed**: 147,354
- **Alignment Rate**: 93.46% (full match) / 95.72% (local match)
- **Variants on Chr3**: 398 total variants
- **A→T SNPs**: 392 variants
- **Indels Detected**: 320 entries

## Tools Used
- **Bowtie2**: Read alignment
- **SAMtools**: File format conversion and mpileup
- **BCFtools**: Variant calling and filtering
- **Docker**: Containerized environment

## Getting Started

### Prerequisites
- Docker installed
- Access to genomic tools container: `gencommand_image:v1.0.0`

### Quick Start
1. Clone this repository
2. Access the Docker container:
   ```bash
   docker run -it --rm -v /path/to/workspace:/workspace gencommand_image:v1.0.0 /bin/bash
   ```
3. Follow the analysis workflow in `project2_analysis.md`

## Analysis Workflow
The complete analysis workflow is documented in `project2_analysis.md` and includes:

1. **Docker Setup**: Environment preparation
2. **Index Creation**: Reference genome indexing
3. **Read Analysis**: Sequence counting and quality assessment
4. **Read Alignment**: Full and local alignment strategies
5. **Variant Calling**: Basic and advanced variant detection
6. **Results Analysis**: Statistical analysis and interpretation

## Contributing
This is an educational project demonstrating genomic analysis techniques. Feel free to use this as a reference for your own genomic analysis workflows.

## License
This project is for educational purposes. Please ensure you have appropriate permissions for any genomic data used in your analysis.

## Contact
For questions about this analysis workflow, please refer to the detailed documentation in `project2_analysis.md`. 