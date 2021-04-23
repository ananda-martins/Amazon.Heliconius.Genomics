##########################################################################################
# Created by  Ananda Martins
# McGill University 
# Edited on December, 22nd, 2020 


# Softwares: 
# fastqc
# multiqc

# Fastqc: runs Quality Control on all the WGS (fastq files) for Heliconius butterflies data
# Multiqc: aggregate results from bioinformatics analyses across many samples into a single report

# Requirements:
# Requires (but not needed if you want to) the "timer.sh" script to calculate the time taken to run analysis

#!/bin/bash

# Prepare variables

# Path to timer (Using to check how long is going to run the analysis)
PROJPATH="/home/anandam/project/anandam/Genomics"

# Source .sh script in bash
source ${PROJPATH}/timer.sh

# Path to sequences
PATHTOWGSSEQUENCES="/home/anandam/project/anandam/Genomics/00_data/raw_data/melpomene"

# Change permission of sequence files
chmod +rwx ${PATHTOWGSSEQUENCES}/*.fastq.gz

# Output path to save fastqc results
OUTQC="/home/anandam/project/anandam/Genomics/1_fastqc/melpomene/fastqc/not_trimmed"

## Necessary to have installed fastqc/0.11.8

# Running FastQC using all WGS files (FastQC on multiple threads)
startimer # Start timer

# FastQC ---------------------------------------------------------------------

fastqc ${PATHTOWGSSEQUENCES}/*.fastq.gz --outdir ${OUTQC}

# Check output, especially the html file
# 1) Per base seq. qual plot is OK (most boxes in the green zone and maybe a couple of towards the end falling below the green)
# 2) Per base seq. content plot has the 4 lines overlapping each other. Sometimes you might see that for the first ~10bases or so, the lines are noisy but from there on they should smoothen out.
# 3) Per seq. GC content plot has single bell shaped hump (more or less) and not two or more. Small shoulders are ok
# 4) If you see adapters in the Adapter content plot, do adapter removal (and maybe trim reads from the end as well for getting rid of low qual bases) and then redo FastQC.

# Prepare variables
MULTIQCPATH="/home/anandam/project/anandam/Genomics/1_fastqc/melpomene/multiqc/not_trimmed"
FASTQCPATH="/home/anandam/project/anandam/Genomics/1_fastqc/melpomene/fastqc/not_trimmed"

# MultiQC ---------------------------------------------------------------------
# Preparing multiQC to run 
# After fastQC ran, you can combine the fastQC output in one report 

## Necessary to have python/3.8

source ~/python_programs/python_multiqc/bin/activate
multiqc --outdir ${MULTIQCPATH} ${FASTQCPATH}
endtimer # End of timer

# View output, especially the html file
# https://www.biostars.org/p/172860/
# Assess if needed to trim
# Assess if quality looks good 
# Assess if there are adapters in the sequences 
