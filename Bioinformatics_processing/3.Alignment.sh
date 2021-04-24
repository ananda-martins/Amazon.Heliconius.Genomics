##########################################################################################
# Created by Ananda Martins
# McGill University 
# Update: January 22nd, 2021 

# Mapping genomes

# One of the first steps we need to take along our pathway to population/speciation
# genomic analyses is mapping our data to a reference genome. Using an alignment 
# software, we essentially find where in the genome our reads originate from and 
# then once these reads are aligned, we are able to either call variants or 
# construct a consensus sequence for our set of aligned reads.

# Softwares: 
# bwa
## More about bwa: http://bio-bwa.sourceforge.net
# samtools
## More about Samtools: http://www.htslib.org/doc/samtools.html
# Picard tools
## More about Picard tools: https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates

#######
## 1 ##  Getting access to the reference genome
#######

# Reference genomes location:

# Download reference genome
# H. erato lativitta - v1 scaffolds (ray phenotype)
wget http://download.lepbase.org/v4/sequence/Heliconius_erato_lativitta_v1_-_scaffolds.fa.gz

# H. melpomeme melpomene - Hmel2.5 - scaffolds (postman phenotype)
wget http://download.lepbase.org/v4/sequence/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa.gz

# Decompress reference genome files
gunzip *scaffolds.fa.gz


#######
## 2 ##  Indexing the reference genome
#######

# Before we can actually perform an alignment, we need to index the reference genome
# This essentially produces an index for rapid searching and aligning. 
# We use the bwa index tool to achieve this. However, the command takes some time 
# to run.

#### Index H. melpomeme melpomene - Hmel2.5 - scaffolds (postman phenotype)
# Needs to have bwa/0.7.17 installed
bwa index Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa

#### Index H. erato lativitta - v1 scaffolds (ray phenotype)
bwa index Heliconius_erato_lativitta_v1_-_scaffolds.fa

#######
## 3 ##  Performing paired end alignments - example with erato
#######

# Important: When bwa aligns reads, it needs access to these files, so they should
# be in the same directory as the reference genome. Then when we actually run the 
# alignment, we tell bwa where the reference is and it does the rest.

# Make a copy of indexed reference genomes and save in  the same directory as the trimmed reads.

#nano 3_mapping_erato_B2.job

#!/bin/bash

#bwa/0.7.17
#samtools/1.3.1

## Path to timer's script - not mandatory, just if you want to check how long takes to run the analysis
PATHTIMER="/home/anandam/project/anandam/Genomics"
# Activate/load timer script
source ${PATHTIMER}/timer.sh

# Set variables for genome refrence and trimmed samples
REFGEN=/home/anandam/project/anandam/Genomics/2_trim/erato/trimmed_seq/1.paired_seq/Batch_2/1_Heliconius_erato_lativitta_v1_-_scaffolds.fa
SAMPLES=/home/anandam/project/anandam/Genomics/2_trim/erato/trimmed_seq/1.paired_seq/Batch_2

# Change permission of sequence files
chmod +rwx ${SAMPLES}/*.fastq.gz

# Making simple list of names
startimer # Start timer

INDS=($(for i in $SAMPLES/*R1_trimmed_paired.fastq.gz;
        do echo $(basename ${i%_R*});
        done))

for IND in ${INDS[@]};
do
# declare variables
FORWARD=/home/anandam/project/anandam/Genomics/2_trim/erato/trimmed_seq/1.paired_seq/Batch_2/${IND}_R1_trimmed_paired.fastq.gz
REVERSE=/home/anandam/project/anandam/Genomics/2_trim/erato/trimmed_seq/1.paired_seq/Batch_2/${IND}_R2_trimmed_paired.fastq.gz
OUTPUT=/home/anandam/project/anandam/Genomics/3_mapping/erato/align/BATCH2/${IND}_sort.bam


# then align and sort
echo "Aligning ${IND} with bwa"
bwa mem -M -t 16 $REFGEN $FORWARD $REVERSE | \
samtools sort -@ 16 -O bam -o $OUTPUT
printf "\n--------------------------------------------------------\n"
done   

endtimer # Finish Timer 14:01:33


# ---------------------------------------------------------------------------------------


#######
## 4 ##  Remove duplicate reads - example with melpomene
#######

#!/bin/bash


# picard/2.20.6 software

## Path to timer script - not mandatory, just if you want to check how long takes to run the analysis
PATHTIMER="/home/anandam/project/anandam/Genomics"
# Activate/load timer script
source ${PATHTIMER}/timer.sh

# Assign variable: path to mapped reads (bam files)
SAMPLES=/home/anandam/project/anandam/Genomics/melpomene/align/BATCH1

# Assign variable for sample names (mapped reads - bam files)
startimer
INDS=($(for i in $SAMPLES/*_sort.bam;
        do echo $(basename ${i%_sort.bam});
        done))

# Run picard MarkDuplicates to each individual
## Remove duplicate reads from the dataset to avoid PCR duplicates and technical 
## duplicates which inflate our sequencing depth and give us false certainty in 
## the genotype calls.
for IND in ${INDS[@]};
do
# Declare variables
Inputfiles=/home/anandam/project/anandam/Genomics/melpomene/align/BATCH1
OutputDir=/home/anandam/project/anandam/Genomics/melpomene/align/BATCH1/duplicated_removed
OUTPUTMetrics=/home/anandam/project/anandam/Genomics/melpomene/align/BATCH1/duplicated_removed/metrics
PATHPICARD=/home/anandam/project/anandam/Genomics

# Run picard - Markduplicates
printf "\n--------------------------------------------------------\n"
echo "Removing duplictes of ${IND}"
printf "\n--------------------------------------------------------\n"
java -jar $PATHPICARD/picard.jar \
MarkDuplicates REMOVE_DUPLICATES=true \
ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
INPUT=${Inputfiles}/${IND}_sort.bam \
OUTPUT=${OutputDir}/${IND}.rmd.bam \
METRICS_FILE=${OUTPUTMetrics}/${IND}.rmd.bam.metrics
done
endtimer


#######
## 5 ##  Index all bam files- example with erato
#######

#!/bin/bash

## Path to timer script - not mandatory, just if you want to check how long takes to run the analysis
PATHTIMER="/home/anandam/project/anandam/Genomics"
# Activate/load timer script
source ${PATHTIMER}/timer.sh

#samtools/1.3.1

startimer
BAMFILES=/home/anandam/project/anandam/Genomics/3_mapping/erato/align/BATCH1/duplicated_removed

INDS=($(for i in $BAMFILES/*.rmd.bam;
        do echo $(basename ${i%.rmd.bam});
        done))

for IND in ${INDS[@]};
do
printf "\n--------------------------------------------------------\n"
echo "INdexing ${IND}"
samtools index ${BAMFILES}/${IND}.rmd.bam
done

endtimer