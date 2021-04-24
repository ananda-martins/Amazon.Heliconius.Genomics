##########################################################################################
# Created by Ananda Martins
# Based on https://speciationgenomics.github.io
# McGill University 
# Update: Feb 14th, 2021 

# Variant Calling

# Call variants from the alignments
# Software: 
# samtools
# bcftools


#######
## 1 ##  Indexing the genome reference … again
#######

#samtools/1.3.1
# Index melpomene refrence genome
samtools faidx 1_Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa  

# This will create a fasta index, denoted by the .fai suffix.

# Check the first lines of the indexed reference genome (with a suffix .fai)
head 1_Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa.fai
## For clarity, the columns are: chromsome name, chromosome length, offset of the 
## first base, fasta line length, fasta line length + 1.


#######
## 2 ##  Variant calling 
#######

# Create txt file with scaffold names
	# Go to reference genome 
cat 1_Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa.fai | cut -f 1 > names_scaffolds.txt

	# Scaffolds Hmel200001o and after are all short, poorly assembled, and possibly junk. 
	# They constitute <1% of the genome.

# Count total of scaffolds (all starts with Hmel)
grep -c "Hmel" names_scaffolds.txt
	# 332

# Extract short, poorly assembled scaffolds to a new txt file
	# starts at Hmel200001o and ends at Hmel200294o
grep "Hmel200" names_scaffolds.txt > junk.scaffolds.txt

# Count number of junk scaffolds in new file
grep -c "Hmel" junk.scaffolds.txt
	# 294

# Extract lines with good scaffolds and save in a new text file
# Check line numbers in file
awk '{ print NR, $1 }' names_scaffolds.txt
	# Needs to remove from line 39 to 332

# Print lines 1 to 38 and save in a new file
awk 'FNR>=1 && FNR<=38' names_scaffolds.txt > scaffolds_sel.txt
	
# Count number of junk scaffolds in new file
awk '{ print NR, $1 }' scaffolds_sel.txt
	# 38
	
# Create one file per chromosome
	# some chromosomes will have more than one scaffold
# Move txt files to new directory
mkdir chr_files
mv *.txt /home/anandam/projects/def-barrett/anandam/Genomics/reference.genomes/melpomene/chr_files
	
# Ex.: chromosome 1
# Print lines 1 and save in a new file
awk 'FNR==1' scaffolds_sel.txt > 1.chr.txt
	# Always check scaffolds_sel.txt document to visualize if the chromosomes have more than 
	# one scaffold and which lines are equivalent to the scaffolds
	awk '{ print NR, $1 }' scaffolds_sel.txt
	

# -------------------------------------------------------------------------------------- #
# ------------------------------------- SUMMARY ---------------------------------------- #
# -------------------------------------------------------------------------------------- #
	# names_scaffolds.txt -> list of all melpomene scaffolds (junk.scaffolds.txt + scaffolds_sel.txt)
	# junk.scaffolds -> list of short, poorly assembled scaffolds
	# scaffolds_sel.txt -> list of selected scaffolds (without junk)
	# For chromosomes that has only one scaffold - use bcftools mpileup -r option
	# For chromosomes that has more than one scaffold - create a txt file with scaffolds
		# and use bcftools mpileup -R option
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #




# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
# ------------------------------ Chromosome 1 ------------------------------------------ #

#!/bin/bash


## Path to timer script - not mandatory, just if you want to check how long takes to run the analysis
PATHTIMER="/home/anandam/projects/def-barrett/anandam/Genomics"
# Activate/load timer script
source ${PATHTIMER}/timer.sh

# Assign variable: path to mapped reads with removed duplicates (bam files)
SAMPLES=/home/anandam/scratch/melpomene_aligned2

# Change permission of sequence files
	#chmod +rwx ${SAMPLES}/*.bam
	#chmod +rwx ${SAMPLES}/*.bam.bai

# Declare variables
## reference genome
REF=/home/anandam/projects/def-barrett/anandam/Genomics/reference.genomes/melpomene/1_Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa
# Change permission of referrence genome
	#chmod +rwx 1_Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.*
		
## Output variant call files
OUTPUTVC=/home/anandam/projects/def-barrett/anandam/Genomics/varcal/melpomene

# Software: bcftools/1.10.2

startimer # start timer

# Run - variant call 
bcftools mpileup -r Hmel201001o -a AD,DP,SP -Ou -f $REF $SAMPLES/*.bam | bcftools call -f GQ,GP -mO z -o $OUTPUTVC/1_chr.vcf.gz
 
endtimer # endtimer 19h16min43sec

# Ps.: Chromosome 1 has only one scaffold: use option mpileup -r


# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
# ------------------------------ Chromosome 3 ------------------------------------------ #

#!/bin/bash

## Path to timer script - not mandatory, just if you want to check how long takes to run the analysis
PATHTIMER="/home/anandam/projects/def-barrett/anandam/Genomics"
# Activate/load timer script
source ${PATHTIMER}/timer.sh

# Assign variable: path to mapped reads with removed duplicates (bam files)
SAMPLES=/home/anandam/scratch/melpomene_aligned2

# Change permission of sequence files
	#chmod +rwx ${SAMPLES}/*.bam
	#chmod +rwx ${SAMPLES}/*.bam.bai

# Declare variables
## reference genome
REF=/home/anandam/projects/def-barrett/anandam/Genomics/reference.genomes/melpomene/1_Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa
# Change permission of referrence genome
	#chmod +rwx 1_Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.*
	
## Txt file with scaffold names - when chromosome has more than one scaffold
CHR=/home/anandam/projects/def-barrett/anandam/Genomics/reference.genomes/melpomene/chr_files/chr3.txt
# Change permission of scaffold names
	chmod +rwx *.txt
		
## Output variant call files
OUTPUTVC=/home/anandam/projects/def-barrett/anandam/Genomics/varcal/melpomene

startimer # start timer

# Run - variant call 
bcftools mpileup -R $CHR -a AD,DP,SP -Ou -f $REF $SAMPLES/*.bam | bcftools call -f GQ,GP -mO z -o $OUTPUTVC/3_chr.vcf.gz
 
endtimer # endtimer


# Ps.: Chromosome 3 has more than one scaffold: create txt containing scaffold names and use option mpileup -R
# Use the same script for all the chromosomes


#######
## 3 ##  Exploring VCF files
#######

# vcf stands for ‘variant call format’ and is a standard format used for 
# variant calling and in population genomics.

# A nice feature of vcf files is that you can access almost any part of 
# the genome you are interested in. To do this though, you need to index 
# the vcf first.

################### Indexing VCF files per chromosome ###################

# Assign directory with vcf files
VCF_files="/home/anandam/projects/def-barrett/anandam/Genomics/varcal/melpomene"

# Assign variable for sample names
CHRS=($(for i in $VCF_files/*.vcf.gz;
do echo $(basename ${i%.vcf.gz});
done))

# Run index for each sample
for CHR in ${CHRS[@]};
do
printf "\n--------------------------------------------------------\n"
echo "Indexing ${CHR}"
bcftools index ${CHR}.vcf.gz
done