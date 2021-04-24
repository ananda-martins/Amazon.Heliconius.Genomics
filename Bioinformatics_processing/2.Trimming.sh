##########################################################################################
# Created by Ananda Martins
# Based on https://speciationgenomics.github.io
# McGill University 
# December 12th, 2020 

# Sometimes Illumina adapter sequences are still present in some reads because 
# and can form adapter dimers or if a DNA fragment is shorter than the read 
# length, the sequencer continues to “read-through” into the adapter at the end
# of the DNA fragment. Trimmomatic is used to trim reads and Illumina adapter 
# sequences.

# Software: Trimmomatic

#!/bin/bash

# Prepare variables
## Path to timer script - not mandatory, just if you want to check how long takes to run the analysis
PATHTIMER="/home/anandam/project/anandam/Genomics"

## Path to raw reads
PATHTOWGSSEQUENCES="/home/anandam/project/anandam/Genomics/00_data/raw_data/melpomene"

## Path to save trimmed sequences
OUTTRIM="/home/anandam/project/anandam/Genomics/2_trim/melpomene/trimmed_seq"

## Path to trimmomatic
TRIMM="/home/anandam/project/anandam/Genomics"

# Change permission of sequence files (needs to do just the first time)
chmod +rwx ${PATHTOWGSSEQUENCES}/*_001.fastq.gz

# Activate/load timer script
source ${PATHTIMER}/timer.sh
            
# Base names - For samples with suffix R1.fastq.gz

startimer

for file in ${PATHTOWGSSEQUENCES}/*_R1_001.fastq.gz
do
printf "\n--------------------------------------------------------\n"
printf "\nNow doing sample\n"
R1="$(basename ${file})"
R2=$(echo $R1| sed 's/_R1_001.fastq.gz/_R2_001.fastq.gz/')
echo $R1
echo $R2
printf "\n"

# Needs trimmomatic/0.36

java -jar $TRIMM/trimmomatic-0.36.jar PE -phred33 \
${PATHTOWGSSEQUENCES}/${R1} \
${PATHTOWGSSEQUENCES}/${R2} \
$OUTTRIM/"${R1%_R1_001.fastq.gz}_R1_trimmed_paired.fastq.gz" \
$OUTTRIM/"${R1%_R1_001.fastq.gz}_R1_trimmed_unpaired.fastq.gz" \
$OUTTRIM/"${R2%_R2_001.fastq.gz}_R2_trimmed_paired.fastq.gz" \
$OUTTRIM/"${R2%_R2_001.fastq.gz}_R2_trimmed_unpaired.fastq.gz" \
ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE-2.fa:3:25:8:2:keepBothReads \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:12
done
endtimer


# Ps.: For trimmomatic parameters check https://www.biostars.org/p/366041/ and
# # http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf



############################ CHANGE TRIMMED FILE NAMES ############################

## Path to trimmed sequences - melpomene
OUTTRIM="/home/anandam/project/anandam/Genomics/2_trim/melpomene/trimmed_seq"

# Change output file names
for file in $OUTTRIM/*_trimmed_paired.fastq.gz
do
mv "$file" "$(echo "$file" | cut -d '.' -f 5,6,7)" # Using "." as separator, will keep parts 5, 6 and 7 of the file name
done
