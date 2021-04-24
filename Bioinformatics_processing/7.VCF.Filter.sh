##########################################################################################
# Created by Ananda Martins
# Based on https://speciationgenomics.github.io
# McGill University 
# Update: Apr 8th, 2021 

# Filter Variants based on VCF stats - Perform the filter using per chromosome vcf files

# Software: 
# vcftools
# bcftools - concatenate per chromosome vcf files and indexing


#######
## 1 ##  Heliconius melpomene filters
#######

# Applying filters to the VCF based on stats performed in subset vcf concatenated file
# Variant Quality => 30
# Variant Mean depth - minimum filter => 9X
# Variant Mean depth - maximum filter => 36X
# Variant missingness => 90% (meaning that we will remove all sites where over 10% of individuals are missing a genotype)
# Minor allele frequency => 0.01
# Remove individuals with < 9X quality mean depth [7IR201, BAR381, BG191, BURIT38_2016, BURIT40_2016, CENT146, CENT147, NOVA230]
# Remove individuals with > 36X quality mean depth [BOQ113, ISA308, MURU246, ZP86]

#!/bin/bash

## Path to timer script - not mandatory, just if you want to check how long takes to run the analysis
PATHTIMER="/home/anandam/projects/def-barrett/anandam/Genomics"
# Activate/load timer script
source ${PATHTIMER}/timer.sh

startimer # start timer

# Assign variables: path to vcf files (per chromosome) and to filtered output directory
VCF_IN=/home/anandam/projects/def-barrett/anandam/Genomics/varcal/melpomene/vcf.updatesampleIDS
VCF_OUT=/home/anandam/projects/def-barrett/anandam/Genomics/varcal/melpomene/vcf.updatesampleIDS/VCF_filtered
REMOVESAMPLES=/home/anandam/projects/def-barrett/anandam/Genomics/varcal/melpomene/jobs/list.samples.to.delete.txt

# Assign variables for filtering parameters
QUAL=30
MIN_DEPTH=9
MAX_DEPTH=36
MISS=0.9
MAF=0.01

# Change permission of files
#chmod +rwx ${VCF_IN}/*.vcf.gz
#chmod +rwx ${VCF_IN}/*.csi

# Assign variable for each chromosome
CHRS=($(for i in $VCF_IN/*.vcf.gz;
        do echo $(basename ${i%.vcf.gz});
        done))	

# Run filter for each chromosome file
for CHR in ${CHRS[@]};
do
printf "\n--------------------------------------------------------\n"
echo "Filtering ${CHR}"
vcftools --gzvcf $VCF_IN/${CHR}.vcf.gz \
--remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL \
--min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
--remove $REMOVESAMPLES \
--minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | bgzip -c > \
$VCF_OUT/${CHR}.filtered.vcf.gz
done

# --remove option includes a file.txt with individuals to be excluded from
# the analysis.

endtimer

# 
#
#

#######
## 2 ##  Concatenate chromosome VCF files
#######

#!/bin/bash

## Path to timer script - not mandatory, just if you want to check how long takes to run the analysis
PATHTIMER="/home/anandam/projects/def-barrett/anandam/Genomics"
# Activate/load timer script
source ${PATHTIMER}/timer.sh

VCF="/home/anandam/projects/def-barrett/anandam/Genomics/varcal/melpomene/updates.vcf/VCF_filtered"
bcftools concat $VCF/*.vcf.gz -Oz -o $VCF/melpomene.vcf.gz

endtimer


#######
## 3 ##  Index vcf filtered file
#######

bcftools index melpomene.vcf.gz


# Check number of variants

bcftools index --nrecords melpomene.vcf.gz # # get total variant count
## 9,436,942

bcftools index --stats melpomene.vcf.gz # # get variant count per chromsome

#.
#.
#.
#.

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------

#######
## 1 ##  Heliconius erato filters
#######

# Applying filters to the VCF based on stats performed in subset vcf concatenated file
# Variant Quality => 30
# Variant Mean depth - minimum filter => 9X
# Variant Mean depth - maximum filter => 30X
# Variant missingness => 90% (meaning that we will remove all sites where over 10% of individuals are missing a genotype)
# Minor allele frequency => 0.01
# Remove individuals with < 9X quality mean depth [CPi227, IG83, IG85, NOVA229, NOVA231, PR118, SAG36, SAG39, TRA263, UTI339]
# Remove individuals with > 30X quality mean depth [BAR380]

#!/bin/bash

## Path to timer script - not mandatory, just if you want to check how long takes to run the analysis
PATHTIMER="/home/anandam/projects/def-barrett/anandam/Genomics"
# Activate/load timer script
source ${PATHTIMER}/timer.sh

startimer # start timer

# Assign variables: path to vcf files (per chromosome) and to filtered output directory
VCF_IN=/home/anandam/projects/def-barrett/anandam/Genomics/varcal/erato/vcf.updatesampleIDS
VCF_OUT=/home/anandam/projects/def-barrett/anandam/Genomics/varcal/erato/vcf.updatesampleIDS/VCF_filtered
REMOVESAMPLES=/home/anandam/projects/def-barrett/anandam/Genomics/varcal/erato/jobs/list.samples.to.delete.txt

# Assign variables for filtering parameters
QUAL=30
MIN_DEPTH=9
MAX_DEPTH=30
MISS=0.9
MAF=0.01

# Change permission of files
#chmod +rwx ${VCF_IN}/*.vcf.gz
#chmod +rwx ${VCF_IN}/*.csi

# Assign variable for each chromosome
CHRS=($(for i in $VCF_IN/*.vcf.gz;
        do echo $(basename ${i%.vcf.gz});
        done))	

# Run filter for each chromosome file
for CHR in ${CHRS[@]};
do
printf "\n--------------------------------------------------------\n"
echo "Filtering ${CHR}"
vcftools --gzvcf $VCF_IN/${CHR}.vcf.gz \
--remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL \
--min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
--remove $REMOVESAMPLES \
--minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | bgzip -c > \
$VCF_OUT/${CHR}.filtered.vcf.gz
done

# --remove option includes a file.txt with individuals to be excluded from
# the analysis.

endtimer

#
#
#


#######
## 2 ##  Concatenate chromosome VCF files
#######

#!/bin/bash

## Path to timer script - not mandatory, just if you want to check how long takes to run the analysis
PATHTIMER="/home/anandam/projects/def-barrett/anandam/Genomics"
# Activate/load timer script
source ${PATHTIMER}/timer.sh

startimer

VCF="/home/anandam/projects/def-barrett/anandam/Genomics/varcal/erato/vcf.updatesampleIDS/VCF_filtered"
OUTPUT="/home/anandam/projects/def-barrett/anandam/Genomics/varcal/erato/vcf.updatesampleIDS/VCF_filtered/concatenated"
bcftools concat $VCF/*.vcf.gz -Oz -o $OUTPUT/erato.vcf.gz

endtimer



#######
## 3 ##  Index vcf filtered file
#######

bcftools index erato.vcf.gz

#.
#.
#.
#.

