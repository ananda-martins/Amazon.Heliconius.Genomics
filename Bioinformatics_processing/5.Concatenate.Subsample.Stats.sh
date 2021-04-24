##########################################################################################
# Created by Ananda Martins
# Based on https://speciationgenomics.github.io
# McGill University 
# Update: Feb 14th, 2021 

# Concatenate vcf files per chromosome and vcf stats


# Software: 
# samtools
# bcftools


# How many unfiltered variants?
#Ex: bcftools view -H 1_chr.vcf.gz | wc -l # Takes a lot of time (not a good idea to use it)
		
# Applying filters: Important to perform some stats analyses on a subset of the vcf files to get an idea of how to set filters.
	

#######
## 1 ##  Concatenate VCF files
#######

#!/bin/bash

## Path to timer script - not mandatory, just if you want to check how long takes to run the analysis
PATHTIMER="/home/anandam/projects/def-barrett/anandam/Genomics"
# Activate/load timer script
source ${PATHTIMER}/timer.sh

startimer

VCF="/home/anandam/projects/def-barrett/anandam/Genomics/varcal/melpomene"
OUTPUT="/home/anandam/projects/def-barrett/anandam/Genomics/varcal/melpomene/concat"
bcftools concat $VCF/*.vcf.gz -Oz -o $OUTPUT/melpomene.vcf.gz

endtimer

#######
## 2 ##  Indexing concatenated VCF 
#######

# Change permission of scaffold names
	#chmod +rwx melpomene.vcf.gz

bcftools index /home/anandam/projects/def-barrett/anandam/Genomics/varcal/melpomene/merged/melpomene.vcf.gz


	# It takes time to perform operations on a large VCF. For this reason it is 
	# a good idea to subsample our variant calls and get an idea of the general 
	# distribution of a few key attributes of the data. It is important to know 
	# first your data and from that determine some sensible thresholds.
	
	# Randomly subsampling a VCF: using vcflib pipeline (https://github.com/vcflib/vcflib)
	# Note that vcfrandomsample cannot handle an uncompressed VCF, so we first 
	# open the file using bcftools and then pipe it to the vcfrandomsample 
	# utility. We set only a single parameter, -r which is a bit confusingly 
	# named for the rate of sampling. This essentially means the fraction of 
	# variants we want to retain. This will give us at least 200 K variants, 
	# depending on the random seed used to start the process. This means that 
	# everyone in the class should get different random subsets - providing us 
	# a nice demonstration for the next step.


#######
## 3 ##  Subsampling concatenated vcf file
#######

#!/bin/bash

## Path to timer script - not mandatory, just if you want to check how long takes to run the analysis
PATHTIMER="/home/anandam/projects/def-barrett/anandam/Genomics"
# Activate/load timer script
source ${PATHTIMER}/timer.sh
	
startimer
	
# Assign directory with vcf files
VCF_files="/home/anandam/projects/def-barrett/anandam/Genomics/varcal/melpomene/concat"

# Change permission of sequence files
#chmod +rwx melpomene.vcf.gz
#chmod +rwx melpomene.vcf.gz.csi

bcftools view $VCF_files/melpomene.vcf.gz | vcfrandomsample -r 0.020 -p 315 > $VCF_files/melpomene.subset.vcf

endtimer

# -p --random-seed N (used random number generator)

#######
## 4 ##  Compress and index our new subset VCF to make it easier to access
#######

# Change permission of sequence files
#chmod +rwx melpomene.subset.vcf

bgzip /home/anandam/projects/def-barrett/anandam/Genomics/varcal/melpomene/merged/melpomene.subset.vcf
bcftools index /home/anandam/projects/def-barrett/anandam/Genomics/varcal/melpomene/merged/melpomene.subset.vcf.gz


#######
## 5 ##  Generating statistics from the subset concatenated VCF 
#######

#!/bin/bash

## Path to timer script - not mandatory, just if you want to check how long takes to run the analysis
PATHTIMER="/home/anandam/projects/def-barrett/anandam/Genomics"
# Activate/load timer script
source ${PATHTIMER}/timer.sh

startimer

# Declaring variables
VCF=/home/anandam/projects/def-barrett/anandam/Genomics/varcal/melpomene/concat
OUT=/home/anandam/projects/def-barrett/anandam/Genomics/varcal/melpomene/concat/stats

# Change permission of sequence files
#chmod +rwx melpomene.subset.vcf.gz
#chmod +rwx melpomene.subset.vcf.gz.csi


# Calculate allele frequency
vcftools --gzvcf ${VCF}/melpomene.subset.vcf.gz --freq2 --out ${OUT} --max-alleles 2

# Calculate mean depth per individual
vcftools --gzvcf ${VCF}/melpomene.subset.vcf.gz --out ${OUT} --depth 

# Calculate mean depth per site
vcftools --gzvcf ${VCF}/melpomene.subset.vcf.gz --out ${OUT} --site-mean-depth 

# Calculate site quality
vcftools --gzvcf ${VCF}/melpomene.subset.vcf.gz --out ${OUT} --site-quality 

# Calculate proportion of missing data per individual
vcftools --gzvcf ${VCF}/melpomene.subset.vcf.gz --out ${OUT} --missing-indv 

# Calculate proportion of missing data per site
vcftools --gzvcf ${VCF}/melpomene.subset.vcf.gz --out ${OUT} --missing-site 

# Calculate heterozygosity and inbreeding coefficient per individual
vcftools --gzvcf ${VCF}/melpomene.subset.vcf.gz --out ${OUT} --het 

endtimer