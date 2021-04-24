##########################################################################################
# Created by Ananda Martins
# McGill University 
# Update: Apr 12th, 2021 

# Rename sample IDs within vcf files
  # Sometimes the sample IDS have too long and/or too complicated names. 
  # It is always good pratice to keep vcf sample names simple. 

# Software: 
# vcftools
# bcftools


# Check names of samples in vcf file
vcf-query -l melpomene.vcf.gz

# Rename sample IDs

# Create a text file with old and new sample ids (Example)
#nano sampleID.txt

#/home/anandam/scratch/melpomene_aligned2/7IR198.rmd.bam 7IR198
#/home/anandam/scratch/melpomene_aligned2/7IR199.rmd.bam 7IR199
#/home/anandam/scratch/melpomene_aligned2/7IR201.rmd.bam 7IR201
#/home/anandam/scratch/melpomene_aligned2/7IR204.rmd.bam 7IR204
#/home/anandam/scratch/melpomene_aligned2/7IR_205.rmd.bam 7IR_205
#/home/anandam/scratch/melpomene_aligned2/7IR208.rmd.bam 7IR208
#/home/anandam/scratch/melpomene_aligned2/7IR_209.rmd.bam 7IR_209
#/home/anandam/scratch/melpomene_aligned2/7IR89_2016.rmd.bam 7IR89_2016
#/home/anandam/scratch/melpomene_aligned2/7IR91_2016.rmd.bam 7IR91_2016
#/home/anandam/scratch/melpomene_aligned2/7IR92_2016.rmd.bam 7IR92_2016
#/home/anandam/scratch/melpomene_aligned2/7IR94_2016.rmd.bam 7IR94_2016
#/home/anandam/scratch/melpomene_aligned2/7IR95_2016.rmd.bam 7IR95_2016
#/home/anandam/scratch/melpomene_aligned2/7IR96_2016.rmd.bam 7IR96_2016
#/home/anandam/scratch/melpomene_aligned2/7IR97_2016.rmd.bam 7IR97_2016

# Assign variables
SAMPLEID="/home/anandam/projects/def-barrett/anandam/Genomics/varcal/melpomene/vcf.updatesampleIDS"
VCFfile="/home/anandam/projects/def-barrett/anandam/Genomics/varcal/melpomene"

# rename based on sampleID txt file
bcftools reheader -s $SAMPLEID/sampleID.txt $VCFfiles/melpomene.vcf.gz > $SAMPLEID/melpomene.vcf.gz

# Indexing VCF file
bcftools index $SAMPLEID/melpomene.vcf.gz

