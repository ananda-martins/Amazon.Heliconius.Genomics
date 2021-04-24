#########################################################################################
# Created by Ananda Martins
# Based on https://speciationgenomics.github.io
# McGill University 
# Update: March 30th, 2021 

#' Stats script for subset of vcf files
#' Helps to decide how to set fliters


# Clear R workspace ---------------------------------------------------------
rm(list=ls())

# Load libraries ------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(gtools)

# Set work directory --------------------------------------------------------
  # Files resulted from 5.Concatenate.Subsample.Stats script
setwd("/Users/anandamartins/VCF_stats/melpomene/subset_stats")


## VCF QUALITY -Phred score --------------------------------------------------
# The first metric we will look at is the (Phred encoded) site quality. 
# This is a measure of how much confidence we have in our variant calls.
# First of all, we read in the site quality report we generated using 
# vcftools. We will use the read_delim command from the readr package 
# (part of the the tidyverse) because it is more efficient for reading 
# in large datafiles. It also allows us to set our own column names.  

var_qual<- read_delim("stats.lqual", delim = "\t",
                      col_names = c("chr", "pos", "qual"), skip = 1)

# Check vcf quality
var_qual

# Plot quality
var_qual_plot <- ggplot(var_qual, aes(qual)) + 
  geom_density(fill = "darkblue", colour = "black", alpha = 0.3) + 
  theme_light() + xlim(0,1200)

# From this we can see that some samples have "low" quality scores - 0-250
# and several have high score ~1000.

# Check plot - better view of samples having 0-300 quality scores
ggplot(var_qual, aes(qual)) + 
  geom_density(fill = "darkblue", colour = "black", alpha = 0.3) + 
  theme_light() + xlim(0,300)

# Check plot - better view of samples having 750-1050 quality scores
ggplot(var_qual, aes(qual)) + 
  geom_density(fill = "darkblue", colour = "black", alpha = 0.3) + 
  theme_light() + xlim(750,1050)

# Remember that a Phred score of 30 represents a 1 in 1000 chance that our 
# SNP call is erroneous. Clearly most sites exceed this - suggesting we have a
# lot of high confidence calls. This is most 
# probably because we have sufficient read depth (as you will see in the 
# next section). However since most sites have a high quality we can 
# see that filtering on this is not really going to be that useful.
# We recommend setting a minimum threshold of 30 and filtering more 
# strongly on other aspects of the data.


#.
#.
#.

## VARIANT MEAN DEPTH -------------------------------------------------------
# Next we will examine the mean depth for each of our variants. 
# This is essentially the number of reads that have mapped to this 
# position. The output we generated with vcftools is the mean of the read 
# depth across all individuals - it is for both alleles at a position and is 
# not partitioned between the reference and the alternative.

## VARIANT MEAN DEPTH -------------------------------------------------------
var_depth <- read_delim("stats.ldepth.mean", delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)


summary(var_depth$mean_depth)
sd(var_depth$mean_depth)

var_depth_plot <- ggplot(var_depth, aes(mean_depth)) + 
  geom_density(fill = "darkblue", colour = "black", alpha = 0.3) + 
  theme_light()

var_depth_plot2 <- ggplot(var_depth, aes(mean_depth)) + 
  geom_density(fill = "darkblue", colour = "black", alpha = 0.3) + 
  theme_light() + xlim(0,50)


# Check plots
var_depth_plot
var_depth_plot2

# This gives a better idea of the distribution. We could set our minimum 
# coverage at the 5 and 95% quantiles but we should keep in mind that 
# the more reads that cover a site, the higher confidence our basecall is.
# 9x as a minimum cutoff for read depth and maximum depth at 36x.

#.
#.
#.

## VARIANT MISSINGNESS ------------------------------------------------------
# Next up we will look at the proportion of missingness at each variant.
# This is a measure of how many individuals lack a genotype at a call 
# site. 
var_miss <- read_delim("stats.lmiss", delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)

summary(var_miss$fmiss)


# Plot 
var_miss_plot <- ggplot(var_miss, aes(fmiss)) + 
  geom_density(fill = "darkblue", colour = "black", alpha = 0.3) + 
  theme_light() + xlim(0,1) + ylim(0,550)


# Check plots
var_miss_plot

# Most sites have almost no missing data(higher density around 0). 
# This means we can be quite conservative when we set our missing data 
# threshold. We will remove all sites where over 10% of 
# individuals are missing a genotype. One thing to note here is that vcftools 
# inverts the direction of missigness, so our 10% threshold means we will 
# tolerate 90% missingness. Typically missingness of 75-95% is used.

#.
#.
#.

## MINOR ALLELE FREQUENCY - This just makes sense when analyzing several samples
# Last of all for our per variant analyses, we will take a look at the 
# distribution of allele frequencies. This will help inform our minor-allele 
# frequency (MAF) thresholds.

var_freq <- read_delim("stats.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)

# find minor allele frequency
var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z)) 

summary(var_freq$maf)

# Plot variant mean depth
mfa_plot <- ggplot(var_freq, aes(maf)) + 
  geom_density(fill = "darkblue", colour = "black", alpha = 0.3) + 
  theme_light()

# Check plots
mfa_plot

# How do we interpret MAF? It is an important measure because low MAF alleles 
# may only occur in one or two individuals. It is possible that some of these 
# low frequency alleles are in fact unreliable base calls - i.e. a source of 
# error.

#.
#.
#.

## MEAN DEPTH PER INDIVIDUAL ------------------------------------------------
ind_depth <- read_delim("stats.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)


# Check mean depth per individual
ind<-as.data.frame(ind_depth$ind)
depths<-as.data.frame(ind_depth$depth)
ind.depths <- cbind(ind,depths)
colnames(ind.depths)[1] <- "ind"
colnames(ind.depths)[2] <- "depth"

# Individuals with low mean depth <9
low.depth.inds<-ind.depths[ind.depths$depth < 9, ] 
# Individuals with high mean depth >36
higher.depth.inds<-ind.depths[ind.depths$depth > 36, ] 

# Plot variant mean depth per individual
ind_depth_plot <- ggplot(ind_depth, aes(depth)) + 
  geom_histogram(fill = "darkblue", colour = "black", alpha = 0.3) + theme_light()

ind_depth_plot
#.
#.
#.

## PROPORTION OF MISSING DATA PER INDIVIDUAL ---------------------------------
ind_miss <- read_delim("stats.imiss", delim = "\t",
                       col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)


# Check proportion of missing data per individual
ind_miss
ind_miss$fmiss
summary(ind_miss$fmiss)


# Plot missing data per individual
ind_miss_plot <- ggplot(ind_miss, aes(fmiss)) + 
  geom_histogram(fill = "darkblue", colour = "black", alpha = 0.3) + theme_light()

ind_miss_plot

#.
#.
#.

## Heterozygosity and inbreeding coefficient per individual
ind_het <- read_delim("stats.het", delim = "\t",
                      col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)


# Check proportion of heterozygosity
ind_het
ind_het$f
summary(ind_het$f)

# Plot heterozigozy and inbreeding
ind_het_plot <- ggplot(ind_het, aes(f)) + 
  geom_histogram(fill = "darkblue", colour = "black", alpha = 0.3) + theme_light()

ind_het_plot