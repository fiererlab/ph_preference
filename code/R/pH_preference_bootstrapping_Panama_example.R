####METHOD TO IDENTIFY ENVIRONMENTAL OPTIMA FROM SEQUENCING DATA 

## Clean/reset environment 
rm(list=ls()) 

## Set working directory
setwd("xxxx")

## Load libraries 
library(dplyr)
library(phyloseq)
library(vegan)
library(ggplot2)
library(permute)
library(lattice)
library(knitr)
library(magrittr)
library(ggpubr)
library("RColorBrewer")

##LOAD DATA
otufile     <- "otu_table_panama.txt"
mapfile     <- "mapfile_panama.txt"

D <- import_qiime(otufilename = otufile, mapfilename = mapfile)
D

#REMOVE MITOCHONDRIA AND CHLOROPLAST
Dt <- D %>% subset_taxa(Family!= "Mitochondria")
Dt <- Dt %>% subset_taxa(Order!= "Chloroplast")
Dt <- Dt %>% subset_taxa(Kingdom!= "Archaea")

#RAREFY DATASET
set.seed(123)
Dt = rarefy_even_depth(Dt, sample.size = min(sample_sums(Dt)), rngseed = T)

Dt <- prune_taxa(taxa_sums(Dt) > 0, Dt) 

##FILTER ASVs THAT OCCUR IN LESS THAN 20 SAMPLES
Dt.filter = filter_taxa(Dt, function(x) sum(x > 0) > 20, TRUE)

Dt.rel <- transform_sample_counts(Dt.filter, function(OTU) OTU/sum(OTU))

new_physeq = Dt.rel

subset.otu = t(otu_table(new_physeq))
subset.ph = as.data.frame(sample_data(new_physeq))$'Soil.pH.water.2017'
subset.otu.ph = as.data.frame(cbind(subset.otu, subset.ph))

##CREATE PH BINS
subset.otu.ph = subset.otu.ph %>% mutate(ph.bin = cut(subset.otu.ph$subset.ph, breaks=c(3.5,4,4.5,5,5.5, 6, 6.5, 7, 7.5, 8, 8.5)))


####BOOTSTRAPPING OF THE RELATIVE ABUNDANCE MAXIMA FOR EACH ASV IN ORDER TO FIND THE PREFERENTIAL ENVIRONMENTAL RANGE. THE METHOD IS VALIDATED AGAINST A RANDOM DISTRIBUTION OF ABUNDANCES TO OBTAIN STATISTICAL SIGNIFICANCE. NON-SIGNIFICANT (P<0.05) BOOTSTRAP TESTS ARE DISCARDED

#1 DEFINE MEDIAN AND BOOTSTRAP FUNCTIONS
samplemax <- function(d, i) {
  return(max(d[i]))
}

myBootFun <- function(d) {
  boot(d, samplemax, R = 1000, sim = "ordinary")
}


#2 BOOTSTRAPPING OF THE RELATIVE ABUNDANCE MAXIMA FOR ALL ASVs

library(boot)

bootstrap.results = list()

for (i in 1:ncol(subset.otu)) {
  
  bootstrap.results[[i]] = myBootFun(subset.otu[,i]) 
  
}  


#3 ESTIMATION OF 95% CONFIDENCE INTERVALS OF THESE MAXIMA FOR EACH ASV (PREFERENTIAL PH RANGE)
ci.results = list()

for (i in 1:length(bootstrap.results)) {
  
  ci.results[[i]] = boot.ci(bootstrap.results[[i]],
                            type = "perc", conf = 0.95)
  }


#4 GET RANGE OF PH BASED ON CONFIDENCE INTERVALS

intervals = list()

for (i in 1:ncol(subset.otu)) {
  
  intervals[[i]] = subset(subset.otu.ph, subset.otu.ph[,i] > ci.results[[i]]$percent[4] & subset.otu.ph[,i] < ci.results[[i]]$percent[5])
  
} 


#SANITY CHECK EXAMPLE
intervals[[91]]$subset.ph
ci.results[[91]]$percent
subset.otu.ph[,91]


#5 CHECK THE RANGE OF PH COVERED BY THESE INTERVALS

ph.intervals = list()

for (i in 1:length(intervals)) {
  
  ph.intervals[[i]] = max(intervals[[i]]$subset.ph) - min(intervals[[i]]$subset.ph)
  
} 


#6 OBTAIN PH PREFERENCE OF EACH ASV

ph.pref = list()

for (i in 1:length(intervals)) {
  
  ph.pref[[i]] = mean(max(intervals[[i]]$subset.ph), min(intervals[[i]]$subset.ph))
  
} 


#7 ASSIGN PH PREFERENCE TO ASVS

ph.intervals.dbb = unlist(ph.intervals)

ph.optima.dbb = data.frame(Range = ph.intervals.dbb, ASV = colnames(subset.otu), tax_table(new_physeq), pH.preference = unlist(ph.pref))

ph.optima.dbb = ph.optima.dbb %>% mutate(ph.bin = cut(ph.optima.dbb$pH.preference, breaks=c(3.5,4,4.5,5,5.5, 6, 6.5, 7, 7.5, 8, 8.5)))


#8 KEEP ONLY ASVs WITH CONFIDENT PH PREFERENCE (PREFERENTIAL PH RANGE =< 0.5) AND EXPORT DATAFRAME

ph.optima.asv.subset.05 = subset(ph.optima.dbb, Range <= 0.5) 

write.csv(ph.optima.asv.subset.05, "xxxx.csv")

