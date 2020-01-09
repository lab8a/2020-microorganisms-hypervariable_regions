# Started on 2019-09-03
# by Rodrigo García-López for Adrian Ochoa-Leyva's Metagenomics and Metatranscriptomics Lab at IBt, UNAM, Cuernavaca, Mexico.
# Under GNU GPLv3 license
# Disclamer: This script was written for use with specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.
# This script was tested with R 3.6.0 - "Planting of a Tree"
# It calculates the alpha value of each column in a contingency table
# As it was created for biom-converted tables, the first line is skipped as it is just a comment (please adjust this if source is different). The first two lines have a # at the start of line.
# No taxonomy is present but can be adjusted with a [-1] over the loaded complete table.
# It first sums the repeated items in the table based on the first column 
# The Rscript bin must to be executed must be the one in my local folder as it has the readr library installed
# Run as follows:
# cat table.tsv|~/bin/R-3.6.0/bin/Rscript calculate_alpha.R
# Tested as: # cat 10_Reduced_sparsity_tables/OTUs_97/silva/06_Final_mix_silva97-table_lvl5.tsv|~/bin/R-3.6.0/bin/Rscript calculate_alpha.R

df <- read.table(file('stdin'), sep="\t",header=T, skip=1, fileEncoding = "UTF-8", comment.char='')
df <- rowsum(df[2:length(df)], group=df[[1]]) # Sum repeated rows
SUM <- apply(df>0,2,sum) #Get the sum of non null observations
library("vegan")
Shannon <- diversity(t(df),index="shannon") # Shannon's entropy index (uncertainty of predicting next draw)
Diversity <-exp(Shannon) # Effective diversity (non-logartithmic version of shannon, true diversity)
Pielou <- Shannon/log(SUM)# Pielou's evenness (distribution of observations among items (less means more abundance variation; sensitive to sample size))
Simpson <- diversity(t(df),index="simpson") # Simpson's measure of dominance (P or drawing two items from same species)
Richness <- t(estimateR(t(df))) # Get the observed species, the chao1 richness, the ACE richness, estimators along with their se
out <- cbind(Richness, Shannon, Diversity, Pielou, Simpson)
write.table(out,"alpha_table.tsv", sep="\t", quote=FALSE, col.names=NA) #Export the table of all alpha indices that were calculated
