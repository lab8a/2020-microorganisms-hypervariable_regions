# Started on 2019-08-06
# by Rodrigo García-López for Adrian Ochoa-Leyva's Metagenomics and Metatranscriptomics Lab at IBt, UNAM, Cuernavaca, Mexico.
# Under GNU GPLv3 license
# Disclamer: This script was written for use with specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.
# This script was tested with R 3.6.0 - "Planting of a Tree"
# It calculates the total items per sample
# As it was created for biom-converted tables, the first line is skipped as it is just a comment (please adjust this if source is different). The first two lines have a # at the start of line.
# No taxonomy is present but can be adjusted with a [-1] over the loaded complete table.
# The Rscript bin must to be executed must be the one in my local folder as it has the readr library installed
# Run as follows:
# cat table.tsv|~/bin/R-3.6.0/bin/Rscript contingency_table-unique_items_per_sample.R

df <- read.table(file('stdin'), sep="\t",header=T, skip=1, fileEncoding = "UTF-8", comment.char='')
df <- rowsum(df[2:length(df)], group=df[[1]])
Totals <- c(colSums(df),"Sums"=sum(colSums(df))) #Get the sum totals for each sample and the grand total
names(Totals) <- paste("Total",names(Totals), sep="_") # Rename the headers to identify totals
df[df>1]=1 # convert to composition (absence/presence) table
Uniques <- c(colSums(df),"Sums"=sum(colSums(df))) #Get the sum of unique items for each sample and the grand total
names(Uniques) <- paste("Uniques",names(Uniques), sep="_") # Rename the headers to identify unique items
library("readr")
writeLines(format_tsv(data.frame(t(c(Uniques,Totals)))), stdout()) #output the error to stdout
