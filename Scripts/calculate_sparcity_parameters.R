# Started on 2019-08-10
# by Rodrigo García-López for Adrian Ochoa-Leyva's Metagenomics and Metatranscriptomics Lab at IBt, UNAM, Cuernavaca, Mexico.
# Under GNU GPLv3 license
# Disclamer: This script was written for use with specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.
# This script was tested with R 3.6.0 - "Planting of a Tree"
# It calculates 2 parameters for filtering tables to reduce sparsity:
# Parameter 1: min_freq: This is is calculated as 0.0001 of the total sum of items
# Parameter 2: min_sample: This is the 25% of the samples
# Both items are rounded to create integers
# The input:
#   As it was created for biom-converted tables, the first line is skipped as it is just a comment (please adjust this if source is different). The first two lines have a # at the start of line.
#   No taxonomy is present but can be adjusted with a [-1] over the loaded complete table.
# Run as follows:
#   cat table.tsv|~/bin/R-3.6.0/bin/Rscript calculate_smallest_sample_size.R

df <- read.table(file('stdin'),sep="\t", header=T, row.names=1, fileEncoding = "UTF-8",comment.char='',skip=1) #This is set to skip the first line
min_freq=round(sum(df)*0.0001)
min_sample=round(ncol(df)*0.25)
out=paste(paste("Freq cutoff:",min_freq, sep="\t"),paste("Sample cutoff:",min_sample,sep="\t"),sep="\n")
write(out,"")
