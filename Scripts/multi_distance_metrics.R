# Started on 2019-09-19
# by Rodrigo García-López for Adrian Ochoa-Leyva's Metagenomics and Metatranscriptomics Lab at IBt, UNAM, Cuernavaca, Mexico.
# Under GNU GPLv3 license
# Disclamer: This script was written for use with specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.
# This script was tested with R 3.6.0 - "Planting of a Tree"
# The script is intended to calculate distance matrices from a contingency table (ideally, normalized, rarefactioned or transformed)
# The input is a square matrix of samples x samples using a metric and has no previous comments in the header (starts in line 1)
# Most metrics have been commented, as they are not normally used. Uncomment them as required

# Run as follows:
# cat table.tsv|~/bin/R-3.6.0/bin/Rscript multi_distance_metrics.R

# Tested with command:
# cat /home/rodrigo/Shrimp_2015-2018/V3-V4_comparison/12_diversity_analyses/02_using_sparsity_reduced_tables/07_transformed_tables/01_Multiregion/06_Final_mix_gg_5/06_Final_mix_gg-lvl5-rar-10000-perm.tsv|~/bin/R-3.6.0/bin/Rscript multi_distance_metrics.R test
# In R:
# df <- read.table("12_diversity_analyses/02_using_sparsity_reduced_tables/07_transformed_tables/01_Multiregion/06_Final_mix_gg_5/06_Final_mix_gg-lvl5-rar-10000-perm.tsv", sep="\t",header=T, skip=0, comment.char='',quote="",fill=F, row.names=1)

args <- commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<1) { # at least, one arguments should be included: <out_name_prefix>
  stop("A minimum of 1 argument is mandatory: cat table.tsv|Rscript multi_distance_metrics.R <output_file>", call.=FALSE)
}
prefix <- as.character(args[1]) # Get a string handle to create output names
df <- read.table(file('stdin'), sep="\t",header=T, skip=0, comment.char='',quote="",fill=F, row.names=1) #Read the table, row names expected in first column, matching colnames in first row

library("vegan")
 ### Calculate distance metrics ### 
jaccard <- vegdist(t(df),method="jaccard",upper=TRUE,diag=TRUE)
write.table(as.matrix(jaccard),paste(prefix,"jaccard.tsv", sep="-"), sep="\t", quote=FALSE, col.names=NA, row.names=TRUE)
bray <- vegdist(t(df),method="bray",upper=TRUE,diag=TRUE)
write.table(as.matrix(bray),paste(prefix,"bray.tsv", sep="-"), sep="\t", quote=FALSE, col.names=NA, row.names=TRUE)
euclidean <- vegdist(t(df),method="euclidean",upper=TRUE,diag=TRUE)
write.table(as.matrix(euclidean),paste(prefix,"euclidean.tsv", sep="-"), sep="\t", quote=FALSE, col.names=NA, row.names=TRUE)
canberra <- vegdist(t(df),method="canberra",upper=TRUE,diag=TRUE)
write.table(as.matrix(canberra),paste(prefix,"canberra.tsv", sep="-"), sep="\t", quote=FALSE, col.names=NA, row.names=TRUE)
manhattan <- vegdist(t(df),method="manhattan",upper=TRUE,diag=TRUE)
write.table(as.matrix(manhattan),paste(prefix,"manhattan.tsv", sep="-"), sep="\t", quote=FALSE, col.names=NA, row.names=TRUE)
chao <- vegdist(t(df),method="chao",upper=TRUE,diag=TRUE)
write.table(as.matrix(chao),paste(prefix,"chao.tsv", sep="-"), sep="\t", quote=FALSE, col.names=NA, row.names=TRUE)
horn <- vegdist(t(df),method="horn",upper=TRUE,diag=TRUE)
write.table(as.matrix(horn),paste(prefix,"horn.tsv", sep="-"), sep="\t", quote=FALSE, col.names=NA, row.names=TRUE)
