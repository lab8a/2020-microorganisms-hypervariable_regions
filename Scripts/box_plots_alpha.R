#!/usr/bin/env /home/rodrigo/bin/R-3.6.0/bin/Rscript
# Started on 2019-09-03
# by Rodrigo García-López for Adrian Ochoa-Leyva's Metagenomics and Metatranscriptomics Lab at IBt, UNAM, Cuernavaca, Mexico.
# Under GNU GPLv3 license
# Disclamer: This script was written for use with specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.
# This script was tested with R 3.6.0 - "Planting of a Tree"
# It was originally intended to compare greengenes and silva results by creating boxplots and compare two populations of alpha indices per index(Chao1 and Shannon [columns 2 and 6] by default in out data)
# This takes two tables, both containing alpha indices calculated for
# Input tables must have samples as rows (order is important as tests are for paired samples), indices as columns (S.obs	S.chao1	se.chao1	S.ACE	se.ACE	Shannon	Diversity	Pielou	Simpson)
# Run as follows:
# ~/bin/R-3.6.0/bin/Rscript box_plots_alpha.R <file1> <file2>
# Tested with command:
# ~/bin/R-3.6.0/bin/Rscript box_plots_alpha.R 12_diversity_analyses/02_using_sparsity_reduced_tables/04_alpha_diversity_R/06_Final_mix_gg_alpha_table_lvl5.tsv 12_diversity_analyses/02_using_sparsity_reduced_tables/04_alpha_diversity_R/06_Final_mix_silva_alpha_table_lvl5.tsv

args <- commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<2) {
  stop("At least two tables must be supplied as arguments: Rscript barplots_cluster_compositions.R <file1> <file2>", call.=FALSE)
}
df1 <- read.table(file=args[1],sep="\t", header=T, row.names=1, fileEncoding = "UTF-8",comment.char='',skip=0)
df2 <- read.table(file=args[2],sep="\t", header=T, row.names=1, fileEncoding = "UTF-8",comment.char='',skip=0)
# First chao1 (2nd column)
t <- round(t.test(df1[,2],df2[,2],paired=TRUE)$p.value,4)
w <- round(wilcox.test(df1[,2],df2[,2],paired=TRUE)$p.value,4)
s1 <- round(shapiro.test(df1[,2])$p.value,4)
s2 <- round(shapiro.test(df2[,2])$p.value,4)
min <- min(df1[,2],df2[,2])
max <- max(df1[,2],df2[,2])
step <- (max-min)/5
pdf("chao1.pdf")
boxplot(df1[,2],df2[,2],border=c("red","blue2"),col=c("pink","lightblue"),lwd=3,xaxt='n',yaxt='n',frame=FALSE)
# axis(1, at=1:2,labels=c("GG","Silva"),cex.axis=2)
axis(1, at=1.5,labels=paste("p-values SW1:",s1,"/ SW2:",s2,"/ t:",t,"/ w:",w),cex.axis=1.5)
axis(2,las=1,at=seq(min,max,step), labels=round(seq(min,max,step),1),cex.axis=3)
dev.off()

# Now for Shannon
t <- round(t.test(df1[,6],df2[,6],paired=TRUE)$p.value,4)
w <- round(wilcox.test(df1[,6],df2[,6],paired=TRUE)$p.value,4) # Carry out a wilcoxon signed rank test (sample order must be the same)
s1 <- round(shapiro.test(df1[,6])$p.value,4)
s2 <- round(shapiro.test(df2[,6])$p.value,4)
min <- min(df1[,6],df2[,6])
max <- max(df1[,6],df2[,6])
step <- (max-min)/5
pdf("shannon.pdf")
boxplot(df1[,6],df2[,6],border=c("chocolate1","chartreuse3"),col=c("bisque","darkolivegreen1"),lwd=3,xaxt='n',yaxt='n',frame=FALSE)
# axis(1, at=1:2,labels=c("GG","Silva"),cex.axis=2)
axis(1, at=1.5,labels=paste("p-values SW1:",s1,"/ SW2:",s2,"/ t:",t,"/ w:",w),cex.axis=1.5)
axis(2,las=1,at=seq(min,max,step), labels=round(seq(min,max,step),digits=1),cex.axis=3)
dev.off()
