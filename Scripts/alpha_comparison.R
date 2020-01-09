# Started on 2019-09-10
# by Rodrigo García-López for Adrian Ochoa-Leyva's Metagenomics and Metatranscriptomics Lab at IBt, UNAM, Cuernavaca, Mexico.
# Under GNU GPLv3 license
# Disclamer: This script was written for use with specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.
# UPDATE: 2019-09-15 - Colored/separated by region (will now only work correctly for the 06_Final_mix sets and also, Pielou's index cannot be calculated for lvl1
# This script was tested with R 3.6.0 - "Planting of a Tree"
# This script reads an otu table (output from converting biom tables from qiime2). Items are in the first column and no taxonomy is found in the last column.
# Alpha diversity indices are calculated for each rarefied table and appended by sample, creating a matrix for each taxonomy index.
# Comparisons between groups are plotted after the generation of the sampled data (by default, organ and region).
# Run as follows:
# cat table.tsv|~/bin/R-3.6.0/bin/Rscript alpha_comparison.R <rar_num> <prefix>

# Tested with command:
# cat /home/rodrigo/Shrimp_2015-2018/V3-V4_comparison/10_Reduced_sparsity_tables/OTUs_97/gg/06_Final_mix_gg97-table_lvl5.tsv|~/bin/R-3.6.0/bin/Rscript alpha_comparison.R 100 test
# In R:
# df <- read.table("10_Reduced_sparsity_tables/OTUs_97/gg/06_Final_mix_gg97-table_lvl5.tsv", sep="\t",header=T, skip=1, comment.char='',quote="",fill=F)

rar_num <- 100 # Set the number of random sampling iterations to carry out (rarefactions)
args <- commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<2) { # at least, three arguments should be included: <min_freq> <min_samples> <#_samples_table_1>
  stop("A minimum of 2 arguments is mandatory: cat table.tsv|Rscript Sparsity_reduction_merged_tables.R <rar_num> <output_file>", call.=FALSE)
}
rar_num <- as.numeric(args[1]) # Get argument one: the total number of rarefactions to carry out
prefix <- as.character(args[2]) # Get a string handle to create output names

df <- read.table(file('stdin'), sep="\t",header=T, skip=1, comment.char='',quote="",fill=F)
df <- rowsum(df[2:length(df)], group=df[[1]]) #sum duplicate items by feature (row) name (if any)
df <- as.matrix(df[order(rowSums(df),decreasing=T),]) # sort input by most abundant first and convert to matrix type object
samples <- colnames(df)

# Rarefy table
library("vegan") # we will create a rarefied set with vegan using the smallest sample size as depth with a certain number of repetitions, to get the average so we can approximate the actual composition more accurately.
min_limit <- min(colSums(df)) # The depth parameter for the rarefaction is set as the depth of the smallest sample
# min_limit <- 7000 # optionally, set the min_limit manually
dfr <- df # Create a copy of current dataframe to get the right dimensions
dfr[dfr!=0] <- 0 # Then empty the newly created matrix
m_shannon <- data.frame(matrix(NA,nrow=rar_num,ncol=ncol(df))) #Initialize empty matrices for each index
colnames(m_shannon) <- samples
m_diversity <- m_shannon # This gets the non-log form of shannon
m_pielou <- m_shannon # This is for the distribution among totals
m_simpson <- m_shannon # This is for the dominance
m_observed <- m_shannon # This is the raw number of unique features
m_chao1 <- m_shannon # This is for the estimated items based on singletons (adjusted if not present)
m_ace <- m_shannon # This is for the estimated items but based on singletons, doubletons, up to 10s, similar to chao

# Now, analyze alpha diversity for each resample
for (i in 1:rar_num){ #Create n tables
	print(i) # print current iteration for tracing purposes when running large number of rarefactions
	rar <- rrarefy(t(df),min_limit)
	m_shannon[i,] <- diversity(rar,index="shannon") # Store Shannon's entropy index (uncertainty of predicting next draw)
	m_diversity[i,] <-exp(m_shannon) # Store the effective diversity (non-logartithmic version of shannon, true diversity)
	SUM <- apply(df>0,2,sum) # Store the sum of non null observations
	m_pielou[i,] <- m_shannon/log(SUM)# Store Pielou's evenness (distribution of observations among items (less means more abundance variation; sensitive to sample size))
	m_simpson[i,] <- diversity(rar,index="simpson") # Store the Simpson's measure of dominance (P or drawing two items from same features)
	Richness <- estimateR(rar) # Get the observed features, the chao1 richness, the ACE richness, estimators along with their se
	m_observed[i,] <- Richness[1,] # Store unique feature counts
	m_chao1[i,] <- Richness[2,] # Store the number of estimated features
	m_ace[i,] <- Richness[2,] # Store the number of estimated features based on the those having 1 through 10
}

# # Output all tables
# write.table(m_shannon,paste(prefix,"shannon",rar_num,"resamples.tsv", sep="_"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE) # and export it as a tsv file
# write.table(m_diversity,paste(prefix,"diversity",rar_num,"resamples.tsv", sep="_"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
# write.table(m_pielou,paste(prefix,"pielou",rar_num,"resamples.tsv", sep="_"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
# write.table(m_simpson,paste(prefix,"simpson",rar_num,"resamples.tsv", sep="_"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
# write.table(m_observed,paste(prefix,"observed",rar_num,"resamples.tsv", sep="_"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
# write.table(m_chao1,paste(prefix,"chao1",rar_num,"resamples.tsv", sep="_"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
# write.table(m_ace,paste(prefix,"ace",rar_num,"resamples.tsv", sep="_"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# Define the groups and create plot parameters
group_V3 <- grep("^V3\\.",samples) # save positions of samples from each group
group_V4 <- grep("^V4\\.",samples)
group_V3V4 <- grep("^V3V4\\.",samples)
colors <- c(rep("coral1",length(group_V3)),rep("cornflowerblue",length(group_V4)),rep("turquoise3",length(group_V3V4)))
colors2 <- c(rep("bisque",length(group_V3)),rep("darkslategray2",length(group_V4)),rep("aquamarine",length(group_V3V4)))
alt_colors <- colors[c(1:4,9:12,17:20,5:8,13:16,21:24)]
alt_colors2 <- colors2[c(1:4,9:12,17:20,5:8,13:16,21:24)]
matrices <- c("m_observed", "m_chao1", "m_ace", "m_shannon", "m_diversity", "m_pielou", "m_simpson")
m_names <- c("Observed features", "Chao1 richness estimator", "ACE richness estimator", "Shannon's entropy index (log)", "Actual diversity (exp from Shannon)", "Pielou's evenness", "Simpson dominance")

collated <- list("H-V3"=1:4,"H-V4"=5:8,"H-V3V4"=9:12,"I-V3"=13:16,"I-V4"=17:20,"I-V3V4"=21:24) # This vector will be used for the evalue comparisons
alphas <- data.frame(matrix(NA,nrow=7,ncol=40)) #Create empty matrix for the results
names(alphas) <- c("Index",colnames(df),apply(t(combn(c("H-V3","H-V4","H-V3V4","I-V3","I-V4","I-V3V4"),2)), 1, function(x) paste(x[1],"vs",x[2],sep="_"))) # Append the names


# Now the boxplots for all samples
for (i in 1:7) {
	# First, ordered by region
	pdf(paste(prefix,matrices[i],rar_num,"resamples.pdf", sep="_"))
	mat <- get(matrices[i])
	name <- m_names[i]
	y <- seq(min(mat),max(mat),round((max(mat)-min(mat))/20,2))
	boxplot(mat,border=as.character(colors),col=as.character(colors2),xaxt='n',yaxt='n',main=name,frame=FALSE,outline=FALSE)
	mtext("Ordered by Region")
	axis(1,at=seq(1,24,1),las=2, cex.axis=0.8, labels=samples)
	axis(2,at=y,las=1,cex.axis=0.8, labels=format(y,nsmall=2))
	abline(v=seq(0.5,100,1),lty=2, col="gray",lwd=0.5)
	abline(v=seq(8.5,20,8),lwd=2)
	abline(v=seq(4.5,24,8),col="darkgray",lwd=2)
	abline(v=seq(4.5,24,8),col="white",lwd=1.75,lty=3)
	
	# Then by organ
	mat <- mat[,c(1:4,9:12,17:20,5:8,13:16,21:24)]
	boxplot(mat,border=as.character(alt_colors),col=as.character(alt_colors2),xaxt='n',yaxt='n',main=name,frame=FALSE,outline=FALSE)
	mtext("Ordered by Organ")
	axis(1,at=seq(1,24,1),las=2, cex.axis=0.8, labels=colnames(mat))
	axis(2,at=y,las=1,cex.axis=0.8, labels=format(y,nsmall=2))
	abline(v=seq(0.5,100,1),lty=2, col="gray",lwd=0.5)
	abline(v=12.5,lwd=2)
	abline(v=c(4.5,8.5,16.5,20.5),col="darkgray",lwd=2)
	abline(v=c(4.5,8.5,16.5,20.5),col="white",lwd=1.75,lty=3)
	
	# Append summarized q.values (sampling-corrected) for each group test (mann-whitney-wilcoxon non-parametric independent singned rank sums)
	mat <- cbind(mat,"H-V3_vs_H-V4"=p.adjust(apply(mat, 1, function(x) wilcox.test(x[collated[[1]]], x[collated[[2]]],paired=FALSE)$p.value),method="fdr"))
	mat <- cbind(mat,"H-V3_vs_H-V3V4"=p.adjust(apply(mat, 1, function(x) wilcox.test(x[collated[[1]]], x[collated[[3]]],paired=FALSE)$p.value),method="fdr"))
	mat <- cbind(mat,"H-V3_vs_I-V3"=p.adjust(apply(mat, 1, function(x) wilcox.test(x[collated[[1]]], x[collated[[4]]],paired=FALSE)$p.value),method="fdr"))
	mat <- cbind(mat,"H-V3_vs_I-V4"=p.adjust(apply(mat, 1, function(x) wilcox.test(x[collated[[1]]], x[collated[[5]]],paired=FALSE)$p.value),method="fdr"))
	mat <- cbind(mat,"H-V3_vs_I-V3V4"=p.adjust(apply(mat, 1, function(x) wilcox.test(x[collated[[1]]], x[collated[[6]]],paired=FALSE)$p.value),method="fdr"))
	mat <- cbind(mat,"H-V4_vs_H-V3V4"=p.adjust(apply(mat, 1, function(x) wilcox.test(x[collated[[2]]], x[collated[[3]]],paired=FALSE)$p.value),method="fdr"))
	mat <- cbind(mat,"H-V4_vs_I-V3"=p.adjust(apply(mat, 1, function(x) wilcox.test(x[collated[[2]]], x[collated[[4]]],paired=FALSE)$p.value),method="fdr"))
	mat <- cbind(mat,"H-V4_vs_I-V4"=p.adjust(apply(mat, 1, function(x) wilcox.test(x[collated[[2]]], x[collated[[5]]],paired=FALSE)$p.value),method="fdr"))
	mat <- cbind(mat,"H-V4_vs_I-V3V4"=p.adjust(apply(mat, 1, function(x) wilcox.test(x[collated[[2]]], x[collated[[6]]],paired=FALSE)$p.value),method="fdr"))
	mat <- cbind(mat,"H-V3V4_vs_I-V3"=p.adjust(apply(mat, 1, function(x) wilcox.test(x[collated[[3]]], x[collated[[4]]],paired=FALSE)$p.value),method="fdr"))
	mat <- cbind(mat,"H-V3V4_vs_I-V4"=p.adjust(apply(mat, 1, function(x) wilcox.test(x[collated[[3]]], x[collated[[5]]],paired=FALSE)$p.value),method="fdr"))
	mat <- cbind(mat,"H-V3V4_vs_I-V3V4"=p.adjust(apply(mat, 1, function(x) wilcox.test(x[collated[[3]]], x[collated[[6]]],paired=FALSE)$p.value),method="fdr"))
	mat <- cbind(mat,"I-V3_vs_I-V4"=p.adjust(apply(mat, 1, function(x) wilcox.test(x[collated[[4]]], x[collated[[5]]],paired=FALSE)$p.value),method="fdr"))
	mat <- cbind(mat,"I-V3_vs_I-V3V4"=p.adjust(apply(mat, 1, function(x) wilcox.test(x[collated[[4]]], x[collated[[6]]],paired=FALSE)$p.value),method="fdr"))
	mat <- cbind(mat,"I-V4_vs_I-V3V4"=p.adjust(apply(mat, 1, function(x) wilcox.test(x[collated[[5]]], x[collated[[6]]],paired=FALSE)$p.value),method="fdr"))
	alphas[i,] <- c(name,colMeans(mat))
# 	header <- c("Index",colnames(mat))
# 	names(alphas) <- header
	# And once more with samples joined by group (region/organ)
	mat <- cbind("H-V3"=apply(mat[,1:4],1,mean),"I-V3"=apply(mat[,13:16],1,mean),"H-V4"=apply(mat[,5:8],1,mean),"I-V4"=apply(mat[,17:20],1,mean),"H-V3V4"=apply(mat[,9:12],1,mean),"I-V3V4"=apply(mat[,21:24],1,mean))
	boxplot(mat,yaxt='n',main=name, border=c("coral1","coral1","cornflowerblue","cornflowerblue","turquoise3","turquoise3"),col=c("bisque","bisque","darkslategray2","darkslategray2","aquamarine","aquamarine"),frame=FALSE,outline=FALSE)
	axis(2,at=y,las=1,cex.axis=0.8, labels=format(y,nsmall=2))
	abline(v=seq(0.5,100,1),lty=2, col="gray",lwd=0.5)
	abline(v=c(2.5,4.5),lwd=2)
	dev.off()
}

write.table(alphas,paste(prefix,rar_num,"alphas.tsv", sep="_"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE) # and export it as a tsv file

# Deprecated extra table for comparing pvalues. No longer used as larger Ns will derive in all p-values being <alpha
# # # # # combo <- t(combn(colnames(mat), 2)) # Create permutations per sample
# # # # # pval <- apply(combo, 1, function(x) wilcox.test(mat[,x[1]], mat[,x[2]])$p.value) # Calculate p.value
# # # # # comp <- data.frame(combo, pval) # combine into long table
# # # # # comp <- comp[order(comp[,3], decreasing=FALSE),] # order the permutations by the p.value
# # # # # write.table(comp, "test.tsv", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
# # # # # 
# # # # # apply(mat, 1, function(x) wilcox.test(x[1:4], x[5:8],paired=FALSE)$p.value)
