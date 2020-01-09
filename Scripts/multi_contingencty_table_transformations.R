# Started on 2019-09-18
# by Rodrigo García-López for Adrian Ochoa-Leyva's Metagenomics and Metatranscriptomics Lab at IBt, UNAM, Cuernavaca, Mexico.
# Under GNU GPLv3 license
# Disclamer: This script was written for use with specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.
# This script was tested with R 3.6.0 - "Planting of a Tree"
# This script reads a contingency table (ideally, absolute abundance) in tsv format, normally the output of converting from biom. Items are in the first column and no taxonomy is found in the last column.
# The objective is to produce from an absolute abundance table; 
# 1.- Total Sum Scaling (TSS) table (in other words, relative abundance table). This does not consider differential sampling depth nor compositionality issues so beware when using it.
# 2.- Rarefied table with <rar_num> resamplings (otherwise called Monte Carlo permutations; rar_num is provided as an argument when running the script). Requires "vegan" library. Reduces bias due to sampling (smaller sampling requires a larger number or permutations). Depth is set to emulate the smallest sample size by default but can be editted manually by changing the value for the object min-limit.
# 3.- Cumulative Sum Scaling table (quantile-scaled based on the cross-table information; quantile passed as argument). Total ratios remain equal but observation are scaled within each sample (alters totals but % remain unaltered). Reduces depth issues and magnifies differential observations. Very sensitive to extremely uneven depth. Requires the "metagenomeSeq" library.
# 4.- Centered-Log-Ratio table (ILR). Geometic mean transformation followed by log transformation (leverages for compositionality issues) requires adjustment of NULL values (0s). Requires library mixOmics from bioconductor
# 5.- Isometric Log Ratio (ILR). Preserves metric properties for geometric elements and avoids the orthogonality bias of the CLR transformation. Difficult to interpret biologically afterwards (mostly for plotting). Requires library mixOmics from bioconductor.
# The last two require a TSS table to function.
# The ilr-transform maps a composition in the D-part Aitchison-simplex isometrically to a D-1 dimensonal euclidian vector. The data can then be analysed in this transformation by all classical multivariate analysis tools. However the interpretation of the results may be difficult, since there is no one-to-one relation between the original parts and the transformed variables.
# The tables can be used for downstream analyses, namely distance matrices calculations to create PCoA and NMDS (for metric and rank-based ordination methods). CLR and ILR are best fitted for PCA (instead of PCoA) as they produce euclidean distances.
# 1 line must be skipped if "# Constructed from biom table" from biom-converted tables is present. set this accordingly at read.table parameters.

# Run as follows:
# cat table.tsv|~/bin/R-3.6.0/bin/Rscript multi_contingencty_table_transformations.R <rar_num> <CSS_quantile> <out_name_prefix>

# Tested with command:
# cat /home/rodrigo/Shrimp_2015-2018/V3-V4_comparison/10_Reduced_sparsity_tables/OTUs_97/12_split_region_tables/06_Final_mix_gg97-table_lvlOTU-V3.tsv|sed "s/onstructed from biom file//"|grep -v "[a-z]"|sed "s/^[V34]*\.//"|~/bin/R-3.6.0/bin/Rscript multi_contingencty_table_transformations.R 10 0.5 test
# In R:
# df <- read.table("10_Reduced_sparsity_tables/OTUs_97/12_split_region_tables/06_Final_mix_gg97-table_lvlOTU-V3.tsv", sep="\t",header=T, skip=1, comment.char='',quote="",fill=F)


rar_num <- 10 # Set the number of random sampling iterations to carry out (rarefactions)
CSS_quantile <- 0.5 # Set the quantiles for the CSS normaliztion (default: median)
prefix <- "test"

args <- commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<3) { # at least, three arguments should be included: <rar_num> <CSS_quantile> <out_name_prefix>
  stop("A minimum of 3 arguments is mandatory: cat table.tsv|Rscript multi_contingencty_table_transformations.R <rar_num> <CSS_quantile> <output_file>", call.=FALSE)
}
rar_num <- as.numeric(args[1]) # Get argument one: the total number of rarefactions to carry out
CSS_quantile <- as.numeric(args[2])
prefix <- as.character(args[3]) # Get a string handle to create output names

df <- read.table(file('stdin'), sep="\t",header=T, skip=1, comment.char='',quote="",fill=F)
df <- rowsum(df[2:length(df)], group=df[[1]]) #sum duplicate items by feature (row) name (if any)
df <- as.matrix(df[order(rowSums(df),decreasing=T),]) # sort input by most abundant first and convert to matrix type object
df <- df[,colSums(df)>0] # NEW filter 2019-09-18: remove those with empty columns (when filtered from large tables)
write.table(df,paste(prefix,"derep_nonempty.tsv", sep="-"), sep="\t", quote=FALSE, col.names=NA) # Output the original matrix with no empty columns nor repeated items

### 1.1.- TSS table (relative)
tss <- t(apply(df,1, function(x) x/colSums(df))) # Create a relative table by dividing each row by the total sample sum
# rel_table <- sweep(df, 2, colSums(df), '/') #Same thing, probably faster
write.table(tss,paste(prefix,"TSS.tsv", sep="-"), sep="\t", quote=FALSE, col.names=NA) # Output the relative frequence matrix

### 1.2.- TSS + CLR table
# CLR and ILR require relative abundance table (any sets summing 1)
library("mixOmics") # This requires mixOmics library found in bioconductor
# Logaritmic transformation is undeterminded for 0 so we need to add an offset. It is common to use 1 but this can only be done if no singletons are presents as otherwise we'd be introducing artificial counts. I agree with  the Peer Bork lab, which has published it is preferable to use a very small value far from the smallest that is actually present in the table.
offset=sort(tss)[sort(tss)>0][1]*0.0000001 #The value is set to a very small fraction of the smallest item. This value is added to the matrix to transform zeroes.
tss_clr <- t(logratio.transfo(t(tss), logratio = "CLR", offset = offset)) # apply CLR transformation (features as columns, samples as rows).
write.table(tss_clr,paste(prefix,"TSS_CLR.tsv", sep="-"), sep="\t", quote=FALSE, col.names=NA) # export the new log-tranformed version of the table

### 1.3.- TSS + ILR table 
offset=sort(tss)[sort(tss)>0][1]*0.0000001
ilr <- t(logratio.transfo(t(tss), logratio = "ILR", offset = offset)) # apply ILR transformation (features as columns, samples as rows).
write.table(ilr,paste(prefix,"TSS_ILR.tsv", sep="-"), sep="\t", quote=FALSE, col.names=NA) # export the new log-transformed table.
# IMPORTANT: The ilr-transform maps a composition in the D-part Aitchison-simplex isometrically to a D-1 dimensonal euclidian vector. The data can then be analysed in this transformation by all classical multivariate analysis tools. However the interpretation of the results may be difficult, since there is no one-to-one relation between the original parts and the transformed variables. Each item in the new table represents a combination and we cannot trace the original samples nor the features. I just included this as many experts in the area have pointed out this may be better than clr, which may be true but no real interpretation seems possible

### 2.1- CSS table
library("metagenomeSeq")
df.metagenomeSeq = newMRexperiment(df, featureData=NULL, libSize=NULL, normFactors=NULL)
df.metagenomeSeq.css <- cumNorm(df.metagenomeSeq, p=CSS_quantile) # Normalize by using a defined quantile (median by default) as scaling parameter
css <- MRcounts(df.metagenomeSeq.css, norm=TRUE, log=F) # Create contingency table based on the new normalized values. This technique partially accounts for compositionality biases too and may be log transformed if required but I will not do that for correlation as this alters linear correlations.
write.table(css,paste(prefix,"CSS",CSS_quantile,"quant.tsv", sep="-"), sep="\t", quote=FALSE, col.names=NA) # export the new normalized version of the table
# Peer Bork tested the original CSS paper for reproducibility and recommends an additional log transformation after the CSS normalization to reach optimal results.
# Doing a TSS after this is not useful as it results in the same table because proportions are kept intact.

### 2.2.- CSS + log
# Using the result of the CSS normalization on the table we can now log-transform but we need to add an offset values because of the zeroes.
css_log <- log(css+1) #For log transformation based on absolute values, we can add 1to all values to avoid zero problems
write.table(css_log,paste(prefix,"CSS_log",CSS_quantile,"quant.tsv", sep="-"), sep="\t", quote=FALSE, col.names=NA) # Save the multi-tranformed table

### 3.1.- Rarefied table - also known as Monte Carlo permutations (total number of iterations set as parameters)
library("vegan")
min_limit <- min(colSums(df)) # The depth parameter for the rarefaction is set as the depth of the smallest sample
# min_limit <- 7000 # optionally, set the min_limit manually
dfr <- df # Create a copy of current dataframe to get the right dimensions
dfr[dfr!=0] <- 0 # Then empty the newly created matrix

for (i in 1:rar_num){ #Create n tables
	print(i) # print current iteration for tracing purposes when running large number of rarefactions
	rar <- rrarefy(t(df),min_limit)
	dfr <- dfr+t(rar)
}
dfr <- dfr/rar_num # From this point forward, we'll work with the mean values for the rarefied table 

write.table(dfr,paste(prefix,"rar",rar_num,"perm.tsv", sep="-"), sep="\t", quote=FALSE, col.names=NA) # Output the mean of all rarefied matrices

### 3.2.- Rarefied + TSS table
dfr_tss <- t(apply(dfr,1, function(x) x/colSums(dfr))) # Using the result of the montecarlo permutations on the table we can now use relative abundances (TSS)
write.table(dfr_tss,paste(prefix,"rar_TSS",rar_num,"perm.tsv", sep="-"), sep="\t", quote=FALSE, col.names=NA) # Output the multi-transformation table

### 3.3.- Rarefied + TSS table + CLR
offset=sort(dfr_tss)[sort(dfr_tss)>0][1]*0.0000001 #The value is set to a very small fraction of the smallest item. This value is added to the matrix to transform zeroes.
dfr_tss_clr <- t(logratio.transfo(t(dfr_tss), logratio = "CLR", offset = offset)) # apply CLR transformation (features as columns, samples as rows).
write.table(dfr_tss,paste(prefix,"rar_TSS_CLR",rar_num,"perm.tsv", sep="-"), sep="\t", quote=FALSE, col.names=NA) # Output the multi-transformation table








