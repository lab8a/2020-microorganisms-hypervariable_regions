# Started on 2019-09-11
# by Rodrigo García-López for Adrian Ochoa-Leyva's Metagenomics and Metatranscriptomics Lab at IBt, UNAM, Cuernavaca, Mexico.
# Under GNU GPLv3 license
# Disclamer: This script was written for use with specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.
# UPDATE: 2019-09-12 - The script was heavily modified to carry out resampling of the original table so that multiple testing could be evaluated for corrected p-values. Trials were carried out with the sparsity-filtered tables. Also, an additional parameter was added for the rarefaction number and the skip first line flag was activated accordingly
# This script was tested with R 3.6.0 - "Planting of a Tree"
# This script reads an otu table in tsv format. Items are in the first column and no taxonomy is found in the last column.
# Beta diversity indices are calculated for creating distance/dissimilarity matrices, then doing a non-metric multidimensional scaling (NMDS)
# 1 line must be skipped if "# Constructed from biom table" from biom-converted tables is present. set this accordingly at read.table parameters.

# Run as follows:
# cat table.tsv|~/bin/R-3.6.0/bin/Rscript beta_diversity.R <rar_num> <out_name_prefix>

# Tested with command:
# cat /home/rodrigo/Shrimp_2015-2018/V3-V4_comparison/10_Reduced_sparsity_tables/OTUs_97/gg/06_Final_mix_gg97-table_lvl5.tsv|~/bin/R-3.6.0/bin/Rscript beta_diversity.R 10 test
# In R:
# df <- read.table("10_Reduced_sparsity_tables/OTUs_97/gg/06_Final_mix_gg97-table_lvl5.tsv", sep="\t",header=T, skip=1, comment.char='',quote="",fill=F)


rar_num <- 10000 # Set the number of random sampling iterations to carry out (rarefactions)
args <- commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<2) { # at least, three arguments should be included: <min_freq> <min_samples> <#_samples_table_1>
  stop("A minimum of 2 arguments is mandatory: cat table.tsv|Rscript beta_diversity.R <rar_num> <output_file>", call.=FALSE)
}
rar_num <- as.numeric(args[1]) # Get argument one: the total number of rarefactions to carry out
prefix <- as.character(args[2]) # Get a string handle to create output names

df <- read.table(file('stdin'), sep="\t",header=T, skip=1, comment.char='',quote="",fill=F)
df <- rowsum(df[2:length(df)], group=df[[1]]) #sum duplicate items by feature (row) name (if any)
df <- as.matrix(df[order(rowSums(df),decreasing=T),]) # sort input by most abundant first and convert to matrix type object

 ### Define groups and set parameters ##
samples <- colnames(df) # also store sample names
unique <- apply(df>0,2,sum) # Get a vector of unique features per sample
reg_org <- region <- organ <- NULL # Initialize empty vector for organ and region
region[grep("^V3\\.", samples)]="coral1" #Populate vectors with the groups
region[grep("^V4\\.", samples)]="cornflowerblue"
region[grep("^V3V4\\.", samples)]="turquoise3"
organ[grep("\\.H", samples)]="darkorchid1"
organ[grep("\\.I", samples)]="chartreuse2"
reg_org[grep("^V3\\.H", samples)]="coral1"
reg_org[grep("^V4\\.H", samples)]="cornflowerblue"
reg_org[grep("^V3V4\\.H", samples)]="turquoise3"
reg_org[grep("^V3\\.I", samples)]="darkred"
reg_org[grep("^V4\\.I", samples)]="darkslategray4"
reg_org[grep("^V3V4\\.I", samples)]="darkolivegreen"

legend_items <- c("coral1","darkred","cornflowerblue","darkslategray4","turquoise3","darkolivegreen")
colors3 <- levels(factor(region))
colors2 <- levels(factor(organ))
colors6 <- levels(factor(reg_org))
# library("betapart") # DEPRECATED
reg_org1 <- region1 <- organ1 <- NULL # Initialize empty vector for organ and region
region1[grep("^V3\\.", samples)]="V3" #Populate vectors with the groups
region1[grep("^V4\\.", samples)]="V4"
region1[grep("^V3V4\\.", samples)]="V3V4"
organ1[grep("\\.H", samples)]="Hepatopancreas"
organ1[grep("\\.I", samples)]="Gut"
reg_org1[grep("^V3\\.H", samples)]="V3-H"
reg_org1[grep("^V4\\.H", samples)]="V4-H"
reg_org1[grep("^V3V4\\.H", samples)]="V3V4-H"
reg_org1[grep("^V3\\.I", samples)]="V3-I"
reg_org1[grep("^V4\\.I", samples)]="V4-I"
reg_org1[grep("^V3V4\\.I", samples)]="V3V4-I"

 ### Set a cross sample test function to analyze all groups
library("vegan") # Distance matrix calculations are included the vegan package
pairwise.adonis <- function(x,factors, sim.method, p.adjust.m) { #Create a function for pairwise comparisons by group; based on code found in https://github.com/pmartinezarbizu/pairwiseAdonis)
	combo = as.matrix(combn(unique(factors),2))
	Comparison = c()
	F.Model =c()
	R2 = c()
	p.value = c()
	for(item in 1:ncol(combo)){
		ad = adonis(x[factors %in% c(as.character(combo[1,item]),as.character(combo[2,item])),] ~
		factors[factors %in% c(as.character(combo[1,item]),as.character(combo[2,item]))] , method =sim.method);
		Comparison = c(Comparison,paste(combo[1,item],'vs',combo[2,item]));
		F.Model =c(F.Model,ad$aov.tab[1,4]);
		R2 = c(R2,ad$aov.tab[1,5]);
		p.value = c(p.value,ad$aov.tab[1,6])
	}
# 	p.adjusted = p.adjust(p.value,method=p.adjust.m)
# 	pairw.res = data.frame(Comparison,F.Model,R2,p.value,p.adjusted)
	pairw.res = data.frame(Comparison,F.Model,R2,p.value)
	return(pairw.res)
}


### Carry out multiple rarefactions also known as Monte Carlo permutations (total number of iterations set as parameters)
min_limit <- min(colSums(df)) # The depth parameter for the rarefaction is set as the depth of the smallest sample
# min_limit <- 7000 # optionally, set the min_limit manually
dfr <- df # Create a copy of current dataframe to get the right dimensions
dfr[dfr!=0] <- 0 # Then empty the newly created matrix
temp=NULL
ad_j_org <- pairwise.adonis(t(df), as.factor(organ1),"jaccard","fdr");ad_j_org[2:4]=0 # Create empty matrix for summary
ad_j_reg <- pairwise.adonis(t(df), as.factor(region1),"jaccard","fdr");ad_j_reg[2:4]=0
ad_j_reg_org <- pairwise.adonis(t(df), as.factor(reg_org1),"jaccard","fdr");ad_j_reg_org[2:4]=0
ad_b_org <- pairwise.adonis(t(df), as.factor(organ1),"bray","fdr");ad_b_org[2:4]=0
ad_b_reg <- pairwise.adonis(t(df), as.factor(region1),"bray","fdr");ad_b_reg[2:4]=0
ad_b_reg_org <- pairwise.adonis(t(df), as.factor(reg_org1),"bray","fdr");ad_b_reg_org[2:4]=0
p_j_org <- data.frame(matrix(NA,nrow=rar_num,ncol=nrow(ad_j_org)));names(p_j_org) <- ad_j_org[,1] # Create a matrix for stored pvalues for each set
p_j_reg <- data.frame(matrix(NA,nrow=rar_num,ncol=nrow(ad_j_reg)));names(p_j_reg) <- ad_j_reg[,1]
p_j_reg_org <- data.frame(matrix(NA,nrow=rar_num,ncol=nrow(ad_j_reg_org)));names(p_j_reg_org) <- ad_j_reg_org[,1]
p_b_org <- data.frame(matrix(NA,nrow=rar_num,ncol=nrow(ad_b_org)));names(p_b_org) <- ad_b_org[,1]
p_b_reg <- data.frame(matrix(NA,nrow=rar_num,ncol=nrow(ad_b_reg)));names(p_b_reg) <- ad_b_reg[,1]
p_b_reg_org <- data.frame(matrix(NA,nrow=rar_num,ncol=nrow(ad_b_reg_org)));names(p_b_reg_org) <- ad_b_reg_org[,1]
all_j_reg <- all_b_reg <- all_j_reg_org <- all_b_reg_org <- c(0,0,0) # These are for storing the total group adonis' results
p_all_j_reg <- p_all_b_reg <- p_all_j_reg_org <- p_all_b_reg_org <- data.frame(matrix(NA,nrow=rar_num,ncol=1)) # This are for the associated pvalues

for (i in 1:rar_num){ #Create n tables
	print(i) # print current iteration for tracing purposes when running large number of rarefactions
	rar <- rrarefy(t(df),min_limit)
	temp <- pairwise.adonis(rar, as.factor(organ1),"jaccard","fdr")[,2:4] # For each pairwise comparison, do a permanova test
	ad_j_org[,2:4] <- ad_j_org[,2:4]+temp #Save the Statistic, R2 and pvalue
	p_j_org[i,] <- unlist(temp[length(temp)]) #and store the p-value
	temp <- pairwise.adonis(rar, as.factor(region1),"jaccard","fdr")[,2:4]
	ad_j_reg[,2:4] <- ad_j_reg[,2:4]+temp
	p_j_reg[i,] <- unlist(temp[length(temp)])
	temp <- pairwise.adonis(rar, as.factor(reg_org1),"jaccard","fdr")[,2:4]
	ad_j_reg_org[,2:4] <- ad_j_reg_org[,2:4]+temp
	p_j_reg_org[i,] <- unlist(temp[length(temp)])
	temp <- pairwise.adonis(rar, as.factor(organ1),"bray","fdr")[,2:4] # For each pa
	ad_b_org[,2:4] <- ad_b_org[,2:4]+temp #Save the Statistic, R2 and pvalue
	p_b_org[i,] <- unlist(temp[length(temp)]) #and store the p-value
	temp <- pairwise.adonis(rar, as.factor(region1),"bray","fdr")[,2:4]
	ad_b_reg[,2:4] <- ad_b_reg[,2:4]+temp
	p_b_reg[i,] <- unlist(temp[length(temp)])
	temp <- pairwise.adonis(rar, as.factor(reg_org1),"bray","fdr")[,2:4]
	ad_b_reg_org[,2:4] <- ad_b_reg_org[,2:4]+temp
	p_b_reg_org[i,] <- unlist(temp[length(temp)])
	jaccard <- vegdist(rar,method="jaccard",upper=TRUE,diag=TRUE) #create distance matrix, and calculate adonis for each group (cross-table)
	temp <- adonis(as.matrix(jaccard) ~ as.factor(region),method=jaccard)$aov.tab[1,4:6]
	all_j_reg <- all_j_reg+temp; p_all_j_reg[i,1] <- temp[3]
	temp <- adonis(as.matrix(jaccard) ~ as.factor(reg_org),method=jaccard)$aov.tab[1,4:6]
	all_j_reg_org <- all_j_reg_org+temp; p_all_j_reg_org[i,1] <- temp[3]
	bray <- vegdist(rar,method="bray",upper=TRUE,diag=TRUE)
	temp <- adonis(as.matrix(bray) ~ as.factor(region),method=bray)$aov.tab[1,4:6]
	all_b_reg <- all_b_reg+temp; p_all_b_reg[i,1] <- temp[3]
	temp <- adonis(as.matrix(bray) ~ as.factor(reg_org),method=bray)$aov.tab[1,4:6]
	all_b_reg_org <- all_b_reg_org+temp; p_all_b_reg_org[i,1] <- temp[3]
	dfr <- dfr+t(rar)
}
df <- dfr/rar_num # From this point forward, we'll work with the mean values for the rarefied table 
ad_j_org[,2:4] <- ad_j_org[2:4]/rar_num # Also get the mean of each statistic and p.value
ad_j_reg[,2:4] <- ad_j_reg[2:4]/rar_num
ad_j_reg_org[,2:4] <- ad_j_reg_org[2:4]/rar_num
ad_b_org[,2:4] <- ad_b_org[2:4]/rar_num
ad_b_reg[,2:4] <- ad_b_reg[2:4]/rar_num
ad_b_reg_org[,2:4] <- ad_b_reg_org[2:4]/rar_num
all_j_reg <- cbind("Region",all_j_reg/rar_num,mean(p.adjust(unlist(p_all_j_reg),"fdr")))
all_j_reg_org <- cbind("Region",all_j_reg_org/rar_num,mean(p.adjust(unlist(p_all_j_reg_org),"fdr")))
all_b_reg <- cbind("Region-organ",all_b_reg/rar_num,mean(p.adjust(unlist(p_all_b_reg),"fdr")))
all_b_reg_org <- cbind("Region-organ",all_b_reg_org/rar_num,mean(p.adjust(unlist(p_all_b_reg_org),"fdr")))

# And append an adjusted pvalue (the mean of the adjusted pvalue vector, obtained for all iterations)
ad_j_org <- cbind(ad_j_org, "q.value"=round(apply(apply(p_j_org, 2, function(x) p.adjust(x,method="fdr")),2,mean),4))
ad_j_reg <- cbind(ad_j_reg, "q.value"=round(apply(apply(p_j_reg, 2, function(x) p.adjust(x,method="fdr")),2,mean),4))
ad_j_reg_org <- cbind(ad_j_reg_org, "q.value"=round(apply(apply(p_j_reg_org, 2, function(x) p.adjust(x,method="fdr")),2,mean),4))
ad_b_org <- cbind(ad_b_org, "q.value"=round(apply(apply(p_b_org, 2, function(x) p.adjust(x,method="fdr")),2,mean),4))
ad_b_reg <- cbind(ad_b_reg, "q.value"=round(apply(apply(p_b_reg, 2, function(x) p.adjust(x,method="fdr")),2,mean),4))
ad_b_reg_org <- cbind(ad_b_reg_org, "q.value"=round(apply(apply(p_b_reg_org, 2, function(x) p.adjust(x,method="fdr")),2,mean),4))

names(all_j_reg) <- names(all_j_reg_org) <- names(all_b_reg) <- names(all_b_reg_org) <- names(ad_j_reg)
rownames(all_j_reg) <- rownames(all_b_reg) <- "Region-All"
rownames(all_j_reg_org) <- rownames(all_j_reg_org) <- "Region/Organ-All"

write.table(df,paste(prefix,"rarefied_matrix.tsv", sep="-"), sep="\t", quote=FALSE, col.names=NA) # Output the mean of all rarefied matrices
write.table(rbind(ad_j_org,all_j_reg,ad_j_reg,all_j_reg_org,ad_j_reg_org),paste(prefix,"jaccard_statistics.tsv", sep="-"), sep="\t", quote=FALSE, row.names=FALSE) # And the mean of all statistics with adjusted p values
write.table(rbind(ad_b_org,all_b_reg,ad_b_reg,all_b_reg_org,ad_b_reg_org),paste(prefix,"bray-curtis_statistics.tsv", sep="-"), sep="\t", quote=FALSE, row.names=FALSE) # Same thing for Bray-Curtis

#### This should be used for inputting the rarefied matrix ###
jaccard <- vegdist(t(df),method="jaccard",upper=TRUE,diag=TRUE)
bray <- vegdist(t(df),method="bray",upper=TRUE,diag=TRUE)
feature_names <- rownames(df) # Store the feature names
rownames(df) <- paste("f",1:nrow(df),sep="") #Replace the name with a features number
write.table(cbind(rownames(df),feature_names),paste(prefix,"feature_list.tsv", sep="-"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

 ########### We'll start beta div calculation with jaccard (composition) ###########
# anosim_jaccard_org <- anosim(jaccard,distance=jaccard,organ,permutations=10000);anosim_jaccard_org
# adonis_jaccard_org <- adonis(as.matrix(jaccard) ~ as.factor(organ),method=jaccard,permutations=10000);adonis_jaccard_org
# anosim_jaccard_reg <- anosim(jaccard,distance=jaccard,region,permutations=10000);anosim_jaccard_reg
# adonis_jaccard_reg <- adonis(as.matrix(jaccard) ~ as.factor(region),method=jaccard,permutations=10000);adonis_jaccard_reg
# anosim_jaccard_reg_org <- anosim(jaccard,distance=jaccard,reg_org,permutations=10000);anosim_jaccard_reg_org
# adonis_jaccard_reg_org <- adonis(as.matrix(jaccard) ~ as.factor(reg_org),method=jaccard,permutations=10000);adonis_jaccard_reg_org

write.table(as.matrix(jaccard),paste(prefix,"jaccard.tsv", sep="-"), sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
# Get the number of dimensions to use for NMDS (too much or too less stress on two axes is uninformative)
stress <- data.frame(matrix(NA,nrow=15,ncol=1)) # Set a list to get the 
for (i in 1:10){
	NMDS <- metaMDS(jaccard,distance=jaccard,k=i,trymax=100,autotransform=FALSE,wascores = FALSE)
	stress[i,1] <- NMDS$stress
}
dimensions <- which(abs(diff(stress[,1]))<0.01)[1] # To set a suitable number of dimensions, use the vector of stress values for each dimension set. The cutoff is set to less than 0.01. Since the difference vector is shifted by -1, we can use the first item not passing (the number actually points to the last item passing).

# NMDS <- metaMDS(jaccard,distance=jaccard,k=dimensions,trymax=1000,autotransform=FALSE,wascores = FALSE) # This uses a previously calculated distance
NMDS <- metaMDS(t(df),distance="jaccard",k=dimensions,trymax=1000,autotransform=FALSE,wascores = TRUE) # This includes species scores and uses the whole table instead
# stressplot(NMDS)
# plot(NMDS, type="n")
# points(NMDS, display=c("sites"), choices=c(1,2), pch=3, col="red")


 ### Color by regions
pdf(paste(prefix,"jaccard_regions.pdf", sep="-"))
# par(mfrow=c(2,2))
# First plot dimensions 1 and 2
ordiplot(NMDS,type="n", main="Jaccard - non-metric multidimensional scaling", choices=c(1,2))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",round(all_j_reg[3],4),"p =",round(all_j_reg[4],4),"q =",round(all_j_reg_org[5],4)))
# # # orditorp(NMDS,display="species",col="gray",air=0.01, cex=0.5, choices=c(1,2))
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=region, choices=c(1,2))
ordihull(NMDS, groups=region, draw="polygon", col=colors3, border=colors3, label=F, choices=c(1,2))
# ordiellipse(NMDS, groups=region, draw="polygon", display="sites", col=colors3, border=colors3, label=F, choices=c(1,2), kind="ehull")
# ordispider(NMDS,groups=region,display="sites",col=colors3,label=F, choices=c(1,2))
# ordibar(NMDS, groups=region, display="sites", col=colors3, label=F, choices=c(1,2))
# ordicluster(NMDS, cluster=hclust(jaccard), choices=c(1,2)) #Hierachical clusering
# Now dimensions 1 and 3
ordiplot(NMDS,type="n", main="Jaccard - non-metric multidimensional scaling", choices=c(1,3))
mtext(paste("Dimensions 1 and 3 - Adonis R2:",round(all_j_reg[3],4),"p =",round(all_j_reg[4],4),"q =",round(all_j_reg_org[5],4)))
# # # orditorp(NMDS,display="species",col="gray",air=0.01, cex=0.5, choices=c(1,3))
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=region, choices=c(1,3))
ordihull(NMDS, groups=region, draw="polygon", col=colors3, border=colors3, label=F, choices=c(1,3))
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main="Jaccard - non-metric multidimensional scaling", choices=c(2,3))
mtext(paste("Dimensions 2 and 3 - Adonis R2:",round(all_j_reg[3],4),"p =",round(all_j_reg[4],4),"q =",round(all_j_reg[5],4)))
# # # orditorp(NMDS,display="species",col="gray",air=0.01, cex=0.5, choices=c(2,3))
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=region, choices=c(2,3))
ordihull(NMDS, groups=region, draw="polygon", col=colors3, border=colors3, label=F, choices=c(2,3))
# ordiplot(NMDS,type="n",xaxt='n',yaxt='n',frame=F,xlab="",ylab="")

# Now the effect of unique items on the clustering dimensions 1 and 2
ordiplot(NMDS,type="n", main="Landscape of unique items - Jaccard", choices=c(1,2))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",round(all_j_reg[3],4),"p =",round(all_j_reg[4],4),"q =",round(all_j_reg[5],4)))
ordisurf(NMDS, unique, col="gray", choices=c(1,2),add=TRUE)
points(NMDS, col=region, pch=19,choices=c(1,2))
legend("topleft", legend=c("V3","V4","V3V4"), pch=19, col=colors3,cex=1)
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main="Landscape of unique items - Jaccard", choices=c(1,3))
mtext(paste("Dimensions 1 and 3 - Adonis R2:",round(all_j_reg[3],4),"p =",round(all_j_reg[4],4),"q =",round(all_j_reg[5],4)))
ordisurf(NMDS, unique, col="gray", choices=c(1,3),add=TRUE)
points(NMDS, col=region, pch=19,choices=c(1,3))
legend("topleft", legend=c("V3","V4","V3V4"), pch=19, col=colors3,cex=1)
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main="Landscape of unique items - Jaccard", choices=c(2,3))
mtext(paste("Dimensions 2 and 3 - Adonis R2:",round(all_j_reg[3],4),"p =",round(all_j_reg[4],4),"q =",round(all_j_reg[5],4)))
ordisurf(NMDS, unique, col="gray", choices=c(2,3),add=TRUE)
points(NMDS, col=region, pch=19,choices=c(2,3))
legend("topleft", legend=c("V3","V4","V3V4"), pch=19, col=colors3,cex=1)

# Now the hierarchical clustering of unique items on the clustering dimensions 1 and 2
ordiplot(NMDS,type="n", main="Hierachical clustering - Jaccard", choices=c(1,2))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",round(all_j_reg[3],4),"p =",round(all_j_reg[4],4),"q =",round(all_j_reg[5],4)))
ordicluster(NMDS, cluster=hclust(jaccard), choices=c(1,2),col="darkgray") #Hierachical clusering
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=region, choices=c(1,2))
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main="Hierachical clustering - Jaccard", choices=c(1,3))
mtext(paste("Dimensions 1 and 3 - Adonis R2:",round(all_j_reg[3],4),"p =",round(all_j_reg[4],4),"q =",round(all_j_reg[5],4)))
ordicluster(NMDS, cluster=hclust(jaccard), choices=c(1,3),col="darkgray") #Hierachical clusering
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=region, choices=c(1,3))
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main="Hierachical clustering - Jaccard", choices=c(2,3))
mtext(paste("Dimensions 2 and 3 - Adonis R2:",round(all_j_reg[3],4),"p =",round(all_j_reg[4],4),"q =",round(all_j_reg[5],4)))
ordicluster(NMDS, cluster=hclust(jaccard), choices=c(2,3),col="darkgray") #Hierachical clusering
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=region, choices=c(2,3))
dev.off()

 ### Color by organs
pdf(paste(prefix,"jaccard_organs.pdf", sep="-"))
# First dimensions 1 and 2
ordiplot(NMDS,type="n", main="Jaccard - non-metric multidimensional scaling", choices=c(1,2))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",round(ad_j_org[3],4),"p =",round(ad_j_org[4],4),"q =",round(ad_j_org[5],4)))
# # # orditorp(NMDS,display="species",col="gray",air=0.01, cex=0.5, choices=c(1,2))
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=organ, choices=c(1,2))
ordihull(NMDS, groups=organ, draw="polygon", col=colors2, border=colors2, label=F, choices=c(1,2))
# Now dimensions 1 and 3
ordiplot(NMDS,type="n", main="Jaccard - non-metric multidimensional scaling", choices=c(1,3))
mtext(paste("Dimensions 1 and 3 - Adonis R2:",round(ad_j_org[3],4),"p =",round(ad_j_org[4],4),"q =",round(ad_j_org[5],4)))
# # # orditorp(NMDS,display="species",col="gray",air=0.01, cex=0.5, choices=c(1,3))
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=organ, choices=c(1,3))
ordihull(NMDS, groups=organ, draw="polygon", col=colors2, border=colors2, label=F, choices=c(1,3))
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main="Jaccard - non-metric multidimensional scaling", choices=c(2,3))
mtext(paste("Dimensions 2 and 3 - Adonis R2:",round(ad_j_org[3],4),"p =",round(ad_j_org[4],4),"q =",round(ad_j_org[5],4)))
# # # orditorp(NMDS,display="species",col="gray",air=0.01, cex=0.5, choices=c(2,3))
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=organ, choices=c(2,3))
ordihull(NMDS, groups=organ, draw="polygon", col=colors2, border=colors2, label=F, choices=c(2,3))
# ordiplot(NMDS,type="n",xaxt='n',yaxt='n',frame=F,xlab="",ylab="")

# Now the effect of unique items on the clustering dimensions 1 and 2
ordiplot(NMDS,type="n", main="Landscape of unique items - Jaccard", choices=c(1,2))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",round(ad_j_org[3],4),"p =",round(ad_j_org[4],4),"q =",round(ad_j_org[5],4)))
ordisurf(NMDS, unique, col="gray", choices=c(1,2),add=TRUE)
points(NMDS, col=organ, pch=19,choices=c(1,2))
legend("topleft", legend=c("V3","V4","V3V4"), pch=19, col=colors2,cex=1)
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main="Landscape of unique items - Jaccard", choices=c(1,3))
mtext(paste("Dimensions 1 and 3 - Adonis R2:",round(ad_j_org[3],4),"p =",round(ad_j_org[4],4),"q =",round(ad_j_org[5],4)))
ordisurf(NMDS, unique, col="gray", choices=c(1,3),add=TRUE)
points(NMDS, col=organ, pch=19,choices=c(1,3))
legend("topleft", legend=c("V3","V4","V3V4"), pch=19, col=colors2,cex=1)
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main="Landscape of unique items - Jaccard", choices=c(2,3))
mtext(paste("Dimensions 2 and 3 - Adonis R2:",round(ad_j_org[3],4),"p =",round(ad_j_org[4],4),"q =",round(ad_j_org[5],4)))
ordisurf(NMDS, unique, col="gray", choices=c(2,3),add=TRUE)
points(NMDS, col=organ, pch=19,choices=c(2,3))
legend("topleft", legend=c("V3","V4","V3V4"), pch=19, col=colors2,cex=1)

# Now the hierarchical clustering of unique items on the clustering dimensions 1 and 2
ordiplot(NMDS,type="n", main="Hierachical clustering - Jaccard", choices=c(1,2))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",round(ad_j_org[3],4),"p =",round(ad_j_org[4],4),"q =",round(ad_j_org[5],4)))
ordicluster(NMDS, cluster=hclust(jaccard), choices=c(1,2),col="darkgray") #Hierachical clusering
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=organ, choices=c(1,2))
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main="Hierachical clustering - Jaccard", choices=c(1,3))
mtext(paste("Dimensions 1 and 3 - Adonis R2:",round(ad_j_org[3],4),"p =",round(ad_j_org[4],4),"q =",round(ad_j_org[5],4)))
ordicluster(NMDS, cluster=hclust(jaccard), choices=c(1,3),col="darkgray") #Hierachical clusering
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=organ, choices=c(1,3))
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main="Hierachical clustering - Jaccard", choices=c(2,3))
mtext(paste("Dimensions 2 and 3 - Adonis R2:",round(ad_j_org[3],4),"p =",round(ad_j_org[4],4),"q =",round(ad_j_org[5],4)))
ordicluster(NMDS, cluster=hclust(jaccard), choices=c(2,3),col="darkgray") #Hierachical clusering
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=organ, choices=c(2,3))
dev.off()

 ### Color by organ/regions
pdf(paste(prefix,"jaccard_reg-org.pdf", sep="-"))
# First dimensions 1 and 2
ordiplot(NMDS,type="n", main="Jaccard - non-metric multidimensional scaling", choices=c(1,2))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",round(all_j_reg_org[3],4),"p =",round(all_j_reg_org[4],4),"q =",round(all_j_reg_org[5],4)))
# # # orditorp(NMDS,display="species",col="gray",air=0.01, cex=0.5, choices=c(1,2))
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=reg_org, choices=c(1,2))
ordihull(NMDS, groups=reg_org, draw="polygon", col=colors6, border=colors6, label=F, choices=c(1,2))
# Now dimensions 1 and 3
ordiplot(NMDS,type="n", main="Jaccard - non-metric multidimensional scaling", choices=c(1,3))
mtext(paste("Dimensions 1 and 3 - Adonis R2:",round(all_j_reg_org[3],4),"p =",round(all_j_reg_org[4],4),"q =",round(all_j_reg_org[5],4)))
# # # orditorp(NMDS,display="species",col="gray",air=0.01, cex=0.5, choices=c(1,3))
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=reg_org, choices=c(1,3))
ordihull(NMDS, groups=reg_org, draw="polygon", col=colors6, border=colors6, label=F, choices=c(1,3))
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main="Jaccard - non-metric multidimensional scaling", choices=c(2,3))
mtext(paste("Dimensions 2 and 3 - Adonis R2:",round(all_j_reg_org[3],4),"p =",round(all_j_reg_org[4],4),"q =",round(all_j_reg_org[5],4)))
# # # orditorp(NMDS,display="species",col="gray",air=0.01, cex=0.5, choices=c(2,3))
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=reg_org, choices=c(2,3))
ordihull(NMDS, groups=reg_org, draw="polygon", col=colors6, border=colors6, label=F, choices=c(2,3))
# ordiplot(NMDS,type="n",xaxt='n',yaxt='n',frame=F,xlab="",ylab="")

# Now the effect of unique items on the clustering dimensions 1 and 2
ordiplot(NMDS,type="n", main="Landscape of unique items - Jaccard", choices=c(1,2))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",round(all_j_reg_org[3],4),"p =",round(all_j_reg_org[4],4),"q =",round(all_j_reg_org[5],4)))
ordisurf(NMDS, unique, col="gray", choices=c(1,2),add=TRUE)
points(NMDS, col=reg_org, pch=19,choices=c(1,2))
legend("topleft", legend=c("V3-I","V3-H","V4-I","V4-H","V3V4-I","V3V4-H"), pch=19, col=legend_items,cex=0.8)
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main="Landscape of unique items - Jaccard", choices=c(1,3))
mtext(paste("Dimensions 1 and 3 - Adonis R2:",round(all_j_reg_org[3],4),"p =",round(all_j_reg_org[4],4),"q =",round(all_j_reg_org[5],4)))
ordisurf(NMDS, unique, col="gray", choices=c(1,3),add=TRUE)
points(NMDS, col=reg_org, pch=19,choices=c(1,3))
legend("topleft", legend=c("V3-I","V3-H","V4-H","V4-I","V3V4-I","V3V4-H"), pch=19, col=legend_items,cex=0.8)
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main="Landscape of unique items - Jaccard", choices=c(2,3))
mtext(paste("Dimensions 2 and 3 - Adonis R2:",round(all_j_reg_org[3],4),"p =",round(all_j_reg_org[4],4),"q =",round(all_j_reg_org[5],4)))
ordisurf(NMDS, unique, col="gray", choices=c(2,3),add=TRUE)
points(NMDS, col=reg_org, pch=19,choices=c(2,3))
legend("topleft", legend=c("V3-I","V3-H","V4-I","V4-H","V3V4-I","V3V4-H"), pch=19, col=legend_items,cex=0.8)

# Now the hierarchical clustering of unique items on the clustering dimensions 1 and 2
ordiplot(NMDS,type="n", main="Hierachical clustering - Jaccard", choices=c(1,2))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",round(all_j_reg_org[3],4),"p =",round(all_j_reg_org[4],4),"q =",round(all_j_reg_org[5],4)))
ordicluster(NMDS, cluster=hclust(jaccard), choices=c(1,2),col="darkgray") #Hierachical clusering
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=reg_org, choices=c(1,2))
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main="Hierachical clustering - Jaccard", choices=c(1,3))
mtext(paste("Dimensions 1 and 3 - Adonis R2:",round(all_j_reg_org[3],4),"p =",round(all_j_reg_org[4],4),"q =",round(all_j_reg_org[5],4)))
ordicluster(NMDS, cluster=hclust(jaccard), choices=c(1,3),col="darkgray") #Hierachical clusering
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=reg_org, choices=c(1,3))
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main="Hierachical clustering - Jaccard", choices=c(2,3))
mtext(paste("Dimensions 2 and 3 - Adonis R2:",round(all_j_reg_org[3],4),"p =",round(all_j_reg_org[4],4),"q =",round(all_j_reg_org[5],4)))
ordicluster(NMDS, cluster=hclust(jaccard), choices=c(2,3),col="darkgray") #Hierachical clusering
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=reg_org, choices=c(2,3))
dev.off()

 ########### We'll continue with bray curtis (abundance) ###########
bray <- vegdist(t(df),method="bray",upper=TRUE,diag=TRUE)
write.table(as.matrix(bray),paste(prefix,"bray_curtis.tsv", sep="-"), sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
# Get the number of dimensions to use for NMDS (too much or too less stress on two axes is uninformative)
stress <- data.frame(matrix(NA,nrow=15,ncol=1)) # Set a list to get the 
for (i in 1:10){
	NMDS <- metaMDS(bray,distance=bray,k=i,trymax=100,autotransform=FALSE,wascores = FALSE)
	stress[i,1] <- NMDS$stress
}
dimensions <- which(abs(diff(stress[,1]))<0.01)[1] # To set a suitable number of dimensions, use the vector of stress values for each dimension set. The cutoff is set to less than 0.01. Since the difference vector is shifted by -1, we can use the first item not passing (the number actually points to the last item passing).

# NMDS <- metaMDS(bray,distance=bray,k=dimensions,trymax=1000,autotransform=FALSE,wascores = FALSE) # This uses a previously calculated distance
NMDS <- metaMDS(t(df),distance="bray",k=dimensions,trymax=1000,autotransform=FALSE,wascores = TRUE) # This includes species scores and uses the whole table instead

 ### Color by regions
pdf(paste(prefix,"bray_regions.pdf", sep="-"))
# par(mfrow=c(2,2))
# First plot dimensions 1 and 2
ordiplot(NMDS,type="n", main="Bray Curtis - non-metric multidimensional scaling", choices=c(1,2))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",round(all_j_reg[3],4),"p =",round(all_b_reg[4],4),"q =",round(all_b_reg[5],4)))
# # # orditorp(NMDS,display="species",col="gray",air=0.01, cex=0.5, choices=c(1,2))
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=region, choices=c(1,2))
ordihull(NMDS, groups=region, draw="polygon", col=colors3, border=colors3, label=F, choices=c(1,2))
# ordiellipse(NMDS, groups=region, draw="polygon", display="sites", col=colors3, border=colors3, label=F, choices=c(1,2), kind="ehull")
# ordispider(NMDS,groups=region,display="sites",col=colors3,label=F, choices=c(1,2))
# ordibar(NMDS, groups=region, display="sites", col=colors3, label=F, choices=c(1,2))
# ordicluster(NMDS, cluster=hclust(bray), choices=c(1,2)) #Hierachical clusering
# Now dimensions 1 and 3
ordiplot(NMDS,type="n", main="Bray Curtis - non-metric multidimensional scaling", choices=c(1,3))
mtext(paste("Dimensions 1 and 3 - Adonis R2:",round(all_b_reg[3],4),"p =",round(all_b_reg[4],4),"q =",round(all_b_reg[5],4)))
# # # orditorp(NMDS,display="species",col="gray",air=0.01, cex=0.5, choices=c(1,3))
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=region, choices=c(1,3))
ordihull(NMDS, groups=region, draw="polygon", col=colors3, border=colors3, label=F, choices=c(1,3))
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main="Bray Curtis - non-metric multidimensional scaling", choices=c(2,3))
mtext(paste("Dimensions 2 and 3 - Adonis R2:",round(all_b_reg[3],4),"p =",round(all_b_reg[4],4),"q =",round(all_b_reg[5],4)))
# # # orditorp(NMDS,display="species",col="gray",air=0.01, cex=0.5, choices=c(2,3))
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=region, choices=c(2,3))
ordihull(NMDS, groups=region, draw="polygon", col=colors3, border=colors3, label=F, choices=c(2,3))
# ordiplot(NMDS,type="n",xaxt='n',yaxt='n',frame=F,xlab="",ylab="")

# Now the effect of unique items on the clustering dimensions 1 and 2
ordiplot(NMDS,type="n", main="Landscape of unique items - Bray Curtis", choices=c(1,2))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",round(all_b_reg[3],4),"p =",round(all_b_reg[4],4),"q =",round(all_b_reg[5],4)))
ordisurf(NMDS, unique, col="gray", choices=c(1,2),add=TRUE)
points(NMDS, col=region, pch=19,choices=c(1,2))
legend("topleft", legend=c("V3","V4","V3V4"), pch=19, col=colors3,cex=1)
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main="Landscape of unique items - Bray Curtis", choices=c(1,3))
mtext(paste("Dimensions 1 and 3 - Adonis R2:",round(all_b_reg[3],4),"p =",round(all_b_reg[4],4),"q =",round(all_b_reg[5],4)))
ordisurf(NMDS, unique, col="gray", choices=c(1,3),add=TRUE)
points(NMDS, col=region, pch=19,choices=c(1,3))
legend("topleft", legend=c("V3","V4","V3V4"), pch=19, col=colors3,cex=1)
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main="Landscape of unique items - Bray Curtis", choices=c(2,3))
mtext(paste("Dimensions 2 and 3 - Adonis R2:",round(all_b_reg[3],4),"p =",round(all_b_reg[4],4),"q =",round(all_b_reg[5],4)))
ordisurf(NMDS, unique, col="gray", choices=c(2,3),add=TRUE)
points(NMDS, col=region, pch=19,choices=c(2,3))
legend("topleft", legend=c("V3","V4","V3V4"), pch=19, col=colors3,cex=1)

# Now the hierarchical clustering of unique items on the clustering dimensions 1 and 2
ordiplot(NMDS,type="n", main="Hierachical clustering - Bray Curtis", choices=c(1,2))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",round(all_b_reg[3],4),"p =",round(all_b_reg[4],4),"q =",round(all_b_reg[5],4)))
ordicluster(NMDS, cluster=hclust(bray), choices=c(1,2),col="darkgray") #Hierachical clusering
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=region, choices=c(1,2))
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main="Hierachical clustering - Bray Curtis", choices=c(1,3))
mtext(paste("Dimensions 1 and 3 - Adonis R2:",round(all_b_reg[3],4),"p =",round(all_b_reg[4],4),"q =",round(all_b_reg[5],4)))
ordicluster(NMDS, cluster=hclust(bray), choices=c(1,3),col="darkgray") #Hierachical clusering
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=region, choices=c(1,3))
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main="Hierachical clustering - Bray Curtis", choices=c(2,3))
mtext(paste("Dimensions 2 and 3 - Adonis R2:",round(all_b_reg[3],4),"p =",round(all_b_reg[4],4),"q =",round(all_b_reg[5],4)))
ordicluster(NMDS, cluster=hclust(bray), choices=c(2,3),col="darkgray") #Hierachical clusering
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=region, choices=c(2,3))
dev.off()

 ### Color by organs
pdf(paste(prefix,"bray_organs.pdf", sep="-"))
# First dimensions 1 and 2
ordiplot(NMDS,type="n", main="Bray Curtis - non-metric multidimensional scaling", choices=c(1,2))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",round(ad_b_org[3],4),"p =",round(ad_b_org[4],4),"q =",round(ad_b_org[5],4)))
# # # orditorp(NMDS,display="species",col="gray",air=0.01, cex=0.5, choices=c(1,2))
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=organ, choices=c(1,2))
ordihull(NMDS, groups=organ, draw="polygon", col=colors2, border=colors2, label=F, choices=c(1,2))
# Now dimensions 1 and 3
ordiplot(NMDS,type="n", main="Bray Curtis - non-metric multidimensional scaling", choices=c(1,3))
mtext(paste("Dimensions 1 and 3 - Adonis R2:",round(ad_b_org[3],4),"p =",round(ad_b_org[4],4),"q =",round(ad_b_org[5],4)))
# # # orditorp(NMDS,display="species",col="gray",air=0.01, cex=0.5, choices=c(1,3))
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=organ, choices=c(1,3))
ordihull(NMDS, groups=organ, draw="polygon", col=colors2, border=colors2, label=F, choices=c(1,3))
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main="Bray Curtis - non-metric multidimensional scaling", choices=c(2,3))
mtext(paste("Dimensions 2 and 3 - Adonis R2:",round(ad_b_org[3],4),"p =",round(ad_b_org[4],4),"q =",round(ad_b_org[5],4)))
# # # orditorp(NMDS,display="species",col="gray",air=0.01, cex=0.5, choices=c(2,3))
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=organ, choices=c(2,3))
ordihull(NMDS, groups=organ, draw="polygon", col=colors2, border=colors2, label=F, choices=c(2,3))
# ordiplot(NMDS,type="n",xaxt='n',yaxt='n',frame=F,xlab="",ylab="")

# Now the effect of unique items on the clustering dimensions 1 and 2
ordiplot(NMDS,type="n", main="Landscape of unique items - Bray Curtis", choices=c(1,2))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",round(ad_b_org[3],4),"p =",round(ad_b_org[4],4),"q =",round(ad_b_org[5],4)))
ordisurf(NMDS, unique, col="gray", choices=c(1,2),add=TRUE)
points(NMDS, col=organ, pch=19,choices=c(1,2))
legend("topleft", legend=c("V3","V4","V3V4"), pch=19, col=colors2,cex=1)
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main="Landscape of unique items - Bray Curtis", choices=c(1,3))
mtext(paste("Dimensions 1 and 3 - Adonis R2:",round(ad_b_org[3],4),"p =",round(ad_b_org[4],4),"q =",round(ad_b_org[5],4)))
ordisurf(NMDS, unique, col="gray", choices=c(1,3),add=TRUE)
points(NMDS, col=organ, pch=19,choices=c(1,3))
legend("topleft", legend=c("V3","V4","V3V4"), pch=19, col=colors2,cex=1)
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main="Landscape of unique items - Bray Curtis", choices=c(2,3))
mtext(paste("Dimensions 2 and 3 - Adonis R2:",round(ad_b_org[3],4),"p =",round(ad_b_org[4],4),"q =",round(ad_b_org[5],4)))
ordisurf(NMDS, unique, col="gray", choices=c(2,3),add=TRUE)
points(NMDS, col=organ, pch=19,choices=c(2,3))
legend("topleft", legend=c("V3","V4","V3V4"), pch=19, col=colors2,cex=1)

# Now the hierarchical clustering of unique items on the clustering dimensions 1 and 2
ordiplot(NMDS,type="n", main="Hierachical clustering - Bray Curtis", choices=c(1,2))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",round(ad_b_org[3],4),"p =",round(ad_b_org[4],4),"q =",round(ad_b_org[5],4)))
ordicluster(NMDS, cluster=hclust(bray), choices=c(1,2),col="darkgray") #Hierachical clusering
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=organ, choices=c(1,2))
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main="Hierachical clustering - Bray Curtis", choices=c(1,3))
mtext(paste("Dimensions 1 and 3 - Adonis R2:",round(ad_b_org[3],4),"p =",round(ad_b_org[4],4),"q =",round(ad_b_org[5],4)))
ordicluster(NMDS, cluster=hclust(bray), choices=c(1,3),col="darkgray") #Hierachical clusering
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=organ, choices=c(1,3))
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main="Hierachical clustering - Bray Curtis", choices=c(2,3))
mtext(paste("Dimensions 2 and 3 - Adonis R2:",round(ad_b_org[3],4),"p =",round(ad_b_org[4],4),"q =",round(ad_b_org[5],4)))
ordicluster(NMDS, cluster=hclust(bray), choices=c(2,3),col="darkgray") #Hierachical clusering
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=organ, choices=c(2,3))
dev.off()

 ### Color by organ/regions
pdf(paste(prefix,"bray_reg-org.pdf", sep="-"))
# First dimensions 1 and 2
ordiplot(NMDS,type="n", main="Bray Curtis - non-metric multidimensional scaling", choices=c(1,2))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",round(all_b_reg_org[3],4),"p =",round(all_b_reg_org[4],4),"q =",round(all_b_reg_org[5],4)))
# # # orditorp(NMDS,display="species",col="gray",air=0.01, cex=0.5, choices=c(1,2))
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=reg_org, choices=c(1,2))
ordihull(NMDS, groups=reg_org, draw="polygon", col=colors6, border=colors6, label=F, choices=c(1,2))
# Now dimensions 1 and 3
ordiplot(NMDS,type="n", main="Bray Curtis - non-metric multidimensional scaling", choices=c(1,3))
mtext(paste("Dimensions 1 and 3 - Adonis R2:",round(all_b_reg_org[3],4),"p =",round(all_b_reg_org[4],4),"q =",round(all_b_reg_org[5],4)))
# # # orditorp(NMDS,display="species",col="gray",air=0.01, cex=0.5, choices=c(1,3))
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=reg_org, choices=c(1,3))
ordihull(NMDS, groups=reg_org, draw="polygon", col=colors6, border=colors6, label=F, choices=c(1,3))
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main="Bray Curtis - non-metric multidimensional scaling", choices=c(2,3))
mtext(paste("Dimensions 2 and 3 - Adonis R2:",round(all_b_reg_org[3],4),"p =",round(all_b_reg_org[4],4),"q =",round(all_b_reg_org[5],4)))
# # # orditorp(NMDS,display="species",col="gray",air=0.01, cex=0.5, choices=c(2,3))
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=reg_org, choices=c(2,3))
ordihull(NMDS, groups=reg_org, draw="polygon", col=colors6, border=colors6, label=F, choices=c(2,3))
# ordiplot(NMDS,type="n",xaxt='n',yaxt='n',frame=F,xlab="",ylab="")

# Now the effect of unique items on the clustering dimensions 1 and 2
ordiplot(NMDS,type="n", main="Landscape of unique items - Bray Curtis", choices=c(1,2))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",round(all_b_reg_org[3],4),"p =",round(all_b_reg_org[4],4),"q =",round(all_b_reg_org[5],4)))
ordisurf(NMDS, unique, col="gray", choices=c(1,2),add=TRUE)
points(NMDS, col=reg_org, pch=19,choices=c(1,2))
legend("topleft", legend=c("V3-I","V3-H","V4-I","V4-H","V3V4-I","V3V4-H"), pch=19, col=legend_items,cex=0.8)
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main="Landscape of unique items - Bray Curtis", choices=c(1,3))
mtext(paste("Dimensions 1 and 3 - Adonis R2:",round(all_b_reg_org[3],4),"p =",round(all_b_reg_org[4],4),"q =",round(all_b_reg_org[5],4)))
ordisurf(NMDS, unique, col="gray", choices=c(1,3),add=TRUE)
points(NMDS, col=reg_org, pch=19,choices=c(1,3))
legend("topleft", legend=c("V3-I","V3-H","V4-H","V4-I","V3V4-I","V3V4-H"), pch=19, col=legend_items,cex=0.8)
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main="Landscape of unique items - Bray Curtis", choices=c(2,3))
mtext(paste("Dimensions 2 and 3 - Adonis R2:",round(all_b_reg_org[3],4),"p =",round(all_b_reg_org[4],4),"q =",round(all_b_reg_org[5],4)))
ordisurf(NMDS, unique, col="gray", choices=c(2,3),add=TRUE)
points(NMDS, col=reg_org, pch=19,choices=c(2,3))
legend("topleft", legend=c("V3-I","V3-H","V4-I","V4-H","V3V4-I","V3V4-H"), pch=19, col=legend_items,cex=0.8)

# Now the hierarchical clustering of unique items on the clustering dimensions 1 and 2
ordiplot(NMDS,type="n", main="Hierachical clustering - Bray Curtis", choices=c(1,2))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",round(all_b_reg_org[3],4),"p =",round(all_b_reg_org[4],4),"q =",round(all_b_reg_org[5],4)))
ordicluster(NMDS, cluster=hclust(bray), choices=c(1,2),col="darkgray") #Hierachical clusering
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=reg_org, choices=c(1,2))
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main="Hierachical clustering - Bray Curtis", choices=c(1,3))
mtext(paste("Dimensions 1 and 3 - Adonis R2:",round(all_b_reg_org[3],4),"p =",round(all_b_reg_org[4],4),"q =",round(all_b_reg_org[5],4)))
ordicluster(NMDS, cluster=hclust(bray), choices=c(1,3),col="darkgray") #Hierachical clusering
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=reg_org, choices=c(1,3))
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main="Hierachical clustering - Bray Curtis", choices=c(2,3))
mtext(paste("Dimensions 2 and 3 - Adonis R2:",round(all_b_reg_org[3],4),"p =",round(all_b_reg_org[4],4),"q =",round(all_b_reg_org[5],4)))
ordicluster(NMDS, cluster=hclust(bray), choices=c(2,3),col="darkgray") #Hierachical clusering
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=reg_org, choices=c(2,3))
dev.off()
