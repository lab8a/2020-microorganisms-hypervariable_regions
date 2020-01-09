# Started on 2019-09-13
# by Rodrigo García-López for Adrian Ochoa-Leyva's Metagenomics and Metatranscriptomics Lab at IBt, UNAM, Cuernavaca, Mexico.
# Under GNU GPLv3 license
# Disclamer: This script was written for use with specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.
# This script was tested with R 3.6.0 - "Planting of a Tree"
# The script is intended to plot a NMDS and PCoA using a distance/similarity matrix and paint by groups
# DEPRECATED: PCA was originally intended to be included but it requires an euclidian matrix (e.g. TSS+CLR transform) so it shouldn't be mixed with methods for distance matrix as input (I have created a separate script for this).
# The input is a square matrix of samples x samples using a metric and has no previous comments in the header

# Run as follows:
# cat table.tsv|~/bin/R-3.6.0/bin/Rscript multi_ordination_methods.R

# Tested with command:
# cat /home/rodrigo/Shrimp_2015-2018/V3-V4_comparison/12_diversity_analyses/02_using_sparsity_reduced_tables/08_distance_matrices/02_split_region/06_Final_mix_gg/06_Final_mix_gg97-tableOTU-V3_unweighted_unifrac.tsv|~/bin/R-3.6.0/bin/Rscript multi_ordination_methods.R test
# In R:
# df <- read.table("12_diversity_analyses/02_using_sparsity_reduced_tables/08_distance_matrices/02_split_region/06_Final_mix_gg/06_Final_mix_gg97-tableOTU-V3_unweighted_unifrac.tsv", sep="\t",header=T, skip=0, comment.char='',quote="",fill=F, row.names=1)

args <- commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<1) { # at least, one arguments should be included: <out_name_prefix>
  stop("A minimum of 1 argument is mandatory: cat table.tsv|Rscript multi_ordination_methods.R <output_file>", call.=FALSE)
}
prefix <- as.character(args[1]) # Get a string handle to create output names
df <- read.table(file('stdin'), sep="\t",header=T, skip=0, comment.char='',quote="",fill=F, row.names=1) #Read the table, row names expected in first column, matching colnames in first row

 ### Define groups and set parameters ##
samples <- colnames(df) # also store sample names
# unique <- apply(df>0,2,sum) # Get a vector of unique features per sample
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
library("vegan") # Distance matrix calculations are included the vegan package

#### This should be used for inputting the rarefied matrix ###
# dist <- vegdist(t(df),method=t(df),upper=TRUE,diag=TRUE)
# bray <- vegdist(t(df),method="bray",upper=TRUE,diag=TRUE)
# feature_names <- rownames(df) # Store the feature names
# rownames(df) <- paste("f",1:nrow(df),sep="") #Replace the name with a features number
# write.table(cbind(rownames(df),feature_names),paste(prefix,"feature_list.tsv", sep="-"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

 ########### We'll start with evaluating the matrix support for groups ###########
# anosim_jaccard_org <- anosim(jaccard,distance=jaccard,organ,permutations=10000);anosim_jaccard_org
# adonis_jaccard_org <- adonis(as.matrix(jaccard) ~ as.factor(organ),method=jaccard,permutations=10000);adonis_jaccard_org
# anosim_jaccard_reg <- anosim(jaccard,distance=jaccard,region,permutations=10000);anosim_jaccard_reg
if (length(levels(as.factor(region))) > 1) { 
	adonis_reg <- adonis(as.matrix(t(df)) ~ as.factor(region),method=t(df),permutations=10000)
	reg_R2 <- round(adonis_reg$aov.tab[1,5],4)
	reg_pval <- round(adonis_reg$aov.tab[1,6],4)
	} else { reg_R2=reg_pval="NA" }
if (length(levels(as.factor(organ))) > 1) { 
	adonis_org <- adonis(as.matrix(t(df)) ~ as.factor(organ),method=t(df),permutations=10000)
	org_R2 <- round(adonis_org$aov.tab[1,5],4)
	org_pval <- round(adonis_org$aov.tab[1,6],4)
	} else { org_R2=org_pval="NA" }
if (length(levels(as.factor(reg_org))) > 1) { 
	adonis_reg_org <- adonis(as.matrix(t(df)) ~ as.factor(reg_org),method=t(df),permutations=10000)
	reg_org_R2 <- round(adonis_reg_org$aov.tab[1,5],4)
	reg_org_pval <- round(adonis_reg_org$aov.tab[1,6],4)
	} else { reg_org_R2=reg_org_pval="NA" }


degf <- length(df)-1 # Get the total degrees of freedom for total communities (n-1)
degf <- ifelse(degf>10,10,degf) #Carry a maximum of 10 tests
# Get the number of dimensions to use for NMDS (too much or too less stress on two axes is uninformative)
converge <- stress <- data.frame(matrix(NA,nrow=degf,ncol=1)) #Set a list to get the stress and another one to see if it converges
for (i in 1:degf){
	NMDS <- metaMDS(t(df),distance=t(df),k=i,trymax=100,autotransform=FALSE,wascores = FALSE)
	stress[i,1] <- NMDS$stress
	converge[i,1] <- NMDS$converged
}
converge[is.na(converge)] <- FALSE # don't consider those that have NAs
stress <- stress[converge[,1],]

dimensions <- which(abs(diff(stress))<0.01)[1]+1 # To set a suitable number of dimensions, use the vector of stress values for each dimension set. The cutoff is set to less than 0.01. Since the difference vector is shifted by 1, we can use the first item not passing plus 1 (the number actually points to the last item passing).
if(is.na(dimensions)){dimensions=length(stress)} # If we have large changes, we may not reach a convergence yet, so use the max number of valid items
dimensions <- ifelse(dimensions<3,3,dimensions) # Force at least 3 dimensions

NMDS <- metaMDS(t(df),distance=t(df),k=dimensions,trymax=1000,autotransform=FALSE,wascores = TRUE) # This includes species scores and uses the whole table instead
tries <- NMDS$tries
converge <- NMDS$converged
# stressplot(NMDS)
# plot(NMDS, type="n")
# points(NMDS, display=c("sites"), choices=c(1,2), pch=3, col="red")


 ### Color by regions
pdf(paste(prefix,"nmds_regions.pdf", sep="-"))
# First plot dimensions 1 and 2
ordiplot(NMDS,type="n", main=paste(prefix,"NMDS-region"), choices=c(1,2))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",reg_R2,"p =",reg_pval))
# # # orditorp(NMDS,display="species",col="gray",air=0.01, cex=0.5, choices=c(1,2))
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=region, choices=c(1,2))
ordihull(NMDS, groups=region, draw="polygon", col=colors3, border=colors3, label=F, choices=c(1,2))
# Now dimensions 1 and 3
ordiplot(NMDS,type="n", main=paste(prefix,"NMDS-region"), choices=c(1,3))
mtext(paste("Dimensions 1 and 3 - Adonis R2:",reg_R2,"p =",reg_pval))
# # # orditorp(NMDS,display="species",col="gray",air=0.01, cex=0.5, choices=c(1,3))
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=region, choices=c(1,3))
ordihull(NMDS, groups=region, draw="polygon", col=colors3, border=colors3, label=F, choices=c(1,3))
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main=paste(prefix,"NMDS-region"), choices=c(2,3))
mtext(paste("Dimensions 2 and 3 - Adonis R2:",reg_R2,"p =",reg_pval))
# # # orditorp(NMDS,display="species",col="gray",air=0.01, cex=0.5, choices=c(2,3))
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=region, choices=c(2,3))
ordihull(NMDS, groups=region, draw="polygon", col=colors3, border=colors3, label=F, choices=c(2,3))
# ordiplot(NMDS,type="n",xaxt='n',yaxt='n',frame=F,xlab="",ylab="")

# Now the effect of unique items on the clustering dimensions 1 and 2
ordiplot(NMDS,type="n", main=paste(prefix,"NMDS-region"), choices=c(1,2))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",reg_R2,"p =",reg_pval))
points(NMDS, col=region, pch=19,choices=c(1,2))
legend("topleft", legend=levels(as.factor(region1)), pch=19, col=colors3,cex=1)
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main=paste(prefix,"NMDS-region"), choices=c(1,3))
mtext(paste("Dimensions 1 and 3 - Adonis R2:",reg_R2,"p =",reg_pval))
points(NMDS, col=region, pch=19,choices=c(1,3))
legend("topleft", legend=levels(as.factor(region1)), pch=19, col=colors3,cex=1)
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main=paste(prefix,"NMDS-region"), choices=c(2,3))
mtext(paste("Dimensions 2 and 3 - Adonis R2:",reg_R2,"p =",reg_pval))
points(NMDS, col=region, pch=19,choices=c(2,3))
legend("topleft", legend=levels(as.factor(region1)), pch=19, col=colors3,cex=1)

# Now the hierarchical clustering of unique items on the clustering dimensions 1 and 2
ordiplot(NMDS,type="n", main=paste(prefix,"Hclust-region"), choices=c(1,2))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",reg_R2,"p =",reg_pval))
ordicluster(NMDS, cluster=hclust(as.dist(df)), choices=c(1,2),col="darkgray") #Hierachical clusering
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=region, choices=c(1,2))
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main=paste(prefix,"Hclust-region"), choices=c(1,3))
mtext(paste("Dimensions 1 and 3 - Adonis R2:",reg_R2,"p =",reg_pval))
ordicluster(NMDS, cluster=hclust(as.dist(df)), choices=c(1,3),col="darkgray") #Hierachical clusering
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=region, choices=c(1,3))
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main=paste(prefix,"Hclust-region"), choices=c(2,3))
mtext(paste("Dimensions 2 and 3 - Adonis R2:",reg_R2,"p =",reg_pval))
ordicluster(NMDS, cluster=hclust(as.dist(df)), choices=c(2,3),col="darkgray") #Hierachical clusering
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=region, choices=c(2,3))
dev.off()

 ### Color by organs
pdf(paste(prefix,"nmds_organs.pdf", sep="-"))
# First plot dimensions 1 and 2
ordiplot(NMDS,type="n", main=paste(prefix,"NMDS-organ"), choices=c(1,2))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",org_R2,"p =",org_pval))
# # # orditorp(NMDS,display="species",col="gray",air=0.01, cex=0.5, choices=c(1,2))
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=organ, choices=c(1,2))
ordihull(NMDS, groups=organ, draw="polygon", col=colors2, border=colors2, label=F, choices=c(1,2))
# Now dimensions 1 and 3
ordiplot(NMDS,type="n", main=paste(prefix,"NMDS-organ"), choices=c(1,3))
mtext(paste("Dimensions 1 and 3 - Adonis R2:",org_R2,"p =",org_pval))
# # # orditorp(NMDS,display="species",col="gray",air=0.01, cex=0.5, choices=c(1,3))
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=organ, choices=c(1,3))
ordihull(NMDS, groups=organ, draw="polygon", col=colors2, border=colors2, label=F, choices=c(1,3))
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main=paste(prefix,"NMDS-organ"), choices=c(2,3))
mtext(paste("Dimensions 2 and 3 - Adonis R2:",org_R2,"p =",org_pval))
# # # orditorp(NMDS,display="species",col="gray",air=0.01, cex=0.5, choices=c(2,3))
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=organ, choices=c(2,3))
ordihull(NMDS, groups=organ, draw="polygon", col=colors2, border=colors2, label=F, choices=c(2,3))
# ordiplot(NMDS,type="n",xaxt='n',yaxt='n',frame=F,xlab="",ylab="")

# Now the items with no polygons dimensions 1 and 2
ordiplot(NMDS,type="n", main=paste(prefix,"NMDS-organ"), choices=c(1,2))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",org_R2,"p =",org_pval))
points(NMDS, col=organ, pch=19,choices=c(1,2))
legend("topleft", legend=levels(as.factor(organ1)), pch=19, col=colors2,cex=1)
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main=paste(prefix,"NMDS-organ"), choices=c(1,3))
mtext(paste("Dimensions 1 and 3 - Adonis R2:",org_R2,"p =",org_pval))
points(NMDS, col=organ, pch=19,choices=c(1,3))
legend("topleft", legend=levels(as.factor(organ1)), pch=19, col=colors2,cex=1)
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main=paste(prefix,"NMDS-organ"), choices=c(2,3))
mtext(paste("Dimensions 2 and 3 - Adonis R2:",org_R2,"p =",org_pval))
points(NMDS, col=organ, pch=19,choices=c(2,3))
legend("topleft", legend=levels(as.factor(organ1)), pch=19, col=colors2,cex=1)

# Now the hierarchical clustering of unique items on the clustering dimensions 1 and 2
ordiplot(NMDS,type="n", main=paste(prefix,"Hclust-organ"), choices=c(1,2))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",org_R2,"p =",org_pval))
ordicluster(NMDS, cluster=hclust(as.dist(df)), choices=c(1,2),col="darkgray") #Hierachical clusering
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=organ, choices=c(1,2))
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main=paste(prefix,"Hclust-organ"), choices=c(1,3))
mtext(paste("Dimensions 1 and 3 - Adonis R2:",org_R2,"p =",org_pval))
ordicluster(NMDS, cluster=hclust(as.dist(df)), choices=c(1,3),col="darkgray") #Hierachical clusering
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=organ, choices=c(1,3))
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main=paste(prefix,"Hclust-organ"), choices=c(2,3))
mtext(paste("Dimensions 2 and 3 - Adonis R2:",org_R2,"p =",org_pval))
ordicluster(NMDS, cluster=hclust(as.dist(df)), choices=c(2,3),col="darkgray") #Hierachical clusering
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=organ, choices=c(2,3))
dev.off()


 ### Color by organ/regions
pdf(paste(prefix,"nmds_reg_org.pdf", sep="-"))
# First plot dimensions 1 and 2
ordiplot(NMDS,type="n", main=paste(prefix,"NMDS-reg-org"), choices=c(1,2))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",reg_org_R2,"p =",reg_org_pval))
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=reg_org, choices=c(1,2))
ordihull(NMDS, groups=reg_org, draw="polygon", col=colors6, border=colors6, label=F, choices=c(1,2))
# Now dimensions 1 and 3
ordiplot(NMDS,type="n", main=paste(prefix,"NMDS-reg-org"), choices=c(1,3))
mtext(paste("Dimensions 1 and 3 - Adonis R2:",reg_org_R2,"p =",reg_org_pval))
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=reg_org, choices=c(1,3))
ordihull(NMDS, groups=reg_org, draw="polygon", col=colors6, border=colors6, label=F, choices=c(1,3))
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main=paste(prefix,"NMDS-reg-org"), choices=c(2,3))
mtext(paste("Dimensions 2 and 3 - Adonis R2:",reg_org_R2,"p =",reg_org_pval))
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=reg_org, choices=c(2,3))
ordihull(NMDS, groups=reg_org, draw="polygon", col=colors6, border=colors6, label=F, choices=c(2,3))

# Now the groups
ordiplot(NMDS,type="n", main=paste(prefix,"NMDS-reg-org"), choices=c(1,2))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",reg_org_R2,"p =",reg_org_pval))
points(NMDS, col=reg_org, pch=19,choices=c(1,2))
legend("topleft", legend=levels(as.factor(reg_org1)), pch=19, col=colors6,cex=1)
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main=paste(prefix,"NMDS-reg-org"), choices=c(1,3))
mtext(paste("Dimensions 1 and 3 - Adonis R2:",reg_org_R2,"p =",reg_org_pval))
points(NMDS, col=reg_org, pch=19,choices=c(1,3))
legend("topleft", legend=levels(as.factor(reg_org1)), pch=19, col=colors6,cex=1)
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main=paste(prefix,"NMDS-reg-org"), choices=c(2,3))
mtext(paste("Dimensions 2 and 4 - Adonis R2:",reg_org_R2,"p =",reg_org_pval))
points(NMDS, col=reg_org, pch=19,choices=c(2,3))
legend("topleft", legend=levels(as.factor(reg_org1)), pch=19, col=colors6,cex=1)

# Now the hierarchical clustering of unique items on the clustering dimensions 1 and 2
ordiplot(NMDS,type="n", main=paste(prefix,"NMDS-Hclust-reg_org"), choices=c(1,2))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",reg_org_R2,"p =",reg_org_pval))
ordicluster(NMDS, cluster=hclust(as.dist(df)), choices=c(1,2),col="darkgray") #Hierachical clusering
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=reg_org, choices=c(1,2))
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main=paste(prefix,"NMDS-Hclust-reg_org"), choices=c(1,3))
mtext(paste("Dimensions 1 and 3 - Adonis R2:",reg_org_R2,"p =",reg_org_pval))
ordicluster(NMDS, cluster=hclust(as.dist(df)), choices=c(1,3),col="darkgray") #Hierachical clusering
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=reg_org, choices=c(1,3))
# Now dimensions 2 and 3
ordiplot(NMDS,type="n", main=paste(prefix,"NMDS-Hclust-reg_org"), choices=c(2,3))
mtext(paste("Dimensions 2 and 3 - Adonis R2:",reg_org_R2,"p =",reg_org_pval))
ordicluster(NMDS, cluster=hclust(as.dist(df)), choices=c(2,3),col="darkgray") #Hierachical clusering
orditorp(NMDS,display="sites",cex=0.8,air=0.01,col=reg_org, choices=c(2,3))
dev.off()


 ##################### Eigenvalue based MDS - Principal Coordinate Analysis #############################
PCOA <- cmdscale(t(df), eig = T, k=length(df)-1)
PCOA.axes <- round(PCOA$eig*100/sum(PCOA$eig),1)

 ### Color by regions
pdf(paste(prefix,"PCoA_regions.pdf", sep="-"))
# First plot dimensions 1 and 2
ordiplot(PCOA,type="n", main=paste(prefix,"PCoA-reg"), choices=c(1,2), xlab=paste("Dim1 - Variation explained:",PCOA.axes[1],"%"), ylab=paste("Dim2 - Variation explained:",PCOA.axes[2],"%"))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",reg_R2,"p =",reg_pval))
orditorp(PCOA,display="sites",cex=0.8,air=0.01,col=region, choices=c(1,2))
ordihull(PCOA, groups=region, draw="polygon", col=colors3, border=colors3, label=F, choices=c(1,2))
# Now dimensions 1 and 3
ordiplot(PCOA,type="n", main=paste(prefix,"PCoA-reg"), choices=c(1,3), xlab=paste("Dim1 - Variation explained:",PCOA.axes[1],"%"), ylab=paste("Dim3 - Variation explained:",PCOA.axes[3],"%"))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",reg_R2,"p =",reg_pval))
orditorp(PCOA,display="sites",cex=0.8,air=0.01,col=region, choices=c(1,3))
ordihull(PCOA, groups=region, draw="polygon", col=colors3, border=colors3, label=F, choices=c(1,3))
# Now dimensions 2 and 3
ordiplot(PCOA,type="n", main=paste(prefix,"PCoA-reg"), choices=c(2,3), xlab=paste("Dim1 - Variation explained:",PCOA.axes[2],"%"), ylab=paste("Dim3 - Variation explained:",PCOA.axes[3],"%"))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",reg_R2,"p =",reg_pval))
orditorp(PCOA,display="sites",cex=0.8,air=0.01,col=region, choices=c(2,3))
ordihull(PCOA, groups=region, draw="polygon", col=colors3, border=colors3, label=F, choices=c(2,3))

# Now the items with no polygons dimensions 1 and 2
ordiplot(PCOA,type="n", main=paste(prefix,"PCoA-reg"), choices=c(1,2), xlab=paste("Dim1 - Variation explained:",PCOA.axes[1],"%"), ylab=paste("Dim2 - Variation explained:",PCOA.axes[2],"%"))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",reg_R2,"p =",reg_pval))
points(PCOA$points[,1], PCOA$points[,2], col=region, pch=19)
legend("topleft", legend=levels(as.factor(region1)), pch=19, col=colors3,cex=1)
# Now dimensions 2 and 3
ordiplot(PCOA,type="n", main=paste(prefix,"PCoA-reg"), choices=c(1,3), xlab=paste("Dim1 - Variation explained:",PCOA.axes[1],"%"), ylab=paste("Dim3 - Variation explained:",PCOA.axes[3],"%"))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",reg_R2,"p =",reg_pval))
points(PCOA$points[,1], PCOA$points[,3], col=region, pch=19)
legend("topleft", legend=levels(as.factor(region1)), pch=19, col=colors3,cex=1)
# Now dimensions 2 and 3
ordiplot(PCOA,type="n", main=paste(prefix,"PCoA-reg"), choices=c(2,3), xlab=paste("Dim1 - Variation explained:",PCOA.axes[2],"%"), ylab=paste("Dim3 - Variation explained:",PCOA.axes[3],"%"))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",reg_R2,"p =",reg_pval))
points(PCOA$points[,2], PCOA$points[,3], col=region, pch=19)
legend("topleft", legend=levels(as.factor(region1)), pch=19, col=colors3,cex=1)

# Now the hierarchical clustering of unique items on the clustering dimensions 1 and 2
ordiplot(PCOA,type="n", main=paste(prefix,"PCoA-Hclust-org"), choices=c(1,2), xlab=paste("Dim1 - Variation explained:",PCOA.axes[1],"%"), ylab=paste("Dim2 - Variation explained:",PCOA.axes[2],"%"))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",reg_R2,"p =",reg_pval))
ordicluster(PCOA, cluster=hclust(as.dist(df)), choices=c(1,2),col="darkgray") #Hierachical clusering
orditorp(PCOA,display="sites",cex=0.8,air=0.01,col=region, choices=c(1,2))
# Now dimensions 2 and 3
ordiplot(PCOA,type="n", main=paste(prefix,"PCoA-Hclust-org"), choices=c(1,3), xlab=paste("Dim1 - Variation explained:",PCOA.axes[1],"%"), ylab=paste("Dim3 - Variation explained:",PCOA.axes[3],"%"))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",reg_R2,"p =",reg_pval))
ordicluster(PCOA, cluster=hclust(as.dist(df)), choices=c(1,3),col="darkgray") #Hierachical clusering
orditorp(PCOA,display="sites",cex=0.8,air=0.01,col=region, choices=c(1,3))
# Now dimensions 2 and 3
ordiplot(PCOA,type="n", main=paste(prefix,"PCoA-Hclust-org"), choices=c(2,3), xlab=paste("Dim1 - Variation explained:",PCOA.axes[2],"%"), ylab=paste("Dim3 - Variation explained:",PCOA.axes[3],"%"))
mtext(paste("Dimensions 2 and 3 - Adonis R2:",reg_R2,"p =",reg_pval))
ordicluster(PCOA, cluster=hclust(as.dist(df)), choices=c(2,3),col="darkgray") #Hierachical clusering
orditorp(PCOA,display="sites",cex=0.8,air=0.01,col=region, choices=c(2,3))
dev.off()


 ### Color by organs
pdf(paste(prefix,"PCoA_organs.pdf", sep="-"))
# First plot dimensions 1 and 2
ordiplot(PCOA,type="n", main=paste(prefix,"PCoA-org"), choices=c(1,2), xlab=paste("Dim1 - Variation explained:",PCOA.axes[1],"%"), ylab=paste("Dim2 - Variation explained:",PCOA.axes[2],"%"))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",org_R2,"p =",org_pval))
orditorp(PCOA,display="sites",cex=0.8,air=0.01,col=organ, choices=c(1,2))
ordihull(PCOA, groups=organ, draw="polygon", col=colors2, border=colors2, label=F, choices=c(1,2))
# Now dimensions 1 and 3
ordiplot(PCOA,type="n", main=paste(prefix,"PCoA-org"), choices=c(1,3), xlab=paste("Dim1 - Variation explained:",PCOA.axes[1],"%"), ylab=paste("Dim3 - Variation explained:",PCOA.axes[3],"%"))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",org_R2,"p =",org_pval))
orditorp(PCOA,display="sites",cex=0.8,air=0.01,col=organ, choices=c(1,3))
ordihull(PCOA, groups=organ, draw="polygon", col=colors2, border=colors2, label=F, choices=c(1,3))
# Now dimensions 2 and 3
ordiplot(PCOA,type="n", main=paste(prefix,"PCoA-org"), choices=c(2,3), xlab=paste("Dim2 - Variation explained:",PCOA.axes[2],"%"), ylab=paste("Dim3 - Variation explained:",PCOA.axes[3],"%"))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",org_R2,"p =",org_pval))
orditorp(PCOA,display="sites",cex=0.8,air=0.01,col=organ, choices=c(2,3))
ordihull(PCOA, groups=organ, draw="polygon", col=colors2, border=colors2, label=F, choices=c(2,3))

# Now the items with no polygons dimensions 1 and 2
ordiplot(PCOA,type="n", main=paste(prefix,"PCoA-org"), choices=c(1,2), xlab=paste("Dim1 - Variation explained:",PCOA.axes[1],"%"), ylab=paste("Dim2 - Variation explained:",PCOA.axes[2],"%"))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",org_R2,"p =",org_pval))
points(PCOA$points[,1], PCOA$points[,2], col=organ, pch=19)
legend("topleft", legend=levels(as.factor(organ1)), pch=19, col=colors2,cex=1)
# Now dimensions 2 and 3
ordiplot(PCOA,type="n", main=paste(prefix,"PCoA-org"), choices=c(1,3), xlab=paste("Dim1 - Variation explained:",PCOA.axes[1],"%"), ylab=paste("Dim3 - Variation explained:",PCOA.axes[3],"%"))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",org_R2,"p =",org_pval))
points(PCOA$points[,1], PCOA$points[,3], col=organ, pch=19)
legend("topleft", legend=levels(as.factor(organ1)), pch=19, col=colors2,cex=1)
# Now dimensions 2 and 3
ordiplot(PCOA,type="n", main=paste(prefix,"PCoA-org"), choices=c(2,3), xlab=paste("Dim2 - Variation explained:",PCOA.axes[2],"%"), ylab=paste("Dim3 - Variation explained:",PCOA.axes[3],"%"))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",org_R2,"p =",org_pval))
points(PCOA$points[,2], PCOA$points[,3], col=organ, pch=19)
legend("topleft", legend=levels(as.factor(organ1)), pch=19, col=colors2,cex=1)

# Now the hierarchical clustering of unique items on the clustering dimensions 1 and 2
ordiplot(PCOA,type="n", main=paste(prefix,"PCoA-Hclust-org"), choices=c(1,2), xlab=paste("Dim1 - Variation explained:",PCOA.axes[1],"%"), ylab=paste("Dim2 - Variation explained:",PCOA.axes[2],"%"))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",org_R2,"p =",org_pval))
ordicluster(PCOA, cluster=hclust(as.dist(df)), choices=c(1,2),col="darkgray") #Hierachical clusering
orditorp(PCOA,display="sites",cex=0.8,air=0.01,col=organ, choices=c(1,2))
# Now dimensions 2 and 3
ordiplot(PCOA,type="n", main=paste(prefix,"PCoA-Hclust-org"), choices=c(1,3), xlab=paste("Dim1 - Variation explained:",PCOA.axes[1],"%"), ylab=paste("Dim3 - Variation explained:",PCOA.axes[3],"%"))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",org_R2,"p =",org_pval))
ordicluster(PCOA, cluster=hclust(as.dist(df)), choices=c(1,3),col="darkgray") #Hierachical clusering
orditorp(PCOA,display="sites",cex=0.8,air=0.01,col=organ, choices=c(1,3))
# Now dimensions 2 and 3
ordiplot(PCOA,type="n", main=paste(prefix,"PCoA-Hclust-org"), choices=c(2,3), xlab=paste("Dim2 - Variation explained:",PCOA.axes[2],"%"), ylab=paste("Dim3 - Variation explained:",PCOA.axes[3],"%"))
mtext(paste("Dimensions 2 and 3 - Adonis R2:",org_R2,"p =",org_pval))
ordicluster(PCOA, cluster=hclust(as.dist(df)), choices=c(2,3),col="darkgray") #Hierachical clusering
orditorp(PCOA,display="sites",cex=0.8,air=0.01,col=organ, choices=c(2,3))
dev.off()

 ### Color by regions/organs
pdf(paste(prefix,"PCoA_reg_org.pdf", sep="-"))
# First plot dimensions 1 and 2
ordiplot(PCOA,type="n", main=paste(prefix,"PCoA-reg/org"), choices=c(1,2), xlab=paste("Dim1 - Variation explained:",PCOA.axes[1],"%"), ylab=paste("Dim2 - Variation explained:",PCOA.axes[2],"%"))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",reg_org_R2,"p =",reg_org_pval))
orditorp(PCOA,display="sites",cex=0.8,air=0.01,col=reg_org, choices=c(1,2))
ordihull(PCOA, groups=reg_org, draw="polygon", col=colors6, border=colors6, label=F, choices=c(1,2))
# Now dimensions 1 and 3
ordiplot(PCOA,type="n", main=paste(prefix,"PCoA-reg/org"), choices=c(1,3), xlab=paste("Dim1 - Variation explained:",PCOA.axes[1],"%"), ylab=paste("Dim3 - Variation explained:",PCOA.axes[3],"%"))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",reg_org_R2,"p =",reg_org_pval))
orditorp(PCOA,display="sites",cex=0.8,air=0.01,col=reg_org, choices=c(1,3))
ordihull(PCOA, groups=reg_org, draw="polygon", col=colors6, border=colors6, label=F, choices=c(1,3))
# Now dimensions 2 and 3
ordiplot(PCOA,type="n", main=paste(prefix,"PCoA-reg/org"), choices=c(2,3), xlab=paste("Dim2 - Variation explained:",PCOA.axes[2],"%"), ylab=paste("Dim3 - Variation explained:",PCOA.axes[3],"%"))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",reg_org_R2,"p =",reg_org_pval))
orditorp(PCOA,display="sites",cex=0.8,air=0.01,col=reg_org, choices=c(2,3))
ordihull(PCOA, groups=reg_org, draw="polygon", col=colors6, border=colors6, label=F, choices=c(2,3))

# Now the items with no polygons dimensions 1 and 2
ordiplot(PCOA,type="n", main=paste(prefix,"PCoA-reg/org"), choices=c(1,2), xlab=paste("Dim1 - Variation explained:",PCOA.axes[1],"%"), ylab=paste("Dim2 - Variation explained:",PCOA.axes[2],"%"))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",reg_org_R2,"p =",reg_org_pval))
points(PCOA$points[,1], PCOA$points[,2], col=reg_org, pch=19)
legend("topleft", legend=levels(as.factor(reg_org1)), pch=19, col=colors6,cex=1)
# Now dimensions 2 and 3
ordiplot(PCOA,type="n", main=paste(prefix,"PCoA-reg/org"), choices=c(1,3), xlab=paste("Dim1 - Variation explained:",PCOA.axes[1],"%"), ylab=paste("Dim3 - Variation explained:",PCOA.axes[3],"%"))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",reg_org_R2,"p =",reg_org_pval))
points(PCOA$points[,1], PCOA$points[,3], col=reg_org, pch=19)
legend("topleft", legend=levels(as.factor(reg_org1)), pch=19, col=colors6,cex=1)
# Now dimensions 2 and 3
ordiplot(PCOA,type="n", main=paste(prefix,"PCoA-reg/org"), choices=c(2,3), xlab=paste("Dim2 - Variation explained:",PCOA.axes[2],"%"), ylab=paste("Dim3 - Variation explained:",PCOA.axes[3],"%"))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",reg_org_R2,"p =",reg_org_pval))
points(PCOA$points[,2], PCOA$points[,3], col=reg_org, pch=19)
legend("topleft", legend=levels(as.factor(reg_org1)), pch=19, col=colors6,cex=1)

# Now the hierarchical clustering of unique items on the clustering dimensions 1 and 2
ordiplot(PCOA,type="n", main=paste(prefix,"PCoA-Hclust-org"), choices=c(1,2), xlab=paste("Dim1 - Variation explained:",PCOA.axes[1],"%"), ylab=paste("Dim2 - Variation explained:",PCOA.axes[2],"%"))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",reg_org_R2,"p =",reg_org_pval))
ordicluster(PCOA, cluster=hclust(as.dist(df)), choices=c(1,2),col="darkgray") #Hierachical clusering
orditorp(PCOA,display="sites",cex=0.8,air=0.01,col=reg_org, choices=c(1,2))
# Now dimensions 2 and 3
ordiplot(PCOA,type="n", main=paste(prefix,"PCoA-Hclust-org"), choices=c(1,3), xlab=paste("Dim1 - Variation explained:",PCOA.axes[1],"%"), ylab=paste("Dim3 - Variation explained:",PCOA.axes[3],"%"))
mtext(paste("Dimensions 1 and 2 - Adonis R2:",reg_org_R2,"p =",reg_org_pval))
ordicluster(PCOA, cluster=hclust(as.dist(df)), choices=c(1,3),col="darkgray") #Hierachical clusering
orditorp(PCOA,display="sites",cex=0.8,air=0.01,col=reg_org, choices=c(1,3))
# Now dimensions 2 and 3
ordiplot(PCOA,type="n", main=paste(prefix,"PCoA-Hclust-org"), choices=c(2,3), xlab=paste("Dim2 - Variation explained:",PCOA.axes[2],"%"), ylab=paste("Dim3 - Variation explained:",PCOA.axes[3],"%"))
mtext(paste("Dimensions 2 and 3 - Adonis R2:",reg_org_R2,"p =",reg_org_pval))
ordicluster(PCOA, cluster=hclust(as.dist(df)), choices=c(2,3),col="darkgray") #Hierachical clusering
orditorp(PCOA,display="sites",cex=0.8,air=0.01,col=reg_org, choices=c(2,3))
dev.off()

#DEPRECATED: PCA uses an euclidian matrix (e.g. a TSS+CLR matrix) and thus shouldn't be mixed with the distance matrices
##################### Euclidian-based MDS - Principal Component Analysis #############################
# # PCA <- rda.default(t(df))# in vegan, givig the function rda() dataframe without predictors runs a PCA very simiar to princomp
# ordiplot(PCA,type="n", main=paste(prefix,"PCoA-Hclust-org"), choices=c(1,2), xlab=paste("Dim1 - Variation explained:",PCOA.axes[1],"%"), ylab=paste("Dim2 - Variation explained:",PCOA.axes[2],"%"))
# orditorp(PCA,display="sites",cex=0.8,air=0.01,col=organ, choices=c(1,2))
