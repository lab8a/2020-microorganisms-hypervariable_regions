# Started on 2019-08-11
# by Rodrigo García-López for Adrian Ochoa-Leyva's Metagenomics and Metatranscriptomics Lab at IBt, UNAM, Cuernavaca, Mexico.
# Under GNU GPLv3 license
# Disclamer: This script was written for use with specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.
# This script was tested with R 3.6.0 - "Planting of a Tree"
# It is intended to run from stdin by passing a contingency table (e.g. OTU table) in tsv format
# As it was created for biom-converted tables, the first line is skipped as it is just a comment (please adjust this if source is different). The first two lines have a # at the start of line.
# First line is expected  to be thet taxonomy or the cluster names
# The output:
# - 1) plot of total items per sample vs quantiles 1-100
# - 2) plot of total items per sample/set for region comparison
# - 3) tsv table with the info to plot 2)
# - 4) recruitment plots of the original table rarefied at a fixed depth, also painted by organ and region
# - 5) sample recruitment by organ and region

# The library "vegan" should be installed to create the recruitment plot.
# Run as follows:
# cat table.tsv|~/bin/R-3.6.0/bin/Rscript create_recruitment_plot_from_table.R

# Tested with df <- read.table("10_Reduced_sparsity_tables/OTUs_97/gg/01_All_gg97-table_lvl7.tsv",sep="\t", header=T, row.names=1, fileEncoding = "UTF-8",comment.char='',skip=1)

df <- read.table(file('stdin'),sep="\t", header=T, row.names=1, fileEncoding = "UTF-8",comment.char='',skip=1)
quantiles_100 <- quantile(colSums(df),probs=seq(0,1,0.01)) #calculate quantiles 1-100 for the total sum per sample
f_gt_10k <- sort(colSums(df))[sort(colSums(df))>10000][1] #first item gt 10k
q10k <- grep(quantiles_100[quantiles_100>f_gt_10k][1],quantiles_100)

pdf("Sample_sum_quantiles.pdf")
plot(quantiles_100, type="l",lwd=2, main="Quantiles of the Sample totals", ylab="Feature sum per sample", xlab="Quantiles 1 - 100",xaxt='n',yaxt='n');axis(1, las=2, at=seq(0,100,5));axis(2, las=1, at=seq(min(quantiles_100),max(quantiles_100),floor((max(quantiles_100)-min(quantiles_100))/20)),cex.axis=0.8);abline(h=c(quantiles_100[6],quantiles_100[11],quantiles_100[q10k]),lty=2:4);abline(v=c(5,10,q10k),lty=2:4);legend("top", lty=2:4,legend=c(paste("Quantile 05 % Freq:",quantiles_100[6]),paste("Quantile 10 % Freq:",quantiles_100[11]),paste("Quantile",q10k,"% Freq:",quantiles_100[q10k])))
dev.off()

q5=quantiles_100[6]
q10=quantiles_100[11]


total_samples <- colSums(df[2:length(df)]) #store the total items
organ=region=NULL # create 2 empty objects
# Now fill each vector with an identifier for type
organ[grep("I", names(total_samples))]="H"
organ[grep("H", names(total_samples))]="I"
region[grep("V3\\.", names(total_samples))]="V3"
region[grep("V4\\.", names(total_samples))]="V4"
region[grep("V3V4\\.", names(total_samples))]="V3V4"
# Now some internal information for plotting
pch=col=pos=NULL # create 2 empty objects
col[grep("I", names(total_samples))]="cornflowerblue"
col[grep("H", names(total_samples))]="aquamarine3"
pch[grep("V3\\.", names(total_samples))]=0
pch[grep("^V4\\.", names(total_samples))]=2
pch[grep("V3V4\\.", names(total_samples))]=1
pos[grep("V3\\.H", names(total_samples))]=0
pos[grep("V3\\.I", names(total_samples))]=10
pos[grep("^V4\\.H", names(total_samples))]=20
pos[grep("^V4\\.I", names(total_samples))]=30
pos[grep("V3V4\\.H", names(total_samples))]=40
pos[grep("V3V4\\.I", names(total_samples))]=50

# Graph the different sample totals grouped by region and organ
df2 <-data.frame("Items"=total_samples,"Organ"=organ,"Region"=region,"pch"=pch,"col"=col,"pos"=pos)
ticks <- seq(min(df2$Items),max(df2$Items),floor((max(df2$Items)-min(df2$Items))/20))
pdf("Sample_sum_per_set.pdf")
plot(df2$pos, df2$Items, col=as.character(df2$col), pch=df2$pch, main="Clustering/Denoising techniques, unique OTUs/ASVs", xlab="", ylab="Total items", cex=2, xlim=c(-5,55),xaxt='n',yaxt='n');abline(v=c(15,35),col="darkgray", lty=9);mtext("Green = Hepatopancreas, Blue = Gut");axis(2, at=ticks,las=1, cex.axis=0.7, labels=format(ticks,scientific=F));axis(1, at=seq(0,50,10), labels=c("V3-H","V3-I","V4-H","V4-I","V3V4-H","V3V4-I"));abline(h=10000,col="darkgray")
dev.off()
write.table(df2, "Sample_sum_per_set-table.tsv",quote=F,sep="\t",col.names=NA)

library("vegan") # For the recruitment plot, we require the vegan library
trans <- t(df[order(rowSums(df),decreasing=T),]) # vegan expects samples as rows, so we transpose the table sorted by decreasing abundance total (clusters)
lentrans <- length(row.names(trans)) #get the total number of items
# Create set_only tables (removing items not in the set) and their lengths


# Now, we'll start with rarefactions of all features 
raremax <- min(rowSums(trans)) # We will first set the raremax to the sample with the smallest number of total items
# Assess how well the min sample items, the quantile 5, and 10, and the 10k items produce rarefactions that resemble the actual species in the set
features <- specnumber(trans) #calculate the total unique items
pdf("rarefaction_cutoff_adjustment.pdf")
par(mfrow = c(2,2))
rarefied <- rarefy(trans,raremax);plot(features, rarefied, xlab = "Total Features", ylab = "Rarefied Features", main = "Total vs Rarefied");mtext(paste("Min sample:",raremax));abline(0, 1)
rarefied <- rarefy(trans,q5);plot(features, rarefied, xlab = "Total Features", ylab = "Rarefied Features", main = "Total vs Rarefied Features");mtext(paste("Quantile 05%:",q5));abline(0, 1)
rarefied <- rarefy(trans,q10);plot(features, rarefied, xlab = "Total Features", ylab = "Rarefied Features", main = "Total vs Rarefied Features");mtext(paste("Quantile 10%:",q10));abline(0, 1)
rarefied <- rarefy(trans,10000);plot(features, rarefied, xlab = "Total Features", ylab = "Rarefied Features", main ="Total vs Rarefied Features");mtext("10000 items");abline(0, 1)
par(mfrow = c(1,1))
dev.off()

# Create several permutations, one combination for a single line (color and width)
lty <- c("solid","dashed","dotted","dotdash","longdash","twodash")
col <- c('chartreuse3', 'cornflowerblue', 'darkgoldenrod1', 'peachpuff3','mediumorchid2', 'turquoise3', 'wheat4','slategray2',"black","burlywood4","aquamarine2","blue2","violetred2","palegreen3","wheat1","magenta1","limegreen","darkorange2","darkgray")
lwd <- c(1)
pars <- expand.grid(col = col, lty = lty, lwd = lwd, stringsAsFactors = FALSE) # This create all permutations, adjust to number of cases required
out <- with(pars[1:nrow(trans), ],rarecurve(trans, step = 20, col = col, lty = lty, label = F))
Nmax <- sapply(out, function(x) max(attr(x, "Subsample"))) #extract vector "Subsample", an attribute of each list in the out df and get the max value (last point in the curve per sample)
Smax <- sapply(out, max) #get the straightforward max (max observed species per sample)

# Now plot everything
pdf("Raw_recruitment-All.pdf")
plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = "Sample Size", ylab = "Unique Features", type = "n",main="Recruitment plot",xaxt='n',yaxt='n') #create empty plot with max rarefaction count and max observed items
mtext("All items shown")
axis(1, at=c(seq(0,max(Nmax),round(max(Nmax)/20))),las=2, cex.axis=0.8, labels=format(seq(0,max(Nmax),round(max(Nmax)/20)),scientific=F))
axis(2, at=c(seq(1,max(Smax),floor(max(Smax)/20))),las=1, cex.axis=0.8, labels=format(seq(1,max(Smax),floor(max(Smax)/20)),scientific=F))
for (i in seq_along(out)) { # advance i for any number of items in out (samples), or else "for each sample"
	N <- attr(out[[i]], "Subsample") #get the rarefaction depth list
	with(pars, lines(N, out[[i]], col = col[i], lty = lty[i], lwd = lwd[i])) #carry out the drawing of each line  using x as the depth, and using each corresponding item in pars (col, lty and lwd item for parameters)
}
abline(v=c(raremax,q5,q10,10000),lty=c(2,3,4,5),col="darkgray") # draw lines at possible cutoffs
legend("topright", legend = c(paste("Min sample:",raremax), paste("Quantile 5:",q5),paste("Quantile 10:",q10),"10000 items"), col = "darkgray", pch = NA, pt.cex=1, cex=.75, bg="white", lty=c(2,3,4,5))
dev.off()

# Now a zoom at 20k
pdf("Raw_recruitment-20K.pdf")
plot(c(1, 20000), c(1, max(Smax)), xlab = "Rarefaction depth", ylab = "Unique Features", type = "n",xaxt='n',yaxt='n',main="Recruitment plot")
mtext("Showing recruitment 1-20000")
axis(1, at=c(seq(0,20000,1000)),las=2, cex.axis=0.8, labels=format(seq(0,20000,1000),scientific=F))
axis(2, at=c(seq(1,max(Smax),floor(max(Smax)/20))),las=1, cex.axis=0.8, labels=format(seq(1,max(Smax),floor(max(Smax)/20)),scientific=F))
for (i in seq_along(out)) { # advance i for any number of items in out (samples), or else "for each sample"
	N <- attr(out[[i]], "Subsample") #get the rarefaction depth list
	with(pars, lines(N, out[[i]], col = col[i], lty = lty[i], lwd = lwd[i])) #carry out the drawing of each line  using x as the depth, and using each corresponding item in pars (col, lty and lwd item for parameters)
}
abline(v=c(raremax,q5,q10,10000),lty=c(2,3,4,5),col="darkgray")
legend("topright", legend = c(paste("Min sample:",raremax), paste("Quantile 5:",q5),paste("Quantile 10:",q10),"10000 items"), col = "darkgray", pch = NA, pt.cex=1, cex=.75, bg="white", lty=c(2,3,4,5))
dev.off()

# Now, we'll change the colors, for this, we need the different groups
# First by region
region <- names(df)
region[grep('^V3\\.',region)]="coral1"
region[grep('^V4\\.',region)]="cornflowerblue"
region[grep('^V3V4\\.',region)]="turquoise3"
region <- as.character(region)

# then by organ
organ <- names(df)
organ[grep('H',organ)]="darkorchid1"
organ[grep('I',organ)]="chartreuse2"
organ <- as.character(organ) 

organt <- names(df)
organt[grep('H',organt)]=1
organt[grep('I',organt)]=2
organt <- as.numeric(organt) 

pars <- data.frame(run = region, org = organ)

### Now run by color
pdf("Recruitment_region-All.pdf")
plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = "Rarefaction depth", ylab = "Unique features", type = "n",xaxt='n',yaxt='n',main="Recruitment plot")
mtext("Colored by region")
axis(1, at=c(seq(0,max(Nmax),round(max(Nmax)/20))),las=2, cex.axis=0.8, labels=format(seq(0,max(Nmax),round(max(Nmax)/20)),scientific=F))
axis(2, at=c(seq(1,max(Smax),floor(max(Smax)/20))),las=1, cex.axis=0.8, labels=format(seq(1,max(Smax),floor(max(Smax)/20)),scientific=F))
for (i in seq_along(out)) { # advance i for any number of items in out (samples), or else "for each sample"
	N <- attr(out[[i]], "Subsample") #get the rarefaction depth list
	with(pars, lines(N, out[[i]], col = as.character(run[i]), lty = organt[i], lwd = 1.5))
}
abline(v=c(raremax,q5,q10,10000),lty=c(2,3,4,5),col="darkgray")
legend("topright", legend = c(paste("Min sample:",raremax), paste("Quantile 5:",q5),paste("Quantile 10:",q10),"10000 items","V3-H","V4-H","V3V4-H","V3-I","V4-I","V3V4-I"), col = c(rep("darkgray",4),"coral1","cornflowerblue","turquoise3","coral1","cornflowerblue","turquoise3"), pch = NA, pt.cex=1, cex=.75, bg="white", lty=c(2,3,4,5,1,1,1,2,2,2),lwd=1.5)
dev.off()

#now a zoom at 20k
pdf("Recruitment_region-20K.pdf")
plot(c(1, 20000), c(1, max(Smax)), xlab = "Rarefaction depth", ylab = "Unique features", type = "n",xaxt='n',yaxt='n',main="Recruitment plot")
mtext("Colored by region - Showing 1-20000")
axis(1, at=c(seq(0,20000,1000)),las=2, cex.axis=0.8, labels=format(seq(0,20000,1000),scientific=F))
axis(2, at=c(seq(1,max(Smax),floor(max(Smax)/20))),las=1, cex.axis=0.8, labels=format(seq(1,max(Smax),floor(max(Smax)/20)),scientific=F))
for (i in seq_along(out)) { # advance i for any number of items in out (samples), or else "for each sample"
	N <- attr(out[[i]], "Subsample") #get the rarefaction depth list
	with(pars, lines(N, out[[i]], col = as.character(run[i]), lty = organt[i], lwd = 1.5))
}
abline(v=c(raremax,q5,q10,10000),lty=c(2,3,4,5),col="darkgray")
legend("topright", legend = c(paste("Min sample:",raremax), paste("Quantile 5:",q5),paste("Quantile 10:",q10),"10000 items","V3-H","V4-H","V3V4-H","V3-I","V4-I","V3V4-I"), col = c(rep("darkgray",4),"coral1","cornflowerblue","turquoise3","coral1","cornflowerblue","turquoise3"), pch = NA, pt.cex=1, cex=.75, bg="white", lty=c(2,3,4,5,1,1,1,2,2,2),lwd=1.5)
dev.off()

#now by organ
pdf("Recruitment_organ-All.pdf")
plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = "Rarefaction depth", ylab = "Unique features", type = "n",xaxt='n',yaxt='n',main="Recruitment plot")
mtext("Colored by organ")
axis(1, at=c(seq(0,max(Nmax),round(max(Nmax)/20))),las=2, cex.axis=0.8, labels=format(seq(0,max(Nmax),round(max(Nmax)/20)),scientific=F))
axis(2, at=c(seq(1,max(Smax),floor(max(Smax)/20))),las=1, cex.axis=0.8, labels=format(seq(1,max(Smax),floor(max(Smax)/20)),scientific=F))
for (i in seq_along(out)) { # advance i for any number of items in out (samples), or else "for each sample"
	N <- attr(out[[i]], "Subsample") #get the rarefaction depth list
	with(pars, lines(N, out[[i]], col = as.character(organ[i]), lty = 1, lwd = 1.5))
}
abline(v=c(raremax,q5,q10,10000),lty=c(2,3,4,5),col="darkgray")
legend("topright", legend = c(paste("Min sample:",raremax), paste("Quantile 5:",q5),paste("Quantile 10:",q10),"10000 items","Hepatopancreas","Gut"), col = c(rep("darkgray",4),"darkorchid1","chartreuse2"), pch = NA, pt.cex=1, cex=.75, bg="white", lty=c(2,3,4,5,1,1), lwd=1.5)
dev.off()

#now a zoom at 20k
pdf("Recruitment_organ-20K.pdf")
plot(c(1, 20000), c(1, max(Smax)), xlab = "Rarefaction depth", ylab = "Unique features", type = "n",xaxt='n',yaxt='n',main="Recruitment plot")
mtext("Colored by organ - Showing 1-20000")
axis(1, at=c(seq(0,20000,1000)),las=2, cex.axis=0.8, labels=format(seq(0,20000,1000),scientific=F))
axis(2, at=c(seq(1,max(Smax),floor(max(Smax)/20))),las=1, cex.axis=0.8, labels=format(seq(1,max(Smax),floor(max(Smax)/20)),scientific=F))
for (i in seq_along(out)) { # advance i for any number of items in out (samples), or else "for each sample"
	N <- attr(out[[i]], "Subsample") #get the rarefaction depth list
	with(pars, lines(N, out[[i]], col = as.character(organ[i]), lty = 1, lwd = 1.5))
}
abline(v=c(raremax,q5,q10,10000),lty=c(2,3,4,5),col="darkgray")
legend("topright", legend = c(paste("Min sample:",raremax), paste("Quantile 5:",q5),paste("Quantile 10:",q10),"10000 items","Hepatopancreas","Gut"), col = c(rep("darkgray",4),"darkorchid1","chartreuse2"), pch = NA, pt.cex=1, cex=.75, bg="white", lty=c(2,3,4,5,1,2), lwd=1.5)
dev.off()

V3set <- trans[grep("^V3\\.", row.names(trans)),];V3set <- V3set[,colSums(V3set)>0]
V4set <- trans[grep("^V4\\.", row.names(trans)),];V4set <- V4set[,colSums(V4set)>0]
V3V4set <- trans[grep("^V3V4\\.", row.names(trans)),];V3V4set <- V3V4set[,colSums(V3V4set)>0]
lenV3 <- nrow(V3set)
lenV4 <- nrow(V4set)
lenV3V4 <- nrow(V3V4set)
# Now get 10k permutations per set
perm=10000 #number of permutations for calculations
item_acc_V3 <-specaccum(V3set,method="random",permutations=perm)
item_acc_V4 <-specaccum(V4set,method="random",permutations=perm)
item_acc_V3V4 <-specaccum(V3V4set,method="random",permutations=perm)
max_OTUs <- max(ncol(V3set),ncol(V4set),ncol(V3V4set)) # get the max number of clusters in any set
pdf("Recruitment_per_sample-region.pdf")
plot(item_acc_V3,ci.type="poly", col="coral1", lwd=2, ci.lty=0, ci.col=rgb(1,0.5,0,alpha=0.3),ylim=c(0,max_OTUs),xaxt='n',yaxt='n',main="Whole sample recruitment on observed features - Per region", xlab="Samples recruited", ylab="Observed features");axis(1, las=2, at=seq(0,lenV3,1),cex.axis=0.8);axis(2, las=1, at=seq(0,max_OTUs,max_OTUs/20), cex.axis=0.8);mtext("Calculated on 10K permutations");legend("bottomright",lty=1, legend=c("V3","V4","V3V4"), col=c("coral1","cornflowerblue","turquoise3"),lwd=2)
lines(item_acc_V4,ci.type="poly", col="cornflowerblue", lwd=2, ci.lty=0, ci.col=rgb(0,0.5,1,alpha=0.3))
lines(item_acc_V3V4,ci.type="poly", col="turquoise3", lwd=2, ci.lty=0, ci.col=rgb(0,1,0.5,alpha=0.3))
dev.off()

# Now with organs
Hset <- trans[grep("H", row.names(trans)),];Hset <- Hset[,colSums(Hset)>0]
Iset <- trans[grep("I", row.names(trans)),];Iset <- Iset[,colSums(Iset)>0]
lenH <- nrow(Hset)
lenI <- nrow(Iset)
# Now get 10k permutations per set
perm=10000 #number of permutations for calculations
item_acc_H <-specaccum(Hset,method="random",permutations=perm)
item_acc_I <-specaccum(Iset,method="random",permutations=perm)
max_OTUs <- max(ncol(Hset),ncol(Iset)) # get the max number of clusters in any set
max_OTUs <- max_OTUs*1.1
pdf("Recruitment_per_sample-organ.pdf")
plot(item_acc_H,ci.type="poly", col="darkorchid1", lwd=2, ci.lty=0, ci.col=rgb(0.5,0,0.5,alpha=0.3),ylim=c(0,max_OTUs),xaxt='n',yaxt='n',main="Whole sample recruitment on observed features - Per organ", xlab="Samples recruited", ylab="Observed features");axis(1, las=2, at=seq(0,lenH,1),cex.axis=0.8);axis(2, las=1, at=seq(0,max_OTUs,max_OTUs/20), cex.axis=0.7);mtext("Calculated on 10K permutations");legend("bottomright",lty=1, legend=c("Hepatopancreas","Gut"), col=c("darkorchid1","chartreuse"),lwd=2)
lines(item_acc_I,ci.type="poly", col="chartreuse", lwd=2, ci.lty=0, ci.col=rgb(0,1,0,alpha=0.3))
dev.off()
