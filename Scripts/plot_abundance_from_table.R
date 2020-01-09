# Started on 2019-08-19
# by Rodrigo García-López for Adrian Ochoa-Leyva's Metagenomics and Metatranscriptomics Lab at IBt, UNAM, Cuernavaca, Mexico.
# Under GNU GPLv3 license
# Disclamer: This script was written for use with specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.
# This script was tested with R 3.6.0 - "Planting of a Tree"
# This script takes otu tables (output from converting biom tables) and calculates different total and unique items at different length frequence cutoffs
# Run as follows:
# cat table.tsv|~/bin/R-3.6.0/bin/Rscript plot_rel_abundance_from_table.R

# Tested with command:
# cat 07.1_pick_otus_open_ref/04_uchime_chimera_removal/05_tables_merged_per_group/03_joined_no_overlap/01_All_gg_vs_silva-all_items.tsv|~/bin/R-3.6.0/bin/Rscript plot_effect_of_freq_cutoff.R

df <- read.table(file('stdin'), sep="\t",header=T, skip=1, fileEncoding = "UTF-8", comment.char='')
df <- rowsum(df[2:length(df)], group=df[[1]]) #sum duplicate items by feature (row) name
df <- as.matrix(df[order(rowSums(df),decreasing=T),]) # sort input by most abundant first
rel <- sweep(df,2,colSums(df),'/') # create relative abundance table
# par(mfrow = c(1,2))
# color = sample(grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]) #create a randomized color vector: from the whole list of R colors, get shades of gray and randomize
# saveRDS(pick_colors, "color_set.rds") #Backup some good color combination
color=readRDS("color_set.rds")
pdf("Rel_abs_abundance.pdf")
max_size <- ifelse((1/ncol(df)*30)>1,1,1/ncol(df)*30)
barplot(rel, col=color, las=2, main="Relative abundance of features per sample",ylab="Feature proportion in sample", border=F,yaxt='n',cex.names=0.4,width=1,xaxt='n') # create the plot for the current table
axis(1, las = 2, at = seq(0.7,ncol(df)*1.2,1.2),labels=colnames(df),cex.axis= max_size) # axes are created separately
axis(2, las = 1, at = seq(0,1,0.05),cex.axis=0.8)
V3_total <- length(grep("^V3\\.",colnames(df))) #count the total items in each set
V4_total <- length(grep("^V4\\.",colnames(df)))
V3V4_total <- length(grep("^V3V4\\.",colnames(df)))
abline(v=c(0.1+V3_total*1.2,0.1+(V3_total+V4_total)*1.2),lwd=2) # create vertical division lines

# Print the legend
max_lab <- ifelse(nrow(df)>50,50,nrow(df)) #get the total number of labels to plot
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1) # create empty canvas
legend("topleft", legend = row.names(df)[1:max_lab], pch=15, pt.cex=0.6, cex=0.4, bty='n', col = color[1:max_lab]) # create the OTU labels for the first 50 items
mtext("Most abundant features", at=0.2, cex=1) #add some title
# and absolute abundance
barplot(df, col=color, las=2, main="Absolute abundance of features per sample",ylab="Feature proportion in sample", border=F,yaxt='n',cex.names=0.4,width=1,xaxt='n') # create the plot for the current table
axis(1, las = 2, at = seq(0.7,ncol(df)*1.2,1.2),labels=colnames(df),cex.axis= max_size) # axes are created separately
axis(2, las = 1, at = seq(0,max(colSums(df)),round(max(colSums(df))/20)),cex.axis=0.8)
abline(v=c(0.1+V3_total*1.2,0.1+(V3_total+V4_total)*1.2),lwd=2) # create vertical division lines
dev.off()

