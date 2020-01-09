# Started on 2019-08-07
# by Rodrigo García-López for Adrian Ochoa-Leyva's Metagenomics and Metatranscriptomics Lab at IBt, UNAM, Cuernavaca, Mexico.
# Under GNU GPLv3 license
# Disclamer: This script was written for use with specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.
# This script was tested with R 3.6.0 - "Planting of a Tree"
# This takes tables the tables containing the total clusters per sample considering all items, only ref-based and only de novo clusters, both for the gg and the silva references.
# The input tables were created with the R script contingency_table-unique_items_per_sample.R
# The Rscript bin must to be executed must be the one in my local folder as it has the readr library installed
# Run as follows:
# cat table.tsv|~/bin/R-3.6.0/bin/Rscript barplots_cluster_compositions.R

# Tested with command:
# cat 07.1_pick_otus_open_ref/04_uchime_chimera_removal/05_tables_merged_per_group/03_joined_no_overlap/01_All_gg_vs_silva-all_items.tsv|~/bin/R-3.6.0/bin/Rscript barplots_cluster_compositions.R

df <- t(read.table(file('stdin'),sep="\t", header=T, row.names=1, fileEncoding = "UTF-8",comment.char='',skip=0)) #Store the complete table transposed
unique_table <- df[grep("^Uniques",rownames(df)),] #get the unique items only
unique_table <- unique_table[-nrow(unique_table),] #remove the sum (last row)
# First, process the unique sets, starting with greengenes
gg <- t(unique_table[,c(3,5)]) #retain only the greengenes table info
V3_total <- length(grep("_V3\\.",colnames(gg))) #count the total items in each set
V4_total <- length(grep("_V4\\.",colnames(gg)))
V3V4_total <- length(grep("_V3V4\\.",colnames(gg)))
pdf("barplots_unique_gg.pdf")
barplot(gg,width=1, space=0, main="Greengenes OTU assignations", col=c("turquoise","coral"),xaxt='n',yaxt='n') #plot using stacked plots 
# axis(1, at=c(V3_total/2,V3_total+V4_total/2,V3_total+V4_total+V3V4_total/2),las=1, cex.axis=2, labels=c("V3","V4","V3V4"))
item_names <- unlist(strsplit(colnames(gg), "_"))[grep("V",unlist(strsplit(colnames(gg), "_")))]
axis(1, at=seq(0.5,length(item_names),1),las=2, cex.axis=1, labels=item_names)
top <- ceiling(max(colSums(gg))/1000)*1000 #get the max height in the plot
axis(2, at=seq(0,top, top/20), las=2) # use it for ploting the y axis
legend("topleft", legend = rownames(gg), pch = 15, col=c("turquoise","coral"), pt.cex=1, cex=0.75, bg="white")
abline(v=c(V3_total,V3_total+V4_total),lwd=3)
dev.off()
# and for silva
gg <- t(unique_table[,c(4,6)]) #retain only the silva table info
V3_total <- length(grep("_V3\\.",colnames(gg)))
V4_total <- length(grep("_V4\\.",colnames(gg)))
V3V4_total <- length(grep("_V3V4\\.",colnames(gg)))
pdf("barplots_unique_silva.pdf")
barplot(gg,width=1, space=0, main="Silva OTU assignations", col=c("turquoise","coral"),xaxt='n',yaxt='n')
# axis(1, at=c(V3_total/2,V3_total+V4_total/2,V3_total+V4_total+V3V4_total/2),las=1, cex.axis=2, labels=c("V3","V4","V3V4"))
item_names <- unlist(strsplit(colnames(gg), "_"))[grep("V",unlist(strsplit(colnames(gg), "_")))]
axis(1, at=seq(0.5,length(item_names),1),las=2, cex.axis=1, labels=item_names)
top <- ceiling(max(colSums(gg))/1000)*1000
axis(2, at=seq(0,top, top/20), las=2)
legend("topleft", legend = rownames(gg), pch = 15, col=c("turquoise","coral"), pt.cex=1, cex=0.75, bg="white")
abline(v=c(V3_total,V3_total+V4_total),lwd=3)
dev.off()

# now the totals in gg
totals_table <- df[grep("^Total",rownames(df)),]
totals_table <- totals_table[-nrow(totals_table),]
gg <- t(totals_table[,c(3,5)])
V3_total <- length(grep("_V3\\.",colnames(gg)))
V4_total <- length(grep("_V4\\.",colnames(gg)))
V3V4_total <- length(grep("_V3V4\\.",colnames(gg)))
pdf("barplots_totals_gg.pdf")
barplot(gg,width=1, space=0, main="Greengenes OTU assignations", col=c("turquoise","coral"),xaxt='n',yaxt='n')
# axis(1, at=c(V3_total/2,V3_total+V4_total/2,V3_total+V4_total+V3V4_total/2),las=1, cex.axis=2, labels=c("V3","V4","V3V4"))
item_names <- unlist(strsplit(colnames(gg), "_"))[grep("V",unlist(strsplit(colnames(gg), "_")))]
axis(1, at=seq(0.5,length(item_names),1),las=2, cex.axis=1, labels=item_names)
top <- ceiling(max(colSums(gg))/1000)*1000 #get the max height in the plot
axis(2, at=seq(0,top, top/20), las=2) # use it for ploting the y axis
abline(v=c(V3_total,V3_total+V4_total),lwd=3)
legend("topright", legend = rownames(gg), pch = 15, col=c("turquoise","coral"), pt.cex=1, cex=0.75, bg="white")
dev.off()
# and for silva
gg <- t(totals_table[,c(4,6)])
V3_total <- length(grep("_V3\\.",colnames(gg)))
V4_total <- length(grep("_V4\\.",colnames(gg)))
V3V4_total <- length(grep("_V3V4\\.",colnames(gg)))
pdf("barplots_totals_silva.pdf")
barplot(gg,width=1, space=0, main="Silva OTU assignations", col=c("turquoise","coral"),xaxt='n',yaxt='n')
# axis(1, at=c(V3_total/2,V3_total+V4_total/2,V3_total+V4_total+V3V4_total/2),las=1, cex.axis=2, labels=c("V3","V4","V3V4"))
item_names <- unlist(strsplit(colnames(gg), "_"))[grep("V",unlist(strsplit(colnames(gg), "_")))]
axis(1, at=seq(0.5,length(item_names),1),las=2, cex.axis=1, labels=item_names)
top <- ceiling(max(colSums(gg))/1000)*1000 #get the max height in the plot
axis(2, at=seq(0,top, top/20), las=2) # use it for ploting the y axis
abline(v=c(V3_total,V3_total+V4_total),lwd=3)
legend("topright", legend = rownames(gg), pch = 15, col=c("turquoise","coral"), pt.cex=1, cex=0.75, bg="white")
dev.off()
