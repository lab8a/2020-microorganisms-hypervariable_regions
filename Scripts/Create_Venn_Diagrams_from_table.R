# Started on 2019-08-06
# by Rodrigo García-López for Adrian Ochoa-Leyva's Metagenomics and Metatranscriptomics Lab at IBt, UNAM, Cuernavaca, Mexico.
# Under GNU GPLv3 license
# Disclamer: This script was written for use with specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.
# UPDATE: 2019-009-02: Some graphical changes for plotting. Added totals
# This script was tested with R 3.6.0 - "Planting of a Tree"
# It is intended to run from stdin by passing a contingency table (e.g. OTU table) in tsv format
# As it was created for biom-converted tables, the first line is skipped as it is just a comment (please adjust this if source is different). The first two lines have a # at the start of line.
# No taxonomy is present but can be adjusted with a [-1] over the loaded complete table.
# V3., V4., and V3V4. name identifiers are used for grouping.
# Run as follows:
# cat table.tsv|~/bin/R-3.6.0/bin/Rscript Create_Venn_Diagrams_from_table.R


df <- read.table(file('stdin'), sep="\t",header=T, skip=1, fileEncoding = "UTF-8", comment.char='')
df <- rowsum(df[2:length(df)], group=df[[1]])
# round(sum(colSums(df))*0.0001)
df <-df[order(rowSums(df),decreasing=T),]#sort by most abundant item (total abundance)
V3 <- grep("^V3\\.",colnames(df)) #get the columns for each set
V4 <- grep("^V4\\.",colnames(df))
V3V4 <- grep("^V3V4\\.",colnames(df))
sumV3 <- rowSums(df[,c(V3)]) #sum the items in each set
sumV4 <- rowSums(df[,c(V4)])
sumV3V4 <- rowSums(df[,c(V3V4)])
Totals_per_set <- cbind(sumV3,sumV4,sumV3V4) #and get a table with this
sumtotals <- rowSums(Totals_per_set) # get the total abundance per cluster
sumV3[sumV3 > 1] <- 7 #Identify those found in each set 7=V3, 17=V4, and 11=V3V4
sumV3V4[sumV3V4 > 1] <- 11
sumV4[sumV4 > 1] <- 17
unique_items <- rowSums(cbind(sumV3,sumV4,sumV3V4)) #now sum them
tab_unique=cbind("V3_only"=ifelse(unique_items==7, 1,0), "V4_only"=ifelse(unique_items==17, 1,0), "V3V4_only"=ifelse(unique_items==11, 1,0), "V3_&_V4"=ifelse(unique_items==24, 1,0), "V3_&_V3V4"=ifelse(unique_items==18, 1,0),"V4_&_V3V4"=ifelse(unique_items==28, 1,0),"All"=ifelse(unique_items==35, 1,0)) #Create a new table identifying each precise combination (determined by the sum in each row)
tab_total <- tab_unique*sumtotals # And create the same table with total abundances
# Create the set sums and intersection sizes for unique items (not used for the eulerr library method)
# # u_V3 <- sum(tab_unique[,"V3_only"])+sum(tab_unique[,"V3_&_V4"])+sum(tab_unique[,"V3_&_V3V4"])+sum(tab_unique[,"All"]) # sum the total unique items in each set
# # u_V4 <- sum(tab_unique[,"V4_only"])+sum(tab_unique[,"V3_&_V4"])+sum(tab_unique[,"V4_&_V3V4"])+sum(tab_unique[,"All"])
# # u_V3V4 <- sum(tab_unique[,"V3V4_only"])+sum(tab_unique[,"V3_&_V3V4"])+sum(tab_unique[,"V4_&_V3V4"])+sum(tab_unique[,"All"])
# # u_V3_n_V4 <- sum(tab_unique[,"V3_&_V4"])+sum(tab_unique[,"All"]) # And store all other intersections (must include the All set
# # u_V3_n_V3V4 <- sum(tab_unique[,"V3_&_V3V4"])+sum(tab_unique[,"All"])
# # u_V4_n_V3V4 <- sum(tab_unique[,"V4_&_V3V4"])+sum(tab_unique[,"All"])
# # u_all <- sum(tab_unique[,"All"])
# # # Now create the set sums and intersection sizes for item abundance
# # c_V3 <- sum(tab_total[,"V3_only"])+sum(tab_total[,"V3_&_V4"])+sum(tab_total[,"V3_&_V3V4"])+sum(tab_total[,"All"]) # sum the total unique items in each set
# # c_V4 <- sum(tab_total[,"V4_only"])+sum(tab_total[,"V3_&_V4"])+sum(tab_total[,"V4_&_V3V4"])+sum(tab_total[,"All"])
# # c_V3V4 <- sum(tab_total[,"V3V4_only"])+sum(tab_total[,"V3_&_V3V4"])+sum(tab_total[,"V4_&_V3V4"])+sum(tab_total[,"All"])
# # c_V3_n_V4 <- sum(tab_total[,"V3_&_V4"])+sum(tab_total[,"All"]) # And all other intersections (including the All set)
# # c_V3_n_V3V4 <- sum(tab_total[,"V3_&_V3V4"])+sum(tab_total[,"All"])
# # c_V4_n_V3V4 <- sum(tab_total[,"V4_&_V3V4"])+sum(tab_total[,"All"])
# # c_all <- sum(tab_total[,"All"])
# # #Test 1 with library VennDiagram
# # library(VennDiagram) # This is not a default library so be sure to install it before
# # overrideTriple = TRUE #This must be set to allow any combination of scaling
# # pdf("Test.pdf")
# # grid.newpage()
# # venn.plot <- draw.triple.venn(area1 = u_V3, area2 = u_V4, area3 = u_V3V4, n12 = u_V3_n_V4, n13 = u_V3_n_V3V4,  n23 = u_V4_n_V3V4, n123 = u_all, category = c('V3', 'V4', 'V5'), fill = c('brown2', 'chartreuse1', 'cyan'), cat.col = c('brown2', 'chartreuse1', 'cyan'), cex = 1, cat.cex = 1, euler.d = TRUE, scaled = TRUE) #First the unique
# # dev.off()
# # There is no exact solution for circular overlaps

# For using other libraries, we need to calculate each section separately for both the unique and total items
u_V3 <- sum(tab_unique[,"V3_only"])
u_V4 <- sum(tab_unique[,"V4_only"])
u_V3V4 <- sum(tab_unique[,"V3V4_only"])
u_V3_n_V4 <- sum(tab_unique[,"V3_&_V4"])
u_V3_n_V3V4 <- sum(tab_unique[,"V3_&_V3V4"])
u_V4_n_V3V4 <- sum(tab_unique[,"V4_&_V3V4"])
u_all <- sum(tab_unique[,"All"])
c_V3 <- sum(tab_total[,"V3_only"])
c_V4 <- sum(tab_total[,"V4_only"])
c_V3V4 <- sum(tab_total[,"V3V4_only"])
c_V3_n_V4 <- sum(tab_total[,"V3_&_V4"])
c_V3_n_V3V4 <- sum(tab_total[,"V3_&_V3V4"])
c_V4_n_V3V4 <- sum(tab_total[,"V4_&_V3V4"])
c_all <- sum(tab_total[,"All"])
uSum_V3=u_V3+u_V3_n_V4+u_V3_n_V3V4+u_all
uSum_V4=u_V4+u_V3_n_V4+u_V4_n_V3V4+u_all
uSum_V3V4=u_V3V4+u_V3_n_V3V4+u_V4_n_V3V4+u_all
cSum_V3=c_V3+c_V3_n_V4+c_V3_n_V3V4+c_all
cSum_V4=c_V4+c_V3_n_V4+c_V4_n_V3V4+c_all
cSum_V3V4=c_V3V4+c_V3_n_V3V4+c_V4_n_V3V4+c_all
#Test 2 with venneuler
# # library(venneuler)  #This library requires the exact non-shared sections instead
# # MyVenn <- venneuler(c(A = u_V3, B = u_V4, C = u_V3V4, "A&B" = u_V3_n_V4, "A&C" = u_V3_n_V3V4, "B&C" = u_V4_n_V3V4, "A&B&C" = u_all))
# # MyVenn$labels <- c("V3","V4","V3V4")
# # pdf("Test_2.pdf")
# # plot(MyVenn)
# # dev.off()
# # The "All" overlap is empty. There is no circular solution

# # #Test 3 with eulerr (allows ellipses for exact overlap)
library("eulerr") #This library requires the exact non-shared sections instead
vd <- euler(c(V3 = u_V3, V4 = u_V4, V3V4 = u_V3V4, "V3&V4" = u_V3_n_V4, "V3&V3V4" = u_V3_n_V3V4, "V4&V3V4" = u_V4_n_V3V4, "V3&V4&V3V4" = u_all),shape = "ellipse") # Ellipses are allowed for exact overlaps
while(sum(vd$residuals)>1){ #it may be necesary to repeat until no residuals are detected (meaning successful plot)
	vd <- euler(c(V3 = u_V3, V4 = u_V4, V3V4 = u_V3V4, "V3&V4" = u_V3_n_V4, "V3&V3V4" = u_V3_n_V3V4, "V4&V3V4" = u_V4_n_V3V4, "V3&V4&V3V4" = u_all),shape = "ellipse")
	}
pdf("Venn_unique.pdf")
plot(vd, key = TRUE, counts = TRUE, main=list(label="Unique items", fontsize=20),fills = list(fill=c('coral1', 'cornflowerblue', 'darkturquoise'),alpha=0.9), edges = TRUE, quantities = list(fontsize = 20,col="black"), col="white", lwd=3, labels=list(col="White", fontsize=40), list(labels=c(uSum_V3,uSum_V4,uSum_V3V4),fontsize=20))
dev.off()# # A good solution is met after several tries, now with totals
vd2 <- euler(c(V3 = c_V3, V4 = c_V4, V3V4 = c_V3V4, "V3&V4" = c_V3_n_V4, "V3&V3V4" = c_V3_n_V3V4, "V4&V3V4" = c_V4_n_V3V4, "V3&V4&V3V4" = c_all),shape = "ellipse") # Ellipses are allowed for exact overlaps
while(sum(vd$residuals)>1){ #it may be necesary to repeat until no residuals are detected (meaning successful plot)
	vd2 <- euler(c(V3 = c_V3, V4 = c_V4, V3V4 = c_V3V4, "V3&V4" = c_V3_n_V4, "V3&V3V4" = c_V3_n_V3V4, "V4&V3V4" = c_V4_n_V3V4, "V3&V4&V3V4" = c_all),shape = "ellipse")
	}
pdf("Venn_totals.pdf")
plot(vd2, key = TRUE, counts = TRUE, main=list(label="Total items", fontsize=20),fills = list(fill=c('coral1', 'cornflowerblue', 'darkturquoise'),alpha=0.9), edges = TRUE, quantities = list(fontsize = 15,col="black"), col="white", lwd=3, labels=list(col="white", fontsize=20), list(labels=c(cSum_V3,cSum_V4,cSum_V3V4),fontsize=20))
dev.off()
