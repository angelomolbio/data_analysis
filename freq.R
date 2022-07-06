#Make this script working and then modify the script to plot the cutoff.val only using only 3 colors.
tmp <- read.csv("cluster_1_singleton_positive_markers_ranked.csv", header=T, row.names=1)
df <- data.frame(tmp$TP, tmp$TN, tmp$cutoff_val)
names(df) <- c("TP", "TN", "cutoff.val");
library(ggplot2)
p <- ggplot(df[which(df[,1] > 0.9 & df[,2] > 0.9),], aes(x=TP, y=TN, color=cutoff.val)) + geom_point()
p + scale_color_gradientn(colours = c("red","blue","yellow"))

#Convert this script in a function

newfunction<- function (saver.annotated, cell.cycle){

#adding CC to cell barcode
cc <- read.table(cell.cycle, sep=",", header=T, row.names=1, stringsAsFactors=F)
rnames.cc <- rownames(cc)
cc <- as.character(unlist(cc))
names(cc) <- rnames.cc
cc <- cc[order(names(cc))]
saver <- read.table(saver.annotated, sep=",", header=T, row.names=1, , stringsAsFactors=F)
saver <- saver[, order(names(saver))]
common <- intersect(names(cc), names(saver))
cc <- cc[which(names(cc) %in% common)]
saver <- saver[,which(names(saver) %in% common)]

if(identical(names(cc), names(saver))){
  names(saver) <- paste(names(saver), cc,sep="_")
  write.table(saver, "saver_cc.csv", sep=",", col.names=NA)
}else{ cat("\nError!\n")}
}

saver.annotated <- "matrix.csv"
cell.cycle <- "cellCycle.csv"
newfunction(saver.annotated,cell.cycle)

head(read.csv("saver_cc.csv",header=TRUE,  sep=",", row.names = 1))

#Write a script to indicate which of the cytobands has 
#the highest number of genes associated and save in a comma 
#delimited file the cytoband and the genes

install.packages("dplyr")

hg38 <- (read.csv("hg38.p13_105.csv",header=TRUE,  sep=","))
table<-data.frame(table(hg38['target']))
head(table)
class(table)
table[order(table$Freq,decreasing = TRUE),]
dim(table)

tabella<-data.frame(hg38)
class(tabella)

tabellafinale <- tabella[which(tabella[,2] == "chr19-p13.3"),]
head (tabellafinale)
write.csv (tabellafinale, "tabellina.csv")
