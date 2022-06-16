library(devtools)
install_github("kendomaniac/docker4seq", ref="master")
library(docker4seq)
downloadContainers(group="docker")
install_github("kendomaniac/rCASC", ref="master")
library(rCASC)
downloadContainers(group="docker")

#Fastqcsetting 
#http://users.path.ox.ac.uk/~pcook/w1/NGS_workshop2.htm
#sudo ln -s /home/angeloscarciglia/Desktop/calogero/FastQC/fastqc /usr/local/bin/fastqc
#chmod + fastqc in FastQC folder downloaded
#./fastqc
# in cd dataset0 cartella del file -> fastqc ra100k.R1.fastq

#MultiQC directly from the terminal with docker container 
#cd folder containing the file 
#docker run -t -v `pwd`:`pwd` -w `pwd` ewels/multiqc

#create image of ubuntu with specific requisiti
#docker pull ubuntu:21.04
#docker run -i -t ubuntu:21.04
#apt-get update/ upgrade/ apt-get install sudo / apt-get install nano
#once done exit and docker commit "istance in which we have done what we need" ubuntu:21.04.angelo (this final point provided the tag included)

#create an R file with the command install.packages('umap', repos='http://cran.us.r-project.org')#saved as umap_command.R in the same folder of Dockerfile
#create the dockerfile and the container from it 
#mkdir DockerFile (generate the folder)
#touch Dockerfile (generate the file with no extension)
#vim Dockerfile (go in the file)
#FROM rocker/r-ubuntu:20.04	
#RUN apt-get -y update
#RUN apt-get -y upgrade
#RUN apt-get -y install sudo 
#RUN apt-get -y install nano
#COPY ./umap_command.R/ home/
#RUN chmod +x /home/umap_command.R
#RUN sudo apt-get update && sudo apt-get -y install libssl-dev

#i to write and esc :wq! to exit saving
###After creating the Dockerfile goes on the terminal and type the command below.
#docker build Dockerfile -t r4:v.0.01 #build the docker 
#docker run -i -t r4:v.0.01 #run the docker
#sudo ls #check that sudo is present and work
#Then go in the home dir
#Rscript ./umap_command.R #run the R file to install umap exit
#docker ps -a #find the instance name ID of the instance previously closed 
#docker commit instaceID (find to be in my case d7a30dac5e11 through docker ps -a) r4:v.0.01 #commit the changes to the docker
#docker run -i -t r4:v.0.01 #run again the docker
#sudo ls #check that the changes were saved
#exit

#LESSON4 CELL RANGER
#docker run r4:v.0.01
#apt-get install -y wget
#in cd /bin paste the wget from cell ranger website 
#gzip -d  cellranger-7.0.0.tar.gz to unzip
#tar xvf cellranger-7.0.0.tar to separate the folders 
#exit and docker commit istance r4:v.0.01
#nano ~/.bashrc enter the editor and add export PATH=/bin/cellranger-7.0.0:$PATH
#to exit editor control+x and then Y and enter
#exit  and commit 
###RUNNING FROM THE INSIDE OF THE DOCKER 
#docker run -i -t -v path/of/the/folder:/home r4:v.0.0.1
#in home create the folder with command mkdir "nomedelladir"
#exit and commit 
###Run cellranger with a shell command to execute mat2csv interactively
#




#PCA function library(docker4seq) pca(
#experiment.table = "./_counts.txt",
#type = c("counts", "FPKM", "TPM"),
#covariatesInNames = FALSE,
#samplesName = TRUE,
#principal.components = c(1, 2),
#legend.position = c("bottom", "bottomleft", "left", "topleft", "top", "topright",
#                    "right", "center"), pdf = TRUE, output.folder = getwd()
#)

install.packages("SAVER")
library(SAVER)

counts2cpm <- function(file, sep = ","){
  tmp <- read.table(file, sep=sep, header=T, row.names=1) 
  col.sum <- apply(tmp, 2, sum)
  tmp1 <- t(tmp)/col.sum
  tmp1 <- t(tmp1)
  tmp1 <- tmp1 * 1000000 
  write.table(tmp1, "cpm.csv", sep=",", col.names=NA) 
  write.table(log2(tmp1 + 1), "log2cpm.csv",sep=",", col.names=NA)
}

getwd()
setwd("/Users/angeloscarciglia/Desktop/esercizio 6")

counts2cpm(file="saver_ctrl.csv",sep =",")
file.rename(from="log2cpm.csv", to ="ctrl_log2cpm.csv")

counts2cpm(file="saver_osi.csv", sep=",")
file.rename(from="log2cpm.csv", to= "osi_log2cpm.csv")

install.packages("umap")
install.packages("ggplot2")
library(ggplot2)
library(umap)

#CTRL samples
ctrl <- read.table("ctrl_log2cpm.csv", sep=",", header=T, row.names=1) 
ctrl.umap <- umap(t(ctrl), random_state=111, n_epochs = 1000) #here you can insert variable min_dist=0.05, n_neighbors= 10 to play with the graph
f=data.frame(x=as.numeric(ctrl.umap$layout[,1]),y=as.numeric(ctrl.umap$layout[,2])) 

#additional info 
#for example select the cells that have value higher than 20 in the x axis 
#length(which(f$x > 20)) # 10
#select the cells that have a value of y axis < -30 and the x axis < -50; it select the rettangle 
#length(intersect(which(f$y < -30), which(f$x < -50)) # 4

#plotting UMAP
sp <- ggplot(f, aes(x=x,y=y)) + geom_point(pch=19, cex=0.3) 
pdf("ctrlUMAPprova2.pdf")
print(sp)
dev.off()

#OSI samples
osi <- read.table("osi_log2cpm.csv", sep=",", header=T, row.names=1) 
osi.umap <- umap(t(osi), random_state=111, n_epochs = 1000) 
f=data.frame(x=as.numeric(osi.umap$layout[,1]),y=as.numeric(osi.umap$layout[,2]))

#plotting UMAP
sp <- ggplot(f, aes(x=x,y=y)) + geom_point(pch=19, cex=0.3)
pdf("osiUMAPprova2.pdf") 
print(sp) 
dev.off()

install.packages("Rtsne")
library(Rtsne)
library(ggplot2)

#CTRL samples
tmp <- read.table("ctrl_log2cpm.csv", sep=",", header=T, row.names=1) 
tmp.labels <- sapply(strsplit(names(tmp), '\\.'), function(x)x[2]) 
cell_line <- as.factor(tmp.labels)
set.seed(111)
tsne.out <- Rtsne(as.matrix(t(tmp)), pca= FALSE, perplexity = 30, theta=0.0) 
f=data.frame (x= as.numeric(tsne.out$Y[,1]),y=tsne.out$Y[,2])
#plotting the tsne
sp <- ggplot(f, aes(x=x,y=y)) + geom_point(pch=19, cex=0.3) 
pdf("ctrl_tsne_noPCAprova3.pdf")
print(sp)
dev.off()
#OSI samples
tmp <- read.table("osi_log2cpm.csv", sep=",", header=T, row.names=1) 
tsne.out <- Rtsne(as.matrix(t(tmp)), pca= FALSE, perplexity = 30, theta=0.0) 
f=data.frame (x= as.numeric(tsne.out$Y[,1]),y=tsne.out$Y[,2])
#plotting the tsne
sp <- ggplot(f, aes(x=x,y=y)) + geom_point(pch=19, cex=0.3) 
pdf("osi_tsne_noPCA.pdf")
print(sp)
dev.off()

if (!require("BiocManager", quietly =TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
browseVignettes("DESeq2")

#Annotation problem: for scRNAseq data the annotation of the genes 
#is ensembleID:symbol while for bulkRNAseq data the annotation is 
#symbol:ensembleID. To solve this problem we have to use the strsplit 
#function to have the same annotation to perform intersection.

getwd()
#start with the bulkRNAseq file p9_wo_w_osi.txt and setwd(/wherethefileis) 
#insert the file as a variable in R
x <- read.table("pc9_wo_w_osi.txt", header=TRUE, sep = "\t", row.names=1)
# with header=TRUE it considers the first line as a header (colum names)
# with row.names=1 it means that the 1st column of the table becomes the rowname of the dataframe #these 2 options are conceptually the same but one for columns and the other for rows
#to check the file to see if it's correct
head(x)
#use strsplit function to create a large list containing the symbol and the ensemble IDs separeted 
y <- strsplit(row.names(x), ":")
#the strsplit f unction take the rownames of the dataframe and separate the name for the ":" character 
#convert the list to a matrix
matrix <- matrix(unlist(y), ncol = 2, byrow = TRUE)
# at this point you have a list that as to be converted in a matrix with the unlist function that use the number of colum=2 and unlist b y row.
#use the paste function to create a vector with the new annotation as in the single cell dataset ensembleID:symbol
z <- paste(matrix[,2], matrix[,1], sep=":")
#the function paste, paste the two colum of the matrix in the order that we want separeted by the character ":" 
#put the new row names in the original dataframe (x) where we had all the data
row.names(x) <- z

#Since you have to perform DE on Acute treatment vs Control and on Chronic treatment 
#vs Control, it is necessary to split the table given in two different matrices. 
#One matrix contains the three ctrl replicates plus three acute treatment replicates 
#and the second matrix contains the three crtl replicates plus three chronic 
#treatment replicates

#to make the two matrix
acute<-x[,c(1,2,3,7,8,9)]
chronic<-x[,c(1,2,3,4,5,6)]
#create the two matrices starting on the dataframe (x) and selecting the coloums of interest and put this in two new variable one contaning
#save the file as txt to be used in future
write.table(acute, file="acute_dataframe.txt", sep="\t", row.names = TRUE, col.names = TRUE) 
write.table(chronic, file="chronic_dataframe.txt", sep="\t", row.names = TRUE, col.names = TRUE)

#Since made from dataframe they are dataframe and you convert them into matrices
#read the files
acute_df<-read.table("acute_dataframe.txt",header = TRUE, sep="\t", row.names = 1) 
chronic_df<-read.table("chronic_dataframe.txt",header = TRUE, sep="\t", row.names = 1)
#convert the dataframe in a matrix
acute_matrix<-as.matrix(acute_df, rownames.force=NA) 
chronic_matrix<-as.matrix(chronic_df, rownames.force=NA)

#create the two vectors containing the samples and the conditions 
samples_acute <- colnames(acute_matrix)
samples_chronic <- colnames(chronic_matrix)
#take the colnames of the two matrix and put it in a vector 
conditions<-c(rep("ctrl",3), rep("treated",3))
#create a new vector that contain three rows contaning ctrl and then three row with treated
#create the dataframe combining the two vectors. The names of the vectors will be considered automatically as header 
coldata_acute <- data.frame(samples_acute, conditions)
coldata_chronic <- data.frame(samples_chronic, conditions)
#put the rownames in the new dataframe taking them from the vectors we created before 
rownames(coldata_acute)<-samples_acute
rownames(coldata_chronic)<-samples_chronic
#factor the conditions column 
coldata_acute$conditions<-factor(coldata_acute$conditions) 
coldata_chronic$conditions<-factor(coldata_chronic$conditions)
#control that effectively the conditions column it is a factor 
class(coldata_acute$conditions) 
class(coldata_chronic$conditions)

# now we are ready to do DESeq on acute since deSEQ wants 3 parameters:
# 1. countData (the matrix with the data)
# 2. colData (our dataframe containing the samples and the conditions)
# 3. design = the parameter for which we have to perform differential expression. Here we put our condition vector with the 2 conditions to
library("DESeq2")# to load the library with all the function needed
#DESeqDatasetFromMatrix generate a Dataset on which the DESeq can be run.
dds_acute <- DESeqDataSetFromMatrix(countData = acute_matrix, colData = coldata_acute, design=~conditions) 
dds_chronic <- DESeqDataSetFromMatrix(countData = chronic_matrix, colData = coldata_chronic, design=~conditions)
# Use the DESeq function to perform the DE analysis 
dds_acute <- DESeq(dds_acute)
dds_chronic <- DESeq(dds_chronic)
# Save the results in a variable 
res_acute <- results(dds_acute) 
res_chronic <- results(dds_chronic)
# Filter the resulted data of the DE analysis for adjusted pValue <=0.05 e |log2FoldChange|>=1. 
res_acute_filter<-res_acute[which(res_acute[,6]<=0.05 & abs(res_acute[,2])>=1),] 
res_chronic_filter<-res_chronic[which(res_chronic[,6]<=0.05 & abs(res_chronic[,2])>=1),]
# This command creates a new variable taking all the rows from the res_* dataframe that have a value in the $padj column <=0.05 and the abs # All the columns of the res_* dataframe are selected and maintained.
# You can see that which is followed by [..., ] the "..." are the rows, the space means all columns.
#Save the filtered data into a file
write.table(res_acute_filter, file="DESeq_results_acute_filter.txt", sep="\t", row.names=TRUE,col.names=TRUE)
write.table(res_chronic_filter, file="DESeq_results_chronic_filter.txt", sep="\t", row.names=TRUE,col.names=TRUE)

#Filter osi_log2cpm.csv using DE of acute RNAseq data and plot UMAP
#to put the single cell data of the osimertinib treatment in a variable using read.table 
osi <- read.table("osi_log2CPM.csv", sep=",", header=T, row.names=1)
#perform the filtering of the single cell data using the DE genes of RNAseq acute treatment 
intersection_osi_acute <- intersect(row.names(res_acute_filter), row.names(osi))
#select the row of the data frame (osi) of the single cell experiment containing the DE genes in common between the single cell data and th 
osiUMAP_acute<-osi[which(row.names(osi)%in%intersection_osi_acute), ]
#save the data in a csv file needed to run the UMAP function 
write.csv(osiUMAP_acute, "osiUMAP_acutelog2cpm.csv", sep=",")
#load the library necessary to run UMAP 
library(umap)
library(ggplot2)
#run UMAP function
osiUMAP_acute <- read.table("osiUMAP_acutelog2cpm.csv", sep=",", header=T, row.names=1) 
osi.umap <- umap(t(osiUMAP_acute), random_state=111, n_epochs = 1000) 
f=data.frame(x=as.numeric(osi.umap$layout[,1]),y=as.numeric(osi.umap$layout[,2]))
#plotting UMAP
sp <- ggplot(f, aes(x=x,y=y)) + geom_point(pch=19, cex=0.3) 
pdf("osiUMAP_acute.pdf")
print(sp)
dev.off()
#Filter osi-log2cpm.csv using DE of chronic RNAseq data and plot UMAP
#perform the filtering of the single cell data using the DE genes of RNAseq chronic treatment 
intersection_osi_chronic <- intersect(row.names(res_chronic_filter), row.names(osi))
#select the row of the data frame (osi) of the single cell experiment containing the DE genes in common between the single cell data and th 
osiUMAP_chronic<-osi[which(row.names(osi)%in%intersection_osi_chronic), ]
#save the data in a csv file needed to run the UMAP function 
write.csv(osiUMAP_chronic,"osiUMAP_chroniclog2cpm.csv", sep=",")
#run UMAP function
osiUMAP_chronic <- read.table("osiUMAP_chroniclog2cpm.csv", sep=",", header=T, row.names=1) 
osi.umap <- umap(t(osiUMAP_chronic), random_state=111, n_epochs = 1000) 
f=data.frame(x=as.numeric(osi.umap$layout[,1]),y=as.numeric(osi.umap$layout[,2]))
#plotting UMAP
sp <- ggplot(f, aes(x=x,y=y)) + geom_point(pch=19, cex=0.3) 
pdf("osiUMAP_chronic.pdf")
print(sp)
dev.off()
?write.csv
#Filter osi-log2cpm.csv using combined DE of acute and chronic RNAseq data and plot UMAP
#use union function to combine and take all the data
#perform first the union on the RNAseq data combining acute and chronic DE genes to have all. 
union_acute_chronic <- union(row.names(res_chronic_filter),row.names(res_acute_filter))
#filter the scRNAseq osi dataset with the genes from chronic and acute united 
osiUMAP_union<-osi[which(row.names(osi)%in%union_acute_chronic), ]
#save the dataset into a csv table
write.csv(osiUMAP_union, "osiUMAP_unionlog2cpm.csv", sep=",")
#run UMAP function
osi_union <- read.table("osiUMAP_unionlog2cpm.csv", sep=",", header=T, row.names=1) 
osi.umap <- umap(t(osi_union), random_state=111, n_epochs = 1000) 
f=data.frame(x=as.numeric(osi.umap$layout[,1]),y=as.numeric(osi.umap$layout[,2]))
#plotting UMAP
sp <- ggplot(f, aes(x=x,y=y)) + geom_point(pch=19, cex=0.3) 
pdf("osiUMAP_union.pdf")
print(sp)
dev.off()




#generate the log2cpm file formt he saver 
counts2cpm <- function(file, sep =","){
  tmp <- read.table(file,sep=sep, header=T, row.names=1)
  col.sum <- apply(tmp,2,sum)
  tmp1 <- t(tmp)/col.sum
  tmp1 <- t(tmp1)
  tmp1 <- tmp1 * 1000000
  write.table(tmp1, "cpm.csv", sep=",", col.names=NA)
  write.table(log2(tmp1+1), "log2cpm.csv", sep=",", col.names=NA)
}

getwd()
setwd("/Users/angeloscarciglia/Desktop")
counts2cpm(file="saver_RNA-5c.csv", sep=",")
file.rename(from="log2cpm.csv", to="RNA-5c_log2cpm.csv")
library(Rtsne)
library(ggplot2)
#to run tsne 
tmp <- read.table("RNA-5c_log2cpm.csv", sep=",", header=T, row.names =1)
tmp.labels <- sapply(strsplit(names(tmp), '\\.'), function(x)x[2])
cell_line <- as.factor(tmp.labels)
set.seed(111)

#we don't want any pca before tsne
tsne_out <- Rtsne(as.matrix(t(tmp)), pca=FALSE, perplexity=30, theta=0.0)
t=data.frame(x=as.numeric(tsne_out$Y[,1]), y=tsne_out$Y[,2])
#plotting tsne
sp <- ggplot (t, aes (x=x, y=y,)) + geom_point(pch=19, cex=0.3)
pdf("tSNE_RNA-5c.pdf")
print(sp)
dev.off()
#running umap                     
library(umap)
library(ggplot2)
umap<- read.table("RNA-5c_log2cpm.csv", sep=",", header=T, row.names=1)
umap.labels <- sapply(strsplit(names(ctrl), '\\.'), function(x)x[2])
cell_line <- as.factor(umap.labels)
sc.umap <- umap (t(ctrl), random_state=111, n_epochs=1000)
u=data.frame (x=as.numeric(sc.umap$layout[,1]), y=as.numeric(sc.umap$layout[,2]))
sp <- ggplot(u,aes(x=x,y=y))+geom_point(pch=19, cex=0.3, color='darkblue')
pdf("UMAP_RNA-5c.pdf")
print(sp)
dev.off()

#to add colours
log2cpm_RNA<- read.table("RNA-5c_log2cpm.csv", sep=",", header=TRUE, row.names=1)
tsne_out <- Rtsne(as.matrix(t(log2cpm_RNA)), pca=FALSE, perplexity=30, theta=0.0)
t=data.frame(x=as.numeric(tsne_out$Y[,1]), y=tsne_out$Y[,2])

