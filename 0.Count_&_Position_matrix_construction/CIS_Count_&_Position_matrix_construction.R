#Clean environment
rm(list=ls(all.names=TRUE))
gc()

#Set the working directory
setwd("/Users/srujansingh/Documents/JeanFan_Lab/Rabb_cold_ischemia")

#load libraries
library(Matrix)

###################### 0 hr COLD ISCHEMIC KIDNEY ################################
#Read in gene expression to create the counts matrix
gexp<-readMM('1_0hr/filtered_feature_bc_matrix/matrix.mtx.gz')
head(gexp)

barcodes<-read.csv('1_0hr/filtered_feature_bc_matrix/barcodes.tsv.gz', header=FALSE, sep="\t")
head(barcodes)

features<-read.csv('1_0hr/filtered_feature_bc_matrix/features.tsv.gz', header=FALSE, sep="\t")
head(features)

rownames(gexp)<-features$V2
colnames(gexp)<-barcodes$V1

#Read in positions of the spots to create the position matrix
pos<-read.csv("1_0hr/spatial/tissue_positions.csv")
pos<-pos[which(pos$in_tissue==1),]
rownames(pos)<-pos$barcode
pos<-pos[colnames(gexp),]

pos<-pos[,-c(1:4)]
colnames(pos)<-c('X','Y')
par(mfrow=c(1,1))
plot(pos$X,pos$Y, main ='0 hr cold ischemia')

#Adding prefix to spot names to distinguish from other time points
rownames(pos)<-paste0("CIS_0h-", rownames(pos))
colnames(gexp)<-paste0("CIS_0h-", colnames(gexp))

#to check if the barcodes of the spots in the count and position matrix are consistent
all(rownames(pos)%in%colnames(gexp))

#Creating a list containing position and gene expression data for 0 hr cold ischemic kidney
CIS_0h=list(pos=pos, gexp=gexp)

###################### 12hr COLD ISCHEMIC KIDNEY ###############################
#Read in gene expression to create the counts matrix
gexp<-readMM('2_12hr/filtered_feature_bc_matrix/matrix.mtx.gz')
head(gexp)

barcodes<-read.csv('2_12hr/filtered_feature_bc_matrix/barcodes.tsv.gz', header=FALSE, sep="\t")
head(barcodes)

features<-read.csv('2_12hr/filtered_feature_bc_matrix/features.tsv.gz', header=FALSE, sep="\t")
head(features)

rownames(gexp)<-features$V2
colnames(gexp)<-barcodes$V1

#Read in positions of the spots to create the position matrix
pos<-read.csv("2_12hr/spatial/tissue_positions.csv")
pos<-pos[which(pos$in_tissue==1),]
rownames(pos)<-pos$barcode
pos<-pos[colnames(gexp),]

pos<-pos[,-c(1:4)]
colnames(pos)<-c('X','Y')
par(mfrow=c(1,1))
plot(pos$X,pos$Y, main ='12 hrs cold ischemia')

#Adding prefix to spot names to distinguish from other time points
rownames(pos)<-paste0("CIS_12h-", rownames(pos))
colnames(gexp)<-paste0("CIS_12h-", colnames(gexp))

#to check if the barcodes of the spots in the count and position matrix are consistent
all(rownames(pos)%in%colnames(gexp))

#Creating a list containing position and gene expression data for 12 hr cold ischemic kidney
CIS_12h=list(pos=pos, gexp=gexp)


###################### 24hrs COLD ISCHEMIC KIDNEY ###############################
#Read in gene expression to create the counts matrix
gexp<-readMM('3_24hr/filtered_feature_bc_matrix/matrix.mtx.gz')
head(gexp)

barcodes<-read.csv('3_24hr/filtered_feature_bc_matrix/barcodes.tsv.gz', header=FALSE, sep="\t")
head(barcodes)

features<-read.csv('3_24hr/filtered_feature_bc_matrix/features.tsv.gz', header=FALSE, sep="\t")
head(features)

rownames(gexp)<-features$V2
colnames(gexp)<-barcodes$V1

#Read in positions of the spots to create the position matrix
pos<-read.csv("3_24hr/spatial/tissue_positions.csv")
pos<-pos[which(pos$in_tissue==1),]
rownames(pos)<-pos$barcode
pos<-pos[colnames(gexp),]

pos<-pos[,-c(1:4)]
colnames(pos)<-c('X','Y')
par(mfrow=c(1,1))
plot(pos$X,pos$Y, main ='24 hrs cold ischemia')

#Adding prefix to spot names to distinguish from other time points
rownames(pos)<-paste0("CIS_24h-", rownames(pos))
colnames(gexp)<-paste0("CIS_24h-", colnames(gexp))

#to check if the barcodes of the spots in the count and position matrix are consistent
all(rownames(pos)%in%colnames(gexp))

#Creating a list containing position and gene expression data for 12 hr cold ischemic kidney
CIS_24h=list(pos=pos, gexp=gexp)


###################### 48hrs COLD ISCHEMIC KIDNEY ###############################
#Read in gene expression to create the counts matrix
gexp<-readMM('4_48hr/filtered_feature_bc_matrix/matrix.mtx.gz')
head(gexp)

barcodes<-read.csv('4_48hr/filtered_feature_bc_matrix/barcodes.tsv.gz', header=FALSE, sep="\t")
head(barcodes)

features<-read.csv('4_48hr/filtered_feature_bc_matrix/features.tsv.gz', header=FALSE, sep="\t")
head(features)

rownames(gexp)<-features$V2
colnames(gexp)<-barcodes$V1

#Read in positions of the spots to create the position matrix.
pos<-read.csv("4_48hr/spatial/tissue_positions.csv")
pos<-pos[which(pos$in_tissue==1),]
rownames(pos)<-pos$barcode
pos<-pos[colnames(gexp),]

pos<-pos[,-c(1:4)]
colnames(pos)<-c('X','Y')
par(mfrow=c(1,1))
plot(pos$X,pos$Y, main ='48 hrs cold ischemia')

#Adding prefix to spot names to distinguish from other time points
rownames(pos)<-paste0("CIS_48h-", rownames(pos))
colnames(gexp)<-paste0("CIS_48h-", colnames(gexp))

#to check if the barcodes of the spots in the count and position matrix are consistent
all(rownames(pos)%in%colnames(gexp))

#Creating a list containing position and gene expression data for 12 hr cold ischemic kidney
CIS_48h=list(pos=pos, gexp=gexp)

#Saving the data (as a list)
save(CIS_0h, CIS_12h, CIS_24h, CIS_48h, file="CIS_data.Rdata")
