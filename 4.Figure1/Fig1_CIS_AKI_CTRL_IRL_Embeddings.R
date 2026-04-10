#Clean environment
rm(list=ls(all.names=TRUE))
gc()

#Set your working directory accordingly 
setwd("/Users/srujansingh/Library/CloudStorage/OneDrive-SharedLibraries-JohnsHopkins/Jean Fan - 2024_spatial_AKI/Manuscript")

## load datasets (located in data folder.)
load("data/AKI_data.RData")
load("data/CIS_data.RData")
load("data/Rabb_ctrl.RData")
load("data/Rabb_irl24h.RData")

#Loading the required library
library(STdeconvolve) #For finding out overdispersed genes
library(harmony) #For implementing harmony to correct for batch effects
library(MUDAN) 
library(openxlsx)
library(RColorBrewer)
library(ggplot2)

################################################################################
################# COMBINED CIS,AKI, CTRL AND IRL DATASET SEGMENTATION ##########

##Removing all mitochondrial genes from the gene experssion table (AKI dataset)
AKI_sham$gexp<-AKI_sham$gexp[-grep("^mt", rownames(AKI_sham$gexp)),]
AKI_4h$gexp<-AKI_4h$gexp[-grep("^mt", rownames(AKI_4h$gexp)),]
AKI_12h$gexp<-AKI_12h$gexp[-grep("^mt", rownames(AKI_12h$gexp)),]
AKI_2d$gexp<-AKI_2d$gexp[-grep("^mt", rownames(AKI_2d$gexp)),]
AKI_6w$gexp<-AKI_6w$gexp[-grep("^mt", rownames(AKI_6w$gexp)),]

##Removing all mitochondrial genes from the gene experssion table (IRL dataset)
# Note IRI dataset is the same as AKI24 dataset in the main article section.
irl1$gexp<-irl1$gexp[-grep("^mt", rownames(irl1$gexp)),]
irl2$gexp<-irl2$gexp[-grep("^mt", rownames(irl2$gexp)),]
irl3$gexp<-irl3$gexp[-grep("^mt", rownames(irl3$gexp)),]
irl4$gexp<-irl4$gexp[-grep("^mt", rownames(irl4$gexp)),]

##Removing all mitochondrial genes from the gene experssion table (ctrl dataset)
ctrl1$gexp<-ctrl1$gexp[-grep("^mt", rownames(ctrl1$gexp)),]
ctrl2$gexp<-ctrl2$gexp[-grep("^mt", rownames(ctrl2$gexp)),]
ctrl3$gexp<-ctrl3$gexp[-grep("^mt", rownames(ctrl3$gexp)),]
ctrl4$gexp<-ctrl4$gexp[-grep("^mt", rownames(ctrl4$gexp)),]

#CIS dataset doesn't have any mitochondrial genes


COD<-as.matrix(read.xlsx("Supplementary_Tables/Fig1_CIS_AKI_CTRL_IRL_Top_Feature_Selected_Genes.xlsx"))

## Impelementing HARMONY
#Getting the data
dataA<-CIS_0h$gexp[COD,]
dataB<-CIS_12h$gexp[COD,]
dataC<-CIS_24h$gexp[COD,]
dataD<-CIS_48h$gexp[COD,]

dataE<-AKI_sham$gexp[COD,]
dataF<-AKI_4h$gexp[COD,]
dataG<-AKI_12h$gexp[COD,]
dataH<-AKI_2d$gexp[COD,]
dataI<-AKI_6w$gexp[COD,]

dataJ<-irl1$gexp[COD,]
dataK<-irl2$gexp[COD,]
dataL<-irl3$gexp[COD,]
dataM<-irl4$gexp[COD,]

dataN<-ctrl1$gexp[COD,]
dataO<-ctrl2$gexp[COD,]
dataP<-ctrl3$gexp[COD,]
dataQ<-ctrl4$gexp[COD,]

#Combine different datasets into a single dataset
cd <- cbind(dataA, dataB, dataC, dataD, 
            dataE, dataF, dataG, dataH, dataI,
            dataJ, dataK, dataL, dataM,
            dataN, dataO, dataP, dataQ)
cd<-as.matrix(cd)
#cd<-cd[,-c(19931,20077)]
cd<-cd[,-c(19931,20077)] #These two spots have no gene expression in them.

# meta data 
meta <- c( rep('CIS_0h', ncol(dataA)), rep('CIS_12h', ncol(dataB)), 
           rep('CIS_24h', ncol(dataC)), rep('CIS_48h', ncol(dataD)), 
           rep('AKI_sham', ncol(dataE)), rep('AKI_4h', ncol(dataF)), rep('AKI_12h', ncol(dataG)),
           rep('AKI_48h', ncol(dataH)), rep('AKI_6w', ncol(dataI)),
           
           rep('IRL_1', ncol(dataJ)), rep('IRL_2', ncol(dataK)), 
           rep('IRL_3', ncol(dataL)), rep('IRL_4', ncol(dataM)),
           rep('CTRL_1', ncol(dataN)), rep('CTRL_2', ncol(dataO)), 
           rep('CTRL_3', ncol(dataP)), rep('CTRL_4', ncol(dataQ)))

names(meta) <- c(colnames(dataA), colnames(dataB), colnames(dataC), colnames(dataD),
                 colnames(dataE), colnames(dataF), colnames(dataG), colnames(dataH), colnames(dataI),
                 colnames(dataJ), colnames(dataK), colnames(dataL), colnames(dataM),
                 colnames(dataN), colnames(dataO), colnames(dataP), colnames(dataQ))
meta<-meta[-c(19931,20077)]

meta <- factor(meta)

#Non-Harmonized Embeddings
emb<-data.frame(readRDS("data/AKI-CIS-irl-ctrl_emb.rds"))
colnames(emb)<-c('tsne1','tsne2')
emb$label<-meta[rownames(emb)]
emb$label<-factor(emb$label, levels=c('CIS_0h','CIS_12h','CIS_24h','CIS_48h',
                                      'AKI_sham','AKI_4h','AKI_12h','AKI_48h','AKI_6w',
                                      'CTRL_1','CTRL_2', 'CTRL_3', 'CTRL_4',
                                      'IRL_1', 'IRL_2', 'IRL_3','IRL_4'))

color_palette<-colorRampPalette(brewer.pal(n=9, name="Set1"))(17)

pdf('Figures/Figure1/pdfs/Non-Harmonized_AKI-CIS-IRL-CTRL_datasets_ggplot.pdf', height=5, width=7)
ggplot(emb, aes(x=tsne1,y=tsne2, col=label)) + geom_point(size=0.3) +
  scale_color_manual(values=color_palette) +
  theme_minimal() +
  theme(panel.grid = element_blank(),          # Remove gridlines
        axis.title = element_blank(),          # Remove axis labels
        axis.text = element_blank(),           # Remove axis text
        axis.ticks = element_blank(),          # Remove axis ticks
        panel.border = element_rect(color="black", fill=NA, linewidth =0.5),  # Add border around plot
        legend.position = "left" 
  ) + guides(color = guide_legend(override.aes = list(size = 2)))
dev.off()

#Harmonized Embeddings
Hemb<-data.frame(readRDS("data/AKI-CIS-irl-ctrl_Harmonized_tSNE_Embeds.rds"))
colnames(Hemb)<-c('tsne1','tsne2')
Hemb$label<-meta[rownames(Hemb)]
Hemb$label<-factor(Hemb$label, levels=c('CIS_0h','CIS_12h','CIS_24h','CIS_48h',
                                      'AKI_sham','AKI_4h','AKI_12h','AKI_48h','AKI_6w',
                                      'CTRL_1','CTRL_2', 'CTRL_3', 'CTRL_4',
                                      'IRL_1', 'IRL_2', 'IRL_3','IRL_4'))

color_palette<-colorRampPalette(brewer.pal(n=9, name="Set1"))(17)

pdf('Figures/Figure1/pdfs/Harmonized_AKI-CIS-IRL-CTRL_datasets_ggplot.pdf', height=5, width=7)
ggplot(Hemb, aes(x=tsne1,y=tsne2, col=label)) + geom_point(size=0.3) +
  scale_color_manual(values=color_palette) +
  theme_minimal() +
  theme(panel.grid = element_blank(),          # Remove gridlines
        axis.title = element_blank(),          # Remove axis labels
        axis.text = element_blank(),           # Remove axis text
        axis.ticks = element_blank(),          # Remove axis ticks
        panel.border = element_rect(color="black", fill=NA, linewidth =0.5),  # Add border around plot
        legend.position = "left" 
  ) + guides(color = guide_legend(override.aes = list(size = 2)))
dev.off()

#Cluster Annotations
#Extracting the compartment specific spots for CIS and AKI:

Clust<-data.frame(readRDS("data/AKI-CIS-irl-ctrl_Harmonized_AllClusters.rds"))
head(Hemb)

Hemb$tissue<-'Other'
Hemb$cluster<-Clust[rownames(Hemb),]
Hemb[which(Hemb$cluster==1),]$tissue<-'I.Medulla'
Hemb[which(Hemb$cluster==7|Hemb$cluster==8),]$tissue<-'O.Medulla'
Hemb[which(Hemb$cluster==2|Hemb$cluster==3|Hemb$cluster==4|
            Hemb$cluster==6|Hemb$cluster==10),]$tissue<-'Cortex'

Hemb$tissue<-factor(Hemb$tissue, levels=c('Cortex','O.Medulla','I.Medulla','Other'))

color_palette<-c("Cortex"="red", "O.Medulla"="blue","I.Medulla"="green","Other"="gray")

pdf('Figures/Figure1/pdfs/Harmonized_AKI-CIS-IRL-CTRL_Clusters_ggplot.pdf', height=5, width=7)
ggplot(Hemb, aes(x=tsne1,y=tsne2, col=tissue)) + geom_point(size=0.3) +
  scale_color_manual(values=color_palette) +
  theme_minimal() +
  theme(panel.grid = element_blank(),          # Remove gridlines
        axis.title = element_blank(),          # Remove axis labels
        axis.text = element_blank(),           # Remove axis text
        axis.ticks = element_blank(),          # Remove axis ticks
        panel.border = element_rect(color="black", fill=NA, linewidth =0.5), # Add border around plot
        legend.position = "left" 
  ) + guides(color = guide_legend(override.aes = list(size = 2)))
dev.off()

## This section is used to generate images in section E of Figure1.

my_theme <- theme_bw() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.text = element_text(),
        legend.position="none"
  )

#CIS Tissue (Cold Ischemia Injury):
CIS_0h.df<-data.frame(CIS_0h$pos)
CIS_0h.df$Compartment<-'Other'
CIS_0h.df$Compartment<-Hemb[rownames(CIS_0h.df),]$tissue

CIS_12h.df<-data.frame(CIS_12h$pos)
CIS_12h.df$Compartment<-'Other'
CIS_12h.df$Compartment<-Hemb[rownames(CIS_12h.df),]$tissue

CIS_24h.df<-data.frame(CIS_24h$pos)
CIS_24h.df$Compartment<-'Other'
CIS_24h.df$Compartment<-Hemb[rownames(CIS_24h.df),]$tissue

CIS_48h.df<-data.frame(CIS_48h$pos)
CIS_48h.df$Compartment<-'Other'
CIS_48h.df$Compartment<-Hemb[rownames(CIS_48h.df),]$tissue

color_palette<-c("Cortex"="red", "O.Medulla"="blue","I.Medulla"="green","Other"="gray")

g1<-ggplot(CIS_0h.df, aes(x=X, y=Y, col=Compartment)) + geom_point(size=0.5) + 
  ggtitle("CIS_0h") + scale_color_manual(values=color_palette) + my_theme

g2<-ggplot(CIS_12h.df, aes(x=X, y=Y, col=Compartment)) + geom_point(size=0.5) + 
  ggtitle("CIS_12h") + scale_color_manual(values=color_palette) + my_theme

g3<-ggplot(CIS_24h.df, aes(x=X, y=Y, col=Compartment)) + geom_point(size=0.5) + 
  ggtitle("CIS_24h") + scale_color_manual(values=color_palette) + my_theme

g4<-ggplot(CIS_48h.df, aes(x=X, y=Y, col=Compartment)) + geom_point(size=0.5) + 
  ggtitle("CIS_48h") + scale_color_manual(values=color_palette) + my_theme

library(gridExtra)
pdf("Figures/Figure1/pdfs/CIS_Compartment_Annotations_ggplot.pdf", height=3, width=12)
grid.arrange(g1,g2, g3, g4, ncol=4)
dev.off()


#AKI Tissue (Warm Ischmia-Reperfusion Injury):
AKI_sham.df<-data.frame(AKI_sham$pos)
AKI_sham.df$Compartment<-'Other'
AKI_sham.df$Compartment<-Hemb[rownames(AKI_sham.df),]$tissue
colnames(AKI_sham.df)<-c('X','Y','Compartment')

AKI_4h.df<-data.frame(AKI_4h$pos)
AKI_4h.df$Compartment<-'Other'
AKI_4h.df$Compartment<-Hemb[rownames(AKI_4h.df),]$tissue
colnames(AKI_4h.df)<-c('X','Y','Compartment')

AKI_12h.df<-data.frame(AKI_12h$pos)
AKI_12h.df$Compartment<-'Other'
AKI_12h.df$Compartment<-Hemb[rownames(AKI_12h.df),]$tissue
colnames(AKI_12h.df)<-c('X','Y','Compartment')

AKI_48h.df<-data.frame(AKI_2d$pos)
AKI_48h.df$Compartment<-'Other'
AKI_48h.df$Compartment<-Hemb[rownames(AKI_48h.df),]$tissue
colnames(AKI_48h.df)<-c('X','Y','Compartment')

color_palette<-c("Cortex"="red", "O.Medulla"="blue","I.Medulla"="green","Other"="gray")

g1<-ggplot(AKI_sham.df, aes(x=X, y=Y, col=Compartment)) + geom_point(size=0.7) + 
  ggtitle("AKI_sham") + scale_color_manual(values=color_palette) + my_theme

g2<-ggplot(AKI_4h.df, aes(x=X, y=Y, col=Compartment)) + geom_point(size=0.7) + 
  ggtitle("AKI_4h") + scale_color_manual(values=color_palette) + my_theme

g3<-ggplot(AKI_12h.df, aes(x=X, y=Y, col=Compartment)) + geom_point(size=0.7) + 
  ggtitle("AKI_12h") + scale_color_manual(values=color_palette) + my_theme

g4<-ggplot(AKI_48h.df, aes(x=X, y=Y, col=Compartment)) + geom_point(size=0.7) + 
  ggtitle("AKI_48h") + scale_color_manual(values=color_palette) + my_theme

library(gridExtra)
pdf("Figures/Figure1/pdfs/AKI_Compartment_Annotations_ggplot.pdf", height=3, width=12)
grid.arrange(g1,g2, g3, g4, ncol=4)
dev.off()


#Similar to CIS and AKI tissue, the different compartments can be visulaized in the
#CTRL and the IRL tissue as well (Section M and N of Supplementary Figure1)

#CTRL Tissue (Native kidney tissues):
CTRL1.df<-data.frame(ctrl1$pos)
CTRL1.df$Compartment<-'Other'
CTRL1.df$Compartment<-Hemb[rownames(CTRL1.df),]$tissue
colnames(CTRL1.df)<-c('X','Y','Compartment')

CTRL2.df<-data.frame(ctrl2$pos)
CTRL2.df$Compartment<-'Other'
CTRL2.df$Compartment<-Hemb[rownames(CTRL2.df),]$tissue
colnames(CTRL2.df)<-c('X','Y','Compartment')

CTRL3.df<-data.frame(ctrl3$pos)
CTRL3.df$Compartment<-'Other'
CTRL3.df$Compartment<-Hemb[rownames(CTRL3.df),]$tissue
colnames(CTRL3.df)<-c('X','Y','Compartment')

CTRL4.df<-data.frame(ctrl4$pos)
CTRL4.df$Compartment<-'Other'
CTRL4.df$Compartment<-Hemb[rownames(CTRL4.df),]$tissue
colnames(CTRL4.df)<-c('X','Y','Compartment')

color_palette<-c("Cortex"="red", "O.Medulla"="blue","I.Medulla"="green","Other"="gray")

g1<-ggplot(CTRL1.df, aes(x=X, y=Y, col=Compartment)) + geom_point(size=0.5) + 
  ggtitle("CTRL1") + scale_color_manual(values=color_palette) + my_theme

g2<-ggplot(CTRL2.df, aes(x=X, y=Y, col=Compartment)) + geom_point(size=0.5) + 
  ggtitle("CTRL2") + scale_color_manual(values=color_palette) + my_theme

g3<-ggplot(CTRL3.df, aes(x=X, y=Y, col=Compartment)) + geom_point(size=0.5) + 
  ggtitle("CTRL3") + scale_color_manual(values=color_palette) + my_theme

g4<-ggplot(CTRL4.df, aes(x=X, y=Y, col=Compartment)) + geom_point(size=0.5) + 
  ggtitle("CTRL4") + scale_color_manual(values=color_palette) + my_theme

library(gridExtra)
pdf("Figures/Figure1/pdfs/CTRL_Compartment_Annotations_ggplot.pdf", height=3, width=12)
grid.arrange(g1,g2, g3, g4, ncol=4)
dev.off()

##IRL Tissue (24 hours warm ischemia-reperfusion injury in male mice):
## Note: IRL is same as AKI24 dataset in the main article section.
IRL1.df<-data.frame(irl1$pos)
IRL1.df$Compartment<-'Other'
IRL1.df$Compartment<-Hemb[rownames(IRL1.df),]$tissue
colnames(IRL1.df)<-c('X','Y','Compartment')

IRL2.df<-data.frame(irl2$pos)
IRL2.df$Compartment<-'Other'
IRL2.df$Compartment<-Hemb[rownames(IRL2.df),]$tissue
colnames(IRL2.df)<-c('X','Y','Compartment')

IRL3.df<-data.frame(irl3$pos)
IRL3.df$Compartment<-'Other'
IRL3.df$Compartment<-Hemb[rownames(IRL3.df),]$tissue
colnames(IRL3.df)<-c('X','Y','Compartment')

IRL4.df<-data.frame(irl4$pos)
IRL4.df$Compartment<-'Other'
IRL4.df$Compartment<-Hemb[rownames(IRL4.df),]$tissue
colnames(IRL4.df)<-c('X','Y','Compartment')

color_palette<-c("Cortex"="red", "O.Medulla"="blue","I.Medulla"="green","Other"="gray")

g1<-ggplot(IRL1.df, aes(x=X, y=Y, col=Compartment)) + geom_point(size=0.5) + 
  ggtitle("IRL1") + scale_color_manual(values=color_palette) + my_theme

g2<-ggplot(IRL2.df, aes(x=X, y=Y, col=Compartment)) + geom_point(size=0.5) + 
  ggtitle("IRL2") + scale_color_manual(values=color_palette) + my_theme

g3<-ggplot(IRL3.df, aes(x=X, y=Y, col=Compartment)) + geom_point(size=0.5) + 
  ggtitle("IRL3") + scale_color_manual(values=color_palette) + my_theme

g4<-ggplot(IRL4.df, aes(x=X, y=Y, col=Compartment)) + geom_point(size=0.5) + 
  ggtitle("IRL4") + scale_color_manual(values=color_palette) + my_theme

library(gridExtra)
pdf("Figures/Figure1/pdfs/IRL_Compartment_Annotations_ggplot.pdf", height=3, width=12)
grid.arrange(g1,g2, g3, g4, ncol=4)
dev.off()
