#Clean environment
rm(list=ls(all.names=TRUE))
gc()

dir<-"/Users/srujansingh/Library/CloudStorage/OneDrive-SharedLibraries-JohnsHopkins/Jean Fan - 2024_spatial_AKI/Revision2_GenomBiol"
setwd(dir)

#Load all libraries:
library(openxlsx)
library(ggplot2)

## load data
load("data/CIS_data.RData")
load("data/AKI_data.RData")

################################################################################
#Extracting the compartment specific spots for CIS and AKI:
#Note: 1.Interface is same as Outer Medulla. 
# 2. Medulla is same as Inner Medulla.

cortex <- readRDS('AKI-CIS-irl-ctrl_Cortex_spots.rds')
interface <- readRDS('AKI-CIS-irl-ctrl_Interface_spots.rds')
medulla <- readRDS('AKI-CIS-irl-ctrl_Medulla_spots.rds')

# ############################## CIS DATASET SPOTS ###############################
#1.All_Cortex_Spots
CIS_cortex.0h<-cortex[grep("^CIS_0h", names(cortex))]
CIS_cortex.12h<-cortex[grep("^CIS_12h", names(cortex))]
CIS_cortex.24h<-cortex[grep("^CIS_24h", names(cortex))]
CIS_cortex.48h<-cortex[grep("^CIS_48h", names(cortex))]

CIS_cortex.0h<-cbind(CIS_0h$pos[names(CIS_cortex.0h),], CIS_cortex.0h)
CIS_cortex.12h<-cbind(CIS_12h$pos[names(CIS_cortex.12h),], CIS_cortex.12h)
CIS_cortex.24h<-cbind(CIS_24h$pos[names(CIS_cortex.24h),], CIS_cortex.24h)
CIS_cortex.48h<-cbind(CIS_48h$pos[names(CIS_cortex.48h),], CIS_cortex.48h)

colnames(CIS_cortex.0h)<-c('X','Y', 'Cluster')
colnames(CIS_cortex.12h)<-c('X','Y', 'Cluster')
colnames(CIS_cortex.24h)<-c('X','Y', 'Cluster')
colnames(CIS_cortex.48h)<-c('X','Y', 'Cluster')

#2.All Interface spots
CIS_interface.0h<-interface[grep("^CIS_0h", names(interface))]
CIS_interface.12h<-interface[grep("^CIS_12h", names(interface))]
CIS_interface.24h<-interface[grep("^CIS_24h", names(interface))]
CIS_interface.48h<-interface[grep("^CIS_48h", names(interface))]

CIS_interface.0h<-cbind(CIS_0h$pos[names(CIS_interface.0h),], CIS_interface.0h)
CIS_interface.12h<-cbind(CIS_12h$pos[names(CIS_interface.12h),], CIS_interface.12h)
CIS_interface.24h<-cbind(CIS_24h$pos[names(CIS_interface.24h),], CIS_interface.24h)
CIS_interface.48h<-cbind(CIS_48h$pos[names(CIS_interface.48h),], CIS_interface.48h)

colnames(CIS_interface.0h)<-c('X','Y', 'Cluster')
colnames(CIS_interface.12h)<-c('X','Y', 'Cluster')
colnames(CIS_interface.24h)<-c('X','Y', 'Cluster')
colnames(CIS_interface.48h)<-c('X','Y', 'Cluster')

#3.All medullary spots
CIS_medulla.0h<-medulla[grep("^CIS_0h", names(medulla))]
CIS_medulla.12h<-medulla[grep("^CIS_12h", names(medulla))]
CIS_medulla.24h<-medulla[grep("^CIS_24h", names(medulla))]
CIS_medulla.48h<-medulla[grep("^CIS_48h", names(medulla))]

CIS_medulla.0h<-cbind(CIS_0h$pos[names(CIS_medulla.0h),], CIS_medulla.0h)
CIS_medulla.12h<-cbind(CIS_12h$pos[names(CIS_medulla.12h),], CIS_medulla.12h)
CIS_medulla.24h<-cbind(CIS_24h$pos[names(CIS_medulla.24h),], CIS_medulla.24h)
CIS_medulla.48h<-cbind(CIS_48h$pos[names(CIS_medulla.48h),], CIS_medulla.48h)

colnames(CIS_medulla.0h)<-c('X','Y', 'Cluster')
colnames(CIS_medulla.12h)<-c('X','Y', 'Cluster')
colnames(CIS_medulla.24h)<-c('X','Y', 'Cluster')
colnames(CIS_medulla.48h)<-c('X','Y', 'Cluster')

############################ normalize #########################################
## easier for plotting
CIS_0h$raw <- CIS_0h$gexp[unique(rownames(CIS_0h$gexp)),]
CIS_12h$raw <- CIS_12h$gexp[unique(rownames(CIS_0h$gexp)),]
CIS_24h$raw <- CIS_24h$gexp[unique(rownames(CIS_0h$gexp)),]
CIS_48h$raw <- CIS_48h$gexp[unique(rownames(CIS_0h$gexp)),]

CIS_0h$mat_notlog <- MERINGUE::normalizeCounts(CIS_0h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
CIS_12h$mat_notlog <- MERINGUE::normalizeCounts(CIS_12h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
CIS_24h$mat_notlog <- MERINGUE::normalizeCounts(CIS_24h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
CIS_48h$mat_notlog <- MERINGUE::normalizeCounts(CIS_48h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)

################################################################################
Cor.HALL<-read.xlsx("Codes2_Tables/CIS_HALLMARK_Pathays_Consistency_R2_Percentile95.xlsx", sheet="Cortex")
Cor.HALL<-data.frame(Cor.HALL)

Int.HALL<-read.xlsx("Codes2_Tables/CIS_HALLMARK_Pathays_Consistency_R2_Percentile95.xlsx", sheet="Outer_Medulla")
Int.HALL<-data.frame(Int.HALL)

Med.HALL<-read.xlsx("Codes2_Tables/CIS_HALLMARK_Pathays_Consistency_R2_Percentile95.xlsx", sheet="Inner_Medulla")
Med.HALL<-data.frame(Med.HALL)

Cor.OXPHOS<-Cor.HALL[Cor.HALL$ID=='HALLMARK_OXIDATIVE_PHOSPHORYLATION',]$core_enrichment #OXPHOS Pathway
Cor.OXPHOS<-trimws(unlist(strsplit(Cor.OXPHOS, ",")))

Int.OXPHOS<-Int.HALL[Int.HALL$ID=='HALLMARK_OXIDATIVE_PHOSPHORYLATION',]$core_enrichment #OXPHOS Pathway
Int.OXPHOS<-trimws(unlist(strsplit(Int.OXPHOS, ",")))

Med.OXPHOS<-Med.HALL[Med.HALL$ID=='HALLMARK_OXIDATIVE_PHOSPHORYLATION',]$core_enrichment #OXPHOS Pathway
Med.OXPHOS<-trimws(unlist(strsplit(Med.OXPHOS, ",")))

mytheme<-theme_minimal() + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add box around plot
        axis.text.x = element_text(size = 14, color="black"), # Increase font size of pathway names
        axis.text.y = element_text(size = 14, color="black"),
        axis.title.x = element_text(size = 16, color="black"),
        axis.title.y = element_text(size = 16, color="black"),
        axis.ticks = element_line(color = "black"),  # Add tick marks
        axis.ticks.length = unit(0.2, "cm"),  # Adjust tick mark length
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold")) 


###########################RAW COUNTS VISUALIZATIONS ###########################
#All Spots:
t0<-colSums(as.matrix(CIS_0h$raw))
t1<-colSums(as.matrix(CIS_12h$raw))
t2<-colSums(as.matrix(CIS_24h$raw))
t3<-colSums(as.matrix(CIS_48h$raw))

df<-data.frame(gexp=c(t0,t1,t2,t3), time=c(rep(0, length(t0)),rep(12, length(t1)),
              rep(24, length(t2)),rep(48, length(t3)) ))

ggplot(df, aes(x=factor(time), y=gexp)) + 
  #geom_point(size=0.5) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.6) +
  ggtitle("Total Counts Whole (Raw)") +
  labs(x="Cold Ischemia Time (hours)", y="Raw Gene Counts")

##xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#CORTEX SPOTS
#TOTAL GENE COUNT (Cortex Spots):
t0_cor<-colSums(as.matrix(CIS_0h$raw[,rownames(CIS_cortex.0h)]))
t1_cor<-colSums(as.matrix(CIS_12h$raw[,rownames(CIS_cortex.12h)]))
t2_cor<-colSums(as.matrix(CIS_24h$raw[,rownames(CIS_cortex.24h)]))
t3_cor<-colSums(as.matrix(CIS_48h$raw[,rownames(CIS_cortex.48h)]))

df<-data.frame(gexp=c(t0_cor,t1_cor,t2_cor,t3_cor), 
               time=c(rep(0, length(t0_cor)),rep(12, length(t1_cor)),
               rep(24, length(t2_cor)),rep(48, length(t3_cor)) ))

res<-lm(df$gexp~df$time)
slope <- coef(res)[2]
intercept <- coef(res)[1]
R2<-summary(res)$r.squared

file<-paste0(dir,"/Codes2_Figures/Figure_S5/pdfs/Cortex_Total_mRNA_Counts_Raw.pdf")
pdf(file, width=7, height=5)
ggplot(df, aes(x=time, y=gexp)) + 
  #geom_point(size=0.5) +
  geom_boxplot(aes(group=factor(time)),outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.6) +
  ggtitle("Total Counts Cortex (Raw)") +
  geom_abline(slope = slope, intercept = intercept, color = "blue", linewidth = 0.5) +
  labs(x="Cold Ischemia Time (hours)", y="Raw Gene Counts") + mytheme +
  annotate("text", label=paste0("Slope = ", signif(slope,2), "\nR2=", 
                                signif(R2,2)), x=0, y=190000, hjust=0, size=5, color="blue")
dev.off()

#OXPHOS GENE RAW COUNTS (Cortex):~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
t0_cor.OXPHOS<-colSums(as.matrix(CIS_0h$raw[Cor.OXPHOS,rownames(CIS_cortex.0h)]))
t1_cor.OXPHOS<-colSums(as.matrix(CIS_12h$raw[Cor.OXPHOS,rownames(CIS_cortex.12h)]))
t2_cor.OXPHOS<-colSums(as.matrix(CIS_24h$raw[Cor.OXPHOS,rownames(CIS_cortex.24h)]))
t3_cor.OXPHOS<-colSums(as.matrix(CIS_48h$raw[Cor.OXPHOS,rownames(CIS_cortex.48h)]))

df<-data.frame(gexp=c(t0_cor.OXPHOS,t1_cor.OXPHOS,t2_cor.OXPHOS,t3_cor.OXPHOS), 
               time=c(rep(0, length(t0_cor.OXPHOS)),rep(12, length(t1_cor.OXPHOS)),
                      rep(24, length(t2_cor.OXPHOS)),rep(48, length(t3_cor.OXPHOS)) ))

res<-lm(df$gexp~df$time)
slope <- coef(res)[2]
intercept <- coef(res)[1]
R2<-summary(res)$r.squared

file<-paste0(dir,"/Codes2_Figures/Figure_S5/pdfs/Cortex_OXPHOS_mRNA_Counts_Raw.pdf")
pdf(file, width=7, height=5)

ggplot(df, aes(x=time, y=gexp)) + 
  #geom_point(size=0.5) + 
  geom_boxplot(aes(group=factor(time)),outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.6) +
  geom_abline(slope = slope, intercept = intercept, color = "blue", linewidth = 0.5) +
  ggtitle("OXPHOS Counts Cortex (Raw)") +
  labs(x="Cold Ischemia Time (hours)", y="Raw Gene Counts") + mytheme +
  annotate("text", label=paste0("Slope = ", signif(slope,2), "\nR2=", 
                                signif(R2,2)), x=0, y=4000, hjust=0, size=5, color="blue")
dev.off()

#OXPHOS Genes (Cortex) PROPORTIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
t0_cor.prop<-t0_cor.OXPHOS/t0_cor
t1_cor.prop<-t1_cor.OXPHOS/t1_cor
t2_cor.prop<-t2_cor.OXPHOS/t2_cor
t3_cor.prop<-t3_cor.OXPHOS/t3_cor

df<-data.frame(gexp=c(t0_cor.prop,t1_cor.prop,t2_cor.prop,t3_cor.prop), 
               time=c(rep(0, length(t0_cor.prop)),rep(12, length(t1_cor.prop)),
                      rep(24, length(t2_cor.prop)),rep(48, length(t3_cor.prop)) ))

res<-lm(df$gexp~df$time)
slope <- coef(res)[2]
intercept <- coef(res)[1]
R2<-summary(res)$r.squared

file<-paste0(dir,"/Codes2_Figures/Figure_S5/pdfs/Cortex_OXPHOS_mRNA_Counts_Proportion.pdf")
pdf(file, width=7, height=5)
ggplot(df, aes(x=time, y=gexp*100)) + 
  #geom_point(size=0.5) + 
  geom_boxplot(aes(group=factor(time)),outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.6) +
  geom_abline(slope = slope*100, intercept = intercept*100, color = "blue", linewidth = 0.5) +
  ggtitle("OXPHOS Proportion Cortex (Raw)") + 
  labs(x="Cold Ischemia Time (hours)", y="Proportion (%)") + mytheme +
  annotate("text", label=paste0("Slope = ", signif(slope*100,2), "\nR2=", 
                                signif(R2,2)), x=0, y=2.8, hjust=0, size=5, color="blue")
dev.off()

##xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#INTERFACE SPOTS
#TOTAL GENE COUNT (Interface Spots):
t0_int<-colSums(as.matrix(CIS_0h$raw[,rownames(CIS_interface.0h)]))
t1_int<-colSums(as.matrix(CIS_12h$raw[,rownames(CIS_interface.12h)]))
t2_int<-colSums(as.matrix(CIS_24h$raw[,rownames(CIS_interface.24h)]))
t3_int<-colSums(as.matrix(CIS_48h$raw[,rownames(CIS_interface.48h)]))

df<-data.frame(gexp=c(t0_int,t1_int,t2_int,t3_int), 
               time=c(rep(0, length(t0_int)),rep(12, length(t1_int)),
                      rep(24, length(t2_int)),rep(48, length(t3_int)) ))

res<-lm(df$gexp~df$time)
slope <- coef(res)[2]
intercept <- coef(res)[1]
R2<-summary(res)$r.squared

file<-paste0(dir,"/Codes2_Figures/Figure_S5/pdfs/Interface_Total_mRNA_Counts_Raw.pdf")
pdf(file, width=7, height=5)
ggplot(df, aes(x=time, y=gexp)) + 
  #geom_point(size=0.5) +
  geom_boxplot(aes(group=factor(time)),outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.6) +
  geom_abline(slope = slope, intercept = intercept, color = "blue", linewidth = 0.5) +
  ggtitle("Total Counts Interface (Raw)") +
  labs(x="Cold Ischemia Time (hours)", y="Raw Gene Counts") + mytheme +
  annotate("text", label=paste0("Slope = ", signif(slope,2), "\nR2=", 
                                signif(R2,2)), x=0, y=110000, hjust=0, size=5, color="blue")
dev.off()

#OXPHOS GENE RAW COUNTS (Interface):~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
t0_int.OXPHOS<-colSums(as.matrix(CIS_0h$raw[Int.OXPHOS,rownames(CIS_interface.0h)]))
t1_int.OXPHOS<-colSums(as.matrix(CIS_12h$raw[Int.OXPHOS,rownames(CIS_interface.12h)]))
t2_int.OXPHOS<-colSums(as.matrix(CIS_24h$raw[Int.OXPHOS,rownames(CIS_interface.24h)]))
t3_int.OXPHOS<-colSums(as.matrix(CIS_48h$raw[Int.OXPHOS,rownames(CIS_interface.48h)]))


df<-data.frame(gexp=c(t0_int.OXPHOS,t1_int.OXPHOS,t2_int.OXPHOS,t3_int.OXPHOS), 
               time=c(rep(0, length(t0_int.OXPHOS)),rep(12, length(t1_int.OXPHOS)),
                      rep(24, length(t2_int.OXPHOS)),rep(48, length(t3_int.OXPHOS)) ))

res<-lm(df$gexp~df$time)
slope <- coef(res)[2]
intercept <- coef(res)[1]
R2<-summary(res)$r.squared

file<-paste0(dir,"/Codes2_Figures/Figure_S5/pdfs/Interface_OXPHOS_mRNA_Counts_Raw.pdf")
pdf(file, width=7, height=5)
ggplot(df, aes(x=time, y=gexp)) + 
  #geom_point(size=0.5) + 
  geom_boxplot(aes(group=factor(time)),outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.6) +
  geom_abline(slope = slope, intercept = intercept, color = "blue", linewidth = 0.5) +
  ggtitle("OXPHOS Counts Interface (Raw)") +
  labs(x="Cold Ischemia Time (hours)", y="Raw Gene Counts") + mytheme +
  annotate("text", label=paste0("Slope = ", signif(slope,2), "\nR2=", 
                                signif(R2,2)), x=0, y=2500, hjust=0, size=5, color="blue")
dev.off()
 
#OXPHOS Genes (Interface) PROPORTIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
t0_int.prop<-t0_int.OXPHOS/t0_int
t1_int.prop<-t1_int.OXPHOS/t1_int
t2_int.prop<-t2_int.OXPHOS/t2_int
t3_int.prop<-t3_int.OXPHOS/t3_int

df<-data.frame(gexp=c(t0_int.prop,t1_int.prop,t2_int.prop,t3_int.prop), 
               time=c(rep(0, length(t0_int.prop)),rep(12, length(t1_int.prop)),
                      rep(24, length(t2_int.prop)),rep(48, length(t3_int.prop)) ))

res<-lm(df$gexp~df$time)
slope <- coef(res)[2]
intercept <- coef(res)[1]
R2<-summary(res)$r.squared

file<-paste0(dir,"/Codes2_Figures/Figure_S5/pdfs/Interface_OXPHOS_mRNA_Counts_Proportion.pdf")
pdf(file, width=7, height=5)
ggplot(df, aes(x=time, y=gexp*100)) + 
  #geom_point(size=0.5) + 
  geom_boxplot(aes(group=factor(time)),outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.6) +
  geom_abline(slope = slope*100, intercept = intercept*100, color = "blue", linewidth = 0.5) +
  ggtitle("OXPHOS Proportion Interface (Raw)") +
  labs(x="Cold Ischemia Time (hours)", y="Proportion (%)") + mytheme +
  annotate("text", label=paste0("Slope = ", signif(slope*100,2), "\nR2=", 
                                signif(R2,2)), x=0, y=3.8, hjust=0, size=5, color="blue")
dev.off()

##xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#Medulla Spots:
t0_med<-colSums(as.matrix(CIS_0h$raw[,rownames(CIS_medulla.0h)]))
t1_med<-colSums(as.matrix(CIS_12h$raw[,rownames(CIS_medulla.12h)]))
t2_med<-colSums(as.matrix(CIS_24h$raw[,rownames(CIS_medulla.24h)]))
t3_med<-colSums(as.matrix(CIS_48h$raw[,rownames(CIS_medulla.48h)]))

df<-data.frame(gexp=c(t0_med,t1_med,t2_med,t3_med), 
               time=c(rep(0, length(t0_med)),rep(12, length(t1_med)),
                      rep(24, length(t2_med)),rep(48, length(t3_med)) ))

res<-lm(df$gexp~df$time)
summary(res)
slope <- coef(res)[2]
intercept <- coef(res)[1]
R2<-summary(res)$r.squared

file<-paste0(dir,"/Codes2_Figures/Figure_S5/pdfs/Medulla_Total_mRNA_Counts_Raw.pdf")
pdf(file, width=7, height=5)
ggplot(df, aes(x=time, y=gexp)) + 
  #geom_point(size=0.5) + 
  geom_boxplot(aes(group=factor(time)),outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.6) +
  geom_abline(slope = slope, intercept = intercept, color = "blue", linewidth = 0.5) +
  ggtitle("Total Counts Medulla (Raw)") +
  labs(x="Cold Ischemia Time (hours)", y="Raw Gene Counts") + mytheme +
  annotate("text", label=paste0("Slope = ", signif(slope,2), "\nR2=", 
                                signif(R2,2)), x=0, y=53000, hjust=0, size=5, color="blue")

dev.off()

#OXPHOS Genes Only (Medulla) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
t0_med.OXPHOS<-colSums(as.matrix(CIS_0h$raw[Med.OXPHOS,rownames(CIS_medulla.0h)]))
t1_med.OXPHOS<-colSums(as.matrix(CIS_12h$raw[Med.OXPHOS,rownames(CIS_medulla.12h)]))
t2_med.OXPHOS<-colSums(as.matrix(CIS_24h$raw[Med.OXPHOS,rownames(CIS_medulla.24h)]))
t3_med.OXPHOS<-colSums(as.matrix(CIS_48h$raw[Med.OXPHOS,rownames(CIS_medulla.48h)]))


df<-data.frame(gexp=c(t0_med.OXPHOS,t1_med.OXPHOS,t2_med.OXPHOS,t3_med.OXPHOS), 
               time=c(rep(0, length(t0_med.OXPHOS)),rep(12, length(t1_med.OXPHOS)),
                      rep(24, length(t2_med.OXPHOS)),rep(48, length(t3_med.OXPHOS)) ))

res<-lm(df$gexp~df$time)
summary(res)
slope <- coef(res)[2]
intercept <- coef(res)[1]
R2<-summary(res)$r.squared

file<-paste0(dir,"/Codes2_Figures/Figure_S5/pdfs/Medulla_OXPHOS_mRNA_Counts_Raw.pdf")
pdf(file, width=7, height=5)
ggplot(df, aes(x=time, y=gexp)) + 
  #geom_point(size=0.5) + 
  geom_boxplot(aes(group = factor(time)),outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.6) +
  geom_abline(slope = slope, intercept = intercept, color = "blue", linewidth = 0.5) +
  ggtitle("OXPHOS Counts Medulla (Raw)") +
  labs(x="Cold Ischemia Time (hours)", y="Raw Gene Counts") + mytheme +
  annotate("text", label=paste0("Slope = ", signif(slope,2), "\nR2=", 
                                signif(R2,2)), x=0, y=1800, hjust=0, size=5, color="blue")
dev.off()

#OXPHOS Genes (Medulla) PROPORTIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
t0_med.prop<-t0_med.OXPHOS/t0_med
t1_med.prop<-t1_med.OXPHOS/t1_med
t2_med.prop<-t2_med.OXPHOS/t2_med
t3_med.prop<-t3_med.OXPHOS/t3_med

df<-data.frame(gexp=c(t0_med.prop,t1_med.prop,t2_med.prop,t3_med.prop), 
               time=c(rep(0, length(t0_med.prop)),rep(12, length(t1_med.prop)),
                      rep(24, length(t2_med.prop)),rep(48, length(t3_med.prop)) ))

res<-lm(df$gexp~df$time)
slope <- coef(res)[2]
intercept <- coef(res)[1]
R2<-summary(res)$r.squared

file<-paste0(dir,"/Codes2_Figures/Figure_S5/pdfs/Medulla_OXPHOS_mRNA_Counts_Proportion.pdf")
pdf(file, width=7, height=5)
#Proportion
ggplot(df, aes(x=time, y=gexp*100)) + 
  #geom_point(size=0.5) + 
  geom_boxplot(aes(group = factor(time)), outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.6) +
  geom_abline(slope = slope*100, intercept = intercept*100, color = "blue", linewidth = 0.5) +
  ggtitle("OXPHOS Proportion Medulla (Raw)") +
  labs(x="Cold Ischemia Time (hours)", y="Proportion (%)") + mytheme +
  annotate("text", label=paste0("Slope = ", signif(slope*100,2), "\nR2=", 
                                signif(R2,2)), x=0, y=6, hjust=0, size=5, color="blue")
dev.off()


