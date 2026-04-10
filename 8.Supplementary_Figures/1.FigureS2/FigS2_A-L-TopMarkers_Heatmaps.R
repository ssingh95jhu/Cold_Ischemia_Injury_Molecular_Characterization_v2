#Clean environment
rm(list=ls(all.names=TRUE))
gc()

setwd("/Users/srujansingh/Library/CloudStorage/OneDrive-SharedLibraries-JohnsHopkins/Jean Fan - 2024_spatial_AKI/Manuscript")

########################## ADDING LIBRARY ######################################
library(openxlsx)
library(ggplot2)
library(pheatmap)

# loading all the 17 spatial transcriptomics datasets belonging to the different
#experimental groups.
#Note: irl is same as AKI24 group as described in the main experimenatl section.
load("data/CIS_data.RData")
load("data/AKI_data.RData")
load("data/Rabb_ctrl.RData")
load("data/Rabb_irl24h.RData")

###################### AKI DATASET #############################################
AKI_2d$gexp<-AKI_2d$gexp[,-1818] #This spot has no genes
AKI_2d$pos<-AKI_2d$pos[colnames(AKI_2d$gexp),]

##Removing all mitochondrial genes from the gene experssion table
AKI_sham$gexp<-AKI_sham$gexp[-grep("^mt", rownames(AKI_sham$gexp)),]
AKI_4h$gexp<-AKI_4h$gexp[-grep("^mt", rownames(AKI_4h$gexp)),]
AKI_12h$gexp<-AKI_12h$gexp[-grep("^mt", rownames(AKI_12h$gexp)),]
AKI_2d$gexp<-AKI_2d$gexp[-grep("^mt", rownames(AKI_2d$gexp)),]
AKI_6w$gexp<-AKI_6w$gexp[-grep("^mt", rownames(AKI_6w$gexp)),]

##Updating the position matrix for AKI
AKI_sham$pos<-AKI_sham$pos[grep("^AKI_sham", colnames(AKI_sham$gexp)),]
AKI_4h$pos<-AKI_4h$pos[grep("^AKI_4h", colnames(AKI_4h$gexp)),]
AKI_12h$pos<-AKI_12h$pos[grep("^AKI_12h", colnames(AKI_12h$gexp)),]
AKI_2d$pos<-AKI_2d$pos[grep("^AKI_2d", colnames(AKI_2d$gexp)),]
AKI_6w$pos<-AKI_6w$pos[grep("^AKI_6w", colnames(AKI_6w$gexp)),]


## Read top marker genes #######################################################
df1<-read.xlsx("Supplementary_Tables/Fig1_TopMarkers_Kidney_Compartments.xlsx", sheet="Cortex")
cortex_genes<-data.frame(Genes=df1$Genes, L2FC=df1$log2FoldChange, pval=df1$pvalue)

df1<-read.xlsx("Supplementary_Tables/Fig1_TopMarkers_Kidney_Compartments.xlsx", sheet="Interface")
interface_genes<-data.frame(Genes=df1$Genes, L2FC=df1$log2FoldChange, pval=df1$pvalue)

df1<-read.xlsx("Supplementary_Tables/Fig1_TopMarkers_Kidney_Compartments.xlsx", sheet="Medulla")
medulla_genes<-data.frame(Genes=df1$Genes, L2FC=df1$log2FoldChange, pval=df1$pvalue)


## Spots #######################################################################
cortex <- readRDS('AKI-CIS-irl-ctrl_Cortex_spots.rds')
interface <- readRDS('AKI-CIS-irl-ctrl_Interface_spots.rds')
medulla <- readRDS('AKI-CIS-irl-ctrl_Medulla_spots.rds')

##### CIS DATASETS HEATMAPS ####################################################
#1.All_Cortex_Spots
CIS_cortex.0h<-names(cortex[grep("^CIS_0h", names(cortex))] )
CIS_cortex.12h<-names(cortex[grep("^CIS_12h", names(cortex))] )
CIS_cortex.24h<-names(cortex[grep("^CIS_24h", names(cortex))] )
CIS_cortex.48h<-names(cortex[grep("^CIS_48h", names(cortex))] )

#2.All Interface spots
CIS_interface.0h<-names(interface[grep("^CIS_0h", names(interface))] )
CIS_interface.12h<-names(interface[grep("^CIS_12h", names(interface))] )
CIS_interface.24h<-names(interface[grep("^CIS_24h", names(interface))] )
CIS_interface.48h<-names(interface[grep("^CIS_48h", names(interface))] )

#3.All medullary spots
CIS_medulla.0h<-names(medulla[grep("^CIS_0h", names(medulla))] )
CIS_medulla.12h<-names(medulla[grep("^CIS_12h", names(medulla))] )
CIS_medulla.24h<-names(medulla[grep("^CIS_24h", names(medulla))] )
CIS_medulla.48h<-names(medulla[grep("^CIS_48h", names(medulla))] )

### COUNT NORMALIZATION OF DATA ################################################
CIS_0h$mat_notlog <- MERINGUE::normalizeCounts(CIS_0h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
CIS_12h$mat_notlog <- MERINGUE::normalizeCounts(CIS_12h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
CIS_24h$mat_notlog <- MERINGUE::normalizeCounts(CIS_24h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
CIS_48h$mat_notlog <- MERINGUE::normalizeCounts(CIS_48h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)

CIS_0h.cortex.cd<-CIS_0h$mat_notlog[,CIS_cortex.0h]
CIS_12h.cortex.cd<-CIS_12h$mat_notlog[,CIS_cortex.12h]
CIS_24h.cortex.cd<-CIS_24h$mat_notlog[,CIS_cortex.24h]
CIS_48h.cortex.cd<-CIS_48h$mat_notlog[,CIS_cortex.48h]

CIS_0h.interface.cd<-CIS_0h$mat_notlog[,CIS_interface.0h]
CIS_12h.interface.cd<-CIS_12h$mat_notlog[,CIS_interface.12h]
CIS_24h.interface.cd<-CIS_24h$mat_notlog[,CIS_interface.24h]
CIS_48h.interface.cd<-CIS_48h$mat_notlog[,CIS_interface.48h]

CIS_0h.medulla.cd<-CIS_0h$mat_notlog[,CIS_medulla.0h]
CIS_12h.medulla.cd<-CIS_12h$mat_notlog[,CIS_medulla.12h]
CIS_24h.medulla.cd<-CIS_24h$mat_notlog[,CIS_medulla.24h]
CIS_48h.medulla.cd<-CIS_48h$mat_notlog[,CIS_medulla.48h]

#HEATMAP CIS (TOP CORTEX GENES)
CIS.cd<-cbind(CIS_0h.cortex.cd, CIS_0h.interface.cd, CIS_0h.medulla.cd, 
              CIS_12h.cortex.cd, CIS_12h.interface.cd, CIS_12h.medulla.cd, 
              CIS_24h.cortex.cd, CIS_24h.interface.cd, CIS_24h.medulla.cd, 
              CIS_48h.cortex.cd, CIS_48h.interface.cd, CIS_48h.medulla.cd)


CIS.cortex.cd<-CIS.cd[cortex_genes$Genes[1:30],]

sample_annotations <- data.frame(
  Timepoint = factor(c(rep("0h", ncol(CIS_0h.cortex.cd)), rep("0h", ncol(CIS_0h.interface.cd)),
                       rep("0h", ncol(CIS_0h.medulla.cd)), rep("12h", ncol(CIS_12h.cortex.cd)),
                       rep("12h", ncol(CIS_12h.interface.cd)), rep("12h", ncol(CIS_12h.medulla.cd)),
                       rep("24h", ncol(CIS_24h.cortex.cd)), rep("24h", ncol(CIS_24h.interface.cd)),
                       rep("24h", ncol(CIS_24h.medulla.cd)), rep("48h", ncol(CIS_48h.cortex.cd)),
                       rep("48h", ncol(CIS_48h.interface.cd)), rep("48h", ncol(CIS_48h.medulla.cd)) ), 
                        levels = c("0h", "12h", "24h", "48h")) ,  # Ensures correct order
  
  TissueType = factor(c(rep("Cortex", ncol(CIS_0h.cortex.cd)), rep("Interface", ncol(CIS_0h.interface.cd)),
                        rep("Medulla", ncol(CIS_0h.medulla.cd)), rep("Cortex", ncol(CIS_12h.cortex.cd)),
                        rep("Interface", ncol(CIS_12h.interface.cd)), rep("Medulla", ncol(CIS_12h.medulla.cd)),
                        rep("Cortex", ncol(CIS_24h.cortex.cd)), rep("Interface", ncol(CIS_24h.interface.cd)),
                        rep("Medulla", ncol(CIS_24h.medulla.cd)), rep("Cortex", ncol(CIS_48h.cortex.cd)),
                        rep("Interface", ncol(CIS_48h.interface.cd)), rep("Medulla", ncol(CIS_48h.medulla.cd)) ), 
                      levels = c("Cortex", "Interface", "Medulla"))  # Example categories
)
rownames(sample_annotations)<-colnames(CIS.cortex.cd)

#Scaling the matrix
mat <- as.matrix(CIS.cortex.cd)
scaled_mat <- t(scale(t(mat)))
  
overall_min<-min(scaled_mat)
overall_max <- quantile(scaled_mat, 0.95)

custom_colors<-colorRampPalette(c("black", "blue", "purple", "yellow"))(100)
breaks <- seq(overall_min, overall_max, length.out = 101)

# Define colors for annotations
ann_colors <- list(
  Timepoint = c("0h" = "lavender", "12h" = "mediumpurple", "24h" = "pink", "48h" = "magenta"),
  TissueType = c("Cortex" = "red", "Interface" = "blue", "Medulla" = "green")
)

# Draw heatmap with annotations
pdf("Figures/Figure1S/pdfs/CIS_Cortex_Marker_Genes_Heatmap.pdf", height=6, width=8)
pheatmap(scaled_mat, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, show_rownames = TRUE,
  color = custom_colors, breaks= breaks,
  annotation_col = sample_annotations,  # Add color-coded sample info
  annotation_colors = ann_colors,  # Apply custom colors
  main = "Top Gene Markers for Cortex (CIS)"
)
dev.off()

##Heatmap for CIS interface:
#Note: interface is same as Outer Medulla as described in the main article section.
CIS.interface.cd<-CIS.cd[interface_genes$Genes[1:30],]

sample_annotations <- data.frame(
  Timepoint = factor(c(rep("0h", ncol(CIS_0h.cortex.cd)), rep("0h", ncol(CIS_0h.interface.cd)),
                       rep("0h", ncol(CIS_0h.medulla.cd)), rep("12h", ncol(CIS_12h.cortex.cd)),
                       rep("12h", ncol(CIS_12h.interface.cd)), rep("12h", ncol(CIS_12h.medulla.cd)),
                       rep("24h", ncol(CIS_24h.cortex.cd)), rep("24h", ncol(CIS_24h.interface.cd)),
                       rep("24h", ncol(CIS_24h.medulla.cd)), rep("48h", ncol(CIS_48h.cortex.cd)),
                       rep("48h", ncol(CIS_48h.interface.cd)), rep("48h", ncol(CIS_48h.medulla.cd)) ), 
                     levels = c("0h", "12h", "24h", "48h")) ,  # Ensures correct order
  
  TissueType = factor(c(rep("Cortex", ncol(CIS_0h.cortex.cd)), rep("Interface", ncol(CIS_0h.interface.cd)),
                        rep("Medulla", ncol(CIS_0h.medulla.cd)), rep("Cortex", ncol(CIS_12h.cortex.cd)),
                        rep("Interface", ncol(CIS_12h.interface.cd)), rep("Medulla", ncol(CIS_12h.medulla.cd)),
                        rep("Cortex", ncol(CIS_24h.cortex.cd)), rep("Interface", ncol(CIS_24h.interface.cd)),
                        rep("Medulla", ncol(CIS_24h.medulla.cd)), rep("Cortex", ncol(CIS_48h.cortex.cd)),
                        rep("Interface", ncol(CIS_48h.interface.cd)), rep("Medulla", ncol(CIS_48h.medulla.cd)) ), 
                      levels = c("Cortex", "Interface", "Medulla"))  # Example categories
)
rownames(sample_annotations)<-colnames(CIS.interface.cd)

mat <- as.matrix(CIS.interface.cd)

#Scaling the matrix
mat <- as.matrix(CIS.interface.cd)
scaled_mat <- t(scale(t(mat)))

overall_min<-min(scaled_mat)
#max(scaled_expr_mat)
overall_max <- quantile(scaled_mat, 0.95)

custom_colors<-colorRampPalette(c("black", "blue", "purple", "yellow"))(100)
breaks <- seq(overall_min, overall_max, length.out = 101)

# Draw heatmap with annotations
pdf("Figures/Figure1S/pdfs/CIS_Interface_Marker_Genes_Heatmap.pdf", height=6, width=8)
pheatmap(scaled_mat, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, show_rownames = TRUE,
         color = custom_colors, breaks= breaks,
         annotation_col = sample_annotations,  # Add color-coded sample info
         annotation_colors = ann_colors,  # Apply custom colors
         main = "Top Gene Markers for Interface (CIS)"
)
dev.off()

#Heatmap for Medulla (CIS)
#Note: medulla is same as Inner Medulla as described in the main article section.
CIS.medulla.cd<-CIS.cd[medulla_genes$Genes[1:30],]

sample_annotations <- data.frame(
  Timepoint = factor(c(rep("0h", ncol(CIS_0h.cortex.cd)), rep("0h", ncol(CIS_0h.interface.cd)),
                       rep("0h", ncol(CIS_0h.medulla.cd)), rep("12h", ncol(CIS_12h.cortex.cd)),
                       rep("12h", ncol(CIS_12h.interface.cd)), rep("12h", ncol(CIS_12h.medulla.cd)),
                       rep("24h", ncol(CIS_24h.cortex.cd)), rep("24h", ncol(CIS_24h.interface.cd)),
                       rep("24h", ncol(CIS_24h.medulla.cd)), rep("48h", ncol(CIS_48h.cortex.cd)),
                       rep("48h", ncol(CIS_48h.interface.cd)), rep("48h", ncol(CIS_48h.medulla.cd)) ), 
                     levels = c("0h", "12h", "24h", "48h")) ,  # Ensures correct order
  
  TissueType = factor(c(rep("Cortex", ncol(CIS_0h.cortex.cd)), rep("Interface", ncol(CIS_0h.interface.cd)),
                        rep("Medulla", ncol(CIS_0h.medulla.cd)), rep("Cortex", ncol(CIS_12h.cortex.cd)),
                        rep("Interface", ncol(CIS_12h.interface.cd)), rep("Medulla", ncol(CIS_12h.medulla.cd)),
                        rep("Cortex", ncol(CIS_24h.cortex.cd)), rep("Interface", ncol(CIS_24h.interface.cd)),
                        rep("Medulla", ncol(CIS_24h.medulla.cd)), rep("Cortex", ncol(CIS_48h.cortex.cd)),
                        rep("Interface", ncol(CIS_48h.interface.cd)), rep("Medulla", ncol(CIS_48h.medulla.cd)) ), 
                      levels = c("Cortex", "Interface", "Medulla"))  # Example categories
)
rownames(sample_annotations)<-colnames(CIS.medulla.cd)

#Scaling the matrix
mat <- as.matrix(CIS.medulla.cd)
scaled_mat <- t(scale(t(mat)))

overall_min<-min(scaled_mat)
overall_max <- quantile(scaled_mat, 0.95)

custom_colors<-colorRampPalette(c("black", "blue", "purple", "yellow"))(100)
breaks <- seq(overall_min, overall_max, length.out = 101)

# Define colors for annotations
ann_colors <- list(
  Timepoint = c("0h" = "lavender", "12h" = "mediumpurple", "24h" = "pink", "48h" = "magenta"),
  TissueType = c("Cortex" = "red", "Interface" = "blue", "Medulla" = "green")
)

# Draw heatmap with annotations
pdf("Figures/Figure1S/pdfs/CIS_Medulla_Marker_Genes_Heatmap.pdf", height=6, width=8)
pheatmap(scaled_mat, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, show_rownames = TRUE,
         color = custom_colors, breaks= breaks,
         annotation_col = sample_annotations,  # Add color-coded sample info
         annotation_colors = ann_colors,  # Apply custom colors
         main = "Top Gene Markers for Medulla (CIS)"
)
dev.off()

##### AKI DATASETS HEATMAPS ####################################################
#1.All_Cortex_Spots
AKI_cortex.sham<-names(cortex[grep("^AKI_sham", names(cortex))] )
AKI_cortex.4h<-names(cortex[grep("^AKI_4h", names(cortex))] )
AKI_cortex.12h<-names(cortex[grep("^AKI_12h", names(cortex))] )
AKI_cortex.2d<-names(cortex[grep("^AKI_2d", names(cortex))] )
AKI_cortex.6w<-names(cortex[grep("^AKI_6w", names(cortex))] )

#2.All Interface spots
AKI_interface.sham<-names(interface[grep("^AKI_sham", names(interface))] )
AKI_interface.4h<-names(interface[grep("^AKI_4h", names(interface))] )
AKI_interface.12h<-names(interface[grep("^AKI_12h", names(interface))] )
AKI_interface.2d<-names(interface[grep("^AKI_2d", names(interface))] )
AKI_interface.6w<-names(interface[grep("^AKI_6w", names(interface))] )

#3.All medullary spots
AKI_medulla.sham<-names(medulla[grep("^AKI_sham", names(medulla))] )
AKI_medulla.4h<-names(medulla[grep("^AKI_4h", names(medulla))] )
AKI_medulla.12h<-names(medulla[grep("^AKI_12h", names(medulla))] )
AKI_medulla.2d<-names(medulla[grep("^AKI_2d", names(medulla))] )
AKI_medulla.6w<-names(medulla[grep("^AKI_6w", names(medulla))] )

### COUNT NORMALIZATION OF DATA ################################################
AKI_sham$mat_notlog <- MERINGUE::normalizeCounts(AKI_sham$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
AKI_4h$mat_notlog <- MERINGUE::normalizeCounts(AKI_4h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
AKI_12h$mat_notlog <- MERINGUE::normalizeCounts(AKI_12h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
AKI_2d$mat_notlog <- MERINGUE::normalizeCounts(AKI_2d$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
#AKI_6w$mat_notlog <- MERINGUE::normalizeCounts(AKI_6w$gexp[unique(rownames(AKI_6w$gexp)),], log=FALSE)

AKI_sham.cortex.cd<-AKI_sham$mat_notlog[,AKI_cortex.sham]
AKI_4h.cortex.cd<-AKI_4h$mat_notlog[,AKI_cortex.4h]
AKI_12h.cortex.cd<-AKI_12h$mat_notlog[,AKI_cortex.12h]
AKI_2d.cortex.cd<-AKI_2d$mat_notlog[,AKI_cortex.2d]
#CIS_6w.cortex.cd<-CIS_48h$mat_notlog[,CIS_cortex.48h]

AKI_sham.interface.cd<-AKI_sham$mat_notlog[,AKI_interface.sham]
AKI_4h.interface.cd<-AKI_4h$mat_notlog[,AKI_interface.4h]
AKI_12h.interface.cd<-AKI_12h$mat_notlog[,AKI_interface.12h]
AKI_2d.interface.cd<-AKI_2d$mat_notlog[,AKI_interface.2d]

AKI_sham.medulla.cd<-AKI_sham$mat_notlog[,AKI_medulla.sham]
AKI_4h.medulla.cd<-AKI_4h$mat_notlog[,AKI_medulla.4h]
AKI_12h.medulla.cd<-AKI_12h$mat_notlog[,AKI_medulla.12h]
AKI_2d.medulla.cd<-AKI_2d$mat_notlog[,AKI_medulla.2d]

#HEATMAP for AKI:
AKI.cd<-cbind(AKI_sham.cortex.cd, AKI_sham.interface.cd, AKI_sham.medulla.cd, 
              AKI_4h.cortex.cd, AKI_4h.interface.cd, AKI_4h.medulla.cd,
              AKI_12h.cortex.cd, AKI_12h.interface.cd, AKI_12h.medulla.cd, 
              AKI_2d.cortex.cd, AKI_2d.interface.cd, AKI_2d.medulla.cd)

#HEATMAP (TOP CORTEX GENES)
AKI.cortex.cd<-AKI.cd[cortex_genes$Genes[1:30],]

sample_annotations <- data.frame(
  Timepoint = factor(c(rep("0h", ncol(AKI_sham.cortex.cd)), rep("0h", ncol(AKI_sham.interface.cd)),
                       rep("0h", ncol(AKI_sham.medulla.cd)), rep("4h", ncol(AKI_4h.cortex.cd)),
                       rep("4h", ncol(AKI_4h.interface.cd)), rep("4h", ncol(AKI_4h.medulla.cd)),
                       rep("12h", ncol(AKI_12h.cortex.cd)), rep("12h", ncol(AKI_12h.interface.cd)),
                       rep("12h", ncol(AKI_12h.medulla.cd)), rep("48h", ncol(AKI_2d.cortex.cd)),
                       rep("48h", ncol(AKI_2d.interface.cd)), rep("48h", ncol(AKI_2d.medulla.cd)) ), 
                     levels = c("0h", "4h", "12h", "48h")) ,  # Ensures correct order
  
  TissueType = factor(c(rep("Cortex", ncol(AKI_sham.cortex.cd)), rep("Interface", ncol(AKI_sham.interface.cd)),
                        rep("Medulla", ncol(AKI_sham.medulla.cd)), rep("Cortex", ncol(AKI_4h.cortex.cd)),
                        rep("Interface", ncol(AKI_4h.interface.cd)), rep("Medulla", ncol(AKI_4h.medulla.cd)),
                        rep("Cortex", ncol(AKI_12h.cortex.cd)), rep("Interface", ncol(AKI_12h.interface.cd)),
                        rep("Medulla", ncol(AKI_12h.medulla.cd)), rep("Cortex", ncol(AKI_2d.cortex.cd)),
                        rep("Interface", ncol(AKI_2d.interface.cd)), rep("Medulla", ncol(AKI_2d.medulla.cd)) ), 
                      levels = c("Cortex", "Interface", "Medulla"))  # Example categories
)
rownames(sample_annotations)<-colnames(AKI.cortex.cd)

#Scaling the matrix
mat <- as.matrix(AKI.cortex.cd)
scaled_mat <- t(scale(t(mat)))

overall_min<-min(scaled_mat)
overall_max <- quantile(scaled_mat, 0.95)

custom_colors<-colorRampPalette(c("black", "blue", "purple", "yellow"))(100)
breaks <- seq(overall_min, overall_max, length.out = 101)

# Define colors for annotations
ann_colors <- list(
  Timepoint = c("0h" = "lavender", "4h" = "mediumpurple", "12h" = "pink", "48h" = "magenta"),
  TissueType = c("Cortex" = "red", "Interface" = "blue", "Medulla" = "green")
)

# Draw heatmap with annotations
pdf("Figures/Figure1S/pdfs/AKI_Cortex_Marker_Genes_Heatmap.pdf", height=6, width=8)
pheatmap(scaled_mat, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, show_rownames = TRUE,
         color = custom_colors, breaks= breaks,
         annotation_col = sample_annotations,  # Add color-coded sample info
         annotation_colors = ann_colors,  # Apply custom colors
         main = "Top Gene Markers for Cortex (AKI)"
)
dev.off()

#Heatmap for AKI (Top Interface Genes)
#Note: interface is same as Outer Medulla as described in the main article section.
AKI.interface.cd<-AKI.cd[interface_genes$Genes[1:30],]

sample_annotations <- data.frame(
  Timepoint = factor(c(rep("0h", ncol(AKI_sham.cortex.cd)), rep("0h", ncol(AKI_sham.interface.cd)),
                       rep("0h", ncol(AKI_sham.medulla.cd)), rep("4h", ncol(AKI_4h.cortex.cd)),
                       rep("4h", ncol(AKI_4h.interface.cd)), rep("4h", ncol(AKI_4h.medulla.cd)),
                       rep("12h", ncol(AKI_12h.cortex.cd)), rep("12h", ncol(AKI_12h.interface.cd)),
                       rep("12h", ncol(AKI_12h.medulla.cd)), rep("48h", ncol(AKI_2d.cortex.cd)),
                       rep("48h", ncol(AKI_2d.interface.cd)), rep("48h", ncol(AKI_2d.medulla.cd)) ), 
                     levels = c("0h", "4h", "12h", "48h")) ,  # Ensures correct order
  
  TissueType = factor(c(rep("Cortex", ncol(AKI_sham.cortex.cd)), rep("Interface", ncol(AKI_sham.interface.cd)),
                        rep("Medulla", ncol(AKI_sham.medulla.cd)), rep("Cortex", ncol(AKI_4h.cortex.cd)),
                        rep("Interface", ncol(AKI_4h.interface.cd)), rep("Medulla", ncol(AKI_4h.medulla.cd)),
                        rep("Cortex", ncol(AKI_12h.cortex.cd)), rep("Interface", ncol(AKI_12h.interface.cd)),
                        rep("Medulla", ncol(AKI_12h.medulla.cd)), rep("Cortex", ncol(AKI_2d.cortex.cd)),
                        rep("Interface", ncol(AKI_2d.interface.cd)), rep("Medulla", ncol(AKI_2d.medulla.cd)) ), 
                      levels = c("Cortex", "Interface", "Medulla"))  # Example categories
)
rownames(sample_annotations)<-colnames(AKI.interface.cd)

#Scaling the matrix
mat <- as.matrix(AKI.interface.cd)
scaled_mat <- t(scale(t(mat)))

overall_min<-min(scaled_mat)
overall_max <- quantile(scaled_mat, 0.95)

custom_colors<-colorRampPalette(c("black", "blue", "purple", "yellow"))(100)
breaks <- seq(overall_min, overall_max, length.out = 101)

# Define colors for annotations
ann_colors <- list(
  Timepoint = c("0h" = "lavender", "4h" = "mediumpurple", "12h" = "pink", "48h" = "magenta"),
  TissueType = c("Cortex" = "red", "Interface" = "blue", "Medulla" = "green")
)

# Draw heatmap with annotations
pdf("Figures/Figure1S/pdfs/AKI_Interface_Marker_Genes_Heatmap.pdf", height=6, width=8)
pheatmap(scaled_mat, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, show_rownames = TRUE,
         color = custom_colors, breaks= breaks,
         annotation_col = sample_annotations,  # Add color-coded sample info
         annotation_colors = ann_colors,  # Apply custom colors
         main = "Top Gene Markers for Interface (AKI)"
)
dev.off()

#Heatmap for AKI (Top medulla genes)
#Note: medulla is same as Inner Medulla as described in the main article section.
AKI.medulla.cd<-AKI.cd[medulla_genes$Genes[1:30],]

sample_annotations <- data.frame(
  Timepoint = factor(c(rep("0h", ncol(AKI_sham.cortex.cd)), rep("0h", ncol(AKI_sham.interface.cd)),
                       rep("0h", ncol(AKI_sham.medulla.cd)), rep("4h", ncol(AKI_4h.cortex.cd)),
                       rep("4h", ncol(AKI_4h.interface.cd)), rep("4h", ncol(AKI_4h.medulla.cd)),
                       rep("12h", ncol(AKI_12h.cortex.cd)), rep("12h", ncol(AKI_12h.interface.cd)),
                       rep("12h", ncol(AKI_12h.medulla.cd)), rep("48h", ncol(AKI_2d.cortex.cd)),
                       rep("48h", ncol(AKI_2d.interface.cd)), rep("48h", ncol(AKI_2d.medulla.cd)) ), 
                     levels = c("0h", "4h", "12h", "48h")) ,  # Ensures correct order
  
  TissueType = factor(c(rep("Cortex", ncol(AKI_sham.cortex.cd)), rep("Interface", ncol(AKI_sham.interface.cd)),
                        rep("Medulla", ncol(AKI_sham.medulla.cd)), rep("Cortex", ncol(AKI_4h.cortex.cd)),
                        rep("Interface", ncol(AKI_4h.interface.cd)), rep("Medulla", ncol(AKI_4h.medulla.cd)),
                        rep("Cortex", ncol(AKI_12h.cortex.cd)), rep("Interface", ncol(AKI_12h.interface.cd)),
                        rep("Medulla", ncol(AKI_12h.medulla.cd)), rep("Cortex", ncol(AKI_2d.cortex.cd)),
                        rep("Interface", ncol(AKI_2d.interface.cd)), rep("Medulla", ncol(AKI_2d.medulla.cd)) ), 
                      levels = c("Cortex", "Interface", "Medulla"))  # Example categories
)
rownames(sample_annotations)<-colnames(AKI.medulla.cd)

#Scaling the matrix
mat <- as.matrix(AKI.medulla.cd)
scaled_mat <- t(scale(t(mat)))

overall_min<-min(scaled_mat)
overall_max <- quantile(scaled_mat, 0.95)

custom_colors<-colorRampPalette(c("black", "blue", "purple", "yellow"))(100)
breaks <- seq(overall_min, overall_max, length.out = 101)


# Define colors for annotations
ann_colors <- list(
  Timepoint = c("0h" = "lavender", "4h" = "mediumpurple", "12h" = "pink", "48h" = "magenta"),
  TissueType = c("Cortex" = "red", "Interface" = "blue", "Medulla" = "green")
)

# Draw heatmap with annotations
pdf("Figures/Figure1S/pdfs/AKI_Medulla_Marker_Genes_Heatmap.pdf", height=6, width=8)
pheatmap(scaled_mat, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, show_rownames = TRUE,
         color = custom_colors, breaks= breaks,
         annotation_col = sample_annotations,  # Add color-coded sample info
         annotation_colors = ann_colors,  # Apply custom colors
         main = "Top Gene Markers for Medulla (AKI)"
)
dev.off()

##################### CTRL DATASETS ############################################
#1.All_Cortex_Spots
ctrl1.cortex<-names(cortex[grep("^CTRL1", names(cortex))] )
ctrl2.cortex<-names(cortex[grep("^CTRL2", names(cortex))] )
ctrl3.cortex<-names(cortex[grep("^CTRL3", names(cortex))] )
ctrl4.cortex<-names(cortex[grep("^CTRL4", names(cortex))] )

#2.All Interface spots
ctrl1.interface<-names(interface[grep("^CTRL1", names(interface))] )
ctrl2.interface<-names(interface[grep("^CTRL2", names(interface))] )
ctrl3.interface<-names(interface[grep("^CTRL3", names(interface))] )
ctrl4.interface<-names(interface[grep("^CTRL4", names(interface))] )

#3.All medullary spots
ctrl1.medulla<-names(medulla[grep("^CTRL1", names(medulla))] )
ctrl2.medulla<-names(medulla[grep("CTRL2", names(medulla))] )
ctrl3.medulla<-names(medulla[grep("^CTRL3", names(medulla))] )
ctrl4.medulla<-names(medulla[grep("^CTRL4", names(medulla))] )

### COUNT NORMALIZATION OF DATA ################################################
ctrl1$mat_notlog <- MERINGUE::normalizeCounts(ctrl1$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
ctrl2$mat_notlog <- MERINGUE::normalizeCounts(ctrl2$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
ctrl3$mat_notlog <- MERINGUE::normalizeCounts(ctrl3$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
ctrl4$mat_notlog <- MERINGUE::normalizeCounts(ctrl4$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
#AKI_6w$mat_notlog <- MERINGUE::normalizeCounts(AKI_6w$gexp[unique(rownames(AKI_6w$gexp)),], log=FALSE)

ctrl1.cortex.cd<-ctrl1$mat_notlog[,ctrl1.cortex]
ctrl2.cortex.cd<-ctrl2$mat_notlog[,ctrl2.cortex]
ctrl3.cortex.cd<-ctrl3$mat_notlog[,ctrl3.cortex]
ctrl4.cortex.cd<-ctrl4$mat_notlog[,ctrl4.cortex]

ctrl1.interface.cd<-ctrl1$mat_notlog[,ctrl1.interface]
ctrl2.interface.cd<-ctrl2$mat_notlog[,ctrl2.interface]
ctrl3.interface.cd<-ctrl3$mat_notlog[,ctrl3.interface]
ctrl4.interface.cd<-ctrl4$mat_notlog[,ctrl4.interface]
#CIS_6w.cortex.cd<-CIS_48h$mat_notlog[,CIS_cortex.48h]

ctrl1.medulla.cd<-ctrl1$mat_notlog[,ctrl1.medulla]
ctrl2.medulla.cd<-ctrl2$mat_notlog[,ctrl2.medulla]
ctrl3.medulla.cd<-ctrl3$mat_notlog[,ctrl3.medulla]
ctrl4.medulla.cd<-ctrl4$mat_notlog[,ctrl4.medulla]

#HEATMAP for CTRL:
CTRL.cd<-cbind(ctrl1.cortex.cd, ctrl1.interface.cd, ctrl1.medulla.cd, 
              ctrl2.cortex.cd, ctrl2.interface.cd, ctrl2.medulla.cd,
              ctrl3.cortex.cd, ctrl3.interface.cd, ctrl3.medulla.cd, 
              ctrl4.cortex.cd, ctrl4.interface.cd, ctrl4.medulla.cd)

#HEATMAP (TOP CORTEX GENES)
CTRL.cortex.cd<-CTRL.cd[cortex_genes$Genes[1:30],]

sample_annotations <- data.frame(
  Specimen = factor(c(rep("Kidney1", ncol(ctrl1.cortex.cd)), rep("Kidney1", ncol(ctrl1.interface.cd)),
                       rep("Kidney1", ncol(ctrl1.medulla.cd)), rep("Kidney2", ncol(ctrl2.cortex.cd)),
                       rep("Kidney2", ncol(ctrl2.interface.cd)), rep("Kidney2", ncol(ctrl2.medulla.cd)),
                       rep("Kidney3", ncol(ctrl3.cortex.cd)), rep("Kidney3", ncol(ctrl3.interface.cd)),
                       rep("Kidney3", ncol(ctrl3.medulla.cd)), rep("Kidney4", ncol(ctrl4.cortex.cd)),
                       rep("Kidney4", ncol(ctrl4.interface.cd)), rep("Kidney4", ncol(ctrl4.medulla.cd)) ), 
                     levels = c("Kidney1", "Kidney2", "Kidney3", "Kidney4")) ,  # Ensures correct order
  
  TissueType = factor(c(rep("Cortex", ncol(ctrl1.cortex.cd)), rep("Interface", ncol(ctrl1.interface.cd)),
                        rep("Medulla", ncol(ctrl1.medulla.cd)), rep("Cortex", ncol(ctrl2.cortex.cd)),
                        rep("Interface", ncol(ctrl2.interface.cd)), rep("Medulla", ncol(ctrl2.medulla.cd)),
                        rep("Cortex", ncol(ctrl3.cortex.cd)), rep("Interface", ncol(ctrl3.interface.cd)),
                        rep("Medulla", ncol(ctrl3.medulla.cd)), rep("Cortex", ncol(ctrl4.cortex.cd)),
                        rep("Interface", ncol(ctrl4.interface.cd)), rep("Medulla", ncol(ctrl4.medulla.cd)) ), 
                      levels = c("Cortex", "Interface", "Medulla"))  # Example categories
)
rownames(sample_annotations)<-colnames(CTRL.cortex.cd)

#Scaling the matrix
mat <- as.matrix(CTRL.cortex.cd)
scaled_mat <- t(scale(t(mat)))

overall_min<-min(scaled_mat)
overall_max <- quantile(scaled_mat, 0.95)

custom_colors<-colorRampPalette(c("black", "blue", "purple", "yellow"))(100)
breaks <- seq(overall_min, overall_max, length.out = 101)

# Define colors for annotations
ann_colors <- list(
  Specimen = c("Kidney1" = "lavender", "Kidney2" = "mediumpurple", "Kidney3" = "pink", "Kidney4" = "magenta"),
  TissueType = c("Cortex" = "red", "Interface" = "blue", "Medulla" = "green")
)

# Draw heatmap with annotations
pdf("Figures/Figure1S/pdfs/CTRL_Cortex_Marker_Genes_Heatmap.pdf", height=6, width=8)
pheatmap(scaled_mat, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, show_rownames = TRUE,
         color = custom_colors, breaks= breaks,
         annotation_col = sample_annotations,  # Add color-coded sample info
         annotation_colors = ann_colors,  # Apply custom colors
         main = "Top Gene Markers for Cortex (CTRL)"
)
dev.off()

#HEATMAP (TOP INTERFACE GENES)
#Note: interface is same as Outer Medulla as described in the main article section.
CTRL.interface.cd<-CTRL.cd[interface_genes$Genes[1:30],]

sample_annotations <- data.frame(
  Specimen = factor(c(rep("Kidney1", ncol(ctrl1.cortex.cd)), rep("Kidney1", ncol(ctrl1.interface.cd)),
                      rep("Kidney1", ncol(ctrl1.medulla.cd)), rep("Kidney2", ncol(ctrl2.cortex.cd)),
                      rep("Kidney2", ncol(ctrl2.interface.cd)), rep("Kidney2", ncol(ctrl2.medulla.cd)),
                      rep("Kidney3", ncol(ctrl3.cortex.cd)), rep("Kidney3", ncol(ctrl3.interface.cd)),
                      rep("Kidney3", ncol(ctrl3.medulla.cd)), rep("Kidney4", ncol(ctrl4.cortex.cd)),
                      rep("Kidney4", ncol(ctrl4.interface.cd)), rep("Kidney4", ncol(ctrl4.medulla.cd)) ), 
                    levels = c("Kidney1", "Kidney2", "Kidney3", "Kidney4")) ,  # Ensures correct order
  
  TissueType = factor(c(rep("Cortex", ncol(ctrl1.cortex.cd)), rep("Interface", ncol(ctrl1.interface.cd)),
                        rep("Medulla", ncol(ctrl1.medulla.cd)), rep("Cortex", ncol(ctrl2.cortex.cd)),
                        rep("Interface", ncol(ctrl2.interface.cd)), rep("Medulla", ncol(ctrl2.medulla.cd)),
                        rep("Cortex", ncol(ctrl3.cortex.cd)), rep("Interface", ncol(ctrl3.interface.cd)),
                        rep("Medulla", ncol(ctrl3.medulla.cd)), rep("Cortex", ncol(ctrl4.cortex.cd)),
                        rep("Interface", ncol(ctrl4.interface.cd)), rep("Medulla", ncol(ctrl4.medulla.cd)) ), 
                      levels = c("Cortex", "Interface", "Medulla"))  # Example categories
)
rownames(sample_annotations)<-colnames(CTRL.interface.cd)

#Scaling the matrix
mat <- as.matrix(CTRL.interface.cd)
scaled_mat <- t(scale(t(mat)))

overall_min<-min(scaled_mat)
overall_max <- quantile(scaled_mat, 0.95)

custom_colors<-colorRampPalette(c("black", "blue", "purple", "yellow"))(100)
breaks <- seq(overall_min, overall_max, length.out = 101)

# Define colors for annotations
ann_colors <- list(
  Specimen = c("Kidney1" = "lavender", "Kidney2" = "mediumpurple", "Kidney3" = "pink", "Kidney4" = "magenta"),
  TissueType = c("Cortex" = "red", "Interface" = "blue", "Medulla" = "green")
)

# Draw heatmap with annotations
pdf("Figures/Figure1S/pdfs/CTRL_Interface_Marker_Genes_Heatmap.pdf", height=6, width=8)
pheatmap(scaled_mat, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, show_rownames = TRUE,
         color = custom_colors, breaks= breaks,
         annotation_col = sample_annotations,  # Add color-coded sample info
         annotation_colors = ann_colors,  # Apply custom colors
         main = "Top Gene Markers for Interface (CTRL)"
)
dev.off()

#Heatmap for CTRL (Top medulla genes)
#Note: medulla is same as Inner Medulla as described in the main article section.
CTRL.medulla.cd<-CTRL.cd[medulla_genes$Genes[1:30],]

sample_annotations <- data.frame(
  Specimen = factor(c(rep("Kidney1", ncol(ctrl1.cortex.cd)), rep("Kidney1", ncol(ctrl1.interface.cd)),
                      rep("Kidney1", ncol(ctrl1.medulla.cd)), rep("Kidney2", ncol(ctrl2.cortex.cd)),
                      rep("Kidney2", ncol(ctrl2.interface.cd)), rep("Kidney2", ncol(ctrl2.medulla.cd)),
                      rep("Kidney3", ncol(ctrl3.cortex.cd)), rep("Kidney3", ncol(ctrl3.interface.cd)),
                      rep("Kidney3", ncol(ctrl3.medulla.cd)), rep("Kidney4", ncol(ctrl4.cortex.cd)),
                      rep("Kidney4", ncol(ctrl4.interface.cd)), rep("Kidney4", ncol(ctrl4.medulla.cd)) ), 
                    levels = c("Kidney1", "Kidney2", "Kidney3", "Kidney4")) ,  # Ensures correct order
  
  TissueType = factor(c(rep("Cortex", ncol(ctrl1.cortex.cd)), rep("Interface", ncol(ctrl1.interface.cd)),
                        rep("Medulla", ncol(ctrl1.medulla.cd)), rep("Cortex", ncol(ctrl2.cortex.cd)),
                        rep("Interface", ncol(ctrl2.interface.cd)), rep("Medulla", ncol(ctrl2.medulla.cd)),
                        rep("Cortex", ncol(ctrl3.cortex.cd)), rep("Interface", ncol(ctrl3.interface.cd)),
                        rep("Medulla", ncol(ctrl3.medulla.cd)), rep("Cortex", ncol(ctrl4.cortex.cd)),
                        rep("Interface", ncol(ctrl4.interface.cd)), rep("Medulla", ncol(ctrl4.medulla.cd)) ), 
                      levels = c("Cortex", "Interface", "Medulla"))  # Example categories
)
rownames(sample_annotations)<-colnames(CTRL.medulla.cd)

#Scaling the matrix
mat <- as.matrix(CTRL.medulla.cd)
scaled_mat <- t(scale(t(mat)))

overall_min<-min(scaled_mat)
overall_max <- quantile(scaled_mat, 0.95)

custom_colors<-colorRampPalette(c("black", "blue", "purple", "yellow"))(100)
breaks <- seq(overall_min, overall_max, length.out = 101)

# Define colors for annotations
ann_colors <- list(
  Specimen = c("Kidney1" = "lavender", "Kidney2" = "mediumpurple", "Kidney3" = "pink", "Kidney4" = "magenta"),
  TissueType = c("Cortex" = "red", "Interface" = "blue", "Medulla" = "green")
)

# Draw heatmap with annotations
pdf("Figures/Figure1S/pdfs/CTRL_Medulla_Marker_Genes_Heatmap.pdf", height=6, width=8)
pheatmap(scaled_mat, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, show_rownames = TRUE,
         color = custom_colors, breaks= breaks,
         annotation_col = sample_annotations,  # Add color-coded sample info
         annotation_colors = ann_colors,  # Apply custom colors
         main = "Top Gene Markers for Medulla (CTRL)"
)
dev.off()

##################### CTRL DATASETS ############################################
#1.All_Cortex_Spots
irl1.cortex<-names(cortex[grep("^IR1", names(cortex))] )
irl2.cortex<-names(cortex[grep("^IR2", names(cortex))] )
irl3.cortex<-names(cortex[grep("^IR3", names(cortex))] )
irl4.cortex<-names(cortex[grep("^IR4", names(cortex))] )

#2.All Interface spots
irl1.interface<-names(interface[grep("^IR1", names(interface))] )
irl2.interface<-names(interface[grep("^IR2", names(interface))] )
irl3.interface<-names(interface[grep("^IR3", names(interface))] )
irl4.interface<-names(interface[grep("^IR4", names(interface))] )

#3.All medullary spots
irl1.medulla<-names(medulla[grep("^IR1", names(medulla))] )
irl2.medulla<-names(medulla[grep("IR2", names(medulla))] )
irl3.medulla<-names(medulla[grep("^IR3", names(medulla))] )
irl4.medulla<-names(medulla[grep("^IR4", names(medulla))] )

### COUNT NORMALIZATION OF DATA ################################################
irl1$mat_notlog <- MERINGUE::normalizeCounts(irl1$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
irl2$mat_notlog <- MERINGUE::normalizeCounts(irl2$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
irl3$mat_notlog <- MERINGUE::normalizeCounts(irl3$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
irl4$mat_notlog <- MERINGUE::normalizeCounts(irl4$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)

irl1.cortex.cd<-irl1$mat_notlog[,irl1.cortex]
irl2.cortex.cd<-irl2$mat_notlog[,irl2.cortex]
irl3.cortex.cd<-irl3$mat_notlog[,irl3.cortex]
irl4.cortex.cd<-irl4$mat_notlog[,irl4.cortex]

irl1.interface.cd<-irl1$mat_notlog[,irl1.interface]
irl2.interface.cd<-irl2$mat_notlog[,irl2.interface]
irl3.interface.cd<-irl3$mat_notlog[,irl3.interface]
irl4.interface.cd<-irl4$mat_notlog[,irl4.interface]
#CIS_6w.cortex.cd<-CIS_48h$mat_notlog[,CIS_cortex.48h]

irl1.medulla.cd<-irl1$mat_notlog[,irl1.medulla]
irl2.medulla.cd<-irl2$mat_notlog[,irl2.medulla]
irl3.medulla.cd<-irl3$mat_notlog[,irl3.medulla]
irl4.medulla.cd<-irl4$mat_notlog[,irl4.medulla]


#Heatmap for IRL:
IRL.cd<-cbind(irl1.cortex.cd, irl1.interface.cd, irl1.medulla.cd, 
               irl2.cortex.cd, irl2.interface.cd, irl2.medulla.cd,
               irl3.cortex.cd, irl3.interface.cd, irl3.medulla.cd, 
               irl4.cortex.cd, irl4.interface.cd, irl4.medulla.cd)

#Heatmap for IRL (Top cortex genes)
IRL.cortex.cd<-IRL.cd[cortex_genes$Genes[1:30],]

sample_annotations <- data.frame(
  Specimen = factor(c(rep("Kidney1", ncol(irl1.cortex.cd)), rep("Kidney1", ncol(irl1.interface.cd)),
                      rep("Kidney1", ncol(irl1.medulla.cd)), rep("Kidney2", ncol(irl2.cortex.cd)),
                      rep("Kidney2", ncol(irl2.interface.cd)), rep("Kidney2", ncol(irl2.medulla.cd)),
                      rep("Kidney3", ncol(irl3.cortex.cd)), rep("Kidney3", ncol(irl3.interface.cd)),
                      rep("Kidney3", ncol(irl3.medulla.cd)), rep("Kidney4", ncol(irl4.cortex.cd)),
                      rep("Kidney4", ncol(irl4.interface.cd)), rep("Kidney4", ncol(irl4.medulla.cd)) ), 
                    levels = c("Kidney1", "Kidney2", "Kidney3", "Kidney4")) ,  # Ensures correct order
  
  TissueType = factor(c(rep("Cortex", ncol(irl1.cortex.cd)), rep("Interface", ncol(irl1.interface.cd)),
                        rep("Medulla", ncol(irl1.medulla.cd)), rep("Cortex", ncol(irl2.cortex.cd)),
                        rep("Interface", ncol(irl2.interface.cd)), rep("Medulla", ncol(irl2.medulla.cd)),
                        rep("Cortex", ncol(irl3.cortex.cd)), rep("Interface", ncol(irl3.interface.cd)),
                        rep("Medulla", ncol(irl3.medulla.cd)), rep("Cortex", ncol(irl4.cortex.cd)),
                        rep("Interface", ncol(irl4.interface.cd)), rep("Medulla", ncol(irl4.medulla.cd)) ), 
                      levels = c("Cortex", "Interface", "Medulla"))  # Example categories
)
rownames(sample_annotations)<-colnames(IRL.cortex.cd)


#Scaling the matrix
mat <- as.matrix(IRL.cortex.cd)
scaled_mat <- t(scale(t(mat)))

overall_min<-min(scaled_mat)
overall_max <- quantile(scaled_mat, 0.95)

custom_colors<-colorRampPalette(c("black", "blue", "purple", "yellow"))(100)
breaks <- seq(overall_min, overall_max, length.out = 101)

# Define colors for annotations
ann_colors <- list(
  Specimen = c("Kidney1" = "lavender", "Kidney2" = "mediumpurple", "Kidney3" = "pink", "Kidney4" = "magenta"),
  TissueType = c("Cortex" = "red", "Interface" = "blue", "Medulla" = "green")
)

# Draw heatmap with annotations
pdf("Figures/Figure1S/pdfs/IRL_Cortex_Marker_Genes_Heatmap.pdf", height=6, width=8)
pheatmap(scaled_mat, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, show_rownames = TRUE,
         color = custom_colors, breaks= breaks,
         annotation_col = sample_annotations,  # Add color-coded sample info
         annotation_colors = ann_colors,  # Apply custom colors
         main = "Top Gene Markers for Cortex (IRL)"
)
dev.off()

#Heatmap for IRL (Top Interface Markers)
#Note: interface is same as Outer Medulla as described in the main article section.
IRL.interface.cd<-IRL.cd[interface_genes$Genes[1:30],]

sample_annotations <- data.frame(
  Specimen = factor(c(rep("Kidney1", ncol(irl1.cortex.cd)), rep("Kidney1", ncol(irl1.interface.cd)),
                      rep("Kidney1", ncol(irl1.medulla.cd)), rep("Kidney2", ncol(irl2.cortex.cd)),
                      rep("Kidney2", ncol(irl2.interface.cd)), rep("Kidney2", ncol(irl2.medulla.cd)),
                      rep("Kidney3", ncol(irl3.cortex.cd)), rep("Kidney3", ncol(irl3.interface.cd)),
                      rep("Kidney3", ncol(irl3.medulla.cd)), rep("Kidney4", ncol(irl4.cortex.cd)),
                      rep("Kidney4", ncol(irl4.interface.cd)), rep("Kidney4", ncol(irl4.medulla.cd)) ), 
                    levels = c("Kidney1", "Kidney2", "Kidney3", "Kidney4")) ,  # Ensures correct order
  
  TissueType = factor(c(rep("Cortex", ncol(irl1.cortex.cd)), rep("Interface", ncol(irl1.interface.cd)),
                        rep("Medulla", ncol(irl1.medulla.cd)), rep("Cortex", ncol(irl2.cortex.cd)),
                        rep("Interface", ncol(irl2.interface.cd)), rep("Medulla", ncol(irl2.medulla.cd)),
                        rep("Cortex", ncol(irl3.cortex.cd)), rep("Interface", ncol(irl3.interface.cd)),
                        rep("Medulla", ncol(irl3.medulla.cd)), rep("Cortex", ncol(irl4.cortex.cd)),
                        rep("Interface", ncol(irl4.interface.cd)), rep("Medulla", ncol(irl4.medulla.cd)) ), 
                      levels = c("Cortex", "Interface", "Medulla"))  # Example categories
)
rownames(sample_annotations)<-colnames(IRL.interface.cd)


#Scaling the matrix
mat <- as.matrix(IRL.interface.cd)
scaled_mat <- t(scale(t(mat)))

overall_min<-min(scaled_mat)
overall_max <- quantile(scaled_mat, 0.95)

custom_colors<-colorRampPalette(c("black", "blue", "purple", "yellow"))(100)
breaks <- seq(overall_min, overall_max, length.out = 101)

# Define colors for annotations
ann_colors <- list(
  Specimen = c("Kidney1" = "lavender", "Kidney2" = "mediumpurple", "Kidney3" = "pink", "Kidney4" = "magenta"),
  TissueType = c("Cortex" = "red", "Interface" = "blue", "Medulla" = "green")
)

# Draw heatmap with annotations
pdf("Figures/Figure1S/pdfs/IRL_Interface_Marker_Genes_Heatmap.pdf", height=6, width=8)
pheatmap(scaled_mat, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, show_rownames = TRUE,
         color = custom_colors, breaks= breaks,
         annotation_col = sample_annotations,  # Add color-coded sample info
         annotation_colors = ann_colors,  # Apply custom colors
         main = "Top Gene Markers for Interface (IRL)"
)
dev.off()

#Heatmap for IRL (Top Medulla Genes)
#Note: medulla is same as Inner Medulla as described in the main article section.
IRL.medulla.cd<-IRL.cd[medulla_genes$Genes[1:30],]

sample_annotations <- data.frame(
  Specimen = factor(c(rep("Kidney1", ncol(irl1.cortex.cd)), rep("Kidney1", ncol(irl1.interface.cd)),
                      rep("Kidney1", ncol(irl1.medulla.cd)), rep("Kidney2", ncol(irl2.cortex.cd)),
                      rep("Kidney2", ncol(irl2.interface.cd)), rep("Kidney2", ncol(irl2.medulla.cd)),
                      rep("Kidney3", ncol(irl3.cortex.cd)), rep("Kidney3", ncol(irl3.interface.cd)),
                      rep("Kidney3", ncol(irl3.medulla.cd)), rep("Kidney4", ncol(irl4.cortex.cd)),
                      rep("Kidney4", ncol(irl4.interface.cd)), rep("Kidney4", ncol(irl4.medulla.cd)) ), 
                    levels = c("Kidney1", "Kidney2", "Kidney3", "Kidney4")) ,  # Ensures correct order
  
  TissueType = factor(c(rep("Cortex", ncol(irl1.cortex.cd)), rep("Interface", ncol(irl1.interface.cd)),
                        rep("Medulla", ncol(irl1.medulla.cd)), rep("Cortex", ncol(irl2.cortex.cd)),
                        rep("Interface", ncol(irl2.interface.cd)), rep("Medulla", ncol(irl2.medulla.cd)),
                        rep("Cortex", ncol(irl3.cortex.cd)), rep("Interface", ncol(irl3.interface.cd)),
                        rep("Medulla", ncol(irl3.medulla.cd)), rep("Cortex", ncol(irl4.cortex.cd)),
                        rep("Interface", ncol(irl4.interface.cd)), rep("Medulla", ncol(irl4.medulla.cd)) ), 
                      levels = c("Cortex", "Interface", "Medulla"))  # Example categories
)
rownames(sample_annotations)<-colnames(IRL.medulla.cd)


#Scaling the matrix
mat <- as.matrix(IRL.medulla.cd)
scaled_mat <- t(scale(t(mat)))

overall_min<-min(scaled_mat)
overall_max <- quantile(scaled_mat, 0.95)

custom_colors<-colorRampPalette(c("black", "blue", "purple", "yellow"))(100)
breaks <- seq(overall_min, overall_max, length.out = 101)

# Define colors for annotations
ann_colors <- list(
  Specimen = c("Kidney1" = "lavender", "Kidney2" = "mediumpurple", "Kidney3" = "pink", "Kidney4" = "magenta"),
  TissueType = c("Cortex" = "red", "Interface" = "blue", "Medulla" = "green")
)

# Draw heatmap with annotations
pdf("Figures/Figure1S/pdfs/IRL_Medulla_Marker_Genes_Heatmap.pdf", height=6, width=8)
pheatmap(scaled_mat, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, show_rownames = TRUE,
         color = custom_colors, breaks= breaks,
         annotation_col = sample_annotations,  # Add color-coded sample info
         annotation_colors = ann_colors,  # Apply custom colors
         main = "Top Gene Markers for Medulla (IRL)"
)
dev.off()