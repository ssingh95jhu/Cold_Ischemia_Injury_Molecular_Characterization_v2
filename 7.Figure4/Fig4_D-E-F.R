#Clean environment
rm(list=ls(all.names=TRUE))
gc()

dir<-"/Users/srujansingh/Library/CloudStorage/OneDrive-SharedLibraries-JohnsHopkins/Jean Fan - 2024_spatial_AKI/Revision2_GenomBiol"
setwd(dir)

########################## ADDING LIBRARY ######################################
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(dplyr)

#################### CREATE WORKBOOKS ##########################################
wb1<-createWorkbook()

######################## 1. PATHWAYS VISUALIZATION ################################
#######Extracting all saved pathways
### 1.CORTEX ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cis_cortex_HALL<-read.xlsx("Codes2_Tables/CIS_HALLMARK_Pathays_Consistency_R2_Percentile95.xlsx", sheet = "Cortex")
rownames(cis_cortex_HALL)<-cis_cortex_HALL$ID
head(cis_cortex_HALL)

aki_cortex_HALL<-read.xlsx("Codes2_Tables/AKI_HALLMARK_Pathways_Consistency_R2_Percentile95.xlsx", sheet = "Cortex")
rownames(aki_cortex_HALL)<-aki_cortex_HALL$ID
head(aki_cortex_HALL)

cis_cortex_KEGG<-read.xlsx("Codes2_Tables/CIS_KEGG_Pathays_Consistency_R2_Percentile95.xlsx", sheet = "Cortex")
rownames(cis_cortex_KEGG)<-cis_cortex_KEGG$ID
head(cis_cortex_KEGG)

aki_cortex_KEGG<-read.xlsx("Codes2_Tables/AKI_KEGG_Pathways_Consistency_R2_Percentile95.xlsx", sheet = "Cortex")
rownames(aki_cortex_KEGG)<-aki_cortex_KEGG$ID
head(aki_cortex_KEGG)

### 2.INTERFACE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cis_interface_HALL<-read.xlsx("Codes2_Tables/CIS_HALLMARK_Pathays_Consistency_R2_Percentile95.xlsx", sheet = "Outer_Medulla")
rownames(cis_interface_HALL)<-cis_interface_HALL$ID
head(cis_interface_HALL)

aki_interface_HALL<-read.xlsx("Codes2_Tables/AKI_HALLMARK_Pathways_Consistency_R2_Percentile95.xlsx", sheet = "Outer_Medulla")
rownames(aki_interface_HALL)<-aki_interface_HALL$ID
head(aki_interface_HALL)

cis_interface_KEGG<-read.xlsx("Codes2_Tables/CIS_KEGG_Pathays_Consistency_R2_Percentile95.xlsx", sheet = "Outer_Medulla")
rownames(cis_interface_KEGG)<-cis_interface_KEGG$ID
head(cis_interface_KEGG)

aki_interface_KEGG<-read.xlsx("Codes2_Tables/AKI_KEGG_Pathways_Consistency_R2_Percentile95.xlsx", sheet = "Outer_Medulla")
rownames(aki_interface_KEGG)<-aki_interface_KEGG$ID
head(aki_interface_KEGG)

### 3.MEDULLA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cis_medulla_HALL<-read.xlsx("Codes2_Tables/CIS_HALLMARK_Pathays_Consistency_R2_Percentile95.xlsx", sheet = "Inner_Medulla")
rownames(cis_medulla_HALL)<-cis_medulla_HALL$ID
head(cis_medulla_HALL)

aki_medulla_HALL<-read.xlsx("Codes2_Tables/AKI_HALLMARK_Pathways_Consistency_R2_Percentile95.xlsx", sheet = "Inner_Medulla")
rownames(aki_medulla_HALL)<-aki_medulla_HALL$ID
head(aki_medulla_HALL)

cis_medulla_KEGG<-read.xlsx("Codes2_Tables/CIS_KEGG_Pathays_Consistency_R2_Percentile95.xlsx", sheet = "Inner_Medulla")
rownames(cis_medulla_KEGG)<-cis_medulla_KEGG$ID
head(cis_medulla_KEGG)

aki_medulla_KEGG<-read.xlsx("Codes2_Tables/AKI_KEGG_Pathways_Consistency_R2_Percentile95.xlsx", sheet = "Inner_Medulla")
rownames(aki_medulla_KEGG)<-aki_medulla_KEGG$ID
head(aki_medulla_KEGG)

######## PLOTTING THE BAR PLOTS ################################################
## CIS CORTEX vs AKI CORTEX (KEGG AND HALLMARK) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
comm_KEGG<-intersect(rownames(cis_cortex_KEGG), rownames(aki_cortex_KEGG) )
df0<-data.frame(cis_cortex_KEGG[comm_KEGG,]$NES, aki_cortex_KEGG[comm_KEGG,]$NES)
df0$des<-cis_cortex_KEGG[comm_KEGG,]$Description
df0$Sdes<-gsub(" - Mus musculus.*", "", df0$des)
df0$Sdes<-paste0("KEGG_",df0$Sdes)
colnames(df0)<-c('CIS_cortex_NES','AKI_cortex_NES', 'CIS_AKI_Path', 'ShortPath')

comm_HALL<-intersect(rownames(cis_cortex_HALL), rownames(aki_cortex_HALL) )
df1<-data.frame(cis_cortex_HALL[comm_HALL,]$NES, aki_cortex_HALL[comm_HALL,]$NES)
df1$des<-cis_cortex_HALL[comm_HALL,]$Description
df1$Sdes<-gsub("HALLMARK_", "HM_", df1$des)
colnames(df1)<-c('CIS_cortex_NES','AKI_cortex_NES', 'CIS_AKI_Path', 'ShortPath')

df<-rbind(df0,df1)
head(df)

df_covary<-df[which((df$CIS_cortex_NES*df$AKI_cortex_NES)>0),]
head(df_covary)

addWorksheet(wb1,"Cortex_Covary")
writeData(wb1, "Cortex_Covary", df_covary)

CIS_Cor<-data.frame(Path=df_covary$ShortPath, CIS_Cortex_NES=df_covary$CIS_cortex_NES, Enrichment="positive")
#CIS_Cor_KEGG_plot<-rbind(head(CIS_Cor_KEGG), tail(CIS_Cor_KEGG))
CIS_Cor_up<-CIS_Cor[which(CIS_Cor$CIS_Cortex_NES>0),]
CIS_Cor_down<-CIS_Cor[which(CIS_Cor$CIS_Cortex_NES<0),]
CIS_Cor_plot<-rbind(CIS_Cor_up, CIS_Cor_down)
CIS_Cor_plot$Enrichment<-ifelse(CIS_Cor_plot$CIS_Cortex_NES<0, "negative", "positive")
CIS_Cor_plot$CIS_Cortex_NES<-ifelse(CIS_Cor_plot$CIS_Cortex_NES<0, -CIS_Cor_plot$CIS_Cortex_NES, CIS_Cor_plot$CIS_Cortex_NES)

AKI_Cor<-data.frame(Path=df_covary$ShortPath, AKI_Cortex_NES=df_covary$AKI_cortex_NES, Enrichment="positive")
# AKI_Cor_KEGG_plot<-rbind(head(AKI_Cor_KEGG), tail(AKI_Cor_KEGG))
AKI_Cor_up<-AKI_Cor[which(AKI_Cor$AKI_Cortex_NES>0),]
AKI_Cor_down<-AKI_Cor[which(AKI_Cor$AKI_Cortex_NES<0),]
AKI_Cor_plot<-rbind(AKI_Cor_up, AKI_Cor_down)
AKI_Cor_plot$Enrichment<-ifelse(AKI_Cor_plot$AKI_Cortex_NES<0, "negative", "positive")
AKI_Cor_plot$AKI_Cortex_NES<-ifelse(AKI_Cor_plot$AKI_Cortex_NES>0, -AKI_Cor_plot$AKI_Cortex_NES, AKI_Cor_plot$AKI_Cortex_NES)

colnames(CIS_Cor_plot)<-c('Pathway', 'NES', 'Enrichment')
colnames(AKI_Cor_plot)<-c('Pathway', 'NES', 'Enrichment')
CIS_AKI_Cortex_plot<-rbind(CIS_Cor_plot,AKI_Cor_plot)
CIS_AKI_Cortex_plot<-CIS_AKI_Cortex_plot[rev(order(CIS_AKI_Cortex_plot$NES)),]

file<-paste0(dir,"/Codes2_Figures/Figure4/pdfs/CIS_vs_AKI_Covary_NES_Cortex_Barplots.pdf")
pdf(file, width=7, height=1.5)
ggplot(CIS_AKI_Cortex_plot, aes(x = rev(factor(Pathway, levels=unique(Pathway))), y = NES , fill = Enrichment)) +
  #geom_bar(stat = "identity", position = "identity", width=0.7, alpha =0.5) + coord_flip() + 
  geom_bar(stat = "identity", position = "identity") + coord_flip() + 
  scale_x_discrete(labels = rev(unique(CIS_AKI_Cortex_plot$Path))) +
  #geom_text(aes(label = KEGG_Path, y=NES*1.8),  color = "black", size = 3, fontface="bold") +
  scale_y_continuous(limits=c(-7,7), breaks = c(seq(-7, 7, by = 2), 0)) +   
  labs(x="",y = "Normalized Enrichment Score", title = "CIS Cortex vs AKI Cortex") +
  scale_fill_manual(values = c("negative" = "deepskyblue", "positive" = "orangered")) + 
  geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.5) +
  #geom_vline(xintercept = 0.5, color = "black", linetype = "solid", size = 1) +
  theme_minimal() + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add box around plot
        #axis.text.x = element_text(size = 14, color="black"), # Increase font size of pathway names
        axis.text.y = element_text(size = 14, color="black"),
        axis.ticks = element_line(color = "black"),  # Add tick marks
        axis.ticks.length = unit(0.2, "cm"),  # Adjust tick mark length
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) 
dev.off()

#Disjoint Pathways ************************************************************
#[NO DISJOINT PATHWAYS]
df_disjoint<-df[which((df$CIS_cortex_NES*df$AKI_cortex_NES)<0),]
head(df_disjoint)

addWorksheet(wb1,"Cortex_Disjoint")
writeData(wb1, "Cortex_Disjoint", df_disjoint)

# 
# CIS_Cor<-data.frame(Path=df_disjoint$ShortPath, CIS_Cortex_NES=df_disjoint$CIS_cortex_NES, Enrichment="positive")
# #CIS_Cor_KEGG_plot<-rbind(head(CIS_Cor_KEGG), tail(CIS_Cor_KEGG))
# CIS_Cor_up<-CIS_Cor[which(CIS_Cor$CIS_Cortex_NES>0),]
# CIS_Cor_down<-CIS_Cor[which(CIS_Cor$CIS_Cortex_NES<0),]
# CIS_Cor_plot<-rbind(CIS_Cor_up, CIS_Cor_down)
# CIS_Cor_plot$Enrichment<-ifelse(CIS_Cor_plot$CIS_Cortex_NES<0, "negative", "positive")
# CIS_Cor_plot$CIS_Cortex_NES<-ifelse(CIS_Cor_plot$CIS_Cortex_NES<0, -CIS_Cor_plot$CIS_Cortex_NES, CIS_Cor_plot$CIS_Cortex_NES)
# 
# AKI_Cor<-data.frame(Path=df_disjoint$ShortPath, AKI_Cortex_NES=df_disjoint$AKI_cortex_NES, Enrichment="positive")
# #AKI_Cor_KEGG_plot<-rbind(head(AKI_Cor_KEGG), tail(AKI_Cor_KEGG))
# AKI_Cor_up<-AKI_Cor[which(AKI_Cor$AKI_Cortex_NES>0),]
# AKI_Cor_down<-AKI_Cor[which(AKI_Cor$AKI_Cortex_NES<0),]
# AKI_Cor_plot<-rbind(AKI_Cor_up, AKI_Cor_down)
# AKI_Cor_plot$Enrichment<-ifelse(AKI_Cor_plot$AKI_Cortex_NES<0, "negative", "positive")
# AKI_Cor_plot$AKI_Cortex_NES<-ifelse(AKI_Cor_plot$AKI_Cortex_NES>0, -AKI_Cor_plot$AKI_Cortex_NES, AKI_Cor_plot$AKI_Cortex_NES)
# 
# colnames(CIS_Cor_plot)<-c('Pathway', 'NES', 'Enrichment')
# colnames(AKI_Cor_plot)<-c('Pathway', 'NES', 'Enrichment')
# CIS_AKI_Cortex_plot<-rbind(CIS_Cor_plot,AKI_Cor_plot)
# #CIS_AKI_Cortex_KEGG_plot<-rbind(AKI_Cor_KEGG_plot,CIS_Cor_KEGG_plot)
# 
# #pdf("Figures/Figure4S/pdfs/CIS_vs_AKI_GSEA_NES_Cortex_Barplots.pdf", width=8, height=16)
# ggplot(CIS_AKI_Cortex_plot, aes(x = factor(Pathway, levels=unique(Pathway)), y = NES , fill = Enrichment)) +
#   #geom_bar(stat = "identity", position = "identity", width=0.7, alpha =0.5) + coord_flip() + 
#   geom_bar(stat = "identity", position = "identity") + coord_flip() + 
#   scale_x_discrete(labels = unique(CIS_AKI_Cortex_KEGG_plot$KEGG_Path)) +
#   #geom_text(aes(label = KEGG_Path, y=NES*1.8),  color = "black", size = 3, fontface="bold") +
#   scale_y_continuous(limits=c(-7,7), breaks = c(seq(-7, 7, by = 2), 0)) +   
#   labs(y = "Normalized Enrichment Score", x="", title = "CISvsAKI_Cortex_Disjoint") +
#   scale_fill_manual(values = c("negative" = "deepskyblue", "positive" = "orangered")) + 
#   geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.5) +
#   #geom_vline(xintercept = 0.5, color = "black", linetype = "solid", size = 1) +
#   theme_minimal() + 
#   theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add box around plot
#         axis.text.y = element_text(size = 14, color="black"), # Increase font size of pathway names
#         axis.ticks = element_line(color = "black"),  # Add tick marks
#         axis.ticks.length = unit(0.2, "cm"),  # Adjust tick mark length
#         axis.line = element_line(color = "black"))
# #dev.off()


# CIS INTERFACE vs AKI INTERFACE (KEGG AND HALLMARK) ~~~~~~~~~~~~~~~~~~~~~~~~~~~
comm_KEGG<-intersect(rownames(cis_interface_KEGG), rownames(aki_interface_KEGG) )
df0<-data.frame(cis_interface_KEGG[comm_KEGG,]$NES, aki_interface_KEGG[comm_KEGG,]$NES)
df0$des<-cis_interface_KEGG[comm_KEGG,]$Description
df0$Sdes<-gsub(" - Mus musculus.*", "", df0$des)
df0$Sdes<-paste0("KEGG_",df0$Sdes)
colnames(df0)<-c('CIS_interface_NES','AKI_interface_NES', 'CIS_AKI_Path', 'ShortPath')

comm_HALL<-intersect(rownames(cis_interface_HALL), rownames(aki_interface_HALL) )
df1<-data.frame(cis_interface_HALL[comm_HALL,]$NES, aki_interface_HALL[comm_HALL,]$NES)
df1$des<-cis_interface_HALL[comm_HALL,]$Description
df1$Sdes<-gsub("HALLMARK_", "HM_", df1$des)
colnames(df1)<-c('CIS_interface_NES','AKI_interface_NES', 'CIS_AKI_Path', 'ShortPath')

df<-rbind(df0,df1)
head(df)

##Covarying Pathways
df_covary<-df[which((df$CIS_interface_NES*df$AKI_interface_NES)>0),]
head(df_covary)

addWorksheet(wb1,"Interface_Covary")
writeData(wb1, "Interface_Covary", df_covary)

CIS_Int<-data.frame(Path=df_covary$ShortPath, CIS_Interface_NES=df_covary$CIS_interface_NES, Enrichment="positive")
#CIS_Med_plot<-rbind(head(CIS_Med), tail(CIS_Med))
CIS_Int_up<-CIS_Int[which(CIS_Int$CIS_Interface_NES>0),]
CIS_Int_down<-CIS_Int[which(CIS_Int$CIS_Interface_NES<0),]
CIS_Int_plot<-rbind(CIS_Int_up, CIS_Int_down)
CIS_Int_plot$Enrichment<-ifelse(CIS_Int_plot$CIS_Interface_NES<0, "negative", "positive")
CIS_Int_plot$CIS_Interface_NES<-ifelse(CIS_Int_plot$CIS_Interface_NES<0, -CIS_Int_plot$CIS_Interface_NES, CIS_Int_plot$CIS_Interface_NES)

AKI_Int<-data.frame(Path=df_covary$ShortPath, AKI_Interface_NES=df_covary$AKI_interface_NES, Enrichment="positive")
# AKI_Int_plot<-rbind(head(AKI_Int), tail(AKI_Int))
AKI_Int_up<-AKI_Int[which(AKI_Int$AKI_Interface_NES>0),]
AKI_Int_down<-AKI_Int[which(AKI_Int$AKI_Interface_NES<0),]
AKI_Int_plot<-rbind(AKI_Int_up, AKI_Int_down)
AKI_Int_plot$Enrichment<-ifelse(AKI_Int_plot$AKI_Interface_NES<0, "negative", "positive")
AKI_Int_plot$AKI_Interface_NES<-ifelse(AKI_Int_plot$AKI_Interface_NES>0, -AKI_Int_plot$AKI_Interface_NES, AKI_Int_plot$AKI_Interface_NES)

colnames(CIS_Int_plot)<-c('Pathway', 'NES', 'Enrichment')
colnames(AKI_Int_plot)<-c('Pathway', 'NES', 'Enrichment')
CIS_AKI_Interface_plot<-rbind(CIS_Int_plot,AKI_Int_plot)
CIS_AKI_Interface_plot<-CIS_AKI_Interface_plot[rev(order(CIS_AKI_Interface_plot$NES)),]

file<-paste0(dir,"/Codes2_Figures/Figure4/pdfs/CIS_vs_AKI_Covary_NES_Interface_Barplots.pdf")
pdf(file, width=7, height=2)
ggplot(CIS_AKI_Interface_plot, aes(x = rev(factor(Pathway, levels=unique(Pathway))), y = NES , fill = Enrichment)) +
  geom_bar(stat = "identity", position = "identity") + coord_flip() + 
  geom_bar(stat = "identity", position = "identity") + coord_flip() + 
  scale_x_discrete(labels = rev(unique(CIS_AKI_Interface_plot$Pathway))) +
  #geom_text(aes(label = Path, y=NES*1.8),  color = "black", size = 3, fontface="bold") +
  scale_y_continuous(limits=c(-7,7), breaks = c(seq(-7, 7, by = 2), 0)) +   
  labs(y = "Normalized Enrichment Score", x="", title = "CIS Outer Medulla vs AKI Outer Medulla") +
  scale_fill_manual(values = c("negative" = "deepskyblue", "positive" = "orangered")) + 
  geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.5) +
  #geom_vline(xintercept = 0.5, color = "black", linetype = "solid", size = 1) +
  theme_minimal() + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add box around plot
        #axis.text.x = element_text(size = 14, color="black"),
        axis.text.y = element_text(size = 14, color="black"), # Increase font size of pathway names
        axis.ticks = element_line(color = "black"),  # Add tick marks
        axis.ticks.length = unit(0.2, "cm"),  # Adjust tick mark length
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
dev.off()

##Disjoint Pathways
#Disjoint Pathways ************************************************************
df_disjoint<-df[which((df$CIS_interface_NES*df$AKI_interface_NES)<0),]
head(df_disjoint)

addWorksheet(wb1,"Interface_Disjoint")
writeData(wb1, "Interface_Disjoint", df_disjoint)

CIS_Int<-data.frame(Path=df_disjoint$ShortPath, CIS_Interface_NES=df_disjoint$CIS_interface_NES, Enrichment="positive")
#CIS_Int_plot<-rbind(head(CIS_Int), tail(CIS_Int))
CIS_Int_up<-CIS_Int[which(CIS_Int$CIS_Interface_NES>0),]
CIS_Int_down<-CIS_Int[which(CIS_Int$CIS_Interface_NES<0),]
CIS_Int_plot<-rbind(CIS_Int_up, CIS_Int_down)
CIS_Int_plot$Enrichment<-ifelse(CIS_Int_plot$CIS_Interface_NES<0, "negative", "positive")
CIS_Int_plot$CIS_Interface_NES<-ifelse(CIS_Int_plot$CIS_Interface_NES<0, -CIS_Int_plot$CIS_Interface_NES, CIS_Int_plot$CIS_Interface_NES)

AKI_Int<-data.frame(Path=df_disjoint$ShortPath, AKI_Interface_NES=df_disjoint$AKI_interface_NES, Enrichment="positive")
#AKI_Int_plot<-rbind(head(AKI_Int), tail(AKI_Int))
AKI_Int_up<-AKI_Int[which(AKI_Int$AKI_Interface_NES>0),]
AKI_Int_down<-AKI_Int[which(AKI_Int$AKI_Interface_NES<0),]
AKI_Int_plot<-rbind(AKI_Int_up, AKI_Int_down)
AKI_Int_plot$Enrichment<-ifelse(AKI_Int_plot$AKI_Interface_NES<0, "negative", "positive")
AKI_Int_plot$AKI_Interface_NES<-ifelse(AKI_Int_plot$AKI_Interface_NES>0, -AKI_Int_plot$AKI_Interface_NES, AKI_Int_plot$AKI_Interface_NES)

colnames(CIS_Int_plot)<-c('Pathway', 'NES', 'Enrichment')
colnames(AKI_Int_plot)<-c('Pathway', 'NES', 'Enrichment')
CIS_AKI_Interface_plot<-rbind(CIS_Int_plot,AKI_Int_plot)
#CIS_AKI_Interface_plot<-rbind(AKI_Int_plot,CIS_Int_plot)
CIS_AKI_Interface_plot<-CIS_AKI_Interface_plot[rev(order(CIS_AKI_Interface_plot$NES)),]

file<-paste0(dir,"/Codes2_Figures/Figure4/pdfs/CIS_vs_AKI_Disjoint_NES_Interface_Barplots.pdf")
pdf(file, width=10, height=10)
ggplot(CIS_AKI_Interface_plot, aes(x = factor(Pathway, levels=unique(Pathway)), y = NES , fill = Enrichment)) +
  #geom_bar(stat = "identity", position = "identity", width=0.7, alpha =0.5) + coord_flip() + 
  geom_bar(stat = "identity", position = "identity") + coord_flip() + 
  scale_x_discrete(labels = unique(CIS_AKI_Interface_plot$Pathway)) +
  #geom_text(aes(label = Path, y=NES*1.8),  color = "black", size = 3, fontface="bold") +
  scale_y_continuous(limits=c(-7,7), breaks = c(seq(-7, 7, by = 2), 0)) +   
  labs(y = "Normalized Enrichment Score", x="", title = "CISvsAKI_Interface_Disjoint") +
  scale_fill_manual(values = c("negative" = "deepskyblue", "positive" = "orangered")) + 
  geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.5) +
  #geom_vline(xintercept = 0.5, color = "black", linetype = "solid", size = 1) +
  theme_minimal() + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add box around plot
        axis.text.x = element_text(size = 16, color="black"),
        axis.text.y = element_text(size = 14, color="black"), # Increase font size of pathway names
        axis.ticks = element_line(color = "black"),  # Add tick marks
        axis.ticks.length = unit(0.2, "cm"),  # Adjust tick mark length
        axis.line = element_line(color = "black"))
dev.off()

# CIS MEDULLA vs AKI MEDULLA (KEGG AND HALLMARK) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Note: MEDULLA is same as Inner Medulla (as described in the main article section)
comm<-intersect(rownames(cis_medulla_KEGG), rownames(aki_medulla_KEGG) )
df0<-data.frame(cis_medulla_KEGG[comm,]$NES, aki_medulla_KEGG[comm,]$NES)
df0$des<-cis_medulla_KEGG[comm,]$Description
df0$Sdes<-gsub(" - Mus musculus.*", "", df0$des)
df0$Sdes<-paste0("KEGG_",df0$Sdes)
colnames(df0)<-c('CIS_medulla_NES','AKI_medulla_NES', 'CIS_AKI_Path', 'ShortPath')

comm_HALL<-intersect(rownames(cis_medulla_HALL), rownames(aki_medulla_HALL) )
df1<-data.frame(cis_medulla_HALL[comm_HALL,]$NES, aki_medulla_HALL[comm_HALL,]$NES)
df1$des<-cis_medulla_HALL[comm_HALL,]$Description
df1$Sdes<-gsub("HALLMARK_", "HM_", df1$des)
colnames(df1)<-c('CIS_medulla_NES','AKI_medulla_NES', 'CIS_AKI_Path', 'ShortPath')

df<-rbind(df0,df1)
head(df)

head(df)

#[NO COVARYING PATHWAY]
# df_covary<-df[which((df$CIS_medulla_NES*df$AKI_medulla_NES)>0),]
# head(df_covary)

addWorksheet(wb1,"Medulla_Covary")
writeData(wb1, "Medulla_Covary", df_covary)
# 
# CIS_Med<-data.frame(Path=df_covary$ShortPath, CIS_Medulla_NES=df_covary$CIS_medulla_NES, Enrichment="positive")
# #CIS_Med_plot<-rbind(head(CIS_Med), tail(CIS_Med))
# CIS_Med_up<-CIS_Med[which(CIS_Med$CIS_Medulla_NES>0),]
# CIS_Med_down<-CIS_Med[which(CIS_Med$CIS_Medulla_NES<0),]
# CIS_Med_plot<-rbind(CIS_Med_up, CIS_Med_down)
# CIS_Med_plot$Enrichment<-ifelse(CIS_Med_plot$CIS_Medulla_NES<0, "negative", "positive")
# CIS_Med_plot$CIS_Medulla_NES<-ifelse(CIS_Med_plot$CIS_Medulla_NES<0, -CIS_Med_plot$CIS_Medulla_NES, CIS_Med_plot$CIS_Medulla_NES)
# 
# AKI_Med<-data.frame(Path=df_covary$ShortPath, AKI_Medulla_NES=df_covary$AKI_medulla_NES, Enrichment="positive")
# # AKI_Med_plot<-rbind(head(AKI_Med), tail(AKI_Med))
# AKI_Med_up<-AKI_Med[which(AKI_Med$AKI_Medulla_NES>0),]
# AKI_Med_down<-AKI_Med[which(AKI_Med$AKI_Medulla_NES<0),]
# AKI_Med_plot<-rbind(AKI_Med_up, AKI_Med_down)
# AKI_Med_plot$Enrichment<-ifelse(AKI_Med_plot$AKI_Medulla_NES<0, "negative", "positive")
# AKI_Med_plot$AKI_Medulla_NES<-ifelse(AKI_Med_plot$AKI_Medulla_NES>0, -AKI_Med_plot$AKI_Medulla_NES, AKI_Med_plot$AKI_Medulla_NES)
# 
# colnames(CIS_Med_plot)<-c('Pathway', 'NES', 'Enrichment')
# colnames(AKI_Med_plot)<-c('Pathway', 'NES', 'Enrichment')
# CIS_AKI_Medulla_plot<-rbind(CIS_Med_plot,AKI_Med_plot)
# 
# #pdf("Figures/FIgure4/pdfs/CIS_vs_AKI_GSEA_NES_Medulla_Barplots.pdf", width=7, height=16)
# ggplot(CIS_AKI_Medulla_plot, aes(x = rev(factor(Pathway, levels=unique(Pathway))), y = NES , fill = Enrichment)) +
#   geom_bar(stat = "identity", position = "identity") + coord_flip() + 
#   scale_x_discrete(labels = rev(unique(CIS_AKI_Medulla_plot$Pathway))) +
#   #geom_text(aes(label = Path, y=NES*1.8),  color = "black", size = 3, fontface="bold") +
#   scale_y_continuous(limits=c(-7,7), breaks = c(seq(-7, 7, by = 2), 0)) +   
#   labs(y = "Normalized Enrichment Score", x="", title = "CIS Inner Medulla vs AKI Inner Medulla") +
#   scale_fill_manual(values = c("negative" = "deepskyblue", "positive" = "orangered")) + 
#   geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.5) +
#   #geom_vline(xintercept = 0.5, color = "black", linetype = "solid", size = 1) +
#   theme_minimal() + 
#   theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add box around plot
#         axis.text.y = element_text(size = 14, color="black"), # Increase font size of pathway names
#         axis.ticks = element_line(color = "black"),  # Add tick marks
#         axis.ticks.length = unit(0.2, "cm"),  # Adjust tick mark length
#         axis.line = element_line(color = "black"),
#         plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
#         legend.position="none")
# #dev.off()

##Disjoint Pathways
#Disjoint Pathways ************************************************************
df_disjoint<-df[which((df$CIS_medulla_NES*df$AKI_medulla_NES)<0),]
head(df_disjoint)

addWorksheet(wb1,"Medulla_Disjoint")
writeData(wb1, "Medulla_Disjoint", df_disjoint)

CIS_Med<-data.frame(Path=df_disjoint$ShortPath, CIS_Medulla_NES=df_disjoint$CIS_medulla_NES, Enrichment="positive")
#CIS_Med_plot<-rbind(head(CIS_Med), tail(CIS_Med))
CIS_Med_up<-CIS_Med[which(CIS_Med$CIS_Medulla_NES>0),]
CIS_Med_down<-CIS_Med[which(CIS_Med$CIS_Medulla_NES<0),]
CIS_Med_plot<-rbind(CIS_Med_up, CIS_Med_down)
CIS_Med_plot$Enrichment<-ifelse(CIS_Med_plot$CIS_Medulla_NES<0, "negative", "positive")
CIS_Med_plot$CIS_Medulla_NES<-ifelse(CIS_Med_plot$CIS_Medulla_NES<0, -CIS_Med_plot$CIS_Medulla_NES, CIS_Med_plot$CIS_Medulla_NES)

AKI_Med<-data.frame(Path=df_disjoint$ShortPath, AKI_Medulla_NES=df_disjoint$AKI_medulla_NES, Enrichment="positive")
#AKI_Med_plot<-rbind(head(AKI_Med), tail(AKI_Med))
AKI_Med_up<-AKI_Med[which(AKI_Med$AKI_Medulla_NES>0),]
AKI_Med_down<-AKI_Med[which(AKI_Med$AKI_Medulla_NES<0),]
AKI_Med_plot<-rbind(AKI_Med_up, AKI_Med_down)
AKI_Med_plot$Enrichment<-ifelse(AKI_Med_plot$AKI_Medulla_NES<0, "negative", "positive")
AKI_Med_plot$AKI_Medulla_NES<-ifelse(AKI_Med_plot$AKI_Medulla_NES>0, -AKI_Med_plot$AKI_Medulla_NES, AKI_Med_plot$AKI_Medulla_NES)

colnames(CIS_Med_plot)<-c('Pathway', 'NES', 'Enrichment')
colnames(AKI_Med_plot)<-c('Pathway', 'NES', 'Enrichment')
CIS_AKI_Medulla_plot<-rbind(CIS_Med_plot,AKI_Med_plot)
#CIS_AKI_Medulla_plot<-rbind(AKI_Med_Plot,CIS_Med_Plot)
CIS_AKI_Medulla_plot<-CIS_AKI_Medulla_plot[rev(order(CIS_AKI_Medulla_plot$NES)),]

file<-paste0(dir,"/Codes2_Figures/Figure4/pdfs/CIS_vs_AKI_Disjoint_NES_Medulla_Barplots.pdf")
pdf(file, width=10, height=10)
ggplot(CIS_AKI_Medulla_plot, aes(x = factor(Pathway, levels=unique(Pathway)), y = NES , fill = Enrichment)) +
  #geom_bar(stat = "identity", position = "identity", width=0.7, alpha =0.5) + coord_flip() + 
  geom_bar(stat = "identity", position = "identity") + coord_flip() + 
  scale_x_discrete(labels = unique(CIS_AKI_Medulla_plot$Pathway)) +
  #geom_text(aes(label = Path, y=NES*1.8),  color = "black", size = 3, fontface="bold") +
  scale_y_continuous(limits=c(-7,7), breaks = c(seq(-7, 7, by = 2), 0)) +   
  labs(y = "Normalized Enrichment Score", x="", title = "CISvsAKI_Medulla_Disjoint") +
  scale_fill_manual(values = c("negative" = "deepskyblue", "positive" = "orangered")) + 
  geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.5) +
  #geom_vline(xintercept = 0.5, color = "black", linetype = "solid", size = 1) +
  theme_minimal() + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add box around plot
        axis.text.x = element_text(size = 16, color="black"),
        axis.text.y = element_text(size = 14, color="black"), # Increase font size of pathway names
        axis.ticks = element_line(color = "black"),  # Add tick marks
        axis.ticks.length = unit(0.2, "cm"),  # Adjust tick mark length
        axis.line = element_line(color = "black"))
dev.off()

saveWorkbook(wb1,"Codes2_Tables/Manuscript/Fig4D-E-F_CISvsAKI_Covary_Disjoint_Pathways.xlsx", overwrite=TRUE)
