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

####### COMPARTMENTAL ENRICHED PATHWAYS BARPLOTS ###############################
#Note: Compartmental enriched pathways are ordered by their normalized enrichment
#scores (NES) obtained from gene set enrichemnt analysis.
#CIS Cortex
CIS_Cortex_KEGG<-read.xlsx("Codes2_Tables/CIS_KEGG_Pathays_Consistency_R2_Percentile95.xlsx", sheet="Cortex")
CIS_Cortex_KEGG$Description<-gsub(" - Mus.*","",CIS_Cortex_KEGG$Description)
CIS_Cortex_KEGG$Description<-paste0("KEGG_", CIS_Cortex_KEGG$Description)
df_Cortex_KEGG<-data.frame(Pathways=CIS_Cortex_KEGG$Description, NES=CIS_Cortex_KEGG$NES, Enrichment="Positive")
df_Cortex_KEGG[which(df_Cortex_KEGG$NES<0),]$Enrichment<-"Negative"
#rownames(df_Cortex)<-CIS_Cortex_KEGG$Description  
CIS_Cortex_HALL<-read.xlsx("Codes2_Tables/CIS_HALLMARK_Pathays_Consistency_R2_Percentile95.xlsx", sheet="Cortex")
CIS_Cortex_HALL$Description<-gsub("HALLMARK_","HM_",CIS_Cortex_HALL$Description)
df_Cortex_HALL<-data.frame(Pathways=CIS_Cortex_HALL$Description, NES=CIS_Cortex_HALL$NES, Enrichment="Positive")
df_Cortex_HALL[which(df_Cortex_HALL$NES<0),]$Enrichment<-"Negative"
#rownames(df_Cortex)<-CIS_Cortex_HALL$Description  

df_Cortex<-rbind(df_Cortex_KEGG, df_Cortex_HALL)
df_Cortex<-df_Cortex[rev(order(df_Cortex$NES)),]
head(df_Cortex)

df_Cortex$Pathways=factor(df_Cortex$Pathways, levels = rev(df_Cortex$Pathways))

pdf("Codes2_Figures/Figure3/pdfs/KEGG_HALL_NES_Cortex_Barplots.pdf", width=8, height=8)
ggplot(df_Cortex, aes(x = Pathways, y = NES, fill=Enrichment)) + geom_bar(stat="identity") +
  scale_fill_manual(values = c("Positive" = "orangered", "Negative" = "deepskyblue")) + 
  labs(x = "", y = "Normalized Enrichment Score", title = "Gene Set Enrichment Analysis (Cortex)") +
  theme_minimal() + scale_y_continuous(limits=c(-5,5), breaks = c(seq(-5, 5, by = 2), 0)) + 
  geom_hline(yintercept=0, linetype="solid", color="black", size=0.5) + coord_flip() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add box around plot
        axis.text.x = element_text(size = 14, color="black"), # Increase font size of pathway names
        axis.text.y = element_text(size = 14, color="black"),
        axis.ticks = element_line(color = "black"),  # Add tick marks
        axis.ticks.length = unit(0.2, "cm"),  # Adjust tick mark length
        axis.line = element_line(color = "black"),
        axis.title.x = element_text(size = 14))
dev.off()

#CIS Interface (Interface is same as Outer Medulla)
CIS_Interface_KEGG<-read.xlsx("Codes2_Tables/CIS_KEGG_Pathays_Consistency_R2_Percentile95.xlsx", sheet="Outer_Medulla")
CIS_Interface_KEGG$Description<-gsub(" - Mus.*","",CIS_Interface_KEGG$Description)
CIS_Interface_KEGG$Description<-paste0("KEGG_", CIS_Interface_KEGG$Description)
df_Interface_KEGG<-data.frame(Pathways=CIS_Interface_KEGG$Description, NES=CIS_Interface_KEGG$NES, Enrichment="Positive")
df_Interface_KEGG[which(df_Interface_KEGG$NES<0),]$Enrichment<-"Negative"
#rownames(df_Interface)<-CIS_Interface_KEGG$Description  

CIS_Interface_HALL<-read.xlsx("Codes2_Tables/CIS_HALLMARK_Pathays_Consistency_R2_Percentile95.xlsx", sheet="Outer_Medulla")
CIS_Interface_HALL$Description<-gsub("HALLMARK_*","HM_",CIS_Interface_HALL$Description)
df_Interface_HALL<-data.frame(Pathways=CIS_Interface_HALL$Description, NES=CIS_Interface_HALL$NES, Enrichment="Positive")
df_Interface_HALL[which(df_Interface_HALL$NES<0),]$Enrichment<-"Negative"

df_Interface<-rbind(df_Interface_KEGG, df_Interface_HALL)
df_Interface<-df_Interface[rev(order(df_Interface$NES)),]
head(df_Interface)

df_Interface$Pathways=factor(df_Interface$Pathways, levels = rev(df_Interface$Pathways))

pdf("Codes2_Figures/Figure3/pdfs/KEGG_HALL_NES_Interface_Barplots.pdf", width=9, height=8)
ggplot(df_Interface, aes(x = Pathways, y = NES, fill=Enrichment)) + geom_bar(stat="identity") +
  scale_fill_manual(values = c("Positive" = "orangered", "Negative" = "deepskyblue")) + 
  labs(x = "", y = "Normalized Enrichment Score", title = "Gene Set Enrichment Analysis (Interface)") +
  theme_minimal() + scale_y_continuous(limits=c(-5,5), breaks = c(seq(-5, 5, by = 2), 0)) + 
  geom_hline(yintercept=0, linetype="solid", color="black", size=0.5) + coord_flip() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add box around plot
        axis.text.x = element_text(size = 14, color="black"), # Increase font size of pathway names
        axis.text.y = element_text(size = 14, color="black"),
        axis.ticks = element_line(color = "black"),  # Add tick marks
        axis.ticks.length = unit(0.2, "cm"),  # Adjust tick mark length
        axis.line = element_line(color = "black"),
        axis.title.x = element_text(size = 14))
dev.off()

#CIS Medulla (Medulla is the same as Inner Medulla)
CIS_Medulla_KEGG<-read.xlsx("Codes2_Tables/CIS_KEGG_Pathays_Consistency_R2_Percentile95.xlsx", sheet="Inner_Medulla")
CIS_Medulla_KEGG$Description<-gsub(" - Mus.*","",CIS_Medulla_KEGG$Description)
CIS_Medulla_KEGG$Description<-paste0("KEGG_", CIS_Medulla_KEGG$Description)
df_Medulla_KEGG<-data.frame(Pathways=CIS_Medulla_KEGG$Description, NES=CIS_Medulla_KEGG$NES, Enrichment="Positive")
#df_Medulla_KEGG[which(df_Medulla_KEGG$NES<0),]$Enrichment<-"Negative"
#rownames(df_Medulla)<-CIS_Medulla_KEGG$Description  
CIS_Medulla_HALL<-read.xlsx("Codes2_Tables/CIS_HALLMARK_Pathays_Consistency_R2_Percentile95.xlsx", sheet="Inner_Medulla")
CIS_Medulla_HALL$Description<-gsub("HALLMARK_*","HM_",CIS_Medulla_HALL$Description)
df_Medulla_HALL<-data.frame(Pathways=CIS_Medulla_HALL$Description, NES=CIS_Medulla_HALL$NES, Enrichment="Positive")
df_Medulla_HALL[which(df_Medulla_HALL$NES<0),]$Enrichment<-"Negative"

df_Medulla<-rbind(df_Medulla_KEGG, df_Medulla_HALL)
df_Medulla<-df_Medulla[rev(order(df_Medulla$NES)),]
head(df_Medulla)

df_Medulla$Pathways=factor(df_Medulla$Pathways, levels = rev(df_Medulla$Pathways))

pdf("Codes2_Figures/Figure3/pdfs/KEGG_HALL_NES_Medulla_Barplots.pdf", width=9, height=8)
ggplot(df_Medulla, aes(x = Pathways, y = NES, fill=Enrichment)) + geom_bar(stat="identity") +
  scale_fill_manual(values = c("Positive" = "orangered", "Negative" = "deepskyblue")) + 
  labs(x = "", y = "Normalized Enrichment Score", title = "Gene Set Enrichment Analysis (Medulla)") +
  theme_minimal() + scale_y_continuous(limits=c(-5,5.1), breaks = c(seq(-5, 5, by = 2), 0)) + 
  geom_hline(yintercept=0, linetype="solid", color="black", size=0.5) + coord_flip() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add box around plot
        axis.text.x = element_text(size = 14, color="black"),# Increase font size of pathway names
        axis.text.y = element_text(size = 14, color="black"),
        axis.ticks = element_line(color = "black"),  # Add tick marks
        axis.ticks.length = unit(0.2, "cm"),  # Adjust tick mark length
        axis.line = element_line(color = "black"),
        axis.title.x = element_text(size = 14))
dev.off()