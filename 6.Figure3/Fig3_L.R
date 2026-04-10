#Clean environment
rm(list=ls(all.names=TRUE))
gc()

setwd("/Users/srujansingh/Library/CloudStorage/OneDrive-SharedLibraries-JohnsHopkins/Jean Fan - 2024_spatial_AKI/Revision2_GenomBiol")

#Load the three prime UTR lengths (Biomart)
load("threeprime_UTR.Rdata")

library(openxlsx)
library(ggplot2)

#Transcript Level 3-prime UTR values from Biomart
Transcript_three.UTR<-threeprime_UTR
Transcript_three.UTR<-Transcript_three.UTR[!is.na(Transcript_three.UTR$utr3_length),]
colnames(Transcript_three.UTR)<-c('Gene','GeneID', 'TranscriptID', 'Start', 'End', 'UTR_length')
head(Transcript_three.UTR)

#Gene Level UTR Length
Gene_three.UTR<-aggregate(UTR_length ~ Gene + GeneID,
                          data = Transcript_three.UTR,
                          FUN = mean)

Gene_three.UTR<-Gene_three.UTR[!duplicated(Gene_three.UTR$Gene),] #remove duplicates
rownames(Gene_three.UTR)<-Gene_three.UTR$Gene
head(Gene_three.UTR)

#Defining the theme:
mytheme<-theme(axis.text.x = element_text(color="black", size=14),  
               axis.text.y = element_text(color="black",size=14),
               axis.title.x=element_text(size=14),
               axis.title.y=element_text(size=14),
               panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
               legend.text = element_text(color="black", size=14),
               panel.grid.major = element_line(color = "grey75", linewidth = 0.5),
               plot.title=element_text(hjust=0.5),
               axis.ticks = element_line(color = "black"), axis.ticks.length = unit(0.1, "cm"))

Med.KEGG<-read.xlsx("Codes2_Tables/CIS_KEGG_Pathays_Consistency_R2_Percentile95.xlsx", sheet="Inner_Medulla")
Med.KEGG<-data.frame(Med.KEGG)

Med.OXPHOS<-Med.KEGG[Med.KEGG$ID=='mmu00190',]$core_enrichment #OXPHOS Pathway
Med.OXPHOS<-trimws(unlist(strsplit(Med.OXPHOS, ",")))

OXPHOS_genes<-intersect(Med.OXPHOS, Gene_three.UTR$Gene)

OXPHOS_Gene_UTR<-Gene_three.UTR[OXPHOS_genes,]
OXPHOS_Gene_UTR$Group<-'OXPHOS'
Rest_Gene_UTR<-Gene_three.UTR[!(rownames(Gene_three.UTR)%in%rownames(OXPHOS_Gene_UTR)), ]
Rest_Gene_UTR$Group<-'Not_OXPHOS'

#OXPHOS vs Up Genes (Inner Medulla)
Med.Up<-data.frame(read.xlsx("Codes2_Tables/Manuscript/ST2_CIS_Up_DEGs_Consistency_R2_Percentile95.xlsx", sheet="Inner_Medulla"))
rownames(Med.Up)<-Med.Up$Gene
Med.Up.genes<-Med.Up$Gene
Rest_Gene_UTR.Up<-intersect(Rest_Gene_UTR$Gene,Med.Up.genes)

#OXPHOS vs Down Genes (Inner Medulla)
Med.Down<-data.frame(read.xlsx("Codes2_Tables/Manuscript/ST1_CIS_Down_DEGs_Consistency_R2_Percentile95.xlsx", sheet="Inner_Medulla"))
Med.Down.genes<-Med.Down$Gene
Rest_Gene_UTR.Down<-intersect(Rest_Gene_UTR$Gene,Med.Down.genes)

##### COMBINE OXPHOS, UP and DONREGULATED GENES ################################
  OXPHOS_Gene_UTR1<-OXPHOS_Gene_UTR
  OXPHOS_Gene_UTR1$Trend<-'OXPHOS'
  
  Rest_Gene_UTR.Up1<-Rest_Gene_UTR[Rest_Gene_UTR.Up,]
  Rest_Gene_UTR.Up1$Trend<-'Up'
  
  Rest_Gene_UTR.Down1<-Rest_Gene_UTR[Rest_Gene_UTR.Down,]
  Rest_Gene_UTR.Down1$Trend<-'Down'
  
  df_Med.Combine<-rbind(OXPHOS_Gene_UTR1,Rest_Gene_UTR.Up1,Rest_Gene_UTR.Down1)
  df_Med.Combine$Group<-factor(df_Med.Combine$Group, levels=c('OXPHOS', 'Not_OXPHOS'))
  df_Med.Combine$Trend<-factor(df_Med.Combine$Trend, levels=c('OXPHOS', 'Up','Down'))
  
  test1<-wilcox.test(OXPHOS_Gene_UTR$UTR_length,Rest_Gene_UTR.Up1$UTR_length)
  test_p1<-paste0("p=", signif(test1$p.value,2))
  test_p1
  
  test2<-wilcox.test(OXPHOS_Gene_UTR$UTR_length,Rest_Gene_UTR.Down1$UTR_length)
  test_p2<-paste0("p=", signif(test2$p.value,2))
  test_p2
  
  test3<-wilcox.test(Rest_Gene_UTR.Up1$UTR_length,Rest_Gene_UTR.Down1$UTR_length)
  test_p3<-paste0("p=", signif(test3$p.value,2))
  test_p3
  
  OXPHOS_median<-paste0("M=",signif(median(OXPHOS_Gene_UTR$UTR_length),3))
  Up_median=paste0("M=", signif(median(Rest_Gene_UTR.Up1$UTR_length),3))
  Down_median=paste0("M=", signif(median(Rest_Gene_UTR.Down1$UTR_length),3))
  
  pdf("Codes2_Figures/Figure3/pdfs/UTR_OXPHOS_vs_Rest_All_DEGs_CIS_InnerMedulla.pdf", width=8, height=6)
  ggplot(df_Med.Combine, aes(x = Trend, y = UTR_length, fill = Trend)) +
    geom_violin(trim = TRUE, alpha = 0.5) +  # Violin
    geom_boxplot(width = 0.05, fill = "white", outlier.shape = NA) +
    ggtitle("Inner_Medulla (OXPHOS vs Rest_Genes)") +
    annotate("text", label=test_p1, x=1, y=6400, hjust=0, size=5, color="black") +
    annotate("text", label=test_p2, x=1.9, y=7600, hjust=0, size=5, color="black") +
    annotate("text", label=test_p3, x=2.3, y=6400, hjust=0, size=5, color="black") +
    annotate("text", label=OXPHOS_median, x=0.5, y=2000, hjust=0, size=5, color="black") +
    annotate("text", label=Up_median, x=1.5, y=2000, hjust=0, size=5, color="black") +
    annotate("text", label=Down_median, x=2.5, y=2000, hjust=0, size=5, color="black") +
    #geom_jitter(width = 0.2, size = 0.5, alpha = 0.7) +  # Points
    theme_classic() + mytheme +
    labs(x = "Gene Group", y = "3' UTR Length (bp)") +
    scale_fill_manual(values = c("OXPHOS" = "orange", "Up" = "red", "Down"='skyblue')) +
    scale_y_continuous(limits = c(0,8000))
  dev.off()
  