#Clean environment
rm(list=ls(all.names=TRUE))
gc()

dir<-"/Users/srujansingh/Library/CloudStorage/OneDrive-SharedLibraries-JohnsHopkins/Jean Fan - 2024_spatial_AKI/Revision2_GenomBiol"
setwd(dir)

########################## ADDING LIBRARY ######################################
library(openxlsx)
library(enrichplot) #For ClusterProfiler()

############ KEGG GSEA PATHWAYS PLOTTING #######################################
########### CIS CORTEX (OXPHOS, Thermogenesis, TCA Cycle) ######################
#Extracting gene set enrichment analysis results (refer to the code 
#Linear_Regression_CIS_AKI.R stored in the the folder 2.Linear_Regression_Modeling)
CIS_Cortex.KEGG<-readRDS("EnrichmentPlots/CIS_Cortex.KEGG_GSEA.rds")
CIS_Cortex.HALL<-readRDS("EnrichmentPlots/CIS_Cortex.HALL.GSEA.rds")

# pathway1<-"Oxidative phosphorylation"
# ID1=which(CIS_Cortex.KEGG$Description==paste0(pathway1, " - Mus musculus (house mouse)") )

pathway2<-"HALLMARK_OXIDATIVE_PHOSPHORYLATION"
ID2=which(CIS_Cortex.HALL$ID==pathway2)

# p1<-gseaplot2(CIS_Cortex.KEGG, geneSetID = c(ID1),  subplots=1:2, title="OXPHOS_Cor_KEGG")
#p2<-gseaplot2(CIS_Cortex.HALL, geneSetID = c(ID2),  subplots=1:2, title="OXPHOS_Cor_HALL")
#plot_list(p1, p2, ncol=1, tag_levels = 'A')

pdf("Codes2_Figures/Figure3/pdfs/OXPHOS_HALL_GSEAPlot_Cortex.pdf", width=10, height=4)
gseaplot2(CIS_Cortex.HALL, geneSetID = c(ID2),  subplots=1:2, title="OXPHOS_Cor_HALL")
dev.off()

# pdf("Figures/Figure3/pdfs/GSEA_Cortex_OXPHOS_TCA_Themogenesis.pdf", height=4, width=8)
# gseaplot2(CIS_Cortex.KEGG, geneSetID = c(ID2),  subplots=1:2, color= c("red" ), title="CIS Cortex") 
# dev.off()


########### CIS INTERFACE (OXPHOS, Thermogenesis, TCA Cycle)
#Extracting gene set enrichment analysis results (refer to the code 
#Linear_Regression_CIS_AKI.R stored in the the folder 2.Linear_Regression_Modeling)
#Note: INTERFACE is same as Outer Medulla (as described in the main article section)
CIS_Interface.KEGG<-readRDS("EnrichmentPlots/CIS_Interface.KEGG_GSEA.rds")
CIS_Interface.HALL<-readRDS("EnrichmentPlots/CIS_Interface.HALL_GSEA.rds")

pathway1<-"Oxidative phosphorylation"
ID1=which(CIS_Interface.KEGG$Description==paste0(pathway1, " - Mus musculus (house mouse)") )
dev.off()

pathway2<-"HALLMARK_OXIDATIVE_PHOSPHORYLATION"
ID2=which(CIS_Interface.HALL$ID==pathway2)

pdf("Codes2_Figures/Figure3/pdfs/OXPHOS_KEGG_GSEAPlot_Interface.pdf", width=10, height=4)
gseaplot2(CIS_Interface.KEGG, geneSetID = c(ID1),  subplots=1:2, title="OXPHOS_Int_KEGG")
dev.off()

pdf("Codes2_Figures/Figure3/pdfs/OXPHOS_HALL_GSEAPlot_Interface.pdf", width=10, height=4)
gseaplot2(CIS_Interface.HALL, geneSetID = c(ID2),  subplots=1:2, title="OXPHOS_Int_HALL")
dev.off()
# plot_list(p1, p2, ncol=1, tag_levels = 'A')

#pdf("Codes2_Figures/Figure3/pdfs/pdfs/OXPHOS_KEGG_HALL.pdf", height=4, width=8)
#gseaplot2(CIS_Interface.KEGG, geneSetID = c(ID1,ID3),  subplots=1:2, color= c("orange","red" ), title="CIS Interface")
#dev.off()

#gseaplot2(CIS_Interface.KEGG, geneSetID = c(ID1),  subplots=1:2, title="")

########### CIS MEDULLA (OXPHOS, Thermogenesis, TCA Cycle)
#Extracting gene set enrichment analysis results (refer to the code 
#Linear_Regression_CIS_AKI.R stored in the the folder 2.Linear_Regression_Modeling)
#Note: MEDULLA is same as Inner Medulla (as described in the main article section)
CIS_Medulla.KEGG<-readRDS("EnrichmentPlots/CIS_Medulla.KEGG_GSEA.rds")
CIS_Medulla.HALL<-readRDS("EnrichmentPlots/CIS_Medulla.HALL_GSEA.rds")

pathway1<-"Oxidative phosphorylation"
ID1=which(CIS_Medulla.KEGG$Description==paste0(pathway1, " - Mus musculus (house mouse)") )

pathway2<-"HALLMARK_OXIDATIVE_PHOSPHORYLATION"
ID2=which(CIS_Medulla.HALL$ID==pathway2)

pdf("Codes2_Figures/Figure3/pdfs/OXPHOS_KEGG_GSEAPlot_Medulla.pdf", width=10, height=4)
gseaplot2(CIS_Medulla.KEGG, geneSetID = c(ID1),  subplots=1:2, title="OXPHOS_Med_KEGG")
dev.off()

pdf("Codes2_Figures/Figure3/pdfs/OXPHOS_HALL_GSEAPlot_Medulla.pdf", width=10, height=4)
gseaplot2(CIS_Medulla.HALL, geneSetID = c(ID2),  subplots=1:2, title="OXPHOS_Med_HALL")
dev.off()

#pdf("Codes2_Figures/Figure3/pdfs/OXPHOS_CIS_Medulla_KEGG_HALL.pdf", height=6, width=8)
#plot_list(p1, p2, ncol=1, tag_levels = 'A')
#dev.off()

#gseaplot2(CIS_Medulla.KEGG, geneSetID = c(ID1),  subplots=1:2, title="")
