#Clean environment
rm(list=ls(all.names=TRUE))
gc()

dir<-"/Users/srujansingh/Library/CloudStorage/OneDrive-SharedLibraries-JohnsHopkins/Jean Fan - 2024_spatial_AKI/Revision2_GenomBiol"
setwd(dir)

#Load all libraries:
library(openxlsx)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(ggplot2)
library(msigdbr) #For HALLMARK Pathway
library(enrichplot) #For ClusterProfiler()
library(pathview)
library(biomaRt)
library(limma)

## load data
load("data/CIS_data.RData")
load("data/AKI_data.RData")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Read the genes: Consistency vs MeanR2 values of all genes.
file<-paste0(dir,"/Codes2_Tables/Consistency_vs_MeanR2_Compartmental.xlsx")

Const_R2_Cor<-read.xlsx(file, sheet="Cortex")
Const_R2_Int<-read.xlsx(file, sheet="Outer_Medulla")
Const_R2_Med<-read.xlsx(file, sheet="Inner_Medulla")

################################################################################
#Extracting the compartment specific spots for CIS and AKI:
#Note: 1.Interface is same as Outer Medulla. 
# 2. Medulla is same as Inner Medulla.

cortex <- readRDS('AKI-CIS-irl-ctrl_Cortex_spots.rds')
interface <- readRDS('AKI-CIS-irl-ctrl_Interface_spots.rds')
medulla <- readRDS('AKI-CIS-irl-ctrl_Medulla_spots.rds')
#all_spots<-readRDS('AKI-CIS-irl-ctrl_All_spots.rds')

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
CIS_0h$mat_notlog <- MERINGUE::normalizeCounts(CIS_0h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
CIS_12h$mat_notlog <- MERINGUE::normalizeCounts(CIS_12h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
CIS_24h$mat_notlog <- MERINGUE::normalizeCounts(CIS_24h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
CIS_48h$mat_notlog <- MERINGUE::normalizeCounts(CIS_48h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)

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
# 
# ################# normalize AKI DATASET (CPM) ##################################

AKI_sham$mat_notlog <- MERINGUE::normalizeCounts(AKI_sham$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
AKI_4h$mat_notlog <- MERINGUE::normalizeCounts(AKI_4h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
AKI_12h$mat_notlog <- MERINGUE::normalizeCounts(AKI_12h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
AKI_2d$mat_notlog <- MERINGUE::normalizeCounts(AKI_2d$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
AKI_6w$mat_notlog <- MERINGUE::normalizeCounts(AKI_6w$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)

######### LINEAR REGRESSION IN A COMPARTMENT-AGNOSTIC MANNER ###################
#### CIS CORTEX  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## subset to just cortex (CIS)

genes.shared <- Reduce(intersect, list(rownames(CIS_0h$mat_notlog), rownames(CIS_12h$mat_notlog),
                                       rownames(CIS_24h$mat_notlog), rownames(CIS_48h$mat_notlog),
                                       rownames(AKI_sham$mat_notlog), rownames(AKI_4h$mat_notlog),
                                       rownames(AKI_12h$mat_notlog), rownames(AKI_2d$mat_notlog)) )

cut_off<-0.026

# #Selecting n=190 random points from the whole CIS tissue cross section
sample_n<-190 #Setting the no. of random spots to be extracted
set.seed(200)
cis_rand_spots.cortex.0h<-sample(rownames(CIS_cortex.0h), sample_n, replace=FALSE)
cis_rand_spots.cortex.12h<-sample(rownames(CIS_cortex.12h), sample_n, replace=FALSE)
cis_rand_spots.cortex.24h<-sample(rownames(CIS_cortex.24h), sample_n, replace=FALSE)
cis_rand_spots.cortex.48h<-sample(rownames(CIS_cortex.48h), sample_n, replace=FALSE)

#Plot to check:
# plot(CIS_0h$pos[cis_rand_spots.cortex.0h,], pch = 21, bg = "black", col = "black", cex = 0.5, main="CIS_Cortex_0h")
# plot(CIS_12h$pos[cis_rand_spots.cortex.12h,], pch = 21, bg = "black", col = "black", cex = 0.5, main="CIS_Cortex_12h")
# plot(CIS_24h$pos[cis_rand_spots.cortex.24h,], pch = 21, bg = "black", col = "black", cex = 0.5, main="CIS_Cortex_24h")
# plot(CIS_48h$pos[cis_rand_spots.cortex.48h,], pch = 21, bg = "black", col = "black", cex = 0.5, main="CIS_Cortex_48h")

CIS_cortex_gexp <- data.frame(cbind(
  as.matrix(CIS_0h$mat_notlog[genes.shared,cis_rand_spots.cortex.0h]),
  as.matrix(CIS_12h$mat_notlog[genes.shared,cis_rand_spots.cortex.12h]),
  as.matrix(CIS_24h$mat_notlog[genes.shared,cis_rand_spots.cortex.24h]),
  as.matrix(CIS_48h$mat_notlog[genes.shared,cis_rand_spots.cortex.48h]) ))

meta <- data.frame(time = c(
  rep(log(0+1), length(cis_rand_spots.cortex.0h)),
  rep(log(12+1), length(cis_rand_spots.cortex.12h)), ## unit of days
  rep(log(24+1), length(cis_rand_spots.cortex.24h)),
  rep(log(48+1), length(cis_rand_spots.cortex.48h))
))  
rownames(meta) <- colnames(CIS_cortex_gexp)
table(meta)

#des = with(meta, model.matrix(~ time)) 
des = model.matrix(~ time, data = meta) #creating the design matrix
fit <- lmFit(CIS_cortex_gexp, design = des) ## fit limma model against time
head(fit$coefficients)
fit <- eBayes(fit)
top <- topTable(fit, number=nrow(CIS_cortex_gexp))
top <- top[top$adj.P.Val < 0.05,]
head(top)

cis_cortex_fc <- top$logFC
names(cis_cortex_fc) <- rownames(top)
cis_cortex_fc<-data.frame(names(cis_cortex_fc) ,cis_cortex_fc)
colnames(cis_cortex_fc)<-c('Gene','Slope')

head(cis_cortex_fc)

# #Calculate R2 value for each gene
fit_des<-fit$coefficients%*%t(des) #Recosntructing the fitted regression line.
Y<-CIS_cortex_gexp
SS_total<-rowSums((Y-rowMeans(Y))^2) # Y - mean(Y); Total Sum of squares per gene
SS_residual<-rowSums((Y-fit_des)^2) #Y - fitted Y, Residual Sum of squares per gene
R2<-1-(SS_residual/SS_total)
length(R2)

Des_R2<-cbind(rownames(Y),R2)
Des_R2<-Des_R2[rownames(top),]

table(rownames(Des_R2)%in%rownames(cis_cortex_fc)) #sanity check

cis_cortex_fc$R2<-as.numeric(Des_R2[,2])

cis_cortex_fc<-cis_cortex_fc[cis_cortex_fc$R2>=cut_off,]
head(cis_cortex_fc)
tail(cis_cortex_fc)

df<-data.frame(cis_cortex_fc)
df$Trend<-'Up'
df[df$Slope<0,]$Trend<-'Down'
df$Trend<-factor(df$Trend, levels=c('Up', 'Down'))
#ggplot(df, aes(x=Slope, y=R2, col=Trend)) + geom_point(size=0.2) + ggtitle("Cortex")

cis_cortex_up<-cis_cortex_fc[cis_cortex_fc$Slope>0,]
cis_cortex_up<-cis_cortex_up[rev(order(cis_cortex_up$Slope)),]
head(cis_cortex_up)
tail(cis_cortex_up)

cis_cortex_down<-cis_cortex_fc[cis_cortex_fc$Slope<0,]
cis_cortex_down<-cis_cortex_down[order(cis_cortex_down$Slope),]
head(cis_cortex_down)
tail(cis_cortex_down)

#Creating a ranked list of genes based on their slope:
cis_cortex_genes<-rbind(cis_cortex_up,cis_cortex_down)
cis_cortex_genes$rank<-c(rank(cis_cortex_up$Slope), -rank(abs(cis_cortex_down$Slope)))
cis_cortex_genes$logrank<-c(log10(rank(cis_cortex_up$Slope))+(1e-10), (log10(rank(abs(cis_cortex_down$Slope)))*-1)-(1e-10) )
colnames(cis_cortex_genes)<-c('Genes', 'Slope','R2', 'Rank', 'LogRank')
cis_cortex_genes<-cis_cortex_genes[rev(order(cis_cortex_genes$Rank)),]
print(cis_cortex_genes)

#Ranked list to be used for GSEA
CIS_Cortex<-cis_cortex_genes$LogRank
names(CIS_Cortex)<-cis_cortex_genes$Genes
head(CIS_Cortex)
tail(CIS_Cortex)

###### GSEA HALLMARK (CIS CORTEX) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Creating background gene list for GSEA 
gene_sets<-msigdbr(species = "mouse", category = "H")
background <- gene_sets[, c("gs_name", "gene_symbol")]

CIS_Cortex.HALL<-GSEA(geneList = CIS_Cortex, TERM2GENE = background, pvalueCutoff = 0.05, pAdjustMethod = "BH",)

CIS_Cortex.HALL.Plot<-CIS_Cortex.HALL
gseaplot(CIS_Cortex.HALL, geneSetID = 1)
saveRDS(CIS_Cortex.HALL.Plot, file="EnrichmentPlots/CIS_Cortex.HALL.GSEA.rds")

if(sign(max(CIS_Cortex.HALL$NES))==1) 
{ CIS_Cortex.HALL<-CIS_Cortex.HALL[rev(order(CIS_Cortex.HALL$NES)),]
} else { CIS_Cortex.HALL<-CIS_Cortex.HALL[order(CIS_Cortex.HALL$NES),] }

CIS_Cortex.HALL$core_enrichment<-gsub("/",",",CIS_Cortex.HALL$core_enrichment)
head(CIS_Cortex.HALL$Description)

###### GSEA KEGG PATHWAY (CIS CORTEX) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Convert to Entrez IDs
gene_ids <- bitr(names(CIS_Cortex), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db, drop=FALSE)
CIS_Cortex.KEGG<-CIS_Cortex
names(CIS_Cortex.KEGG)<-gene_ids$ENTREZID
head(CIS_Cortex.KEGG)

# Run KEGG GSEA
CIS_Cortex.KEGG <- gseKEGG(
  geneList = CIS_Cortex.KEGG,          # Ranked gene list with Entrez IDs
  organism = "mmu",                 # "hsa" for human; use "mmu" for mouse, etc.
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,              # P-value threshold
  verbose = FALSE,                   # Suppress verbose output
)

CIS_Cortex.KEGG.Plot<-CIS_Cortex.KEGG
gseaplot(CIS_Cortex.KEGG, geneSetID = 1)
saveRDS(CIS_Cortex.KEGG.Plot, file="EnrichmentPlots/CIS_Cortex.KEGG_GSEA.rds")

if(sign(max(CIS_Cortex.KEGG$NES))==1) 
{ CIS_Cortex.KEGG<-CIS_Cortex.KEGG[rev(order(CIS_Cortex.KEGG$NES)),]
} else { CIS_Cortex.KEGG<-CIS_Cortex.KEGG[order(CIS_Cortex.KEGG$NES),] }

head(CIS_Cortex.KEGG$Description)

#To convert the entrezID back to gene symbols
entrez_ids<-CIS_Cortex.KEGG$core_enrichment

split_ids <- strsplit(entrez_ids, "/") #To separate the entrezIDs 

gene_symbol<-lapply(split_ids, function(entrez_sublist) {
  bitr(entrez_sublist, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Mm.eg.db)
})

symbols_list <- lapply(gene_symbol, function(df) df$SYMBOL)
head(symbols_list)

xyz<-lapply(seq_along(CIS_Cortex.KEGG$core_enrichment), function(i) {
  CIS_Cortex.KEGG$core_enrichment[i]<-symbols_list[i]
  return(CIS_Cortex.KEGG$core_enrichment[i])
})

head(xyz)

CIS_Cortex.modKEGG<-CIS_Cortex.KEGG
CIS_Cortex.modKEGG$core_enrichment<-xyz

head(CIS_Cortex.modKEGG$core_enrichment)

### CIS INTERFACE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##Note: INTERFACE is same as Outer Medulla.
genes.shared <- Reduce(intersect, list(rownames(CIS_0h$mat_notlog), rownames(CIS_12h$mat_notlog),
                                       rownames(CIS_24h$mat_notlog), rownames(CIS_48h$mat_notlog),
                                       rownames(AKI_sham$mat_notlog), rownames(AKI_4h$mat_notlog),
                                       rownames(AKI_12h$mat_notlog), rownames(AKI_2d$mat_notlog)) )

set.seed(300)
cis_rand_spots.interface.0h<-sample(rownames(CIS_interface.0h), sample_n, replace=FALSE)
cis_rand_spots.interface.12h<-sample(rownames(CIS_interface.12h), sample_n, replace=FALSE)
cis_rand_spots.interface.24h<-sample(rownames(CIS_interface.24h), sample_n, replace=FALSE)
cis_rand_spots.interface.48h<-sample(rownames(CIS_interface.48h), sample_n, replace=FALSE)

CIS_interface_gexp <- data.frame(cbind(
  as.matrix(CIS_0h$mat_notlog[genes.shared,cis_rand_spots.interface.0h]),
  as.matrix(CIS_12h$mat_notlog[genes.shared,cis_rand_spots.interface.12h]),
  as.matrix(CIS_24h$mat_notlog[genes.shared,cis_rand_spots.interface.24h]),
  as.matrix(CIS_48h$mat_notlog[genes.shared,cis_rand_spots.interface.48h]) ))

meta <- data.frame(time = c(
  rep(log(0+1), length(cis_rand_spots.interface.0h)),
  rep(log(12+1), length(cis_rand_spots.interface.12h)), ## unit of days
  rep(log(24+1), length(cis_rand_spots.interface.24h)),
  rep(log(48+1), length(cis_rand_spots.interface.48h))
))  
rownames(meta) <- colnames(CIS_interface_gexp)
table(meta)

des = model.matrix(~ time, data = meta) #creating the design matrix
fit <- lmFit(CIS_interface_gexp, design = des) ## fit limma model against time
head(fit$coefficients)
fit <- eBayes(fit)
top <- topTable(fit, number=nrow(CIS_interface_gexp))
top <- top[top$adj.P.Val < 0.05,]
head(top)

cis_interface_fc <- top$logFC
names(cis_interface_fc) <- rownames(top)
cis_interface_fc<-data.frame(names(cis_interface_fc) ,cis_interface_fc)
colnames(cis_interface_fc)<-c('Gene','Slope')

head(cis_interface_fc)

# #Calculate R2 value for each gene
fit_des<-fit$coefficients%*%t(des) #Recosntructing the fitted regression line.
Y<-CIS_interface_gexp
SS_total<-rowSums((Y-rowMeans(Y))^2) # Y - mean(Y); Total Sum of squares per gene
SS_residual<-rowSums((Y-fit_des)^2) #Y - fitted Y, Residual Sum of squares per gene
R2<-1-(SS_residual/SS_total)
length(R2)

Des_R2<-cbind(rownames(Y),R2)
Des_R2<-Des_R2[rownames(top),]

table(rownames(Des_R2)%in%rownames(cis_interface_fc)) #sanity check

cis_interface_fc$R2<-as.numeric(Des_R2[,2])

cis_interface_fc<-cis_interface_fc[cis_interface_fc$R2>=cut_off,]
head(cis_interface_fc)
tail(cis_interface_fc)

df<-data.frame(cis_interface_fc)
df$Trend<-'Up'
df[df$Slope<0,]$Trend<-'Down'
df$Trend<-factor(df$Trend, levels=c('Up', 'Down'))

cis_interface_up<-cis_interface_fc[cis_interface_fc$Slope>=0,]
cis_interface_up<-cis_interface_up[rev(order(cis_interface_up$Slope)),]
head(cis_interface_up)
tail(cis_interface_up)

cis_interface_down<-cis_interface_fc[cis_interface_fc$Slope<0,]
cis_interface_down<-cis_interface_down[order(cis_interface_down$Slope),]
head(cis_interface_down)
tail(cis_interface_down)

#Creating a ranked list of genes based on their slope
cis_interface_genes<-rbind(cis_interface_up,cis_interface_down)
cis_interface_genes$rank<-c(rank(cis_interface_up$Slope), -rank(abs(cis_interface_down$Slope)))
cis_interface_genes$logrank<-c(log10(rank(cis_interface_up$Slope))+(1e-10), (log10(rank(abs(cis_interface_down$Slope)))*-1)-(1e-10) )
colnames(cis_interface_genes)<-c('Genes', 'Slope', 'R2', 'Rank', 'LogRank')
cis_interface_genes<-cis_interface_genes[rev(order(cis_interface_genes$Rank)),]
print(cis_interface_genes)

#Ranked list to be used for GSEA
CIS_Interface<-cis_interface_genes$LogRank
names(CIS_Interface)<-cis_interface_genes$Genes
head(CIS_Interface)
tail(CIS_Interface)

###### GSEA HALLMARK (CIS INTERFACE) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Creating background gene list for GSEA 
gene_sets<-msigdbr(species = "mouse", category = "H")
background <- gene_sets[, c("gs_name", "gene_symbol")]

CIS_Interface.HALL<-GSEA(geneList = CIS_Interface, TERM2GENE = background, pvalueCutoff = 0.05, pAdjustMethod = "BH",)

CIS_Interface.HALL.Plot<-CIS_Interface.HALL
gseaplot(CIS_Interface.HALL.Plot, geneSetID = 1)
saveRDS(CIS_Interface.HALL.Plot, file="EnrichmentPlots/CIS_Interface.HALL_GSEA.rds")

if(sign(max(CIS_Interface.HALL$NES))==1) 
{ CIS_Interface.HALL<-CIS_Interface.HALL[rev(order(CIS_Interface.HALL$NES)),]
} else { CIS_Interface.HALL<-CIS_Interface.HALL[order(CIS_Interface.HALL$NES),] }

CIS_Interface.HALL$core_enrichment<-gsub("/",",",CIS_Interface.HALL$core_enrichment)
head(CIS_Interface.HALL$Description)

###### GSEA KEGG PATHWAY (CIS INTERFACE) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Convert to Entrez IDs
gene_ids <- bitr(names(CIS_Interface), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db, drop=FALSE)
CIS_Interface.KEGG<-CIS_Interface
names(CIS_Interface.KEGG)<-gene_ids$ENTREZID
head(CIS_Interface.KEGG)

# Run KEGG GSEA
CIS_Interface.KEGG <- gseKEGG(
  geneList = CIS_Interface.KEGG,          # Ranked gene list with Entrez IDs
  organism = "mmu",                 # "hsa" for human; use "mmu" for mouse, etc.
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,              # P-value threshold
  verbose = FALSE                   # Suppress verbose output
)

CIS_Interface.KEGG.Plot<-CIS_Interface.KEGG
gseaplot(CIS_Interface.KEGG.Plot, geneSetID = 1)
saveRDS(CIS_Interface.KEGG.Plot, file="EnrichmentPlots/CIS_Interface.KEGG_GSEA.rds")

if(sign(max(CIS_Interface.KEGG$NES))==1) 
{ CIS_Interface.KEGG<-CIS_Interface.KEGG[rev(order(CIS_Interface.KEGG$NES)),]
} else { CIS_Interface.KEGG<-CIS_Interface.KEGG[order(CIS_Interface.KEGG$NES),] }

head(CIS_Interface.KEGG$Description)

#To convert the entrezID back to gene symbols
entrez_ids<-CIS_Interface.KEGG$core_enrichment

split_ids <- strsplit(entrez_ids, "/") #To separate the entrezIDs 

gene_symbol<-lapply(split_ids, function(entrez_sublist) {
  bitr(entrez_sublist, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Mm.eg.db)
})

symbols_list <- lapply(gene_symbol, function(df) df$SYMBOL)
head(symbols_list)

xyz<-lapply(seq_along(CIS_Interface.KEGG$core_enrichment), function(i) {
  CIS_Interface.KEGG$core_enrichment[i]<-symbols_list[i]
  return(CIS_Interface.KEGG$core_enrichment[i])
})

head(xyz)

CIS_Interface.modKEGG<-CIS_Interface.KEGG
CIS_Interface.modKEGG$core_enrichment<-xyz

#### CIS MEDULLA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Note: MEDULLA is same as Inner Medulla.

genes.shared <- Reduce(intersect, list(rownames(CIS_0h$mat_notlog), rownames(CIS_12h$mat_notlog),
                                       rownames(CIS_24h$mat_notlog), rownames(CIS_48h$mat_notlog),
                                       rownames(AKI_sham$mat_notlog), rownames(AKI_4h$mat_notlog),
                                       rownames(AKI_12h$mat_notlog), rownames(AKI_2d$mat_notlog)) )

sample_n<-190
#Selecting 190 random points from the whole CIS medulla cross section
set.seed(400)
cis_rand_spots.medulla.0h<-sample(rownames(CIS_medulla.0h), sample_n, replace=FALSE)
cis_rand_spots.medulla.12h<-sample(rownames(CIS_medulla.12h), sample_n, replace=FALSE)
cis_rand_spots.medulla.24h<-sample(rownames(CIS_medulla.24h), sample_n, replace=FALSE)
cis_rand_spots.medulla.48h<-sample(rownames(CIS_medulla.48h), sample_n, replace=FALSE)

CIS_medulla_gexp <- data.frame(cbind(
  as.matrix(CIS_0h$mat_notlog[genes.shared,cis_rand_spots.medulla.0h]),
  as.matrix(CIS_12h$mat_notlog[genes.shared,cis_rand_spots.medulla.12h]),
  as.matrix(CIS_24h$mat_notlog[genes.shared,cis_rand_spots.medulla.24h]),
  as.matrix(CIS_48h$mat_notlog[genes.shared,cis_rand_spots.medulla.48h]) ))

meta <- data.frame(time = c(
  rep(log(0+1), length(cis_rand_spots.medulla.0h)),
  rep(log(12+1), length(cis_rand_spots.medulla.12h)), ## unit of days
  rep(log(24+1), length(cis_rand_spots.medulla.24h)),
  rep(log(48+1), length(cis_rand_spots.medulla.48h))
))  
rownames(meta) <- colnames(CIS_medulla_gexp)
table(meta)

des = model.matrix(~ time, data = meta) #creating the design matrix
fit <- lmFit(CIS_medulla_gexp, design = des) ## fit limma model against time
head(fit$coefficients)
fit <- eBayes(fit)
top <- topTable(fit, number=nrow(CIS_medulla_gexp))
top <- top[top$adj.P.Val < 0.05,]
head(top)

cis_medulla_fc <- top$logFC
names(cis_medulla_fc) <- rownames(top)
cis_medulla_fc<-data.frame(names(cis_medulla_fc) ,cis_medulla_fc)
colnames(cis_medulla_fc)<-c('Gene','Slope')

head(cis_medulla_fc)

# #Calculate R2 value for each gene
fit_des<-fit$coefficients%*%t(des) #Recosntructing the fitted regression line.
Y<-CIS_medulla_gexp
SS_total<-rowSums((Y-rowMeans(Y))^2) # Y - mean(Y); Total Sum of squares per gene
SS_residual<-rowSums((Y-fit_des)^2) #Y - fitted Y, Residual Sum of squares per gene
R2<-1-(SS_residual/SS_total)
length(R2)

Des_R2<-cbind(rownames(Y),R2)
Des_R2<-Des_R2[rownames(top),]

table(rownames(Des_R2)%in%rownames(cis_medulla_fc)) #sanity check

cis_medulla_fc$R2<-as.numeric(Des_R2[,2])

cis_medulla_fc<-cis_medulla_fc[cis_medulla_fc$R2>=cut_off,]
head(cis_medulla_fc)
tail(cis_medulla_fc)

df<-data.frame(cis_medulla_fc)
df$Trend<-'Up'
df[df$Slope<0,]$Trend<-'Down'
df$Trend<-factor(df$Trend, levels=c('Up', 'Down'))

cis_medulla_up<-cis_medulla_fc[cis_medulla_fc$Slope>=0,]
cis_medulla_up<-cis_medulla_up[rev(order(cis_medulla_up$Slope)),]
head(cis_medulla_up)
tail(cis_medulla_up)

cis_medulla_down<-cis_medulla_fc[cis_medulla_fc$Slope<0,]
cis_medulla_down<-cis_medulla_down[order(cis_medulla_down$Slope),]
head(cis_medulla_down)
tail(cis_medulla_down)

#Creating a ranked list of genes based on their slope
cis_medulla_genes<-rbind(cis_medulla_up,cis_medulla_down)
cis_medulla_genes$rank<-c(rank(cis_medulla_up$Slope), -rank(abs(cis_medulla_down$Slope)))
cis_medulla_genes$logrank<-c(log10(rank(cis_medulla_up$Slope))+(1e-10), (log10(rank(abs(cis_medulla_down$Slope)))*-1)-(1e-10) )
colnames(cis_medulla_genes)<-c('Genes', 'Slope', 'R2', 'Rank', 'LogRank')
cis_medulla_genes<-cis_medulla_genes[rev(order(cis_medulla_genes$Rank)),]
print(cis_medulla_genes)

#Ranked list to be used for GSEA
CIS_Medulla<-cis_medulla_genes$LogRank
names(CIS_Medulla)<-cis_medulla_genes$Genes
head(CIS_Medulla)
tail(CIS_Medulla)

###### GSEA HALLMARK (CIS MEDULLA) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Creating background gene list for GSEA 
gene_sets<-msigdbr(species = "mouse", category = "H")
background <- gene_sets[, c("gs_name", "gene_symbol")]

CIS_Medulla.HALL<-GSEA(geneList = CIS_Medulla, TERM2GENE = background, pvalueCutoff = 0.05, pAdjustMethod = "BH",)

CIS_Medulla.HALL.Plot<-CIS_Medulla.HALL
gseaplot(CIS_Medulla.HALL.Plot, geneSetID = 1)
saveRDS(CIS_Medulla.HALL.Plot, file="EnrichmentPlots/CIS_Medulla.HALL_GSEA.rds")

if(sign(max(CIS_Medulla.HALL$NES))==1) 
{ CIS_Medulla.HALL<-CIS_Medulla.HALL[rev(order(CIS_Medulla.HALL$NES)),]
} else { CIS_Medulla.HALL<-CIS_Medulla.HALL[order(CIS_Medulla.HALL$NES),] }

CIS_Medulla.HALL$core_enrichment<-gsub("/",",",CIS_Medulla.HALL$core_enrichment)
head(CIS_Medulla.HALL$Description)

# Run KEGG GSEA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Convert to Entrez IDs
gene_ids <- bitr(names(CIS_Medulla), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db, drop=FALSE)
CIS_Medulla.KEGG<-CIS_Medulla
names(CIS_Medulla.KEGG)<-gene_ids$ENTREZID
head(CIS_Medulla.KEGG)

CIS_Medulla.KEGG <- gseKEGG(
  geneList = CIS_Medulla.KEGG,          # Ranked gene list with Entrez IDs
  organism = "mmu",                 # "hsa" for human; use "mmu" for mouse, etc.
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,              # P-value threshold
  verbose = FALSE                   # Suppress verbose output
)


CIS_Medulla.KEGG.Plot<-CIS_Medulla.KEGG
gseaplot(CIS_Medulla.KEGG.Plot, geneSetID = 1)
saveRDS(CIS_Medulla.KEGG.Plot, file="EnrichmentPlots/CIS_Medulla.KEGG_GSEA.rds")

if(sign(max(CIS_Medulla.KEGG$NES))==1) 
{ CIS_Medulla.KEGG<-CIS_Medulla.KEGG[rev(order(CIS_Medulla.KEGG$NES)),]
} else { CIS_Medulla.KEGG<-CIS_Medulla.KEGG[order(CIS_Medulla.KEGG$NES),] }

head(CIS_Medulla.KEGG$Description)

#To convert the entrezID back to gene symbols
entrez_ids<-CIS_Medulla.KEGG$core_enrichment

split_ids <- strsplit(entrez_ids, "/") #To separate the entrezIDs 

gene_symbol<-lapply(split_ids, function(entrez_sublist) {
  bitr(entrez_sublist, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Mm.eg.db)
})

symbols_list <- lapply(gene_symbol, function(df) df$SYMBOL)
head(symbols_list)

xyz<-lapply(seq_along(CIS_Medulla.KEGG$core_enrichment), function(i) {
  CIS_Medulla.KEGG$core_enrichment[i]<-symbols_list[i]
  return(CIS_Medulla.KEGG$core_enrichment[i])
})

head(xyz)

CIS_Medulla.modKEGG<-CIS_Medulla.KEGG
CIS_Medulla.modKEGG$core_enrichment<-xyz

head(CIS_Medulla.modKEGG$core_enrichment)

################################################################################
############################ AKI SPOTS #########################################
#1.All_Cortex_Spots
AKI_cortex.sham<-cortex[grep("^AKI_sham", names(cortex))]
AKI_cortex.4h<-cortex[grep("^AKI_4h", names(cortex))]
AKI_cortex.12h<-cortex[grep("^AKI_12h", names(cortex))]
AKI_cortex.2d<-cortex[grep("^AKI_2d", names(cortex))]

AKI_cortex.sham<-cbind(AKI_sham$pos[names(AKI_cortex.sham),], AKI_cortex.sham)
AKI_cortex.4h<-cbind(AKI_4h$pos[names(AKI_cortex.4h),], AKI_cortex.4h)
AKI_cortex.12h<-cbind(AKI_12h$pos[names(AKI_cortex.12h),], AKI_cortex.12h)
AKI_cortex.2d<-cbind(AKI_2d$pos[names(AKI_cortex.2d),], AKI_cortex.2d)

colnames(AKI_cortex.sham)<-c('X','Y', 'Cluster')
colnames(AKI_cortex.4h)<-c('X','Y', 'Cluster')
colnames(AKI_cortex.12h)<-c('X','Y', 'Cluster')
colnames(AKI_cortex.2d)<-c('X','Y', 'Cluster')

#To plot and check/verify the spots are right
ggplot(data.frame(AKI_cortex.sham), aes(x=X, y=Y, col=Cluster)) + geom_point(size=0.5) + ggtitle("AKI_Cortex_sham")
ggplot(data.frame(AKI_cortex.4h), aes(x=X, y=Y, col=Cluster)) + geom_point(size=0.5) + ggtitle("AKI_Cortex_4h")
ggplot(data.frame(AKI_cortex.12h), aes(x=X, y=Y, col=Cluster)) + geom_point(size=0.5) + ggtitle("AKI_Cortex_12h")
ggplot(data.frame(AKI_cortex.2d), aes(x=X, y=Y, col=Cluster)) + geom_point(size=0.5) + ggtitle("AKI_Cortex_2d")

#2.All Interface spots
AKI_interface.sham<-interface[grep("^AKI_sham", names(interface))]
AKI_interface.4h<-interface[grep("^AKI_4h", names(interface))]
AKI_interface.12h<-interface[grep("^AKI_12h", names(interface))]
AKI_interface.2d<-interface[grep("^AKI_2d", names(interface))]

AKI_interface.sham<-cbind(AKI_sham$pos[names(AKI_interface.sham),], AKI_interface.sham)
AKI_interface.4h<-cbind(AKI_4h$pos[names(AKI_interface.4h),], AKI_interface.4h)
AKI_interface.12h<-cbind(AKI_12h$pos[names(AKI_interface.12h),], AKI_interface.12h)
AKI_interface.2d<-cbind(AKI_2d$pos[names(AKI_interface.2d),], AKI_interface.2d)

colnames(AKI_interface.sham)<-c('X','Y', 'Cluster')
colnames(AKI_interface.4h)<-c('X','Y', 'Cluster')
colnames(AKI_interface.12h)<-c('X','Y', 'Cluster')
colnames(AKI_interface.2d)<-c('X','Y', 'Cluster')

#To plot and check/verify the spots are right
ggplot(data.frame(AKI_interface.sham), aes(x=X, y=Y, col=Cluster)) + geom_point(size=0.5) + ggtitle("AKI_Interface_sham")
ggplot(data.frame(AKI_interface.4h), aes(x=X, y=Y, col=Cluster)) + geom_point(size=0.5) + ggtitle("AKI_Interface_4h")
ggplot(data.frame(AKI_interface.12h), aes(x=X, y=Y, col=Cluster)) + geom_point(size=0.5) + ggtitle("AKI_Interface_12h")
ggplot(data.frame(AKI_interface.2d), aes(x=X, y=Y, col=Cluster)) + geom_point(size=0.5) + ggtitle("AKI_Interface_2d")


#3.All medullary spots
AKI_medulla.sham<-medulla[grep("^AKI_sham", names(medulla))]
AKI_medulla.4h<-medulla[grep("^AKI_4h", names(medulla))]
AKI_medulla.12h<-medulla[grep("^AKI_12h", names(medulla))]
AKI_medulla.2d<-medulla[grep("^AKI_2d", names(medulla))]

AKI_medulla.sham<-cbind(AKI_sham$pos[names(AKI_medulla.sham),], AKI_medulla.sham)
AKI_medulla.4h<-cbind(AKI_4h$pos[names(AKI_medulla.4h),], AKI_medulla.4h)
AKI_medulla.12h<-cbind(AKI_12h$pos[names(AKI_medulla.12h),], AKI_medulla.12h)
AKI_medulla.2d<-cbind(AKI_2d$pos[names(AKI_medulla.2d),], AKI_medulla.2d)

colnames(AKI_medulla.sham)<-c('X','Y', 'Cluster')
colnames(AKI_medulla.4h)<-c('X','Y', 'Cluster')
colnames(AKI_medulla.12h)<-c('X','Y', 'Cluster')
colnames(AKI_medulla.2d)<-c('X','Y', 'Cluster')

# #To plot and check/verify the spots are right
# ggplot(data.frame(AKI_medulla.sham), aes(x=X, y=Y, col=Cluster)) + geom_point(size=0.5) + ggtitle("AKI_Medulla_sham")
# ggplot(data.frame(AKI_medulla.4h), aes(x=X, y=Y, col=Cluster)) + geom_point(size=0.5) + ggtitle("AKI_Medulla_4h")
# ggplot(data.frame(AKI_medulla.12h), aes(x=X, y=Y, col=Cluster)) + geom_point(size=0.5) + ggtitle("AKI_Medulla_12h")
# ggplot(data.frame(AKI_medulla.2d), aes(x=X, y=Y, col=Cluster)) + geom_point(size=0.5) + ggtitle("AKI_Medulla_2d")

################################################################################
#### LINEAR REGRESSION ANALYSIS FOR WARM ISCHMEMIA-REPERFUSION (AKI) DATASET####

#Read the genes: Consistency vs MeanR2 values of all genes.
file<-paste0(dir,"/Codes2_Tables/AKI_Consistency_vs_MeanR2_Compartmental.xlsx")

AKI_Const_R2_Cor<-read.xlsx(file, sheet="Cortex")
AKI_Const_R2_Int<-read.xlsx(file, sheet="Outer_Medulla")
AKI_Const_R2_Med<-read.xlsx(file, sheet="Inner_Medulla")

#### AKI CORTEX ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## subset to just cortex (AKI)
genes.shared <- Reduce(intersect, list(rownames(CIS_0h$mat_notlog), rownames(CIS_12h$mat_notlog),
                                       rownames(CIS_24h$mat_notlog), rownames(CIS_48h$mat_notlog),
                                       rownames(AKI_sham$mat_notlog), rownames(AKI_4h$mat_notlog),
                                       rownames(AKI_12h$mat_notlog), rownames(AKI_2d$mat_notlog)) )

#cutoff<-quantile(AKI_Const_R2_Cor$mean_R2, 0.95, na.rm=TRUE)
cut_off<-0.093

#Selecting 190 random points from the whole AKI cortex cross section
set.seed(600)
aki_rand_spots.cortex.sham<-sample(rownames(AKI_cortex.sham), sample_n, replace=FALSE)
aki_rand_spots.cortex.4h<-sample(rownames(AKI_cortex.4h), sample_n, replace=FALSE)
aki_rand_spots.cortex.12h<-sample(rownames(AKI_cortex.12h), sample_n, replace=FALSE)
aki_rand_spots.cortex.2d<-sample(rownames(AKI_cortex.2d), sample_n, replace=FALSE)

#Plot to check:
plot(AKI_sham$pos[aki_rand_spots.cortex.sham,], pch = 21, bg = "black", col = "black", cex = 0.5, main="AKI_Cortex_sham")
plot(AKI_4h$pos[aki_rand_spots.cortex.4h,], pch = 21, bg = "black", col = "black", cex = 0.5, main="AKI_Cortex_4h")
plot(AKI_12h$pos[aki_rand_spots.cortex.12h,], pch = 21, bg = "black", col = "black", cex = 0.5, main="AKI_Cortex_12h")
plot(AKI_2d$pos[aki_rand_spots.cortex.2d,], pch = 21, bg = "black", col = "black", cex = 0.5, main="AKI_Cortex_2d")

AKI_cortex_gexp <- data.frame(cbind(
  as.matrix(AKI_sham$mat_notlog[genes.shared,aki_rand_spots.cortex.sham]),
  as.matrix(AKI_4h$mat_notlog[genes.shared,aki_rand_spots.cortex.4h]),
  as.matrix(AKI_12h$mat_notlog[genes.shared,aki_rand_spots.cortex.12h]),
  as.matrix(AKI_2d$mat_notlog[genes.shared,aki_rand_spots.cortex.2d]) ))

meta <- data.frame(time = c(
  rep(log(0+1), length(aki_rand_spots.cortex.sham)),
  rep(log(4+1), length(aki_rand_spots.cortex.4h)), ## unit of days
  rep(log(12+1), length(aki_rand_spots.cortex.12h)),
  rep(log(48+1), length(aki_rand_spots.cortex.2d))
))  
rownames(meta) <- colnames(AKI_cortex_gexp)
table(meta)

des = model.matrix(~ time, data = meta) #creating the design matrix
fit <- lmFit(AKI_cortex_gexp, design = des) ## fit limma model against time

head(fit$coefficients)
fit <- eBayes(fit)
top <- topTable(fit, number=nrow(AKI_cortex_gexp))
top <- top[top$adj.P.Val < 0.05,]
head(top)

aki_cortex_fc <- top$logFC
names(aki_cortex_fc) <- rownames(top)
aki_cortex_fc<-data.frame(names(aki_cortex_fc) ,aki_cortex_fc)
colnames(aki_cortex_fc)<-c('Gene','Slope')

head(aki_cortex_fc)

# #Calculate R2 value for each gene
fit_des<-fit$coefficients%*%t(des) #Recosntructing the fitted regression line.

Y<-AKI_cortex_gexp
SS_total<-rowSums((Y-rowMeans(Y))^2) # Y - mean(Y); Total Sum of squares per gene
SS_residual<-rowSums((Y-fit_des)^2) #Y - fitted Y, Residual Sum of squares per gene
R2<-1-(SS_residual/SS_total)
length(R2)

Des_R2<-cbind(rownames(Y),R2)
Des_R2<-Des_R2[rownames(top),]

table(rownames(Des_R2)%in%rownames(aki_cortex_fc)) #sanity check

aki_cortex_fc$R2<-as.numeric(Des_R2[,2])

aki_cortex_fc<-aki_cortex_fc[aki_cortex_fc$R2>=cut_off,]
head(aki_cortex_fc)
tail(aki_cortex_fc)

df<-data.frame(aki_cortex_fc)
df$Trend<-'Up'
df[df$Slope<0,]$Trend<-'Down'
df$Trend<-factor(df$Trend, levels=c('Up', 'Down'))

aki_cortex_up<-aki_cortex_fc[aki_cortex_fc$Slope>=0,]
aki_cortex_up<-aki_cortex_up[rev(order(aki_cortex_up$Slope)),]
head(aki_cortex_up)
tail(aki_cortex_up)

aki_cortex_down<-aki_cortex_fc[aki_cortex_fc$Slope<0,]
aki_cortex_down<-aki_cortex_down[order(aki_cortex_down$Slope),]
head(aki_cortex_down)
tail(aki_cortex_down)

#Creating a ranked list of genes based on their slope
aki_cortex_genes<-rbind(aki_cortex_up,aki_cortex_down)
print(aki_cortex_genes)
aki_cortex_genes$rank<-c(rank(aki_cortex_up$Slope), -rank(abs(aki_cortex_down$Slope)))
aki_cortex_genes$logrank<-c(log10(rank(aki_cortex_up$Slope))+(1e-10), (log10(rank(abs(aki_cortex_down$Slope)))*-1)-(1e-10) )
colnames(aki_cortex_genes)<-c('Genes', 'Slope', 'R2', 'Rank', 'LogRank')
aki_cortex_genes<-aki_cortex_genes[rev(order(aki_cortex_genes$Rank)),]
print(aki_cortex_genes)

#Ranked list to be used for GSEA
AKI_Cortex<-aki_cortex_genes$LogRank
names(AKI_Cortex)<-aki_cortex_genes$Genes
head(AKI_Cortex)
tail(AKI_Cortex)

###### GSEA HALLMARK (AKI CORTEX) **********************************************
#Creating background gene list for GSEA 
gene_sets<-msigdbr(species = "mouse", category = "H")
background <- gene_sets[, c("gs_name", "gene_symbol")]

AKI_Cortex.HALL<-GSEA(geneList = AKI_Cortex, TERM2GENE = background, pvalueCutoff = 0.05, pAdjustMethod = "BH",)

AKI_Cortex.HALL.Plot<-AKI_Cortex.HALL
gseaplot(AKI_Cortex.HALL, geneSetID = 1)
saveRDS(AKI_Cortex.HALL.Plot, file="EnrichmentPlots/AKI_Cortex.HALL_GSEA.rds")

if(sign(max(AKI_Cortex.HALL$NES))==1) 
{ AKI_Cortex.HALL<-AKI_Cortex.HALL[rev(order(AKI_Cortex.HALL$NES)),]
} else { AKI_Cortex.HALL<-AKI_Cortex.HALL[order(AKI_Cortex.HALL$NES),] }

AKI_Cortex.HALL$core_enrichment<-gsub("/",",",AKI_Cortex.HALL$core_enrichment)
head(AKI_Cortex.HALL$Description)

###### GSEA KEGG PATHWAY (AKI CORTEX) ******************************************
# Convert to Entrez IDs
gene_ids <- bitr(names(AKI_Cortex), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db, drop=FALSE)
AKI_Cortex.KEGG<-AKI_Cortex
names(AKI_Cortex.KEGG)<-gene_ids$ENTREZID
head(AKI_Cortex.KEGG)

# Run KEGG GSEA
AKI_Cortex.KEGG <- gseKEGG(
  geneList = AKI_Cortex.KEGG,          # Ranked gene list with Entrez IDs
  organism = "mmu",                 # "hsa" for human; use "mmu" for mouse, etc.
  pvalueCutoff = 0.05,              # P-value threshold
  verbose = FALSE                   # Suppress verbose output
)

AKI_Cortex.KEGG.Plot<-AKI_Cortex.KEGG
gseaplot(AKI_Cortex.KEGG, geneSetID = 1)
saveRDS(AKI_Cortex.KEGG.Plot, file="EnrichmentPlots/AKI_Cortex.KEGG_GSEA.rds")

if(sign(max(AKI_Cortex.KEGG$NES))==1) 
{ AKI_Cortex.KEGG<-AKI_Cortex.KEGG[rev(order(AKI_Cortex.KEGG$NES)),]
} else { AKI_Cortex.KEGG<-AKI_Cortex.KEGG[order(AKI_Cortex.KEGG$NES),] }

head(AKI_Cortex.KEGG$Description)

#To convert the entrezID back to gene symbols
entrez_ids<-AKI_Cortex.KEGG$core_enrichment

split_ids <- strsplit(entrez_ids, "/") #To separate the entrezIDs 

gene_symbol<-lapply(split_ids, function(entrez_sublist) {
  bitr(entrez_sublist, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Mm.eg.db)
})

symbols_list <- lapply(gene_symbol, function(df) df$SYMBOL)
head(symbols_list)

xyz<-lapply(seq_along(AKI_Cortex.KEGG$core_enrichment), function(i) {
  #print(AKI_Cortex.KEGG$core_enrichment[i])
  #print(symbols_list[[i]])
  AKI_Cortex.KEGG$core_enrichment[i]<-symbols_list[i]
  #print(AKI_Cortex.KEGG$core_enrichment[i])
  return(AKI_Cortex.KEGG$core_enrichment[i])
})

head(xyz)

AKI_Cortex.modKEGG<-AKI_Cortex.KEGG
AKI_Cortex.modKEGG$core_enrichment<-xyz

#head(CIS_Global.KEGG$core_enrichment)
head(AKI_Cortex.modKEGG$core_enrichment)

#### AKI INTERFACE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## subset to just interface (AKI)
#Note: Interface is same as Outer medulla.

#Selecting 190 random points from the whole AKI interface cross section
set.seed(700)
aki_rand_spots.interface.sham<-sample(rownames(AKI_interface.sham), sample_n, replace=FALSE)
aki_rand_spots.interface.4h<-sample(rownames(AKI_interface.4h), sample_n, replace=FALSE)
aki_rand_spots.interface.12h<-sample(rownames(AKI_interface.12h), sample_n, replace=FALSE)
aki_rand_spots.interface.2d<-sample(rownames(AKI_interface.2d), sample_n, replace=FALSE)

#Plot to check:
plot(AKI_sham$pos[aki_rand_spots.interface.sham,], pch = 21, bg = "black", col = "black", cex = 0.5, main="AKI_Interface_sham")
plot(AKI_4h$pos[aki_rand_spots.interface.4h,], pch = 21, bg = "black", col = "black", cex = 0.5, main="AKI_Interface_4h")
plot(AKI_12h$pos[aki_rand_spots.interface.12h,], pch = 21, bg = "black", col = "black", cex = 0.5, main="AKI_Interface_12h")
plot(AKI_2d$pos[aki_rand_spots.interface.2d,], pch = 21, bg = "black", col = "black", cex = 0.5, main="AKI_Interface_2d")

AKI_interface_gexp <- data.frame(cbind(
  as.matrix(AKI_sham$mat_notlog[genes.shared,aki_rand_spots.interface.sham]),
  as.matrix(AKI_4h$mat_notlog[genes.shared,aki_rand_spots.interface.4h]),
  as.matrix(AKI_12h$mat_notlog[genes.shared,aki_rand_spots.interface.12h]),
  as.matrix(AKI_2d$mat_notlog[genes.shared,aki_rand_spots.interface.2d]) ))

meta <- data.frame(time = c(
  rep(log(0+1), length(aki_rand_spots.interface.sham)),
  rep(log(4+1), length(aki_rand_spots.interface.4h)), ## unit of days
  rep(log(12+1), length(aki_rand_spots.interface.12h)),
  rep(log(48+1), length(aki_rand_spots.interface.2d))
))  
rownames(meta) <- colnames(AKI_interface_gexp)
table(meta)

des = model.matrix(~ time, data = meta) #creating the design matrix
fit <- lmFit(AKI_interface_gexp, design = des) ## fit limma model against time

head(fit$coefficients)
fit <- eBayes(fit)
top <- topTable(fit, number=nrow(AKI_interface_gexp))
top <- top[top$adj.P.Val < 0.05,]
head(top)

aki_interface_fc <- top$logFC
names(aki_interface_fc) <- rownames(top)
aki_interface_fc<-data.frame(names(aki_interface_fc) ,aki_interface_fc)
colnames(aki_interface_fc)<-c('Gene','Slope')

head(aki_interface_fc)

# #Calculate R2 value for each gene
fit_des<-fit$coefficients%*%t(des) #Recosntructing the fitted regression line.

Y<-AKI_interface_gexp
SS_total<-rowSums((Y-rowMeans(Y))^2) # Y - mean(Y); Total Sum of squares per gene
SS_residual<-rowSums((Y-fit_des)^2) #Y - fitted Y, Residual Sum of squares per gene
R2<-1-(SS_residual/SS_total)
length(R2)

Des_R2<-cbind(rownames(Y),R2)
Des_R2<-Des_R2[rownames(top),]

table(rownames(Des_R2)%in%rownames(aki_interface_fc)) #sanity check

aki_interface_fc$R2<-as.numeric(Des_R2[,2])

#cutoff<-quantile(AKI_Const_R2_Int$mean_R2, 0.95, na.rm=TRUE)

aki_interface_fc<-aki_interface_fc[aki_interface_fc$R2>=cut_off,]
head(aki_interface_fc)
tail(aki_interface_fc)

df<-data.frame(aki_interface_fc)
df$Trend<-'Up'
df[df$Slope<0,]$Trend<-'Down'
df$Trend<-factor(df$Trend, levels=c('Up', 'Down'))

aki_interface_up<-aki_interface_fc[aki_interface_fc$Slope>=0,]
aki_interface_up<-aki_interface_up[rev(order(aki_interface_up$Slope)),]
head(aki_interface_up)
tail(aki_interface_up)

aki_interface_down<-aki_interface_fc[aki_interface_fc$Slope<0,]
aki_interface_down<-aki_interface_down[order(aki_interface_down$Slope),]
head(aki_interface_down)
tail(aki_interface_down)

#Creating a ranked list of genes based on their slope
aki_interface_genes<-rbind(aki_interface_up,aki_interface_down)
print(aki_interface_genes)
aki_interface_genes$rank<-c(rank(aki_interface_up$Slope), -rank(abs(aki_interface_down$Slope)))
aki_interface_genes$logrank<-c(log10(rank(aki_interface_up$Slope))+(1e-10), (log10(rank(abs(aki_interface_down$Slope)))*-1)-(1e-10) )
colnames(aki_interface_genes)<-c('Genes', 'Slope', 'R2', 'Rank', 'LogRank')
aki_interface_genes<-aki_interface_genes[rev(order(aki_interface_genes$Rank)),]
print(aki_interface_genes)


#Ranked list to be used for GSEA
AKI_Interface<-aki_interface_genes$LogRank
names(AKI_Interface)<-aki_interface_genes$Genes
head(AKI_Interface)
tail(AKI_Interface)

###### GSEA HALLMARK (AKI INTERFACE) **********************************************
#Creating background gene list for GSEA 
gene_sets<-msigdbr(species = "mouse", category = "H")
background <- gene_sets[, c("gs_name", "gene_symbol")]

AKI_Interface.HALL<-GSEA(geneList = AKI_Interface, TERM2GENE = background, pvalueCutoff = 0.05, pAdjustMethod = "BH",)

AKI_Interface.HALL.Plot<-AKI_Interface.HALL
gseaplot(AKI_Interface.HALL, geneSetID = 1)
saveRDS(AKI_Interface.HALL.Plot, file="EnrichmentPlots/AKI_Interface.HALL_GSEA.rds")

if(sign(max(AKI_Interface.HALL$NES))==1) 
{ AKI_Interface.HALL<-AKI_Interface.HALL[rev(order(AKI_Interface.HALL$NES)),]
} else { AKI_Interface.HALL<-AKI_Interface.HALL[order(AKI_Interface.HALL$NES),] }

AKI_Interface.HALL$core_enrichment<-gsub("/",",",AKI_Interface.HALL$core_enrichment)
head(AKI_Interface.HALL$Description)

###### GSEA KEGG PATHWAY (AKI INTERFACE) ******************************************
# Convert to Entrez IDs
gene_ids <- bitr(names(AKI_Interface), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db, drop=FALSE)
AKI_Interface.KEGG<-AKI_Interface
names(AKI_Interface.KEGG)<-gene_ids$ENTREZID
head(AKI_Interface.KEGG)

# Run KEGG GSEA
AKI_Interface.KEGG <- gseKEGG(
  geneList = AKI_Interface.KEGG,          # Ranked gene list with Entrez IDs
  organism = "mmu",                 # "hsa" for human; use "mmu" for mouse, etc.
  pvalueCutoff = 0.05,              # P-value threshold
  verbose = FALSE                   # Suppress verbose output
)

AKI_Interface.KEGG.Plot<-AKI_Interface.KEGG
gseaplot(AKI_Interface.KEGG, geneSetID = 1)
saveRDS(AKI_Interface.KEGG.Plot, file="EnrichmentPlots/AKI_Interface.KEGG_GSEA.rds")

if(sign(max(AKI_Interface.KEGG$NES))==1) 
{ AKI_Interface.KEGG<-AKI_Interface.KEGG[rev(order(AKI_Interface.KEGG$NES)),]
} else { AKI_Interface.KEGG<-AKI_Interface.KEGG[order(AKI_Interface.KEGG$NES),] }

head(AKI_Interface.KEGG$Description)

#To convert the entrezID back to gene symbols
entrez_ids<-AKI_Interface.KEGG$core_enrichment

split_ids <- strsplit(entrez_ids, "/") #To separate the entrezIDs 

gene_symbol<-lapply(split_ids, function(entrez_sublist) {
  bitr(entrez_sublist, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Mm.eg.db)
})

symbols_list <- lapply(gene_symbol, function(df) df$SYMBOL)
head(symbols_list)

xyz<-lapply(seq_along(AKI_Interface.KEGG$core_enrichment), function(i) {
  #print(AKI_Interface.KEGG$core_enrichment[i])
  #print(symbols_list[[i]])
  AKI_Interface.KEGG$core_enrichment[i]<-symbols_list[i]
  #print(AKI_Interface.KEGG$core_enrichment[i])
  return(AKI_Interface.KEGG$core_enrichment[i])
})

head(xyz)

AKI_Interface.modKEGG<-AKI_Interface.KEGG
AKI_Interface.modKEGG$core_enrichment<-xyz

#head(CIS_Global.KEGG$core_enrichment)
head(AKI_Interface.modKEGG$core_enrichment)

#### AKI MEDULLA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## subset to just medulla (AKI)
#Note: medulla is same as Inner Medulla

#Selecting 190 random points from the whole AKI medulla cross section
set.seed(800)
aki_rand_spots.medulla.sham<-sample(rownames(AKI_medulla.sham), sample_n, replace=FALSE)
aki_rand_spots.medulla.4h<-sample(rownames(AKI_medulla.4h), sample_n, replace=FALSE)
aki_rand_spots.medulla.12h<-sample(rownames(AKI_medulla.12h), sample_n, replace=FALSE)
aki_rand_spots.medulla.2d<-sample(rownames(AKI_medulla.2d), sample_n, replace=FALSE)

#Plot to check:
plot(AKI_sham$pos[aki_rand_spots.medulla.sham,], pch = 21, bg = "black", col = "black", cex = 0.5, main="AKI_Medulla_sham")
plot(AKI_4h$pos[aki_rand_spots.medulla.4h,], pch = 21, bg = "black", col = "black", cex = 0.5, main="AKI_Medulla_4h")
plot(AKI_12h$pos[aki_rand_spots.medulla.12h,], pch = 21, bg = "black", col = "black", cex = 0.5, main="AKI_Medulla_12h")
plot(AKI_2d$pos[aki_rand_spots.medulla.2d,], pch = 21, bg = "black", col = "black", cex = 0.5, main="AKI_Medulla_2d")

AKI_medulla_gexp <- data.frame(cbind(
  as.matrix(AKI_sham$mat_notlog[genes.shared,aki_rand_spots.medulla.sham]),
  as.matrix(AKI_4h$mat_notlog[genes.shared,aki_rand_spots.medulla.4h]),
  as.matrix(AKI_12h$mat_notlog[genes.shared,aki_rand_spots.medulla.12h]),
  as.matrix(AKI_2d$mat_notlog[genes.shared,aki_rand_spots.medulla.2d]) ))

meta <- data.frame(time = c(
  rep(log(0+1), length(aki_rand_spots.medulla.sham)),
  rep(log(4+1), length(aki_rand_spots.medulla.4h)), ## unit of days
  rep(log(12+1), length(aki_rand_spots.medulla.12h)),
  rep(log(48+1), length(aki_rand_spots.medulla.2d))
))  
rownames(meta) <- colnames(AKI_medulla_gexp)
table(meta)

des = model.matrix(~ time, data = meta) #creating the design matrix
fit <- lmFit(AKI_medulla_gexp, design = des) ## fit limma model against time

head(fit$coefficients)
fit <- eBayes(fit)
top <- topTable(fit, number=nrow(AKI_medulla_gexp))
top <- top[top$adj.P.Val < 0.05,]
head(top)

aki_medulla_fc <- top$logFC
names(aki_medulla_fc) <- rownames(top)
aki_medulla_fc<-data.frame(names(aki_medulla_fc) ,aki_medulla_fc)
colnames(aki_medulla_fc)<-c('Gene','Slope')

head(aki_medulla_fc)

# #Calculate R2 value for each gene
fit_des<-fit$coefficients%*%t(des) #Recosntructing the fitted regression line.

Y<-AKI_medulla_gexp
SS_total<-rowSums((Y-rowMeans(Y))^2) # Y - mean(Y); Total Sum of squares per gene
SS_residual<-rowSums((Y-fit_des)^2) #Y - fitted Y, Residual Sum of squares per gene
R2<-1-(SS_residual/SS_total)
length(R2)

Des_R2<-cbind(rownames(Y),R2)
Des_R2<-Des_R2[rownames(top),]

table(rownames(Des_R2)%in%rownames(aki_medulla_fc)) #sanity check

aki_medulla_fc$R2<-as.numeric(Des_R2[,2])

#cutoff<-quantile(AKI_Const_R2_Med$mean_R2, 0.95, na.rm=TRUE)

aki_medulla_fc<-aki_medulla_fc[aki_medulla_fc$R2>=cut_off,]
head(aki_medulla_fc)
tail(aki_medulla_fc)

df<-data.frame(aki_medulla_fc)
df$Trend<-'Up'
df[df$Slope<0,]$Trend<-'Down'
df$Trend<-factor(df$Trend, levels=c('Up', 'Down'))

aki_medulla_up<-aki_medulla_fc[aki_medulla_fc$Slope>=0,]
aki_medulla_up<-aki_medulla_up[rev(order(aki_medulla_up$Slope)),]
head(aki_medulla_up)
tail(aki_medulla_up)

aki_medulla_down<-aki_medulla_fc[aki_medulla_fc$Slope<0,]
aki_medulla_down<-aki_medulla_down[order(aki_medulla_down$Slope),]
head(aki_medulla_down)
tail(aki_medulla_down)

#Creating a ranked list of genes based on their slope
aki_medulla_genes<-rbind(aki_medulla_up,aki_medulla_down)
print(aki_medulla_genes)
aki_medulla_genes$rank<-c(rank(aki_medulla_up$Slope), -rank(abs(aki_medulla_down$Slope)))
aki_medulla_genes$logrank<-c(log10(rank(aki_medulla_up$Slope))+(1e-10), (log10(rank(abs(aki_medulla_down$Slope)))*-1)-(1e-10) )
colnames(aki_medulla_genes)<-c('Genes', 'Slope', 'R2', 'Rank', 'LogRank')
aki_medulla_genes<-aki_medulla_genes[rev(order(aki_medulla_genes$Rank)),]
print(aki_medulla_genes)

#Ranked list to be used for GSEA
AKI_Medulla<-aki_medulla_genes$LogRank
names(AKI_Medulla)<-aki_medulla_genes$Genes
head(AKI_Medulla)
tail(AKI_Medulla)

###### GSEA HALLMARK (AKI MEDULLA) **********************************************
#Creating background gene list for GSEA 
gene_sets<-msigdbr(species = "mouse", category = "H")
background <- gene_sets[, c("gs_name", "gene_symbol")]

AKI_Medulla.HALL<-GSEA(geneList = AKI_Medulla, TERM2GENE = background, pvalueCutoff = 0.05, pAdjustMethod = "BH",)

AKI_Medulla.HALL.Plot<-AKI_Medulla.HALL
gseaplot(AKI_Medulla.HALL, geneSetID = 1)
saveRDS(AKI_Medulla.HALL.Plot, file="EnrichmentPlots/AKI_Medulla.HALL_GSEA.rds")

if(sign(max(AKI_Medulla.HALL$NES))==1) 
{ AKI_Medulla.HALL<-AKI_Medulla.HALL[rev(order(AKI_Medulla.HALL$NES)),]
} else { AKI_Medulla.HALL<-AKI_Medulla.HALL[order(AKI_Medulla.HALL$NES),] }

AKI_Medulla.HALL$core_enrichment<-gsub("/",",",AKI_Medulla.HALL$core_enrichment)
head(AKI_Medulla.HALL$Description)

###### GSEA KEGG PATHWAY (AKI MEDULLA) ******************************************
# Convert to Entrez IDs
gene_ids <- bitr(names(AKI_Medulla), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db, drop=FALSE)
AKI_Medulla.KEGG<-AKI_Medulla
names(AKI_Medulla.KEGG)<-gene_ids$ENTREZID
head(AKI_Medulla.KEGG)

# Run KEGG GSEA
AKI_Medulla.KEGG <- gseKEGG(
  geneList = AKI_Medulla.KEGG,          # Ranked gene list with Entrez IDs
  organism = "mmu",                 # "hsa" for human; use "mmu" for mouse, etc.
  pvalueCutoff = 0.05,              # P-value threshold
  verbose = FALSE                   # Suppress verbose output
)

AKI_Medulla.KEGG.Plot<-AKI_Medulla.KEGG
gseaplot(AKI_Medulla.KEGG, geneSetID = 1)
saveRDS(AKI_Medulla.KEGG.Plot, file="EnrichmentPlots/AKI_Medulla.KEGG_GSEA.rds")

if(sign(max(AKI_Medulla.KEGG$NES))==1) 
{ AKI_Medulla.KEGG<-AKI_Medulla.KEGG[rev(order(AKI_Medulla.KEGG$NES)),]
} else { AKI_Medulla.KEGG<-AKI_Medulla.KEGG[order(AKI_Medulla.KEGG$NES),] }

head(AKI_Medulla.KEGG$Description)

#To convert the entrezID back to gene symbols
entrez_ids<-AKI_Medulla.KEGG$core_enrichment

split_ids <- strsplit(entrez_ids, "/") #To separate the entrezIDs 

gene_symbol<-lapply(split_ids, function(entrez_sublist) {
  bitr(entrez_sublist, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Mm.eg.db)
})

symbols_list <- lapply(gene_symbol, function(df) df$SYMBOL)
head(symbols_list)

xyz<-lapply(seq_along(AKI_Medulla.KEGG$core_enrichment), function(i) {

  AKI_Medulla.KEGG$core_enrichment[i]<-symbols_list[i]
  return(AKI_Medulla.KEGG$core_enrichment[i])
})

head(xyz)

AKI_Medulla.modKEGG<-AKI_Medulla.KEGG
AKI_Medulla.modKEGG$core_enrichment<-xyz

#head(CIS_Global.KEGG$core_enrichment)
head(AKI_Medulla.modKEGG$core_enrichment)


################################################################################
#Save Workbooks (CIS Dataset):
wb0<-createWorkbook()

addWorksheet(wb0,"Inner_Medulla")
writeData(wb0, "Inner_Medulla", CIS_Medulla.HALL)

addWorksheet(wb0,"Outer_Medulla")
writeData(wb0, "Outer_Medulla", CIS_Interface.HALL)

addWorksheet(wb0,"Cortex")
writeData(wb0, "Cortex", CIS_Cortex.HALL)

#dir<-"/Users/srujansingh/Library/CloudStorage/OneDrive-SharedLibraries-JohnsHopkins/Jean Fan - 2024_spatial_AKI/Revision2_GenomBiol"
file<-paste0(dir,"/Codes2_Tables/CIS_HALLMARK_Pathays_Consistency_R2_Percentile95.xlsx")
saveWorkbook(wb0, file, overwrite = TRUE)


wb<-createWorkbook()

addWorksheet(wb,"Inner_Medulla")
writeData(wb, "Inner_Medulla", CIS_Medulla.modKEGG)

addWorksheet(wb,"Outer_Medulla")
writeData(wb, "Outer_Medulla", CIS_Interface.modKEGG)

addWorksheet(wb,"Cortex")
writeData(wb, "Cortex", CIS_Cortex.modKEGG)

file<-paste0(dir,"/Codes2_Tables/CIS_KEGG_Pathays_Consistency_R2_Percentile95.xlsx")
saveWorkbook(wb, file, overwrite = TRUE)


wb1<-createWorkbook()

addWorksheet(wb1,"Inner_Medulla")
writeData(wb1, "Inner_Medulla", cis_medulla_up)

addWorksheet(wb1,"Outer_Medulla")
writeData(wb1, "Outer_Medulla", cis_interface_up)

addWorksheet(wb1,"Cortex")
writeData(wb1, "Cortex", cis_cortex_up)

file<-paste0(dir,"/Codes2_Tables/CIS_Up_DEGs_Consistency_R2_Percentile95.xlsx")
saveWorkbook(wb1, file, overwrite = TRUE)


wb2<-createWorkbook()

addWorksheet(wb2,"Inner_Medulla")
writeData(wb2, "Inner_Medulla", cis_medulla_down)

addWorksheet(wb2,"Outer_Medulla")
writeData(wb2, "Outer_Medulla", cis_interface_down)

addWorksheet(wb2,"Cortex")
writeData(wb2, "Cortex", cis_cortex_down)

file<-paste0(dir,"/Codes2_Tables/CIS_Down_DEGs_Consistency_R2_Percentile95.xlsx")
saveWorkbook(wb2, file, overwrite = TRUE)

wb3<-createWorkbook()

addWorksheet(wb3,"Inner_Medulla")
writeData(wb3, "Inner_Medulla", cis_medulla_genes)

addWorksheet(wb3,"Outer_Medulla")
writeData(wb3, "Outer_Medulla", cis_interface_genes)

addWorksheet(wb3,"Cortex")
writeData(wb3, "Cortex", cis_cortex_genes)

file<-paste0(dir,"/Codes2_Tables/CIS_All_DEGs_Rank_Consistency_R2_Percentile95.xlsx")
saveWorkbook(wb3, file, overwrite = TRUE)


################################################################################
#Save Workbooks (AKI Dataset):
wb0<-createWorkbook()

addWorksheet(wb0,"Inner_Medulla")
writeData(wb0, "Inner_Medulla", AKI_Medulla.HALL)

addWorksheet(wb0,"Outer_Medulla")
writeData(wb0, "Outer_Medulla", AKI_Interface.HALL)

addWorksheet(wb0,"Cortex")
writeData(wb0, "Cortex", AKI_Cortex.HALL)

#dir<-"/Users/srujansingh/Library/CloudStorage/OneDrive-SharedLibraries-JohnsHopkins/Jean Fan - 2024_spatial_AKI/Revision2_GenomBiol"
file<-paste0(dir,"/Codes2_Tables/AKI_HALLMARK_Pathways_Consistency_R2_Percentile95.xlsx")
saveWorkbook(wb0, file, overwrite = TRUE)


wb<-createWorkbook()

addWorksheet(wb,"Inner_Medulla")
writeData(wb, "Inner_Medulla", AKI_Medulla.modKEGG)

addWorksheet(wb,"Outer_Medulla")
writeData(wb, "Outer_Medulla", AKI_Interface.modKEGG)

addWorksheet(wb,"Cortex")
writeData(wb, "Cortex", AKI_Cortex.modKEGG)

file<-paste0(dir,"/Codes2_Tables/AKI_KEGG_Pathways_Consistency_R2_Percentile95.xlsx")
saveWorkbook(wb, file, overwrite = TRUE)


wb1<-createWorkbook()

addWorksheet(wb1,"Inner_Medulla")
writeData(wb1, "Inner_Medulla", aki_medulla_up)

addWorksheet(wb1,"Outer_Medulla")
writeData(wb1, "Outer_Medulla", aki_interface_up)

addWorksheet(wb1,"Cortex")
writeData(wb1, "Cortex", aki_cortex_up)

file<-paste0(dir,"/Codes2_Tables/AKI_Up_DEGs_Consistency_R2_Percentile95.xlsx")
saveWorkbook(wb1, file, overwrite = TRUE)


wb2<-createWorkbook()

addWorksheet(wb2,"Inner_Medulla")
writeData(wb2, "Inner_Medulla", aki_medulla_down)

addWorksheet(wb2,"Outer_Medulla")
writeData(wb2, "Outer_Medulla", aki_interface_down)

addWorksheet(wb2,"Cortex")
writeData(wb2, "Cortex", aki_cortex_down)

file<-paste0(dir,"/Codes2_Tables/AKI_Down_DEGs_Consistency_R2_Percentile95.xlsx")
saveWorkbook(wb2, file, overwrite = TRUE)

wb3<-createWorkbook()

addWorksheet(wb3,"Inner_Medulla")
writeData(wb3, "Inner_Medulla", aki_medulla_genes)

addWorksheet(wb3,"Outer_Medulla")
writeData(wb3, "Outer_Medulla", aki_interface_genes)

addWorksheet(wb3,"Cortex")
writeData(wb3, "Cortex", aki_cortex_genes)

file<-paste0(dir,"/Codes2_Tables/AKI_All_DEGs_Rank_Consistency_R2_Percentile95.xlsx")
saveWorkbook(wb3, file, overwrite = TRUE)
