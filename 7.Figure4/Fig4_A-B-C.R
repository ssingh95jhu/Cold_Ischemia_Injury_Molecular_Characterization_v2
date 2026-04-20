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

################################################################################
#### Extracting data from the stored tables:
## DEGS DETECTED IN DIFFERENT COMPARTMENTS OF THE KIDNEY (CIS & AKI DATASETS)
cis_cortex_up <- read.xlsx("Codes2_Tables/Manuscript/ST2_CIS_Up_DEGs_Consistency_R2_Percentile95.xlsx", sheet = "Cortex")
cis_cortex_down<-read.xlsx("Codes2_Tables/Manuscript/ST1_CIS_Down_DEGs_Consistency_R2_Percentile95.xlsx", sheet = "Cortex")
cis_cortex_degs<-rbind(cis_cortex_up, cis_cortex_down)
rownames(cis_cortex_degs)<-cis_cortex_degs$Gene
head(cis_cortex_degs)
tail(cis_cortex_degs)

cis_interface_up <- read.xlsx("Codes2_Tables/Manuscript/ST2_CIS_Up_DEGs_Consistency_R2_Percentile95.xlsx", sheet = "Outer_Medulla")
cis_interface_down<-read.xlsx("Codes2_Tables/Manuscript/ST1_CIS_Down_DEGs_Consistency_R2_Percentile95.xlsx", sheet = "Outer_Medulla")
cis_interface_degs<-rbind(cis_interface_up, cis_interface_down)
rownames(cis_interface_degs)<-cis_interface_degs$Gene
head(cis_interface_degs)
tail(cis_interface_degs)

cis_medulla_up <- read.xlsx("Codes2_Tables/Manuscript/ST2_CIS_Up_DEGs_Consistency_R2_Percentile95.xlsx", sheet = "Inner_Medulla")
cis_medulla_down<-read.xlsx("Codes2_Tables/Manuscript/ST1_CIS_Down_DEGs_Consistency_R2_Percentile95.xlsx", sheet = "Inner_Medulla")
cis_medulla_degs<-rbind(cis_medulla_up, cis_medulla_down)
rownames(cis_medulla_degs)<-cis_medulla_degs$Gene
head(cis_medulla_degs)
tail(cis_medulla_degs)

aki_cortex_up<-read.xlsx("Codes2_Tables/AKI_Up_DEGs_Consistency_R2_Percentile95.xlsx", sheet = "Cortex")
aki_cortex_down<-read.xlsx("Codes2_Tables/AKI_Down_DEGs_Consistency_R2_Percentile95.xlsx", sheet = "Cortex")
aki_cortex_degs<-rbind(aki_cortex_up, aki_cortex_down)
rownames(aki_cortex_degs)<-aki_cortex_degs$Gene
head(aki_cortex_degs)
tail(aki_cortex_degs)

aki_interface_up<-read.xlsx("Codes2_Tables/AKI_Up_DEGs_Consistency_R2_Percentile95.xlsx", sheet = "Outer_Medulla")
aki_interface_down<-read.xlsx("Codes2_Tables/AKI_Down_DEGs_Consistency_R2_Percentile95.xlsx", sheet = "Outer_Medulla")
aki_interface_degs<-rbind(aki_interface_up, aki_interface_down)
rownames(aki_interface_degs)<-aki_interface_degs$Gene
head(aki_interface_degs)
tail(aki_interface_degs)

aki_medulla_up<-read.xlsx("Codes2_Tables/AKI_Up_DEGs_Consistency_R2_Percentile95.xlsx", sheet = "Inner_Medulla")
aki_medulla_down<-read.xlsx("Codes2_Tables/AKI_Down_DEGs_Consistency_R2_Percentile95.xlsx", sheet = "Inner_Medulla")
aki_medulla_degs<-rbind(aki_medulla_up, aki_medulla_down)
rownames(aki_medulla_degs)<-aki_medulla_degs$Gene
head(aki_medulla_degs)
tail(aki_medulla_degs)

wb1<-createWorkbook()

## 1.CIS CORTEX vs AKI CORTEX ****************************************************
comm_genes<-intersect(rownames(cis_cortex_degs), rownames(aki_cortex_degs) )
df<-data.frame(CIS_CORTEX_DEGS=cis_cortex_degs[comm_genes,]$Slope, AKI_CORTEX_DEGS=aki_cortex_degs[comm_genes,]$Slope)
df$gene <- comm_genes
head(df)

#1.a.Pearson Correlation:
CIS_AKI_Cor_Corr<-cor.test(x=df$CIS_CORTEX_DEGS, y=df$AKI_CORTEX_DEGS, method="pearson" )
CIS_AKI_Cor_Corr <- data.frame(
  Statistic = CIS_AKI_Cor_Corr$statistic,
  Correlation = CIS_AKI_Cor_Corr$estimate,
  p_value = CIS_AKI_Cor_Corr$p.value,
  Conf_Lower = CIS_AKI_Cor_Corr$conf.int[1],
  Conf_Upper = CIS_AKI_Cor_Corr$conf.int[2]
)

#1.b.Goodness of fit:
model <- lm(AKI_CORTEX_DEGS ~ CIS_CORTEX_DEGS, data = df)
summary<-summary(model)
R2<-summary$r.squared # Look for R-squared
p<-summary$fstatistic
pvalue<-pf(p[1], p[2], p[3], lower.tail = FALSE)
print(pvalue)

CIS_AKI_Cor_Fit<-data.frame(R2=R2, pvalue=pvalue)
print(CIS_AKI_Cor_Fit)

intercept<-summary$coefficients[1,1]
slope<-summary$coefficients[2,1]

highlight_list<-c("Spp1")

df$quadrant <- ifelse((df$CIS_CORTEX_DEGS > 0 & df$AKI_CORTEX_DEGS > 0) | 
                        (df$CIS_CORTEX_DEGS < 0 & df$AKI_CORTEX_DEGS < 0), 
                      "green", "darkviolet")

file<-paste0(dir,"/Codes2_Figures/Figure4/pdfs/CISvsAKI_DEGs_Cortex_9000x32000.pdf")
pdf(file, height=3.5, width=5)
ggplot(df, aes(x = CIS_CORTEX_DEGS, y = AKI_CORTEX_DEGS)) +
  geom_point(aes(color=quadrant),alpha=0.8, shape=16) +
  geom_rect(aes(xmin=-250, xmax=250, ymin=-480, ymax=480), fill=NA, color='red', linetype='solid', size=0.5) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  geom_abline(slope = slope, intercept = intercept, color = "blue", linewidth = 0.5) + # Pearson correlation line
  geom_text_repel(
    data = \(d) d %>% filter(abs(CIS_CORTEX_DEGS) > 500 | abs(AKI_CORTEX_DEGS) > 500),
    mapping = aes(label = gene, color = ifelse(gene %in% highlight_list, "red", "black")), 
    max.overlaps = 30, force=10, box.padding=0.5) + 
  scale_x_continuous(limits = c(-6000, 3000)) +
  scale_y_continuous(limits = c(-21000, 11000)) +
  scale_color_identity() + theme_minimal() +
  theme(axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        axis.ticks = element_line(color = "black"), axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) +
  labs(x = "Regression Slopes (CIS Cortex)",  y = "Regression Slopes (AKI Cortex)", 
       title = "CIS Cortex vs AKI Cortex")
dev.off()

file<-paste0(dir,"/Codes2_Figures/Figure4/pdfs/CISvsAKI_DEGs_Cortex_500x960.pdf")
pdf(file, height=5, width=5)
ggplot(df, aes(x = CIS_CORTEX_DEGS, y = AKI_CORTEX_DEGS)) +
  #geom_point(alpha=0.2) +
  geom_point(aes(color=quadrant),alpha=0.8, shape=16) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  geom_abline(slope = slope, intercept = intercept, color = "blue", linewidth = 0.5) +  # Pearson correlation line
  geom_text_repel(
    data = \(d) d %>% filter(abs(CIS_CORTEX_DEGS) >50 | abs(AKI_CORTEX_DEGS) >50),
    mapping = aes(label = gene), max.overlaps = 15, box.padding=0.5) + 
  scale_x_continuous(limits = c(-250, 250)) +
  scale_y_continuous(limits = c(-480, 480)) + 
  scale_color_identity() + theme_minimal() +
  theme(axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        axis.ticks = element_line(color = "black"), axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(color = "red", fill = NA, size = 1)) +
  labs(x = "Regression Slopes (CIS Cortex)", y = "Regression Slopes (AKI Cortex)")
dev.off()

cortex_covary<-df[which((df$CIS_CORTEX_DEGS*df$AKI_CORTEX_DEGS)>0),]
cortex_covary<-cbind(cortex_covary[,3],cortex_covary[,1:2])
colnames(cortex_covary)<-c('Gene', 'CIS_Slope', 'AKI_Slope')
head(cortex_covary)

cortex_disjoint<-df[which((df$CIS_CORTEX_DEGS*df$AKI_CORTEX_DEGS)<0),]
cortex_disjoint<-cbind(cortex_disjoint[,3],cortex_disjoint[,1:2])
colnames(cortex_disjoint)<-c('Gene', 'CIS_Slope', 'AKI_Slope')
head(cortex_disjoint)

addWorksheet(wb1,"Cortex_Covary")
addWorksheet(wb1,"Cortex_Disjoint")
writeData(wb1, "Cortex_Covary", cortex_covary)
writeData(wb1, "Cortex_Disjoint", cortex_disjoint)

addWorksheet(wb1, "Corr_CISvsAKI_Cortex")
writeData(wb1, "Corr_CISvsAKI_Cortex", CIS_AKI_Cor_Corr)

addWorksheet(wb1, "R2_CISvsAKI_Cortex")
writeData(wb1, "R2_CISvsAKI_Cortex", CIS_AKI_Cor_Fit)

## 2.CIS INTERFACE vs AKI INTERFACE **********************************************
comm_genes<-intersect(rownames(cis_interface_degs), rownames(aki_interface_degs) )
df<-data.frame(CIS_INTERFACE_DEGS=cis_interface_degs[comm_genes,]$Slope, AKI_INTERFACE_DEGS=aki_interface_degs[comm_genes,]$Slope)
df$gene <- comm_genes
head(df)

#2.a Pearson Correlation:
CIS_AKI_Int_Corr<-cor.test(x=df$CIS_INTERFACE_DEGS, y=df$AKI_INTERFACE_DEGS, method="pearson" )
CIS_AKI_Int_Corr <- data.frame(
  Statistic = CIS_AKI_Int_Corr$statistic,
  Correlation = CIS_AKI_Int_Corr$estimate,
  p_value = CIS_AKI_Int_Corr$p.value,
  Conf_Lower = CIS_AKI_Int_Corr$conf.int[1],
  Conf_Upper = CIS_AKI_Int_Corr$conf.int[2]
)

#2.b Goodness of Fit:
model <- lm(AKI_INTERFACE_DEGS ~ CIS_INTERFACE_DEGS, data = df)
summary<-summary(model)
R2<-summary$r.squared # Look for R-squared
p<-summary$fstatistic
pvalue<-pf(p[1], p[2], p[3], lower.tail = FALSE)
print(pvalue)

CIS_AKI_Int_Fit<-data.frame(R2=R2, pvalue=pvalue)
print(CIS_AKI_Int_Fit)

intercept<-summary$coefficients[1,1]
slope<-summary$coefficients[2,1]

df$quadrant <- ifelse((df$CIS_INTERFACE_DEGS > 0 & df$AKI_INTERFACE_DEGS > 0) | 
                        (df$CIS_INTERFACE_DEGS < 0 & df$AKI_INTERFACE_DEGS < 0), 
                      "green", "darkviolet")

file<-paste0(dir,"/Codes2_Figures/Figure4/pdfs/CISvsAKI_DEGs_Interface_9000x32000.pdf")
pdf(file, height=3.5, width=5)
ggplot(df, aes(x = CIS_INTERFACE_DEGS, y = AKI_INTERFACE_DEGS)) +
  #geom_point(alpha=0.2) +
  geom_point(aes(color=quadrant),alpha=0.8, shape=16) +
  geom_rect(aes(xmin=-250, xmax=250, ymin=-480, ymax=480), fill=NA, color='red', linetype='solid', size=0.5) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  geom_abline(slope = slope, intercept = intercept, color = "blue", linewidth = 0.5) +  # Pearson correlation line
  geom_text_repel(
    data = \(d) d %>% filter(abs(CIS_INTERFACE_DEGS) > 500 | abs(AKI_INTERFACE_DEGS) > 500),
    mapping = aes(label = gene, color = ifelse(gene %in% highlight_list, "red", "black")), 
    max.overlaps = 30, force=10, box.padding=0.5) + 
  scale_x_continuous(limits = c(-6000, 3000)) +
  scale_y_continuous(limits = c(-21000, 11000)) + 
  scale_color_identity() + theme_minimal() +
  theme(axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        axis.ticks = element_line(color = "black"), axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) +
  labs(x = "Regression Slopes (CIS Outer Medulla)",  y = "Regression Slopes (AKI Outer Medulla)", 
       title = "CIS Outer Medulla vs AKI Outer Medulla")
dev.off()

file<-paste0(dir,"/Codes2_Figures/Figure4/pdfs/CISvsAKI_DEGs_Interface_500x960.pdf")
pdf(file, height=5, width=5)
ggplot(df, aes(x = CIS_INTERFACE_DEGS, y = AKI_INTERFACE_DEGS)) +
  #geom_point(alpha=0.2) +
  geom_point(aes(color=quadrant),alpha=0.8, shape=16) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  geom_abline(slope = slope, intercept = intercept, color = "blue", linewidth = 0.5) +
  geom_text_repel(
    data = \(d) d %>% filter(abs(CIS_INTERFACE_DEGS) >50 | abs(AKI_INTERFACE_DEGS) >50),
    mapping = aes(label = gene), max.overlaps = 15, box.padding=0.5) + 
  scale_x_continuous(limits = c(-250, 250)) +
  scale_y_continuous(limits = c(-480, 480)) + 
  scale_color_identity() + theme_minimal() +
  theme(axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        axis.ticks = element_line(color = "black"), axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(color = "red", fill = NA, size = 1))+
  labs(x = "Regression Slopes (CIS Outer Medulla)",  y = "Regression Slopes (AKI Outer Medulla)")
dev.off()

interface_covary<-df[which((df$CIS_INTERFACE_DEGS*df$AKI_INTERFACE_DEGS)>0),]
interface_covary<-cbind(interface_covary[,3],interface_covary[,1:2])
colnames(interface_covary)<-c('Gene', 'CIS_Slope', 'AKI_Slope')
head(interface_covary)

interface_disjoint<-df[which((df$CIS_INTERFACE_DEGS*df$AKI_INTERFACE_DEGS)<0),]
interface_disjoint<-cbind(interface_disjoint[,3],interface_disjoint[,1:2])
colnames(interface_disjoint)<-c('Gene', 'CIS_Slope', 'AKI_Slope')
head(interface_disjoint)

addWorksheet(wb1,"Interface_Covary")
addWorksheet(wb1,"Interface_Disjoint")
writeData(wb1, "Interface_Covary", interface_covary)
writeData(wb1, "Interface_Disjoint", interface_disjoint)

addWorksheet(wb1, "Corr_CISvsAKI_Interface")
writeData(wb1, "Corr_CISvsAKI_Interface", CIS_AKI_Int_Corr)

addWorksheet(wb1, "R2_CISvsAKI_Interface")
writeData(wb1, "R2_CISvsAKI_Interface", CIS_AKI_Int_Fit)

## 3.CIS MEDULLA vs AKI MEDULLA **************************************************
comm_genes<-intersect(rownames(cis_medulla_degs), rownames(aki_medulla_degs) )
df<-data.frame(CIS_MEDULLA_DEGS=cis_medulla_degs[comm_genes,]$Slope, AKI_MEDULLA_DEGS=aki_medulla_degs[comm_genes,]$Slope)
df$gene <- comm_genes
head(df)

#3.a. Pearson Correlation
CIS_AKI_Med_Corr<-cor.test(x=df$CIS_MEDULLA_DEGS, y=df$AKI_MEDULLA_DEGS, method="pearson" )
CIS_AKI_Med_Corr <- data.frame(
  Statistic = CIS_AKI_Med_Corr$statistic,
  Correlation = CIS_AKI_Med_Corr$estimate,
  p_value = CIS_AKI_Med_Corr$p.value,
  Conf_Lower = CIS_AKI_Med_Corr$conf.int[1],
  Conf_Upper = CIS_AKI_Med_Corr$conf.int[2]
)

#3.b Goodness of Fit:
model <- lm(AKI_MEDULLA_DEGS ~ CIS_MEDULLA_DEGS, data = df)
summary<-summary(model)
R2<-summary$r.squared # Look for R-squared
p<-summary$fstatistic
pvalue<-pf(p[1], p[2], p[3], lower.tail = FALSE)
print(pvalue)

CIS_AKI_Med_Fit<-data.frame(R2=R2, pvalue=pvalue)
print(CIS_AKI_Med_Fit)

intercept<-summary$coefficients[1,1]
slope<-summary$coefficients[2,1]

df$quadrant <- ifelse((df$CIS_MEDULLA_DEGS > 0 & df$AKI_MEDULLA_DEGS > 0) | 
                        (df$CIS_MEDULLA_DEGS < 0 & df$AKI_MEDULLA_DEGS < 0), 
                      "green", "darkviolet")

file<-paste0(dir,"/Codes2_Figures/Figure4/pdfs/CISvsAKI_DEGs_Medulla_9000x32000.pdf")
pdf(file, height=3.5, width=5)
ggplot(df, aes(x = CIS_MEDULLA_DEGS, y = AKI_MEDULLA_DEGS)) +
  #geom_point(alpha=0.2) +
  geom_point(aes(color=quadrant),alpha=0.8, shape=16) +
  geom_rect(aes(xmin=-250, xmax=250, ymin=-480, ymax=480), fill=NA, color='red', linetype='solid', size=0.5) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  geom_abline(slope = slope, intercept = intercept, color = "blue", linewidth = 0.5) +
  geom_text_repel(
    data = \(d) d %>% filter(abs(CIS_MEDULLA_DEGS) >500 | abs(AKI_MEDULLA_DEGS) > 500),
    mapping = aes(label = gene, color = ifelse(gene %in% highlight_list, "red", "black")), 
    max.overlaps = 40, force=10, box.padding=0.5) + 
  scale_x_continuous(limits = c(-6000, 3000)) +
  scale_y_continuous(limits = c(-21000, 11000)) + 
  scale_color_identity() + theme_minimal() +
  theme(axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        axis.ticks = element_line(color = "black"), axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) +
  labs(x = "Regression Slopes (CIS Inner Medulla)",  y = "Regression Slopes (AKI Inner Medulla)", 
       title = "CIS Inner Medulla vs AKI Inner Medulla")
dev.off()

file<-paste0(dir,"/Codes2_Figures/Figure4/pdfs/CISvsAKI_DEGs_Medulla_500x960.pdf")
pdf(file, height=5, width=5)
ggplot(df, aes(x = CIS_MEDULLA_DEGS, y = AKI_MEDULLA_DEGS)) +
  #geom_point(alpha=0.2) +
  geom_point(aes(color=quadrant),alpha=0.8, shape=16) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  geom_abline(slope = slope, intercept = intercept, color = "blue", linewidth = 0.5) +
  geom_text_repel(
    data = \(d) d %>% filter(abs(CIS_MEDULLA_DEGS) >50 | abs(AKI_MEDULLA_DEGS) >50),
    mapping = aes(label = gene), max.overlaps = 20, box.padding=0.5) + 
  scale_x_continuous(limits = c(-250, 250)) +
  scale_y_continuous(limits = c(-480, 480)) + 
  scale_color_identity() + theme_minimal() +
  theme(axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        axis.ticks = element_line(color = "black"), axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(color = "red", fill = NA, size = 1)) +
  labs(x = "Regression Slopes (CIS Inner Medulla)",  y = "Regression Slopes (AKI Inner Medulla)")
dev.off()

medulla_covary<-df[which((df$CIS_MEDULLA_DEGS*df$AKI_MEDULLA_DEGS)>0),]
medulla_covary<-cbind(medulla_covary[,3],medulla_covary[,1:2])
colnames(medulla_covary)<-c('Gene', 'CIS_Slope', 'AKI_Slope')
head(medulla_covary)

medulla_disjoint<-df[which((df$CIS_MEDULLA_DEGS*df$AKI_MEDULLA_DEGS)<0),]
medulla_disjoint<-cbind(medulla_disjoint[,3],medulla_disjoint[,1:2])
colnames(medulla_disjoint)<-c('Gene', 'CIS_Slope', 'AKI_Slope')
head(medulla_disjoint)

addWorksheet(wb1,"Medulla_Covary")
addWorksheet(wb1,"Medulla_Disjoint")
writeData(wb1, "Medulla_Covary", medulla_covary)
writeData(wb1, "Medulla_Disjoint", medulla_disjoint)

addWorksheet(wb1, "Corr_CISvsAKI_Medulla")
writeData(wb1, "Corr_CISvsAKI_Medulla", CIS_AKI_Med_Corr)

addWorksheet(wb1, "R2_CISvsAKI_Medulla")
writeData(wb1, "R2_CISvsAKI_Medulla", CIS_AKI_Med_Fit)

saveWorkbook(wb1,"Codes2_Tables/Manuscript/Fig4_CISvsAKI_Covary_Disjoint_DEGs_Correlation_R2.xlsx", overwrite=TRUE)

