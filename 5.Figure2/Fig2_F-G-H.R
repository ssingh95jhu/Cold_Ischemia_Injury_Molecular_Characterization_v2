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
#### Extracting linear regression results from the stored tables:
## DEGS DETECTED IN DIFFERENT COMPARTMENTS OF THE KIDNEY
cis_cortex_up <- read.xlsx("Codes2_Tables/CIS_Up_DEGs_Consistency_R2_Percentile95.xlsx", sheet = "Cortex")
cis_cortex_down<-read.xlsx("Codes2_Tables/CIS_Down_DEGs_Consistency_R2_Percentile95.xlsx", sheet = "Cortex")
cis_cortex_degs<-rbind(cis_cortex_up, cis_cortex_down)
rownames(cis_cortex_degs)<-cis_cortex_degs$Gene
head(cis_cortex_degs)
tail(cis_cortex_degs)

cis_interface_up <- read.xlsx("Codes2_Tables/CIS_Up_DEGs_Consistency_R2_Percentile95.xlsx", sheet = "Outer_Medulla")
cis_interface_down<-read.xlsx("Codes2_Tables/CIS_Down_DEGs_Consistency_R2_Percentile95.xlsx", sheet = "Outer_Medulla")
cis_interface_degs<-rbind(cis_interface_up, cis_interface_down)
rownames(cis_interface_degs)<-cis_interface_degs$Gene
head(cis_interface_degs)
tail(cis_interface_degs)

cis_medulla_up <- read.xlsx("Codes2_Tables/CIS_Up_DEGs_Consistency_R2_Percentile95.xlsx", sheet = "Inner_Medulla")
cis_medulla_down<-read.xlsx("Codes2_Tables/CIS_Down_DEGs_Consistency_R2_Percentile95.xlsx", sheet = "Inner_Medulla")
cis_medulla_degs<-rbind(cis_medulla_up, cis_medulla_down)
rownames(cis_medulla_degs)<-cis_medulla_degs$Gene
head(cis_medulla_degs)
tail(cis_medulla_degs)


######## VISUALIZING COMMON GENES BETWEEN CIS COMPARTMENTS #####################
wb1<-createWorkbook()

####### 1.CIS CORTEX GENES vs CIS INTERFACE GENES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
comm_genes<-intersect(rownames(cis_cortex_degs), rownames(cis_interface_degs) )
df<-data.frame(CIS_CORTEX_DEGS=cis_cortex_degs[comm_genes,]$Slope, CIS_INTERFACE_DEGS=cis_interface_degs[comm_genes,]$Slope)
df$gene <- comm_genes
head(df)


df$quadrant <- "green"
df[which(df$CIS_CORTEX_DEGS*df$CIS_INTERFACE_DEGS<0),]$quadrant<-"darkviolet"

#1.a. Pearson Correlation
Cor_Int_Corr<-cor.test(x=df$CIS_CORTEX_DEGS, y=df$CIS_INTERFACE_DEGS, method="pearson" )
Cor_Int_Corr <- data.frame(
  Statistic = Cor_Int_Corr$statistic,
  Correlation = Cor_Int_Corr$estimate,
  p_value = Cor_Int_Corr$p.value,
  Conf_Lower = Cor_Int_Corr$conf.int[1],
  Conf_Upper = Cor_Int_Corr$conf.int[2]
)

df_covary<-df[which(df$CIS_CORTEX_DEGS*df$CIS_INTERFACE_DEGS>0),]
df_covary<-data.frame(Gene=df_covary$gene, Cortex_Slope=df_covary$CIS_CORTEX_DEGS, Interface_Slope=df_covary$CIS_INTERFACE_DEGS)
head(df_covary)

df_disjoint<-df[which(df$CIS_CORTEX_DEGS*df$CIS_INTERFACE_DEGS<0),]
df_disjoint<-data.frame(Gene=df_disjoint$gene, Cortex_Slope=df_disjoint$CIS_CORTEX_DEGS, Interface_Slope=df_disjoint$CIS_INTERFACE_DEGS)
head(df_disjoint)

addWorksheet(wb1, "Cov_Cor_Int")
writeData(wb1, "Cov_Cor_Int", df_covary)

addWorksheet(wb1, "Dis_Cor_Int")
writeData(wb1, "Dis_Cor_Int", df_disjoint)

addWorksheet(wb1, "Corr_Cor_Int")
writeData(wb1, "Corr_Cor_Int", Cor_Int_Corr)

highlight_list<-c("Fth1", "Umod", "Igfbp7","Spp1", "Kap")

model <- lm(CIS_INTERFACE_DEGS ~ CIS_CORTEX_DEGS, data = df)
summary<-summary(model)
R2<-summary$r.squared # Look for R-squared
p<-summary$fstatistic
pvalue<-pf(p[1], p[2], p[3], lower.tail = FALSE)
print(pvalue)

CIS_Cor_Int_Fit<-data.frame(R2=R2, pvalue=pvalue)
print(CIS_Cor_Int_Fit)

addWorksheet(wb1, "R2_Cor_Int")
writeData(wb1, "R2_Cor_Int", CIS_Cor_Int_Fit)

slope <- coef(model)[2]  # Extract slope
intercept <- coef(model)[1]  # Extract intercept

file<-paste0(dir,"/Codes2_Figures/Figure2/pdfs/CIS_DEGs_Cortex_vs_Interface_9000x9000.pdf")
pdf(file, height=3, width=5)
ggplot(df, aes(x = CIS_CORTEX_DEGS, y = CIS_INTERFACE_DEGS)) +
  # geom_point(alpha = 0.2) +
  geom_point(aes(color=quadrant),alpha=0.8, shape=16) +
  geom_rect(aes(xmin = -200, xmax = 200, ymin = -200, ymax = 200), 
            fill = NA, color = 'red', linetype = 'solid', size = 0.5) +
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) +
  geom_abline(slope = slope, intercept = intercept, color = "blue", linewidth = 0.5) +  # Pearson correlation line
  geom_text_repel(
    data = \(d) d %>% filter(abs(CIS_CORTEX_DEGS) > 500 | abs(CIS_INTERFACE_DEGS) > 500),
    mapping = aes(label = gene, color = ifelse(gene %in% highlight_list, "red", "black")), 
    max.overlaps = 30, force = 10, box.padding = 0.5) +  
  scale_x_continuous(limits = c(-6000, 3000)) +
  scale_y_continuous(limits = c(-6000, 3000)) + 
  scale_color_identity() + theme_minimal() +
  theme(axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        axis.ticks = element_line(color = "black"), axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(color = "black", fill = NA, size = 1))
dev.off()

#Inset plot for the above plot
file<-paste0(dir,"/Codes2_Figures/Figure2/pdfs/CIS_DEGs_Cortex_vs_Interface_400x400.pdf")
pdf(file, height=5, width=5)
ggplot(df, aes(x = CIS_CORTEX_DEGS, y = CIS_INTERFACE_DEGS)) +
  # geom_point(alpha=0.2) +
  geom_point(aes(color=quadrant),alpha=0.8, shape=16) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  geom_abline(slope = slope, intercept = intercept, color = "blue", linewidth = 0.5) +  # Pearson correlation line
  geom_text_repel(
    data = \(d) d %>% filter((CIS_CORTEX_DEGS) <100 | (CIS_INTERFACE_DEGS) <100),
    mapping = aes(label = gene), max.overlaps = 30, box.padding=0.5) + 
  scale_x_continuous(limits = c(-200, 200)) +
  scale_y_continuous(limits = c(-200, 200))  + 
  scale_color_identity() + theme_minimal() +
  theme(axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        axis.ticks = element_line(color = "black"), axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(color = "red", fill = NA, size = 1))
dev.off()

####### CIS CORTEX GENES vs CIS MEDULLA GENES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
comm_genes<-intersect(rownames(cis_cortex_degs), rownames(cis_medulla_degs) )
df<-data.frame(CIS_CORTEX_DEGS=cis_cortex_degs[comm_genes,]$Slope, CIS_MEDULLA_DEGS=cis_medulla_degs[comm_genes,]$Slope)
df$gene <- comm_genes
head(df)

df$quadrant <- "green"
df[which(df$CIS_CORTEX_DEGS*df$CIS_MEDULLA_DEGS<0),]$quadrant<-"darkviolet"

Cor_Med_Corr<-cor.test(x=df$CIS_CORTEX_DEGS, y=df$CIS_MEDULLA_DEGS, method="pearson" )
Cor_Med_Corr<- data.frame(
  Statistic = Cor_Med_Corr$statistic,
  Correlation = Cor_Med_Corr$estimate,
  p_value = Cor_Med_Corr$p.value,
  Conf_Lower = Cor_Med_Corr$conf.int[1],
  Conf_Upper = Cor_Med_Corr$conf.int[2]
)

df_covary<-df[which(df$CIS_CORTEX_DEGS*df$CIS_MEDULLA_DEGS>0),]
df_covary<-data.frame(Gene=df_covary$gene, Cortex_Slope=df_covary$CIS_CORTEX_DEGS, Medulla_Slope=df_covary$CIS_MEDULLA_DEGS)
head(df_covary)

df_disjoint<-df[which(df$CIS_CORTEX_DEGS*df$CIS_MEDULLA_DEGS<0),]
df_disjoint<-data.frame(Gene=df_disjoint$gene, Cortex_Slope=df_disjoint$CIS_CORTEX_DEGS, Medulla_Slope=df_disjoint$CIS_MEDULLA_DEGS)
head(df_disjoint)

addWorksheet(wb1, "Cov_Cor_Med")
writeData(wb1, "Cov_Cor_Med", df_covary)

addWorksheet(wb1, "Dis_Cor_Med")
writeData(wb1, "Dis_Cor_Med", df_disjoint)

addWorksheet(wb1, "Corr_Cor_Med")
writeData(wb1, "Corr_Cor_Med", Cor_Med_Corr)

highlight_list<-c("Fth1", "Umod", "Igfbp7","Spp1", "Kap")

model <- lm(CIS_MEDULLA_DEGS ~ CIS_CORTEX_DEGS, data = df)
summary<-summary(model)
R2<-summary$r.squared # Look for R-squared
p<-summary$fstatistic
pvalue<-pf(p[1], p[2], p[3], lower.tail = FALSE)
print(pvalue)

CIS_Cor_Med_Fit<-data.frame(R2=R2, pvalue=pvalue)
print(CIS_Cor_Med_Fit)

addWorksheet(wb1, "R2_Cor_Med")
writeData(wb1, "R2_Cor_Med", CIS_Cor_Med_Fit)

slope <- coef(model)[2]  # Extract slope
intercept <- coef(model)[1]  # Extract intercept

file<-paste0(dir,"/Codes2_Figures/Figure2/pdfs/CIS_DEGs_Cortex_vs_Medulla_9000x9000.pdf")
pdf(file, height=3, width=5)
ggplot(df, aes(x = CIS_CORTEX_DEGS, y = CIS_MEDULLA_DEGS)) +
  # geom_point(alpha = 0.2) +
  geom_point(aes(color=quadrant),alpha=0.8, shape=16) +
  geom_rect(aes(xmin = -200, xmax = 200, ymin = -200, ymax = 200), 
            fill = NA, color = 'red', linetype = 'solid', size = 0.5) +
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) +
  geom_abline(slope = slope, intercept = intercept, color = "blue", linewidth = 0.5) +  # Pearson correlation line
  geom_text_repel(
    data = \(d) d %>% filter(abs(CIS_CORTEX_DEGS) > 500 | abs(CIS_MEDULLA_DEGS) > 500),
    mapping = aes(label = gene, color = ifelse(gene %in% highlight_list, "red", "black")), 
    max.overlaps = 30, force = 10, box.padding = 0.5) +  
  scale_x_continuous(limits = c(-6000, 3000)) +
  scale_y_continuous(limits = c(-6000, 3000)) + 
  scale_color_identity() + theme_minimal() +
  theme(axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        axis.ticks = element_line(color = "black"), axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(color = "black", fill = NA, size = 1))
dev.off()

file<-paste0(dir,"/Codes2_Figures/Figure2/pdfs/CIS_DEGs_Cortex_vs_Medulla_400x400.pdf")
pdf(file, height=5, width=5)
ggplot(df, aes(x = CIS_CORTEX_DEGS, y = CIS_MEDULLA_DEGS)) +
  # geom_point(alpha=0.2) +
  geom_point(aes(color=quadrant),alpha=0.8, shape=16) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  geom_abline(slope = slope, intercept = intercept, color = "blue", linewidth = 0.5) +  # Pearson correlation line
  geom_text_repel(
    data = \(d) d %>% filter((CIS_CORTEX_DEGS) <100 | (CIS_MEDULLA_DEGS) <100),
    mapping = aes(label = gene), max.overlaps = 20, box.padding=0.5) + 
  scale_x_continuous(limits = c(-200, 200)) +
  scale_y_continuous(limits = c(-200, 200)) + 
  scale_color_identity() + theme_minimal() +
  theme(axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        axis.ticks = element_line(color = "black"), axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(color = "red", fill = NA, size = 1))
dev.off()

####### CIS INTERFACE GENES vs CIS MEDULLA GENES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
comm_genes<-intersect(rownames(cis_interface_degs), rownames(cis_medulla_degs) )
df<-data.frame(CIS_INTERFACE_DEGS=cis_interface_degs[comm_genes,]$Slope, CIS_MEDULLA_DEGS=cis_medulla_degs[comm_genes,]$Slope)
df$gene <- comm_genes
head(df)

df$quadrant <- "green"
df[which(df$CIS_INTERFACE_DEGS*df$CIS_MEDULLA_DEGS<0),]$quadrant<-"darkviolet"

Int_Med_Corr<-cor.test(x=df$CIS_INTERFACE_DEGS, y=df$CIS_MEDULLA_DEGS, method="pearson" )
Int_Med_Corr<- data.frame(
  Statistic = Int_Med_Corr$statistic,
  Correlation = Int_Med_Corr$estimate,
  p_value = Int_Med_Corr$p.value,
  Conf_Lower = Int_Med_Corr$conf.int[1],
  Conf_Upper = Int_Med_Corr$conf.int[2]
)

df_covary<-df[which(df$CIS_INTERFACE_DEGS*df$CIS_MEDULLA_DEGS>0),]
df_covary<-data.frame(Gene=df_covary$gene, Interface_Slope=df_covary$CIS_INTERFACE_DEGS, Medulla_Slope=df_covary$CIS_MEDULLA_DEGS)
head(df_covary)

df_disjoint<-df[which(df$CIS_INTERFACE_DEGS*df$CIS_MEDULLA_DEGS<0),]
df_disjoint<-data.frame(Gene=df_disjoint$gene, Interface_Slope=df_disjoint$CIS_INTERFACE_DEGS, Medulla_Slope=df_disjoint$CIS_MEDULLA_DEGS)
head(df_disjoint)

addWorksheet(wb1, "Cov_Int_Med")
writeData(wb1, "Cov_Int_Med", df_covary)

addWorksheet(wb1, "Dis_Int_Med")
writeData(wb1, "Dis_Int_Med", df_disjoint)

addWorksheet(wb1, "Corr_Int_Med")
writeData(wb1, "Corr_Int_Med", Int_Med_Corr)


highlight_list<-c("Fth1", "Umod", "Igfbp7","Spp1", "Kap")

model <- lm(CIS_MEDULLA_DEGS ~ CIS_INTERFACE_DEGS, data = df)
summary<-summary(model)
R2<-summary$r.squared # Look for R-squared
p<-summary$fstatistic
pvalue<-pf(p[1], p[2], p[3], lower.tail = FALSE)
print(pvalue)

CIS_Int_Med_Fit<-data.frame(R2=R2, pvalue=pvalue)
print(CIS_Int_Med_Fit)

addWorksheet(wb1, "R2_Int_Med")
writeData(wb1, "R2_Int_Med", CIS_Int_Med_Fit)

slope <- coef(model)[2]  # Extract slope
intercept <- coef(model)[1]  # Extract intercept

file<-paste0(dir,"/Codes2_Figures/Figure2/pdfs/CIS_DEGs_Interface_vs_Medulla_9000x9000.pdf")
pdf(file, height=3, width=5)
ggplot(df, aes(x = CIS_INTERFACE_DEGS, y = CIS_MEDULLA_DEGS)) +
  #geom_point(alpha = 0.2) +
  geom_point(aes(color=quadrant),alpha=0.8, shape=16) +
  geom_rect(aes(xmin = -200, xmax = 200, ymin = -200, ymax = 200), 
            fill = NA, color = 'red', linetype = 'solid', size = 0.5) +
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) +
  geom_abline(slope = slope, intercept = intercept, color = "blue", linewidth = 0.5) +  # Pearson correlation line
  geom_text_repel(
    data = \(d) d %>% filter(abs(CIS_INTERFACE_DEGS) > 500 | abs(CIS_MEDULLA_DEGS) > 500),
    mapping = aes(label = gene, color = ifelse(gene %in% highlight_list, "red", "black")), 
    max.overlaps = 30, force = 10, box.padding = 0.5
  ) +  
  scale_x_continuous(limits = c(-6000, 3000)) +
  scale_y_continuous(limits = c(-6000, 3000)) + 
  scale_color_identity() + 
  scale_color_identity() + theme_minimal() +
  theme(axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        axis.ticks = element_line(color = "black"), axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(color = "black", fill = NA, size = 1))
dev.off()

file<-paste0(dir,"/Codes2_Figures/Figure2/pdfs/CIS_DEGs_Interface_vs_Medulla_400x400.pdf")
pdf(file, height=5, width=5)
ggplot(df, aes(x = CIS_INTERFACE_DEGS, y = CIS_MEDULLA_DEGS)) +
  #geom_point(alpha=0.2) +
  geom_point(aes(color=quadrant),alpha=0.8, shape=16) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  geom_abline(slope = slope, intercept = intercept, color = "blue", linewidth = 0.5) +  # Pearson correlation line
  geom_text_repel(
    data = \(d) d %>% filter((CIS_INTERFACE_DEGS) <100 | (CIS_MEDULLA_DEGS) <100),
    mapping = aes(label = gene), max.overlaps = 20, box.padding=0.5) + 
  scale_x_continuous(limits = c(-200, 200)) +
  scale_y_continuous(limits = c(-200, 200)) + 
  scale_color_identity() + theme_minimal() +
  theme(axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        axis.ticks = element_line(color = "black"), axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(color = "red", fill = NA, size = 1))
dev.off()

saveWorkbook(wb1,"Codes2_Tables/Manuscript/Fig2_CIS_Compartmental_Covary_Disjoint_DEGs.xlsx", overwrite=TRUE)
