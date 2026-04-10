#Clean environment
rm(list=ls(all.names=TRUE))
gc()

setwd("/Users/srujansingh/Library/CloudStorage/OneDrive-SharedLibraries-JohnsHopkins/Jean Fan - 2024_spatial_AKI/Manuscript")

########################## ADDING LIBRARY ######################################
library(openxlsx)
library(ggplot2)

# load coled ischemia kidney data
load("data/CIS_data.RData")

## Cold Ischemia Dataset Quality Control #######################################
#cutsom theme
my_theme<-theme_minimal() + 
  theme(
        axis.text.x = element_text(size = 12, color = "black"),   # Adjust x-axis label size & color
        axis.text.y = element_text(size = 12, color = "black", angle=0),
        axis.title.x = element_text(size = 14, color = "black"),  # Adjust x-axis title size & color
        axis.title.y = element_text(size = 14, color = "black"),   # Adjust y-axis title size & color
        axis.ticks.length = unit(0.2, "cm"), 
        axis.ticks = element_line(size = 0.5),
        plot.title=element_text(hjust=0.5),
        panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add a border
  )

########## CIS Transcripts per Spot  QC ########################################
#CIS 0 hours
df1<-data.frame(Genes_per_Spot=log10(colSums(as.matrix(CIS_0h$gexp)) +1))
g1<-ggplot(df1, aes(x = Genes_per_Spot)) +
  geom_histogram(binwidth = 0.05, fill = "gray", color = "black") +  # Adjust bin width as needed
  labs(title = "CIS_0h", x = "Log10 (Transcripts per spot)", y = "Frequency") +
  scale_x_continuous(breaks = seq(1.5, 5.5, by = 0.25), limits = c(1.5, 5.5)) +  
  scale_y_continuous(breaks = seq(0, 400, by = 100), limits = c(0, 400)) + 
  my_theme

#CIS 12 hours
df2<-data.frame(Genes_per_Spot=log10(colSums(as.matrix(CIS_12h$gexp)) +1))
g2<-ggplot(df2, aes(x = Genes_per_Spot)) +
  geom_histogram(binwidth = 0.05, fill = "gray", color = "black") +  # Adjust bin width as needed
  labs(title = "CIS_12h", x = "Log10 (Transcripts per spot)", y = "Frequency") +
  scale_x_continuous(breaks = seq(1.5, 5.5, by = 0.25), limits = c(1.5, 5.5)) +  
  scale_y_continuous(breaks = seq(0, 400, by = 100), limits = c(0, 400)) + 
  my_theme

#CIS 24 hours
df3<-data.frame(Genes_per_Spot=log10(colSums(as.matrix(CIS_24h$gexp)) +1))
g3<-ggplot(df3, aes(x = Genes_per_Spot)) +
  geom_histogram(binwidth = 0.05, fill = "gray", color = "black") +  # Adjust bin width as needed
  labs(title = "CIS_24h", x = "Log10 (Transcripts per spot)", y = "Frequency") +
  scale_x_continuous(breaks = seq(1.5, 5.5, by = 0.25), limits = c(1.5, 5.5)) +  
  scale_y_continuous(breaks = seq(0, 400, by = 100), limits = c(0, 400)) + 
  my_theme

#CIS 48 hours
df4<-data.frame(Genes_per_Spot=log10(colSums(as.matrix(CIS_48h$gexp)) +1))
g4<-ggplot(df4, aes(x = Genes_per_Spot)) +
  geom_histogram(binwidth = 0.05, fill = "gray", color = "black") +  # Adjust bin width as needed
  labs(title = "CIS_48h", x = "Log10 (Transcripts per spot)", y = "Frequency") +
  scale_x_continuous(breaks = seq(1.5, 5.5, by = 0.25), limits = c(1.5, 5.5)) +  
  scale_y_continuous(breaks = seq(0, 400, by = 100), limits = c(0, 400)) + 
  my_theme

df<-c(df1,df2,df3,df4)
range(df)

library(gridExtra)
pdf("Figures/Figure1S/pdfs/CIS_QC_Transcripts_per_Spot.pdf", height=10, width=7.5)
grid.arrange(g1,g2,g3,g4, ncol=1)
dev.off()


########## CIS Spots per Gene QC ###############################################
#CIS 0 hours
df5<-data.frame(Spots_per_Gene=log10(rowSums(as.matrix(CIS_0h$gexp)) +1))
g5<-ggplot(df5, aes(x = Spots_per_Gene)) +
  geom_histogram(binwidth = 0.05, fill = "gray", color = "black") +  # Adjust bin width as needed
  labs(title = "CIS_0h", x = "Log10 (Spots per gene)", y = "Frequency") +
  scale_x_continuous(breaks = seq(0, 7, by = 1), limits = c(0, 7)) +  
  scale_y_continuous(breaks = seq(0, 800, by = 200), limits = c(0, 800)) + 
  my_theme

#CIS 12 hours
df6<-data.frame(Spots_per_Gene=log10(rowSums(as.matrix(CIS_12h$gexp)) +1))
g6<-ggplot(df6, aes(x = Spots_per_Gene)) +
  geom_histogram(binwidth = 0.05, fill = "gray", color = "black") +  # Adjust bin width as needed
  labs(title = "CIS_12h", x = "Log10 (Spots per gene)", y = "Frequency") +
  scale_x_continuous(breaks = seq(0, 7, by = 1), limits = c(0, 7)) +  
  scale_y_continuous(breaks = seq(0, 800, by = 200), limits = c(0, 800)) + 
  my_theme

#CIS 24 hours
df7<-data.frame(Spots_per_Gene=log10(rowSums(as.matrix(CIS_24h$gexp)) +1))
g7<-ggplot(df7, aes(x = Spots_per_Gene)) +
  geom_histogram(binwidth = 0.05, fill = "gray", color = "black") +  # Adjust bin width as needed
  labs(title = "CIS_24h", x = "Log10 (Spots per gene)", y = "Frequency") +
  scale_x_continuous(breaks = seq(0, 7, by = 1), limits = c(0, 7)) +  
  scale_y_continuous(breaks = seq(0, 800, by = 200), limits = c(0, 800)) + 
  my_theme

#CIS 48 hours
df8<-data.frame(Spots_per_Gene=log10(rowSums(as.matrix(CIS_48h$gexp)) +1))
g8<-ggplot(df8, aes(x = Spots_per_Gene)) +
  geom_histogram(binwidth = 0.05, fill = "gray", color = "black") +  # Adjust bin width as needed
  labs(title = "CIS_48h", x = "Log10 (Spots per gene)", y = "Frequency") +
  scale_x_continuous(breaks = seq(0, 7, by = 1), limits = c(0, 7)) +  
  scale_y_continuous(breaks = seq(0, 800, by = 200), limits = c(0, 800)) + 
  my_theme

library(gridExtra)
pdf("Figures/Figure1S/pdfs/CIS_QC_Spots_per_Gene.pdf", height=10, width=7.5)
grid.arrange(g5,g6,g7,g8, ncol=1)
dev.off()

