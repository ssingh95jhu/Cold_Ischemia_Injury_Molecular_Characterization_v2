library(ggplot2)
library(ggsci)

#Uqcr11 (Oxidative Phosphorylation)
setwd("/Users/srujansingh/Library/CloudStorage/OneDrive-SharedLibraries-JohnsHopkins/Jean Fan - 2024_spatial_AKI/Revision1/qPCR")
gene1<-'Uqcr11'
filename1<-paste0("/Users/srujansingh/Library/CloudStorage/OneDrive-SharedLibraries-JohnsHopkins/Jean Fan - 2024_spatial_AKI/Revision1/qPCR/",gene1,".csv")
data <- read.csv(filename1)

head(data)

data1<-data[,c(2,3,13)]
head(data1)


model_cortex<-wilcox.test(RQ2 ~ Time, data = subset(data1, Region=="Cortex"))
model_medulla<-wilcox.test(RQ2 ~ Time, data = subset(data1, Region=="Medulla"))

p.value=data.frame(Region=c("Cortex", "Medulla"),
                   pval=c(signif(model_cortex$p.value,2), signif(model_medulla$p.value,2)) )

custom_colors<-c(Cortex="red", Medulla="magenta")

pdf("/Users/srujansingh/Library/CloudStorage/OneDrive-SharedLibraries-JohnsHopkins/Jean Fan - 2024_spatial_AKI/Revision1/Figures/3/Uqcr11_qPCR.pdf", width=6, height=5)
ggplot(data1, aes(x = Time, y = RQ2, fill = Region)) +
  geom_boxplot(alpha = 1, width = 0.6, outlier.shape=NA) +
  scale_fill_manual(values=custom_colors) +
  geom_jitter(alpha=1, position = position_jitter(width = 0.1)) +
  facet_wrap(~Region) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14),     
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),     
        axis.text.y = element_text(size = 12)) +
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=c(0,4))+
  labs(x="Cold Ischemia Time (hours)", y="Relative mRNA Fold Change", title = "Uqcr11") +
geom_text(
  data = p.value,
  aes(x = 1.5, y = Inf, label = paste0("p=",pval)),   
  vjust = 1.5,                            
  size = 4
)
dev.off()

#Cox6a1 (Oxidative Phosphorylation)
gene2<-'Cox6a1'
filename2<-paste0("/Users/srujansingh/Library/CloudStorage/OneDrive-SharedLibraries-JohnsHopkins/Jean Fan - 2024_spatial_AKI/Revision1/qPCR/",gene2,".csv")
data <- read.csv(filename2)

data<-data[which((data$Time=="0h")|(data$Time=="48h")),]
data$Time<-as.numeric(sub( "h","", data$Time))

data1<-data[,c(2,3,13)]
data1$Time<-factor(data1$Time, levels=c(0,48))
head(data1)

model_cortex<-wilcox.test(RQ ~ Time, data = subset(data1, Region=="Cortex"))
model_medulla<-wilcox.test(RQ ~ Time, data = subset(data1, Region=="Medulla"))

p.value=data.frame(Region=c("Cortex", "Medulla"),
                   pval=c(signif(model_cortex$p.value,2), signif(model_medulla$p.value,2)) )

custom_colors<-c(Cortex="red", Medulla="magenta")

pdf("/Users/srujansingh/Library/CloudStorage/OneDrive-SharedLibraries-JohnsHopkins/Jean Fan - 2024_spatial_AKI/Revision1/Figures/3/Cox6a1_qPCR.pdf", width=6, height=5)
ggplot(data1, aes(x = Time, y = RQ, fill = Region)) +
  geom_boxplot(alpha = 1, width = 0.6, outlier.shape=NA) +
  scale_fill_manual(values=custom_colors) +
  geom_jitter(alpha=1, position = position_jitter(width = 0.1)) +
  facet_wrap(~Region) +
  theme_bw() +
  theme (axis.title.x = element_text(size = 14),     
         axis.title.y = element_text(size = 14),
         axis.text.x = element_text(size = 12),     
         axis.text.y = element_text(size = 12)) +
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=c(0,2))+
  labs(x="Cold Ischemia Time (hours)", y="Relative mRNA Fold Change", title = "Cox6a1") +
  geom_text(
    data = p.value,
    aes(x = 1.5, y = Inf, label = paste0("p=", pval)),   # x = middle, y = top
    vjust = 1.5,                            # move slightly below top
    size = 4
  )
dev.off()