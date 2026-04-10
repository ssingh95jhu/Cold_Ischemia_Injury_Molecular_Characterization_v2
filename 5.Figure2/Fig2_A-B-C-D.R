#Clean environment
rm(list=ls(all.names=TRUE))
gc()

dir<-"/Users/srujansingh/Library/CloudStorage/OneDrive-SharedLibraries-JohnsHopkins/Jean Fan - 2024_spatial_AKI/Revision2_GenomBiol"
setwd(dir)

## load data
load("data/CIS_data.RData")

########################## ADDING LIBRARY ######################################
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(patchwork)
library(cowplot)
library(viridis) #color blind friendly
library(RColorBrewer) #for RColorbrewer

################# NORMALIZE CIS DATASET ########################################
#CPM Normalization
CIS_0h$mat_notlog <- MERINGUE::normalizeCounts(CIS_0h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
CIS_12h$mat_notlog <- MERINGUE::normalizeCounts(CIS_12h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
CIS_24h$mat_notlog <- MERINGUE::normalizeCounts(CIS_24h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
CIS_48h$mat_notlog <- MERINGUE::normalizeCounts(CIS_48h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)

# ##############################################################################
# #Extracting the compartment specific spots for all the datasets:
cortex <- readRDS('AKI-CIS-irl-ctrl_Cortex_spots.rds')
interface <- readRDS('AKI-CIS-irl-ctrl_Interface_spots.rds')
medulla <- readRDS('AKI-CIS-irl-ctrl_Medulla_spots.rds')

#1.All_Cortex_Spots (CIS)
CIS_cortex.0h<-names(cortex[grep("^CIS_0h", names(cortex))])
CIS_cortex.12h<-names(cortex[grep("^CIS_12h", names(cortex))])
CIS_cortex.24h<-names(cortex[grep("^CIS_24h", names(cortex))])
CIS_cortex.48h<-names(cortex[grep("^CIS_48h", names(cortex))])

#2.All_Interface_Spots (CIS)
CIS_interface.0h<-names(interface[grep("^CIS_0h", names(interface))])
CIS_interface.12h<-names(interface[grep("^CIS_12h", names(interface))])
CIS_interface.24h<-names(interface[grep("^CIS_24h", names(interface))])
CIS_interface.48h<-names(interface[grep("^CIS_48h", names(interface))])

#3.All_Medulla_Spots (CIS)
CIS_medulla.0h<-names(medulla[grep("^CIS_0h", names(medulla))])
CIS_medulla.12h<-names(medulla[grep("^CIS_12h", names(medulla))])
CIS_medulla.24h<-names(medulla[grep("^CIS_24h", names(medulla))])
CIS_medulla.48h<-names(medulla[grep("^CIS_48h", names(medulla))])


####################### PLOTTING AND VISUALIZATION #############################
## FOR COLD ISCHEMIA TISSUE SPECIMEN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Igfbp7
g<-'Igfbp7' #Insert any gene of interest, for e.g., Igfbp7
overall_min <- min(CIS_0h$mat_notlog[g,], CIS_12h$mat_notlog[g,], 
                   CIS_24h$mat_notlog[g,], CIS_48h$mat_notlog[g,])

overall_max <- quantile(c(CIS_0h$mat_notlog[g,], CIS_12h$mat_notlog[g,], 
                          CIS_24h$mat_notlog[g,], CIS_48h$mat_notlog[g,]), 0.95)

#Set up a custom color scale (as per your choice)
color_scale <- scale_color_gradientn(colors=viridis(5),
                                     #colors=c("grey","magenta","yellow","red"),
                                     values = c(0, 0.8, 1),
                                     limits = c(overall_min, overall_max), oob = scales::squish)

#Set up a custom theme for plotting
my_theme <- theme_bw() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.text = element_text(),
        legend.position="none"
  )

#Spatiotemporal Heatmaps
df <- data.frame(CIS_0h$pos, gexp=CIS_0h$mat_notlog[g,])
g1 <- ggplot(df, aes(x=X, y=Y, col=gexp)) + geom_point(size=0.5, alpha=1) + color_scale + my_theme  + ggtitle('CIS 0h') 
df <- data.frame(CIS_12h$pos, gexp=CIS_12h$mat_notlog[g,])
g2 <- ggplot(df, aes(x=X, y=Y, col=gexp)) + geom_point(size=0.5, alpha=1)  + color_scale  + my_theme + ggtitle('CIS 12h')
df <- data.frame(CIS_24h$pos, gexp=CIS_24h$mat_notlog[g,])
g3 <- ggplot(df, aes(x=X, y=Y, col=gexp)) + geom_point(size=0.5, alpha=1) + color_scale + my_theme + ggtitle('CIS 24h')
df <- data.frame(CIS_48h$pos, gexp=CIS_48h$mat_notlog[g,])
g4 <- ggplot(df, aes(x=X, y=Y, col=gexp)) + geom_point(size=0.5, alpha=1) + color_scale  + theme(legend.text = element_text(size=10), legend.ticks = element_line(color = "black")) + ggtitle('CIS 48h')

# Extract the legend from one plot
legend <- get_legend(g4)
g4 <- ggplot(df, aes(x=X, y=Y, col=gexp)) + geom_point(size=0.5, alpha=1) + color_scale  + my_theme + ggtitle('CIS 48h')

comb_plot<-plot_grid(g1, g2, g3, g4, legend, ncol = 5, rel_widths = c(1, 1, 1, 1, 0.25))

#RibbonPlots
df1 <- data.frame(time = c(rep(0, length(CIS_cortex.0h)), rep(12, length(CIS_cortex.12h)), 
                        rep(24, length(CIS_cortex.24h)), rep(48, length(CIS_cortex.48h))),
                  gexp = c(CIS_0h$mat_notlog[g,CIS_cortex.0h], CIS_12h$mat_notlog[g,CIS_cortex.12h], 
                           CIS_24h$mat_notlog[g,CIS_cortex.24h], CIS_48h$mat_notlog[g,CIS_cortex.48h]),
                  
                  group = 'Cortex')
head(df1)

df2 <- data.frame(time = c(rep(0, length(CIS_interface.0h)), rep(12, length(CIS_interface.12h)), 
                        rep(24, length(CIS_interface.24h)), rep(48, length(CIS_interface.48h))),
                  gexp = c(CIS_0h$mat_notlog[g,CIS_interface.0h], CIS_12h$mat_notlog[g,CIS_interface.12h], 
                           CIS_24h$mat_notlog[g,CIS_interface.24h], CIS_48h$mat_notlog[g,CIS_interface.48h]), group = 'O.Medulla')
head(df2)

df3 <- data.frame(time = c(rep(0, length(CIS_medulla.0h)), rep(12, length(CIS_medulla.12h)), 
                        rep(24, length(CIS_medulla.24h)), rep(48, length(CIS_medulla.48h))),
                  gexp = c(CIS_0h$mat_notlog[g,CIS_medulla.0h], CIS_12h$mat_notlog[g,CIS_medulla.12h], 
                           CIS_24h$mat_notlog[g,CIS_medulla.24h], CIS_48h$mat_notlog[g,CIS_medulla.48h]), group = 'I.Medulla')
head(df3)

df4 <- data.frame(time = c(rep(0, ncol(CIS_0h$mat_notlog)), rep(12, ncol(CIS_12h$mat_notlog)), 
                        rep(24, ncol(CIS_24h$mat_notlog)), rep(48, ncol(CIS_48h$mat_notlog))),
                  gexp = c(CIS_0h$mat_notlog[g,], CIS_12h$mat_notlog[g,], 
                           CIS_24h$mat_notlog[g,], CIS_48h$mat_notlog[g,]), group = 'Global')
head(df4)

df<-rbind(df1, df2, df3, df4)

df_comb <- aggregate(gexp ~ time + group, data = df,
                        FUN = function(x) c(mean = mean(x),
                                            max=max(x),
                                            min=min(x),
                                            sd = sd(x),
                                            se = sd(x)/sqrt(length(x))))

df_comb <- cbind(df_comb[,1:2], df_comb$gexp)
colnames(df_comb) <- c("time","group","mean","max", "min","sd","se")
df_comb$group<-factor(df_comb$group, levels=c("Global", "Cortex", "O.Medulla", "I.Medulla"))


###Ribbon Plots: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 glineSD<-ggplot(df_comb,aes(x=time, y=mean, color=group, fill=group, group=group)) +
   geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd),alpha=0.25, colour=NA) +
  geom_line(size=0.5) + geom_point(size=0.7) + 
  scale_color_manual(values=c("Global"="black","Cortex"="red", 
                              "O.Medulla"="blue", "I.Medulla"="green")) +
  scale_fill_manual(values=c("Global"="black","Cortex"="red", 
                             "O.Medulla"="blue","I.Medulla"="green")) +
   #ylim(overall_min0, overall_max0) + 
   theme_light() +
   theme(axis.text.x = element_text(color="black", size=11),  
         axis.text.y = element_text(color="black",size=11),
         axis.title.x=element_text(size=12),
         axis.title.y=element_text(size=12),
         panel.border = element_rect(color = "black", fill = NA, size = 1),
         legend.text = element_text(color="black", size=12),
         panel.grid.minor = element_blank(),
         axis.ticks = element_line(color = "black"), axis.ticks.length = unit(0.1, "cm")) +
   xlab("time (hours)") + ylab("Avg. gene exp. (CPM)")
 glineSD

file<-paste0(dir,"/Codes2_Figures/Figure2/pdfs/Igfbp7_CIS_AllTimepoints.pdf")  
pdf(file, height=3, width=15)
gridExtra::grid.arrange(comb_plot,glineSD, ncol=2, widths=c(16, 5), top=g)
dev.off()

#Ranbp3l ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
g<-'Ranbp3l' #Insert any gene of interest, for e.g., Igfbp7
overall_min <- min(CIS_0h$mat_notlog[g,], CIS_12h$mat_notlog[g,], 
                   CIS_24h$mat_notlog[g,], CIS_48h$mat_notlog[g,])

overall_max <- quantile(c(CIS_0h$mat_notlog[g,], CIS_12h$mat_notlog[g,], 
                          CIS_24h$mat_notlog[g,], CIS_48h$mat_notlog[g,]), 0.95)

#Set up a custom color scale (as per your choice)
color_scale <- scale_color_gradientn(colors=viridis(5),
                                     #colors=c("grey","magenta","yellow","red"),
                                     values = c(0, 0.8, 1),
                                     limits = c(overall_min, overall_max), oob = scales::squish)

#Set up a custom theme for plotting
my_theme <- theme_bw() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.text = element_text(),
        legend.position="none"
  )

#Spatiotemporal Heatmaps
df <- data.frame(CIS_0h$pos, gexp=CIS_0h$mat_notlog[g,])
g1 <- ggplot(df, aes(x=X, y=Y, col=gexp)) + geom_point(size=0.5, alpha=1) + color_scale + my_theme  + ggtitle('CIS 0h') 
df <- data.frame(CIS_12h$pos, gexp=CIS_12h$mat_notlog[g,])
g2 <- ggplot(df, aes(x=X, y=Y, col=gexp)) + geom_point(size=0.5, alpha=1)  + color_scale  + my_theme + ggtitle('CIS 12h')
df <- data.frame(CIS_24h$pos, gexp=CIS_24h$mat_notlog[g,])
g3 <- ggplot(df, aes(x=X, y=Y, col=gexp)) + geom_point(size=0.5, alpha=1) + color_scale + my_theme + ggtitle('CIS 24h')
df <- data.frame(CIS_48h$pos, gexp=CIS_48h$mat_notlog[g,])
g4 <- ggplot(df, aes(x=X, y=Y, col=gexp)) + geom_point(size=0.5, alpha=1) + color_scale  + theme(legend.text = element_text(size=10), legend.ticks = element_line(color = "black")) + ggtitle('CIS 48h')

# Extract the legend from one plot
legend <- get_legend(g4)
g4 <- ggplot(df, aes(x=X, y=Y, col=gexp)) + geom_point(size=0.5, alpha=1) + color_scale  + my_theme + ggtitle('CIS 48h')

comb_plot<-plot_grid(g1, g2, g3, g4, legend, ncol = 5, rel_widths = c(1, 1, 1, 1, 0.25))

#RibbonPlots
df1 <- data.frame(time = c(rep(0, length(CIS_cortex.0h)), rep(12, length(CIS_cortex.12h)), 
                           rep(24, length(CIS_cortex.24h)), rep(48, length(CIS_cortex.48h))),
                  gexp = c(CIS_0h$mat_notlog[g,CIS_cortex.0h], CIS_12h$mat_notlog[g,CIS_cortex.12h], 
                           CIS_24h$mat_notlog[g,CIS_cortex.24h], CIS_48h$mat_notlog[g,CIS_cortex.48h]),
                  
                  group = 'Cortex')
head(df1)

df2 <- data.frame(time = c(rep(0, length(CIS_interface.0h)), rep(12, length(CIS_interface.12h)), 
                           rep(24, length(CIS_interface.24h)), rep(48, length(CIS_interface.48h))),
                  gexp = c(CIS_0h$mat_notlog[g,CIS_interface.0h], CIS_12h$mat_notlog[g,CIS_interface.12h], 
                           CIS_24h$mat_notlog[g,CIS_interface.24h], CIS_48h$mat_notlog[g,CIS_interface.48h]), group = 'O.Medulla')
head(df2)

df3 <- data.frame(time = c(rep(0, length(CIS_medulla.0h)), rep(12, length(CIS_medulla.12h)), 
                           rep(24, length(CIS_medulla.24h)), rep(48, length(CIS_medulla.48h))),
                  gexp = c(CIS_0h$mat_notlog[g,CIS_medulla.0h], CIS_12h$mat_notlog[g,CIS_medulla.12h], 
                           CIS_24h$mat_notlog[g,CIS_medulla.24h], CIS_48h$mat_notlog[g,CIS_medulla.48h]), group = 'I.Medulla')
head(df3)

df4 <- data.frame(time = c(rep(0, ncol(CIS_0h$mat_notlog)), rep(12, ncol(CIS_12h$mat_notlog)), 
                           rep(24, ncol(CIS_24h$mat_notlog)), rep(48, ncol(CIS_48h$mat_notlog))),
                  gexp = c(CIS_0h$mat_notlog[g,], CIS_12h$mat_notlog[g,], 
                           CIS_24h$mat_notlog[g,], CIS_48h$mat_notlog[g,]), group = 'Global')
head(df4)

df<-rbind(df1, df2, df3, df4)

df_comb <- aggregate(gexp ~ time + group, data = df,
                     FUN = function(x) c(mean = mean(x),
                                         max=max(x),
                                         min=min(x),
                                         sd = sd(x),
                                         se = sd(x)/sqrt(length(x))))

df_comb <- cbind(df_comb[,1:2], df_comb$gexp)
colnames(df_comb) <- c("time","group","mean","max", "min","sd","se")
df_comb$group<-factor(df_comb$group, levels=c("Global", "Cortex", "O.Medulla", "I.Medulla"))


###Ribbon Plots: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
glineSD<-ggplot(df_comb,aes(x=time, y=mean, color=group, fill=group, group=group)) +
  geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd),alpha=0.25, colour=NA) +
  geom_line(size=0.5) + geom_point(size=0.7) + 
  scale_color_manual(values=c("Global"="black","Cortex"="red", 
                              "O.Medulla"="blue", "I.Medulla"="green")) +
  scale_fill_manual(values=c("Global"="black","Cortex"="red", 
                             "O.Medulla"="blue","I.Medulla"="green")) +
  #ylim(overall_min0, overall_max0) + 
  theme_light() +
  theme(axis.text.x = element_text(color="black", size=11),  
        axis.text.y = element_text(color="black",size=11),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.text = element_text(color="black", size=12),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"), axis.ticks.length = unit(0.1, "cm")) +
  xlab("time (hours)") + ylab("Avg. gene exp. (CPM)")
glineSD

file<-paste0(dir,"/Codes2_Figures/Figure2/pdfs/Ranbp3l_CIS_AllTimepoints.pdf")  
pdf(file, height=3, width=15)
gridExtra::grid.arrange(comb_plot,glineSD, ncol=2, widths=c(16, 5), top=g)
dev.off()

#Ttc36 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
g<-'Ttc36' #Insert any gene of interest, for e.g., Igfbp7
overall_min <- min(CIS_0h$mat_notlog[g,], CIS_12h$mat_notlog[g,], 
                   CIS_24h$mat_notlog[g,], CIS_48h$mat_notlog[g,])

overall_max <- quantile(c(CIS_0h$mat_notlog[g,], CIS_12h$mat_notlog[g,], 
                          CIS_24h$mat_notlog[g,], CIS_48h$mat_notlog[g,]), 0.95)

#Set up a custom color scale (as per your choice)
color_scale <- scale_color_gradientn(colors=viridis(5),
                                     #colors=c("grey","magenta","yellow","red"),
                                     values = c(0, 0.8, 1),
                                     limits = c(overall_min, overall_max), oob = scales::squish)

#Set up a custom theme for plotting
my_theme <- theme_bw() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.text = element_text(),
        legend.position="none"
  )

#Spatiotemporal Heatmaps
df <- data.frame(CIS_0h$pos, gexp=CIS_0h$mat_notlog[g,])
g1 <- ggplot(df, aes(x=X, y=Y, col=gexp)) + geom_point(size=0.5, alpha=1) + color_scale + my_theme  + ggtitle('CIS 0h') 
df <- data.frame(CIS_12h$pos, gexp=CIS_12h$mat_notlog[g,])
g2 <- ggplot(df, aes(x=X, y=Y, col=gexp)) + geom_point(size=0.5, alpha=1)  + color_scale  + my_theme + ggtitle('CIS 12h')
df <- data.frame(CIS_24h$pos, gexp=CIS_24h$mat_notlog[g,])
g3 <- ggplot(df, aes(x=X, y=Y, col=gexp)) + geom_point(size=0.5, alpha=1) + color_scale + my_theme + ggtitle('CIS 24h')
df <- data.frame(CIS_48h$pos, gexp=CIS_48h$mat_notlog[g,])
g4 <- ggplot(df, aes(x=X, y=Y, col=gexp)) + geom_point(size=0.5, alpha=1) + color_scale  + theme(legend.text = element_text(size=10), legend.ticks = element_line(color = "black")) + ggtitle('CIS 48h')

# Extract the legend from one plot
legend <- get_legend(g4)
g4 <- ggplot(df, aes(x=X, y=Y, col=gexp)) + geom_point(size=0.5, alpha=1) + color_scale  + my_theme + ggtitle('CIS 48h')

comb_plot<-plot_grid(g1, g2, g3, g4, legend, ncol = 5, rel_widths = c(1, 1, 1, 1, 0.25))

#RibbonPlots
df1 <- data.frame(time = c(rep(0, length(CIS_cortex.0h)), rep(12, length(CIS_cortex.12h)), 
                           rep(24, length(CIS_cortex.24h)), rep(48, length(CIS_cortex.48h))),
                  gexp = c(CIS_0h$mat_notlog[g,CIS_cortex.0h], CIS_12h$mat_notlog[g,CIS_cortex.12h], 
                           CIS_24h$mat_notlog[g,CIS_cortex.24h], CIS_48h$mat_notlog[g,CIS_cortex.48h]),
                  
                  group = 'Cortex')
head(df1)

df2 <- data.frame(time = c(rep(0, length(CIS_interface.0h)), rep(12, length(CIS_interface.12h)), 
                           rep(24, length(CIS_interface.24h)), rep(48, length(CIS_interface.48h))),
                  gexp = c(CIS_0h$mat_notlog[g,CIS_interface.0h], CIS_12h$mat_notlog[g,CIS_interface.12h], 
                           CIS_24h$mat_notlog[g,CIS_interface.24h], CIS_48h$mat_notlog[g,CIS_interface.48h]), group = 'O.Medulla')
head(df2)

df3 <- data.frame(time = c(rep(0, length(CIS_medulla.0h)), rep(12, length(CIS_medulla.12h)), 
                           rep(24, length(CIS_medulla.24h)), rep(48, length(CIS_medulla.48h))),
                  gexp = c(CIS_0h$mat_notlog[g,CIS_medulla.0h], CIS_12h$mat_notlog[g,CIS_medulla.12h], 
                           CIS_24h$mat_notlog[g,CIS_medulla.24h], CIS_48h$mat_notlog[g,CIS_medulla.48h]), group = 'I.Medulla')
head(df3)

df4 <- data.frame(time = c(rep(0, ncol(CIS_0h$mat_notlog)), rep(12, ncol(CIS_12h$mat_notlog)), 
                           rep(24, ncol(CIS_24h$mat_notlog)), rep(48, ncol(CIS_48h$mat_notlog))),
                  gexp = c(CIS_0h$mat_notlog[g,], CIS_12h$mat_notlog[g,], 
                           CIS_24h$mat_notlog[g,], CIS_48h$mat_notlog[g,]), group = 'Global')
head(df4)

df<-rbind(df1, df2, df3, df4)

df_comb <- aggregate(gexp ~ time + group, data = df,
                     FUN = function(x) c(mean = mean(x),
                                         max=max(x),
                                         min=min(x),
                                         sd = sd(x),
                                         se = sd(x)/sqrt(length(x))))

df_comb <- cbind(df_comb[,1:2], df_comb$gexp)
colnames(df_comb) <- c("time","group","mean","max", "min","sd","se")
df_comb$group<-factor(df_comb$group, levels=c("Global", "Cortex", "O.Medulla", "I.Medulla"))


###Ribbon Plots: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
glineSD<-ggplot(df_comb,aes(x=time, y=mean, color=group, fill=group, group=group)) +
  geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd),alpha=0.25, colour=NA) +
  geom_line(size=0.5) + geom_point(size=0.7) + 
  scale_color_manual(values=c("Global"="black","Cortex"="red", 
                              "O.Medulla"="blue", "I.Medulla"="green")) +
  scale_fill_manual(values=c("Global"="black","Cortex"="red", 
                             "O.Medulla"="blue","I.Medulla"="green")) +
  #ylim(overall_min0, overall_max0) + 
  theme_light() +
  theme(axis.text.x = element_text(color="black", size=11),  
        axis.text.y = element_text(color="black",size=11),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.text = element_text(color="black", size=12),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"), axis.ticks.length = unit(0.1, "cm")) +
  xlab("time (hours)") + ylab("Avg. gene exp. (CPM)")
glineSD

file<-paste0(dir,"/Codes2_Figures/Figure2/pdfs/Ttc36_CIS_AllTimepoints.pdf")  
pdf(file, height=3, width=15)
gridExtra::grid.arrange(comb_plot,glineSD, ncol=2, widths=c(16, 5), top=g)
dev.off()

#Slc34a1~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
g<-'Slc34a1' #Insert any gene of interest, for e.g., Igfbp7
overall_min <- min(CIS_0h$mat_notlog[g,], CIS_12h$mat_notlog[g,], 
                   CIS_24h$mat_notlog[g,], CIS_48h$mat_notlog[g,])

overall_max <- quantile(c(CIS_0h$mat_notlog[g,], CIS_12h$mat_notlog[g,], 
                          CIS_24h$mat_notlog[g,], CIS_48h$mat_notlog[g,]), 0.95)

#Set up a custom color scale (as per your choice)
color_scale <- scale_color_gradientn(colors=viridis(5),
                                     #colors=c("grey","magenta","yellow","red"),
                                     values = c(0, 0.8, 1),
                                     limits = c(overall_min, overall_max), oob = scales::squish)

#Set up a custom theme for plotting
my_theme <- theme_bw() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.text = element_text(),
        legend.position="none"
  )

#Spatiotemporal Heatmaps
df <- data.frame(CIS_0h$pos, gexp=CIS_0h$mat_notlog[g,])
g1 <- ggplot(df, aes(x=X, y=Y, col=gexp)) + geom_point(size=0.5, alpha=1) + color_scale + my_theme  + ggtitle('CIS 0h') 
df <- data.frame(CIS_12h$pos, gexp=CIS_12h$mat_notlog[g,])
g2 <- ggplot(df, aes(x=X, y=Y, col=gexp)) + geom_point(size=0.5, alpha=1)  + color_scale  + my_theme + ggtitle('CIS 12h')
df <- data.frame(CIS_24h$pos, gexp=CIS_24h$mat_notlog[g,])
g3 <- ggplot(df, aes(x=X, y=Y, col=gexp)) + geom_point(size=0.5, alpha=1) + color_scale + my_theme + ggtitle('CIS 24h')
df <- data.frame(CIS_48h$pos, gexp=CIS_48h$mat_notlog[g,])
g4 <- ggplot(df, aes(x=X, y=Y, col=gexp)) + geom_point(size=0.5, alpha=1) + color_scale  + theme(legend.text = element_text(size=10), legend.ticks = element_line(color = "black")) + ggtitle('CIS 48h')

# Extract the legend from one plot
legend <- get_legend(g4)
g4 <- ggplot(df, aes(x=X, y=Y, col=gexp)) + geom_point(size=0.5, alpha=1) + color_scale  + my_theme + ggtitle('CIS 48h')

comb_plot<-plot_grid(g1, g2, g3, g4, legend, ncol = 5, rel_widths = c(1, 1, 1, 1, 0.25))

#RibbonPlots
df1 <- data.frame(time = c(rep(0, length(CIS_cortex.0h)), rep(12, length(CIS_cortex.12h)), 
                           rep(24, length(CIS_cortex.24h)), rep(48, length(CIS_cortex.48h))),
                  gexp = c(CIS_0h$mat_notlog[g,CIS_cortex.0h], CIS_12h$mat_notlog[g,CIS_cortex.12h], 
                           CIS_24h$mat_notlog[g,CIS_cortex.24h], CIS_48h$mat_notlog[g,CIS_cortex.48h]),
                  
                  group = 'Cortex')
head(df1)

df2 <- data.frame(time = c(rep(0, length(CIS_interface.0h)), rep(12, length(CIS_interface.12h)), 
                           rep(24, length(CIS_interface.24h)), rep(48, length(CIS_interface.48h))),
                  gexp = c(CIS_0h$mat_notlog[g,CIS_interface.0h], CIS_12h$mat_notlog[g,CIS_interface.12h], 
                           CIS_24h$mat_notlog[g,CIS_interface.24h], CIS_48h$mat_notlog[g,CIS_interface.48h]), group = 'O.Medulla')
head(df2)

df3 <- data.frame(time = c(rep(0, length(CIS_medulla.0h)), rep(12, length(CIS_medulla.12h)), 
                           rep(24, length(CIS_medulla.24h)), rep(48, length(CIS_medulla.48h))),
                  gexp = c(CIS_0h$mat_notlog[g,CIS_medulla.0h], CIS_12h$mat_notlog[g,CIS_medulla.12h], 
                           CIS_24h$mat_notlog[g,CIS_medulla.24h], CIS_48h$mat_notlog[g,CIS_medulla.48h]), group = 'I.Medulla')
head(df3)

df4 <- data.frame(time = c(rep(0, ncol(CIS_0h$mat_notlog)), rep(12, ncol(CIS_12h$mat_notlog)), 
                           rep(24, ncol(CIS_24h$mat_notlog)), rep(48, ncol(CIS_48h$mat_notlog))),
                  gexp = c(CIS_0h$mat_notlog[g,], CIS_12h$mat_notlog[g,], 
                           CIS_24h$mat_notlog[g,], CIS_48h$mat_notlog[g,]), group = 'Global')
head(df4)

df<-rbind(df1, df2, df3, df4)

df_comb <- aggregate(gexp ~ time + group, data = df,
                     FUN = function(x) c(mean = mean(x),
                                         max=max(x),
                                         min=min(x),
                                         sd = sd(x),
                                         se = sd(x)/sqrt(length(x))))

df_comb <- cbind(df_comb[,1:2], df_comb$gexp)
colnames(df_comb) <- c("time","group","mean","max", "min","sd","se")
df_comb$group<-factor(df_comb$group, levels=c("Global", "Cortex", "O.Medulla", "I.Medulla"))


###Ribbon Plots: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
glineSD<-ggplot(df_comb,aes(x=time, y=mean, color=group, fill=group, group=group)) +
  geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd),alpha=0.25, colour=NA) +
  geom_line(size=0.5) + geom_point(size=0.7) + 
  scale_color_manual(values=c("Global"="black","Cortex"="red", 
                              "O.Medulla"="blue", "I.Medulla"="green")) +
  scale_fill_manual(values=c("Global"="black","Cortex"="red", 
                             "O.Medulla"="blue","I.Medulla"="green")) +
  #ylim(overall_min0, overall_max0) + 
  theme_light() +
  theme(axis.text.x = element_text(color="black", size=11),  
        axis.text.y = element_text(color="black",size=11),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.text = element_text(color="black", size=12),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"), axis.ticks.length = unit(0.1, "cm")) +
  xlab("time (hours)") + ylab("Avg. gene exp. (CPM)")
glineSD

file<-paste0(dir,"/Codes2_Figures/Figure2/pdfs/Slc34a1_CIS_AllTimepoints.pdf")  
pdf(file, height=3, width=15)
gridExtra::grid.arrange(comb_plot,glineSD, ncol=2, widths=c(16, 5), top=g)
dev.off()