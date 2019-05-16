# load relevant libraries
library(R.matlab)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)

# Change working directory to where file is (only works in R studio)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# import data
data_frommat <- readMat("multiplecycles_out.mat",header=T) # read in as list
df_frommat <- as.data.frame(data_frommat) # convert list to data frame
names(df_frommat) <- c('t','bal','plasmaT','tumorM','tumorT')

g1 <- ggplot(df_frommat) + 
  geom_line(aes(x=t,y=plasmaT,color="Regimen 1")) +
  scale_x_continuous("Time (Months)",limits=c(0,6)) +
  scale_y_continuous("TMZ in Plasma (M)") +
  scale_color_manual(name='',values=c('Regimen 1'='#CC0033')) +
  theme_bw()

g2 <- ggplot(df_frommat) + 
  geom_line(aes(x=t,y=tumorT,color="Regimen 1")) +
  scale_x_continuous("Time (Months)",limits=c(0,6)) +
  scale_y_continuous("TMZ in Tumor (M)") +
  scale_color_manual(name='',values=c('Regimen 1'='#CC0033')) +
  theme_bw()

g3 <- ggplot(df_frommat) + 
  geom_line(aes(x=t,y=tumorM,color="Regimen 1")) +
  scale_x_continuous("Time (Months)",limits=c(0,6)) +
  scale_y_continuous("MTIC in Tumor (M)") +
  scale_color_manual(name='',values=c('Regimen 1'='#CC0033')) +
  theme_bw()

g4 <- ggplot(df_frommat) + 
  geom_line(aes(x=t,y=bal,color="Regimen 1")) +
  scale_x_continuous("Time (Months)",limits=c(0,6)) +
  scale_y_continuous("Mass Balance") +
  scale_color_manual(name='',values=c('Regimen 1'='#CC0033')) +
  theme_bw()

p1 <- plot_grid(g1 + theme(legend.position="none"),
                g2 + theme(legend.position="none"),
                nrow=1)
p2 <- plot_grid(g3 + theme(legend.position="none"),
                g4 + theme(legend.position="none"),
                nrow=1)

title = textGrob("Concentration Profiles for Multiple Cycles",gp=gpar(fontsize=15,font=3))
out = grid.arrange(p1,p2,nrow=2,top=title)
ggsave(file="MultipleCycles.pdf",plot=out,width=10,height=6) # Save the figure 