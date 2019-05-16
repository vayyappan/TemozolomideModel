# load relevant libraries
library(ggplot2) 
library(R.matlab)
library(grid)
library(gridExtra)
library(reshape)
library(cowplot)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# import data
data_frommat <- readMat("tumorPD_5on23off.mat",header=T) # read in as list
data_frommat2 <- readMat("tumorPD_7on8off.mat",header=T) 
data_frommat3 <- readMat("tumorPD_28on.mat",header=T) 
data_frommatM0 <- readMat("tumorPDM0_5on23off.mat",header=T)
data_frommatM02 <- readMat("tumorPDM0_7on8off.mat",header=T) 
data_frommatM03 <- readMat("tumorPDM0_28on.mat",header=T) 

df1 <- data.frame(data_frommat) # convert list to data frame
df2 <- data.frame(data_frommat2) 
df3 <- data.frame(data_frommat3) 
dfM0 <- data.frame(data_frommatM0)
dfM02 <- data.frame(data_frommatM02)
dfM03 <- data.frame(data_frommatM03)

df1[2] <- 1e9*df1[2]
df2[2] <- 1e9*df2[2]
df3[2] <- 1e9*df3[2]
dfM0[2] <- 1e9*dfM0[2]
dfM02[2] <- 1e9*dfM02[2]
dfM03[2] <- 1e9*dfM03[2]

names(df1) <- c('t','adductT','resT','susT','totalT')
names(df2) <- c('t2','adductT2','resT2','susT2','totalT2')
names(df3) <- c('t3','adductT3','resT3','susT3','totalT3')
names(dfM0) <- c('tM0','adductTM0','resTM0','susTM0','totalTM0')
names(dfM02) <- c('tM02','adductTM02','resTM02','susTM02','totalTM02')
names(dfM03) <- c('tM03','adductTM03','resTM03','susTM03','totalTM03')

g1 <- ggplot() + 
  geom_line(data=df1,aes(x=t,y=adductT,color="Regimen 1")) +
  geom_line(data=df2,aes(x=t2,y=adductT2,color="Regimen 2")) +
  geom_line(data=df3,aes(x=t3,y=adductT3,color="Regimen 3")) +
  scale_x_continuous("Time (Months)",limits=c(0,6)) +
  scale_y_continuous("MGMT Expressed") +
  scale_color_manual(name='',values=c('Regimen 1'='#CC0033','Regimen 2'='#00CCFF','Regimen 3'='#9933FF')) +
  theme_bw(base_size = 18) +
  theme(plot.title=element_text(hjust=0.5)) +
  ggtitle("[Tumor DNA Adduct] (nM)")

g2 <- ggplot() +
  geom_line(data=df1,aes(x=t,y=totalT,color="Regimen 1")) +
  geom_line(data=df2,aes(x=t2,y=totalT2,color="Regimen 2")) +
  geom_line(data=df3,aes(x=t3,y=totalT3,color="Regimen 3")) +
  scale_x_continuous("Time (Months)",limits=c(0,6)) +
  scale_y_continuous("") +
  scale_color_manual(name='',values=c('Regimen 1'='#CC0033','Regimen 2'='#00CCFF','Regimen 3'='#9933FF')) +
  theme_bw(base_size = 18) +
  theme(plot.title=element_text(hjust=0.5)) +
  ggtitle("Tumor Growth (mL)")

g3 <- ggplot() + 
  geom_line(data=dfM0,aes(x=tM0,y=adductTM0,color="Regimen 1")) +
  geom_line(data=dfM02,aes(x=tM02,y=adductTM02,color="Regimen 2")) +
  geom_line(data=dfM03,aes(x=tM03,y=adductTM03,color="Regimen 3")) +
  scale_x_continuous("Time (Months)",limits=c(0,6)) +
  scale_y_continuous("MGMT Silenced") +
  scale_color_manual(name='',values=c('Regimen 1'='#CC0033','Regimen 2'='#00CCFF','Regimen 3'='#9933FF')) +
  theme_bw(base_size = 18) +
  ggtitle("")

g4 <- ggplot() +
  geom_line(data=dfM0,aes(x=tM0,y=totalTM0,color="Regimen 1")) +
  geom_line(data=dfM02,aes(x=tM02,y=totalTM02,color="Regimen 2")) +
  geom_line(data=dfM03,aes(x=tM03,y=totalTM03,color="Regimen 3")) +
  scale_x_continuous("Time (Months)",limits=c(0,6)) +
  scale_y_continuous("") +
  scale_color_manual(name='',values=c('Regimen 1'='#CC0033','Regimen 2'='#00CCFF','Regimen 3'='#9933FF')) +
  theme_bw(base_size = 18) +
  ggtitle("")

p1 <- plot_grid(g1 + theme(legend.position="none"),
                g2 + theme(legend.position="none"),
                nrow=1)
p2 <- plot_grid(g3 + theme(legend.position="none"),
                g4 + theme(legend.position="none"),
                nrow=1)
figs <- grid.arrange(p1,p2,nrow=2)
legend <- get_legend(g1 + theme(legend.position="bottom"))

#title = textGrob("Tumor pharmacodynamics and treatment efficacy",gp=gpar(fontsize=15,font=3))
out <- grid.arrange(figs,legend,nrow=2,heights=c(6,1))
ggsave(file="Fig5_TumorDiffRegimens.pdf",plot=out,width=10,height=7) # Save the figure