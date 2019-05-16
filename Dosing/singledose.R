# load relevant libraries
library(R.matlab)
library(ggplot2)
library(grid)
library(gridExtra)
# install.packages('cowplot') may be necessary
library(cowplot)

# Change working directory to where file is (only works in R studio)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# import data
data_frommat <- readMat("onedose_out.mat",header=T) # read in as list
data_frommat2 <- readMat("onedose_out2.mat",header=T) 
data_frommat3 <- readMat("onedose_out3.mat",header=T)

df1 <- as.data.frame(data_frommat) # convert list to data frame
df2 <- as.data.frame(data_frommat2)
df3 <- as.data.frame(data_frommat3)

names(df1) <- c('tumorM','t','plasmaT','tumorT','bal')
names(df2) <- c('tumorM2','t2','plasmaT2','tumorT2','bal2')
names(df3) <- c('tumorM3','t3','plasmaT3','tumorT3','bal3')

g1 <- ggplot() + 
  geom_line(data=df1,aes(x=t,y=plasmaT,color="Regimen 1")) +
  geom_line(data=df2,aes(x=t2,y=plasmaT2,color="Regimen 2")) +
  geom_line(data=df3,aes(x=t3,y=plasmaT3,color="Regimen 3")) +
  scale_x_continuous("Time (hours)") +
  scale_y_continuous("TMZ in Plasma (M)") +
  scale_color_manual(name='',values=c('Regimen 1'='#CC0033','Regimen 2'='#00CCFF','Regimen 3'='#9933FF')) +
  theme(legend.position=c(0.8,0.8),legend.background=element_rect(fill='transparent'),text=element_text(size=14)) +
  theme_bw()
  
g2 <- ggplot() + 
  geom_line(data=df1,aes(x=t,y=tumorT,color="Regimen 1")) +
  geom_line(data=df2,aes(x=t2,y=tumorT2,color="Regimen 2")) +
  geom_line(data=df3,aes(x=t3,y=tumorT3,color="Regimen 3")) +
  scale_x_continuous("Time (hours)") +
  scale_y_continuous("TMZ in Tumor (M)") +
  scale_color_manual(name='',values=c('Regimen 1'='#CC0033','Regimen 2'='#00CCFF','Regimen 3'='#9933FF')) +
  theme(legend.position=c(0.8,0.8),legend.background=element_rect(fill='transparent'),text=element_text(size=14)) +
  theme_bw()

g3 <- ggplot() + 
  geom_line(data=df1,aes(x=t,y=tumorM,color="Regimen 1")) +
  geom_line(data=df2,aes(x=t2,y=tumorM2,color="Regimen 2")) +
  geom_line(data=df3,aes(x=t3,y=tumorM3,color="Regimen 3")) +
  scale_x_continuous("Time (hours)") +
  scale_y_continuous("MTIC in Tumor (M)") +
  scale_color_manual(name='',values=c('Regimen 1'='#CC0033','Regimen 2'='#00CCFF','Regimen 3'='#9933FF')) +
  theme_bw()

g4 <- ggplot() + 
  geom_line(data=df1,aes(x=t,y=bal,color="Regimen 1")) +
  geom_line(data=df2,aes(x=t2,y=bal2,color="Regimen 2")) +
  geom_line(data=df3,aes(x=t3,y=bal3,color="Regimen 3")) +
  scale_x_continuous("Time (hours)") +
  scale_y_continuous("Mass Balance") +
  scale_color_manual(name='',values=c('Regimen 1'='#CC0033','Regimen 2'='#00CCFF','Regimen 3'='#9933FF')) +
  theme_bw()

p1 <- plot_grid(g1 + theme(legend.position="none"),
                g2 + theme(legend.position="none"),
                nrow=1)
p2 <- plot_grid(g3 + theme(legend.position="none"),
                g4 + theme(legend.position="none"),
                nrow=1)
legend <- get_legend(g1 + theme(legend.position="bottom"))


title = textGrob("Concentration profiles for a single dose",gp=gpar(fontsize=15,font=3))
out <- grid.arrange(p1,p2,legend,nrow=3,heights=c(2.75,2.75,0.5))
ggsave(file="SingleDose.pdf",plot=out,width=10,height=6) # Save the figure 
