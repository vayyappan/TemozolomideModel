# load relevant libraries
library(R.matlab)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(stringr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



# import data
data_frommat <- readMat("singlecycle_reg1.mat",header=T) # read in as list
data_frommat2 <- readMat("singlecycle_reg2.mat",header=T)
data_frommat3 <- readMat("singlecycle_reg3.mat",header=T)

df1 <- data.frame(data_frommat) # convert list to data frame
df2 <- data.frame(data_frommat2)
df3 <- data.frame(data_frommat3)

# Log scale
x = c(1,3,4,5)
df1[x] <- 1e6*(df1)[x]
df2[x] <- 1e6*(df2)[x]
df3[x] <- 1e6*(df3)[x]

names(df1) <- c('tumorM','t','plasmaT','tumorT','bal')
names(df2) <- c('tumorM2','t2','plasmaT2','tumorT2','bal2')
names(df3) <- c('tumorM3','t3','plasmaT3','tumorT3','bal3')

g1 <- ggplot() + 
  geom_line(data=df1,aes(x=t,y=plasmaT,color="Regimen 1")) +
  geom_line(data=df2,aes(x=t2,y=plasmaT2,color="Regimen 2")) +
  geom_line(data=df3,aes(x=t3,y=plasmaT3,color="Regimen 3")) +
  scale_x_continuous("Time (days)",limits=c(0,28)) +
  scale_y_continuous("[Plasma TMZ] (uM)") +
  scale_color_manual(name='',values=c('Regimen 1'='#CC0033','Regimen 2'='#00CCFF','Regimen 3'='#9933FF')) +
  theme_bw(base_size=18) +
  ggtitle("")

g2 <- ggplot() +
  geom_line(data=df1,aes(x=t,y=tumorT,color="Regimen 1")) +
  geom_line(data=df2,aes(x=t2,y=tumorT2,color="Regimen 2")) +
  geom_line(data=df3,aes(x=t3,y=tumorT3,color="Regimen 3")) +
  scale_x_continuous("Time (days)",limits=c(0,28)) +
  scale_y_continuous("[Tumor TMZ] (uM)") +
  scale_color_manual(name='',values=c('Regimen 1'='#CC0033','Regimen 2'='#00CCFF','Regimen 3'='#9933FF')) +
  theme_bw(base_size=18) +
  ggtitle("")

g3 <- ggplot() +
  geom_line(data=df1,aes(x=t,y=tumorM,color="Regimen 1")) +
  geom_line(data=df2,aes(x=t2,y=tumorM2,color="Regimen 2")) +
  geom_line(data=df3,aes(x=t3,y=tumorM3,color="Regimen 3")) +
  scale_x_continuous("Time (days)",limits=c(0,28)) +
  scale_y_continuous("[Tumor MTIC] (uM)") +
  scale_color_manual(name='',values=c('Regimen 1'='#CC0033','Regimen 2'='#00CCFF','Regimen 3'='#9933FF')) +
  theme_bw(base_size=18) +
  ggtitle("")

g4 <- ggplot() +
  geom_line(data=df1,aes(x=t,y=bal,color="Regimen 1")) +
  geom_line(data=df2,aes(x=t2,y=bal2,color="Regimen 2")) +
  geom_line(data=df3,aes(x=t3,y=bal3,color="Regimen 3")) +
  scale_x_continuous("Time (days)",limits=c(0,28)) +
  scale_y_continuous("Mass Balance (umol)") +
  scale_color_manual(name='',values=c('Regimen 1'='#CC0033','Regimen 2'='#00CCFF','Regimen 3'='#9933FF')) +
  theme_bw(base_size=18) +
  ggtitle("")

p1 <- plot_grid(g1 + theme(legend.position="none"),
               g2 + theme(legend.position="none"),
               nrow=1)
p2 <- plot_grid(g3 + theme(legend.position="none"),
                g4 + theme(legend.position="none"),
                nrow=1)
figs <- grid.arrange(p1,p2,nrow=2)
legend <- get_legend(g1 + theme(legend.position="bottom"))

#title = textGrob("Single cycle concentration profiles for different dosing regimens",gp=gpar(fontsize=15,font=3))
out <- grid.arrange(figs,legend,nrow=2,heights=c(6,1))
ggsave(file="Fig2_DiffRegimens.pdf",plot=out,width=10,height=7) # Save the figure 