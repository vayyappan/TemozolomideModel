# load relevant libraries
library(ggplot2) 
library(R.matlab)
library(grid)
library(gridExtra)
library(reshape)
library(cowplot)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# part i
Data_frommat1 <- readMat("tumVolM.mat", header=T)
df1 <- as.data.frame(Data_frommat1)
names(df1) <- c('T', 'Missed', 'Delay-05.6d', 'Delay-11.2d', 'Delay-16.8d', 'Delay-22.4d', 'Double Dose After Skip', 'Control')
Data_frommat2 <- readMat("tumVolS.mat", header=T)
df2 <- as.data.frame(Data_frommat2)
names(df2) <- c('T', 'Missed', 'Delay-05.6d', 'Delay-11.2d', 'Delay-16.8d', 'Delay-22.4d', 'Double Dose After Skip', 'Control')
Data_frommat3 <- readMat("plasTM.mat", header=T) #actually tumor TMZ, just misnamed
df3 <- as.data.frame(Data_frommat3)
names(df3) <- c('T', 'Missed', 'Delay-05.6d', 'Delay-11.2d', 'Delay-16.8d', 'Delay-22.4d', 'Double Dose After Skip', 'Control')
Data_frommat4 <- readMat("plasTS.mat", header=T) #actually tumor TMZ, just misnamed
df4 <- as.data.frame(Data_frommat4)
names(df4) <- c('T', 'Missed', 'Delay-05.6d', 'Delay-11.2d', 'Delay-16.8d', 'Delay-22.4d', 'Double Dose After Skip', 'Control')
Data_frommat5 <- readMat("tumDM.mat", header=T)
df5 <- as.data.frame(Data_frommat5)*1e9
names(df5) <- c('T', 'Missed', 'Delay-05.6d', 'Delay-11.2d', 'Delay-16.8d', 'Delay-22.4d', 'Double Dose After Skip', 'Control')
Data_frommat6 <- readMat("tumDS.mat", header=T)
df6 <- as.data.frame(Data_frommat6)*1e9
names(df6) <- c('T', 'Missed', 'Delay-05.6d', 'Delay-11.2d', 'Delay-16.8d', 'Delay-22.4d', 'Double Dose After Skip', 'Control')

g <- ggplot(df1) +
  geom_line(aes(x = df1$T, y=df1$`Delay-05.6d`,color='Delay-05.6d'),size=1) +
  geom_line(aes(x = df1$T, y=df1$`Delay-11.2d`,color='Delay-11.2d'),size=1) +
  geom_line(aes(x = df1$T, y=df1$`Delay-16.8d`,color='Delay-16.8d'),size=1) +
  geom_line(aes(x = df1$T, y=df1$`Delay-22.4d`,color='Delay-22.4d'),size=1) +
  geom_line(aes(x = df1$T, y=df1$`Double Dose After Skip`,color='Double Dose After Miss'),size=1) +
  geom_line(aes(x = df1$T, y=df1$`Missed`,color='Missed Dose'),size=1) + 
  geom_line(aes(x = df1$T, y=df1$`Control`,color='Control'),size=1) +
  scale_x_continuous("Time (months)") +
  scale_y_continuous("MGMT Expressed)") +
  scale_color_manual(name='',values=c('Control'='#CC0033','Delay-05.6d'='#FF9933','Delay-11.2d'='#EFC411','Delay-16.8d'='#00CC00','Delay-22.4d'='#00CCFF','Double Dose After Miss'='#334CFF','Missed Dose'='#9933FF')) +
  theme_bw(base_size=18) +
  theme(plot.title=element_text(hjust=0.5)) +
  ggtitle("","Tumor Volume (mL)")

g1 <- ggplot(df2) +
  geom_line(aes(x = df2$T, y=df2$`Delay-05.6d`,color='Delay-05.6d'),size=1) +
  geom_line(aes(x = df2$T, y=df2$`Delay-11.2d`,color='Delay-11.2d'),size=1) +
  geom_line(aes(x = df2$T, y=df2$`Delay-16.8d`,color='Delay-16.8d'),size=1) +
  geom_line(aes(x = df2$T, y=df2$`Delay-22.4d`,color='Delay-22.4d'),size=1) +
  geom_line(aes(x = df2$T, y=df2$`Double Dose After Skip`,color='Double Dose After Miss'),size=1) +
  geom_line(aes(x = df2$T, y=df2$`Missed`,color='Missed Dose'),size=1) + 
  geom_line(aes(x = df2$T, y=df2$`Control`,color='Control'),size=1) +
  scale_x_continuous("Time (months)") +
  scale_y_continuous("MGMT Silenced") +
  scale_color_manual(name='',values=c('Control'='#CC0033','Delay-05.6d'='#FF9933','Delay-11.2d'='#EFC411','Delay-16.8d'='#00CC00','Delay-22.4d'='#00CCFF','Double Dose After Miss'='#334CFF','Missed Dose'='#9933FF')) +
  theme_bw(base_size=18) +
  ggtitle("")

g2 <- ggplot(df3) +
  geom_line(aes(x = df3$T, y=df3$`Delay-05.6d`,color='Delay-05.6d'),size=1) +
  geom_line(aes(x = df3$T, y=df3$`Delay-11.2d`,color='Delay-11.2d'),size=1) +
  geom_line(aes(x = df3$T, y=df3$`Delay-16.8d`,color='Delay-16.8d'),size=1) +
  geom_line(aes(x = df3$T, y=df3$`Delay-22.4d`,color='Delay-22.4d'),size=1) +
  geom_line(aes(x = df3$T, y=df3$`Double Dose After Skip`,color='Double Dose After Miss'),size=1) +
  geom_line(aes(x = df3$T, y=df3$`Missed`,color='Missed Dose'),size=1) + 
  geom_line(aes(x = df3$T, y=df3$`Control`,color='Control'),size=1) +
  scale_x_continuous("Time (months)") +
  scale_y_continuous("") +
  scale_color_manual(name='',values=c('Control'='#CC0033','Delay-05.6d'='#FF9933','Delay-11.2d'='#EFC411','Delay-16.8d'='#00CC00','Delay-22.4d'='#00CCFF','Double Dose After Miss'='#334CFF','Missed Dose'='#9933FF')) +
  theme_bw(base_size=18) +
  theme(plot.title=element_text(hjust=0.5)) +
  ggtitle("","[Tumor TMZ] (uM)")

g3 <- ggplot(df4) +
  geom_line(aes(x = df4$T, y=df4$`Delay-05.6d`,color='Delay-05.6d'),size=1) +
  geom_line(aes(x = df4$T, y=df4$`Delay-11.2d`,color='Delay-11.2d'),size=1) +
  geom_line(aes(x = df4$T, y=df4$`Delay-16.8d`,color='Delay-16.8d'),size=1) +
  geom_line(aes(x = df4$T, y=df4$`Delay-22.4d`,color='Delay-22.4d'),size=1) +
  geom_line(aes(x = df4$T, y=df4$`Double Dose After Skip`,color='Double Dose After Miss'),size=1) +
  geom_line(aes(x = df4$T, y=df4$`Missed`,color='Missed Dose'),size=1) + 
  geom_line(aes(x = df4$T, y=df4$`Control`,color='Control'),size=1) +
  scale_x_continuous("Time (months)") +
  scale_y_continuous("") +
  scale_color_manual(name='',values=c('Control'='#CC0033','Delay-05.6d'='#FF9933','Delay-11.2d'='#EFC411','Delay-16.8d'='#00CC00','Delay-22.4d'='#00CCFF','Double Dose After Miss'='#334CFF','Missed Dose'='#9933FF')) +
  theme_bw(base_size=18) +
  ggtitle("")

g4 <- ggplot(df5) +
  geom_line(aes(x = df5$T, y=df5$`Delay-05.6d`,color='Delay-05.6d'),size=1) +
  geom_line(aes(x = df5$T, y=df5$`Delay-11.2d`,color='Delay-11.2d'),size=1) +
  geom_line(aes(x = df5$T, y=df5$`Delay-16.8d`,color='Delay-16.8d'),size=1) +
  geom_line(aes(x = df5$T, y=df5$`Delay-22.4d`,color='Delay-22.4d'),size=1) +
  geom_line(aes(x = df5$T, y=df5$`Double Dose After Skip`,color='Double Dose After Miss'),size=1) +
  geom_line(aes(x = df5$T, y=df5$`Missed`,color='Missed Dose'),size=1) + 
  geom_line(aes(x = df5$T, y=df5$`Control`,color='Control'),size=1) +
  scale_x_continuous("Time (months)") +
  scale_y_continuous("") +
  scale_color_manual(name='',values=c('Control'='#CC0033','Delay-05.6d'='#FF9933','Delay-11.2d'='#EFC411','Delay-16.8d'='#00CC00','Delay-22.4d'='#00CCFF','Double Dose After Miss'='#334CFF','Missed Dose'='#9933FF')) +
  theme_bw(base_size=18) +
  theme(plot.title=element_text(hjust=0.5)) +
  ggtitle("","[Tumor DNA Adduct] (nM)")

g5 <- ggplot(df6) +
  geom_line(aes(x = df6$T, y=df6$`Delay-05.6d`,color='Delay-05.6d'),size=1) +
  geom_line(aes(x = df6$T, y=df6$`Delay-11.2d`,color='Delay-11.2d'),size=1) +
  geom_line(aes(x = df6$T, y=df6$`Delay-16.8d`,color='Delay-16.8d'),size=1) +
  geom_line(aes(x = df6$T, y=df6$`Delay-22.4d`,color='Delay-22.4d'),size=1) +
  geom_line(aes(x = df6$T, y=df6$`Double Dose After Skip`,color='Double Dose After Miss'),size=1) +
  geom_line(aes(x = df6$T, y=df6$`Missed`,color='Missed Dose'),size=1) + 
  geom_line(aes(x = df6$T, y=df6$`Control`,color='Control'),size=1) +
  scale_x_continuous("Time (monts)") +
  scale_y_continuous("") +
  scale_color_manual(name='',values=c('Control'='#CC0033','Delay-05.6d'='#FF9933','Delay-11.2d'='#EFC411','Delay-16.8d'='#00CC00','Delay-22.4d'='#00CCFF','Double Dose After Miss'='#334CFF','Missed Dose'='#9933FF')) +
  theme_bw(base_size=18) +
  ggtitle("")

p1 <- plot_grid(g + theme(legend.position="none"),
                g2 + theme(legend.position="none"),
                g4 + theme(legend.position="none"),
                nrow=1)
p2 <- plot_grid(g1 + theme(legend.position="none"),
                g3 + theme(legend.position="none"),
                g5 + theme(legend.position="none"),
                nrow=1)

figs <- grid.arrange(p1,p2,nrow=2)
legend <- get_legend(g + theme(legend.position="bottom"))

out <- grid.arrange(figs,legend,nrow=2,heights=c(7,1))
ggsave(file="Fig6_DoseCurves.pdf",plot
       =out,width=16,height=8) 
