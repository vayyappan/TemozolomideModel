# we want boxplots
# we are reading in plasma TMZ AUC for the mgmt cohort

library(ggplot2)
library(R.matlab)
library(gridExtra)
library(reshape2)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# part iii
Data_frommat5 <- readMat("missDose_F.mat", header=T)
df5 <- as.data.frame(Data_frommat5)
df5 <- abs(df5)

missed <- subset(df5, select=c('personIndex','finalVol.total.1', 'finalVol.total.2'))
meltedBox1 <- melt(missed, id.vars="personIndex")
meltedBox1$variable <- factor(meltedBox1$variable, levels = c('finalVol.total.1', 'finalVol.total.2'), ordered=TRUE)
labz <- c('MGMT Expressing ', 'MGMT Silenced')
box13 <- ggplot(meltedBox1, aes(x=meltedBox1$variable, y=meltedBox1$value)) + geom_boxplot(fill=c("indianred2", "#CCCCCC")) + theme_bw() + theme(text = element_text(size=40)) + ggtitle("Missed \n") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=labz) + theme(axis.text.x=element_text(angle=45, hjust=1)) + scale_y_continuous("Final Tumor Volume (mL)") + xlab(" ") 

missed2 <- subset(df5, select=c('personIndex','finalVol.total.3', 'finalVol.total.4'))
meltedBox2 <- melt(missed2, id.vars="personIndex")
meltedBox2$variable <- factor(meltedBox2$variable, levels = c('finalVol.total.3', 'finalVol.total.4'), ordered=TRUE)
labz <- c('MGMT Expressing ', 'MGMT Silenced')
box14 <- ggplot(meltedBox2, aes(x=meltedBox2$variable, y=meltedBox2$value)) + geom_boxplot(fill=c("indianred2", "#CCCCCC")) + theme_bw() + theme(text = element_text(size=40)) + ggtitle("Delay \n 5.6d") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=labz) + theme(axis.text.x=element_text(angle=45, hjust=1)) + scale_y_continuous("Final Tumor Volume (mL)") + xlab(" ") 

missed3 <- subset(df5, select=c('personIndex','finalVol.total.6', 'finalVol.total.5'))
meltedBox3 <- melt(missed3, id.vars="personIndex")
meltedBox3$variable <- factor(meltedBox3$variable, levels = c('finalVol.AUC.total.6', 'finalVol.total.5'), ordered=TRUE)
labz <- c('MGMT Expressing ', 'MGMT Silenced')
box15 <- ggplot(meltedBox3, aes(x=meltedBox3$variable, y=meltedBox3$value)) + geom_boxplot(fill=c("indianred2", "#CCCCCC")) + theme_bw() + theme(text = element_text(size=40)) + ggtitle("Delay \n 11.2d") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=labz) + theme(axis.text.x=element_text(angle=45, hjust=1)) + scale_y_continuous("Final Tumor Volume (mL)") + xlab(" ") 

missed4 <- subset(df5, select=c('personIndex','finalVol.total.7', 'finalVol.total.8'))
meltedBox4 <- melt(missed4, id.vars="personIndex")
meltedBox4$variable <- factor(meltedBox4$variable, levels = c('finalVol.total.7', 'finalVol.total.8'), ordered=TRUE)
labz <- c('MGMT Expressing ', 'MGMT Silenced')
box16 <- ggplot(meltedBox4, aes(x=meltedBox4$variable, y=meltedBox4$value)) + geom_boxplot(fill=c("indianred2", "#CCCCCC")) + theme_bw() + theme(text = element_text(size=40)) + ggtitle("Delay \n 16.8d") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=labz) + theme(axis.text.x=element_text(angle=45, hjust=1)) + scale_y_continuous("Final Tumor Volume (mL)") + xlab(" ") 

missed5 <- subset(df5, select=c('personIndex','finalVol.total.9', 'finalVol.total.10'))
meltedBox5 <- melt(missed5, id.vars="personIndex")
meltedBox5$variable <- factor(meltedBox5$variable, levels = c('finalVol.total.9', 'finalVol.total.10'), ordered=TRUE)
labz <- c('MGMT Expressing ', 'MGMT Silenced')
box17 <- ggplot(meltedBox5, aes(x=meltedBox5$variable, y=meltedBox5$value)) + geom_boxplot(fill=c("indianred2", "#CCCCCC")) + theme_bw() + theme(text = element_text(size=40)) + ggtitle("Delay \n 22.4d") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=labz) + theme(axis.text.x=element_text(angle=45, hjust=1)) + scale_y_continuous("Final Tumor Volume (mL)") + xlab(" ") 

missed6 <- subset(df5, select=c('personIndex','finalVol.total.11', 'finalVol.total.12'))
meltedBox6 <- melt(missed6, id.vars="personIndex")
meltedBox6$variable <- factor(meltedBox6$variable, levels = c('finalVol.total.11', 'finalVol.total.12'), ordered=TRUE)
labz <- c('MGMT Expressing ', 'MGMT Silenced')
box18 <- ggplot(meltedBox6, aes(x=meltedBox6$variable, y=meltedBox6$value)) + geom_boxplot(fill=c("indianred2", "#CCCCCC")) + theme_bw() + theme(text = element_text(size=40)) + ggtitle("Double Dose \n After Skip") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=labz) + theme(axis.text.x=element_text(angle=45, hjust=1)) + scale_y_continuous("Final Tumor Volume (mL)") + xlab(" ") 

missed7 <- subset(df5, select=c('personIndex','finalVol.total.13', 'finalVol.total.14'))
meltedBox7 <- melt(missed7, id.vars="personIndex")
meltedBox7$variable <- factor(meltedBox7$variable, levels = c('finalVol.total.13', 'finalVol.total.14'), ordered=TRUE)
labz <- c('MGMT Expressing ', 'MGMT Silenced')
box19 <- ggplot(meltedBox7, aes(x=meltedBox7$variable, y=meltedBox7$value)) + geom_boxplot(fill=c("indianred2", "#CCCCCC")) + theme_bw() + theme(text = element_text(size=40)) + ggtitle("Control \n") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=labz) + theme(axis.text.x=element_text(angle=45, hjust=1)) + scale_y_continuous("Final Tumor Volume (mL)") + xlab(" ") 

out <- grid.arrange(box19, box13, box14, box15, box16, box17, box18, ncol=7)
ggsave(file="SupFig6_DoseAnalysis3.pdf",plot
       =out,width=49,height=16) 