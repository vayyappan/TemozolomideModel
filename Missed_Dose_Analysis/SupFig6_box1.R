# we want boxplots
# we are reading in plasma TMZ AUC for the mgmt cohort

library(ggplot2)
library(R.matlab)
library(gridExtra)
library(reshape2)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# part i
Data_frommat3 <- readMat("missDose_TAA.mat", header=T)
df3 <- as.data.frame(Data_frommat3)
df3 <- log10(df3)

missed <- subset(df3, select=c('personIndex','tumorAdduct.AUC.total.1', 'tumorAdduct.AUC.total.2'))
meltedBox1 <- melt(missed, id.vars="personIndex")
meltedBox1$variable <- factor(meltedBox1$variable, levels = c('tumorAdduct.AUC.total.1', 'tumorAdduct.AUC.total.2'), ordered=TRUE)
labz <- c('MGMT Expressing ', 'MGMT Silenced')
box1 <- ggplot(meltedBox1, aes(x=meltedBox1$variable, y=meltedBox1$value)) + geom_boxplot(fill=c("indianred2", "#CCCCCC")) + theme_bw() + theme(text = element_text(size=40)) + ggtitle("Missed \n") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=labz) + theme(axis.text.x=element_text(angle=45, hjust=1)) + scale_y_continuous("log10(Tumor Adduct AUC)") + xlab(" ") 

missed2 <- subset(df3, select=c('personIndex','tumorAdduct.AUC.total.3', 'tumorAdduct.AUC.total.4'))
meltedBox2 <- melt(missed2, id.vars="personIndex")
meltedBox2$variable <- factor(meltedBox2$variable, levels = c('tumorAdduct.AUC.total.3', 'tumorAdduct.AUC.total.4'), ordered=TRUE)
labz <- c('MGMT Expressing ', 'MGMT Silenced')
box2 <- ggplot(meltedBox2, aes(x=meltedBox2$variable, y=meltedBox2$value)) + geom_boxplot(fill=c("indianred2", "#CCCCCC")) + theme_bw() + theme(text = element_text(size=40)) + ggtitle("Delay \n 5.6d") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=labz) + theme(axis.text.x=element_text(angle=45, hjust=1)) + scale_y_continuous("log10(Tumor Adduct AUC)") + xlab(" ") 

missed3 <- subset(df3, select=c('personIndex','tumorAdduct.AUC.total.5', 'tumorAdduct.AUC.total.6'))
meltedBox3 <- melt(missed3, id.vars="personIndex")
meltedBox3$variable <- factor(meltedBox3$variable, levels = c('tumorAdduct.AUC.total.5', 'tumorAdduct.AUC.total.6'), ordered=TRUE)
labz <- c('MGMT Expressing ', 'MGMT Silenced')
box3 <- ggplot(meltedBox3, aes(x=meltedBox3$variable, y=meltedBox3$value)) + geom_boxplot(fill=c("indianred2", "#CCCCCC")) + theme_bw() + theme(text = element_text(size=40)) + ggtitle("Delay \n 11.2d") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=labz) + theme(axis.text.x=element_text(angle=45, hjust=1)) + scale_y_continuous("log10(Tumor Adduct AUC)") + xlab(" ") 

missed4 <- subset(df3, select=c('personIndex','tumorAdduct.AUC.total.7', 'tumorAdduct.AUC.total.8'))
meltedBox4 <- melt(missed4, id.vars="personIndex")
meltedBox4$variable <- factor(meltedBox4$variable, levels = c('tumorAdduct.AUC.total.7', 'tumorAdduct.AUC.total.8'), ordered=TRUE)
labz <- c('MGMT Expressing ', 'MGMT Silenced')
box4 <- ggplot(meltedBox4, aes(x=meltedBox4$variable, y=meltedBox4$value)) + geom_boxplot(fill=c("indianred2", "#CCCCCC")) + theme_bw() + theme(text = element_text(size=40)) + ggtitle("Delay \n 16.8d") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=labz) + theme(axis.text.x=element_text(angle=45, hjust=1)) + scale_y_continuous("log10(Tumor Adduct AUC)") + xlab(" ") 

missed5 <- subset(df3, select=c('personIndex','tumorAdduct.AUC.total.9', 'tumorAdduct.AUC.total.10'))
meltedBox5 <- melt(missed5, id.vars="personIndex")
meltedBox5$variable <- factor(meltedBox5$variable, levels = c('tumorAdduct.AUC.total.9', 'tumorAdduct.AUC.total.10'), ordered=TRUE)
labz <- c('MGMT Expressing ', 'MGMT Silenced')
box5 <- ggplot(meltedBox5, aes(x=meltedBox5$variable, y=meltedBox5$value)) + geom_boxplot(fill=c("indianred2", "#CCCCCC")) + theme_bw() + theme(text = element_text(size=40)) + ggtitle("Delay \n 22.4d") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=labz) + theme(axis.text.x=element_text(angle=45, hjust=1)) + scale_y_continuous("log10(Tumor Adduct AUC)") + xlab(" ") 

missed6 <- subset(df3, select=c('personIndex','tumorAdduct.AUC.total.11', 'tumorAdduct.AUC.total.12'))
meltedBox6 <- melt(missed6, id.vars="personIndex")
meltedBox6$variable <- factor(meltedBox6$variable, levels = c('tumorAdduct.AUC.total.11', 'tumorAdduct.AUC.total.12'), ordered=TRUE)
labz <- c('MGMT Expressing ', 'MGMT Silenced')
box6 <- ggplot(meltedBox6, aes(x=meltedBox6$variable, y=meltedBox6$value)) + geom_boxplot(fill=c("indianred2", "#CCCCCC")) + theme_bw() + theme(text = element_text(size=40)) + ggtitle("Double Dose \n After Skip") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=labz) + theme(axis.text.x=element_text(angle=45, hjust=1)) + scale_y_continuous("log10(Tumor Adduct AUC)") + xlab(" ") 

missed7 <- subset(df3, select=c('personIndex','tumorAdduct.AUC.total.13', 'tumorAdduct.AUC.total.14'))
meltedBox7 <- melt(missed7, id.vars="personIndex")
meltedBox7$variable <- factor(meltedBox7$variable, levels = c('tumorAdduct.AUC.total.13', 'tumorAdduct.AUC.total.14'), ordered=TRUE)
labz <- c('MGMT Expressing ', 'MGMT Silenced')
box7 <- ggplot(meltedBox7, aes(x=meltedBox7$variable, y=meltedBox7$value)) + geom_boxplot(fill=c("indianred2", "#CCCCCC")) + theme_bw() + theme(text = element_text(size=40)) + ggtitle("Control \n") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=labz) + theme(axis.text.x=element_text(angle=45, hjust=1)) + scale_y_continuous("log10(Tumor Adduct AUC)") + xlab(" ") 

out <- grid.arrange(box7, box1, box2, box3, box4, box5, box6, ncol=7)
ggsave(file="SupFig6_DoseAnalysis1.pdf",plot
       =out,width=49,height=16) 
