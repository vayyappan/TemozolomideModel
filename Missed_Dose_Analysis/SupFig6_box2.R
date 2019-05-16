# we want boxplots
# we are reading in plasma TMZ AUC for the mgmt cohort

library(ggplot2)
library(R.matlab)
library(gridExtra)
library(reshape2)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# part ii
Data_frommat4 <- readMat("missDose_CA.mat", header=T)
df4 <- as.data.frame(Data_frommat4)
df4 <- abs(df4)
df4 <- log10(df4)

missed <- subset(df4, select=c('personIndex','CtroughAdduct.total.1', 'CtroughAdduct.total.2'))
meltedBox1 <- melt(missed, id.vars="personIndex")
meltedBox1$variable <- factor(meltedBox1$variable, levels = c('CtroughAdduct.total.1', 'CtroughAdduct.total.2'), ordered=TRUE)
labz0 <- c('MGMT Expressing ', 'MGMT Silenced')
box7 <- ggplot(meltedBox1, aes(x=meltedBox1$variable, y=meltedBox1$value)) + geom_boxplot(fill=c("indianred2", "#CCCCCC")) + theme_bw() + theme(text = element_text(size=40)) + ggtitle("Missed \n") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=labz) + theme(axis.text.x=element_text(angle=45, hjust=1)) + scale_y_continuous("log10(Ctrough Adduct)") + xlab(" ") 

missed2 <- subset(df4, select=c('personIndex','CtroughAdduct.total.3', 'CtroughAdduct.total.4'))
meltedBox2 <- melt(missed2, id.vars="personIndex")
meltedBox2$variable <- factor(meltedBox2$variable, levels = c('CtroughAdduct.total.3', 'CtroughAdduct.total.4'), ordered=TRUE)
labz <- c('MGMT Expressing ', 'MGMT Silenced')
box8 <- ggplot(meltedBox2, aes(x=meltedBox2$variable, y=meltedBox2$value)) + geom_boxplot(fill=c("indianred2", "#CCCCCC")) + theme_bw() + theme(text = element_text(size=40)) + ggtitle("Delay \n 5.6d") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=labz) + theme(axis.text.x=element_text(angle=45, hjust=1)) + scale_y_continuous("log10(Ctrough Adduct)") + xlab(" ") 

missed3 <- subset(df4, select=c('personIndex','CtroughAdduct.total.5', 'CtroughAdduct.total.6'))
meltedBox3 <- melt(missed3, id.vars="personIndex")
meltedBox3$variable <- factor(meltedBox3$variable, levels = c('CtroughAdduct.total.5', 'CtroughAdduct.total.6'), ordered=TRUE)
labz <- c('MGMT Expressing ', 'MGMT Silenced')
box9 <- ggplot(meltedBox3, aes(x=meltedBox3$variable, y=meltedBox3$value)) + geom_boxplot(fill=c("indianred2", "#CCCCCC")) + theme_bw() + theme(text = element_text(size=40)) + ggtitle("Delay \n 11.2d") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=labz) + theme(axis.text.x=element_text(angle=45, hjust=1)) + scale_y_continuous("log10(Ctrough Adduct)") + xlab(" ") 

missed4 <- subset(df4, select=c('personIndex','CtroughAdduct.total.7', 'CtroughAdduct.total.8'))
meltedBox4 <- melt(missed4, id.vars="personIndex")
meltedBox4$variable <- factor(meltedBox4$variable, levels = c('CtroughAdduct.total.7', 'CtroughAdduct.total.8'), ordered=TRUE)
labz <- c('MGMT Expressing ', 'MGMT Silenced')
box10 <- ggplot(meltedBox4, aes(x=meltedBox4$variable, y=meltedBox4$value)) + geom_boxplot(fill=c("indianred2", "#CCCCCC")) + theme_bw() + theme(text = element_text(size=40)) + ggtitle("Delay \n 16.8d") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=labz) + theme(axis.text.x=element_text(angle=45, hjust=1)) + scale_y_continuous("log10(Ctrough Adduct)") + xlab(" ") 

missed5 <- subset(df4, select=c('personIndex','CtroughAdduct.total.9', 'CtroughAdduct.total.10'))
meltedBox5 <- melt(missed5, id.vars="personIndex")
meltedBox5$variable <- factor(meltedBox5$variable, levels = c('CtroughAdduct.total.9', 'CtroughAdduct.total.10'), ordered=TRUE)
labz <- c('MGMT Expressing ', 'MGMT Silenced')
box11 <- ggplot(meltedBox5, aes(x=meltedBox5$variable, y=meltedBox5$value)) + geom_boxplot(fill=c("indianred2", "#CCCCCC")) + theme_bw() + theme(text = element_text(size=40)) + ggtitle("Delay \n 22.4d") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=labz) + theme(axis.text.x=element_text(angle=45, hjust=1)) + scale_y_continuous("log10(Ctrough Adduct)") + xlab(" ") 

missed6 <- subset(df4, select=c('personIndex','CtroughAdduct.total.11', 'CtroughAdduct.total.12'))
meltedBox6 <- melt(missed6, id.vars="personIndex")
meltedBox6$variable <- factor(meltedBox6$variable, levels = c('CtroughAdduct.total.11', 'CtroughAdduct.total.12'), ordered=TRUE)
labz <- c('MGMT Expressing ', 'MGMT Silenced')
box12 <- ggplot(meltedBox6, aes(x=meltedBox6$variable, y=meltedBox6$value)) + geom_boxplot(fill=c("indianred2", "#CCCCCC")) + theme_bw() + theme(text = element_text(size=40)) + ggtitle("Double Dose \n After Skip") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=labz) + theme(axis.text.x=element_text(angle=45, hjust=1)) + scale_y_continuous("log10(Ctrough Adduct)") + xlab(" ") 

missed7 <- subset(df4, select=c('personIndex','CtroughAdduct.total.13', 'CtroughAdduct.total.14'))
meltedBox7 <- melt(missed7, id.vars="personIndex")
meltedBox7$variable <- factor(meltedBox7$variable, levels = c('CtroughAdduct.total.13', 'CtroughAdduct.total.14'), ordered=TRUE)
labz <- c('MGMT Expressing ', 'MGMT Silenced')
box13 <- ggplot(meltedBox7, aes(x=meltedBox7$variable, y=meltedBox7$value)) + geom_boxplot(fill=c("indianred2", "#CCCCCC")) + theme_bw() + theme(text = element_text(size=40)) + ggtitle("Control \n") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=labz) + theme(axis.text.x=element_text(angle=45, hjust=1)) + scale_y_continuous("log10(Ctrough Adduct)") + xlab(" ") 

out <- grid.arrange(box13, box7, box8, box9, box10, box11, box12, ncol=7)
ggsave(file="SupFig6_DoseAnalysis2.pdf",plot
       =out,width=49,height=16) 
