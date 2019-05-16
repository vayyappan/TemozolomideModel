# Vinay, Kayla, and Claire
# April 20, 2019
library(ggplot2) 
library(R.matlab)
library(reshape)
library(grid)
library(gridExtra)

# Supplementary Figure 5 Visualization - Initial Tumor Size Variability
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Import data.

# MGMT silenced
data_Plasma_TMZ_AUC_MGMT_silenced <- readMat("size_plasT_AUC_silence.mat",header=T,numerals='no.loss') # from mat file
data_Tumor_TMZ_AUC_MGMT_silenced <- readMat("size_tumT_AUC_silence.mat",header=T,numerals='no.loss')
data_Tumor_Adduct_AUC_MGMT_silenced <- readMat("size_tumorAdd_AUC_silence.mat",header=T,numerals='no.loss')
data_Ctrough_MGMT_silenced <- readMat("size_CtroughAdduct_silence.mat",header=T,numerals='no.loss')
data_Tvol_MGMT_silenced <- readMat("size_finalVol_silence.mat",header=T,numerals='no.loss')

# MGMT not silenced
data_Plasma_TMZ_AUC_MGMT <- readMat("size_plasT_AUC_mgmt.mat",header=T,numerals='no.loss') 
data_Tumor_TMZ_AUC_MGMT <- readMat("size_tumT_AUC_mgmt.mat",header=T,numerals='no.loss')
data_Tumor_Adduct_AUC_MGMT <- readMat("size_tumorAdd_AUC_mgmt.mat",header=T,numerals='no.loss')
data_Ctrough_MGMT <- readMat("size_CtroughAdduct_mgmt.mat",header=T,numerals='no.loss') 
data_Tvol_MGMT <- readMat("size_finalVol_mgmt.mat",header=T,numerals='no.loss')


# Transform into data frames
# MGMT silenced
Plasma_TMZ_AUC_MGMT_silenced <- log10(as.data.frame(data_Plasma_TMZ_AUC_MGMT_silenced))
Tumor_TMZ_AUC_MGMT_silenced <- log10(as.data.frame(data_Tumor_TMZ_AUC_MGMT_silenced))
Tumor_Adduct_AUC_MGMT_silenced <- log10(as.data.frame(data_Tumor_Adduct_AUC_MGMT_silenced))
Ctrough_MGMT_silenced <- log10(abs(as.data.frame(data_Ctrough_MGMT_silenced)))
Tvol_MGMT_silenced<- as.data.frame(data_Tvol_MGMT_silenced)

# MGMT not silenced
Plasma_TMZ_AUC_MGMT <- log10(as.data.frame(data_Plasma_TMZ_AUC_MGMT))
Tumor_TMZ_AUC_MGMT <- log10(as.data.frame(data_Tumor_TMZ_AUC_MGMT))
Tumor_Adduct_AUC_MGMT <- log10(as.data.frame(data_Tumor_Adduct_AUC_MGMT))
Ctrough_MGMT <- log10(abs(as.data.frame(data_Ctrough_MGMT)))
Tvol_MGMT <- as.data.frame(data_Tvol_MGMT)

# Melt them together
pop_ts_Plasma_TMZ_AUC <- data.frame(Plasma_TMZ_AUC_MGMT,Plasma_TMZ_AUC_MGMT_silenced)
colnames(pop_ts_Plasma_TMZ_AUC) <- c("Patient","MGMT Expressed","Patient","MGMT Silenced")
pop_ts_Plasma_TMZ_AUC <- melt(pop_ts_Plasma_TMZ_AUC,id.vars="Patient")

pop_ts_Tumor_TMZ_AUC <- data.frame(Tumor_TMZ_AUC_MGMT,Tumor_TMZ_AUC_MGMT_silenced)
colnames(pop_ts_Tumor_TMZ_AUC) <- c("Patient","MGMT Expressed","Patient","MGMT Silenced")
pop_ts_Tumor_TMZ_AUC <- melt(pop_ts_Tumor_TMZ_AUC,id.vars="Patient")

pop_ts_Tumor_Adduct_AUC <- data.frame(Tumor_Adduct_AUC_MGMT,Tumor_Adduct_AUC_MGMT_silenced)
colnames(pop_ts_Tumor_Adduct_AUC) <- c("Patient","MGMT Expressed","Patient","MGMT Silenced")
pop_ts_Tumor_Adduct_AUC <- melt(pop_ts_Tumor_Adduct_AUC,id.vars="Patient")

pop_ts_Ctrough <- data.frame(Ctrough_MGMT,Ctrough_MGMT_silenced)
colnames(pop_ts_Ctrough) <- c("MGMT Expressed","Patient","MGMT Silenced","Patient")
pop_ts_Ctrough <- melt(pop_ts_Ctrough,id.vars="Patient")

pop_ts_Tvol <- data.frame(Tvol_MGMT,Tvol_MGMT_silenced)
colnames(pop_ts_Tvol) <- c("MGMT Expressed","Patient","MGMT Silenced","Patient")
pop_ts_Tvol <- melt(pop_ts_Tvol,id.vars="Patient")


box_5a <- ggplot(pop_ts_Plasma_TMZ_AUC,aes(x=variable,y=value,fill=variable)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#CC0033","#CCCCCC")) +
  theme(text=element_text(size=8)) +
  xlab("") +
  ylab("log10(Plasma TMZ AUC)") +
  scale_y_continuous(limits=c(-5.978057,-5.978043),breaks=c(-5.978057,-5.978055,-5.978053,-5.978051,-5.978049,-5.978047,-5.978045,-5.978043)) +
  theme_bw() +
  guides(fill=FALSE)

box_5b <- ggplot(pop_ts_Tumor_TMZ_AUC,aes(x=variable,y=value,fill=variable)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#CC0033","#CCCCCC")) +
  theme(text=element_text(size=8)) +
  xlab("") +
  ylab("log10(Tumor TMZ AUC)") +
  scale_y_continuous(limits=c(-7.5,-2.0),breaks=c(-7.0,-6.0,-5.0,-4.0,-3.0,-2.0)) +
  theme_bw() +
  guides(fill=FALSE)

box_5c <- ggplot(pop_ts_Tumor_Adduct_AUC,aes(x=variable,y=value,fill=variable)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#CC0033","#CCCCCC")) +
  theme(text=element_text(size=8)) +
  xlab("") +
  ylab("log10(Tumor Adduct AUC)") +
  scale_y_continuous(limits=c(-14,-2),breaks=c(-14,-12,-10,-8,-6,-4,-2)) +
  theme_bw() +
  guides(fill=FALSE)

box_5d <- ggplot(pop_ts_Ctrough,aes(x=variable,y=value,fill=variable)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#CC0033","#CCCCCC")) +
  theme(text=element_text(size=8)) +
  xlab("") +
  ylab("log10(Ctrough)") +
  scale_y_continuous(limits=c(-25,-3),breaks=c(-25,-20,-15,-10,-5)) +
  theme_bw() +
  guides(fill=FALSE)

box_5e <- ggplot(pop_ts_Tvol,aes(x=variable,y=value,fill=variable)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#CC0033","#CCCCCC")) +
  theme(text=element_text(size=8)) +
  xlab("") +
  ylab("Final Tumor Volume (mL)") +
  scale_y_continuous(limits=c(0,350),breaks=c(0,50,100,150,200,250,300,350)) +
  theme_bw() +
  guides(fill=FALSE)

out_5 <- grid.arrange(box_5a,box_5b,box_5c,box_5d,box_5e,nrow=2,top=textGrob("Initial Tumor Size Variability",gp=gpar(fontsize=15,font=3)))
ggsave(file="SupFig5_TumorSize_Variability.pdf",plot=out_5,width=10,height=6) # Save the figure  
