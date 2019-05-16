# Vinay, Kayla, and Claire
# April 20, 2019
library(ggplot2) 
library(R.matlab)
library(reshape)
library(grid)
library(gridExtra)

# Figure 4 Full PK/PD Variability
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# MGMT silenced
data_Plasma_TMZ_AUC_MGMT_silenced <- readMat("full_plasT_AUC_silence.mat",header=T,numerals='no.loss') # from mat file
data_Tumor_TMZ_AUC_MGMT_silenced <- readMat("full_tumT_AUC_silence.mat",header=T,numerals='no.loss')
data_Tumor_Adduct_AUC_MGMT_silenced <- readMat("full_tumorAdd_AUC_silence.mat",header=T,numerals='no.loss')
data_Ctrough_MGMT_silenced <- readMat("full_CtroughAdduct_silence.mat",header=T,numerals='no.loss')
data_Tvol_MGMT_silenced <- readMat("full_finalVol_silence.mat",header=T,numerals='no.loss')

# MGMT not silenced
data_Plasma_TMZ_AUC_MGMT <- readMat("full_plasT_AUC_mgmt.mat",header=T,numerals='no.loss') 
data_Tumor_TMZ_AUC_MGMT <- readMat("full_tumT_AUC_mgmt.mat",header=T,numerals='no.loss')
data_Tumor_Adduct_AUC_MGMT <- readMat("full_tumorAdd_AUC_mgmt.mat",header=T,numerals='no.loss')
data_Ctrough_MGMT <- readMat("full_CtroughAdduct_mgmt.mat",header=T,numerals='no.loss') 
data_Tvol_MGMT <- readMat("full_finalVol_mgmt.mat",header=T,numerals='no.loss')


# Transform into data frames
# MGMT silenced
Plasma_TMZ_AUC_MGMT_silenced <- log10(10e6*abs(as.data.frame(data_Plasma_TMZ_AUC_MGMT_silenced)))
Tumor_TMZ_AUC_MGMT_silenced <- log10(abs(as.data.frame(data_Tumor_TMZ_AUC_MGMT_silenced)))
Tumor_Adduct_AUC_MGMT_silenced <- log10(abs(as.data.frame(data_Tumor_Adduct_AUC_MGMT_silenced)))
Ctrough_MGMT_silenced <- log10(abs(as.data.frame(data_Ctrough_MGMT_silenced)))
Tvol_MGMT_silenced<- as.data.frame(data_Tvol_MGMT_silenced)

# MGMT not silenced
Plasma_TMZ_AUC_MGMT <- log10(abs(as.data.frame(data_Plasma_TMZ_AUC_MGMT)))
Tumor_TMZ_AUC_MGMT <- log10(abs(as.data.frame(data_Tumor_TMZ_AUC_MGMT)))
Tumor_Adduct_AUC_MGMT <- log10(abs(as.data.frame(data_Tumor_Adduct_AUC_MGMT)))
Ctrough_MGMT <- log10(abs(as.data.frame(data_Ctrough_MGMT)))
Tvol_MGMT <- as.data.frame(data_Tvol_MGMT)

# Melt them together
pop_a_Plasma_TMZ_AUC <- data.frame(Plasma_TMZ_AUC_MGMT,Plasma_TMZ_AUC_MGMT_silenced)
colnames(pop_a_Plasma_TMZ_AUC) <- c("Patient","MGMT Expressed","Patient","MGMT Silenced")
pop_a_Plasma_TMZ_AUC <- melt(pop_a_Plasma_TMZ_AUC,id.vars="Patient")

pop_a_Tumor_TMZ_AUC <- data.frame(Tumor_TMZ_AUC_MGMT,Tumor_TMZ_AUC_MGMT_silenced)
colnames(pop_a_Tumor_TMZ_AUC) <- c("Patient","MGMT Expressed","Patient","MGMT Silenced")
pop_a_Tumor_TMZ_AUC <- melt(pop_a_Tumor_TMZ_AUC,id.vars="Patient")

pop_a_Tumor_Adduct_AUC <- data.frame(Tumor_Adduct_AUC_MGMT,Tumor_Adduct_AUC_MGMT_silenced)
colnames(pop_a_Tumor_Adduct_AUC) <- c("Patient","MGMT Expressed","Patient","MGMT Silenced")
pop_a_Tumor_Adduct_AUC <- melt(pop_a_Tumor_Adduct_AUC,id.vars="Patient")

pop_a_Ctrough <- data.frame(Ctrough_MGMT,Ctrough_MGMT_silenced)
colnames(pop_a_Ctrough) <- c("MGMT Expressed","Patient","MGMT Silenced","Patient")
pop_a_Ctrough <- melt(pop_a_Ctrough,id.vars="Patient")

pop_a_Tvol <- data.frame(Tvol_MGMT,Tvol_MGMT_silenced)
colnames(pop_a_Tvol) <- c("MGMT Expressed","Patient","MGMT Silenced","Patient")
pop_a_Tvol <- melt(pop_a_Tvol,id.vars="Patient")


box_4a <- ggplot(pop_a_Plasma_TMZ_AUC,aes(x=variable,y=value,fill=variable)) +
  geom_boxplot() +
  scale_fill_manual(values=c("indianred2","#CCCCCC")) +
  theme(text=element_text(size=8)) +
  xlab("") +
  ylab("log10(Plasma TMZ Exposure)") +
  scale_y_continuous(limits=c(-5.979,-5.974),breaks=c(-5.979,-5.978,-5.977,-5.976,-5.975,-5.974)) +
  theme_bw() +
  guides(fill=FALSE)

box_4b <- ggplot(pop_a_Tumor_TMZ_AUC,aes(x=variable,y=value,fill=variable)) +
  geom_boxplot() +
  scale_fill_manual(values=c("indianred2","#CCCCCC")) +
  theme(text=element_text(size=8)) +
  xlab("") +
  ylab("log10(Tumor TMZ Exposure)") +
  scale_y_continuous(limits=c(-7.5,-2.0),breaks=c(-7.0,-6.0,-5.0,-4.0,-3.0,-2.0)) +
  theme_bw() +
  guides(fill=FALSE)

box_4c <- ggplot(pop_a_Tumor_Adduct_AUC,aes(x=variable,y=value,fill=variable)) +
  geom_boxplot() +
  scale_fill_manual(values=c("indianred2","#CCCCCC")) +
  theme(text=element_text(size=8)) +
  xlab("") +
  ylab("log10(Tumor DNA Adduct Exposure)") +
  scale_y_continuous(limits=c(-14.0,-2.0),breaks=c(-14.0,-12.0,-10.0,-8.0,-6.0,-4.0,-2.0)) +
  theme_bw() +
  guides(fill=FALSE)

box_4d <- ggplot(pop_a_Ctrough,aes(x=variable,y=value,fill=variable)) +
  geom_boxplot() +
  scale_fill_manual(values=c("indianred2","#CCCCCC")) +
  theme(text=element_text(size=8)) +
  xlab("") +
  ylab("log10(Ctrough of Tumor DNA Adduct)") +
  scale_y_continuous(limits=c(-25,-3),breaks=c(-25,-20,-15,-10,-5)) +
  theme_bw() +
  guides(fill=FALSE)

box_4e <- ggplot(pop_a_Tvol,aes(x=variable,y=value,fill=variable)) +
  geom_boxplot() +
  scale_fill_manual(values=c("indianred2","#CCCCCC")) +
  theme(text=element_text(size=8)) +
  xlab("") +
  ylab("Final Tumor Volume (mL)") +
  scale_y_continuous(limits=c(0,350),breaks=c(0,50,100,150,200,250,300,350)) +
  theme_bw() +
  guides(fill=FALSE)

out_6 <- grid.arrange(box_4a,box_4b,box_4c,box_4d,box_4e,nrow=2,top=textGrob("Full Parameter Variability",gp=gpar(fontsize=15,font=3)))
ggsave(file="Fig4_FullVariability.pdf",plot=out_6,width=10,height=7) # Save the figure   

