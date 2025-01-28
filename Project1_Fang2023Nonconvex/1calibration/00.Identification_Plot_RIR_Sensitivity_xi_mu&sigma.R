########################################################
#SetUP !!!
########################################################
rm(list=ls())
setwd("~/Dropbox/02.Projects/2018.10.04_JMP/03.Model/FHANK-Ss_short_paper/1calibration")

#Pkg
library(ggplot2)   # for general plotting
library(dplyr)
library(plyr)      # for fortifying shapefiles
library(viridis)
library(latex2exp)


########################################################
# MuXI Elasiticity variation Plot
########################################################
#Load data
data <- read.csv(file = "02.Identification_results/Identification_df_Auto2.csv") 
muxi <- c("phi_k","muxi","els")
data1 <- data[muxi]
data1 <- mutate(data1, ID = group_indices(data1,phi_k))
data1$els = -data1$els
data1$ID <- as.character(data1$ID)
typeof(data1$ID)

#For slides
data1_benchmark <- data1[ which(data1$ID==1), ]


#Color Legends
cols <- c("steelblue")
types <- c("solid")
group_labels <- c("Variance-Fixed Model")


# Plot for Benchmark Only
p1 <- ggplot(data1_benchmark) +
  geom_smooth(aes(x=muxi, y=els, group=ID, color=ID),size=1.6) + 
  geom_hline(yintercept=-5.0, linetype="dashed", color = "red", size=1.0) +
  geom_vline(xintercept=0.3, linetype='dashed') +  
  ylab("PE real interest rate elasticity of inv.") + 
  xlab(TeX('$\\mu_\\xi}$')) +
  #coord_fixed(ratio.values / ratio.display) +
  #ylim(0.0,15) +
  #xlim(0.0,1.0) +
  theme_linedraw(base_size=20) +
  theme(legend.position = c(0.7,0.25)) +
  theme(legend.text = element_text( size = 20)) +
  scale_colour_manual(name = "",labels = group_labels,values = cols) 
p1  
ggsave(file="00.Output/Identification_RIR_sensitivity_benchmark_mu.eps", device = cairo_ps,
       height=8,width=8,units="in",dpi=200)




########################################################
# SigmaXi Elasiticity variation Plot
########################################################
#Load data
data <- read.csv(file = "02.Identification_results/Identification_df_Auto3.csv") 
muxi <- c("phi_k","sigmaxi","els")
data1 <- data[muxi]
data1 <- mutate(data1, ID = group_indices(data1,phi_k))
data1$els = -data1$els
data1$ID <- as.character(data1$ID)
typeof(data1$ID)

#For slides
data1_benchmark <- data1[ which(data1$ID==1), ]


#Color Legends
cols <- c("#FC4E07")
types <- c("solid")
group_labels <- c("Mean-Fixed Model")


# Plot for Benchmark Only
p1 <- ggplot(data1_benchmark) +
  geom_smooth(aes(x=sigmaxi, y=els, group=ID, color=ID),size=1.6) + 
  geom_hline(yintercept=-5.0, linetype="dashed", color = "red", size=1.0) +
  geom_vline(xintercept=0.17, linetype='dashed') +  
  ylab("PE real interest rate elasticity of inv.") + 
  xlab(TeX('$\\sigma_\\xi}$')) +
  #coord_fixed(ratio.values / ratio.display) +
  #ylim(0.0,15) +
  #xlim(0.0,1.0) +
  theme_linedraw(base_size=20) +
  theme(legend.position = c(0.7,0.25)) +
  theme(legend.text = element_text( size = 20)) +
  scale_colour_manual(name = "",labels = group_labels,values = cols) 
p1  
ggsave(file="00.Output/Identification_RIR_sensitivity_benchmark_sigma.eps", device = cairo_ps,
       height=8,width=8,units="in",dpi=200)




