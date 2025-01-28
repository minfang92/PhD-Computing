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
# XIbar Elasiticity variation Plot
########################################################
#Load data
data <- read.csv(file = "02.Identification_results/Identification_df_Auto1_wage.csv") 
xibar <- c("phi_k","xi_bar","els")
data1 <- data[xibar]
data1 <- mutate(data1, ID = group_indices(data1,phi_k))
data1$els = data1$els
data1$ID <- as.character(data1$ID)
typeof(data1$ID)

#For slides
data1_benchmark <- data1[ which(data1$ID==1), ]


#Color Legends
cols <- c( "blue")
types <- c( "solid" )
group_labels <- c( "Bundled Model" )


# Plot with others
p1 <- ggplot(data1_benchmark, aes(x=xi_bar, y=els, group=ID)) +
  geom_smooth(aes(linetype=ID, color=ID),size=1.6) + 
  geom_hline(yintercept=-5.0, linetype="dashed", color = "red", size=1.0) +
  geom_vline(xintercept=0.6, linetype='dashed') + 
  ylab("PE wage elasticity of inv.") + 
  xlab(TeX('$\\bar{\\xi}$')) +
  theme_linedraw(base_size=20) +
  theme(legend.position = c(0.7,0.25)) +
  theme(legend.text = element_text( size = 20)) +
  scale_colour_manual(name = "",labels = group_labels,values = cols) +
  scale_linetype_manual(name = "",labels = group_labels,values = types)
p1  
ggsave(file="00.Output/Identification_wage_sensitivity_xibar.eps", device = cairo_ps,
       height=8,width=8,units="in",dpi=200)






