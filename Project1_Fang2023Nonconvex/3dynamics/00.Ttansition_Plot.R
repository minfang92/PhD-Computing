########################################################
#SetUP !!!
########################################################
rm(list=ls())
setwd("~/Dropbox/02.Projects/2018.10.04_JMP/03.Model/FHANK-Ss_short_paper/3dynamics")

#Pkg
library(ggplot2)
library(tidyr)
library(tibble)
library(hrbrthemes)
library(dplyr)
library(latex2exp)
library(RColorBrewer)

#Color Legends
cols <- c("#FC4E07", "#FC4E07", "steelblue", "steelblue")
types <- c( "solid", "dotdash", "solid", "dotdash")
group_labels <- c("Recession, Large", "Recession, Small",  "Boom, Large", "Boom, Small")

########################################################
# bundled
########################################################

#Load data
data <- read.csv(file = "transition_bundled.csv") 

inv <- data %>%
  # Data wrangling
  as_tibble() %>%
  rowid_to_column(var="t") %>%
  gather(key="ID", value="value", -1) %>%
  subset(.,t<=9) %>%
  subset(.,t>1)

#-------------------------------------------------
# Plot
p <- ggplot(inv, aes(x=t-1, y=value, group=ID, shape=ID)) +
  geom_line(aes(linetype=ID, color=ID),size=1.6) + 
  geom_point(aes(color=ID),size=4) + 
  theme_linedraw(base_size=20) +
  ylab("Impulse Response of Inv.(absolute value in %)") + 
  xlab("Quarters") +
  theme(legend.position = c(0.7,0.8)) +
  theme(legend.text = element_text( size = 12)) +
  theme(legend.key.size =  unit(0.5, "in"),legend.key = element_blank()) +
  scale_colour_manual(name = "",labels = group_labels,values = cols) +
  scale_linetype_manual(name = "",labels = group_labels,values = types) +
  scale_shape_manual(name = "",labels = group_labels,values = c(15,16,17,18))
p  

#export figure
ggsave(p,filename="Transition_bundled.eps", path="00.Output", height=8,width=5,units="in",dpi=200)


########################################################
# zero mu
########################################################

#Load data
data <- read.csv(file = "transition_zero_mu.csv") 

inv <- data %>%
  # Data wrangling
  as_tibble() %>%
  rowid_to_column(var="t") %>%
  gather(key="ID", value="value", -1) %>%
  subset(.,t<=9) %>%
  subset(.,t>1)

#-------------------------------------------------
# Plot
p <- ggplot(inv, aes(x=t-1, y=value, group=ID, shape=ID)) +
  geom_line(aes(linetype=ID, color=ID),size=1.6) + 
  geom_point(aes(color=ID),size=4) + 
  theme_linedraw(base_size=20) +
  ylab("Impulse Response of Inv.(absolute value in %)") + 
  xlab("Quarters") +
  theme(legend.position = c(0.7,0.8)) +
  theme(legend.text = element_text( size = 12)) +
  theme(legend.key.size =  unit(0.5, "in"),legend.key = element_blank()) +
  scale_colour_manual(name = "",labels = group_labels,values = cols) +
  scale_linetype_manual(name = "",labels = group_labels,values = types) +
  scale_shape_manual(name = "",labels = group_labels,values = c(15,16,17,18))
p  

#export figure
ggsave(p,filename="Transition_zero_mu.eps", path="00.Output", height=8,width=5,units="in",dpi=200)


########################################################
# zero sigma
########################################################

#Load data
data <- read.csv(file = "transition_zero_sigma.csv") 

inv <- data %>%
  # Data wrangling
  as_tibble() %>%
  rowid_to_column(var="t") %>%
  gather(key="ID", value="value", -1) %>%
  subset(.,t<=9) %>%
  subset(.,t>1)

#-------------------------------------------------
# Plot
p <- ggplot(inv, aes(x=t-1, y=value, group=ID, shape=ID)) +
  geom_line(aes(linetype=ID, color=ID),size=1.6) + 
  geom_point(aes(color=ID),size=4) + 
  theme_linedraw(base_size=20) +
  ylab("Impulse Response of Inv.(absolute value in %)") + 
  xlab("Quarters") +
  theme(legend.position = c(0.7,0.8)) +
  theme(legend.text = element_text( size = 12)) +
  theme(legend.key.size =  unit(0.5, "in"),legend.key = element_blank()) +
  scale_colour_manual(name = "",labels = group_labels,values = cols) +
  scale_linetype_manual(name = "",labels = group_labels,values = types) +
  scale_shape_manual(name = "",labels = group_labels,values = c(15,16,17,18))
p  

#export figure
ggsave(p,filename="Transition_zero_sigma.eps", path="00.Output", height=8,width=5,units="in",dpi=200)

