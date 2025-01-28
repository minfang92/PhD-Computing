########################################################
#SetUP !!!
########################################################
rm(list=ls())
setwd("~/Dropbox/02.Projects/2018.10.04_JMP/03.Model/FHANK-Ss_short_paper/2steadystate")
#Pkg
library(ggplot2)
library(tidyr)
library(tibble)
library(hrbrthemes)
library(dplyr)
library(latex2exp)
library(RColorBrewer)


########################################################
# Investment Rate Plot
########################################################
#Load data
inv <- read.csv(file = "InvestmentRate_bundled.csv") 

#Wide to Long
inv <- inv %>%
  # Data wrangling
  as_tibble() %>%
  rowid_to_column(var="k") %>%
  gather(key="z", value="value", -1) %>%
  mutate(z=as.numeric(gsub("x","",z))) 

#Factor
summary(inv$value)
inv <- inv %>%
  # create a new variable from value
  mutate(valuefactor=cut(value,breaks=c(min(value,na.rm=T)-0.1,-0.01,0.01,0.10,0.20,0.30,0.40,max(value,na.rm=T)+0.1),
                         labels=c("<0","~0%","0-10%","10-20%","20-30%","30-40%",">40%"))) %>%
  # change level order
  mutate(cvaluefactor=factor(as.character(valuefactor),levels=rev(levels(valuefactor))))

#Plot
textcol <- "grey40"
p <- ggplot(inv,aes(x=z,y=k,fill=valuefactor)) + 
  geom_tile(colour="white",size=0.1)+
  guides(fill=guide_legend(title="Investment\nRate"))+
  labs(x="Productivity Grid",y="Capital Stock Grid",title="Investment Rate at Steady State (Bundled Model)")+
  scale_fill_manual(values=c("#00AFBB","#ddf1da","#fdae61","#f46d43","#FC4E07","#d53e4f","#CC79A7"))+
  #coord_fixed()+
  theme_grey(base_size=15)+
  theme(legend.position="right",legend.direction="vertical",
        legend.title=element_text(colour=textcol),
        legend.margin=margin(grid::unit(0,"cm")),
        legend.text=element_text(colour=textcol,size=7,face="bold"),
        legend.key.height=grid::unit(0.8,"cm"),
        legend.key.width=grid::unit(0.2,"cm"),
        axis.text.x=element_text(size=10,colour=textcol),
        axis.text.y=element_text(vjust=0.2,colour=textcol),
        axis.ticks=element_line(size=0.4),
        plot.background=element_blank(),
        panel.border=element_blank(),
        plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
        plot.title=element_text(colour=textcol,hjust=0,size=14,face="bold"))
p
#export figure
ggsave(p,filename="SS_Investment_Policy_over_Distribution_bundled.eps",path = "00.Output", height=5.5,width=8.8,units="in",dpi=200)





########################################################
# Prob. of Ajustment Plot
########################################################
#Load data
prob <- read.csv(file = "Prob_of_Adjustment_bundled.csv") 

#Wide to Long
prob <- prob %>%
  # Data wrangling
  as_tibble() %>%
  rowid_to_column(var="k") %>%
  gather(key="z", value="value", -1) %>%
  mutate(z=as.numeric(gsub("x","",z))) 

#Factor
summary(prob$value)
prob <- prob %>%
  # create a new variable from value
  mutate(valuefactor=cut(value,breaks=c(min(value,na.rm=T)-0.1,0.011,0.15,0.30,0.45,0.60,0.75,max(value,na.rm=T)+0.1),
                         labels=c("~0%","0-15%","15-30%","30-45%","45-60%","60-75%",">75%"))) %>%
  # change level order
  mutate(cvaluefactor=factor(as.character(valuefactor),levels=rev(levels(valuefactor))))

#Plot
textcol <- "grey40"
p <- ggplot(prob,aes(x=z,y=k,fill=valuefactor)) + 
  geom_tile(colour="white",size=0.1)+
  guides(fill=guide_legend(title="Adjustment\nProbability"))+
  labs(x="Productivity Grid",y="Capital Stock Grid",title="Adjustment Probability at Steady State (Bundled Model)")+
  scale_fill_manual(values=(brewer.pal(7,"OrRd"))) +
  #coord_fixed()+
  theme_grey(base_size=15)+
  theme(legend.position="right",legend.direction="vertical",
        legend.title=element_text(colour=textcol),
        legend.margin=margin(grid::unit(0,"cm")),
        legend.text=element_text(colour=textcol,size=7,face="bold"),
        legend.key.height=grid::unit(0.8,"cm"),
        legend.key.width=grid::unit(0.2,"cm"),
        axis.text.x=element_text(size=10,colour=textcol),
        axis.text.y=element_text(vjust=0.2,colour=textcol),
        axis.ticks=element_line(size=0.4),
        plot.background=element_blank(),
        panel.border=element_blank(),
        plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
        plot.title=element_text(colour=textcol,hjust=0,size=14,face="bold"))
p
#export figure
ggsave(p,filename="SS_Adjustment_Policy_over_Distribution_bundled.eps",path = "00.Output", height=5.5,width=8.8,units="in",dpi=200)







