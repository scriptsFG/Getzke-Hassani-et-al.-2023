
library(tidyverse)
library(ggpubr)
library(rstatix)
#######



#presence of DAPG and Pyoverdine 
rm(list = ls())


#read meta data
read_tsv("S3_data.txt")->molecules

#summerize most abundant BGCs by Class
filter(molecules,Genomes!="Root401")->molecules
filter(molecules,Culture_collection=="AtRoot")->AtRoot
filter(AtRoot,genus=="Pseudomonas")->AtRoot
AtRoot %>% 
  select_if(is_numeric) %>% 
  select_if(colSums(.) > 0) %>% 
  bind_cols("Genomes"=AtRoot$Genomes, .)->AtRoot
pivot_longer(AtRoot, cols = 2:20, names_to ="BGC", values_to = "n")->AtRoot

ggballoonplot(data=AtRoot, x="BGC", y="Genomes",size="n", size.range = c(0, 20))+
  ylab("")+
  xlab("")+
  #ggtitle("Halo width of R401(-mutant) on different lawn bacteria")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.border = element_rect(colour = 'black', size=1, fill = NA),
        axis.text=element_text(size=12, color= "black"),
        axis.text.x = element_text(size=12, color= "black", angle = 90, hjust=0.99),
        axis.title=element_text(size=12),
        legend.position = "none", 
        strip.background = element_rect(color="black", fill=NA, size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black"))+
  scale_y_discrete(limits=rev(c("Root569", "Root9","Root562", "Root401_Pacbio", "Root329", "Root71", "Root68")),labels=rev(c("R569", "R9","R562", "R401", "R329", "R71", "R68")))+
  #scale_x_discrete(limits=c("56", "DD", "56DD"), labels=c("dNRPS", "dDAPG",  "dNRPSdDAPG"))+
  #scale_x_discrete(limits=c("54","56", "DD", "54DD","56DD"), labels=c("dAcT","dNRPS", "dDAPG", "dAcTdDAPG", "dNRPSdDAPG"))+
  #scale_fill_manual(values=colors)+
  NULL
#ggsave("Ps_balloon.pdf")
#ggsave("Ps_balloon_full.pdf")




