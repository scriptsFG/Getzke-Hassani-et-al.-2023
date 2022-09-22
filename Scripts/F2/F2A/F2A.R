
library(tidyverse)
library(ggpubr)
library(rstatix)
#######

rm(list = ls())

colors<-c("Agaricomycetes"="#bfbfbf", "Ascomycota"="#8B008B", "Dothideomycetes"="#DB7093", "Eurotiomycetes"="#bfbfbf", "Glomeromycetes"="#bfbfbf", "Lecanoromycetes"="#bfbfbf", "Leotiomycetes"="#bfbfbf", "Microbotryomycetes"="#bfbfbf", "Mucoromycotina"="#bfbfbf", "Pezizomycetes"="#bfbfbf", "Pezizomycotina"="#bfbfbf", "Saccharomycetes"="#bfbfbf", "Sordariomycetes"="#8B008B", "Tremellomycetes"="#bfbfbf", "[Saprospirae]"="#bfbfbf", "Acidimicrobiia"="#bfbfbf",  "Actinobacteria"="#B22222",  "Alphaproteobacteria"="#32CD32", "Bacilli"="#FFA500",  "Bacteroidia"="#bfbfbf",  "Betaproteobacteria"="#228B22",  "Chlamydiia"="#bfbfbf",  "Chloroflexi"="#bfbfbf", 
           "Clostridia"="#bfbfbf",  "Cytophagia"="#bfbfbf",  "DA052"="#bfbfbf",  "Deltaproteobacteria"="#bfbfbf", 
           "Fibrobacteria"="#bfbfbf",  "Flavobacteriia"="#4682B4",  "Gammaproteobacteria"="#006400",  "Gemm-1"="#bfbfbf",  "Nitrospira"="#bfbfbf", 
           "Rubrobacteria"="#bfbfbf",  "SC3"="#bfbfbf",  "Sphingobacteriia"="#bfbfbf",  "Thermoleophilia"="#bfbfbf",  "TK10"="#bfbfbf",  "TM7-3"="#bfbfbf",  "VHS-B5-50"="#bfbfbf")


#read meta data
read_tsv("F2A_data.txt")->nBGCs
read_tsv("interaction_producers.txt")->producers
producers[,c(1,5)]->producers2

#summerize most abundant BGCs by Class
filter(nBGCs,Genomes!="Root401")->nBGCs
filter(nBGCs,Culture_collection=="AtRoot")->AtRoot

AtRoot %>% 
  select_if(is_numeric) %>% 
  select_if(colSums(.) > 19.5) %>% 
  bind_cols("class"=AtRoot$class, .) %>% 
  group_by(class) %>% 
  summarise(across(where(is.numeric), ~ mean(.)))->AtRoot_class
AtRoot_class$class->AtRoot_class$Genomes
AtRoot[, c(1, 3:56)]->AtRoot2
left_join(producers2, AtRoot2)->producers3
bind_rows(AtRoot_class, producers3)->AtRoot_class
AtRoot_class %>%
  select_if(~ !any(is.na(.)))->AtRoot_class
pivot_longer(AtRoot_class, cols = 2:21, names_to ="BGC", values_to = "n")->AtRoot_class


read_tsv("At_R_isolates.txt")->subset
left_join(subset, AtRoot)->subset
subset[, 3:56]->subset
subset %>% 
  select_if(colSums(.) > 0)->subset
subset %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))->subset
rowSums(subset[1,])->subset


filter(nBGCs,Culture_collection=="AtLeaf")->AtLeaf
AtLeaf[, 3:56]->AtLeaf
AtLeaf %>% 
  select_if(colSums(.) > 0)->AtLeaf
AtLeaf %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))->AtLeaf
rowSums(AtLeaf[1,])->AtLeaf



ggballoonplot(data=AtRoot_class, x="BGC", y="Genomes", fill="class",size="n", size.range = c(0.00000001, 12))+
  ylab("")+
  xlab("")+
  #ggtitle("Halo width of R401(-mutant) on different lawn bacteria")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.border = element_rect(colour = 'black', size=1, fill = NA),
        axis.text=element_text(size=12, color= "black"),
        axis.text.x = element_text(size=12, color= "black", angle = 45, hjust=0.99),
        axis.title=element_text(size=12),
        legend.position = "none", 
        strip.background = element_rect(color="black", fill=NA, size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black"))+
  scale_y_discrete(limits=rev(c( "Root63", "Root1310","Root920", "Root342", "Root690", "Root569", "Root71", "Root68","Root401_Pacbio", "Root562", "Actinobacteria", "Bacilli","Flavobacteriia", "Alphaproteobacteria", "Betaproteobacteria", "Gammaproteobacteria")),labels=rev(c( "R63", "R1310", "R920", "R342","R690", "R569", "R71", "R68","R401", "R562", "Actinobacteria", "Bacilli","Flavobacteriia", "Alphaproteobacteria", "Betaproteobacteria", "Gammaproteobacteria")))+
  #scale_x_discrete(limits=c("56", "DD", "56DD"), labels=c("dNRPS", "dDAPG",  "dNRPSdDAPG"))+
  #scale_x_discrete(limits=c("54","56", "DD", "54DD","56DD"), labels=c("dAcT","dNRPS", "dDAPG", "dAcTdDAPG", "dNRPSdDAPG"))+
  scale_fill_manual(values=colors)+
  NULL
#ggsave("BGC_balloon.pdf")
#ggsave("BGC_balloon_full.pdf")


AtRoot_class <- transform(AtRoot_class, Group = ifelse(str_detect(Genomes, "Root"), "Strain", "Class"))
AtRoot_class %>% 
  group_by(BGC) %>% 
  dunn_test(n~Group, p.adjust.method ="BH") %>% 
  select(BGC, p.adj)->test
test


#compute nBGC and nBGC families
bind_cols("Genomes"=AtRoot$Genomes, "class"=AtRoot$class,"nBGCs"=rowSums(AtRoot[ , 3:56], na.rm=TRUE),"nBGC_Families"=rowSums(AtRoot[ , 3:56] > 0))->AtRoot_sums
AtRoot_sums %>% 
  group_by(class) %>% 
  summarise("nBGCs"=mean(nBGCs),"nBGC_Families"=mean(nBGC_Families))->AtRoot_sums2
AtRoot_sums2$class->AtRoot_sums2$Genomes
left_join(producers2, AtRoot_sums)->producers4
bind_rows(AtRoot_sums2, producers4)->AtRoot_sums3
pivot_longer(AtRoot_sums3, cols = 2:3, names_to ="BGC", values_to = "n")->AtRoot_sums3


AtRoot_sums3 <- transform(AtRoot_sums3, Group = ifelse(str_detect(Genomes, "Root"), "Strain", "Class"))
AtRoot_sums3 %>% 
  group_by(BGC) %>% 
  dunn_test(n~Group, p.adjust.method ="BH") %>% 
  select(BGC, p.adj)->test
test


ggballoonplot(data=AtRoot_sums3, x="BGC", y="Genomes", fill="grey",size="n", size.range = c(0, 15))+
  ylab("")+
  xlab("")+
  #ggtitle("Halo width of R401(-mutant) on different lawn bacteria")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.border = element_rect(colour = 'black', size=1, fill = NA),
        axis.text=element_text(size=12, color= "black"),
        axis.text.x = element_text(size=12, color= "black", angle = 45, hjust=0.99),
        axis.title=element_text(size=12),
        legend.position = "none", 
        strip.background = element_rect(color="black", fill=NA, size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black"))+
  scale_y_discrete(limits=rev(c( "Root63", "Root1310","Root920", "Root342", "Root690", "Root569", "Root71", "Root68","Root401_Pacbio", "Root562", "Actinobacteria", "Bacilli","Flavobacteriia", "Alphaproteobacteria", "Betaproteobacteria", "Gammaproteobacteria")),labels=rev(c( "R63", "R1310", "R920", "R342","R690", "R569", "R71", "R68","R401", "R562", "Actinobacteria", "Bacilli","Flavobacteriia", "Alphaproteobacteria", "Betaproteobacteria", "Gammaproteobacteria")))+
  #scale_x_discrete(limits=c("56", "DD", "56DD"), labels=c("dNRPS", "dDAPG",  "dNRPSdDAPG"))+
  #scale_x_discrete(limits=c("54","56", "DD", "54DD","56DD"), labels=c("dAcT","dNRPS", "dDAPG", "dAcTdDAPG", "dNRPSdDAPG"))+
  #scale_fill_manual(values=colors)+
  NULL
#ggsave("BGC_balloon2.pdf")
#ggsave("BGC_balloon_full2.pdf")








