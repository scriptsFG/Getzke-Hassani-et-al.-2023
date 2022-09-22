library(tidyverse)
library(ggpubr)
library(rstatix)

#######
#######

#comparison across all culture collections - S6D
rm(list = ls())

colors<-c("Leaf"="#bfbfbf", "Root"="#343a40")
#read meta data
read_tsv("summary_BGC_Families_for_felix.txt")->nBGCs
filter(nBGCs,family=="Pseudomonadaceae")->nBGCs
read_tsv("pseudomonas_bgc_family.txt")->Karasov
filter(nBGCs,Genomes!="Root401")->nBGCs


filter(nBGCs, Culture_collection=="AtLeaf")->AtLeaf
filter(nBGCs, Culture_collection!="AtLeaf")->Roots

rep("Root", length(Roots$Culture_collection))->Roots$Organ
rep("Leaf", length(AtLeaf$Culture_collection))->AtLeaf$Organ
rep("Leaf", length(Karasov$Genomes))->Karasov$Organ
rep("Karasov", length(Karasov$Genomes))->Karasov$Culture_collection
bind_rows(Roots, AtLeaf, Karasov)->full


full %>% group_by(Culture_collection) %>%tally()->n_entries



#compute nBGC and nBGC families
bind_cols("Genomes"=full$Genomes, "Organ"=full$Organ, "Culture_collection"=full$Culture_collection, "nBGCs"=rowSums(full[ , 3:56], na.rm=TRUE),"nBGC_Families"=rowSums(full[ , 3:56] > 0, na.rm=TRUE))->full3
full3 %>% 
  group_by(Organ, Culture_collection) %>% 
  summarise("nBGCs"=sum(nBGCs),"nBGC_Families"=sum(nBGC_Families))->Chi
full3 %>% 
  group_by(Organ, Culture_collection) %>% 
  summarise("nBGCs"=mean(nBGCs),"nBGC_Families"=mean(nBGC_Families))->full3

pivot_longer(full3, cols = 3:4, names_to ="BGC", values_to = "n")->full4


ggballoonplot(data=full4, x="BGC", y="Culture_collection", fill='Organ',size="n", size.range = c(0, 15))+
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
        #legend.position = "none", 
        strip.background = element_rect(color="black", fill=NA, size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black"))+
  scale_y_discrete(limits=rev(c("AtLeaf","Karasov","AtRoot","NC", "LjRoot")))+
  #scale_x_discrete(limits=c("56", "DD", "56DD"), labels=c("dNRPS", "dDAPG",  "dNRPSdDAPG"))+
  #scale_x_discrete(limits=c("54","56", "DD", "54DD","56DD"), labels=c("dAcT","dNRPS", "dDAPG", "dAcTdDAPG", "dNRPSdDAPG"))+
  scale_fill_manual(values=colors)+
  NULL
#ggsave("total_BGC_all.pdf")
#ggsave("total_BGC_all_full.pdf")


#Chi Square test

#nBGCs
M <- as.table(rbind(c(33.05, 15.95), c(755.8, 762.2)))
dimnames(M) <- list(Soil = c("Root", "Leaf"),
                    Rescue_category = c("Compound","NoCompound"))
(Xsq <- chisq.test(M))  # Prints test summary
Xsq$observed   # observed counts (same as M)
Xsq$expected   # expected counts under the null
Xsq$residuals  # Pearson residuals
Xsq$stdres     # standardized residuals

#nBGC_families
M <- as.table(rbind(c(21.15, 27.85), c(466.45, 1051.55)))
dimnames(M) <- list(Soil = c("Root", "Leaf"),
                    Rescue_category = c("Compound","NoCompound"))
(Xsq <- chisq.test(M))  # Prints test summary
Xsq$observed   # observed counts (same as M)
Xsq$expected   # expected counts under the null
Xsq$residuals  # Pearson residuals
Xsq$stdres     # standardized residuals








#DAPG and Pyoverdine for all strains - S6E
rm(list = ls())

colors<-c("Leaf"="#bfbfbf", "Root"="#343a40")
#read meta data
read_tsv("summary_molecules_prediction_for_felix.txt")->nBGCs
filter(nBGCs,family=="Pseudomonadaceae")->nBGCs
read_tsv("pseudomonas_bgc_molecules_prediction.txt")->Karasov
filter(nBGCs,Genomes!="Root401")->nBGCs



filter(nBGCs, Culture_collection=="AtLeaf")->AtLeaf
filter(nBGCs, Culture_collection!="AtLeaf")->Roots

rep("Root", length(Roots$Culture_collection))->Roots$Organ
rep("Leaf", length(AtLeaf$Culture_collection))->AtLeaf$Organ
rep("Leaf", length(Karasov$Genomes))->Karasov$Organ
rep("Karasov", length(Karasov$Genomes))->Karasov$Culture_collection
bind_rows(Roots, AtLeaf, Karasov)->full

full %>% group_by(Culture_collection) %>%tally()->n_entries
full %>% 
  select(Genomes, Organ, Culture_collection, class, family, genus, `2,4-diacetylphloroglucinol`, pyoverdin)->full2
full2$`2,4-diacetylphloroglucinol`[is.na(full2$`2,4-diacetylphloroglucinol`)]<- 0

full2 %>% 
  group_by(Organ, Culture_collection) %>% 
  summarise("DAPG"=sum(`2,4-diacetylphloroglucinol`),"pyoverdin"=sum(pyoverdin))->chi

full2 %>% 
  group_by(Organ, Culture_collection) %>% 
  summarise("DAPG"=mean(`2,4-diacetylphloroglucinol`),"pyoverdin"=mean(pyoverdin))->full2

pivot_longer(full2, cols = 3:4, names_to ="BGC", values_to = "n")->full3


full3 %>% 
  group_by(Organ, BGC) %>% 
  summarise(avg=mean(n))->full4

#Root/Leaf
2.3822751/1.5120545

ggballoonplot(data=full3, x="BGC", y="Culture_collection", fill='Organ',size="n", size.range = c(0, 15))+
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
        #legend.position = "none", 
        strip.background = element_rect(color="black", fill=NA, size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black"))+
  scale_y_discrete(limits=rev(c("AtLeaf","Karasov","AtRoot","NC", "LjRoot")))+
  #scale_x_discrete(limits=c("56", "DD", "56DD"), labels=c("dNRPS", "dDAPG",  "dNRPSdDAPG"))+
  #scale_x_discrete(limits=c("54","56", "DD", "54DD","56DD"), labels=c("dAcT","dNRPS", "dDAPG", "dAcTdDAPG", "dNRPSdDAPG"))+
  scale_fill_manual(values=colors)+
  NULL
#ggsave("DAPG_pyoverdin_all.pdf")
#ggsave("DAPG_pyoverdin_all_full.pdf")


#DAPG
M <- as.table(rbind(c(27, 22), c(0, 1518)))
dimnames(M) <- list(Soil = c("Root", "Leaf"),
                    Rescue_category = c("Compound","NoCompound"))
(Xsq <- chisq.test(M))  # Prints test summary
Xsq$observed   # observed counts (same as M)
Xsq$expected   # expected counts under the null
Xsq$residuals  # Pearson residuals
Xsq$stdres     # standardized residuals
#pyoverdin
M <- as.table(rbind(c(38.66666667, 10.33333333), c(735.3333333, 782.6666667)))
dimnames(M) <- list(Soil = c("Root", "Leaf"),
                    Rescue_category = c("Compound","NoCompound"))
(Xsq <- chisq.test(M))  # Prints test summary
Xsq$observed   # observed counts (same as M)
Xsq$expected   # expected counts under the null
Xsq$residuals  # Pearson residuals
Xsq$stdres     # standardized residuals

