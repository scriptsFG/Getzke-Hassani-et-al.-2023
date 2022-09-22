library(tidyverse)
library(PMCMRplus)
library("phyloseq")
#######

rm(list = ls())

#####
# 54 = ∆pvdy
# 56 = ∆pvdl
# DD = ∆phld

#colors for plotting
colors<-c("HK"="#C2C2C1", "WT" = "#1F1F1F", "54"="#F94144","56"="#577590","DD"="#F9C74F", "54DD"="#F8961E", "56DD"="#90BE6D")
shapes<-c("I"=21, "II"=22, "III"=24)



#read 16s and ITS1 table
read_tsv(file = "ASV_table.txt")->counts
column_to_rownames(counts, var="...1")->counts



#load mapping file
read_tsv("mapping_file.txt")->map
str_c(map$Compartment, map$Experiment, sep = "_")->map$Compartment_Experiment

#add strain meta data (taxonomy)
read_tsv("tax_file.txt")->tax


#generate phyloseq data format
#otu_tables
otu_table(counts, taxa_are_rows=T)->countsPS


#tax_tables
#tax[is.na(tax)] <- "undetermined"
column_to_rownames(tax, "Isolate")->tax
tax_table(as.matrix(tax))->taxPS

#sample_table
filter(map, map$Sample_Name!="H2O")->map
column_to_rownames(map, "Sample_Name")->mapPS
sample_data(mapPS)->mapPS

#generate phylo_seq objects
phyloseq(countsPS, taxPS, mapPS)->dataPS



#####plot Shannon indeces
prune_samples(dataPS@sam_data$Compartment=="Root", dataPS)->dataPS_R
prune_samples(dataPS@sam_data$Compartment=="Soil", dataPS)->dataPS_S


dataPS_R <- prune_taxa(taxa_sums(dataPS_R) > 0, dataPS_R)
dataPS_S <- prune_taxa(taxa_sums(dataPS_S) > 0, dataPS_S)

plot_richness(dataPS_R, x="R401", measures=c("Shannon"))->ShannonR
plot_richness(dataPS_S, x="R401", measures=c("Shannon"))->ShannonS

ShannonR[["data"]]->ShannonR
ShannonS[["data"]]->ShannonS


ggplot(data=ShannonR, aes(x=R401, y=value, color=R401))+
  #specify the plot
  geom_jitter(aes(fill=R401, color=R401, shape=Experiment), width = 0.2, alpha=0.8, shape=16)+
  geom_boxplot(outlier.shape = NA, lwd=0.9, fill=NA)+
  facet_wrap(~Compartment)+
  xlab("")+
  ylab("Shannon Index")+
  #ggtitle("Halo width of R401(-mutant) on different lawn bacteria")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.border = element_rect(colour = 'black', size=1, fill = NA),
        axis.text=element_text(size=12, color= "black"),
        axis.text.x = element_text(size=12, color= "black", angle = 90, hjust=1),
        axis.title=element_text(size=12), 
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.position="none",
        strip.background = element_rect(color="black", fill=NA, size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black"))+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  scale_shape_manual(values=shapes)+
  scale_x_discrete(limits=c("HK", "WT", "54", "56", "DD", "54DD", "56DD"), labels=c("HK","WT", "dAcT", "dNRPS", "dDAPG", "dAcTdDAPG", "dNRPSdDAPG"))+
  #scale_y_continuous(limits=c(0, 5))+
  NULL
#ggsave("Shannon_boxplot_root.pdf")

dunnR <- kwAllPairsDunnTest(value ~ as.factor(R401), data=ShannonR, p.adjust.method ="BH")
dunnR

ggplot(data=ShannonS, aes(x=R401, y=value, color=R401))+
  #specify the plot
  geom_jitter(aes(fill=R401, color=R401, shape=Experiment), width = 0.2, alpha=0.8, shape=16)+
  geom_boxplot(outlier.shape = NA, lwd=0.9, fill=NA)+
  facet_wrap(~Compartment)+
  xlab("")+
  ylab("Shannon Index")+
  #ggtitle("Halo width of R401(-mutant) on different lawn bacteria")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.border = element_rect(colour = 'black', size=1, fill = NA),
        axis.text=element_text(size=12, color= "black"),
        axis.text.x = element_text(size=12, color= "black", angle = 90, hjust=1),
        axis.title=element_text(size=12), 
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.position="none",
        strip.background = element_rect(color="black", fill=NA, size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black"))+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  scale_shape_manual(values=shapes)+
  scale_x_discrete(limits=c("HK", "WT", "54", "56", "DD", "54DD", "56DD"), labels=c("HK","WT", "dAcT", "dNRPS", "dDAPG", "dAcTdDAPG", "dNRPSdDAPG"))+
  #scale_y_continuous(limits=c(0, 5))+
  NULL
#ggsave("Shannon_boxplot_soil.pdf")

dunnS <- kwAllPairsDunnTest(value ~ as.factor(R401), data=ShannonS, p.adjust.method ="BH")
dunnS


