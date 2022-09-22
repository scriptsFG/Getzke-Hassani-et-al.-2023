library(tidyverse)
library(rstatix)
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


#subset Root samples from PS data for plotting
prune_samples(dataPS@sam_data$Compartment=="Root" | dataPS@sam_data$Compartment=="Soil", dataPS)->dataPS_RS


#compute RA
dataPS_RS  = transform_sample_counts(dataPS_RS, function(x) x / sum(x) )

prune_samples(dataPS_RS@sam_data$Compartment=="Root", dataPS_RS)->dataPS_R
prune_samples(dataPS_RS@sam_data$Compartment=="Soil", dataPS_RS)->dataPS_S

#melt for permanova
psmelt(dataPS_RS)->data

data %>% 
  filter(., Compartment == "Root") %>% 
  filter(., R401 != "HK") %>% 
  filter(., OTU != "Root401")->stats
data %>% 
  filter(., Compartment == "Root") %>% 
  filter(., R401 != "HK") %>% 
  filter(., OTU != "Root401") %>% 
  group_by(OTU, R401) %>%
  summarise(avgRA=median(Abundance))->data
filter(data, R401=="WT")->dataWT
dataWT$R401<-NULL
colnames(dataWT)[2]<-c("avgRAWT")
left_join(data, dataWT)->data
((data$avgRA/data$avgRAWT)*100)->data$normRA
data.frame(1:17, c("Root918", "Root930", "Soil728", "Soil736", "Root4", "Root265", "Root935", "Root105", "Root241", "Root483D1", "Root483D2", "Root131", "Root335", "Rs_GMI1000", "Root29", "Root627", "Root667"))->order
colnames(order)<-c("order", "OTU")
inner_join(data, order)->data
data[is.na(data)] <- 100
replace(data$normRA, data$normRA<100,100)->data$normRA
replace(data$normRA, data$normRA==Inf,100)->data$normRA
data$normRA-100->data$normRA


ggballoonplot(data=data, y="R401", x="order", fill='grey', size="normRA")+
  ylab("")+
  xlab("")+
  #ggtitle("Halo width of R401(-mutant) on different lawn bacteria")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.border = element_rect(colour = 'black', size=1, fill = NA),
        axis.text=element_text(size=12, color= "black"),
        axis.text.x = element_text(size=12, color= "black", angle = 90, hjust=1),
        axis.title=element_text(size=12), 
        strip.background = element_rect(color="black", fill=NA, size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black"))+
  scale_x_discrete(limits=c(1:17),labels=c("Root918", "Root930", "Soil728", "Soil736", "Root4", "Root265", "Root935", "Root105", "Root241", "Root483D1", "Root483D2", "Root131", "Root335", "Rs_GMI1000", "Root29", "Root627", "Root667"))+
  #scale_x_discrete(limits=c("56", "DD", "56DD"), labels=c("dNRPS", "dDAPG",  "dNRPSdDAPG"))+
  scale_y_discrete(limits=c("54","56", "DD", "54DD","56DD"), labels=c("dAcT","dNRPS", "dDAPG", "dAcTdDAPG", "dNRPSdDAPG"))+
  #scale_y_continuous(limits=c(0, 200))+
  scale_fill_viridis_c(option = "B")+
  NULL
#ggsave("Root_C_balloon_full.pdf", width=11, height=20, unit="cm")


stats %>% 
  group_by(OTU) %>% 
  dunn_test(Abundance~R401, p.adjust.method ="BH") %>% 
  filter(., group2=="WT") %>% 
  filter(., p.adj<=0.05)->test


