library(vegan)
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



###halo assay data
read_tsv("halo_for_F5CF.txt")->halo
data.frame(c("WT","54", "56","DD", "54DD", "56DD"), c("WT","dAcT", "dNRPS","dDAPG", "dDAPGdAcT", "dDAPGdNRPS"))->labels
colnames(labels)<-c("New", "R401")
left_join(halo, labels)->halo
str_replace(string=halo$lawn_bacterium, pattern = "R",  replacement = "Root")->halo$lawn_bacterium
str_replace(string=halo$lawn_bacterium, pattern = "S",  replacement = "Soil")->halo$lawn_bacterium
str_replace(string=halo$lawn_bacterium, pattern = "Roots GMI1000",  replacement = "Rs_GMI1000")->halo$lawn_bacterium
filter(halo, halo$R401!="cAcT")->halo
halo$R401<-NULL
colnames(halo)[1]<-c("OTU")
colnames(halo)[4]<-c("R401")
(((1/(halo$avgnormhalowidth2))*10000)-100)*(-1)->halo$avgnormhalowidth2
halo %>%
  group_by(R401) %>%
  summarise(avgnormhalowidth2=mean(avgnormhalowidth2))->halo



##### Shannon indeces

dataPS_RS <- prune_taxa(taxa_sums(dataPS_RS) > 0, dataPS_RS)
plot_richness(dataPS_RS, x="R401", measures=c("Shannon"))->Shannon
#plot_richness(dataPS_S, x="R401", measures=c("Shannon"))->Shannon
#plot_richness(dataPS_R, x="R401")->Shannon
Shannon[["data"]]->Shannon



ggplot(data=Shannon, aes(x=R401, y=value, color=R401))+
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


Shannon %>% 
  group_by(R401, Compartment) %>% 
  summarise(Shannon=mean(value))->Shannon



inner_join(Shannon, halo)->halo
filter(halo, Compartment=="Root")->haloR
filter(halo, Compartment=="Soil")->haloS

ggplot(halo, aes(x=avgnormhalowidth2, y=Shannon))+
  geom_point(aes(size=0.1, color=R401, fill=R401))+
  geom_smooth(method='lm', se=T, aes(color=Compartment))+
  #stat_summary(fun.data=)
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.border = element_rect(colour = 'black', size=1, fill = NA),
        axis.text=element_text(size=12, color= "black"),
        axis.title=element_text(size=12), 
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.position="none")+
  #scale_shape(solid = F)+
  #scale_size_manual(values=c(5, 5))+
  #scale_shape_manual(values=c(15, 17))+
  #scale_color_manual(values=correlations$color)+
  #scale_alpha_manual(values=c(0.8, 1))+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  #scale_size_manual(values=c(1, 2))+
  #scale_x_discrete(limits=factor(c(1,2)), labels=c("Soil", "Root"))+
  #scale_x_continuous(limits=c(0, 1), breaks = c(0,0.2, 0.4, 0.6, 0.8, 1))+
  #geom_vline(xintercept = 0, size=0.5)+
  ylab("Average Shannon diversity")+
  xlab("Average relative decrease in halo width [%]")+
#facet_wrap(~OTU)+
NULL
#ggsave("Correlation_halo_Shannon.pdf")


cor.test(x=haloR$avgnormhalowidth2, y=haloR$Shannon, method=c("spearman"), exact = F)->test
test$p.value
cor.test(x=haloS$avgnormhalowidth2, y=haloS$Shannon, method=c("spearman"), exact = F)->test
test$p.value

#regression analysis
summary(lm(avgnormhalowidth2 ~ Shannon, data=haloR)) 
summary(lm(avgnormhalowidth2 ~ Shannon, data=haloS)) 
