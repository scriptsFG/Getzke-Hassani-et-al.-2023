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

#####ordination
dataPS_R_ord <- ordinate(dataPS_R, "MDS", "bray")
as.data.frame(dataPS_R_ord[[4]][,1:2])->ordR
rownames_to_column(ordR, var = "Sample_Name")->ordR
left_join(ordR, map)->ordR


ggplot(ordR, aes(x = Axis.1, y = Axis.2))+
  #geom_polygon(aes(color=R401, alpha=0.5, fill=R401))+
  geom_point(size=3,aes(color=R401, fill=R401, alpha=0.95))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.border = element_rect(colour = 'black', size=0.5, fill = NA),
        axis.text=element_text(size=10, color= "black"),
        axis.text.x =element_text(size=10, color= "black", angle = 90, hjust = 1),
        axis.title=element_text(size=10), 
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.position="none",
        strip.background = element_rect(color="black", fill=NA, size=0.5, linetype="solid"),
        strip.text.x = element_text(size = 8, color = "black"))+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  #scale_shape_manual(values=shapes)+
  #stat_ellipse(aes(color=R401), level = 0.95) +
  NULL
#ggsave("PCoA_RA_Root.pdf", width=12, height=12, unit="cm")


#####ordination
dataPS_S_ord <- ordinate(dataPS_S, "MDS", "bray")
as.data.frame(dataPS_S_ord[[4]][,1:2])->ordS
rownames_to_column(ordS, var = "Sample_Name")->ordS
left_join(ordS, map)->ordS


ggplot(ordS, aes(x = Axis.1, y = Axis.2))+
  #geom_polygon(aes(color=R401, alpha=0.5, fill=R401))+
  geom_point(size=3,aes(color=R401, fill=R401, alpha=0.95))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.border = element_rect(colour = 'black', size=0.5, fill = NA),
        axis.text=element_text(size=10, color= "black"),
        axis.text.x =element_text(size=10, color= "black", angle = 90, hjust = 1),
        axis.title=element_text(size=10), 
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.position="none",
        strip.background = element_rect(color="black", fill=NA, size=0.5, linetype="solid"),
        strip.text.x = element_text(size = 8, color = "black"))+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  #scale_shape_manual(values=shapes)+
  #stat_ellipse(aes(color=R401), level = 0.95) +
  NULL
#ggsave("PCoA_RA_Soil.pdf", width=12, height=12, unit="cm")





#permanova 


#######permanova on RA


#Root
# remove "#" in order to test different comparisons
data %>% 
  #filter(., R401!= "DD") %>% 
  #filter(., R401!= "54DD") %>%
  #filter(., R401!= "56DD") %>% 
  #filter(., R401!= "54") %>% 
  #filter(., R401!= "56") %>% 
  #filter(., R401!= "HK") %>% 
  filter(., Compartment== "Root")->perm
perm->mapAdonis
#generate otu table
#select only needed cols
perm %>% select(Sample, OTU, Abundance)->perm


#spread into OTU table
spread(perm, key=Sample, value=Abundance)->perm
column_to_rownames(perm, var="OTU")->perm
perm %>% arrange(row.names(perm))->perm
#generate mapping file
mapAdonis %>% select(Sample, R401)->mapAdonis
distinct(mapAdonis)->mapAdonis
row.names(mapAdonis)<-NULL
mapAdonis->mapAdonis2
rownames(mapAdonis)<-mapAdonis$Sample
mapAdonis$Sample<-NULL
mapAdonis %>% arrange(row.names(mapAdonis))->mapAdonis


#perform permanova analysis
dist_perm=vegdist(t(perm), method="bray")
mod=adonis(dist_perm~R401, data=mapAdonis)
res_perm<-matrix(ncol=4,nrow=5)
for (e in 1:5){
  
  res_perm[e,1]<-mod$aov.tab$Df[e]
  res_perm[e,2]<-mod$aov.tab$F.Model[e]
  res_perm[e,3]<-mod$aov.tab$`Pr(>F)`[e]
  res_perm[e,4]<-mod$aov.tab$R2[e]
}
colnames(res_perm)<-c("df","F","P","R2")
rownames(res_perm)<-c("R401","Residuals","Total", "", "")
res_perm


#Soil
# remove "#" in order to test different comparisons
data %>% 
  #filter(., R401!= "DD") %>%
  #filter(., R401!= "54DD") %>%
  #filter(., R401!= "56DD") %>%
  #filter(., R401!= "54") %>%
  #filter(., R401!= "56") %>%
  #filter(., R401!= "HK") %>% 
  filter(., Compartment== "Soil")->perm
perm->mapAdonis
#generate otu table
#select only needed cols
perm %>% select(Sample, OTU, Abundance)->perm


#spread into OTU table
spread(perm, key=Sample, value=Abundance)->perm
column_to_rownames(perm, var="OTU")->perm
perm %>% arrange(row.names(perm))->perm
#generate mapping file
mapAdonis %>% select(Sample, R401)->mapAdonis
distinct(mapAdonis)->mapAdonis
row.names(mapAdonis)<-NULL
mapAdonis->mapAdonis2
rownames(mapAdonis)<-mapAdonis$Sample
mapAdonis$Sample<-NULL
mapAdonis %>% arrange(row.names(mapAdonis))->mapAdonis


#perform permanova analysis
dist_perm=vegdist(t(perm), method="bray")
mod=adonis(dist_perm~R401, data=mapAdonis)
res_perm<-matrix(ncol=4,nrow=5)
for (e in 1:5){
  
  res_perm[e,1]<-mod$aov.tab$Df[e]
  res_perm[e,2]<-mod$aov.tab$F.Model[e]
  res_perm[e,3]<-mod$aov.tab$`Pr(>F)`[e]
  res_perm[e,4]<-mod$aov.tab$R2[e]
}
colnames(res_perm)<-c("df","F","P","R2")
rownames(res_perm)<-c("R401","Residuals","Total", "", "")
res_perm

