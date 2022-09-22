

## load libraries
library(tidyverse)
library(PMCMRplus)
##clean environment
rm(list = ls())

colors<-c("HK"="#C2C2C1", "WT" = "#1F1F1F", "dAcT"="#F94144","dNRPS"="#577590","dDAPG"="#F9C74F", "dDAPGdAcT"="#F8961E", "dDAPGdNRPS"="#90BE6D", "cAcT"="#402425", "dAcTdNRPS"="#83677B")
colors2<-c("Agaricomycetes"="#bfbfbf", "Ascomycota"="#8B008B", "Dothideomycetes"="#DB7093", "Eurotiomycetes"="#bfbfbf", "Glomeromycetes"="#bfbfbf", "Lecanoromycetes"="#bfbfbf", "Leotiomycetes"="#bfbfbf", "Microbotryomycetes"="#bfbfbf", "Mucoromycotina"="#bfbfbf", "Pezizomycetes"="#bfbfbf", "Pezizomycotina"="#bfbfbf", "Saccharomycetes"="#bfbfbf", "Sordariomycetes"="#8B008B", "Tremellomycetes"="#bfbfbf", "[Saprospirae]"="#bfbfbf", "Acidimicrobiia"="#bfbfbf",  "Actinobacteria"="#B22222",  "Alphaproteobacteria"="#32CD32", "Bacilli"="#FFA500",  "Bacteroidia"="#bfbfbf",  "Betaproteobacteria"="#228B22",  "Chlamydiia"="#bfbfbf",  "Chloroflexi"="#bfbfbf", 
           "Clostridia"="#bfbfbf",  "Cytophagia"="#bfbfbf",  "DA052"="#bfbfbf",  "Deltaproteobacteria"="#bfbfbf", 
           "Fibrobacteria"="#bfbfbf",  "Flavobacteriia"="#4682B4",  "Gammaproteobacteria"="#006400",  "Gemm-1"="#bfbfbf",  "Nitrospira"="#bfbfbf", 
           "Rubrobacteria"="#bfbfbf",  "SC3"="#bfbfbf",  "Sphingobacteriia"="#bfbfbf",  "Thermoleophilia"="#bfbfbf",  "TK10"="#bfbfbf",  "TM7-3"="#bfbfbf",  "VHS-B5-50"="#bfbfbf")



##load data file
read_tsv("F4G_data.txt")->data


#load colony information
read_tsv("F4G_setup.txt")->colonies
colonies[1:7,1:2]->colonies
left_join(data, colonies)->data


##aggregate
aggregate(data$`Colony_diameter`, by=list(data$`Lawn_bacterium`, data$Plate, data$`Mutant`),FUN=mean)->data2
colnames(data2)<-c("lawn_bacterium","plate","R401","colony_size")
aggregate(data$`Halo_diameter`, by=list(data$`Lawn_bacterium`, data$Plate, data$`Mutant`),FUN=mean)->data3
colnames(data3)<-c("lawn_bacterium", "plate", "R401","halo_diameter")
left_join(data2, data3)->data

#compute halo_width
(data$halo_diameter-data$colony_size)/2->data$halo_width




#fuse lawn bacteria tax with data frame
#load lawn information
read_tsv("tax_file.txt")->tax
colnames(tax)[2]<-c("lawn_bacterium")
left_join(data, tax, by="lawn_bacterium")->data
data[is.na(data)] <- 0




#compute average susceptibility towards R401
filter(data, data$R401=="WT")->R401
R401 %>%
  group_by(lawn_bacterium) %>%
  summarise(avghalowidth=mean(halo_width))->avgR401
rep(1, length(avgR401$lawn_bacterium))->avgR401$position

#compute remaining antimicrobial activity in dDAPGdAcT for every strain
filter(data, data$R401=="dDAPGdAcT")->res
res %>%
  group_by(lawn_bacterium) %>%
  summarise(res_halowidth=mean(halo_width))->res
left_join(res, tax)->res

#compute lack of halo production compared to WT
colnames(R401)[6]<-c("WT_halo_width")
R401 %>% select(lawn_bacterium, plate, WT_halo_width)->R401
left_join(data, R401)->R401
replace(R401$halo_width, R401$halo_width<=0,0.00001)->R401$halo_width
1/R401$halo_width->R401$halo_width2
1/R401$WT_halo_width->R401$WT_halo_width2
(R401$halo_width2/R401$WT_halo_width2)*100->R401$norm_halo_width2
(R401$halo_width/R401$WT_halo_width)*100->R401$norm_halo_width
R401 %>%
  group_by(lawn_bacterium, R401) %>%
  summarise(avgnormhalowidth=mean(norm_halo_width), avgnormhalowidth2=mean(norm_halo_width2))->avgnormR401
write.table(avgnormR401, "halo_for_F5CF.txt",sep="\t",row.names=FALSE)

#compute average for each mutant
avgnormR401 %>%
  group_by(R401) %>%
  summarise(avgavgnormhalowidth=mean(avgnormhalowidth))->avgavgnormR401
rep(1, length(avgavgnormR401$R401))->avgavgnormR401$position




ggplot(data = avgnormR401, aes(x = lawn_bacterium, y = R401)) +
  geom_tile(aes(fill = avgnormhalowidth))+
  #scale_fill_distiller(palette = "YlGnBu") +
  #scale_x_discrete(limits=c("S782", "S736", "S728", "R85", "R930", "R918", "R22", "R553", "R4", "S805", "R265", "R135", "R935", "R901", "R667", "R559", "R627", "R179", "Rs", "R29", "R402", "R418", "R189", "R239", "R131", "R241", "R436", "R105", "R50", "R1334", "R483D2", "R670", "R483D1"))+
  scale_x_discrete(limits=c("R22", "R930", "S728", "R85", "S805", "S736", "S782", "R4", "R553", "R265", "R29", "R402", "Rs", "R179", "R667",  "R483D1", "R670", "R483D2", "R1334", "R50", "R105"))+
  scale_y_discrete(limits=c("dDAPGdAcT","dDAPGdNRPS","dDAPG","dAcT","dNRPS", "cAcT","WT"))+
  labs(title = "", y = "Root401 genotype", x= "Target bacteria")+
  #scale_fill_discrete(low = "blue", high = "red")+
  #scale_color_gradient2('avgnormhalowidth', limits=c(-10, 5), low = 'blue', mid = 'white', high = 'red')
  scale_fill_gradient2('avgnormhalowidth',  limits=c(0, 100), midpoint=c(1), mid='white', high = "#1d2d44", na.value="black")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = 'grey'), axis.text.x = element_text(angle = 90, hjust = 1))
#ggsave("heatmap.pdf")

ggplot(data = avgR401, aes(x = lawn_bacterium, y = position)) +
  geom_tile(aes(fill = avghalowidth))+
  #scale_fill_distiller(palette = "YlGnBu") +
  #scale_x_discrete(limits=c("S782", "S736", "S728", "R85", "R930", "R918", "R22", "R553", "R4", "S805", "R265", "R135", "R935", "R901", "R667", "R559", "R627", "R179", "Rs", "R29", "R402", "R418", "R189", "R239", "R131", "R241", "R436", "R105", "R50", "R1334", "R483D2", "R670", "R483D1"))+
  scale_x_discrete(limits=c("R22", "R930", "S728", "R85", "S805", "S736", "S782", "R4", "R553", "R265", "R29", "R402", "Rs", "R179", "R667",  "R483D1", "R670", "R483D2", "R1334", "R50", "R105"))+
  #scale_y_discrete(limits=c("WT", "cAcT","dNRPS","dAcT","dDAPG", "dDAPGdNRPS", "dDAPGdAcT"))+
  labs(title = "", y = "Target microbe", x= "Effector microbe")+
  #scale_fill_discrete(low = "blue", high = "red")+
  #scale_color_gradient2('avgnormhalowidth', limits=c(-10, 5), low = 'blue', mid = 'white', high = 'red')
  scale_fill_gradient2('avgnormhalowidth', low = "white",  high = "black")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))
#ggsave("Target_panel.pdf")

ggplot(data = avgavgnormR401, aes(x = position, y = R401)) +
  geom_tile(aes(fill = avgavgnormhalowidth))+
  #scale_fill_distiller(palette = "YlGnBu") +
  #scale_x_discrete(limits=c("S782", "S736", "S728", "R85", "R930", "R918", "R22", "R553", "R4", "S805", "R265", "R135", "R935", "R901", "R667", "R559", "R627", "R179", "Rs", "R29", "R402", "R418", "R189", "R239", "R131", "R241", "R436", "R105", "R50", "R1334", "R483D2", "R670", "R483D1"))+
  #scale_x_discrete(limits=c("S782", "S736", "S728", "R85", "R930", "R22", "R553", "R4", "S805", "R265", "R667", "R179", "Rs", "R29", "R402",  "R105", "R50", "R1334", "R483D2", "R670", "R483D1"))+
  scale_y_discrete(limits=c("dDAPGdAcT","dDAPGdNRPS","dDAPG","dAcT","dNRPS", "cAcT","WT"))+
  labs(title = "", y = "Target microbe", x= "Effector microbe")+
  #scale_fill_discrete(low = "blue", high = "red")+
  #scale_color_gradient2('avgnormhalowidth', limits=c(-10, 5), low = 'blue', mid = 'white', high = 'red')
  scale_fill_gradient2('avgnormhalowidth',  limits=c(0, 100), midpoint=c(1), mid='white', high = "#1d2d44", na.value="black")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))
#ggsave("R401_panel.pdf")

###color code for lawn bacteria
rep(1, length(tax$lawn_bacterium))->tax$position
ggplot(data = tax, aes(x = lawn_bacterium, y = position)) +
  geom_tile(aes(fill = Class))+
  #scale_fill_distiller(palette = "YlGnBu") +
  #scale_x_discrete(limits=c("S782", "S736", "S728", "R85", "R930", "R918", "R22", "R553", "R4", "S805", "R265", "R135", "R935", "R901", "R667", "R559", "R627", "R179", "Rs", "R29", "R402", "R418", "R189", "R239", "R131", "R241", "R436", "R105", "R50", "R1334", "R483D2", "R670", "R483D1"))+
  scale_x_discrete(limits=c("R22", "R930", "S728", "R85", "S805", "S736", "S782", "R4", "R553", "R265", "R29", "R402", "Rs", "R179", "R667",  "R483D1", "R670", "R483D2", "R1334", "R50", "R105"))+
  #scale_y_discrete(limits=c("WT", "cAcT","dNRPS","dAcT","dDAPG", "dDAPGdNRPS", "dDAPGdAcT"))+
  labs(title = "", y = "Target microbe", x= "Effector microbe")+
  #scale_fill_discrete(low = "blue", high = "red")+
  #scale_color_gradient2('avgnormhalowidth', limits=c(-10, 5), low = 'blue', mid = 'white', high = 'red')
  #scale_fill_gradient2('avgnormhalowidth', low = "white",  high = "black")+
  scale_fill_manual(values=colors2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))
#ggsave("Tax_panel.pdf")




