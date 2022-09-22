library(tidyverse)
#######


rm(list = ls())

#colors for plotting
colors<-c("Rs"="#006400", "L182"="#FFA500", "L394"="#4682B4",  "Agaricomycetes"="#bfbfbf", "Ascomycota"="#8B008B", "Dothideomycetes"="#DB7093", "Eurotiomycetes"="#bfbfbf", "Glomeromycetes"="#bfbfbf", "Lecanoromycetes"="#bfbfbf", "Leotiomycetes"="#bfbfbf", "Microbotryomycetes"="#bfbfbf", "Mucoromycotina"="#bfbfbf", "Pezizomycetes"="#bfbfbf", "Pezizomycotina"="#bfbfbf", "Saccharomycetes"="#bfbfbf", "Sordariomycetes"="#8B008B", "Tremellomycetes"="#bfbfbf", "[Saprospirae]"="#bfbfbf", "Acidimicrobiia"="#bfbfbf",  "Actinobacteria"="#B22222",  "Alphaproteobacteria"="#32CD32", "Bacilli"="#FFA500",  "Bacteroidia"="#bfbfbf",  "Betaproteobacteria"="#228B22",  "Chlamydiia"="#bfbfbf",  "Chloroflexi"="#bfbfbf", 
           "Clostridia"="#bfbfbf",  "Cytophagia"="#bfbfbf",  "DA052"="#bfbfbf",  "Deltaproteobacteria"="#bfbfbf", 
           "Fibrobacteria"="#bfbfbf",  "Flavobacteriia"="#4682B4",  "Gammaproteobacteria"="#006400",  "Gemm-1"="#bfbfbf",  "Nitrospira"="#bfbfbf", 
           "Rubrobacteria"="#bfbfbf",  "SC3"="#bfbfbf",  "Sphingobacteriia"="#bfbfbf",  "Thermoleophilia"="#bfbfbf",  "TK10"="#bfbfbf",  "TM7-3"="#bfbfbf",  "VHS-B5-50"="#bfbfbf")


#read tables
read_tsv(file = "F3B_data.txt")->all
data.frame("Class"="Rs", "halo"=5.840)->Rs




all %>% 
  filter(R401>0) %>% 
  group_by(Class) %>% 
  summarise(halo=mean(R401), sum(R401>0))->all
#bind_rows(all, leaf)->all
bind_rows(all, data.frame("Class"="yaxis","halo"=0))->all
all$halo*10->all$halo

bind_rows(all, Rs)->all

ggplot(all, aes(Class, halo, fill = Class)) + 
  geom_col(width = 0.8, color="#000000") + 
  geom_hline(yintercept = 0, color = "#000000", size = 0.5) +
  #geom_hline(yintercept = 0.1, color = "#000000", size = 0.1) +
  geom_hline(yintercept = 2, color = "#000000", size = 0.1) +
  #geom_hline(yintercept = 0.3, color = "#000000", size = 0.1) +
  geom_hline(yintercept = 4, color = "#000000", size = 0.1) +
  #geom_hline(yintercept = 0.5, color = "#000000", size = 0.1) +
  geom_hline(yintercept = 6, color = "#000000", size = 0.1) +
  #geom_hline(yintercept = 0.7, color = "#000000", size = 0.1) +
  geom_hline(yintercept = 8, color = "#000000", size = 0.1) +
  geom_vline(xintercept = 0.5, color = "#000000", size = 0.5)+
  scale_y_continuous(limits = c(-2, 8)) +
  #scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  scale_x_discrete(limits=c("yaxis", "Actinobacteria", "Bacilli",  "Flavobacteriia",  "Alphaproteobacteria", "Betaproteobacteria","Gammaproteobacteria","yaxis", "Rs"))+
  #geom_text(label = "Gross\nValue Added\nper Hour",
  #          x = 0, y = 0, color = "white", vjust = 1.7, size = 5,
  #          check_overlap = TRUE, fontface = "bold") +
  coord_polar()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.border = element_rect(colour = 'black', size=0.5, fill = NA),
        axis.text=element_text(size=10, color= "black"),
        axis.text.x =element_text(size=NA, color= NA),
        axis.title=element_text(size=10), 
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.position="none",
        strip.background = element_rect(color="black", fill=NA, size=0.5, linetype="solid"),
        strip.text.x = element_text(size = 8, color = "black"))+
  NULL
  

#ggsave("R401_targets_new.pdf")




