library(tidyverse)
#######

rm(list = ls())

#colors for plotting
colors<-c("Agaricomycetes"="#bfbfbf", "Ascomycota"="#8B008B", "Dothideomycetes"="#DB7093", "Eurotiomycetes"="#bfbfbf", "Glomeromycetes"="#bfbfbf", "Lecanoromycetes"="#bfbfbf", "Leotiomycetes"="#bfbfbf", "Microbotryomycetes"="#bfbfbf", "Mucoromycotina"="#bfbfbf", "Pezizomycetes"="#bfbfbf", "Pezizomycotina"="#bfbfbf", "Saccharomycetes"="#bfbfbf", "Sordariomycetes"="#8B008B", "Tremellomycetes"="#bfbfbf", "[Saprospirae]"="#bfbfbf", "Acidimicrobiia"="#bfbfbf",  "Actinobacteria"="#B22222",  "Alphaproteobacteria"="#32CD32", "Bacilli"="#FFA500",  "Bacteroidia"="#bfbfbf",  "Betaproteobacteria"="#228B22",  "Chlamydiia"="#bfbfbf",  "Chloroflexi"="#bfbfbf", 
           "Clostridia"="#bfbfbf",  "Cytophagia"="#bfbfbf",  "DA052"="#bfbfbf",  "Deltaproteobacteria"="#bfbfbf", 
           "Fibrobacteria"="#bfbfbf",  "Flavobacteriia"="#4682B4",  "Gammaproteobacteria"="#006400",  "Gemm-1"="#bfbfbf",  "Nitrospira"="#bfbfbf", 
           "Rubrobacteria"="#bfbfbf",  "SC3"="#bfbfbf",  "Sphingobacteriia"="#bfbfbf",  "Thermoleophilia"="#bfbfbf",  "TK10"="#bfbfbf",  "TM7-3"="#bfbfbf",  "VHS-B5-50"="#bfbfbf")



#mean_ant vs n_ant colored by class
read_tsv(file = "F3A_data.txt")->meta
meta$mean_ant*100->meta$mean_ant
meta$mean_sus*100->meta$mean_sus
#plot correlations of ant and pim score
ggplot(meta, aes( x=mean_ant, y=n_ant, color=Class))+
  #geom_density()+
  geom_jitter(size=3, alpha=1)+
  #stat_summary(fun.data=)
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.border = element_rect(colour = 'black', size=1, fill = NA),
        axis.text=element_text(size=14, color= "black"),
        axis.title=element_text(size=14), 
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
  #scale_y_continuous(limits=c(0, 0.8))+
  #geom_vline(xintercept = 0, size=0.5)+
  xlab("Average antagonistic activity [mm]")+
  ylab("Number of observed antagonstic interactions")+
  #geom_smooth(method='lm', color='black', se=F)+
  NULL
#ggsave("ant_ant.pdf", width =13, height=13, unit="cm")



