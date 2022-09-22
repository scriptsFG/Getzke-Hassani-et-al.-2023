

## load libraries
library(tidyverse)
library(PMCMRplus)
##clean environment
rm(list = ls())

colors<-c("Fe" = "#1F1F1F", "WT" = "#1F1F1F", "dpvdy"="#F94144","dpvdl"="#577590","dpvdyPVDY"="#402425", "dpvdydpvdl"="#83677B", "tn5pvdy"="#F94144", "dpvds"="#8CA5BA", "dphld"="#F9C74F", "dpvdydphld"="#F8961E", "dpvdldphld"="#90BE6D")


##load data file
read_tsv("S4D_data.txt", col_select = c(1:8))->data


#load colony information
read_tsv("S4D_setup.txt")->colonies

inner_join(data, colonies)->data




##aggregate
data %>% 
  group_by(Lawn_bacterium, Iron, Plate, Mutant, Producer, Mutant_setup) %>%
  summarise(halo=mean(Halo_size_mm))->data
filter(data, Iron=="Fe")->Fe
filter(Fe, Mutant=="WT")->Fe
Fe$Iron->Fe$Mutant
filter(data, Iron!="Fe")->data
bind_rows(data, Fe)->data
filter(data, Mutant_setup=="III")->data3









ggplot(data=data3, aes(x=Mutant, y=halo))+
  #specify the plot
  geom_jitter(aes(color=Mutant, fill=Mutant), width = 0.2, alpha=0.3)+
  scale_shape_manual(values = c(1:4))+
  geom_boxplot(outlier.shape = NA, lwd=0.9, fill=NA, aes(color=Mutant))+
  #facet_wrap(Iron~Lawn_bacterium)+
  xlab("")+
  ylab("Halo Width (mm)")+
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
  scale_x_discrete(limits=c("WT", "dpvdy","dpvdl","Fe"))+
  #scale_x_discrete(limits=c("WT", "dpvdy","dpvdl","Fe"), labels=c())+
  geom_vline(xintercept = 3.5)+
  scale_y_continuous(limits = c(0,7))+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  ggtitle("")
#ggsave("Halo_Rs_setup3.pdf", width=12, height=14, unit="cm")


dunn <- kwAllPairsDunnTest(halo ~ as.factor(Mutant), data=data3, p.adjust.method ="BH")
dunn
plot(dunn)



