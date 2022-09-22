################Rs halo assay######################


## load libraries
library(tidyverse)
library(PMCMRplus)
##clean environment
rm(list = ls())

colors<-c("WT" = "#1F1F1F", "dphld"="#F9C74F")


##load data file
read_tsv("F3E_data.txt", col_select = c(1:6))->data


#load colony information
read_tsv("F3E_setup.txt")->colonies

inner_join(data, colonies)->data




##aggregate
data %>% 
  group_by(Lawn_bacterium, Plate, Mutant) %>%
  summarise(halo=mean(`Halo_diameter cm`), colony=mean(`Colony_diameter cm`))->data
((data$halo)-data$colony)*10->data$size_mm
filter(data, Lawn_bacterium=="Rs4")->data

data %>% 
  group_by(Mutant) %>% 
  summarise(ave=mean(size_mm))->test
##make boxplots
ggplot(data=data, aes(x=Mutant, y=size_mm))+
  #specify the plot
  geom_jitter(aes(color=Mutant, fill=Mutant), width = 0.2, alpha=0.3)+
  scale_shape_manual(values = c(1:4))+
  geom_boxplot(outlier.shape = NA, lwd=0.9, fill=NA, aes(color=Mutant))+
  facet_wrap(~Lawn_bacterium)+
  xlab("")+
  ylab("Halo size [mm]")+
  #ggtitle("Halo width of R401(-mutant) on different lawn bacteria")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.border = element_rect(colour = 'black', size=1, fill = NA),
        axis.text=element_text(size=14, color= "black"),
        axis.text.x = element_text(size=14, color= "black", angle = 90, hjust=1),
        axis.title=element_text(size=14), 
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.position="none",
        strip.background = element_rect(color="black", fill=NA, size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black"))+
  scale_x_discrete(limits=c("WT", "dphld"))+
  #geom_vline(xintercept = 7.5)+
  scale_y_continuous(limits=c(0,6))+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  #ggtitle("Rs")+
  NULL
#ggsave("Halo_pvd_Rs.pdf", height=9, width=5, unit="cm")


dunn <- kwAllPairsDunnTest(size_mm ~ as.factor(Mutant), data=data, p.adjust.method ="BH")
dunn
plot(dunn)

