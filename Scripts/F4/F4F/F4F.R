
## load libraries
library(tidyverse)

##clean environment
rm(list = ls())

colors<-c("R401" = "#1F1F1F", "dpvdl"="#577590")




##load data file
read_tsv("F4F_data.txt")->data
read_tsv("F4F_setup.txt")->setup
read_tsv("F4F_ODs.txt")->ODs

str_c(setup$Row, setup$Column)->setup$Well
left_join(data, setup)->data
left_join(data, ODs)->data

filter(data, Strain=="empty")->empty

empty %>% 
  group_by(Time) %>% 
  summarise(empty=mean(Mean))->empty


data %>% 
  group_by(Time, Strain, OD, Bio, Dilution) %>% 
  summarise(mean=mean(Mean))->data


left_join(data, empty)->data


#normalize 
#	2 * (%siderophore_activity)/((input_volume)/(OD600))	== Fe_binding capacity (nmol/ml/OD)

# %siderophore_activity = (value- blank) / blank 
#input volume (2)
(2*((data$empty-data$mean)/data$empty))/((0.025*data$Dilution)/data$OD)->data$Sid_activity


###over time
data %>% 
  group_by(Time, Strain, Dilution) %>% 
  filter(Strain=="R401"|Strain=="dpvdl") %>% 
  summarise(mean=mean(Sid_activity), sd=sd(Sid_activity))->data_kin


filter(data_kin, Dilution == 1)->data_kin1
filter(data_kin, Dilution == 5)->data_kin5
filter(data_kin, Dilution == 10)->data_kin10

ggplot(data_kin1) + 
  #geom_point(alpha = 0.4)+
  geom_line(aes(x = Time, y = mean, color=Strain))+
  geom_ribbon(aes(x = Time, y = mean,  fill=Strain, ymin = mean - sd/sqrt(3), ymax = mean + sd/sqrt(3)), alpha=0.1) +
  #stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), aes())+
  theme(panel.background = element_rect(F))+
  theme(axis.line = element_line(T))+
  #geom_hline(yintercept = S782_halo_mean)+
  #scale_y_continuous(breaks = c(seq(0,20,5)),minor_breaks=NULL)+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  theme(axis.text.x = element_text(color="Black"))+
  #theme(legend.position = "none")+ 
  #scale_color_manual(values=colors)+
  #scale_fill_manual(values=colors)+
  ylab("Fe binding capacity [nmol/ml/OD]")+
  xlab("Time [min]")
#ggsave("Sid_assay_hamD_overtime.pdf")








data %>% 
  filter(Time == "2400.3") %>% 
  filter(Strain=="R401"|Strain=="dpvdl") %>% 
  filter(Dilution==1)->data2
ggplot(data=data2, aes(x=Strain, y=Sid_activity))+
  #specify the plot
  geom_jitter(aes(color=Strain, fill=Strain), width = 0.2, alpha=0.3)+
  #  scale_shape_manual(values = c(1:4))+
  geom_boxplot(outlier.shape = NA, lwd=0.9, fill=NA, aes(color=Strain))+
  #  facet_wrap(vars(Dilution, Time), scales="free")+
  xlab("")+
  ylab("Fe binding capacity [nmol/ml/OD]")+
  #ggtitle("Halo width of R401(-mutant) on different lawn bacteria")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.border = element_rect(colour = 'black', size=1, fill = NA),
        axis.text=element_text(size=12, color= "black"),
        axis.text.x = element_text( size=12, color= "black"),
        axis.title=element_text(size=12), 
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.position="none",
        strip.background = element_rect(color="black", fill=NA, size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black"))+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  scale_x_discrete(limits=c("R401", "dpvdl"))+
  #scale_y_continuous(limits = c(0, 120))+
  NULL  

#ggsave("Siderophore_assay3_box_new.pdf", width=8, height=6, unit="cm")  

dunn <- kwAllPairsDunnTest(Sid_activity ~ as.factor(Strain), data=data2, p.adjust.method ="BH")
dunn
plot(dunn)

