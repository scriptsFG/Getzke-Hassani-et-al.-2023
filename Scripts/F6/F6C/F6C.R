# load libraries
library(tidyverse)
library(PMCMRplus)

##clean environment
rm(list = ls())

colors<-c("HK"="#C2C2C1", "WT" = "#1F1F1F", "dAcT"="#F94144","dNRPS"="#577590","dDAPG"="#F9C74F", "dDAPGdAcT"="#F8961E", "dDAPGdNRPS"="#90BE6D", "cAcT"="#402425", "dAcTdNRPS"="#83677B")



##load data file
read_tsv("F6C_data.txt", col_select = c(1:6))->data
data[is.na(data)] <- 0


pivot_longer(data, cols = 3:6, names_to ="dilution", values_to = "raw_counts")->data
filter(data, data$raw_counts!=9999)->data
filter(data, data$raw_counts!=0)->data


#load colony information
read_tsv("F6C_setup.txt")->setup

#diution factor
data.frame("dilution"=c("Undiluted", "Dilution_1","Dilution_2", "Dilution_3"), "factor"=c(1,10,100, 1000))->dilution_factor

#fuse data frames
left_join(data, dilution_factor)->data
left_join(data, setup)->data


#take average of dilutions
data$raw_counts*data$factor->data$counts

data %>%
  group_by(Sample, Bio, Tech, Treatment, Root_FW) %>%
  summarise(counts=median(counts))->data


data$counts/data$Root_FW->data$counts


##make boxplots
ggplot(data=data, aes(x=Treatment, y=counts))+
  #specify the plot
  geom_jitter(aes(color=Treatment, fill=Treatment), width = 0.2, alpha=0.3)+
  scale_shape_manual(values = c(1:4))+
  geom_boxplot(outlier.shape = NA, lwd=0.9, fill=NA, aes(color=Treatment))+
  xlab("")+
  ylab("CFU/mg root")+
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
  scale_x_discrete(limits=c("HK", "WT", "dAcT","dNRPS","dDAPG",  "dDAPGdAcT", "dDAPGdNRPS"))+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)
#ggsave("colony_counts.pdf", width = 8, height=13, unit="cm")



aov <- aov(counts ~ Treatment, data = data)
summary(aov)
TukeyHSD(aov)->tukey
tukey

