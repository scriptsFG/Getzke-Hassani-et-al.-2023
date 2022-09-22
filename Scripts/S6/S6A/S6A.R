## load libraries
library(tidyverse)
library(PMCMRplus)

##clean environment
rm(list = ls())

colors<-c("HK"="#C2C2C1", "WT" = "#1F1F1F", "dAcT"="#F94144","dNRPS"="#577590","dDAPG"="#F9C74F", "dDAPGdAcT"="#F8961E", "dDAPGdNRPS"="#90BE6D", "cAcT"="#402425", "dAcTdNRPS"="#83677B")


##load data file
read_tsv("S6A_data.txt")->data



##make boxplots
ggplot(data=data, aes(x=R401, y=counts))+
  #specify the plot
  geom_jitter(aes(color=R401, fill=R401), width = 0.2, alpha=0.3)+
  #scale_shape_manual(values = c(1:4))+
  geom_boxplot(outlier.shape = NA, lwd=0.9, fill=NA, aes(color=R401))+
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
  scale_x_discrete(limits=c("WT", "dDAPGdAcT"), labels=c("WT", "dpvdydphld"))+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  NULL
#ggsave("colony_counts_FP.pdf", width=3, height=6, unit="cm")

dunn <- kwAllPairsDunnTest(counts ~ as.factor(R401), data=data, p.adjust.method ="BH")
dunn


aov <- aov(counts ~ R401, data = data)
summary(aov)
TukeyHSD(aov)->tukey
tukey

