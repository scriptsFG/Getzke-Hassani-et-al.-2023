library(tidyverse)
library(PMCMRplus)
#######

rm(list = ls())

#colors for plotting
colors<-c("HK"="#C2C2C1", "WT" = "#1F1F1F", "54"="#F94144","56"="#577590","DD"="#F9C74F", "54DD"="#F8961E", "56DD"="#90BE6D")
shapes<-c("I"=21, "II"=22, "III"=24)


#read 16s and ITS1 table
read_tsv(file = "S5C_data.txt")->freshweights
log2(freshweights$Freshweight)->freshweights$log2Freshweight

#calculate mean of HK
filter(freshweights, freshweights$R401=="HK")->HK
subset(HK, HK$Freshweight!=0)->HK
median(HK$log2Freshweight)->median_HK


ggplot(freshweights, aes(x=R401, y=log2Freshweight))+
  geom_jitter(aes(x=R401, y=log2Freshweight, color=R401), width = 0.3, alpha=0.6)+
  #scale_shape_manual(values = c(1:4))+
  geom_boxplot(outlier.shape = NA, lwd=1, fill=NA, aes(color=R401))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.border = element_rect(colour = 'black', size=1, fill = NA),
        axis.text.y =element_text(size=12, color= "black", vjust = 0.5),
        axis.text.x =element_text(size=12, color= "black", angle = 90, hjust = 1, vjust = 0.5),
        axis.title=element_text(size=12), 
        legend.background=element_blank(),
        legend.key=element_blank(),
        plot.background = element_blank(),
        legend.position="none",
        strip.background = element_rect(color="black", fill=NA, size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black"))+
  #scale_shape(solid = F)+
  #scale_size_manual(values=c(5, 5))+
  #scale_shape_manual(values=c(15, 17))+
  #scale_color_manual(values=correlations$color)+
  #scale_alpha_manual(values=c(0.8, 1))+
  scale_color_manual(values=colors)+
  #scale_fill_manual(values=colors)+
  geom_hline(yintercept = median_HK, size=1)+
  #scale_size_manual(values=c(1, 2))+
  scale_x_discrete(limits=c("HK", "WT", "54", "56", "DD", "54DD", "56DD"), labels=c("HK", "WT", "dAcT", "dNRPS", "dDAPG", "dAcTdDAPG", "dNRPSdDAPG"))+
  #scale_y_continuous(limits=c(-2.5,6.5))+
  #geom_vline(xintercept = 4.5, size=0.5)+
  ylab("log2(Shoot Freshweight)")+
  xlab("")+
  #facet_wrap(~x)+
  NULL
#ggsave("R401_freshweights.pdf")


dunnR <- kwAllPairsDunnTest(log2Freshweight ~ as.factor(R401), data=freshweights, p.adjust.method ="BH")
dunnR
