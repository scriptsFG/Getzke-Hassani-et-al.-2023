library(tidyverse)
library(PMCMRplus)
#######

rm(list = ls())

#####
# 54 = ∆pvdy
# 56 = ∆pvdl
# DD = ∆phld

#colors for plotting
colors<-c("HK"="#C2C2C1", "WT" = "#1F1F1F", "54"="#F94144","56"="#577590","DD"="#F9C74F", "54DD"="#F8961E", "56DD"="#90BE6D")




###halo assay data
read_tsv("halo_for_F5CF.txt")->halo
data.frame(c("WT","54", "56","DD", "54DD", "56DD"), c("WT","dAcT", "dNRPS","dDAPG", "dDAPGdAcT", "dDAPGdNRPS"))->labels
colnames(labels)<-c("New", "R401")
left_join(halo, labels)->halo
str_replace(string=halo$lawn_bacterium, pattern = "R",  replacement = "Root")->halo$lawn_bacterium
str_replace(string=halo$lawn_bacterium, pattern = "S",  replacement = "Soil")->halo$lawn_bacterium
str_replace(string=halo$lawn_bacterium, pattern = "Roots GMI1000",  replacement = "Rs_GMI1000")->halo$lawn_bacterium
filter(halo, halo$R401!="cAcT")->halo
halo$R401<-NULL
colnames(halo)[1]<-c("OTU")
colnames(halo)[4]<-c("R401")
(((1/(halo$avgnormhalowidth2))*10000)-100)*(-1)->halo$avgnormhalowidth2
halo %>%
  group_by(R401) %>%
  summarise(avgnormhalowidth2=mean(avgnormhalowidth2))->halo



#read pre calculated %variance per mutant and compartment
read_tsv("perc_var_explained.txt", col_select = c(1,2,4))->var_exp

left_join(halo, var_exp)->halo
halo$Var_exp*100->halo$Var_exp




ggplot(halo, aes(x=avgnormhalowidth2, y=Var_exp))+
  geom_point(aes(size=0.1, color=R401, fill=R401))+
  geom_smooth(method='lm', se=T, aes(color=Compartment))+
  #stat_summary(fun.data=)
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.border = element_rect(colour = 'black', size=1, fill = NA),
        axis.text=element_text(size=12, color= "black"),
        axis.title=element_text(size=12), 
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
  #scale_x_discrete(limits=factor(c(1,2)), labels=c("Soil", "Root"))+
  #scale_x_continuous(limits=c(0, 1), breaks = c(0,0.2, 0.4, 0.6, 0.8, 1))+
  #geom_vline(xintercept = 0, size=0.5)+
  ylab("Explained variance [%]")+
  xlab("Average relative decrease in halo width [%]")+
  ggtitle("p=0.0103, R2=0.7984 ")
#facet_wrap(~OTU)+
NULL
#ggsave("Correlation_halo_Var_exp.pdf")


filter(halo, Compartment=="Root")->haloR
filter(halo, Compartment=="Soil")->haloS

cor.test(x=haloR$avgnormhalowidth2, y=haloR$Var_exp, method=c("spearman"), exact = F)->test
test$p.value
cor.test(x=haloS$avgnormhalowidth2, y=haloS$Var_exp, method=c("spearman"), exact = F)->test
test$p.value

#regression analysis
summary(lm(avgnormhalowidth2 ~ Var_exp, data=haloR)) 
summary(lm(avgnormhalowidth2 ~ Var_exp, data=haloS)) 





