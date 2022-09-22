library(tidyverse)
#######

rm(list = ls())

#colors for plotting
colors<-c("WT" = "#1F1F1F", "dAcT"="#F94144","dNRPS"="#577590","dDAPG"="#F9C74F", "dDAPGdAcT"="#F8961E", "dDAPGdNRPS"="#90BE6D")


#read data file for ODS
read_tsv(file = "S6B_data.txt")->are
read_tsv(file = "S6C_data.txt")->fe



#read setup
read_tsv(file = "S6BC_setup.txt")->setup
str_c(setup$Row, setup$Col)->setup$Well

rep("ARE", length(are$Mean))->are$Iron
rep("Fe", length(fe$Mean))->fe$Iron




rbind(are, fe)->od
round(od$`Time [h]`, digits = 3)->od$`Time [h]`
left_join(od, setup)->od
od %>% 
  group_by(`Time [h]`, Iron, R401, Well) %>%
  summarise(OD=mean(Mean))->od
colnames(od)[1]<-c("Time")
filter(od, R401!="X")->od
od %>% 
  filter(R401!="XXX") %>% 
  filter(R401!="none")->od
od %>%
  group_by(R401, Time, Iron) %>%
  summarise(OD_mean=mean(OD), OD_sd=sd(OD))->od_mean

#kinetic
ggplot(od_mean) + 
  #geom_point(alpha = 0.4)+
  geom_ribbon(aes(x = Time, y = OD_mean,  fill=R401, ymin = OD_mean - OD_sd, ymax = OD_mean + OD_sd), alpha=0.5) +
  geom_line(aes(x = Time, y = OD_mean,  color=R401)) +
  #stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), aes())+
  theme(panel.background = element_rect(F))+
  theme(axis.line = element_line(T))+
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
        #legend.position="none",
        strip.background = element_rect(color="black", fill=NA, size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black"))+
  #geom_hline(yintercept = S782_halo_mean)+
  #scale_y_continuous(breaks = c(seq(0,20,5)),minor_breaks=NULL)+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  #theme(axis.text.x = element_text(color="Black", hjust = 1))+
  ylab("OD600")+
  xlab("Time [h]")+
  facet_wrap(~Iron)+
  NULL

#ggsave("Growth_speed_ARE_OD_curves.pdf",  width = 30, height=12, units = "cm")

