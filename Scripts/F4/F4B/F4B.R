#load libraries
library(tidyverse)


rm(list = ls())

colors<-c("n"="#90B597", "y"="#DAE6DC")


#####load all files 

#read all tables individually
read_tsv(file = "Measurements21.1.20.tsv")->df21_1
subset(df21_1, select = -c(Measurement_time))->df21_1
read_tsv(file = "Measurements30.1.20.tsv")->df30_1
subset(df30_1, select = -c(Measurement_time))->df30_1
read_tsv(file = "Measurements7.2.20.tsv")->df7_2
subset(df7_2, select = -c(Measurement_time))->df7_2
read_tsv(file = "Measurements14.2.20.tsv")->df14_2
subset(df14_2, select = -c(Measurement_time))->df14_2
read_tsv(file = "Measurements21.2.20.tsv.tsv")->df21_2
subset(df21_2, select = -c(Measurement_time))->df21_2
#generate one df
rbind(df21_1, df30_1, df7_2, df14_2)->df_parsed 

###remove plates that need to be repeated
filter(df_parsed, Filename!="P55.xlsx")->df_parsed
filter(df_parsed, Filename!="P72.xlsx")->df_parsed
filter(df_parsed, Filename!="P80.xlsx")->df_parsed
filter(df_parsed, Filename!="P55_Rs.xlsx")->df_parsed
filter(df_parsed, Filename!="P72_Rs.xlsx")->df_parsed
filter(df_parsed, Filename!="P80_Rs.xlsx")->df_parsed

#generate one df
rbind(df_parsed, df21_2)->df_parsed 
###modify labeling
#remove .xlxs
df_parsed$Filename %>%
  str_remove_all(".xlsx")->df_parsed$Filename
#split P and Rs
separate(df_parsed, Filename, c("plate", "cond"), "_")->df_parsed
df_parsed$cond[is.na(df_parsed$cond)] <- "R401"
###spread the data frame
#spread all measurements into different cols
subset(df_parsed, select = -c(Measurement_mode))->df_parsed
spread(df_parsed, Measurement_wavelength, Value, fill = NA, convert = T, drop = TRUE, sep = NULL)->df_parsed
(df_parsed)
#spread condition
filter(df_parsed, cond=="Rs")->df_Rs
colnames(df_Rs)<-c("x", "y", "Rs", "Rs_475nm/510nm", "Rs_485nm/535nm", "Rs_600")
filter(df_parsed, cond=="R401")->df_R401
colnames(df_R401)<-c("well", "plate", "R401", "R401_475nm/510nm", "R401_485nm/535nm", "R401_600")
cbind(df_R401, df_Rs)->df_parsed
subset(df_parsed, select = -c(x, y, R401, Rs))->df_parsed
#split wells into rows and cols
substr(df_parsed$well, 1, 1)->df_parsed$row
substr(df_parsed$well, 2, 3)->df_parsed$col
###add boarder labeling
arrange(df_parsed, plate)->df_parsed
##########
##########
##########
##########modify if more data is added
rep(c("y", "y", "y", "y", "y", "y", "y", "y", "y", "y", "y", "y", "y", "n", "n", "y", "n", "n", "n", "n", "n", "n", "n", "n", "y", "n", "n", "y", "n", "n", "n", "n", "n", "n", "n", "n", "y", "n", "n", "y", "n", "n", "n", "n", "n", "n", "n", "n", "y", "n", "n", "y", "n", "n", "n", "n", "n", "n", "n", "n", "y", "n", "n", "y", "n", "n", "n", "n", "n", "n", "n", "n", "y", "n", "n", "y", "n", "n", "n", "n", "n", "n", "n", "n", "y", "y", "y", "y", "y", "y", "y", "y", "y", "y", "y", "y"), 91)->borders
as.data.frame(borders)->borders
cbind(df_parsed, borders)->df_parsed
#reposition cols
df_parsed %>% select(borders, everything())->df_parsed
df_parsed %>% select(col, everything())->df_parsed
df_parsed %>% select(row, everything())->df_parsed
df_parsed %>% select(well, everything())->df_parsed

###load other data with different format
read.table(file = "plate84_105to120.txt", header=T)->df_rest
as.integer(df_parsed$col)->df_parsed$col
colnames(df_parsed)->colnames(df_rest)
bind_rows(df_parsed, df_rest)->df_parsed


#take snapshot of df for plotting
df_parsed->df_plot
###split everything by borders
#fuse plate and well info
df_parsed %>% 
  unite(plate_well, plate, well, sep = "_", remove = FALSE)->df_parsed
#generate df for borders
filter(df_parsed, borders=="y")->df_borders
colnames(df_borders)<-c("plate_well", "well", "row", "col", "borders", "plate", "B_R401_475nm/510nm", "B_R401_485nm/535nm", "B_R401_600", "B_Rs_475nm/510nm", "B_Rs_485nm/535nm", "B_Rs_600")
#generate df for central wells
filter(df_parsed, borders=="n")->df_center
colnames(df_center)<-c("plate_well", "well", "row", "col", "borders", "plate", "C_R401_475nm/510nm", "C_R401_485nm/535nm", "C_R401_600", "C_Rs_475nm/510nm", "C_Rs_485nm/535nm", "C_Rs_600")
#join both df back together
bind_rows(df_center, df_borders)->df_parsed


############################extracting outliers using 3x absolute deviation around median (https://www.sciencedirect.com/science/article/pii/S0022103113000668#f0005)
###extract all "normal growers"
#caluclate median per plate
df_parsed %>%
  group_by(plate) %>%
  summarise_at(vars(`C_R401_600`), funs(median(., na.rm=TRUE)))->`C_R401_600_MEDIAN`
colnames(`C_R401_600_MEDIAN`)<-c("plate", "C_R401_600_MEDIAN")
#caluclate mad per plate
df_parsed %>%
  group_by(plate) %>%
  summarise_at(vars(`C_R401_600`), funs(mad(., na.rm=TRUE)))->`C_R401_600_MAD`
colnames(`C_R401_600_MAD`)<-c("plate", "C_R401_600_MAD")
#calculate Threshhold 
`C_R401_600_MEDIAN`$C_R401_600_MEDIAN-3*`C_R401_600_MAD`$C_R401_600_MAD->`C_R401_600_MAD`$`C_R401_600_TH`
#join TH with df
left_join(df_parsed, C_R401_600_MAD)->df_parsed

#extract all values in which C_R401_600 exceeds threshhold
filter(df_parsed, df_parsed$`C_R401_600`> df_parsed$`C_R401_600_TH`)->df_parsed


###extract all Fluorescence values that exceed the 
#caluclate median per plate
df_parsed %>%
  group_by(plate) %>%
  summarise_at(vars(`C_Rs_475nm/510nm`), funs(median(., na.rm=TRUE)))->`C_Rs_475nm/510nm_MEDIAN`
colnames(`C_Rs_475nm/510nm_MEDIAN`)<-c("plate", "C_Rs_475nm/510nm_MEDIAN")
#caluclate mad per plate
df_parsed %>%
  group_by(plate) %>%
  summarise_at(vars(`C_Rs_475nm/510nm`), funs(mad(., na.rm=TRUE)))->`C_Rs_475nm/510nm_MAD`
colnames(`C_Rs_475nm/510nm_MAD`)<-c("plate", "C_Rs_475nm/510nm_MAD")
#calculate Threshhold 
`C_Rs_475nm/510nm_MEDIAN`$`C_Rs_475nm/510nm_MEDIAN`+3*`C_Rs_475nm/510nm_MAD`$`C_Rs_475nm/510nm_MAD`->`C_Rs_475nm/510nm_MAD`$`C_Rs_475nm/510nm_TH`
#join TH with df
left_join(df_parsed, `C_Rs_475nm/510nm_MAD`)->df_parsed

#extract all values in which C_Rs_475nm/510nm exceeds TH
filter(df_parsed, df_parsed$`C_Rs_475nm/510nm`> df_parsed$`C_Rs_475nm/510nm_TH`)->df_parsed


######add information about potentially contaminated wells:
read_tsv(file = "contaminated_wells.txt")->cont_wells
fill(cont_wells, plate)->cont_wells
separate(cont_wells, plate, c("plate", "cond"), "_")->cont_wells
cont_wells$cond[is.na(cont_wells$cond)] <- "R401"
cont_wells$cont<-c(rep("y", 271))
left_join(cont_wells, df_plot)->cont_wells
select(cont_wells, plate, well, cont)->cont_wells
left_join(df_plot, cont_wells)->df_plot
df_plot[is.na(df_plot)] <- "n"
filter(df_plot, cont=="n")->df_plot

#plot by plate

#extract candidate information
#write_tsv(df_parsed, "candidates.txt")


# #caluclate mean per plate
# df_parsed %>%
#   group_by(plate) %>%
#   summarise_at(vars(`C_Rs_475nm/510nm`), funs(mean(., na.rm=TRUE)))->`C_Rs_475nm/510nm_Mean`
# colnames(`C_Rs_475nm/510nm_Mean`)<-c("plate", "C_Rs_475nm/510nm_Mean")
# left_join(df_parsed, `C_Rs_475nm/510nm_Mean`)->df_parsed
# #extract all values in which C_Rs_475nm/510nm surpases the plates 97th percentile
# filter(df_parsed, df_parsed$`C_Rs_475nm/510nm`> 2*df_parsed$`C_Rs_475nm/510nm_Mean`)->df_parsed

############################extracting outliers using quantiles
###extract all "normal growers"
#caluclate percentiles per plate
# df_parsed %>% 
#   group_by(plate) %>% 
#   summarize(quants = quantile(`C_R401_600`, probs = c(0.10), na.rm = T))->`C_R401_600_50Q`
# colnames(`C_R401_600_50Q`)<-c("plate", "C_R401_600_50Q")
# left_join(df_parsed, `C_R401_600_50Q`)->df_parsed
# #extract all values in which C_Rs_475nm/510nm surpases the plates 97th percentile
# filter(df_parsed, df_parsed$`C_R401_600`>df_parsed$`C_R401_600_50Q`)->df_parsed
# 
# ###extract all Fluorescence values that exceed the 95th percentile
# #caluclate percentiles per plate
# df_parsed %>% 
#   group_by(plate) %>% 
#   summarize(quants = quantile(`C_Rs_475nm/510nm`, probs = c(0.95), na.rm = T))->`C_Rs_475nm/510nm_97Q`
# colnames(`C_Rs_475nm/510nm_97Q`)<-c("plate", "C_Rs_475nm/510nm_97Q")
# left_join(df_parsed, `C_Rs_475nm/510nm_97Q`)->df_parsed
# #extract all values in which C_Rs_475nm/510nm surpases the plates 97th percentile
# filter(df_parsed, df_parsed$`C_Rs_475nm/510nm` > df_parsed$`C_Rs_475nm/510nm_97Q`)->df_parsed
# 


#plot distribution
ggplot(df_plot, aes(`Rs_475nm/510nm`)) +
  geom_histogram()


#plot data against perfect normal distribution
ggplot(df_plot, aes(sample = `Rs_475nm/510nm`)) +
  geom_qq() +
  geom_qq_line(col = "red")

#modify below to plot individual plate
# a<-"P108"
# filter(df_plot, plate==a)->df_plot2
# filter(df_parsed, plate==a)->df_parsed2
# 
# filter(`C_Rs_475nm/510nm_MAD`, plate==a)->MAD2
# filter(C_R401_600_MAD, plate==a)->C_R401_600_MAD2




ggplot()+
  geom_point(df_plot, mapping=aes(x=`R401_600`, y=`Rs_475nm/510nm`, alpha=0.5, color=borders, fill=borders))+
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
  #geom_hline(yintercept = 3000, size=1)+
  #geom_vline(xintercept = 0.75, size=1)+
  #scale_size_manual(values=c(1, 2))+
  #scale_y_discrete(limits=c(0, -2, -4, -6, -8))+
  #scale_x_continuous(limits=c(0, 1), breaks = c(0,0.2, 0.4, 0.6, 0.8, 1))+
  #geom_vline(xintercept = 0, size=0.5)+
  ylab("Fluo_Rs+R401")+
  xlab("OD600_R401")+
  geom_point(data=df_parsed, mapping=aes(x=C_R401_600, y=`C_Rs_475nm/510nm`))+
  geom_hline(data=`C_Rs_475nm/510nm_MAD`, mapping=aes(yintercept=`C_Rs_475nm/510nm_TH`))+
  geom_vline(data=`C_R401_600_MAD`, mapping=aes(xintercept=`C_R401_600_TH`))+
  #geom_point(data=cont_wells, mapping=aes(x=R401_600, y=`Rs_475nm/510nm`))+
  facet_wrap(.~plate, ncol = 12)
  
  

#ggsave("R401_Tn5_all_plates.pdf", width = 50, height = 50, limitsize = F, units = c("cm"))

#ggsave("R401_Tn5_P108.pdf")

#####normalize each plate by the average fluorescence and OD
filter(df_plot, borders=="n")->normFluo
aggregate(normFluo$`Rs_475nm/510nm`, by=list(normFluo$plate), FUN=median)->normFluo
colnames(normFluo)<-c("plate",  "normFluo")
aggregate(df_plot$R401_600, by=list(df_plot$plate), FUN=median)->normOD
colnames(normOD)<-c("plate", "normOD")
left_join(df_plot, normFluo)->df_plot
left_join(df_plot, normOD)->df_plot
df_plot$`Rs_475nm/510nm`/df_plot$normFluo->df_plot$newFluo
df_plot$R401_600/df_plot$normOD->df_plot$newOD
df_plot %>% 
  filter(well!="A1" & well!="A12" & well!="H1" & well!="H12")->df_plot  


#generate new dataframe to highlight candidates
select(df_parsed, plate, well)->df_highlights
inner_join(df_highlights, df_plot)->df_highlights


read_tsv("candidates_names.txt")->candidates38
inner_join(candidates38, df_plot)->candidates38

#P1G10
df_plot %>%
  filter(plate=="P37") %>%
  filter(well=="G7")->P1G10



ggplot()+
  geom_point(df_plot, mapping=aes(x=R401_600, y=newFluo, color=borders, fill=borders), alpha=0.8,size=2)+
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
  #geom_hline(yintercept = 3000, size=1)+
  #geom_vline(xintercept = 0.75, size=1)+
  #scale_size_manual(values=c(1, 2))+
  scale_y_continuous(limits=c(0,6))+
  #scale_x_continuous(limits=c(0, 1), breaks = c(0,0.2, 0.4, 0.6, 0.8, 1))+
  #geom_vline(xintercept = 0, size=0.5)+
  ylab("GFP-fluorescence [norm. RFU]")+
  xlab("OD600")+
  geom_point(data=df_highlights, mapping=aes(x=R401_600, y=newFluo), size=2, color='#44664B')+
  geom_point(data=candidates38, mapping=aes(x=R401_600, y=newFluo), size=2, color='#1E2D21')+
  #geom_point(data=P1G10, mapping=aes(x=R401_600, y=newFluo), size=2, color="#FFB703")+
  #geom_hline(data=`C_Rs_475nm/510nm_MAD`, mapping=aes(yintercept=`C_Rs_475nm/510nm_TH`))+
  #geom_vline(data=`C_R401_600_MAD`, mapping=aes(xintercept=`C_R401_600_TH`))
  #geom_point(data=cont_wells, mapping=aes(x=R401_600, y=`Rs_475nm/510nm`))+
  #facet_wrap(.~plate, ncol = 12)
  NULL
#ggsave("R401_Tn5_all_plates_norm.pdf")

df_plot %>% 
  filter(borders=="n") %>% 
  filter(R401_600<=1) ->df_plot2


ggplot(df_plot2, aes(x=R401_600, y=newFluo, alpha=0.5, color=borders, fill=borders))+
  geom_point(aes(x=R401_600, y=newFluo, alpha=0.5, color=borders, fill=borders))+
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
  #geom_hline(yintercept = 3000, size=1)+
  #geom_vline(xintercept = 0.75, size=1)+
  #scale_size_manual(values=c(1, 2))+
  #scale_y_discrete(limits=c(0, -2, -4, -6, -8))+
  #scale_x_continuous(limits=c(0, 1), breaks = c(0,0.2, 0.4, 0.6, 0.8, 1))+
  #geom_vline(xintercept = 0, size=0.5)+
  ylab("GFP-fluorescence [norm. RFU]")+
  xlab("OD600")+
  ggtitle("p1.155269e-280; R0.6423")+
  geom_smooth(method='lm', color='black', se=F)+
  NULL

#ggsave("growth_speed_effect.pdf")

##stats
cor.test(x=df_plot2$R401_600, y=df_plot2$newFluo, method=c("pearson"), exact = F)->stats
stats$p.value

summary(lm(df_plot2$R401_600~df_plot2$newFluo))











