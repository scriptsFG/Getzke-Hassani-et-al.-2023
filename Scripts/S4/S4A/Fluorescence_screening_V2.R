#Fluorescence-based screening of R569 mTn5 library for PVD biosynthesis mutants
#Milena Malisic

library(ggplot2)
library(tidyverse)
library(readxl)


getwd()

Theme <- theme(strip.text = element_text(size= 12),
               aspect.ratio = 1,
               plot.background= element_blank(), 
               panel.background = element_blank(),
               axis.line = element_line(size = 0.2),
               axis.title = element_text(size = 9), 
               axis.ticks.y = element_line(size = 0.2),
               axis.text.x = element_text(size = 7, angle = 0, hjust = 0.5, vjust=0.5),
               legend.position = "bottom",
               axis.text.y = element_text(size = 7, angle = 0),
               legend.text = element_text())

#Data

Fluor_Screen <-  read_xlsx("R569_Fluor_Screen.xlsx")

## Calculate median absolute deviation across all screening plates ## 
# Abs600 and Fluorescence was BOTH measured in 2-fold dilution, therefore was not taken into consideration for calculations # 

median_Abs600 <- median(Fluor_Screen$Abs600_raw)

mad_Abs600 <- mad(Fluor_Screen$Abs600_raw)

median_Fluor <- median(Fluor_Screen$Fluorescence_raw)
mad_Fluor <- mad(Fluor_Screen$Fluorescence_raw)

threshold_Abs600 <- median_Abs600-1*mad_Abs600 # 1x absolute deviation from median for Abs600
threshold_Fluor_strict <- median_Fluor-6*mad_Fluor # 6x absolute devation from median for Fluor

Candidates <- (data=subset(Fluor_Screen, Abs600_raw > 0.3994 & Fluorescence_raw < 514.3868))
write.csv(Candidates, "PVD_candidates_6foldmad.csv")


p1 <- ggplot(Fluor_Screen, aes(x=Abs600_raw, y=log10(Fluorescence_raw)))+
  geom_point()+
  xlab("Abs600") +
  ylab("log10_Fluorescence")+ 
  geom_point(data=subset(Fluor_Screen, Abs600_raw > 0.3994 & Fluorescence_raw < 2721.7195), color = "#636363")+
  geom_point(data=subset(Fluor_Screen, Abs600_raw > 0.3994 & Fluorescence_raw < 514.3868), color = "#bdbdbd")+
  geom_point(data=subset(Fluor_Screen, Well_Plate=="17.1_D8"), color = "#577590")+
  geom_point(data=subset(Fluor_Screen, Well_Plate=="19.1_D8"), color = "#F94144")+
  ylab("log10_Fluorescence")+
  Theme
(p1)

ggsave("R569_PVD_Screen.pdf",p1, width = 4, height = 4)