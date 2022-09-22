#load libraries
library(tidyverse)
library(PMCMRplus)


rm(list = ls())



read_tsv(file = "ant_sus_scores.txt")->scores


dunn <- kwAllPairsDunnTest(BGCs ~ as.factor(Fraction), data=scores, p.adjust.method ="BH")
dunn

