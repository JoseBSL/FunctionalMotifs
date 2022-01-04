
#Load libraries
library(tidyverse)
library(kableExtra)

#Read data
networks <- read_csv("Data/Csv/data_for_motifs_analysis_1.csv")

#Subset by interaction greater than 0
d <- subset(networks, Interaction>0)
#Create dataframe
plant_species <- as.data.frame(unique(d$Plant_species))
#Divide dataset in two
x <- plant_species[1:252,]
y <- plant_species[253:503,]
n <- max(length(x), length(y))
length(x) <- n                      
length(y) <- n
#cbind with max length
d <- as.data.frame(cbind(x, y))
#Convert NA to blank spaces
d[is.na(d)] <- ""
#Create table
kbl(d,booktabs = T,col.names = c("Plant species", "Plant species")) %>%
column_spec(1,italic = T) %>%   column_spec(2,italic = T)
  