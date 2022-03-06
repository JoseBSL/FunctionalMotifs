library(ggplot2)
library(scatterpie)
library(RColorBrewer)
library(tidyverse)

plant_means_reordered <- read_csv("Data/Csv/plant_abs_freq_means.csv")
pollinator_means_reordered <- read_csv("Data/Csv/pollinator_abs_freq_means.csv")
#Filter out in order to have just 6 FG'S
pollinator_means_reordered <- pollinator_means_reordered %>% filter(!Node_FG %in% c("Birds", "Lizards", "Other_insects"))


# Pollinator--------
pollinator_trial <- pollinator_means_reordered %>% select(-mean, -SE)# %>% 
  #spread(Node_FG,mean_natural_units)

#Select potions of interest just for 3 species
s <- c(seq(from=3, to=16))
#filter dataset
d <- pollinator_trial %>% filter(position %in% s) %>% filter(Node_FG=="Bee")
d$indirect <- NA
#Add number of indirect interactions manually
#Faster for now, do not know how to code it easily

#3species
d$indirect[d$position=="3"] <- "1"
d$indirect[d$position=="4"] <- "0"
d$indirect[d$position=="5"] <- "0"
d$indirect[d$position=="6"] <- "1"
#4species
d$indirect[d$position=="7"] <- "0"
d$indirect[d$position=="8"] <- "2"
d$indirect[d$position=="9"] <- "1"
d$indirect[d$position=="10"] <- "1"
d$indirect[d$position=="11"] <- "1"
d$indirect[d$position=="12"] <- "1"
d$indirect[d$position=="13"] <- "1"
d$indirect[d$position=="14"] <- "1"
d$indirect[d$position=="15"] <- "2"
d$indirect[d$position=="16"] <- "0"
#5species
d$indirect[d$position=="17"] <- "0"
d$indirect[d$position=="18"] <- "0"
d$indirect[d$position=="19"] <- "0"
d$indirect[d$position=="20"] <- "0"
d$indirect[d$position=="21"] <- "0"
d$indirect[d$position=="22"] <- "0"
d$indirect[d$position=="23"] <- "0"
d$indirect[d$position=="24"] <- "0"
d$indirect[d$position=="25"] <- "0"
d$indirect[d$position=="26"] <- "0"
d$indirect[d$position=="27"] <- "0"
d$indirect[d$position=="28"] <- "0"
d$indirect[d$position=="29"] <- "0"
d$indirect[d$position=="30"] <- "0"
d$indirect[d$position=="31"] <- "0"
d$indirect[d$position=="32"] <- "0"
d$indirect[d$position=="33"] <- "0"
d$indirect[d$position=="34"] <- "0"
d$indirect[d$position=="35"] <- "0"
d$indirect[d$position=="36"] <- "0"
d$indirect[d$position=="37"] <- "0"
d$indirect[d$position=="38"] <- "0"
d$indirect[d$position=="39"] <- "0"
d$indirect[d$position=="40"] <- "0"
d$indirect[d$position=="41"] <- "0"
d$indirect[d$position=="42"] <- "0"
d$indirect[d$position=="43"] <- "0"
d$indirect[d$position=="44"] <- "0"
d$indirect[d$position=="45"] <- "0"




d$position <- as.factor(d$position)
ggplot(d,aes(position,mean_natural_units))+
  geom_col(position = "dodge")+
  facet_grid(.~indirect,scale='free_x', space = "free_x") + ggtitle("Bees")
