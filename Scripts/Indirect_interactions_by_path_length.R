
#Visualize indirect interactions of motif classes by path length 
library(ggplot2)
library(scatterpie)
library(RColorBrewer)
library(tidyverse)


direct_indirect <- read_csv("Data/Data_processing/Motifs_connections/motif_interactions.csv")
colnames(direct_indirect)
direct_indirect <- direct_indirect %>% select(motif_id,position,direct_interactions, indirect_interactions,ind_same_type_over_directed_int,type)


#Add path length classification
direct_indirect$Broad_categories <- NA
direct_indirect$Broad_categories[direct_indirect$motif_id==2 | direct_indirect$motif_id==7 | 
                                   direct_indirect$motif_id==17 | direct_indirect$motif_id==3 | 
                                   direct_indirect$motif_id==4 | direct_indirect$motif_id==8] <- "Fan" 

direct_indirect$Broad_categories[direct_indirect$motif_id==15 | direct_indirect$motif_id==11] <- "Medium-weak" 

direct_indirect$Broad_categories[direct_indirect$motif_id==6 | direct_indirect$motif_id==16 | 
                                   direct_indirect$motif_id==12] <- "Strong" 

direct_indirect$Broad_categories[direct_indirect$motif_id==5 | direct_indirect$motif_id==9 |  direct_indirect$motif_id==14|  direct_indirect$motif_id==13|  direct_indirect$motif_id==10] <- "Weak" 


#Filter motif 1 and convert to factor
colnames(direct_indirect)
direct_indirect$Broad_categories <- factor(direct_indirect$Broad_categories, levels=c("Strong", "Fan", "Medium-weak", "Weak"))
direct_indirect <- direct_indirect %>% filter(!motif_id==1)

#Convert to factor
direct_indirect$position <- as.factor(direct_indirect$position)
direct_indirect$motif_id <- as.factor(direct_indirect$motif_id)

direct_indirect$ind_same_type_over_directed_int <- as.factor(direct_indirect$ind_same_type_over_directed_int)
 ggplot(direct_indirect %>% filter(type == "pollinator"),aes(motif_id,indirect_interactions))+
  geom_col(position = "dodge")+
  facet_grid(.~Broad_categories,scale='free_x', space = "free_x") + ggtitle("Pollinator")+
  ylab("Number of indirect interactions") + xlab("Motif number")

