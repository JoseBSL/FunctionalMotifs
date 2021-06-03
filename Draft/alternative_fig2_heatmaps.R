
library(lme4)
library(nlme)
library(tidyverse)
library(stringi)
source("../Scripts/Scripts_Alfonso/add_study_id.R")

#############################################
# EXTRACT MEAN VALUES
#############################################


all_position_percentiles <- read_csv("../Data/Csv/Motifs_positions_and_null_models/GF_positions_frequency_percentile.csv")

all_position_percentiles$position <- stri_extract_first_regex(all_position_percentiles$position,
                                                              "[0-9]+") %>% as.numeric()

all_position_percentiles <- add_study_id(all_position_percentiles)
#str(all_position_percentiles)
#head(all_position_percentiles)

list_plants_FG <- as.character(1:10)

plant_position_percentiles_filtered <- all_position_percentiles %>% 
  filter(Node_FG %in% list_plants_FG, !position %in% c(1,2),
         !is.na(expected_freq_its_GF)) #ME DA ERROR
pollinator_position_percentiles_filtered <- all_position_percentiles %>% 
  filter(!Node_FG %in% list_plants_FG, !position %in% c(1,2),
         !is.na(expected_freq_its_GF))


# Plants-------------

# Mean percentile and SE taking into account the study system
plant_codes <- plant_position_percentiles_filtered$position %>% unique()
plant_means <- plant_position_percentiles_filtered %>% select(position,Node_FG) %>% unique()
plant_means$mean <- NA
plant_means$SE <- NA


for(i.pos in 1:length(plant_codes)){
  
  m <- lmer(percentil_its_GF ~ Node_FG + (1|study_id),
            data = plant_position_percentiles_filtered %>%
              filter(position == plant_codes[i.pos]))  
  
  temp <- fixed.effects(m) %>% unname()
  plant_means$mean[plant_means$position == plant_codes[i.pos]] <- c(temp[1],
                                                                    temp[1]+temp[2],
                                                                    temp[1]+temp[3],
                                                                    temp[1]+temp[4],
                                                                    temp[1]+temp[5]) ## to get real means per FG
  temp <- sqrt(diag(vcov(m)))
  plant_means$SE[plant_means$position == plant_codes[i.pos]] <- c(temp[1], 
                                                                  temp[1]+temp[2],
                                                                  temp[1]+temp[3],
                                                                  temp[1]+temp[4],
                                                                  temp[1]+temp[5])
  
}


# Pollinators--------

# Mean percentile and SE taking into account the study system
pollinator_codes <- pollinator_position_percentiles_filtered$position %>% unique()
pollinator_means <- pollinator_position_percentiles_filtered %>% select(position,Node_FG) %>% unique()
pollinator_means$mean <- NA
pollinator_means$SE <- NA
pollinator_means_reordered <- NULL

for(i.pos in 1:length(pollinator_codes)){
  
  m <- lmer(percentil_its_GF ~ Node_FG + (1|study_id), #HERE WE NEED STUDY SYSTEM ONLY!! 
            data = pollinator_position_percentiles_filtered %>%
              filter(position == pollinator_codes[i.pos]),
            control = lmerControl(optimizer ="Nelder_Mead")) # I changed the optimizer because motif 8 (i.pos = 3) caused convergence problems, max gradient = 0.0022 > 0.002.
  
  temp <- fixed.effects(m) %>% unname()
  pollinator_means_aux <- pollinator_means[pollinator_means$position == pollinator_codes[i.pos],] %>% arrange(Node_FG)
  
  pollinator_means_aux$mean <- c(temp[1],
                                 temp[1]+temp[2],
                                 temp[1]+temp[3],
                                 temp[1]+temp[4],
                                 temp[1]+temp[5],
                                 temp[1]+temp[6],
                                 temp[1]+temp[7],
                                 temp[1]+temp[8],
                                 temp[1]+temp[9]) ## to get real means per FG
  
  
  temp <- sqrt(diag(vcov(m)))
  pollinator_means_aux$SE <- c(temp[1],
                               temp[1]+temp[2],
                               temp[1]+temp[3],
                               temp[1]+temp[4],
                               temp[1]+temp[5],
                               temp[1]+temp[6],
                               temp[1]+temp[7],
                               temp[1]+temp[8],
                               temp[1]+temp[9])
  
  pollinator_means_reordered <- bind_rows(pollinator_means_reordered,pollinator_means_aux)
  
}

# There are lower limits (mean - SE) smaller than 0 (see Other insects)
# We set the lower limit of those error bars to zero

pollinator_means_reordered <- pollinator_means_reordered %>% mutate(lower = mean-SE)
pollinator_means_reordered$lower[pollinator_means_reordered$lower < 0] <- 0

#############################################
# ADDING BROAD MOTIF CATEGORIES TO POSITIONS
#############################################

motifs_raw_data <- read_csv("../Data/Data_processing/Motifs_connections/motif_pattern_connections.csv")

motifs_raw_data$Broad_categories <- NA

motifs_raw_data$Broad_categories[motifs_raw_data$motif_id %in% c(5,9,14,13,10)] <- "Core-peripherical"
motifs_raw_data$Broad_categories[motifs_raw_data$motif_id %in% c(6,16,12)] <- "Complete"
motifs_raw_data$Broad_categories[motifs_raw_data$motif_id %in% c(2,3,7,4,17,8)] <- "Fan"
motifs_raw_data$Broad_categories[motifs_raw_data$motif_id %in% c(15,11)] <- "Asymmetric complete"

motifs_raw_data$plant <- gsub("[^0-9.-]", "", motifs_raw_data$plant)
motifs_raw_data$pollinator <- gsub("[^0-9.-]", "", motifs_raw_data$pollinator)
motifs_raw_data$plant <- as.numeric(motifs_raw_data$plant)
motifs_raw_data$pollinator <- as.numeric(motifs_raw_data$pollinator)
motifs_raw_data <- motifs_raw_data %>% filter(!is.na(Broad_categories))

fan_pollinator <- motifs_raw_data$pollinator[motifs_raw_data$Broad_categories == "Fan"] %>% unique()
core_per_pollinator <- motifs_raw_data$pollinator[motifs_raw_data$Broad_categories == "Core-peripherical"] %>% unique()
complete_pollinator <- motifs_raw_data$pollinator[motifs_raw_data$Broad_categories == "Complete"] %>% unique()
asymm_pollinator <- motifs_raw_data$pollinator[motifs_raw_data$Broad_categories == "Asymmetric complete"] %>% unique()

pollinator_means_reordered$Broad_categories <- NA

pollinator_means_reordered$Broad_categories[pollinator_means_reordered$position %in% fan_pollinator] <- "Fan" 
pollinator_means_reordered$Broad_categories[pollinator_means_reordered$position %in% asymm_pollinator] <- "Medium-weak" 
pollinator_means_reordered$Broad_categories[pollinator_means_reordered$position %in% complete_pollinator] <- "Strong" 
pollinator_means_reordered$Broad_categories[pollinator_means_reordered$position %in% core_per_pollinator] <- "Weak" 

ggplot(data = pollinator_means_reordered) + 
  geom_bar(aes(x=as.factor(Node_FG), y=mean, 
                 fill = as.factor(Broad_categories)),stat="identity")+
  geom_errorbar(aes(x = as.factor(Node_FG), ymin=lower, ymax=mean+SE), 
                width=1.0,size=1)+
  facet_wrap(~position,ncol = 5)+
  geom_point(aes(x = as.factor(Node_FG), y=mean), 
             size=2)+
  labs(x=NULL, y = "Percentile",fill =NULL)+
  theme_bw()+
  theme(legend.position = "bottom")+#,axis.text.x=element_text(angle=90,vjust=0.5, hjust=1))+
  coord_flip()



fan_plant <- motifs_raw_data$plant[motifs_raw_data$Broad_categories == "Fan"] %>% unique()
core_per_plant <- motifs_raw_data$plant[motifs_raw_data$Broad_categories == "Core-peripherical"] %>% unique()
complete_plant <- motifs_raw_data$plant[motifs_raw_data$Broad_categories == "Complete"] %>% unique()
asymm_plant <- motifs_raw_data$plant[motifs_raw_data$Broad_categories == "Asymmetric complete"] %>% unique()

plant_means$Broad_categories <- NA

plant_means$Broad_categories[plant_means$position %in% fan_plant] <- "Fan" 
plant_means$Broad_categories[plant_means$position %in% asymm_plant] <- "Medium-weak" 
plant_means$Broad_categories[plant_means$position %in% complete_plant] <- "Strong" 
plant_means$Broad_categories[plant_means$position %in% core_per_plant] <- "Weak" 


ggplot(data = plant_means) + 
  geom_bar(aes(x=as.factor(Node_FG), y=mean, 
               fill = as.factor(Broad_categories)),stat="identity")+
  geom_errorbar(aes(x = as.factor(Node_FG), ymin=mean-SE, ymax=mean+SE), 
                width=1.0,size=1)+
  facet_wrap(~position,ncol = 5)+
  geom_point(aes(x = as.factor(Node_FG), y=mean), 
             size=2)+
  labs(x=NULL, y = "Percentile",fill = NULL)+
  theme_bw()+
  theme(legend.position = "bottom")+#,axis.text.x=element_text(angle=90,vjust=0.5, hjust=1))+
  coord_flip()

########################################
# CREATE HEATMAPS
########################################

library(ComplexHeatmap)

# PLANTS-----

# Prepare data for heatmap
heat_plants_aux <- plant_means %>% spread(Node_FG,mean) %>% select(-SE) 

heat_plants <- heat_plants_aux %>%
  group_by(position,Broad_categories) %>% 
  summarise_all(sum,na.rm = T)

heat_plants_fan <- heat_plants %>% filter(Broad_categories == "Fan") %>% select(-Broad_categories) %>% 
  as.matrix()

rownames(heat_plants_fan) <- heat_plants_fan[,1]

heat_plants_Weak <- heat_plants %>% filter(Broad_categories == "Weak") %>% select(-Broad_categories) %>% 
  as.matrix()

rownames(heat_plants_Weak) <- heat_plants_Weak[,1]

heat_plants_Strong <- heat_plants %>% filter(Broad_categories == "Strong") %>% select(-Broad_categories) %>% 
  as.matrix()

rownames(heat_plants_Strong) <- heat_plants_Strong[,1]

heat_plants_Medium_weak <- heat_plants %>% filter(Broad_categories == "Medium-weak") %>% select(-Broad_categories) %>% 
  as.matrix()

rownames(heat_plants_Medium_weak) <- heat_plants_Medium_weak[,1]

par(mfrow=c(2,2))

hp1 <- Heatmap(heat_plants_fan[,c(2:ncol(heat_plants_fan))], 
        name = "mean\npercentile", #title of legend
        column_title = "Functional\ Groups", row_title = "Fan",
        row_names_gp = gpar(fontsize = 12) # Text size for row names
)

hp2 <- Heatmap(heat_plants_Weak[,c(2:ncol(heat_plants_Weak))], 
        name = "mean\npercentile", #title of legend
        column_title = "Functional\ Groups", row_title = "Weak",
        row_names_gp = gpar(fontsize = 12) # Text size for row names
)

hp3 <- Heatmap(heat_plants_Strong[,c(2:ncol(heat_plants_Strong))], 
        name = "mean\npercentile", #title of legend
        column_title = "Functional\ Groups", row_title = "Strong",
        row_names_gp = gpar(fontsize = 12) # Text size for row names
)

hp4 <- Heatmap(heat_plants_Medium_weak[,c(2:ncol(heat_plants_Medium_weak))], 
        name = "mean\npercentile", #title of legend
        column_title = "Functional\ Groups", row_title = "Medium-weak",
        row_names_gp = gpar(fontsize = 12) # Text size for row names
)

hp1 %v% hp2 %v% hp3 %v% hp4

# POLLINATORS-----------     

filter_rare_groups <- T

if(filter_rare_groups == T){
  heat_pollinators_aux <- pollinator_means_reordered %>% filter(!Node_FG %in% c("Birds",
                                                                                "Lizards",
                                                                                "Other_insects")) %>% 
    spread(Node_FG,mean) %>% select(-SE,-lower)
  
  colnames(heat_pollinators_aux) <- c("position","Broad_categories","Bee",            
                                      "Coleoptera","Lepidoptera","Non-bee\nHymenoptera", 
                                      "Non-syrphids\ndiptera","Syrphids")
}else{
  heat_pollinators_aux <- pollinator_means_reordered %>% spread(Node_FG,mean) %>% select(-SE,-lower)
  
  colnames(heat_pollinators_aux) <- c("position","Broad_categories","Bee","Birds",            
                                      "Coleoptera","Lepidoptera","Lizards","Non-bee\nHymenoptera", 
                                      "Non-syrphids\ndiptera","Other insects","Syrphids")
}


heat_pollinators <- heat_pollinators_aux %>%
  group_by(position,Broad_categories) %>% 
  summarise_all(sum,na.rm = T)

heat_pollinators_fan <- heat_pollinators %>% filter(Broad_categories == "Fan") %>% select(-Broad_categories) %>% 
  as.matrix()

rownames(heat_pollinators_fan) <- heat_pollinators_fan[,1]

heat_pollinators_Weak <- heat_pollinators %>% filter(Broad_categories == "Weak") %>% select(-Broad_categories) %>% 
  as.matrix()

rownames(heat_pollinators_Weak) <- heat_pollinators_Weak[,1]

heat_pollinators_Strong <- heat_pollinators %>% filter(Broad_categories == "Strong") %>% select(-Broad_categories) %>% 
  as.matrix()

rownames(heat_pollinators_Strong) <- heat_pollinators_Strong[,1]

heat_pollinators_Medium_weak <- heat_pollinators %>% filter(Broad_categories == "Medium-weak") %>% select(-Broad_categories) %>% 
  as.matrix()

rownames(heat_pollinators_Medium_weak) <- heat_pollinators_Medium_weak[,1]


h1 <- Heatmap(heat_pollinators_fan[,c(2:ncol(heat_pollinators_fan))], 
              name = "mean\npercentile", #title of legend
              column_title = "Functional\ Groups", row_title = "Fan",
              row_names_gp = gpar(fontsize = 12) # Text size for row names
)

h2 <- Heatmap(heat_pollinators_Weak[,c(2:ncol(heat_pollinators_Weak))], 
              name = "mean\npercentile", #title of legend
              column_title = "Functional\ Groups", row_title = "Weak",
              row_names_gp = gpar(fontsize = 12) # Text size for row names
)

h3 <- Heatmap(heat_pollinators_Strong[,c(2:ncol(heat_pollinators_Strong))], 
              name = "mean\npercentile", #title of legend
              column_title = "Functional\ Groups", row_title = "Strong",
              row_names_gp = gpar(fontsize = 12) # Text size for row names
)

h4 <- Heatmap(heat_pollinators_Medium_weak[,c(2:ncol(heat_pollinators_Medium_weak))], 
              name = "mean\npercentile", #title of legend
              column_title = "Functional\ Groups", row_title = "Medium-weak",
              row_names_gp = gpar(fontsize = 12) # Text size for row names
)

h1 %v% h2 %v% h3 %v% h4
