
# Load libraries

library(lme4)
library(nlme)
library(tidyverse)
library(stringi)
source("Scripts/Scripts_Alfonso/add_study_id.R")


pollinator_positions_codes <- c("np2","np4","np6","np8","np11","np12","np14","np16",
                                "np18","np21","np22","np24","np25","np28","np29","np31",
                                "np34","np35","np38","np41","np42","np44","np46")
plant_positions_codes <- c("np1","np3","np5","np7","np9","np10","np13","np15",
                                "np17","np19","np20","np23","np26","np27","np30","np32",
                                "np33","np36","np37","np39","np40","np43","np45")


#############################################
# EXTRACT MEAN VALUES
#############################################
all_position_percentiles <- read_csv("Data/Csv/Motifs_positions_and_null_models/GF_positions_frequency_percentile.csv")

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

motifs_raw_data <- read_csv("Data/Data_processing/Motifs_connections/motif_pattern_connections.csv")

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


fan_plant <- motifs_raw_data$plant[motifs_raw_data$Broad_categories == "Fan"] %>% unique()
core_per_plant <- motifs_raw_data$plant[motifs_raw_data$Broad_categories == "Core-peripherical"] %>% unique()
complete_plant <- motifs_raw_data$plant[motifs_raw_data$Broad_categories == "Complete"] %>% unique()
asymm_plant <- motifs_raw_data$plant[motifs_raw_data$Broad_categories == "Asymmetric complete"] %>% unique()

plant_means$Broad_categories <- NA

plant_means$Broad_categories[plant_means$position %in% fan_plant] <- "Fan" 
plant_means$Broad_categories[plant_means$position %in% asymm_plant] <- "Medium-weak" 
plant_means$Broad_categories[plant_means$position %in% complete_plant] <- "Strong" 
plant_means$Broad_categories[plant_means$position %in% core_per_plant] <- "Weak" 

###########################
# ESTIMATING SPECIFICITY
###########################
# Reference: 10.1371/journal.pone.0114674

motif_pattern_connections <- read_csv("Data/Data_processing/Motifs_connections/motif_pattern_connections.csv")

number_plant <- motif_pattern_connections %>% group_by(motif_id,plant) %>% 
  count() %>% group_by(motif_id) %>% count() %>% rename(number_plants = n)

number_pollinator <- motif_pattern_connections %>% group_by(motif_id,pollinator) %>% 
  count() %>% group_by(motif_id) %>% count() %>% rename(number_pollinators = n)

specificity_plant_aux <- motif_pattern_connections %>% 
  group_by(motif_id,plant) %>% count()

specificity_pollinator_aux <- motif_pattern_connections %>% 
  group_by(motif_id,pollinator) %>% count()

specificity_plant_aux$position <- stri_extract_first_regex(specificity_plant_aux$plant,
                                                              "[0-9]+") %>% as.numeric()
specificity_pollinator_aux$position <- stri_extract_first_regex(specificity_pollinator_aux$pollinator,
                                                  "[0-9]+") %>% as.numeric()

specificity_plant_aux2 <- specificity_plant_aux %>% 
  left_join(number_pollinator, by = "motif_id")

specificity_pollinator_aux2 <- specificity_pollinator_aux %>% 
  left_join(number_plant, by = "motif_id")

specificity_plant <- specificity_plant_aux2 %>% ungroup %>%
  mutate(s=(number_pollinators-n)/(number_pollinators-1)) %>%
  select(position,s) %>% unique()

specificity_pollinator <- specificity_pollinator_aux2 %>% ungroup %>%
  mutate(s=(number_plants-n)/(number_plants-1)) %>%
  select(position,s) %>% unique()


plant_means_s <- plant_means %>% left_join(specificity_plant, by = "position")

pollinator_means_reordered_s <- pollinator_means_reordered %>% 
  left_join(specificity_pollinator, by = "position")


# We replace NaN of some specialist positions by one
plant_means_s$s[is.nan(plant_means_s$s)] <- 1.0
pollinator_means_reordered_s$s[is.nan(pollinator_means_reordered_s$s)] <- 1.0

library(patchwork)

#select poll functional groups of interest
levels(factor(pollinator_means_reordered_s$Node_FG))
vars <- c("Birds","Lizards", "Other_insects")
pollinator_means_reordered_s_1 <- pollinator_means_reordered_s %>% filter(!Node_FG %in% vars)
  
p1 <- ggplot(pollinator_means_reordered_s_1,aes(x=s,y=mean))+
  geom_point(alpha=0.5, color="black")+
  geom_smooth(method = "lm",color="black")+
  facet_wrap(~Node_FG)+
  labs(x="Specificty",y="Average frequency per position") +
  theme_bw() + ggtitle("Floral visitors") +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- ggplot(plant_means_s,aes(x=s,y=mean))+
  geom_point(alpha=0.5, color="black")+
  geom_smooth(method = "lm", color="black")+
  facet_wrap(~Node_FG)+
  labs(x="Specificty",y="Average frequency per position")+
  theme_bw() + ggtitle("Plants") +
  theme(plot.title = element_text(hjust = 0.5))

p1/p2

#save data 
write.csv(pollinator_means_reordered_s_1, "Data/Csv/pollinator_means_reordered_s_1.csv")
write.csv(plant_means_s, "Data/Csv/plant_means_s.csv")


