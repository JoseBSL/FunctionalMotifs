
# Load libraries
library(nlme)
library(tidyverse)
library(stringi)
source("Scripts/Scripts_Alfonso/add_study_id.R")
library(glmmTMB)
library(lme4)
library(visreg)


############################################
# RECALCULATING RESULTS WITH FEWER NETWORKS
############################################


percentage_studies <- .60
repetitions <- 1000

################################################################################
# Initializing variables

pollinator_positions_codes <- c("np2","np4","np6","np8","np11","np12","np14","np16",
                                "np18","np21","np22","np24","np25","np28","np29","np31",
                                "np34","np35","np38","np41","np42","np44","np46")
plant_positions_codes <- c("np1","np3","np5","np7","np9","np10","np13","np15",
                                "np17","np19","np20","np23","np26","np27","np30","np32",
                                "np33","np36","np37","np39","np40","np43","np45")


shade_data_pollinator_sp_plot <- NULL
mean_line_data_pollinator_sp_plot <- NULL
shade_data_pollinator_int_plot <- NULL
mean_line_data_pollinator_int_plot <- NULL
shade_data_plant_sp_plot <- NULL
mean_line_data_plant_sp_plot <- NULL
shade_data_plant_int_plot <- NULL
mean_line_data_plant_int_plot <- NULL

#############################################
# RANDOM CALCULATIONS
#############################################

set.seed(1234)

for (rep_i in 1:repetitions) {
  
  #############################################
  # SELECT STUDIES
  #############################################
  
  all_position_percentiles_raw <- read_csv("Data/Csv/Motifs_positions_and_null_models/GF_positions_frequency_percentile.csv")
  
  all_network_ids <- unique(all_position_percentiles_raw$Network_id)
  
  selected_network_indexes <- sample.int(length(all_network_ids), 
                                         round(percentage_studies*
                                                 length(all_network_ids)),
                                         replace = FALSE) %>% sort()
  
  all_position_percentiles <- all_position_percentiles_raw %>% 
    filter(Network_id %in% all_network_ids[selected_network_indexes])
  
  #############################################
  # EXTRACT MEAN VALUES
  #############################################
  
  all_position_percentiles$position <- stri_extract_first_regex(all_position_percentiles$position,
                                                                "[0-9]+") %>% as.numeric()
  
  all_position_percentiles <- add_study_id(all_position_percentiles)
  
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
  
  
  
  pollinator_freq_position_s <- pollinator_position_percentiles_filtered %>% 
    left_join(specificity_pollinator, by = "position")
  
  pollinator_freq_position_s$s[is.nan(pollinator_freq_position_s$s)] <- 1.0
  
  ###########################
  # ADDING INDIRECT INT DATA
  ###########################
  
  motif_interactions <- read_csv("Data/Data_processing/Motifs_connections/motif_interactions.csv") %>%
    select(position,indirect_interactions) %>% unique()
  
  pollinator_freq_position_s_ind <- pollinator_freq_position_s %>% 
    left_join(motif_interactions, by = "position") %>% 
    filter(!Node_FG %in% c("Birds","Lizards","Other_insects")) %>%
    rename(Percentil=percentil_its_GF,Specificity=s,Indirect_interactions=indirect_interactions,Group=Node_FG)
  
  
  pollinator_freq_position_s_ind$position <- as.factor(pollinator_freq_position_s_ind$position)
  
  model_freq_specif_ind <- glmmTMB(Percentil ~ Specificity*Group +
                                     scale(Indirect_interactions)*Group + 
                                     (1|study_id/position),
                                   pollinator_freq_position_s_ind)
  
  
  pollinators_sp_plot <- visreg(model_freq_specif_ind, "Specificity", 
                                by="Group", 
                                scale = "response",
                                strip.names=c("Bee",
                                              "Coleoptera",
                                              "Lepidoptera",
                                              "Non-bee-Hymenoptera",
                                              "Non-syrphids-diptera",
                                              "Syrphids"),gg = TRUE) + theme_bw()
  
  pollinators_ind_int_plot <- visreg(model_freq_specif_ind, "Indirect_interactions", 
                                     by="Group", 
                                     scale = "response",
                                     strip.names=c("Bee",
                                                   "Coleoptera",
                                                   "Lepidoptera",
                                                   "Non-bee-Hymenoptera",
                                                   "Non-syrphids-diptera",
                                                   "Syrphids"),gg = TRUE)+ theme_bw()+xlab("Indirect interactions")
  

  
  ##############################################################
  # EXTRACT DATA FROM POLLINATOR'S
  ##############################################################
  
  df_from_pollinators_sp_plot <- ggplot_build(pollinators_sp_plot)
  df_from_pollinators_int_plot <- ggplot_build(pollinators_ind_int_plot)
  
  # Extract data from pollinator's specificity plot--------
  
  
  
  shade_data_x_pollinator_sp_plot <- df_from_pollinators_sp_plot$data[[1]] %>% select(x) %>% unique() %>% pull()
  
  for (x_i in 1:length(shade_data_x_pollinator_sp_plot)) {
    shade_data_aux_pollinator_sp_plot <- df_from_pollinators_sp_plot$data[[1]] %>% filter(x == shade_data_x_pollinator_sp_plot[x_i]) %>% 
      group_by(PANEL) %>% summarise(ymin = min(y), ymax = max(y)) %>%
      mutate(x = shade_data_x_pollinator_sp_plot[x_i],
             repetition = rep_i)
    shade_data_pollinator_sp_plot <- bind_rows(shade_data_pollinator_sp_plot,
                                               shade_data_aux_pollinator_sp_plot)
  }
  
  mean_line_data_pollinator_sp_plot <- bind_rows(mean_line_data_pollinator_sp_plot,
                                                 df_from_pollinators_sp_plot$data[[2]] %>%
                                                   mutate(repetition = rep_i))
  
  # Extract data from pollinator's ind. interactions' plot--------
  
  
  shade_data_x_pollinator_int_plot <- df_from_pollinators_int_plot$data[[1]] %>% select(x) %>% unique() %>% pull()
  
  for (x_i in 1:length(shade_data_x_pollinator_int_plot)) {
    shade_data_aux_pollinator_int_plot <- df_from_pollinators_int_plot$data[[1]] %>% filter(x == shade_data_x_pollinator_int_plot[x_i]) %>% 
      group_by(PANEL) %>% summarise(ymin = min(y), ymax = max(y)) %>%
      mutate(x = shade_data_x_pollinator_int_plot[x_i],
             repetition = rep_i)
    shade_data_pollinator_int_plot <- bind_rows(shade_data_pollinator_int_plot,
                                                shade_data_aux_pollinator_int_plot)
  }
  
  mean_line_data_pollinator_int_plot <- bind_rows(mean_line_data_pollinator_int_plot,
                                                  df_from_pollinators_int_plot$data[[2]] %>%
                                                    mutate(repetition = rep_i))
  
  
  
  #---------------------------------
  # Percentil analises with LMMs for Plants
  
  
  plant_freq_position_s <- plant_position_percentiles_filtered %>% 
    left_join(specificity_plant, by = "position")
  
  plant_freq_position_s$s[is.nan(plant_freq_position_s$s)] <- 1.0
  
  # Indirect interaction data
  motif_interactions <- read_csv("Data/Data_processing/Motifs_connections/motif_interactions.csv") %>%
    select(position,indirect_interactions) %>% unique()
  
  plant_freq_position_s_ind <- plant_freq_position_s %>% 
    left_join(motif_interactions, by = "position") %>% 
    filter(!Node_FG %in% c("Birds","Lizards","Other_insects")) %>%
    rename(Percentil=percentil_its_GF,Specificity=s,Indirect_interactions=indirect_interactions,Group=Node_FG)
  
  
  plant_freq_position_s_ind$position <- as.factor(plant_freq_position_s_ind$position)
  
  model_freq_specif_ind_plant <- glmmTMB(Percentil ~ scale(Specificity)*Group +
                                           scale(Indirect_interactions)*Group + 
                                           (1|study_id/position),
                                         plant_freq_position_s_ind)
  
  # Plant Model REAL OUTPUTS!!!
  
  plants_sp_plot <- visreg(model_freq_specif_ind_plant, "Specificity", 
                           by="Group", scale = "response",
                           strip.names=c("Selfing herbs",
                                         "Small outcrossing perennials",
                                         "Self-incompatible perennials\nwith large flowers",
                                         "Tall plants with small\nunisexual flowers",
                                         "Short-lived outcrossers with\nlong zygomorphic flowers"),gg = TRUE)+
    theme_bw()
  
  plants_ind_int_plot <- visreg(model_freq_specif_ind_plant, "Indirect_interactions", 
                                by="Group", scale = "response",
                                strip.names=c("Selfing herbs",
                                              "Small outcrossing perennials",
                                              "Self-incompatible perennials\nwith large flowers",
                                              "Tall plants with small\nunisexual flowers",
                                              "Short-lived outcrossers with\nlong zygomorphic flowers"),gg = TRUE)+
    theme_bw()+xlab("Indirect interactions")
  

  
  
  ##############################################################
  # EXTRACT DATA FROM POLLINATOR'S
  ##############################################################
  
  df_from_plants_sp_plot <- ggplot_build(plants_sp_plot)
  df_from_plants_int_plot <- ggplot_build(plants_ind_int_plot)
  
  # Extract data from plants' specificity plot--------
  
  shade_data_x_plant_sp_plot <- df_from_plants_sp_plot$data[[1]] %>% select(x) %>% unique() %>% pull()
  
  for (x_i in 1:length(shade_data_x_plant_sp_plot)) {
    shade_data_aux_plant_sp_plot <- df_from_plants_sp_plot$data[[1]] %>% filter(x == shade_data_x_plant_sp_plot[x_i]) %>% 
      group_by(PANEL) %>% summarise(ymin = min(y), ymax = max(y)) %>%
      mutate(x = shade_data_x_plant_sp_plot[x_i],
             repetition = rep_i)
    shade_data_plant_sp_plot <- bind_rows(shade_data_plant_sp_plot,
                                          shade_data_aux_plant_sp_plot)
  }
  
  mean_line_data_plant_sp_plot <- bind_rows(mean_line_data_plant_sp_plot,
                                            df_from_plants_sp_plot$data[[2]] %>%
                                              mutate(repetition = rep_i))
  
  # Extract data from plant's ind. interactions' plot--------
  
  
  shade_data_x_plant_int_plot <- df_from_plants_int_plot$data[[1]] %>% select(x) %>% unique() %>% pull()
  
  for (x_i in 1:length(shade_data_x_plant_int_plot)) {
    shade_data_aux_plant_int_plot <- df_from_plants_int_plot$data[[1]] %>% filter(x == shade_data_x_plant_int_plot[x_i]) %>% 
      group_by(PANEL) %>% summarise(ymin = min(y), ymax = max(y)) %>%
      mutate(x = shade_data_x_plant_int_plot[x_i],
             repetition = rep_i)
    shade_data_plant_int_plot <- bind_rows(shade_data_plant_int_plot,
                                           shade_data_aux_plant_int_plot)
  }
  
  mean_line_data_plant_int_plot <- bind_rows(mean_line_data_plant_int_plot,
                                             df_from_plants_int_plot$data[[2]] %>%
                                               mutate(repetition = rep_i))
  
  file_suffix <- paste0("pollinator_sp_plot_",repetitions,"_repetitions_",
                        round(100*percentage_studies),"_percent.csv")
  
  write_csv(shade_data_pollinator_sp_plot,
            paste0("Data/Csv/shade_data_",file_suffix))
  write_csv(mean_line_data_pollinator_sp_plot,
            paste0("Data/Csv/mean_line_data_",file_suffix))
  
  file_suffix <- paste0("pollinator_int_plot_",repetitions,"_repetitions_",
                        round(100*percentage_studies),"_percent.csv")
  
  write_csv(shade_data_pollinator_int_plot,
            paste0("Data/Csv/shade_data_",file_suffix))
  write_csv(mean_line_data_pollinator_int_plot,
            paste0("Data/Csv/mean_line_data_",file_suffix))
  
  file_suffix <- paste0("plant_sp_plot_",repetitions,"_repetitions_",
                        round(100*percentage_studies),"_percent.csv")
  
  write_csv(shade_data_plant_sp_plot,
            paste0("Data/Csv/shade_data_",file_suffix))
  write_csv(mean_line_data_plant_sp_plot,
            paste0("Data/Csv/mean_line_data_",file_suffix))
 
  file_suffix <- paste0("plant_int_plot_",repetitions,"_repetitions_",
                        round(100*percentage_studies),"_percent.csv")
  
  write_csv(shade_data_plant_int_plot,
            paste0("Data/Csv/shade_data_",file_suffix))
  write_csv(mean_line_data_plant_int_plot,
            paste0("Data/Csv/mean_line_data_",file_suffix))

  
}


#######################################
# PLOT GRAPHS
#######################################

library(tidyverse)

# Load data for all datasets------------------------

percentage_studies_load <- 1.0
repetitions_load <- 1

file_suffix <- paste0("pollinator_sp_plot_",repetitions_load,"_repetitions_",
                      round(100*percentage_studies_load),"_percent.csv")

shade_data_pollinator_sp_plot_all <- read_csv(paste0("Data/Csv/shade_data_",file_suffix))
mean_line_data_pollinator_sp_plot_all <- read_csv(paste0("Data/Csv/mean_line_data_",file_suffix))

file_suffix <- paste0("pollinator_int_plot_",repetitions_load,"_repetitions_",
                      round(100*percentage_studies_load),"_percent.csv")

shade_data_pollinator_int_plot_all <- read_csv(paste0("Data/Csv/shade_data_",file_suffix))
mean_line_data_pollinator_int_plot_all <- read_csv(paste0("Data/Csv/mean_line_data_",file_suffix))

file_suffix <- paste0("plant_sp_plot_",repetitions_load,"_repetitions_",
                      round(100*percentage_studies_load),"_percent.csv")

shade_data_plant_sp_plot_all <- read_csv(paste0("Data/Csv/shade_data_",file_suffix))
mean_line_data_plant_sp_plot_all <- read_csv(paste0("Data/Csv/mean_line_data_",file_suffix))

file_suffix <- paste0("plant_int_plot_",repetitions_load,"_repetitions_",
                      round(100*percentage_studies_load),"_percent.csv")

shade_data_plant_int_plot_all <- read_csv(paste0("Data/Csv/shade_data_",file_suffix))
mean_line_data_plant_int_plot_all <- read_csv(paste0("Data/Csv/mean_line_data_",file_suffix))

# Load data for a fraction of the total datasets-----------------

percentage_studies_load <- .60
repetitions_load <- 1000

file_suffix <- paste0("pollinator_sp_plot_",repetitions_load,"_repetitions_",
                      round(100*percentage_studies_load),"_percent.csv")

shade_data_pollinator_sp_plot <- read_csv(paste0("Data/Csv/shade_data_",file_suffix))
mean_line_data_pollinator_sp_plot <- read_csv(paste0("Data/Csv/mean_line_data_",file_suffix))

file_suffix <- paste0("pollinator_int_plot_",repetitions_load,"_repetitions_",
                      round(100*percentage_studies_load),"_percent.csv")

shade_data_pollinator_int_plot <- read_csv(paste0("Data/Csv/shade_data_",file_suffix))
mean_line_data_pollinator_int_plot <- read_csv(paste0("Data/Csv/mean_line_data_",file_suffix))

file_suffix <- paste0("plant_sp_plot_",repetitions_load,"_repetitions_",
                      round(100*percentage_studies_load),"_percent.csv")

shade_data_plant_sp_plot <- read_csv(paste0("Data/Csv/shade_data_",file_suffix))
mean_line_data_plant_sp_plot <- read_csv(paste0("Data/Csv/mean_line_data_",file_suffix))

file_suffix <- paste0("plant_int_plot_",repetitions_load,"_repetitions_",
                      round(100*percentage_studies_load),"_percent.csv")

shade_data_plant_int_plot <- read_csv(paste0("Data/Csv/shade_data_",file_suffix))
mean_line_data_plant_int_plot <- read_csv(paste0("Data/Csv/mean_line_data_",file_suffix))

pollinator_names <- list(
  '1'="Bee",
  '2'="Coleoptera",
  '3'="Lepidoptera",
  '4'="Non-bee-Hymenoptera",
  '5'="Non-syrphids-diptera",
  '6'="Syrphids"
)

plant_names <- list(
  '1'="Selfing herbs",
  '2'="Small outcrossing perennials",
  '3'="Self-incompatible perennials\nwith large flowers",
  '4'="Tall plants with small\nunisexual flowers",
  '5'="Short-lived outcrossers with\nlong zygomorphic flowers")

pollinator_labeller <- function(variable,value){
  return(pollinator_names[value])
}

plant_labeller <- function(variable,value){
  return(plant_names[value])
}


final_pollinator_sp_plot <- ggplot(shade_data_pollinator_sp_plot)+
  geom_ribbon(aes(x=x, ymax=ymax, ymin= ymin, fil = as.factor(repetition)), alpha=10/repetitions_load)+
  facet_grid( ~ PANEL, labeller=pollinator_labeller)+
  geom_line(data = mean_line_data_pollinator_sp_plot, 
            aes(x = x, y = y, co = as.factor(repetition)), alpha=0.2)+
  geom_line(data = shade_data_pollinator_sp_plot_all,
              aes(x=x, y=ymax, color ="red"), size = 1.2,linetype = "dashed")+
  geom_line(data = shade_data_pollinator_sp_plot_all,
            aes(x=x, y= ymin, color ="red"), size = 1.2, linetype = "dashed")+
  geom_line(data = mean_line_data_pollinator_sp_plot_all, 
            aes(x = x, y = y, color ="red"), size = 1.2)+
  theme_bw()+
  ggtitle("Floral visitors") +
  labs(y="Percentile", x =NULL)+
  theme(legend.position = "none")


final_pollinator_int_plot <- ggplot(shade_data_pollinator_int_plot)+
  geom_ribbon(aes(x=x, ymax=ymax, ymin= ymin, fil = as.factor(repetition)), alpha=10/repetitions_load)+
  facet_grid( ~ PANEL, labeller=pollinator_labeller)+
  geom_line(data = mean_line_data_pollinator_int_plot, 
            aes(x = x, y = y, co = as.factor(repetition)), alpha=0.2)+
  geom_line(data = shade_data_pollinator_int_plot_all,
            aes(x=x, y=ymax, color ="red"), size = 1.2,linetype = "dashed")+
  geom_line(data = shade_data_pollinator_int_plot_all,
            aes(x=x, y= ymin, color ="red"), size = 1.2, linetype = "dashed")+
  geom_line(data = mean_line_data_pollinator_int_plot_all, 
            aes(x = x, y = y, color ="red"), size = 1.2)+
  theme_bw()+
  ggtitle("Floral visitors") +
  labs(y="Percentile", x =NULL)+
  theme(legend.position = "none")


final_plant_sp_plot <- ggplot(shade_data_plant_sp_plot)+
  geom_ribbon(aes(x=x, ymax=ymax, ymin= ymin, fil = as.factor(repetition)), alpha=10/repetitions_load)+
  facet_grid( ~ PANEL, labeller=plant_labeller)+
  geom_line(data = mean_line_data_plant_sp_plot, 
            aes(x = x, y = y, co = as.factor(repetition)), alpha=0.2)+
  geom_line(data = shade_data_plant_sp_plot_all,
            aes(x=x, y=ymax, color ="red"), size = 1.2,linetype = "dashed")+
  geom_line(data = shade_data_plant_sp_plot_all,
            aes(x=x, y= ymin, color ="red"), size = 1.2, linetype = "dashed")+
  geom_line(data = mean_line_data_plant_sp_plot_all, 
            aes(x = x, y = y, color ="red"), size = 1.2)+
  theme_bw()+
  ggtitle("Plants")+
  labs(x = "Specialization", y="Percentile")+
  theme(legend.position = "none")


final_plant_int_plot <- ggplot(shade_data_plant_int_plot)+
  geom_ribbon(aes(x=x, ymax=ymax, ymin= ymin, fil = as.factor(repetition)), alpha=10/repetitions_load)+
  facet_grid( ~ PANEL, labeller=plant_labeller)+
  geom_line(data = mean_line_data_plant_int_plot, 
            aes(x = x, y = y, co = as.factor(repetition)), alpha=0.2)+
  geom_line(data = shade_data_plant_int_plot_all,
            aes(x=x, y=ymax, color ="red"), size = 1.2,linetype = "dashed")+
  geom_line(data = shade_data_plant_int_plot_all,
            aes(x=x, y= ymin, color ="red"), size = 1.2, linetype = "dashed")+
  geom_line(data = mean_line_data_plant_int_plot_all, 
            aes(x = x, y = y, color ="red"), size = 1.2)+
  theme_bw()+
  ggtitle("Plants")+
  labs(y="Percentile", x ="Number of indirect interactions")+
  theme(legend.position = "none")

library(patchwork)

final_pollinator_sp_plot/final_plant_sp_plot & theme(axis.text=element_text(size=10),
                                                    axis.title=element_text(size=16)) & 
  plot_annotation(tag_levels = c('A'),  tag_suffix = ')') &  theme(plot.tag = element_text(size = 22))

final_pollinator_int_plot/final_plant_int_plot & theme(axis.text=element_text(size=10),
                                                       axis.title=element_text(size=16)) & 
  plot_annotation(tag_levels = c('A'),  tag_suffix = ')') &  theme(plot.tag = element_text(size = 22))
