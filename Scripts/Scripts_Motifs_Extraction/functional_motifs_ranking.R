
# Load libraries
library(tidyverse)

# Load raw data with plant and pollinators
networks_raw <- read_csv("Data/Csv/data_for_motifs_analysis.csv")

# Load raw data with plant and pollinators
frequency_raw <- read_csv("Data/Csv/network_frequency_motifs.csv")

# List the files with summary information about the amount of functional motifs
folder_motif_data <- "Data/Csv/Motifs links"

motif_files <- list.files(folder_motif_data) 
motif_summary_files <-  motif_files[grep("Summary_functional_motifs_",motif_files)]

all_networks_scores <- NULL

# Extract positions and save the results
for (i in 1:length(motif_summary_files)){
  
  file_i <- paste0(folder_motif_data, "/", motif_summary_files[i])
  summary_i <- read_csv(file_i)
  
  # Add number of nodes in each type of motif
  summary_i$nodes <- 2
  summary_i$nodes[summary_i$Motif_pattern_id %in% c(2:3)] <- 3
  summary_i$nodes[summary_i$Motif_pattern_id %in% c(4:7)] <- 4
  summary_i$nodes[summary_i$Motif_pattern_id %in% c(8:17)] <- 5
  
  # Add number of plant nodes in each type of motif
  summary_i$Plant_nodes <- 4
  summary_i$Plant_nodes[summary_i$Motif_pattern_id %in% 
                          c(1,3,4,8)] <- 1
  summary_i$Plant_nodes[summary_i$Motif_pattern_id %in% 
                          c(2,5,6,9,10,11,12)] <- 2
  summary_i$Plant_nodes[summary_i$Motif_pattern_id %in% 
                          c(7,13,14,15,16)] <- 3
  
  # Add number of pollinator nodes in each type of motif
  summary_i$Pollinator_nodes <- 4
  summary_i$Pollinator_nodes[summary_i$Motif_pattern_id %in% 
                               c(1,2,7,17)] <- 1
  summary_i$Pollinator_nodes[summary_i$Motif_pattern_id %in% 
                               c(3,5,6,13,14,15,16)] <- 2
  summary_i$Pollinator_nodes[summary_i$Motif_pattern_id %in% 
                               c(4,9,10,11,12)] <- 3
  
  # Number of plant/pollinator species
  species_i <- networks_raw %>%
    filter(Network_id == unique(summary_i$Network_id)) %>%
    dplyr::select(Pollinator_species,Plant_species)
  
  summary_i$Total_plants <- unique(species_i$Plant_species) %>% length()
  summary_i$Total_pollinators <- unique(species_i$Pollinator_species) %>% length()
    
  
  # Calculate normalize_sum
  summary_i$Total_motifs <- sum(summary_i$Counts)
  summary_i <- summary_i %>% mutate(normalize_sum = Counts/Total_motifs)
  
  # Calculate normalize_sizeclass
  frequency_i <- frequency_raw %>% 
    filter(Network_id == unique(summary_i$Network_id)) %>%
    dplyr::select(nodes,frequency) %>% group_by(nodes) %>% 
    count(wt = frequency) %>% rename(Total_sizeclass = n)
  
  summary_i <- summary_i %>% left_join(frequency_i, by = "nodes") %>%
    mutate(normalize_sizeclass = Counts/Total_sizeclass)
  
  # Calculate normalize_nodesets
  summary_i <- summary_i %>%
    mutate(Total_nodesets = choose(Total_pollinators,Pollinator_nodes)*
             choose(Total_plants,Plant_nodes),
           normalize_nodesets = Counts/Total_nodesets)
  
  
  all_networks_scores <- bind_rows(all_networks_scores,summary_i)
  
}

# Save percentages of functional motifs

write_csv(all_networks_scores,"Data/Csv/functional_motifs_percentages.csv")

# Extract mean values and sd for the previous percentages

functional_motif_results <- all_networks_scores %>%
  dplyr::select(Motif_pattern_id,Motif_functional_ID) %>% unique()

functional_motif_results$Mean_normalize_sum <- NA
functional_motif_results$Mean_normalize_sizeclass <- NA
functional_motif_results$Mean_normalize_nodesets <- NA
functional_motif_results$Sd_normalize_sum <- NA
functional_motif_results$Sd_normalize_sizeclass <- NA
functional_motif_results$Sd_normalize_nodesets <- NA

number_network <- all_networks_scores$Network_id %>% unique() %>% length()

for(i in 1:nrow(functional_motif_results)){
  
  aux_i <- all_networks_scores %>% 
    filter(Motif_pattern_id == Motif_pattern_id[i],
           Motif_functional_ID == Motif_functional_ID[i]) %>% 
    dplyr::select(normalize_sum,normalize_sizeclass,normalize_nodesets)
  
  # We add zeros for each network that does not have the functional group
  
  if (nrow(aux_i) < number_network){
    aux_i[c((nrow(aux_i)+1):number_network),] <- 0
  }
  
  functional_motif_results$Mean_normalize_sum[i] <- mean(aux_i$normalize_sum)
  functional_motif_results$Mean_normalize_sizeclass[i] <- mean(aux_i$normalize_sizeclass)
  functional_motif_results$Mean_normalize_nodesets[i] <- mean(aux_i$normalize_nodesets)
  functional_motif_results$Sd_normalize_sum[i] <- sd(aux_i$normalize_sum)
  functional_motif_results$Sd_normalize_sizeclass[i] <- sd(aux_i$normalize_sizeclass)
  functional_motif_results$Sd_normalize_nodesets[i] <- sd(aux_i$normalize_nodesets)
}

# Save percentages of functional motifs

write_csv(functional_motif_results,"Data/Csv/functional_motif_mean_values.csv")

# Prepare ranking plots
functional_motif_results_label <- functional_motif_results
functional_motif_results_label$Label <- paste0("Motif ",
                                               functional_motif_results_label$Motif_pattern_id,
                                               "\n",
                                               functional_motif_results_label$Motif_functional_ID)


x <- functional_motif_results_label %>% 
  arrange(desc(Mean_normalize_sum)) %>% 
  head(10) %>% select(Label,Motif_functional_ID,Mean_normalize_sum,Sd_normalize_sum)

x_gr <- x %>% gather(key="mean",value = "value",-c(Label,Motif_functional_ID))
x_gr$mean[x_gr$mean=="Mean_normalize_sum"] <- "Mean"
x_gr$mean[x_gr$mean=="Sd_normalize_sum"] <- "Sd"

ggplot(arrange(x_gr,mean))+
geom_col(aes(x = 100*value, y = reorder(Label,value),fill = mean),
         position = position_stack(reverse = TRUE))+
  labs(x="Average percentage (%)",
       y = NULL,
       fill = NULL,
       title ="Normalization: Total number of motifs")+
  theme_bw()

x2 <- functional_motif_results_label %>% 
  arrange(desc(Mean_normalize_sizeclass)) %>% 
  head(10) %>% select(Label,Motif_functional_ID,Mean_normalize_sizeclass,
                      Sd_normalize_sizeclass)

x2_gr <- x2 %>% gather(key="mean",value = "value",-c(Label,Motif_functional_ID))
x2_gr$mean[x2_gr$mean=="Mean_normalize_sizeclass"] <- "Mean"
x2_gr$mean[x2_gr$mean=="Sd_normalize_sizeclass"] <- "Sd"

ggplot(arrange(x2_gr,mean))+
  geom_col(aes(x = 100*value, y = reorder(Label,value),fill = mean),
           position = position_stack(reverse = TRUE))+
  labs(x="Average percentage (%)",
       y = NULL,
       fill = NULL,
       title ="Normalization: Total number of motifs in a given size class")+
  theme_bw()



x3 <- functional_motif_results_label %>% 
  arrange(desc(Mean_normalize_nodesets)) %>% 
  head(10) %>% select(Label,Motif_functional_ID,Mean_normalize_nodesets,
                      Sd_normalize_nodesets)

x3_gr <- x3 %>% gather(key="mean",value = "value",-c(Label,Motif_functional_ID))
x3_gr$mean[x3_gr$mean=="Mean_normalize_nodesets"] <- "Mean"
x3_gr$mean[x3_gr$mean=="Sd_normalize_nodesets"] <- "Sd"

ggplot(arrange(x3_gr,mean))+
  geom_col(aes(x = 100*value, y = reorder(Label,value),fill = mean),
           position = position_stack(reverse = TRUE))+
  labs(x="Average percentage (%)",
       y = NULL,
       fill = NULL,
       title ="Normalization: Total number of node sets")+
  theme_bw()
