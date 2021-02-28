
# This script generate a file with the motif positions of those nodes in 
# "Data/Csv/Motifs" files (which contain the links of each motif that is
# present in a given network).
# The resulting full info is also storaged at "Data/Csv/Motifs": "Motifs_positions_" files
# A summary with the counts for each functional motif is storaged at "Data/Csv/Motifs":
# "Summary_functional_motifs_" files

# Load libraries
library(tidyverse)

# Load raw data with functional groups
networks_raw <- read_csv("Data/Csv/data_for_motifs_analysis.csv")

# List the files with information about the links of the motifs
folder_motif_data <- "Data/Csv/Motifs links"

motif_files <- list.files(folder_motif_data) 
motif_link_files <-  motif_files[grep("Motifs_links_",motif_files)]

# Extract positions and save the results

for (i in 1:length(motif_link_files)){

  
  start_time <- Sys.time()
  print(motif_link_files[i])
  
  # Open link file and rearrange it
  file_i <- paste0(folder_motif_data, "/", motif_link_files[i])
  motif_links_i <- read_csv(file_i) %>% 
    gather("Type","Species",-c(Network_id, Motif_pattern_id,
                               Nodes, Motif_number)) %>%
    arrange(Motif_pattern_id, Nodes,Motif_number) %>%
    group_by(Network_id, Motif_pattern_id,
             Nodes,Motif_number, Type,Species) %>% count() %>%
    rename(Degree = n,Node_id = Species)
  
  # Classify nodes: Plant - Pollinator
  motif_links_i$Type[motif_links_i$Type == "Plant_species"] <- "Plant"
  motif_links_i$Type[motif_links_i$Type == "Pollinator_species"] <- "Pollinator"
  
  # Assign position depending on the motif pattern, type and degree
  motif_links_i$Position <- NA
  
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 1) &
                  (motif_links_i$Type == "Plant")] <- 1
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 1) &
                           (motif_links_i$Type == "Pollinator")] <- 2
  
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 2) &
                           (motif_links_i$Type == "Plant")] <- 3
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 2) &
                           (motif_links_i$Type == "Pollinator")] <- 4
  
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 3) &
                           (motif_links_i$Type == "Plant")] <- 5
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 3) &
                           (motif_links_i$Type == "Pollinator")] <- 6
  
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 4) &
                           (motif_links_i$Type == "Plant")] <- 7
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 4) &
                           (motif_links_i$Type == "Pollinator")] <- 8
  
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 5) &
                           (motif_links_i$Type == "Plant") &
                           (motif_links_i$Degree == 1)] <- 9
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 5) &
                           (motif_links_i$Type == "Plant") &
                           (motif_links_i$Degree == 2)] <- 10
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 5) &
                           (motif_links_i$Type == "Pollinator")&
                           (motif_links_i$Degree == 1)] <- 11
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 5) &
                           (motif_links_i$Type == "Pollinator")&
                           (motif_links_i$Degree == 2)] <- 12
  
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 6) &
                           (motif_links_i$Type == "Plant")] <- 13
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 6) &
                           (motif_links_i$Type == "Pollinator")] <- 14
  
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 7) &
                           (motif_links_i$Type == "Plant")] <- 15
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 7) &
                           (motif_links_i$Type == "Pollinator")] <- 16
  
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 8) &
                           (motif_links_i$Type == "Plant")] <- 17
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 8) &
                           (motif_links_i$Type == "Pollinator")] <- 18
  
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 9) &
                           (motif_links_i$Type == "Plant") &
                           (motif_links_i$Degree == 1)] <- 19
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 9) &
                           (motif_links_i$Type == "Plant") &
                           (motif_links_i$Degree == 3)] <- 20
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 9) &
                           (motif_links_i$Type == "Pollinator")&
                           (motif_links_i$Degree == 1)] <- 21
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 9) &
                           (motif_links_i$Type == "Pollinator")&
                           (motif_links_i$Degree == 2)] <- 22
  
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 10) &
                           (motif_links_i$Type == "Plant")] <- 23
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 10) &
                           (motif_links_i$Type == "Pollinator")&
                           (motif_links_i$Degree == 1)] <- 24
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 10) &
                           (motif_links_i$Type == "Pollinator")&
                           (motif_links_i$Degree == 2)] <- 25
  
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 11) &
                           (motif_links_i$Type == "Plant") &
                           (motif_links_i$Degree == 2)] <- 26
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 11) &
                           (motif_links_i$Type == "Plant") &
                           (motif_links_i$Degree == 3)] <- 27
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 11) &
                           (motif_links_i$Type == "Pollinator")&
                           (motif_links_i$Degree == 1)] <- 28
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 11) &
                           (motif_links_i$Type == "Pollinator")&
                           (motif_links_i$Degree == 2)] <- 29
  
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 12) &
                           (motif_links_i$Type == "Plant")] <- 30
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 12) &
                           (motif_links_i$Type == "Pollinator")] <- 31
  
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 13) &
                           (motif_links_i$Type == "Plant") &
                           (motif_links_i$Degree == 1)] <- 32
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 13) &
                           (motif_links_i$Type == "Plant") &
                           (motif_links_i$Degree == 2)] <- 33
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 13) &
                           (motif_links_i$Type == "Pollinator")&
                           (motif_links_i$Degree == 1)] <- 34
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 13) &
                           (motif_links_i$Type == "Pollinator")&
                           (motif_links_i$Degree == 3)] <- 35
  
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 14) &
                           (motif_links_i$Type == "Plant") &
                           (motif_links_i$Degree == 1)] <- 36
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 14) &
                           (motif_links_i$Type == "Plant") &
                           (motif_links_i$Degree == 2)] <- 37
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 14) &
                           (motif_links_i$Type == "Pollinator")] <- 38
  
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 15) &
                           (motif_links_i$Type == "Plant") &
                           (motif_links_i$Degree == 1)] <- 39
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 15) &
                           (motif_links_i$Type == "Plant") &
                           (motif_links_i$Degree == 2)] <- 40
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 15) &
                           (motif_links_i$Type == "Pollinator")&
                           (motif_links_i$Degree == 2)] <- 41
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 15) &
                           (motif_links_i$Type == "Pollinator")&
                           (motif_links_i$Degree == 3)] <- 42
  
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 16) &
                           (motif_links_i$Type == "Plant")] <- 43
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 16) &
                           (motif_links_i$Type == "Pollinator")] <- 44
  
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 17) &
                           (motif_links_i$Type == "Plant")] <- 45
  motif_links_i$Position[(motif_links_i$Motif_pattern_id == 17) &
                           (motif_links_i$Type == "Pollinator")] <- 46
  
  
  # Add functional groups
  network_i <- motif_links_i$Network_id %>% unique()
  
  pollinator_FG <- networks_raw %>% filter(Network_id == network_i) %>% 
    dplyr::select(Pollinator_species,Pollinator_functional_group) %>% unique() %>%
    rename(Node_id = Pollinator_species, Functional_group = Pollinator_functional_group)
  
  plant_FG <- networks_raw %>% filter(Network_id == network_i) %>% 
    dplyr::select(Plant_species,Plant_functional_groups) %>% unique() %>%
    rename(Node_id = Plant_species, Functional_group = Plant_functional_groups)
  plant_FG$Functional_group <- as.character(plant_FG$Functional_group)
  
  all_FG <- bind_rows(pollinator_FG, plant_FG)
    
  motif_links_i_FG <- motif_links_i %>% left_join(all_FG, by = "Node_id")
  
  # Rearrange data frame and display equal node positions alphabetically
  
  motif_links_i_FG <- motif_links_i_FG %>% 
    arrange(Motif_pattern_id,Motif_number,Position,Node_id)
  
  # Add motif functional ID
  
  motif_links_i_FG$Motif_functional_ID <- NA
  
  amount_motifs_patterns  <- motif_links_i_FG %>% ungroup() %>%
    dplyr::select(Motif_pattern_id,Motif_number) %>% unique() %>% 
    group_by(Motif_pattern_id) %>% count()
  
  motif_patterns <- amount_motifs_patterns$Motif_pattern_id
  
  for(i.pattern in 1:length(motif_patterns)){
    
    for(j.pattern in 1:amount_motifs_patterns$n[i.pattern]){
      
      positions_FG <- motif_links_i_FG %>% ungroup() %>%
        filter(Motif_pattern_id == motif_patterns[i.pattern],
               Motif_number ==  j.pattern) %>% select(Functional_group) %>% pull()
      
      motif_links_i_FG$Motif_functional_ID[(motif_links_i_FG$Motif_pattern_id == 
                                              motif_patterns[i.pattern]) &
                                           (motif_links_i_FG$Motif_number ==
                                              j.pattern)] <- paste(positions_FG,
                                                                   collapse = "_")
      
    }

  }
  
  # Create a summary file
  
  summary_i <- motif_links_i_FG %>% ungroup() %>% 
    dplyr::select(Network_id,Motif_pattern_id,Motif_number,Motif_functional_ID) %>%
    unique() %>% group_by(Network_id,Motif_pattern_id,Motif_functional_ID) %>% count() %>%
    rename(Counts = n) %>% arrange(desc(Counts))
  
  # save file with the full info of motifs positions and the summary of functional motifs
  
  new_file_i <- paste0(folder_motif_data, "/Motifs_positions_", network_i,".csv")
  write_csv(motif_links_i_FG,new_file_i)
  
  summary_file_i <- paste0(folder_motif_data, "/Summary_functional_motifs_", network_i, ".csv")
  write_csv(summary_i,summary_file_i)
  
  #Print time consumed
  end_time <- Sys.time()
  print(end_time-start_time)
  
}
