
# This script extracts the interactions of (up to 5 nodes) motifs and save them in the
# folder "Data/Csv/Motifs links"

# Load libraries
library(tidyverse)
library(igraph)
source("Scripts/Scripts_Motifs_Extraction/aux_functions.R")

# CLEANING----------
# Load data networks raw data and select those interactions whose Interaction frequency is 
# greater than a given threshold

int.threshold <- 1

networks <- read_csv("Data/Csv/data_for_motifs_analysis.csv") %>% select(-X1) %>%
  filter(Interaction >= int.threshold)

patterns <- read_csv("Data/Data_processing/Motifs_connections/motif_pattern_connections.csv")

# MOTIF MINING --------

# Get the list of networks and motif patterns that are available to iterate
list_Network_id <- networks$Network_id %>% unique()
list_parterns <- patterns$motif_id %>% unique()

# Select networks with fewer motifs

min_num_motifs <- 55001
max_num_motifs <- 110000

small_network_motifs <- read_csv("Data/Csv/network_frequency_motifs.csv") %>% filter(nodes<=5) %>%
  group_by(Network_id) %>% count(wt = frequency) %>% filter(n >= min_num_motifs,
                                                            n <= max_num_motifs) %>%
  dplyr::select(Network_id) %>% pull()

small_network_index <- which(list_Network_id %in% small_network_motifs)

# Extract motifs' links for selected networks

for (i.network in small_network_index[2:length(small_network_index)]){#1:length(list_Network_id)){
  
  start_time <- Sys.time()
  print(list_Network_id[i.network])
  print(i.network)
  
  # Initialize variables for storing purposes
  motifs_connections <- NULL # Store the connections of each motif
  df_motif_8_17 <- NULL # Auxiliary variable that storages info about patterns 8 and 17
  
  
  # Create a graph for the i-th network
  networks_i <- networks %>% filter(Network_id == list_Network_id[i.network])
  edge_list_i <- networks_i[,c("Pollinator_species","Plant_species")]
  
  g_i <- igraph::graph_from_edgelist(as.matrix(edge_list_i), directed = F)
  
  # To create a two mode network we define the types
  V(g_i)$type <- bipartite_mapping(g_i)$type  ## Add the "type" attribute to the network.
  
  
  for(i.pattern_index in 1:length(list_parterns)){
    
    i.pattern <- list_parterns[i.pattern_index]
    
    if(i.pattern == 8){
      
      results_8_17 <- connections_pattern_i(motifs_connections,patterns,
                                            list_parterns,i.pattern,df_motif_8_17)
      
      motifs_connections <- results_8_17[[1]]
      df_motif_8_17 <- results_8_17[[2]]
      
    }else{
      
      motifs_connections <- connections_pattern_i(motifs_connections,patterns,
                                                  list_parterns,i.pattern,df_motif_8_17)
    }
    
    
  }
  
  # Save results in folder
  folder_motifs <- "Data/Csv/Motifs links"
  file_motifs <- paste0(folder_motifs,"/Motifs_links_",list_Network_id[i.network],".csv") 
  write_csv(motifs_connections,file_motifs)
  
  #Print time consumed
  end_time <- Sys.time()
  print(end_time-start_time)
}
