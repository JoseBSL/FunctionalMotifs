
# This script cleans the raw data file for motif analyses. Then, we extract the following
# information for a given networks:
# 1) The frequency each motif types
# 2) The frequency of motif positions for each node in each network

# Load libraries
library(tidyverse)
library(bmotif)
library(igraph)

# CLEANING----------
# Load data networks raw data and select those interactions whose Interaction frequency is 
# greater than a given threshold

int.threshold <- 1

networks <- read_csv("Data/Csv/data_for_motifs_analysis.csv") %>% select(-X1) %>%
  filter(Interaction >= int.threshold)

# MOTIF MINING --------


list_Network_id <- networks$Network_id %>% unique()

freq_motifs <- NULL #Store the frequency each motif types
node_positions <- NULL # Store the frequency of motif positions for each node

# We are not using weigths to estimate motif info
mean_weight_i <-  F
weights_method_i <-  'none'

# Get info

for (i.network in 1:length(list_Network_id)){
  
  networks_i <- networks %>% filter(Network_id == list_Network_id[i.network])
  edge_list_i <- networks_i[,c("Pollinator_species","Plant_species")]
  
  g_i <- igraph::graph_from_edgelist(as.matrix(edge_list_i), directed = FALSE)
  
  # To create a two mode network we define the types
  V(g_i)$type <- bipartite_mapping(g_i)$type  ## Add the "type" attribute to the network.
  
  # Get incidence matrix
  inc_matrix_i <- as_incidence_matrix(g_i, types = NULL, attr = NULL, names = TRUE,
                                      sparse = FALSE)
  
  # Frequency of different motifs --------
  
  # Since we will compare multiple networks, we normalise motif frequencies.
  
  # ‘normalise_sum’ expresses the frequency of each motif as a proportion of the total number
  # of motifs in the network.
  
  # ‘normalise_sizeclass’ expresses the frequency of each motif as a proportion of the total 
  # number of motifs within its size class (i.e., 2 nodes, 3 nodes...).
  
  # ‘normalise_nodesets’ expresses the frequency of each motif as the number of species
  # combinations that occur in a motif as a proportion of the number of species combinations that
  # could occur in that motif.
  
  freq_motifs_i <- mcount(inc_matrix_i, normalisation = T,
                          mean_weight = mean_weight_i, standard_dev = F)
  
  # Nodes position in motifs --------
  # Check Normalisation information
  node_positions_i <- node_positions(inc_matrix_i, level = "all",
                                     weights_method = weights_method_i)
  
  
  # Add network id
  
  col_names_freq_i <- colnames(freq_motifs_i)
  col_names_posit_i <- colnames(node_positions_i)
  
  freq_motifs_i$Network_id  <- list_Network_id[i.network]
  node_positions_i$Network_id  <- list_Network_id[i.network]
  
  # Rearrange columns
  
  freq_motifs_i <- freq_motifs_i[,c("Network_id",col_names_freq_i)]
  node_positions_i <- node_positions_i[,c("Network_id",col_names_posit_i)]
  
  
  # Store results----
  
  freq_motifs <-bind_rows(freq_motifs,freq_motifs_i)
  node_positions <- bind_rows(node_positions,node_positions_i)
  
}

# Save results
write_csv(freq_motifs,"Data/Csv/network_frequency_motifs.csv")
write_csv(node_positions,"Data/Csv/node_frequency_positions_motifs.csv")
