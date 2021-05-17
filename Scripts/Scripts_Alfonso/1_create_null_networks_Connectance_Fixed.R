
# Extract motif positions for the nodes in the networks in
# "Data/Csv/data_for_motifs_analysis_1.csv"
# To do so, the script calculates the incidence matrices of
# the BINARY networks

# Besides, the script storage the position of 1000 null models where
# the connectance is conserved. We used 'vaznull', which is described
# as The algorithm was described as follows: "The algorithm randomized 
# the total number of individual interactions observed in the original 
# interaction matrix, F. To this end, the algorithm first created a 
# binary matrix, assigning interspecific interactions according to 
# species-specific probabilities, requiring that each species had at 
# least one interaction. As in Vazquez et al. (2005b), the 
# species-specific probabilities were proportional to species' relative 
# abundances (probabilities are in fact approximately proportional and 
# not equal to relative abundances because of the requirement that each 
# species receives at least one interaction; this requirement causes 
# probabilities to deviate from relative abundances, especially for rare 
# species). Once the number of filled cells in the original matrix was 
# reached, the remaining interactions were distributed among the filled 
# cells, so that connectance in the original and randomized matrices was 
# the same." (Vazquez et al. 2007, page 1122-1123).

# Load libraries
library(tidyverse)
library(igraph)
library(bmotif)


# CLEANING----------
# Load data networks raw data and select those interactions whose Interaction frequency is 
# greater than a given threshold
  
int.threshold <- 1

networks <- read_csv("Data/Csv/data_for_motifs_analysis_1.csv") %>% select(-X1) %>%
  filter(Interaction >= int.threshold)

# MOTIF MINING --------

# Get the list of networks and motif patterns that are available to iterate
list_Network_id <- networks$Network_id %>% unique()


# CONFIGURATION OF BMOTIF PACKAGE----
# We are not using weigths to estimate motif info
mean_weight_i <-  F
weights_method_i <-  'none'
normalisation_method_i = "none" #Here we dont use normalization. We will normalize once we
# sum the frequecies of the nodes that belong to a given functional group by sizeclass
# (in 2_processing_positions_null_Connectance_Fixed.R)

# sizeclass: divides the position measure for each node by the total number of times that 
# node appears in any position within the same motif size class (the number of nodes a motif 
# contains).
# "sizeclass_NAzero": same as ’sizeclass’ but replaces all NA values with 0. If a species
# does not occur in any motifs in a given size class, ’sizeclass’ normalisation will return
# NAs. ’sizeclass_NAzero’ avoids this by replacing NAs with zero.


six_node_conf = F # only up to 5 nodes motifs

set.seed(123)
for (i.network in 1:length(list_Network_id)){
  
  print(list_Network_id[i.network])
  print(i.network)
  
  # Create a graph for the i-th network
  networks_i <- networks %>% filter(Network_id == list_Network_id[i.network])
  
  networks_i <- networks_i %>% 
    group_by(Pollinator_species,Plant_species,
             Pollinator_functional_group,Plant_functional_groups) %>%
    count(wt = Interaction) %>% rename(Interaction = n)
  
  edge_list_i <- networks_i[,c("Pollinator_species","Plant_species")]
  
  g_i <- igraph::graph_from_edgelist(as.matrix(edge_list_i), directed = F)
  
  # To create a two mode network we define the types
  V(g_i)$type <- bipartite_mapping(g_i)$type  ## Add the "type" attribute to the network.
  
  # Create incidence matrix
  incidence_mat_i <- as_incidence_matrix(g_i)
  
  # # Add weights to the incidence matrix
  # 
  # el <- cbind(a=networks_i$Pollinator_species, 
  #             b=networks_i$Plant_species)
  # 
  # incidence_mat_i[el[,1:2]] <- networks_i$Interaction
  
  # Extract node positions
  node_positions_i <- bmotif::node_positions(incidence_mat_i, six_node = six_node_conf,
                                             level = "all", weights_method = weights_method_i,
                                             normalisation = normalisation_method_i)
  
  # Extract motifs' frequencies
  freq_motifs_i <- mcount(incidence_mat_i, six_node = six_node_conf, normalisation = T,
                               mean_weight = mean_weight_i, standard_dev = F)
  
  Node_id  <- row.names(node_positions_i)
  
  # Add network id
  col_names_posit_i <- colnames(node_positions_i)
  col_names_freq_i <- colnames(freq_motifs_i)
  
  node_positions_i$Network_id  <- list_Network_id[i.network]
  freq_motifs_i$Network_id  <- list_Network_id[i.network]
  
  node_positions_i$Node_id  <- Node_id
  node_positions_i$Type_network  <- "observed"
  freq_motifs_i$Type_network  <- "observed"
  
  
  # Add FG
  
  species_FG <- networks_i %>% ungroup() %>% 
    dplyr::select(Pollinator_species,Plant_species, Pollinator_functional_group,
                  Plant_functional_groups)
  
  list_FG <- tibble(Node_id = c(species_FG$Pollinator_species,
                                species_FG$Plant_species),
                    Node_FG = c(species_FG$Pollinator_functional_group,
                                species_FG$Plant_functional_groups)) %>%
    distinct()
  
  
  node_positions_i <- node_positions_i %>% 
    left_join(list_FG, by = "Node_id")
  
  # Rearrange columns
  node_positions_i <- node_positions_i[,c("Network_id","Type_network",
                                          "Node_id","Node_FG",
                                          col_names_posit_i)]
  
  freq_motifs_i <- freq_motifs_i[,c("Network_id","Type_network",col_names_freq_i)]
  
  # Generate null models for a given incidence matrix------
  
  null_mod <- bipartite::nullmodel(incidence_mat_i, N = 1000, method = "vaznull")
  
  for(i.null in 1:length(null_mod)){
    
    if (sum(null_mod[[i.null]] > 1)){print("ALERT!! NON-BINARY MATRIX")}
    
    rownames(null_mod[[i.null]]) <- rownames(incidence_mat_i)
    colnames(null_mod[[i.null]]) <- colnames(incidence_mat_i)
    
    # Extract node positions
    node_positions_i_null <- bmotif::node_positions(null_mod[[i.null]], six_node = six_node_conf,
                                                    level = "all", 
                                                    weights_method = weights_method_i,
                                                    normalisation = normalisation_method_i)
    
    # Extract motifs' frequencies
    freq_motifs_i_null <- mcount(null_mod[[i.null]], six_node = six_node_conf, normalisation = T,
                            mean_weight = mean_weight_i, standard_dev = F)
    
    
    # Add network id
    col_names_posit_i <- colnames(node_positions_i_null)
    col_names_freq_i <- colnames(freq_motifs_i_null)
    
    node_positions_i_null$Network_id  <- list_Network_id[i.network]
    freq_motifs_i_null$Network_id  <- list_Network_id[i.network]
    
    node_positions_i_null$Node_id  <- Node_id
    node_positions_i_null$Type_network  <- paste0("random ",i.null)
    freq_motifs_i_null$Type_network  <- paste0("random ",i.null)
    
    # Add FG
    
    node_positions_i_null <- node_positions_i_null %>% 
      left_join(list_FG, by = "Node_id")
    
    # Rearrange columns
    node_positions_i_null <- node_positions_i_null[,c("Network_id","Type_network",
                                                      "Node_id","Node_FG",
                                                      col_names_posit_i)]
    
    freq_motifs_i_null <- freq_motifs_i_null[,c("Network_id","Type_network",
                                                col_names_freq_i)]
    
    node_positions_i <- bind_rows(node_positions_i,node_positions_i_null)
    freq_motifs_i <- bind_rows(freq_motifs_i,freq_motifs_i_null)
    
  }
  
  # Save results in folder
  folder_motifs_pos <- "Data/Csv/Motifs positions and null models"
  file_motifs_pos <- paste0(folder_motifs_pos,"/Binary_Mat_Sp_Motifs_positions_null_Conn_FiX_",list_Network_id[i.network],".csv") 
  write_csv(node_positions_i,file_motifs_pos)
  
  folder_motifs_fre <- "Data/Csv/Motifs frequencies and null models"
  file_motifs_fre <- paste0(folder_motifs_fre,"/Binary_Mat_Sp_Motifs_frequencies_null_Conn_FiX_",list_Network_id[i.network],".csv") 
  write_csv(freq_motifs_i,file_motifs_fre)
  
}
