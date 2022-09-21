
library(igraph)
library(tidyverse)
library(DescTools)


# Load data networks raw data and select those interactions whose Interaction frequency is 
# greater than a given threshold

int.threshold <- 1

networks <- read_csv("Data/Csv/data_for_motifs_analysis_1.csv") %>% select(-X1) %>%
  filter(Interaction >= int.threshold)

# MOTIF MINING --------

# Get the list of networks and motif patterns that are available to iterate
list_Network_id <- networks$Network_id %>% unique()

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
  
  
  degree_sequence <- igraph::degree(g_i) %>% as.numeric()
  observed_gini <- DescTools::Gini(degree_sequence, unbiased=FALSE)
  
  # Create incidence matrix
  incidence_mat_i <- as_incidence_matrix(g_i)
  
  # Generate null models for a given incidence matrix------
  
  null_mod <- bipartite::nullmodel(incidence_mat_i, N = 1000, method = "vaznull")
  
  null_ginis_network_i <- NULL
  
  for(i.null in 1:length(null_mod)){
    
    if (sum(null_mod[[i.null]] > 1)){print("ALERT!! NON-BINARY MATRIX")}
    
    null_g_i <- igraph::graph_from_incidence_matrix(null_mod[[i.null]])
    null_degree_sequence <- igraph::degree(null_g_i) %>% as.numeric()
    null_gini <- DescTools::Gini(null_degree_sequence, unbiased=FALSE)
    
    null_ginis_network_i <- c(null_ginis_network_i, null_gini)
    
  }
  
  percentil_gini <- sum(null_ginis_network_i<observed_gini)/length(null_ginis_network_i)
  
  network_i_data <- tibble(Network_id=list_Network_id[i.network],
                           type=c("observed",rep("null",1000)),
                           gini=c(observed_gini,null_ginis_network_i),
                           percentile_observed_gini=percentil_gini)
  
  # Save results in folder
  folder_gini_pos <- "Data/Csv/gini_coefficients"
  file_gini <- paste0(folder_gini_pos,"/gini_null_Conn_FiX_",list_Network_id[i.network],".csv") 
  write_csv(network_i_data,file_gini)
  
}
