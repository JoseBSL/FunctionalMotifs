
# This script extracts the interactions of (up to 5 nodes) motifs and save them in the
# folder "Data/Csv/Motifs links"

# Load libraries
library(tidyverse)
library(igraph)

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

min_num_motifs <- 0
max_num_motifs <- 10000

small_network_motifs <- read_csv("Data/Csv/network_frequency_motifs.csv") %>% filter(nodes<=5) %>%
  group_by(Network_id) %>% count(wt = frequency) %>% filter(n >= min_num_motifs,
                                                            n <= max_num_motifs) %>%
  dplyr::select(Network_id) %>% pull()

small_network_index <- which(list_Network_id %in% small_network_motifs)

# Extract motifs' links for selected networks

for (i.network in small_network_index[5:length(small_network_index)]){#1:length(list_Network_id)){
  
  start_time <- Sys.time()
  print(list_Network_id[i.network])
  print(i.network)
  
  # Initialize variables for storing purposes
  motifs_connections <- NULL #Store the connections of each motif
  
  # Create a graph for the i-th network
  networks_i <- networks %>% filter(Network_id == list_Network_id[i.network])
  edge_list_i <- networks_i[,c("Pollinator_species","Plant_species")]
  
  g_i <- igraph::graph_from_edgelist(as.matrix(edge_list_i), directed = F)
  
  # To create a two mode network we define the types
  V(g_i)$type <- bipartite_mapping(g_i)$type  ## Add the "type" attribute to the network.
  
  for(i.pattern_index in 1:length(list_parterns)){
    
    i.pattern <- list_parterns[i.pattern_index]
    
    # Create a graph for the i-th pattern
    
    pattern_i <- patterns %>% filter(motif_id == list_parterns[i.pattern])
    edge_pattern_i  <- pattern_i [,c("pollinator","plant")]
    p_i <- igraph::graph_from_edgelist(as.matrix(edge_pattern_i), directed = FALSE)
    
    # To create a two mode network we define the types
    V(p_i)$type <- bipartite_mapping(p_i)$type
    #plot.igraph(p_i)
    #degree_distribution(p_i)
    
    iso <- subgraph_isomorphisms(p_i, g_i)
    
    # WARNING subgraph_isomorphisms works with directed graphs and
    # at the end of the day it duplicates the motifs in our list
    # check manual: https://igraph.org/r/doc/subgraph_isomorphisms.html
    
    motifs <- lapply(iso, function (x) { induced_subgraph(g_i, x) })
    
    # Filter non-correct motifs-----
    
    # Remove elements with wrong number of pollinators in the pattern
    # or incorrect edge distribution
    
    pollinators_pattern_i <- pattern_i$pollinator %>% unique() %>% length()
    
    if(length(motifs)>0){
     
      # Create a df to store node information
      
      df_motif_i <- data.frame(matrix(ncol = 2, nrow = length(motifs)))
      colnames(df_motif_i) <- c("nodes_list","use")
      df_motif_i$use <- T
      df_motif_i$motif_i_ID <- rownames(df_motif_i)
      
      # fil_init_time <- Sys.time()
      
      for (i.motifs in 1:length(motifs)){
        
        nodes_i <- V(motifs[[i.motifs]])$name %>% sort()
        
        df_motif_i$nodes_list[i.motifs] <- I(list(nodes_i))
        
        pollinators_motif_i <- sum(nodes_i %in% networks_i$Pollinator_species)
        
        if(pollinators_motif_i!=pollinators_pattern_i){
          
          df_motif_i$use[i.motifs] <- F
          
        }else if(all(degree_distribution(p_i) != degree_distribution(motifs[[i.motifs]]))){
          
          df_motif_i$use[i.motifs] <- F
          
        }
        
      }
      
      # fil_end_time <- Sys.time()
      # fil_end_time - fil_init_time # Time difference of 10.71666 mins for 94008 motifs
      
      # Remove wrong motifs 
      df_motif_i <- df_motif_i %>% filter(use==T)
      
      if(nrow(df_motif_i)>0){ # If there are
        
        #repetitions <- df_motif_i %>% group_by(nodes_list) %>% count()
        
        # Get the index of each unique set of nodes
        species_index <- df_motif_i %>% group_by(nodes_list) %>%
          slice(which.min(motif_i_ID)) %>% data.frame() %>% 
          select(nodes_list,motif_i_ID)
        
        
        for(i.motifs.index in 1:nrow(species_index)){
          
          i.motifs <- as.numeric(species_index$motif_i_ID[i.motifs.index])
          
          # Get edge list for each motif
          
          nodes_i <- V(motifs[[i.motifs]])$name
          
          edge_list_motif_i <- get.edges(motifs[[i.motifs]],es = E(motifs[[i.motifs]])) %>%
            as.tibble() %>%
            mutate(Plant_species = nodes_i[V1],Pollinator_species = nodes_i[V2]) %>% select(-V1,-V2)
          
          # Check and correct node roles: plant/pollinator
          
          for(i.edge in 1:nrow(edge_list_motif_i)){
            
            if(!edge_list_motif_i$Pollinator_species[i.edge] %in% edge_list_i$Pollinator_species){
              poll_aux <- edge_list_motif_i$Plant_species[i.edge]
              edge_list_motif_i$Plant_species[i.edge] <- edge_list_motif_i$Pollinator_species[i.edge]
              edge_list_motif_i$Pollinator_species[i.edge] <- poll_aux
            }
          }
          
          # Add network id, motif pattern and motif number to our edge list
          
          col_names_edge_list_motif_i <- colnames(edge_list_motif_i)
          
          edge_list_motif_i$Network_id  <- list_Network_id[i.network]
          edge_list_motif_i$Motif_pattern_id <- i.pattern
          edge_list_motif_i$Nodes <- length(nodes_i)
          edge_list_motif_i$Motif_number <- i.motifs.index
          
          # Rearrange columns
          
          edge_list_motif_i <- edge_list_motif_i[,c("Network_id","Motif_pattern_id","Nodes",
                                                    "Motif_number",col_names_edge_list_motif_i)]
          
          # Store edge_list
          
          motifs_connections <- bind_rows(motifs_connections, edge_list_motif_i)
        }
        
      }
       
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
