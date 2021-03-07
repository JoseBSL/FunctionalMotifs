

connections_pattern_i <- function(g_i,motifs_connections,patterns,
                                  list_parterns,i.pattern,df_motif_8_17,df_motif_9_13,
                                  df_motif_10_14,df_motif_11_15,df_motif_12_16){
  
  if (!i.pattern %in% c(8,17,9,13,10,14,11,15,12,16)){
    
    # Create a graph for the i-th pattern
    
    pattern_i <- patterns %>% filter(motif_id == list_parterns[i.pattern])
    edge_pattern_i  <- pattern_i[,c("pollinator","plant")]
    p_i <- igraph::graph_from_edgelist(as.matrix(edge_pattern_i), directed = FALSE)
    
    # To create a two mode network we define the types
    V(p_i)$type <- bipartite_mapping(p_i)$type
    
    # Calculate subgraph isomorphisms by using the patern of a given motif
    
    # WARNING subgraph_isomorphisms works with directed graphs and
    # at the end of the day it duplicates the motifs in our list
    # check manual: https://igraph.org/r/doc/subgraph_isomorphisms.html
    
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
        
        pollinators_motif_i <- sum(nodes_i %in% networks_i$Pollinator_functional_group)
        
        if(pollinators_motif_i!=pollinators_pattern_i){
          
          df_motif_i$use[i.motifs] <- F
          
        }else if(any(degree_distribution(p_i) != degree_distribution(motifs[[i.motifs]]))){
          
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
            mutate(Plant_functional_groups = nodes_i[V1],
                   Pollinator_functional_group = nodes_i[V2]) %>% select(-V1,-V2)
          
          # Check and correct node roles: plant/pollinator
          
          for(i.edge in 1:nrow(edge_list_motif_i)){
            
            if(!edge_list_motif_i$Pollinator_functional_group[i.edge] %in% edge_list_i$Pollinator_functional_group){
              poll_aux <- edge_list_motif_i$Plant_functional_groups[i.edge]
              edge_list_motif_i$Plant_functional_groups[i.edge] <- edge_list_motif_i$Pollinator_functional_group[i.edge]
              edge_list_motif_i$Pollinator_functional_group[i.edge] <- poll_aux
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
    
    return(motifs_connections)
    
  }else if(i.pattern == 8){
    
    # Create a graph for the i-th pattern
    
    pattern_i <- patterns %>% filter(motif_id == list_parterns[i.pattern])
    edge_pattern_i  <- pattern_i[,c("pollinator","plant")]
    p_i <- igraph::graph_from_edgelist(as.matrix(edge_pattern_i), directed = FALSE)
    
    # To create a two mode network we define the types
    V(p_i)$type <- bipartite_mapping(p_i)$type
    
    # Calculate subgraph isomorphisms by using the patern of a given motif
    
    # WARNING subgraph_isomorphisms works with directed graphs and
    # at the end of the day it duplicates the motifs in our list
    # check manual: https://igraph.org/r/doc/subgraph_isomorphisms.html
    
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
      df_motif_i$use <- 8
      df_motif_i$motif_i_ID <- rownames(df_motif_i)
      
      # fil_init_time <- Sys.time()
      
      for (i.motifs in 1:length(motifs)){
        
        nodes_i <- V(motifs[[i.motifs]])$name %>% sort()
        
        df_motif_i$nodes_list[i.motifs] <- I(list(nodes_i))
        
        pollinators_motif_i <- sum(nodes_i %in% networks_i$Pollinator_functional_group)
        
        if(any(degree_distribution(p_i) != degree_distribution(motifs[[i.motifs]]))){
          
          df_motif_i$use[i.motifs] <- 0
          
        }else if(pollinators_motif_i!=pollinators_pattern_i){
          
          df_motif_i$use[i.motifs] <- 17
          
        }
        
      }
      
      # Store the resulting df
      
      df_motif_8_17 <- df_motif_i
      
      # Remove wrong motifs 
      df_motif_i <- df_motif_i %>% filter(use==8)
      
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
            mutate(Plant_functional_groups = nodes_i[V1],Pollinator_functional_group = nodes_i[V2]) %>% select(-V1,-V2)
          
          # Check and correct node roles: plant/pollinator
          
          for(i.edge in 1:nrow(edge_list_motif_i)){
            
            if(!edge_list_motif_i$Pollinator_functional_group[i.edge] %in% edge_list_i$Pollinator_functional_group){
              poll_aux <- edge_list_motif_i$Plant_functional_groups[i.edge]
              edge_list_motif_i$Plant_functional_groups[i.edge] <- edge_list_motif_i$Pollinator_functional_group[i.edge]
              edge_list_motif_i$Pollinator_functional_group[i.edge] <- poll_aux
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
    
    return(list(motifs_connections,df_motif_8_17))
    
  }else if(i.pattern == 9){
    
    # Create a graph for the i-th pattern
    
    pattern_i <- patterns %>% filter(motif_id == list_parterns[i.pattern])
    edge_pattern_i  <- pattern_i[,c("pollinator","plant")]
    p_i <- igraph::graph_from_edgelist(as.matrix(edge_pattern_i), directed = FALSE)
    
    # To create a two mode network we define the types
    V(p_i)$type <- bipartite_mapping(p_i)$type
    
    # Calculate subgraph isomorphisms by using the patern of a given motif
    
    # WARNING subgraph_isomorphisms works with directed graphs and
    # at the end of the day it duplicates the motifs in our list
    # check manual: https://igraph.org/r/doc/subgraph_isomorphisms.html
    
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
      df_motif_i$use <- 9
      df_motif_i$motif_i_ID <- rownames(df_motif_i)
      
      # fil_init_time <- Sys.time()
      
      for (i.motifs in 1:length(motifs)){
        
        nodes_i <- V(motifs[[i.motifs]])$name %>% sort()
        
        df_motif_i$nodes_list[i.motifs] <- I(list(nodes_i))
        
        pollinators_motif_i <- sum(nodes_i %in% networks_i$Pollinator_functional_group)
        
        if(any(degree_distribution(p_i) != degree_distribution(motifs[[i.motifs]]))){
          
          df_motif_i$use[i.motifs] <- 0
          
        }else if(pollinators_motif_i!=pollinators_pattern_i){
          
          df_motif_i$use[i.motifs] <- 13
          
        }
        
      }
      
      # Store the resulting df
      
      df_motif_9_13 <- df_motif_i
      
      # Remove wrong motifs 
      df_motif_i <- df_motif_i %>% filter(use==9)
      
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
            mutate(Plant_functional_groups = nodes_i[V1],Pollinator_functional_group = nodes_i[V2]) %>% select(-V1,-V2)
          
          # Check and correct node roles: plant/pollinator
          
          for(i.edge in 1:nrow(edge_list_motif_i)){
            
            if(!edge_list_motif_i$Pollinator_functional_group[i.edge] %in% edge_list_i$Pollinator_functional_group){
              poll_aux <- edge_list_motif_i$Plant_functional_groups[i.edge]
              edge_list_motif_i$Plant_functional_groups[i.edge] <- edge_list_motif_i$Pollinator_functional_group[i.edge]
              edge_list_motif_i$Pollinator_functional_group[i.edge] <- poll_aux
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
    
    return(list(motifs_connections,df_motif_9_13))
    
  }else if(i.pattern == 10){
    
    # Create a graph for the i-th pattern
    
    pattern_i <- patterns %>% filter(motif_id == list_parterns[i.pattern])
    edge_pattern_i  <- pattern_i[,c("pollinator","plant")]
    p_i <- igraph::graph_from_edgelist(as.matrix(edge_pattern_i), directed = FALSE)
    
    # To create a two mode network we define the types
    V(p_i)$type <- bipartite_mapping(p_i)$type
    
    # Calculate subgraph isomorphisms by using the patern of a given motif
    
    # WARNING subgraph_isomorphisms works with directed graphs and
    # at the end of the day it duplicates the motifs in our list
    # check manual: https://igraph.org/r/doc/subgraph_isomorphisms.html
    
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
      df_motif_i$use <- 10
      df_motif_i$motif_i_ID <- rownames(df_motif_i)
      
      # fil_init_time <- Sys.time()
      
      for (i.motifs in 1:length(motifs)){
        
        nodes_i <- V(motifs[[i.motifs]])$name %>% sort()
        
        df_motif_i$nodes_list[i.motifs] <- I(list(nodes_i))
        
        pollinators_motif_i <- sum(nodes_i %in% networks_i$Pollinator_functional_group)
        
        if(any(degree_distribution(p_i) != degree_distribution(motifs[[i.motifs]]))){
          
          df_motif_i$use[i.motifs] <- 0
          
        }else if(pollinators_motif_i!=pollinators_pattern_i){
          
          df_motif_i$use[i.motifs] <- 14
          
        }
        
      }
      
      # Store the resulting df
      
      df_motif_10_14 <- df_motif_i
      
      # Remove wrong motifs 
      df_motif_i <- df_motif_i %>% filter(use==10)
      
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
            mutate(Plant_functional_groups = nodes_i[V1],Pollinator_functional_group = nodes_i[V2]) %>% select(-V1,-V2)
          
          # Check and correct node roles: plant/pollinator
          
          for(i.edge in 1:nrow(edge_list_motif_i)){
            
            if(!edge_list_motif_i$Pollinator_functional_group[i.edge] %in% edge_list_i$Pollinator_functional_group){
              poll_aux <- edge_list_motif_i$Plant_functional_groups[i.edge]
              edge_list_motif_i$Plant_functional_groups[i.edge] <- edge_list_motif_i$Pollinator_functional_group[i.edge]
              edge_list_motif_i$Pollinator_functional_group[i.edge] <- poll_aux
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
    
    return(list(motifs_connections,df_motif_10_14))
    
  }else if(i.pattern == 11){
    
    # Create a graph for the i-th pattern
    
    pattern_i <- patterns %>% filter(motif_id == list_parterns[i.pattern])
    edge_pattern_i  <- pattern_i[,c("pollinator","plant")]
    p_i <- igraph::graph_from_edgelist(as.matrix(edge_pattern_i), directed = FALSE)
    
    # To create a two mode network we define the types
    V(p_i)$type <- bipartite_mapping(p_i)$type
    
    # Calculate subgraph isomorphisms by using the patern of a given motif
    
    # WARNING subgraph_isomorphisms works with directed graphs and
    # at the end of the day it duplicates the motifs in our list
    # check manual: https://igraph.org/r/doc/subgraph_isomorphisms.html
    
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
      df_motif_i$use <- 11
      df_motif_i$motif_i_ID <- rownames(df_motif_i)
      
      # fil_init_time <- Sys.time()
      
      for (i.motifs in 1:length(motifs)){
        
        nodes_i <- V(motifs[[i.motifs]])$name %>% sort()
        
        df_motif_i$nodes_list[i.motifs] <- I(list(nodes_i))
        
        pollinators_motif_i <- sum(nodes_i %in% networks_i$Pollinator_functional_group)
        
        if(any(degree_distribution(p_i) != degree_distribution(motifs[[i.motifs]]))){
          
          df_motif_i$use[i.motifs] <- 0
          
        }else if(pollinators_motif_i!=pollinators_pattern_i){
          
          df_motif_i$use[i.motifs] <- 15
          
        }
        
      }
      
      # Store the resulting df
      
      df_motif_11_15 <- df_motif_i
      
      # Remove wrong motifs 
      df_motif_i <- df_motif_i %>% filter(use==11)
      
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
            mutate(Plant_functional_groups = nodes_i[V1],Pollinator_functional_group = nodes_i[V2]) %>% select(-V1,-V2)
          
          # Check and correct node roles: plant/pollinator
          
          for(i.edge in 1:nrow(edge_list_motif_i)){
            
            if(!edge_list_motif_i$Pollinator_functional_group[i.edge] %in% edge_list_i$Pollinator_functional_group){
              poll_aux <- edge_list_motif_i$Plant_functional_groups[i.edge]
              edge_list_motif_i$Plant_functional_groups[i.edge] <- edge_list_motif_i$Pollinator_functional_group[i.edge]
              edge_list_motif_i$Pollinator_functional_group[i.edge] <- poll_aux
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
    
    return(list(motifs_connections,df_motif_11_15))
    
  }else if(i.pattern == 12){
    
    # Create a graph for the i-th pattern
    
    pattern_i <- patterns %>% filter(motif_id == list_parterns[i.pattern])
    edge_pattern_i  <- pattern_i[,c("pollinator","plant")]
    p_i <- igraph::graph_from_edgelist(as.matrix(edge_pattern_i), directed = FALSE)
    
    # To create a two mode network we define the types
    V(p_i)$type <- bipartite_mapping(p_i)$type
    
    # Calculate subgraph isomorphisms by using the patern of a given motif
    
    # WARNING subgraph_isomorphisms works with directed graphs and
    # at the end of the day it duplicates the motifs in our list
    # check manual: https://igraph.org/r/doc/subgraph_isomorphisms.html
    
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
      df_motif_i$use <- 12
      df_motif_i$motif_i_ID <- rownames(df_motif_i)
      
      # fil_init_time <- Sys.time()
      
      for (i.motifs in 1:length(motifs)){
        
        nodes_i <- V(motifs[[i.motifs]])$name %>% sort()
        
        df_motif_i$nodes_list[i.motifs] <- I(list(nodes_i))
        
        pollinators_motif_i <- sum(nodes_i %in% networks_i$Pollinator_functional_group)
        
        if(any(degree_distribution(p_i) != degree_distribution(motifs[[i.motifs]]))){
          
          df_motif_i$use[i.motifs] <- 0
          
        }else if(pollinators_motif_i!=pollinators_pattern_i){
          
          df_motif_i$use[i.motifs] <- 16
          
        }
        
      }
      
      # Store the resulting df
      
      df_motif_12_16 <- df_motif_i
      
      # Remove wrong motifs 
      df_motif_i <- df_motif_i %>% filter(use==12)
      
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
            mutate(Plant_functional_groups = nodes_i[V1],Pollinator_functional_group = nodes_i[V2]) %>% select(-V1,-V2)
          
          # Check and correct node roles: plant/pollinator
          
          for(i.edge in 1:nrow(edge_list_motif_i)){
            
            if(!edge_list_motif_i$Pollinator_functional_group[i.edge] %in% edge_list_i$Pollinator_functional_group){
              poll_aux <- edge_list_motif_i$Plant_functional_groups[i.edge]
              edge_list_motif_i$Plant_functional_groups[i.edge] <- edge_list_motif_i$Pollinator_functional_group[i.edge]
              edge_list_motif_i$Pollinator_functional_group[i.edge] <- poll_aux
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
    
    return(list(motifs_connections,df_motif_12_16))
    
  }else if(i.pattern == 13){
    
    # Since we will use our previous calculations for the 8th pattern, we calculated
    # the motifs that are based on them.
    
    i.pattern_ini <- 9
    
    # Create a graph for the i-th pattern
    
    pattern_i <- patterns %>% filter(motif_id == list_parterns[i.pattern_ini])
    edge_pattern_i  <- pattern_i[,c("pollinator","plant")]
    p_i <- igraph::graph_from_edgelist(as.matrix(edge_pattern_i), directed = FALSE)
    
    # To create a two mode network we define the types
    V(p_i)$type <- bipartite_mapping(p_i)$type
    
    # Calculate subgraph isomorphisms by using the patern of a given motif
    
    # WARNING subgraph_isomorphisms works with directed graphs and
    # at the end of the day it duplicates the motifs in our list
    # check manual: https://igraph.org/r/doc/subgraph_isomorphisms.html
    
    iso <- subgraph_isomorphisms(p_i, g_i)
    
    # WARNING subgraph_isomorphisms works with directed graphs and
    # at the end of the day it duplicates the motifs in our list
    # check manual: https://igraph.org/r/doc/subgraph_isomorphisms.html
    
    motifs <- lapply(iso, function (x) { induced_subgraph(g_i, x) })
    
    if(!is.null(df_motif_9_13)){
      
      # Remove wrong motifs 
      df_motif_i <- df_motif_9_13 %>% filter(use==13)
      
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
            mutate(Plant_functional_groups = nodes_i[V1],Pollinator_functional_group = nodes_i[V2]) %>% select(-V1,-V2)
          
          # Check and correct node roles: plant/pollinator
          
          for(i.edge in 1:nrow(edge_list_motif_i)){
            
            if(!edge_list_motif_i$Pollinator_functional_group[i.edge] %in% edge_list_i$Pollinator_functional_group){
              poll_aux <- edge_list_motif_i$Plant_functional_groups[i.edge]
              edge_list_motif_i$Plant_functional_groups[i.edge] <- edge_list_motif_i$Pollinator_functional_group[i.edge]
              edge_list_motif_i$Pollinator_functional_group[i.edge] <- poll_aux
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
    
    return(motifs_connections)
    
  }else if(i.pattern == 14){
    
    # Since we will use our previous calculations for the 8th pattern, we calculated
    # the motifs that are based on them.
    
    i.pattern_ini <- 10
    
    # Create a graph for the i-th pattern
    
    pattern_i <- patterns %>% filter(motif_id == list_parterns[i.pattern_ini])
    edge_pattern_i  <- pattern_i[,c("pollinator","plant")]
    p_i <- igraph::graph_from_edgelist(as.matrix(edge_pattern_i), directed = FALSE)
    
    # To create a two mode network we define the types
    V(p_i)$type <- bipartite_mapping(p_i)$type
    
    # Calculate subgraph isomorphisms by using the patern of a given motif
    
    # WARNING subgraph_isomorphisms works with directed graphs and
    # at the end of the day it duplicates the motifs in our list
    # check manual: https://igraph.org/r/doc/subgraph_isomorphisms.html
    
    iso <- subgraph_isomorphisms(p_i, g_i)
    
    # WARNING subgraph_isomorphisms works with directed graphs and
    # at the end of the day it duplicates the motifs in our list
    # check manual: https://igraph.org/r/doc/subgraph_isomorphisms.html
    
    motifs <- lapply(iso, function (x) { induced_subgraph(g_i, x) })
    
    if(!is.null(df_motif_10_14)){
      
      # Remove wrong motifs 
      df_motif_i <- df_motif_10_14 %>% filter(use==14)
      
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
            mutate(Plant_functional_groups = nodes_i[V1],Pollinator_functional_group = nodes_i[V2]) %>% select(-V1,-V2)
          
          # Check and correct node roles: plant/pollinator
          
          for(i.edge in 1:nrow(edge_list_motif_i)){
            
            if(!edge_list_motif_i$Pollinator_functional_group[i.edge] %in% edge_list_i$Pollinator_functional_group){
              poll_aux <- edge_list_motif_i$Plant_functional_groups[i.edge]
              edge_list_motif_i$Plant_functional_groups[i.edge] <- edge_list_motif_i$Pollinator_functional_group[i.edge]
              edge_list_motif_i$Pollinator_functional_group[i.edge] <- poll_aux
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
    
    return(motifs_connections)
    
  }else if(i.pattern == 15){
    
    # Since we will use our previous calculations for the 8th pattern, we calculated
    # the motifs that are based on them.
    
    i.pattern_ini <- 11
    
    # Create a graph for the i-th pattern
    
    pattern_i <- patterns %>% filter(motif_id == list_parterns[i.pattern_ini])
    edge_pattern_i  <- pattern_i[,c("pollinator","plant")]
    p_i <- igraph::graph_from_edgelist(as.matrix(edge_pattern_i), directed = FALSE)
    
    # To create a two mode network we define the types
    V(p_i)$type <- bipartite_mapping(p_i)$type
    
    # Calculate subgraph isomorphisms by using the patern of a given motif
    
    # WARNING subgraph_isomorphisms works with directed graphs and
    # at the end of the day it duplicates the motifs in our list
    # check manual: https://igraph.org/r/doc/subgraph_isomorphisms.html
    
    iso <- subgraph_isomorphisms(p_i, g_i)
    
    # WARNING subgraph_isomorphisms works with directed graphs and
    # at the end of the day it duplicates the motifs in our list
    # check manual: https://igraph.org/r/doc/subgraph_isomorphisms.html
    
    motifs <- lapply(iso, function (x) { induced_subgraph(g_i, x) })
    
    if(!is.null(df_motif_11_15)){
      
      # Remove wrong motifs 
      df_motif_i <- df_motif_11_15 %>% filter(use==15)
      
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
            mutate(Plant_functional_groups = nodes_i[V1],Pollinator_functional_group = nodes_i[V2]) %>% select(-V1,-V2)
          
          # Check and correct node roles: plant/pollinator
          
          for(i.edge in 1:nrow(edge_list_motif_i)){
            
            if(!edge_list_motif_i$Pollinator_functional_group[i.edge] %in% edge_list_i$Pollinator_functional_group){
              poll_aux <- edge_list_motif_i$Plant_functional_groups[i.edge]
              edge_list_motif_i$Plant_functional_groups[i.edge] <- edge_list_motif_i$Pollinator_functional_group[i.edge]
              edge_list_motif_i$Pollinator_functional_group[i.edge] <- poll_aux
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
    
    return(motifs_connections)
    
  }else if(i.pattern == 16){
    
    # Since we will use our previous calculations for the 8th pattern, we calculated
    # the motifs that are based on them.
    
    i.pattern_ini <- 12
    
    # Create a graph for the i-th pattern
    
    pattern_i <- patterns %>% filter(motif_id == list_parterns[i.pattern_ini])
    edge_pattern_i  <- pattern_i[,c("pollinator","plant")]
    p_i <- igraph::graph_from_edgelist(as.matrix(edge_pattern_i), directed = FALSE)
    
    # To create a two mode network we define the types
    V(p_i)$type <- bipartite_mapping(p_i)$type
    
    # Calculate subgraph isomorphisms by using the patern of a given motif
    
    # WARNING subgraph_isomorphisms works with directed graphs and
    # at the end of the day it duplicates the motifs in our list
    # check manual: https://igraph.org/r/doc/subgraph_isomorphisms.html
    
    iso <- subgraph_isomorphisms(p_i, g_i)
    
    # WARNING subgraph_isomorphisms works with directed graphs and
    # at the end of the day it duplicates the motifs in our list
    # check manual: https://igraph.org/r/doc/subgraph_isomorphisms.html
    
    motifs <- lapply(iso, function (x) { induced_subgraph(g_i, x) })
    
    if(!is.null(df_motif_12_16)){
      
      # Remove wrong motifs 
      df_motif_i <- df_motif_12_16 %>% filter(use==16)
      
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
            mutate(Plant_functional_groups = nodes_i[V1],Pollinator_functional_group = nodes_i[V2]) %>% select(-V1,-V2)
          
          # Check and correct node roles: plant/pollinator
          
          for(i.edge in 1:nrow(edge_list_motif_i)){
            
            if(!edge_list_motif_i$Pollinator_functional_group[i.edge] %in% edge_list_i$Pollinator_functional_group){
              poll_aux <- edge_list_motif_i$Plant_functional_groups[i.edge]
              edge_list_motif_i$Plant_functional_groups[i.edge] <- edge_list_motif_i$Pollinator_functional_group[i.edge]
              edge_list_motif_i$Pollinator_functional_group[i.edge] <- poll_aux
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
    
    return(motifs_connections)
    
  }else if(i.pattern == 17){
    
    # Since we will use our previous calculations for the 8th pattern, we calculated
    # the motifs that are based on them.
    
    i.pattern_ini <- 8
    
    # Create a graph for the i-th pattern
    
    pattern_i <- patterns %>% filter(motif_id == list_parterns[i.pattern_ini])
    edge_pattern_i  <- pattern_i[,c("pollinator","plant")]
    p_i <- igraph::graph_from_edgelist(as.matrix(edge_pattern_i), directed = FALSE)
    
    # To create a two mode network we define the types
    V(p_i)$type <- bipartite_mapping(p_i)$type
    
    # Calculate subgraph isomorphisms by using the patern of a given motif
    
    # WARNING subgraph_isomorphisms works with directed graphs and
    # at the end of the day it duplicates the motifs in our list
    # check manual: https://igraph.org/r/doc/subgraph_isomorphisms.html
    
    iso <- subgraph_isomorphisms(p_i, g_i)
    
    # WARNING subgraph_isomorphisms works with directed graphs and
    # at the end of the day it duplicates the motifs in our list
    # check manual: https://igraph.org/r/doc/subgraph_isomorphisms.html
    
    motifs <- lapply(iso, function (x) { induced_subgraph(g_i, x) })

    if(!is.null(df_motif_8_17)){
      
      # Remove wrong motifs 
      df_motif_i <- df_motif_8_17 %>% filter(use==17)
      
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
            mutate(Plant_functional_groups = nodes_i[V1],Pollinator_functional_group = nodes_i[V2]) %>% select(-V1,-V2)
          
          # Check and correct node roles: plant/pollinator
          
          for(i.edge in 1:nrow(edge_list_motif_i)){
            
            if(!edge_list_motif_i$Pollinator_functional_group[i.edge] %in% edge_list_i$Pollinator_functional_group){
              poll_aux <- edge_list_motif_i$Plant_functional_groups[i.edge]
              edge_list_motif_i$Plant_functional_groups[i.edge] <- edge_list_motif_i$Pollinator_functional_group[i.edge]
              edge_list_motif_i$Pollinator_functional_group[i.edge] <- poll_aux
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
    
    return(motifs_connections)
    
  }
  
  
}

  
  