
get_motif_probabilities <- function(motif_id_i_probability){
  
  motif_id_i_probability$motif_functional_ID <- NA
  
  for(i in 1:nrow(motif_id_i_probability)){
    
    motif_ID_aux <- motif_id_i_probability[i,c(1:(length(motif_id_i_probability)-2))]
    motif_FG_ID_aux <- motif_ID_aux %>% gather("position","Node_FG") %>% 
      mutate(position = gsub("[^0-9.-]", "", position)) %>% 
      arrange(position,Node_FG)
    
    motif_fucntional_ID <- paste(motif_id_i,
                                 paste(motif_FG_ID_aux$Node_FG,collapse = "_"),
                                 sep = "_")
    motif_id_i_probability$motif_functional_ID[i] <- motif_fucntional_ID
    
  }
  
  # Filter repeated codes
  motif_i_probability <- motif_id_i_probability[1,]
  
  for(i in 2:nrow(motif_id_i_probability)){
    
    if(!motif_id_i_probability$motif_functional_ID[i] %in% 
       motif_i_probability$motif_functional_ID){
      motif_i_probability <- bind_rows(motif_i_probability,
                                       motif_id_i_probability[i,])
    }
    
  }
  
  # Add motif combination number for motif of type i
  motif_i_probability$motif_combination <- as.numeric(rownames(motif_i_probability))
  
  # Add node probability
  motif_i_probability_final <- motif_i_probability %>% gather(position,Node_FG,-motif,
                                                              -motif_functional_ID,
                                                              -motif_combination) %>%
    mutate(position = gsub("[^0-9.-]", "", position),
           position = as.numeric(position)) %>%
    left_join(node_prob, by = c("position", "Node_FG")) %>% 
    rename(node_probability = probability)
  
  # Add motif probability
  motif_i_probability_final$motif_probability <- NA
  
  for(i in 1:max(motif_i_probability_final$motif_combination)){
    motif_combo_positions <- which(motif_i_probability_final$motif_combination == i)
    motif_combi_aux <- motif_i_probability_final %>% filter(motif_combination == i)
    motif_i_probability_final$motif_probability[motif_combo_positions] <- 
      prod(motif_combi_aux$node_probability)
  }
  
  return(motif_i_probability_final)
  
}
