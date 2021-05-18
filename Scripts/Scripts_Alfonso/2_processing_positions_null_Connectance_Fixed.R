# We wnat to know:
#   1) Is the network structure significant?
#   2) Are FG linked to positions?

# NOTES
# 
# para los cálculos de cuantiles he decidido asignarle un NA a las "posiciones de polinizadores"
# cuando analizo GF de plantas, y vice versa (a las posiciones de plantas, cuando analizo
#                                             polinizadores).
# El motivo es que la frecuencia con la que una posición de polinizadores aparece en un nodo planta es cero tanto en las redes observadas como en el modelo nulo, y ello genera este tipo de resultados:
#   x <- rep(0,7)
# quantile(x, type=1)
# 0%  25%  50%  75% 100%
# 0     0   0     0   0
# Es decir, el cero observado se corresponde con todos los percentiles de la distribución random.
# 
# Para poder reportar un único valor de percentil y distinguir los 0% y 100% que realmente
# son significativos de los que aparecen como un artefacto de una distribución uniforme,
# la asignación de los NAs es la única manera que se me ha ocurrido.
# 
# Si no me equivoco, una posible consecuencia de la introducción de los NAs es,
# por ejemplo, la necesidad de ajustar 2 GLMMs para las posiciones: uno de polinizadores
# y otro para las plantas.


# Load libraries
library(tidyverse)
source("Scripts/Scripts_Alfonso/sizeclass_normalization.R") # To normalize frequencies of FGs

# List the files with information about the positions
folder_motif_data <- "Data/Csv/Motifs_positions_and_null_models"

motif_files <- list.files(folder_motif_data) 
motif_link_files <-  motif_files[grep("Binary_Mat_Sp_Motifs_positions_null_Conn_FiX_",motif_files)]
pollinator_positions_codes <- c("np2","np4","np6","np8","np11","np12","np14","np16",
                                "np18","np21","np22","np24","np25","np28","np29","np31",
                                "np34","np35","np38","np41","np42","np44","np46")
plant_positions_codes <- c("np1","np3","np5","np7","np9","np10","np13","np15",
                                "np17","np19","np20","np23","np26","np27","np30","np32",
                                "np33","np36","np37","np39","np40","np43","np45")


for(i in 1:length(motif_link_files)){
  
  start_time <- Sys.time()
  print(motif_link_files[i])
  
  # Open link file and rearrange it
  file_i <- paste0(folder_motif_data, "/", motif_link_files[i])
  positions_i <- read_csv(file_i)
  positions_i_observed <- positions_i %>% filter(Type_network=="observed")
  # 
  # # Plot matrix
  # library('plot.matrix')
  # positions_i_observed_sort <- positions_i_observed %>% arrange(Node_FG)
  # matrix_obs_i <- as.matrix(log10(positions_i_observed_sort[,5:ncol(positions_i_observed_sort)]))
  # row.names(matrix_obs_i) <- paste0(positions_i_observed_sort$Node_id," (",
  #                                   positions_i_observed_sort$Node_FG,")")
  # par(mar=c(4, 22, 2, 3.5)) # adapt margins
  # plot(matrix_obs_i ,axis.row =list(side=2, las=2),
  #      main=unique(positions_i_observed_sort$Network_id),xlab='',ylab='',
  #      col = brewer.pal(9, "RdYlBu"))#cm.colors)
  
  ##############
  # Aggregating nodes and 
  
  positions_i_agg_raw <- positions_i %>% select(-Node_id) %>%
    group_by(Network_id,Type_network,Node_FG) %>% summarise_all(sum)
  
  ##############
  # Normalize frequencies by sizeclass:
  positions_i_agg <- sizeclass_normalization(positions_i_agg_raw)
  
  ###############
  # Initialize the variables to store results
  positions_i_observed_agg <- positions_i_agg %>% filter(Type_network=="observed")
  positions_i_signif_agg_singleFG <- positions_i_observed_agg # percentil single FG
  positions_i_signif_agg <- positions_i_observed_agg # percentil all FG
  positions_i_average_freq_agg_singleFG <- positions_i_observed_agg # percentil single FG
  positions_i_average_freq_agg <- positions_i_observed_agg # percentil all FG
  
  list_positions_agg <- colnames(positions_i_observed_agg)[4:ncol(positions_i_observed_agg)]
  list_nodes_agg <- positions_i_observed_agg$Node_FG
  
  list_plants_FG <- as.character(1:10)
  
  for(i.FG in 1:length(list_nodes_agg)){
    
    # For comparison with this very FG
    random_pos_FG_i_agg_singleFG <- positions_i_agg %>% filter(Type_network!="observed",
                                                               Node_FG == list_nodes_agg[i.FG])
    
    row_index_agg_singleFG <- which(positions_i_observed_agg$Node_FG == list_nodes_agg[i.FG])
    
    # For comparison with all the FG
    if (list_nodes_agg[i.FG] %in% list_plants_FG){ # For comparison with all the FG: Plants
      
      # For comparison with all FG
      random_pos_FG_i_agg <- positions_i_agg %>% filter(Type_network!="observed",
                                                        Node_FG %in% list_plants_FG)
      
      row_index_agg <- which(positions_i_observed_agg$Node_FG %in% list_plants_FG)
      
    }else{ # For comparison with all FG: Animals
      
      random_pos_FG_i_agg <- positions_i_agg %>% filter(Type_network!="observed",
                                                        !Node_FG %in% list_plants_FG)
      
      row_index_agg <- which(!positions_i_observed_agg$Node_FG %in% list_plants_FG)
    }
    
    
    for(i.pos in 1:length(list_positions_agg)){
      
      col_index_agg <- which(colnames(positions_i_observed_agg)==list_positions_agg[i.pos])
      
      values_pos_FG_i_agg_singleFG <- random_pos_FG_i_agg_singleFG[,col_index_agg] %>% pull()
      values_pos_FG_i_agg <- random_pos_FG_i_agg[,col_index_agg] %>% pull()
      
      # For comparison with this very FG
      threshold_singleFG <- positions_i_observed_agg[row_index_agg_singleFG,col_index_agg] %>%
        pull()
      
      lower_less_singleFG <- sum(values_pos_FG_i_agg_singleFG < threshold_singleFG) /
        length(values_pos_FG_i_agg_singleFG)
      positions_i_signif_agg_singleFG[row_index_agg_singleFG,col_index_agg] <- lower_less_singleFG
      
      # Add av. freq for this very FG
      
      positions_i_average_freq_agg_singleFG[row_index_agg_singleFG,col_index_agg] <- 
        mean(values_pos_FG_i_agg_singleFG)
      
      
      for (row_index_agg_type in row_index_agg){
        
        # For comparison with all FG
        threshold <- positions_i_observed_agg[row_index_agg_type,col_index_agg] %>% pull()
        
        lower_less <- sum(values_pos_FG_i_agg < threshold) / length(values_pos_FG_i_agg)
        positions_i_signif_agg[row_index_agg_type,col_index_agg] <- lower_less
        
        # Add av. freq for all FG
        positions_i_average_freq_agg[row_index_agg_singleFG,col_index_agg] <- 
          mean(values_pos_FG_i_agg)
        
        # # Estimate z-score
        # 
        # threshold <- positions_i_observed_agg[row_index_agg_type,col_index_agg] %>% pull()
        # positions_i_signif_agg[row_index_agg_type,col_index_agg] <- 
        #   (threshold - mean(values_pos_FG_i_agg)) / sd(values_pos_FG_i_agg)
        
      }
      
    }
  }
  
  # Set to NA percentiles of pollinators in plant positions
  erase.row <- which(positions_i_signif_agg$Node_FG %in% list_plants_FG)
  erase.col <- which(colnames(positions_i_signif_agg) %in% pollinator_positions_codes)
  
  positions_i_signif_agg_singleFG[erase.row,erase.col] <- NA
  positions_i_signif_agg[erase.row,erase.col] <- NA
  positions_i_average_freq_agg_singleFG[erase.row,erase.col] <- NA
  positions_i_average_freq_agg[erase.row,erase.col] <- NA
  
  # Set to NA percentiles of plants in pollinator positions
  erase.row <- which(!positions_i_signif_agg$Node_FG %in% list_plants_FG)
  erase.col <- which(colnames(positions_i_signif_agg) %in% plant_positions_codes)
  
  positions_i_signif_agg_singleFG[erase.row,erase.col] <- NA
  positions_i_signif_agg[erase.row,erase.col] <- NA
  positions_i_average_freq_agg_singleFG[erase.row,erase.col] <- NA
  positions_i_average_freq_agg[erase.row,erase.col] <- NA
  
  
  #######################
  # Joint all tables
  
  tibl_obs_freq <- positions_i_observed_agg %>% select(-Type_network) %>%
    gather("position", "observed_freq",-Network_id,-Node_FG)
  tibl_exp_freq <- positions_i_average_freq_agg_singleFG %>% select(-Type_network) %>%
    gather("position", "expected_freq_its_GF",-Network_id,-Node_FG)
  tibl_exp_freq_allGF <- positions_i_average_freq_agg %>% select(-Type_network) %>%
    gather("position", "expected_freq_all_GF",-Network_id,-Node_FG)
  tibl_exp_perc <- positions_i_signif_agg_singleFG %>% select(-Type_network) %>%
    gather("position", "percentil_its_GF",-Network_id,-Node_FG)
  tibl_exp_perc_allGF <- positions_i_signif_agg %>% select(-Type_network) %>%
    gather("position", "percentil_all_GF",-Network_id,-Node_FG)
  
  
  results <- tibl_obs_freq %>% left_join(tibl_exp_freq, by = c("Network_id",
                                                               "Node_FG","position")) %>%
    left_join(tibl_exp_freq_allGF, by = c("Network_id","Node_FG","position")) %>%
    left_join(tibl_exp_perc, by = c("Network_id","Node_FG","position")) %>%
    left_join(tibl_exp_perc_allGF, by = c("Network_id","Node_FG","position"))
  
  # Save results
  
  # save file with the full info of motifs positions and the summary of functional motifs
  network_i <- positions_i_signif_agg$Network_id %>% unique()
  new_file_i <- paste0(folder_motif_data, "/Position_quantiles_Conn_Fixed_", network_i,".csv")
  write_csv(results,new_file_i)

  end_time <- Sys.time()
  print(end_time-start_time)
  
}


####################
# Save every Position_quantiles_Conn_Fixed in a single file

# List the files with information about the positions
folder_motif_data <- "Data/Csv/Motifs_positions_and null_models"

motif_files <- list.files(folder_motif_data) 
quant_link_files <-  motif_files[grep("Position_quantiles_Conn_Fixed_",motif_files)]

GF_positions_frequency_percentile <- NULL

for(i in 1:length(quant_link_files)){
  
  file_i <- paste0(folder_motif_data, "/", quant_link_files[i])
  quantiles_i <- read_csv(file_i)
  
  GF_positions_frequency_percentile <- bind_rows(GF_positions_frequency_percentile,quantiles_i)
  
}

final_file_i <- paste0(folder_motif_data, "/GF_positions_frequency_percentile.csv")
write_csv(GF_positions_frequency_percentile,final_file_i)
