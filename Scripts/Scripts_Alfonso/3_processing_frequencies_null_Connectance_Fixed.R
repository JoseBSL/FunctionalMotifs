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

# List the files with information about the frequencies of motifs
folder_motif_data <- "Data/Csv/Motifs frequencies and null models"

motif_files <- list.files(folder_motif_data) 
motif_freq_files <-  motif_files[grep("Binary_Mat_Sp_Motifs_frequencies_null_Conn_FiX_",
                                      motif_files)]

for(i in 1:length(motif_freq_files)){
  
  start_time <- Sys.time()
  print(motif_freq_files[i])
  
  # Open link file and rearrange it
  file_i <- paste0(folder_motif_data, "/", motif_freq_files[i])
  freq_i <- read_csv(file_i)
  freq_i_observed <- freq_i %>% filter(Type_network=="observed") %>% 
    select(Network_id,motif,nodes,normalise_sizeclass) %>% 
    rename(observed_freq_sizeclass = normalise_sizeclass)
  
  
  exp_freq_i <- freq_i %>% filter(Type_network!="observed") %>% 
    select(Network_id,motif,nodes,normalise_sizeclass) %>% 
    rename(expected_freq_sizeclass = normalise_sizeclass) %>%
    group_by(Network_id,motif,nodes) %>% summarize(m = mean(expected_freq_sizeclass)) %>%
    rename(expected_freq_sizeclass = m) %>% ungroup()

  perc_freq_i_observed <- freq_i_observed %>% rename(percentil_sizeclass = observed_freq_sizeclass)
  
  for(j in 1:nrow(perc_freq_i_observed)){
    
    motif_j <- perc_freq_i_observed$motif[j]
    threshold_j <- perc_freq_i_observed$percentil_sizeclass[j]
    random_freq_i <- freq_i %>% filter(Type_network!="observed",motif == motif_j) %>%
      select(normalise_sizeclass) %>% pull()
    lower_less <- sum(random_freq_i < threshold_j) / length(random_freq_i)
    perc_freq_i_observed$percentil_sizeclass[j] <- lower_less
    
  }
   
  # Create resulting table
  
  result <- freq_i_observed %>% left_join(exp_freq_i, by = c("Network_id","motif","nodes")) %>%
    left_join(perc_freq_i_observed, by = c("Network_id","motif","nodes"))
  
  # Save results
  
  # save file with the full info of motifs positions and the summary of functional motifs
  network_i <- freq_i$Network_id %>% unique()
  new_file_i <- paste0(folder_motif_data, "/Motif_quantiles_Conn_Fixed_", network_i,".csv")
  write_csv(result,new_file_i)

  end_time <- Sys.time()
  print(end_time-start_time)
  
}

#############################

# Save every Motif_quantiles_Conn_Fixed_ in a single file

# List the files with information about the positions
folder_motif_data <- "Data/Csv/Motifs frequencies and null models"

motif_files <- list.files(folder_motif_data) 
quant_link_files <-  motif_files[grep("Motif_quantiles_Conn_Fixed_",motif_files)]

Motifs_frequency_percentile <- NULL

for(i in 1:length(quant_link_files)){
  
  file_i <- paste0(folder_motif_data, "/", quant_link_files[i])
  quantiles_i <- read_csv(file_i)
  
  Motifs_frequency_percentile <- bind_rows(Motifs_frequency_percentile,quantiles_i)
  
}

final_file_i <- paste0(folder_motif_data, "/Motifs_frequency_percentile.csv")
write_csv(Motifs_frequency_percentile,final_file_i)
