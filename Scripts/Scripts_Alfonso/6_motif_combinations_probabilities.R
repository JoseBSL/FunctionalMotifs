
library(tidyverse)
source("Scripts/Scripts_Alfonso/get_motif_probabilities_function.R")

pollinator_means <- read_csv("Data/Csv/pollinator_abs_freq_means.csv")
plant_means <- read_csv("Data/Csv/plant_abs_freq_means.csv")

pollinator_node_prob_aux <- pollinator_means %>% select(-mean, -SE) %>% 
  spread(Node_FG,mean_natural_units)

pollinator_node_prob_aux$total <- 
  rowSums(pollinator_node_prob_aux[,c(2:ncol(pollinator_node_prob_aux))])

pollinator_node_prob <- pollinator_node_prob_aux %>% 
  gather("Node_FG","frequency",-position,-total) %>%
  mutate(probability = frequency/total) %>% select(-total,-frequency)

plant_node_prob_aux <- plant_means %>% select(-mean, -SE) %>% 
  spread(Node_FG,mean_natural_units)

plant_node_prob_aux$total <- 
  rowSums(plant_node_prob_aux[,c(2:ncol(plant_node_prob_aux))])

plant_node_prob <- plant_node_prob_aux %>% 
  gather("Node_FG","frequency",-position,-total) %>%
  mutate(probability = frequency/total) %>% select(-total,-frequency)

node_prob <- bind_rows(plant_node_prob,pollinator_node_prob) %>%
  mutate(position = as.numeric(position))

########################################
# CREATE OF POSSIBLE MOTIF COMBINATIONS
########################################

list_pollinator_FG <- pollinator_node_prob$Node_FG %>% unique()
list_plant_FG <- plant_node_prob$Node_FG %>% unique()

print_positions <- function(motif_id_i){
  motif_connections <- read_csv("Data/Data_processing/Motifs_connections/motif_pattern_connections.csv")
  motif_i <- motif_connections %>% filter(motif_id == motif_id_i )
  
  plant_i <- motif_i$plant %>% unique()
  cat("plant: ",plant_i,"\n")
  pollinator_i <- motif_i$pollinator %>% unique()
  cat("pollinator: ",pollinator_i)
}

motif_combinations_prob <- NULL

#-----------------
motif_id_i <- 1
print_positions(motif_id_i)
motif_id_i_probability <- crossing(p1 = list_plant_FG,
                                   p2 = list_pollinator_FG) %>% mutate(motif = motif_id_i)

motif_combinations_prob <- bind_rows(motif_combinations_prob,
                                     get_motif_probabilities(motif_id_i_probability))
#-----------------
motif_id_i <- 2
print_positions(motif_id_i)
motif_id_i_probability <- crossing(p3a = list_plant_FG,
                                   p3b = list_plant_FG,
                                   p4 = list_pollinator_FG) %>% mutate(motif = motif_id_i)

motif_combinations_prob <- bind_rows(motif_combinations_prob,
                                     get_motif_probabilities(motif_id_i_probability))
#-----------------
motif_id_i <- 3
print_positions(motif_id_i)
motif_id_i_probability <- crossing(p5 = list_plant_FG,
                                   p6a = list_pollinator_FG,
                                   p6b = list_pollinator_FG) %>% mutate(motif = motif_id_i)

motif_combinations_prob <- bind_rows(motif_combinations_prob,
                                     get_motif_probabilities(motif_id_i_probability))
#-----------------
motif_id_i <- 4
print_positions(motif_id_i)
motif_id_i_probability <- crossing(p7 = list_plant_FG,
                                   p8a = list_pollinator_FG,
                                   p8b = list_pollinator_FG,
                                   p8c = list_pollinator_FG) %>% mutate(motif = motif_id_i)

motif_combinations_prob <- bind_rows(motif_combinations_prob,
                                     get_motif_probabilities(motif_id_i_probability))

#-----------------
motif_id_i <- 5
print_positions(motif_id_i)
motif_id_i_probability <- crossing(p9 = list_plant_FG,
                                   p10 = list_plant_FG,
                                   p12 = list_pollinator_FG,
                                   p11 = list_pollinator_FG) %>% mutate(motif = motif_id_i)

motif_combinations_prob <- bind_rows(motif_combinations_prob,
                                     get_motif_probabilities(motif_id_i_probability))

#-----------------
motif_id_i <- 6
print_positions(motif_id_i)
motif_id_i_probability <- crossing(p13a = list_plant_FG,
                                   p13b = list_plant_FG,
                                   p14a = list_pollinator_FG,
                                   p14b = list_pollinator_FG) %>% mutate(motif = motif_id_i)

motif_combinations_prob <- bind_rows(motif_combinations_prob,
                                     get_motif_probabilities(motif_id_i_probability))

#-----------------
motif_id_i <- 7
print_positions(motif_id_i)
motif_id_i_probability <- crossing(p15a = list_plant_FG,
                                   p15b = list_plant_FG,
                                   p15c = list_plant_FG,
                                   p16 = list_pollinator_FG) %>% mutate(motif = motif_id_i)

motif_combinations_prob <- bind_rows(motif_combinations_prob,
                                     get_motif_probabilities(motif_id_i_probability))

#-----------------
motif_id_i <- 8
print_positions(motif_id_i)
motif_id_i_probability <- crossing(p17 = list_plant_FG,
                                   p18a = list_pollinator_FG,
                                   p18b = list_pollinator_FG,
                                   p18c = list_pollinator_FG,
                                   p18d = list_pollinator_FG) %>% mutate(motif = motif_id_i)

motif_combinations_prob <- bind_rows(motif_combinations_prob,
                                     get_motif_probabilities(motif_id_i_probability))

#-----------------
motif_id_i <- 9
print_positions(motif_id_i)
motif_id_i_probability <- crossing(p20 = list_plant_FG,
                                   p19 = list_plant_FG,
                                   p21a = list_pollinator_FG,
                                   p21b = list_pollinator_FG,
                                   p22 = list_pollinator_FG) %>% mutate(motif = motif_id_i)

motif_combinations_prob <- bind_rows(motif_combinations_prob,
                                     get_motif_probabilities(motif_id_i_probability))

#-----------------
motif_id_i <- 10
print_positions(motif_id_i)
motif_id_i_probability <- crossing(p23a = list_plant_FG,
                                   p23b = list_plant_FG,
                                   p24a = list_pollinator_FG,
                                   p24b = list_pollinator_FG,
                                   p25 = list_pollinator_FG) %>% mutate(motif = motif_id_i)

motif_combinations_prob <- bind_rows(motif_combinations_prob,
                                     get_motif_probabilities(motif_id_i_probability))

#-----------------
motif_id_i <- 11
print_positions(motif_id_i)
motif_id_i_probability <- crossing(p27 = list_plant_FG,
                                   p26 = list_plant_FG,
                                   p28 = list_pollinator_FG,
                                   p29a = list_pollinator_FG,
                                   p29b = list_pollinator_FG) %>% mutate(motif = motif_id_i)

motif_combinations_prob <- bind_rows(motif_combinations_prob,
                                     get_motif_probabilities(motif_id_i_probability))

#-----------------
motif_id_i <- 12
print_positions(motif_id_i)
motif_id_i_probability <- crossing(p30a = list_plant_FG,
                                   p30b = list_plant_FG,
                                   p31a = list_pollinator_FG,
                                   p31b = list_pollinator_FG,
                                   p31c = list_pollinator_FG) %>% mutate(motif = motif_id_i)

motif_combinations_prob <- bind_rows(motif_combinations_prob,
                                     get_motif_probabilities(motif_id_i_probability))

#-----------------
motif_id_i <- 13
print_positions(motif_id_i)
motif_id_i_probability <- crossing(p33 = list_plant_FG,
                                   p32a = list_plant_FG,
                                   p32b = list_plant_FG,
                                   p34 = list_pollinator_FG,
                                   p35 = list_pollinator_FG) %>% mutate(motif = motif_id_i)

motif_combinations_prob <- bind_rows(motif_combinations_prob,
                                     get_motif_probabilities(motif_id_i_probability))
#-----------------
motif_id_i <- 14
print_positions(motif_id_i)
motif_id_i_probability <- crossing(p36a = list_plant_FG,
                                   p36b = list_plant_FG,
                                   p37 = list_plant_FG,
                                   p38a = list_pollinator_FG,
                                   p38b = list_pollinator_FG) %>% mutate(motif = motif_id_i)

motif_combinations_prob <- bind_rows(motif_combinations_prob,
                                     get_motif_probabilities(motif_id_i_probability))

#-----------------
motif_id_i <- 15
print_positions(motif_id_i)
motif_id_i_probability <- crossing(p39 = list_plant_FG,
                                   p40a = list_plant_FG,
                                   p40b = list_plant_FG,
                                   p41 = list_pollinator_FG,
                                   p42 = list_pollinator_FG) %>% mutate(motif = motif_id_i)

motif_combinations_prob <- bind_rows(motif_combinations_prob,
                                     get_motif_probabilities(motif_id_i_probability))

#-----------------
motif_id_i <- 16
print_positions(motif_id_i)
motif_id_i_probability <- crossing(p43a = list_plant_FG,
                                   p43b = list_plant_FG,
                                   p43c = list_plant_FG,
                                   p44a = list_pollinator_FG,
                                   p44b = list_pollinator_FG) %>% mutate(motif = motif_id_i)

motif_combinations_prob <- bind_rows(motif_combinations_prob,
                                     get_motif_probabilities(motif_id_i_probability))

#-----------------
motif_id_i <- 17
print_positions(motif_id_i)
motif_id_i_probability <- crossing(p45a = list_plant_FG,
                                   p45b = list_plant_FG,
                                   p45c = list_plant_FG,
                                   p45d = list_plant_FG,
                                   p46 = list_pollinator_FG) %>% mutate(motif = motif_id_i)

motif_combinations_prob <- bind_rows(motif_combinations_prob,
                                     get_motif_probabilities(motif_id_i_probability))

# Save results
write_csv(motif_combinations_prob,"Data/Csv/node_motifs_theoretical_probability.csv")

# Number of motif combinations
motif_combinations_prob %>% 
  select(motif,motif_combination,motif_functional_ID,motif_probability) %>%
  unique()
  count() %>% filter(n>1)