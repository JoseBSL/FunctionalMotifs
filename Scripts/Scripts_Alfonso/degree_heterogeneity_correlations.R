
########################################################################################################################################################
#1) LOAD NETWORK DATA
########################################################################################################################################################
#Set working directory to read files
dir_ini <- getwd()
setwd("Data/Data_networks_quantitative")
temp <- list.files(pattern="*.csv")
my.list <- list(for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i])))
my_files <- list.files(pattern = "\\.csv$")
my_data <- lapply(my_files, function(i){read.csv(i, check.names=FALSE, row.names = 1)})
my_network_ids <- gsub(".csv", "", my_files)

library(igraph)
library(tidyverse)
library(DescTools)

setwd(dir_ini)
int.threshold <- 1
networks <- read_csv("Data/Csv/data_for_motifs_analysis_1.csv") %>% select(-X1) %>%
  filter(Interaction >= int.threshold)

degree_sequence_from_edgelist <- function(network_name){
  
  networks_i <- networks %>% filter(Network_id == network_name)
  
  networks_i <- networks_i %>% 
    group_by(Pollinator_species,Plant_species,
             Pollinator_functional_group,Plant_functional_groups) %>%
    count(wt = Interaction) %>% rename(Interaction = n)
  
  edge_list_i <- networks_i[,c("Pollinator_species","Plant_species")]
  
  g_i <- igraph::graph_from_edgelist(as.matrix(edge_list_i), directed = F)
  
  degree_sequence <- igraph::degree(g_i) %>% as.numeric()
  
  return(degree_sequence)
}

gini_coeff_fun <- function(network_name){
  
  degree_sequence <- degree_sequence_from_edgelist(network_name)
  Gini_network <- DescTools::Gini(degree_sequence, unbiased=FALSE)
  
  return(Gini_network)
  
}

randic_index_fun <- function(network_name){
  
  degree_sequence <- degree_sequence_from_edgelist(network_name)
  output <- 0
  
  for (i in 1:length(degree_sequence)) {
    
    for (j in 1:length(degree_sequence)) {
      
      output <- output + (degree_sequence[i]*degree_sequence[j])^(-0.5)
      
    }
    
  }
  
  return(output)
  
}

#Calculate metrics
connectance = lapply(my_data, function(i) {networklevel(as.matrix(i), index="connectance")})
degree_dist = lapply(my_data, function(i) {networklevel(as.matrix(i), index="links per species")})
nestedness = lapply(my_data, function(i) {networklevel(as.matrix(i), index="nestedness")})
gini_coefficients = lapply(my_network_ids, function(i) {gini_coeff_fun(i)})
randic_index = lapply(my_network_ids, function(i) {randic_index_fun(i)})

#Connectance/nestedness
connectance_d = bind_rows(connectance)
nestedness_d = bind_rows(nestedness)
degree_dist_d = bind_rows(degree_dist)
gini_dist_d = as.numeric(gini_coefficients)
randic_index_d = as.numeric(randic_index)

cor.test(connectance_d$connectance, nestedness_d$nestedness, method = "spearman")
cor.test(connectance_d$connectance, degree_dist_d$`links per species`, method = "spearman")
cor.test(connectance_d$connectance, nestedness_d$nestedness, method = "spearman")
cor.test(connectance_d$connectance, gini_dist_d, method = "spearman")
cor.test(connectance_d$connectance, randic_index_d, method = "spearman")
