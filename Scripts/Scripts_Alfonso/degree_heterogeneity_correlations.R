
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
library(bipartite)
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

VAR_fun <- function(network_name){
  
  degree_sequence <- degree_sequence_from_edgelist(network_name)
  
  return(var(degree_sequence))
  
}

heterogeneity_index_fun <- function(network_name){
  
  networks_i <- networks %>% filter(Network_id == network_name)
  
  networks_i <- networks_i %>% 
    group_by(Pollinator_species,Plant_species,
             Pollinator_functional_group,Plant_functional_groups) %>%
    count(wt = Interaction) %>% rename(Interaction = n)
  
  edge_list_i <- networks_i[,c("Pollinator_species","Plant_species")]
  
  plant_degree <- edge_list_i %>% ungroup() %>% dplyr::select(Plant_species) %>% 
    dplyr::group_by(Plant_species) %>%
    count() %>% rename(degree_Plant = n)
  
  plant_nodes <- nrow(plant_degree)
  
  pollinator_degree <- edge_list_i %>% ungroup() %>% dplyr::select(Pollinator_species) %>% 
    dplyr::group_by(Pollinator_species) %>%
    count() %>% rename(degree_Pollinator = n)
  
  pollinator_nodes <- nrow(pollinator_degree)
  
  total_nodes <- plant_nodes + pollinator_nodes
  
  edge_list_i_final <- edge_list_i %>% 
    left_join(pollinator_degree, by = "Pollinator_species") %>%
    left_join(plant_degree, by = "Plant_species") %>%
    mutate(fun_degr_Poll = 1/sqrt(degree_Pollinator),
           fun_degr_Plant = 1/sqrt(degree_Plant),
           irregularity = ( fun_degr_Poll - fun_degr_Plant )^2)

  return(sum(edge_list_i_final$irregularity)/(total_nodes-2*sqrt(total_nodes-1)))
  
}


#Calculate metrics
connectance = lapply(my_data, function(i) {networklevel(as.matrix(i), index="connectance")})
degree_dist = lapply(my_data, function(i) {networklevel(as.matrix(i), index="links per species")})
nestedness = lapply(my_data, function(i) {networklevel(as.matrix(i), index="nestedness")})
gini_coefficients = lapply(my_network_ids, function(i) {gini_coeff_fun(i)})
heterogeneity_index = lapply(my_network_ids, function(i) {heterogeneity_index_fun(i)})
VAR = lapply(my_network_ids, function(i) {VAR_fun(i)})

#Connectance/nestedness
connectance_d = bind_rows(connectance)
nestedness_d = bind_rows(nestedness)
degree_dist_d = bind_rows(degree_dist)

gini_dist_d = as.numeric(gini_coefficients)
heterogeneity_index_dist_d = as.numeric(heterogeneity_index)
VAR_dist_d = as.numeric(VAR)

cor.test(connectance_d$connectance, nestedness_d$nestedness, method = "spearman")
cor.test(connectance_d$connectance, degree_dist_d$`links per species`, method = "spearman")
cor.test(connectance_d$connectance, nestedness_d$nestedness, method = "spearman")
cor.test(connectance_d$connectance, gini_dist_d, method = "spearman")
cor.test(connectance_d$connectance, heterogeneity_index_dist_d, method = "spearman")
cor.test(connectance_d$connectance, VAR_dist_d, method = "spearman")


#Bind cols
#1st create datframe of interest
gini_dist_d = data.frame(gini_index=as.numeric(gini_coefficients))
#Now bind cols in one 
d = bind_cols(connectance_d,nestedness_d, gini_dist_d)

d %>% cor(method ="spearman") %>%
  ggcorrplot::ggcorrplot(lab = TRUE, type = "lower") +
  ggtitle("Network metric correlations")

library(ggstatsplot)

ggcorrmat(d, type="nonparametric")+
  ggtitle("Network metric correlations")

