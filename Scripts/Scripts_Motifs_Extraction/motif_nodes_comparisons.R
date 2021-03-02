
# Explore the effect of network size (number of nodes) in the number of motifs


# Load libraries
library(tidyverse)
require(scales)

# Load raw data 

# Networks
int.threshold <- 1

networks <- read_csv("Data/Csv/data_for_motifs_analysis.csv") %>% select(-X1) %>%
  filter(Interaction >= int.threshold)

list_Network_id <- networks$Network_id %>% unique()

topology <- tibble(Network_id = list_Network_id)
topology$Nodes <- NA
topology$Plants <- NA
topology$Pollinators <- NA
topology$Connectance <- NA

for (i.network in 1:length(list_Network_id)){#1:length(list_Network_id)){

  # Create a graph for the i-th network
  networks_i <- networks %>% filter(Network_id == list_Network_id[i.network])
  edge_list_i <- networks_i[,c("Pollinator_species","Plant_species")]
  
  topology$Plants[i.network] <- edge_list_i$Plant_species %>% unique() %>% length()
  topology$Pollinators[i.network] <- edge_list_i$Pollinator_species %>% unique() %>% length()
  topology$Nodes[i.network] <- topology$Plants[i.network] + topology$Pollinators[i.network]
  topology$Connectance[i.network] <- 
    nrow(edge_list_i) / (topology$Plants[i.network] * topology$Pollinators[i.network])
}

ggplot(topology, aes(x = Connectance, y = Nodes))+
  geom_point() + geom_smooth(method = "loess")

ggplot(topology, aes(x = Plants, y = Nodes))+
  geom_point() + geom_smooth(method = "loess")

ggplot(topology, aes(x = Pollinators, y = Nodes))+
  geom_point() + geom_smooth(method = "loess")

ggplot(topology, aes(x = Pollinators, y = Connectance))+
  geom_point() + geom_smooth(method = "loess")

# Motif frequencies
motifs_raw <- read_csv("Data/Csv/network_frequency_motifs.csv")

total_motifs <- motifs_raw %>% group_by(Network_id) %>% count(wt = frequency) %>%
  rename(total_motifs = n)

motifs_3_nodes <- motifs_raw %>% filter(nodes==3) %>% group_by(Network_id) %>% 
  count(wt = frequency) %>%
  rename(motifs_3_nodes = n)

motifs_4_nodes <- motifs_raw %>% filter(nodes==4) %>% group_by(Network_id) %>% 
  count(wt = frequency) %>%
  rename(motifs_4_nodes = n)

motifs_5_nodes <- motifs_raw %>% filter(nodes==5) %>% group_by(Network_id) %>% 
  count(wt = frequency) %>%
  rename(motifs_5_nodes = n)


topology_motifs <- topology %>% left_join(total_motifs, by = "Network_id") %>%
  left_join(motifs_3_nodes, by = "Network_id") %>%
  left_join(motifs_4_nodes, by = "Network_id") %>%
  left_join(motifs_5_nodes, by = "Network_id")

ggplot(topology_motifs %>% filter(Connectance > 0),
       aes(x = Pollinators*Plants, y = total_motifs, color = Connectance))+
  geom_point() + geom_smooth(method = "lm")+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  xlab("Number of pollinators x Number of plants")+
  ylab("Number of (up to 5) motifs")+
  theme_bw()

m1 <- lm(log10(Pollinators*Plants) ~ log10(total_motifs),
         topology_motifs %>% filter(Connectance > 0))
summary(m1)

ggplot(topology_motifs %>% filter(Nodes < 75),
       aes(x = Nodes, y = log10(total_motifs)))+
  geom_point() + geom_smooth(method = "loess")

ggplot(topology_motifs,
       aes(x = Connectance, y = log10(total_motifs)))+
  geom_point() + geom_smooth(method = "loess")

ggplot(topology_motifs,
       aes(x = Nodes, y = log10(motifs_3_nodes)))+
  geom_point() + geom_smooth(method = "loess")

ggplot(topology_motifs,
       aes(x = Nodes, y = log10(motifs_4_nodes)))+
  geom_point() + geom_smooth(method = "loess")

ggplot(topology_motifs,
       aes(x = Nodes, y = log10(motifs_3_nodes)))+
  geom_point() + geom_smooth(method = "loess")
