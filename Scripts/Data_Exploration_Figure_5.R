

library(tidyverse)
library(RColorBrewer)

# Data with pollinator/plant labels

plant_means_reordered <- read_csv("Data/Csv/plant_abs_freq_means.csv")
pollinator_means_reordered <- read_csv("Data/Csv/pollinator_abs_freq_means.csv")
#Filter out in order to have just 6 FG'S
pollinator_means_reordered <- pollinator_means_reordered %>% filter(!Node_FG %in% c("Birds", "Lizards", "Other_insects"))


####################################################################################################
#Read data to check top 5 motif of each category of observed propabilities

#Under, No signi, Over-rep
####################################################################################################
# Check top Z-scores of each group (Under, No diff, Over)

#Read data
data <- read.csv("Data/Csv/motifs_observed_probability_SIMUL_CI.csv")

#Find critical value of Z-score
p <- 0.05 #cutoff probability 95% confidence
critical_value <- qnorm(p/2) #double tail probability divide by 2

data_1 <- data %>%
  mutate(infra_over_represented = case_when(
    z_score < -abs(critical_value) ~ "infra",
    between(z_score, -abs(critical_value), abs(critical_value)) ~ "no_diff",
    z_score > abs(critical_value) ~ "over"
  ))


#data_2 <- filter(data_1, abs(z_score)<2000) 
#str(data_2)


#Add number of nodes of each motif
data_1$N_node
data_1$N_node[data_1$motif=="1"] <- 2 # 2 SPECIES
data_1$N_node[data_1$motif=="2"] <- 3 # 3 SPECIES
data_1$N_node[data_1$motif=="3"] <- 3
data_1$N_node[data_1$motif=="4"] <- 4 # 4 SPECIES
data_1$N_node[data_1$motif=="5"] <- 4
data_1$N_node[data_1$motif=="6"] <- 4
data_1$N_node[data_1$motif=="7"] <- 4
data_1$N_node[data_1$motif=="8"] <- 5 # 5 SPECIES 
data_1$N_node[data_1$motif=="9"] <- 5
data_1$N_node[data_1$motif=="10"] <- 5
data_1$N_node[data_1$motif=="11"] <- 5
data_1$N_node[data_1$motif=="12"] <- 5
data_1$N_node[data_1$motif=="13"] <- 5
data_1$N_node[data_1$motif=="14"] <- 5
data_1$N_node[data_1$motif=="15"] <- 5
data_1$N_node[data_1$motif=="16"] <- 5
data_1$N_node[data_1$motif=="17"] <- 5

data_1 %>%                              
  group_by(infra_over_represented) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

#Select under-represented by the ones with highest porbability of being observed (top10)
under_data <- filter(data_1, infra_over_represented=="infra") 

under_top_5 <- filter(under_data, row_number(desc(abs(motif_observed_probability)))<= 5) %>%
  select(motif_functional_ID,motif_observed_probability,counts_observed,z_score)

#Select no statistical difference by the ones with highest porbability of being observed (top5)
no_diff_data <- filter(data_1, infra_over_represented=="no_diff") 

no_diff_top_5 <- filter(no_diff_data, row_number(desc(abs(motif_observed_probability)))<= 5) %>%
  select(motif_functional_ID,motif_observed_probability,counts_observed,z_score)

#Select over-represenetd by the ones with highest porbability of being observed (top5)
over_data <- filter(data_1, infra_over_represented=="over") 

over_top_5 <- filter(over_data, row_number(desc(abs(motif_observed_probability)))<= 5) %>%
  select(motif_functional_ID,motif_observed_probability,counts_observed,z_score)


#Now check the percentage of ach 2/3/4/5 species nodes on each of the 3 categories

under_data %>%                              
  group_by(N_node) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

no_diff_data %>%                              
  group_by(N_node) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

over_data %>%                              
  group_by(N_node) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
