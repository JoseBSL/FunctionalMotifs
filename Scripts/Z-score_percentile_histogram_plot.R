
#load library
library(dplyr)
library(ggplot2)
#Read data
data <- read.csv("Data/Csv/motifs_observed_probability_SIMUL_CI.csv")

#Check NA's per column
data %>%
  summarise_all(funs(sum(is.na(.))))

# At the moment this is just -2 and 2 but check how to find critical value of z-score
data_1 <- data %>%
  mutate(infra_over_represented = case_when(
    z_score < -2 ~ "infra",
    between(z_score, -2, 2) ~ "no_diff",
    z_score > 2 ~ "over"
  ))
     
#check levels
levels(factor(data_1$infra_over_represented))

#Check proportion of infra/over and no statistical difference
data_1 %>% 
  group_by(infra_over_represented) %>%
  summarise(no_rows = length(infra_over_represented)) %>%
  mutate(proportion = no_rows / sum(no_rows))

data_2 <- filter(data_1, z_score<3000) 

#Plot histogram of Z-scores
ggplot(data_2, aes(x=z_score)) + 
  geom_histogram(color="black", fill="white")



