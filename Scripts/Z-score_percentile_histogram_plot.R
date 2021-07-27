
#load library
library(dplyr)
library(ggplot2)
library(tidyverse)
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
str(data_2)
#Plot histogram of Z-scores
ggplot(data_2, aes(x=z_score, color=infra_over_represented, fill=infra_over_represented)) + 
  geom_histogram(bins = 50, alpha = 0.5, position = "identity",lwd = 0.25)+
  geom_vline(xintercept = -2)+
  geom_vline(xintercept = 2)+
  xlim(-20,20) + ylab("Frequency")+  theme_bw() +
  scale_fill_manual(name="Motif frequencies" ,values=c("coral2", "palegreen3", "cyan3"), labels=c("Under-represented",
                    "No statistical difference", "Over-represented")) +
  scale_color_manual(name="Motif frequencies" ,values=c("coral2", "palegreen3", "cyan3"), labels=c("Under-represented",
                     "No statistical difference", "Over-represented"))


###############
# PERCENTILES
###############

data_perc <- as.tibble(data)

# There are 827 motifs with upper and lower CI equal to zero
data_perc %>% filter(lower_CI == upper_CI)
data_perc %>% filter(lower_CI == upper_CI, lower_CI>0)

data_perc$infra_over_represented <- "no_diff"
data_perc$infra_over_represented[data_perc$round_motif_observed_probability<data_perc$lower_CI] <- "infra"
data_perc$infra_over_represented[data_perc$round_motif_observed_probability>data_perc$upper_CI] <- "over"

# There are 19514 motifs whose observed probability (0) is expected and equal to the lower bound (0)
data_perc %>% filter(percentil_observed == 0, infra_over_represented == "no_diff")

min(data_perc$percentil_observed[data_perc$infra_over_represented=="infra"])
max(data_perc$percentil_observed[data_perc$infra_over_represented=="infra"])

min(data_perc$percentil_observed[data_perc$infra_over_represented=="over"])
max(data_perc$percentil_observed[data_perc$infra_over_represented=="over"])

min(data_perc$percentil_observed[data_perc$infra_over_represented=="no_diff"])
max(data_perc$percentil_observed[data_perc$infra_over_represented=="no_diff"])

data_perc %>% group_by(infra_over_represented) %>% count()

#Plot histogram of percentil
ggplot(data, aes(x=percentil_observed)) + 
  geom_histogram(color="black", fill="white",bins = 100)+
  geom_vline(xintercept = 0.025)+
  geom_vline(xintercept = 0.975)

str(data)
