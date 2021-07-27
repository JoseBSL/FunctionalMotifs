
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

data_2 <- filter(data_1, abs(z_score)<2000) 
str(data_2)
#Plot histogram of Z-scores of those motifs with observed probability > 0

# Justificación: si son observados 0 (o incluso 1 vez) no es muy interesante lo que les pase.
# Podrian estar underrepresented, pero como con la Chi, si el expected es super bajo, 
# diferenciar si esta under-represented sera dificil de detectar .

ggplot(data_2 %>% filter(round_motif_observed_probability>0), aes(x=z_score, color=infra_over_represented, fill=infra_over_represented)) + 
  geom_histogram(bins = 51, alpha = 0.5, position = "identity",lwd = 0.25)+
  geom_vline(xintercept = -2)+
  geom_vline(xintercept = 2)+
  xlim(-20,20) + ylab("Frequency")+  theme_bw() +
  scale_fill_manual(name="Motif frequencies" ,values=c("coral2", "palegreen3", "cyan3"), labels=c("Under-represented",
                    "No statistical difference", "Over-represented")) +
  scale_color_manual(name="Motif frequencies" ,values=c("coral2", "palegreen3", "cyan3"), labels=c("Under-represented",
                     "No statistical difference", "Over-represented"))

