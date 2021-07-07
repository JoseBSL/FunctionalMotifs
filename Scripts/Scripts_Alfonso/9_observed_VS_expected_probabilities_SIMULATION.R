
library(tidyverse)

motifs_observed_probability <- read_csv("Data/Csv/motifs_observed_probability.csv")
motifs_expected_probability <- read_csv("Data/Csv/node_motifs_theoretical_probability.csv") %>%
  select(motif,motif_functional_ID,motif_probability) %>% unique() %>%
  rename(motif_expected_probability = motif_probability)

motifs_probability <- motifs_expected_probability %>% 
  left_join(motifs_observed_probability, by = c("motif","motif_functional_ID"))

motifs_probability[is.na(motifs_probability)] <- 0

motifs_probability <- motifs_probability %>% arrange(desc(motif_observed_probability))

#################################################
# RANDOM SIMULATION OF THE EXPECTED PROBABILITIES
#################################################

size = 1e7

simulations <- 1000

samples_random_probability <- NULL

for(i in 1:round(simulations,0)){
  
  random_i <- sample(motifs_probability$motif_functional_ID, 
                     size, replace = TRUE, 
                     prob = motifs_probability$motif_expected_probability) %>%
    as_tibble() %>% group_by(value) %>% count() %>% ungroup %>%
    mutate(random_probability = n/size,sample = i) %>% rename(motif_functional_ID = value) %>%
    dplyr::select(motif_functional_ID,sample,random_probability)
  
  samples_random_probability <- bind_rows(samples_random_probability,random_i)
}

write_csv(samples_random_probability,"Data/Csv/motifs_samples_random_probability_SIMU.csv")

#################################################
# EXTRACT CIs FROM SIMULATED PROBABILITIES
#################################################

# Variable to storage CI 95%
motifs_observed_probability_CI <- motifs_probability
motifs_observed_probability_CI$round_motif_observed_probability <- 
  round(motifs_observed_probability_CI$motif_observed_probability, 7) #1e-7 is the maximum resolution we can observe
motifs_observed_probability_CI$percentil_observed <- NA
motifs_observed_probability_CI$lower_CI <- NA
motifs_observed_probability_CI$upper_CI <- NA

motifs_ID_list <- motifs_probability$motif_functional_ID

# CI bounds from simulated probabilities and percentile of observed probabilities

for(i.motif in 1:length(motifs_ID_list)){
  
  random_samples_motif <- samples_random_probability %>% ungroup() %>%
    filter(motif_functional_ID == motifs_ID_list[i.motif]) %>% 
    dplyr::select(random_probability) %>% pull()
  
  if(length(random_samples_motif) < simulations){
    
    random_samples_motif <- c(random_samples_motif,
                              rep(0, (simulations-length(random_samples_motif))))
    
  }
  
  CI <- quantile(random_samples_motif, prob=c(0.025,0.975)) %>% unname()
  
  motifs_observed_probability_CI$percentil_observed[i.motif] <- 
    sum(random_samples_motif < motifs_observed_probability_CI$round_motif_observed_probability[i.motif])/length(random_samples_motif)
  motifs_observed_probability_CI$lower_CI[i.motif] <- CI[1]
  motifs_observed_probability_CI$upper_CI[i.motif] <- CI[2]
  
}

write_csv(motifs_observed_probability_CI,"Data/Csv/motifs_observed_probability_SIMUL_CI.csv")

underrepresented_motifs <-  motifs_observed_probability_CI %>% filter(round_motif_observed_probability<lower_CI)
overrepresented_motifs <-  motifs_observed_probability_CI %>% filter(round_motif_observed_probability>upper_CI)

write_csv(underrepresented_motifs,"Data/Csv/underrepresented_motifs.csv")
write_csv(overrepresented_motifs,"Data/Csv/overrepresented_motifs.csv")
