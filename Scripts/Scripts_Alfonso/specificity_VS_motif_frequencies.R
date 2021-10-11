
library(lme4)
library(nlme)
library(stringi)
library(tidyverse)
source("Scripts/Scripts_Alfonso/add_study_id.R")
d <- read_csv("Data/Csv/Motifs_frequencies_and_null_models/Motifs_frequency_percentile.csv")
d <- add_study_id(d)


# Mean percentile and SE taking into account the study system
motif_codes <- d %>% filter(motif != 1) %>% # We remove links (motif code: 1)
  select(motif) %>% unique() %>% pull()
motif_means <- tibble(motif = motif_codes)
motif_means$mean <- NA
motif_means$SE <- NA

for(i.motif in 1:length(motif_codes)){
  
  m <- lmer(percentil_sizeclass ~ 1+(1|study_id), #HERE WE NEED STUDY SYSTEM ONLY!! 
            data = subset(d, motif == motif_codes[i.motif]))  
  
  motif_means$mean[i.motif] <- fixed.effects(m) %>% unname()
  motif_means$SE[i.motif] <- sqrt(diag(vcov(m)))
  
}

###########################
# ESTIMATING SPECIFICITY
###########################
# Reference: 10.1371/journal.pone.0114674

motif_pattern_connections <- read_csv("Data/Data_processing/Motifs_connections/motif_pattern_connections.csv")

number_plant <- motif_pattern_connections %>% group_by(motif_id,plant) %>% 
  count() %>% group_by(motif_id) %>% count() %>% rename(number_plants = n)

number_pollinator <- motif_pattern_connections %>% group_by(motif_id,pollinator) %>% 
  count() %>% group_by(motif_id) %>% count() %>% rename(number_pollinators = n)

specificity_plant_aux <- motif_pattern_connections %>% 
  group_by(motif_id,plant) %>% count()

specificity_pollinator_aux <- motif_pattern_connections %>% 
  group_by(motif_id,pollinator) %>% count()

specificity_plant_aux$position <- stri_extract_first_regex(specificity_plant_aux$plant,
                                                           "[0-9]+") %>% as.numeric()
specificity_pollinator_aux$position <- stri_extract_first_regex(specificity_pollinator_aux$pollinator,
                                                                "[0-9]+") %>% as.numeric()

specificity_plant_aux2 <- specificity_plant_aux %>% 
  left_join(number_pollinator, by = "motif_id")

specificity_pollinator_aux2 <- specificity_pollinator_aux %>% 
  left_join(number_plant, by = "motif_id")

specificity_plant <- specificity_plant_aux2 %>% ungroup %>%
  mutate(s=(number_pollinators-n)/(number_pollinators-1)) %>%
  select(position,s) %>% unique()

specificity_pollinator <- specificity_pollinator_aux2 %>% ungroup %>%
  mutate(s=(number_plants-n)/(number_plants-1)) %>%
  select(position,s) %>% unique()

# We replace NaN of some specialist positions by one
specificity_plant$s[is.nan(specificity_plant$s)] <- 1.0
specificity_pollinator$s[is.nan(specificity_pollinator$s)] <- 1.0

plant_motif_pattern_connections <- motif_pattern_connections %>%
  select(motif_id,plant) %>% unique()

plant_motif_pattern_connections$position <- stri_extract_first_regex(plant_motif_pattern_connections$plant,
                                                                     "[0-9]+") %>% as.numeric()

pollinator_motif_pattern_connections <- motif_pattern_connections %>%
  select(motif_id,pollinator) %>% unique()

pollinator_motif_pattern_connections$position <- 
  stri_extract_first_regex(pollinator_motif_pattern_connections$pollinator,
                           "[0-9]+") %>% as.numeric()

plant_motif_s <- plant_motif_pattern_connections %>%
  left_join(specificity_plant, by = "position")

mean_plant_motif_s <- plant_motif_s %>% ungroup() %>% group_by(motif_id) %>%
  summarise(mean_s_plant = mean(s))

pollinator_motif_s <- pollinator_motif_pattern_connections %>%
  left_join(specificity_pollinator, by = "position")

mean_pollinator_motif_s <- pollinator_motif_s %>% ungroup() %>% 
  group_by(motif_id) %>%
  summarise(mean_s_pollinator = mean(s))


motif_means_s <- motif_means %>% rename(motif_id = motif) %>%
  left_join(mean_plant_motif_s, by = "motif_id") %>%
  left_join(mean_pollinator_motif_s, by = "motif_id") %>%
  rename(mean_frequency = mean)

# library(plotly)
# plot_ly(x=motif_means_s$mean_s_plant, 
#         y=motif_means_s$mean_s_pollinator, 
#         z=motif_means_s$mean_frequency, type="scatter3d", mode="markers")

library(glmmTMB)
m_lm <- glmmTMB(mean_frequency~mean_s_plant*mean_s_pollinator, data = motif_means_s, family=beta_family())
#m_lm <- lm(mean_frequency~mean_s_plant*mean_s_pollinator,motif_means_s)
summary(m_lm)

library(betareg)
m_lm2 <- betareg(mean_frequency~mean_s_plant*mean_s_pollinator, data = motif_means_s, link = "logit")
#m_lm <- lm(mean_frequency~mean_s_plant*mean_s_pollinator,motif_means_s)
summary(m_lm2)

library(performance)
r2_efron(m_lm)
r2(m_lm2)

#0.67

library(visreg)
visreg2d(m_lm, "mean_s_plant", "mean_s_pollinator",scale ="response")

motif_means_s$predicted_freq <- predict(m_lm,motif_means_s,type ="response")

ggplot(motif_means_s,aes(x=mean_frequency, 
                        y = predicted_freq))+
  #geom_point()+
  geom_abline(slope=1,intercept = 0)+
  geom_text(label=motif_means_s$motif_id)+
  theme_bw()+
  xlim(0,1)+ylim(0,1)
