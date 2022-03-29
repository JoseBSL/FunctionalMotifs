
# Load libraries

library(lme4)
library(nlme)
library(tidyverse)
library(stringi)
source("Scripts/Scripts_Alfonso/add_study_id.R")


pollinator_positions_codes <- c("np2","np4","np6","np8","np11","np12","np14","np16",
                                "np18","np21","np22","np24","np25","np28","np29","np31",
                                "np34","np35","np38","np41","np42","np44","np46")
plant_positions_codes <- c("np1","np3","np5","np7","np9","np10","np13","np15",
                                "np17","np19","np20","np23","np26","np27","np30","np32",
                                "np33","np36","np37","np39","np40","np43","np45")


#############################################
# EXTRACT MEAN VALUES
#############################################
all_position_percentiles <- read_csv("Data/Csv/Motifs_positions_and_null_models/GF_positions_frequency_percentile.csv")

all_position_percentiles$position <- stri_extract_first_regex(all_position_percentiles$position,
                                                              "[0-9]+") %>% as.numeric()

all_position_percentiles <- add_study_id(all_position_percentiles)
#str(all_position_percentiles)
#head(all_position_percentiles)

list_plants_FG <- as.character(1:10)

plant_position_percentiles_filtered <- all_position_percentiles %>% 
  filter(Node_FG %in% list_plants_FG, !position %in% c(1,2),
         !is.na(expected_freq_its_GF)) #ME DA ERROR
pollinator_position_percentiles_filtered <- all_position_percentiles %>% 
  filter(!Node_FG %in% list_plants_FG, !position %in% c(1,2),
         !is.na(expected_freq_its_GF))


# Plants-------------

# Mean percentile and SE taking into account the study system
plant_codes <- plant_position_percentiles_filtered$position %>% unique()
plant_means <- plant_position_percentiles_filtered %>% select(position,Node_FG) %>% unique()
plant_means$mean <- NA
plant_means$SE <- NA


for(i.pos in 1:length(plant_codes)){
  
  m <- lmer(percentil_its_GF ~ Node_FG + (1|study_id),
            data = plant_position_percentiles_filtered %>%
              filter(position == plant_codes[i.pos]))  
  
  temp <- fixed.effects(m) %>% unname()
  plant_means$mean[plant_means$position == plant_codes[i.pos]] <- c(temp[1],
                                                                    temp[1]+temp[2],
                                                                    temp[1]+temp[3],
                                                                    temp[1]+temp[4],
                                                                    temp[1]+temp[5]) ## to get real means per FG
  temp <- sqrt(diag(vcov(m)))
  plant_means$SE[plant_means$position == plant_codes[i.pos]] <- c(temp[1], 
                                                                  temp[1]+temp[2],
                                                                  temp[1]+temp[3],
                                                                  temp[1]+temp[4],
                                                                  temp[1]+temp[5])
  
}


# Pollinators--------

# Mean percentile and SE taking into account the study system
pollinator_codes <- pollinator_position_percentiles_filtered$position %>% unique()
pollinator_means <- pollinator_position_percentiles_filtered %>% select(position,Node_FG) %>% unique()
pollinator_means$mean <- NA
pollinator_means$SE <- NA
pollinator_means_reordered <- NULL

for(i.pos in 1:length(pollinator_codes)){
  
  m <- lmer(percentil_its_GF ~ Node_FG + (1|study_id), #HERE WE NEED STUDY SYSTEM ONLY!! 
            data = pollinator_position_percentiles_filtered %>%
              filter(position == pollinator_codes[i.pos]),
            control = lmerControl(optimizer ="Nelder_Mead")) # I changed the optimizer because motif 8 (i.pos = 3) caused convergence problems, max gradient = 0.0022 > 0.002.
  
  temp <- fixed.effects(m) %>% unname()
  pollinator_means_aux <- pollinator_means[pollinator_means$position == pollinator_codes[i.pos],] %>% arrange(Node_FG)
  
  pollinator_means_aux$mean <- c(temp[1],
                                 temp[1]+temp[2],
                                 temp[1]+temp[3],
                                 temp[1]+temp[4],
                                 temp[1]+temp[5],
                                 temp[1]+temp[6],
                                 temp[1]+temp[7],
                                 temp[1]+temp[8],
                                 temp[1]+temp[9]) ## to get real means per FG
  
  
  temp <- sqrt(diag(vcov(m)))
  pollinator_means_aux$SE <- c(temp[1],
                               temp[1]+temp[2],
                               temp[1]+temp[3],
                               temp[1]+temp[4],
                               temp[1]+temp[5],
                               temp[1]+temp[6],
                               temp[1]+temp[7],
                               temp[1]+temp[8],
                               temp[1]+temp[9])
  
  pollinator_means_reordered <- bind_rows(pollinator_means_reordered,pollinator_means_aux)
  
}

# There are lower limits (mean - SE) smaller than 0 (see Other insects)
# We set the lower limit of those error bars to zero

pollinator_means_reordered <- pollinator_means_reordered %>% mutate(lower = mean-SE)
pollinator_means_reordered$lower[pollinator_means_reordered$lower < 0] <- 0

#############################################
# ADDING BROAD MOTIF CATEGORIES TO POSITIONS
#############################################

motifs_raw_data <- read_csv("Data/Data_processing/Motifs_connections/motif_pattern_connections.csv")

motifs_raw_data$Broad_categories <- NA

motifs_raw_data$Broad_categories[motifs_raw_data$motif_id %in% c(5,9,14,13,10)] <- "Core-peripherical"
motifs_raw_data$Broad_categories[motifs_raw_data$motif_id %in% c(6,16,12)] <- "Complete"
motifs_raw_data$Broad_categories[motifs_raw_data$motif_id %in% c(2,3,7,4,17,8)] <- "Fan"
motifs_raw_data$Broad_categories[motifs_raw_data$motif_id %in% c(15,11)] <- "Asymmetric complete"

motifs_raw_data$plant <- gsub("[^0-9.-]", "", motifs_raw_data$plant)
motifs_raw_data$pollinator <- gsub("[^0-9.-]", "", motifs_raw_data$pollinator)
motifs_raw_data$plant <- as.numeric(motifs_raw_data$plant)
motifs_raw_data$pollinator <- as.numeric(motifs_raw_data$pollinator)
motifs_raw_data <- motifs_raw_data %>% filter(!is.na(Broad_categories))

fan_pollinator <- motifs_raw_data$pollinator[motifs_raw_data$Broad_categories == "Fan"] %>% unique()
core_per_pollinator <- motifs_raw_data$pollinator[motifs_raw_data$Broad_categories == "Core-peripherical"] %>% unique()
complete_pollinator <- motifs_raw_data$pollinator[motifs_raw_data$Broad_categories == "Complete"] %>% unique()
asymm_pollinator <- motifs_raw_data$pollinator[motifs_raw_data$Broad_categories == "Asymmetric complete"] %>% unique()

pollinator_means_reordered$Broad_categories <- NA

pollinator_means_reordered$Broad_categories[pollinator_means_reordered$position %in% fan_pollinator] <- "Fan" 
pollinator_means_reordered$Broad_categories[pollinator_means_reordered$position %in% asymm_pollinator] <- "Medium-weak" 
pollinator_means_reordered$Broad_categories[pollinator_means_reordered$position %in% complete_pollinator] <- "Strong" 
pollinator_means_reordered$Broad_categories[pollinator_means_reordered$position %in% core_per_pollinator] <- "Weak" 


fan_plant <- motifs_raw_data$plant[motifs_raw_data$Broad_categories == "Fan"] %>% unique()
core_per_plant <- motifs_raw_data$plant[motifs_raw_data$Broad_categories == "Core-peripherical"] %>% unique()
complete_plant <- motifs_raw_data$plant[motifs_raw_data$Broad_categories == "Complete"] %>% unique()
asymm_plant <- motifs_raw_data$plant[motifs_raw_data$Broad_categories == "Asymmetric complete"] %>% unique()

plant_means$Broad_categories <- NA

plant_means$Broad_categories[plant_means$position %in% fan_plant] <- "Fan" 
plant_means$Broad_categories[plant_means$position %in% asymm_plant] <- "Medium-weak" 
plant_means$Broad_categories[plant_means$position %in% complete_plant] <- "Strong" 
plant_means$Broad_categories[plant_means$position %in% core_per_plant] <- "Weak" 

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


plant_means_s <- plant_means %>% left_join(specificity_plant, by = "position")

pollinator_means_reordered_s <- pollinator_means_reordered %>% 
  left_join(specificity_pollinator, by = "position")


# We replace NaN of some specialist positions by one
plant_means_s$s[is.nan(plant_means_s$s)] <- 1.0
pollinator_means_reordered_s$s[is.nan(pollinator_means_reordered_s$s)] <- 1.0

library(patchwork)

#select poll functional groups of interest
levels(factor(pollinator_means_reordered_s$Node_FG))
vars <- c("Birds","Lizards", "Other_insects")
pollinator_means_reordered_s_1 <- pollinator_means_reordered_s %>% filter(!Node_FG %in% vars)
  
p1 <- ggplot(pollinator_means_reordered_s_1,aes(x=s,y=mean))+
  geom_point(alpha=0.5, color="black")+
  geom_smooth(method = "lm",color="black")+
  facet_wrap(~Node_FG)+
  labs(x="Specificty",y="Average frequency per position") +
  theme_bw() + ggtitle("Floral visitors") +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- ggplot(plant_means_s,aes(x=s,y=mean))+
  geom_point(alpha=0.5, color="black")+
  geom_smooth(method = "lm", color="black")+
  facet_wrap(~Node_FG)+
  labs(x="Specificty",y="Average frequency per position")+
  theme_bw() + ggtitle("Plants") +
  theme(plot.title = element_text(hjust = 0.5))

p1/p2

#save data 
write.csv(pollinator_means_reordered_s_1, "Data/Csv/pollinator_means_reordered_s_1.csv")
write.csv(plant_means_s, "Data/Csv/plant_means_s.csv")

library(glmmTMB)
library(performance)
m_lm_pos_poll <- glmmTMB(mean~s*Node_FG, data = pollinator_means_reordered_s_1,
                         family = beta_family())
summary(m_lm_pos_poll)
r2(m_lm_pos_poll)

m_lm_pos_poll_rd <- glmmTMB(mean ~ s + (1|Node_FG), data = pollinator_means_reordered_s_1,
                            family = beta_family())

summary(m_lm_pos_poll_rd)
r2(m_lm_pos_poll_rd)


library(visreg)
visreg2d(m_lm_pos_poll, "Node_FG","s",scale ="response")
visreg(m_lm_pos_poll,"s",scale ="response")
visreg(m_lm_pos_poll,"Node_FG",scale ="response")

visreg2d(m_lm_pos_poll_rd, "Node_FG","s",scale ="response")
visreg(m_lm_pos_poll_rd,"s",scale ="response")
visreg(m_lm_pos_poll_rd,"Node_FG",scale ="response")


#---------------------------------
# Percentil analises with LMMs for Pollinators


library(lme4)
library(glmmTMB)
library(performance)

pollinator_freq_position_s <- pollinator_position_percentiles_filtered %>% 
  left_join(specificity_pollinator, by = "position")

pollinator_freq_position_s$s[is.nan(pollinator_freq_position_s$s)] <- 1.0

# Indirect interaction data
motif_interactions <- read_csv("Data/Data_processing/Motifs_connections/motif_interactions.csv") %>%
  select(position,indirect_interactions) %>% unique()

pollinator_freq_position_s_ind <- pollinator_freq_position_s %>% 
  left_join(motif_interactions, by = "position") %>% 
  filter(!Node_FG %in% c("Birds","Lizards","Other_insects")) %>%
  rename(Percentil=percentil_its_GF,Specificity=s,Indirect_interactions=indirect_interactions,Group=Node_FG)

# pollinator_freq_position_s_ind$Percentil[pollinator_freq_position_s_ind$Percentil==0] <- 
#   1e-10
# 
# pollinator_freq_position_s_ind$Percentil[pollinator_freq_position_s_ind$Percentil==1] <- 
#   1-1e-10

pollinator_freq_position_s_ind$position <- as.factor(pollinator_freq_position_s_ind$position)

model_freq_specif_ind <- lmer(Percentil ~ Specificity*Group +
                                   scale(Indirect_interactions)*Group + 
                                   (1|study_id/position),
                           pollinator_freq_position_s_ind)


model_freq_specif_ind <- glmmTMB(Percentil ~ Specificity*Group +
                       scale(Indirect_interactions)*Group + 
                       (1|study_id/position),
                     pollinator_freq_position_s_ind)


#Save model data
saveRDS(model_freq_specif_ind, "Data/RData/model_freq_specif_ind.RData")
saveRDS(pollinator_freq_position_s_ind, "Data/RData/pollinator_freq_position_s_ind.RData")


summary(model_freq_specif_ind)
r2(model_freq_specif_ind)

library(visreg)
visreg2d(model_freq_specif_ind, "Group","Specificity",scale ="response")
visreg2d(model_freq_specif_ind, "Group","Indirect_interactions",scale ="response")

# Pollinator Model REAL OUTPUTS!!!
model_freq_specif_ind <- readRDS("Data/RData/model_freq_specif_ind.RData")

pollinators_sp_plot <- visreg(model_freq_specif_ind, "Specificity", by="Group",
                         strip.names=c("Bee",
                                       "Coleoptera",
                                       "Lepidoptera",
                                       "Non-bee-Hymenoptera",
                                       "Non-syrphids-diptera",
                                       "Syrphids"),gg = TRUE) + theme_bw()

pollinators_ind_int_plot <- visreg(model_freq_specif_ind, "Indirect_interactions", by="Group",
                              strip.names=c("Bee",
                                            "Coleoptera",
                                            "Lepidoptera",
                                            "Non-bee-Hymenoptera",
                                            "Non-syrphids-diptera",
                                            "Syrphids"),gg = TRUE)+ theme_bw()+xlab("Indirect interactions")



library(patchwork)
pollinators_sp_plot/pollinators_ind_int_plot

library("DHARMa")
pollinator_freq_position_s_ind$scale_Indirect_interactions <- scale(pollinator_freq_position_s_ind$Indirect_interactions)
model_freq_specif_ind2 <- lmer(Percentil ~ Specificity*Group +
                                 scale_Indirect_interactions*Group + 
                                (1|study_id/position),
                              pollinator_freq_position_s_ind)
simulationOutput <- simulateResiduals(fittedModel = model_freq_specif_ind2, n = 1050, use.u = T)
plot(simulationOutput)
testDispersion(simulationOutput)
plotResiduals(simulationOutput, form = pollinator_freq_position_s_ind$Group)
plotResiduals(simulationOutput, form = pollinator_freq_position_s_ind$Specificity)
plotResiduals(simulationOutput, form = pollinator_freq_position_s_ind$Indirect_interactions)
########################################3
#


############################
#Visualization of slopes by using visreg: This is not OK because I will use individual fits

library(ggeffects)
library(scales)
Bee_freq_position_s_ind <- pollinator_freq_position_s_ind %>% 
  filter(Group =="Bee")

model_Bee <- glmmTMB(Percentil ~ Specificity + 
                       Indirect_interactions + 
                       (1|study_id), Bee_freq_position_s_ind,
                     family = beta_family())

coleoptera_freq_position_s_ind <- pollinator_freq_position_s_ind %>% 
  filter(Group =="Coleoptera")

model_coleoptera <- glmmTMB(Percentil ~ Specificity + 
                       Indirect_interactions + 
                       (1|study_id), coleoptera_freq_position_s_ind,
                     family = beta_family())

Lepidoptera_freq_position_s_ind <- pollinator_freq_position_s_ind %>% 
  filter(Group =="Lepidoptera")
model_Lepidoptera <- glmmTMB(Percentil ~ Specificity + 
                       Indirect_interactions + 
                       (1|study_id), Lepidoptera_freq_position_s_ind,
                     family = beta_family())

NonbeeHymenoptera_freq_position_s_ind <- pollinator_freq_position_s_ind %>% 
  filter(Group =="Non-bee-Hymenoptera")
model_NonbeeHymenoptera <- glmmTMB(Percentil ~ Specificity + 
                               Indirect_interactions + 
                               (1|study_id), NonbeeHymenoptera_freq_position_s_ind,
                             family = beta_family())
Nonsyrphidsdiptera_freq_position_s_ind <- pollinator_freq_position_s_ind %>% 
  filter(Group =="Non-syrphids-diptera")
model_Nonsyrphidsdiptera <- glmmTMB(Percentil ~ Specificity + 
                                     Indirect_interactions + 
                                     (1|study_id), Nonsyrphidsdiptera_freq_position_s_ind,
                                   family = beta_family())

Syrphids_freq_position_s_ind <- pollinator_freq_position_s_ind %>% 
  filter(Group =="Syrphids")
model_Syrphids <- glmmTMB(Percentil ~ Specificity + 
                                      Indirect_interactions + 
                                      (1|study_id), Syrphids_freq_position_s_ind,
                                    family = beta_family())



# #####################################
# dev.off()
# par(mfrow = c(4,3),mar=c(3.95,4,1,1))
# 
# visreg(model_Bee,"Specificity",xlab="Specificity",ylab="Percentil",
#        main="Bee",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
# points(Percentil ~ Specificity, data = Bee_freq_position_s_ind, 
#        col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
#        pch = 20) 
# 
# visreg(model_coleoptera,"Specificity",xlab="Specificity",ylab="Percentil",
#        main="Coleoptera",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
# points(Percentil ~ Specificity, data = coleoptera_freq_position_s_ind, 
#        col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
#        pch = 20) 
# 
# visreg(model_Lepidoptera,"Specificity",xlab="Specificity",ylab="Percentil",
#        main="Lepidoptera",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
# points(Percentil ~ Specificity, data = Lepidoptera_freq_position_s_ind, 
#        col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
#        pch = 20) 
# 
# visreg(model_NonbeeHymenoptera,"Specificity",xlab="Specificity",ylab="Percentil",
#        main="Non-bee-Hymenoptera",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
# points(Percentil ~ Specificity, data = NonbeeHymenoptera_freq_position_s_ind, 
#        col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
#        pch = 20) 
# 
# visreg(model_Nonsyrphidsdiptera,"Specificity",xlab="Specificity",ylab="Percentil",
#        main="Non-syrphids-diptera",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
# points(Percentil ~ Specificity, data = Nonsyrphidsdiptera_freq_position_s_ind, 
#        col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
#        pch = 20) 
# 
# visreg(model_Syrphids,"Specificity",xlab="Specificity",ylab="Percentil",
#        main="Syrphids",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
# points(Percentil ~ Specificity, data = Syrphids_freq_position_s_ind, 
#        col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
#        pch = 20) 
# 
# 
# ##############
# 
# visreg(model_Bee,"Indirect_interactions",xlab="I. interactions",ylab="Percentil",
#        main="Bee",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
# points(Percentil ~ Indirect_interactions, data = Bee_freq_position_s_ind, 
#        col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
#        pch = 20) 
# 
# visreg(model_coleoptera,"Indirect_interactions",xlab="I. interactions",ylab="Percentil",
#        main="Coleoptera",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
# points(Percentil ~ Indirect_interactions, data = coleoptera_freq_position_s_ind, 
#        col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
#        pch = 20) 
# 
# visreg(model_Lepidoptera,"Specificity",xlab="Specificity",ylab="Percentil",
#        main="Lepidoptera",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
# points(Percentil ~ Indirect_interactions, data = Lepidoptera_freq_position_s_ind, 
#        col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
#        pch = 20) 
# 
# visreg(model_NonbeeHymenoptera,"Indirect_interactions",xlab="I. interactions",ylab="Percentil",
#        main="Non-bee-Hymenoptera",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
# points(Percentil ~ Indirect_interactions, data = NonbeeHymenoptera_freq_position_s_ind, 
#        col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
#        pch = 20) 
# 
# visreg(model_Nonsyrphidsdiptera,"Indirect_interactions",xlab="I. interactions",ylab="Percentil",
#        main="Non-syrphids-diptera",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
# points(Percentil ~ Indirect_interactions, data = Nonsyrphidsdiptera_freq_position_s_ind, 
#        col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
#        pch = 20) 
# 
# visreg(model_Syrphids,"Indirect_interactions",xlab="I. interactions",ylab="Percentil",
#        main="Syrphids",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
# points(Percentil ~ Indirect_interactions, data = Syrphids_freq_position_s_ind, 
#        col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
#        pch = 20) 
# 
# # save
# #
# #
# #
# #
# #
# dev.off()



#---------------------------------
# Percentil analises with LMMs for Plants


library(lme4)
library(performance)

plant_freq_position_s <- plant_position_percentiles_filtered %>% 
  left_join(specificity_plant, by = "position")

plant_freq_position_s$s[is.nan(plant_freq_position_s$s)] <- 1.0

# Indirect interaction data
motif_interactions <- read_csv("Data/Data_processing/Motifs_connections/motif_interactions.csv") %>%
  select(position,indirect_interactions) %>% unique()

plant_freq_position_s_ind <- plant_freq_position_s %>% 
  left_join(motif_interactions, by = "position") %>% 
  filter(!Node_FG %in% c("Birds","Lizards","Other_insects")) %>%
  rename(Percentil=percentil_its_GF,Specificity=s,Indirect_interactions=indirect_interactions,Group=Node_FG)

# plant_freq_position_s_ind$Percentil[plant_freq_position_s_ind$Percentil==0] <- 
#   1e-10
# 
# plant_freq_position_s_ind$Percentil[plant_freq_position_s_ind$Percentil==1] <- 
#   1-1e-10

plant_freq_position_s_ind$position <- as.factor(plant_freq_position_s_ind$position)

model_freq_specif_ind_plant <- glmmTMB(Percentil ~ scale(Specificity)*Group +
                                         scale(Indirect_interactions)*Group + 
                                   (1|study_id/position),
                                 plant_freq_position_s_ind)


summary(model_freq_specif_ind_plant)
r2(model_freq_specif_ind_plant)

library(visreg)
dev.off()
visreg2d(model_freq_specif_ind_plant, "Group","Specificity",scale ="response")
visreg2d(model_freq_specif_ind_plant, "Group","Indirect_interactions",scale ="response")

# Plant Model REAL OUTPUTS!!!

plants_sp_plot <- visreg(model_freq_specif_ind_plant, "Specificity", by="Group",
       strip.names=c("Selfing herbs",
                     "Small outcrossing perennials",
                     "Self-incompatible perennials\nwith large flowers",
                     "Tall plants with small\nunisexual flowers",
                     "Short-lived outcrossers with\nlong zygomorphic flowers"),gg = TRUE)+
  theme_bw()

plants_ind_int_plot <- visreg(model_freq_specif_ind_plant, "Indirect_interactions", by="Group",
                         strip.names=c("Selfing herbs",
                                       "Small outcrossing perennials",
                                       "Self-incompatible perennials\nwith large flowers",
                                       "Tall plants with small\nunisexual flowers",
                                       "Short-lived outcrossers with\nlong zygomorphic flowers"),gg = TRUE)+
  theme_bw()+xlab("Indirect interactions")

model_freq_specif_ind_plant$fit

library(patchwork)
plants_sp_plot/plants_ind_int_plot

#Save model data
saveRDS(model_freq_specif_ind_plant, "Data/RData/model_freq_specif_ind_plant.RData")
saveRDS(plant_freq_position_s_ind, "Data/RData/plant_freq_position_s_ind.RData")


library("DHARMa")
plant_freq_position_s_ind$scale_Indirect_interactions <- scale(plant_freq_position_s_ind$Indirect_interactions)
model_freq_specif_ind_plant2 <- lmer(Percentil ~ Specificity*Group +
                                 scale_Indirect_interactions*Group + 
                                 (1|study_id/position),
                               plant_freq_position_s_ind)
simulationOutput_plants <- simulateResiduals(fittedModel = model_freq_specif_ind_plant2, n = 1050, use.u=T)
plot(simulationOutput_plants)
testDispersion(simulationOutput_plants)
plotResiduals(simulationOutput_plants, form = plant_freq_position_s_ind$Group)
plotResiduals(simulationOutput_plants, form = plant_freq_position_s_ind$Specificity)
plotResiduals(simulationOutput_plants, form = plant_freq_position_s_ind$Indirect_interactions)
plotResiduals(simulationOutput_plants, form = plant_freq_position_s_ind$position)

hist(plant_freq_position_s_ind$Percentil, xlab = "Response", main = "")
########################################3
#


# ############################
# #Visualization of slopes by using visreg
# 
# library(ggeffects)
# library(scales)
# G1_freq_position_s_ind <- plant_freq_position_s_ind %>% 
#   filter(Group =="1")
# 
# model_G1 <- glmmTMB(Percentil ~ Specificity + 
#                        Indirect_interactions + 
#                        (1|study_id), G1_freq_position_s_ind,
#                      family = beta_family())
# 
# G2_freq_position_s_ind <- plant_freq_position_s_ind %>% 
#   filter(Group =="2")
# 
# model_G2 <- glmmTMB(Percentil ~ Specificity + 
#                               Indirect_interactions + 
#                               (1|study_id), G2_freq_position_s_ind,
#                             family = beta_family())
# 
# G3_freq_position_s_ind <- plant_freq_position_s_ind %>% 
#   filter(Group =="3")
# model_G3<- glmmTMB(Percentil ~ Specificity + 
#                                Indirect_interactions + 
#                                (1|study_id), G3_freq_position_s_ind,
#                              family = beta_family())
# 
# G4_freq_position_s_ind <- plant_freq_position_s_ind %>% 
#   filter(Group =="4")
# model_G4 <- glmmTMB(Percentil ~ Specificity + 
#                                      Indirect_interactions + 
#                                      (1|study_id), G4_freq_position_s_ind,
#                                    family = beta_family())
# G5_freq_position_s_ind <- plant_freq_position_s_ind %>% 
#   filter(Group =="5")
# model_G5 <- glmmTMB(Percentil ~ Specificity + 
#                                       Indirect_interactions + 
#                                       (1|study_id), G5_freq_position_s_ind,
#                                     family = beta_family())
# 
# 
# 
# 
# #####################################
# dev.off()
# par(mfrow = c(5,2),mar=c(3.95,4,1,1))
# 
# visreg(model_G1,"Specificity",xlab="Specificity",ylab="Percentil",
#        main="Selfing herbs",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
# points(Percentil ~ Specificity, data = G1_freq_position_s_ind, 
#        col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
#        pch = 20) 
# 
# visreg(model_G1,"Indirect_interactions",xlab="I. interactions",ylab="Percentil",
#        main="Selfing herbs",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
# points(Percentil ~ Indirect_interactions, data = G1_freq_position_s_ind, 
#        col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
#        pch = 20) 
# 
# 
# visreg(model_G2,"Specificity",xlab="Specificity",ylab="Percentil",
#        main="Small outcrossing perennials",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
# points(Percentil ~ Specificity, data = G2_freq_position_s_ind, 
#        col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
#        pch = 20) 
# 
# visreg(model_G2,"Indirect_interactions",xlab="I. interactions",ylab="Percentil",
#        main="Small outcrossing perennials",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
# points(Percentil ~ Indirect_interactions, data = G2_freq_position_s_ind, 
#        col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
#        pch = 20) 
# 
# 
# 
# visreg(model_G3,"Specificity",xlab="Specificity",ylab="Percentil",
#        main="Self-incompatible perennials with large flowers",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
# points(Percentil ~ Specificity, data = G3_freq_position_s_ind, 
#        col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
#        pch = 20) 
# 
# visreg(model_G3,"Indirect_interactions",xlab="I. interactions",ylab="Percentil",
#        main="Self-incompatible perennials with large flowers",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
# points(Percentil ~ Indirect_interactions, data = G3_freq_position_s_ind, 
#        col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
#        pch = 20) 
# 
# 
# 
# visreg(model_G4,"Specificity",xlab="Specificity",ylab="Percentil",
#        main="Tall plants with small unisexual flowers",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
# points(Percentil ~ Specificity, data = G4_freq_position_s_ind, 
#        col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
#        pch = 20) 
# 
# visreg(model_G4,"Indirect_interactions",xlab="I. interactions",ylab="Percentil",
#        main="Tall plants with small unisexual flowers",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
# points(Percentil ~ Indirect_interactions, data = G4_freq_position_s_ind, 
#        col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
#        pch = 20) 
# 
# 
# 
# visreg(model_G5,"Specificity",xlab="Specificity",ylab="Percentil",
#        main="Short-lived outcrossers with long zygomorphic flowers",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
# points(Percentil ~ Specificity, data = G5_freq_position_s_ind, 
#        col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
#        pch = 20) 
# 
# visreg(model_G5,"Indirect_interactions",xlab="I. interactions",ylab="Percentil",
#        main="Short-lived outcrossers with long zygomorphic flowers",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
# points(Percentil ~ Indirect_interactions, data = G5_freq_position_s_ind, 
#        col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
#        pch = 20) 
# 
# 
# 
# # save
# #
# #
# #
# #
# #
# dev.off()
