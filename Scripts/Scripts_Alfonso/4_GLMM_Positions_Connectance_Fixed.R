
library(tidyverse)
library(lme4)
library(nlme)

# List the files with information about the positions
file_pos <- "Data/Csv/Motifs positions and null models/GF_positions_frequency_percentile.csv"
all_position_percentiles <- read_csv(file_pos)

# Create different dataframes for plant and pollinators, and remove data for postions np1 and np2,
# and leave only those positions that can be occupied either by plants or poll. 
list_plants_FG <- as.character(1:10)

plant_position_percentiles_filtered <- all_position_percentiles %>% 
  filter(Node_FG %in% list_plants_FG, !position %in% c("np1","np2"),
         !is.na(expected_freq_its_GF))
pollinator_position_percentiles_filtered <- all_position_percentiles %>% 
  filter(!Node_FG %in% list_plants_FG, !position %in% c("np1","np2"),
         !is.na(expected_freq_its_GF))

  
# Plot positions

ggplot(plant_position_percentiles_filtered, 
       aes(x = position, y = percentil_its_GF, color = Node_FG))+
  geom_boxplot(alpha=0.3)+
  facet_wrap(~Node_FG)

ggplot(plant_position_percentiles_filtered, 
       aes(x = position, y = percentil_all_GF, color = Node_FG))+
  geom_boxplot(alpha=0.3)+
  facet_wrap(~Node_FG)

ggplot(plant_position_percentiles_filtered, 
       aes(x = percentil_its_GF, y = percentil_all_GF, color = Node_FG))+
  geom_point(alpha=0.3)+
  facet_wrap(~position)

ggplot(pollinator_position_percentiles_filtered, 
       aes(x = position, y = percentil_its_GF, color = Node_FG))+
  geom_boxplot(alpha=0.3)+
  facet_wrap(~Node_FG)

ggplot(pollinator_position_percentiles_filtered, 
       aes(x = position, y = percentil_all_GF, color = Node_FG))+
  geom_boxplot(alpha=0.3)+
  facet_wrap(~Node_FG)

ggplot(pollinator_position_percentiles_filtered, 
       aes(x = percentil_its_GF, y = percentil_all_GF, color = Node_FG))+
  geom_point(alpha=0.3)+
  facet_wrap(~position)

#############################################
# EXTRACT 95% CI for Pollinators percentiles
#############################################

pollinator_positions_codes <- pollinator_position_percentiles_filtered$position %>% unique()
pollinator_positions_CI <- tibble(position = pollinator_positions_codes)
pollinator_positions_CI$mean <- NA
pollinator_positions_CI$lower_CI <- NA
pollinator_positions_CI$upper_CI <- NA

for(i.pos in 1:length(pollinator_positions_codes)){
  pollinator_m <- lmer(percentil_all_GF ~ 1+(1|Node_FG),
                       data = pollinator_position_percentiles_filtered %>% 
                         filter(position == pollinator_positions_codes[i.pos]))
  
  CI_m <- confint(pollinator_m)
  intercept_row <- which(row.names(CI_m)=="(Intercept)")
  pollinator_positions_CI$mean[i.pos] <- fixed.effects(pollinator_m) %>% unname()
  pollinator_positions_CI$lower_CI[i.pos] <- CI_m[intercept_row,1]
  pollinator_positions_CI$upper_CI[i.pos] <- CI_m[intercept_row,2]
}

results_pollinator_position <- pollinator_position_percentiles_filtered %>% 
  left_join(pollinator_positions_CI, by = "position")

results_pollinator_position$evaluation <- "[2.5,97.5]"
results_pollinator_position$evaluation[results_pollinator_position$percentil_all_GF > 
                                         results_pollinator_position$upper_CI ] <- "(97.5,100]"
results_pollinator_position$evaluation[results_pollinator_position$percentil_all_GF < 
                                         results_pollinator_position$lower_CI ] <- "[0,2.5)"

ggplot(results_pollinator_position, 
       aes(x = position, y = percentil_all_GF, color = evaluation))+
  geom_point(alpha=0.3)+
  facet_wrap(~Node_FG)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1))+#,legend.position="bottom")+
  labs(x=NULL, color = "Percentile's C.I.")


library(RColorBrewer)
ggplot(results_pollinator_position %>%
         mutate(evaluation = factor(evaluation, 
                                    levels=c("(97.5,100]", "[2.5,97.5]", "[0,2.5)")),
                position = factor(position),
                counts = 1), 
       aes(x = position, y = counts, fill = evaluation))+
  geom_bar(position="stack", stat="identity")+
  scale_fill_brewer(palette = "RdYlBu") +
  facet_wrap(~Node_FG)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1))+#,legend.position="bottom")+
  labs(x=NULL, y = "Number of observations", fill = "Percentile's C.I.")


ggplot(results_pollinator_position %>% 
         mutate(evaluation = factor(evaluation, 
                                    levels=c("(97.5,100]", "[2.5,97.5]", "[0,2.5)")),
                position = factor(position),
                counts = 1), 
       aes(x = position, y = counts, fill = evaluation))+
  geom_bar(position="fill", stat="identity")+
  scale_fill_brewer(palette = "RdYlBu") +
  facet_wrap(~Node_FG)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1))+#,legend.position="bottom")+
  labs(x=NULL, y = "Number of observations", fill = "Percentile's C.I.")


#############################################
# EXTRACT 95% CI for Plant percentiles
#############################################

plant_positions_codes <- plant_position_percentiles_filtered$position %>% unique()
plant_positions_CI <- tibble(position = plant_positions_codes)
plant_positions_CI$mean <- NA
plant_positions_CI$lower_CI <- NA
plant_positions_CI$upper_CI <- NA

for(i.pos in 1:length(plant_positions_codes)){
  plant_m <- lmer(percentil_all_GF ~ 1+(1|Node_FG),
                       data = plant_position_percentiles_filtered %>% 
                         filter(position == plant_positions_codes[i.pos]))
  
  CI_m <- confint(plant_m)
  intercept_row <- which(row.names(CI_m)=="(Intercept)")
  plant_positions_CI$mean[i.pos] <- fixed.effects(plant_m) %>% unname()
  plant_positions_CI$lower_CI[i.pos] <- CI_m[intercept_row,1]
  plant_positions_CI$upper_CI[i.pos] <- CI_m[intercept_row,2]
}

results_plant_position <- plant_position_percentiles_filtered %>% 
  left_join(plant_positions_CI, by = "position")

results_plant_position$evaluation <- "[2.5,97.5]"
results_plant_position$evaluation[results_plant_position$percentil_all_GF > 
                                         results_plant_position$upper_CI ] <- "(97.5,100]"
results_plant_position$evaluation[results_plant_position$percentil_all_GF < 
                                         results_plant_position$lower_CI ] <- "[0,2.5)"

ggplot(results_plant_position, 
       aes(x = position, y = percentil_all_GF, color = evaluation))+
  geom_point(alpha=0.3)+
  facet_wrap(~Node_FG)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1))+#,legend.position="bottom")+
  labs(x=NULL, color = "Percentile's C.I.")

library(RColorBrewer)
ggplot(results_plant_position %>% 
         mutate(evaluation = factor(evaluation, 
                                    levels=c("(97.5,100]", "[2.5,97.5]", "[0,2.5)")),
                position = factor(position),
                counts = 1), 
       aes(x = position, y = counts, fill = evaluation))+
  geom_bar(position="stack", stat="identity")+
  scale_fill_brewer(palette = "RdYlBu") +
  facet_wrap(~Node_FG)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1))+#,legend.position="bottom")+
  labs(x=NULL, y = "Number of observations", fill = "Percentile's C.I.")

ggplot(results_plant_position %>% 
         mutate(evaluation = factor(evaluation, 
                                    levels=c("(97.5,100]", "[2.5,97.5]", "[0,2.5)")),
                position = factor(position),
                counts = 1), 
       aes(x = position, y = counts, fill = evaluation))+
  geom_bar(position="fill", stat="identity")+
  scale_fill_brewer(palette = "RdYlBu") +
  facet_wrap(~Node_FG)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1))+#,legend.position="bottom")+
  labs(x=NULL, y = "Number of observations", fill = "Percentile's C.I.")



library(ggplot2)
library(scatterpie)

plant_scatter_data <- results_plant_position %>% 
  mutate(evaluation = factor(evaluation, 
                             levels=c("(97.5,100]", "[2.5,97.5]", "[0,2.5)")),
         position = factor(position,
                           levels = 
                             colnames(pollinatior_position_percentiles[,plant_positions])))

plant_scatter_data_total_counts <- results_plant_position %>% group_by(Node_FG,position) %>%
  count() %>% rename(total_counts = n)

positions_data = tibble(position = paste0("np",c(1:46)),num_position = c(1:46))
positions_data$degree <- 
  c(1,1,1,2,2,1,3,1,1,2,1,2,2,2,1,3,4,
    1,1,3,1,2,2,1,2,2,3,1,2,3,2,1,2,1,
    3,1,2,2,1,2,2,3,2,3,1,4)

positions_data$Motif_nodes <-
  c(2,2,3,3,3,3,rep(4,10),rep(5,30))

plant_scatter_data_FG <- results_plant_position %>% group_by(Node_FG,position,evaluation) %>%
  count() %>% rename(counts = n) %>% spread(evaluation,counts) %>%
  left_join(positions_data, by = "position") %>%
  left_join(plant_scatter_data_total_counts, by = c("Node_FG","position")) %>%
  rename(High=`(97.5,100]`,Low=`[0,2.5)`,Normal=`[2.5,97.5]`)


ggplot() + geom_scatterpie(aes(x=num_position, y=degree, group=Node_FG), data = plant_scatter_data_FG,
                           cols=as.character(c("High","Normal","Low")))+coord_equal()+
  facet_wrap(~Node_FG)+
  theme_bw()
