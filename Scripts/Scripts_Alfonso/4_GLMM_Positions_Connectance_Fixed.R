
library(tidyverse)
library(lme4)
library(nlme)

# List the files with information about the positions
folder_motif_data <- "Data/Csv/Motifs positions and null models"

position_files <- list.files(folder_motif_data) 
position_percentile_files <-  position_files[grep("Position_quantiles_Conn_Fixed_",
                                                  position_files)]

all_position_percentiles <- NULL

for(i in 1:length(position_percentile_files)){
  
  # Open link file and rearrange it
  file_i <- paste0(folder_motif_data, "/", position_percentile_files[i])
  all_position_percentiles_i <- read_csv(file_i)
  
  all_position_percentiles <- bind_rows(all_position_percentiles, all_position_percentiles_i)
  
}

list_plants_FG <- as.character(1:10)

plant_position_percentiles <- all_position_percentiles %>% 
  filter(Node_FG %in% list_plants_FG)
pollinatior_position_percentiles <- all_position_percentiles %>% 
  filter(!Node_FG %in% list_plants_FG)

# Remove plant positions from  the data for pollinators' FG

pollinator_positions <- which(is.na(plant_position_percentiles[1,]))
pollinator_position_percentiles_filtered <- 
  pollinatior_position_percentiles[,c(1,3,pollinator_positions)] %>%
  gather("position", "percentile",-Network_id,-Node_FG)

# Remove pollinator positions from  the data for plant' FG
plant_positions <- which(is.na(pollinatior_position_percentiles[1,]))
plant_position_percentiles_filtered <- plant_position_percentiles[,c(1,3,plant_positions)] %>%
  gather("position", "percentile",-Network_id,-Node_FG)



  
# Plot positions

ggplot(plant_position_percentiles_filtered, 
       aes(x = position, y = percentile, color = Node_FG))+
  geom_boxplot(alpha=0.3)+
  facet_wrap(~Node_FG)

ggplot(pollinator_position_percentiles_filtered, 
       aes(x = position, y = percentile, color = Node_FG))+
  geom_boxplot(alpha=0.3)+
  facet_wrap(~Node_FG)

#############################################
# EXTRACT 95% CI for Pollinators percentiles
#############################################

pollinator_positions_codes <- pollinator_position_percentiles_filtered$position %>% unique()
pollinator_positions_CI <- tibble(position = pollinator_positions_codes)
pollinator_positions_CI$mean <- NA
pollinator_positions_CI$lower_CI <- NA
pollinator_positions_CI$upper_CI <- NA

for(i.pos in 1:length(pollinator_positions_codes)){
  pollinator_m <- lmer(percentile ~ 1+(1|Node_FG),
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
results_pollinator_position$evaluation[results_pollinator_position$percentile > 
                                         results_pollinator_position$upper_CI ] <- "(97.5,100]"
results_pollinator_position$evaluation[results_pollinator_position$percentile < 
                                         results_pollinator_position$lower_CI ] <- "[0,2.5)"

ggplot(results_pollinator_position, 
       aes(x = position, y = percentile, color = evaluation))+
  geom_point(alpha=0.3)+
  facet_wrap(~Node_FG)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1))+#,legend.position="bottom")+
  labs(x=NULL, color = "Percentile's C.I.")

library(RColorBrewer)
ggplot(results_pollinator_position %>% 
         mutate(evaluation = factor(evaluation, 
                                    levels=c("(97.5,100]", "[2.5,97.5]", "[0,2.5)")),
                position = factor(position,
                                  levels = 
                                    colnames(pollinatior_position_percentiles[,pollinator_positions]))), 
       aes(x = position, y = percentile, fill = evaluation))+
  geom_bar(position="stack", stat="identity")+
  scale_fill_brewer(palette = "RdYlBu") +
  facet_wrap(~Node_FG)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1))+#,legend.position="bottom")+
  labs(x=NULL, y = "Number of observations", fill = "Percentile's C.I.")

ggplot(results_pollinator_position %>% 
         mutate(evaluation = factor(evaluation, 
                                    levels=c("(97.5,100]", "[2.5,97.5]", "[0,2.5)")),
                position = factor(position,
                                  levels = 
                                    colnames(pollinatior_position_percentiles[,pollinator_positions]))), 
       aes(x = position, y = percentile, fill = evaluation))+
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
  plant_m <- lmer(percentile ~ 1+(1|Node_FG),
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
results_plant_position$evaluation[results_plant_position$percentile > 
                                         results_plant_position$upper_CI ] <- "(97.5,100]"
results_plant_position$evaluation[results_plant_position$percentile < 
                                         results_plant_position$lower_CI ] <- "[0,2.5)"

ggplot(results_plant_position, 
       aes(x = position, y = percentile, color = evaluation))+
  geom_point(alpha=0.3)+
  facet_wrap(~Node_FG)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1))+#,legend.position="bottom")+
  labs(x=NULL, color = "Percentile's C.I.")

library(RColorBrewer)
ggplot(results_plant_position %>% 
         mutate(evaluation = factor(evaluation, 
                                    levels=c("(97.5,100]", "[2.5,97.5]", "[0,2.5)")),
                position = factor(position,
                                  levels = 
                                    colnames(pollinatior_position_percentiles[,plant_positions]))), 
       aes(x = position, y = percentile, fill = evaluation))+
  geom_bar(position="stack", stat="identity")+
  scale_fill_brewer(palette = "RdYlBu") +
  facet_wrap(~Node_FG)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1))+#,legend.position="bottom")+
  labs(x=NULL, y = "Number of observations", fill = "Percentile's C.I.")

ggplot(results_plant_position %>% 
         mutate(evaluation = factor(evaluation, 
                                    levels=c("(97.5,100]", "[2.5,97.5]", "[0,2.5)")),
                position = factor(position,
                                  levels = 
                                    colnames(pollinatior_position_percentiles[,plant_positions]))), 
       aes(x = position, y = percentile, fill = evaluation))+
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
