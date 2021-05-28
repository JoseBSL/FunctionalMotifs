library(lme4)
library(nlme)
library(tidyverse)
library(stringi)
library(ggpubr)
source("Scripts/Scripts_Alfonso/add_study_id.R")

all_position_percentiles <- read_csv("Data/Csv/Motifs_positions_and_null_models/GF_positions_frequency_percentile.csv")

all_position_percentiles$position <- stri_extract_first_regex(all_position_percentiles$position,
                                                              "[0-9]+") %>% as.numeric()

all_position_percentiles <- add_study_id(all_position_percentiles)
str(all_position_percentiles)
head(all_position_percentiles)

list_plants_FG <- as.character(1:10)

plant_position_percentiles_filtered <- all_position_percentiles %>% 
  filter(Node_FG %in% list_plants_FG, !position %in% c(1,2),
         !is.na(expected_freq_its_GF)) #ME DA ERROR
pollinator_position_percentiles_filtered <- all_position_percentiles %>% 
  filter(!Node_FG %in% list_plants_FG, !position %in% c(1,2),
         !is.na(expected_freq_its_GF))

ggplot(plant_position_percentiles_filtered, aes(percentil_its_GF))+
  geom_histogram(color="black", fill="white")+
  facet_wrap(~position)+
  theme_bw()+
  labs(x="Sizeclass percentile", title = "Histograms: Percentiles for plant positions")

ggplot(pollinator_position_percentiles_filtered, aes(percentil_its_GF))+
  geom_histogram(color="black", fill="white")+
  facet_wrap(~position)+
  theme_bw()+
  labs(x="Sizeclass percentile",title = "Histograms: Percentiles for pollinator positions")



# estimated mean taking into account the random str is 0.XX +- 0.XX (look at the (Intercept) estimate)


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

GF_plant_pos <- ggplot(NULL) + 
  geom_point(data = plant_position_percentiles_filtered,
             aes(x=as.factor(position), y=percentil_its_GF, 
                 color = as.factor(position)),
             position = "jitter",alpha=0.5)+
  geom_errorbar(data = plant_means,aes(x = as.factor(position), ymin=mean-SE, ymax=mean+SE), 
                width=1.0,size=1)+
  facet_wrap(~Node_FG,ncol = 5)+
  geom_point(data = plant_means,aes(x = as.factor(position), y=mean), 
             size=2)+
  labs(x="Position", y = "Percentile",title = "Plants")+
  theme_bw()+
  theme(legend.position = "none")+#,axis.text.x=element_text(angle=90,vjust=0.5, hjust=1))+
  coord_flip()


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

GF_poll_pos <- ggplot(NULL) + 
  geom_point(data = pollinator_position_percentiles_filtered,
             aes(x=as.factor(position), y=percentil_its_GF, 
                 color = as.factor(position)),
             position = "jitter",alpha=0.5)+
  geom_errorbar(data = pollinator_means_reordered, aes(x = as.factor(position), ymin=lower,
                                                       ymax=mean+SE), 
                width=1.0,size=1)+
  facet_wrap(~Node_FG,ncol = 5)+
  geom_point(data = pollinator_means_reordered,aes(x = as.factor(position), y=mean), 
             size=2)+
  labs(x="Position", y = "Percentile", title = "Pollinators")+
  theme_bw()+
  theme(legend.position = "none")+#,axis.text.x=element_text(angle=90,vjust=0.5, hjust=1))+
  coord_flip()


# Plant and Pollinators--------

ggarrange(GF_plant_pos,GF_poll_pos,
          ncol = 1, nrow = 2,heights = c(.90, 1.4))



str(pollinator_position_percentiles_filtered)


library(tidyHeatmap)

tidy_heatmap(data_exprs,
             rows = external_gene_name,
             columns = sample,
             values = expression,
             scale = "row",
             annotation_col = c(sample_type, condition, group),
             annotation_row = c(is_immune_gene, direction),
             gaps_row = direction,
             gaps_col = group
)

