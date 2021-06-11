# Load libraries
library(tidyverse)
library(stringi)
library(glmmTMB)
source("Scripts/Scripts_Alfonso/add_study_id.R")

# List the files with information about the positions
folder_motif_data <- "Data/Csv/Motifs_positions_and_null_models"

motif_files <- list.files(folder_motif_data) 
motif_link_files <-  motif_files[grep("Binary_Mat_Sp_Motifs_positions_null_Conn_FiX_",motif_files)]
pollinator_positions_codes <- c("np2","np4","np6","np8","np11","np12","np14","np16",
                                "np18","np21","np22","np24","np25","np28","np29","np31",
                                "np34","np35","np38","np41","np42","np44","np46")
plant_positions_codes <- c("np1","np3","np5","np7","np9","np10","np13","np15",
                           "np17","np19","np20","np23","np26","np27","np30","np32",
                           "np33","np36","np37","np39","np40","np43","np45")

pollinator_absolut_frequencies <- NULL
plant_absolut_frequencies <- NULL

for(i in 1:length(motif_link_files)){

  print(motif_link_files[i])
  
  # Open link file and rearrange it
  file_i <- paste0(folder_motif_data, "/", motif_link_files[i])
  positions_i <- read_csv(file_i)
  
  positions_i_observed <- positions_i %>% filter(Type_network=="observed") %>%
    select(-Node_id)
  
  positions_i_observed_agg <- positions_i_observed %>% 
    group_by(Network_id,Type_network,Node_FG) %>% summarise_all(sum)
  
  
  ##############
  # Break down the data for plants and pollinator
  
  poll_position_col <- which(colnames(positions_i_observed_agg) %in% pollinator_positions_codes)
  plant_position_col <- which(colnames(positions_i_observed_agg) %in% plant_positions_codes)
  
  pollinator_frequencies_raw_i <- positions_i_observed_agg %>% filter(!Node_FG %in% c(1:5))
  
  pollinator_frequencies_agg_i <-  pollinator_frequencies_raw_i[,-plant_position_col]
  
  
  plant_frequencies_raw_i <- positions_i_observed_agg %>% filter(Node_FG %in% c(1:5))
  
  plant_frequencies_agg_i <-  plant_frequencies_raw_i[,-poll_position_col]

  
  pollinator_absolut_frequencies <- bind_rows(pollinator_absolut_frequencies,
                                               pollinator_frequencies_agg_i)
  
  plant_absolut_frequencies <- bind_rows(plant_absolut_frequencies,
                                          plant_frequencies_agg_i)
  
}

pollinator_absolut_frequencies <- pollinator_absolut_frequencies %>% 
  ungroup %>% select(-Type_network) %>%
  gather(position,frequency,-Network_id,-Node_FG)

pollinator_absolut_frequencies$position <- stri_extract_first_regex(
  pollinator_absolut_frequencies$position,"[0-9]+") %>% as.numeric()

pollinator_absolut_frequencies <- add_study_id(pollinator_absolut_frequencies)

plant_absolut_frequencies <- plant_absolut_frequencies %>% 
  ungroup %>% select(-Type_network) %>%
  gather(position,frequency,-Network_id,-Node_FG)

plant_absolut_frequencies$position <- stri_extract_first_regex(
  plant_absolut_frequencies$position,"[0-9]+") %>% as.numeric()

plant_absolut_frequencies <- add_study_id(plant_absolut_frequencies)

############################################
# EXTRACT MEAN FREQUENCIES PER FG, POSITION
############################################

# Pollinators--------

# Mean percentile and SE taking into account the study system
pollinator_codes <- pollinator_absolut_frequencies$position %>% unique()
pollinator_means <- pollinator_absolut_frequencies %>% select(position,Node_FG) %>% unique()
pollinator_means$mean <- NA
pollinator_means$SE <- NA
pollinator_means_reordered <- NULL

for(i.pos in 1:length(pollinator_codes)){
  
  pollinator_absolut_frequencies_aux <- pollinator_absolut_frequencies %>%
    filter(position == pollinator_codes[i.pos])
  pollinator_absolut_frequencies_aux$study_id <- as.factor(pollinator_absolut_frequencies_aux$study_id)
  pollinator_absolut_frequencies_aux$Node_FG <- as.factor(pollinator_absolut_frequencies_aux$Node_FG)
  pollinator_absolut_frequencies_aux$Network_id <- as.factor(pollinator_absolut_frequencies_aux$Network_id)
  
  m <- glmmTMB(frequency ~ Node_FG + (1|study_id),
               family = nbinom2(),
               data = pollinator_absolut_frequencies_aux)
  
  #summary(m)
 
  m_fixef <- fixef(m)
  
  temp <- m_fixef$cond %>% unname()
  
  pollinator_means_aux <- pollinator_means[pollinator_means$position == pollinator_codes[i.pos],] %>% arrange(Node_FG)
  
  pollinator_means_aux$mean <- c(temp[1],
                                 temp[1]+temp[2],
                                 temp[1]+temp[3],
                                 temp[1]+temp[4],
                                 temp[1]+temp[5],
                                 temp[1]+temp[6],
                                 temp[1]+temp[7],
                                 temp[1]+temp[8],
                                 temp[1]+temp[9])
  
  
  temp <- sqrt(diag(vcov(m)$cond)) %>% unname()
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

pollinator_means_reordered$mean_natural_units <- exp(pollinator_means_reordered$mean)


# Pollinators--------

# Mean percentile and SE taking into account the study system
plant_codes <- plant_absolut_frequencies$position %>% unique()
plant_means <- plant_absolut_frequencies %>% select(position,Node_FG) %>% unique()
plant_means$mean <- NA
plant_means$SE <- NA
plant_means_reordered <- NULL

for(i.pos in 1:length(plant_codes)){
  
  plant_absolut_frequencies_aux <- plant_absolut_frequencies %>%
    filter(position == plant_codes[i.pos])
  plant_absolut_frequencies_aux$study_id <- as.factor(plant_absolut_frequencies_aux$study_id)
  plant_absolut_frequencies_aux$Node_FG <- as.factor(plant_absolut_frequencies_aux$Node_FG)
  plant_absolut_frequencies_aux$Network_id <- as.factor(plant_absolut_frequencies_aux$Network_id)
  
  m <- glmmTMB(frequency ~ Node_FG + (1|study_id),
               family = nbinom2(),
               data = plant_absolut_frequencies_aux)
  
  #summary(m)
  
  m_fixef <- fixef(m)
  
  temp <- m_fixef$cond %>% unname()
  
  plant_means_aux <- plant_means[plant_means$position == plant_codes[i.pos],] %>% arrange(Node_FG)
  
  plant_means_aux$mean <- c(temp[1],
                                 temp[1]+temp[2],
                                 temp[1]+temp[3],
                                 temp[1]+temp[4],
                                 temp[1]+temp[5])
  
  
  temp <- sqrt(diag(vcov(m)$cond)) %>% unname()
  plant_means_aux$SE <- c(temp[1],
                               temp[1]+temp[2],
                               temp[1]+temp[3],
                               temp[1]+temp[4],
                               temp[1]+temp[5])
  
  plant_means_reordered <- bind_rows(plant_means_reordered,plant_means_aux)
  
}

plant_means_reordered$mean_natural_units <- exp(plant_means_reordered$mean)

############################
# PIE CHARTS
############################

library(ggplot2)
library(scatterpie)
library(RColorBrewer)

# Pollinator--------
pollinator_pies <- pollinator_means_reordered %>% select(-mean, -SE) %>% 
  spread(Node_FG,mean_natural_units)

pollinator_pies$x <- 0
pollinator_pies$y <- 0

ggplot() + geom_scatterpie(aes(x=x, y=y, group=position,r=5), 
                           data = pollinator_pies,
                           cols=as.character(c("Bee","Birds","Coleoptera","Lepidoptera",
                                               "Lizards","Non-bee-Hymenoptera",
                                               "Non-syrphids-diptera","Other_insects",
                                               "Syrphids")))+
  coord_equal()+
  scale_fill_brewer(palette = "Paired")+
  facet_wrap(~position,nrow = 3)+
  theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title.x=element_blank(),axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),axis.title.y=element_blank(),
          axis.text.y=element_blank(),axis.ticks.y=element_blank(),
          legend.position="bottom")+
  labs(title = "Pollinator positions",fill=NULL)
    
# plant--------
plant_pies <- plant_means_reordered %>% select(-mean, -SE) %>% 
  spread(Node_FG,mean_natural_units)

plant_pies$x <- 0
plant_pies$y <- 0

ggplot() + geom_scatterpie(aes(x=x, y=y, group=position,r=5), 
                           data = plant_pies,
                           cols=as.character(c("1","2","3","4","5")))+
  coord_equal()+
  scale_fill_brewer(palette = "Dark2")+
  facet_wrap(~position,nrow = 3)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank(),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),axis.title.y=element_blank(),
        axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        legend.position="bottom")+
  labs(title = "Plant positions",fill=NULL)


#################################
# CREATE PIE-MOTIFS
#################################

equidistribution <- function(number_nodes){
  
  x_width <- 5.5/(number_nodes+1)
  
  x_result <- rep(-1.5,number_nodes)
  
  for(i in 1:number_nodes){
    
    x_result[i] <- x_result[i] + i*x_width
    
  }
  
  return(x_result)
  
}

plot_pie_motif_positions <- function(motif_number){
  motifs_raw_data <- read_csv("Data/Data_processing/Motifs_connections/motif_pattern_connections.csv")
  
  motif_data <- motifs_raw_data %>% filter(motif_id == motif_number)
  
  
  # Create dataframe with node positions
  
  labdown = motif_data$plant %>% unique()
  
  
  xdown = equidistribution(length(labdown))
  
  ydown = rep(0,length(labdown))
  
  labup = motif_data$pollinator %>% unique()
  
  xup = equidistribution(length(labup))
  
  yup = rep(1,length(labup))
  
  d_nodes=data.frame(x=c(xdown,xup),
                     y=c(ydown,yup),
                     label = c(labdown,labup))
  
  
  # Create a dataframe with links
  
  plant_positions <- data.frame(plant = labdown,
                                xdown = xdown,
                                ydown = ydown)
  
  pollinator_positions <- data.frame(pollinator = labup,
                                     xup = xup,
                                     yup = yup)
  
  link_df <- motif_data %>% left_join(plant_positions, by = "plant") %>% 
    left_join(pollinator_positions, by = "pollinator")
  
  # Before plotting the motif, we remove letters from the node position id
  
  d_nodes$label <- gsub("[^0-9.-]", "", d_nodes$label)
  plant_pies_aux <- plant_pies %>% select(-x,-y) %>% mutate(position = as.character(position))
  
  d_nodes_plant <- d_nodes %>% rename(position = label) %>% 
    left_join(plant_pies_aux, by = "position") %>% as_tibble() %>%
    filter(!is.na(`1`))
  
  pollinator_pies_aux <- pollinator_pies %>% select(-x,-y) %>% mutate(position = as.character(position))
  
  d_nodes_pollinator <- d_nodes %>% rename(position = label) %>% 
    left_join(pollinator_pies_aux, by = "position") %>% as_tibble() %>%
    filter(!is.na(Bee))
  
  # Define the number of colors you want
  nb.cols <- 5
  mycolors <- c(brewer.pal(5, "Dark2"),brewer.pal(9, "Paired"))
  
  ggplot() +
    scale_shape_identity() +
    geom_segment(data = link_df, aes(x=xdown, y=ydown, xend = xup, yend = yup), 
                 size = 1, color = "grey20")+
    geom_scatterpie(aes(x=x, y=y, group=position,r=0.3), 
                    data = d_nodes_plant,
                    cols=as.character(c("1","2","3","4","5")))+
    geom_scatterpie(aes(x=x, y=y, group=position,r=.3), 
                    data = d_nodes_pollinator,
                    cols=as.character(c("Bee","Birds","Coleoptera","Lepidoptera",
                                        "Lizards","Non-bee-Hymenoptera",
                                        "Non-syrphids-diptera","Other_insects",
                                        "Syrphids")))+
    coord_equal()+
    scale_fill_manual(values = mycolors)+
    xlim(-1.5,4)+ ylim(-0.5,1.5)+
    theme(panel.background = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+
    labs(x=NULL,y=NULL)+
    guides(color = FALSE)+
    ggtitle(paste0("Motif ",motif_number))+
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")
}

motif_number <- 1
p1 <- plot_pie_motif_positions(motif_number)

motif_number <- 2
p2 <- plot_pie_motif_positions(motif_number)

motif_number <- 3
p3 <- plot_pie_motif_positions(motif_number)

motif_number <- 4
p4 <- plot_pie_motif_positions(motif_number)

motif_number <- 5
p5 <- plot_pie_motif_positions(motif_number)

motif_number <- 6
p6 <- plot_pie_motif_positions(motif_number)

motif_number <- 7
p7 <- plot_pie_motif_positions(motif_number)

motif_number <- 8   
p8 <- plot_pie_motif_positions(motif_number)

motif_number <- 9
p9 <- plot_pie_motif_positions(motif_number)

motif_number <- 10
p10 <- plot_pie_motif_positions(motif_number)

motif_number <- 11
p11 <- plot_pie_motif_positions(motif_number)

motif_number <- 12
p12 <- plot_pie_motif_positions(motif_number)

motif_number <- 13
p13 <- plot_pie_motif_positions(motif_number)

motif_number <- 14
p14 <- plot_pie_motif_positions(motif_number)

motif_number <- 15
p15 <- plot_pie_motif_positions(motif_number)

motif_number <- 16
p16 <- plot_pie_motif_positions(motif_number)

motif_number <- 17
p17 <- plot_pie_motif_positions(motif_number)


library(patchwork)

(p1 & ylab("2 species")| p2 & ylab("3 species") | p3 )/ (p4 & ylab("4 species")| p5 | p6 | p7)/ (p8  & ylab("5 species")| p9 | p10 | p11 | p12) / (p13 & ylab("5 species")| p14 | p15 | p16 | p17) & theme(axis.title.y = element_text(color="black", size=14, face="bold")) #+
#plot_annotation(tag_levels = '1') 

