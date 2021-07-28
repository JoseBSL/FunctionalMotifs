library(tidyverse)
library(RColorBrewer)

equidistribution <- function(number_nodes){
  
  x_width <- 5.5/(number_nodes+1)
  
  x_result <- rep(-1.5,number_nodes)
  
  for(i in 1:number_nodes){
    
    x_result[i] <- x_result[i] + i*x_width
    
  }
  
  return(x_result)
  
}

plot_motif_label <- function(label){
  
  label_split <- str_split(label,"_")[[1]]
  motif_number <- as.numeric(label_split[1])
  nodes <- label_split[2:length(label_split)]
  plant_nodes <- nodes[nodes %in% plant_means_reordered$Node_FG]
  pollinator_nodes <- nodes[nodes %in% pollinator_means_reordered$Node_FG]
  
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
  
  
  # Add FG IDs to nodes
  FG_IDs <- tibble(label = c(sort(labdown),sort(labup)),
                   FG = c(plant_nodes,pollinator_nodes))
  
  d_nodes_FG <- d_nodes %>% left_join(FG_IDs, by = "label")
  levels_FG <- c("Bee","Birds","Coleoptera","Lepidoptera",
                 "Lizards","Non-bee-Hymenoptera", "Non-syrphids-diptera","Other_insects",
                 "Syrphids","1","2","3","4","5")
  
  d_nodes_FG$FG <- factor(d_nodes_FG$FG, levels = levels_FG)
  
  
  levels(d_nodes_FG$FG) <- c("Bee","Birds","Coleoptera","Lepidoptera",
                             "Lizards","Non-bee-Hymenoptera", "Non-syrphids-diptera","Other_insects",
                             "Syrphids","1","2","3","4","5")
  
  # Before plotting the motif, we remove letters from the node position id
  
  d_nodes_FG$label <- gsub("[^0-9.-]", "", d_nodes_FG$label)
  
  
  # Define the number of colors you want
  mycolors <- c(brewer.pal(9, "Paired"),brewer.pal(5, "Dark2"))
  
  ggplot(data=d_nodes_FG, aes(x, y)) +
    scale_shape_identity() +
    geom_segment(data = link_df, aes(x=xdown, y=ydown, xend = xup, yend = yup), size = 1, color = "grey20")+
    geom_point(aes(color = FG), size=10)+
    geom_text(aes(x, y , label = label), size=3)+
    scale_color_manual(drop = FALSE, values = mycolors,
                       breaks=c("Bee","Birds","Coleoptera","Lepidoptera",
                                "Lizards","Non-bee-Hymenoptera", "Non-syrphids-diptera","Other_insects",
                                "Syrphids","1","2","3","4","5"), 
                       labels=c("Bee","Birds","Coleoptera","Lepidoptera",
                                "Lizards","Non-bee-Hymenoptera", "Non-syrphids-diptera",
                                "Other insects", "Syrphids", 
                                "Selfing herbs","Small outcrossing perennials",
                                "Self-incomp. perennials\nwith large flowers",
                                "Tall plants with small\nunisexual flowers",
                                "Short-lived outcrossers with\nlong zygomorphic flowers"))+
    xlim(-1.5,4)+ ylim(-0.5,1.5)+
    theme(panel.background = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+
    labs(x=NULL,y=NULL)+
    ggtitle(paste0(label))+
    theme(plot.title = element_text(hjust = 0.5))+
    guides(color=guide_legend(title="Functional groups"))
  
  
}

# Data with pollinator/plant labels

plant_means_reordered <- read_csv("Data/Csv/plant_abs_freq_means.csv")
pollinator_means_reordered <- read_csv("Data/Csv/pollinator_abs_freq_means.csv")



# Examples

label <- "9_2_4_Coleoptera_Non-syrphids-diptera_Bee"
p1 <- plot_motif_label(label)
label <- "17_2_2_2_4_Bee"
p2 <- plot_motif_label(label)
label <- "3_5_Bee_Birds"
p3 <- plot_motif_label(label)

library(patchwork)
(p1|p2|p3) + plot_layout(guides = 'collect') &
  theme(legend.position = "bottom")


# Check top Z-scores of each group (Under, No diff, Over)

#Read data
data <- read.csv("Data/Csv/motifs_observed_probability_SIMUL_CI.csv")

#Find critical value of Z-score
p <- 0.05 #cutoff probability 95% confidence
critical_value <- qnorm(p/2) #double tail probability divide by 2

data_1 <- data %>%
  mutate(infra_over_represented = case_when(
    z_score < -abs(critical_value) ~ "infra",
    between(z_score, -abs(critical_value), abs(critical_value)) ~ "no_diff",
    z_score > abs(critical_value) ~ "over"
  ))


#data_2 <- filter(data_1, abs(z_score)<2000) 
#str(data_2)

#Select under-represented by the ones with highest porbability of being observed (top10)
under_data <- filter(data_1, infra_over_represented=="infra") 
str(under_data)

under_top_10 <- filter(under_data, row_number(desc(abs(motif_observed_probability)))<= 10) %>%
  select(motif_functional_ID,motif_observed_probability,counts_observed,z_score)

#Select no statistical difference by the ones with highest porbability of being observed (top10)
no_diff_data <- filter(data_1, infra_over_represented=="no_diff") 
str(no_diff_data)

no_diff_top_10 <- filter(no_diff_data, row_number(desc(abs(motif_observed_probability)))<= 10) %>%
  select(motif_functional_ID,motif_observed_probability,counts_observed,z_score)

#Select over-represenetd by the ones with highest porbability of being observed (top10)
over_data <- filter(data_1, infra_over_represented=="over") 
str(over_data)

over_top_10 <- filter(over_data, row_number(desc(abs(motif_observed_probability)))<= 20) %>%
  select(motif_functional_ID,motif_observed_probability,counts_observed,z_score)



#Trying a random plot// Aggregate of motifs and counts observed

sum_counts <- aggregate(counts_observed ~ motif, data=data, sum)

ggplot(data, aes(x = motif, y = counts_observed)) +
  geom_jitter()+ scale_x_continuous("motif", labels = as.character(motif), breaks = motif)
