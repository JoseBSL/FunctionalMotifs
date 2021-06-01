
plot_motif_positions <- function(motif_number){
  
  motifs_raw_data <- read_csv("Data/Data_processing/Motifs_connections/motif_pattern_connections.csv")
  
  motif_data <- motifs_raw_data %>% filter(motif_id == motif_number)
  
  
  # Create dataframe with node positions
  
  labdown = motif_data$plant %>% unique()
  
  if(length(labdown)==2){
    
    xdown = c(0.75,2.25)
    
  }else if(length(labdown)==3){
    
    xdown = c(0,1.5,3)
    
  }else{
    
    xdown = c(1.5)
  }
  
  ydown = rep(0,length(labdown))
  
  labup = motif_data$pollinator %>% unique()
  
  if(length(labup)==2){
    
    xup = c(0.75,2.25)
    
  }else if(length(labup)==3){
    
    xup = c(0,1.5,3)
    
  }else{
    
    xup = c(1.5)
  }
  
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
  
  ggplot(data=d_nodes, aes(x, y)) +
    scale_shape_identity() +
    geom_segment(data = link_df, aes(x=xdown, y=ydown, xend = xup, yend = yup), size = 2, color = "grey")+
    geom_point(aes(color = as.factor(label)),size=20)+
    geom_text(aes(x, y , label = label),size=6)+
    scale_color_brewer(palette = "Set2")+
    xlim(-1,4)+ ylim(-0.5,1.5)+
    theme(panel.background = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+
    labs(x=NULL,y=NULL)+
    guides(color = FALSE)
  
}