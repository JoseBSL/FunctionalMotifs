########################################################################################################################################################
#SCRIPT TO MERGE LONG DATA WITH Z-SCORES, TRAIT DATA AND FUNCTIONAL GROUPS

#1) LOAD DATA 

#3) MERGE WITH TRAIT DATA

#4) MERGE WITH FUNCTIONAL GROUPS

#5) SAVE DATA
########################################################################################################################################################

#LOAD LIBRARIES
library(data.table)
library(readxl) #read trait data
library(dplyr)
########################################################################################################################################################
#1) LOAD LONG FORMAT DATA WITH Z-SCORES, TRAIT DATA AND FUNCTIONAL GROUP DATA
########################################################################################################################################################
long_networks <- read.csv("Data/Csv/long_format_quantitative_networks.csv")
head(long_networks)
t_data <- read_excel("Data/Trait_data_raw/Trait_data_final.xlsx")
#head(t_data)
#Read cluster data
hclust_d_5 <- read.csv("Data/Csv/imputed_trait_data_hclust_5_clusters_forest_data.csv") #5 clusters

########################################################################################################################################################
#3) MERGE WITH TRAIT DATA
########################################################################################################################################################
#process trait data all species
t_data <- t_data[1:1712,] #reading to max length of the excel sheet with data
#filter data, select species with flower level info and capitulum
t_data_filtered <- filter(t_data, Info_level == "flower" |  Info_level == "capitulum")
t_data_filtered <- as.data.frame(t_data_filtered)
colnames(t_data_filtered)[1] <- "Plant_species"
#select columns of interest
traits <- t_data_filtered %>% select(Species_geonet,Order_all,Family_all,Genus_all,Species_all,Breeding_system,IMPUTED_Compatibility,Autonomous_selfing_level,Autonomous_selfing_level_data_type,Autonomous_selfing_level_fruit_set,Flower_morphology,Flower_symmetry,Flowers_per_plant,Flowers_per_inflorescence,Floral_unit_width,Corolla_diameter_mean,Corolla_length_mean,STYLE_IMPUTED,OVULES_IMPUTED,life_form,lifespan,IMPUTED_plant_height_mean_m,Nectar_presence_absence)
str(t_data_filtered)

#Remove duplicated species
t <- traits[!duplicated(traits$Species_all), ]
colnames(t)[1] <- "Plant_species"

#MERGE NETWORK AND TRAIT DATA
long_format_trait_data <- merge(long_networks, t, by = "Plant_species", all.x =T)
str(long_format_trait_data)
########################################################################################################################################################
#4) MERGE WITH FUNCTIONAL GROUPS
########################################################################################################################################################
#5 CLUSTERS
#Some data processing before merging
#colnames to merge by same id
colnames(hclust_d_5)[1] <-"Species_all"
#convert to chracters
hclust_d_5$Species_all <- as.character(hclust_d_5$Species_all)
long_format_trait_data$Species_all <- as.character(long_format_trait_data$Species_all)

#convert character NA's to NA's
long_format_trait_data$Species_all[long_format_trait_data$Species_all=="NA"]<- NA
hclust_d_5$Species_all[hclust_d_5$Species_all=="NA"]<- NA

#remove NA'S
long_format_trait_data <- long_format_trait_data[!is.na(long_format_trait_data$Species_all),]
hclust_d_5 <- hclust_d_5[!is.na(hclust_d_5$Species_all),]

#merge
final_d <- merge(long_format_trait_data,hclust_d_5, by="Species_all", all.x  = T) #now data is ready for analysis
head(final_d)


#Select columns of interest
final_d_1 <- final_d[,c("Species_all","Id","Pollinator_species","Interaction", "order","family","genus","guild","Order_all","Family_all","Genus_all","Clusters")]
head(final_d_1)

#change colnames to more informative ones
colnames(final_d_1) <- c("Plant_species", "Network_id", "Pollinator_species", "Interaction", "Pollinator_order", "Pollinator_family", "Pollinator_genus", 
                         "Pollinator_functional_group","Plant_order", "Plant_family", "Plant_genus", "Plant_functional_groups")


#Exclude rows without info to species level
final_d_1 <- final_d_1[ grep("Alternanthera sp", final_d_1$Plant_species, invert = TRUE) , ]
final_d_1 <- final_d_1[ grep("Arrabidaea sp", final_d_1$Plant_species, invert = TRUE) , ]
final_d_1 <- final_d_1[ grep("Carduus sp.", final_d_1$Plant_species, invert = TRUE) , ]
final_d_1 <- final_d_1[ grep("Croton sp. 1", final_d_1$Plant_species, invert = TRUE) , ]
final_d_1 <- final_d_1[ grep("Croton sp. 2", final_d_1$Plant_species, invert = TRUE) , ]
final_d_1 <- final_d_1[ grep("Croton sp. 2", final_d_1$Plant_species, invert = TRUE) , ]
final_d_1 <- final_d_1[ grep("Eupatorium sp", final_d_1$Plant_species, invert = TRUE) , ]
final_d_1 <- final_d_1[ grep("Gymnocalycium sp", final_d_1$Plant_species, invert = TRUE) , ]
final_d_1 <- final_d_1[ grep("Jatropha sp", final_d_1$Plant_species, invert = TRUE) , ]
final_d_1 <- final_d_1[ grep("Leontodon sp.", final_d_1$Plant_species, invert = TRUE) , ]
final_d_1 <- final_d_1[ grep("Linum sp.", final_d_1$Plant_species, invert = TRUE) , ]
final_d_1 <- final_d_1[ grep("Malpighiaceae sp 1", final_d_1$Plant_species, invert = TRUE) , ]
final_d_1 <- final_d_1[ grep("Myrtaceae sp 1", final_d_1$Plant_species, invert = TRUE) , ]
final_d_1 <- final_d_1[ grep("Oxalis sp", final_d_1$Plant_species, invert = TRUE) , ]
final_d_1 <- final_d_1[ grep("Pavonia sp", final_d_1$Plant_species, invert = TRUE) , ]
final_d_1 <- final_d_1[ grep("Pectis sp", final_d_1$Plant_species, invert = TRUE) , ]
final_d_1 <- final_d_1[ grep("Piriqueta sp", final_d_1$Plant_species, invert = TRUE) , ]
final_d_1 <- final_d_1[ grep("Ranunculus sp.", final_d_1$Plant_species, invert = TRUE) , ]
final_d_1 <- final_d_1[ grep("Retama sp.", final_d_1$Plant_species, invert = TRUE) , ]
final_d_1 <- final_d_1[ grep("Senegalia sp", final_d_1$Plant_species, invert = TRUE) , ]
final_d_1 <- final_d_1[ grep("Sida sp.", final_d_1$Plant_species, invert = TRUE) , ]
final_d_1 <- final_d_1[ grep("Tradescantia sp", final_d_1$Plant_species, invert = TRUE) , ]
final_d_1 <- final_d_1[ grep("Trifolium sp.", final_d_1$Plant_species, invert = TRUE) , ]
final_d_1 <- final_d_1[ grep("Vernonia sp", final_d_1$Plant_species, invert = TRUE) , ]
final_d_1 <- final_d_1[ grep("Xyris sp.", final_d_1$Plant_species, invert = TRUE) , ]



#remove .csv from Network id column
final_d_1$Network_id <- gsub(".csv", "", final_d_1$Network_id)

na <- final_d_1[is.na(final_d_1$Plant_functional_groups),]
levels(factor(na$Pollinator_functional_group))
#No NA's in funcional groups of pollinators and plants

########################################################################################################################################################
#5) SAVE DATA
########################################################################################################################################################
write.csv(final_d_1, "Data/Csv/data_for_motifs_analysis_1.csv")
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################


