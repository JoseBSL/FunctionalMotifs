########################################################################################################################################################
#SCRIPT FOR DATA PREPARATION->MERGE LONG FORMAT WITH SPP NAMES (JUST POLLINATORS) AND ADD GUILDS

#1) READ NETWORK DATA

#2) ADD POLLINATOR GUILD (Output of taxsize and manually -non found searches-)

#3) SAVE DATA
########################################################################################################################################################

#LOAD LIBRARIES
library(brms)
library(ape)
library(tidybayes)
library(ggplot2)
library(reshape2)
library(data.table)
library(dplyr)
library(bipartite)
library(readxl)
library(rtrees)
library(ape)
library(dplyr)
library(tidyverse)

########################################################################################################################################################
#1) LOAD NETWORK DATA
########################################################################################################################################################
#Set working directory to read files
setwd("~/R_Projects/Reproductive traits/Data/Data_networks_quantitative") 

temp <- list.files(pattern="*.csv")
my.list <- list(for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i])))
my_files <- list.files(pattern = "\\.csv$")
my_data <- lapply(my_files, function(i){read.csv(i, check.names=FALSE)})
#Add "id" to the list to the long format data
data_id_list <- lapply(seq_along(my_data), 
                       function(x) cbind(my_data[[x]], unique.id=my_files[x]))


#For loop to melt each data frame and merge
i <- NULL
all_long <- NULL

for (i in data_id_list){
  i <- melt(i)
  all_long <- rbind(all_long, i)
}

#Renaming columns
colnames(all_long) <- c("Plant_species", "Id", "Pollinator_species", "Interaction") 

########################################################################################################################################################
#2) ADD POLLINATOR GUILD
########################################################################################################################################################

#Load poll guild data (FILLED MANUALLY AND WITH TAXIZE PREVIOUSLY, I STILL HAVE TO ADD SOME FAM/GENUS/ORDER
#BECAUSE OF SOME NETWORKS WERE ADDED LATE)
#PLANT SPECIES NAMES HAS BEEN ALSO BEEN SEARCHED BUT PREVIOUSLY WITH TAXSIZE

setwd("~/R_Projects/Reproductive traits") 
poll_names <- read.csv("Data/Data_processing/pollinator_species_names/poll_spp_names_corrected.csv")
#select unique word in order to meger
all_long$genus_old <- word(all_long$Pollinator_species)
#merge data by genus with the long format dataframe of the networks
all_long_poll_names <- merge(all_long, poll_names, by= "genus_old", all.x  = T)

#Now I add manually some species that are unfilled
#convert columns to characters
all_long_poll_names$order <- as.character(all_long_poll_names$order)
all_long_poll_names$family <- as.character(all_long_poll_names$family)
all_long_poll_names$genus <- as.character(all_long_poll_names$genus)
#Acentron lentifera
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Acentron lentifera"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Acentron lentifera"] <- "Megachilidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Acentron lentifera"] <- "Megachile"
#Alepidosceles imitatrix
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Alepidosceles imitatrix"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Alepidosceles imitatrix"] <- "Apidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Alepidosceles imitatrix"] <- "Alepidosceles"
#Alloscirtetica brethesi
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Alloscirtetica brethesi"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Alloscirtetica brethesi"] <- "Apidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Alloscirtetica brethesi"] <- "Alloscirtetica"
#Anthanassa frisia
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Anthanassa frisia"] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Anthanassa frisia"] <- "Nymphalidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Anthanassa frisia"] <- "Anthanassa"
#Aphidoidea sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Aphidoidea sp.1"] <- "Hemiptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Aphidoidea sp.1"] <- NA
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Aphidoidea sp.1"] <- NA
#Apodemia
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Apodemia sp."] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Apodemia sp."] <- "Riodinidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Apodemia sp."] <- "Apodemia"
#Arhysosage melanica
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Arhysosage melanica"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Arhysosage melanica"] <- "Andrenidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Arhysosage melanica"] <- "Arhysosage"
#Aricoris sp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Aricoris sp."] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Aricoris sp."] <- "Riodinidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Aricoris sp."] <- "Aricoris"
#Ascia monuste
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Ascia monuste"] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Ascia monuste"] <- "Pieridae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Ascia monuste"] <- "Ascia"
#Atrichopogon yamabukiensis
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Atrichopogon yamabukiensis"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Atrichopogon yamabukiensis"] <- "Ceratopogonidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Atrichopogon yamabukiensis"] <- "Atrichopogon"
#Aricoris sp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Aricoris sp."] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Aricoris sp."] <- "Riodinidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Aricoris sp."] <- "Ariconius"
#Augochlorella aurata
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Augochlorella aurata"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Augochlorella aurata"] <- "Halictidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Augochlorella aurata"] <- "Augochlorella"
#Babiohaltica sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Babiohaltica sp.1"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Babiohaltica sp.1"] <- "Chrysomelidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Babiohaltica sp.1"] <- "Babiohaltica"
#Babiohaltica sp.2
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Babiohaltica sp.2"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Babiohaltica sp.2"] <- "Chrysomelidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Babiohaltica sp.2"] <- "Babiohaltica"
#Beetle sp.15
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Beetle sp.15"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Beetle sp.15"] <- NA
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Beetle sp.15"] <- NA
#Beetle sp.8
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Beetle sp.8"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Beetle sp.8"] <- NA
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Beetle sp.8"] <- NA
#Bombilius sp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Bombilius sp."] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Bombilius sp."] <- "Bombyliidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Bombilius sp."] <- "Bombylius"
#Bruchini sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Bruchini sp.1"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Bruchini sp.1"] <- "Chrysomelidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Bruchini sp.1"] <- NA
#Caenohalictus sp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Caenohalictus sp."] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Caenohalictus sp."] <- "Halictidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Caenohalictus sp."] <- "Caenohalictus"
#Camptodes vittata
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Camptodes vittata"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Camptodes vittata"] <- "Nitidulidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Camptodes vittata"] <- "Camptodes"
#Cardiophorus sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Cardiophorus sp.1"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Cardiophorus sp.1"] <- "Elateridae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Cardiophorus sp.1"] <- "Cardiophorus"
#Cecidomie sp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Cecidomie sp."] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Cecidomie sp."] <- "Cecidomyiidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Cecidomie sp."] <- NA
#Cecidomyiidae sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Cecidomyiidae sp.1"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Cecidomyiidae sp.1"] <- "Cecidomyiidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Cecidomyiidae sp.1"] <- NA
#Celama sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Celama sp.1"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Celama sp.1"] <- "Nolidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Celama sp.1"] <- "Celama"
#Cenophodes tamsi
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Cenophodes tamsi"] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Cenophodes tamsi"] <- "Sphingidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Cenophodes tamsi"] <- "Cenophodes"
#Centris sp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Centris sp."] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Centris sp."] <- "Apidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Centris sp."] <- "Centris"
#Centris tarsata
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Centris tarsata"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Centris tarsata"] <- "Apidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Centris tarsata"] <- NA
#Centris varia
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Centris varia"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Centris varia"] <- "Apidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Centris varia"] <- "Centris"
#Centris fuscata
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Centris fuscata"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Centris fuscata"] <- "Apidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Centris fuscata"] <- "Centris"
#Centris brethesi
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Centris brethesi"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Centris brethesi"] <- "Apidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Centris brethesi"] <- "Centris"
#Ceratalictus sp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Ceratalictus sp."] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Ceratalictus sp."] <- "Halictidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Ceratalictus sp."] <- "Ceratalictus"
#Chaetonerius alluaudii
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Chaetonerius alluaudii"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Chaetonerius alluaudii"] <- "Neriidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Chaetonerius alluaudii"] <- "Chaetonerius"
#Chalcidoidea sp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Chalcidoidea sp."] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Chalcidoidea sp."] <- NA
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Chalcidoidea sp."] <- NA
#Chalcidoidea sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Chalcidoidea sp.1"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Chalcidoidea sp.1"] <- NA
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Chalcidoidea sp.1"] <- NA
#Chelostoma philadelphi
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Chelostoma philadelphi"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Chelostoma philadelphi"] <- "Megachilidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Chelostoma philadelphi"] <- "Chelostoma"
#Chioides catillus
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Chioides catillus"] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Chioides catillus"] <- "Hesperiidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Chioides catillus"] <- "Chioides"
#Chloroprocta idioidea
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Chloroprocta idioidea"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Chloroprocta idioidea"] <- "Calliphoridae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Chloroprocta idioidea"] <- "Chloroprocta"
#Chlorostilbon lucidus
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Chlorostilbon lucidus"] <- "Apodiformes"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Chlorostilbon lucidus"] <- "Trochilidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Chlorostilbon lucidus"] <- "Chlorostilbon"
#Chryptocephalus sp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Chryptocephalus sp."] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Chryptocephalus sp."] <- "Chrysomelidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Chryptocephalus sp."] <- "Cryptocephalus"
#Chrysauginae sp1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Chrysauginae sp1"] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Chrysauginae sp1"] <- "Pyralidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Chrysauginae sp1"] <- NA
#Cixiidae sp
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Cixiidae sp"] <- "Hemiptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Cixiidae sp"] <- NA
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Cixiidae sp"] <- NA
#Coccinelidae sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Coccinelidae sp.1"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Coccinelidae sp.1"] <- "Coccinellidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Coccinelidae sp.1"] <- NA
#Coleoptera sp.2
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Coleoptera sp.2"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Coleoptera sp.2"] <- NA
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Coleoptera sp.2"] <- NA
#Coleoptera sp.7
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Coleoptera sp.7"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Coleoptera sp.7"] <- NA
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Coleoptera sp.7"] <- NA
#Coleoptera sp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Coleoptera sp."] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Coleoptera sp."] <- NA
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Coleoptera sp."] <- NA
#Coleoptera sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Coleoptera sp.1"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Coleoptera sp.1"] <- NA
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Coleoptera sp.1"] <- NA
#Coleoptera sp.4
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Coleoptera sp.4"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Coleoptera sp.4"] <- NA
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Coleoptera sp.4"] <- NA
#Coleoptera sp.3
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Coleoptera sp.3"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Coleoptera sp.3"] <- NA
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Coleoptera sp.3"] <- NA
#Coleoptera sp2
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Coleoptera sp2"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Coleoptera sp2"] <- NA
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Coleoptera sp2"] <- NA
#Coleoptera sp.5
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Coleoptera sp.5"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Coleoptera sp.5"] <- NA
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Coleoptera sp.5"] <- NA
#Compsomyiops fulvicrura
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Compsomyiops fulvicrura"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Compsomyiops fulvicrura"] <- "Calliphoridae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Compsomyiops fulvicrura"] <- "Compsomyiops"
#Copestilum aricia
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Copestilum aricia"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Copestilum aricia"] <- "Syrphidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Copestilum aricia"] <- "Copestylum"
#Copestilum sp.2
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Copestilum sp.2"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Copestilum sp.2"] <- "Syrphidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Copestilum sp.2"] <- "Copestylum"
#Copestilum sp.3
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Copestilum sp.3"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Copestilum sp.3"] <- "Syrphidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Copestilum sp.3"] <- "Copestylum"
#Coranarta cordigera
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Coranarta cordigera"] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Coranarta cordigera"] <- "Noctuidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Coranarta cordigera"] <- "Coronarta"
#Cosmosatyrus chilensis
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Cosmosatyrus chilensis"] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Cosmosatyrus chilensis"] <- "Nymphalidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Cosmosatyrus chilensis"] <- "Cosmosatyrus"
#Cossonini sp.2
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Cossonini sp.2"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Cossonini sp.2"] <- "Curculionidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Cossonini sp.2"] <- NA
#Cossonini sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Cossonini sp.1"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Cossonini sp.1"] <- "Curculionidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Cossonini sp.1"] <- NA
#Cossonini sp1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Cossonini sp1"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Cossonini sp1"] <- "Curculionidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Cossonini sp1"] <- NA
#Cryptophagidae sp.2
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Cryptophagidae sp.2"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Cryptophagidae sp.2"] <- "Cryptophagidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Cryptophagidae sp.2"] <- NA
#Cryptophagidae sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Cryptophagidae sp.1"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Cryptophagidae sp.1"] <- "Cryptophagidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Cryptophagidae sp.1"] <- NA
#Cryptophagidae sp2
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Cryptophagidae sp2"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Cryptophagidae sp2"] <- "Cryptophagidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Cryptophagidae sp2"] <- NA
#Cryptophleps nigrihalteratus
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Cryptophleps nigrihalteratus"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Cryptophleps nigrihalteratus"] <- "Dolichopodidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Cryptophleps nigrihalteratus"] <- "Cryptophleps"
#Cryptorhynchidius graniger
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Cryptorhynchidius graniger"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Cryptorhynchidius graniger"] <- "Curculionidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Cryptorhynchidius graniger"] <- "Cryptorhynchidius"
#Dasyhelea tamsi
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Dasyhelea tamsi"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Dasyhelea tamsi"] <- "Ceratopogonidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Dasyhelea tamsi"] <- "Dasyhelea"
#Delta alluaudi
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Delta alluaudi"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Delta alluaudi"] <- "Eumenidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Delta alluaudi"] <- "Delta"
#Dermestidae sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Dermestidae sp.1"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Dermestidae sp.1"] <- "Dermestidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Dermestidae sp.1"] <- NA
#Dermestidae sp.2
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Dermestidae sp.1"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Dermestidae sp.1"] <- "Dermestidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Dermestidae sp.1"] <- NA
#Dermestidae sp1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Dermestidae sp1"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Dermestidae sp1"] <- "Dermestidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Dermestidae sp1"] <- NA
#Dixidae sp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Dixidae sp."] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Dixidae sp."] <- "Dixidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Dixidae sp."] <- NA
#Dixidae sp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Dixidae sp."] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Dixidae sp."] <- "Formicidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Dixidae sp."] <-"Dolichoderus"
#Elateridae sp.2
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Elateridae sp.2"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Elateridae sp.2"] <- "Elateridae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Elateridae sp.2"] <- NA
#Elateridae sp.3
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Elateridae sp.3"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Elateridae sp.3"] <- "Elateridae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Elateridae sp.3"] <- NA
#Elateroidea sp.2
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Elateroidea sp.2"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Elateroidea sp.2"] <- "Elateridae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Elateroidea sp.2"] <- NA
#Elateroidea sp.3
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Elateroidea sp.3"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Elateroidea sp.3"] <- "Elateridae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Elateroidea sp.3"] <- NA
#Ematurga atomaria
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Ematurga atomaria"] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Ematurga atomaria"] <- "Geometridae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Ematurga atomaria"] <- "Ematurga"
#Endotricha mesenterialis
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Endotricha mesenterialis"] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Endotricha mesenterialis"] <- "Pyralidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Endotricha mesenterialis"] <- "Endotricha"
#Epanthidium bicoloratum
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Epanthidium bicoloratum"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Epanthidium bicoloratum"] <- "Megachilidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Epanthidium bicoloratum"] <- "Epanthidium"
#Epitragus sp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Epitragus sp."] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Epitragus sp."] <- "Ephydridae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Epitragus sp."] <- NA
#Ephydridae sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Ephydridae sp.1"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Ephydridae sp.1"] <- "Ephydridae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Ephydridae sp.1"] <- NA
#Epitragus sp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Epitragus sp."] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Epitragus sp."] <- "Tenebrionidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Epitragus sp."] <- "Epitragus"
#Eudiagogini sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Eudiagogini sp.1"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Eudiagogini sp.1"] <- "Curculionidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Eudiagogini sp.1"] <- "Eudiagogini"
#Eudiagogini sp.3
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Eudiagogini sp.3"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Eudiagogini sp.3"] <- "Curculionidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Eudiagogini sp.3"] <- "Eudiagogini"
#Eurycratus sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Eurycratus sp.1"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Eurycratus sp.1"] <- "Halictidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Eurycratus sp.1"] <- "Eurycratus"
#Eupetersia scotti
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Eupetersia scotti"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Eupetersia scotti"] <- "Halictidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Eupetersia scotti"] <- "Eupetersia"
#Exapion ulicis
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Exapion ulicis"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Exapion ulicis"] <- "Curculionidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Exapion ulicis"] <- "Exapion"
#Figitidae sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Figitidae sp.1"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Figitidae sp.1"] <- "Figitidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Figitidae sp.1"] <- NA
#Figitidae sp.2
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Figitidae sp.2"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Figitidae sp.2"] <- "Figitidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Figitidae sp.2"] <- NA
#Geraeus sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Geraeus sp.1"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Geraeus sp.1"] <- "Curculionidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Geraeus sp.1"] <- "Geraeus"
#Glutophrissa drusilla
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Glutophrissa drusilla"] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Glutophrissa drusilla"] <- "Pieridae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Glutophrissa drusilla"] <- "Appias"
#Gracilodes sp1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Gracilodes sp1"] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Gracilodes sp1"] <- "Erebidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Gracilodes sp1"] <- "Gracilodes"
#Hapalips sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Hapalips sp.1"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Hapalips sp.1"] <- "Erotylidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Hapalips sp.1"] <- "Hapalips"
#Hapalips sp1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Hapalips sp1"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Hapalips sp1"] <- "Erotylidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Hapalips sp1"] <- "Hapalips"
#Hasinamelissa sp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Hasinamelissa sp."] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Hasinamelissa sp."] <- "Apidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Hasinamelissa sp."] <- "Hasinamelissa"
#Hasinamelissa sp
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Hasinamelissa sp"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Hasinamelissa sp"] <- "Apidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Hasinamelissa sp"] <- "Hasinamelissa"
#Helorus coruscus
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Helorus coruscus"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Helorus coruscus"] <- "Heloridae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Helorus coruscus"] <- "Helorus"
#Hemidactylus frenatus
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Hemidactylus frenatus"] <- "Squamata"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Hemidactylus frenatus"] <- NA
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Hemidactylus frenatus"] <- NA
#Hemilucilia segmentaria
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Hemilucilia segmentaria"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Hemilucilia segmentaria"] <- "Calliphoridae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Hemilucilia segmentaria"] <- "Hemilucilia"
#Herpetogramma sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Herpetogramma sp.1"] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Herpetogramma sp.1"] <- "Crambidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Herpetogramma sp.1"] <- "Herpetogramma"
#Hesp.eriidae sp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Hesp.eriidae sp."] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Hesp.eriidae sp."] <- "Hesperiidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Hesp.eriidae sp."] <- NA
#Hesperiinae sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Hesperiinae sp.1"] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Hesperiinae sp.1"] <- "Hesperiidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Hesperiinae sp.1"] <- NA
#Hesperiinae sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Hesperiinae sp.1"] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Hesperiinae sp.1"] <- "Hesperiidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Hesperiinae sp.1"] <- NA
#Heterocera sp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Heterocera sp."] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Heterocera sp."] <- NA
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Heterocera sp."] <- NA
#Heterocera
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Heterocera"] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Heterocera"] <- NA
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Heterocera"] <- NA
#Hirmoneura sp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Hirmoneura sp."] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Hirmoneura sp."] <- "Nemestrinidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Hirmoneura sp."] <- "Hirmoneura"
#Hymenoptera sp.2
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Hymenoptera sp.2"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Hymenoptera sp.2"] <- NA
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Hymenoptera sp.2"] <- NA
#Hymenoptera sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Hymenoptera sp.1"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Hymenoptera sp.1"] <- NA
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Hymenoptera sp.1"] <- NA
#Hymenoptera sp.3
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Hymenoptera sp.3"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Hymenoptera sp.3"] <- NA
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Hymenoptera sp.3"] <- NA
#Hymenoptera sp.4
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Hymenoptera sp.4"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Hymenoptera sp.4"] <- NA
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Hymenoptera sp.4"] <- NA
#Hymenoptera sp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Hymenoptera sp."] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Hymenoptera sp."] <- NA
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Hymenoptera sp."] <- NA
#Hypenodinae sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Hypenodinae sp.1"] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Hypenodinae sp.1"] <- "Erebidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Hypenodinae sp.1"] <- NA
#Hypolimnas misippus
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Hypolimnas misippus"] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Hypolimnas misippus"] <- "Nymphalidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Hypolimnas misippus"] <- "Hypolimnas"
#Hypsipetes crassirostris
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Hypsipetes crassirostris"] <- "Passeriformes"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Hypsipetes crassirostris"] <- NA
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Hypsipetes crassirostris"] <- NA
#Hyrdotaea nigrisquama
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Hyrdotaea nigrisquama"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Hyrdotaea nigrisquama"] <- "Muscidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Hyrdotaea nigrisquama"] <- "Hyrdotaea"
#Inostemma sp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Inostemma sp."] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Inostemma sp."] <- "Platygastridae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Inostemma sp."] <- "Inostemma"
#Itame brunneata
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Itame brunneata"] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Itame brunneata"] <- "Geometridae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Itame brunneata"] <- "Macaria"
#Larocanthidium nigrilum
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Larocanthidium nigrilum"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Larocanthidium nigrilum"] <- "Megachilidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Larocanthidium nigrilum"] <- "Larocanthidium"
#Lassiopogon sp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Lassiopogon sp."] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Lassiopogon sp."] <- "Asilidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Lassiopogon sp."] <- "Lasiopogon"
#Leptometriella nigra
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Leptometriella nigra"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Leptometriella nigra"] <- "Apidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Leptometriella nigra"] <- "Leptometriella"
#Liosarcophaga spilargyra
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Liosarcophaga spilargyra"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Liosarcophaga spilargyra"] <- "Sarcophagagidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Liosarcophaga spilargyra"] <- "Liosarcophaga"
#Lobosciara bilotata
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Lobosciara bilotata"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Lobosciara bilotata"] <- "Sciaridae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Lobosciara bilotata"] <- "Lobosciara"
#Lonchoptera furcata
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Lonchoptera furcata"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Lonchoptera furcata"] <- "Lonchopteridae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Lonchoptera furcata"] <- "Lonchoptera"
#Melectoides cockerelli
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Melectoides cockerelli"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Melectoides cockerelli"] <- "Apidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Melectoides cockerelli"] <- "Melectoides"
#Mesonychium jenseni
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Mesonychium jenseni"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Mesonychium jenseni"] <- "Apidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Mesonychium jenseni"] <- "Mesonychium"
#Microlepidoptera
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Microlepidoptera"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Microlepidoptera"] <- NA
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Microlepidoptera"] <- NA
#Mischocyttarus sp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Mischocyttarus sp."] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Mischocyttarus sp."] <- "Vespidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Mischocyttarus sp."] <- "Mischocyttarus"
#near Eurycratus sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="near Eurycratus sp.1"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="near Eurycratus sp.1"] <- "Salpingidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="near Eurycratus sp.1"] <- NA
#near Eurycratus sp1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="near Eurycratus sp1"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="near Eurycratus sp1"] <- "Salpingidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="near Eurycratus sp1"] <- NA
#Nola sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Nola sp.1"] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Nola sp.1"] <- "Nolidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Nola sp1"] <- NA
#Nola sp1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Nola sp1"] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Nola sp1"] <- "Nolidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Nola sp1"] <- NA
#Orthonevra geniculata
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Orthonevra geniculata"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Orthonevra geniculata"] <- "Syrphidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Orthonevra geniculata"] <- "Othonevra"
#Oscinosoma speighti
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Oscinosoma speighti"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Oscinosoma speighti"] <- "Chloropidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Oscinosoma speighti"] <- "Oscinosoma"
#Pelecorhychus rubidus
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Pelecorhychus rubidus"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Pelecorhychus rubidus"] <- "Tabanoidea"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Pelecorhychus rubidus"] <- "Pelecorhychus"
#Physiphoa azurea
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Physiphoa azurea"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Physiphoa azurea"] <- "Ulidiidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Physiphoa azurea"] <- "Physiphora"
#Phystis simois variegata
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Phystis simois variegata"] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Phystis simois variegata"] <- "Nymphalidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Phystis simois variegata"] <- "Phystis"
#Pipunculidae
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Pipunculidae"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Pipunculidae"] <- "Pipunculidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Pipunculidae"] <- NA
#Platygastridae sp.2
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Platygastridae sp.2"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Platygastridae sp.2"] <- "Platygastridae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Platygastridae sp.2"] <- NA
#Platygastridae sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Platygastridae sp.1"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Platygastridae sp.1"] <- "Platygastridae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Platygastridae sp.1"] <- NA
#Plebeia catamarcensis
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Plebeia catamarcensis"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Plebeia catamarcensis"] <- "Apidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Plebeia catamarcensis"] <- "Plebeia"
#Poecilognathus sp.2
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Poecilognathus sp.2"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Poecilognathus sp.2"] <- "Bombyliidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Poecilognathus sp.2"] <- "Poecilognathus"
#Poecilognathus sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Poecilognathus sp.1"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Poecilognathus sp.1"] <- "Bombyliidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Poecilognathus sp.1"] <- "Poecilognathus"
#Poecilognathus sp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Poecilognathus sp."] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Poecilognathus sp."] <- "Bombyliidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Poecilognathus sp."] <- "Poecilognathus"
#Polemistus sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Polemistus sp.1"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Polemistus sp.1"] <- "Cabronidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Polemistus sp.1"] <- "Polemistus"
#Pompilidae sp.2
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Pompilidae sp.2"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Pompilidae sp.2"] <- "Andrenidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Pompilidae sp.2"] <- NA
#Psaenythia
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Pompilidae sp.2"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Pompilidae sp.2"] <- "Andrenidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Pompilidae sp.2"] <- "Psaenythia"
#Psaenythia rufipes
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Psaenythia rufipes"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Psaenythia rufipes"] <- "Andrenidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Psaenythia rufipes"] <- "Psaenythia"
#Pseudocentron curvipes
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Pseudocentron curvipes"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Pseudocentron curvipes"] <- "Megachilidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Pseudocentron curvipes"] <- "Pseudocentron"
#Pseudocentron sp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Pseudocentron sp."] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Pseudocentron sp."] <- "Megachilidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Pseudocentron sp."] <- "Pseudocentron"
#Pseudolycoriella setigera
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Pseudolycoriella setigera"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Pseudolycoriella setigera"] <- "Sciaridae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Pseudolycoriella setigera"] <- "Pseudolycoriella"
#Pseudolycoriella microcteniuni
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Pseudolycoriella microcteniuni"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Pseudolycoriella microcteniuni"] <- "Sciaridae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Pseudolycoriella microcteniuni"] <- "Pseudolycoriella"
#Pseudomyrmex gracilis
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Pseudomyrmex gracilis"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Pseudomyrmex gracilis"] <- "Formicidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Pseudomyrmex gracilis"] <- "Pseudomyrmex"
#Psylotrix viridicoerulea
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Psylotrix viridicoerulea"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Psylotrix viridicoerulea"] <- "Dasytidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Psylotrix viridicoerulea"] <- "Psylotrix"
#Ptiloglossa willinki
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Ptiloglossa willinki"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Ptiloglossa willinki"] <- "Colletidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Ptiloglossa willinki"] <- "Ptiloglossa"
#Pyrginae sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Pyrginae sp.1"] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Pyrginae sp.1"] <- "Hesperiidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Pyrginae sp.1"] <- NA
#Pyrisitia leuce leuce
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Pyrisitia leuce leuce"] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Pyrisitia leuce leuce"] <- "Pieridae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Pyrisitia leuce leuce"] <- "Pyristia"
#Pyrisitia nise
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Pyrisitia nise"] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Pyrisitia nise"] <- "Pieridae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Pyrisitia nise"] <- "Pyristia"
#Pyronota festiva
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Pyronota festiva"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Pyronota festiva"] <- "Scarabaeidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Pyronota festiva"] <- "Pyronota"
#Rhinia apicalis
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Rhinia apicalis"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Rhinia apicalis"] <- "Rhiniidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Rhinia apicalis"] <- "Rhinia"
#Riodinidae sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Riodinidae sp.1"] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Riodinidae sp.1"] <- "Riodinidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Riodinidae sp.1"] <- NA
#Sciomyzidae
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Sciomyzidae"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Sciomyzidae"] <- "Sciomyzidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Sciomyzidae"] <- NA
#Sciomyzidae sp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Sciomyzidae sp."] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Sciomyzidae sp."] <- "Sciomyzidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Sciomyzidae sp."] <- NA
#Sciomyzidae sp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Scolytidae sp.4"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Scolytidae sp.4"] <- "Curculionidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Scolytidae sp.4"] <- NA
#Scolytidae sp4
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Scolytidae sp4"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Scolytidae sp4"] <- "Curculionidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Scolytidae sp4"] <- NA
#Scolytidae sp.2
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Scolytidae sp.2"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Scolytidae sp.2"] <- "Curculionidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Scolytidae sp.2"] <- NA
#Scolytidae sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Scolytidae sp.1"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Scolytidae sp.1"] <- "Curculionidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Scolytidae sp.1"] <- NA
#Scolytidae sp.3
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Scolytidae sp.3"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Scolytidae sp.3"] <- "Curculionidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Scolytidae sp.3"] <- NA
#Simuliidae sp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Simuliidae sp."] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Simuliidae sp."] <- "Simuliidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Simuliidae sp."] <- "Simuliidae sp."
#Siricidae sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Siricidae sp.1"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Siricidae sp.1"] <- "Siricidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Siricidae sp.1"] <- NA
#Sirphidae sp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Sirphidae sp."] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Sirphidae sp."] <- "Syrphidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Sirphidae sp."] <- NA
#Smitta mahensis
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Smitta mahensis"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Smitta mahensis"] <- "Chironomidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Smitta mahensis"] <- "Smitta"
#Spintherophytasp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Spintherophytasp."] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Spintherophytasp."] <- "Chrysomelidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Spintherophytasp."] <- "Spintherophyta"
#Staphylinidae sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Staphylinidae sp.1"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Staphylinidae sp.1"] <- "Staphylinidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Staphylinidae sp.1"] <- NA
#Staphylinidae sp.2
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Staphylinidae sp.2"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Staphylinidae sp.2"] <- "Staphylinidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Staphylinidae sp.2"] <- NA
#Staphylinidae sp1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Staphylinidae sp1"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Staphylinidae sp1"] <- "Staphylinidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Staphylinidae sp1"] <- NA
#Stenomorda near disparilis
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Stenomorda near disparilis"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Stenomorda near disparilis"] <- "Mordellidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Stenomorda near disparilis"] <- "Stenomorda"
#Stenopterus rufus
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Stenopterus rufus"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Stenopterus rufus"] <- "Cerambycidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Stenopterus rufus"] <- "Stenopturus"
#Stenopterus rufus
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Stenotrupis sp.1"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Stenotrupis sp.1"] <- "Curculionidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Stenotrupis sp.1"] <- "Stenotrupis"
#Stenopterus rufus
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Stenurella melanura"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Stenurella melanura"] <- "Cerambycidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Stenurella melanura"] <- "Stenurella"
#Stictopterinae sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Stictopterinae sp.1"] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Stictopterinae sp.1"] <- "Noctuidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Stictopterinae sp.1"] <-  NA
#Tachipompillus sp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Tachipompillus sp."] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Tachipompillus sp."] <- "Pompilidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Tachipompillus sp."] <-  "Tachipompillus"
#Tatochila autodice
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Tatochila autodice"] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Tatochila autodice"] <- "Pieridae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Tatochila autodice"] <-  "Tatochila"
#Tatochila autodice
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Tetragonisca angustula"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Tetragonisca angustula"] <- "Apidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Tetragonisca angustula"] <-  "Tetragonisca"
#Tetragonoschema missionarium
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Tetragonoschema missionarium"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Tetragonoschema missionarium"] <- "Buprestidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Tetragonoschema missionarium"] <-  "Tetragonoschema"
#Thripomorpha rufithorax
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Thripomorpha rufithorax"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Thripomorpha rufithorax"] <- "Scatopsidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Thripomorpha rufithorax"] <-  "Thripomorpha"
#Thynnidae sp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Thynnidae sp."] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Thynnidae sp."] <- "Thynnidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Thynnidae sp."] <-  NA
#Toxophora aurea
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Toxophora aurea"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Toxophora aurea"] <- "Bombyliidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Toxophora aurea"] <-  "Toxophora"
#Trialeurodes sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Trialeurodes sp.1"] <- "Hemiptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Trialeurodes sp.1"] <- "Aleyrodidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Trialeurodes sp.1"] <-  "Trialeurodes"
#Trupanea alboapicata
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Trupanea alboapicata"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Trupanea alboapicata"] <- "Tephritidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Trupanea alboapicata"] <-  "Trupanea"
#Vesp.idae sp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Vesp.idae sp."] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Vesp.idae sp."] <- "Vespidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Vesp.idae sp."] <-  NA
#Vesp.idae sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Vesp.idae sp.1"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Vesp.idae sp.1"] <- "Vespidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Vesp.idae sp.1"] <-  NA
#Xenocalliphora hortona
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Xenocalliphora hortona"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Xenocalliphora hortona"] <- "Calliphoridae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Xenocalliphora hortona"] <- "Xenocalliphora"
#Xenosciomyza subalpina
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Xenosciomyza subalpina"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Xenosciomyza subalpina"] <- "Helosciomyzidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Xenosciomyza subalpina"] <- "Xenosciomyza"
#Yphthimoides celmis
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Yphthimoides celmis"] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Yphthimoides celmis"] <- "Nymphalidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Yphthimoides celmis"] <- "Yphthimoides"
#Zeta sp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Zeta sp."] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Zeta sp."] <- "Vespidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Zeta sp."] <- "Zeta"
#Zizina sp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Zizina sp."] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Zizina sp."] <- "Lycaenidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Zizina sp."] <- "Zizinia"

#Acromyrmex lobicornis
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Acromyrmex lobicornis"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Acromyrmex lobicornis"] <- "Formicidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Acromyrmex lobicornis"] <- "Acromyrmex"

#Audre sp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Audre sp."] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Audre sp."] <- "Papilionidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Audre sp."] <- "Audre"

#Chaetophthalamus dorsalis
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Chaetophthalamus dorsalis"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Chaetophthalamus dorsalis"] <- NA
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Chaetophthalamus dorsalis"] <- NA

#Danaidae sp..
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Danaidae sp."] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Danaidae sp."] <- "Nymphalidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Danaidae sp."] <- NA

#Dermestidae sp.2
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Dermestidae sp.2"] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Dermestidae sp.2"] <- "Dermestidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Dermestidae sp.2"] <- NA

#Dolichoderus bispinosus
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Dolichoderus bispinosus"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Dolichoderus bispinosus"] <- "Formicidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Dolichoderus bispinosus"] <- "Dolichoderus"

#Ephydridae
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Ephydridae"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Ephydridae"] <- "Ephydridae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Ephydridae"] <- NA

#Eudaeus sp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Eudaeus sp."] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Eudaeus sp."] <- "Curculionidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Eudaeus sp."] <- NA

#Exonera sp.1
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Exonera sp.1"] <- "Hymenoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Exonera sp.1"] <- NA
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Exonera sp.1"] <- NA

#near Eurycratus sp1 
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="near Eurycratus sp1 "] <- "Coleoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="near Eurycratus sp1 "] <- NA
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="near Eurycratus sp1 "] <- NA

#Physiphora azurea
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Physiphora azurea"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Physiphora azurea"] <- "Ulidiidae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Physiphora azurea"] <- "Physiphora"

#Pieris sp.
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Pieris sp."] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Pieris sp."] <- "Pieridae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Pieris sp."] <- "Pieris"

#Pieris rapae
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Pieris rapae"] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Pieris rapae"] <- "Pieridae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Pieris rapae"] <- "Pieris"

#Pieris rapae
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Pieris napi"] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Pieris napi"] <- "Pieridae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Pieris napi"] <- "Pieris"

#Pieris rapae
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Pieris sp.1"] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Pieris sp.1"] <- "Pieridae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Pieris sp.1"] <- "Pieris"

#Pieris brassicae
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Pieris brassicae"] <- "Lepidoptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Pieris brassicae"] <- "Pieridae"
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Pieris brassicae"] <- "Pieris"

#Prohardya carinata
all_long_poll_names$order[all_long_poll_names$Pollinator_species=="Prohardya carinata"] <- "Diptera"
all_long_poll_names$family[all_long_poll_names$Pollinator_species=="Prohardya carinata"] <- NA
all_long_poll_names$genus[all_long_poll_names$Pollinator_species=="Prohardya carinata"] <- NA





#remove network "mistakes"
all_long_poll_names <- all_long_poll_names[ grep("Agrostis atritegulator", all_long_poll_names$Pollinator_species, invert = TRUE) , ]
all_long_poll_names <- all_long_poll_names[ grep("Alophera sp.1", all_long_poll_names$Pollinator_species, invert = TRUE) , ]
all_long_poll_names <- all_long_poll_names[ grep("Anosymia sp.", all_long_poll_names$Pollinator_species, invert = TRUE) , ]
all_long_poll_names <- all_long_poll_names[ grep("Mydeae sp.", all_long_poll_names$Pollinator_species, invert = TRUE) , ]
#removing earthworms
all_long_poll_names <- all_long_poll_names[ grep("Panara naso", all_long_poll_names$Pollinator_species, invert = TRUE) , ]
#removing Unidentified
all_long_poll_names <- all_long_poll_names[ grep("sp.1", all_long_poll_names$Pollinator_species, invert = TRUE) , ]
all_long_poll_names <- all_long_poll_names[ grep("sp.2", all_long_poll_names$Pollinator_species, invert = TRUE) , ]
all_long_poll_names <- all_long_poll_names[ grep("sp.3", all_long_poll_names$Pollinator_species, invert = TRUE) , ]
all_long_poll_names <- all_long_poll_names[ grep("sp.4", all_long_poll_names$Pollinator_species, invert = TRUE) , ]
all_long_poll_names <- all_long_poll_names[ grep("sp.5", all_long_poll_names$Pollinator_species, invert = TRUE) , ]
all_long_poll_names <- all_long_poll_names[ grep("Unientified sp.1", all_long_poll_names$Pollinator_species, invert = TRUE) , ]
all_long_poll_names <- all_long_poll_names[ grep("Unidentified sp.1", all_long_poll_names$Pollinator_species, invert = TRUE) , ]
all_long_poll_names <- all_long_poll_names[ grep("Unidentified sp.5", all_long_poll_names$Pollinator_species, invert = TRUE) , ]
all_long_poll_names <- all_long_poll_names[ grep("Unidentified sp.8", all_long_poll_names$Pollinator_species, invert = TRUE) , ]
all_long_poll_names <- all_long_poll_names[ grep("Unientified sp.2", all_long_poll_names$Pollinator_species, invert = TRUE) , ]
all_long_poll_names <- all_long_poll_names[ grep("Unientified sp.3", all_long_poll_names$Pollinator_species, invert = TRUE) , ]
all_long_poll_names <- all_long_poll_names[ grep("Unidentified sp.6", all_long_poll_names$Pollinator_species, invert = TRUE) , ]
all_long_poll_names <- all_long_poll_names[ grep("Unidentified sp.7", all_long_poll_names$Pollinator_species, invert = TRUE) , ]
all_long_poll_names <- all_long_poll_names[ grep("Unidentified sp.2", all_long_poll_names$Pollinator_species, invert = TRUE) , ]
all_long_poll_names <- all_long_poll_names[ grep("Unidentified sp.3", all_long_poll_names$Pollinator_species, invert = TRUE) , ]
all_long_poll_names <- all_long_poll_names[ grep("Unidentified sp.4", all_long_poll_names$Pollinator_species, invert = TRUE) , ]
all_long_poll_names <- all_long_poll_names[ grep("Unientified sp.4", all_long_poll_names$Pollinator_species, invert = TRUE) , ]
all_long_poll_names <- all_long_poll_names[ grep("Unidentified sp.9", all_long_poll_names$Pollinator_species, invert = TRUE) , ]


all_long_poll_names$guild <- NA

#Add pollinator guilds for analysis
all_long_poll_names$guild[all_long_poll_names$order=="Lepidoptera"] <- "Lepidoptera"
all_long_poll_names$guild[all_long_poll_names$order=="Coleoptera"] <- "Coleoptera"
all_long_poll_names$guild[all_long_poll_names$order=="Diptera"] <- "Non-syrphids-diptera"
all_long_poll_names$guild[all_long_poll_names$family=="Syrphidae"] <- "Syrphids"
all_long_poll_names$guild[all_long_poll_names$order=="Hymenoptera"] <- "Non-bee-Hymenoptera"
all_long_poll_names$guild[all_long_poll_names$family=="Apidae"] <- "Bee"
all_long_poll_names$guild[all_long_poll_names$family=="Megachilidae"] <- "Bee"
all_long_poll_names$guild[all_long_poll_names$family=="Halictidae"] <- "Bee"
all_long_poll_names$guild[all_long_poll_names$family=="Andrenidae"] <- "Bee"
all_long_poll_names$guild[all_long_poll_names$family=="Colletidae"] <- "Bee"
all_long_poll_names$guild[all_long_poll_names$family=="Melittidae"] <- "Bee"
all_long_poll_names$guild[all_long_poll_names$family=="Stenotritidae"] <- "Bee"
#Now the other guilds

all_long_poll_names$guild[all_long_poll_names$Pollinator_species=="Acarina sp.1"] <- "Other_insects"
all_long_poll_names$guild[all_long_poll_names$order=="Neuroptera"] <- "Other_insects"
all_long_poll_names$guild[all_long_poll_names$order=="Neuroptera"] <- "Other_insects"
all_long_poll_names$guild[all_long_poll_names$order=="Neuroptera"] <- "Other_insects"
all_long_poll_names$guild[all_long_poll_names$order=="Hemiptera"] <- "Other_insects"
all_long_poll_names$guild[all_long_poll_names$order=="Orthoptera"] <- "Other_insects"
all_long_poll_names$guild[all_long_poll_names$order=="Dermaptera"] <- "Other_insects"
all_long_poll_names$guild[all_long_poll_names$order=="Blattodea"] <- "Other_insects"
all_long_poll_names$guild[all_long_poll_names$order=="Isopoda"] <- "Other_insects"
all_long_poll_names$guild[all_long_poll_names$order=="Araneae"] <- "Other_insects"
all_long_poll_names$guild[all_long_poll_names$order=="Thysanoptera"] <- "Other_insects"


all_long_poll_names$guild[all_long_poll_names$order=="Apodiformes"] <- "Birds"
all_long_poll_names$guild[all_long_poll_names$Pollinator_species=="Cinnyris dussumieri"] <- "Birds"
all_long_poll_names$guild[all_long_poll_names$order=="Passeriformes"] <- "Birds"

all_long_poll_names$guild[all_long_poll_names$order=="Squamata"] <- "Lizards"


na <- all_long_poll_names[is.na(all_long_poll_names$guild),]
#check for na

#Aggregate by poll guild
#all_poll <- reshape2::dcast(Plant_species + guild +Id ~ "Interaction", value.var = "Interaction", fun.aggregate = sum, data = all_long_poll_names, na.rm= TRUE)
#head(all_poll)
#b <- subset(all_long_poll_names, guild=="Bee")
#table(b$genus)

########################################################################################################################################################
#3)SAVE DATA
########################################################################################################################################################
#save poll guild data
write.csv(all_long_poll_names, "Data/Csv/long_format_quantitative_networks.csv")
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################



