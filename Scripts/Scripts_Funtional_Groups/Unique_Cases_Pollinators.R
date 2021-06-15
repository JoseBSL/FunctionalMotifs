
library(dplyr)

final_d_1 <- read.csv("Data/Csv/data_for_motifs_analysis_1.csv")

distinct_poll <- final_d_1 %>% distinct(Pollinator_species, .keep_all = T)

poll_info <- distinct_poll %>% select(Pollinator_order, Pollinator_family, Pollinator_genus,Pollinator_species, Pollinator_functional_group)

write.csv(poll_info, "Data/Csv/unique_cases_poll_species.csv")

str(poll_info)
