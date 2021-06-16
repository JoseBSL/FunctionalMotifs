


library("FactoMineR")
library("factoextra")

data(housetasks)

dt <- as.table(as.matrix(housetasks))

chisq <- chisq.test(housetasks)


res.ca <- CA(housetasks, graph = FALSE)

fviz_ca_biplot(res.ca, repel = TRUE)


fviz_ca_biplot(res.ca, 
               map ="rowprincipal", arrow = c(TRUE, TRUE),
               repel = TRUE)



######READ NOW MY DATA######
library(tidyverse)
library(igraph)
library(bmotif)
networks <- read_csv("Data/Csv/data_for_motifs_analysis_1.csv") %>% select(-X1) %>%
  filter(Interaction >= int.threshold)

str(networks)


d = networks %>% select(Pollinator_functional_group,Plant_functional_groups)

library(reshape2)

d_tab <- dcast(data=d,
      Pollinator_functional_group ~ Plant_functional_groups,
      fun.aggregate = length,
      value.var = "Plant_functional_groups")

rownames(d_tab) <- d_tab$Pollinator_functional_group

d_tab = d_tab %>% select(-Pollinator_functional_group)

d_mat <- as.table(as.matrix(d_tab))

chisq <- chisq.test(d_mat)

res.ca <- CA(d_mat, graph = FALSE)

fviz_ca_biplot(res.ca, 
               map ="rowprincipal", arrow = c(TRUE, TRUE),
               repel = TRUE)
fviz_ca_biplot(res.ca, map ="colgreen", arrow = c(TRUE, FALSE),
               repel = TRUE)

fviz_ca_biplot(res.ca, repel = TRUE)
