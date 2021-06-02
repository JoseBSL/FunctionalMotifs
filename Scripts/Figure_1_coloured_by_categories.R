library(lme4)
library(nlme)
library(tidyverse)
source("Scripts/Scripts_Alfonso/add_study_id.R")
d <- read_csv("Data/Csv/Motifs_frequencies_and_null_models/Motifs_frequency_percentile.csv")
d <- add_study_id(d)
#str(d)
#head(d)

#ggplot(d %>% filter(motif != 1), aes(percentil_sizeclass))+
#  geom_histogram(color="black", fill="white")+
#  facet_wrap(~motif)+
#  theme_bw()+
#  labs(x="Sizeclass percentile")
#
# Motifs tend to be over- and under-represensented in real networks

# Mean percentile and SE taking into account the study system
motif_codes <- d %>% filter(motif != 1) %>% # We remove links (motif code: 1)
  select(motif) %>% unique() %>% pull()
motif_means <- tibble(motif = motif_codes)
motif_means$mean <- NA
motif_means$SE <- NA

for(i.motif in 1:length(motif_codes)){
  
  m <- lmer(percentil_sizeclass ~ 1+(1|study_id), #HERE WE NEED STUDY SYSTEM ONLY!! 
            data = subset(d, motif == motif_codes[i.motif]))  
  
  motif_means$mean[i.motif] <- fixed.effects(m) %>% unname()
  motif_means$SE[i.motif] <- sqrt(diag(vcov(m)))
  
}


d$Broad_categories 

d$Broad_categories[d$motif==2 | d$motif==7 | d$motif==17 | d$motif==3|d$motif==4|d$motif==8] <- "Fan" 

d$Broad_categories[d$motif==15 | d$motif==11] <- "Medium-weak" 

d$Broad_categories[d$motif==6 | d$motif==16 |  d$motif==12] <- "Strong" 

d$Broad_categories[d$motif==5 | d$motif==9 |  d$motif==14|  d$motif==13|  d$motif==10] <- "Weak" 


ggplot(NULL) + 
  geom_point(data = d %>% filter(motif != 1),
             aes(y=as.factor(motif), x=percentil_sizeclass, 
                 color = as.factor(Broad_categories)),
             position = "jitter",alpha=0.5)+
  geom_errorbar(data = motif_means,aes(y = as.factor(motif), xmin=mean-SE, xmax=mean+SE), 
                width=1.0,size=1)+
  geom_point(data = motif_means,aes(y = as.factor(motif), x=mean), 
             size=2)+
  labs(y="Motif", x = "Percentile")+
  theme_bw()

#Fig 1 Should be a dot plot (or forest plot) of the mean +- SE of the 15 motifs.
#NICE! Can we put motifs in Y axes, and % in X?