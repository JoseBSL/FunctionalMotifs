
library(tidyverse)

# https://bartomeuslab.com/2014/12/17/preferring-a-preference-index/

chi_pref <- function(obs, exp, alpha = 0.05){
  chi <- chisq.test(obs, p = exp, rescale.p = TRUE, simulate.p.value = F, B = 2000)
  print(chi) #tells you if there is an overall preference. (sig = pref)
  res <- chi$residuals
  #res <- (obs-exp)/sqrt(exp) #hand calculation, same result.
  #calculate bonferoni Z-statistic for each plant.
  alpha <- alpha
  k <- length(obs)
  n <- sum(obs)
  p_obs <- obs/n
  ak <- alpha/(2*k)
  Zak <- abs(qnorm(ak))
  low_interval <- p_obs - (Zak*(sqrt(p_obs*(1-p_obs)/n)))
  upper_interval <- p_obs + (Zak*(sqrt(p_obs*(1-p_obs)/n)))
  p_exp <- exp/sum(exp)
  sig <- ifelse(p_exp >= low_interval & p_exp <= upper_interval, "ns", "sig")
  plot(c(0,k+1), c(min(low_interval),max(upper_interval)), type = "n", 
       ylab = "Preference", xlab = "items", las = 1)
  arrows(x0 = c(1:k), y0 = low_interval, x1 = c(1:k), y1 = upper_interval, code = 3
         ,angle = 90)
  points(p_exp, col = "red")
  out <- data.frame(chi_test_p = rep(chi$p.value, length(res)), 
                    chi_residuals = res, sig = sig)
  out
}


motifs_observed_probability <- read_csv("Data/Csv/motifs_observed_probability.csv")
motifs_expected_probability <- read_csv("Data/Csv/node_motifs_theoretical_probability.csv") %>%
  select(motif,motif_functional_ID,motif_probability) %>% unique() %>%
  rename(motif_expected_probability = motif_probability)

motifs_probability <- motifs_expected_probability %>% 
  left_join(motifs_observed_probability, by = c("motif","motif_functional_ID"))

motifs_probability[is.na(motifs_probability)] <- 0

motifs_probability <- motifs_probability %>% arrange(desc(motif_observed_probability))

obs <- motifs_probability$counts_observed
exp <- motifs_probability$motif_expected_probability

chi_test <- chi_pref(obs, exp, alpha = 0.05)

res <- chisq.test(obs, p = exp, rescale.p = TRUE, simulate.p.value = T, B = 2000)
res 
res$expected %>% as.vector()
res$observed %>% as.vector()

# Null hypothesis (H0): There is no significant difference between the observed and the 
# expected value.
# The p-value of the test is less than the significance level alpha = 0.05. 
# We can conclude that the motifs are significantly not distributed as expected 
# with a p-value = 2.2e-16.

motifs_probability$chi_test_p <- chi_test$chi_test_p
motifs_probability$chi_residuals <- chi_test$chi_residuals
motifs_probability$sig <- chi_test$sig

motifs_probability <- motifs_probability %>% arrange(desc(sig),desc(abs(chi_residuals)))

# expected values
chi_pref(motifs_probability$counts_observed,
         motifs_probability$motif_expected_probability,
         alpha = 0.05)

hist(motifs_probability$counts_observed[motifs_probability$counts_observed>1000],300)


chi_pref(obs = c(0,25,200), exp = c(0.00003,50,100),alpha = 0.05)


Convictions <- matrix(c(2, 10, 15, 3), nrow = 2,
                      dimnames =
                        list(c("Dizygotic", "Monozygotic"),
                             c("Convicted", "Not convicted")))


#######################
# Note that, the chi-square test should be used only when all calculated expected 
# values are greater than 5. In our case there are several categories equal to zero

expected <- sum(motifs_probability$counts_observed)*
  motifs_probability$motif_expected_probability
expected[expected < 5] 
