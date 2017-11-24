##http://onlinelibrary.wiley.com/doi/10.1111/j.1558-5646.2012.01633.x/full

## Test on worms data (local alice) :
Joelle_data <- read.csv("../../../EvolutionFinalData.csv")

###############################
# parasite prevalence (proportion of infected hosts among all the hosts 
# examined) and mean load (the mean number of parasites in all hosts, including
# the zero values of uninfected hosts) with 95% confidence limits over the
# complete dataset (Rózsa et al. 2000). Confidence intervals for prevalence and 
# load were estimated using Sterne's Exact method (Reiczigel 2003) and bootstrap
# (1000 replicates), respectively.

# Using QP3.0, we modeled the aggregation of parasites within hosts using the 
# negative binomial distribution, calculating maximum likelihood estimates 
# (MLEs) of parameter k (Bliss and Fisher 1953). QP3.0 also calculates a χ2 
# goodness of fit between expected and observed frequencies under the negative
# binomial model. 

library(exactci)

fit_neg_bin <- function(x) {
  x = as.vector(na.omit(x))
  pos = sum(table(x)[names(table(x)) != 0])
  n = sum(table(x))
  mean_intensity = mean(na.omit(x))
  prev = round(pos/n * 100,2)
  CI_low = binom.exact(pos, n, p = pos/n,
                       alternative = c("two.sided", "less", "greater"),
                       tsmethod = c("minlike"),
                       conf.level = 0.95)$conf.int[1]*100
  CI_up = binom.exact(pos, n, p = pos/n,
                      alternative = c("two.sided", "less", "greater"),
                      tsmethod = c("minlike"),
                      conf.level = 0.95)$conf.int[2]*100
  fit = fitdistr(x, densfun = "negative binomial")
  # Goodness of fit with the chi squared test  
  df = as.data.frame(table(x))
  names(df) = c("values", "observed_freq")
  df$values = as.numeric(as.character(df$values))
  df$expected_freq = dnbinom(df$values, size = fit$estimate[1],
                             mu = fit$estimate[2]) * sum(df$observed_freq)
  Chi2 = chisq.test(x = df$observed_freq, y = df$expected_freq)
  # Plot the histogram & fitted distribution
  hist(x, breaks = 50, freq = F)
  fitD = dnbinom(0:length(x), size = fit$estimate[1], mu = fit$estimate[2])
  lines(fitD, lwd="3", col="red")
  return(list(pos = pos, n = n, mean_intensity = mean_intensity,
              prev = prev, CI_low = CI_low, CI_up = CI_up, fit.neg.bin = fit,
              Chi2 = Chi2))
}

# caecal nematode pinworms (A. tetraptera and S. obvelata): 70.9% of mice were
# infected with an average load of 39.2 pinworms per mouse 
fit_neg_bin(Joelle_data$Aspiculuris.Syphacia)

# second most widespread infection was due to the caecal nematode Trichuris 
# muris (whipworm) found in 21.1% of mice with an average load of 3.9 
# individuals per mouse (Fig. 2D).
fit_neg_bin(Joelle_data$Trichuris)

# The larval stage of Taenia taeniaeformis (tapeworm) was the most prevalent 
# cestode, found in the liver of 10.7% of individuals (Fig. 2E). 
fit_neg_bin(Joelle_data$Taenia)

# The fourth most prevalent parasite was Mastophorus muris, found in the stomach
# of 9.4% of individuals (Fig. 2F).
fit_neg_bin(Joelle_data$Mastophorus)

