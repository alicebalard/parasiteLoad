##http://onlinelibrary.wiley.com/doi/10.1111/j.1558-5646.2012.01633.x/full

## Test on worms data (local alice) :
Joelle_data <- na.omit(read.csv("../../../EvolutionFinalData.csv"))

###############################
## Prevalence of the bigger groups 
prev_and_mean <- function(x) {
  prev = sum(table(x)[names(table(x)) != 0]) /
    sum(table(x)) * 100
  mean = mean(x)
  return(paste0("Prevalence of ", round(prev,2), "%, mean of ", round(mean,2), " per mouse"))
}

# caecal nematode pinworms (A. tetraptera and S. obvelata): 70.9% of mice were infected
# with an average load of 39.2 pinworms per mouse 
prev_and_mean(Joelle_data$Aspiculuris.Syphacia)

# second most widespread infection was due to the caecal nematode Trichuris muris (whipworm)
#found in 21.1% of mice with an average load of 3.9 individuals per mouse (Fig. 2D).
prev_and_mean(Joelle_data$Trichuris)

#The larval stage of Taenia taeniaeformis (tapeworm) was the most prevalent cestode,
# found in the liver of 10.7% of individuals (Fig. 2E). 
prev_and_mean(Joelle_data$Taenia)

# The fourth most prevalent parasite was Mastophorus muris, 
#found in the stomach of 9.4% of individuals (Fig. 2F).
prev_and_mean(Joelle_data$Mastophorus)

###############################
## Which distribution does fit?
hist(Joelle_data$Aspiculuris.Syphacia)
