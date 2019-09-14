library(parasiteLoad)

toyData <- read.csv("data/BALdata.csv") # update with final published table later
toyData <-toyData[toyData$Aspiculuris_Syphacia != 0,]

# Export
library(optimx)

## needed
# FitnegBin
# gtest
m <- mlHyb(Aspiculuris_Syphacia~., data = toyData, model = "negbin")
m
