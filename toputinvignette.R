toyData <- read.csv("https://raw.githubusercontent.com/alicebalard/Article_IntensityEimeriaHMHZ/master/data/cleanedData.csv")
toyData <-toyData[toyData$Aspiculuris_Syphacia != 0,]

m <- mlHyb(Aspiculuris_Syphacia~., data = toyData, model = "negbin", hybridEffect = TRUE)
m

m0 <- mlHyb(Aspiculuris_Syphacia~., data = toyData, model = "negbin", hybridEffect = FALSE)
m0

lrt.mlHyb(m)
##### to do

# groups
mlHyb(Aspiculuris_Syphacia~Sex, data = toyData, model = "negbin", hybridEffect = TRUE)
mlHyb(Aspiculuris_Syphacia~., data = toyData, model = "negbin", hybridEffect = TRUE)

# fix sides?

# print table 4 hypothesis

# plots
# plot(lm(Aspiculuris_Syphacia~Sex, data = toyData)) Ca ira


