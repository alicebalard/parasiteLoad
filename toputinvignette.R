toyData <- read.csv("https://raw.githubusercontent.com/alicebalard/Article_IntensityEimeriaHMHZ/master/data/cleanedData.csv")
toyData <-toyData[toyData$Aspiculuris_Syphacia != 0,]

m <- mlHyb(Aspiculuris_Syphacia~., data = toyData, model = "negbin", hybridEffect = TRUE)
m

m0 <- mlHyb(Aspiculuris_Syphacia~., data = toyData, model = "negbin",
                          hybridEffect = FALSE)
m0
