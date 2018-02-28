source("MLE_hybrid_functions.R")

## Import data WATWM
Joelle_data <- read.csv("../examples/Reproduction_WATWM/EvolutionFinalData.csv")
# to check
Joelle_data <- Joelle_data[complete.cases(Joelle_data$HI),]
# pinworms (A. tetraptera and S. obvelata)
# Trichuris muris (whipworm)
dataTrichuris <- Joelle_data[complete.cases(Joelle_data$Trichuris),]
dataTrichuris_F <- dataTrichuris[dataTrichuris$Sex == "F",]
dataTrichuris_M <- dataTrichuris[dataTrichuris$Sex == "M",]
dataTrichuris_F$Sex <- droplevels(dataTrichuris_F$Sex)
dataTrichuris_M$Sex <- droplevels(dataTrichuris_M$Sex)

# Taenia taeniaeformis (tapeworm) 

##### Input end #####
paramBounds <- c(L1start = 10, L1LB = 0, L1UB = 700, 
                 L2start = 10, L2LB = 0, L2UB = 700, 
                 alphaStart = 0, alphaLB = -5, alphaUB = 5,
                 A1start = 10, A1LB = 0, A1UB = 1000, 
                 A2start = 10, A2LB = 0, A2UB = 1000, 
                 Zstart = 0, ZLB = -5, ZUB = 5)

TrichurisFit <- fit(data = dataTrichuris, response = "Trichuris",
                    hybridIndex = HI, paramBounds = paramBounds)

TrichurisFit_F <- fit(data = dataTrichuris_F, response = "Trichuris",
                      hybridIndex = HI, paramBounds = paramBounds)

TrichurisFit_M <- fit(data = dataTrichuris_M, response = "Trichuris",
                      hybridIndex = HI, paramBounds = paramBounds)

####### Is alpha significant for each hypothesis?

# H0: the expected load for the subspecies and between sexes is the same
isAlphaSignif(TrichurisFit$fitAlphaBasic, TrichurisFit$fitNoAlphaBasic)
H0 <- TrichurisFit$fitAlphaBasic

# H1: the mean load across sexes is the same, but can differ across subspecies
isAlphaSignif(TrichurisFit$fitAlphaAdvanced, TrichurisFit$fitNoAlphaAdvanced)
H1 <- TrichurisFit$fitAlphaAdvanced

# H2: the mean load across subspecies is the same, but can differ between the sexes
isAlphaSignif(TrichurisFit_F$fitAlphaBasic, TrichurisFit_F$fitNoAlphaBasic)
isAlphaSignif(TrichurisFit_M$fitAlphaBasic, TrichurisFit_M$fitNoAlphaBasic)
H2 <- list(female = TrichurisFit_F$fitAlphaBasic,
           male = TrichurisFit_M$fitAlphaBasic)

# H3: the mean load can differ both across subspecies and between sexes
isAlphaSignif(TrichurisFit_F$fitAlphaAdvanced, TrichurisFit_F$fitNoAlphaAdvanced)
isAlphaSignif(TrichurisFit_M$fitAlphaAdvanced, TrichurisFit_M$fitNoAlphaAdvanced)
H3 <- list(female = TrichurisFit_F$fitAlphaAdvanced,
           male = TrichurisFit_M$fitAlphaAdvanced)

####### Compare the hypotheses with G-tests 
# H1 vs H0
Gtest(model0 = H0, model1 = H1)

# H2 vs H0
Gtest(model0 = H0, model1 = H2)

# H3 vs H1
Gtest(model0 = H1, model1 = H3)

# H3 vs H2
Gtest(model0 = H2, model1 = H3)
