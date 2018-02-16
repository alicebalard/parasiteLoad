source("../R/WithMLE2.R")

## Import data WATWM
Joelle_data <- read.csv("../examples/Reproduction_WATWM/EvolutionFinalData.csv")
# to check
Joelle_data[is.na(Joelle_data)] <- 0
# pinworms (A. tetraptera and S. obvelata)
# Trichuris muris (whipworm)
# Taenia taeniaeformis (tapeworm) 

################## Data analysis ################## 
system.time(ModelPinworm <- myFun(data = Joelle_data, 
                                  hybridIndex = HI, 
                                  response = "Aspiculuris.Syphacia"))

system.time(ModelWhipworm <- myFun(data = Joelle_data, 
                                   hybridIndex = HI, 
                                   response = "Trichuris"))

system.time(ModelTapeworm <- myFun(data = Joelle_data, 
                                   hybridIndex = HI, 
                                   response = "Taenia"))

system.time(ModelMasto <- myFun(data = Joelle_data, 
                                hybridIndex = HI, 
                                response = "Mastophorus"))


finalFun(ModelPinworm)

################## Plotting ################## 

## profile investigates behavior of objective function near the MLE
system.time(myProf <- profile(ModelPinworm$H1$fitAlpha))

## Marginal confidence interval
myConfInt <- confint(myProf)

## Marginal confidence interval for alpha
alphaCILB <- myConfInt[rownames(myConfInt) == "alpha"][1]
alphaCIUB <- myConfInt[rownames(myConfInt) == "alpha"][2]

## Draw the line for the parameters at their MLE, alpha varying
DF <- data.frame(HI = seq(0,1,0.01),
                 loadMLE = MeanLoad(L1 = coef(ModelPinworm$H1$fitAlpha)[names(coef(ModelPinworm$H1$fitAlpha)) == "L1"],
                                    L2 =  coef(ModelPinworm$H1$fitAlpha)[names(coef(ModelPinworm$H1$fitAlpha)) == "L2"],
                                    alpha =  coef(ModelPinworm$H1$fitAlpha)[names(coef(ModelPinworm$H1$fitAlpha)) == "alpha"], 
                                    hybridIndex = seq(0,1,0.01)),
                 loadMLEnoAlpha = MeanLoad(L1 = coef(ModelPinworm$H1$fitNoAlpha)[names(coef(ModelPinworm$H1$fitNoAlpha)) == "L1"],
                                           L2 =  coef(ModelPinworm$H1$fitNoAlpha)[names(coef(ModelPinworm$H1$fitNoAlpha)) == "L2"],
                                           alpha =  coef(ModelPinworm$H1$fitNoAlpha)[names(coef(ModelPinworm$H1$fitNoAlpha)) == "alpha"], 
                                           hybridIndex = seq(0,1,0.01)),
                 loadMLEAlphaLB = MeanLoad(L1 = coef(ModelPinworm$H1$fitAlpha)[names(coef(ModelPinworm$H1$fitAlpha)) == "L1"],
                                           L2 =  coef(ModelPinworm$H1$fitAlpha)[names(coef(ModelPinworm$H1$fitAlpha)) == "L2"],
                                           alpha =  alphaCILB,
                                           hybridIndex = seq(0,1,0.01)),
                 loadMLEAlphaUB = MeanLoad(L1 = coef(ModelPinworm$H1$fitAlpha)[names(coef(ModelPinworm$H1$fitAlpha)) == "L1"],
                                           L2 =  coef(ModelPinworm$H1$fitAlpha)[names(coef(ModelPinworm$H1$fitAlpha)) == "L2"],
                                           alpha =  alphaCIUB,
                                           hybridIndex = seq(0,1,0.01)))

ggplot() +
  geom_point(data = Joelle_data, aes(x = HI, y = log10(Aspiculuris.Syphacia + 1), color = Sex)) +
  geom_ribbon(aes(x = DF$HI, 
                  ymin = log10(DF$loadMLEAlphaUB + 1), 
                  ymax = log10(DF$loadMLEAlphaLB + 1)),
              fill = "grey", alpha = .5) +
  geom_line(aes(x = DF$HI, y = log10(DF$loadMLE + 1))) +
  geom_line(aes(x = DF$HI, y = log10(DF$loadMLEnoAlpha + 1)), linetype="dotted") +
  scale_color_manual(values = c("red", "darkblue")) +
  theme_linedraw()

