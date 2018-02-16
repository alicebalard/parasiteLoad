source("../R/WithMLE2.R")

## Import data WATWM
NewWorms <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/raw_data/MiceTable_2014to2017.csv")
NewWorms$pinworms <- NewWorms$Aspiculuris_tetraptera + NewWorms$Syphacia_obvelata

# Correct mistakes
NewWorms$HI[!is.na(NewWorms$HI) & NewWorms$HI > 1] <- NewWorms$HI[!is.na(NewWorms$HI) & NewWorms$HI > 1]/1000

NewWorms$Sex[!is.na(NewWorms$Sex) & NewWorms$Sex == "male (young)"] <- "M"
NewWorms$Sex[!is.na(NewWorms$Sex) & NewWorms$Sex == "female (pregnant)"] <- "F"
NewWorms$Sex[!is.na(NewWorms$Sex) & NewWorms$Sex == "young"] <- NA
NewWorms$Sex[!is.na(NewWorms$Sex) & NewWorms$Sex == "female (young)"] <- "F"
NewWorms$Sex <- droplevels(NewWorms$Sex)

# remove unusable lines
NewWorms <- NewWorms[complete.cases(NewWorms[c("HI", "Sex", "pinworms")]),]

################## Data analysis ################## 
system.time(ModelNewAspi <- myFun(data = NewWorms, 
                                  hybridIndex = HI, 
                                  response = "pinworms"))

finalFun(ModelNewAspi)

ggplot() +
  geom_point(data = aspi, aes(x = HI, y = log10(pinworms + 1), color = Sex)) 
################## Plotting ################## 

## profile investigates behavior of objective function near the MLE
system.time(myProf <- profile(ModelNewAspi$H3$fitAlpha))

## Marginal confidence interval
myConfInt <- confint(myProf)
myConfInt

## Marginal confidence interval for alpha
alphaCILB_M <- myConfInt[rownames(myConfInt) == "alpha.(Intercept)"][1] +
  myConfInt[rownames(myConfInt) == "alpha.SexM"][1]
alphaCILB_F <- myConfInt[rownames(myConfInt) == "alpha.(Intercept)"][1] 
alphaCIUB_M <- myConfInt[rownames(myConfInt) == "alpha.(Intercept)"][2] +
  myConfInt[rownames(myConfInt) == "alpha.SexM"][2]
alphaCIUB_F <- myConfInt[rownames(myConfInt) == "alpha.(Intercept)"][2] 

## Draw the line for the parameters at their MLE, alpha varying
DF <- data.frame(HI = seq(0,1,0.01),
                 loadMLE_M = MeanLoad(L1 = coef(ModelPinworm$H3$fitAlpha)[names(coef(ModelPinworm$H3$fitAlpha)) == "L1.(Intercept)"] +
                                        coef(ModelPinworm$H3$fitAlpha)[names(coef(ModelPinworm$H3$fitAlpha)) == "L1.SexM"],
                                      L2 = coef(ModelPinworm$H3$fitAlpha)[names(coef(ModelPinworm$H3$fitAlpha)) == "L2.(Intercept)"] +
                                        coef(ModelPinworm$H3$fitAlpha)[names(coef(ModelPinworm$H3$fitAlpha)) == "L2.SexM"],
                                      alpha =   coef(ModelPinworm$H3$fitAlpha)[names(coef(ModelPinworm$H3$fitAlpha)) == "alpha.(Intercept)"] +
                                        coef(ModelPinworm$H3$fitAlpha)[names(coef(ModelPinworm$H3$fitAlpha)) == "alpha.SexM"],
                                      hybridIndex = seq(0,1,0.01)),
                 loadMLE_F = MeanLoad(L1 = coef(ModelPinworm$H3$fitAlpha)[names(coef(ModelPinworm$H3$fitAlpha)) == "L1.(Intercept)"],
                                      L2 = coef(ModelPinworm$H3$fitAlpha)[names(coef(ModelPinworm$H3$fitAlpha)) == "L2.(Intercept)"], 
                                      alpha =   coef(ModelPinworm$H3$fitAlpha)[names(coef(ModelPinworm$H3$fitAlpha)) == "alpha.(Intercept)"],
                                      hybridIndex = seq(0,1,0.01)),
                 loadMLE_M_noalpha = MeanLoad(L1 = coef(ModelPinworm$H3$fitNoAlpha)[names(coef(ModelPinworm$H3$fitNoAlpha)) == "L1.(Intercept)"] +
                                                coef(ModelPinworm$H3$fitNoAlpha)[names(coef(ModelPinworm$H3$fitNoAlpha)) == "L1.SexM"],
                                              L2 = coef(ModelPinworm$H3$fitNoAlpha)[names(coef(ModelPinworm$H3$fitNoAlpha)) == "L2.(Intercept)"] +
                                                coef(ModelPinworm$H3$fitNoAlpha)[names(coef(ModelPinworm$H3$fitNoAlpha)) == "L2.SexM"],
                                              alpha = 0,
                                              hybridIndex = seq(0,1,0.01)),
                 loadMLE_F_noalpha = MeanLoad(L1 = coef(ModelPinworm$H3$fitNoAlpha)[names(coef(ModelPinworm$H3$fitNoAlpha)) == "L1.(Intercept)"],
                                              L2 = coef(ModelPinworm$H3$fitNoAlpha)[names(coef(ModelPinworm$H3$fitNoAlpha)) == "L2.(Intercept)"], 
                                              alpha = 0,
                                              hybridIndex = seq(0,1,0.01)),
                 loadMLEAlphaLB_M = MeanLoad(L1 = coef(ModelPinworm$H1$fitAlpha)[names(coef(ModelPinworm$H1$fitAlpha)) == "L1"],
                                             L2 =  coef(ModelPinworm$H1$fitAlpha)[names(coef(ModelPinworm$H1$fitAlpha)) == "L2"],
                                             alpha =  alphaCIUB_M,
                                             hybridIndex = seq(0,1,0.01)),
                 loadMLEAlphaLB_F = MeanLoad(L1 = coef(ModelPinworm$H1$fitAlpha)[names(coef(ModelPinworm$H1$fitAlpha)) == "L1"],
                                             L2 =  coef(ModelPinworm$H1$fitAlpha)[names(coef(ModelPinworm$H1$fitAlpha)) == "L2"],
                                             alpha =  alphaCIUB_F,
                                             hybridIndex = seq(0,1,0.01)),
                 loadMLEAlphaUB_M = MeanLoad(L1 = coef(ModelPinworm$H1$fitAlpha)[names(coef(ModelPinworm$H1$fitAlpha)) == "L1"],
                                             L2 =  coef(ModelPinworm$H1$fitAlpha)[names(coef(ModelPinworm$H1$fitAlpha)) == "L2"],
                                             alpha =  alphaCILB_M,
                                             hybridIndex = seq(0,1,0.01)),
                 loadMLEAlphaUB_F = MeanLoad(L1 = coef(ModelPinworm$H1$fitAlpha)[names(coef(ModelPinworm$H1$fitAlpha)) == "L1"],
                                             L2 =  coef(ModelPinworm$H1$fitAlpha)[names(coef(ModelPinworm$H1$fitAlpha)) == "L2"],
                                             alpha =  alphaCILB_F,
                                             hybridIndex = seq(0,1,0.01)))
ggplot() +
  geom_point(data = NewWorms, aes(x = HI, y = log10(pinworms + 1), color = Sex)) +
  geom_ribbon(aes(x = DF$HI, 
                  ymin = log10(DF$loadMLEAlphaLB_M + 1),
                  ymax = log10(DF$loadMLEAlphaUB_M + 1)), fill = "blue", alpha = .5) +
  geom_ribbon(aes(x = DF$HI, 
                  ymin = log10(DF$loadMLEAlphaLB_F + 1),
                  ymax = log10(DF$loadMLEAlphaUB_F + 1)), fill = "red", alpha = .5) +
  geom_line(aes(x = DF$HI, y = log10(DF$loadMLE_M + 1)), col = "blue") +
  geom_line(aes(x = DF$HI, y = log10(DF$loadMLE_F + 1)), col = "red") +
  geom_line(aes(x = DF$HI, y = log10(DF$loadMLE_M_noalpha + 1)), col = "blue", linetype="dotted") +
  geom_line(aes(x = DF$HI, y = log10(DF$loadMLE_F_noalpha + 1)), col = "red", linetype="dotted") +
  theme_linedraw()

