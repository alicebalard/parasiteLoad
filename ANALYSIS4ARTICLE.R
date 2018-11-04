# Data analysis : test of hybrid immune vigor in response to parasite infection
# author: Alice Balard

## Installation
library(parasiteLoad)
library(bbmle)
require(optimx) # for bbmle it needs to be required(?)
library(ggplot2)
library("VennDiagram")
# library(MASS)
library(fitdistrplus) # evaluate distribution
# install.packages("epiR")
library(epiR) # Sterne's exact method
# install.packages("simpleboot")
library(simpleboot) # BS
library(boot) # BS
library(ggmap)
# devtools::install_github("3wen/legendMap")
library(legendMap)
library(gridExtra) # several plots in one panel

##### Functions used for data analysis
plotMap <- function(df, margin = 1, zoom = 8,
                    source = "stamen", maptype = "toner-lite"){
  area <- get_map(location =
                    c(min(df$Longitude - margin),
                      min(df$Latitude - margin),
                      max(df$Longitude + margin),
                      max(df$Latitude +  margin)),
                  source = source, maptype= maptype,
                  zoom = zoom)
  ggmap(area) +
    geom_point(data = df, shape = 21, size = 2,
               aes(Longitude, Latitude, fill = HI), alpha = .4) + # set up the points
    scale_fill_gradient("Hybrid\nindex", high="red",low="blue") +
    geom_rect(xmin = 11, xmax = 12.7, ymin = 50, ymax = 51.3, fill = "white") +
    scale_bar(lon = 11.8, lat = 50.9, arrow_length = 10, arrow_distance = 20,
              distance_lon = 20, distance_lat = 5, distance_legend = 10,
              dist_unit = "km", orientation = TRUE, legend_size = 2, arrow_north_size = 3)
}

##### Prepare datasets for each analysis
WATWMdata <- read.csv("data/WATWMdata.csv")
BALdata <- read.csv("data/BALdata.csv")

# Keep individuals with hybrid index and sex
WATWMdata <- WATWMdata[!is.na(WATWMdata$HI) & !is.na(WATWMdata$Sex),]
BALdata <- BALdata[!is.na(BALdata$HI) & !is.na(BALdata$Sex),]

# pinworms balard
pinwormsdata_bal <- BALdata[!is.na(BALdata$Aspiculuris_Syphacia),]
pinwormsdata_bal$`Aspiculuris.Syphacia+1` <-
  pinwormsdata_bal$Aspiculuris_Syphacia + 1

# pinworms "where are the wormy mice"
pinwormsdata_watwm <- WATWMdata[!is.na(WATWMdata$Aspiculuris.Syphacia),]
pinwormsdata_watwm$`Aspiculuris.Syphacia+1` <-
  pinwormsdata_watwm$Aspiculuris.Syphacia + 1

# Eimeria qpcr balard
qpcrdata <- BALdata[!is.na(BALdata$delta_ct_cewe_MminusE) |
                     !is.na(BALdata$delta_ct_ilwe_MminusE),]
df <- qpcrdata[, c("delta_ct_cewe_MminusE", "delta_ct_ilwe_MminusE")]
qpcrdata$delta_ct_max_MminusE <- apply(df, 1, function(x){max(x, na.rm = T)})
rm(df)
# threshold of detection by qPCR = -5. Then we add -5 to all to have null or positive values
qpcrdata$delta_ct_max_MminusE[qpcrdata$delta_ct_max_MminusE <= -5] <- -5
qpcrdata$delta_ct_max_MminusEPOS <- qpcrdata$delta_ct_max_MminusE + 5

# Eimeria oocysts flotation balard
flotationdata <- BALdata[!is.na(BALdata$OPG),]
flotationdata$OPG <- round(flotationdata$OPG)
flotationdata$`OPG+1` <- flotationdata$OPG + 1

# Body condition index data TBC...

## A. General parasitology

myQuantitativeParasitology <- function(x){

  intensity <- round(mean(x[x>0]),1)
  abundance <- round(mean(x), 1)
  max <- max(x)
  Ni <- length(x)
  NiPos <- length(x[x>0])

  # Confidence intervals for prevalence calculated with Sterne's exact method
  sternetest <- epiR::epi.prev(pos = length(x[x > 0]), tested = length(x),
                               se = 1, sp=1, conf.level = .95, method = "sterne")
  cilow <- sternetest$ap["lower"]
  cihigh <- sternetest$ap["upper"]
  prevalence <- sternetest$ap["est"]

  # "Confidence intervals based on normal theory perform poorly
  # for skewed distributions, in particular if the sample is small,
  # but BCa (bias-corrected and accelerated) bootstrap confidence
  # intervals proposed by Efron and Tibshirani (1993) offer a solution
  # even in such cases"http://zoologia.hu/list/quant_large.pdf

  set.seed(1234) # to make it replicable

  bA <- one.boot(x, mean, R=1000)
  myBSA <- boot.ci(bA, type=c("bca"))
  confinfA <- round(myBSA$bca[4],1)
  confsupA <- round(myBSA$bca[5],1)

  bI <- one.boot(x[x>0], mean, R=1000)
  myBSI <- boot.ci(bI, type=c("bca"))
  confinfI <- round(myBSI$bca[4],1)
  confsupI <- round(myBSI$bca[5],1)

  ## Printout results
  Result <- cat(paste0("Prevalence % [CI 95%] (N infected hosts/ N hosts)\n",
                       round(prevalence,1), " [", round(cilow,1), "-", round(cihigh,1), "]",
                       " (", NiPos, "/", Ni,  ")\n",
                       "Abundance [CI 95%] (Max parasite load)\n",
                       round(abundance,1), " [", confinfA, "-", confsupA, "]",
                       " (", max, ")\n",
                       "Intensity [CI 95%] (Max parasite load)\n",
                       round(intensity,1), " [", confinfI, "-", confsupI, "]",
                       " (", max, ")"))

  return(Result)
}

myQuantitativeParasitology(pinwormsdata_bal$Aspiculuris_Syphacia)
myQuantitativeParasitology(qpcrdata$delta_ct_max_MminusEPOS)
myQuantitativeParasitology(flotationdata$OPG)

listWorms <- c("Hymenolepis", "Taenia", "Aspiculuris_Syphacia", "Trichuris",
               "Heterakis", "Mastophorus")
apply(BALdata[,listWorms], 2, table)

## B. Maps

# Plot all samples : pinworms
# pdf(file = "../figures/mapAllworms_2014to2017", width = 7, height = 6)
plotMap(df = pinwormsdata_bal)
# dev.off()

# Plot infected samples : pinworms
# pdf(file = "../figures/mapAllworms_2014to2017", width = 7, height = 6)
plotMap(df = pinwormsdata_bal[pinwormsdata_bal$Aspiculuris_Syphacia >=1,])
# dev.off()

# Plot all samples : eimeria (tissues)
# pdf(file = "../figures/mapAllworms_2014to2017", width = 7, height = 6)
plotMap(df = qpcrdata)
# dev.off()

# Plot infected samples : eimeria (tissues)
# pdf(file = "../figures/mapAllworms_2014to2017", width = 7, height = 6)
plotMap(df = qpcrdata[qpcrdata$delta_ct_max_MminusEPOS > 0,])
# dev.off()

## C. Choose distributions for each dataset

# macro parasite counts -> negbin (ref)
# qpcr -> explore best distribution

## MLE of parameter k negbin
hist(x, prob=T, breaks = max(x))
# fit the negative binomial distribution
fit <- fitdist(x, "nbinom")
# get the fitted densities. mu and size from fit.
fitD <- dnbinom(0:max(x), size=fit$estimate[1], mu=fit$estimate[2])
# add fitted line (blue) to histogram
lines(fitD, lwd="3", col="blue")
# Goodness of fit with the chi squared test
# get the frequency table
t <- table(x)
# convert to dataframe
df <- as.data.frame(t)
df$x <- as.numeric(as.character(df$x))
df$expp <- pnbinom(q = df$x, size=fit$estimate[1], mu = fit$estimate[2])
df$obsp <- cumsum(df$Freq / sum(df$Freq))
# perform the chi-squared test
chisq.test(x = df$expp, y = df$obsp)

# Pearson's Chi-squared test
#
# data:  df$expp and df$obsp
# X-squared = 9120, df = 9025, p-value = 0.239

# Observed and expected frequencies do not differ significantly (P =
# 0.05) so there is no statistical evidence to reject H 0

size=fit$estimate[1]
size
# The negative binomial distribution seems to describe parasite load well for all parasites

### qPCR balard

### flotation balard

## D. Compare old transect / new transect

## E. Apply & plot our model

## F. Supplementary

# Choice of qPCR
adjrsq <- summary(lm(qpcrdata$delta_ct_max_MminusE ~ qpcrdata$OPG))$adj.r.squared
F1 <- ggplot(qpcrdata, aes(delta_ct_max_MminusE, OPG +1)) +
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = 1, ymax = min(qpcrdata$OPG[qpcrdata$OPG >0], na.rm = T) +1,
           alpha = 0.2) +
  geom_point() +
  scale_y_log10() +
  stat_smooth(method = "lm", fill = "pink", col = "red") +
  xlab("Eimeria detection in tissues\n(Eimeria-host DNA log-ratio)") +
  ylab("Eimeria detection in feces\n(Eimeria oocysts per gram of feces)") +
  annotate("text", x = 1, y = 1300, label = paste0("R = ", round(adjrsq, 2))) +
  geom_density(data = qpcrdata[!is.na(qpcrdata$OPG) & qpcrdata$OPG == 0,],
               aes(x = delta_ct_max_MminusE)) +
  theme_classic()
F1

ggplot() +
geom_density(data = qpcrdata[!is.na(qpcrdata$OPG) & qpcrdata$OPG == 0,],
             aes(x = delta_ct_max_MminusE, stat(count)))

## bad correlation + grey zone for OPG + not same samples tested

## + Venn diagram OPG/qPCR

area1 = nrow(subset(qpcrdata, delta_ct_max_MminusE > -5))
area2 = nrow(subset(qpcrdata, !is.na(OPG) & OPG > 0))
## areas of 2-group overlap
cross.area = nrow(subset(qpcrdata, delta_ct_max_MminusE > -5 & !is.na(OPG) & OPG > 0))
grid.newpage()
F2 <- draw.pairwise.venn(area1 = area1,
                     area2 = area2,
                     cross.area = cross.area,
                     category = c("Eimeria detection\nin tissues", "Eimeria detection\nin feces"),
                     col = "transparent",
                     fill = c("purple","yellow"),
                     alpha = 0.40,
                     cex = 1.5, cat.cex = 1.3, fontfamily = "serif", fontface = "bold",
                     cat.pos = c(-20,33), cat.dist = c(-.09, -.1),
                     cat.fontfamily = "serif")

grid.arrange(gTree(children=F2), F1, ncol = 2)


### pinworms (A. tetraptera and S. obvelata (WATWMdata$Aspiculuris.Syphacia))

```{r runfitJo, fig.width=7, fig.height=4}
fit_Joelle_negbin <- analyse(WATWMdata, "Aspiculuris.Syphacia+1",
                             model = "negbin", group = "Sex")
fit_Joelle_negbin

plot_Joelle_negbin <- bananaPlots(mod = fit_Joelle_negbin$H3,
                                  data = WATWMdata,
                                  response = "Aspiculuris.Syphacia+1",
                                  islog10 = TRUE, group = "Sex")
plot_Joelle_negbin
```

### And Jennys now:

```{r runfitJe,  fig.width=7, fig.height=4}
fit_pinworms_negbin <- analyse(pinworms_data, response = "Aspiculuris.Syphacia+1",
                               model = "negbin", group = "Sex")
fit_pinworms_negbin

plot_pinworms_negbin <- bananaPlots(mod = fit_pinworms_negbin$H3,
                                    data = pinworms_data,
                                    response = "Aspiculuris.Syphacia+1",
                                    islog10 = TRUE, group = "Sex")
plot_pinworms_negbin

pdf(file = "../figures/part1.pdf", width=7, height=4)
plot_pinworms_negbin
dev.off()
```

To compare both dataset hybrid effect (alpha) we run 2 models with our data:

  * model A with fixed alpha = 1.44 (males) alpha = 1.34 (females) (Bairds et al. 2012 value for pinworms)

* model B with variable alpha (model B being fit_pinworms_negbin)

then we compare via G-test.

```{r}
defaultConfig <- list(optimizer = "optimx",
                      method = c("L-BFGS-B", "bobyqa"),
                      control = list(follow.on = TRUE))

## Female: alpha fixed at 1.34
model = "negbin"
data = pinworms_data[pinworms_data$Sex =="F",]
response = "Aspiculuris.Syphacia+1"
paramBounds <- getParamBounds(model, data, response)
paramBounds[["alphaStart"]] <- 1.34
paramBounds[["alphaUB"]] <- 1.35
paramBounds[["alphaLB"]] <- 1.33

fitF <- FitAdvancedAlphaNegbin(data, response, hybridIndex = "HI",
                               paramBounds, config = defaultConfig)

## Male: alpha fixed at 1.44
model = "negbin"
data = pinworms_data[pinworms_data$Sex =="M",]
response = "Aspiculuris.Syphacia+1"
paramBounds <- getParamBounds(model, data, response)
paramBounds[["alphaStart"]] <- 1.44
paramBounds[["alphaUB"]] <- 1.45
paramBounds[["alphaLB"]] <- 1.43

fitM <- FitAdvancedAlphaNegbin(data, response, hybridIndex = "HI",
                               paramBounds, config = defaultConfig)

# Gtest
LL1 <- logLik(fitF) + logLik(fitF)
LL2 <- logLik(fit_pinworms_negbin$H3$groupA) +
  logLik(fit_pinworms_negbin$H3$groupB)
dLL <- LL1 - LL2
dDF <- 2 # 2 alpha fixed in one case
1 - pchisq(2*0.66, df=2)
```
