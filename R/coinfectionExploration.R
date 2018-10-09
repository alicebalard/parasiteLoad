## Installation
library(parasiteLoad)
library(bbmle)
require(optimx) # for bbmle it needs to be required(?)
library(ggplot2)
library(MASS)


## Prepare dataset for each analysis

# add 1 for worms/oocysts count!!
WATWMdata$`Aspiculuris.Syphacia+1` <- WATWMdata$Aspiculuris.Syphacia + 1
BALdata$`Aspiculuris.Syphacia+1` <- BALdata$Aspiculuris_Syphacia + 1
BALdata$`OPG+1` <- BALdata$OPG + 1

# and select variables for each analysis
qpcr_data <- BALdata[!is.na(BALdata$HI) &
                       !is.na(BALdata$Sex) &
                       !is.na(BALdata$delta_ct_MminusEAllPos), ]

flotation_data <- BALdata[!is.na(BALdata$`OPG+1`) &
                            !is.na(BALdata$HI) &
                            !is.na(BALdata$Sex), ]


pinworms_data <- BALdata[!is.na(BALdata$`Aspiculuris.Syphacia+1`) &
                           !is.na(BALdata$HI) &
                           !is.na(BALdata$Sex),]

body_data <- BALdata[!is.na(BALdata$Body_weight) &
                       !is.na(BALdata$Body_length) &
                       !is.na(BALdata$HI) &
                       !is.na(BALdata$Sex) &
                       !is.na(BALdata$qPCRstatus) ,]

WATWMdata <- WATWMdata[!is.na(WATWMdata$`Aspiculuris.Syphacia+1`) &
                         !is.na(WATWMdata$HI) &
                         !is.na(WATWMdata$Sex), ]

# all worms together for comparison
d1 <- WATWMdata[c("Sex", "Aspiculuris.Syphacia+1", "HI")]
d2 <- pinworms_data[c("Sex", "Aspiculuris.Syphacia+1", "HI")]
d1$batch <- "WATWM"
d2$batch <- "JENNY"
allWorms <- rbind(d1,d2)
allWorms$batch <- as.factor(allWorms$batch)


## Test coinfection
BALdata$coinfStatus[!is.na(BALdata$Aspiculuris_Syphacia) &
                      BALdata$delta_ct_MminusEAllPos > 0]

                      BALdata$Aspiculuris_Syphacia > 0 &
                      !is.na(BALdata$delta_ct_MminusEAllPos) &
                      BALdata$delta_ct_MminusEAllPos > 0]

ggplot(BALdata, aes(x = delta_ct_MminusEAllPos, y = Aspiculuris_Syphacia)) +
  geom_point(pch = 21, size = 3, aes(fill = HI)) +
  scale_fill_continuous(high = "red", low = "blue") +
  theme_classic()

