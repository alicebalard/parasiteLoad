## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>"
)

## ----install-------------------------------------------------------------
library(parasiteLoad)
library(bbmle)
require(optimx) # for bbmle it needs to be required(?)
library(ggplot2)
library(MASS)
library(fitdistrplus) # evaluate distribution
library(epiR) # Sterne's exact method
library(simpleboot) # BS
library(boot) # BS

