## ----setup, echo = FALSE, results = "hide"------------------------------------
options(signif = 3, digits = 3)
knitr::opts_chunk$set(tidy = FALSE, cache = TRUE, autodep = TRUE, fig.height = 5.5,
                      message = FALSE, error = FALSE, warning = TRUE)
set.seed(0xdada)

## -----------------------------------------------------------------------------
library(CAGEr)
data(FANTOM5humanSamples)
head(FANTOM5humanSamples)
nrow(FANTOM5humanSamples)

## -----------------------------------------------------------------------------
astrocyteSamples <- FANTOM5humanSamples[grep("Astrocyte", 
				FANTOM5humanSamples[,"description"]),]
astrocyteSamples

## -----------------------------------------------------------------------------
data(FANTOM5mouseSamples)
head(FANTOM5mouseSamples)
nrow(FANTOM5mouseSamples)

## -----------------------------------------------------------------------------
astrocyteSamples[,"sample"]

## -----------------------------------------------------------------------------
library(FANTOM3and4CAGE)
data(FANTOMhumanSamples)
head(FANTOMhumanSamples)
data(FANTOMmouseSamples)
head(FANTOMmouseSamples)

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

