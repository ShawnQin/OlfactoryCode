#*************************************************************
# this program analyze data set from Saito et al, 2009 
# test if the binidng affinity follows log normal distribution
# noticed that in that paper EC50 is used as sensitivity
# they tested 219 mouse receptor and 245 human receptor with 93 odorant
# and identified only 10+52 agnoist odorants
# the date set contian a matrix with 62x63 elements, with most of the 
# elements 0
#
# last revised on 5/19/2017
# ************************************************************

#load essential package
library(xlsx)
library(RColorBrewer) #this package is used for color settings
library(ggplot2)
library(fitdistrplus) #statistical test

# source self written scripts for distribution fitting
source("testDistriOSNOr.R")   #plot all single OSN or Odor
source("testDistAll.R")     # plot overall response distribtuion

#load the data
dataFile <- "/Users/shan/Documents/GoogleDrive/olfactoryCoding/data/Saito2009s3.xls"
rawData<-read.xlsx(dataFile,1)
Data <- rawData[,-1]       # the first column is non-number
allResp <- Data[Data !=0]  #all non-zero elements
linearResp <- exp(allResp) #linearize

#statistical property of these sensitivity
hist(allResp,20)

#test for lognormal distribution

