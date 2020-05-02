rm(list = ls(all = TRUE))
## load packages ####
library(saezero, lib.loc = "/home/lyux/R/x86_64-pc-linux-gnu-library/")
library(dplyr)
library(lme4)
source("simulation/utility-functions.R")

args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
  seed <- 2016
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

## parameter configuration ####
true.para <- list(beta = c(-13, 2), sig2lu = 0.22, sig2le = 1.23,
                  alpha = c(-20, 5), sig2lb = 0.52)

## parametric bootstrap ####
## each takes about 6~7.5 hours
sim_mc(seed = seed, nsim = 200, true.para = true.para,
       rho = 0.9, link = "logit", alts = FALSE, bootstrap = TRUE)