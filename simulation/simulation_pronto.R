rm(list = ls(all = TRUE))
## load packages ####
library(saezero, lib.loc = "/home/lyux/R/x86_64-pc-linux-gnu-library/")
library(dplyr)
library(lme4)
source("simulation/utility-functions.R")

## parameter configuration ####
true.para <- list(beta = c(-13, 2), sig2lu = 0.22, sig2le = 1.23,
                  alpha = c(-20, 5), sig2lb = 0.52)

## compare with alternative predictors ####
## takes about 3~3.5 hours
system.time(
  sim_mc(nsim = 1000, true.para = true.para, rho = 0.9, link = "logit",
         alts = TRUE, bootstrap = FALSE, seed = 2020)
)

## compare with EB(0) predictor as rho changes ####
rhos <- c(-0.9, -0.6, -0.3, 0, 0.3, 0.6)
## each takes about 3~3.5 hours
system.time(
  mclapply(rhos, sim_mc, nsim = 1000, true.para = true.para,
           link = "logit", alts = FALSE, bootstrap = FALSE,
           seed = 2020, mc.cores = 6)
)