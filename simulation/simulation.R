rm(list = ls(all = TRUE))
## load packages
if(!require(saezero)) devtools::intall_github("XiaodanLyu/saezero")
library(saezero)
library(dplyr)
library(lme4)
source("simulation/utility-functions.R")

## parameter configuration ####
true.para <- list(beta = c(-13, 2), sig2lu = 0.22, sig2le = 1.23,
                  alpha = c(-20, 5), sig2lb = 0.52)

# set.seed(2020)
# system.time(
#   sim_mc(nsim = 3, true.para, rho = 0.9, link = "logit",
#          alts = TRUE, bootstrap = TRUE)
# )

library(parallel)
set.seed(200405)
system.time(
  mclapply(c(-0.9, -0.6, -0.3, 0, 0.3, 0.6), sim_mc,
           nsim = 1000, true.para = true.para, link = "logit",
           alts = FALSE, bootstrap = FALSE, mc.cores = 6)
)
