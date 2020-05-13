rm(list = ls(all = TRUE))
## load packages ####
if(!require(saezero)) devtools::intall_github("XiaodanLyu/saezero")
library(saezero)
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
for (rho in rhos){
  system.time(
    sim_mc(rho = rho, nsim = 1000, true.para = true.para, link = "logit", 
           alts = FALSE, bootstrap = FALSE, seed = 2020)
  )
}

## parametric bootstrap ####
## each takes about 6~7.5 hours
for (seed in 2015:2019){
  system.time(
    sim_mc(seed = seed, nsim = 200, true.para = true.para,
           rho = 0.9, link = "logit", alts = FALSE, bootstrap = TRUE)
  )
}

## tables and figures ####
setwd("simulation")
## ---- tb1
load("intermediate_results/Cor0.9_D60_nsim1000_seed2020_Result.RData")
mcmse.eb <- tapply(colMeans((pred.eb.store-YbarNis.store)^2), nis, mean)
## MC MSE of the proposed EB predictor
sprintf("%.2f", mcmse.eb*10^5)
## MSE differences between the alternative predictors and the EB predictor
dmse.alt <- function(pred.store){
  mcmse <- apply((pred.store-YbarNis.store)^2, 1, tapply, nis, mean)
  dmse <- mcmse-as.vector(mcmse.eb)
  diff <- rowMeans(dmse)*1e+5
  err <- 1.96*sqrt(apply(dmse, 1, var)/1000)*1e+5
  paste0(sprintf("%.2f", diff), " (", sprintf("%.2f", err), ")")
}
dmse.result <- do.call(
  "cbind",
  lapply(list(eb0 = pred.eb0.store,
              pi = pred.pi.store,
              zi = pred.zi.store,
              si = pred.si.store),
         dmse.alt)) 
dmse.result

## ---- tb2
## MSE differences between the EB(0) predictor and the EB predictor
dmse.store <- c()
rhos <- c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
for (rho in rhos){
  load(sprintf("intermediate_results/Cor%s_D60_nsim1000_seed2020_Result.RData", rho))
  mcmse.eb <- apply((pred.eb.store-YbarNis.store)^2, 1, tapply, nis, mean)
  mcmse.eb0 <- apply((pred.eb0.store-YbarNis.store)^2, 1, tapply, nis, mean)
  dmse <- mcmse.eb0-mcmse.eb
  diff <- rowMeans(dmse)*1e+5
  err <- 1.96*sqrt(apply(dmse, 1, var)/1000)*1e+5
  dmse.store <- cbind(dmse.store, paste0(sprintf("%.2f", diff), " (", sprintf("%.2f", err), ")"))
}
colnames(dmse.store) <- rhos
dmse.store

## ---- tb3
## assemble simulation results
seeds <- 2015:2019
m12boot <- m1biasboot <- m2boot <- mseboot <- c()
ybar <- predeb <- onestep <- c()
for (seed in seeds){
  load(sprintf("intermediate_results/bootstrap/Cor0.9_D60_nsim200_seed%.0f_Result.RData", seed))
  m12boot <- rbind(m12boot, M12boot.store)
  m1biasboot <- rbind(m1biasboot, M1biasboot.store)
  m2boot <- rbind(m2boot, M2boot.store)
  mseboot <- rbind(mseboot, MSEboot.store)
  ybar <- rbind(ybar, YbarNis.store)
  predeb <- rbind(predeb, pred.eb.store)
  onestep <- rbind(onestep, mse.eb.store)
}
## bias and coverage for the MSE estimators
mcmse <- colMeans((predeb - ybar)^2)
evalmse <- function(estmse){
  bmse <- apply(estmse, 1, function(x) tapply(x-mcmse, nis, mean))
  cp <- apply(abs(predeb - ybar) <= 1.96*sqrt(estmse), 1, tapply, nis, mean)
  bmse.val <- rowMeans(bmse)*1e+5
  bmse.err <- 1.96*sqrt(apply(bmse, 1, var)/1000)*1e+5
  cp.val <- rowMeans(cp)*100
  cp.err <- 1.96*sqrt(apply(cp, 1, var)/1000)*100
  data.frame(
    bmse = paste0(sprintf("%.2f", bmse.val), " (", sprintf("%.2f", bmse.err), ")"),
    cp = paste0(sprintf("%.2f", cp.val), " (", sprintf("%.2f", cp.err), ")"))
}
output <- do.call("cbind", lapply(list(
  onestep = onestep, boot = mseboot, semiboot = onestep + m2boot),
  evalmse)) 
output

## ---- tbs2
## composition of the three MSE components
bootres <- do.call("cbind", lapply(
  list(m2boot/mseboot, m1biasboot/mseboot, m12boot/mseboot), colMeans))
bootres <- apply(bootres, 2, tapply, nis, mean)
round(bootres*100, 2)
