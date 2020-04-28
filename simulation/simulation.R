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
mcmse.eb <- colMeans((pred.eb.store-YbarNis.store)^2)
dmse.alt <- function(pred.store){
  mcmse <- apply((pred.store-YbarNis.store)^2, 1, tapply, nis, mean)
  dmse <- mcmse-mcmse.eb
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
dmse.store <- c()
rhos <- c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
for (rho in rhos){
  load(sprintf("intermediate_results/Cor%s_D60_nsim1000_seed2020_Result.RData", rho))
  mcmse.eb <- apply((pred.eb.store-YbarNis.store)^2, 1, tapply, nis, mean)
  mcmse.eb0 <- apply((pred.eb0.store-YbarNis.store)^2, 1, tapply, nis, mean)
  dmse <- mcmse.eb0-mcmse.eb
  diff <- rowMeans(dmse)*1e+5
  err <- 1.96*sqrt(apply(dmse, 1, var)/1000)*1e+5
  # rej <- ifelse(diff>err, TRUE, FALSE)
  dmse.store <- cbind(dmse.store, paste0(sprintf("%.2f", diff), " (", sprintf("%.2f", err), ")"))
}
colnames(dmse.store) <- rhos
dmse.store

## ---- tb3
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
bootres <- do.call("cbind", lapply(
  list(m2boot/mseboot, m1biasboot/mseboot, m12boot/mseboot), colMeans))
bootres <- apply(bootres, 2, tapply, nis, mean)
round(bootres*100, 2)

## ---- figs2
library(dplyr)
library(ggplot2)
library(patchwork)
ebeb0.store <- c()
for (rho in rhos){
  load(sprintf("intermediate_results/Cor%s_D60_nsim1000_seed2020_Result.RData", rho))
  mcmse.eb <- colMeans((pred.eb.store-YbarNis.store)^2)
  mcmse.eb0 <- colMeans((pred.eb0.store-YbarNis.store)^2)
  rbmse <- colMeans(mse.eb.store)/mcmse.eb-1
  rbmse0 <- colMeans(mse.eb0.store)/mcmse.eb0-1
  cp <- colMeans(abs(pred.eb.store - YbarNis.store) <= 1.96*sqrt(mse.eb.store))
  cp0 <- colMeans(abs(pred.eb0.store - YbarNis.store) <= 1.96*sqrt(mse.eb0.store))
  ebeb0.store <- rbind(ebeb0.store, cbind(rho = rho, nis, rbmse, rbmse0, cp, cp0))
}
tb <- ebeb0.store %>% as_tibble() %>% 
  tidyr::gather(metric, value, rbmse:cp0) %>% 
  group_by(rho, nis, metric) %>% summarise(value = mean(value)) %>% 
  ungroup() %>% mutate(nis = forcats::fct_relevel(factor(nis), "5")) 
g1 <- tb %>% filter(grepl("rbmse", metric)) %>%
  rename(RBMSE = value) %>% 
  ggplot(aes(x = rho, y = RBMSE, color = metric)) + 
  geom_point(aes(shape = metric), size = rel(2)) + geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(labels = function(x) sprintf("%.2f", x)) +
  scale_x_continuous(breaks = rhos) + 
  guides(color = FALSE, shape = FALSE) +
  theme_bw(base_size = 15) +
  facet_grid(~nis)
g2 <- tb %>% filter(grepl("cp", metric)) %>% 
  rename(CP = value) %>% 
  mutate(metric = recode_factor(metric, "cp" = "EB", "cp0" = "EB0")) %>% 
  ggplot(aes(x = rho, y = CP, color = metric)) + 
  geom_point(aes(shape = metric), size = rel(2)) + geom_line() +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  scale_x_continuous(breaks = rhos) + 
  theme_bw(base_size = 15) + labs(color = "", shape = "") + 
  facet_grid(~nis) + theme(legend.position = "bottom") 
g1 / g2
