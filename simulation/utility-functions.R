sim_mc <- function(nsim, true.para, rho = 0.9, link = "logit",
                   alts = FALSE, bootstrap = FALSE){

  N <- 10000
  D <- 60
  nis_uni <- c(5, 10, 20)
  nis <- rep(nis_uni, times = D/length(nis_uni))
  Nis <- round(N*nis/sum(nis), 0)
  N <- sum(Nis)
  
  mulx <- 4.45
  sig2lx <- 0.055
  lxN <- mulx + rnorm(N)*sqrt(sig2lx)
  areapop <- rep(1:length(Nis), times = Nis)
  Xs <- data.frame(x = lxN, area = areapop)
  
  fit <- list(fixed = list(p1 = true.para$beta, p0 = true.para$alpha),
              errorvar = true.para$sig2le, refvar1 = true.para$sig2lu,
              refvar0 = true.para$sig2lb, refcor = rho)
  
  iter <- 0
  ## true area means
  YbarNis.store <- c()
  ## EB and One-step MSE estimator
  pred.eb.store <- mse.eb.store <- c()
  ## EB(0) estimator
  pred.eb0.store <- mse.eb0.store <- c()
  ## alternate methods: plug-in, zero-ignored MMSE, shifted MMSE
  pred.pi.store <- pred.zi.store <- pred.si.store <- c()
  ## bootstrap terms
  M2boot.store <- MSEboot.store <- c()
  M1biasboot.store <- M12boot.store <- c()
  
  repeat{
    iter <- iter + 1
    yN <- simLBH(fit, Xs, f_pos = ~x, f_zero = ~x, f_area = ~area, link = link)
    
    YbarNis <- tapply(yN, areapop, mean)
    YbarNis.store <- rbind(YbarNis.store, YbarNis)
    
    popmc <- data.frame(area = areapop, y = yN, x = lxN)
    smc <- sapply(1:D, function(i) sample((1:N)[areapop==i], size = nis[i], replace = F))
    smc <- unlist(smc)
    srsmc <- popmc[smc,]
    Xoosmc <- popmc[-smc, -2]
    
    ## EB prediction
    srsmc_2p <- as.2pdata(f_pos = y~x, f_zero = ~x, f_area = ~area, data = srsmc)
    fitmc <- mleLBH(srsmc_2p, link)
    ebpred <- ebLBH(Xaux = Xoosmc, data_2p = srsmc_2p, fit = fitmc)
    pred.eb.store <- rbind(pred.eb.store, ebpred$eb)
    mse.eb.store <- rbind(mse.eb.store, ebpred$mse)
    
    ## EB(0) prediction
    fit0mc <- fitmc$fit0
    attr(fit0mc, "link") <- link
    eb0pred <- ebLBH(Xaux = Xoosmc, data_2p = srsmc_2p, fit = fit0mc)
    pred.eb0.store <- rbind(pred.eb0.store, eb0pred$eb)
    mse.eb0.store <- rbind(mse.eb0.store, eb0pred$mse)
    
    ## alternate methods ####
    if(alts){
      oosmc <- which(!(1:N)%in%smc)
      hij <- ebLogNormal(Xoosmc, f_pos = y~x, f_area = ~area, data = srsmc)
      ## PI estimator
      mzero <- glmer(delta~x+(1|area), family = binomial(link), nAGQ = 10,
                     data = srsmc %>% mutate(delta = y>0))
      phatpi <- predict(mzero, newdata = Xoosmc, type = "response")
      pipred <- popmc %>% mutate(y = replace(y, oosmc, hij*phatpi)) %>% 
        group_by(area) %>% summarise(y = mean(y)) %>% pull(y)
      pred.pi.store <- rbind(pred.pi.store, pipred)
      ## ZI estimator
      zipred <- popmc %>% mutate(y = replace(y, oosmc, hij)) %>% 
        filter(y>0) %>% group_by(area) %>% summarise(y = mean(y)) %>% pull(y)
      pred.zi.store <- rbind(pred.zi.store, zipred)
      ## SI estimator
      eps <- min(srsmc$y[srsmc$y>0], na.rm = TRUE)
      hijeps <- ebLogNormal(Xoosmc, f_pos = y~x, f_area = ~area,
                            data = srsmc %>% mutate(y = y + eps))
      sipred <- popmc %>% mutate(y = replace(y, oosmc, pmax(0, hijeps-eps))) %>% 
        group_by(area) %>% summarise(y = mean(y)) %>% pull(y)
      pred.si.store <- rbind(pred.si.store, sipred)
    }
    
    ## bootstrap terms ####
    if(bootstrap){
      boot.res <- pbmseLBH(Xs, srsmc_2p, smc, fitmc, link)
      M2boot.store <- rbind(M2boot.store, boot.res$EBM2)
      MSEboot.store <- rbind(MSEboot.store, boot.res$EBMSE)
      M1biasboot.store <- rbind(M1biasboot.store, boot.res$EBM1hat - ebpred$mse)
      M12boot.store <- rbind(M12boot.store, boot.res$EBM12)
    }
    
    print(paste("iter", iter))
    if(iter >= nsim){break}
  }
  
  f.name <- sprintf("simulation/intermediate_results/Cor%s_D%.0f_nsim%.0f_Result.RData",
                    rho, length(nis), nsim)
  save(list = c("nis", grep(".store", ls(all = T), value = TRUE)),
       file = f.name, envir = environment())
  
  return(NULL)
}

ebLogNormal <- function(Xaux, f_pos, f_area, data){
  
  area <- data[, all.vars(f_area)]
  area_oos <- Xaux[, all.vars(f_area)]
  nis <- tapply(area, area, length)
  mis <- tapply(area_oos, area_oos, length)
  Gs <- saezero:::Gmat(nis)
  Gns <- saezero:::Gmat(mis)
  
  alldata <- model.frame(f_pos, data)
  deltas <- ifelse(alldata[,1]>0, 1, 0)
  posdata <- alldata[which(deltas==1),]
  lys <- log(posdata[,1])
  Xs1  <- model.matrix(f_pos, posdata)
  
  areapos <- area[deltas==1]
  fit_pos <- lmer(lys~Xs1-1+(1|areapos))
  beta <- matrix(fixef(fit_pos), ncol=1)
  sig2le <- sigma(fit_pos)^2
  sig2lu <- unlist(VarCorr(fit_pos))
  
  nti <- tapply(deltas, area, sum)
  rt <- rep(0, length(area))
  rt[deltas==1] <- lys - Xs1%*%beta
  
  gamma <- sig2lu/(sig2lu+sig2le/nti)
  rbar <- drop(t(Gs)%*%rt)/nti
  rbar[nti==0] <- 0
  mu_u <- gamma*rbar
  mu_u[nti==0] <- 0
  sig2_u <- gamma*sig2le/nti
  sig2_u[nti==0] <- sig2lu
  Xs1_oos <- model.matrix(f_pos[-2], data = Xaux)
  yhat <- exp(Xs1_oos%*%beta + Gns%*%(mu_u+sig2_u/2) + sig2le/2)
  
  return(yhat)
}

pbmseLBH <- function(Xpop, sample_2p, smc, fit, link = "logit", B = 1){

  b <- 0
  pop_boot.store <- eb_boot.store <- mmse_boot.store <- m1_boot.store <- c()
  repeat{
    b <- b + 1
    ys <- simLBH(fit, Xpop, f_pos = ~x, f_zero = ~x, f_area = ~area)
    pop_boot <- tapply(ys, Xpop$area, mean)
    sample_boot <- Xpop %>% mutate(y = ys) %>% slice(smc)
    sample_boot_2p <- as.2pdata(f_pos = y~x, f_zero = ~x,
                                f_area = ~area, data = sample_boot)
    fit_boot <- mleLBH(sample_boot_2p, link = link)
    eb_boot <- ebLBH(Xpop[-smc], data_2p = sample_boot_2p, fit = fit_boot)$eb
    mmse_boot <- ebLBH(Xpop[-smc], data_2p = sample_boot_2p, fit = fit)$eb
    m1_boot <- ebLBH(Xpop[-smc], data_2p = sample_2p, fit = fit_boot)$mse
    pop_boot.store <- rbind(pop_boot.store, pop_boot)
    eb_boot.store <- rbind(eb_boot.store, eb_boot)
    mmse_boot.store <- rbind(mmse_boot.store, mmse_boot)
    m1_boot.store <- rbind(m1_boot.store, m1_boot)
    if (b>=B) break
  }
  EBM2 <- colMeans((eb_boot.store - mmse_boot.store)^2)
  EBMSE <- colMeans((pop_boot.store - eb_boot.store)^2)
  EBM1hat <- colMeans(m1_boot.store)
  EBM12 <- colMeans((mmse_boot.store-pop_boot.store)*(eb_boot.store-mmse_boot.store))
  
  return(data.frame(EBM2, EBMSE, EBM12, EBM1hat))
}
