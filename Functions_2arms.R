########################################################################
########################################################################
##Functions used to run each method/analyze methods in for 2 arms sims:
######NOTE: THE PSCL APPROACH IS NOT APPLICABLE HERE  as it was developed
#######ONLY for incorporating information on the control arm
########################################################################
########################################################################
##on my computer: 
#wd <- "~/Desktop/Dissertation/Paper1/Simulations_Code/Simulations_2Arms"

##on MSI: 
wd <- "/home/murra484/haine108/Paper1/Simulations_Code/Simulations_2Arms"

#setwd(wd)
########################################################################
########################################################################
##inputs: 
########################################################################
########################################################################
#all functions below take the following inputs: 
#dat = list of data frames with the first as the trial
#the data frames contain the following named columns: PS, logitPS, Y, Trt, trial
##M (for MCMC approaches) and nburnin (MEM)
##M1 and M2 and nburnin1 and nburnin2 for the (SS-MIX-MEM approach) 
##K for the PSCL and the SS-MIX-MEM approach (if capped)

########################################################################
########################################################################
##No Borrowing Function:
########################################################################
NB <- function(dat, M, nburnin){
      Mall <- M+nburnin
      library(dplyr); library(tidyr)
      ###############functions to be used: 
      ###################################
      ####mean and variance of a beta with specified a and b:
      mean_beta <- function(a, b){
          mean <- a/(a+b)
          return(mean)
      }
      ####variance of a beta with specified a and b:
      var_beta <- function(a, b){
        var<- (a*b)/(((a+b)^2)*(a+b+1))
        return(var)
      }
      ####odds ratio:
      odds_ratio <- function(p1, p2){
        OR <- ( p1/(1-p1) ) / (p2/(1-p2))
        return(OR)
      }
      ###set up data: ##########
      RCT <- dat[[1]]
      treated <- filter(RCT, Trt == 1)
      untreated <- filter(RCT, Trt == 0)
      
      ###control arm draws: #####
      y_ctrl <- sum(untreated$Y); n_ctrl <- nrow(untreated)
      a_ctrl <- 1+y_ctrl; b_ctrl <- 1+n_ctrl-y_ctrl
      ctrl_draws <- rbeta(Mall, a_ctrl, b_ctrl)
      #print(head(ctrl_draws))
      ###treat arm draws: #####
      y_trt <- sum(treated$Y); n_trt <- nrow(treated)
      a_trt <- 1+y_trt; b_trt <- 1+n_trt-y_trt
      trt_draws <- rbeta(Mall, a_trt, b_trt)
      #print(head(trt_draws))
      ###odds ratio draws: ####
      #jac <- c(1/est_ptrt+(1/(1-est_ptrt)), 0, 0, -(1/(1-est_pct))-(1/est_pct))
      OR_draws <- data.frame( cbind(ctrl_draws, trt_draws) )
      colnames(OR_draws) = c("pctrl", "ptrt" )
      #print(head(OR_draws))
      logOR <- OR_draws %>% 
        mutate(logOR = log(ptrt)-log(1-ptrt)+log(1-pctrl)-log(pctrl) )
      #print(head(logOR))
      logOR <- logOR[-c(1:nburnin), ]
      #print(head(logOR))
      
      ####mean and var of OR: ####
      est <- mean(logOR$logOR)
      var_est <- var(logOR$logOR)
      CI <- quantile(logOR$logOR, probs = c(0.025, 0.975))
      ###means and variances of pct and ptrt: 
      est_pct <- mean_beta(a_ctrl, b_ctrl); est_ptrt = mean_beta(a_trt, b_trt)
      var_pct <- var_beta(a_ctrl, b_ctrl); var_ptrt = var_beta(a_trt, b_trt)
      CI_pct <- quantile(ctrl_draws, probs = c(0.025, 0.975)); CI_ptrt <- quantile(trt_draws, probs = c(0.025, 0.975))
      ###mean and var of log OR delta method: 
      # library(msm)
      # est_logOR_delta <- log(est_ptrt)-(log(1-est_ptrt))+log(1-est_pct)-log(est_pct)
      # 
      # var_logOR_delta <- ( deltamethod(~ log((x1*(1-x2))/( x2*(1-x1)) ),
      #                                c(est_ptrt, est_pct),
      #                                matrix(c( var_ptrt, 0, 0, var_pct), ncol = 2 )) )^2

      
      ret.list = list(logOR = list(est = est, 
                                   var = var_est, 
                                   ci = CI), 
                      pc = list(est = est_pct, 
                                 var = var_pct, 
                                ci = CI_pct), 
                      pt = list(est = est_ptrt, 
                                  var = var_ptrt, 
                                ci = CI_ptrt))
      
      return(ret.list)
      
}

# NB_check <- NB(dat, M, nburnin)

########################################################################
########################################################################
##MEM Function:
########################################################################
##source the Cpp function: (for speed)
library(Rcpp)
Rcpp::sourceCpp(paste0(wd, "/MEM_cpp.cpp") )
#############function for capped prior version: 
###solves for the pip for the specified external data set:
solve_prior = function(nc, no, yo, ac, bc, ao, bo, K){
  capped_MEM.vals = capped_MEM(nc, no, yo, ac, bc, K)
  C.vect = capped_MEM.vals$x.val
  #####parameters under maximal borrowing UNEXCHANGEABLE
  yc = nc*(yo/no)
  ac0 = ac+yc; bc0 = bc+nc-yc
  ao0 = ao+yo; bo0 = bo+no-yo
  
  #####parameters under maximal borrowing EXCHANGEABLE
  ac1 = ac+yc+yo; bc1 = bc+nc-yc+no-yo
  
  ####the marginal likelihood under the non-exchangeable assumption: (beta())/()*()/()
  P0 = exp(lbeta(ac0,bc0)-lbeta(ac, bc)+lbeta(ao0,bo0)-lbeta(ao,bo)  )
  P1 = exp(lbeta(ac1,bc1)-lbeta(ac,bc) )
  #print(P0)
  #print(P1)
  
  if( ((P0 == P1) & P0 == 0) || (is.na(P0) )|| is.na(P1)) {
    P0 = .1; P1 = .9
    #print("NOTE TO SELF: TRY BIC")
  }
  
  if(length(C.vect) == 1){
    C = C.vect
    num = P0*(1-C)
    den =  C*P1-C*P0+P0
    pi1 = num/den
  }
  if(length(C.vect) == 2){
    pi1.v = c()
    for (i in 1:2){
      C = C.vect[i]
      num = P0*(1-C)
      den =  C*P1-C*P0+P0
      pi1.v[i] = num/den
    }
    pi1 = min(pi1.v)
  }
  #print(pi1)
  if(is.na(pi1)){pi1 = 0}
  pip = min(pi1, 0.5)
  #print(pip)
  return(pip)
  
}
capped_MEM = function(nc, no, yo, ac, bc, K){
  mean_beta <- function(a, b){
    mean <- a/(a+b)
    return(mean)
  }
  var_beta <- function(a, b){
    var<- (a*b)/(((a+b)^2)*(a+b+1))
    return(var)
  }
  K = min(no, K)
  ###1 is the noborrowing approach in the SS-MIX/MEM 
  ac1 = ac+(nc*(yo/no)); bc1 = bc+nc-(nc*(yo/no))
  #print(ac1)
  ac2 = ac+yo+(nc*(yo/no)); bc2 = bc+nc+no-yo-(nc*(yo/no))
  #print(ac2)
  mu1 = mean_beta(ac1, bc1); #sigma1 = var(ac1, bc1)
  mu2 = mean_beta(ac2, bc2); #sigma2 = var(ac2,bc2)
  
  # print(mu1)
  # print(mu2)
  sigma1 = var_beta(ac1, bc1)
  sigma2 = var_beta(ac2, bc2)
  #print(sigma1)
  #print(sigma2)
  sigmaNB = sigma1
  A = (2*mu1*mu2-mu2^2-mu1^2)
  B = (sigma1+mu1^2-sigma2-mu2^2-2*mu1*mu2+2*mu2^2)
  C = sigma2+mu2^2-mu2^2-(sigmaNB/(1+(K/nc)))
  # print(A)
  # print(B)
  # print(C)
  x1 = (-B+sqrt(B^2-4*A*C))/(2*A)
  x2 = (-B-sqrt(B^2-4*A*C))/(2*A)
  
  x.pos = c(x1, x2)
  x.val = ifelse( (x1 <= 1 && x1 >= 0) && (x2 <= 1 && x2 >= 0), x.pos, 
                  ifelse(x1 <= 1 && x1 >= 0, x1,  
                         ifelse(x2 <= 1 && x2 >= 0, x2, 1)) )
  ###check:
  #print(A*x1^2+B*x1+C)
  #print(A*x2^2+B*x2+C)
  ret.list = list(pos = x.pos, 
                  x.val = x.val)
  return(ret.list)
}
#######mem function: 
MEM <- function(dat, M, nburnin, prior.text = "FB"){
  library(dplyr); library(tidyr)
  ###set up data: ##########################
  #####trial data:
  RCT <- dat[[1]]
  treated <- filter(RCT, Trt == 1)
  untreated <- filter(RCT, Trt == 0)
  ####treated external:
  ext <- dat[[2]]
  trt_ext <- filter(ext, Trt == 1)
  untrt_ext <- filter(ext, Trt == 0)
  
  ####get sufficient statistics for the trial, and the external: 
  r_trt_t = sum(treated$Y); n_trt_t = nrow(treated)
  r_ctrl_t = sum(untreated$Y); n_ctrl_t = nrow(untreated)
  if(nrow(trt_ext) == 0){
    r_trt_ext = n_trt_ext = 0
  }
  else{
    r_trt_ext = sum(trt_ext$Y); n_trt_ext = nrow(trt_ext)
  }
  if(nrow(untrt_ext) == 0){
    r_untrt_ext = n_untrt_ext = 0
  }
  else{
    r_untrt_ext = sum(untrt_ext$Y); n_untrt_ext = nrow(untrt_ext)
  }

  ####get the prior values: 
  if(prior.text == "FB"){ pip_ctrl = pip_trt = 0.50 }
  if(prior.text == "capped"){ 
    pip_ctrl = solve_prior(n_ctrl_t, n_untrt_ext, r_untrt_ext, 1, 1, 1, 1, n_ctrl_t)  
    pip_trt = solve_prior(n_trt_t, n_trt_ext, r_trt_ext, 1, 1, 1, 1, n_trt_t) 
  }
  

  ######get all the draws using Rcpp version of MEM function:###############
  ##all hyperparams are 1 for the a, b of beta distributions:
  list_draws = MEM_binary(r_ctrl = r_ctrl_t, n_ctrl = n_ctrl_t,
                          r_trt = r_trt_t, n_trt = n_trt_t, 
                          ext_r_ctrl = r_untrt_ext, ext_N_ctrl = n_untrt_ext,
                          ext_r_trt = r_trt_ext, ext_N_trt = n_trt_ext,
                          nburnin, M, pip_ctrl, pip_trt, 1,1,1,1,1,1,1,1)
  pc_res = list(est = mean(list_draws$pc[-c(1:nburnin)]), 
                var = var(list_draws$pc[-c(1:nburnin)]), 
                ci = quantile(list_draws$pc[-c(1:nburnin)], probs = c(0.025, 0.975)))
  pt_res = list(est = mean(list_draws$pt[-c(1:nburnin)]), 
                var = var(list_draws$pt[-c(1:nburnin)]),
                ci = quantile(list_draws$pt[-c(1:nburnin)], probs = c(0.025, 0.975)))
  logOR_res = list(est = mean(list_draws$logOR[-c(1:nburnin)]), 
                var = var(list_draws$logOR[-c(1:nburnin)]),
                ci = quantile(list_draws$logOR[-c(1:nburnin)], probs = c(0.025, 0.975)))
  ind_mem_ctrl = mean(list_draws$mem_ctrl[-c(1:nburnin)])
  ind_mem_trt = mean(list_draws$mem_trt[-c(1:nburnin)])
  weights_ctrl = mean( list_draws$weights_ctrl[2] ) 
  weights_trt = mean( list_draws$weights_trt[2] ) 
  results = list(pc = pc_res, 
                 pt = pt_res, 
                 logOR = logOR_res, 
                 mem_ctrl = ind_mem_ctrl, 
                 mem_trt = ind_mem_trt, 
                 weights_ctrl = weights_ctrl, 
                 weights_trt= weights_trt,
                 pip_ctrl = pip_ctrl,
                 pip_trt = pip_trt)
  return(results)
  
  
}
########################################################################
########################################################################
##SS-MIX-MEM Function:
########################################################################
##source the Cpp functions: 
Rcpp::sourceCpp(paste0(wd, "/SS_MIX_cpp.cpp") )
SS_MIX_MEM = function(dat, M1, M2, nburnin1, nburnin2, prior.text){
  library(dplyr); library(tidyr)
  ####get dat set up for functions: ######
  RCT <- dat[[1]]; PS_rct <- RCT$logitPS
  treated <- filter(RCT, Trt == 1)
  untreated <- filter(RCT, Trt == 0)
  
  ext <- dat[[2]]
  ext_trt <- filter(ext, Trt == 1)
  ps_trt_ext <- ext_trt$logitPS; y_trt_ext <- ext_trt$Y
  ext_ctrl <- filter(ext, Trt == 0)
  ps_ctrl_ext <- ext_ctrl$logitPS; y_ctrl_ext <- ext_ctrl$Y
  
  ##step 1: run the SS-MIX step for each data frame:#######
  Mall1 = M1+nburnin1
  ##step 2: run the MEM step for each draw from each data frame:(so get M1*M2 total draws)#######
  ####get sufficient statistics for the trial:
  r_trt_t = sum(treated$Y); n_trt_t = nrow(treated)
  r_ctrl_t = sum(untreated$Y); n_ctrl_t = nrow(untreated)
  
  if(nrow(ext_trt) == 0){borrow_text = "control"}
  if(nrow(ext_trt) > 0){borrow_text = "both"}
  print(borrow_text)
  #####get the sufficient statistics and the pip for the external data:
  if(borrow_text == "both"){
    SS_ctrl_all <- SS_MIX(PS_rct, ps_ctrl_ext,  y_ctrl_ext, Mall1,
                          0.00001, .001,.001, 0.00001, .001,  .001,  1, 1)
    SS_ctrl <- cbind(SS_ctrl_all$yobsIN, SS_ctrl_all$nobsIN)[-c(1:nburnin1), ]
    SS_trt_all <- SS_MIX(PS_rct, ps_trt_ext,  y_trt_ext, Mall1,
                         0.00001, .001,.001, 0.00001, .001,  .001,  1, 1)
    SS_trt <- cbind(SS_trt_all$yobsIN, SS_trt_all$nobsIN)[-c(1:nburnin1), ]
    w_ctrl <- c(SS_ctrl_all$w)[-c(1:nburnin1)]; w_trt = c(SS_trt_all$w)[-c(1:nburnin1)]
    mix_ctrl = mean(SS_ctrl_all$w[-c(1:nburnin1)]); mix_trt = mean(SS_trt_all$w[-c(1:nburnin1)])
    n_in_ctrl = mean(SS_ctrl[,2]); y_in_ctrl = mean(SS_ctrl[,1])
    n_in_trt = mean(SS_trt[,2]); y_in_trt = mean(SS_trt[,1])
    ####if prior.text = "capped" 
    if(prior.text == "FB") {pip_ctrl <- pip_trt <- 0.5}
    if(prior.text == "capped"){
      pip_ctrl <- mean(do.call(rbind, 
                               lapply(1:M1, function(x){ solve_prior(n_ctrl_t, SS_ctrl[x, 2], SS_ctrl[x, 1], 1, 1, 1, 1, n_ctrl_t)}  ) ))
      pip_trt <- mean(do.call(rbind, 
                               lapply(1:M1, function(x){ solve_prior(n_trt_t, SS_trt[x, 2], SS_trt[x, 1], 1, 1, 1, 1, n_trt_t)}  ) ))
    }
    #####step 2 a): if prior.text == "capped" then solve for PIP for control and treated arms: 
  }
  if(borrow_text == "control"){
    SS_trt_all <- NA
    SS_ctrl_all <- SS_MIX(PS_rct, ps_ctrl_ext,  y_ctrl_ext, Mall1,
                          0.00001, .001,.001, 0.00001, .001,  .001,  1, 1)
    SS_ctrl <- cbind(SS_ctrl_all$yobsIN, SS_ctrl_all$nobsIN)[-c(1:nburnin1), ]
    SS_trt <- cbind( rep(0, M1), rep(0, M1) )
    mu_ctrl <- c(SS_ctrl_all$muct)[-c(1:nburnin1)]
    w_ctrl <- c(SS_ctrl_all$w)[-c(1:nburnin1)]
    mix_trt = NA; mix_ctrl = mean(SS_ctrl_all$w[-c(1:nburnin1)] )
    n_in_ctrl = mean(SS_ctrl[,2]); y_in_ctrl = mean(SS_ctrl[,1])
    n_in_trt = mean(SS_trt[,2]); y_in_trt = mean(SS_trt[,1])
    if(prior.text == "FB") {pip_ctrl <- pip_trt <- 0.5}
    if(prior.text == "capped"){
      pip_ctrl <- mean(do.call(rbind, 
                               lapply(1:M1, function(x){ solve_prior(n_ctrl_t, 
                                                                     SS_ctrl[x, 2], 
                                                                     SS_ctrl[x, 1], 1, 1, 1, 1, 
                                                                     n_ctrl_t)}  ) ))
      pip_trt <- 0.5
    }
  }
  print(paste0("got past SS-MIX Step:"))
  ######Run the MEM step:
  list_draws = lapply(1:M1, function(x) { MEM_binary(r_ctrl = r_ctrl_t, n_ctrl = n_ctrl_t,
                          r_trt = r_trt_t, n_trt = n_trt_t, 
                          ext_r_ctrl = SS_ctrl[x, 1] , ext_N_ctrl = SS_ctrl[x, 2],
                          ext_r_trt = SS_trt[x, 1], ext_N_trt = SS_trt[x, 2],
                          nburnin2, M2, pip_ctrl, pip_trt, 1,1,1,1,1,1,1,1) } )
  pc_draws = c( sapply(1:M1, function(x) {list_draws[[x]]$pc[-c(1:nburnin2)]}) )
  pc_res = list(est = mean(pc_draws), 
                var = var(pc_draws), 
                ci = quantile(pc_draws, probs = c(0.025, 0.975)))
  pt_draws = c( sapply(1:M1, function(x) {list_draws[[x]]$pt[-c(1:nburnin2)]}) )
  pt_res = list(est = mean(pt_draws), 
                var = var(pt_draws),
                ci = quantile(pt_draws, probs = c(0.025, 0.975)))
  logOR_draws = c( sapply(1:M1, function(x) {list_draws[[x]]$logOR[-c(1:nburnin2)]}) )
  logOR_res = list(est = mean(logOR_draws), 
                   var = var(logOR_draws),
                   ci = quantile(logOR_draws, probs = c(0.025, 0.975)))
  ind_mem_ctrl_draws = c( sapply(1:M1, function(x) {list_draws[[x]]$mem_ctrl[-c(1:nburnin2)]}) )
  ind_mem_ctrl = mean(ind_mem_ctrl_draws)
  ind_mem_trt_draws = c( sapply(1:M1, function(x) {list_draws[[x]]$mem_trt[-c(1:nburnin2)]}) )
  ind_mem_trt = mean(ind_mem_trt_draws)
  weights_ctrl = do.call(rbind, lapply(1:M1, function(x) {list_draws[[x]]$weights_ctrl[2]}) )
  weights_trt = do.call(rbind, lapply(1:M1, function(x) {list_draws[[x]]$weights_trt[2]}) ) 
  print(paste0("got past MEMs Step:"))
  results = list(pc = pc_res, 
                 pt = pt_res, 
                 logOR = logOR_res, 
                 mem_ctrl = ind_mem_ctrl, 
                 mem_trt = ind_mem_trt,
                 weights_ctrl = mean( weights_ctrl ), 
                 weights_trt= mean (weights_trt ), 
                 mix_param_ctrl = mix_ctrl, 
                 mix_param_trt = mix_trt, 
                 n_in_ctrl = n_in_ctrl, 
                 y_in_ctrl = y_in_ctrl, 
                 n_in_trt = n_in_trt, 
                 y_in_trt = y_in_trt, 
                 pip_ctrl = pip_ctrl, 
                 pip_trt = pip_trt, 
                 # draws_mem_ctrl = ind_mem_ctrl_draws, 
                 # weights_mem_ctrl = weights_ctrl,
                 # draws_w = w_ctrl,
                 # draws_muct = mu_ctrl)
                 Z_vec_ctrl = SS_ctrl_all$Zvec[-c(1:nburnin1), ], 
                 Z_vec_trt = SS_trt_all$Zvec[-c(1:nburnin1), ])  
  return(results)
  
  
}
#######################################################################
##generate data for example for each of the functions above: 
#######################################################################
#######################################################################
# library(boot)
# PS_1 = rnorm(200, 0, 1);
# comps = rbinom(1000, 1, 0.50); means = c(0, -4); sds = c(1, 1)
# PS_2 = rnorm(1000, means[comps+1], sds[comps+1])
# probs = c(0.30, 0.60); Y_trt = rbinom(100, 1, prob = 0.7)
# Y_ctrl = rbinom(100, 1, prob = 0.30); Y_2 = rbinom(1000, 1, prob = probs[comps+1])
# 
# M1 = 10000; nburnin1 = 1000; all = M1+nburnin1
# M = 10000; nburnin = 1000
# rct.dat = data.frame(logitPS = PS_1,
#                      Y = c(Y_ctrl, Y_trt),
#                      Trt = c(rep(0, length(Y_ctrl)), rep(1, length(Y_trt))) ) %>% mutate( PS = inv.logit(logitPS))
# 
# ext.dat = data.frame(logitPS = PS_2,
#                      Y = Y_2,
#                      Trt = rep(0, 1000) ) %>% mutate( PS = inv.logit(logitPS))
# dat = list(); dat[[1]] = rct.dat; dat[[2]] = ext.dat
# ######################################################################
# #example run for each method:
# ######################################################################
# MEM_FB <- MEM(dat, M, nburnin, "FB")
# MEM_capped <- MEM(dat, M, nburnin, "capped")
# SS_MIX_MEM_CHECK_FB = SS_MIX_MEM(dat, 100, 100, 10, 10, "FB", "control")
# SS_MIX_MEM_CHECK_capped = SS_MIX_MEM(dat, 100, 100, 10, 10, "capped", "control")
# NB_check <- NB(dat, M, nburnin)
# PSCL_check <- PSCL(dat)
# 
# 
