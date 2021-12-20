########################################################################
########################################################################
##Simulation Analysis Code:
########################################################################
########################################################################
##on my computer: 
#wd <- "~/Desktop/Dissertation/Paper1/Simulations_Code/"
.libPaths("/home/ssafo/haine108/Dissertation/Rcode/RLibs")

nsim = 1000
######MSI parameters: 
wd <- "/home/murra484/haine108/Paper1/Simulations_Code/Simulations_1Arm"
#######get job number/row of grid we're on dependent on scenario from MSI: 
ncores = 24
job = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
args=(commandArgs(T))

print(args[1])
print(args[2])



method = gsub("\\method=", "", args[1])
scenario = gsub("scenario=", "", args[2], fixed = T)


print(method)

print(scenario)

#setwd(wd)

####set parameter grid based on scenario: ##############
if(scenario == 1){
  ######scenario = 1: Bias = 0, w = 1 (ext is from same pop as RCT)
  w = c(0.90, 1); Bias = 0
  #####go across many different OR options: 
  logOR = seq(-2.5, 2.5, 0.5); OR = exp(logOR)
  grid_options = data.frame(expand.grid(w, Bias, logOR))
  colnames(grid_options) = c("w", "Bias", "logOR")
  print(dim(grid_options))
}
if(scenario == 2){
  ######scenario = 2: Bias = 0, w = varies (breaking MEMs)
  w = seq(0, 1, 0.2); Bias = 0
  #####go across many different OR options: 
  logOR = seq(-1, 1, 0.5); OR = exp(logOR)
  grid_options = data.frame(expand.grid(w, Bias, logOR))
  colnames(grid_options) = c("w", "Bias", "logOR")
  print(dim(grid_options))
}
if(scenario == 3){
  ######scenario = 3: Bias = varies, w = varies a little bit (breaking PSCL approach)
  w = seq(0.5, 1, 0.5); Bias = seq(-2, 2, 0.5)
  #####go across many different OR options: 
  logOR = seq(0, 1, 0.5); OR = exp(logOR)
  grid_options = data.frame(expand.grid(w, Bias, logOR))
  colnames(grid_options) = c("w", "Bias", "logOR")
  print(dim(grid_options))
}

####set current parameters: 
cur_params = grid_options[job, ]
############source the functions to be used: 
source(paste0(wd, "/Functions.R") )

########################################################################
########################################################################
##data generation function:
########################################################################
find_beta0_tau = function(pc, pt, nMCMC, no_covs){
  library(MASS)
  ##set up values:
  mean.vect = rep(0, no_covs);
  rho = 0.1; offdiag <- rho*1
  Sigma = diag(1, nrow = no_covs, ncol = no_covs)
  Sigma[lower.tri(Sigma)] <- offdiag
  Sigma[upper.tri(Sigma)] <-offdiag
  X <- mvrnorm(n = nMCMC,
               mean.vect,
               Sigma)
  beta = rep(1, no_covs)
  library(boot)
  
  ####find beta0 #########
  mean.b <- mean(X%*%beta) 
  beta0.vect <- seq(-3+logit( pc) - mean.b, 3+logit( pc) - mean.b, .01)
  mean_inv.logit <- sapply(1:length( beta0.vect), function(x) abs( pc - mean( inv.logit(beta0.vect[x]+X%*%beta)) ) )
  mean_inv.logit_ind <- which.min ( mean_inv.logit )
  beta0 = beta0.vect[mean_inv.logit_ind]
  
  X <- mvrnorm(n = nMCMC,
               mean.vect,
               Sigma)
  
  Y.ctrl  <- rbinom(n = nrow(X), size = 1,
                    prob = inv.logit(beta0+X%*%beta) )
  
  X <- mvrnorm(n = nMCMC,
               mean.vect,
               Sigma)
  
  ###find tau:########
  mean.b <- beta0+mean(X%*%beta)    
  tau.vect <- seq(-3+logit( pt) - mean.b, 3+logit( pt) - mean.b, .01)
  mean_inv.logit <- sapply(1:length( tau.vect), function(x) abs( pt - mean( inv.logit(tau.vect[x]+beta0+X%*%beta)) ) )
  mean_inv.logit_ind <- which.min ( mean_inv.logit )
  tau = tau.vect[mean_inv.logit_ind]
  
  X <- mvrnorm(n = nMCMC,
               mean.vect,
               Sigma)
  
  Y.trt  <- rbinom(n = nrow(X), size = 1,
                    prob = inv.logit(beta0+tau+X%*%beta) )
  
  ret.list = list(beta0 = beta0, tau = tau, check.ctrl = mean(Y.ctrl), check.trt = mean(Y.trt), pc = pc, pt = pt)
  
  return(ret.list)
  
}
####useful functions for parameter definitions: ###########
log_oddsratio <- function(pc, pt){
  logOR = log((pt/(1-pt))/(pc/(1-pc)))
  return(logOR)
}
odds <- function(prob){
  return(prob/(1-prob))
}
get_pt_fromOR <- function(pc, OR){
  pc_odds = odds(pc)
  num = pc_odds*OR
  den = 1+pc_odds*OR
  pt = num/den
  return(pt)
}

######function to generate the data:##############################
###inputs:
######mix_prop: from 0 to 1 if 1 that means all are from component 1 or in the RCT
data_generation <- function(mix_prop, logOddsRatio, bias_ext){
  #####across all scenarios we have: 
  pctrl = 0.30; num_covs = 5; n_trial = 200; n_ext = 1000
  ext_dat_mean = c(0, -5)
  ###find the beta0 and the tau for these parameters: #######
  OR <- exp(logOddsRatio) 
  ptrt <- get_pt_fromOR(pctrl, OR )
  beta <- rep(1, num_covs)
  beta0_tau <- find_beta0_tau(pctrl, ptrt, 10000, num_covs)
  
  #####generate X matrices: #####
      library(MASS)
      #########trial X matrix: #####
        mean.vect <- rep(0, num_covs);
        rho <- 0.1; offdiag <- rho*1
        Sigma <- diag(1, nrow = num_covs, ncol = num_covs)
        Sigma[lower.tri(Sigma)] <- offdiag
        Sigma[upper.tri(Sigma)] <-offdiag
        X_trial <- mvrnorm(n = n_trial,
               mean.vect,
               Sigma)
      #########external data X matrix: #####
      components <- sample(1:2,prob=c(mix_prop, 1-mix_prop),size=n_ext,replace=TRUE)
      mus <- c(0, -5)
      X_ext <- do.call(rbind, lapply(components, 
                                          function(x){
                                            mvrnorm(n=1,
                                                    mu=rep(mus[x], num_covs),
                                                    Sigma) } ) ) 
      
  #####generate outcomes: #####
    library(boot)
    ########trial choose treated and control arm outcomes:#####
      Y_ctrl <- rbinom(n_trial/2, 1, prob = inv.logit(beta0_tau$beta0+X_trial%*%beta) )
      Y_trt <- rbinom(n_trial/2, 1, prob = inv.logit(beta0_tau$beta0+beta0_tau$tau+X_trial%*%beta) )
    ########external RWD outcomes: #####
      Y_ext <- rbinom(n_ext, 1, prob = inv.logit(beta0_tau$beta0+X_ext%*%beta+bias_ext) )
      
  #####generate propensity scores: #####
      ps_trial_dat <- data.frame(cbind( data.frame(Y = c(Y_ctrl, Y_trt), 
                                                   Trt = c( rep(0, n_trial/2), rep(1, n_trial/2)), 
                                                   Trial = 1), X_trial) )
      ps_ext_dat <- data.frame(cbind( data.frame(Y = Y_ext, 
                                                   Trt = c( rep(0, n_ext) ), 
                                                   Trial = 0), X_ext) )
      allps_dat <- rbind(ps_trial_dat, ps_ext_dat)
      PS_all <- glm(as.formula( paste0("Trial~", paste0("X", seq(1, num_covs), collapse = "+") ) ), 
                      data = allps_dat, family = binomial(link = "logit"))
      #########trial PS: #####
      PS_trial = predict(PS_all, newdata = ps_trial_dat, type = "response")
      
      #########external RWD PS: #####
      PS_ext = predict(PS_all, newdata = ps_ext_dat, type = "response")
  ########save as a list with each element of the list specified as 
  #################### a data frame with the following named columns: Y, PS, logitPS, Trt
      library(dplyr); library(tidyr)
      trial <- data.frame(cbind( data.frame(Y = c(Y_ctrl, Y_trt), 
                                            PS = PS_trial,
                                            Trt = c( rep(0, n_trial/2), rep(1, n_trial/2)), 
                                            Trial = 1), X_trial) ) %>% mutate(logitPS = logit(PS))
      ext <- data.frame(cbind( data.frame(Y = Y_ext, 
                                          Trt = c( rep(0, n_ext) ), 
                                          Trial = 0, 
                                          trueZ = components, 
                                          PS = PS_ext), X_ext) ) %>% mutate(logitPS = logit(PS)) 
      dat_gen = list()
      dat_gen[[1]] = trial; dat_gen[[2]] = ext
  
  return(dat_gen)
}
########################################################################
########################################################################
##Code to analyse the simulation results: 
######inputs are: true parameters list with the following: (pc, pt, logOR)
######method: string with options of: "SS_MIX_MEM_Capped", "MEM_Capped", "PSCL", etc. 
######NB_results: list from the run_sims done for the no borrowing method
######n_rct_ctrl: number of individuals in the control arm of the rct
########################################################################
analyse_sims <- function(parameters, sim_results, NB_results, n_rct_ctrl){
  ESSS<- function(var_NB, var_B, N){
    esss = ((var_NB/var_B)-1)*N
    return(esss)
  }
  # print(parameters$logOR)
  # print(parameters$pt)
  # print(parameters$pc)
  
  ##get the Bias, Variance, Relative MSE, and ESSS from the simulations
  ##for the log(OR) term values####
  truelogOR <- parameters$logOR
  bias_logOR <- sim_results$logOR$est - truelogOR
  var_logOR <- sim_results$logOR$var
  
  MSE_logOR <- (bias_logOR^2+var_logOR)
  ESSS_logOR <- ESSS(NB_results$logOR$var, var_logOR, n_rct_ctrl)
  logOR <- list(bias = bias_logOR, 
                var = var_logOR, 
                MSE = MSE_logOR, 
                ESSS = ESSS_logOR)
  
  ##for the pc term values####
  truepc <- parameters$pc
  bias_pc <- sim_results$pc$est - truepc
  var_pc <- sim_results$pc$var
  MSE_pc <- (bias_pc^2+var_pc)
  ESSS_pc <- ESSS(NB_results$pc$var, var_pc, n_rct_ctrl)
  pc <- list(bias = bias_pc, 
             var = var_pc, 
             MSE = MSE_pc, 
             ESSS = ESSS_pc)
  ##for the pt term values####
  truept <- parameters$pt
  bias_pt <- sim_results$pt$est - truept
  var_pt <- sim_results$pt$var
  MSE_pt <- (bias_pt^2+var_pt)
  ESSS_pt <- ESSS(NB_results$pt$var, var_pt, n_rct_ctrl)
  pt <- list(bias = bias_pt, 
             var = var_pt, 
             MSE = MSE_pt, 
             ESSS = ESSS_pt)
  
  ret.list = list(pc = pc, pt = pt, logOR = logOR)
  return(ret.list)
}
########################################################################
########################################################################
##Code to run the simulation: 
########################################################################
run_sims <- function(nsim, method, parameters, scenario, job, n_rct_ctrl = 50, wd, no_cores){
  library(parallel, quietly = TRUE)
  print(paste0( "job = ", job) )
  ####for NB approach generate all the data and save all the NB results for each dataset:
  new_wd = paste0(wd, "/Scenario", scenario)
  
  ###set seed: 
  set.seed(1993+1990+job)
  
  ####get current parameter values: 
  print(paste0("got to current parameter values:"))
  cur_mix = parameters$w; logOR_cur = parameters$logOR; bias_cur = parameters$Bias
  trueparams = list(logOR = logOR_cur, 
                    pc = 0.30, 
                    pt = get_pt_fromOR(0.30, exp(logOR_cur)))
  if(method == "NB"){
    print(paste0("got to data generation"))
    dat_cur <- lapply(1:nsim, function(x) data_generation(cur_mix, logOR_cur, bias_cur))
    print(paste0("got to NB_run"))
    NB_cur <- mclapply(1:nsim, function(x) NB(dat_cur[[x]], 10000, 1000), mc.cores = no_cores)
    ####save results for other methods to use: 
    print(paste0("got to saving data:"))
    dat_name = paste0(new_wd, "/data/scen_", scenario, "job_", job, ".RDS" )
    saveRDS(dat_cur, file = dat_name)
    print(paste0("got to saving NB res:"))
    NB_name = paste0(new_wd, "/NB/scen_", scenario, "job_", job, ".RDS" )
    saveRDS(NB_cur, file = NB_name)
    print(paste0("got to analyse sims:"))
    print((NB_cur[[1]]))
    print(NB_cur[[nsim]])
    analysis = lapply(1:nsim, function(x) analyse_sims(trueparams, NB_cur[[x]], NB_cur[[x]], n_rct_ctrl)) 
  }
  if(method != "NB"){
    dat_cur <- readRDS( paste0(new_wd, "/data/scen_", scenario, "job_", job, ".RDS" ) )
    print(paste0("read in data"))
    NB_res <- readRDS( paste0(new_wd, "/NB/scen_", scenario, "job_", job, ".RDS" ) )
    print(paste0("got NB res"))
    if(method == "PSCL"){
      run_PSCL = function(dat, i){
        #print(i)
        temp <- list()
        r <- NULL
        attempt <- 1
        while( is.null(r) && attempt <= 20 ) {
          print(paste0("attempt: ", attempt) )
          attempt <- attempt + 1
          try(
            temp <- PSCL(dat)
          )
        }
        return(temp)
      }
      res = mclapply(1:nsim, function(x) {run_PSCL(dat_cur[[x]],x)}, mc.cores = no_cores-5)
      }
       
      
    if(method == "SS_MIX_MEM_Capped"){res = mclapply(1:nsim, function(x) {
      SS_MIX_MEM(dat_cur[[x]], 3000, 8000, 3000, 1000, "capped", "control")}, mc.cores = no_cores ) }
    if(method == "SS_MIX_MEM_FB"){res = mclapply(1:nsim, function(x) {
     SS_MIX_MEM(dat_cur[[x]], 3000, 8000, 3000, 1000, "FB", "control")} , mc.cores = no_cores) }
    if(method == "MEM_Capped"){res = mclapply(1:nsim, function(x) {MEM(dat_cur[[x]], 20000, 10000, "capped") } , mc.cores = no_cores) }
    if(method == "MEM_FB"){res = mclapply(1:nsim, function(x) {MEM(dat_cur[[x]], 20000, 10000, "FB") } , mc.cores = no_cores) }
    
    print(paste0("got to analyse sims:"))
    print((res[[1]]))
    print(res[[nsim]])
    analysis = lapply(1:nsim, function(x) {analyse_sims(trueparams, res[[x]], NB_res[[x]], n_rct_ctrl)})
  }
  
  ######get means of analysis results for each term: 
  library(dplyr); library(tidyr)
  print(paste0("got to pcres:"))
  pc_res =  colMeans( do.call(rbind, lapply(1:nsim, function(x) as.numeric (analysis[[x]]$pc) )), na.rm = T  )  
  pc_res = data.frame(bias = pc_res[1], var = pc_res[2], MSE = pc_res[3], ESSS = pc_res[4],  
                       method = method, 
                       metric = "pc", 
                       scenario = scenario, 
                       true_w = parameters$w, 
                       true_logOR = logOR_cur, 
                       true_bias = bias_cur )
  print(paste0("got to pt res:"))
  pt_res =  colMeans( do.call(rbind, lapply(1:nsim, function(x) as.numeric (analysis[[x]]$pt) )), na.rm = T)  
  pt_res = data.frame(bias = pt_res[1], var = pt_res[2], MSE = pt_res[3], ESSS = pt_res[4],  
                      method = method, 
                      metric = "pt", 
                      scenario = scenario, 
                      true_w = parameters$w, 
                      true_logOR = logOR_cur, 
                      true_bias = bias_cur )
  print(paste0("got to logOR:"))
  logOR_res =  colMeans( do.call(rbind, lapply(1:nsim, function(x) as.numeric (analysis[[x]]$logOR) )), na.rm = T  )  
  logOR_res = data.frame(bias = logOR_res[1], var = logOR_res[2], MSE = logOR_res[3], ESSS = logOR_res[4],  
                      method = method, 
                      metric = "logOR", 
                      scenario = scenario, 
                      true_w = parameters$w, 
                      true_logOR = logOR_cur, 
                      true_bias = bias_cur )
  print(paste0("got to saving all res:"))
  all_res = data.frame(rbind(pc_res, pt_res, logOR_res))
  if(method == "PSCL"){
    saveRDS(res, paste0(new_wd, "/PSCL/scen_", scenario, "job_", job, "method_", method, ".RDS" ) )}
  if(method %in% c("SS_MIX_MEM_Capped", "SS_MIX_MEM_FB") ){
    weights = mean( do.call(rbind, lapply(1:nsim, function(x) as.numeric (res[[x]]$weights_ctrl) ))  ) 
    pip = mean( do.call(rbind, lapply(1:nsim, function(x) as.numeric (res[[x]]$pip_ctrl) ))  ) 
    mix = mean( do.call(rbind, lapply(1:nsim, function(x) as.numeric (res[[x]]$mix_param_ctrl) ))  ) 
    mem = mean( do.call(rbind, lapply(1:nsim, function(x) as.numeric (res[[x]]$mem_ctrl) ))  ) 
    mix_mem = data.frame("mix" = mix, "mem" = mem, "weights" = weights)
    write.table(mix_mem, paste0(new_wd, "/", method,"/scen_", scenario, "job_", job, "_mix_mems_info_", method, ".txt" ) )}
  if(method %in% c("MEM_Capped", "MEM_FB") ){
    weights = mean( do.call(rbind, lapply(1:nsim, function(x) as.numeric (res[[x]]$weights_ctrl) ))  ) 
    pip = mean( do.call(rbind, lapply(1:nsim, function(x) as.numeric (res[[x]]$pip_ctrl) ))  ) 
    mem = mean( do.call(rbind, lapply(1:nsim, function(x) as.numeric (res[[x]]$mem_ctrl) ))  ) 
    mix_mem = data.frame("mem" = mem, "weights" = weights, "pip" = pip)
    write.table(mix_mem, paste0(new_wd, "/", method,"/scen_", scenario, "job_", job, "_mems_info_", method, ".txt" ) )}
  
  write.table(all_res, paste0(new_wd, "/", method, "/scen_", scenario, "_job_", job, "_method_", method, ".txt" ))
  
  # ret.list = list(pc = pc_res, pt = pt_res, logOR = logOR_res)
  # 
  # return(ret.list)
  
}

run_sims(nsim, method, cur_params, scenario, job, n_rct_ctrl = 100, wd, ncores)

# # ###test with 2 nsim: 
# nsim = 2; cur_params = grid_options[job, ]
# NB_check = run_sims(nsim, "NB", cur_params, 1, 1, 50, wd, 2)
# ####check all methods:
# methods.vect = c("NB","PSCL", "MEM_Capped", "MEM_FB", "SS_MIX_MEM_Capped", "SS_MIX_MEM_FB")
# 
# system.time( {all_methods = lapply(1:length( methods.vect), function(x) 
#   run_sims(nsim, methods.vect[x], parameters, 1, 1, 50, wd, 2) )} )
