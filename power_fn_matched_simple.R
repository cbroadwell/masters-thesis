##  Power function for matched pairs setting
#   Programmer: Carly Broadwell
#   Last updated: 06-25-2018
#                 07-13-2018  analyses: matched t-test, 2 types glmm, conditional gee, matched perm.test
#                 07-31-2018  removed random cluster effect from stratified GLMM (left it in metaanalysis as specified in Thompson paper)

source("/Users/carlyb/Desktop/Thesis/HCV_Simulation_Study/gen_matched_simple.R")
source("/Users/carlyb/Desktop/Thesis/HCV_Simulation_Study/matched_perm_test.R")

library(lme4)
library(geepack)
library(survival)

power_fn_matched_simple <- function(nreps,nc,c.size,u,s,te,ta,nperms){
  c_ttest <- 0
  c_glmm_strat <- 0
  c_glmm_meta <- 0
  c_glmm_meta_rand <- 0
  c_gee <- 0
  c_perm <- 0
  
  for (i in 1:nreps){
    #generate data and compute summary data
    dat <- gen_matched_simple(nc,c.size,u,s,te,ta)
    summary.dat <- as.data.frame(as.matrix(aggregate(y~MatchID + ClusterID + x, data = dat, FUN = function(z) c(mn = mean(z),v = abs(var(z))))))

    #matched pairs t-test
    if (as.numeric(t.test(y.mn ~ x, data = summary.dat, paired = TRUE)$p.value) < 0.05) {
      c_ttest <- c_ttest + 1
    }
    
    #stratified GLMM
    if (as.numeric(summary(glmer(y ~ x + strata(MatchID), data = dat, family = binomial))$coefficients[2,4]) < 0.05) {
      c_glmm_strat <- c_glmm_strat + 1
    }
    
    #GLMM based on metaanalysis approach
    if (as.numeric(summary(glmer(y ~ x + MatchID + (1 | ClusterID), data = dat, family =  binomial))$coefficients[2,4]) < 0.05) {
      c_glmm_meta <- c_glmm_meta + 1
    }
    
    #GLMM based on metaanalysis approach, but with MatchID as random effect 
    if (as.numeric(summary(glmer(y ~ x + (1 | MatchID) + (1 | ClusterID), data = dat, family = binomial))$coefficients[2,4]) < 0.05) {
      c_glmm_meta_rand <- c_glmm_meta_rand + 1
    }
    
    #stratified GEE
    if (as.numeric(summary(geese(y ~ x + strata(MatchID), data = dat, id = ClusterID, family = binomial, corstr = "exchangeable"))$mean[2,4]) < 0.05) {
      c_gee <- c_gee + 1
    } 
    
    #matched pairs permutation test
    if (matched.perm.test(summary.dat, nperms) < 0.05) {
      c_perm <- c_perm + 1
    }
  }
  
  c.vec <- c(c_ttest,c_glmm_strat,c_glmm_meta,c_glmm_meta_rand,c_gee,c_perm)
  power.vec <- c.vec / nreps
  
  return(power.vec)
}




