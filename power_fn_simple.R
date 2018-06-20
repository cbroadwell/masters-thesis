# Writing power function for simple case
# Last updated: 2018-06-12
# Programmer: Carly Broadwell
# Purpose:  Write a function to calculate the power for the 4 tests of interest
# Parameters that can be varied:
#   nreps: the number of datasets to generate and analyze/collect p-values from
#   Parameters from "simple" data generation:
#     Nclust: the number of clusters
#     csizes: the sizes of the clusters
#     u: odds of outcome at baseline
#     sj: standard deviation of cluster intercepts
#     trteff: ODDS RATIO - TREATMENT EFFECT - v important param to vary
#   Parameters from tests:
#     nperm: number of permutations for permutation test

source("/Users/carlyb/Desktop/Thesis/HCV_Simulation_Study/perm_test.R")
source("/Users/carlyb/Desktop/Thesis/HCV_Simulation_Study/gen_sim_data_simple.R")
# library(lme4)
# library(gee)
# library(geepack)
library(weights)

#Updated 2018-06-12 becase data generation updated to be a matrix, not dataframe
#Also updated because RW and I decided to do 2-sided p-values, not 1-sided as previously coded
#Updated 2018-06-14 to use summary data, and weighted t-tests and permutation tests

#library(tidyverse)

# nc <- 11
# cs <- c(643,271,208,159,141,189,220,246,236,115,204)
# m_u <- log(.55/.45)
# s <- sqrt(0.01011)
# te <- log(1.2)
# pt <- 0.5
# test.mat <- generate_one_meas(nc,cs,m_u,s,te,pt)
# 
# clust.summary <- as.matrix(aggregate(y~ClusterID + x,data = test.mat, FUN = function(x) c(mn = mean(x), variance = abs(var(x)))))
# 
# clust.summary[,3][clust.summary[,2]==0]
# clust.summary[,3][clust.summary[,2]==1]
# wtt <- wtd.t.test(clust.summary[,3][clust.summary[,2]==1] , clust.summary[,3][clust.summary[,2]==0], weight = (1/clust.summary[,4][clust.summary[,2]==1]), weighty = (1/clust.summary[,4][clust.summary[,2]==0]),samedata = FALSE)
# 
# as.numeric(wtt$coefficients[3])

power_fn_simple <- function(nreps,Nclust,csizes,u,sj,trteff,ptrt,nperms){
  
  #set up counters to collect pvals
  c_ttest <- 0
  c_wtd_ttest <- 0
  c_glmm <- 0
  c_gee <- 0
  p_permtest <- 1:nreps #function written to have 1-sided pval
  
  #generate "nreps" datasets, carry out analyses, collect pvals in vectors
  for (z in 1:nreps){
    dat <- generate_one_meas(N.clust = Nclust, clust.sizes = csizes, mu = u, tau.sd = sj, TE = trteff, p.trt = ptrt)
    df.dat <- as.data.frame(dat)
    cl.summary <- as.matrix(aggregate(y~ClusterID + x,data = dat, FUN = function(z) c(mn = mean(z), variance = abs(var(z)))))

    if (t.test(df.dat$ y~ df.dat$x)$p.value < 0.05){
      c_ttest <- c_ttest + 1
    }
    
    if (as.numeric(wtd.t.test(x = cl.summary[,3][cl.summary[,2]==1] , y = cl.summary[,3][cl.summary[,2]==0], weight = (1/cl.summary[,4][cl.summary[,2]==1]), weighty = (1/cl.summary[,4][cl.summary[,2]==0]), samedata = FALSE)$coefficients[3]) < 0.05){
      c_wtd_ttest <- c_wtd_ttest + 1
    }
    
    if (as.numeric(summary(glmer(y ~ x + (1 | ClusterID), data = df.dat, family = binomial))$coefficients[2,4]) < 0.05){
      c_glmm <- c_glmm + 1
    }
    
    if (as.numeric(summary(geeglm(y ~ x, family=binomial, data = df.dat, id = ClusterID))$coefficients[2,4]) < 0.05 ){
      c_gee <- c_gee + 1
    }
    
    p_permtest[z] <- as.numeric(perm_test(dat, "ClusterID", "y", "x", nperms)[5])
    
  }
  #make vector of counts
  counts.vec <- c(c_ttest,c_wtd_ttest,c_glmm,c_gee)
  power.vec <- 1:5

  power.vec[1:4] <- counts.vec / nreps

  
  #add the permutation test pvals
  power.vec[5] <- mean(as.numeric((p_permtest < 0.05)))
  
  #name the columns of the power vector
  #names(power.vec) <- c("t-test","weighted t-test","GLMM","GEE","Permutation test")
  return(power.vec)
}

 #null_test <- power_fn_simple(10,nc,cs,m_u,s,te,pt,1000)
 #weak_test <- power_fn_simple(10,nc,cs,m_u,s,log(1.2),pt,1000)
 #strong_test <- power_fn_simple(10,nc,cs,m_u,s,log(1.5),pt,1000)

####OLD CODE

# power_fn_simple <- function(nreps,Nclust,csizes,u,sj,trteff,ptrt,nperms){
#   
#   #set up vectors to collect pvals
#   p_ttest <- rep(0,nreps)
#   p_glmm <- rep(0,nreps)
#   p_gee <- rep(0,nreps)
#   p_permtest <- rep(0,nreps)
#   
#   #generate "nreps" datasets, carry out analyses, collect pvals in vectors
#   for (z in 1:nreps){
#     dat <- generate_one_meas(N.clust = Nclust, clust.sizes = csizes, mu = u, tau.sd = sj, TE = trteff, p.trt = ptrt)
#     p_ttest[z] <- t.test(dat$y ~ dat$x)$p.value
#     p_glmm[z] <- as.numeric(summary(glmer(y ~ x + (1 | ClusterID), data = dat, family = binomial))$coefficients[2,4])
#     p_gee[z] <- as.numeric(summary(geeglm(y ~ x, family=binomial, data = dat, id = ClusterID))$coefficients[2,4])
#     p_permtest[z] <- as.numeric(perm_test(dat, "ClusterID", "y", "x", nperms)[6])
#   }
#   #make data frame of pvals
#   p_vals.dat <- data.frame(p_ttest,p_glmm,p_gee,p_permtest)
#   
#   #set up vector to collect power
#   power.vec <- rep(0,ncol(p_vals.dat))
#   
#   #compute the pctage of time that a significant result is acheived
#   for (yy in 1:ncol(p_vals.dat)){
#     power.vec[yy] <- mean(as.numeric((p_vals.dat[,yy] < 0.05)))
#   }
# 
#   #name the columns of the power vector
#   #names(power.vec) <- c("t-test","GLMM","GEE","Permutation test")
#   return(power.vec)
# }

#test1 <- power_fn_simple(100,10,rep(100,10),0.5,0.25,1.75,0.5,10000)
#test1

#?glmer

#summary(glmer(df[,"y"] ~ df[,"x"] + (1 | df[,"ClusterID"]), family = binomial))

