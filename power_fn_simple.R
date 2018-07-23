# Writing power function for simple case
# Last updated: 2018-07-22
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
source("/Users/carlyb/Desktop/Thesis/HCV_Simulation_Study/weighted_perm_test.R")
source("/Users/carlyb/Desktop/Thesis/HCV_Simulation_Study/getICC.R")
# library(lme4)
# library(gee)
# library(geepack)
library(weights)

#Updated 2018-06-12 becase data generation updated to be a matrix, not dataframe
#Also updated because RW and I decided to do 2-sided p-values, not 1-sided as previously coded
#Updated 2018-06-14 to use summary data, and weighted t-tests and permutation tests
#Updated 2018-07-02 - switched back to dataframe and not matrix
#                   - reorganized order of calculations
#                   - better distinguish individual and summary-level analyses
#                   - change calculation of weights for t-test
#Updated 2018-07-22 - Remove individual t-test
#                   - Add unweighted cluster-level t-test
#                   - Fix weighted cluster-level t-test (ICC by treatment arm)
#                   - Add weighted perm test

# nc <- 11
# cs <- c(643,271,208,159,141,189,220,246,236,115,204)
# m_u <- log(.55/.45)
# s <- sqrt(0.01011)
# te <- log(1.2)
# pt <- 0.5
# test.dat <- generate_one_meas(nc,cs,m_u,s,te,pt)


power_fn_simple <- function(nreps,Nclust,csizes,u,sj,trteff,ptrt,nperms){
  
  #set up counters to collect pvals
  c_glmm <- 0
  c_gee <- 0
  c_ttest <- 0
  c_wtd_ttest <- 0
  p_permtest <- 1:nreps #function written to have 1-sided pval
  c_wtd_perm_test <- 0
  
  #generate "nreps" datasets, carry out analyses, collect pvals in vectors
  for (z in 1:nreps){
    #create individual-level dataset, and cluster-level dataset
    df.dat <- generate_one_meas(N.clust = Nclust, clust.sizes = csizes, mu = u, tau.sd = sj, TE = trteff, p.trt = ptrt)
    cl.summary.dat <- as.data.frame(as.matrix(aggregate(y~ClusterID + x,data = df.dat, FUN = function(z) c(mn = mean(z), variance = abs(var(z)), size = length(z)))))
    
    #calculate ICC for weighted t-test - using function
    #need to know sizes of clusters with x = 0 and x = 1, separate data by treatment arm
    csizes_x0 <- cl.summary.dat$y.size[cl.summary.dat$x == 0]
    csizes_x1 <- cl.summary.dat$y.size[cl.summary.dat$x == 1]
    x0.dat <- cl.summary.dat[which(cl.summary.dat$x == 0),]
    x1.dat <- cl.summary.dat[which(cl.summary.dat$x == 1),]
    
    #compute ICC by treatment arm
    icc_est_x0 <- getICC(x0.dat,csizes_x0,"y.mn")[2]
    icc_est_x1 <- getICC(x1.dat,csizes_x1,"y.mn")[2]
    
    icc_est <- 1:Nclust
    icc_est[cl.summary.dat$x == 0] <- icc_est_x0
    icc_est[cl.summary.dat$x == 1] <- icc_est_x1
    
    #create weights from summary data and ICC
    wts <- cl.summary.dat$y.size / (1 + icc_est*(cl.summary.dat$y.size - 1))
    
    
    #individual-level analyses
      
      #GLMM
    if (as.numeric(summary(glmer(y ~ x + (1 | ClusterID), data = df.dat, family = binomial))$coefficients[2,4]) < 0.05){
      c_glmm <- c_glmm + 1
    }
      #GEE
    if (as.numeric(summary(geeglm(y ~ x, family=binomial, data = df.dat, id = ClusterID))$coefficients[2,4]) < 0.05 ){
      c_gee <- c_gee + 1
    }
      #weighted permutation test
    if (as.numeric(perm.gee(df.dat$y ~ 1, df.dat$ClusterID, df.dat$x, binomial)$pval.perm) < 0.05) {
      c_wtd_perm_test = c_wtd_perm_test + 1
    }
    
    #cluster-level analyses  
      #unweighted t-test
    if (as.numeric(t.test(y.mn ~ x, data = cl.summary.dat)$p.value < 0.05)){
      c_ttest <- c_ttest + 1
    }
      #weighted t-test
    if (as.numeric(wtd.t.test(x = cl.summary.dat$y.mn[cl.summary.dat$x == 1], y = cl.summary.dat$y.mn[cl.summary.dat$x == 0], weight = (wts[cl.summary.dat$x == 1]), weighty = (wts[cl.summary.dat$x == 0]), samedata = FALSE)$coefficients[3]) < 0.05){
      c_wtd_ttest <- c_wtd_ttest + 1
    }
    
    p_permtest[z] <- as.numeric(perm_test(df.dat, "ClusterID", "y", "x", nperms)[5])
    
  }
  #make vector of counts
  counts.vec <- c(c_glmm, c_gee, c_ttest, c_wtd_ttest, c_wtd_perm_test)
  power.vec <- 1:6
  
  power.vec[1:5] <- counts.vec / nreps
  
  
  #add the permutation test pvals
  power.vec[6] <- mean(as.numeric((p_permtest < 0.05)))
  
  #name the columns of the power vector
  #names(power.vec) <- c("GLMM","GEE","tTest","WtdtTest","WtdPermTest","PermTest")
  
  return(power.vec)
}

 nc <- 14
 cs <- c(643,271,208,159,141,189,220,246,236,115,204,1320,812,195)
 m_u <- log(.53/.47)
 s <- sqrt(0.1728)
 pt <- 0.5
 nperms <- 1000
 testdat <- generate_one_meas(nc, cs, m_u, s, log(1.01), pt)
 testsum <- as.data.frame(as.matrix(aggregate(y~ClusterID + x,data = testdat, FUN = function(z) c(mn = mean(z), variance = abs(var(z)), size = length(z)))))
 #calculate ICC for weighted t-test - using function
 #need to know sizes of clusters with x = 0 and x = 1, separate data by treatment arm
 csizesx0 <- testsum$y.size[testsum$x == 0]
 csizesx1 <- testsum$y.size[testsum$x == 1]
 x0dat <- testsum[which(testsum$x == 0),]
 x1dat <- testsum[which(testsum$x == 1),]
 
 #compute ICC by treatment arm
 icc_est_x0 <- getICC(x0dat,csizesx0,"y.mn")[2]
 icc_est_x1 <- getICC(x1.dat,csizes_x1)[2]
 
 null_test <- power_fn_simple(10,nc,cs,m_u,s,log(1.01),pt,1000)
 weak_test <- power_fn_simple(10,nc,cs,m_u,s,log(1.2),pt,1000)
 strong_test <- power_fn_simple(10,nc,cs,m_u,s,log(1.5),pt,1000)

 # 
# 
# df.test <- generate_one_meas(N.clust = Nclust, clust.sizes = csizes, mu = m_u, tau.sd = sj, TE = log(1.5), p.trt = ptrt)
# cl.summary.test <- as.data.frame(as.matrix(aggregate(y~ClusterID + x,data = df.test, FUN = function(z) c(mn = mean(z), variance = abs(var(z)), size = length(z)))))
# t.test(df.test$y ~ df.test$x)
# rho <- suppressWarnings(as.numeric(iccbin(cid = ClusterID, y = y, data = df.test, method = "aov", ci.type = "aov")$estimates[2]))
# summary(glmer(y ~ x + (1 | ClusterID), data = df.test, family = binomial))
# t.test(y.mn ~ x, data = cl.summary.test)

# cl.summary.df <- as.data.frame(cl.summary.test)
# 
# rho <- as.numeric(iccbin(cid = ClusterID, y = y, data = df.test, method = "aov", ci.type = "aov")$estimates[2])
# wts <- cl.summary.df$y.size / (1 + rho*(cl.summary.df$y.size - 1))
# cl.summary.df$y.mn[cl.summary.df$x == 1]
# wts[cl.summary.df$x == 1]
# 
# wtd.t.test(x = cl.summary.df$y.mn[cl.summary.df$x == 1], y = cl.summary.df$y.mn[cl.summary.df$x == 0], weight = (wts[cl.summary.df$x == 1]), weighty = (wts[cl.summary.df$x == 0]), samedata = FALSE)
# system.time(null_test <- power_fn_simple(10,Nclust,csizes,m_u,sj,log(1.02),0.5,1000))

#gee.fit <- geeglm(y ~ x, family=binomial, data = test.df, id = ClusterID)
#summary(gee.fit)

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

