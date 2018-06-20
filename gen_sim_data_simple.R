
# Generating Simulated Data for simple case
# Last updated: 2018-06-15
# Programmer: Carly Broadwell
# Purpose:  Write a function to generate cluster data which will be analyzed under four different approaches. Potential 
#           parameters to be varied: 
#           I (number of clusters) 
#           n_i (number of observations per cluster, not necessarily fixed)
#           ICC - how does this impact simulation?
#           J (number of measurement occasions) - one baseline, and one or two post-baseline
#           Will also want to vary the treatment effect - acheived in varying beta
library(lme4)
library(gee)
library(geepack)

#N.clust is the number of clusters
#clust.sizes is a vector of length N.clust containing the sizes of each cluster
#mu is the odds of the outcome at baseline
#tau.sd is the standard deviation of the log of the cluster-intercepts
#p.trt is the proportion of the clusters assigned to the intervention.
#UPDATED 2018-06-15: TE argument is now the log of the OR (previously was OR) and mu argument is now the log odds at baseline

generate_one_meas = function(N.clust,clust.sizes,mu,tau.sd,TE,p.trt) {
  #calculate total number of obs
  N.total <- sum(clust.sizes)
  
  #allocate treatment by cluster
  x_short.vec <- 1:N.clust
  x_short.vec[1:round(N.clust*p.trt)] <- 0
  x_short.vec[(round(N.clust*p.trt)+1):length(x_short.vec)] <- 1
  
  #generate N.clust random intercepts
  #let sd.tau be a parameter for now - need to work out how it fits with the ICC
  #also make ClusterID variable
  gamma.vec <- rnorm(N.clust,0,tau.sd)
  
  #Clust.ID.vec will be a vector of length N.total identifying the Cluster ID for each individual (same for all ind. in a cluster)
  #multi_gamma.vec w/b a vector of length N.total identifying the random intercept for each individual (same for all ind. in a cluster)
  #x.vec w/b a vector of length N.total identifying the treatment assignment for each individual (same for all ind in a cluster)
  
  #create empty vectors 
  multi_gamma.vec <- 1:N.total
  Clust.ID.vec <- 1:N.total
  x.vec <- 1:N.total
  
  #initialize count to keep track of size
  count <- 0 
  
  #assign values to individuals by cluster
  for (i in 1:N.clust){
    this.N <- clust.sizes[i]
    multi_gamma.vec[(count + 1):(count + this.N)] <- gamma.vec[i]
    Clust.ID.vec[(count + 1):(count + this.N)] <- i
    x.vec[(count + 1):(count + this.N)] <- x_short.vec[i]
    count <- count + this.N
  }
  
  #generate vector of probabilities
  p.vec <- plogis(mu + multi_gamma.vec + TE*x.vec)
  
  #use as probability of Y=1
  y.vec <- rbinom(N.total,1,p.vec)
  
  #bind all vectors into a dataframe #UPDATED 2018-06-12 to be a matrix for EFFICIENCY GAINS
  sim.data <- as.matrix(cbind(Clust.ID.vec,x.vec,y.vec))
  colnames(sim.data) <- c("ClusterID","x","y")
  
  return(sim.data)
}

#null_test <- generate_one_meas(11,as.integer(abs(rnorm(11,200,5))),log(0.55/0.45),sqrt(0.06),log(1),0.5)
# weak_test <- generate_one_meas(11,as.integer(abs(rnorm(11,200,5))),log(0.55/0.45),sqrt(0.06),log(1.2),0.5)
# strong_test <- generate_one_meas(11,as.integer(abs(rnorm(11,200,5))),log(0.55/0.45),sqrt(0.06),log(1.5),0.5)
# 
# null_test.df <- as.data.frame(null_test)
# summary(glmer(y ~ x + (1 | ClusterID), data = null_test.df, family = binomial))$coefficients[2,4]
# summary(geeglm(y ~ x, family = binomial, data = null_test.df, id = ClusterID))$coefficients[2,4]
# 
# weak_test.df <- as.data.frame(weak_test)
# summary(glmer(y ~ x + (1 | ClusterID), data = weak_test.df, family = binomial))$coefficients[2,4]
# summary(geeglm(y ~ x, family = binomial, data = weak_test.df, id = ClusterID))$coefficients[2,4]
# 
# strong_test.df <- as.data.frame(strong_test)
# summary(glmer(y ~ x + (1 | ClusterID), data = strong_test.df, family = binomial))$coefficients[2,4]
# summary(geeglm(y ~ x, family = binomial, data = strong_test.df, id = ClusterID))$coefficients[2,4]

#Nclust=12
#csize = as.integer(rnorm(Nclust,1500,30))
#test <- generate_one_meas(Nclust,csize,0.5,0.06,1.5,0.5)
#model <- glmer(y ~ x + (1 | ClusterID), data=test, family = binomial)

#Nclust = 100
#csize = as.integer(rnorm(100,100,6))
#test1 <- generate_one_meas(N.clust=Nclust,clust.sizes=csize,mu=2,tau.sd=log(2),TE=1.5,p.trt=0.5)
#model <- glmer(y ~ x + (1 | ClusterID), data=test1, family = binomial)
#exp(as.numeric(fixef(model)[2])) 


# #t-test 
# ttest <- t.test(test1$y~test1$x)
# p_ttest <- ttest$p.value
# 
# #glmm with random cluster intercept
# glmm_fit <- glmer(y ~ x + (1 | ClusterID), data=test1, family = binomial)
# p_glmm <- as.numeric(summary(glmm_fit)$coefficients[2,4])
# 
# #gee model
# gee_fit <- geeglm(y ~ x, family = binomial, data = test1, id=ClusterID)
# p_gee <- summary(gee_fit)$coefficients[2,4]
# 
# #permutation test - code by hand
# #smaller functions:   pre_perm_test to get cluster-level summaries
# #                     one.perm to do one permutation
# 
# pre_perm_test <- function(dat,groupvar,yvar,xvar){
#   Nclust <- length(unique(dat[[groupvar]]))
#   clustmeans <- rep(0,Nclust)
#   clusttrt <- rep(0,Nclust)
#   for (j in 1:Nclust){
#     clustdata <- dat[ which(dat[[groupvar]]== j),]
#     clustmeans[j] <- mean(clustdata[[yvar]])
#     clusttrt[j] <- clustdata[[xvar]][1]
#   }
#   return(cbind(clustmeans,clusttrt))
# }
# one.perm <- function(x,y){
#   x_perm <- sample(x)
#   return(mean(y[x_perm==0])-mean(y[x_perm==1]))
# }
# perm_test <-function(dat,groupvar,yvar,xvar,nperms){
#   pre_perm_test(dat,groupvar,yvar,xvar)
#   clustmeans <- pre_perm_test(test1,"ClusterID","y","x")[,1]
#   clusttrt <- pre_perm_test(test1,"ClusterID","y","x")[,2]
#   #calculate observed value
#   t_obs <- mean(clustmeans[clusttrt==0]) - mean(clustmeans[clusttrt==1])
#   #permute trt allocation nperms times
#   mult.perms <- replicate(nperms,one.perm(clusttrt,clustmeans))
#   #calculate 2sided pval
#   pval.2sided <- mean(abs(mult.perms)>abs(t_obs))
#   #calculate 1sided pval
#   pval.1sided <- mean((t_obs > mult.perms))
#   result <- list(groupvar,yvar,xvar,nperms,pval.2sided,pval.1sided)
#   names(result) <- c("Grouping Variable","Y Variable","X Variable","Number of Permutations","2-sided p-value","1-sided p-value")
#   return(result)
# }
# permtest <- perm_test(test1,"ClusterID","y","x",10000)
# p_permtest <- as.numeric(permtest[6])
# 
# #structure all together
# pval.vec <- c(p_ttest,p_glmm,p_gee,p_permtest)
# pval.vec
