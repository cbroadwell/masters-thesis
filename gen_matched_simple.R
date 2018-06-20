# Generating matched cluster data for simple case
# Last updated: 2018-06-10
# Programmer: Carly Broadwell

## A note: at some point it might be useful to be able to match "as closely as possible" but not perfectly
## i.e., pair clusters that are generally close in size but not the exact same size


# N.clust is the total number of clusters
# clust.sizes is a vector of length N.clust containing the cluster sizes
# mu is estimated log odds of outcome at baseline 
# tau.sd is estimated between-cluster standard deviation
# TE is log of treatment effect
# t_alloc is treatment allocation ratio. If t_alloc = 1, then each experimental cluster has one 
#   matched control cluster. If t_alloc = 2, then each experimental cluster will have two.... etc.
library(reshape2)

gen_matched_simple = function(N.clust,clust.sizes,mu,tau.sd,TE,t_alloc){
  
  #create matching variable based on size
  clust.sizes <- sort(clust.sizes)
  m_id.vec <- vector('numeric',length(clust.sizes))
  m_id.vec[1:length(c.size)] <- 1:(length(clust.sizes)/(t_alloc + 1))
  m_id.vec <- sort(m_id.vec)
  
  #create treatment assignment vector at cluster level
  x_short.vec <- vector('numeric',length(clust.sizes))
  x_short.vec[1:length(clust.sizes)] <- c(0,1)
  
  #create gammas based on between cluster variance at cluster level
  gamma.vec <- rnorm(length(clust.sizes),0,tau.sd)
  
  #calculate total number of obs in data
  N.total <- sum(clust.sizes)
  
  #initialized long (individual) vectors
  multi_gamma.vec <- vector('numeric',N.total)
  ClustID.vec <- vector('numeric',N.total)
  MatchID.vec <- vector('numeric',N.total)
  x.vec <- vector('numeric',N.total)
  
  #replace with values within each cluster
  count <- 0
  for (i in 1:N.clust){
    this.N <- clust.sizes[i]
    multi_gamma.vec[(count + 1):(count + this.N)] <- gamma.vec[i]
    ClustID.vec[(count + 1):(count + this.N)] <- i
    x.vec[(count + 1):(count + this.N)] <- x_short.vec[i]
    MatchID.vec[(count + 1):(count + this.N)] <- m_id.vec[i]
    count <- count + this.N
  }

  p.vec <- plogis(mu + multi_gamma.vec + TE*x.vec)
  y.vec <- rbinom(N.total,1,p.vec)
  
  result.mat <- matrix(c(ClustID.vec, MatchID.vec, x.vec, y.vec), nrow = N.total, ncol = 4)
  colnames(result.mat) <- c("ClusterID","MatchID","x","y")
  
  sum.result.mat <- aggregate(y ~ ClusterID + x + MatchID, data = result.mat, FUN = function(z) c(mean(z), var(z)))
  return(list(result.mat, sum.result.mat))
  
}

nc <- 14
c.size <- c(643,271,208,159,141,189,220,246,236,115,204,1320,812,195)
u <- log(.55/.45)
s<- sqrt(0.01011)
te <- log(1.5)
ta <- 1
# 
out_test <- gen_matched_simple(nc,c.size,u,s,te,ta)
head(out_test[[1]])
head(out_test[[2]])

summary.dat<-reshape(data=out_test[[2]],varying = NULL,v.names = c("y"), timevar = "x", idvar = "MatchID", drop = "ClusterID", direction = "wide")
summary.dat


t.test(summary.dat$y.1, summary.dat$y.0, paired = TRUE)

