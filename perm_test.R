# Coding for a permutation test (unmatched)
# Last updated: 2018-06-12
# Programmer: Carly Broadwell
# Purpose: write code to conduct a permutation test 
# Parameters: 
#   dat = dataset
#   groupvar = grouping variable (e.g., Cluster ID) with which to subset the data and calculate group-level summaries from individual data
#   yvar = response variable. held fixed under permutation test
#   xvar = exposure/treatment assignment variable. this is what will be permuted.
#   null hypothesis is that we can permute the exposure randomly and have similar results in the difference in proportions/means in the untreated and treated as we did in the observed data 


#Updated 2018-06-12 because data generation updated to result in matrix, instead of data frame, for efficiency.
pre_perm_test <- function(dat,groupvar,yvar,xvar){
  Nclust <- length(unique(dat[,groupvar]))
  clustmeans <- rep(0,Nclust)
  clusttrt <- rep(0,Nclust)
  for (j in 1:Nclust){
    clustdata <- dat[ which(dat[,groupvar]== j),]
    clustmeans[j] <- mean(clustdata[,yvar])
    clusttrt[j] <- clustdata[1,xvar]
  }
  return(cbind(clustmeans,clusttrt))
}
one.perm <- function(x,y){
  x_perm <- sample(x)
  return(mean(y[x_perm==0])-mean(y[x_perm==1]))
}
perm_test <-function(dat,groupvar,yvar,xvar,nperms){
  pre_perm_test(dat,groupvar,yvar,xvar)
  clustmeans <- pre_perm_test(dat,"ClusterID","y","x")[,1]
  clusttrt <- pre_perm_test(dat,"ClusterID","y","x")[,2]
  #calculate observed value
  t_obs <- mean(clustmeans[clusttrt==0]) - mean(clustmeans[clusttrt==1])
  #permute trt allocation nperms times
  mult.perms <- replicate(nperms,one.perm(clusttrt,clustmeans))
  #calculate 2sided pval
  pval.2sided <- mean(abs(mult.perms)>abs(t_obs))
  #calculate 1sided pval
  pval.1sided <- mean((t_obs > mult.perms))
  result <- list(groupvar,yvar,xvar,nperms,pval.2sided,pval.1sided)
  names(result) <- c("Grouping Variable","Y Variable","X Variable","Number of Permutations","2-sided p-value","1-sided p-value")
  return(result)
}

