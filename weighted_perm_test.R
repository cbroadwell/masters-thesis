#weighted permutation test code
#https://cdsweb01.fhcrc.org/statmethods/documents/permutation_test.htm
#Tom Braun at Biostatistics Department, U. of Michigan
#Ziding Feng at the Fred Hutchinson Cancer Center.

#Reference:
  
  #Optimal Permutation Tests for the Analysisi of Group Randomized Trials

#Braun T and Feng Z.

#JASA 2001 96: 1124-1132

perm.gee <- function(formula, clus, rand, fam, rho)
{
  #This function applies GEE to a set of clustered data and then computes the
  #permutation distribution of the EE.
  
  #formula = model formula, i.e. y ~ age + gender
  #fam = likelihood for outcomes (use either gaussian, poisson, or binomial)
  #clus = variable defining level of clustering
  #rand = variable indicating treatment arm
  #If outcome is normal, need estimate of scale parameter.
  mod <- glm(formula, family = fam)
  p <- length(mod$coef)
  scale <- mod$deviance/mod$df.residual
  
  #Define all necessary parts of score which vary based on "family"
  varfn <- fam()$variance	#Variance function
  invlink <- fam()$linkinv	#Inverse link function
  fam.ind <- match(fam()$family, c("gaussian", "poisson", "binomial"))
  
  #Compute pieces of score function UNDER NULL HYPOTHESIS
  X <- model.matrix(formula)
  b <- mod$coef
  eta <- X[,-p] %*% as.matrix(b[-p])
  mu <- invlink(eta)
  deta.dmu <- switch(fam.ind, 1, 1/mu, 1/(mu*(1-mu)))
  if (length(deta.dmu)==1) deta.dmu <- rep(deta.dmu,length(mu))
  y <- model.extract(model.frame(formula), response)
  e <- as.numeric(y - mu)
  v.mu <- varfn(mu)
  if (length(v.mu)==1) v.mu <- rep(v.mu,length(mu))
  
  #Estimate the intra-cluster correlation based on the residuals above. 
  #CBB ELIMINATING - computed previously in function, no need to redo
  #Inserting computed value as variable "rho"
  # if (p==2) rho <- get.rho.aov(dat, clust, yvar,)
  # if (p!=2) rho <- get.rho.aov(the.data,T)
  
  #Compute scores for each cluster without regard to treatment assignment
  cat("Computing scores:	1...")
  scores <- NULL
  
  for(i in unique(clus)) 
  {
    cat(i)
    j <- clus == i
    m <- sum(j)
    if (p==1)
    {
      a <- -rho/(1+(m-1)*rho)
      V.inv <- matrix(a, nrow=m, ncol=m)
      diag(V.inv) <- 1+a
      V.inv <- V.inv / (1-rho) / v.mu[j]
    }
    else
    {
      v1 <- matrix(rho, nrow=m, ncol=m)
      diag(v1) <- 1
      v2 <- sqrt(as.matrix(v.mu[j]) %*% 
                   t(as.matrix(v.mu[j])))
      V <- v1 * v2 * scale
      V.inv <- solve(V)
    }
    
    U <- t(1/deta.dmu[j]) %*% V.inv %*% e[j]
    
    scores <- cbind(scores, U)
    #	if(i==max(clus)) cat("\n")
  }
  
  scores1 <- scores
  #Compute scores again, but use cluster sizes as weights.
  cat("3...")
  scores2 <- NULL
  for(i in unique(clus)) 
  {
    j <- clus == i
    m <- sum(j)
    U <- sum(e[j])
    
    scores2 <- cbind(scores2, U)
  }
  #Compute scores one more time, but use no weights
  cat("4...","\n")
  scores3 <- NULL
  for(i in unique(clus)) 
  {
    j <- clus == i
    m <- sum(j)
    U <- sum(e[j]) / m
    
    scores3 <- cbind(scores3, U)
  }
  
  k <- length(unique(clus))/2	#Number of pairs
  
  #Create matrix of all 2^k possible treatment permutations if number of
  #clusters is small enough
  if (k <=13)
  { 
    cat("Computing perm distns...\n")
    perm.trt <- allperms(k)	#Permutations for treated clusters
    perm.trt <- perm.trt*2 - 1	#Change (T, F) to (+1, -1)
    
    perm.ctrl <- -perm.trt	#Permutations for untreated clusters
    
    #Now combine permutations pairwise into one matrix
    perms <- NULL
    for (i in 1:k) perms <- rbind(perms, perm.trt[i,], perm.ctrl[i,])
  }
  
  #If there are too many clusters, do a subsample of the permutations
  if (k > 13)
  {
    cat("Computing (approx) perm distns...\n")
    perms <- get.subset(2000, k)
  }
  
  #Compute the OBSERVED sum-of-scores, with and without the correct weights.
  trt.obs <- as.matrix(tapply(rand,clus,min))
  obs.ss <- as.numeric(scores %*% trt.obs)
  obs.ss1 <- as.numeric(scores1 %*% trt.obs)
  obs.ss2 <- as.numeric(scores2 %*% trt.obs)
  obs.ss3 <- as.numeric(scores3 %*% trt.obs)
  
  #Compute permutation distribution of obs.ss and obs.ss2
  perm.ss <- sort(scores %*% perms)
  perm.ss1 <- sort(scores1 %*% perms)
  perm.ss2 <- sort(scores2 %*% perms)
  perm.ss3 <- sort(scores3 %*% perms)
  
  #Compute p-values
  pval <- sum(perm.ss >= obs.ss)/(ncol(perms))
  pval1 <- sum(perm.ss1 >= obs.ss1)/(ncol(perms))
  pval2 <- sum(perm.ss2 >= obs.ss2)/(ncol(perms))
  pval3 <- sum(perm.ss3 >= obs.ss3)/(ncol(perms))
  
  
  list(pval.perm = pval, pval.perm.exact = pval1, pval.perm.samp = pval2,
       pval.perm.indep = pval3)
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

allperms <- function(n)
{
  #Creates matrix of dim n x 2^n, each column of which is one of the 2^n
  #permutations of a vector of 0's and 1`s.
  
  k <- 2^n
  powers <- 2^(1:n)
  all <- sapply(0:(k - 1), function(i, powers)
    i %% powers < powers/2, powers = powers)
  all
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

get.subset <- function(N, k)
{
  #Will find a subset of N treatment allocations for k pairs of clusters.
  
  perms <- matrix(nrow=2*k, ncol=N)
  
  for (i in 1:k)
  {
    perms[2*i-1,] <- sample(c(-1,1), replace=T, size=N) 
    perms[2*i,] <- -perms[2*i-1,]
  }
  perms
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
