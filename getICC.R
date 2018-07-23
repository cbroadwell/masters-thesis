## ICC estimation: 
#     E(s^2) = pi(1-pi)/(harmonic mean of cluster sizes) + sigma^2between
#     so sigma^2B = s^2 - p(1-p)/mH

## A function - modified to handle individual-level data

getICC <- function(dat,c.sizes,yvar){
  HarMean <- 1/mean(1/c.sizes)
  S.est <- var(dat[[yvar]])
  p.est <- mean(dat[[yvar]])
  Sb.est <- S.est - p.est*(1-p.est)/HarMean
  ICC.est <- Sb.est / p.est*(1-p.est)
  return(c(Sb.est,ICC.est))
}


### Helper function for Bruan/Feng permutation code
get.rho.aov <- function(my.data, adj=F)
{
  #Computes intra-cluster correlation for binary data using ANOVA method.
  ni <- table(my.data$cluster) 
  N <- length(my.data$frtveg)
  k <- max(cluster)
  n0 <- 1/(k-1)*(N - 1/N*sum(ni^2))
  
  if (adj==T) 
  { 
    mod <- anova(lm(frtveg~white+male+age+educ+as.factor(cluster), data=my.data))
    ms.b <- mod[5,3]
    ms.w <- mod[6,3]
  }
  
  if (adj==F)
  {
    mod <- anova(lm(frtveg~as.factor(cluster), data=my.data))
    ms.b <- mod[1,3]
    ms.w <- mod[2,3]
  }
  rho.aov <- max((ms.b-ms.w)/(ms.b+(n0-1)*ms.w),0)
  rho.aov
}
