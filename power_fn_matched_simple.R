# Writing analysis/power function for matched pairs data in simple case
# Last updated: 2018-06-10
# Programmer: Carly Broadwell

source("/Users/carlyb/Desktop/Thesis/HCV_Simulation_Study/gen_matched_simple.R")
nc <- 12
c.size <- abs(as.integer(rnorm(6,100,5)))
u <- 0.5
s<- sqrt(0.06)
te <- 1
ta <- 1
test.dat <- gen_matched_simple(nc,c.size,u,s,te,ta)

#Matched t-test - based on cluster summaries
c_summaries <- aggregate(test.dat$y,by = list(test.dat$ClusterID), FUN = mean)
str(c_summaries)

t.test(test.dat$y ~ test.dat$x, paired = TRUE)
