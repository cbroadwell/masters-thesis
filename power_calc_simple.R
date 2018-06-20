#source("/Users/carlyb/Desktop/Thesis/HCV_Simulation_Study/gen_sim_data_simple.R")
#source("/Users/carlyb/Desktop/Thesis/HCV_Simulation_Study/perm_test.R")
source("/Users/carlyb/Desktop/Thesis/HCV_Simulation_Study/power_fn_simple.R")

library(ggplot2)
library(reshape2)

power_plot <- function(nreps,Nclust,csizes,u,sj,TE.vec,ptrt,nperms){
  result <- matrix(0,length(TE.vec),6)
  for (k in 1:length(TE.vec)){
    power.vec <- power_fn_simple(nreps,Nclust,csizes,u,sj,TE.vec[k],ptrt,nperms)
    result[k,1] <- TE.vec[k]
    result[k,2:6] <- power.vec
  }
  result.dat <- data.frame(result)
  names(result.dat) <- c("Treatment Effect","t-test","Weighted t-test","GLMM","GEE","Permutation test")
  return(result.dat)
}

#######
# Power calc conditions:
#   11 clusters (England only)
#   mu = log(.54/.46)
#   tau = sqrt(0.01011)
#   1:1 trt allocation
#   1 post-baseline measurement
#   No additional covariates
#   Unmatched design
set.seed(061018)

nreps <- 500
Nclust <- 11
csizes <- c(643,271,208,159,141,189,220,246,236,115,204)
m_u <- log(.54/.46)
sj <- sqrt(0.01011)
TE.vec <- log(seq(1,1.5,.01))
ptrt <- 0.5
nperms <- 1000

for_plot_a1 <- power_plot(nreps,Nclust,csizes,m_u,sj,TE.vec,ptrt,nperms)

#write.csv(for_plot_a1,file="/Users/carlyb/Desktop/Thesis/HCV_Simulation_Study/power_calc_a1.csv")
#save(for_plot_a1,file="/Users/carlyb/Desktop/Thesis/HCV_Simulation_study/power_calc_a1.Rdata")

for_plot_long_a1 <- melt(for_plot_a1,id="Treatment Effect")
p1 <- ggplot(for_plot_long_a1, aes(x=for_plot_long_a1[["Treatment Effect"]],y=value,color=variable))
p2 <- p1 + geom_line()

avgcsize <- round(mean(csizes),0)
label_txt1 <- paste0("mu",":",round(m_u,2))
label_txt2 <- paste0("tau",":",round(sj,2))

p3 <- p2 + xlab("Treatment Effect (log scale)") + ylab("Power") + ggtitle("Power",subtitle = paste0(nreps," reps,",Nclust," clusters, avg cluster size  ", avgcsize))
p4 <- p3 + geom_text(aes(x=0.01,y=0.65,label=label_txt1),parse = T,color="black")
p5 <- p4 + geom_text(aes(x=0.01,y=0.6,label=label_txt2),parse = T,color="black")
p6 <- p5 + guides(color=guide_legend(title="Test"))
#png('power_curve_a1.png')
#p6
#dev.off()


#######
# Power calc conditions:
#   14 clusters (England and Scotland)
#   mu = log(.53/.47)
#   tau = sqrt(0.01)
#   1:1 trt allocation
#   1 post-baseline measurement
#   No additional covariates
#   Unmatched design

Nclust <- 14
csizes <- c(643,271,208,159,141,189,220,246,236,115,204,1320,812,195)
m_u <- log(.53/.47)
sj <- sqrt(0.01)

for_plot_a2 <- power_plot(nreps,Nclust,csizes,m_u,sj,TE.vec,ptrt,nperms)

#write.csv(for_plot_a2,file="/Users/carlyb/Desktop/Thesis/HCV_Simulation_Study/power_calc_a2.csv")
#save(for_plot_a2,file="/Users/carlyb/Desktop/Thesis/HCV_Simulation_study/power_calc_a2.Rdata")

for_plot_long_a2 <- melt(for_plot_a2,id="Treatment Effect")
p7 <- ggplot(for_plot_long_a2, aes(x=for_plot_long_a2[["Treatment Effect"]],y=value,color=variable))
p8 <- p7 + geom_line()

avgcsize <- round(mean(csizes),0)
label_txt1 <- paste0("mu",":",round(u,2))
label_txt2 <- paste0("tau",":",round(sj,2))

p9 <- p8 + xlab("Treatment Effect (log scale)") + ylab("Power") + ggtitle("Power",subtitle = paste0(nreps," reps,",Nclust," clusters, avg cluster size  ", avgcsize))
p10 <- p9 + geom_text(aes(x=0.9,y=0.65,label=label_txt1),parse = T,color="black")
p11 <- p10 + geom_text(aes(x=0.9,y=0.6,label=label_txt2),parse = T,color="black")
p12 <- p11 + guides(color=guide_legend(title="Test"))
p12
