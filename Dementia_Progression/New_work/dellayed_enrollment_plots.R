library(rstan)
library(bayesplot)
library(ggplot2)
library(coda)
library(circular)
library(here)
library(expm) 
library(msm)
library(fGarch)
library(splines2)
library(splines)
library(gridExtra)
library(ggpubr)


# Load no bias outputs
load(file="~/Documents/Research Comp/Own Work/NoBiasOutput/full_msm.rda")
load(file="~/Documents/Research Comp/Own Work/NoBiasOutput/bayes_models.rda")

### Intercepts ###

# Compute coverages
covs <- matrix(0,22,50)

# betaMat is used to generate the data in Simulate_cav.r.
trueBetaMat <- matrix(c(-2.54,  0.11, -0.56,
                        -2.94, -0.24,  0.15,
                        -1.10, -0.15, -0.03,
                        -3.92,  0.23,  0.21,
                        -2.12,  0.08,  1.17), nrow=5, byrow=TRUE)


cmap <- matrix(c( 1, 4, 7,10,13,  16,17,18,19,20,21,
                  2, 5, 8,11,14,   0, 0, 0, 0, 0, 0,
                  3, 6, 9,12,15,   0, 0, 0, 0, 0, 0), 3, 11, byrow=TRUE)

index <- cmap[1,1:5] # Indices from cmap.
truth <- t(trueBetaMat)[index]

# Bayes
for(i in 1:50){
  
  covs[, i] <- ifelse(
    truth >= summary(models[[i]])$summary[1:22 , 4] & 
      truth <= summary(models[[i]])$summary[1:22 , 8], 1, 0)
  
  b_0s <- rstan::extract(models[[i]], pars=c("b_0"))
  b_1s <- rstan::extract(models[[i]], pars=c("b_1"))
  # Compute the intercepts for un-centered data
  for(j in 1:5){
    covs[7+j , i] <- ifelse(
      truth[7+j] >= quantile(b_0s[[1]][,j] - 8*b_1s[[1]][,j], 0.025) &
        truth[7+j] <= quantile(b_0s[[1]][,j] - 8*b_1s[[1]][,j], 0.975), 
      1, 0)
  }
}

# msm
covg_msm <- c(0,0,0,0,0)
for(l in 1:5){
  covg_msm[l] <- round( (sum( lower_msm[,index[l]] - 8*trueBetaMat[l,2]  <= truth[l] & 
                                truth[l] <= upper_msm[,index[l]] - 8*trueBetaMat[l,2]))/ 50, 3)
}

labels <- c('State 1 (well)   --->   State 2 (mild)',
            'State 1 (well)   --->   State 4 (dead)',
            'State 2 (mild)   --->   State 3 (severe)',
            'State 2 (mild)   --->   State 4 (dead)',
            'State 3 (severe)   --->   State 4 (dead)')

# Box plots
pdf(paste0("./New_work/",'NB_cav_intercept.pdf')) 
for(l in 1:5){
  par( mar=c(5, 2.5, 2.5, 2.5))
  b_cov <- round(mean(covs[7+l,1:61]), 2)
  intercepts_msm <- full_msm[[1]][7+l, ]
  
  boxplot(output_bayes[7+l, ], intercepts_msm, 
          names=c('Bayes','MLE (msm)'),
          ylab=NA, 
          main=bquote(hat(beta)[0]^{(.(l))} ~ (.(labels[l]))), 
          cex =1.2)
  abline(h=truth[l],col='green',lwd=4,lty=5)
  mtext(paste0('.95 coverage = ',toString(b_cov),' Bayes, ',toString(covg_msm[l]),' MLE (msm)'), 
        side=1, line=2.5, cex=1.2)
}
dev.off()


### Beta 1 ###

# Compute coverages
index <- cmap[2,1:5]
truth <- trueBetaMat[,2]

# Bayes
for(i in 1:50){
  b_0s <- rstan::extract(models[[i]], pars=c("b_0"))
  b_1s <- rstan::extract(models[[i]], pars=c("b_1"))
  for(j in 1:5){
    covs[12+j , i] <- ifelse(truth[j] >= quantile(b_1s[[1]][,j], 0.025) &
                              truth[j] <= quantile(b_1s[[1]][,j], 0.975), 1, 0)
  }
}

# msm
covg_msm <- c(0,0,0,0,0)
for(l in 1:5){
  covg_msm[l] <- round( (sum( full_msm[[2]][, index[l]]  <= truth[l] & 
                                truth[l] <= full_msm[[3]][, index[l]]))/ 50, 3)
}


# Box plots
pdf(paste0("./New_work/",'NB_cav_b1.pdf')) 
for(l in 1:5){
  par( mar=c(5, 2.5, 2.5, 2.5))
  b_cov <- round(mean(covs[12+l,1:50]), 2)
  intercepts_msm <- full_msm[[1]][12+l, ]
  
  boxplot(output_bayes[12+l, ], intercepts_msm, 
          names=c('Bayes','MLE (msm)'),
          ylab=NA, 
          main=bquote(hat(beta)[1]^{(.(l))} ~ (.(labels[l]))), 
          cex =1.2)
  abline(h=truth[l],col='green',lwd=4,lty=5)
  mtext(paste0('.95 coverage = ',toString(b_cov),' Bayes, ',toString(covg_msm[l]),' MLE (msm)'), 
        side=1, line=2.5, cex=1.2)
}
dev.off()


### Beta 2 ###

# Compute coverages
index <- cmap[3,1:5]
truth <- trueBetaMat[,3]

# Bayes
for(i in 1:50){
  b_2s <- rstan::extract(models[[i]], pars=c("b_2"))
  for(j in 1:5){
    covs[17+j , i] <- ifelse(truth[j] >= quantile(b_2s[[1]][,j], 0.025) &
                               truth[j] <= quantile(b_2s[[1]][,j], 0.975), 1, 0)
  }
}

# msm
covg_msm <- c(0,0,0,0,0)
for(l in 1:5){
  covg_msm[l] <- round( (sum( full_msm[[2]][, index[l]]  <= truth[l] & 
                                truth[l] <= full_msm[[3]][, index[l]]))/ 50, 3)
}


# Box plots
pdf(paste0("./New_work/",'NB_cav_b2.pdf')) 
for(l in 1:5){
  par( mar=c(5, 2.5, 2.5, 2.5))
  b_cov <- round(mean(covs[17+l,1:50]), 2)
  intercepts_msm <- full_msm[[1]][17+l, ]
  
  boxplot(output_bayes[17+l, ], intercepts_msm, 
          names=c('Bayes','MLE (msm)'),
          ylab=NA, 
          main=bquote(hat(beta)[2]^{(.(l))} ~ (.(labels[l]))), 
          cex =1.2)
  abline(h=truth[l],col='green',lwd=4,lty=5)
  mtext(paste0('.95 coverage = ',toString(b_cov),' Bayes, ',toString(covg_msm[l]),' MLE (msm)'), 
        side=1, line=2.5, cex=1.2)
}
dev.off()



### Miscclassification probabilities ###

labels <- c('P( observed state 2 | true state 1 )',
            'P( observed state 1 | true state 2 )',
            'P( observed state 3 | true state 2 )',
            'P( observed state 2 | true state 3 )')

# Compute coverages
index <- cmap[1,6:9]
truth <- c( .01, .24, .06, .11)

# Bayes
for(i in 1:50){
  ps <- rstan::extract(models[[i]], pars=c("p"))
  for(j in 1:4){
    covs[j , i] <- ifelse(truth[j] >= quantile(ps$p[,j], 0.025) &
                               truth[j] <= quantile(ps$p[,j], 0.975), 1, 0)
  }
}

# msm
covg_msm <- c(0,0,0,0)
for(l in 1:4){
  covg_msm[l] <- round( (sum( full_msm[[2]][, index[l]]  <= truth[l] & 
                                truth[l] <= full_msm[[3]][, index[l]]))/ 50, 3)
}


# Box plots
pdf(paste0("./New_work/",'NB_cav_p.pdf')) 
par(mfrow=c(2, 2), mar=c(5, 2.5, 2.5, 2.5))
for(l in 1:4){
  par( mar=c(5, 2.5, 2.5, 2.5))
  b_cov <- round(mean(covs[l,1:50]), 2)
  intercepts_msm <- full_msm[[1]][l, ]
  
  boxplot(output_bayes[l, ], intercepts_msm, 
          names=c('Bayes','MLE (msm)'),
          ylab=NA, 
          main=bquote(hat(p)[.(l)] ~ (.(labels[l]))), 
          cex =1.2)
  abline(h=truth[l],col='green',lwd=4,lty=5)
  mtext(paste0('.95 coverage = ',toString(b_cov),' Bayes, ',toString(covg_msm[l]),' MLE (msm)'), 
        side=1, line=2.5, cex=1)
}
dev.off()


# Load bias outputs
load(file="./New_work/Output/full_msm.rda")
load(file="./New_work/Output/models.rda")
load(file="./New_work/Output/Output_mle.rda")

labels <- c('State 1 (well)   --->   State 2 (mild)',
            'State 1 (well)   --->   State 4 (dead)',
            'State 2 (mild)   --->   State 3 (severe)',
            'State 2 (mild)   --->   State 4 (dead)',
            'State 3 (severe)   --->   State 4 (dead)')

### Intercepts ###

# Compute coverages
covs <- matrix(0,22,50)

index <- cmap[1,1:5] # Indices from cmap.
truth <- t(trueBetaMat)[index]

# Bayes
for(i in 1:50){
  
  covs[, i] <- ifelse(
    truth >= summary(models[[i]])$summary[1:22 , 4] & 
      truth <= summary(models[[i]])$summary[1:22 , 8], 1, 0)
  
  b_0s <- rstan::extract(models[[i]], pars=c("b_0"))
  b_1s <- rstan::extract(models[[i]], pars=c("b_1"))
  
  for(j in 1:5){
    covs[7+j , i] <- ifelse(truth[7+j] >= quantile(b_0s[[1]][,j]-8*b_1s[[1]][,j], 0.025) &
                              truth[7+j] <= quantile(b_0s[[1]][,j]-8*b_1s[[1]][,j], 0.975), 1, 0)
  }
}

# msm
covg_msm <- c(0,0,0,0,0)
for(l in 1:5){
  covg_msm[l] <- round( (sum( lower_msm[,index[l]] - 8*trueBetaMat[l,2]  <= truth[l] & 
                                truth[l] <= upper_msm[,index[l]] - 8*trueBetaMat[l,2]))/ 50, 3)
}

# Box plots
pdf(paste0("./New_work/",'cav_intercept.pdf')) 
for(l in 1:5){
  par( mar=c(5, 2.5, 2.5, 2.5))
  b_cov <- round(mean(covs[12+l,1:50]), 2)
  intercepts_msm <- full_msm[[1]][7+l, ]
  
  boxplot(output_bayes[7+l, ], mle_res[7+l,], intercepts_msm, 
          names=c('Bayes', 'MLE','MLE (msm)'),
          ylab=NA, 
          main=bquote(hat(beta)[1]^{(.(l))} ~ (.(labels[l]))), 
          cex =1.2)
  abline(h=truth[l],col='green',lwd=4,lty=5)
  mtext(paste0('.95 coverage = ',toString(b_cov),' Bayes, ',toString(covg_msm[l]),' MLE (msm)'), 
        side=1, line=2.5, cex=1.2)
}
dev.off()


### Beta 1 ###

# Compute coverages
index <- cmap[2,1:5]
truth <- trueBetaMat[,2]

# Bayes
for(i in 1:50){
  b_0s <- rstan::extract(models[[i]], pars=c("b_0"))
  b_1s <- rstan::extract(models[[i]], pars=c("b_1"))
  for(j in 1:5){
    covs[12+j , i] <- ifelse(truth[j] >= quantile(b_1s[[1]][,j], 0.025) &
                               truth[j] <= quantile(b_1s[[1]][,j], 0.975), 1, 0)
  }
}

# msm
covg_msm <- c(0,0,0,0,0)
for(l in 1:5){
  covg_msm[l] <- round( (sum( full_msm[[2]][, index[l]]  <= truth[l] & 
                                truth[l] <= full_msm[[3]][, index[l]]))/ 50, 3)
}

# Box plots
pdf(paste0("./New_work/",'cav_b1.pdf')) 
for(l in 1:5){
  par( mar=c(5, 2.5, 2.5, 2.5))
  b_cov <- round(mean(covs[12+l,1:50]), 2)
  intercepts_msm <- full_msm[[1]][12+l, ]
  
  boxplot(output_bayes[12+l, ], mle_res[12+l,], intercepts_msm, 
          names=c('Bayes', 'MLE','MLE (msm)'),
          ylab=NA, 
          main=bquote(hat(beta)[1]^{(.(l))} ~ (.(labels[l]))), 
          cex =1.2)
  abline(h=truth[l],col='green',lwd=4,lty=5)
  mtext(paste0('.95 coverage = ',toString(b_cov),' Bayes, ',toString(covg_msm[l]),' MLE (msm)'), 
        side=1, line=2.5, cex=1.2)
}
dev.off()



### Beta 2 ###

# Compute coverages
index <- cmap[3,1:5]
truth <- trueBetaMat[,3]

# Bayes
for(i in 1:50){
  b_2s <- rstan::extract(models[[i]], pars=c("b_2"))
  for(j in 1:5){
    covs[17+j , i] <- ifelse(truth[j] >= quantile(b_2s[[1]][,j], 0.025) &
                               truth[j] <= quantile(b_2s[[1]][,j], 0.975), 1, 0)
  }
}

# msm
covg_msm <- c(0,0,0,0,0)
for(l in 1:5){
  covg_msm[l] <- round( (sum( full_msm[[2]][, index[l]]  <= truth[l] & 
                                truth[l] <= full_msm[[3]][, index[l]]))/ 50, 3)
}


# Box plots
pdf(paste0("./New_work/",'cav_b2.pdf')) 
for(l in 1:5){
  par( mar=c(5, 2.5, 2.5, 2.5))
  b_cov <- round(mean(covs[17+l,1:50]), 2)
  intercepts_msm <- full_msm[[1]][17+l, ]
  
  boxplot(output_bayes[17+l, ], mle_res[17+l,], intercepts_msm, 
          names=c('Bayes', 'MLE','MLE (msm)'),
          ylab=NA, 
          main=bquote(hat(beta)[2]^{(.(l))} ~ (.(labels[l]))), 
          cex =1.2)
  abline(h=truth[l],col='green',lwd=4,lty=5)
  mtext(paste0('.95 coverage = ',toString(b_cov),' Bayes, ',toString(covg_msm[l]),' MLE (msm)'), 
        side=1, line=2.5, cex=1.2)
}
dev.off()



### Miscclassification probabilities ###

labels <- c('P( observed state 2 | true state 1 )',
            'P( observed state 1 | true state 2 )',
            'P( observed state 3 | true state 2 )',
            'P( observed state 2 | true state 3 )')

# Compute coverages
index <- cmap[1,6:9]
truth <- c( .01, .24, .06, .11)

# Bayes
for(i in 1:50){
  ps <- rstan::extract(models[[i]], pars=c("p"))
  for(j in 1:4){
    covs[j , i] <- ifelse(truth[j] >= quantile(ps$p[,j], 0.025) &
                            truth[j] <= quantile(ps$p[,j], 0.975), 1, 0)
  }
}

# msm
covg_msm <- c(0,0,0,0)
for(l in 1:4){
  covg_msm[l] <- round( (sum( full_msm[[2]][, index[l]]  <= truth[l] & 
                                truth[l] <= full_msm[[3]][, index[l]]))/ 50, 3)
}


# Box plots
pdf(paste0("./New_work/",'cav_p.pdf')) 
par(mfrow=c(2, 2), mar=c(5, 2.5, 2.5, 2.5))
for(l in 1:4){
  par( mar=c(5, 2.5, 2.5, 2.5))
  b_cov <- round(mean(covs[l,1:50]), 2)
  intercepts_msm <- full_msm[[1]][l, ]
  
  boxplot(output_bayes[l, ], mle_res[l,], intercepts_msm, 
          names=c('Bayes', 'MLE','MLE (msm)'),
          ylab=NA, 
          main=bquote(hat(p)[.(l)] ~ (.(labels[l]))), 
          cex =1.2)
  abline(h=truth[l],col='green',lwd=4,lty=5)
  mtext(paste0('.95 coverage = ',toString(b_cov),' Bayes, ',toString(covg_msm[l]),' MLE (msm)'), 
        side=1, line=2.5, cex=1)
}
dev.off()

