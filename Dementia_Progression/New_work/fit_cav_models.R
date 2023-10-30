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
library(matrixStats)
library(ggmcmc)


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
Dir <- "./New_work/NoBiasData/" #/Data for data with enrollment bias

# To store bayesian results
models <- list()

# To store msm results
msm_res <- matrix(0, 22, 50)
upper_msm <- NULL
lower_msm <- NULL

# To store corrected mle results
mle_res <- matrix(0, 22, 50)

for(n in 1:50){
  
  ### Bayesian ###
  
  load(file=paste0(Dir,'cavData', n, '.rda'))
  cavData <- cavData[cavData$obstrue == 0,]
  cav.data <- list(N = 4,
                   R = 5,
                   T = nrow(cavData),
                   ID = cavData$ptnum,
                   Ts = cavData$years,
                   y = cavData$state,
                   sex = cavData$sex,
                   years = cavData$years,
                   years_i = cavData$iyears
  )
  # Fit the model in stan 
  fit_cav_sim <- stan(file="./New_work/cav.stan", data=cav.data, refresh=5000, iter = 1250, 
                      chains = 2, seed = n, cores = 2) 
  # Store model 
  models <- append(models, fit_cav_sim)
  
  
  
  ###### Fit msm ###### 
  # Code adapted from original paper (Williams et al. (2020))
  upper_msm <- NULL
  lower_msm <- NULL
  parEst_msm <- NULL
  qmat <- matrix(c( 0,exp(-2),      0,exp(-2),
                    0,      0,exp(-2),exp(-2),
                    0,      0,      0,exp(-2),
                    0,      0,      0,      0), ncol=4, byrow=TRUE)
  dimnames(qmat) <- list( c('Well', 'Mild','Severe','Death'), c('Well', 'Mild','Severe','Death'))
  
  emat = matrix(c(1, exp(-3), 0, 0,
                  exp(-3), 1, exp(-3), 0,
                        0, exp(-3), 1, 0,
                        0, 0, 0, 1), ncol=4, byrow=TRUE)
  emat = emat / rowSums(emat)
  dimnames(emat) <- list( c('Well','Mild','Severe','Death'), c('Well','Mild','Severe','Death'))
  
  
  Output_msm <- msm(state ~ years, subject=ptnum, data=cavData, qmatrix=qmat, covariates= ~ 1 + iyears + sex,
                    center=FALSE, covinits=list(iyears=c(0,0,0,0,0), sex=c(0,0,0,0,0)), obstrue=obstrue,
                    ematrix=emat, initprobs=c(1, exp(-4.5), exp(-5), 0), est.initprobs=TRUE, deathexact=4,
                    censor=99, censor.states=1:3, method='BFGS', control=list(fnscale=4000, maxit=10000))
  # msm CIs
  for(r in c(5,13,10,14,15)){  
    for(l in 1:3){  
      temp <- c( temp, Output_msm$Qmatrices[[l]][r])  
      temp_upper_msm <- c( temp_upper_msm, 
                           Output_msm$Qmatrices[[l]][r] + 1.96*Output_msm$QmatricesSE[[l]][r])
      temp_lower_msm <- c( temp_lower_msm, 
                           Output_msm$Qmatrices[[l]][r] - 1.96*Output_msm$QmatricesSE[[l]][r])
    }  
  }
  
  for(r in c(5,2,10,7)){
    temp <- c( temp, Output_msm$Ematrices[[1]][r])  
    temp_upper_msm <- c( temp_upper_msm, Output_msm$EmatricesU[[1]][r])
    temp_lower_msm <- c( temp_lower_msm, Output_msm$EmatricesL[[1]][r])
  }
  temp <- c( temp, Output_msm$opt$par[20:21])
  temp_upper_msm <- c( temp_upper_msm, Output_msm$ci[36:37,2])
  temp_lower_msm <- c( temp_lower_msm, Output_msm$ci[36:37,1])
  
  upper_msm <- rbind( upper_msm, temp_upper_msm)
  lower_msm <- rbind( lower_msm, temp_lower_msm)
  
  # Parameter estimates, giving same structure as Stan
  msm_res[1:4, n] <- c(ematrix.msm(Output_msm)[[1]][1,2],
                       ematrix.msm(Output_msm)[[1]][2,1],
                       ematrix.msm(Output_msm)[[1]][2,3],
                       ematrix.msm(Output_msm)[[1]][3,2] )
  
  msm_res[5:7, n] <- Output_msm$estimates.t[35:37]
  
  msm_res[8:12, n] <-c(Output_msm$Qmatrices[[1]][1,2],
                       Output_msm$Qmatrices[[1]][1,4],
                       Output_msm$Qmatrices[[1]][2,3],
                       Output_msm$Qmatrices[[1]][2,4],
                       Output_msm$Qmatrices[[1]][3,4] )
  
  msm_res[13:17, n] <-c(Output_msm$Qmatrices[[2]][1,2],
                        Output_msm$Qmatrices[[2]][1,4],
                        Output_msm$Qmatrices[[2]][2,3],
                        Output_msm$Qmatrices[[2]][2,4],
                        Output_msm$Qmatrices[[2]][3,4] )
  
  msm_res[18:22, n] <-c(Output_msm$Qmatrices[[3]][1,2],
                        Output_msm$Qmatrices[[3]][1,4],
                        Output_msm$Qmatrices[[3]][2,3],
                        Output_msm$Qmatrices[[3]][2,4],
                        Output_msm$Qmatrices[[3]][3,4] )
  
  
  
  ### Fit MLE using the authors' bias correction ###
  # Code adapted from original paper (Williams et al. (2020))
  
  cmap <- matrix(c( 1, 4, 7,10,13,  16,17,18,19,20,21,
                    2, 5, 8,11,14,   0, 0, 0, 0, 0, 0,
                    3, 6, 9,12,15,   0, 0, 0, 0, 0, 0), 3, 11, byrow=TRUE)
  
  rcoef <- data.frame(response = c(1,1,1,1),
                      lp	   	 = c(1,2,3,4),
                      term     = c('(Intercept)','(Intercept)','(Intercept)','(Intercept)'),
                      coef 	 = c(1,2,3,4),
                      init 	 = rep(-3,4), stringsAsFactors=FALSE ) 
  
  pcoef <- data.frame(lp   = c(1,2),		 
                      term = c('(Intercept)','(Intercept)'),
                      coef = c(1,2),
                      init = c(-4.5,-5), stringsAsFactors=FALSE )
  
  qmat <- matrix(c( 0,exp(-2),      0,exp(-2),
                    0,      0,exp(-2),exp(-2),
                    0,      0,      0,exp(-2),
                    0,      0,      0,      0), ncol=4, byrow=TRUE)
  dimnames(qmat) <- list( c('Well', 'Mild','Severe','Death'), c('Well', 'Mild','Severe','Death'))
  
  qcoef <- data.frame(state1 = c(  rep(1,3),  rep(1,3),  rep(2,3), rep(2,3), rep(3,3)),
                      state2 = c(  rep(2,3),  rep(4,3),  rep(3,3), rep(4,3), rep(4,3)),
                      term   = c('(Intercept)','iyears','sex', 
                                 '(Intercept)','iyears','sex',
                                 '(Intercept)','iyears','sex',
                                 '(Intercept)','iyears','sex',
                                 '(Intercept)','iyears','sex'),
                      coef   = c( 1, 2, 3,
                                  4, 5, 6,
                                  7, 8, 9,
                                  10,11,12,
                                  13,14,15),
                      init   = c(-2, 0, 0,
                                 -2, 0, 0,
                                 -2, 0, 0,
                                 -2, 0, 0,
                                 -2, 0, 0), stringsAsFactors=FALSE) 
  
  steps  = 15000
  burnin = 10000
  
  Output_mle <- hmm(hbind(years, state) ~ 1 + iyears + sex, data=cavData, mc.cores=2, entry=c(1,1,1,0), 
                    death=4, id=ptnum, rfun=list(misclassResp), rcoef=rcoef, pfun=pfun, pcoef=pcoef, qmatrix=qmat, 
                    qcoef=qcoef, otype=obstype_hmm, scale=FALSE, mfun=hmmscore, 
                    mpar=list(gr="hmmboth", iter=10000))
  
  temp_upper_mle <- NULL
  temp_lower_mle <- NULL
  
  parEst_mle <- Output_mle$fit$coef	
  
  meanYears <- round( mean(cavData$years), 0)
  center <- meanYears
  N <- length(unique(cavData$ptnum))
  meanSampleSize <- c(meanSampleSize, N)
  index <- cmap[1,1:5] # Indices from cmap.
  truth <- t(trueBetaMat)[index]
  
  for(l in 1:5){
    
    mle_res[7+l, n] <- parEst_mle[index[l]] #- center*parEst_mle[index[l] + 1] # when data are centered
    mle_res[12+l, n] <- parEst_mle[index[l] + 1]
    mle_res[17+l, n] <- parEst_mle[index[l] + 2]
    
  }
  
  index <- cmap[1,6:9] # Indices from cmap.
  # logistic transforms of the parameter etimates.
  mle_res[1:4, n] <- c( exp(parEst_mle[index[1]]) / sum( c( 1, exp(parEst_mle[index[1]])) ),
                        exp(parEst_mle[index[2:3]]) / sum( c( 1, exp(parEst_mle[index[2:3]])) ),
                        exp(parEst_mle[index[4]]) / sum( c( 1, exp(parEst_mle[index[4]])) ) )
  
  index <- cmap[1,10:11] # Indices from cmap.
  
  # logistic transforms of the parameter etimates.
  mle_res[6:7, n] <- exp(parEst_mle[index]) / sum( c( 1, exp(parEst_mle[index]), 0) )
  mle_res[5,n] <- 1-sum(mle_res[6:7, k])
}
  
# save model outputs
save(models, file=paste0("./New_work/NoBiasOutput/models.rda")) # Use /Output/ when enrollment bias
full_msm <- list(msm_res, lower_msm, upper_msm)
save(full_msm, file="./New_work/NoBiasOutput/full_msm.rda")
# Only used for enrollment bias, if no bias, is equivalent to msm
save(mle_res, file=paste0("./New_work/Output/Output_mle.rda"))



# Produce trace plots 
chs12 <- ggs(models[[1]], inc_warmup = TRUE) 
chs34 <- ggs(models[[6]], inc_warmup = TRUE) 
chs34$Chain <- ifelse(chs34$Chain == 1, 3, 4)

chains <- rbind(chs12, chs34)

# Misclassification probs
chains %>% 
  filter(Parameter %in% c("p[1]", "p[2]", "p[3]", "p[4]")) %>%
  mutate(chain = factor(Chain)) %>% 
  ggplot(aes(x = Iteration, y = value)) +
  # this marks off the warmups
  annotate(geom = "rect",
           xmin = 0, xmax = 1250, ymin = -Inf, ymax = Inf,
           alpha = 1/6, linewidth = 0) +
  geom_line(aes(color = chain),
            linewidth = .15)+ 
  theme_bw(base_size = 13) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(plot.title = element_text(face = "bold"), 
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.97, .45))  +
  facet_wrap(~ Parameter, scales = "free_y") 

# Intercepts
chains %>% 
  filter(Parameter %in% c("b_0[1]", "b_0[2]", "b_0[3]", "b_0[4]", "b_0[5]")) %>%
  mutate(chain = factor(Chain)) %>% 
  ggplot(aes(x = Iteration, y = value)) +
  # this marks off the warmups
  annotate(geom = "rect",
           xmin = 0, xmax = 1250, ymin = -Inf, ymax = Inf,
           alpha = 1/6, linewidth = 0) +
  geom_line(aes(color = chain),
            linewidth = .15)+ 
  theme_bw(base_size = 13) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(plot.title = element_text(face = "bold"), 
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.75, .35))  +
  facet_wrap(~ Parameter, scales = "free_y") 

# Beta 1
chains %>% 
  filter(Parameter %in% c("b_1[1]", "b_1[2]", "b_1[3]", "b_1[4]", "b_1[5]")) %>%
  mutate(chain = factor(Chain)) %>% 
  ggplot(aes(x = Iteration, y = value)) +
  # this marks off the warmups
  annotate(geom = "rect",
           xmin = 0, xmax = 1250, ymin = -Inf, ymax = Inf,
           alpha = 1/6, linewidth = 0) +
  geom_line(aes(color = chain),
            linewidth = .15) + 
  theme_bw(base_size = 13) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(plot.title = element_text(face = "bold"), 
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.75, .35))  +
  facet_wrap(~ Parameter, scales = "free_y") 

# Beta 2
chains %>% 
  filter(Parameter %in% c("b_2[1]", "b_2[2]", "b_2[3]", "b_2[4]", "b_2[5]")) %>%
  mutate(chain = factor(Chain)) %>% 
  ggplot(aes(x = Iteration, y = value)) +
  annotate(geom = "rect",
           xmin = 0, xmax = 1250, ymin = -Inf, ymax = Inf,
           alpha = 1/6, linewidth = 0) +
  geom_line(aes(color = chain),
            linewidth = .15) + 
  theme_bw(base_size = 13) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(plot.title = element_text(face = "bold"), 
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.75, .35))  +
  facet_wrap(~ Parameter, scales = "free_y") 

