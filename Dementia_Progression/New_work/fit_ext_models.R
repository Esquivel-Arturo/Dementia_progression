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

# First the data must be simulated using the "Simulate_cav.r" (same folder as this script)

# Remove un-observed cases
cavData <- cavData[cavData$obstrue == 0,]

# Plots showing the data and the real state-dependent distributions
# Only state with response observed (1-3)
obs <- cavData[cavData$obstrue==0 & cavData$obstype_hmm !=2, ]
# Sorted response values observed
x <- sort(obs$obs)

p1 <- ggplot(data.frame(x = c(0, 1))) +
  xlim(c(0, 26)) + ylim(c(0,0.25))+
  stat_function(
    fun = dgamma, args=list(6, 0.95), aes(colour="No Cav")) +
  stat_function(
    fun = dgamma, args=list(48, 4), aes(colour="Mild")) +
  stat_function(
    fun = dsnorm, args=list(18, 3, 0.7), aes(colour="Severe")) +
  scale_color_manual(values = c("No Cav" = "#003f5c", "Mild" = "#bc5090", "Severe" = "#ffa600"), 
                     name = 'State',
                     labels = c("No Cav", "Mild", "Severe")
  ) +
  labs(x = "F-V", y = "Density", title = "A: Sampling Distributions")  +
  theme_bw(base_size = 13)  +
  theme(plot.title = element_text(face = "bold"))

# Real state sequence (generated at the same time as the data)
s1_p <- sum(states[obs$obstrue==0 & obs$obstype_hmm !=2] == 1)/nrow(obs)
s2_p <- sum(states[obs$obstrue==0 & obs$obstype_hmm !=2] == 2)/nrow(obs)
s3_p <- sum(states[obs$obstrue==0 & obs$obstype_hmm !=2] == 3)/nrow(obs)

fx1 <- dgamma(x, 6, 0.95)
fx2 <- dgamma(x, 48, 4)
fx3 <- dsnorm(x, 18, 3, 0.7)
sam_dist <- data.frame(x = x, 
                       fx1 = fx1,
                       fx2 = fx2,
                       fx3 = fx3)


p2 <- ggplot(sam_dist, aes(x = x)) +
  geom_histogram(aes(y = after_stat(density)), fill = "lightgrey") +
  geom_line(aes(y= s1_p*fx1, colour = "No Cav")) + 
  geom_line(aes(y= s2_p*fx2, colour = "Mild")) +
  geom_line(aes(y= s3_p*fx3, colour = "Severe")) +
  geom_line(aes(y= fx1*s1_p + fx2*s2_p + fx3*s3_p, colour = "Mixture"))  +
  scale_color_manual(values = c("No Cav" = "#003f5c", "Mild" = "#bc5090", "Severe" = "#ffa600", "Mixture" = "black"), 
                     name = 'State',
                     labels = c("No Cav" = "No Cav", "Mixture", "Severe", "Mixture")
  ) +
  labs(x = "F-V", y = "Density", title = "B: Empirical vs Sampling Distribution")  +
  theme_bw(base_size = 13)  +
  theme(plot.title = element_text(face = "bold"))

ggarrange(p1, p2, ncol=1, nrow=2, common.legend = TRUE, legend="right") 



### Fit the models ###

### Gaussian Formulation ###
cav.data <- list(N = 4,
                 R = 5,
                 T = nrow(cavData),
                 ID = cavData$ptnum,
                 Ts = cavData$years,
                 y = cavData$state,
                 sex = cavData$sex,
                 years = cavData$years,
                 years_i = cavData$iyears,
                 obs = cavData$obs
)

fit_cav_n <- stan(file="./New_work/cav_normal.stan", data=cav.data, refresh=300, iter = 1500, 
                  chains = 2, seed = 16, cores = 2)




### P-splines ###

# calculation of the B-spline design matrix 
ord<-4
degree<-ord-1
dr <- range(obs$obs)
# Get regularly spaced knots
knots <- c(min(dr) - 0.001, 
           quantile(x = dr, probs = seq(0, 15, length.out = 15)[2:14] / 15), 
           max(dr) + 0.001)
# Basis design matrix
B0<-spline.des(
  knots,
  seq(knots[1], knots[length(knots)], length=100000),
  degree+1,
  outer.ok=T
  )$design 
nb<-ncol(B0)
K<-(nb-1)/2
ws<-rep(NA,nb)
# Regularise the bases
for (k in 1:nb){
  # this computes the integrals of the B-spline basis functions (which are then standardized below)
  ws[k] <- (diff(c(knots[1],knots[length(knots)]))/100000*sum(B0[,k]))^(-1) 
}
# Get the standarised splines for each observed response 
B3<-t(t(spline.des(knots, obs$obs, degree+1, outer.ok=T)$design)*ws) 
# Matrix with same nrow as data entries
B <- matrix(0 , nrow(cavData), ncol(B3))
# Fill entries with observed responses with the corresponding basis value
B[cavData$obstrue==0 & cavData$obstype_hmm !=2, ] <- B3

# fit the model
cav.data <- list(N = 4,
                 R = 5,
                 T = nrow(cavData),
                 K = 13,
                 ID = cavData$ptnum,
                 Ts = cavData$years,
                 y = cavData$state,
                 sex = cavData$sex,
                 years = cavData$years,
                 years_i = cavData$iyears,
                 B1 = B
)

fit_cav_sp <- stan(file="./New_work/cav_splines.stan", data=cav.data, refresh=300, iter = 1500, 
                   chains = 2, seed = 16, cores = 2) 



### Produce comparison plots ###

# Parameters for the normal formulation 
n_fx1 <- dnorm(x, summary(fit_cav_n)$summary[23, 1], summary(fit_cav_n)$summary[26, 1])
n_fx2 <- dnorm(x, summary(fit_cav_n)$summary[24, 1], summary(fit_cav_n)$summary[26, 1])
n_fx3 <- dnorm(x, summary(fit_cav_n)$summary[25, 1], summary(fit_cav_n)$summary[26, 1])
# Credible Intervals
ln_fx1 <- dnorm(x, summary(fit_cav_n)$summary[23, 4], summary(fit_cav_n)$summary[26, 4])
ln_fx2 <- dnorm(x, summary(fit_cav_n)$summary[24, 4], summary(fit_cav_n)$summary[26, 4])
ln_fx3 <- dnorm(x, summary(fit_cav_n)$summary[25, 4], summary(fit_cav_n)$summary[26, 4])

un_fx1 <- dnorm(x, summary(fit_cav_n)$summary[23, 8], summary(fit_cav_n)$summary[26, 8])
un_fx2 <- dnorm(x, summary(fit_cav_n)$summary[24, 8], summary(fit_cav_n)$summary[26, 8])
un_fx3 <- dnorm(x, summary(fit_cav_n)$summary[25, 8], summary(fit_cav_n)$summary[26, 8])
# Gaussian estimated distributions
n_dist <- data.frame(x = x, 
                     fx1 = fx1,
                     fx2 = fx2,
                     fx3 = fx3,
                     es_f1 = n_fx1,
                     es_f2 = n_fx2,
                     es_f3 = n_fx3,
                     l_f1 = ln_fx1,
                     l_f2 = ln_fx2,
                     l_f3 = ln_fx3,
                     u_f1 = un_fx1,
                     u_f2 = un_fx2,
                     u_f3 = un_fx3)

p1 <- ggplot(n_dist, aes(x = x)) +
  geom_histogram(aes(y = after_stat(density)), fill = "lightgrey") +
  geom_line(aes(y= s1_p*fx1, colour = "No Cav", linetype = "True")) + 
  geom_line(aes(y= s2_p*fx2, colour = "Mild", linetype = "True")) +
  geom_line(aes(y= s3_p*fx3, colour = "Severe", linetype = "True")) +
  geom_line(aes(y= fx1*s1_p + fx2*s2_p + fx3*s3_p, colour = "Mixture", linetype = "True")) +
  geom_line(aes(y= s1_p*es_f1, colour = "No Cav", linetype = "Estimated")) + 
  geom_line(aes(y= s2_p*es_f2, colour = "Mild", linetype = "Estimated")) +
  geom_line(aes(y= s3_p*es_f3, colour = "Severe", linetype = "Estimated")) +
  geom_ribbon(aes(ymin=s1_p*l_f1, ymax=s1_p*u_f1), fill="#003f5c", alpha=0.45) +
  geom_ribbon(aes(ymin=s2_p*l_f2, ymax=s2_p*u_f2), fill="#bc5090", alpha=0.45) +
  geom_ribbon(aes(ymin=s3_p*l_f3, ymax=s3_p*u_f3), fill="#ffa600", alpha=0.45) +
  geom_line(aes(y= es_f1*s1_p + es_f2*s2_p + es_f3*s3_p, colour = "Mixture", linetype = "Estimated")) +
  geom_ribbon(aes(ymin=l_f1*s1_p + l_f2*s2_p + l_f3*s3_p,
                  ymax=u_f1*s1_p + u_f2*s2_p + u_f3*s3_p), fill="black", alpha=0.45) +
  scale_linetype_manual(name = "Distribution", values=c("dashed", "solid")) +
  scale_color_manual(values = c("No Cav" = "#003f5c", "Mild" = "#bc5090", "Severe" = "#ffa600", "Mixture" = "black"), 
                     name = 'State',
                     labels = c("No Cav" = "No Cav", "Mixture", "Severe", "Mixture")
  ) +
  labs(x = "F-V", y = "Density", title = "A: Gaussian")  +
  theme_bw(base_size = 13)  +
  theme(plot.title = element_text(face = "bold"))

# Estimated spline coefficients
a1 <- summary(fit_cav_sp)$summary[23:33, 1]
a2 <- summary(fit_cav_sp)$summary[34:44, 1]
a3 <- summary(fit_cav_sp)$summary[45:55, 1]
# Credible intervals
la1 <- summary(fit_cav_sp)$summary[23:33, 4]
la2 <- summary(fit_cav_sp)$summary[34:44, 4]
la3 <- summary(fit_cav_sp)$summary[45:55, 4]

ua1 <- summary(fit_cav_sp)$summary[23:33, 8]
ua2 <- summary(fit_cav_sp)$summary[34:44, 8]
ua3 <- summary(fit_cav_sp)$summary[45:55, 8]

# Weighted to sum to one 
w1 <- exp(a1)/sum(exp(a1))
sp_fx1 <- rowSums(sweep(as.matrix(Bp), MARGIN=2, w1, `*`))
lw1 <- exp(la1)/sum(exp(la1))
lsp_fx1 <- rowSums(sweep(as.matrix(Bp), MARGIN=2, lw1, `*`))
uw1 <- exp(ua1)/sum(exp(ua1))
usp_fx1 <- rowSums(sweep(as.matrix(Bp), MARGIN=2, uw1, `*`))

w2 <- exp(a2)/sum(exp(a2))
sp_fx2 <- rowSums(sweep(as.matrix(Bp), MARGIN=2, w2, `*`))
lw2 <- exp(la2)/sum(exp(la2))
lsp_fx2 <- rowSums(sweep(as.matrix(Bp), MARGIN=2, lw2, `*`))
uw2 <- exp(ua2)/sum(exp(ua2))
usp_fx2 <- rowSums(sweep(as.matrix(Bp), MARGIN=2, uw2, `*`))

w3 <- exp(a3)/sum(exp(a3))
sp_fx3 <- rowSums(sweep(as.matrix(Bp), MARGIN=2, w3, `*`))
lw3 <- exp(la3)/sum(exp(la3))
lsp_fx3 <- rowSums(sweep(as.matrix(Bp), MARGIN=2, lw3, `*`))
uw3 <- exp(ua3)/sum(exp(ua3))
usp_fx3 <- rowSums(sweep(as.matrix(Bp), MARGIN=2, uw3, `*`))
# Splines estimated distributions
sp_dist <- data.frame(x = x, 
                      fx1 = fx1,
                      fx2 = fx2,
                      fx3 = fx3,
                      es_f1 = sp_fx1,
                      es_f2 = sp_fx2,
                      es_f3 = sp_fx3,
                      l_f1 = lsp_fx1,
                      l_f2 = lsp_fx2,
                      l_f3 = lsp_fx3,
                      u_f1 = usp_fx1,
                      u_f2 = usp_fx2,
                      u_f3 = usp_fx3)

p2 <- ggplot(sp_dist, aes(x = x)) +
  geom_histogram(aes(y = after_stat(density)), fill = "lightgrey") +
  geom_line(aes(y= s1_p*fx1, colour = "No Cav", linetype = "True")) + 
  geom_line(aes(y= s2_p*fx2, colour = "Mild", linetype = "True")) +
  geom_line(aes(y= s3_p*fx3, colour = "Severe", linetype = "True")) +
  geom_line(aes(y= fx1*s1_p + fx2*s2_p + fx3*s3_p, colour = "Mixture", linetype = "True")) +
  geom_line(aes(y= s1_p*es_f1, colour = "No Cav", linetype = "Estimated")) + 
  geom_line(aes(y= s2_p*es_f2, colour = "Mild", linetype = "Estimated")) +
  geom_line(aes(y= s3_p*es_f3, colour = "Severe", linetype = "Estimated")) +
  geom_ribbon(aes(ymin=s1_p*l_f1, ymax=s1_p*u_f1), fill="#003f5c", alpha=0.45) +
  geom_ribbon(aes(ymin=s2_p*l_f2, ymax=s2_p*u_f2), fill="#bc5090", alpha=0.45) +
  geom_ribbon(aes(ymin=s3_p*l_f3, ymax=s3_p*u_f3), fill="#ffa600", alpha=0.45) +
  geom_line(aes(y= es_f1*s1_p + es_f2*s2_p + es_f3*s3_p, colour = "Mixture", linetype = "Estimated")) +
  geom_ribbon(aes(ymin=l_f1*s1_p + l_f2*s2_p + l_f3*s3_p,
                  ymax=u_f1*s1_p + u_f2*s2_p + u_f3*s3_p), fill="black", alpha=0.45) +
  scale_linetype_manual(name = "Distribution", values=c("dashed", "solid")) +
  scale_color_manual(values = c("No Cav" = "#003f5c", "Mild" = "#bc5090", "Severe" = "#ffa600", "Mixture" = "black"), 
                     name = 'State',
                     labels = c("No Cav" = "No Cav", "Mixture", "Severe", "Mixture")
  ) +
  labs(x = "F-V", y = "Density", title = "B: P-splines")  +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggarrange(p1, p2, ncol=1, nrow=2, common.legend = TRUE, legend="right") 




