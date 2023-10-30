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

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# First a synthetic data set must be simulated using "Simulate.r" 
load(file="./TimeCapsule_Code/Simulation_output/demData1.rda")

# Use 1 and zero to encode dem diagnosis
demData$DemStatus <- ifelse(demData$DemStatus == 1, 0, 1)

# Remove unobserved cases 
data <- demData[demData$obstype != 0, ]

# mark missing observations and give them a value
data$lpib <- ifelse(is.na(data$lpib), -100, data$lpib)
data$misspib <- ifelse(data$lpib == -100, 0, 1)

data$thickness <- ifelse(is.na(data$thickness), -100, data$thickness)
data$missthick <- ifelse(data$thickness == -100, 0, 1)

data$mmse <- ifelse(is.na(data$mmse), -100, data$mmse)
data$DemStatus <- ifelse(is.na(data$DemStatus), -100, data$DemStatus)

# compute integer years enrolled 
entry_ages <- data |> group_by(ptnum) |>
  summarise(min(age)) |> as.data.frame()

data <- data |> left_join(entry_ages, by=c("ptnum"="ptnum"))

data$iyears <- floor(data$age - data$`min(age)`)

# Basis for state 1 -> state 2 rate cubic splines
Basis <- bs( 50:120, knots=c(55,65,75,90), degree=3, intercept=TRUE)

# It only has been fitted to a small portion of the data
mcsa.data <- list(N = 7,
                 R = 13,
                 T = nrow(data[1:800,]),
                 ID = data[1:800,]$ptnum,
                 pib = data[1:800,]$lpib,
                 thick = data[1:800,]$thickness,
                 mmse = data[1:800,]$mmse,
                 dem = data[1:800,]$DemStatus,
                 sex = data[1:800,]$male,
                 age = data[1:800,]$age,
                 age_i = data[1:800,]$iage,
                 educ = data[1:800,]$educ,
                 apoe = data[1:800,]$apoe4,
                 ntest = data[1:800,]$ntests,
                 iyears = data[1:800,]$iyears,
                 B = Basis,
                 obstype = data[1:800,]$obstype,
                 missp = data[1:800,]$misspib,
                 misst = data[1:800,]$missthick
)

# Fit stan model
fit_mcsa <- stan(file="./New_work/mcsa.stan", data=mcsa.data, refresh=300, iter = 1250, 
                    chains = 3, seed = 29, cores = 3) 
 


