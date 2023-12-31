#
# Fit the same model with msm and hmm, to verify that I have the
#  right likelihood.
library(msm)
library(hmm)
# source("../trial/loadall.R") # for testing without loading
# load('../data/test1.rda')  #simple data set with 4 subjects

# The age intervals are about a year, so make transitions around 5-15% per
# year.  This makes the loglik far from 1.
sname <- levels(test1$state)
qmat <- matrix(0, 6, 6, dimnames=list(from= sname, to=sname))
qmat[1,2] <- qmat[1,3] <- .05
qmat[2,4] <- .06
qmat[3,4] <- .07
qmat[3,5] <- .05
qmat[4,5] <- .15
qmat[1:5,6]  <- c(.05, .05, .05, .07, .15)

# first, a non HMM
test1$istate <- as.numeric(test1$state)
# msm does not use the subset argument
mfit1 <- msm(istate ~ age, data=test1, subject=clinic, 
             qmatrix = qmat, fixedpar=TRUE)

#
# one can fool an HMM into doing non-HMM by starting with p0= 1 for
#  all the possible starting states, and then making the
#  response function  "perfection".
#
init4 <- function(nstate, ...) {
    init <- rep(0, nstate)
    init[1:4] <- 1
    init
}
    
hfit1 <- hmm(hbind(age, state) ~ 1, data=test1, mc.cores=1,
             id = clinic, qmatrix = qmat, rfun=hmm_noerror,
             pfun=init4, mfun=hmmtest, mpar=list(fn="hmmloglik"))

all.equal(-2*hfit1$loglik, mfit1$minus2loglik)

# do the computation by hand
phat <- matrix(1, 4, 4)
idlist <- unique(test1$clinic)
rmat <- qmat
diag(rmat) <- diag(rmat) - rowSums(rmat)
for (i in 1:4) {
    tdata <- subset(test1, clinic==idlist[i])
    delta <- diff(tdata$age)
    oldstate <- tdata$istate[-nrow(tdata)]
    newstate <- tdata$istate[-1]

    for (j in 1:length(delta)) {
        P <- expm(rmat * delta[j])
        phat[i,j] <- P[oldstate[j], newstate[j]]
    }
}
p2 <- apply(phat, 1, prod)
all.equal(sum(log(p2)), hfit1$loglik)


# Add covariates
mfit1b <- msm(istate ~ age, data=test1, subject=clinic, 
              qmatrix = qmat, fixedpar=TRUE,
              covariates= list("1-3"= ~educ,  "1-6"= ~male, "2-6"= ~male),
              covinits= list(educ=.1, male=c(.2, .3)))

q1b <- data.frame(state1=c(1,1,2), state2=c(3,6,6), term=c(1,2,2), coef=1:3,
                  init=c(.1, .2, .3))

# To get the exact same answer we have to use the exact same centering as
#  msm uses
center <- attr(mfit1b$data$mm.cov, "means")
tscale <- list(center=c(0, center), scale=rep(1, 1+length(center)))
hfit1b <- hmm(hbind(age, state) ~ educ + male, data=test1, mc.cores=1,
              id = clinic, qmatrix = qmat, rfun=hmm_noerror, qcoef=q1b,
              pfun=init4, mfun=hmmtest, mpar=list(fn="hmmloglik"),
              scale=tscale)
all.equal(-2*hfit1b$loglik, mfit1b$minus2loglik)  

# do this second computation by hand
# Note that msm uses centered covariates, and hmm has an option for 
#  centering and/or scaled. 
# msm and hsm has somewhat differnent estimates of a "mean".
phat2 <- matrix(1, 4, 4)
for (i in 1:4) {
    tdata <- subset(test1, clinic==idlist[i])
    delta <- diff(tdata$age)
    oldstate <- tdata$istate[-nrow(tdata)]
    newstate <- tdata$istate[-1]

    for (j in 1:length(delta)) {
        rmat <- qmat
        rmat[1,3] <- rmat[1,3] *exp(q1b$init[1]*(tdata$educ[j]-center[1]))
        rmat[1,6] <- rmat[1,6] *exp(q1b$init[2]*(tdata$male[j]-center[2]))
        rmat[2,6] <- rmat[2,6] *exp(q1b$init[3]*(tdata$male[j]-center[2]))
        diag(rmat) <- diag(rmat) - rowSums(rmat)
        P <- expm(rmat * delta[j])
        phat2[i,j] <- P[oldstate[j], newstate[j]]
    }
}
p2 <- apply(phat2, 1, prod)
all.equal(sum(log(p2)), hfit1b$loglik) # succeeds

# Repeat, and let the routines know that 6=death is exact
mfit2 <- msm(istate ~ age, data=test1, subject= clinic, 
             qmatrix = qmat, fixedpar=TRUE, death=6)
otype <- 1 + 1*(test1$istate==6)
hfit2 <-   hmm(hbind(age, state) ~ 1, data=test1, mc.cores=3,
               id = clinic, qmatrix = qmat, rfun= hmm_noerror,
               pfun=init4, mfun=hmmtest, mpar=list(fn="hmmloglik"),
               otype= otype, death=6)
all.equal(-2*hfit2$loglik, mfit2$minus2loglik)

# Repeat again, with misclassification probabilities
# Here is a miss function for 6 states and fixed probs
e1 <- .12  # an A- as A+ or vice versa
e2 <- .2   # an N- as N+ or vice versa
temp <- outer(c("Acorrect"= (1-e1), "Afalse"= e1), 
              c("Ncorrect"=(1-e2), "Nfalse"= e2), '*')

missmat <- rbind(temp[c(1,2,3,4)], temp[c(2, 1, 4,3)],
                 temp[c(3,4,1,2)], temp[c(4, 3, 2, 1)])
missmat <- cbind(missmat, 0, 0)
missmat <- rbind(missmat, c(0,0,0,.1,.9,0), c(0,0,0,0,0,1))

hmiss <- function(y, ...) {
    missmat[,y] 
}

mfit3 <- msm(istate ~ age, data=test1, subject= clinic, 
             qmatrix = qmat, fixedpar=TRUE, death=6,
             ematrix=missmat, initprob=c(1,1,1,1,0,0)/4)

init6 <- function(nstate, ...) {
    c(1,1,1,1,0,0)/4
}
hfit3 <-  hmm(hbind(age, state) ~ 1, data=test1, mc.cores=3,
               id = clinic, qmatrix = qmat, rfun=hmiss,
               pfun=init6, mfun=hmmtest, mpar=list(fn="hmmloglik"),
               otype= otype, death=6)
all.equal(-2*hfit3$loglik, mfit3$minus2loglik)
 

# These models use the cav data set from the msm package
Qm <- rbind(c(0, .148, 0, .0171),
            c(0,  0,  .202, .081),
            c(0,  0,   0,  .126),
            c(0,  0,   0,   0))  #page 38, msm manual

ematrix <- rbind(c(.8, .2, 0, 0),
                 c(.1, .8, .1, 0),
                 c(0, 0.3, .7, 0),
                 c(0, 0,    0, 1))


mfit4a <- msm(state ~ years, subject=PTNUM, data=cav,
            qmatrix=Qm, ematrix=ematrix, death=4,
            obstrue= firstobs, fixedpar=TRUE)

cinit <- function(nstate, ...) {
    init <- rep(0.0, nstate)
    init[1] <- 1   # everyone starts in state 1
    init
}

otype <- with(cav, ifelse(state==4, 2, 1))
first <- !duplicated(cav$PTNUM)
otype[first] <- 0   # the first state is exact, so no response function.

efun <- function(y, nstate, eta, ...) {
    ptemp <- exp(eta[1,])
    emat <- diag(4)
    emat[1,2] <- ptemp[2]    # fill in by column
    emat[2,1] <- ptemp[1]
    emat[2,3] <- ptemp[4]
    emat[3,2] <- ptemp[3]
    emat <- emat/rowSums(emat)
    emat[,y, drop=FALSE]
}

rcoef <- data.frame(lp=1:4, term=0, coef=1:4, 
                    init= log(c(1/8, 2/8, 3/7, 1/8)))

# Verify that I have it set up correctly
# On the linear predictor scale each element is log(e[i,j]/e[i,i])
all.equal(ematrix, efun(1:4, 4, matrix(rcoef$init, 1)))


hfit4a <- hmm(cbind(years, state) ~ 1, data=cav, mc.cores=1,
            id = PTNUM, qmatrix=Qm, rfun=efun, pfun=cinit,
            death=4, otype=otype,  rcoef=rcoef, mfun= hmmtest)
all.equal(-2*hfit4a$loglik, mfit4a$minus2loglik)


# Now with covariates
mfit4b  <- msm(state ~ years, subject=PTNUM, data=cav,
               qmatrix=Qm, ematrix=ematrix, death=4, 
               obstrue= firstobs, covariates=list("1-2"= ~sex, "1-4"= ~sex),
               fixedpar=TRUE, covinits=list(sex=c(1.1, 2.1)))

mscale <- list(center=c(0, attr(mfit4b$data$mm.cov, "means")),
               scale=rep(1.0, 2)) 
qcoef <- data.frame(state1=c(1,1), state2=c(2,4), 
                    term =1, coef=1:2, init=c(1.1, 2.1))
hfit4b <- hmm(cbind(years, state) ~ sex, data=cav, mc.cores=1,
              id = PTNUM, qmatrix=Qm, rfun=efun, pfun=cinit,
              death=4, otype=otype,  rcoef=rcoef, mfun= hmmtest,
              qcoef=qcoef, scale=mscale)
all.equal(-2*hfit4b$loglik, mfit4b$minus2loglik)

