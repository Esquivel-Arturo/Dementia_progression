#
# The same fits as jackson.r, but now check gradients
#
# source('../trial/loadall.R')
# load("../data/test1.rda")
library(hmm)

sname <- levels(test1$state)
qmat <- matrix(0, 6, 6, dimnames=list(from= sname, to=sname))
qmat[1,2] <- qmat[1,3] <- .01
qmat[2,4] <- .03
qmat[3,4] <- .04
qmat[3,5] <- .03
qmat[4,5] <- .1
qmat[,6]  <- c(.02, .02, .02, .02, .06, 0)

# first, a non HMM
test1$istate <- as.numeric(test1$state)
test1$id <- match(test1$clinic, unique(test1$clinic))  # easier to use
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
    
# Do a numeric derivative for all parameters, on a subset of the
#  data
eps <- 1e-8
pp <- which(qmat >0)
dvec <- double(length(pp))
lines <- 1:17   #all rows
for (j in 1:length(pp)) {
    qtemp <- qmat
    qtemp[pp[j]] <- exp(log(qmat[pp[j]]) + eps)  # change beta
    tfit <-  hmm(hbind(age, state) ~ 1, data=test1, mc.cores=1,
                 id = id, qmatrix = qtemp, rfun= hmm_noerror, 
                 pfun=init4, mfun=hmmtest,
                 subset=lines)
    dvec[j] <- tfit$loglik
}
hfit <- hmm(hbind(age, state) ~ 1, data=test1, mc.cores=1,
                 id = id, qmatrix = qmat, rfun= hmm_noerror, 
                 pfun=init4, mfun=hmmtest, mpar=list(fn="hmmboth"),
            subset=lines)
d2 <- (dvec- hfit$loglik)/eps
all.equal(d2, hfit$fit$deriv, tolerance=sqrt(eps))
#rbind(d2, hfit$fit$deriv)


# Round 2: let the routines know that 6=death is exact
otype <- 1 + 1*(test1$istate==6)
hfit2 <-   hmm(hbind(age, state) ~ 1, data=test1, mc.cores=1,
               id =  id, qmatrix = qmat, rfun=hmm_noerror,
               pfun=init4, mfun=hmmtest, mpar=list(fn="hmmboth"),
               otype= otype, death=6)
eps <- 1e-8
pp <- which(qmat >0)
dvec <- double(length(pp))
for (j in 1:length(pp)) {
    qtemp <- qmat
    qtemp[pp[j]] <- exp(log(qmat[pp[j]]) + eps)  # change beta
    tfit <-  hmm(hbind(age, state) ~ 1, data=test1, mc.cores=1,
                 id =  id, qmatrix = qtemp, rfun= hmm_noerror,  
                 pfun=init4, mfun=hmmtest,
                 otype=otype, death=6)
    dvec[j] <- tfit$loglik
}
d2 <- (dvec- hfit2$loglik)/eps
all.equal(d2, hfit2$fit$deriv, tolerance=sqrt(eps))


# Round 3: initial probabilites of the vectors are added
pdat <- data.frame(lp=1:3, term=0, coef=1:3, init=c(-1, -2, -3))

hfit3 <-   hmm(hbind(age, state) ~ 1, data=test1, mc.cores=1,
               id =  id, qmatrix = qmat, rfun=hmm_noerror,
               pfun=hmminit, pcoef=pdat, 
               mfun=hmmtest, mpar=list(fn="hmmboth"),
               otype= otype, death=6)

dvec <- double(length(pp) + 3)
for (j in 1:length(pp)) {
    qtemp <- qmat
    qtemp[pp[j]] <- exp(log(qmat[pp[j]]) + eps)  # change beta
    tfit <-  hmm(hbind(age, state) ~ 1, data=test1, mc.cores=1,
                 id =  id, qmatrix = qtemp, rfun= hmm_noerror,  
                 pfun=hmminit, pcoef=pdat, mfun=hmmtest,
                 otype=otype, death=6)
    dvec[j] <- tfit$loglik
}
for (j in 1:3) {
    ptemp <- pdat
    ptemp$init[j] <- ptemp$init[j] + eps
    tfit <-  hmm(hbind(age, state) ~ 1, data=test1, mc.cores=1,
                 id =  id, qmatrix = qmat, rfun= hmm_noerror,  
                 pfun=hmminit, pcoef=ptemp, mfun=hmmtest,
                 otype=otype, death=6)
    dvec[j+length(pp)] <- tfit$loglik
}
d3 <- (dvec- hfit3$loglik)/eps
all.equal(d3, hfit3$fit$deriv, tolerance=sqrt(eps))


# Round 4, add an entry state
otype2 <- otype
otype2[1:3] <- c(0,0,3)
otype2[10:11] <- c(0,3)

evec <- c(1, .3, .4, .5, .1, 0)  # a range of values, on purpose
hfit4 <-   hmm(hbind(age, state) ~ 1, data=test1, mc.cores=1,
               id =  id, qmatrix = qmat, rfun=hmm_noerror,
               pfun=hmminit, pcoef=pdat, entry=evec,
               mfun=hmmtest, mpar=list(fn="hmmboth"),
               otype= otype2, death=6)

dvec <- double(length(pp) + 3)
for (j in 1:length(pp)) {
    qtemp <- qmat
    qtemp[pp[j]] <- exp(log(qmat[pp[j]]) + eps)  # change beta
    tfit <-  hmm(hbind(age, state) ~ 1, data=test1, mc.cores=1,
                 id =  id, qmatrix = qtemp, rfun= hmm_noerror,  
                 pfun=hmminit, pcoef=pdat, mfun=hmmtest,
                 otype=otype2, death=6, entry=evec)
    dvec[j] <- tfit$loglik
}
for (j in 1:3) {
    ptemp <- pdat
    ptemp$init[j] <- ptemp$init[j] + eps
    tfit <-  hmm(hbind(age, state) ~ 1, data=test1, mc.cores=1,
                 id =  id, qmatrix = qmat, rfun= hmm_noerror,  
                 pfun=hmminit, pcoef=ptemp, mfun=hmmtest,
                 otype=otype2, death=6, entry= evec)
    dvec[j+length(pp)] <- tfit$loglik
}
d4 <- (dvec- hfit4$loglik)/eps
all.equal(d4, hfit4$fit$deriv, tolerance=sqrt(eps))
#rbind(d4, hfit4$fit$deriv)

   
# Now add misclassification probabilities into the mix,
#  death and dementia are perfect. The next bit is copied
#  from multi.R
astate <- c(1,2,1,2, NA, NA)[test1$state]
nstate <- c(1,1,2,2, NA, NA)[test1$state]
cstate <- c(NA, NA, NA, NA, 1,2)[test1$state]

err1 <- function(y, nstate, eta, gradient) {
    statemap <- matrix(c(1,2,1,2,0,0,2,1,2,1,0,0), ncol=2)
    hmulti(y, nstate, eta, gradient, statemap)
}

err2 <- function(y, nstate, eta, gradient) {
    statemap <- matrix(c(1,1,2,2,0,0,2,2,1,1,0,0), ncol=2)
    hmulti(y, nstate, eta, gradient, statemap)
}

# death/dementia is perfect
err3 <- function(y, nstate, ...) {
     temp <- matrix(0L, nrow=nstate, ncol=length(y))
     temp[5,] <- ifelse(y==1, 1, 0)
     temp[6,] <- ifelse(y==2, 1, 0)
     temp
}

# hfit4b: the response rates are fixed.  Make sure we don't mess up
#  the other derivatives.
rcoef <- data.frame(response=1:2, lp=1:2, term=0, coef=0,
                    init=log(c(.12/.88, .2/.8)))

hfit4b <-  hmm(hbind(age, astate, nstate, cstate) ~ 1, data=test1, mc.cores=1,
              id = clinic, qmatrix = qmat, 
              rfun = list(err1, err2, err3), rcoef=rcoef,
              pfun=hmminit, pcoef=pdat, mfun=hmmtest, mpar=list(fn="hmmboth"),
              otype= otype2, death=6, entry=evec)

dvec4 <- double(length(pp) + 3)
for (j in 1:length(pp)) {
    qtemp <- qmat
    qtemp[pp[j]] <- exp(log(qmat[pp[j]]) + eps)  # change beta
    tfit <-  hmm(hbind(age, astate, nstate, cstate) ~ 1, data=test1, mc.cores=1,
                 id =  id, qmatrix = qtemp, 
                 rfun = list(err1, err2, err3), rcoef=rcoef,
                 pfun=hmminit, pcoef=pdat, mfun=hmmtest,
                 otype=otype2, death=6, entry=evec)
    dvec4[j] <- tfit$loglik
}
for (j in 1:3) {
    ptemp <- pdat
    ptemp$init[j] <- ptemp$init[j] + eps
    tfit <-  hmm(hbind(age, astate, nstate, cstate) ~ 1, data=test1, mc.cores=1,
                 id =  id, qmatrix = qmat,
                 rfun = list(err1, err2, err3), rcoef=rcoef,
                 pfun=hmminit, pcoef=ptemp, mfun=hmmtest,
                 otype=otype2, death=6, entry= evec)
    dvec4[j+length(pp)] <- tfit$loglik
}
d4b <- (dvec4- hfit4b$loglik)/eps
all.equal(d4b, hfit4b$fit$deriv, tolerance=sqrt(eps))

#
# get derivatives wrt the rates
#
rcoef <- data.frame(response=1:2, lp=1:2, term=0, coef=1:2,
                    init=log(c(.12/.88, .2/.8)))

hfit5 <-  hmm(hbind(age, astate, nstate, cstate) ~ 1, data=test1, mc.cores=1,
              id = clinic, qmatrix = qmat, 
              rfun = list(err1, err2, err3), rcoef=rcoef,
              pfun=hmminit, pcoef=pdat, mfun=hmmtest, mpar=list(fn="hmmboth"),
              otype= otype2, death=6, entry=evec)

all.equal(hfit5$loglik, hfit4b$loglik) # solution is unchanged

dvec5 <- double(length(pp) + 5)
for (j in 1:length(pp)) { #rates
    qtemp <- qmat
    qtemp[pp[j]] <- exp(log(qmat[pp[j]]) + eps)  # change beta
    tfit <-  hmm(hbind(age, astate, nstate, cstate) ~ 1, data=test1,mc.cores=1,
                 id =  id, qmatrix = qtemp, 
                 rfun = list(err1, err2, err3), rcoef=rcoef,
                 pfun=hmminit, pcoef=pdat, mfun=hmmtest,
                 otype=otype2, death=6, entry=evec)
    dvec5[j] <- tfit$loglik
}
for (j in 1:3) { #initial
    ptemp <- pdat
    ptemp$init[j] <- ptemp$init[j] + eps
    tfit <-  hmm(hbind(age, astate, nstate, cstate) ~ 1, data=test1,mc.cores=1,
                 id =  id, qmatrix = qmat, 
                 rfun = list(err1, err2, err3), rcoef=rcoef,
                 pfun=hmminit, pcoef=ptemp, mfun=hmmtest,
                 otype=otype2, death=6, entry= evec)
    dvec5[j+length(pp)+2] <- tfit$loglik
}
for (j in 1:2) {  #response
    rtemp <- rcoef
    rtemp$init[j] <- rtemp$init[j] + eps
    tfit <-  hmm(hbind(age, astate, nstate, cstate) ~ 1, data=test1,mc.cores=1,
                 id =  clinic, qmatrix = qmat, 
                 rfun = list(err1, err2, err3), rcoef=rtemp,
                pfun=hmminit, pcoef=pdat, mfun=hmmtest,
                 otype=otype2, death=6, entry= evec)
    dvec5[j + length(pp)] <- tfit$loglik
}
    
d5 <- (dvec5- hfit5$loglik)/eps
all.equal(d5, hfit5$fit$deriv, tolerance=sqrt(eps))

# Make some of the parameters fixed.  The routine needs to skip over those
#  linear predictors.
qfix <- data.frame(state1=2, state2=4, term=0, coef=0, init=log(.03))
rfix <- rcoef
rfix$coef[1] <- 0
pfix <- pdat
pfix$coef[2] <- 0

hfit6 <-  hmm(hbind(age, astate, nstate, cstate) ~ 1, data=test1, mc.cores=1,
              id = clinic, qmatrix = qmat, qcoef=qfix,
              rfun = list(err1, err2, err3), rcoef=rfix,
              pfun=hmminit, pcoef=pfix, mfun=hmmtest, mpar=list(fn="hmmboth"),
              otype= otype2, death=6, entry=evec)
all.equal(hfit6$fit$deriv, hfit5$fit$deriv[-c(3, 12,15)])
