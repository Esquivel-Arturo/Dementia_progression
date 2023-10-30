# A run with covariates.  This uses the cav data from the msm package
library(msm)
library(hmm)

Qm <- rbind(c(0, .148, 0, .0171),
            c(0,  0,  .202, .081),
            c(0,  0,   0,  .126),
            c(0,  0,   0,   0))  #page 38, msm manual

ematrix <- rbind(c(.8, .2, 0, 0),
                 c(.1, .8, .1, 0),
                 c(0, 0.3, .7, 0),
                 c(0, 0,    0, 1))

mfit4a  <- msm(state ~ years, subject=PTNUM, data=cav,
               qmatrix=Qm, ematrix=ematrix, death=4, fixedpar=TRUE,
               obstrue= firstobs, covariates=list("1-2"= ~sex, "1-4"= ~sex),
               covinits=list(sex=c(-.5, .2)))

cinit <- function(nstate, ...) {
    init <- rep(0.0, nstate)
    init[1] <- 1   # everyone starts in state 1
    init
}
 
# the first state is exact, so no response function.
otype <- ifelse(cav$firstobs, 0, ifelse(cav$state==4, 2,1))

rcoef <- data.frame(lp=1:4, term=0, coef=1:4, 
                    init= log(c(2/8, 1/8, 1/8, 3/7)))
errfun <- function(y, nstate, eta, gradient) {
    # true states 1, 2, and 3 have separate linear predictors
    temp1 <- hmulti(y, nstate, eta[,1], gradient,
                    statemap= rbind(1:2, 0,0,0))
    temp2 <- hmulti(y, nstate, eta[,2:3], gradient,
                    statemap= rbind(0, c(2,1,3), 0, 0))
    temp3 <- hmulti(y, nstate, eta[,4], gradient,
                    statemap= rbind(0, 0, c(3,2), 0))
    pmat <- rbind(temp1[1,], temp2[2,], temp3[3,], ifelse(y==4,1,0))
    if (gradient) {
        gmat <- array(c(attr(temp1, 'gradient'), attr(temp2, 'gradient'),
                        attr(temp3, 'gradient')), dim=c(nstate, length(y), 4))
        attr(pmat, "gradient") <- gmat
    }
    pmat
}
        
# Verify that I have it set up correctly
# On the linear predictor scale each element is log(e[i,j]/e[i,i])
eta.temp <- outer(c(1,1,1,1), rcoef$init, '*')
temp <- errfun(1:4, 4, eta.temp, FALSE)
all.equal(ematrix, temp)

# check out the derivatives wrt the 4 linear predictors
eps <- 1e-8
etemp <- errfun(1:4, 4, eta.temp, TRUE)

for (i in 1:4) {
    tinit <- rcoef$init
    tinit[i] <- tinit[i] + eps
    dtemp <- (errfun(1:4, 4, outer(c(1,1,1,1), tinit, '*'), FALSE) -temp)/ eps
    print(all.equal(dtemp, attr(etemp, "gradient")[,,i], tol=sqrt(eps)))
}

q4 <- data.frame(state1= c(1,1,2,2,3,1, 1),
                 state2= c(2,4,3,4,4,2, 4),
                 term  = c(0,0,0,0,0,1, 1),
                 coef  = 1:7,
                 init  = mfit4a$estimates[1:7]) 

# Now do derivatives.
# To reprise the msm answers use the same centering
hscale <- list(center=c(0, attr(mfit4a$data$mm.cov, "means")), scale=c(1,1))

hfit4a <- hmm(cbind(years, state) ~ sex, data=cav, mc.cores=5,
              id = PTNUM, qmatrix=Qm, rfun=errfun, pfun=cinit,
              death=4, otype=otype,  rcoef=rcoef, qcoef= q4, 
              mfun= "hmmtest", mpar=list(fn="hmmboth"),
              scale=hscale) 

dvec4 <- double(length(hfit4a$coef))
for (i in 1:7) {
    qtemp <- q4
    qtemp$init[i] <- qtemp$init[i] + eps
    tfit <- hmm(cbind(years, state) ~ sex, data=cav, mc.cores=10,
                id = PTNUM, qmatrix=Qm, rfun=errfun, pfun=cinit,
                death=4, otype=otype,  rcoef=rcoef, qcoef= qtemp, 
                mfun= "hmmtest", scale=hscale) 
    dvec4[i] <- tfit$loglik
}
for (i in 1:4) {
    rtemp <- rcoef
    rtemp$init[i] <- rcoef$init[i] + eps
    tfit <- hmm(cbind(years, state) ~ sex, data=cav, mc.cores=10,
                id = PTNUM, qmatrix=Qm, rfun=errfun, pfun=cinit,
                death=4, otype=otype,  rcoef=rtemp, qcoef= q4, 
                mfun= "hmmtest", scale=hscale) 
    dvec4[i+7] <- tfit$loglik
}    
d4 <- (dvec4 - hfit4a$loglik)/eps
all.equal(d4, hfit4a$fit$deriv, tolerance=sqrt(eps))
# msm thinks of -2*loglik as the function of interest
all.equal(hfit4a$fit$deriv, mfit4a$deriv[-(8:10)]/-2)

# check the hmmgrad function
#  the derivative ends up labeled as "loglik" (no user would ever do this)
hfit4c <- hmm(cbind(years, state) ~ sex, data=cav, mc.cores=5,
              id = PTNUM, qmatrix=Qm, rfun=errfun, pfun=cinit,
              death=4, otype=otype,  rcoef=rcoef, qcoef= q4, 
              mfun= "hmmtest", mpar=list(fn="hmmgrad"),
              scale=hscale) 
all.equal(hfit4a$fit$deriv, hfit4c$loglik)

mfit4  <- msm(state ~ years, subject=PTNUM, data=cav,
              qmatrix=Qm, ematrix=ematrix, death=4,
              obstrue= firstobs, covariates=list("1-2"= ~sex, "1-4"= ~sex),
              covinits=list(sex=c(-.5, .2)))
hfit4 <- hmm(cbind(years, state) ~ sex, data=cav, mc.cores=6,
            id = PTNUM, qmatrix=Qm, rfun=errfun, pfun=cinit,
            death=4, otype=otype,  rcoef=rcoef, qcoef= q4, 
            mfun= optim,  
            mpar=list(control=list(fnscale= -1), method='BFGS', gr= "hmmgrad"))

all.equal( -2*hfit4$loglik, mfit4$minus2loglik)

# The coefs from msm have some extra bits: slots for coefficients we didn't
#  use (sex effect on the 1-3 transition for instance), number of categories
#  of y, etc.  Due to different optim defaults we don't exactly match values.
indx <- c(1:7, which(names(mfit4$estimates)=="p"))
rbind(hmm=hfit4$coefficients, msm=mfit4$estimates[indx])

# Check the icoef option
hfit4x <- hmm(cbind(years, state) ~ sex, data=cav, mc.cores=6,
              id = PTNUM, qmatrix=Qm, rfun=errfun, pfun=cinit,
              death=4, otype=otype,  rcoef=rcoef, qcoef= q4, 
              mfun="hmmtest", mpar=list(fn="hmmboth"), 
              icoef = hfit4$coef)
all.equal(hfit4x$loglik, hfit4$loglik)

# This last is surprising, given that we start at the same place.
#   msm must have a tighter convergence criteria
rbind(hmm= hfit4$fit$counts, msm=mfit4$paramdata$opt$counts)

# derivatives are closer to zero as well
range(abs(hfit4x$fit$deriv + mfit4$deriv/2))

