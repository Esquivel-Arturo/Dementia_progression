# Automatically generated from the hmmcode
mlogit <- function(eta, gradient=FALSE) {
    m <- ncol(eta)
    denom <- 1 + rowSums(exp(eta))
    pi    <- cbind(1, exp(eta))/denom
    
    if (gradient) {
        if (nrow(eta) ==1) {  # only one subject
            pi2 <- drop(pi)
            #temp <- diag(pi) - outer(pi, pi) # next line is a touch faster
            temp <- diag(pi2) -  rep(pi2, m+1) * rep(pi2, each= m+1)
            dmat <- temp[,-1,drop=FALSE]
            dim(dmat) <- c(1, dim(dmat))  # make it a "row" per subject
        }
        else {
            dmat <- array(0., dim=c(nrow(eta), m+1, m))
            dmat[,1,] <- -pi[,-1]/denom
            for (j in 1:m) { # for each column of eta
                dmat[,j+1,] <- -pi[,j+1]* pi[,-1]
                dmat[,j+1, j] <- pi[,j+1]*(1-pi[,j+1])
            }
        }
        attr(pi, "gradient") <- dmat
    }
    pi
}
hmminit <- function(nstate, eta, gradient=FALSE) {
    if (is.vector(eta)) eta <- matrix(eta, nrow=1)
    extra <- nstate - (1 + ncol(eta))
    if (extra < 0) 
        stop("too many linear predictors for the number of states")
    temp <- mlogit(eta, gradient)

    if (extra ==0) temp[1,]
    else {
        if (gradient) {
            d2 <- matrix(0, nstate, ncol(eta))
            dmat <- attr(temp, "gradient")[1,,]
            d2[1:nrow(dmat),] <- dmat

            p <- c(as.vector(temp), rep(0., extra))
            attr(p, "gradient") <- d2
            p
        }
        else c(temp, rep(0., extra))
    }
}    
hmmesetup <- function(emat) {
    if (any(is.na(emat)) || any(emat != floor(emat)) || any(emat < -1))
        stop("emat must be an array of integers")
    if (all(emat <1)) stop("no linear predictors are specified")
    temp <- table(row(emat), factor(emat, c(-1, 1:max(emat))),
                  useNA="no")
    if (any(temp[,1] != 1))
        stop("each row of emat must have exactly one reference state")
    if (any(colSums(temp)==0)) stop("there are unused linear predictors")
    index <- tapply(emat, row(emat), function(x) x[x>0])
    
    new  <- matrix(0L, nrow=nrow(emat), ncol=ncol(emat))
    for (i in 1:nrow(emat)) {
        used <- match(emat[i,], index[[i]], nomatch=0)
        new[i, emat[i,]==-1] <- 1  # the reference cell
        new[i, used>0] <- 1 + used[used>0]
    }
    result <- list(index = index, emat= new)
    class(result) <- "hmmesetup"
    result
}
hmmemat <- function(y, nstate, eta, gradient=FALSE, statemap, setup) {
    if (!is.matrix(eta)) eta <- matrix(eta, nrow=length(y))
    if (!inherits(setup, "hmmesetup")) 
        stop("setup must be the result of hmmesetup")
    if (length(statemap) != nstate) stop("wrong length for statemap")
    if (any(statemap <1 | statemap > nrow(setup$emat)))
        stop("statemap does not match setup")
    if (any(y<1 | y > ncol(setup$emat))) stop("y does not match setup")
    
    ny <- length(y)
    rmat <- matrix(0., nrow=nstate, ncol=ny)
    if (gradient) gmat <- array(0., dim=c(nstate, ny, ncol(eta)))
    for (i in 1:length(setup$index)) {
        sindx <- which(statemap ==i)
        if (length(setup$index[[i]]) ==0) {
            # no errors for this state, e.g., death
            indx2 <- which(setup$emat[i,] ==1) # there should be only 1
            rmat[sindx, y %in% indx2] <- 1
            # gradient is zero for these cells
        }
        else {
            mtemp <- mlogit(eta[, setup$index[[i]], drop=FALSE], gradient)
            indx2 <- setup$emat[i,]
            # say that indx2 = 0 2 1 0.  
            #  mtemp will have one row per subject, rmat one col per subject  
            #  rmat[, y==1] is left alone, also rmat[,y==4]; these outcomes
            #      by definition can't happen, so rmat stays at 0
            #  rmat[i, y==2] = mtemp[y==2, indx2[2]] for each i in sindx
            #  rmat[i, y==3] = mtemp[y==3, indx2[3]] for each i in sindx
            #
            #  indx2[y] is the column of mtemp to grab for each y, but mind
            #    the zeros, rows go from 1 to ny
            #  As we march down mtemp from 1 to ny (across rmat) the element
            #    we want is "none" if indx2[y]==0, and 1:ny + ny*(indx2[y]-1)
            #    for the others = standard R indexing.
            #  The last bit of trickery is that rmat will have repeated rows
            #    for elements of sindx, and since the first index varies
            #    fastest we can simply "stutter" the index.
            indx2 <- indx2[y]
            indx3 <- ifelse(indx2==0, 0, 1:ny + ny*(indx2-1)) # matrix index
            indx3 <- rep(indx3, each=length(sindx))  # stutter
            rmat[sindx, indx2>0] <- mtemp[indx3]
            if (gradient) {
                indx4 <- setup$index[[i]]
                for (k in 1:length(indx4)) 
                    gmat[sindx, indx2 >0, indx4[k]] <-
                         (attr(mtemp, "gradient")[,,k])[indx3]
            }
        }
    }
    if (gradient) attr(rmat, "gradient") <- gmat
    rmat
    }
# a multinomial response function
hmulti <- function(y, nstate, eta, gradient=FALSE, statemap) {
    if (!is.matrix(eta)) eta <- matrix(eta, nrow=length(y))
    temp <- mlogit(eta, gradient)
    if (nrow(temp) != length(y)) 
        stop("nrow(eta) != length(y)")
    stopifnot(is.matrix(statemap), nrow(statemap)== nstate,
              ncol(statemap)== (ncol(eta) +1))
    if (gradient) {
        gmat <- array(0., dim=c(nstate, length(y), ncol(eta)))
        gtemp <- attr(temp, "gradient")
    }
    rmat <- matrix(0., nrow= nstate, ncol =length(y))
        
    # This function is called once per subject: both y and nstate are short
    #  A loop is simpler than fancy indexing
    # (It looks like you could do it all at once, but gtemp[vector, vector, ]
    #  will grab too much.)
    for (i in 1:nstate) {
        indx1 <- match(y, statemap[i,], nomatch=0)
        for (obs in which(indx1>0)) {
            rmat[i, obs] <- temp[obs, indx1[obs]]
            if (gradient) gmat[i,obs,] <- gtemp[obs, indx1[obs],]
        }
    }
    if (gradient) attr(rmat, "gradient") <- gmat
    rmat
}
hmmncut <- function(y, nstate, eta, gradient=FALSE, cuts, statemap) {
    xd <- function(x, sd)  # x * dnorm(x)
        ifelse(is.finite(x), x * dnorm(x, 0, sd), 0)
        
    # This only accepts a single parameter
    if (is.matrix(eta) && ncol(eta) != 1) 
        stop("only a single linear predictor is allowed")
    estd <- exp(as.vector(eta))
    if (length(statemap) != nstate) stop("wrong length for statemap")
    if (!is.matrix(cuts) || ncol(cuts) !=2) 
        stop("cuts must be a 2 column matrix")
    if (any(statemap!= floor(statemap) | statemap <0))
        stop ("statemap must be integers >=0")
    if (any(statemap > nrow(cuts)))
        stop("statemap and cuts do not agree")
    yprob <- matrix(0., nstate, length(y))
    if (gradient) ygrad <- array(0., dim=c(nstate, length(y), 1))

    for (i in which(statemap>0)) {
        j <- statemap[i]
        yprob[i,] <- pnorm(cuts[j,2] -y, 0, estd) - pnorm(cuts[j,1] -y, 0, estd)
        if (gradient) 
            ygrad[i,,1] <- xd(cuts[j,1]-y, estd) - xd(cuts[j,2]-y, estd)
    }

    if (gradient) attr(yprob, "gradient") <- ygrad
    yprob
}
        
hmmlcut <- function(y, nstate, eta, gradient=FALSE, cuts, statemap) {
    # The canonical form of the distribution has variance pi^2/3
    #  We want eta=0 to be a variance of 1, so the variance parameter
    #  is rescaled.  
    
    xd <- function(x, sd)  # x * dnorm(x)
        ifelse(is.finite(x), x * dlogis(x, 0, sd), 0)
        
    # This only accepts a single parameter
    if (is.matrix(eta) && ncol(eta) != 1) 
        stop("only a single linear predictor is allowed")
    estd <- exp(as.vector(eta))* sqrt(3)/pi
    if (length(statemap) != nstate) stop("wrong length for statemap")
    if (!is.matrix(cuts) || ncol(cuts) !=2) 
        stop("cuts must be a 2 column matrix")
    if (any(statemap!= floor(statemap) | statemap <0))
        stop ("statemap must be integers >=0")
    if (any(statemap > nrow(cuts)))
        stop("statemap and cuts do not agree")
    yprob <- matrix(0., nstate, length(y))
    if (gradient) ygrad <- array(0., dim=c(nstate, length(y), 1))

    for (i in which(statemap >0)) {
        j <- statemap[i]
        yprob[i,] <- plogis(cuts[j,2] -y, 0, estd) - 
                     plogis(cuts[j,1] -y, 0, estd)
        if (gradient) 
            ygrad[i,,1] <- xd(cuts[j,1]-y, estd) - xd(cuts[j,2]-y, estd)
    }

    if (gradient) attr(yprob, "gradient") <- ygrad
    yprob
}
