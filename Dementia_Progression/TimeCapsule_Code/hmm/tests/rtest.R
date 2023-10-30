library(hmm)
#
# Test the continuous response functions
#
f1 <- function(y, eta, gradient=FALSE) {
    cuts <- matrix(c(-Inf, 1.3, 1.5, 1.3, 1.5, Inf), 3, 2)
    hmmncut(y, nstate=4, eta=eta, gradient=gradient, cuts, c(1,2,3,0))
}

# Check out the function
yy <- seq(1, 2, by=.25)
test1 <- f1(yy, -2)
ny <- length(yy)
tmat <- matrix(0, 4, ny)
estd <- exp(-2)

tmat[1,] <- pnorm(1.3-yy, 0, estd)
tmat[2,] <- pnorm(1.5-yy, 0, estd) - pnorm(1.3-yy, 0, estd)
tmat[3,] <- 1 - pnorm(1.5-yy, 0, estd)
all.equal(tmat, test1)

# Now the derivatives wit respect to eta
tfun <- function(y, eta, c1, c2) {
    s <- exp(eta)
    pnorm(c2-y, 0, s) - pnorm(c1-y, 0, s)
}
eps= 1e-8
test2 <- attr(f1(yy, -2, gradient=TRUE), "gradient")

tmat[1,] <- (tfun(yy, eps-2, -Inf, 1.3) - tfun(yy, -2, -Inf, 1.3))/eps
tmat[2,] <- (tfun(yy, eps-2,  1.3, 1.5) - tfun(yy, -2,  1.3, 1.5))/eps
tmat[3,] <- (tfun(yy, eps-2,  1.5, Inf) - tfun(yy, -2,  1.5, Inf))/eps

all.equal(tmat, test2[,,1], tolerance=sqrt(eps))

# Repeat with the logistic
f2 <- function(y, eta, gradient=FALSE) {
    cuts <- matrix(c(-Inf, 1.3, 1.5, 1.3, 1.5, Inf), 3, 2)
    hmmlcut(y, nstate=4, eta=eta, gradient=gradient, cuts, c(1,2,3,0))
}

test1 <- f2(yy, -2)
estd <- exp(-2)*sqrt(3)/pi

tmat[1,] <- plogis(1.3-yy, 0, estd)
tmat[2,] <- plogis(1.5-yy, 0, estd) - plogis(1.3-yy, 0, estd)
tmat[3,] <- 1 - plogis(1.5-yy, 0, estd)
all.equal(tmat, test1)

# Now the derivatives wit respect to eta
tfun <- function(y, eta, c1, c2) {
    s <- exp(eta) * sqrt(3)/pi
    plogis(c2-y, 0, s) - plogis(c1-y, 0, s)
}
eps= 1e-8
test2 <- attr(f2(yy, -2, gradient=TRUE), "gradient")

tmat[1,] <- (tfun(yy, eps-2, -Inf, 1.3) - tfun(yy, -2, -Inf, 1.3))/eps
tmat[2,] <- (tfun(yy, eps-2,  1.3, 1.5) - tfun(yy, -2,  1.3, 1.5))/eps
tmat[3,] <- (tfun(yy, eps-2,  1.5, Inf) - tfun(yy, -2,  1.5, Inf))/eps

all.equal(tmat, test2[,,1], tolerance=sqrt(eps))
