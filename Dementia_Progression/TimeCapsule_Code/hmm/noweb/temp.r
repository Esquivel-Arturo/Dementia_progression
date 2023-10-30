# Test out my Gaussian derivation
#

d1 <- function(z, s= .1) 
    pnorm(1.5-z, 0, s) - pnorm(1.3-z, 0, s)

d2 <- function(z, s=.1) 
      ((1.3-z)*dnorm(1.3-z, 0,s) - (1.5-z)*dnorm(1.5-z, 0, s))/s


d3 <- function(x, s=.1, eps=1e-6) 
    (d1(x, s+eps) - d1(x, s))/eps

xx <- seq(1, 2, length=51)
plot(xx, d1(xx))
abline(v=c(1.3, 1.5), col=2)

plot(xx, d3(xx))
plot(xx, d2(xx))

#
# do it harder
#
ifun <- function(z, s=.1) {
    tfun <- function(x) dnorm(x, 0, s)* x^2/s^3
    n <- length(z)
    t1 <- double(n)
    for (i in 1:n) t1[i] <- integrate(tfun, 1.3-z[i], 1.5-z[i])$value
    t2 <- d1(z,s)/s
    t1 - t2
}

d1b <- function(z, s=.1) {
    tfun <- function(x) dnorm(x, 0, s)
    integrate(tfun, 1.3-z, 1.5-z)$value
}

matplot(xx, cbind(d3(xx), ifun(xx), d2(xx)), type='l')

itest <- function(z, s=.1) {
    tfun <- function(x) dnorm(x, 0, s)* x^2/s^3
    integrate(tfun, 1.3-z, 1.5-z)$value
}

erf <- function(x) 2*pnorm(x*sqrt(2)) -1
test2 <- function(x, s=.1) erf(x/(s * sqrt(2))/(2*s) - x*dnorm(x,0,s)/s^2)

etest <- function(z, a=50) {
    tfun <- function(x) x^2* exp(-a * x*x)
    t1 <- integrate(tfun, -4, 3)
    tf2 <- function(x) .25*sqrt(pi/a^3)*erf(x*sqrt(a)) - x*exp(-a*x*x)/(2*a)
    c(t1$value, tf2(3) - tf2(-4))
}
e2 <- function(z, s=.1) {
    a <- 1/(2*s*s)
    const <- 1/(s * sqrt(2*pi))
    tfun <- function(x) const * x^2* exp(-a * x*x)/ s^3
    t1 <- integrate(tfun, 1.3-z, 1.5-z)$value
    tf2 <- function(x) c(.25*sqrt(pi/a^3)*erf(x*sqrt(a)) -x*exp(-a*x*x)/(2*a))

    tf3 <- function(x) c(erf(x/(s*sqrt(2)))/(2*s) -x*dnorm(x,0,s)/ s)
    tf4 <- function(x) pnorm(x,0,s)/s - x*dnorm(x,0,s)/s

    c(t1, const *(tf2(1.5-z) - tf2(1.3-z))/ s^3, tf3(1.5-z) - tf3(1.3-z),
      tf4(1.5-z) - tf4(1.3-z))
}

e3 <-function(z, s=.1) 
    c(pnorm(1.5-z ,0,s)/s, -(1.5-z)*dnorm(1.5-z,0,s)/s,
      pnorm(1.3-z ,0,s)/s, -(1.3-z)*dnorm(1.3-z,0,s)/s, d1(z,s)/s)

d2x <- function(z, s=.1) 
    c(d1(z,s)/s^2, d1(z,s)/s, 
      ((1.5-z)*dnorm(1.5-z, 0,s) - (1.3-z)*dnorm(1.3-z, 0, s))/s)


