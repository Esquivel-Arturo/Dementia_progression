### R code from vignette source 'hmmuser.Rnw'

###################################################
### code chunk number 1: model (eval = FALSE)
###################################################
## fit <- hmm(hbind(age, y) ~ iage + ns(iage,3) + male, data=data2,
##            subset=1:50, id=clinic, otype=otype,
##            qmatrix= qmat, qcoef=qcon, rcoef=rcon, pcoef=pcoef)


###################################################
### code chunk number 2: qmat
###################################################
qmat <- matrix(c(0, 1, 1, 0, 0, 1,
                 0, 0, 0, 1, 0, 1,
                 0, 0, 0, 1, 1, 1,
                 0, 0, 0, 0, 1, 1, 
                 0, 0, 0, 0, 1, 1,
                 0, 0, 0, 0, 0, 1), ncol=6, byrow=TRUE)
states <- c("A-N-", "A+N-", "A-N+", "A+N+", "demented", "dead")
dimnames(qmat) <- list(states, states)
qmat


###################################################
### code chunk number 3: qcon
###################################################
states <- c("A-N-", "A+N-", "A-N+", "A+N+", "demented", "dead")
qcon1 <- data.frame(state1 = states[c(1,3,1,2,3,4,5)],
                    state2 = states[c(2,4,3,4,5,5,6)],
                    term = rep(c("ns(iage, 3)", "iage"), c(1,6)),
                    coef = 1:7)
# force equal age coeffiecients for death
qcon2 <- data.frame(state1 = states[rep(1:4, 3)],
                    state2 = rep("dead", 12),
                    term   = rep(c("(Intercept)", "iage","male"), c(4,4,4)),
                    coef = rep(12:14, c(4,4,4)))
rbind(qcon1, qcon2)


