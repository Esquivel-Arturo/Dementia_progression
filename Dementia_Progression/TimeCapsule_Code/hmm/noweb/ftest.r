tfun <- function(formula, data) {
    mf <- model.frame(formula, data=data)
    y <- model.response(mf)
    x <- model.matrix(delete.response(terms(mf)), mf)
    browser()
    y
}

library(survival)
test <- ovarian
test$z <- rep(letters[1:3], length=26)

tfun(Surv(futime, fustat) ~ rx, ovarian)

tfun(cbind(futime, z) ~ rx +age,  data = test)

test$z2 <- factor(test$z)
tfun(cbind(futime, z2) ~ rx, test)

tfun(data.frame(futime, z) ~ rx, test)
