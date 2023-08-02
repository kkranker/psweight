loud = FALSE

library(haven)
library(dplyr)
cattaneo2 <- read_dta("C:/Users/kkranker/Code/Stata/Regression-Adjusted-Means/data/cattaneo2.dta")
if(loud) View(cattaneo2)

library(stats)
lm.out1 <- lm(bweight ~ mbsmoke + mmarried + mage + fbaby + medu, data = cattaneo2)
summary(lm.out1)
lm.out1$coefficients[2]

library(lmw)
lmw.out1 <- lmw( ~ mbsmoke + mmarried + mage + fbaby + medu, data = cattaneo2,
                estimand = "ATE", method = "URI",
                treat = "mbsmoke")
if(loud) print(lmw.out1)
summary(lmw.out1)

lmw.fit1 <- lmw_est(lmw.out1, outcome = "bweight")
if(loud) print(lmw.fit1)
summary(lmw.fit1)


## unpack

X <- cattaneo2 %>%
  mutate(one=1) %>%
  select(one, mbsmoke, mmarried, mage, fbaby, medu) %>%
  as.matrix()
w <- rep(1, nrow(X))
rw <- sqrt(w)
#qr_X <- qr(rw*X)
qr_X <- qr(X)
qr_X[["qr"]]
p <- qr_X$rank
XtX1 <- chol2inv(qr_X$qr[1:p, 1:p, drop = FALSE])
XtX1