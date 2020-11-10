## mle_gompertz_simu_for_working_team.R

## josh g. (Nov 9, 2020)

## Overview: we demonstrate that the truncated gompertz is possible to fit successfully to
## simulated gompertz data.

## As extra-credit, we do with a covariate

## And as extra-extra credit, we estimate SE of our estimates using Hessian.

## Notes:
## * We use mode M parameterization of Gompertz
## * Not sure which optimizer works best but I just use default
## * Problems with optimization often are coding problems with the objective function.


## Contents:
## 1. Gompertz functions
## 2. Simulate data (without covariate) with truncation
## 3. ML estimates
## 4. Simulate data with a covariate (with trunc) and re-estimate
## 5. Hessian example for estimated standard errors


library(txtplot) ## for ascii plots within terminal


## 1. Gompertz functions

library(flexsurv)
source("gompertz_functions.R") ## has functions like rgomp_mode()

## 2. Simulate data (without covariate) with truncation
set.seed(314)
n = 1000
M.true = 88
beta.true = 1/10
u = 90 ## upper bound
l = 80 ## lower bound
y = rgomp_mode(n = n, M = M.true, beta = beta.true) ## untruncated ages at death
yt = y[l < y & y < u] ## truncated ages at death


## 3. MLE estimation with known truncation

minusLogLik_1 =  function(yt, p, u, l)
{
    M = exp(p[1])
    beta = exp(p[2])
    n = length(yt)
    ## denom : F(u) - F(l)
    denom = pgomp_mode(u, M = M, beta = beta) - pgomp_mode(l, M = M, beta = beta)
    ## L = f1/denom * f2/denom ... = f1*f2*f3 / denom^n
    ## LL = sum(log(fi)) - n * log(denom)
    LL = sum(dgomp_mode(yt, M = M, beta = beta, log = TRUE)) - n* log(denom)
    mLL = -LL
    return(mLL)
}


## Maximize the likelihood (by minimizing the minusLogLik)

## starting values
p.start = c(log.M = log(M.true * .8), ## kind of close to true values
            log.beta = log(beta.true * 1.2))

## optimizer
fit1 = optim(par = p.start, fn = minusLogLik_1, yt = yt, u = u, l = l)

## get out estimates
log.est = fit1$par
est = round(exp(log.est),3)
names(est) = c("M.hat", "beta.hat")
true <- c("M.true" = M.true, "beta.true" = beta.true)
print(cbind(est, true))
##             est true
## M.hat    87.313 88.0
## beta.hat  0.088  0.1

## pretty good match (if you increase sample size to 10000, probably
## do even better) or could simulate this smalelr sample many times.

## Note: you can use Hessian form optimizer to get estimated variance (and SE) of estimates



## 4. Simulate data with a covariate (with trunc) and re-estimate

## Note: we use proportional hazards gompertz model

## h(x)[i] = alpha * exp(beta * x) * exp(b*x[i])
## where exp(b * x[i]) is the proportional effect of x on hazards



## Simulate with n obs
## set.seed(23)
n = 10000 ## bigger sample
x = rnorm(n) ## just let the covariate be normal(0,1)
b.true = .6 ## the coeficient (effect) we will try to estimate
alpha = 10^-4.5
rate.vec = alpha * exp(x*b.true) ## proportional hazard parameters for our simulated sample
M.vec <- getMode(alpha = rate.vec, beta = beta.true) ## in terms of mode
## now simulate, taking advantage of ability to use vectors as arguments for parameter values
y = rgomp_mode(n = n, M = M.vec, beta = beta.true)
## truncated death ages
yt = y[l < y & y < u] ## truncated ages at death
## corresponding truncated coviareate
xt = x[l < y & y < u]


minusLogLik_2 =  function(yt, xt, p, u, l)
{
    M = exp(p[1])
    beta = exp(p[2])
    b = p[3] ## no exp()
    alpha = getAlpha(M = M, beta = beta)
    rate.vec = alpha * exp(xt * b)
    M.vec <- getMode(alpha = rate.vec, beta = beta)
    ## denom : F(u) - F(l)
    denom = pgomp_mode(u, M = M.vec, beta = beta) - pgomp_mode(l, M = M.vec, beta = beta)
    ## replace 0s with a very small number for numerical stability (before taking logs)
    eps = 10^-6
    denom[denom == 0] <- eps
    ## L = f1/denom * f2/denom ... = f1*f2*f3 / denom^n
    ## LL = sum(log(fi)) - sum( log(denom))
    LL = sum(dgomp_mode(yt, M = M.vec, beta = beta, log = TRUE)) - sum(log(denom))
    mLL = -LL
    return(mLL)
}

## starting values
p.start = c(log.M = log(M.true * .9), ## kind of close to true values
            log.beta = log(beta.true * 1.2),
            b = b.true * .5)

## optimizer
fit2 = optim(par = p.start, fn = minusLogLik_2, yt = yt, xt = xt, u = u, l = l)

## get out estimates
log.est = fit2$par
est = round(c(exp(log.est[1:2]), log.est[3]),3)
names(est) = c("M.hat", "beta.hat", "b.hat")
true <- c("M.true" = M.true, "beta.true" = beta.true, b.true = b.true)
print(cbind(est, true))
##             est true
## M.hat    78.982 88.0
## beta.hat  0.075  0.1
## b.hat     0.728  0.6
## ok, not terrible for b.hat (note: doesn't work so well with smaller sample)


## 5. Hessian example for estimated standard errors
set.seed(123)
fit2.hess = optim(par = p.start, fn = minusLogLik_2, yt = yt, xt = xt, u = u, l = l,
                  hessian = TRUE)
##https://stats.stackexchange.com/questions/27033/in-r-given-an-output-from-optim-with-a-hessian-matrix-how-to-calculate-paramet
fit<- fit2.hess
## fisher_info<-solve(-fit$hessian)
H = fit$hessian
fisher_info = solve(H)
sigma.hat<-sqrt(diag(fisher_info))
upper<-fit$par+1.96*sigma.hat
lower<-fit$par-1.96*sigma.hat
interval<-data.frame(value=fit$par, upper=upper, lower=lower)
##               value      upper      lower
## log.M     4.3730357  4.4242775  4.3217938
## log.beta -2.3946092 -2.0045801 -2.7846383
## b         0.6719671  0.8268553  0.5170789

exp(interval[1:2,])
##                value      upper       lower
## log.M    79.28394719 83.4524931 75.32362485
## log.beta  0.09120831  0.1347169  0.06175142

## ok, this all looks reasonable but I'm not 100% sure.
## And I doubt that exp(interval) is really exactly correct.
