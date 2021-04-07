## Fixed version of Serge's original 'broken' file

## Description of fixes:
## (1) fixed bug in extracting "B" the effect of covariate
## (2) made sure covar_effect.vec was the right dimensions
## (3) worked with whole file, not just sample (still not very big)
## (4) jacked up maxit on optim() so it is 5,000 and converges.
## (5) used "edu" in years, rather than log income

setwd("~/GDrive/demographer/github/censoc_mort")
library(data.table)
library(tidyverse)
source("./gompertz_functions.R")


## 1. Write the mll function



mll.gomp.multi.cohort.cov <- function(p,  ## as many alphas as cohorts, and 1 beta
                                      A)  ## data matrix with y, u, l, and covariates, including cohort
{
    ## (1) get parameters alpha, beta, and B
    ## we'll convert later to one element per individual
    cohorts = names(table(A$cohort))  ## small vector of cohorts, e.g, 1915-1919.
    k <- length(cohorts)
    M = exp(p[1:k]); names(M) = cohorts
    beta = exp(p[k+1]); names(beta) = "beta" ## slope of gompertz
    ## convert M to alpha
    alpha = getAlpha(M, beta); names(alpha) = cohorts

    ## new version just has one covariate
    ## e.g., educ in years

    b = p[length(p)] ## covariate coef is last one

    ## (2) build rate.vec, which has the combined effect of cohort and covariates, 1 element per individual
    alpha.vec <- alpha[paste(A$cohort)] ## this now has one alpha for each individual
    ## note: effects are multiplicative of form haz_i = base * exp(covar_i * b)
    covar_effect.vec <- exp(A$covar * b) ## this gives one covar_effect for each individual
    if (length(alpha.vec) != length(covar_effect.vec))
        print("warning: length(alpha.vec) != length(covar_effect.vec)")
    rate.vec = alpha.vec * covar_effect.vec

    y = A$y
    u.vec = A$u
    l.vec = A$l

    ## (3) obtain likelihood Numerator has difference in cumulative
    ## probabilities, 1 year apart.  This works with data that is
    ## exact to the year of death like we have in CenSoc.
    num = pgompertz(y + 1, shape = beta, rate = rate.vec)-
        pgompertz(y, shape = beta, rate = rate.vec)
    denom = pgompertz(u.vec, shape = beta, rate = rate.vec)-
        pgompertz(l.vec, shape = beta, rate = rate.vec)
    R = num/denom
    minusloglik = -sum(log(R))
    return(minusloglik)
}


## 3. Build data matrix A using CenSoc data

dt = fread("./censoc_numident_demo_v1.csv")  ## demo Numident

cohorts = 1915:1919
## let's try with continuous variable for education
dt[educd == "1 year of college", edu := 13]
dt[educd == "2 years of college", edu := 14]
dt[educd == "3 years of college", edu := 15]
dt[educd == "4 years of college", edu := 16]
dt[educd == "5+ years of college", edu := 18]
dt[grepl("Grade ", educd), edu := gsub("Grade ", "", educd)]
dt[educd == "No schooling completed", edu := NA]
## to check, look at: dt[, table(educd, edu)]

## other covariates
dt$gender <- ifelse(dt$sex == "Male", 0, 1)
dt$rural <- ifelse(dt$urban == "Urban", 0, 1)
#dt$nonwage <- ifelse(is.na(dt$incnonwg), 1, 0)
dt$nonwage <- ifelse(dt$incnonwg == "N/A", 0, 1)

## create sample
my.dt <- dt[dyear %in% 1988:2005 & byear %in% cohorts &
            !is.na(edu)]





## now create the needed format for optimization function

A = my.dt[, .(y = death_age,
              u = 2005 - byear,
              l = 1988 - byear,
              cohort = byear,
              covar = gender)]  ## change out the covar as desired
At = A[l < y & y < u]
## check densities
library(txtplot)
txtdensity(At[cohort == 1917]$y)
txtdensity(At[cohort == 1919]$y)
txtdensity(At[cohort == 1915]$y)
## ok, it looks like we don't have end-of-interval issues


## starting values, just let all the modes be the same ("80"
M.start = rep(80, length(cohorts))
names(M.start) = paste0("coh", cohorts)


## and let slope of gompertz start at .1
## and the effect of 1 year of educ to lower mortality by 10%
p.start = c("log.of" = log(M.start),
            "log.beta" = log(.1),
            "b" = -.1) ## note we use "b" in original scale, not logged

## run optimizer
## big fix is making sure maxit is a big number!
fit = optim(par = p.start, fn = mll.gomp.multi.cohort.cov,
            A = At,
            hessian = TRUE,
            control = list(maxit = 5000))
## let's get standard errors from Hessian
hess <- fit$hess
## fisher_info <- -solve(hess)
fisher_info <- solve(hess)
se.vec <- sqrt(diag(fisher_info))
se.vec

## extract results
p = fit$par
p.upper = p + 2*se.vec
p.lower = p - 2*se.vec
hat = c(exp(p[-length(p)]), p[length(p)])
hat.upper = c(exp(p.upper[-length(p)]), p.upper[length(p)])
hat.lower = c(exp(p.lower[-length(p)]), p.lower[length(p)])
out = cbind(hat, hat.lower, hat.upper)
rownames(out) = c(paste0("coh", cohorts), "beta", "b")
coefs <- print(round(out,2))

# write results to working directory
sink("coefficients.txt", append = TRUE)
coefs
sink()

##           hat hat.lower hat.upper
## coh1915 80.78     79.19     82.40
## coh1916 79.63     78.02     81.27
## coh1917 79.15     77.54     80.80
## coh1918 77.97     76.37     79.59
## coh1919 78.09     76.48     79.73
## beta     0.13      0.12      0.14
## b       -0.04     -0.06     -0.02

## compare to OLS
lm1 <- lm(y ~ covar + as.factor(cohort), data = A)

# write summary to working directory
sink("ols.summary.txt", append = TRUE)
summary(lm1)
sink()

## -------------- edu -----------
##                       Estimate Std. Error t value Pr(>|t|)
## (Intercept)           79.83078    0.27172 293.800  < 2e-16 ***
## covar                  0.12047    0.02209   5.454 5.09e-08 ***
## as.factor(cohort)1916 -0.99186    0.20004  -4.958 7.28e-07 ***
## as.factor(cohort)1917 -1.92877    0.19819  -9.732  < 2e-16 ***
## as.factor(cohort)1918 -2.99152    0.19572 -15.284  < 2e-16 ***
## as.factor(cohort)1919 -3.74162    0.19626 -19.064  < 2e-16 ***

## b = -.04  vs ols 0.12

## -4% on mortality should change life expectancy by about .04*.1 * 80 = [1] 0.32 years
## here we're getting about 1/2 to 1/3 of this.
## truncation effect should be about

A[, txtdensity(y)]
## mode 80 with beta = 1/10, and truncation at about 68 and 88

my.alpha = getAlpha(M = 80, beta = .1)
get.trunc.mean.gomp <- function(alpha, beta, l, u)
{
    x = 0:110
    hx = alpha * exp(beta*x)
    Hx = cumsum(hx)
    lx = c(1, exp(-Hx))
    sum(lx)
    dx = -diff(lx)
    s <- x %in% l:u ## approximate truncation for these cohorts
    sum((dx*x)[s])/sum(dx[s])
}
get.trunc.mean.gomp(my.alpha, .1, 68, 88)
## [1] 78.3002
get.trunc.mean.gomp(my.alpha*(.96), .1, 68, 88)
## [1] 78.41795
## so truncated mean goes up by .11
## very close to OLS


