## Leslie's questions
## 1. Should you read Monica's dissertation? Where should we start?

## 2. 1 cohort or a few cohorts, are we agreed on that?
## (josh's answer, maybe add a just sex for 1918 for even simpler ...)

## Maria's questions
## 1. What kind of stan should we use rstanr etc? (for our work, and for sharing with the world?)

## Jordan
## 1. What's the exact set-up of this experiment? (specs, please)

## Nathan
## 1. Global Q. : other users, what are they doing?
## (lit search, also just reaching out by email ..?)

## Josh's response to Leslie2 and Jordan1

## Let's start small and work up:

## (0**) Just male-female for cohort of 1918 without weights (no weights, for now).
## (1) male-female and income for 1918
## (2) male-female and income for all of the cohorts below (1915:1920 -- I think)

## One thing that I think is worth doing, is reparameterizing Gompertz in terms of modal age at death.
## Monica has a blog post on this



## Some useful common resources:

## https://www.monicaalexander.com/posts/2018-02-15-gompertz/

## library(flexsurv)

## we should also have our own functions for parameterizing with
## M instead of alpha ... josh has written these and will share upon being reminded.

## josh should share MLE stuff.










## ------- old code below is from discussion with Achim ----

## Quick extraction and regression example using
## censoc_numident_demo_v1.csv
## data downloaded from censoc.berkeley.edu

## josh and achim (nov 5, 2020)


library(data.table) ## just for some data manipulation (different approach than tidyverse, but unfortunately not base-R)

dt = fread("~/Downloads/censoc_numident_demo_v1/censoc_numident_demo_v1.csv") ## read in the data easily
## recode annual wage income in 1939 for missing values
dt[, inc := incwage]
dt[incwage == 0, inc := NA]             # 0 doesn't mean no earnings, it means missing for many diff. reasons
dt[incwage == 999999, inc := NA]        # 999999 means not eligible for question, for many diff. reasons
library(txtplot)                        # ascii plotting (old school)
## check to make sure it looks right
txtdensity(dt[!is.na(inc)]$inc)         # unimodal, skew right
txtdensity(log(dt[!is.na(inc)]$inc))    # with logs, fairly symmetric and "normal"


## now we run some regressions (ignoring truncation in the estimation,
## although we do a bit to compensate, such as putting in byear fixed
## effects -- and making the truncation explicit in the subset
## command for "dyear %in%"

m <- lm(death_age ~ inc + sex + as.factor(byear),
        data = dt,
        subset = dyear %in% 1988:2005 & byear %in% 1915:1920)
## summary(m)
##                        Estimate Std. Error t value Pr(>|t|)
## (Intercept)          81.9930481  0.1699769 482.378  < 2e-16 ***
## inc                   0.0003816  0.0001238   3.082  0.00206 **
## sexMale              -1.2464907  0.1082242 -11.518  < 2e-16 ***
## as.factor(byear)1916 -1.0395974  0.1798306  -5.781 7.66e-09 ***
## as.factor(byear)1917 -1.7444560  0.1788659  -9.753  < 2e-16 ***
## as.factor(byear)1918 -2.7287784  0.1809039 -15.084  < 2e-16 ***
## as.factor(byear)1919 -3.6390555  0.1856717 -19.599  < 2e-16 ***
## as.factor(byear)1920 -4.3745215  0.1923473 -22.743  < 2e-16 ***
## Note: the sample size is about 10,000 people


## focus on one birth cohort
m.1918 <- lm(death_age ~ inc + sex,
        data = dt,
        subset = dyear %in% 1988:2005 & byear %in% 1918)
m.1918.sex.only <- lm(death_age ~ sex,
        data = dt,
        subset = dyear %in% 1988:2005 & byear %in% 1918)


library(lfe)                            # fixed effects felm() package
m.fe <- felm(death_age ~ inc + sex | as.factor(byear),
        data = dt,
        subset = dyear %in% 1988:2005 & byear %in% 1915:1918)


library(stargazer)                      # nice display of regression results package

stargazer(m, m.fe, m.1918, type = "text")

## =========================================================================================
##                                              Dependent variable:
##                      --------------------------------------------------------------------
##                                                   death_age
##                                 OLS                  felm                  OLS
##                                 (1)                   (2)                  (3)
## -----------------------------------------------------------------------------------------
## inc                          0.0004***             0.0003**              -0.0001
##                              (0.0001)              (0.0001)              (0.0003)

## sexMale                      -1.246***             -1.402***            -1.575***
##                               (0.108)               (0.131)              (0.257)

## as.factor(byear)1916         -1.040***
##                               (0.180)

## as.factor(byear)1917         -1.744***
##                               (0.179)

## as.factor(byear)1918         -2.729***
##                               (0.181)

## as.factor(byear)1919         -3.639***
##                               (0.186)

## as.factor(byear)1920         -4.375***
##                               (0.192)

## Constant                     81.993***                                  79.751***
##                               (0.170)                                    (0.261)

## -----------------------------------------------------------------------------------------
## Observations                   9,364                 6,442                1,624
## R2                             0.094                 0.054                0.023
## Adjusted R2                    0.093                 0.054                0.022
## Residual Std. Error      5.038 (df = 9356)     5.024 (df = 6436)    5.041 (df = 1621)
## F Statistic          138.040*** (df = 7; 9356)                   19.130*** (df = 2; 1621)
## =========================================================================================
## Note:                                                         *p<0.1; **p<0.05; ***p<0.01


## include sex-only model for our MLE experiments
#
stargazer(m, m.1918, m.1918.sex.only, type = "text")
## ================================================================================================
##                                                  Dependent variable:
##                      ---------------------------------------------------------------------------
##                                                       death_age
##                                 (1)                      (2)                      (3)
## ------------------------------------------------------------------------------------------------
## inc                          0.0004***                 -0.0001
##                              (0.0001)                  (0.0003)

## sexMale                      -1.246***                -1.575***       ---->>>   -1.396***  <<--
##                               (0.108)                  (0.257)                  (0.189)

## as.factor(byear)1916         -1.040***
##                               (0.180)

## as.factor(byear)1917         -1.744***
##                               (0.179)

## as.factor(byear)1918         -2.729***
##                               (0.181)

## as.factor(byear)1919         -3.639***
##                               (0.186)

## as.factor(byear)1920         -4.375***
##                               (0.192)

## Constant                     81.993***                79.751***                79.516***
##                               (0.170)                  (0.261)                  (0.136)

## ------------------------------------------------------------------------------------------------
## Observations                   9,364                    1,624                    2,790
## R2                             0.094                    0.023                    0.019
## Adjusted R2                    0.093                    0.022                    0.019
## Residual Std. Error      5.038 (df = 9356)        5.041 (df = 1621)        5.000 (df = 2788)
## F Statistic          138.040*** (df = 7; 9356) 19.130*** (df = 2; 1621) 54.362*** (df = 1; 2788)
## ================================================================================================
## Note:                                                                *p<0.1; **p<0.05; ***p<0.01

## for birth cohort of 1918
## OLS is estimating: beta = d(Y | Y %in% 70:87) / d X
## Achim-method : beta_a = d(Y | Y <= 87) / dX

## could we invert the problem? Not completely, because we're dealing with survival ... and that affects the way we think about it.

## Achim-method-inverted? : beta_a = d(Y | Y >= 70) / dX

## What do we really want? Something like

## beta_true = d(Y | Y >= 65) / dX

## for example, in parametric MLE

## L = f(y | alpha(X), beta(X)) / [F(u | alpha(X), beta(X)) - F(l | alpha(X), beta(X))]

## then with alpha.hat and beta.hat and X, I would estimate

## d(Y(X) | Y >= 65) / dX






## Last: look at density of deaths for cohort of 1918


dt[dyear %in% 1988:2005 & byear %in% 1918,  txtdensity(death_age)]
dt[sex == "Male" & dyear %in% 1988:2005 & byear %in% 1918,  txtdensity(death_age)]
dt[sex == "Female" & dyear %in% 1988:2005 & byear %in% 1918,  txtdensity(death_age)]

## dt[dyear %in% 1988:2005 & byear %in% 1918,  txtdensity(death_age)]
## 0.07 +-----+-------------+--------------+-------------+--------+
##      |                                     ***********         |
##      |                                   ***         ***       |
## 0.06 +                               *****             **      +
##      |                      **********                  **     |
##      |                  *****                            **    |
## 0.05 +             ******                                 *    +
##      |          ****                                      **   |
## 0.04 +      *****                                          *   +
##      |     **                                              **  |
##      |    **                                                *  |
## 0.03 +   **                                                 *  +
##      |  **                                                     |
##      |  *                                                      |
## 0.02 +-----+-------------+--------------+-------------+--------+
##           70            75             80            85
## NULL
## > dt[sex == "Male" & dyear %in% 1988:2005 & byear %in% 1918,  txtdensity(death_age)]
## 0.06 +-----+-------------+--------------+-------------+--------+
##      |                           *******************           |
##      |           ******** ********                 ****        |
##      |         ***      ***                           **       |
## 0.05 +        **                                       **      +
##      |       **                                         *      |
##      |      **                                           *     |
##      |     **                                            **    |
## 0.04 +    **                                              *    +
##      |    *                                               *    |
##      |   **                                                *   |
##      |   *                                                 *   |
## 0.03 +  *                                                   *  +
##      |  *                                                   *  |
##      +-----+-------------+--------------+-------------+--------+
##           70            75             80            85
## NULL
## > dt[sex == "Female" & dyear %in% 1988:2005 & byear %in% 1918,  txtdensity(death_age)]
## 0.08 +-----+-------------+--------------+-------------+--------+
##      |                                     ************        |
## 0.07 +                                    **          ***      +
##      |                                  ***             **     |
## 0.06 +                                ***                *     +
##      |                         ********                  **    |
## 0.05 +                   *******                          **   +
##      |                 ***                                 *   |
## 0.04 +               ***                                   **  +
##      |             ***                                      *  |
## 0.03 +     ********                                            +
##      |    **                                                   |
## 0.02 +   **                                                    +
##      |  **                                                     |
##      +-----+-------------+--------------+-------------+--------+
##           70            75             80            85
## NULL
## >
