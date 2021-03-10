#setwd("~/GDrive/demographer/github/censoc_mort")
setwd("~/Documents/MLE")
library(data.table)
source("./gompertz_functions.R")

## 1. Write the mll function 

mll.gomp.multicohort.alpha.fe <- function(p, # as many alphas as cohorts, and 1 beta
                                          A) # data matrix with y, u, l, and covariates including cohort
{
  ## (1) get parameters alpha, beta, and B
  ## These will small in number, we'll convert later to one element per individual.
  cohorts = names(table(A$cohort))  ## in data matrix A, each birth year is a name is a cohort
  k <- length(cohorts)
  M = exp(p[1:k]); names(M) = cohorts
  beta = exp(p[k+1])
  ## convert M to alpha
  alpha = getAlpha(M, beta); names(alpha) = cohorts
  
  if(length(p) > k) ## if covariates exist
    B = p[(k+1):length(k)]
  
  ## (2) build rate.vec which has the combined effect of cohort and covariates, 1 element for each indiv
  alpha.vec <- alpha[paste(A$cohort)] ## this is now has one alpha for each individual
  ## (if covariates exist, could do something like next line ...)
  ##     covar_effect.vec <- exp(B * log_wages)
  covar_effect.vec <- 1 ## assume no covariates
  rate.vec = alpha.vec * covar_effect.vec
  y = A$y
  u.vec = A$u
  l.vec = A$l
  
  ## (3) likelihood
  ## Numerator has difference in cumulative probabilities, 1 year apart. This works with data that is
  ## exact to the year of death like we have in CenSoc.
  num = pgompertz(y+1, shape = beta, rate = rate.vec)-
    pgompertz(y, shape = beta, rate = rate.vec)
  denom = pgompertz(u.vec, shape = beta, rate = rate.vec)-
    pgompertz(l.vec, shape = beta, rate = rate.vec)
  R = num/denom
  minusloglik = -sum(log(R))
  return(minusloglik)
}


## 2. Assign values required for mll.gomp.multicohort.alpha.fe function

## scalars
beta_true = 1/10

## the following variables just have 5 elements
cohorts = 1915:1919
first_dyear = 1988
last_dyear = 2005
uppers <- last_dyear - cohorts
lowers <- first_dyear - cohorts
modes_true <- c(86, 85, 84, 83, 82) ## listed in reverse order compared to Josh's example; derived from
## modal values found using Mode function (see bottom of script) on
## each birth year cohort of the demo data
names(modes_true) <- cohorts


## 3. Serge's kludgy way to build dataframe "A" using censoc_demo data

#dt = fread("./censoc_numident_demo_v1/censoc_numident_demo_v1.csv")
#dt <- fread("/data/josh/CenSoc/censoc_data/censoc_dmf_v1/censoc_dmf_v1.csv")
dt <- fread("/data/josh/CenSoc/censoc_data/censoc_numident_v2/censoc_numident_v2.csv")

getEstimates <- function()
{    
  cohrt <- dt[byear == 1915]  ## subset a birth cohort
  #Mode(cohrt$death_age)  ## mode = 88, based on full subset (currently commented out; see Mode function at line 165)
  cohrt1 <- cohrt[sample(nrow(cohrt), 10000), ]  ## randomly choose cases
  cohrt1 <- cohrt1[, c("byear","death_age")]  ## subset all rows in just two columns
  #Mode(cohrt1$death_age)  ## mode = 82, based on random sample of 1,000
  cohrt1$u = 90  ## figure from five lines above 
  cohrt1$l = 72  ## figure from five lines above
  cohrt1$M = 86
  names(cohrt1)[2] = "y"
  names(cohrt1)[1] = "cohort"
  
  cohrt <- dt[byear == 1916]
  #Mode(cohrt$death_age)  ## mode = 83
  cohrt2 <- cohrt[sample(nrow(cohrt), 10000), ]
  cohrt2 <- cohrt2[, c("byear","death_age")]
  #Mode(cohrt2$death_age)  ## mode = 82
  cohrt2$u = 89
  cohrt2$l = 71
  cohrt2$M = 85
  names(cohrt2)[2] = "y"
  names(cohrt2)[1] = "cohort"
  
  cohrt <- dt[byear == 1917]
  #Mode(cohrt$death_age)  ## mode = 84
  cohrt3 <- cohrt[sample(nrow(cohrt), 10000), ]
  cohrt3 <- cohrt3[, c("byear","death_age")]
  #Mode(cohrt3$death_age)  ## mode = 82
  cohrt3$u = 88
  cohrt3$l = 70
  cohrt3$M = 84
  names(cohrt3)[2] = "y"
  names(cohrt3)[1] = "cohort"
  
  cohrt <- dt[byear == 1918]
  #Mode(cohrt$death_age)  ## mode = 83
  cohrt4 <- cohrt[sample(nrow(cohrt), 10000), ]
  cohrt4 <- cohrt4[, c("byear","death_age")]
  #Mode(cohrt4$death_age)  ## mode = 82
  cohrt4$u = 87
  cohrt4$l = 69
  cohrt4$M = 83
  names(cohrt4)[2] = "y"
  names(cohrt4)[1] = "cohort"
  
  cohrt <- dt[byear == 1919]
  #Mode(cohrt$death_age)  ## mode = 82
  cohrt5 <- cohrt[sample(nrow(cohrt), 10000), ]
  cohrt5 <- cohrt5[, c("byear","death_age")]
  #Mode(cohrt5$death_age)  ## mode = 82 
  cohrt5$u = 86
  cohrt5$l = 68
  cohrt5$M = 82
  names(cohrt5)[2] = "y"
  names(cohrt5)[1] = "cohort"
  
  
  ## combine five cohorts into a single dataframe "A"
  A <- rbind(cohrt1, cohrt2, cohrt3, cohrt4, cohrt5)
  
  ## check to see that mean rises with cohort
  A[, mean(y), by = cohort]
  
  ## truncate the cohort to upper and lower truncation points
  At = A[l < y & y < u]
  
  ## floor so it is like "death_age", so 76.89 becomes 76
  At[, y := floor(y)]
  
  ## starting values
  p.start = log(c(rep(80, length(cohorts)),
                  beta_true * .8))
  names(p.start) <- c(paste0("coh", cohorts), "beta")
  
  ## run optimizer
  fit = optim(par = p.start, fn = mll.gomp.multicohort.alpha.fe,
              A = At,
              control = list(maxit = 10000))
  
  hat <- round(exp(fit$par), digits = 3)
  true <- c(modes_true, beta_true)
  
  print(cbind(hat, true))
  myDF <- cbind(hat, true)
  
  ## run the line below just once at the beginning to create the CSV file to which the results from
  ## subsequent replications of the getEstimates function will be written to
  #write.csv(myDF, file = "myDF.csv")
  
  write.table(myDF, "myDF.csv", sep = ",", col.names = !file.exists("myDF.csv"), append = T)
}


replicate(10, getEstimates())




# ## functions to obtain modal values

# ## option 1
# Mode <- function(x) {
#     ux <- unique(x)
#     ux[which.max(tabulate(match(x, ux)))]
# }
# 
# ## option 2
# Mode <- function(x, na.rm = FALSE) {
#     if(na.rm){
#         x = x[!is.na(x)]
#     }
#     ux <- unique(x)
#     return(ux[which.max(tabulate(match(x, ux)))])
# }
