## parameterize gompertz with M (mode) and beta (slope)

library(flexsurv)

getMode <- function(alpha, beta)
{
    M = (1/beta) * log(beta/alpha)
    return(M)
}
getAlpha <- function(M, beta)
{
    alpha = beta/ exp(M*beta)
    return(alpha)
}

dgomp_mode <- function(x, M, beta, ...)
{
    alpha = getAlpha(M, beta)
    dgompertz(x, shape = beta, rate = alpha, ...)
}

pgomp_mode <- function(q, M, beta, ...)
{
    alpha = getAlpha(M, beta)
    pgompertz(q, shape = beta, rate = alpha, ...)
}

rgomp_mode <- function(n, M, beta)
{
    alpha = getAlpha(M, beta)
    rgompertz(n, shape = beta, rate = alpha)
}

