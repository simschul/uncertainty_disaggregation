#'
#'
#'
#' @author Simon Schulte
#' Date: 2021-11-29 10:18:16
#'
#' Content:
#'


############################################################################## #
##### load packages ############################################################
############################################################################## #

library(data.table)
library(tidyverse)
library(units)
library(ggforce)
library(my.utils)

############################################################################## #
##### settings #################################################################
############################################################################## #
options("datatable.print.class" = TRUE)
theme_set(theme_bw())

############################################################################## #
##### load data #############################################################
############################################################################## #

## 1. the standard dirichlet (from gtools) =====================================
## see: https://stat.ethz.ch/pipermail/r-help/2000-December/009561.html


logD <- function(a) {
  sum(lgamma(a)) - lgamma(sum(a))
}


######################################################
ddirichlet <- function(x, alpha)
  ## probability density for the Dirichlet function, where x=vector of
  ## probabilities
  ## and (alpha-1)=vector of observed samples of each type
  ## ddirichlet(c(p,1-p),c(x1,x2)) == dbeta(p,x1,x2)
{
  s <- sum((alpha - 1) * log(x))
  exp(sum(s) - logD(alpha))
}


######################################################
n <- 1000
a <- c(1, 1, 1)

rdirichlet <- function(n, a)
  ## pick n random deviates from the Dirichlet function with shape
  ## parameters a
{
  l <- length(a)
  
  x <- matrix(rgamma(l * n, a), ncol = l, byrow = TRUE)
  
  sm <- x %*% rep(1, l)
  
  x / as.vector(sm)
  
}

rdirichlet(1000, c(1, 1, 1)) %>% boxplot


# 2. Gernerlised dirichlet with MGML package ===============
library(MGML)

rgdirmn

# ===================


#
# Generation of `n` random variates.
#
rGDirichlet <- function(n, alpha, beta) {
  d <- max(length(alpha), length(beta))
  Y <-
    matrix(rbeta(d * n, alpha, beta), d) # Each sample is in a column of Y
  q <- rep(0, n)
  for (i in seq_len(d)) {
    # Convert `Y` to `X` (in place)
    Y[i,] <- Y[i,] * (1 - q)
    q <- q + Y[i,]
  }
  return(Y)
}
#
# Multivariate moments.
#
mGDirichlet <- function(r, a, b, use.log = TRUE) {
  d <- sum(r) - cumsum(r)
  z <- sum(lgamma(a + b) + lgamma(a + r) + lgamma(b + d) -
             (lgamma(a) + lgamma(b) + lgamma(a + b + r + d)))
  if (!isTRUE(use.log))
    z <- exp(z) # The `r` moment of a GDirichlet(a, b) variate
  z
}



# 3. from stackexchange ========================================================
# from : https://stats.stackexchange.com/questions/534411/how-to-generate-data-from-a-generalized-dirichlet-distribution

# Multivariate moments.
#
rGDirichlet <- function(r, a, b, use.log = TRUE) {
  d <- sum(r) - cumsum(r)
  z <- sum(lgamma(a + b) + lgamma(a + r) + lgamma(b + d) -
             (lgamma(a) + lgamma(b) + lgamma(a + b + r + d)))
  if (!isTRUE(use.log))
    z <- exp(z) # The `r` moment of a GDirichlet(a, b) variate
  z
}

# Random draws.
#
set.seed(17)
alpha <- 0.1 + round(rexp(4, 1 / 3), 1)            # alpha parameters
beta <- 0.1 + round(rexp(length(alpha), 1 / 3), 1) # beta parameters
X <- rGDirichlet(1e4, alpha, beta)
#
# Check some low moments for agreement with theory.
#
m <- 3 # Maximum total order of the moments
R <-
  as.matrix(do.call(expand.grid, lapply(seq_along(alpha), function(i)
    0:m)))
R <- R[rowSums(R) <= m,]
rownames(R) <- apply(R, 1, paste0, collapse = "")

Log.Moments <- apply(R, 1, function(r) {
  c(Theoretical = mGDirichlet(r, alpha, beta),
    Empirical = log(mean(apply(X, 2, function(x)
      prod(x ^ r)))))
})
plot(
  t(Log.Moments),
  pch = 21,
  bg = hsv(rowSums(R) / (m + 1), .8, 1, .25),
  cex = 1.5,
  main = paste0("Logs of Moments Up to Total Order ", m)
)
abline(0:1, col = "Gray", lwd = 2)
mtext(paste0(
  "alpha = (",
  paste(signif(alpha, 2), collapse = ", "),
  ")  beta = (",
  paste(signif(beta, 2), collapse = ", "),
  ")"
), side = 3)



# 4. from Wong 1998 ===============================================
# from wong 1998
# pseudocode:
# xl = rbeta(cq, ill)
# sum = Xl
# f o r ( j in 2 : k )
# {
#   xj. = rbeta(~j, flj) • (1 -sum = sum + Xj
# }



rgdirichlet <- function(n, alpha, beta) {
  k <- length(alpha)
  x <- matrix(0, nrow = n, ncol = k)
  x[, 1] <- rbeta(n, alpha[1], beta[1])
  sum <- x[, 1]
  for (i in 2:k) {
    x[, i] <- rbeta(n, alpha[i], beta[i]) * (1 - sum)
    sum <- sum + x[, i]
  }
  return(x)
}

alpha <- c(0.3, 0.5, 0.2)
beta <- c(0.3, 0.04, 0.01)
n <- 100
set.seed(12)

test <- rgdirichlet(n, alpha, beta)
apply(test, 1, sum)

plot(test[1, ], type = 'lines', ylim = c(0, 1))
for (i in 2:nrow(test)) {
  lines(test[i, ], col = i)
}

boxplot(test)


# 5. from Plessis 2010 ========================================================
# see apendix 2

## The way to go!!!
mu <- c(0.3, 0.5, 0.2)
u <- c(0.3, 0.04, 0.01)
n <- 100
set.seed(12)

rgdirichlet <- function(n, mu, u) {
  alpha <- (mu / u) ^ 2
  beta <- mu / (u) ^ 2
  k <- length(alpha)
  x <- matrix(0, nrow = n, ncol = k)
  for (i in 1:k) {
    x[, i] <- rgamma(n, shape = alpha[i], rate = beta[i])
  }
  return(x / rowSums(x)) 
}

boxplot(x)
plot_dirichlet(x, mu = mu, u = u, type = 2, xlab = '', ylab = '')

test <- rgdirichlet(n, mu, u)
apply(test, 1, sum)


plot_dirichlet <- function(x, mu = NULL, u = NULL, type, ...) {
  
  # type line plot ====================================== 
  if (type == 1) {
    plot(x[1, ],
         type = 'lines',
         ylim = c(0, 1),
         col = 'grey80', ...)
    
    for (i in 2:nrow(x)) {
      lines(x[i, ], col = 'grey80')
    }
    
    if (!is.null(mu)) lines(mu, col = 'red', lwd = 3)
    if (!is.null(u)) {
      lines(mu + u,
            col = 'red',
            lwd = 3,
            lty = 2)
      lines(mu - u,
            col = 'red',
            lwd = 3,
            lty = 2)  
    }  
  }
  
  # type 95 CI =========================================
  if (type == 2) {
    percent <- x %>%
      apply(2, quantile, probs = c(0.025, 0.5, 0.975))
    
    plot(
      percent[1, ],
      type = 'point',
      ylim = c(0, 1),
      col = 'grey30',
      pch = 25, ...
    )
    points(percent[2, ], col = 'grey30', pch = 8)
    points(percent[3, ], col = 'grey30', pch = 24)
    
    if (!is.null(mu)) points(mu, col = 'red', pch = 8)
    
    if (!is.null(u)) {
      points(mu + u,
             col = 'red',
             pch = 24)
      points(mu - u,
             col = 'red',
             pch = 25)  
    }
  }
  
  # type boxplot =====================================
  if (type == 3) {
    boxplot(x, ...)
    
    if (!is.null(mu)) points(mu, col = 'red', lwd = 3)
    
    if (!is.null(u)) {
      points(mu + u,
             col = 'red',
             lwd = 3,
             pch = 4)
      points(mu - u,
             col = 'red',
             lwd = 3,
             pch = 4)  
    } 
  }
}

plot_dirichlet(test, mu = mu, u = u, type = 2, xlab = '', ylab = '')



# _b) test function -----------------------------------------------------------

# __i. random uncertainties for disaggreagte =======================================================

library(gtools)
N <- 25 # number of simulations 
n <- 1000 # number of draws per simulation
k <- 5 # number of shares

# draw shares
mumat <- rdirichlet(N, rep(1, k))

# draw uncertainty related to shares (between 0 and 100 %)
umat <- mumat * matrix(runif(k*N), ncol = k)

# draw n times from generalised dirichlet
draws <- array(dim = c(n, k, N))

for (i in 1:N) {
  draws[,,i] <- rgdirichlet(n, mu = mumat[i,], u = umat[i,])
}

# plot results
par(mfrow = c(5,5), mar = c(2,0.2,0.2,2))
for (i in 1:N) {
  plot_dirichlet(draws[,,i], mu = mumat[i,], u = umat[i,], type = 1, 
                 xlab = '', ylab = '')
}

# __ii. fixed uncertainties for all dissagreagete values ===================================================================

N <- 25 # number of simulations 
n <- 200 # number of draws per simulation
k <- 5 # number of shares
u <- 0.05 # uncertainty of shares/dissagregate values

# draw shares
mumat <- rdirichlet(N, rep(1, k))

# used fixed (relative) uncertainty of shares 
umat <- mumat * u

# draw n times from generalised dirichlet
draws <- array(dim = c(n, k, N))

for (i in 1:N) {
  draws[,,i] <- rgdirichlet(n, mu = mumat[i,], u = umat[i,])
}

# plot results
par(mfrow = c(5,5), mar = c(2,0.2,0.2,2))
for (i in 1:N) {
  plot_dirichlet(draws[,,i], mu = mumat[i,], u = umat[i,], type = 2, 
                 xlab = '', ylab = '')
}


# __iii. mix of diff fixed uncertainties disaggregate values ===================================================================

N <- 25 # number of simulations 
n <- 200 # number of draws per simulation
k <- 5 # number of shares
u <- c(0.01, 0.3) # uncertainty of shares/dissagregate values

# draw shares
mumat <- rdirichlet(N, rep(1, k))

# used fixed (relative) uncertainty of shares 
umat <- mumat * sample(u, size = length(mumat), 
                       prob = c(0.5, 0.5), replace = TRUE)

# draw n times from generalised dirichlet
draws <- array(dim = c(n, k, N))

for (i in 1:N) {
  draws[,,i] <- rgdirichlet(n, mu = mumat[i,], u = umat[i,])
}

# plot results
par(mfrow = c(5,5), mar = c(2,0.2,0.2,2))
for (i in 1:N) {
  plot_dirichlet(draws[,,i], mu = mumat[i,], u = umat[i,], type = 2, 
                 xlab = '', ylab = '')
}




# 6. fitting the parameters using the modcmfitr package =========================
# code from: https://cran.r-project.org/web/packages/modcmfitr/vignettes/modcmfitrOverview.pdf
# does not seem to make better fits

library(modcmfitr)
library(nloptr)
set.seed(12345)
Outcomes <- c("A","B","C")
RawData <- matrix(data = c(0.43, 0.55, 0.65,
                           0.16, 0.27, 0.46,
                           0.03, 0.18, 0.23
),ncol=3,byrow=TRUE)


mu <- c(0.3, 0.5, 0.2)
u <- c(0.3, 0.04, 0.01)

RawData <- matrix(data = c(0, 0.3, 0.6,
                           0.46, 0.5, 0.54,
                           0.19, 0.2, 0.21
),ncol=3,byrow=TRUE)

SearchParams <- c(10000,1000) # number of iterations, max number of searches
# (set to 100 here; you probably want more).
ModCMorCM <- 0 # flag to determine whether fitting mCM or CM distribution
Quantiles <- c(0.025,0.5,0.975) # example here is 95% credibility limits and median.
mCM <- fitModCM(Outcomes, RawData, SearchParams, ModCMorCM, Quantiles)




Z <- mCM[1:2,1:4]
rownames(Z) <- c("Z1","Z2")
print(Z)

mCMSamples <- rModCM(100,Z)
colnames(mCMSamples) <- c("A", "B", "C")
print(mCMSamples)

boxplot(mCMSamples, ylim = c(0,1))
points(RawData[,1], col = 'red')
points(RawData[,2], col = 'red')
points(RawData[,3], col = 'red')





