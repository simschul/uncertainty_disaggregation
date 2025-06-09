### distributions.R:  distribution and parameter estimation functions for the lcmix package

### UNIVARIATE DISTRIBUTION DEFINITIONS

## The categorical distribution.

# density function
dcategorical <- function(x, prob, log=FALSE)
  if(log) log(prob)[x] else prob[x]

# random generation function
rcategorical <- function(n, prob)
  sample.int(length(prob), n, TRUE, prob)

## The univariate Pearson Type VII distribution (alpha, 1) parameterization

# probability density function
dpvii <- function(x, mean, scale, shape, log=FALSE)
{
  shape.prime = shape + 0.5
  Delta = (x - mean)^2
  delta = Delta / scale
  scale.prime = 1 + 0.5*delta

  retn = (- 0.5*LC_LOG2PI - 0.5*log(scale) - lgamma(shape)
          + lgamma(shape.prime) - shape.prime*log(scale.prime) )

  if(log) retn else exp(retn)
}

# random generation function
rpvii <- function(n, mean, scale, shape)
  rnorm(n, mean, sqrt(scale / rgamma(n, shape, 1)))

## The Weibull, shape-decay parameterization ("WeiSD") distribution.
## If X ~ WeiSD(shape, decay) then:
## f(x) = shape * decay * x^(shape-1) * exp(-decay * x^shape)
## F(x) = 1 - exp(-decay * x^shape)
## S(x) = exp(-decay * x^shape)
## H(x) = decay * x^shape
## h(x) = decay * shape * x^(shape-1)
## This implies that the "shape" parameters for WSR and the standard R
## parameterization of the Weibull distribution are the same, while
## using the "scale" parameter in the standard R parameterization we
## have decay = scale^(-shape), or equivalently, scale = decay^(-1/shape).

# probability density function
dweisd <- function(x, shape, decay, log=FALSE)
  dweibull(x, shape=shape, scale=decay^(-1/shape), log=log)

# cumulative density function
pweisd <- function(q, shape, decay, lower.tail=TRUE, log.p=FALSE)
  pweibull(q, shape=shape, scale=decay^(-1/shape),
           lower.tail=lower.tail, log.p=log.p)

# quantile function
qweisd <- function(p, shape, decay, lower.tail=TRUE, log.p=FALSE)
  qweibull(p, shape=shape, scale=decay^(-1/shape),
           lower.tail=lower.tail, log.p=log.p)

# random generation function
rweisd <- function(n, shape, decay)
  rweibull(n, shape=shape, scale=decay^(-1/shape))

# expected value
eweisd <- function(shape, decay) decay^(-1/shape) * gamma(1 + 1/shape)

# variance
vweisd <- function(shape, decay)
  decay^(-2/shape) * gamma(1 + 2/shape) - eweisd(shape, decay)^2

# median
medweisd <- function(shape, decay) decay^(-1/shape) * log(2)^(1/shape)

# for X ~ WeiSD(shape, decay), return the rth raw moment of x, i.e. E(X^r)
rmomweisd <- function(r, shape, decay) decay^(-r/shape) * gamma(1 + r/shape)

### MULTIVARIATE DISTRIBUTION DEFINITIONS

## Multivariate exponential distribution using normal copula

# density function
dmvexp <- function(x, rate, corr=diag(ncol(x)), log=FALSE)
{
  # munge data PRN
  if(is.data.frame(x))
    x = as.matrix(x)
  else if(is.vector(x))
    x = matrix(x, nrow=1)
  ## gather statistics

  D = ncol(x)

  cdf = sapply(1:D, function(d) pexp(x[,d], rate[d]))
  normcdf = qnorm(cdf)
  normcdf = forceRange(normcdf, LC_MINSTDNORM, LC_MAXSTDNORM)
  # deal with outliers so small (large) that their effective CDF is 0 (1); this represents a robustification approach

  ## calculate copula contribution to log-PDF

  res.copula.add = dmvnorm(normcdf, cov=corr, log=TRUE)
  res.copula.sub = rowSums(dnorm(normcdf, log=TRUE))
  res.copula = res.copula.add - res.copula.sub

  ## calculate marginal contribution to log-PDF

  res.data = rowSums(sapply(1:D, function(d)
    dexp(x[,d], rate[d], log=TRUE)))

  ## final calculations and return

  retn = res.copula + res.data

  if(log) retn else exp(retn)
}

# random generation function
rmvexp <- function(n, rate=1, corr=diag(length(rate)))
{
  ## extract parameters, do sanity checks, deal with univariate case

  if(!is.matrix(corr) || !isSymmetric(corr))
    stop("'corr' must be a symmetric matrix")
  D = ncol(corr)

  Dr = length(rate)
  if(Dr > D)
    warning("'rate' longer than width of 'corr', truncating to fit")
  if(Dr != D)
    rate = rep(rate, length.out=D)

  if(D == 1) rexp(n, rate)

  ## generate standard multivariate normal matrix, convert to CDF

  Z = rmvnorm(n, cov=corr)
  cdf = pnorm(Z)

  ## convert to exp, return

  sapply(1:D, function(d) qexp(cdf[,d], rate[d]))
}

## Multivariate gamma distribution using normal copula (TO DO:  figure out an elegant way to allow for scale parameter; it may be as simple as throwing the "scale=1/rate" default into the function calls and then handing scale off to the dgamma() and pgamma() functions along with shape and rate ...)

# density function
dmvgamma <- function(x, shape, rate, corr=diag(ncol(x)), log=FALSE)
{
  # munge data PRN
  if(is.data.frame(x))
    x = as.matrix(x)
  else if(is.vector(x))
    x = matrix(x, nrow=1)

  ## gather statistics

  D = ncol(x)

  cdf = sapply(1:D, function(d) pgamma(x[,d], shape[d], rate[d]))
  normcdf = qnorm(cdf)
  normcdf = forceRange(normcdf, LC_MINSTDNORM, LC_MAXSTDNORM)
  # deal with outliers so small (large) that their effective CDF is 0 (1); this represents a robustification approach

  ## calculate copula contribution to log-PDF

  res.copula.add = dmvnorm(normcdf, cov=corr, log=TRUE)
  res.copula.sub = rowSums(dnorm(normcdf, log=TRUE))
  res.copula = res.copula.add - res.copula.sub

  ## calculate marginal contribution to log-PDF

  res.data = rowSums(sapply(1:D, function(d)
    dgamma(x[,d], shape[d], rate[d], log=TRUE)))

  ## final calculations and return

  retn = res.copula + res.data

  if(log) retn else exp(retn)
}

# random generation function
rmvgamma <- function(n, shape=1, rate=1, corr=diag(length(shape)))
{
  ## extract parameters, do sanity checks, deal with univariate case

  if(!is.matrix(corr) || !isSymmetric(corr))
    stop("'corr' must be a symmetric matrix")
  D = ncol(corr)

  Ds = length(shape)
  if(Ds > D)
    warning("'shape' longer than width of 'corr', truncating to fit")
  if(Ds != D)
    shape = rep(shape, length.out=D)

  Dr = length(rate)
  if(Dr > D)
    warning("'rate' longer than width of 'corr', truncating to fit")
  if(Dr != D)
    rate = rep(rate, length.out=D)

  if(D == 1) rgamma(n, shape, rate)

  ## generate standard multivariate normal matrix, convert to CDF

  Z = rmvnorm(n, cov=corr)
  #Z = mvtnorm::rmvnorm(n, sigma=corr)

  cdf = pnorm(Z)

  ## convert to gamma, return

  sapply(1:D, function(d) qgamma(cdf[,d], shape[d], rate[d]))
}

## Multivariate normal distribution

# density function
# TO DO:  can we speed this up (and still retain numerical stability) by using the more familiar LU-decomposition-based expressions such as solve() and determinant() for the terms of the log-PDF?	 Something to experiment with, once we're sure that everything else is working ...
dmvnorm <- function(x, mean=rep(0, ncol(x)), cov=diag(ncol(x)), log=FALSE)
{
  ## preprocessing

  # munge data PRN
  if(is.data.frame(x))
    x = as.matrix(x)
  else if(is.vector(x))
    x = matrix(x, nrow=1)

  # gather statistics
  N = nrow(x)
  D = ncol(x)

  # adjust data by mean
  x = t(t(x) - mean)

  ## build terms of log-PDF

  qr.cov = qr(cov)
  determinant.term = sum(log(abs(diag(qr.cov$qr))))

  constant.term = D * LC_LOG2PI

  distance.term = if(N == 1)
    as.vector(x %*% qr.solve(qr.cov, t(x)))
  else
    rowSums((x %*% qr.solve(qr.cov)) * x)

  ## calculate and return (log-)density

  res = -0.5 * (determinant.term + constant.term + distance.term)

  if(log) res else exp(res)
}

# random generation function
rmvnorm <- function(n, mean=NULL, cov=NULL)
{
  ## munge parameters PRN and deal with the simplest univariate case

  if(is.null(mean))
    if(is.null(cov))
      return(rnorm(n))
  else
    mean = rep(0, nrow(cov))
  else if (is.null(cov))
    cov = diag(length(mean))

  ## gather statistics, do sanity checks

  D = length(mean)
  if (D != nrow(cov) || D != ncol(cov))
    stop("length of mean must equal nrow and ncol of cov")

  E = eigen(cov, symmetric=TRUE)
  if (any(E$val < 0))
    stop("Numerically negative definite covariance matrix")

  ## generate values and return

  mean.term = mean
  covariance.term = E$vec %*% (t(E$vec) * sqrt(E$val))
  independent.term = matrix(rnorm(n*D), nrow=D)

  drop(t(mean.term + (covariance.term %*% independent.term)))
}

## Multivariate Pearson Type VII distribution

# probability density function
dmvpvii <- function(x, mean, scale, shape, log=FALSE)
{
  # munge data PRN
  if(is.data.frame(x))
    x = as.matrix(x)
  else if(is.vector(x))
    x = matrix(x, nrow=1)

  # gather statistics and named parameters
  D = ncol(x)
  x.prime = t(t(x) - mean)
  Sigma = scale
  delta = rowSums(x.prime %*% .fast_inv_sym(scale) * x.prime)
  # squared Mahalanobis distance
  lambda.prime = 1 + 0.5*delta
  alpha = shape
  alpha.prime = alpha + 0.5*D

  # calculate log density
  retn = ( -0.5*D*LC_LOG2PI + lgamma(alpha.prime) - lgamma(alpha)
           - 0.5*.fast_labs_det(Sigma) - alpha.prime*log(lambda.prime) )

  # send it back
  if(log) retn else {
    attributes(retn) = NULL # strip "logarithm = TRUE" attribute
    exp(retn)
  }
}

# random generation function
rmvpvii <- function(n, mean, scale, shape)
  t(mean + t(rmvnorm(n, cov=scale) / sqrt(rgamma(n, shape, 1))))

## Multivariate Weibull (WeiSD) distribution using normal copula

# density function (x is a N-by-D matrix, shape and decay are D-length vectors, corr is a D-by-D matrix)
dmvweisd <- function(x, shape, decay, corr=diag(ncol(x)), log=FALSE)
{
  ## munge data PRN

  if(is.data.frame(x))
    x = as.matrix(x)
  else if(is.vector(x))
    x = matrix(x, nrow=1)

  ## gather statistics

  D = ncol(x)

  cdf = sapply(1:D, function(d) pweisd(x[,d], shape[d], decay[d]))
  normcdf = qnorm(cdf)
  normcdf = forceRange(normcdf, LC_MINSTDNORM, LC_MAXSTDNORM)
  # deal with outliers so small (large) that their effective CDF is 0 (1); this represents a robustification approach

  ## calculate copula contribution to log-PDF

  res.copula.add = dmvnorm(normcdf, cov=corr, log=TRUE)
  res.copula.sub = rowSums(dnorm(normcdf, log=TRUE))
  res.copula = res.copula.add - res.copula.sub

  ## calculate marginal contribution to log-PDF

  res.data = rowSums(sapply(1:D, function(d)
    dweisd(x[,d], shape[d], decay[d], log=TRUE)))

  ## final calculations and return

  retn = res.copula + res.data

  if(log) retn else exp(retn)
}

# random generation function
rmvweisd <- function(n, shape=1, decay=1, corr=diag(length(shape)))
{
  ## extract parameters, do sanity checks, deal with univariate case

  if(!is.matrix(corr) || !isSymmetric(corr))
    stop("'corr' must be a symmetric matrix")
  D = ncol(corr)

  Ds = length(shape)
  if(Ds > D)
    warning("'shape' longer than width of 'corr', truncating to fit")
  if(Ds != D)
    shape = rep(shape, length.out=D)

  Dd = length(decay)
  if(Dd > D)
    warning("'decay' longer than width of 'corr', truncating to fit")
  if(Dd != D)
    decay = rep(decay, length.out=D)

  if(D == 1) rweisd(n, shape, decay)

  ## generate standard multivariate normal matrix, convert to CDF

  Z = rmvnorm(n, cov=corr)
  cdf = pnorm(Z)

  ## convert to Weibull (WeiSD), return

  sapply(1:D, function(d) qweisd(cdf[,d], shape[d], decay[d]))
}

### PARAMETER ESTIMATION FUNCTIONS, WITH OPTIONAL WEIGHTS

## exponential distribution

thetahat.exp <- function(x, w=1, aslist=TRUE)
{
  ## sanity check and gather statistics

  mx = min(x)
  if(mx < 0)
    stop("data must be non-negative")
  if(mx < LC_EPS)
    x = forceRange(x, LC_EPS) # prevent instability

  N = length(x)
  Nw = length(w)
  if(Nw > N)
    warning("weights longer than data, truncating to fit")
  if(Nw != N)
    w = rep(w, length.out=N)

  ## MLE

  rate = sum(w) / sum(w * x)

  ## package it up, send it back

  if(aslist)
    list(rate=rate)
  else
    c(rate=rate)
}

thetahat.mvexp <- function(x, w=1, aslist=TRUE)
{
  ## munge data PRN, gather statistics, do sanity checks

  if(!is.matrix(x)) x = as.matrix(x)

  N = nrow(x)
  D = ncol(x)
  Nw = length(w)
  if(Nw > N)
    warning("weights longer than data, truncating to fit")
  if(Nw != N)
    w = rep(w, length.out=N)

  ## marginal parameter estimation

  res = sapply(1:D, function(d) thetahat.exp(x[,d], w, FALSE))
  names(res) = colnames(x)

  retn = list(rate=res)

  ## copula parameter estimation

  cdf = sapply(1:D, function(d) pexp(x[,d], retn$rate[d]))
  normcdf = qnorm(cdf)
  normcdf = forceRange(normcdf, LC_MINSTDNORM, LC_MAXSTDNORM)
  # deal with outliers so small (large) that their effective CDF is 0 (1); this represents a robustification approach

  retn$corr = .rhohat(normcdf, w)

  ## package it up, send it back

  if(aslist)
    return(retn)
  else
    unlist(retn)
}

## gamma distribution

thetahat.gamma <- function(x, w=1, aslist=TRUE, shape.guess=NULL)
{
  ## sanity check and gather statistics

  mx = min(x)
  if(mx < 0)
    stop("data must be non-negative")
  if(mx < LC_EPS)
    x = forceRange(x, LC_EPS) # prevent instability

  N = length(x)
  Nw = length(w)
  if(Nw > N)
    warning("weights longer than data, truncating to fit")
  if(Nw != N)
    w = rep(w, length.out=N)

  sw = sum(w)
  lx = log(x)

  ## take advantage of existing multivariate normal estimators to do initial method-of-moments parameter estimation

  if(is.null(shape.guess)) {
    xmoments = thetahat.norm(x, w)
    shape.guess = xmoments$mean^2 / xmoments$var
  }

  ## profile likelihood parameter estimation

  toptim <- function(shape)
  {
    rate = shape*sw / sum(w*x)

    term1 = shape * log(rate) * sw
    term2 = lgamma(shape) * sw
    term3 = (shape-1) * sum(w * lx)
    term4 = shape * sw # = rate * sum(w*x)

    -(term1 - term2 + term3 - term4)
  }

  shape = nlminb(shape.guess, toptim, lower=LC_EPS)$par
  rate = shape*sw / sum(w*x)
  scale = 1/rate

  ## package it up, send it back

  if(aslist)
    list(shape=shape, rate=rate, scale=scale)
  else
    c(shape=shape, rate=rate, scale=scale)
}

thetahat.mvgamma <- function(x, w=1, aslist=TRUE, shape.guess=NULL)
{
  ## munge data PRN, gather statistics, do sanity checks

  if(!is.matrix(x)) x = as.matrix(x)

  N = nrow(x)
  D = ncol(x)
  Nw = length(w)
  if(Nw > N)
    warning("weights longer than data, truncating to fit")
  if(Nw != N)
    w = rep(w, length.out=N)

  ## marginal parameter estimation
  res = sapply(1:D, function(d)
    thetahat.gamma(x[,d], w, FALSE, shape.guess[[d]]))
  colnames(res) = colnames(x)

  retn = list(shape=res["shape",], rate=res["rate",])

  ## copula parameter estimation

  cdf = sapply(1:D, function(d) pgamma(x[,d], retn$shape[d], retn$rate[d]))
  normcdf = qnorm(cdf)
  normcdf = forceRange(normcdf, LC_MINSTDNORM, LC_MAXSTDNORM)
  # deal with outliers so small (large) that their effective CDF is 0 (1); this represents a robustification approach
  colnames(normcdf) = colnames(x)

  retn$corr = .rhohat(normcdf, w)

  ## package it up, send it back

  if(aslist)
    return(retn)
  else
    unlist(retn)
}

## normal distribution

thetahat.norm <- function(x, w=1, aslist=TRUE)
{
  ## gather statistics, do sanity checks

  N = length(x)
  Nw = length(w)
  if(Nw > N)
    warning("weights longer than data, truncating to fit")
  if(Nw != N)
    w = rep(w, length.out=N)

  ## parameter estimation

  sw = sum(w)

  mean = sum(w * x) / sw
  var = sum(w * x^2) / sw - mean^2
  sd = sqrt(var)

  ## package it up, send it back

  if(aslist)
    list(mean=mean, var=var, sd=sd)
  else
    c(mean=mean, var=var, sd=sd)
}

thetahat.mvnorm <- function(x, w=1, aslist=TRUE)
{
  ## munge data PRN, gather statistics, do sanity checks

  if(!is.matrix(x)) x = as.matrix(x)

  N = nrow(x)
  Nw = length(w)
  if(Nw > N)
    warning("weights longer than data, truncating to fit")
  if(Nw != N)
    w = rep(w, length.out=N)

  ## parameter estimation

  sw = sum(w)
  wX = w * x

  mean = colSums(wX) / sw
  cov = t(x) %*% (wX) / sw - tcrossprod(mean)
  # "tcrossprod(mean)" is equivalent to "mean %*% t(mean)" or "outer(mean, mean)" but should be slightly faster

  ## package it up, send it back

  if(aslist)
    list(mean=mean, cov=cov)
  else
    rbind(mean, cov)
}

## Pearson Type VII (PVII) distribution

# If iter.max > 0, a complete EM estimation will be carried out, otherwise only 1 step.  If not NULL, "theta" should contain current parameter estimates (elements $mean, $shape, and $rate.)
thetahat.pvii <- function(x, w=1, aslist=TRUE, iter.max=LC_ITER_MAX,
                          theta=NULL)
{
  ## sanity check and gather statistics

  N = length(x)
  Nw = length(w)
  if(Nw > N)
    warning("weights longer than data, truncating to fit")
  if(Nw != N)
    w = rep(w, length.out=N)

  sw = sum(w)

  # TO DO:  you should also sanity check theta if not NULL

  ## initialize PRN and do first iteration

  if(is.null(theta)) theta = .pvii_init_params(x, w)

  theta = .pvii_emstep(x, w, sw, theta)
  iter = 1

  ## EM iteration PRN

  if(iter.max) {

    obj.old = (theta$obj - LC_EPS) / 2
    # dummy value to ensure iteration (subtract LC_EPS before dividing in the incredibly unlikely event that theta$obj is very very near 0)

    while(abs(theta$obj/obj.old - 1) > LC_ITER_TOL && iter < iter.max)
    {
      obj.old = theta$obj
      iter = iter + 1
      theta = .pvii_emstep(x, w, sw, theta)
    }
  }

  ## package it up, send it back

  theta$obj = NULL # we don't want this
  retn = if(aslist) theta else unlist(theta)
  attr(retn, "iter") = iter
  return(retn)
}
thetahat.mvpvii <- function(x, w=1, aslist=TRUE, iter.max=LC_ITER_MAX,
                            theta=NULL)
{
  ## sanity check and gather statistics

  N = nrow(x)
  Nw = length(w)
  if(Nw > N)
    warning("weights longer than data, truncating to fit")
  if(Nw != N)
    w = rep(w, length.out=N)

  sw = sum(w)

  # TO DO:  you should also sanity check theta if not NULL

  ## initialize PRN and do first iteration

  if(is.null(theta)) theta = .mvpvii_init_params(x, w)

  tX = t(x)
  theta = .mvpvii_emstep(x, tX, w, sw, theta)
  iter = 1

  ## EM iteration PRN

  if(iter.max) {

    obj.old = (theta$obj - LC_EPS) / 2
    # dummy value to ensure iteration (subtract LC_EPS before dividing in the incredibly unlikely event that theta$obj is very very near 0)

    while(abs(theta$obj/obj.old - 1) > LC_ITER_TOL && iter < iter.max)
    {
      obj.old = theta$obj
      iter = iter + 1
      theta = .mvpvii_emstep(x, tX, w, sw, theta)
    }
  }

  ## package it up, send it back

  theta$obj = NULL # we don't want this
  retn = if(aslist) theta else unlist(theta)
  attr(retn, "iter") = iter
  return(retn)
}

## Weibull (WeiSD) distribution

thetahat.weisd <- function(x, w=1, aslist=TRUE, shape.guess=NULL)
{
  ## sanity check and gather statistics

  mx = min(x)
  if(mx < 0)
    stop("data must be non-negative")
  if(mx < LC_EPS)
    x = forceRange(x, LC_EPS) # prevent instability

  N = length(x)
  Nw = length(w)
  if(Nw > N)
    warning("weights longer than data, truncating to fit")
  if(Nw != N)
    w = rep(w, length.out=N)

  lx = log(x)
  sw = sum(w)
  swlx = sum(w * lx)

  ## initial method-of-moments parameter estimation

  if(is.null(shape.guess))
    shape.guess = pi / sqrt(6 * thetahat.norm(lx, w)$var)
  # see Lawless (1982) pp. 18-19

  ## profile likelihood parameter estimation

  toptim <- function(shape)
    -log(shape)*sw + log(sum(w * x^shape))*sw - (shape-1)*swlx

  shape = nlminb(shape.guess, toptim, lower=LC_EPS)$par
  decay = sw / sum(w * x^shape)

  ## package it up, send it back

  if(aslist)
    list(shape=shape, decay=decay)
  else
    c(shape=shape, decay=decay)
}

thetahat.mvweisd <- function(x, w=1, aslist=TRUE, shape.guess=NULL)
{
  ## munge data PRN, gather statistics, do sanity checks

  if(!is.matrix(x)) x = as.matrix(x)

  N = nrow(x)
  D = ncol(x)
  Nw = length(w)
  if(Nw > N)
    warning("weights longer than data, truncating to fit")
  if(Nw != N)
    w = rep(w, length.out=N)

  ## marginal parameter estimation

  res = sapply(1:D, function(d)
    thetahat.weisd(x[,d], w, FALSE, shape.guess[[d]]))
  colnames(res) = colnames(x)

  retn = list(shape=res["shape",], decay=res["decay",])

  ## copula parameter estimation

  cdf = sapply(1:D, function(d) pweisd(x[,d], retn$shape[d], retn$decay[d]))
  normcdf = qnorm(cdf)
  normcdf = forceRange(normcdf, LC_MINSTDNORM, LC_MAXSTDNORM)
  # deal with outliers so small (large) that their effective CDF is 0 (1); this represents a robustification approach
  colnames(normcdf) = colnames(x)

  retn$corr = .rhohat(normcdf, w)

  ## package it up, send it back

  if(aslist)
    return(retn)
  else
    unlist(retn)
}


