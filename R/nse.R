#' @name nse
#' @docType package
#' @title Computation of Numerical Standard Errors in R
#' @description \code{nse} is an \R package for computing the numerical standard errors, an estimate 
#' of the standard deviation of a simulation result, if the simulation experiment were to be repeated 
#' many times. The package currently implements more than thirty estimators, including batch means 
#' estimators (Geyer, 1992, Section 3.2), initial sequence estimators Geyer (1992, Equation 3.3), 
#' spectrum at zero estimators (Heidelberger and Welch, 1981; Flegal and Jones, 2010), heteroskedasticity 
#' and autocorrelation consistent (HAC) kernel estimators (Newey and West, 1987; Andrews, 1991; Andrews and 
#' Monahan, 1992; Newey and West, 1994; Hirukawa, 2010), and bootstrap estimators Politis and 
#' Romano (1992, 1994); Politis and White (2004). The full set of estimators are described in 
#' Ardia et al. (2016). 
#' 
#' @section Functions:
#' \itemize{
#' \item \code{\link{nse.geyer}}: Geyer NSE estimator.
#' \item \code{\link{nse.spec0}}: Spectral density at zero NSE estimator.
#' \item \code{\link{nse.nw}}: Newey-West NSE estimator.
#' \item \code{\link{nse.andrews}}: Andrews NSE estimator.
#' \item \code{\link{nse.hiruk}}: Hirukawa NSE estimator.
#' \item \code{\link{nse.boot}}: Bootstrap NSE estimator.
#' }
#' 
#' @section Update:
#' The latest version of the package is available at \url{https://github.com/keblu/nse}.
#' @author David Ardia and Keven Bluteau
#' @references 
#' Andrews, D.W.K. (1991). 
#' Heteroskedasticity and autocorrelation consistent covariance matrix estimation. 
#' \emph{Econometrica} \bold{59}(3), pp.817--858. \doi{10.2307/2938229}.
#' 
#' Andrews, D.W.K, Monahan, J.C. (1992).
#' An improved heteroskedasticity and autocorrelation consistent covariance matrix estimator. 
#' \emph{Econometrica} \bold{60}(4), pp.953--966. \doi{10.2307/2951574}.
#' 
#' Ardia, D., Bluteau, K., Hoogerheide, L. (2016).
#' Comparison of multiple methods for computing numerical standard errors: An extensive Monte Carlo study. 
#' Working paper. \url{https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2741587}.
#' 
#' Flegal, J.M., Hughes, J., Vats D. (2010). 
#' Batch means and spectral variance estimators in Markov chain Monte Carlo.
#' \emph{Annals of Statistics} \bold{38}(2), pp.1034--1070. \doi{10.1214/09-aos735}.
#' 
#' Geyer, C.J. (1992). 
#' Practical Markov chain Monte Carlo. 
#' \emph{Statistical Science} \bold{7}(4), pp.473--483.
#' 
#' Heidelberger, P., Welch, Peter D. (1981).
#' A spectral method for confidence interval generation and run length control in simulations. 
#' \emph{Communications of the ACM} \bold{24}(4), pp.233--245. \doi{10.1145/358598.358630}.
#' 
#' Hirukawa, M. (2010). 
#' A two-stage plug-in bandwidth selection and its implementation for covariance estimation.
#' \emph{Econometric Theory}, \bold{26}(3), pp.710--743. \doi{10.1017/s0266466609990089}.
#' 
#' Newey, W.K., West, K.D. (1987).
#' A simple, positive semi-definite, heteroskedasticity and autocorrelationconsistent covariance matrix. 
#' \emph{Econometrica} \bold{55}(3), pp.703--708. \doi{10.2307/1913610}.
#' 
#' Newey, W.K., West, K.D. (1994) .
#' Automatic lag selection in covariance matrix estimation.
#' \emph{Review of Economic Studies} \bold{61}(4), pp.631--653. \doi{doi: 10.3386/t0144}
#' 
#' Politis, D.N., Romano, and J.P. (1992).
#' A circular block-resampling procedure for stationary data.
#' In \emph{Exploring the limits of bootstrap}, John Wiley & Sons, pp.263--270.
#' 
#' Politis, D.N., Romano, and J.P. (1994). 
#' The stationary bootstrap.
#' \emph{Journal of the American Statistical Association} \bold{89}(428), pp.1303--1313. \doi{10.2307/2290993}.
#' 
#' Politis, D.N., White, H. (2004).
#' Automatic block-length selection for the dependent bootstrap.
#' \emph{Econometric Reviews} \bold{23}(1), pp.53--70. \doi{10.1081/etc-120028836}.
NULL

#' @name nse.geyer
#' @title Geyer NSE estimator
#' @description Calculate the variance of sample mean with the method of Geyer (1992).
#' @note  The type "iseq" is a wrapper around \code{\link[mcmc]{initseq}} from the 
#' MCMC package and gives the positive intial sequence estimator.
#'  The type "bm" is the batch mean estimator.
#'  The type "obm" is the overlapping batch mean estimator.
#'  The type "iseq.bm" is a combinaison of the two.
#' @param x A numeric vector or a matrix(only for type "bm").
#' @param type The type c("iseq","bm","obm","iseq.bm").
#' @param nbatch An optional parameter for the type m and iseq.bm.
#' @param iseq.type constraints on function, ("pos") nonnegative, ("dec") nonnegative 
#' and nonincreasing, and ("con") nonnegative, nonincreasing, and convex. The default value is "pos".
#' @import mcmc
#' @author David Ardia and Keven Bluteau
#' @references 
#' Geyer, C.J. (1992). 
#' Practical Markov chain Monte Carlo. 
#' \emph{Statistical Science} \bold{7}(4), pp.473--483.
#' @return  The variance estimator in the univariate case or the variance-covariance 
#' matrix estimator in the multivariate case.
#' @examples
#' n    = 1000
#' ar   = 0.9
#' mean = 1
#' sd   = 1
# 
#' set.seed(1234)
#' x = as.vector(arima.sim(n = n, list(ar = ar), sd = sd) + mean)
# 
#' nse.geyer(x = x, type = "bm", nbatch = 30)
#' nse.geyer(x = x, type = "obm", nbatch = 30)
#' nse.geyer(x = x, type = "iseq", iseq.type = "pos")
#' nse.geyer(x = x, type = "iseq.bm", iseq.type = "con")
#'  
#'@export
nse.geyer = function(x, type = c("iseq", "bm", "obm", "iseq.bm"), 
                     nbatch = 30, iseq.type = c("pos", "dec", "con")) {
  if (is.vector(x)) {
    x = matrix(x,ncol = 1)
  }
  n = dim(x)[1]
  
  type = type[1]  
  if (type == "iseq") {
    
    f.error.multivariate(x) 
    
    if (iseq.type == "pos") {
      var.iseq = mcmc::initseq(x = x)$var.pos
    } else if(iseq.type == "dec") {
      var.iseq = mcmc::initseq(x = x)$var.dec
    } else if(iseq.type == "con") {
      var.iseq = mcmc::initseq(x = x)$var.con
    } else {
      stop("Invalid iseq.type : must be one of c('pos','dec','con')")
    }
    
    iseq = var.iseq / n # Intial sqequence Geyer (1992)
    out = iseq
    
  } else if (type == "bm") {
    
    ncol  = dim(x)[2]
    x     = as.data.frame(x)
    batch = matrix(unlist(lapply(split(x, ceiling(seq_along(x[,1]) / (n / nbatch))), 
                                 FUN = function(x) colMeans(x))), ncol = ncol,byrow = TRUE)
    out   = stats::var(x = batch) / (nbatch - 1)
    
    if (is.matrix(out) && dim(out) == c(1,1)) {
      out = as.vector(out)
    }
    
  } else if(type == "obm") {
    
    out = as.vector(mcmcse::mcse(x, method = "obm")$se^2)
    
  } else if(type == "iseq.bm"){
    
    f.error.multivariate(x)
    batch = unlist(lapply(split(x, ceiling(seq_along(x) / (n / nbatch))), FUN = mean))
    
    iseq.type = iseq.type[1]
    if(iseq.type == "pos") {
      var.iseq = mcmc::initseq(batch)$var.pos
    } else if(iseq.type == "dec") {
      var.iseq = mcmc::initseq(x = batch)$var.dec
    } else if(iseq.type == "con") {
      var.iseq = mcmc::initseq(x = batch)$var.con
    } else {
      stop("Invalid iseq.type : must be one of c('pos','dec','con')")
    }
    
    iseq.bm = var.iseq / nbatch
    out = iseq.bm
    
  } else {
    stop("Invalid type : must be of type c('iseq','bm','iseq.bm')")
  }
  out = unname(out)
  return(out)
}

#' @name nse.spec0
#' @title Spectral density at zero NSE estimator
#' @description Calculate the variance of the mean with the spectrum at zero estimator.
#' @note  This is a wrapper around \code{\link[coda]{spectrum0.ar}} from the CODA package 
#' and \code{\link[mcmcse]{mcse}} from the mcmcse package.
#' @param x A numeric vector.
#' @param type A character string denoting the method to use in estimating the spectral density function.
#' @param lag.prewhite Prewhite the serie before analysis (integer or NULL, i.e. automatic selection).
#' @return The variance estimator.
#' @references 
#' Flegal, J.M., Hughes, J., Vats D. (2010). 
#' Batch means and spectral variance estimators in Markov chain Monte Carlo.
#' \emph{Annals of Statistics} \bold{38}(2), pp.1034--1070. \doi{10.1214/09-aos735}.
#' @author David Ardia and Keven Bluteau
#' @import coda sapa mcmcse
#' @examples 
#' n    = 1000
#' ar   = 0.9
#' mean = 1
#' sd   = 1
#' set.seed(1234)   
#' x = as.vector(arima.sim(n = n, list(ar = ar), sd = sd) + mean)
#'  
#' nse.spec0(x = x, type = "ar", lag.prewhite = 0)
#' nse.spec0(x = x, type = "ar", lag.prewhite = 1)
#' nse.spec0(x = x, type = "ar", lag.prewhite = NULL)
#' 
#' nse.spec0(x = x, type = "glm", lag.prewhite = 0)
#' nse.spec0(x = x, type = "glm", lag.prewhite = 1)
#' nse.spec0(x = x, type = "glm", lag.prewhite = NULL)
#' 
#' nse.spec0(x = x, type = "wosa", lag.prewhite = 0)
#' nse.spec0(x = x, type = "wosa", lag.prewhite = 1)
#' nse.spec0(x = x, type = "wosa", lag.prewhite = NULL)
#' 
#' nse.spec0(x = x, type = "bartlett", lag.prewhite = 0)
#' nse.spec0(x = x, type = "bartlett", lag.prewhite = 1)
#' nse.spec0(x = x, type = "bartlett", lag.prewhite = NULL)
#' 
#' nse.spec0(x = x, type = "tukey", lag.prewhite = 0)
#' nse.spec0(x = x, type = "tukey", lag.prewhite = 1)
#' nse.spec0(x = x, type = "tukey", lag.prewhite = NULL)
#'  
#'@export
nse.spec0 = function(x, type = c("ar", "glm", "wosa", "bartlett", "tukey"), lag.prewhite = 0) {
  
  scale = 1
  if (is.vector(x)){
    x = matrix(x,ncol = 1)
  }
  f.error.multivariate(x)
  
  n = dim(x)[1]
  tmp   = f.prewhite(x, ar.order = lag.prewhite) 
  x     = tmp$ar.resid
  scale = tmp$scale
  
  type = type[1]
  if (type == "ar") {
    spec0 = coda::spectrum0.ar(x)$spec
  }else if (type == "glm") {
    spec0 = coda::spectrum0(x)$spec
  }else if (type == "wosa") {
    spec0 = sapa::SDF(x, method = "wosa", single.sided = TRUE)[1]
  }else if (type == "tukey") {
    out = as.vector((mcmcse::mcse(x, method = "tukey")$se)^2) * scale
    return(out)
  }else if (type == "bartlett") {
    out = as.vector((mcmcse::mcse(x, method = "bartlett")$se)^2) * scale
    return(out)
  }else {
    stop("Invalid type : must be one of c('ar','bartlett','wosa','tukey')")
  }
  spec0 = spec0 * scale
  out   = spec0 / n
  out   = unname(out)
  return(out)
}

#' @name nse.nw
#' @title Newey-West NSE estimator
#' @description Calculate the variance of the mean with the Newey West (1987, 1994) HAC estimator.
#' @note This is a wrapper around \code{\link[sandwich]{lrvar}} from the sandwich package.
#' @param x A numeric vector or matrix.
#' @param lag.prewhite Prewhite the serie before analysis (integer or NULL, i.e. automatic selection).
#' @return The variance estimator in the univariate case or the variance-covariance matrix 
#' estimator in the multivariate case.
#' @author David Ardia and Keven Bluteau
#' @references 
#' Newey, W.K., West, K.D. (1987).
#' A simple, positive semi-definite, heteroskedasticity and autocorrelationconsistent covariance matrix. 
#' \emph{Econometrica} \bold{55}(3), pp.703--708. \doi{10.2307/1913610}.
#' 
#' Newey, W.K., West, K.D. (1994) .
#' Automatic lag selection in covariance matrix estimation.
#' \emph{Review of Economic Studies} \bold{61}(4), pp.631--653. \doi{doi: 10.3386/t0144}
#' @import sandwich
#' @examples 
#' n    = 1000
#' ar   = 0.9
#' mean = 1
#' sd   = 1
#'
#' set.seed(1234)   
#' x = as.vector(arima.sim(n = n, list(ar = ar), sd = sd) + mean)
#'  
#' nse.nw(x = x, lag.prewhite = 0)
#' nse.nw(x = x, lag.prewhite = 1)
#' nse.nw(x = x, lag.prewhite = NULL)
#'@export
nse.nw <- function(x, lag.prewhite = 0) {
  tmp = f.prewhite(x, ar.order = lag.prewhite) 
  lag.prewhite = tmp$ar.order
  out = sandwich::lrvar(x = x, type = "Newey-West", prewhite = lag.prewhite, adjust = TRUE)
  out = unname(out)
  return(out)
}

#' @name nse.andrews
#' @title Andrews NSE estimator
#' @description Calculate the variance of the mean with the kernel based variance estimator indtroduced by Andrews (1991).
#' @note  This is a wrapper around \code{\link[sandwich]{lrvar}} from the sandwich package and use Andrews (1991) automatic bandwidth estimator.
#' @param x A numeric vector or matrix.
#' @param type The type of kernel used c("bartlett","parzen","qs","trunc","tukey").
#' @param lag.prewhite Prewhite the serie before analysis (integer or NULL, i.e. automatic selection).
#' @param approx Andrews approximation c("AR(1)", "ARMA(1,1)")
#' @return The variance estimator in the univariate case or the variance-covariance matrix estimator in the multivariate case.
#' @references 
#' Andrews, D.W.K. (1991).
#' Heteroskedasticity and autocorrelation consistent covariance matrix estimation. 
#' \emph{Econometrica} \bold{59}(3), pp.817--858. \doi{10.2307/2938229}.
#' 
#' Andrews, D.W.K, Monahan, J.C. (1992).
#' An improved heteroskedasticity and autocorrelation consistent covariance matrix estimator. 
#' \emph{Econometrica} \bold{60}(4), pp.953--966. \doi{10.2307/2951574}.
#' Newey, W.K., West, K.D. (1987).
#' A simple, positive semi-definite, heteroskedasticity and autocorrelationconsistent covariance matrix. 
#' \emph{Econometrica} \bold{55}(3), pp.703--708. \doi{10.2307/1913610}.
#' 
#' Newey, W.K., West, K.D. (1994) .
#' Automatic lag selection in covariance matrix estimation.
#' \emph{Review of Economic Studies} \bold{61}(4), pp.631--653. \doi{doi: 10.3386/t0144}
#' @author David Ardia and Keven Bluteau
#' @import sandwich
#' @examples 
#' n    = 1000
#' ar   = 0.9
#' mean = 1
#' sd   = 1
#'  
#' set.seed(1234)   
#' x = as.vector(arima.sim(n = n, list(ar = ar), sd = sd) + mean)
#'  
#'nse.andrews(x = x, type = "bartlett", lag.prewhite = 0)
#'nse.andrews(x = x, type = "bartlett", lag.prewhite = 1)
#'nse.andrews(x = x, type = "bartlett", lag.prewhite = NULL)
#'
#'nse.andrews(x = x, type = "parzen", lag.prewhite = 0)
#'nse.andrews(x = x, type = "parzen", lag.prewhite = 1)
#'nse.andrews(x = x, type = "parzen", lag.prewhite = NULL)
#'
#'nse.andrews(x = x, type = "tukey", lag.prewhite = 0)
#'nse.andrews(x = x, type = "tukey", lag.prewhite = 1)
#'nse.andrews(x = x, type = "tukey", lag.prewhite = NULL)
#'
#'nse.andrews(x = x, type = "qs", lag.prewhite = 0)
#'nse.andrews(x = x, type = "qs", lag.prewhite = 1)
#'nse.andrews(x = x, type = "qs", lag.prewhite = NULL)
#'
#'nse.andrews(x = x, type = "trunc", lag.prewhite = 0)
#'nse.andrews(x = x, type = "trunc", lag.prewhite = 1)
#'nse.andrews(x = x, type = "trunc", lag.prewhite = NULL)
#'@export
nse.andrews <- function(x, type = c("bartlett", "parzen", "tukey", "qs", "trunc"), lag.prewhite = 0, approx = c("AR(1)", "ARMA(1,1)")) {
  tmp = f.prewhite(x, ar.order = lag.prewhite) 
  lag.prewhite = tmp$ar.order
  type.sandwich = f.type.sandwich(type.in = type)
  out = sandwich::lrvar(x = x, type = "Andrews", prewhite = lag.prewhite, adjust = TRUE, kernel = type.sandwich, approx = approx)
  out = unname(out)
  return(out)
}

#' @name nse.hiruk
#' @title Hirukawa NSE estimator
#' @description Calculate the variance of the mean with the kernel based variance estimator 
#' by Andrews (1991) using Hirukawa (2010) automatic bandwidth estimator.
#' @details This is a wrapper around \code{\link[sandwich]{lrvar}} from the sandwich package 
#' and use Hirukawa (2010) automatic bandwidth estimator.
#' @param x A numeric vector.
#' @param lag.prewhite Prewhite the serie before analysis (integer or NULL, i.e. automatic selection).
#' @param type The type of kernel used c("Bartlett","Parzen").
#' @author David Ardia and Keven Bluteau
#' @references 
#' Hirukawa, M. (2010). 
#' A two-stage plug-in bandwidth selection and its implementation for covariance estimation.
#' \emph{Econometric Theory}, \bold{26}(3), pp.710--743. \doi{10.1017/s0266466609990089}.
#' @import sandwich
#' @return The variance estimator.
#' @examples
#' n    = 1000
#' ar   = 0.9
#' mean = 1
#' sd   = 1
#' 
#' set.seed(1234) 
#' x = as.vector(arima.sim(n = n, list(ar = ar), sd = sd) + mean)
#'  
#' nse.hiruk(x = x, type = "bartlett", lag.prewhite = 0)
#' nse.hiruk(x = x, type = "bartlett", lag.prewhite = 1)
#' nse.hiruk(x = x, type = "bartlett", lag.prewhite = NULL)
#' 
#' nse.hiruk(x = x, type = "parzen", lag.prewhite = 0)
#' nse.hiruk(x = x, type = "parzen", lag.prewhite = 1)
#' nse.hiruk(x = x, type = "parzen", lag.prewhite = NULL)
#' 
#'@export
nse.hiruk <- function(x, type = c("bartlett", "parzen"), lag.prewhite = 0) {
  f.error.multivariate(x)
  tmp = f.prewhite(x, ar.order = lag.prewhite) 
  lag.prewhite = tmp$ar.order
  bandwidth = f.hiruk.bandwidth.solve(x, type = type, lag.prewhite = lag.prewhite)
  type.sandwich = f.type.sandwich(type.in = type)
  out = sandwich::lrvar(x = x, type = "Andrews", prewhite = lag.prewhite, adjust = TRUE, kernel = type.sandwich, bw = bandwidth)
  out = unname(out)
  return(out)
}

#' @name nse.boot
#' @title Bootstrap NSE estimator
#' @description Calculate the variance of the mean with a bootstrap variance estimator.
#' @note  Use the automatic blocksize in \code{\link[np]{b.star}} from th np package which is based 
#' on Politis and White (2004) and Patton and al (2009). 
#' Two bootstrap schemes are available; The stationary bootstrap of Politis and Romano  (1994)
#' and the circular bootstrap of Politis and Romano (1992).
#' @param x  A numeric vector or a matrix.
#' @param nb The number of bootstrap replications.
#' @param type The bootstrap schemes c("stationary","circular").
#' @param b the block length. If NULL automatic block length selection.
#' @param lag.prewhite Prewhite the serie before analysis (integer or NULL, i.e. automatic selection).
#' @return The variance estimator in the univariate case or the variance-covariance matrix 
#' estimator in the multivariate case.
#' @author David Ardia and Keven Bluteau
#' @references 
#' Politis, D.N., Romano, and J.P. (1992).
#' A circular block-resampling procedure for stationary data.
#' In \emph{Exploring the limits of bootstrap}, John Wiley & Sons, pp.263--270.
#' 
#' Politis, D.N., Romano, and J.P. (1994). 
#' The stationary bootstrap.
#' \emph{Journal of the American Statistical Association} \bold{89}(428), pp.1303--1313. \doi{10.2307/2290993}.
#' 
#' Politis, D.N., White, H. (2004).
#' Automatic block-length selection for the dependent bootstrap.
#' \emph{Econometric Reviews} \bold{23}(1), pp.53--70. \doi{10.1081/etc-120028836}.
#' @import np Rcpp stats
#' @useDynLib nse
#' @importFrom Rcpp evalCpp   
#' @examples  
#' n    = 1000
#' ar   = 0.9
#' mean = 1
#' sd   = 1
#' 
#' set.seed(1234) 
#' x = as.vector(arima.sim(n = n, list(ar = ar), sd = sd) + mean)
#'  
#' set.seed(1234)
#' nse.boot(x = x, nb = 1000, type = "stationary", b = NULL, lag.prewhite = 0)
#' nse.boot(x = x, nb = 1000, type = "stationary", b = NULL, lag.prewhite = 1)
#' nse.boot(x = x, nb = 1000, type = "stationary", b = NULL, lag.prewhite = NULL)
#' 
#' nse.boot(x = x, nb = 1000, type = "stationary", b = 10, lag.prewhite = 0)
#' nse.boot(x = x, nb = 1000, type = "stationary", b = 10, lag.prewhite = 1)
#' nse.boot(x = x, nb = 1000, type = "stationary", b = 10, lag.prewhite = NULL)
#' 
#' nse.boot(x = x, nb = 1000, type = "circular", b = NULL, lag.prewhite = 0)
#' nse.boot(x = x, nb = 1000, type = "circular", b = NULL, lag.prewhite = 1)
#' nse.boot(x = x, nb = 1000, type = "circular", b = NULL, lag.prewhite = NULL)
#' 
#' nse.boot(x = x, nb = 1000, type = "circular", b = 10, lag.prewhite = 0)
#' nse.boot(x = x, nb = 1000, type = "circular", b = 10, lag.prewhite = 1)
#' nse.boot(x = x, nb = 1000, type = "circular", b = 10, lag.prewhite = NULL)
#'@export
nse.boot <- function(x, nb, type = c("stationary", "circular"), b = NULL, lag.prewhite = 0){
  
  x = as.vector(x)
  
  # prewhiteneing
  tmp   = f.prewhite(x, ar.order = lag.prewhite) 
  x     = tmp$ar.resid
  scale = tmp$scale
  
  # optimal block length selection
  if (is.null(b)){
    blockSize = np::b.star(data = x, round = TRUE)
    if(type == "stationary"){
      b = as.numeric(blockSize[1,1])
    } else if(type == "circular"){
      b = as.numeric(blockSize[1,2])
    }
  }
  
  out = scale * stats::var(f.bootstrap(x = x, nb = nb, statistic = colMeans, b = b, type = type)$statistic)
  out = as.vector(out)
  out = unname(out)
  return(out)
}
