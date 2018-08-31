---
title: 'nse: Computation of numerical standard errors in R'
date: "30 August 2018"
tags:
- nse
- numerical standard error
- simulation
- MCMC
- HAC
- bootstrap
authors:
- affiliation: Institute of Financial Analysis - University of Neuchatel
  name: David Ardia
  orcid: 0000-0003-2823-782X
- affiliation: Institute of Financial Analysis - University of Neuchatel
  name: Keven Bluteau
  orcid: 0000-0003-2990-4807
---

[![status](http://joss.theoj.org/papers/9ad602c309d8cbd10ecdc546fdb25462/status.svg)](http://joss.theoj.org/papers/9ad602c309d8cbd10ecdc546fdb25462)
[![CRAN](http://www.r-pkg.org/badges/version/nse)](https://cran.r-project.org/package=nse) 
[![Downloads](http://cranlogs.r-pkg.org/badges/nse?color=brightgreen)](http://www.r-pkg.org/pkg/nse)
[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/nse?color=brightgreen)](http://www.r-pkg.org/pkg/nse)
[![Pending Pull-Requests](http://githubbadges.herokuapp.com/keblu/nse/pulls.svg?style=flat)](https://github.com/keblu/nse/pulls)
[![Github Issues](http://githubbadges.herokuapp.com/keblu/nse/issues.svg)](https://github.com/keblu/nse/issues)

# Summary
`nse` ([Ardia and Bluteau, 2017](http://dx.doi.org/10.21105/joss.00172)) is an R package for computing the numerical standard error (NSE), an estimate 
of the standard deviation of a simulation result, if the simulation experiment were to be repeated 
many times. The package provides a set of wrappers around several R packages, which give access to 
more than thirty NSE estimators, including batch means estimators, initial sequence 
estimators, spectrum at zero estimators, heteroskedasticity and autocorrelation 
consistent (HAC) kernel estimators and bootstrap estimators. See [Ardia and Bluteau (2017)](http://dx.doi.org/10.21105/joss.00172) for details. The full set of methods available in `nse` is summarized in [Ardia et al. (2018)](https://doi.org/10.1515/jtse-2017-0011) together with several examples of applications in econometrics and finance.

The latest stable version of `nse` is available at [https://cran.r-project.org/package=nse](https://cran.r-project.org/package=nse).

The latest development version of `nse` is available at [https://github.com/keblu/nse)](https://github.com/keblu/nse).

# Installation
To install the latest stable version of `nse`, run the following commands in R:

    R> install.packages("nse")

To install the development version of `nse`, run the following commands in R:

    R> install.packages("devtools")

    R> devtools::install_github("keblu/nse", dependencies = TRUE)

Then check the help of the various files and run the examples:

    R> library("nse")

    R> ?nse
    
Please cite `nse` in publications:

    R> citation("nse")
    
    
# References
Ardia, D., Bluteau, K., Hoogerheide, L.F. (2018).      
Methods for computing numerical standard errors: Review and application to Value-at-Risk estimation.        
_Journal of Time Series Econometrics_ **10**(2) pp 1-9.    
[https://doi.org/10.1515/jtse-2017-0011](https://doi.org/10.1515/jtse-2017-0011)      
[http://dx.doi.org/10.2139/ssrn.2741587](http://dx.doi.org/10.2139/ssrn.2741587)

Ardia, D., Bluteau, K. (2017).      
nse: Computation of numerical standard errors in R.       
_Journal of Open Source Software_ **10**(2).      
[http://dx.doi.org/10.21105/joss.00172](http://dx.doi.org/10.21105/joss.00172)  
