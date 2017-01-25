---
title: 'nse: Computation of Numerical Standard Errors in R'
date: "24 January 2017"
tags:
- nse
- numerical standard error
- simulation
- MCMC
- HAC
- bootstrap
authors:
- affiliation: Institute of Financial Analysis - University of Neuchâtel
  name: David Ardia
  orcid: 0000-0003-2823-782X
- affiliation: Institute of Financial Analysis - University of Neuchâtel
  name: Keven Bluteau
  orcid: null
---

[![CRAN](http://www.r-pkg.org/badges/version/nse)](https://cran.r-project.org/package=nse) 
[![Downloads](http://cranlogs.r-pkg.org/badges/nse?color=brightgreen)](http://www.r-pkg.org/pkg/nse)
[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/nse?color=brightgreen)](http://www.r-pkg.org/pkg/nse)

# Introduction

`nse` is an R package for computing the numerical standard error (NSE), an estimate 
of the standard deviation of a simulation result, if the simulation experiment were to be repeated 
many times. The package provides a set of wrappers around several R packages, which give access to 
more than thirty NSE estimators, including batch means estimators, initial sequence 
estimators, spectrum at zero estimators, heteroskedasticity and autocorrelation 
consistent (HAC) kernel estimators and bootstrap estimators. 

The full set of methods available in the package is summarized in [Ardia et al. (2016)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2741587) together with several examples of applications of NSE in econometrics and finance.

The latest stable version of the package is available at 'https://cran.r-project.org/package=nse'.

The latest development version of the package is available at 'https://github.com/keblu/nse'.

# Installation
To install the latest stable version of the package, run the following commands in R:

    R> install.packages("nse")

To install the development version of the package, run the following commands in R:

    R> install.packages("devtools")

    R> devtools::install_github("keblu/nse", dependencies = TRUE)

Then check the help of the various files, and run the examples:

    R> library("nse")

    R> ?nse
