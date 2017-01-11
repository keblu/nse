---
title: 'nse: Computation of Numerical Standard Errors in R'
bibliography: paper.bib
date: "10 January 2017"
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

# Summary

`nse` is an R package (@R) for computing the numerical standard error (NSE), an estimate of the
standard deviation of a simulation result, if the simulation experiment were to be repeated
many times. The package provides a set of wrappers around several R packages, which give access to more than thirty estimators, including 
batch means estimators (@Geyer1992[Section 3.2]), initial sequence estimators (@Geyer1992[Equation 3.3]), spectrum at zero estimators (@HeidelbergerWelch1981,@FlegalJones2010), heteroskedasticity and autocorrelation 
consistent (HAC) kernel estimators (@NeweyWest1987,@Andrews1991,@AndrewsMonahan1992,@NeweyWest1994,@Hirukawa2010), and bootstrap estimators @PolitisRomano1992,@PolitisRomano1994,@PolitisWhite2004). The full set of methods available is 
presented in @ArdiaEtAl2016 together with several examples of applications of NSE in econometrics and finance. The latest version of the package is available at 'https://github.com/keblu/nse'.

# References