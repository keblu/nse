---
title: 'nse: Computation of Numerical Standard Errors in R'
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

# Introduction

`nse` is an R package for computing the numerical standard errors, an estimate of the
standard deviation of a simulation result, if the simulation experiment were to be repeated
many times. The package currently implements more than thirty estimators, including 
batch means estimators, initial sequence estimators, spectrum at zero 
estimators, heteroskedasticity and autocorrelation 
consistent (HAC) kernel estimators and bootstrap estimators. 

The full set of methods available is presented in Ardia et al. (2016). 

The latest version of the package is available at 'https://github.com/keblu/nse'.

# Installation
To install the package, run the following commands in R:

R> install.packages("devtools")

R> devtools::install_github("keblu/nse", dependencies = TRUE)

Then check the help of the various files, and run the examples:

R> library("nse")

R> ?nse
