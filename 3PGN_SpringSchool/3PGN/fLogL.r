#'---
#'title: "Some likelihood functions"
#'author: Francesco Minunno
#'email: francesco.minunno@helsinki.fi
#'date: 29 April 2015
#'---
#' 

#' Synopsis: In this R script likelihood functions are defined.
#'

#' Gaussian function
#' ===============================
#'
#' This function gives the likelihood value in logarithmic scale for data that are normally distributed.
#'
#' Arguments are:
#'
#' * sims = numeric vector of simulated data;
#' * data = numeric vector of observed data;
#' * data_s = numeric vector of standard deviations of the observed data.
#'
fLogL_Gaus <- function(sims,data,data_s)
{ 
Ri <- (sims -  data)/data_s
logLi <- -0.5*(Ri^2) - 0.5*log(2*pi) - log(  data_s)
sum(logLi)
}

#' Likelihood from Sivia 2006
#' ===============================
#'
#' This is a heavy-tailed likelihood function, so is more suitable for datasets where outliers might be present. 
#' Likelihood values are in logarithmic scale
#' see Sivia, D.S., 2006. Data Analysis: a Bayesian Tutorial, second ed. O. U. Press, 260 pp. for more details.
#' 
#' Arguments are:
#'
#' * sims = numeric vector of simulated data;
#' * data = numeric vector of observed data;
#' * data_s = numeric vector of standard deviations of the observed data.
#' 
fLogL_Siv <- function(sims,data,data_s)
{ 
Ri         <- (sims - data) / data_s
i0         <- which( abs(Ri)<1.e-08 )

logLi      <- log(1-exp(-0.5*Ri^2)) - log(Ri^2) - 0.5*log(2*pi) - log(data_s)
logLi[i0]  <- -0.5*log(2*pi) - log(2*data_s[i0])
sum(logLi)
}




