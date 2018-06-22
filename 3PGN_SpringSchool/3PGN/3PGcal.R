setwd("C://Users//andrey//Documents//3PGN//3PGN")

library(lhs)
library(coda)

par= read.delim("~/3PGN/3PGN/input/parameters.txt", header=T)

dyn.load("3PG.dll")
#'
.Fortran('init_model')
#' 
#' Define the function to run 3PGN in R
model_3PGN <- function(pValues){
  .Fortran('changepars', pValues)  
  y <- array(1,dim=c(nMonths,8,nSites))
  output <- .Fortran('model',0,y)[[2]]
}

#' Setting model for runs
nMonths <-156   # number of months for wich running the model
nSites <- 5 	# number of sites
varnames <- c('standAge','NEP','dbh','H','WF','WR','WS','StandVol') # output variables considered
load('pSet.rdata') #load parameters for Eucalypus Globulus in Portugal
load('data.rdata')
parNam = as.character(par[,1])
inddata <- which(!is.na(data[,3,1]))
std=c(0,0.5,1,1,0.3,3,3,2)

data[,2,1]=data[,2,1]+rnorm(156,0,std[2]*abs(data[,2,1])*0.1)

for(i in 3:8){
data[inddata,i,]=abs(data[inddata,i,]+rnorm(length(inddata),0,std[i]*0.1*data[inddata,i,]))
}

pSet=par[,2]
#' Likelihood function and Prior
source('fLogL.r') # Read likelihood script
fLogL <- fLogL_Gaus # Choose likelihood: fLogL_Gaus, flogL_Siv
fLogPrior <- function(pSet,parmin,parmax){sum(dunif(pSet, min=parmin, max=parmax, log=T))} # Define Prior distribution

#' Define the function for calibration. The only argument to be passed is a parameter vector
#'
fLogPost <- function(pValues){
  pSet <- pValues
  logprior  <- fLogPrior(pSet,par$min,par$max)
  if (logprior == -Inf) {logpost=-Inf}else{
    out_3PGN <- model_3PGN(pSet)
    
    loglikelihood_NEP= fLogL(out_3PGN[,2,1],data[,2,1],abs(0.1*data[,2,1]))
    loglikelihood_DBH= fLogL(out_3PGN[inddata,3,],data[inddata,3,],abs(0.1*data[inddata,3,]))
    loglikelihood_H= fLogL(out_3PGN[inddata,4,],data[inddata,4,],abs(0.1*data[inddata,4,]))
    loglikelihood_WF= fLogL(out_3PGN[inddata,5,],data[inddata,5,],abs(0.1*data[inddata,5,]))  
    loglikelihood_WR= fLogL(out_3PGN[inddata,6,],data[inddata,6,],abs(0.1*data[inddata,6,])) 
    loglikelihood_WS= fLogL(out_3PGN[inddata,7,],data[inddata,7,],abs(0.1*data[inddata,7,])) 
    loglikelihood_StandVol= fLogL(out_3PGN[inddata,8,],data[inddata,8,],abs(0.1*data[inddata,8,]))
    loglikelihood= loglikelihood_NEP+loglikelihood_DBH+ 
    loglikelihood_H+loglikelihood_WF+loglikelihood_WR+loglikelihood_WS+
    loglikelihood_StandVol
    logpost <- loglikelihood+logprior}
  return((logpost)) 
}

npar=51

source('DEMC.ZS.r')

# parameters of MCMC 
niter <- 3e5 		#total number of iterations
nBI <- 0.5*niter		#burn-in length 
eps <- 0.001 
Npop <- 3
m0 <- 10*npar
CR <- 1.	# crossover probability
F <- 2.38	# differential weighting factor 


Z <- randomLHS(n=npar, k=m0) #use Latin Hypercube to sample the initial population
for (i in 1:npar) Z[,i] <- qunif(Z[,i], min=par$min, max=par$max,log=FALSE) 
#'
BC1_DEMCzs <- DE_MC.ZS(Npop = Npop, Z=Z, FUN=fLogPost, X= matrix(Z[,1:Npop], ncol = Npop), CR= 1.0, F = 2.38, pSnooker= 0.1, pGamma1 = 0.1, n.generation = floor(niter/Npop), n.thin = max(floor(niter/100000),1), n.burnin = floor(nBI/Npop), eps.mult =0.2,eps.add = 0) 

Z <- randomLHS(n=npar, k=m0) #use Latin Hypercube to sample the initial population
for (i in 1:npar) Z[,i] <- qunif(Z[,i], min=par$min, max=par$max,log=FALSE) 
#'
BC2_DEMCzs <- DE_MC.ZS(Npop = Npop, Z=Z, FUN=fLogPost, X= matrix(Z[,1:Npop], ncol = Npop), CR= 1.0, F = 2.38, pSnooker= 0.1, pGamma1 = 0.1, n.generation = floor(niter/Npop), n.thin = max(floor(niter/100000),1), n.burnin = floor(nBI/Npop), eps.mult =0.2,eps.add = 0) 
chain1 <- as.mcmc(t(BC1_DEMCzs$Draws))
chain2 <- as.mcmc(t(BC2_DEMCzs$Draws))
convergence_test <- gelman.diag(chains, confidence = 0.95,autoburnin=FALSE,multivariate=TRUE)

#'
#' Analysis of the results
#' ===============================
#'
#' Check convergence using the method proposed by Gelman and Rubin 1992 
rownames(BC1_DEMCzs$Draws) <- par$nam[parind]
rownames(BC2_DEMCzs$Draws) <- par$nam[parind]
chain1 <- as.mcmc(t(BC1_DEMCzs$Draws))
chain2 <- as.mcmc(t(BC2_DEMCzs$Draws))
chains <- mcmc.list(chain1,chain2)
convergence_test <- gelman.diag(chains, confidence = 0.95,autoburnin=FALSE,multivariate=TRUE)
convergence_test
print(convergence_test)
gelman.plot(chains,autoburnin=FALSE)

 

#'
#' Plot traceplot and marginal posterior distribution
plot(chain1)



nSample <- 1000
ind <- floor(seq(1,(dim(BC1_DEMCzs$Draws)[2]), length.out=nSample))
par(mfrow=c(2,1),oma=c(0,0,2,0))

plot(data[,2,1],col=2,xlab='',ylab='tDW/ha', main='NEP 1')
p <-pSet
for (i in 1:nSample){
  pSet <- BC1_DEMCzs$Draws[,ind[i]]
  out_3PGN <- model_3PGN(pSet)
  if (i==1) out_3PGN_mean <- out_3PGN[,2,1]
  out_3PGN_mean <- (out_3PGN_mean*(i-1)+out_3PGN[,2,1])/i
  lines(out_3PGN[,2,1], col='light green')
}
lines(out_3PGN_mean, col='dark green')
legend('topleft', c('obs', 'sim', 'mod unc'),lty=c(0,1,1),col=c(2,'dark green','light green'),pch=c(1,NA,NA))

varunits=c('Years','tDW/ha','cm','m','tDW/ha','tDW/ha','tDW/ha','m³/ha')
for (j in 1:5){
  for (k in 3:8){
    paste(varnames[k],j,sep=" ")
    nSample <- 1000
    ind <- floor(seq(1,(dim(BC1_DEMCzs$Draws)[2]), length.out=nSample))
    par(mfrow=c(2,1),oma=c(0,0,2,0))
    
    plot(data[,k,j],col=2,xlab='',ylab=varunits[k], main=paste(varnames[k],j,sep=" "))
    p <-pSet
    for (i in 1:nSample){
      pSet <- BC1_DEMCzs$Draws[,ind[i]]
      out_3PGN <- model_3PGN(pSet)
      if (i==1) out_3PGN_mean <- out_3PGN[,k,j]
      out_3PGN_mean <- (out_3PGN_mean*(i-1)+out_3PGN[,k,j])/i
      lines(out_3PGN[,k,j], col='light green')
     
    }
    lines(out_3PGN_mean, col='dark green')
    legend('topleft', c('obs', 'sim', 'mod unc'),lty=c(0,1,1),col=c(2,'dark green','light green'),pch=c(1,NA,NA))
  
  
  }
}



nSample <- 1000
ind <- floor(seq(1,(dim(BC1_DEMCzs$Draws)[2]), length.out=nSample))

set.seed(2)
NEPerror <- NEPsim <- matrix(0, nrow = nSample, ncol = 156)
for (i in 1:nSample) {
  pSet <- BC1_DEMCzs$Draws[,ind[i]]
  out_3PGN <- model_3PGN(pSet)
  NEPsim[i, ] <-  out_3PGN[,2,1]
  NEPerror[i, ] <- rnorm(156, mean = 0, sd = abs(0.1*data[,2,1]))
}

#' ## Plot of GPP over time: model prediction incl. error
plot(data[,2,1],col=2,xlab='',ylab='tDW/ha', main='NEP 1')
for(i in 1:nrow(NEPsim)) {
  lines(NEPsim[i, ] + NEPerror[i, ], col = rgb(0.8, 0.8, 0, 0.05))
}
for(i in 1:nrow(NEPsim)) {
  lines(NEPsim[i, ], col = rgb(0, 0.8, 0, 0.05))
}
lines(colMeans(NEPsim), col = "blue")
points(out_3PGN[,2,1], col=2, cex = 0.4, pch = 19)
legend('topleft', 
       c('obs', 'sim', 'mod unc', 'pred unc'),
       lty=c(0,1,1,1),
       col=c(2,'blue',rgb(0, 0.8, 0, 1), rgb(0.8, 0.8, 0, 1)),
       pch=c(1,NA,NA,NA), bty = "n")



##Other predictions
varunits=c('Years','tDW/ha','cm','m','tDW/ha','tDW/ha','tDW/ha','m³/ha')
for (j in 1:5){
  for (k in 3:8){
    paste(varnames[k],j,sep=" ")
    
    nSample <- 1000
    ind <- floor(seq(1,(dim(BC1_DEMCzs$Draws)[2]), length.out=nSample))
    
    set.seed(2)
    NEPerror <- NEPsim <- matrix(0, nrow = nSample, ncol = 156)
    for (i in 1:nSample) {
      pSet <- BC1_DEMCzs$Draws[,ind[i]]
      out_3PGN <- model_3PGN(pSet)
      NEPsim[i, ] <-  out_3PGN[,k,j]
      NEPerror[i, ] <- rnorm(ncol(NEPerror), mean = 0, sd = abs(0.1*data[inddata,k,j]))
    }
    
    #' ## Plot of GPP over time: model prediction incl. error
    plot(data[,k,j],col=2,xlab='',ylab=varunits[k], main=paste(varnames[k],j,sep=" "))
    for(i in 1:nrow(NEPsim)) {
      lines(NEPsim[i, ] + NEPerror[i, ], col = rgb(0.8, 0.8, 0, 0.05))
    }
    for(i in 1:nrow(NEPsim)) {
      lines(NEPsim[i, ], col = rgb(0, 0.8, 0, 0.05))
    }
    lines(colMeans(NEPsim), col = "blue")
    points(out_3PGN[,k,j], col=2, cex = 0.4, pch = 19)
    legend('topleft', 
           c('obs', 'sim', 'mod unc', 'pred unc'),
           lty=c(0,1,1,1),
           col=c(2,'blue',rgb(0, 0.8, 0, 1), rgb(0.8, 0.8, 0, 1)),
           pch=c(1,NA,NA,NA), bty = "n")
    
    med_resid <- apply(data[,k,j] - t(NEPsim[,k,j]), 1, median)
    plot(med_pred, med_resid)
  }
}


