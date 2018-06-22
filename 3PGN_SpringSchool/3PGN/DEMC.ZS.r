#Differential Evolution Markov Chain with fewer chains and snooker updater
#(ter Braak & Vrugt 2008)

DE_MC.ZS <- function(Npop = 3, Z, FUN, X= matrix(Z[,1:Npop], ncol = Npop), CR= 1.0, F = 2.38, pSnooker= 0.1, pGamma1 = 0.1, n.generation = 10, n.thin = 1, n.burnin = 0, eps.mult =0.2,eps.add = 0, ...){ 
# Differential Evolution Markov Chain applied to X with logposterior specified by FUN 
# Z is the initial population: a matrix of number of parameters by number of individuals (d x m0) 
# X is the initial active population that will be evolved: a matrix of number of parameters by number of individuals (d x N) 
# value 
# FUN(theta ,...) with theta a k vector of parameter values and ... are other arguments to FUN (e.g. NULL, data, ..) 
# 
# eps.mult >0 gives d-dimensional variation around gamma. It adds scaled uncorrelated noise to the proposal. 
#             Its advantage over eps.add is that its effect scales with the differences of vectors in the population 
#             whereas eps.add does not. if the variance of a dimenstion is close to 0, eps.mult gives smaller changes. 
# eps.add >0  is needed to garantee that all positions in the space can be reached. 
#             For targets without gaps, it can set small or even to 0. 
# 
# Value 
# $Draws d x (Npop*n) array  with n the number of retained simulations  [post process with monitor.DE.MC] 
#         matrix(Draws[j,],nrow= Npop) is an Npop x n array (for each j in 1:d) 
# $ accept.prob 
# $ X.final 
# $ logfitness.X.final 
# 
#  ter  Braak C. J. F., and Vrugt J. A. (2008). Differential Evolution Markov Chain 
#  with snooker updater and fewer chains. Statistics and Computing 
#    http://dx.doi.org/10.1007/s11222-008-9104-9 . 
# 
#    identifier   symbol in paper 
#     Npop          N 
#     Npar          d 
#     gamma_par     gamma 
#     gamma_snooker gamma 
#     M0            M0 
#     mZ            M 
#     n.thin        K 
# 
# see also 
# ter Braak, C. J. F. (2006). A Markov Chain Monte Carlo version of 
#    the genetic algorithm Differential Evolution: easy Bayesian computing 
#    for real parameter spaces. Statistics and Computing, 16, 239-249. 
# 
# 
#  Cajo ter Braak, Biometris, Wageningen UR, cajo.terbraak@wur.nl   23 October 2008 
#  
# 
M0 = mZ = ncol(Z) 
#Npop = ncol(X) 
Npar = nrow(X) 
Npar12  =(Npar - 1)/2  # factor for Metropolis ratio DE Snooker update 
F2 = F/sqrt(2*Npar) 
F1 = 1.0 
accept = rep(NA,n.generation) 
iseq = 1:Npop 
rr = NULL 
r_extra = 0 
logfitness_X = apply (X, 2, FUN, ...) 
for (iter in 1:n.generation) { 
   if (iter%%(n.generation/10) == 0) print(c(iter*Npop,logfitness_X)) 
   accepti = 0 
   for (i in iseq){ 
     # select to random different individuals from Z in rr, a 2-vector 
     if ( runif(1)< pSnooker ) {  
     # DE-Snooker update 
       # if (Npop >1) { z =  X[,sample(iseq[-i],1)]} else { # no real advantage and precludes parallel computing 
        rr = sample(1:mZ, 3, replace = FALSE) 
        z = Z[,rr[3]] 
       # }                                                  # no real advantage and precludes parallel computing 
        x_z = X[,i] - z 
        D2 = max(sum(x_z*x_z),1.0e-300) 
        gamma_snooker = runif(1, min=1.2,max=2.2) 
        #gamma_snooker =1.7 
        projdiff = sum((Z[,rr[1]] -Z[,rr[2]]) *x_z)/D2  # inner_product of difference with x_z / squared norm x_z 
        x_prop = X[,i] + (gamma_snooker * projdiff) * x_z 
        x_z = x_prop - z 
        D2prop = max(sum(x_z*x_z), 1.0e-30) 
        r_extra = Npar12 * (log(D2prop) - log(D2))   # Npar12  =(Npar - 1)/2  # extra term in logr for accept - reject ratio 
    } else { 
    # DE-parallel direction update 
       if ( runif(1)< pGamma1 ) { gamma_par = F1 # to be able to jump between modes 
        } else { 
        gamma_par = F2 * runif(Npar, min=1-eps.mult, max=1+eps.mult)    # multiplicative error to be applied to the difference 
         # gamma_par = F2 
       } 
       rr = sample(1:mZ, 2, replace = FALSE) 
       if (eps.add ==0) {  # avoid generating normal random variates if possible 
         x_prop = X[,i] + gamma_par * (Z[,rr[1]]-Z[,rr[2]]) } else { 
         x_prop = X[,i] + gamma_par * (Z[,rr[1]]-Z[,rr[2]])  +  eps.add*rnorm(Npar,0,1) 
       } 
       r_extra = 0 
      } 
     logfitness_x_prop = FUN(x_prop,  ...) 
     logr =  logfitness_x_prop - logfitness_X[i] 
    # print(c(logfitness_X[i], logfitness_x_prop ,logr,r_extra)) 
     if (!is.na(logr) & (logr + r_extra)> log(runif(1)) ){ 
        accepti = accepti+1 
        X[,i] = x_prop 
        logfitness_X[i] = logfitness_x_prop 
     } 
  } # i loop 
  accept[iter] = accepti 
  if (!(iter%%n.thin) ){   
   Z = cbind(Z,X) 
   mZ = ncol(Z) 
   } 

} # n.generation 
 list(Draws= Z[,-(1:(M0 + Npop* floor(n.burnin/n.thin)))] , accept.prob= accept/Npop, X.final = X, logfitness.X.final = logfitness_X) 
} 
