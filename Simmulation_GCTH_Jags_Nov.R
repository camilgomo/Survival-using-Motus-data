## GCTH Simulation CJS time variation 
rm(list=ls())
getwd()

library('dclone')
library('R2jags')
library('rjags')
library('runjags')
library('doParallel')
library('random')
n.cores <- 3


# Define parameter values
n.occasions <- 5                   # Number of capture occasions
marked <- rep(500, n.occasions-1)   # Annual number of newly marked individuals
phi <- c(0.814,0.881,0.84,0.877) #Falts hacerlo con solo Phi promedio yp variable con diferentes valores 
#phi <- rep(0.59, n.occasions−1) 
p <- c(0.25,0.387,0.05,0.351)
#p <- rep(0.4, n.occasions−1)

# Define matrices with survival and recapture probabilities
PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked))
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked))

# Define function to simulate a capture-history (CH) matrix
simul.cjs <- function(PHI, P, marked){
  n.occasions <- dim(PHI)[2] + 1
  CH <- matrix(0, ncol = n.occasions, nrow = sum(marked))
  # Define a vector with the occasion of marking
  mark.occ <- rep(1:length(marked), marked[1:length(marked)])
  # Fill the CH matrix
  for (i in 1:sum(marked)){
    CH[i, mark.occ[i]] <- 1       # Write an 1 at the release occasion
    if (mark.occ[i]==n.occasions) next
    for (t in (mark.occ[i]+1):n.occasions){
      # Bernoulli trial: does individual survive occasion?
      sur <- rbinom(1, 1, PHI[i,t-1])
      if (sur==0) break		# If dead, move to next individual 
      # Bernoulli trial: is individual recaptured? 
      rp <- rbinom(1, 1, P[i,t-1])
      if (rp==1) CH[i,t] <- 1
    } #t
  } #i
  return(CH)
}

# Execute function
CH <- simul.cjs(PHI, P, marked)

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
first <- apply(CH, 1, get.first)
time<-seq(1:5)
ntime<-length(time)
nind<-nrow(CH)

sink("GCTH_SIM.jags")
cat("
model{
for(t in 1:(ntime)){ # uninformative priors for annual survival
  phi0[t] ~ dunif(0,1)
  lphi0[t] <- log(phi0[t]/(1-phi0[t]))
  lp0[t] <- log(p0[t]/(1-p0[t]))
}
for(t in 1:(ntime-1)){ 
  p0[t] ~ dunif(0,1)
}
p0[5]~dunif(0.2,0.4) # constrained uniform prior for last interval to provide better estimate of survival 

mig.surv1<-phi0[2]*phi0[3]*phi0[4]*phi0[5]# calculate survival across all stages
mig.surv2<-phi0[2]*phi0[3]*phi0[4] # survival across all but final stage 

# model likelihood expression:
for(i in 1:nind){
    for(t in 1:first[i]){
    	z[i,t] ~ dbern(1)
    	}
    for(t in (first[i]+1):ntime){
    	logit(p[i,t]) <- lp0[t]
    	logit(phi[i,t]) <- lphi0[t]
    mu1[i,t] <- p[i,t]*z[i,t]
    y[i,t] ~ dbern(mu1[i,t])
    mu2[i,t] <- z[i,t-1]*phi[i,t]
    z[i,t] ~ dbern(mu2[i,t])
    }
    }
#
}

",fill=TRUE)
sink()

# Bundle data 
jags.data<-list(y=CH, nind=nind, ntime=ntime, first=first, time=time) 


# function to estimate when the alive state of an individual is known, not essential but helps
known.state.cjs <- function(ch){
  state <- ch
  for (i in 1:dim(ch)[1]){
    n1 <- min(which(ch[i,]==1))
    n2 <- max(which(ch[i,]==1))
    state[i,n1:n2] <- 1
    state[i,n1] <- NA
  }
  state[state==0] <- NA
  return(state)
}

# Initial values for all variables with priors - phi, p and z state
inits <- function(){list(phi0=runif((ntime),0,1), p0=c((runif((ntime-1),0,1)),(runif(1,0.2,0.4))),z=known.state.cjs(CH))}
# 
# Parameters to be estimate - annual survival (phi0), annual recapture (p0) and the two estimates for survival across migratory stages
parameters<-c("phi0","p0","mig.surv1","mig.surv2")

# MCMC settings 
ni <- 50000 # iterations
nt <- 5 # thin
nb <- 30000 # burnin period
nc <- 3 # chains

# Call JAGS from R 
GCTH_sim<- jags(jags.data, inits, parameters, "GCTH_SIM.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, DIC=TRUE)

# Summarize posteriors
print(GCTH_sim, digits = 3)