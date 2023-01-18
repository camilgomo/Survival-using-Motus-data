## GCTH Motus Survival Analysis 
rm(list=ls())
getwd()

library('dclone')
library('R2jags')
library('rjags')
library('runjags')
library('doParallel')
library('random')
n.cores <- 3

gcthsurv<-read.csv("C:\\Users\\cg478\\Documents\\Camila\\Publicaciones\\Survival Migration MOTUS\\scripts\\SurvivalMOTUS\\gcth.motus.jags.csv",header=TRUE)
gcthsurv<-as.matrix(gcthsurv[,6:10])
y<-gcthsurv
time<-seq(1:5)
ntime<-length(time)
nind<-nrow(y)

# for each individual, this code below identifies when the time interval when each individual was marked and when last seen, this is used for a function below
last<-first<-rep(NA,nind)
Zst<-y
for (i in 1:nind){
  h<-as.vector(y[i,])# treat each individual i as a separate vector
  first[i] <- min((1:5)[!is.na(h) & h==1]) # for each individual find occasion of first marking
  last[i]  <- max((1:5)[!is.na(h) & h==1]) # for each individual find occasion of last observation
  Zst[i,first[i]:last[i]]<-1 # true state is alive between first and last obs
  Zst[i,1:first[i]]<-0 # 0s before first capture
}

sink("GCTH.jags")
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
jags.data<-list(y=y, nind=nind, ntime=ntime, first=first, time=time)

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
inits <- function(){list(phi0=runif((ntime),0,1), p0=c((runif((ntime-1),0,1)),(runif(1,0.2,0.4))),z=known.state.cjs(y))}
# 
# Parameters to be estimate - annual survival (phi0), annual recapture (p0) and the two estimates for survival across migratory stages
parameters<-c("phi0","p0","mig.surv1","mig.surv2")

# MCMC settings 
ni <- 50000 # iterations
nt <- 5 # thin
nb <- 30000 # burnin period
nc <- 3 # chains

# Call JAGS from R 
GCTH<- jags(jags.data, inits, parameters, "GCTH.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, DIC=TRUE)

# Summarize posteriors

print(GCTH, digits = 3)

plot(density(GCTH$BUGSoutput$sims.list$phi0[,1]), xlim = c(0, 1), ylim = c(0, 5), main = "", xlab = expression(phi[1]), ylab = "Density", frame = FALSE, lwd = 2)
abline(h = 1, lty = 2, lwd = 2)

plot(density(GCTH$BUGSoutput$sims.list$phi0[,2]), xlim = c(0, 1), ylim = c(0, 5), main = "", xlab = expression(phi[2]), ylab = "Density", frame = FALSE, lwd = 2)
abline(h = 1, lty = 2, lwd = 2)

plot(density(GCTH$BUGSoutput$sims.list$phi0[,3]), xlim = c(0, 1), ylim = c(0, 5), main = "", xlab = expression(phi[3]), ylab = "Density", frame = FALSE, lwd = 2)
abline(h = 1, lty = 2, lwd = 2)

plot(density(GCTH$BUGSoutput$sims.list$phi0[,4]), xlim = c(0, 1), ylim = c(0, 5), main = "", xlab = expression(phi[4]), ylab = "Density", frame = FALSE, lwd = 2)
abline(h = 1, lty = 2, lwd = 2)

plot(density(GCTH$BUGSoutput$sims.list$phi0[,5]), xlim = c(0, 1), ylim = c(0, 5), main = "", xlab = expression(phi[5]), ylab = "Density", frame = FALSE, lwd = 2)
abline(h = 1, lty = 2, lwd = 2)


plot(density(GCTH$BUGSoutput$sims.list$p0[,1]), xlim = c(0, 1), ylim = c(0, 5), main = "", xlab = expression(p[1]), ylab = "Density", frame = FALSE, lwd = 2)
abline(h = 1, lty = 2, lwd = 2)

plot(density(GCTH$BUGSoutput$sims.list$p0[,2]), xlim = c(0, 1), ylim = c(0, 5), main = "", xlab = expression(p[2]), ylab = "Density", frame = FALSE, lwd = 2)
abline(h = 1, lty = 2, lwd = 2)

plot(density(GCTH$BUGSoutput$sims.list$p0[,3]), xlim = c(0, 1), ylim = c(0, 5), main = "", xlab = expression(p[3]), ylab = "Density", frame = FALSE, lwd = 2)
abline(h = 1, lty = 2, lwd = 2)

plot(density(GCTH$BUGSoutput$sims.list$p0[,4]), xlim = c(0, 1), ylim = c(0, 5), main = "", xlab = expression(p[4]), ylab = "Density", frame = FALSE, lwd = 2)
abline(h = 1, lty = 2, lwd = 2)

plot(density(GCTH$BUGSoutput$sims.list$p0[,5]), xlim = c(0, 1), ylim = c(0, 5), main = "", xlab = expression(p[5]), ylab = "Density", frame = FALSE, lwd = 2)
abline(h = 1, lty = 2, lwd = 2)

############################
sink("GCTH2.jags")
cat("
model{
for(t in 1:(ntime)){ # uninformative priors for annual survival
  phi0[t] ~ dunif(0,1)
  lphi0[t] <- log(phi0[t]/(1-phi0[t]))
  lp0[t] <- log(p0[t]/(1-p0[t]))
}
  p0[1] ~ dnorm(0.21,0.788)    #priors for p based on detectability of survivors
   p0[2] ~ dnorm(0.21,0.37)
    p0[3] ~ dnorm(0.151,0.255)
     p0[4] ~ dnorm(0.277,0.399)
      p0[5] ~ dnorm(0.035,0.101)

mig.surv1<-phi0[2]*phi0[3]*phi0[4]*phi0[5] # calculate survival across all stages
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
jags.data<-list(y=y, nind=nind, ntime=ntime, first=first, time=time)

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
inits <- function(){list(phi0=runif((ntime),0,1), p0=c((runif((ntime-1),0,1)),(runif(1,0,1))),z=known.state.cjs(y))}
# 
# Parameters to be estimate - annual survival (phi0), annual recapture (p0) and the two estimates for survival across migratory stages
parameters<-c("phi0","p0","mig.surv1","mig.surv2")

# MCMC settings 
ni <- 50000 # iterations
nt <- 5 # thin
nb <- 30000 # burnin period
nc <- 3 # chains

# Call JAGS from R 
GCTH2<- jags(jags.data, inits, parameters, "GCTH2.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, DIC=TRUE)

# Summarize posteriors
print(GCTH2, digits = 3)


