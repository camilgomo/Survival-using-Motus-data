### Survival analyses code for: Gómez et al.Estimating apparent survival along hemispheric migratory routes using Motus tracking. Movement Ecology.

rm(list=ls())
getwd()

library('dclone')
library('R2jags')
library('rjags')
library('runjags')
library('doParallel')
library('random')
n.cores <- 3

surv <- surv<-read.csv("allSpp.motus.jags.csv",header=TRUE)

head(surv)

y<-surv[,c('T1', 'T2', 'T3', 'T4', 'T5', 'T6')]
time<-seq(1:6)
n.occasions <-length(time)
nind<-nrow(y)
group <- surv[, "Species"]
group <- ifelse(surv[,"Species"]=="GCTH",1,ifelse(surv[,"Species"]=="BLPW",2,ifelse(surv[,"Species"]=="SWTH",3,4)))
n.group <- length(unique(group)) # number of species
CH <- as.matrix(y)

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Function to create a matrix with information about known latent state z
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

# Function to create a matrix of initial values for latent state z
cjs.init.z <- function(ch,f){
  for (i in 1:dim(ch)[1]){
    if (sum(ch[i,])==1) next
    n2 <- max(which(ch[i,]==1))
    ch[i,f[i]:n2] <- NA
  }
  for (i in 1:dim(ch)[1]){
    ch[i,1:f[i]] <- NA
  }
  return(ch)
}

# Specify model in BUGS language
sink("cjs-int.jags")
cat("
model {

# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      phi[i,t] <- alpha[group[i],t] # assuming full time dependence in survival
      p[i,t] <- beta[group[i],t]  # assuming full time dependence in recapture
      } #t
   } #i

for (g in 1:n.group){
	for (t in 1:(n.occasions-1)){
		alpha[g,t] ~ dunif(0,1)
	}
	for (t in 1:(n.occasions-2)){
		beta[g,t] ~ dunif(0,1)
	}
  beta[g,5] ~ dunif(0.1,0.5)
 
  }

### Derived parameters

# Migratory period survival

GCTH.migsurv1<-alpha[1,2]*alpha[1,3]*alpha[1,4] # survival during migration
GCTH.migsurv2<-alpha[1,2]*alpha[1,3]*alpha[1,4]*alpha[1,5] # survival during migration plus return
BLPW.migsurv1<-alpha[2,1]*alpha[2,2]*alpha[2,3]*alpha[2,4] 
BLPW.migsurv2<-alpha[2,1]*alpha[2,2]*alpha[2,3]*alpha[2,4]*alpha[2,5] 
SWTH.migsurv1<-alpha[3,1]*alpha[3,2]*alpha[3,3]*alpha[3,4] 
SWTH.migsurv2<-alpha[3,1]*alpha[3,2]*alpha[3,3]*alpha[3,4]*alpha[3,5] 

# Comparison of species differences for each time interval

for (t in 1:(n.occasions-1)){
 	GCTH.BLPW[t]<-alpha[1,t] - alpha[2,t] # ignore first interval for comparisons with GCTH
 	GCTH.SWTH[t]<-alpha[1,t] - alpha[3,t]
	BLPW.SWTH[t]<-alpha[2,t] - alpha[3,t]
	}

# Comparison of entire migration for SWTH vs BLPW 

BLPW.SWTH.mig<-BLPW.migsurv1-SWTH.migsurv1

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- 1
   for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
      } #t
   } #i
}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(y = CH,f = f,nind = dim(CH)[1],n.occasions = dim(CH)[2],z = known.state.cjs(CH),n.group=n.group,g = length(unique(group)),group = group)

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f))} 

# Parameters monitored
parameters <- c("alpha","beta","GCTH.migsurv1","GCTH.migsurv2","BLPW.migsurv1","BLPW.migsurv2",
"SWTH.migsurv1","SWTH.migsurv2","GCTH.BLPW","GCTH.SWTH","BLPW.SWTH","BLPW.SWTH.mig")

# MCMC settings
ni <- 25000
nt <- 5
nb <- 10000
nc <- 3

# Call JAGS from R (BRT 7 min)
cjs.int <- jags(jags.data, inits, parameters, "cjs-int.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(cjs.int, digits = 3)   
