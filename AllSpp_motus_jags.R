## All Species together Motus Survival Analysis with group effect
rm(list=ls())
getwd()

library('dclone')
library('R2jags')
library('rjags')
library('runjags')
library('doParallel')
library('random')
n.cores <- 3

surv<-read.csv("C:\\Users\\camil\\OneDrive\\Documents\\Camila\\Publicaciones\\Survival Migration MOTUS\\scripts\\SurvivalMOTUS\\allSpp.motus.jags.csv",header=TRUE)
head()

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

#Fixed group and time effects

# Specify model in BUGS language
sink("cjs-add.jags")
cat("
model {

# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      logit(phi[i,t]) <- beta.phi[group[i]] + gamma.phi[t]
      logit(p[i,t]) <- beta.p[group[i]] + gamma.p[t]
      } #t
   } #i
# for survival parameters
for (t in 1:(n.occasions-1)){
   gamma.phi[t] ~ dnorm(0, 0.01)I(-10, 10)          # Priors for time effects
   phi.g1[t] <- 1 / (1+exp(-gamma.phi[t]))          # Back-transformed survival of group 1
   phi.g2[t] <- 1 / (1+exp(-gamma.phi[t]-beta.phi[2]))  # Back-transformed survival of group 2
   phi.g3[t] <- 1 / (1+exp(-gamma.phi[t]-beta.phi[2]-beta.phi[3]))  # Back-transformed survival of group 3
   }
beta.phi[1] <- 0                            # Corner constraint
beta.phi[2] ~ dnorm(0, 0.01)I(-10, 10)      # Prior for difference in group survival
beta.phi[3] ~ dnorm(0, 0.01)I(-10, 10)      # Prior for difference in group survival


# for recapture parameters
for (t in 1:(n.occasions-1)){
   gamma.p[t] ~ dnorm(0, 0.01)I(-10, 10)          # Priors for time effects
   p.g1[t] <- 1 / (1+exp(-gamma.p[t]))          # Back-transformed recapture of group 1
   p.g2[t] <- 1 / (1+exp(-gamma.p[t]-beta.p[2]))  # Back-transformed recapture of group 2
   p.g3[t] <- 1 / (1+exp(-gamma.p[t]-beta.p[2]-beta.p[3]))  # Back-transformed recapture of group 3
   }
beta.p[1] <- 0                            # Corner constraint
beta.p[2] ~ dnorm(0, 0.01)I(-10, 10)      # Prior for difference in group recapture
beta.p[3] ~ dnorm(0, 0.01)I(-10, 10)  


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
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH), g = length(unique(group)), group = group)

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), beta.phi = c(NA, rnorm(1), rnorm(1)), beta.p = c(NA, rnorm(1), rnorm(1)))}  

# Parameters monitored
parameters <- c("phi.g1", "phi.g2","phi.g3", "p.g1", "p.g2", "p.g3", "beta.phi", "beta.p")

# MCMC settings
ni <- 200000
nt <- 5
nb <- 150000
nc <- 3

# Call JAGS from R (BRT 7 min)
cjs.add <- jags(jags.data, inits, parameters, "cjs-add.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(cjs.add, digits = 3)   

# Figure of species survival
lower.f <- upper.f <- lower.m <- upper.m <- numeric()
for (t in 1:(n.occasions-1)){
  lower.f[t] <- quantile(cjs.add$BUGSoutput$sims.list$phi.g1[,t], 0.025)
  upper.f[t] <- quantile(cjs.add$BUGSoutput$sims.list$phi.g1[,t], 0.975)
  lower.m[t] <- quantile(cjs.add$BUGSoutput$sims.list$phi.g2[,t], 0.025)
  upper.m[t] <- quantile(cjs.add$BUGSoutput$sims.list$phi.g2[,t], 0.975)
  lower.m[t] <- quantile(cjs.add$BUGSoutput$sims.list$phi.g3[,t], 0.025)
  upper.m[t] <- quantile(cjs.add$BUGSoutput$sims.list$phi.g3[,t], 0.975)
}
windows()
plot(x=(1:(n.occasions-1))-0.1, y = cjs.add$BUGSoutput$mean$phi.g1, type = "b", pch = 16, ylim = c(0.2, 1), ylab = "Survival probability", xlab = "Latitude", bty = "n", cex = 1.5, axes = FALSE)
axis(1, at = 1:6, labels = rep(NA,6), tcl = -0.25)
axis(1, at = seq(1,6,1), labels = c("4","11","29","39","49", ">59"))
axis(2, at = seq(0.2, 1, 0.1), labels = c("0.2", NA, "0.4", NA, "0.6", NA, "0.8", NA, "1.0"), las = 1)
segments((1:(n.occasions-1))-0.1, lower.f, (1:(n.occasions-1))-0.1, upper.f)
points(x = (1:(n.occasions-1))+0.1, y = cjs.add$BUGSoutput$mean$phi.g2, type = "b", pch = 1, lty = 2, cex = 1.5)
segments((1:(n.occasions-1))+0.1, lower.m, (1:(n.occasions-1))+0.1, upper.m)
points(x = (1:(n.occasions-1))-0.05, y = cjs.add$BUGSoutput$mean$phi.g3, type = "b", pch = 2, lty = 2, cex = 1.5)
segments((1:(n.occasions-1))-0.05, lower.m, (1:(n.occasions-1))-0.05, upper.m)


##### Interaction species x time (independent variation of survival in time)

#Fixed group and time effects

# Specify model in BUGS language
sink("cjs-add.jags")
cat("
model {

# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      phi[i,t] <- eta.phi[group[i],t]
      p[i,t] <- p.g[group[i],t]
      } #t
   } #i
# for survival parameters
for (u in 1:g){
   for (t in 1:(n.occasions-1)){
      eta.phi[u,t] ~ dunif(0, 1)     # Prior for time and group-spec. survival
      } #t
   } #g
# for recapture parameters
for (u in 1:g){
for (t in 1:(n.occasions-1)){
   p.g[u,t] ~ dunif(0, 1)              # Priors for group-spec. recapture
} #t
   } #g

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
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH), g = length(unique(group)), group = group)

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), beta.phi = c(NA, rnorm(1), rnorm(1)), beta.p = c(NA, rnorm(1), rnorm(1)))}  

# Parameters monitored
parameters <- c("phi.g1", "phi.g2","phi.g3", "p.g1", "p.g2", "p.g3", "eta.phi", "p.g")

# MCMC settings
ni <- 2000
nt <- 5
nb <- 15000
nc <- 3

# Call JAGS from R (BRT 7 min)
cjs.add <- jags(jags.data, inits, parameters, "cjs-add.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(cjs.add, digits = 3)   

# Figure of species survival interaction model
lower.f <- upper.f <- lower.m <- upper.m <- numeric()
for (t in 1:(n.occasions-1)){
  lower.f[t] <- quantile(cjs.add$BUGSoutput$sims.list$phi.g1[,t], 0.025)
  upper.f[t] <- quantile(cjs.add$BUGSoutput$sims.list$phi.g1[,t], 0.975)
  lower.m[t] <- quantile(cjs.add$BUGSoutput$sims.list$phi.g2[,t], 0.025)
  upper.m[t] <- quantile(cjs.add$BUGSoutput$sims.list$phi.g2[,t], 0.975)
  lower.m[t] <- quantile(cjs.add$BUGSoutput$sims.list$phi.g3[,t], 0.025)
  upper.m[t] <- quantile(cjs.add$BUGSoutput$sims.list$phi.g3[,t], 0.975)
}
windows()
plot(x=(1:(n.occasions-1))-0.1, y = cjs.add$BUGSoutput$mean$phi.g1, type = "b", pch = 16, ylim = c(0.2, 1), ylab = "Survival probability", xlab = "Latitude", bty = "n", cex = 1.5, axes = FALSE)
axis(1, at = 1:6, labels = rep(NA,6), tcl = -0.25)
axis(1, at = seq(1,6,1), labels = c("4","11","29","39","49", ">59"))
axis(2, at = seq(0.2, 1, 0.1), labels = c("0.2", NA, "0.4", NA, "0.6", NA, "0.8", NA, "1.0"), las = 1)
segments((1:(n.occasions-1))-0.1, lower.f, (1:(n.occasions-1))-0.1, upper.f)
points(x = (1:(n.occasions-1))+0.1, y = cjs.add$BUGSoutput$mean$phi.g2, type = "b", pch = 1, lty = 2, cex = 1.5)
segments((1:(n.occasions-1))+0.1, lower.m, (1:(n.occasions-1))+0.1, upper.m)
points(x = (1:(n.occasions-1))-0.05, y = cjs.add$BUGSoutput$mean$phi.g3, type = "b", pch = 2, lty = 2, cex = 1.5)
segments((1:(n.occasions-1))-0.05, lower.m, (1:(n.occasions-1))-0.05, upper.m)
