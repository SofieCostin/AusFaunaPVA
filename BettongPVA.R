# Sofie Costin, Corey Bradshaw, Frédérik Saltré
# Global Ecology, Flinders University globalecologyflinders.com
# Brush-tailed bettong PVA
# requires library - Plotly
# updated 1/3/22

## remove everything
rm(list = ls())

# libraries
library(plotly)
options(scipen = 1000)

## functions
# beta distribution shape parameter estimator function
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

## source/matrix operators
source("matrixOperators.r")


# Leslie matrix ----------------------------------------------------

age.max = 6 # maximum age of females in wild from Threatened Species Listing Advice (FIND SOURCE)

# Stages max values from Thompson et al 2015

# stage 1 (pouch young) = 110/365 days = 0.30 of 1 year
# Stage 2 (independent young) =  70/365 days  = 0.19 of 1 year
# stage 3 (sexually active) = 185/365 days = 0.51 of 1 year

# survival and mortality from Pacioni et al 2017 (Supplementary 1)

# Create vectors: fertility -----------------------------------------------

pr.breed <- 0.893
pr.breed.yr1 <- 0.893*0.51
no.offspring <- 3
no.Foffspring <- no.offspring/2
m.vec <- c(pr.breed.yr1*no.Foffspring, rep(pr.breed*no.Foffspring, age.max)) ## female offspring produced each year
plot(0:6,m.vec,pch=19,type="b")

# fertility errors
m.sd.vec <- c(rep(0.26, age.max+1)) #mean and standard deviations vector, juvenile and adult fertility 

# Create vectors: Survival ------------------------------------------------

#survival
s.vec <- c(((0.956*0.3)+(0.85*0.19)+(0.961*0.51)), rep(0.961, age.max))  ## female  survival # probability of surviving from one year to the next. e.g surviving fourth year of life

# survival errors
s.sd.vec <- c(rep(0.15, age.max+1)) #mean and standard deviations vector, juvenile and adult survival

# Create matrix -----------------------------------------------------------

popmat <- matrix(data = 0, nrow=age.max+1, ncol=age.max+1) # creates a matrix where data = 0 and dimensions age.max x age.max
diag(popmat[2:(age.max+1),]) <- s.vec[1:(age.max)] #  populates popmat diagonally with s.vec
# popmat[age.max+1,age.max+1] <- s.vec[1:(age.max)] # would fill 7,7 with s.vec, but not needed because pr(breed) at age.max = 0
popmat[1,] <- m.vec[age.max+1] ## row 1 of popmat populated with m.vec
popmat.orig <- popmat ## save original matrix
popmat

## matrix properties
max.lambda(popmat) ## 1-yr lambda
max.r(popmat) # rate of population change, 1-yr
ssd <- stable.stage.dist(popmat) ## stable stage distribution
R.val(popmat, age.max) # reproductive value
gen.l <- G.val(popmat, age.max) # mean generation length


# Initial population vector -----------------------------------------------


pop.found <- 500 # change this to change mvp 
init.vec <- ssd * pop.found #initial population vector
plot(0:6,ssd,pch=19,type="b")


# Project -----------------------------------------------------------------


## set time limit for projection in 1-yr increments

yr.now <- 2020 # 
yr.end <- yr.now + round(40*gen.l,0) # set projection end date to 40 generations (IUCN)
t <- (yr.end - yr.now) #time frame
yr.vec <- seq(yr.now,yr.end) # year vector

tot.F <- sum(popmat.orig[1,])
popmat <- popmat.orig # reset matrix 

## set population storage matrices
n.mat <- matrix(0, nrow=age.max+1,ncol=(t+1))  #empty matrix
n.mat[,1] <- init.vec # fill first matrix column with initial population vector

## set up projection loop
for (i in 1:t) {
  n.mat[,i+1] <- popmat %*% n.mat[,i]
}

n.pred <- colSums(n.mat) # number of animals through time period, no density reduction treatment, no carry capacity
yrs <- seq(yr.now, yr.end, 1)
plot(yrs, n.pred,type="b",lty=2,pch=19,xlab="year",ylab="N")

## compensatory density feedback
## population rate of increase relative to carrying capacity. Larger distance between population and K = faster population growth
K.max <- pop.found*2
K.vec <- c(1,K.max*0.5,0.75*K.max,0.85*K.max,0.95*K.max)
plot(K.vec)

red.vec <- c(1,0.98,0.85,0.75,0.6)
plot(K.vec,red.vec,pch=19,type="b")
Kred.dat <- data.frame(K.vec,red.vec)

## logistic power function a/(1+(x/b)^c) 
## fits logistic power function to population relative to carrying capacity, K
param.init <- c(1, K.max, 2) #(1, 15000, 2.5)
fit.lp <- nls(red.vec ~ a/(1+(K.vec/b)^c), 
              data = Kred.dat,
              algorithm = "port",
              start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
              trace = TRUE,      
              nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
fit.lp.summ <- summary(fit.lp)
plot(K.vec,red.vec,pch=19,xlab="N",ylab="reduction factor")
K.vec.cont <- seq(1,2*pop.found,1)
pred.lp.fx <- coef(fit.lp)[1]/(1+(K.vec.cont/coef(fit.lp)[2])^coef(fit.lp)[3])
lines(K.vec.cont,pred.lp.fx,lty=2,lwd=3,col="red")

a.lp <- coef(fit.lp)[1]
b.lp <- coef(fit.lp)[2]
c.lp <- coef(fit.lp)[3]

## compensatory density-feedback deterministic model
## set population storage matrices
n.mat <- matrix(0, nrow=age.max+1,ncol=(t+1))
n.mat[,1] <- init.vec
popmat <- popmat.orig

## set up projection loop
for (i in 1:t) {
  totN.i <- sum(n.mat[,i])
  pred.red <- as.numeric(a.lp/(1+(totN.i/b.lp)^c.lp))
  diag(popmat[2:(age.max+1),]) <- s.vec[1:(age.max)]*pred.red
  # popmat[age.max+1,age.max+1] <- s.vec[age.max+1]
  popmat[1,] <- m.vec
  n.mat[,i+1] <- popmat %*% n.mat[,i]
}

n.pred <- colSums(n.mat)
plot(yrs, n.pred,type="l",lty=2,pch=19,xlab="year",ylab="N", ylim=c(0,4 * pop.found)) # untreated population increases, rate of increase relative to K, no stochastic sampling
abline(h=K.max, lty=2, col="red") # carrying capacity


# ## stochastic projection with density feedback
# 
# ## set storage matrices & vectors
# iter <- 1000 # number of iterations
# itdiv <- iter/10
# 
# n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
# s.arr <- m.arr <- array(data=NA, dim=c(t+1, age.max+1, iter))
# 
# for (e in 1:iter) {
#   popmat <- popmat.orig
#   
#   n.mat <- matrix(0, nrow=age.max+1,ncol=(t+1))
#   n.mat[,1] <- init.vec
#   
#   for (i in 1:t) {
#     # stochastic survival values
#     s.alpha <- estBetaParams(s.vec.upd, s.sd.vec^2)$alpha
#     s.beta <- estBetaParams(s.vec.upd, s.sd.vec^2)$beta
#     s.stoch <- rbeta(length(s.alpha), s.alpha, s.beta)
#     
#     # stochastic fertility sampler (Gaussian)
#     fert.stch <- rnorm(length(popmat[,1]), m.vec, m.sd.vec)
#     m.arr[i,,e] <- ifelse(fert.stch < 0, 0, fert.stch)
#     
#     totN.i <- sum(n.mat[,i], na.rm=T)
#     pred.red <- as.numeric(a.lp/(1+(totN.i/b.lp)^c.lp))
#     
#     diag(popmat[2:(age.max+1),]) <- (s.stoch[-(age.max+1)])*pred.red
#     popmat[age.max+1,age.max+1] <- (s.stoch[age.max+1])*pred.red
#     popmat[1,] <- m.arr[i,,e]
#     n.mat[,i+1] <- popmat %*% n.mat[,i]
#     
#     s.arr[i,,e] <- s.stoch * pred.red
#     
#   } # end i loop
#   
#   n.sums.mat[e,] <- ((as.vector(colSums(n.mat))))
#   
#   if (e %% itdiv==0) print(e) 
#   
# } # end e loop
# 
# n.md <- apply(n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
# n.up <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
# n.lo <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
# 
# # plot probability of N, survival and fertility 
# 
# par(mfrow=c(1,3))
# plot(yrs,n.md,type="l", main = "", xlab="year", ylab="pN1", lwd=2, ylim=c(0.95*min(n.lo),1.05*max(n.up)))
# lines(yrs,n.lo,lty=2,col="red",lwd=1.5)
# lines(yrs,n.up,lty=2,col="red",lwd=1.5)
# 
# s.add <- m.add  <- rep(0, age.max+1)
# for (m in 1:iter) {
#   s.add <- rbind(s.add, s.arr[ceiling(gen.l):(t+1),,m])
#   m.add <- rbind(m.add, m.arr[ceiling(gen.l):(t+1),,m])
# }
# s.add <- s.add[-1,]
# m.add <- m.add[-1,]
# 
# s.md <- apply(s.add, MARGIN=2, median, na.rm=T) # mean s over all iterations
# s.up <- apply(s.add, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
# s.lo <- apply(s.add, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
# 
# plot(age.vec,s.md,type="l", main = "", xlab="age", ylab="s", lwd=2, ylim=c(0.95*min(s.lo),1.05*max(s.up)))
# lines(age.vec,s.lo,lty=2,col="red",lwd=1.5)
# lines(age.vec,s.up,lty=2,col="red",lwd=1.5)
# 
# m.md <- apply(m.add, MARGIN=2, median, na.rm=T) # mean s over all iterations
# m.up <- apply(m.add, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
# m.lo <- apply(m.add, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
# 
# plot(age.vec,m.md,type="l", main = "", xlab="age", ylab="m", lwd=2, ylim=c(0.95*min(m.lo),1.05*max(m.up)))
# lines(age.vec,m.lo,lty=2,col="red",lwd=1.5)
# lines(age.vec,m.up,lty=2,col="red",lwd=1.5)
# par(mfrow=c(1,1))



# Decrement founding population size and determine MVP & Pr(Qext) -------------------------------------------------


Qthresh <- 50 # quasiextinction threshold. Set at 50 (50f Ne to avoid inbreeding depression), then at whatever cons. managers want 

# sequence vector for founding N
pop.found.vec <- seq(pop.found, 10, -5) # change the increments to get a smoother line


# Set iterations ----------------------------------------------------------


iter <- 10000 # iterations to run for each projection loop change to 10,000 for final run
itdiv <- iter/1000 #final model rate at iter/1000

## set time limit for projection in 1-yr increments
yr.now <- 2020
#************************
yr.end <- yr.now + round(40*gen.l,0) # set projection end date
#************************
t <- (yr.end - yr.now) #timeframe
yr.vec <- seq(yr.now,yr.end)

# p loop storage
PrQExt <- minNmd <- minNlo <- minNup <- rep(NA, length(pop.found.vec))

for (p in 1:length(pop.found.vec)) {
  
  ## initial population vector
  ssd <- stable.stage.dist(popmat)
  init.vec <- ssd * pop.found.vec[p] #initial population vector
  
  ## stochastic projection with density feedback
  ## set storage matrices & vectors
  n.sums.mat <- qExt.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
  
  for (e in 1:iter) {
    
    popmat <- popmat.orig
    
    n.mat <- matrix(0, nrow=age.max+1,ncol=(t+1))
    n.mat[,1] <- init.vec
    
    for (i in 1:t) {
      # stochastic survival values
      s.alpha <- estBetaParams(s.vec, s.sd.vec^2)$alpha
      s.beta <- estBetaParams(s.vec, s.sd.vec^2)$beta
      s.stoch <- rbeta(length(s.alpha), s.alpha, s.beta)
      
      if (rbinom(1, 1, 0.14/gen.l) == 1) { # catastrophe
        cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
        cat.beta <- estBetaParams(0.5, 0.05^2)$beta
        s.stoch <- s.stoch * (rbeta(1, cat.alpha, cat.beta))
      }
      
      # stochastic fertility sampler (Gaussian)
      fert.stch <- rnorm(length(popmat[,1]), m.vec, m.sd.vec)
      fert.stoch <- ifelse(fert.stch < 0, 0, fert.stch)
      
      totN.i <- sum(n.mat[,i], na.rm=T)
      pred.red <- as.numeric(a.lp/(1+(totN.i/b.lp)^c.lp))
      
      diag(popmat[2:(age.max+1),]) <- (s.stoch[-(age.max+1)])*pred.red
      popmat[age.max+1,age.max+1] <- (s.stoch[age.max+1])*pred.red
      popmat[1,] <- fert.stoch
      n.mat[,i+1] <- popmat %*% n.mat[,i]
      
    } # end i loop
    
    n.sums.mat[e,] <- ((as.vector(colSums(n.mat))))
    qExt.mat[e,] <- ifelse(n.sums.mat[e,] < Qthresh, 1, 0)
    
    if (e %% itdiv==0) print(e) 
    
  } # end e loop
  
  n.md <- apply(n.sums.mat, MARGIN=2, median, na.rm=T) # median over all iterations
  n.up <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  n.lo <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
  Qmax <- apply(qExt.mat, MARGIN=1, max, na.rm=T)
  minN <- apply(n.sums.mat, MARGIN=1, min, na.rm=T) # takes the min from n.sums.mat 

  
  minNmd[p] <- median(minN)
  minNlo[p] <- quantile(minN, probs=0.025)
  minNup[p] <- quantile(minN, probs=0.975)
  
  PrQExt[p] <- sum(Qmax)/iter
  
  print(" ", quote = FALSE)
  print("*********************************************", quote = FALSE)
  print(paste("founding N = ", pop.found.vec[p]), quote = FALSE)
  print("*********************************************", quote = FALSE)
  print(" ", quote = FALSE)
  
} # end p loop

par(mfrow=c(2,1))
plot(pop.found.vec, PrQExt, type="l", xlab="founding N", ylab="Pr(quasi-extinction)")

plot(pop.found.vec, minNmd, type="l", xlab="founding N", ylab="minimum N", ylim=c(min(minNlo), max(minNup)))
lines(pop.found.vec, minNlo, lty=2, col="red")
lines(pop.found.vec, minNup, lty=2, col="red")
par(mfrow=c(1,1))

