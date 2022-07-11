# Sofie Costin, Corey Bradshaw, Frederik Saltre
# Global Ecology, Flinders University globalecologyflinders.com
# Brush-tailed bettong Bettongia penicillata (BP) PVA
# requires library - Plotly
# update 03/07/2022
References:
Pacioni C, Williams MR, Lacy RC, Spencer PB, Wayne AF. 2017. Predators and genetic fitness: key threatening factors for the conservation of a bettong species. Pacific Conservation Biology 23:200-212.
Thompson CK, Wayne AF, Godfrey SS, Thompson RA. 2015. Survival, age estimation and sexual maturity of pouch young of the brush-tailed bettong (Bettongia penicillata) in captivity. Australian Mammalogy 37:29-38.


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

# Leslie matrix -------------------------------------

age.max = 6 # maximum age of females in wild 

age.vec <- 0:age.max

# Create vectors: fertility -----------------------------------------------

# Stages max values from Thompson et al 2015

# stage 1 (pouch young) = 110/365 days = 0.30 of 1 year
# Stage 2 (independent young) =  70/365 days  = 0.19 of 1 year
# stage 3 (sexually active) = 185/365 days = 0.51 of 1 year

# survival and mortality from Pacioni et al 2017 (Supplementary 1)

m.vec <- c((0.51*0.893), rep(0.893, age.max)) ## female offspring produced each year
plot(0:5,m.vec,pch=19,type="b")

# fertility errors
m.sd.vec <- c(rep(0.35, age.max+1)) #mean and standard deviations vector, juvenile and adult fertility 

# Create vectors: Survival ------------------------------------------------

#survival
s.vec <- c((((0.956*0.3)+(0.85*0.19)+(0.961*0.51))*0.52), rep(0.961*0.52, age.max))  ## female  survival # probability of surviving from one year to the next. e.g surviving fourth year of life
#s.vec <- c(0.4, rep(0.8, age.max-1), 0.1)  ## female  survival # probability of surviving from one year to the next. e.g surviving fourth year of life

# survival errors
s.sd.vec <- c(rep(0.15, age.max+1)) #mean and standard deviations vector, juvenile and adult survival

# Create matrix -----------------------------------------------------------

popmat <- matrix(data = 0, nrow=age.max+1, ncol=age.max+1) # creates a matrix where data = 0 and dimensions age.max x age.max
diag(popmat[2:(age.max+1),]) <- s.vec[1:(age.max)] #  populates popmat diagonally with s.vec
# popmat[age.max+1,age.max+1] <- s.vec[1:(age.max)] # would fill 7,7 with s.vec, but not needed because pr(breed) at age.max = 0
popmat[1,] <- m.vec ## row 1 of popmat populated with m.vec
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
plot(0:age.max,ssd,pch=19,type="b")


# Project 40 generations-----------------------------------------------------------------


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
K.vec <- c(1,K.max*0.8,0.9*K.max,0.95*K.max,0.99*K.max)
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


# Decrement founding population size and determine MVP & Pr(Qext) -------------------------------------------------

Qthresh <- 25 # quasiextinction threshold. Set at 25 (25f Ne to avoid inbreeding depression), then at whatever cons. managers want 

# sequence vector for founding N
pop.found.vec <- seq(pop.found, 0, -5) # change the increments down to get a smoother line


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
s.stor.mat <- matrix(data = NA, nrow = iter, ncol = t)
fert.stor.mat <- matrix(data = NA, nrow = iter, ncol = t)

for (p in 1:length(pop.found.vec)) {
  
  ## initial population vector
  popmat <- popmat.orig
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
      s.stoch1 <- rbeta(length(s.alpha), s.alpha, s.beta)
      s.stoch <- ifelse(is.na(s.stoch1) == T, 0, s.stoch1)
      
      if (rbinom(1, 1, 0.14/gen.l) == 1) { # catastrophe
        s.stoch <- 0.5 * s.stoch
      }
      s.stor.mat[e,i] <- s.stoch[4]
      
      # print(s.stoch)
      # stochastic fertility sampler (Gaussian)
      fert.stch <- rnorm(length(popmat[,1]), m.vec, m.sd.vec)
      fert.stoch <- ifelse(fert.stch < 0, 0, fert.stch)
      fert.stor.mat[e,i] <- fert.stoch[4]
      
      # print(fert.stoch)
      totN.i <- sum(n.mat[,i], na.rm=T)
      pred.red <- as.numeric(a.lp/(1+(totN.i/b.lp)^c.lp))
      diag(popmat[2:(age.max+1),]) <- (s.stoch[-(age.max+1)])*pred.red
      popmat[age.max+1,age.max+1] <- (s.stoch[age.max+1])*pred.red
      popmat[1,] <- fert.stoch
      n.mat[,i+1] <- popmat %*% n.mat[,i]
    } # end i loop
    
    n.sums.mat[e,] <- ((as.vector(colSums(n.mat, na.rm=T))))
    qExt.mat[e,] <- ifelse(n.sums.mat[e,] < Qthresh, 1, 0)
    
    if (e %% itdiv==0) print(e)
  } # end e loop
  
  n.md <- apply(n.sums.mat, MARGIN=2, median, na.rm=T) # median over all iterations
  n.up <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  n.lo <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
  plot(yr.vec,n.md,type="l")
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


# Save output -------------------------------------------------------------

# write.csv(s.stor.mat, file = "BettongSStorMat4.csv")
# write.csv(fert.stor.mat, file = "BettongFertStorMat4.csv")

alldata <- data.frame(pop.found.vec, minNmd, minNlo, minNup, PrQExt)
write.csv(alldata, file = "BP40Galldata.csv")
 
 # final plots -------------------------------------------------------------

 ## plot in ggplot 
Bettongplot <- ggplot(data = alldata, mapping = aes(x=pop.found.vec, y=PrQExt)) +
 geom_line(aes(y=PrQExt), color = "black") + 
 geom_vline(xintercept = 185, linetype = 2, color = "red") +
 scale_x_continuous(limits = c(0,250), breaks = seq(0,250, by = 10),expand = c(0,0.7)) +
 scale_y_continuous(limits = c(0,1), breaks = seq(0,1, by = 0.1), expand = c(0,0))+
 theme_bw() +
 labs(x = "Founding N", y = "Pr(quasi-ext)")
Bettongplot
ggsave("BP40Gplot.png")

BettongplotminN <- ggplot(data = alldata, mapping = aes(x=pop.found.vec, y=minNmd)) +
 geom_line(aes(y=minNmd), color = "black") + 
 geom_line(aes(y=minNlo), color = "red", linetype = 2) + 
 geom_line(aes(y=minNup), color = "red", linetype = 2) +
 scale_x_continuous(limits = c(0,520), breaks = seq(0,500, by = 50),expand = c(0,0)) +
 scale_y_continuous(limits = c(0,500), breaks = seq(0,500, by = 50))+
 theme_bw() +
 labs(x = "Founding N", y = "lowest N")
BettongplotminN
ggsave("BP40GplotminN.png")

alldata

# Project 100 years -----------------------------------------------------------------


## set time limit for projection in 1-yr increments

yr.now <- 2020 # 
yr.end <- yr.now + 100 # set projection end date to 40 generations (IUCN)
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
K.vec <- c(1,K.max*0.8,0.9*K.max,0.95*K.max,0.99*K.max)
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


# Decrement founding population size and determine MVP & Pr(Qext) -------------------------------------------------

Qthresh <- 25 # quasiextinction threshold. Set at 25 (25f Ne to avoid inbreeding depression), then at whatever cons. managers want 

# sequence vector for founding N
pop.found.vec <- seq(pop.found, 0, -5) # change the increments down to get a smoother line


# Set iterations ----------------------------------------------------------


iter <- 10000 # iterations to run for each projection loop change to 10,000 for final run
itdiv <- iter/1000 #final model rate at iter/1000


# p loop storage
PrQExt <- minNmd <- minNlo <- minNup <- rep(NA, length(pop.found.vec))
s.stor.mat <- matrix(data = NA, nrow = iter, ncol = t)
fert.stor.mat <- matrix(data = NA, nrow = iter, ncol = t)

for (p in 1:length(pop.found.vec)) {
  
  ## initial population vector
  popmat <- popmat.orig
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
      s.stoch1 <- rbeta(length(s.alpha), s.alpha, s.beta)
      s.stoch <- ifelse(is.na(s.stoch1) == T, 0, s.stoch1)
      
      if (rbinom(1, 1, 0.14/gen.l) == 1) { # catastrophe
        s.stoch <- 0.5 * s.stoch
      }
      s.stor.mat[e,i] <- s.stoch[4]
      
      # print(s.stoch)
      # stochastic fertility sampler (Gaussian)
      fert.stch <- rnorm(length(popmat[,1]), m.vec, m.sd.vec)
      fert.stoch <- ifelse(fert.stch < 0, 0, fert.stch)
      fert.stor.mat[e,i] <- fert.stoch[4]
      
      # print(fert.stoch)
      totN.i <- sum(n.mat[,i], na.rm=T)
      pred.red <- as.numeric(a.lp/(1+(totN.i/b.lp)^c.lp))
      diag(popmat[2:(age.max+1),]) <- (s.stoch[-(age.max+1)])*pred.red
      popmat[age.max+1,age.max+1] <- (s.stoch[age.max+1])*pred.red
      popmat[1,] <- fert.stoch
      n.mat[,i+1] <- popmat %*% n.mat[,i]
    } # end i loop
    
    n.sums.mat[e,] <- ((as.vector(colSums(n.mat, na.rm=T))))
    qExt.mat[e,] <- ifelse(n.sums.mat[e,] < Qthresh, 1, 0)
    
    if (e %% itdiv==0) print(e)
  } # end e loop
  
  n.md <- apply(n.sums.mat, MARGIN=2, median, na.rm=T) # median over all iterations
  n.up <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  n.lo <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
  plot(yr.vec,n.md,type="l")
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


# Save output -------------------------------------------------------------

# write.csv(s.stor.mat, file = "BettongSStorMat4.csv")
# write.csv(fert.stor.mat, file = "BettongFertStorMat4.csv")

alldata <- data.frame(pop.found.vec, minNmd, minNlo, minNup, PrQExt)
write.csv(alldata, file = "BP100Yalldata.csv")

# final plots -------------------------------------------------------------

## plot in ggplot 
Bettongplot <- ggplot(data = alldata, mapping = aes(x=pop.found.vec, y=PrQExt)) +
  geom_line(aes(y=PrQExt), color = "black") + 
  geom_vline(xintercept = 185, linetype = 2, color = "red") +
  scale_x_continuous(limits = c(0,250), breaks = seq(0,250, by = 10),expand = c(0,0.7)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1, by = 0.1), expand = c(0,0))+
  theme_bw() +
  labs(x = "Founding N", y = "Pr(quasi-ext)")
Bettongplot
ggsave("BP100Yplot.png")

BettongplotminN <- ggplot(data = alldata, mapping = aes(x=pop.found.vec, y=minNmd)) +
  geom_line(aes(y=minNmd), color = "black") + 
  geom_line(aes(y=minNlo), color = "red", linetype = 2) + 
  geom_line(aes(y=minNup), color = "red", linetype = 2) +
  scale_x_continuous(limits = c(0,520), breaks = seq(0,500, by = 50),expand = c(0,0)) +
  scale_y_continuous(limits = c(0,500), breaks = seq(0,500, by = 50))+
  theme_bw() +
  labs(x = "Founding N", y = "lowest N")
BettongplotminN
ggsave("BP100YplotminN.png")

alldata

