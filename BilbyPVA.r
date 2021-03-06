# Sofie Costin, Corey Bradshaw, Frederik Saltre
# Global Ecology, Flinders University globalecologyflinders.com
# Bilby Macrotis Lagotis (ML) PVA
# requires library - Plotly
# updated 03/07/2022

# References:
# Southgate RI, Christie P, Bellchambers K. 2000. Breeding biology of captive, reintroduced and wild greater bilbies, Macrotis lagotis (Marsupialia : Peramelidae). Wildlife Research 27:621.
# Southgate R, Possingham H. 1995. Modelling the reintroduction of the greater bilby Macrotis lagotis using the metapopulation model analysis of the likelihood of extinction (ALEX). Biological Conservation 73:151-160.


## remove everything
#rm(list = ls())

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

age.max <- 5 # maximum age of females in wild from (Southgate, Christie & Bellchambers 2000; Southgate & Possingham 1995)

age.vec <- 0:age.max

# Create vectors: fertility -----------------------------------------------

# (Southgate, Christie, & Bellchambers 2000; Southgate & Possingham 1995)
# Stage 1 + < 3 months (Juvenile) = 3/12 months = 0.25 of first year
# Stage 2 4-6 months (Immature) = 3/12 months = 0.25 of first year
# Stage 3 7-9 months (Subadult) = 3 / 12 months = 0.25 of first year
# Stage 3 10 months + (Adult) = 3/12 months = 0.25 of first year

# probability of breeding 
juv.fert <- 0 
imm.fert <- 0 
sAd.fert <- 0.5 
Ad.fert <- 0.89

pr.breed.yr1 <- (juv.fert*0.25) + (imm.fert*0.25) + (sAd.fert*0.25) + (Ad.fert*0.25)
pr.breed <- 0.89

no.young <- ((0.19*1)+(0.61*2)+(0.08*3))
no.Fyoung <- no.young/2

m.vec <- c((pr.breed.yr1*no.Fyoung), rep((pr.breed*no.Fyoung), age.max)) # female offspring produced per year
m.vec

# fertility errors

pr.breed.sd <- 0.2 #based on mean of other spp.
m.sd.vec <- c(rep(pr.breed.sd, age.max+1)) #mean and standard deviations vector, juvenile and adult fertility 
m.sd.vec

# Create vectors: Survival ------------------------------------------------
# from Southgate, Christie, & Bellchambers 2000

juv.surv <- 0.92
imm.surv <- 0.91
sAd.surv <- 0.94
ad.surv <- 0.92


s.vec <- c((((juv.surv*0.25) + (imm.surv*0.25) + (sAd.surv*0.25) + (ad.surv*0.25))*0.6), rep(ad.surv*0.6, age.max))
s.vec ## female  survival # probability of surviving from one year to the next. e.g surviving fourth year of life

# survival errors

s.sd.vec <- c(rep(0.122, age.max+1)) #mean and standard deviations vector, juvenile and adult survival from mean of other spp.

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


# pop.found ---------------------------------------------------------------
pop.found <- 500 # change this to change mvp 
init.vec <- ssd * pop.found #initial population vector
plot(0:age.max,ssd,pch=19,type="b")


# Project 40 generations -----------------------------------------------------------------


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

# compensatory density feedback

# population rate of increase relative to carrying capacity. Larger distance between population and K = faster population growth
K.max <- pop.found*2
K.vec <- c(1,K.max*0.83,0.92*K.max,0.98*K.max,0.99*K.max)
plot(K.vec)


# red.vec -----------------------------------------------------------------
red.vec <- c(1,0.988,0.855,0.71,0.6)
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

# compensatory density-feedback deterministic model
# set population storage matrices
n.mat <- matrix(0, nrow=age.max+1,ncol=(t+1))
n.mat[,1] <- init.vec
popmat <- popmat.orig

# set up projection loop
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


# Qthresh -----------------------------------------------------------------
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


alldata <- data.frame(pop.found.vec, minNmd, minNlo, minNup, PrQExt)
alldata
write.csv(alldata, file = "ML40Galldata.csv")

# final plots -------------------------------------------------------------

## plot in ggplot
Bilbyplot <- ggplot(data = alldata, mapping = aes(x=pop.found.vec, y=PrQExt)) +
  geom_line(aes(y=PrQExt), color = "black") +
  geom_vline(xintercept = 125, linetype = 2, color = "red") +
  scale_x_continuous(limits = c(0,150), breaks = seq(0,150, by = 10),expand = c(0,0.5)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1, by = 0.1), expand = c(0,0))+
  theme_bw() +
  labs(x = "Founding N", y = "Pr(quasi-ext)")
Bilbyplot
ggsave("ML40Gplot.png")

BilbyplotminN <- ggplot(data = alldata, mapping = aes(x=pop.found.vec, y=minNmd)) +
  geom_line(aes(y=minNmd), color = "black") +
  geom_line(aes(y=minNlo), color = "red", linetype = 2) +
  geom_line(aes(y=minNup), color = "red", linetype = 2) +
  scale_x_continuous(limits = c(0,520), breaks = seq(0,500, by = 50),expand = c(0,0)) +
  scale_y_continuous(limits = c(0,500), breaks = seq(0,500, by = 50))+
  theme_bw() +
  labs(x = "Founding N", y = "lowest N")
BilbyplotminN
ggsave("ML40GplotminN.png")


# Project 100 years -----------------------------------------------------------------


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

# compensatory density feedback

# population rate of increase relative to carrying capacity. Larger distance between population and K = faster population growth
K.max <- pop.found*2
K.vec <- c(1,K.max*0.83,0.92*K.max,0.98*K.max,0.99*K.max)
plot(K.vec)


# red.vec -----------------------------------------------------------------
red.vec <- c(1,0.988,0.855,0.71,0.6)
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

# compensatory density-feedback deterministic model
# set population storage matrices
n.mat <- matrix(0, nrow=age.max+1,ncol=(t+1))
n.mat[,1] <- init.vec
popmat <- popmat.orig

# set up projection loop
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


# Qthresh -----------------------------------------------------------------
Qthresh <- 25 # quasiextinction threshold. Set at 25 (25f Ne to avoid inbreeding depression), then at whatever cons. managers want 

# sequence vector for founding N
pop.found.vec <- seq(pop.found, 0, -5) # change the increments down to get a smoother line


# Set iterations ----------------------------------------------------------
iter <- 10000 # iterations to run for each projection loop change to 10,000 for final run
itdiv <- iter/1000 #final model rate at iter/1000

## set time limit for projection in 1-yr increments
yr.now <- 2020
#************************
yr.end <- yr.now + round(100,0) # set projection end date
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

alldata <- data.frame(pop.found.vec, minNmd, minNlo, minNup, PrQExt)
alldata
write.csv(alldata, file = "ML100Yalldata5.csv")

# final plots -------------------------------------------------------------

## plot in ggplot
Bilbyplot <- ggplot(data = alldata, mapping = aes(x=pop.found.vec, y=PrQExt)) +
  geom_line(aes(y=PrQExt), color = "black") +
  geom_vline(xintercept = 125, linetype = 2, color = "red") +
  scale_x_continuous(limits = c(0,150), breaks = seq(0,150, by = 10),expand = c(0,0.5)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1, by = 0.1), expand = c(0,0))+
  theme_bw() +
  labs(x = "Founding N", y = "Pr(quasi-ext)")
Bilbyplot
ggsave("ML100Yplot5.png")

BilbyplotminN <- ggplot(data = alldata, mapping = aes(x=pop.found.vec, y=minNmd)) +
  geom_line(aes(y=minNmd), color = "black") +
  geom_line(aes(y=minNlo), color = "red", linetype = 2) +
  geom_line(aes(y=minNup), color = "red", linetype = 2) +
  scale_x_continuous(limits = c(0,520), breaks = seq(0,500, by = 50),expand = c(0,0)) +
  scale_y_continuous(limits = c(0,500), breaks = seq(0,500, by = 50))+
  theme_bw() +
  labs(x = "Founding N", y = "lowest N")
BilbyplotminN
ggsave("ML100YplotminN.png")

