# Sofie Costin, Corey Bradshaw, Frederik Saltre
# Global Ecology, Flinders University globalecologyflinders.com
# Chuditch Dasyurus geoffroii (DG) PVA
# requires library - Plotly
# updated 04/07/2022

# Morris K, Johnson B, Orell P, Gaikhorst G, Wayne A, Moro D. 2003. Recovery of the threatened chuditch (Dasyurus geoffroii): a case study. Predators with pouches: the biology of carnivorous marsupials:435-451.
# Stead-Richardson, E.J., Bradshaw, S.D., Bradshaw, F.J. and Gaikhorst, G., 2001. Monitoring the oestrous cycle of the chuditch (Dasyurus geoffroii)(Marsupialia: Dasyuridae): non-invasive analysis of faecal oestradiol-17b. Australian Journal of Zoology, 49(2), pp.183-193.)


## remove everything
# rm(list = ls())


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

age.max <- 3 # maximum age of females in wild from (Fauna Profile DCBA)

age.vec <- 0:age.max

# Create vectors: fertility -----------------------------------------------


# Stage 1 = < 8 months (juvenile) = 8/12 months =  0.67 of first year
# Stage 2 = 8 + months (adult) = 4/12 months = 0.33 of first year

fert <- c((3/6), (3/5), (1/2), (3/4), (3/3))  # proportion of breeding females reported by Morris et al. 2003

# probability of breeding 
juv.fert <- 0
Ad.fert <- mean(fert)

pr.breed.yr1 <- (juv.fert * 0.67) +  (Ad.fert * 0.33)
pr.breed <- Ad.fert

# stochastic young sampler
young.stoch <- rpois(age.max+1, 5) # poisson distribution around mean of 5 young.
no.young <- (young.stoch)
no.Fyoung <- no.young/2 # number of female young assuming 1:1 sex ratio

# m.vec <- c(pr.breed.yr1*no.Fyoung[1], pr.breed*no.Fyoung[2:4]) # female offspring produced per year
m.vec <- c(0.33165, 2.68, 2.01, 1.675)
m.vec

# fertility errors
pr.breed.sd <- sd(fert) 
m.sd.vec <- c(rep(pr.breed.sd, age.max+1)) #mean and standard deviations vector, juvenile and adult fertility 
m.sd.vec


# Create vectors: Survival ------------------------------------------------
# based on survival data for D. hallucatus from Rew-Duffy et al. 2020

# mass
DM.mass <- 1.68 ## 1.7 kg (Glen 2008-AJZ)

# surv.yr1a <- ((((41+35)/76) * (10/12)) + (((40+32)/76) * (2/12)))
surv.yr2a <- ((((40+32)/76) * (3/12)) + (((33+31)/76) * (5/12) + (((20+29)/76)*(4/12))))
surv.yr3a <- ((((20+29)/76) * (1/12)) + (((16+23)/76) * (5/12)) + (((9+19)/76) * (5/12)) + (((3+6)/76) * (1/12)))
surv.yr4a <- ((((3+6)/76) * (2/12)) + (((1/76)*(10/12))))

surv.yr1 <- surv.yr2a*0.25
surv.yr2 <- surv.yr2a*0.7
surv.yr3 <- surv.yr3a*0.7
surv.yr4 <- surv.yr4a*0.7

s.vec <- c(surv.yr1, surv.yr2, surv.yr3, surv.yr4)
s.vec ## female  survival # probability of surviving from one year to the next. e.g surviving fourth year of life


survivors.low<- c(0.41,0.41,0.40,0.33,0.20,0.16,0.09,0.03,0,0,0)
survivors.high<- c(0.35,0.35,0.32,0.31,0.29,0.23,0.19,0.06,0.01,0.01,0)
survivors <-cbind.data.frame(survivors.high, survivors.low)
survivors$mean <- (survivors.high+survivors.low)/2
survivors$low <- abs((survivors$mean-survivors.high))
survivors$sd <- (survivors$low / 2) #based on the range rule

# convert to years

surv.sd.yr1 <- ((0.015 * (10/12)) + (0.02 * (2/12)))
surv.sd.yr2 <- ((0.02 * (3/12)) + (0.005 * (5/12)) + (0.0225*(4/12)))
surv.sd.yr3 <- ((0.0225 * (1/12)) + (0.0175 * (5/12)) + (0.025 * (5/12)) + (0.0075 * (1/12)))
surv.sd.yr4 <- ((0.0075 * (2/12)) + ((0.0025*(10/12))))

# <- c(rep(0.122, age.max+1))
s.sd.vec <- c(surv.sd.yr1, surv.sd.yr2, surv.sd.yr3, surv.sd.yr4) #standard deviations vector, juvenile and adult survival

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


# Project 40 Generations  -----------------------------------------------------------------


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

red.vec <- c(1,0.999,0.995,0.98,0.95)
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
pop.found.vec <- seq(pop.found, 0, -10) # change the increments down to get a smoother line


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

# write.csv(s.stor.mat, file = "ChuditchSStorMat4.csv")
# write.csv(fert.stor.mat, file = "ChuditchFertStorMat4.csv")

par(mfrow=c(2,1))
plot(pop.found.vec, PrQExt, type="l", xlab="founding N", ylab="Pr(quasi-extinction)")
plot(pop.found.vec, minNmd, type="l", xlab="founding N", ylab="minimum N", ylim=c(min(minNlo), max(minNup)))
lines(pop.found.vec, minNlo, lty=2, col="red")
lines(pop.found.vec, minNup, lty=2, col="red")
par(mfrow=c(1,1))


# Save output -------------------------------------------------------------

alldata <- data.frame(pop.found.vec, minNmd, minNlo, minNup, PrQExt)
write.csv(alldata, file = "DG40Galldata.csv")

# final plots -------------------------------------------------------------

## plot in ggplot 
Chuditchplot <- ggplot(data = alldata, mapping = aes(x=pop.found.vec, y=PrQExt)) +
  geom_line(aes(y=PrQExt), color = "black") + 
  geom_vline(xintercept = 100, linetype = 2, color = "red") +
  scale_x_continuous(limits = c(0,150), breaks = seq(0,150, by = 10),expand = c(0,0.7)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1, by = 0.1))+
  theme_bw() +
  labs(x = "Founding N", y = "Pr(quasi-ext)")
Chuditchplot
ggsave("DG40Gplot.png")

ChuditchplotminN <- ggplot(data = alldata, mapping = aes(x=pop.found.vec, y=minNmd)) +
  geom_line(aes(y=minNmd), color = "black") + 
  geom_line(aes(y=minNlo), color = "red", linetype = 2) + 
  geom_line(aes(y=minNup), color = "red", linetype = 2) +
  scale_x_continuous(limits = c(0,520), breaks = seq(0,500, by = 50),expand = c(0,0)) +
  scale_y_continuous(limits = c(0,500), breaks = seq(0,500, by = 50))+
  theme_bw() +
  labs(x = "Founding N", y = "lowest N")
ChuditchplotminN
ggsave("DG40GplotminN.png")

alldata

# Project 100 years  -----------------------------------------------------------------


## set time limit for projection in 1-yr increments

yr.now <- 2020 # 
yr.end <- yr.now + round(100,0) # set projection end date to 100 years
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

red.vec <- c(1,0.999,0.995,0.98,0.95)
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
pop.found.vec <- seq(pop.found, 0, -10) # change the increments down to get a smoother line


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
write.csv(alldata, file = "DG100Yalldata.csv")

# final plots -------------------------------------------------------------

## plot in ggplot 
Chuditchplot <- ggplot(data = alldata, mapping = aes(x=pop.found.vec, y=PrQExt)) +
  geom_line(aes(y=PrQExt), color = "black") + 
  geom_vline(xintercept = 50, linetype = 2, color = "red") +
  scale_x_continuous(limits = c(0,150), breaks = seq(0,150, by = 10),expand = c(0,0.7)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1, by = 0.1))+
  theme_bw() +
  labs(x = "Founding N", y = "Pr(quasi-ext)")
Chuditchplot
ggsave("DG100Yplot.png")

ChuditchplotminN <- ggplot(data = alldata, mapping = aes(x=pop.found.vec, y=minNmd)) +
  geom_line(aes(y=minNmd), color = "black") + 
  geom_line(aes(y=minNlo), color = "red", linetype = 2) + 
  geom_line(aes(y=minNup), color = "red", linetype = 2) +
  scale_x_continuous(limits = c(0,520), breaks = seq(0,500, by = 50),expand = c(0,0)) +
  scale_y_continuous(limits = c(0,500), breaks = seq(0,500, by = 50))+
  theme_bw() +
  labs(x = "Founding N", y = "lowest N")
ChuditchplotminN
ggsave("DG100YplotminN.png")


