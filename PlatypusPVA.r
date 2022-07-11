# Sofie Costin, Corey Bradshaw, Frederik Saltre
# Global Ecology, Flinders University globalecologyflinders.com
# Platypus Ornithorhynchus anatinus (OA) PVA
# requires library - Plotly
# updated 09/07/2022
# life history data from Bino, G., Grant, T. R. & Kingsford, R. T. 2015. Life history and dynamics of a platypus (Ornithorhynchus anatinus) population: four decades of mark-recapture surveys. Scientific Reports, 5:1, 16073.

## remove everything
rm(list = ls())

# libraries
library(plotly)
library(ggplot2)


# functions ---------------------------------------------------------------

# beta distribution shape parameter estimator function
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

## source/matrix operators
source("matrixOperators.r")


# create Leslie matrix ----------------------------------------------------


age.max <- 21 # maximum age of females in wild

age.vec <- 0:age.max


# create vectors: fertility -----------------------------------------------

 # 5 / 70 with 3 eggs
# 54 / 70 with 2 eggs
# 11 / 70 with 1 egg
# 300 of 812 lactating at any one time (0.37)

pr.breed <- 0.62
pr.breed.sd <- 0.44/10
no.egg <- ((5*3) + (54*2) + (11*1))/70
no.Fegg <- no.egg/2
m.vec <- c(0,0,rep(pr.breed*no.Fegg, age.max-1))
m.vec

breed.stoch <- rbeta(1, estBetaParams(pr.breed, pr.breed.sd^2)$alpha, estBetaParams(pr.breed, pr.breed.sd^2)$beta)
egg.stoch <- sample(c(1:3),1,replace=T,prob=c(11/70,54/70,5/70)) # stochastic egg sampler

m.stoch <- breed.stoch*egg.stoch
m.stoch

plot(0:age.max,m.vec,pch=19,type="b")

# fertility errors
m.sd.vec <- c(0, 0, rep(0.15, age.max-1)) #mean and standard deviations vector, juvenile and adult fertility 

# estimated max r
OA.mass <- .8641 # overall mean adult female weight from 1973 - 2014
rm.allom.pred <- 10^(0.6914 - (0.2622*log10(OA.mass*1000)))
lm.pred <- exp(rm.allom.pred)
lm.pred

# create vectors: survival ----------------------------------------------------------------

j.suv <- 0.27/0.6 # juvenile female survival 
a.suv <- 0.76 # adult female survival
s.vec <- c(j.suv, j.suv, rep(a.suv,age.max-1))
plot(0:age.max,s.vec,pch=19,type="b")
jsupd <- j.suv+(0.154*j.suv) # incorporates juvenile dispersal 
ssupd <- a.suv+(0.145*a.suv) # incorporates adult dispersal 
s.vec.upd <- c(jsupd, jsupd, rep(ssupd,age.max-1))
plot(0:age.max,s.vec.upd,pch=19,type="b")

# survival errors
s.sd.vec <- c(0.04, 0.04, rep(0.05,age.max-1)) #mean and standard deviations vector, juvenile and adult survival


# create matrix -----------------------------------------------------------

popmat <- matrix(data = 0, nrow=age.max+1, ncol=age.max+1) # creates a matrix where data = 0 and dimensions age.max x age.max
diag(popmat[2:(age.max+1),]) <- s.vec.upd[1:(age.max)] #  populates popmat diagonally with s.vec
popmat[age.max+1,age.max+1] <- s.vec.upd[age.max+1] 
popmat[1,] <- m.vec ## row 1 of popmat populated with m.vec
popmat.orig <- popmat ## save original matrix
popmat

## matrix properties. Functions from "matrixOperators.R"
max.lambda(popmat) ## 1-yr lambda
max.r(popmat) # rate of population change, 1-yr
ssd <- stable.stage.dist(popmat) ## stable stage distribution

plot(0:age.max,ssd,pch=19,type="b")

R.val(popmat, age.max) # reproductive value
gen.l <- G.val(popmat, age.max) # mean generation length


## matrix properties
max.lambda(popmat) ## 1-yr lambda
max.r(popmat) # rate of population change, 1-yr
stable.stage.dist(popmat) ## stable stage distribution
R.val(popmat, age.max) # reproductive value
gen.l <- G.val(popmat, age.max) # mean generation length


# initial population vector -----------------------------------------------

pop.found <- 500 # change this to change initial population size 
ssd <- stable.stage.dist(popmat) # stable stage distribution
init.vec <- ssd * pop.found #initial population vector


# PROJECT 40 GENERATIONS -----------------------------------------------------------------


## set time limit for projection in 1-yr increments
yr.now <- 2020 # 
yr.end <- yr.now + round(40*gen.l,0) # set projection end date to 40 generations
t <- (yr.end - yr.now) #time frame

tot.F <- sum(popmat.orig[1,])
popmat <- popmat.orig # reset matrix 
yr.vec <- seq(yr.now,yr.end) # year vector

## set population storage matrices
n.mat <- matrix(0, nrow=age.max+1,ncol=(t+1))  #empty matrix
n.mat[,1] <- init.vec # fill first matrix column with initial population vector

## set up projection loop
for (i in 1:t) {
  n.mat[,i+1] <- popmat %*% n.mat[,i]
}

n.pred <- colSums(n.mat) # number of animals through time period, no density reduction treatment, no carrying capacity
yrs <- seq(yr.now, yr.end, 1)
plot(yrs, n.pred,type="b",lty=2,pch=19,xlab="year",ylab="N")

# compensatory density feedback
# population rate of increase relative to carrying capacity. Larger distance between population and K = faster population growth
K.max <- pop.found
K.vec <- c(1,K.max*0.5,K.max*0.83,K.max*0.90,K.max)
plot(K.vec)

red.vec <- c(1,0.999,0.99,0.985,0.976) # reduction vector
plot(K.vec,red.vec,pch=19,type="b")
Kred.dat <- data.frame(K.vec,red.vec)

# logistic power function a/(1+(x/b)^c) 
# fits logistic power function to population relative to carrying capacity, K
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
  diag(popmat[2:(age.max+1),]) <- s.vec.upd[1:(age.max)]*pred.red
  popmat[age.max+1,age.max+1] <- s.vec.upd[age.max+1]
  popmat[1,] <- m.vec
  n.mat[,i+1] <- popmat %*% n.mat[,i]
}

n.pred <- colSums(n.mat)
plot(yrs, n.pred,type="l",lty=2,pch=19,xlab="year",ylab="N") # untreated population increases, rate of increase relative to K, no stochastic sampling
abline(h=K.max, lty=2, col="red") # carrying capacity


# decrement founding population size and determine MVP and Pr(Qext) -------

# quasi-extinction threshold -----------------------------------------------
Qthresh <- 25 # quasiextinction threshold. Set at 25 (25f Ne to avoid inbreeding depression according to 50/500 rule), then at whatever cons. managers want 

# sequence vector for founding N
pop.found.vec <- seq(pop.found, 0, -5) # change the increments down to get a smoother line

# set iterations ----------------------------------------------------------

iter <- 10000 # iterations to run for each projection loop. Change to 10,000 for final run
itdiv <- iter/1000 #final model rate at iter/1000

## set time limit for projection in 1-yr increments
yr.now <- 2020 # 
yr.end <- yr.now + round(40*gen.l,0) # set projection end date to 40 generations
t <- (yr.end - yr.now) #time frame

# p loop storage ----------------------------------------------------------

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
      s.alpha <- estBetaParams(s.vec.upd, s.sd.vec^2)$alpha
      s.beta <- estBetaParams(s.vec.upd, s.sd.vec^2)$beta
      s.stoch <- rbeta(length(s.alpha), s.alpha, s.beta)
      
      if (rbinom(1, 1, 0.14/gen.l) == 1) { # catastrophe
        s.stoch <- 0.5 * s.stoch
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
  
  n.md <- apply(n.sums.mat, MARGIN=2, median, na.rm=T) # median n over all iterations
  n.up <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper n over all iterations
  n.lo <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower n over all iterations
  
  Qmax <- apply(qExt.mat, MARGIN=1, max, na.rm=T)
  minN <- apply(n.sums.mat, MARGIN=1, min, na.rm=T)
  
  
  minNmd[p] <- median(minN)
  minNlo[p] <- quantile(minN, probs=0.025)
  minNup[p] <- quantile(minN, probs=0.975)
  
  PrQExt[p] <- sum(Qmax)/iter
  
  print("", quote = FALSE)
  print("################", quote = FALSE)
  print(paste("founding N = ", pop.found.vec[p]), quote = FALSE)
  print("################", quote = FALSE)
  print("", quote = FALSE)
  
} # end p loop


# create plots -------------------------------------------------------------------

par(mfrow=c(2,1))
plot(pop.found.vec, PrQExt, type="l", xlab="founding N", ylab="Pr(quasi-extinction)")
plot(pop.found.vec, minNmd, type="l", xlab="founding N", ylab="lowest N", ylim=c(min(minNlo), max(minNup)))
lines(pop.found.vec, minNlo, lty=2, col="red")
lines(pop.found.vec, minNup, lty=2, col="red")
par(mfrow=c(1,1))


# save output -------------------------------------------------------------
alldata <- data.frame(pop.found.vec, minNmd, minNlo, minNup, PrQExt)
write.csv(alldata, file = "OA40Galldata.csv")

# final plots -------------------------------------------------------------

## plot in ggplot 
platyplot <- ggplot(data = alldata, mapping = aes(x=pop.found.vec, y=PrQExt)) +
  geom_line(aes(y=PrQExt), color = "black") + 
  geom_vline(xintercept = 160, linetype = 2, color = "red") +
  scale_x_continuous(limits = c(0,200), breaks = seq(0,200, by = 10),expand = c(0,0.5)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1, by = 0.1), expand = c(0,0))+
  theme_bw() +
  labs(x = "Founding N", y = "Pr(quasi-ext)")
platyplot
ggsave("OA40Gplot.png")

platyplotplotminN <- ggplot(data = alldata, mapping = aes(x=pop.found.vec, y=minNmd)) +
  geom_line(aes(y=minNmd), color = "black") + 
  geom_line(aes(y=minNlo), color = "red", linetype = 2) + 
  geom_line(aes(y=minNup), color = "red", linetype = 2) +
  scale_x_continuous(limits = c(0,pop.found+10), breaks = seq(0,pop.found, by =50),expand = c(0,0)) +
  scale_y_continuous(limits = c(0,minNup+10), expand = c(0,0))+
  theme_bw() +
  labs(x = "Founding N", y = "lowest N")
platyplotplotminN
ggsave("OA40GplotplotminN.png")


# PROJECT 100 YEARS -----------------------------------------------------------------

## set time limit for projection in 1-yr increments
yr.now <- 2020 # 
yr.end <- yr.now + round(100,0) # set projection end date to 100 years (IUCN)
t <- (yr.end - yr.now) #time frame

tot.F <- sum(popmat.orig[1,])
popmat <- popmat.orig # reset matrix 
yr.vec <- seq(yr.now,yr.end) # year vector

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
K.max <- pop.found
K.vec <- c(1,K.max*0.5,K.max*0.83,K.max*0.90,K.max)
plot(K.vec)

red.vec <- c(1,0.999,0.99,0.985,0.976)
plot(K.vec,red.vec,pch=19,type="b")
Kred.dat <- data.frame(K.vec,red.vec)

# logistic power function a/(1+(x/b)^c) 
# fits logistic power function to population relative to carrying capacity, K
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
  diag(popmat[2:(age.max+1),]) <- s.vec.upd[1:(age.max)]*pred.red
  popmat[age.max+1,age.max+1] <- s.vec.upd[age.max+1]
  popmat[1,] <- m.vec
  n.mat[,i+1] <- popmat %*% n.mat[,i]
}

n.pred <- colSums(n.mat)
plot(yrs, n.pred,type="l",lty=2,pch=19,xlab="year",ylab="N") # untreated population increases, rate of increase relative to K, no stochastic sampling
abline(h=K.max, lty=2, col="red") # carrying capacity

# decrement founding population size and determine MVP and Pr(Qext) -------

# quasiextinction threshold -----------------------------------------------
Qthresh <- 25 # quasi-extinction threshold. Set at 25 (25f Ne to avoid inbreeding depression), then at whatever cons. managers want 

# sequence vector for founding N
pop.found.vec <- seq(pop.found, 0, -5) # change the increments down to get a smoother line

# set iterations ----------------------------------------------------------

iter <- 10000 # iterations to run for each projection loop change to 10,000 for final run
itdiv <- iter/1000 #final model rate at iter/1000


# p loop storage ----------------------------------------------------------

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
      s.alpha <- estBetaParams(s.vec.upd, s.sd.vec^2)$alpha
      s.beta <- estBetaParams(s.vec.upd, s.sd.vec^2)$beta
      s.stoch <- rbeta(length(s.alpha), s.alpha, s.beta)
      
      if (rbinom(1, 1, 0.14/gen.l) == 1) { # catastrophe
        s.stoch <- 0.5 * s.stoch
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
  
  n.md <- apply(n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
  n.up <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  n.lo <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
  Qmax <- apply(qExt.mat, MARGIN=1, max, na.rm=T)
  minN <- apply(n.sums.mat, MARGIN=1, min, na.rm=T)
  
  
  minNmd[p] <- median(minN)
  minNlo[p] <- quantile(minN, probs=0.025)
  minNup[p] <- quantile(minN, probs=0.975)
  
  PrQExt[p] <- sum(Qmax)/iter
  
  print("", quote = FALSE)
  print("################", quote = FALSE)
  print(paste("founding N = ", pop.found.vec[p]), quote = FALSE)
  print("################", quote = FALSE)
  print("", quote = FALSE)
  
} # end p loop


# plots -------------------------------------------------------------------

par(mfrow=c(2,1))
plot(pop.found.vec, PrQExt, type="l", xlab="founding N", ylab="Pr(quasi-extinction)")

plot(pop.found.vec, minNmd, type="l", xlab="founding N", ylab="lowest N", ylim=c(min(minNlo), max(minNup)))
lines(pop.found.vec, minNlo, lty=2, col="red")
lines(pop.found.vec, minNup, lty=2, col="red")
par(mfrow=c(1,1))


# save output -------------------------------------------------------------

alldata <- data.frame(pop.found.vec, minNmd, minNlo, minNup, PrQExt)
write.csv(alldata, file = "OA100Yalldata.csv")

# final plots -------------------------------------------------------------

## plot in ggplot 
OAplot <- ggplot(data = alldata, mapping = aes(x=pop.found.vec, y=PrQExt)) +
  geom_line(aes(y=PrQExt), color = "black") + 
  geom_vline(xintercept = 160, linetype = 2, color = "red") +
  scale_x_continuous(limits = c(0,200), breaks = seq(0,200, by = 10),expand = c(0,0.5)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1, by = 0.1), expand = c(0,0))+
  theme_bw() +
  labs(x = "Founding N", y = "Pr(quasi-ext)")
OAplot
ggsave("OA100Yplot.png")

OAplotminN <- ggplot(data = alldata, mapping = aes(x=pop.found.vec, y=minNmd)) +
  geom_line(aes(y=minNmd), color = "black") + 
  geom_line(aes(y=minNlo), color = "red", linetype = 2) + 
  geom_line(aes(y=minNup), color = "red", linetype = 2) +
  scale_x_continuous(limits = c(0,pop.found+10), breaks = seq(0,pop.found, by =50),expand = c(0,0)) +
  scale_y_continuous(limits = c(0,minNup+10), expand = c(0,0))+
  theme_bw() +
  labs(x = "Founding N", y = "lowest N")
OAplotminN
ggsave("OA100YplotminN.png")
