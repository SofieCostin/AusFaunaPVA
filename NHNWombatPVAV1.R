# Sofie Costin, Corey Bradshaw, Frederik Saltre
# Global Ecology, Flinders University globalecologyflinders.com
# Northern hairy-nosed wombat Lasiorhinus krefftii (LK) PVA
# requires library - Plotly
# update 23/06/2022

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

LK.age.max = 30 # maximum age of females in wild from Johnson 1991

LK.age.vec <- 0:LK.age.max

LK.mass <- 32.5 #(Martin & Carver 2021)

# Create vectors: fertility -----------------------------------------------
# male to female ratio 2.25:1 (Banks et al 2003)
# females breed 2 in every 3 years (Crossman 1994)
# 50 - 75% of females breeding any one year
# one offspring per year (Johnson 1991)
# sexual maturity at 3 years (Johnson 1991)

LK.pr.breed <- 0.6666
LK.breed.f <- LK.pr.breed/5 # females. sex ration 1:1 historically


LK.m.vec <- c(rep(0, 3), rep(LK.breed.f, LK.age.max-2)) ## female offspring produced each year
plot(0:LK.age.max,LK.m.vec,pch=19,type="b")

# fertility errors
LK.m.sd.vec <- c(rep(0.2, LK.age.max+1)) #mean and standard deviations vector, juvenile and adult fertility #mean SD?

# Create vectors: Survival ------------------------------------------------

LK.s.vec <- c(0.43, rep(0.95, LK.age.max))
LK.s.sd.vec <- c(rep(0.14, LK.age.max+1))

# Create matrix -----------------------------------------------------------

popmat <- matrix(data = 0, nrow=LK.age.max+1, ncol=LK.age.max+1) # creates a matrix where data = 0 and dimensions age.max x age.max
diag(popmat[2:(LK.age.max+1),]) <- LK.s.vec[1:(LK.age.max)] #  populates popmat diagonally with s.vec
popmat[1,] <- LK.m.vec ## row 1 of popmat populated with m.vec
LK.popmat.orig <- popmat ## save original matrix
popmat

## matrix properties
max.lambda(popmat) ## 1-yr lambda
max.r(popmat) # rate of population change, 1-yr
ssd <- stable.stage.dist(popmat) ## stable stage distribution
R.val(popmat, LK.age.max) # reproductive value
gen.l <- G.val(popmat, LK.age.max) # mean generation length


# Initial population vector -----------------------------------------------


LK.pop.found <- 500 # change this to change mvp 
LK.init.vec <- ssd * LK.pop.found #initial population vector
plot(0:LK.age.max,ssd,pch=19,type="b")


# Project -----------------------------------------------------------------


## set time limit for projection in 1-yr increments

yr.now <- 2020 # 
yr.end <- yr.now + round(40*gen.l,0) # set projection end date to 40 generations (IUCN)
t <- (yr.end - yr.now) #time frame
yr.vec <- seq(yr.now,yr.end) # year vector

tot.F <- sum(LK.popmat.orig[1,])
popmat <- LK.popmat.orig # reset matrix 

## set population storage matrices
n.mat <- matrix(0, nrow=LK.age.max+1,ncol=(t+1))  #empty matrix
n.mat[,1] <- LK.init.vec # fill first matrix column with initial population vector

## set up projection loop
for (i in 1:t) {
  n.mat[,i+1] <- popmat %*% n.mat[,i]
}

n.pred <- colSums(n.mat) # number of animals through time period, no density reduction treatment, no carry capacity
yrs <- seq(yr.now, yr.end, 1)
plot(yrs, n.pred,type="b",lty=2,pch=19,xlab="year",ylab="N")

# compensatory density feedback
LK.K.max <- 2*LK.pop.found
LK.K.vec <- c(1, LK.K.max/2, 0.75*LK.K.max, LK.K.max) 
LK.red.vec <- c(1,0.985,0.951,0.873)
plot(LK.K.vec, LK.red.vec,pch=19,type="b")
LK.Kred.dat <- data.frame(LK.K.vec, LK.red.vec)

# logistic power function a/(1+(x/b)^c)
LK.param.init <- c(1, 2*LK.K.max, 2)
LK.fit.lp <- nls(LK.red.vec ~ a/(1+(LK.K.vec/b)^c), 
                 data = LK.Kred.dat,
                 algorithm = "port",
                 start = c(a = LK.param.init[1], b = LK.param.init[2], c = LK.param.init[3]),
                 trace = TRUE,      
                 nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
LK.fit.lp.summ <- summary(LK.fit.lp)
plot(LK.K.vec, LK.red.vec, pch=19,xlab="N",ylab="reduction factor")
LK.K.vec.cont <- seq(1,2*LK.pop.found,1)
LK.pred.lp.fx <- coef(LK.fit.lp)[1]/(1+(LK.K.vec.cont/coef(LK.fit.lp)[2])^coef(LK.fit.lp)[3])
lines(LK.K.vec.cont, LK.pred.lp.fx, lty=3,lwd=3,col="red")

LK.a.lp <- coef(LK.fit.lp)[1]
LK.b.lp <- coef(LK.fit.lp)[2]
LK.c.lp <- coef(LK.fit.lp)[3]

## compensatory density-feedback deterministic model
## set population storage matrices
LK.n.mat <- matrix(0, nrow=LK.age.max+1, ncol=(t+1))
LK.n.mat[,1] <- LK.init.vec
LK.popmat <- LK.popmat.orig

## set up projection loop
for (i in 1:t) {
  LK.totN.i <- sum(LK.n.mat[,i])
  LK.pred.red <- as.numeric(LK.a.lp/(1+(LK.totN.i/LK.b.lp)^LK.c.lp))
  diag(LK.popmat[2:(LK.age.max+1),]) <- (LK.s.vec[-(LK.age.max+1)])*LK.pred.red
  LK.popmat[LK.age.max+1,LK.age.max+1] <- 0
  LK.popmat[1,] <- LK.m.vec
  LK.n.mat[,i+1] <- LK.popmat %*% LK.n.mat[,i]
}

LK.n.pred <- colSums(LK.n.mat)
plot(yrs, LK.n.pred, type="l",lty=2,pch=19,xlab="year",ylab="N")
abline(h=LK.pop.found, lty=2, col="red", lwd=2)


# Decrement founding population size and determine MVP & Pr(Qext) -------------------------------------------------

Qthresh <- 25 # quasiextinction threshold. Set at 50 (50f Ne to avoid inbreeding depression), then at whatever cons. managers want 

# sequence vector for founding N
pop.found.vec <- seq(LK.pop.found, 0, -5) # change the increments down to get a smoother line


# Set iterations ----------------------------------------------------------


iter <- 100 # iterations to run for each projection loop change to 10,000 for final run
itdiv <- iter/10 #final model rate at iter/1000

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
  popmat <- LK.popmat.orig
  init.vec <- ssd * pop.found.vec[p] #initial population vector
  
  ## stochastic projection with density feedback
  ## set storage matrices & vectors
  n.sums.mat <- qExt.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
  
  for (e in 1:iter) {
    popmat <- LK.popmat.orig
    n.mat <- matrix(0, nrow=LK.age.max+1,ncol=(t+1))
    n.mat[,1] <- init.vec
    for (i in 1:t) {
      
      # stochastic survival values
      s.alpha <- estBetaParams(LK.s.vec, LK.s.sd.vec^2)$alpha
      s.beta <- estBetaParams(LK.s.vec, LK.s.sd.vec^2)$beta
      s.stoch1 <- rbeta(length(s.alpha), s.alpha, s.beta)
      s.stoch <- ifelse(is.na(s.stoch1) == T, 0, s.stoch1)
      
      if (rbinom(1, 1, 0.14/gen.l) == 1) { # catastrophe
        s.stoch <- 0.5 * s.stoch
      }
      s.stor.mat[e,i] <- s.stoch[4]
      
      # print(s.stoch)
      # stochastic fertility sampler (Gaussian)
      fert.stch <- rnorm(length(popmat[,1]), LK.m.vec, LK.m.sd.vec)
      fert.stoch <- ifelse(fert.stch < 0, 0, fert.stch)
      fert.stor.mat[e,i] <- fert.stoch[4]
      
      # print(fert.stoch)
      totN.i <- sum(n.mat[,i], na.rm=T)
      pred.red <- as.numeric(LK.a.lp/(1+(totN.i/LK.b.lp)^LK.c.lp))
      diag(popmat[2:(LK.age.max+1),]) <- (s.stoch[-(LK.age.max+1)])*pred.red
      popmat[LK.age.max+1,LK.age.max+1] <- (s.stoch[LK.age.max+1])*pred.red
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
plot(pop.found.vec, minNmd, type="l", xlab="founding N", ylab="lowest N", ylim=c(min(minNlo), max(minNup)))
lines(pop.found.vec, minNlo, lty=2, col="red")
lines(pop.found.vec, minNup, lty=2, col="red")
par(mfrow=c(1,1))


# Save output -------------------------------------------------------------

setwd("C:/Users/sofie/OneDrive - Flinders/Research projects/Honours/NHN wombat")

write.csv(s.stor.mat, file = "NHNWombatSStorMat4.csv")
write.csv(fert.stor.mat, file = "NHNWombatFertStorMat4.csv")

alldata <- data.frame(pop.found.vec, minNmd, minNlo, minNup, PrQExt)
write.csv(alldata, file = "NHNWombatalldata.csv")

# final plots -------------------------------------------------------------

## plot in ggplot 
NHNWombatplot <- ggplot(data = alldata, mapping = aes(x=pop.found.vec, y=PrQExt)) +
  geom_line(aes(y=PrQExt), color = "black") + 
  geom_vline(xintercept = 75, linetype = 2, color = "red") +
  scale_x_continuous(limits = c(0,150), breaks = seq(0,150, by = 10),expand = c(0,0.7)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1, by = 0.1), expand = c(0,0))+
  theme_bw() +
  labs(x = "Founding N", y = "Pr(quasi-ext)")
NHNWombatplot
ggsave("NHNWombatplot.png")

NHNWombatplotminN <- ggplot(data = alldata, mapping = aes(x=pop.found.vec, y=minNmd)) +
  geom_line(aes(y=minNmd), color = "black") + 
  geom_line(aes(y=minNlo), color = "red", linetype = 2) + 
  geom_line(aes(y=minNup), color = "red", linetype = 2) +
  scale_x_continuous(limits = c(0,520), breaks = seq(0,500, by = 50),expand = c(0,0)) +
  scale_y_continuous(limits = c(0,500), breaks = seq(0,500, by = 50))+
  theme_bw() +
  labs(x = "Founding N", y = "lowest N")
NHNWombatplotminN
ggsave("NHNWombatplotminN.png")

alldata

