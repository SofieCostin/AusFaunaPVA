# Sofie Costin, Corey Bradshaw, Frédérik Saltré
# Global Ecology, Flinders University globalecologyflinders.com
# Brush-tailed bettong PVA
# requires library - Plotly
# updated 15/2/22

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

# create Leslie matrix
age.max = 6 # maximum age of females in wild from Threatened Species Listing Advice (FIND SOURCE)

# Stages max values from Thompson et al 2015

# stage 1 (pouch young) = 110/365 days = 0.30 of 1 year
# Stage 2 (independent young) =  70/365 days  = 0.19 of 1 year
# stage 3 (sexually active) = 185/365 days = 0.51 of 1 year

# survival and mortality from Pacioni et al 2017 (Supplementary 1)

## create vectors 
#fertility 
m.vec <- c((0.51*89.3), rep(89.3, age.max-1)) ## female offspring produced each year
plot(0:5,m.vec,pch=19,type="b")

# fertility errors
m.sd.vec <- c(rep(3.5, age.max)) #mean and standard deviations vector, juvenile and adult fertility 

#survival
s.vec <- c(((0.956*0.3)+(0.85*0.19)+(0.961*0.51)), rep(0.961, age.max-1))  ## female  survival # probability of surviving from one year to the next. e.g surviving fourth year of life

# survival errors
s.sd.vec <- c(rep(0.15, age.max)) #mean and standard deviations vector, juvenile and adult survival

                                    
# create matrix
popmat <- matrix(data = 0, nrow=age.max, ncol=age.max) ##creates a matrix where data = 0 and dimensions age.max x age.max; 16 x 16
diag(popmat[2:age.max,]) <- s.vec[1:5] ## diagonally in popmat from row 2, col 1 to row 15, col 5 populated with s.vec
popmat[age.max,age.max] <- s.vec[6] ## position [6,6] 
popmat[1,] <- m.vec ## row 1 of popmat populated with m.vec
popmat.orig <- popmat ## save original matrix
popmat

## matrix properties. Functions from "matrixOperators.R"
max.lambda(popmat) ## 1-yr lambda. 
max.r(popmat) # rate of population change, 1-yr
ssd <- stable.stage.dist(popmat) ## stable stage distribution

plot(0:5,ssd,pch=19,type="b")

R.val(popmat, age.max) # reproductive value
gen.l <- G.val(popmat, age.max) # mean generation length

################################################################################################################################
# returns number of female offspring produced per female during its lifetime as 2214.879, mean generation time 1.996848
# this is an issue, max should be 3 young x 6 years = ~18 young
################################################################################################################################

## matrix properties
max.lambda(popmat) ## 1-yr lambda
max.r(popmat) # rate of population change, 1-yr
stable.stage.dist(popmat) ## stable stage distribution
R.val(popmat, age.max) # reproductive value
gen.l <- G.val(popmat, age.max) # mean generation length

## initial population vector
pop.found <- 50 # change this to change mvp 
ssd <- stable.stage.dist(popmat)
init.vec <- ssd * pop.found #initial population vector

#################
## project
## set time limit for projection in 1-yr increments
yr.now <- 2020 # update if more data available post-2010
#************************
yr.end <- 2030 # set projection end date
#************************
t <- (yr.end - yr.now) #timeframe

tot.F <- sum(popmat.orig[1,])
popmat <- popmat.orig #resets matrix 
yr.vec <- seq(yr.now,yr.end) #year vector, 2020, 2021, 2022...

## set population storage matrices
n.mat <- matrix(0, nrow=age.max,ncol=(t+1)) #empty matrix
n.mat[,1] <- init.vec #fill first matrix column with initial population vector

## set up projection loop
for (i in 1:t) {
  n.mat[,i+1] <- popmat %*% n.mat[,i]
}

n.pred <- colSums(n.mat) #number of predators - cats - through time period, no density reduction treatment, no carry capacity
yrs <- seq(yr.now, yr.end, 1)
plot(yrs, n.pred,type="b",lty=2,pch=19,xlab="year",ylab="N")

# compensatory density feedback # K = carry capacity
#population rate of increase relative to carry capacity. Larger distance between population and K = faster population growth
K.max <- 2*pop.found
# K.min <- 1 #not used
K.vec <- c(1,pop.found/2,pop.found,0.75*K.max,K.max) #1= K.min, .75 = red.thresh??
# red.thresh <- 0.75 #not used
red.vec <- c(1,0.965,0.89,0.79,0.71)
# red.vec <- c(1,0.99,0.97,0.95,0.93)
plot(K.vec,red.vec,pch=19,type="b")
Kred.dat <- data.frame(K.vec,red.vec)

# logistic power function a/(1+(x/b)^c) #fits logistic power function to population relative to carry capacity, K
param.init <- c(1, 100, 2.5)
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
n.mat <- matrix(0, nrow=age.max,ncol=(t+1))
n.mat[,1] <- init.vec
popmat <- popmat.orig

## set up projection loop
for (i in 1:t) {
  totN.i <- sum(n.mat[,i])
  pred.red <- a.lp/(1+(totN.i/b.lp)^c.lp)
  diag(popmat[2:age.max,]) <- s.vec[1:5]*pred.red
  popmat[age.max,age.max] <- s.vec[5]*pred.red
  n.mat[,i+1] <- popmat %*% n.mat[,i]
}



n.pred <- colSums(n.mat)
plot(yrs, n.pred,type="b",lty=2,pch=19,xlab="year",ylab="N",ylim=c(0,1.05*K.max)) #untreated population increases, rate of increase relative to K, no stochastic sampling
abline(h=K.max, lty=2, col="red") #carry capacity

#################################################### 
## iterations and quasi ext for each following model
####################################################
iter <- 10000 #final model run at 10 000
itdiv <- iter/1000 #final model rate at iter/1000

################################################################################################################
## untreated population
###############################################################################################################
## stochatic projection with density feedback
## set storage matrices & vectors

n.sums.mat <- matrix(data = 0, nrow = iter, ncol = (t+1)) #storage matrix

for (e in 1:iter) {
  popmat <- popmat.orig
  
  n.mat <- matrix(0, nrow=age.max,ncol=(t+1))
  n.mat[,1] <- init.vec
  
  for (i in 1:t) {
    # stochastic survival values
    s.alpha <- estBetaParams(s.vec, s.sd.vec^2)$alpha
    s.beta <- estBetaParams(s.vec, s.sd.vec^2)$beta
    s.stoch <- rbeta(length(s.alpha), s.alpha, s.beta)
    
    # stochastic fertilty sampler (gaussian)
    fert.stch <- rnorm(length(popmat[,1]), popmat[1,], m.sd.vec)
    fert.stoch <- ifelse(fert.stch < 0, 0, fert.stch)
    
    totN.i <- sum(n.mat[,i])
    pred.red <- a.lp/(1+(totN.i/b.lp)^c.lp)
    
    popmat[1,] <- fert.stoch
    diag(popmat[2:age.max,]) <- s.stoch*pred.red
    #popmat[age.max,age.max] <- 0
    
    n.mat[,i+1] <- popmat %*% n.mat[,i]
    
  } # end i loop
  
  n.sums.mat[e,] <- ((as.vector(colSums(n.mat))/pop.found))
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

n.md <- apply(n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
n.up <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
n.lo <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

plot(yrs,n.md,type="l", main = "Min N with SD for untr pop", xlab="year", ylab="Minimum population", lwd=2, ylim=c(0.95*min(n.lo),1.05*max(n.up)))
lines(yrs,n.lo,lty=2,col="red",lwd=1.5)
lines(yrs,n.up,lty=2,col="red",lwd=1.5)

untreated <- data.frame(yrs, n.md, n.lo, n.up)

###############################################################################################################################
## constant proportional yearly harvest
###############################################################################################################################

# harvest rate
harv.prop.consist <- seq(0.2,0.99,0.05) #sequence harvest/culling quotas, e.g remove 0.5-.99 porportion of founding pop

# define our quasi-extinction probability storage vector
min.med.n <- min.lo.n <- min.up.n <- rep(0,length(harv.prop.consist))

for (s in 1:length(harv.prop.consist)) {
  
  # set storage matrices & vectors
  n.sums.mat <- matrix(data = 0, nrow = iter, ncol = (t+1))
  
  for (e in 1:iter) {
    popmat <- popmat.orig
    
    n.mat <- matrix(0, nrow=age.max,ncol=(t+1))
    n.mat[,1] <- init.vec
    
    for (i in 1:t) {
      # stochastic survival values
      s.alpha <- estBetaParams(s.vec, s.sd.vec^2)$alpha
      s.beta <- estBetaParams(s.vec, s.sd.vec^2)$beta
      s.stoch <- rbeta(length(s.alpha), s.alpha, s.beta)
      
      # stochastic fertilty sampler (gaussian)
      fert.stch <- rnorm(length(popmat[,1]), popmat[1,], m.sd.vec)
      fert.stoch <- ifelse(fert.stch < 0, 0, fert.stch)
      
      totN.i <- sum(n.mat[,i])
      pred.red <- a.lp/(1+(totN.i/b.lp)^c.lp)
      
      popmat[1,] <- fert.stoch
      diag(popmat[2:age.max,]) <- s.stoch*pred.red
      #popmat[age.max,age.max] <- 0
      
      n.mat[,i+1] <- popmat %*% n.mat[,i]
      
      ## harvest things here
      n.mat[,i+1] <- n.mat[,i+1] - round(stable.stage.dist(popmat) * round(sum(n.mat[,i+1])*harv.prop.consist[s], 0), 0)
      
      
      if (length(which(n.mat[,i+1] < 0)) > 0) {
        n.mat[which(n.mat[,i+1] < 0), i+1] <- 0
      }
      
    } # end i loop
    
    n.sums.mat[e,] <- as.vector((colSums(n.mat))/pop.found) # / pop.mat for min proportion remaining population 
    
    if (e %% itdiv==0) print(e) 
    
  } # end e loop
  
  # calculate minimum population size
  
  min.pop.vec <- apply(n.sums.mat, MARGIN=1, min)
  min.med.n[s] <- median(min.pop.vec, na.rm=T)
  min.lo.n[s] <- quantile(min.pop.vec, probs=0.025, na.rm=T)
  min.up.n[s] <- quantile(min.pop.vec, probs=0.975, na.rm=T)
  
  n.md <- apply((n.sums.mat), MARGIN=2, mean, na.rm=T) # minimum over all iterations
  n.up <- apply((n.sums.mat), MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  n.lo <- apply((n.sums.mat), MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
  plot(yrs,n.md,type="l",xlab="year", ylab="minimum N", lwd=2, ylim=c(0.95*min(n.lo),1.05*max(n.up)))
  lines(yrs,n.lo,lty=2,col="red",lwd=1.5)
  lines(yrs,n.up,lty=2,col="red",lwd=1.5)
  
  print("##############")
  print(paste("harvest proportion = ", harv.prop.consist[s], sep=""))
  print("##############")
  
} # ends S loop

plot(harv.prop.consist, min.med.n, type="l", pch=19, xlab="harvest proportion", ylab="min N", ylim=c(min(min.lo.n),max(min.up.n)))
lines(harv.prop.consist, min.lo.n, col="red", lty=2)
lines(harv.prop.consist, min.up.n, col="red", lty=2)

minn.prop.pop <- data.frame(harv.prop.consist, min.med.n, min.lo.n, min.up.n)


##################################################################################################################################################
## high harvest for first 2 years, constant proportional harvest in remaining years
#####################################################################################################################################################

# harvest rate
harv.prop.init <- seq(0.5,0.9,0.05)
harv.prop.maint <- seq(0.1,0.5,0.05)

# storage
minn.med.mat <- minn.lo.mat <- minn.up.mat <- pmin.med.mat <- pmin.lo.mat <- pmin.up.mat <- matrix(data=NA, ncol=length(harv.prop.maint), nrow=length(harv.prop.init)) #storage matrices

for (m in 1:length(harv.prop.maint)) {
  
  for (n in 1:length(harv.prop.init)) {
    
    # storage
    n.sums.mat <- p.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
    
    for (e in 1:iter) {
      
      popmat <- popmat.orig
      
      n.mat <- matrix(0, nrow=age.max,ncol=(t+1))
      n.mat[,1] <- init.vec
      
      for (i in 1:t) {
        # stochastic survival values
        s.alpha <- estBetaParams(s.vec, s.sd.vec^2)$alpha
        s.beta <- estBetaParams(s.vec, s.sd.vec^2)$beta
        s.stoch <- rbeta(length(s.alpha), s.alpha, s.beta)
        
        # stochastic fertilty sampler (gaussian)
        fert.stch <- rnorm(length(popmat[,1]), popmat[1,], m.sd.vec)
        fert.stoch <- ifelse(fert.stch < 0, 0, fert.stch)
        
        totN.i <- sum(n.mat[,i])
        pred.red <- a.lp/(1+(totN.i/b.lp)^c.lp)
        
        popmat[1,] <- fert.stoch
        diag(popmat[2:age.max,]) <- s.stoch*pred.red
        #popmat[age.max,age.max] <- 0
        
        n.mat[,i+1] <- popmat %*% n.mat[,i]
        
        ## harvest 
        if (i < 3) {
          n.mat[,i+1] <- n.mat[,i+1] - round(stable.stage.dist(popmat) * round(sum(n.mat[,i+1])*harv.prop.init[n], 0), 0)
        } else {
          n.mat[,i+1] <- n.mat[,i+1] - round(stable.stage.dist(popmat) * round(sum(n.mat[,i+1])*harv.prop.maint[m], 0), 0)
        }
        
        if (length(which(n.mat[,i+1] < 0)) > 0) {
          n.mat[which(n.mat[,i+1] < 0), i+1] <- 0
        }
        
      } # end i loop
      
      n.sums.mat[e,] <- as.vector(colSums(n.mat))
      p.sums.mat[e,] <- n.sums.mat[e,] / pop.found
      
      if (e %% itdiv==0) print(e) 
    } # end e loop (stochastic iterations)
    
    min.pop.vec <- apply(n.sums.mat, MARGIN=1, min, na.rm=T)
    min.ppop.vec <- apply(p.sums.mat, MARGIN=1, min, na.rm=T)
    
    # median, lower & upper minimum population sizes
    minn.med.mat[n, m] <- median(min.pop.vec, na.rm=T) 
    minn.lo.mat[n, m] <- quantile(min.pop.vec, probs=0.025, na.rm=T)
    minn.up.mat[n, m] <- quantile(min.pop.vec, probs=0.975, na.rm=T)
    
    # median, lower & upper minimum proportional population sizes
    pmin.med.mat[n, m] <- median(min.ppop.vec, na.rm=T)
    pmin.lo.mat[n, m] <- quantile(min.ppop.vec, probs=0.025, na.rm=T) 
    pmin.up.mat[n, m] <- quantile(min.ppop.vec, probs=0.975, na.rm=T)
    
    
    print("##############################")
    print(paste("init harvest proportion = ", harv.prop.init[n], sep=""))
    print("##############################")
    
  } # end n loop (initial harvest rate)
  
  print("##############################")
  print(paste("maint harvest proportion = ", harv.prop.maint[m], sep=""))
  print("##############################")
  
} # end m loop (maintenance harvest rate)

## plot 3D surfaces
f1 <- list(
  family = "Avenir Light",
  size = 26,
  color = "black"
)
f2 <- list(
  family = "Avenir Light",
  size = 18,
  color = "black"
)
f3 <- list(
  family = "Avenir Light",
  size = 16,
  color = "black"
)

# minimum proportional population size (median)
par(mar=c(5,5,2,8))
pminmed3d <- plot_ly(z = ~pmin.med.mat, autocontour=F, type="contour", line = list(smoothing = 0.90), contours = list(start=0.01, end=0.32, size=0.025, showlabels = TRUE, labelfont=list(
  size=18, family="Avenir Light", face="bold", color="white"))) %>%
  colorbar(title = "med min pN1", titlefont=f2, tickfont=f2) %>%
  layout(
    xaxis = list(title="maintenance cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.1,0.5,0.1)), tickvals=seq(0,8,2)),
    yaxis = list(title="initial cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.5,0.9,0.1)), tickvals=seq(0,8,2)))
pminmed3d

twophase.med <- data.frame(pmin.med.mat)
colnames(twophase.med) <- harv.prop.maint
rownames(twophase.med) <- harv.prop.init

twophase.lo <- data.frame(pmin.lo.mat)
colnames(twophase.lo) <- harv.prop.maint
rownames(twophase.lo) <- harv.prop.init

twophase.up <- data.frame(pmin.up.mat)
colnames(twophase.up) <- harv.prop.maint
rownames(twophase.up) <- harv.prop.init

pmin3d <- plot_ly(showscale = FALSE) %>% 
  add_surface(z = ~pmin.med.mat) %>%
  add_surface(z = ~pmin.lo.mat, opacity = 0.55) %>%
  add_surface(z = ~pmin.up.mat, opacity = 0.55) %>%
  layout(scene = list(
    xaxis = list(title="maintenance cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.1,0.5,0.1)), tickvals=seq(0,8,2)),
    yaxis = list(title="initial cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.5,0.9,0.1)), tickvals=seq(0,8,2)),
    zaxis = list(title="min pN1", tickfont=f3, titlefont=f1)))
pmin3d

# quasi ext (median)
par(mar=c(5,5,2,8))
minmed3d <- plot_ly(z = ~qext.mat, autocontour=T, type="contour", line = list(smoothing = 0.90), contours = list(showlabels = TRUE, labelfont=list(
  size=18, family="Avenir Light", face="bold", color="white"))) %>%
  colorbar(title = "qE", titlefont=f2, tickfont=f2) %>%
  layout(
    xaxis = list(title="maintenance cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.1,0.5,0.1)), tickvals=seq(0,8,2)),
    yaxis = list(title="initial cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.5,0.9,0.1)), tickvals=seq(0,8,2)))
minmed3d

min3d <- plot_ly(showscale = FALSE) %>% 
  add_surface(z = ~minn.med.mat) %>%
  add_surface(z = ~minn.lo.mat, opacity = 0.55) %>%
  add_surface(z = ~minn.up.mat, opacity = 0.55) %>%
  layout(scene = list(
    xaxis = list(title="maintenance cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.1,0.5,0.1)), tickvals=seq(0,8,2)),
    yaxis = list(title="initial cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.5,0.9,0.1)), tickvals=seq(0,8,2)),
    zaxis = list(title="min N1", tickfont=f3, titlefont=f1)))
min3d


#######################################################################################################################################################################################################
######################################################################################### COSTS #################################################################################################
########################################################################################################################################################################################################
## high harvest for initial 2 years, consistent harvest remaining years
############################################################################################################################################################################################################
## contributed by CJA Bradshaw
###########################################################################################

Dudley.area <- 375*100 #ha
KI.area <- 4405*100  #ha

## cost parameters
felixer.unit <- 13000 # AU$ #cost from felix vs felixer report 
trap.unit <- c(157,297) # AU$ # cost per trap from traps.com.au
shoot.ph <- 518.54/20 # ammo & labour (total AU$ over 20 hours) # Holmes et al 2015
only.bait.unit <- (2.07 + 0.2) # From Curiosity correspondence. $2.07 per bait + $0.20  + $250 administration fee per order, + freight fee. 500 baits per pack

# Felixers data from Moseby et al (2020)
# 20 felixers for "felixer paddock", Arid Recovery (n1 = 48), killed 31 cats (n2 = 17), over 41 days 
num.felixer <- round((20/48) * pop.found, 0)
pfelixer.killr <- (31/48 * (1/(41/365)))


felixer.area <- 26*100 #ha; density from Arid Recovery trial, "felixer paddock" = 26 km^2
felixer.dens <- 20/felixer.area #20 felixer traps over the area, average density 0.77 felixers/km^2
KI.felixer.num <- round(KI.area * felixer.dens, 0) # number of felixers needed if same density was applied throughout Kangaroo Island
KI.felixer.num # not neccessarily reflective of the use of felixers as they are used in targeted areas and spread sporadically, as opposed to systematicaly placed like traps or baits

# traps
# 40 traps killed 21 over 148 days Hodgens 
ptrap.killr <- (21/262 * (1/(148/365)))
trap.dens <- 40/Dudley.area
KI.trap.num <- round(KI.area * trap.dens, 0)
KI.trap.num

# shooting
# 14725 person-hours killed 872 (+ 172 from wounds) cats (Marion) Parkes et al. 2014 & Bloomer & Bester 1992
# assume cats not killed by Felixers & traps shot by hunters
cats.pph <- (872+172)/14725


# baiting 
# 943 baits killed 11 cats over 18.86 km^2. Pre-baiting dens = 1.18 cats/km^2, post-baiting = 0.58 cats/km^2. Ref, 'Dudley peninsula feral cat eradication operations plan: summary may 2020 - mid 2023"
# KI uses Curiosity (PAPP)
## Kangaroo Island area - 4,405 km^2  or 440 500 ha 
# can't bait built-up areas, need 500m buffer zone around towns, built up areas 362 ha, how many built up areas? 
## parndana (second largest town) approx area as circle - 2km^2 (??), + 500m buffer = area 10km^2. 5 'main towns' KI. 5*2 = 10km^2 or 1000 ha, :- approx 1000 ha can't bait urban
# can't bait beaches. KI 540 km coastline, arbitary 100m buffer around coastline = 594 km can't bait + buffer zone. 540 * 1.1 = area no bait
## Dirk Hartog Island, 15 cats collared, average density 0.701 cats/km^2 (average area = 10.515 km^2 (A = 15*0.701), 50 baits per km, baits = 50*10.515 = 525.75), 14 died following bait consumption ... 525.74/14 = 37.55 baits/cat
# Dirk Hartog Island, used eradicat (1080)
nobaitfarm <- (2303 - (2303*.94))*100 #ha; can't be baited 
nobaitcoast <- (540 * 1.1)*100 #ha; dist around costline, *1.1 for the 100m buffer around coast 
nobaittown <- 1000 #ha; can't bait town 
nobaitarea <- nobaitfarm + nobaitcoast + nobaittown # total area can't be baited
baitareaKI <- KI.area - nobaitarea #ha; area eligible for baiting 
baitdens <- 50/100 # 50 baits per km^2 converted to baits per ha
baitnum <- (baitareaKI * baitdens) #number of baits need for entire Island

baitadminfee <- 250 # administration fee, once off for baits, or twice off for two years 
baitdrop.time <- (30/60/60) * (baitareaKI/100) # 50 baits drops every 30 seconds or 1 km^2, with plane speed 240km/h - bait area/100 to convert back to km^2
trips <- baitnum/3500 #can only take 3500 baits per trip 
upandback <- seq(1,32,0.5) #ha; dist from airport to start of each bait transect
averagedist <- (sum(upandback))/(length(upandback)) #ha; average dist from airport to transect
baitreloadtime <- ((averagedist*trips)*2)/240 #52 trips needed to drop all baits, *2 for to and from airport, plan speed 240km
planehph <- 750  #cost per hour plane hire when dropping baits, inc wages of 2x pilots (1x loading and dropping baits)
planeferrycost <- (3*600)*2 # 3 hour flight William Creek to Kangaroo Island (*2 for return), $600 p/h for charter
baitreloadcost <- baitreloadtime*planehph #average extra cost for reloading baits
baitdropcost <- baitdrop.time*planehph
planecost <- planeferrycost + baitreloadcost +  baitdropcost
cost.total.bait <- planecost + baitadminfee + (only.bait.unit * baitnum) #cost of total baiting to cover entire island 

  
pbait.killr <- 14/525.74 # 14 cats killed by 525.74 baits



###########################################################################################

###########################################################################################
## Type III functional response (reduction in capture efficiency with decreasing density)
max.eff <- 1 # max efficiency 
min.eff <- 0 #min efficiency
max.pN <- 1 #max population proportion
min.pN <- 0 #min population proportion
infl.eff <- 0.5

pN.vec <- c(min.pN, 0.2, 0.4, 0.5, 0.7, 0.8, max.pN) 
eff.vec <- c(min.eff, 0.05, 0.3, infl.eff, 0.85, 0.95, max.eff)
plot(pN.vec, eff.vec, type="b", pch=19)
eff.dat <- data.frame(pN.vec, eff.vec)
colnames(eff.dat) <- c("pN", "eff")

# a/(1 + b*e^(-cx)) (logistic)
param.init <- c(1, 85, 8.9)
fit.eff <- nls(eff ~ a/(1+(b*exp(-c*pN))), 
               data = eff.dat,
               algorithm = "port",
               start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
               trace = TRUE,      
               nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
fit.eff.summ <- summary(fit.eff)
plot(pN.vec,eff.vec,pch=19,xlab="pN",ylab="efficiency")
pN.vec.cont <- seq(0,1,0.01)
pred.eff.fx <- coef(fit.eff)[1]/(1+(coef(fit.eff)[2]*exp(-coef(fit.eff)[3]*pN.vec.cont)))
lines(pN.vec.cont,pred.eff.fx,lty=2,lwd=3,col="red")

a.eff <- coef(fit.eff)[1]
b.eff <- coef(fit.eff)[2]
c.eff <- coef(fit.eff)[3]

###########################################################################################


#################################################### 
## iterations and quasi ext for each following model
####################################################
iter <- 10000 #final model run at 10 000
itdiv <- iter/1000 #final model rate at iter/1000

## run choices
## make up shortfall in kill by ...
#shortfall.method <- "F" # adding Felixer units
shortfall.method <- "T" # adding traps
#shortfall.method <- "H" # increasing hunting pressure

# harvest rate
harv.prop.init <- seq(0.5,0.9,0.05)
harv.prop.maint <- seq(0.1,0.5,0.05)
q.ext <- 20

# storage
qext.mat <- minn.med.mat <- minn.lo.mat <- minn.up.mat <- pmin.med.mat <- pmin.lo.mat <- pmin.up.mat <- totcost.med <- totcost.lo <- totcost.up <- matrix(data=NA, ncol=length(harv.prop.maint), nrow=length(harv.prop.init))

for (m in 1:length(harv.prop.maint)) {
  
  for (n in 1:length(harv.prop.init)) {
    
    # storage
   init.k.sums.mat <- k.sums.mat <- n.sums.mat <- p.sums.mat <- totalcost.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
    
    for (e in 1:iter) {
      
      popmat <- popmat.orig
      
      init.k.mat <- n.mat <- k.mat <- matrix(0, nrow=age.max,ncol=(t+1))
      n.mat[,1] <- init.vec
      
      for (i in 1:t) {
        # stochastic survival values
        s.alpha <- estBetaParams(s.vec, s.sd.vec^2)$alpha
        s.beta <- estBetaParams(s.vec, s.sd.vec^2)$beta
        s.stoch <- rbeta(length(s.alpha), s.alpha, s.beta)
        
        # stochastic fertilty sampler (gaussian)
        fert.stch <- rnorm(length(popmat[,1]), popmat[1,], m.sd.vec)
        fert.stoch <- ifelse(fert.stch < 0, 0, fert.stch)
        
        totN.i <- sum(n.mat[,i])
        pred.red <- a.lp/(1+(totN.i/b.lp)^c.lp)
        
        popmat[1,] <- fert.stoch  # add new stochastically resampled fertilities
        diag(popmat[2:age.max,]) <- s.stoch*pred.red
        #popmat[age.max,age.max] <- s.stoch[age.max]*pred.red
        
        n.mat[,i+1] <- popmat %*% n.mat[,i]
        
        ## harvest things here
        if (i < 3) {
          k.mat[,i+1] <- round(stable.stage.dist(popmat) * round(sum(n.mat[,i+1])*harv.prop.init[n], 0), 0)
          n.mat[,i+1] <- n.mat[,i+1] - k.mat[,i+1]
          init.k.mat[,i+1] <- n.mat[,i+1]
        } else {
          k.mat[,i+1] <- round(stable.stage.dist(popmat) * round(sum(n.mat[,i+1])*harv.prop.maint[m], 0), 0)
          n.mat[,i+1] <- n.mat[,i+1] - k.mat[,i+1]
        }
        
        if (length(which(n.mat[,i+1] < 0)) > 0) {
          n.mat[which(n.mat[,i+1] < 0), i+1] <- 0
        }
        if (length(which(k.mat[,i+1] < 0)) > 0) {
          k.mat[which(k.mat[,i+1] < 0), i+1] <- 0
        }
        
      } # end i loop
      
      init.k.sums.mat[e,] <- as.vector(colSums(init.k.mat))
      k.sums.mat[e,] <- as.vector(colSums(k.mat))
      n.sums.mat[e,] <- as.vector(colSums(n.mat))
      p.sums.mat[e,] <- n.sums.mat[e,] / pop.found
      
      # cost of cats killed here
      eff.vec.iter <- a.eff/(1+(b.eff*exp(-c.eff*p.sums.mat[e,]))) # efficiency this iteration
      
      # calculate numbers killed per year using baiting and trapping first two years
      bait.kill.base <- round(init.k.sums.mat[e,] * (eff.vec.iter*pbait.killr), 0)
      trap.kill.base <- round(k.sums.mat[e,] * (eff.vec.iter*ptrap.killr), 0)
      bt.kill.base <- trap.kill.base + bait.kill.base
      shortfall <- k.sums.mat[e,] - bt.kill.base # how many cats not being killed by these methods?
      
      #base cost
      base.cost <- (cost.total.bait*2) + (KI.trap.num*runif(1,min=trap.unit[1],max=trap.unit[2])) # at initial roll-out numbers
      
      # make up shortfall
      if (shortfall.method == "H") {
        makeup.iter <- shoot.ph*(shortfall / (cats.pph*eff.vec.iter)) # how many person-hours required to make up shortfall?
      }
      if (shortfall.method == "F") {
        makeup.iter <- felixer.unit*(shortfall / (pfelixer.killr*eff.vec.iter)) # how many person-hours required to make up shortfall?
      }
      if (shortfall.method == "T") {
        makeup.iter <- (runif(1,min=trap.unit[1],max=trap.unit[2]))*(shortfall / (ptrap.killr*eff.vec.iter)) # how many person-hours required to make up shortfall?
      }
      
      totalcost.mat[e,] <- base.cost + makeup.iter 
      
      if (e %% itdiv==0) print(e) 
    } # end e loop (stochastic iterations)
    
    min.pop.vec <- apply(n.sums.mat, MARGIN=1, min, na.rm=T)
    min.ppop.vec <- apply(p.sums.mat, MARGIN=1, min, na.rm=T)
    
    # median, lower & upper minimum population sizes
    minn.med.mat[n, m] <- median(min.pop.vec, na.rm=T) 
    minn.lo.mat[n, m] <- quantile(min.pop.vec, probs=0.025, na.rm=T)
    minn.up.mat[n, m] <- quantile(min.pop.vec, probs=0.975, na.rm=T)
    
    # median, lower & upper minimum proportional population sizes
    pmin.med.mat[n, m] <- median(min.ppop.vec, na.rm=T)  
    pmin.lo.mat[n, m] <- quantile(min.ppop.vec, probs=0.025, na.rm=T) 
    pmin.up.mat[n, m] <- quantile(min.ppop.vec, probs=0.975, na.rm=T)
    
    # quasi-extinction
    qext.mat[n, m] <- (sum(ifelse(round(min.pop.vec, 0) < q.ext, 1, 0)) / iter)
    
    ## costs
    totcost.vec <- apply(totalcost.mat, MARGIN=1, sum, na.rm=T)
    totcost.med[n, m] <- median(totcost.vec, na.rm=T)
    colnames(totcost.med) <- harv.prop.maint
    rownames(totcost.med) <- harv.prop.init
    
    totcost.lo[n, m] <- quantile(totcost.vec, probs=0.025, na.rm=T)
    colnames(totcost.lo) <- harv.prop.maint
    rownames(totcost.lo) <- harv.prop.init
    
    totcost.up[n, m] <- quantile(totcost.vec, probs=0.975, na.rm=T)
    colnames(totcost.up) <- harv.prop.maint
    rownames(totcost.up) <- harv.prop.init
    
    print("##############################")
    print(paste("init harvest proportion = ", harv.prop.init[n], sep=""))
    print("##############################")
    
  } # end n loop (initial harvest rate)
  
  print("##############################")
  print(paste("maint harvest proportion = ", harv.prop.maint[m], sep=""))
  print("##############################")
  
} # end m loop (maintenance harvest rate)

## plot 3D surfaces
f1 <- list(
  family = "Avenir Light",
  size = 26,
  color = "black"
)
f2 <- list(
  family = "Avenir Light",
  size = 18,
  color = "black"
)
f3 <- list(
  family = "Avenir Light",
  size = 16,
  color = "black"
)

# total cost (median)
par(mar=c(5,5,2,8))
costcontmed3d <- plot_ly(z = ~totcost.med, autocontour=T, type="contour", line = list(smoothing = 0.90), contours = list(showlabels = TRUE, labelfont=list(
  size=18, family="Avenir Light", face="bold", color="white"))) %>%
  colorbar(title = "tot $", titlefont=f2, tickfont=f2) %>%
  layout(
    xaxis = list(title="maintenance cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.1,0.5,0.1)), tickvals=seq(0,8,2)),
    yaxis = list(title="initial cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.5,0.9,0.1)), tickvals=seq(0,8,2)))
costcontmed3d

cost3d <- plot_ly(showscale = FALSE) %>% 
  add_surface(z = ~totcost.med) %>%
  add_surface(z = ~totcost.lo, opacity = 0.55) %>%
  add_surface(z = ~totcost.up, opacity = 0.55) %>%
  layout(scene = list(
    xaxis = list(title="maintenance cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.1,0.5,0.1)), tickvals=seq(0,8,2)),
    yaxis = list(title="initial cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.5,0.9,0.1)), tickvals=seq(0,8,2)),
    zaxis = list(title="tot $", tickfont=f3, titlefont=f1)))
cost3d


################################################################################################################################################################################################
#################################################################################################### trap-neuter-release #######################################################################
#################################################################################################################################################################################################
### SAME METHODS AS ABOVE, ALTERED FERTILITY INSTEAD OF SURVIVAL

TNR <- seq(.01,.9,.01)

min.med.n <- min.lo.n <- min.up.n <- rep(0,length(TNR))

for (s in 1:length(TNR)) {
  
  #storage matrix
  n.sums.mat <- matrix(0, nrow = iter, ncol = (t+1))
  
  for (e in 1:iter){
    popmat <- popmat.orig
    
    n.mat <- matrix(0, nrow = age.max, ncol = (t+1))
    n.mat[,1] <- init.vec
    
    for (i in 1:t) {
      # stochastic survival values
      s.alpha <- estBetaParams(s.vec, s.sd.vec^2)$alpha
      s.beta <- estBetaParams(s.vec, s.sd.vec^2)$beta
      s.stoch <- rbeta(length(s.alpha), s.alpha, s.beta)
      
      # stochastic fertilty sampler (gaussian)
      fert.stch <- rnorm(length(popmat[,1]), popmat[1,], m.sd.vec)
      fert.stoch <- ifelse(fert.stch < 0, 0, fert.stch)
      
      totN.i <- sum(n.mat[,i])
      pred.red <- a.lp/(1+(totN.i/b.lp)^c.lp)
      
      popmat[1,] <- fert.stoch  # add new stochastically resampled fertilities
      diag(popmat[2:age.max,]) <- s.stoch*pred.red # add new stochastically resampled survivals
      #popmat[age.max,age.max] <- 0 # add new stochastically resampled survivals
      
      #fertility reduction 
      popmat[1,] <- popmat[1,]*TNR[s]
      
      # project
      n.mat[,i+1] <- popmat %*% n.mat[,i]
      
      if (length(which(n.mat[,i+1] < 0)) > 0) {
        n.mat[which(n.mat[,i+1] < 0), i+1] <- 0
      }
      
    } #end i loop
    
    n.sums.mat[e,] <- as.vector(colSums(n.mat))
    
    if (e %% itdiv==0) print(e) 
    
  } #end e loop
  
  min.pop.vec <- apply(n.sums.mat, MARGIN=1, min, na.rm=T)
  min.med.n[s] <- median(min.pop.vec, na.rm=T)
  min.lo.n[s] <- quantile(min.pop.vec, probs=0.025, na.rm=T)
  min.up.n[s] <- quantile(min.pop.vec, probs=0.975, na.rm=T)
  
  n.md <- apply(n.sums.mat, MARGIN=2, mean, na.rm=T) # minimum over all iterations
  n.up <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  n.lo <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
  plot(yrs,n.md,type="l",xlab="year", ylab="minimum N", lwd=2, ylim=c(0.95*min(n.lo),1.05*max(n.up)))
  lines(yrs,n.lo,lty=2,col="red",lwd=1.5)
  lines(yrs,n.up,lty=2,col="red",lwd=1.5)
  
  print("##############")
  print(paste("TNR = ", TNR[s], sep=""))
  print("##############")
  
} #end s loop

plot(1-TNR, min.med.n/pop.found, type="l", pch=19, xlab="proportion spayed each year", ylab="proportion of N1", ylim=c(min(min.lo.n/pop.found),max(min.up.n/pop.found)))
lines(1-TNR, min.lo.n/pop.found, col="red", lty=2)
lines(1-TNR, min.up.n/pop.found, col="red", lty=2)

spay.out <- data.frame(1-TNR, min.med.n/pop.found, min.up.n/pop.found, min.lo.n/pop.found)
colnames(spay.out) <- c("pSpay","pNmed","pNup","pNlo")


TNR.pop <- data.frame(1-TNR, min.med.n, min.up.n, min.up.n)
colnames(TNR.pop) <- c('proportion spayed', 'Median', 'Upper', 'Lower')
