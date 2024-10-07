## Loading Packages ##
library(pracma)
library(primePCA)
library(earth)
library(doParallel)
library(foreach)
library(VariableScreening)
library(dr)
library(MAVE)
library(readxl)
source('RPEDR.R')
source("RPEDR_dim_est.R")
source("RPEDR_double.R")


music.data <- read.table("Real_Data/default_plus_chromatic_features_1059_tracks.txt", sep = ",")

## setting d0 ##
d0 <- 3
set.seed(200)

n <- nrow(music.data)
s <- sample(n, 2*n/3)
p <- ncol(music.data) - 2

train.x <- as.matrix(subset(music.data[s,], select = -c(`V117`, `V118`)))
train.y <- music.data[s,]$V117

test.x <- as.matrix(subset(music.data[-s,], select = -c(`V117`, `V118`)))
test.y <- music.data[-s,]$V117


## MARS ##
mars.train.data <- cbind.data.frame(train.y, train.x)
names(mars.train.data)[1] <- "Y"

smp=seq(min(p,4))
gcv.mars=Inf
for (df in smp) {
  fiti=earth(Y ~., data=mars.train.data, degree=df)
  if(fiti$gcv<gcv.mars){
    fit.mars=fiti
    gcv.mars=fiti$gcv
    df=df
  }
}

mars.est.reg <- predict(fit.mars, test.x)
rmse.mars <- sqrt(mean((mars.est.reg - test.y)^2))
cat("The RMSE of MARS is", rmse.mars, '\n')


## SIR ##
m.sir <- dr(train.y ~ train.x, method="sir")
m.sir.dir <- m.sir$evectors[,1:d0]
v.n <- as.numeric(gsub("[^0-9]", "", rownames(m.sir.dir)))
m.sir.dir <- pracma::gramSchmidt(as.matrix(m.sir.dir))$Q


Z <- train.x[,v.n] %*% m.sir.dir
proj.data <- cbind.data.frame(train.y, Z)
names(proj.data)[1] <- "Y"
smp=seq(min(p,4))
gcv.sir=Inf
for (df in smp) {
  fiti=earth(Y ~., data=proj.data, degree=df)
  if(fiti$gcv<gcv.sir){
    fit.sir=fiti
    gcv.sir=fiti$gcv
    df=df
  }
}
sir.est.reg = predict(fit.sir, test.x[,v.n] %*% m.sir.dir)
rmse.sir <- sqrt(mean((sir.est.reg - test.y)^2))
cat("The RMSE of SIR is", rmse.sir, '\n')



## pHd ##
m.phd <- dr(train.y ~ train.x, method="phdres")
m.phd.dir <- m.phd$evectors[,1:d0]
v.n <- as.numeric(gsub("[^0-9]", "", rownames(m.phd.dir)))
m.phd.dir <- pracma::gramSchmidt(as.matrix(m.phd.dir))$Q

Z <- train.x[,v.n] %*% m.phd.dir
proj.data <- cbind.data.frame(train.y, Z)
names(proj.data)[1] <- "Y"
smp=seq(min(p,4))
gcv.phd=Inf
for (df in smp) {
  fiti=earth(Y ~., data=proj.data, degree=df)
  if(fiti$gcv<gcv.phd){
    fit.phd=fiti
    gcv.phd=fiti$gcv
    df=df
  }
}
phd.est.reg = predict(fit.phd, test.x[,v.n] %*% m.phd.dir)
rmse.phd <- sqrt(mean((phd.est.reg - test.y)^2))
cat("The RMSE of pHd is", rmse.phd, '\n')



## MAVE ##
dr.mave <- mave(train.y ~ train.x, method="meanMAVE")
mave.est.reg <- predict(dr.mave, test.x, dim=d0)
rmse.mave <- sqrt(mean((mave.est.reg - test.y)^2))
cat("The RMSE of MAVE is", rmse.mave, '\n')


### For reproducing results for gKDR and drMARS, see https://github.com/liuyu-star/drMARS ###


## RPEDR ##

L <- 200
M <- 10*p
d <- ceiling(sqrt(p))
m.RPEDR <- RPEDR(X=train.x, Y=train.y, 
                 d=d, L=L, M=M, version = "mars", dist="mixed", n.cores=20)
RP.dir <- m.RPEDR$U[,1:d0]

Z.RP <- train.x %*% RP.dir
proj.data <- cbind.data.frame(train.y, Z.RP)
names(proj.data)[1] <- "Y"
smp=seq(min(p,4))
gcv.RP=Inf
for (df in smp) {
  fiti=earth(Y ~., data=proj.data, degree=df)
  if(fiti$gcv<gcv.RP){
    fit.RP=fiti
    gcv.RP=fiti$gcv
    df.RP=df
  }
}
RP.est.reg = predict(fit.RP, test.x %*% RP.dir)
rmse.rpe <- sqrt(mean((RP.est.reg - test.y)^2))
cat("The RMSE of RPE is", rmse.rpe, '\n')



### RPEDR double ###
m.RPEDR.double.dir <- RPEDR_double(X=train.x, Y=train.y, 
                                   d=d, L=L, M=M, d0=d0, version = "mars", dist="mixed", n.cores=20)
Z.RP <- train.x %*% m.RPEDR.double.dir
proj.data <- cbind.data.frame(train.y, Z.RP)
names(proj.data)[1] <- "Y"
smp=seq(min(p,4))
gcv.RP2=Inf
for (df in smp) {
  fiti=earth(Y ~., data=proj.data, degree=df)
  if(fiti$gcv<gcv.RP2){
    fit.RP2=fiti
    gcv.RP2=fiti$gcv
    df.RP=df
  }
}
RP.d.est.reg = predict(fit.RP2, test.x %*% m.RPEDR.double.dir)
rmse.rpe2 <- sqrt(mean((RP.d.est.reg - test.y)^2))
cat("The RMSE of RPE2 is", rmse.rpe2, '\n')