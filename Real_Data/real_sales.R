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

residential.data <- read_excel("Real_Data/Residential-Building-Data-Set.xlsx", range = "F2:DE374")
n <- nrow(residential.data)
p <- ncol(residential.data) - 2

set.seed(101)
s <- sample(n, n/3)
d0 <- 10 # Setting d0

sales.data <- subset(residential.data, select=-`V-10`)
names(sales.data) <- c(paste0("X", 1:p), "Y")
sales.data[, -(p+1)] <- scale(sales.data[, -(p+1)])
sales.train.X <- as.matrix(subset(sales.data[-s,], select=-`Y`))
sales.train.Y <- sales.data[-s,]$Y
sales.test.X <- as.matrix(subset(sales.data[s,], select=-`Y`))
sales.test.Y <- sales.data[s,]$Y


## MARS ##
mars.train.data <- cbind.data.frame(sales.train.Y, sales.train.X)
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

mars.est.reg <- predict(fit.mars, sales.test.X)
rmse.mars <- sqrt(mean((mars.est.reg - sales.test.Y)^2))
cat("The RMSE of MARS is", rmse.mars, '\n')




## SIR ##
s1 <- Sys.time()
m.sir <- dr(sales.train.Y ~ sales.train.X, method="sir")
m.sir.dir <- m.sir$evectors[,1:d0]
v.n <- as.numeric(gsub("[^0-9]", "", rownames(m.sir.dir)))
m.sir.dir <- pracma::gramSchmidt(as.matrix(m.sir.dir))$Q


Z <- sales.train.X[,v.n] %*% m.sir.dir
proj.data <- cbind.data.frame(sales.train.Y, Z)
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
sir.est.reg = predict(fit.sir, sales.test.X[,v.n] %*% m.sir.dir)
rmse.sir <- sqrt(mean((sir.est.reg - sales.test.Y)^2))
cat("The RMSE of SIR is", rmse.sir, '\n')



## pHd ##
m.phd <- dr(sales.train.Y ~ sales.train.X, method="phdres")
m.phd.dir <- m.phd$evectors[,1:d0]
v.n <- as.numeric(gsub("[^0-9]", "", rownames(m.phd.dir)))
m.phd.dir <- pracma::gramSchmidt(as.matrix(m.phd.dir))$Q

Z <- sales.train.X[,v.n] %*% m.phd.dir
proj.data <- cbind.data.frame(sales.train.Y, Z)
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
phd.est.reg = predict(fit.phd, sales.test.X[,v.n] %*% m.phd.dir)
rmse.phd <- sqrt(mean((phd.est.reg - sales.test.Y)^2))
cat("The RMSE of pHd is", rmse.phd, '\n')



## MAVE ##
dr.mave <- mave(sales.train.Y ~ sales.train.X, method="meanMAVE")
mave.est.reg <- predict(dr.mave, sales.test.X, dim=d0)
rmse.mave <- sqrt(mean((mave.est.reg - sales.test.Y)^2))
cat("The RMSE of MAVE is", rmse.mave, '\n')


### For reproducing results for gKDR and drMARS, see https://github.com/liuyu-star/drMARS ###


## RPEDR ##

L <- 200
M <- 10*p
d <- ceiling(sqrt(p))
m.RPEDR <- RPEDR(X=construction.train.X, Y=construction.train.Y, 
                 d=d, L=L, M=M, version = "mars", dist="mixed", n.cores=20)
RP.dir <- m.RPEDR$U[,1:d0]

Z.RP <- sales.train.X %*% RP.dir
proj.data <- cbind.data.frame(sales.train.Y, Z.RP)
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
RP.est.reg = predict(fit.RP, sales.test.X %*% RP.dir)
rmse.rpe <- sqrt(mean((RP.est.reg - sales.test.Y)^2))
cat("The RMSE of RPE is", rmse.rpe, '\n')




### RPEDR double ###
m.RPEDR.double.dir <- RPEDR_double(X=sales.train.X, Y=sales.train.Y, 
                                   d=d, L=L, M=M, d0=d0, version = "mars", dist="mixed", n.cores=20)
Z.RP <- sales.train.X %*% m.RPEDR.double.dir
proj.data <- cbind.data.frame(sales.train.Y, Z.RP)
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
RP.d.est.reg = predict(fit.RP2, sales.test.X %*% m.RPEDR.double.dir)
rmse.rpe2 <- sqrt(mean((RP.d.est.reg - sales.test.Y)^2))
cat("The RMSE of RPE2 is", rmse.rpe2, '\n')