source('RPEDR.R')
source("RPEDR_dim_est.R")
source("RPEDR_double.R")

### single simulation test ###
n <- 200
#N <- n
p <- 20
L <- 200
M <- 10*p
d <- ceiling(sqrt(p))

X <- matrix(rnorm(n*p),n,p)
## Model 1a ##
fx <- (X[,1] + X[,2])^2
d0 <- 1
V0 <- matrix(0,p,d0)
V0[,1] <- c(rep(1/sqrt(2),2), rep(0,p-2))

x = X
y = fx + rnorm(n)/2


# d0 prespecified, dimension reduction directions provided by RPEDR
m.RPEDR <- RPEDR(X=x, Y=y, d=d, L=L, M=M, version = "mars", dist="mixed", n.cores=10)
m.RPEDR.dir <- m.RPEDR$U[,seq(d0)]

# No preknowledge of d0, the dimension estimated by Algorithm 2
est.dim <- RPEDR_dim_est(m.RPEDR, L=L, R=10000, dist="mixed", n.cores=10)

# Dimension reduction directions provided by double RPEDR
m.RPEDR.double.dir <- RPEDR_double(X=x, Y=y, d=d, L=L, M=M, d0=d0, version = "mars", dist="mixed", n.cores=10)
