source('RPEDR.R')

### single simulation test ###
n <- 200
N <- n
p <- 10
L <- 200
M <- 10*p
d <- ceiling(sqrt(p))
D <- ceiling(sqrt(p))

X <- matrix(rnorm((n+N)*p),n+N,p)
## Model 1a ##
fx <- (X[,1] + X[,2])^2
d0 <- 1
V0 <- matrix(0,p,d0)
V0[,1] <- c(rep(1/sqrt(2),2), rep(0,p-2))

x=X[seq(n),]
snr=2
y = fx[seq(n)] + rnorm(n)/snr

m.RPEDR <- RPEDR(X=x, Y=y, d=d, L=L, M=M, D=D, version = "mars", dist="mixed", n.cores=10)
m.RPEDR.dir <- m.RPEDR$U[,seq(d0)]