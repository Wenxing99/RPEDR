## Required packages ##
library(earth)
library(doParallel)
library(foreach)
library(regpro)
source("RPEDR.R")
source("RPEDR_dim_est.R")

### RPEDR Function ###
RPEDR_double <- function(X, Y, d, L, M, d0, ortho = FALSE, 
                  dist = "mixed", ..., version = "mars", 
                  kernel = "normal", n.cores) {
  m.RPEDR <- RPEDR(X=X, Y=Y, d=d, L=L, M=M, version = "mars", dist="mixed", n.cores=n.cores, ...)
  est.dim <- RPEDR_dim_est(m.RPEDR, L=L, dist="mixed", n.cores=n.cores, ...)
  
  if (est.dim > d0) {
    X.double <- X%*%m.RPEDR$U[,1:est.dim]
    p.double <- est.dim
    M.double <- p.double*10
    d.double <- ceiling(sqrt(p.double))
    
    m.RPEDR.d <- RPEDR(X=X.double, Y=Y, d=d.double, L=L, M=M.double, version = "mars", dist="mixed", n.cores=n.cores, ...)
    m.RPEDR.d.dir <- (m.RPEDR$U[,1:est.dim]%*%m.RPEDR.d$U)[,seq(d0)]
  }else{
    m.RPEDR.d.dir <- m.RPEDR$U[,seq(d0)]
  }
  
  return(m.RPEDR.d.dir)
}