## Required packages ##
library(earth)
library(doParallel)
library(foreach)
library(regpro)

### RPEDR_dim_est Function ###
RPEDR_dim_est <- function(model_RPEDR, L, R = 10000, ortho = FALSE, 
                         dist = "mixed", ..., n.cores) {
  
  ### Helper function for random matrix generation ###
  generate_random_matrix <- function(p, d, L, dist, ortho, ...) {
    if (ortho) {
      A <- matrix(, nrow = p, ncol = 0)
      for (l in 1:L) {
        temp <- matrix(rt(p*p, 1), p, p)
        temp.Q <- qr.Q(qr(temp))
        A <- cbind(A, temp.Q[, 1:d])
      }
    } else {
      if (dist == "t") {
        A0 <- matrix(rt(L*d*p, ...), p, L*d)
      } else if (dist == "norm") {
        A0 <- matrix(rnorm(L*d*p, ...), p, L*d)
      } else if (dist == "cauchy") {
        A0 <- matrix(rcauchy(L*d*p, ...), p, L*d)
      } else if (dist == "mixed") {
        A0 <- matrix(0, p, L*d)
        unif.n <- runif(L)
        norm.ind <- unif.n > 0.5
        A0.norm <- matrix(rnorm(sum(norm.ind) * d * p), p, round(sum(norm.ind) * d))
        A0.cauchy <- matrix(rt(sum(!norm.ind) * d * p, 1), p, round(sum(!norm.ind) * d))
        idx.col.norm <- sort(as.vector(outer(which(norm.ind==1), 1:d, FUN = function(i, j) (i - 1) * d + j)))
        idx.col.cauchy <- setdiff(1:(L * d), idx.col.norm)
        A0[, idx.col.norm] <- A0.norm
        A0[, idx.col.cauchy] <- A0.cauchy
      }
      A <- t(t(A0) / sqrt(diag(t(A0) %*% A0)))
    }
    return(A)
  }
  
  p <- length(model_RPEDR$D)
  d <- ceiling(sqrt(p))
  
  n.cores <- min(n.cores, parallel::detectCores() - 2)
  cl <- makeCluster(n.cores)
  registerDoParallel(cl)
  
  D.list <- foreach(r=1:R) %dopar% { 
    A <- generate_random_matrix(p, d, L, dist, ortho, ...)
    proj <- matrix(0, p, p)
    for (l in seq_len(L)) {
      idx <- (l-1)*d + 1:d
      proj <-  proj +  1 / L * A[,idx] %*% t(A[,idx]) 
    }
    D.out <- svd(proj)$d   
    return(D.out)
  }
  
  stopCluster(cl)
  
  D.matrix <- matrix(, nrow=R, ncol=p)
  for (r in seq_len(R)) {
    D.matrix[r,] <- D.list[[r]]
  }
  
  D.median <- apply(D.matrix, 2, quantile, probs=0.5)
  
  est.d <- which.max(model_RPEDR$D<D.median)-1
  
  return(est.d)
}