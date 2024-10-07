## Required packages ##
library(earth)
library(doParallel)
library(foreach)
library(regpro)

### RPEDR Function ###
RPEDR <- function(X, Y, d, L, M, ortho = FALSE, dist = "mixed", ..., version = "mars", 
                  kernel = "normal", n.cores) {
  
  ### Helper function for random matrix generation ###
  generate_random_matrix <- function(p, d, M, dist, ortho, ...) {
    if (ortho) {
      A <- matrix(, nrow = p, ncol = 0)
      for (m in 1:M) {
        temp <- matrix(rt(p*p, 1), p, p)
        temp.Q <- qr.Q(qr(temp))
        A <- cbind(A, temp.Q[, 1:d])
      }
    } else {
      if (dist == "t") {
        A0 <- matrix(rt(M*d*p, ...), p, M*d)
      } else if (dist == "norm") {
        A0 <- matrix(rnorm(M*d*p, ...), p, M*d)
      } else if (dist == "cauchy") {
        A0 <- matrix(rcauchy(M*d*p, ...), p, M*d)
      } else if (dist == "mixed") {
        A0 <- matrix(0, p, M*d)
        unif.n <- runif(M)
        norm.ind <- unif.n > 0.5
        A0.norm <- matrix(rnorm(sum(norm.ind) * d * p), p, round(sum(norm.ind) * d))
        A0.cauchy <- matrix(rt(sum(!norm.ind) * d * p, 1), p, round(sum(!norm.ind) * d))
        idx.col.norm <- sort(as.vector(outer(which(norm.ind==1), 1:d, FUN = function(i, j) (i - 1) * d + j)))
        idx.col.cauchy <- setdiff(1:(M * d), idx.col.norm)
        A0[, idx.col.norm] <- A0.norm
        A0[, idx.col.cauchy] <- A0.cauchy
      }
      A <- t(t(A0) / sqrt(diag(t(A0) %*% A0)))
    }
    return(A)
  }
  
  ### Helper function for selecting best projection within group and return its outer product###
  best_proj_outer <- function(X, Y, A, s, idx, version, kernel) {
    Z <- X %*% A
    Out <- rep(0, length(idx))
    
    for (m in 1:length(idx)) {
      current_idx <- idx[[m]]
      Z_proj <- Z[, current_idx]
      
      if (version == "kernel") {
        fitted <- if (kernel == "normal") {
          apply(Z_proj[s, ], 1, function(z) {
            regpro::kernesti.regr(z, Z_proj[-s, ], Y[-s], h = 0.1, kernel = "gauss")
          })
        } else if (kernel == "loclin-normal") {
          apply(Z_proj[s, ], 1, function(z) {
            regpro::loclin(z, Z_proj[-s, ], Y[-s], h = 1, kernel = "gauss")
          })
        } else if (kernel == "loclin-uniform") {
          apply(Z_proj[s, ], 1, function(z) {
            regpro::loclin(z, Z_proj[-s, ], Y[-s], h = 0.1, kernel = "uniform")
          })
        } else {
          kreg.estimator <- kreg(x = Z_proj[-s, ], y = Y[-s], grid = Z_proj[s, ], bandwidth = 0.01, kernel = kernel)
          kreg.estimator$y
        }
      } else if (version == "poly") {
        train.data <- as.data.frame(cbind(Y = Y[-s], Z_proj[-s, ]))
        colnames(train.data)[-1] <- paste0("X", 1:(dim(train.data)[2]-1))
        formula <- as.formula(paste0("Y ~ .^2 + ", 
                                     paste0("I(", colnames(train.data)[-1], "^2)", collapse = "+")))
        fit <- lm(formula, data=train.data)
        newdata <- as.data.frame(Z_proj[s, ])
        colnames(newdata) <- paste0("X", 1:(dim(train.data)[2]-1))
        fitted <- predict(fit, newdata)
      } else if (version == "linear") {
        fit <- lm(Y ~ ., data = cbind(Y = Y[-s], as.data.frame(Z_proj[-s, ])))
        fitted <- predict(fit, as.data.frame(Z_proj[s, ]))
      } else if (version == "mars") {
        fit <- earth::earth(Y ~ ., data = cbind(Y = Y[-s], as.data.frame(Z_proj[-s, ])), degree = 3)
        fitted <- predict(fit, Z_proj[s, ])
      }
      
      Out[m] <- mean((fitted - Y[s])^2)
    }
    
    min_idx <- which.min(Out)
    return(1/L*(A[ ,idx[[min_idx]]])%*%t(A[ ,idx[[min_idx]]]))
  }
  
  n <- nrow(X)
  p <- ncol(X)
  
  Pistar <- matrix(0, p, p)
  
  n.cores <- min(n.cores, parallel::detectCores() - 2)
  cl <- makeCluster(n.cores)
  registerDoParallel(cl)
  
  idx_list <- lapply(1:M, function(m) ((m - 1) * d + 1):(m * d))
  
  P.list <- foreach(l = 1:L) %dopar% {
    A <- generate_random_matrix(p, d, M, dist, ortho, ...)
    s <- sample(n, n / 3)
    best_proj_outer(X, Y, A, s, idx_list, version, kernel)
  }
  
  stopCluster(cl)
  
  Pistar <- Reduce("+", P.list)
  SVD <- svd(Pistar)
  
  list(U = SVD$u, D = SVD$d)
}