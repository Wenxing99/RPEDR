## Required packages ##
library(earth)
library(doParallel)
library(foreach)
library(regpro)
library(doRNG)

### RPEDR Function ###
RPEDR <- function(X, Y, d, L, M,
                              ortho    = FALSE,
                              dist     = "mixed",
                              ...,
                              n1       = ceiling(2*nrow(X)/3),
                              version  = "mars",
                              kernel   = "normal",
                              n.cores  = 10,
                              rng_seed = 123) {
  
  n <- nrow(X)
  p <- ncol(X)
  
  Pistar <- matrix(0, p, p)
  
  generate_single_projection <- function(p, d, dist, ortho, ...) {
    if (ortho) {
      temp   <- matrix(rt(p * p, 1), p, p)
      temp.Q <- qr.Q(qr(temp))
      A_proj <- temp.Q[, 1:d, drop = FALSE]
    } else {
      if (dist == "t") {
        A0 <- matrix(rt(d * p, ...), p, d)
      } else if (dist == "norm") {
        A0 <- matrix(rnorm(d * p, ...), p, d)
      } else if (dist == "cauchy") {
        A0 <- matrix(rcauchy(d * p, ...), p, d)
      } else if (dist == "mixed") {
        if (runif(1) > 0.5) {
          A0 <- matrix(rnorm(d * p, ...), p, d)
        } else {
          A0 <- matrix(rt(d * p, 1, ...), p, d)
        }
      } else {
        stop("Unknown dist: ", dist)
      }
      cs <- sqrt(colSums(A0^2))
      cs[cs == 0] <- 1
      A_proj <- sweep(A0, 2, cs, "/")
    }
    A_proj
  }
  
  
  compute_mse_for_Z <- function(Z_proj, Y, s, version, kernel) {
    ## Z_proj: n × d
    if (version == "kernel") {
      
      if (kernel == "normal") {
        fitted <- apply(Z_proj[s, , drop = FALSE], 1, function(z) {
          regpro::kernesti.regr(z,
                                Z_proj[-s, , drop = FALSE],
                                Y[-s],
                                h = 0.1,
                                kernel = "gauss")
        })
      } else if (kernel == "loclin-normal") {
        fitted <- apply(Z_proj[s, , drop = FALSE], 1, function(z) {
          regpro::loclin(z,
                         Z_proj[-s, , drop = FALSE],
                         Y[-s],
                         h = 1,
                         kernel = "gauss")
        })
      } else if (kernel == "loclin-uniform") {
        fitted <- apply(Z_proj[s, , drop = FALSE], 1, function(z) {
          regpro::loclin(z,
                         Z_proj[-s, , drop = FALSE],
                         Y[-s],
                         h = 0.1,
                         kernel = "uniform")
        })
      } else {
        kreg.estimator <- regpro::kreg(
          x         = Z_proj[-s, , drop = FALSE],
          y         = Y[-s],
          grid      = Z_proj[s, , drop = FALSE],
          bandwidth = 0.01,
          kernel    = kernel
        )
        fitted <- kreg.estimator$y
      }
      
    } else if (version == "poly") {
      
      train.data <- as.data.frame(cbind(Y = Y[-s], Z_proj[-s, , drop = FALSE]))
      colnames(train.data)[-1] <- paste0("X", 1:(ncol(train.data) - 1))
      formula <- as.formula(paste0(
        "Y ~ .^2 + ",
        paste0("I(", colnames(train.data)[-1], "^2)", collapse = "+")
      ))
      fit <- lm(formula, data = train.data)
      newdata <- as.data.frame(Z_proj[s, , drop = FALSE])
      colnames(newdata) <- paste0("X", 1:(ncol(train.data) - 1))
      fitted <- predict(fit, newdata)
      
    } else if (version == "linear") {
      
      train.data <- cbind.data.frame(Y[-s], Z_proj[-s, , drop = FALSE])
      colnames(train.data)[1]  <- "Y"
      colnames(train.data)[-1] <- paste0("X", 1:(ncol(train.data) - 1))
      fit <- lm(Y ~ ., data = train.data)
      newdata <- as.data.frame(Z_proj[s, , drop = FALSE])
      colnames(newdata) <- paste0("X", 1:(ncol(train.data) - 1))
      fitted <- predict(fit, newdata)
      
    } else if (version == "mars") {
      
      train.data <- cbind.data.frame(Y[-s], Z_proj[-s, , drop = FALSE])
      colnames(train.data)[1] <- "Y"
      fit <- earth::earth(Y ~ ., data = train.data, degree = 3)
      fitted <- predict(fit, Z_proj[s, , drop = FALSE])
      
    } else {
      stop("Unknown version: ", version)
    }
    
    mean((fitted - Y[s])^2)
  }
  
  
  # n.cores <- n.cores
  cl <- makeCluster(n.cores)
  on.exit({
    try(stopCluster(cl), silent = TRUE)
  }, add = TRUE)
  registerDoParallel(cl)
  
  parallel::clusterEvalQ(cl, {
    suppressMessages(library(earth))
    # suppressMessages(library(regpro))
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      RhpcBLASctl::blas_set_num_threads(1)
      RhpcBLASctl::omp_set_num_threads(1)
    } else {
      Sys.setenv(OMP_NUM_THREADS      = "1",
                 MKL_NUM_THREADS      = "1",
                 OPENBLAS_NUM_THREADS = "1")
    }
    NULL
  })
  
  P.list <- foreach(l = 1:L, .options.RNG = rng_seed) %dorng% {
    
    s <- sample(n, n-n1)      
    
    Out      <- rep(NA_real_, M)
    A_best   <- matrix(0, p, d)
    best_mse <- Inf
    
    for (m in 1:M) {
      A_proj <- generate_single_projection(p, d, dist, ortho, ...)
      Z_proj <- X %*% A_proj
      
      mse_m <- compute_mse_for_Z(Z_proj, Y, s, version, kernel)
      Out[m] <- mse_m
      
      if (mse_m < best_mse) {
        best_mse <- mse_m
        A_best   <- A_proj
      }
    }
    
    (A_best %*% t(A_best)) / L
  }
  
  stopCluster(cl)
  
  Pistar <- Reduce("+", P.list)
  SVD <- svd(Pistar)
  
  list(U = SVD$u, D = SVD$d)
}