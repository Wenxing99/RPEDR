## Required packages ##
library(earth)
library(doParallel)
library(foreach)

### RPEDR ###
RPEDR <- function(X, Y, d, L, M, D = 5, ortho = FALSE, dist="mixed", ..., version = "kernel", 
                  kernel = "normal", n.cores) {
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  Pistar <- matrix(0,p,p)
  
  if (n.cores > (parallel::detectCores()-2)) {
    n.cores <- parallel::detectCores()-2
  }
    
  cl <- makeCluster(n.cores)  
  registerDoParallel(cl)
  
  P.list <- foreach(l=1:L) %dopar% { 
    
    if (ortho) {
      #################################################
      ### New random matrices generating code block ###
      
      A <- matrix(, nrow = p, ncol = 0)
      
      for (m in 1:M) {
        temp <- matrix(rt(p*p,1),p,p)
        temp.Q <- qr.Q(qr(temp))
        A <- cbind(A, temp.Q[,1:d])
      }
      #################################################
    }
    
    else {
      ######################################################
      ### Original random matrices generating code block ###
      if (dist=="t") {
        A0 <- matrix(rt(M*d*p, ...), p, M*d)
      }
      else if (dist=="norm") {
        A0 <- matrix(rnorm(M*d*p, ...), p, M*d)
      }
      else if (dist=="cauchy") {
        A0 <- matrix(rcauchy(M*d*p, ...), p, M*d)
      }
      else if (dist=="mixed") {
        A0 <- matrix(0,p,M*d)
        unif.n <- runif(M)
        norm.ind <- unif.n > 0.5
        cauchy.ind <- unif.n <= 0.5
        norm.idx <- which(norm.ind==1)
        cauchy.idx <- which(cauchy.ind==1)
        num.norm.proj <- sum(norm.ind)
        num.cauchy.proj <- sum(cauchy.ind)
        A0.norm <- matrix(rnorm(num.norm.proj*d*p), p, num.norm.proj*d)
        A0.cauchy <- matrix(rt(num.cauchy.proj*d*p, 1), p, num.cauchy.proj*d)
        idx.col.norm <- c()
        idx.col.cauchy <- c()
        for (i in seq_len(d)) {
          idx.col.norm <- sort(union(idx.col.norm, (norm.idx-1)*d+i))
          idx.col.cauchy <- sort(union(idx.col.cauchy, (cauchy.idx-1)*d+i))
        }
        A0[,idx.col.norm] <- A0.norm
        A0[,idx.col.cauchy] <- A0.cauchy
      }
      
      A <- t(t(A0)/sqrt(diag(t(A0)%*%A0))) 
      
      ######################################################
    }
    
    Out <- rep(0,M)
    s <- sample(n,n/3)
    
    if (version == "kernel") {
      
      Z <- X%*%A # random projection, Z: n*(M*d)
      
      for (m in 1:M) {
        idx <- (m-1)*d + 1:d
        if (kernel == "normal") {
          fitted <- apply(Z[s,idx], 1, function(z){regpro::kernesti.regr(z, Z[-s,idx] ,Y[-s], h = 0.1,
                                                                         kernel="gauss")})
        }
        else if (kernel == "loclin-normal") {
          fitted <- apply(Z[s,idx], 1, function(z){loclin(z, Z[-s,idx] ,Y[-s], h = 1,
                                                          kernel="gauss")})
        }
        else if (kernel == "loclin-uniform") {
          fitted <- apply(Z[s,idx], 1, function(z){loclin(z, Z[-s,idx] ,Y[-s], h = 0.1,
                                                          kernel="uniform")})
        }
        else {
          kreg.estimator <- kreg(x=Z[-s,idx],y=Y[-s],grid=Z[s,idx], bandwidth = 0.01, kernel=kernel)
          fitted <- kreg.estimator$y
        }
        Out[m] <- mean((fitted-Y[s])^2)
      }
    }
    
    else if (version == "poly") {
      
      Z <- X%*%A # random projection, Z: n*(M*d)
      
      for (m in 1:M) {
        idx <- (m-1)*d + 1:d
        train.data <- cbind.data.frame(Y[-s], Z[-s,idx])
        colnames(train.data)[1] <- "Y"
        colnames(train.data)[-1] <- paste0("X", 1:(dim(train.data)[2]-1))
        formula <- as.formula(paste0("Y ~ .^2 + ", 
                                     paste0("I(", colnames(train.data)[-1], "^2)", collapse = "+")))
        print(formula)
        fit.RP <- lm(formula, data=train.data)
        newdata <- as.data.frame(Z[s,idx])
        colnames(newdata) <- paste0("X", 1:(dim(train.data)[2]-1))
        print(colnames(newdata))
        fitted <- predict(fit.RP, newdata)
        Out[m] <- mean((fitted-Y[s])^2)
      }
      
    }
    
    else if (version == "linear") {
      
      Z <- X%*%A # random projection, Z: n*(M*d)
      
      for (m in 1:M) {
        idx <- (m-1)*d + 1:d
        train.data <- cbind.data.frame(Y[-s], Z[-s,idx])
        colnames(train.data)[1] <- "Y"
        colnames(train.data)[-1] <- paste0("X", 1:(dim(train.data)[2]-1))
        fit.RP <- lm(Y ~ ., data=train.data)
        newdata <- as.data.frame(Z[s,idx])
        colnames(newdata) <- paste0("X", 1:(dim(train.data)[2]-1))
        fitted <- predict(fit.RP, newdata)
        Out[m] <- mean((fitted-Y[s])^2)
      }
    }
    
    else if (version == "mars") {
      
      Z <- X%*%A # random projection, Z: n*(M*d)
      
      for (m in 1:M) {
        idx <- (m-1)*d + 1:d
        train.data <- cbind.data.frame(Y[-s], Z[-s,idx])
        colnames(train.data)[1] <- "Y"
        fit.RP <- earth::earth(Y ~., data=train.data, degree=3)
        
        fitted <- predict(fit.RP, Z[s,idx])
        Out[m] <- mean((fitted-Y[s])^2)
      }
    }
    
    
    min_idx <- (which.min(Out)-1)*d + 1:d
    return(1/L*(A[ ,min_idx])%*%t(A[ ,min_idx]))
  }
  
  stopCluster(cl)
  
  for (l in 1:L) {
    Pistar <- Pistar + P.list[[l]]
  }
  
  SVD <- svd(Pistar)
  
  list(U=SVD$u, D=SVD$d)
  
}