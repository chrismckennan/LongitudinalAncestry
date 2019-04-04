require(parallel)
require(mvtnorm)

##Estimate posterior probabilities

PostProbs <- function(SXX, SXY, Post, var.vec, rho.vec) {
  names.cols <- c("P(b0, b7 > 0)", "P(b0, b7 < 0)", "P(b0 > 0)", "P(b0 < 0)", "P(b7 > 0)", "P(b7 < 0)", "P(b0 > 0 & b7 <= 0)", "P(b0 < 0 & b7 >= 0)", "P(b0 <= 0, b7 > 0)", "P(b0 >= 0, b7 < 0)", "E(b_0 | Data)", "E(b_7 | Data)")
  
  p <- length(SXX)
  k.var <- length(var.vec)
  k.rho  <- length(rho.vec)
  K <- 1 + (2+k.rho)*k.var
  
  cl <- makeCluster(max(detectCores() - 1, 1))
  clusterExport(cl=cl, varlist=c("SXX", "SXY", "Post", "k.var", "k.rho", "K", "var.vec", "rho.vec"), envir=environment())
  clusterEvalQ(cl=cl, expr={library(mvtnorm)})
  PPS <- t(parSapply(cl=cl, X=1:p, FUN=Est.pps))
  stopCluster(cl=cl)
  
  colnames(PPS) <- names.cols
  return(PPS)
}


####Functions to parallelize####
#returns ( P(b0, b7 > 0); P(b0, b7 < 0); P(b0 > 0); P(b0 < 0); P(b7 > 0); P(b7 < 0); P(b0 > 0 & b7 <= 0); P(b0 < 0 & b7 >= 0); P(b0 <= 0, b7 > 0); P(b0 >= 0, b7 < 0); E(b_0 | Data); E(b_7 | Data) )
#gamma = beta_7 - rho*beta_0
Est.pps <- function(g) {
  pi.g <- Post[,g]
  sxx.g <- SXX[[g]]
  sxy.g <- SXY[[g]]
  
  out <- rep(0, 12)
  pi.rhos <- pi.g[(K-k.var*k.rho+1):K]  #( (var.vec, rho[1]), ..., (var.vec, rho[k.rho]) )
  pi.10 <- pi.g[2:(k.var+1)]   #(var.vec, 0)
  pi.01 <- pi.g[(k.var+2):(2*k.var+1)]   #(0, var.vec)
  
  count <- 1
  for (s in 1:k.rho) {
    rho.s <- rho.vec[s]
    if (rho.s < 1) {Omega.s <- diag(2); Omega.s[1,2] <- -rho.s; Omega.s[2,1] <- -rho.s; Omega.s <- 1/(1-rho.s^2)*Omega.s}
    for (k in 1:k.var) {
      tau2.k <- var.vec[k]
      if (rho.s < 1) {
        Sigma.sk <- solve(sxx.g + 1/tau2.k * Omega.s)
        mu.sk <- as.vector(Sigma.sk %*% sxy.g)
        out[11] <- out[11] + pi.rhos[count] * mu.sk[1]; out[12] <- out[12] + pi.rhos[count] * mu.sk[2]
      } else {
        sigma2.sk <- 1/(sum(sxx.g) + 1/tau2.k)
        mu.sk <- sum(sxy.g) * sigma2.sk
        out[11:12] <- out[11:12] + rep(pi.rhos[count] * mu.sk, length=2)
      }
      
      ##All of P(b0, b7 > 0) and P(b0, b7 < 0), indices 1:2##
      if (rho.s < 1) {
        out[1] <- out[1] + pi.rhos[count] * pmvnorm(lower=c(0,0), upper=c(Inf, Inf), mean=mu.sk, sigma=Sigma.sk, algorithm=TVPACK())
        out[2] <- out[2] + pi.rhos[count] * pmvnorm(lower=c(-Inf,-Inf), upper=c(0, 0), mean=mu.sk, sigma=Sigma.sk, algorithm=TVPACK())
      } else {
        p0 <- pnorm(q = 0, mean = mu.sk, sd = sqrt(sigma2.sk))
        out[1] <- out[1] + pi.rhos[count] * (1-p0)
        out[2] <- out[2] + pi.rhos[count] * p0
      }
      
      ##Part of P(b0 > 0); P(b0 < 0); P(b7 > 0); P(b7 < 0), indices 3:6##
      if (rho.s < 1) {
        p0.b0 <- pnorm(q = 0, mean = mu.sk[1], sd = sqrt(Sigma.sk[1,1]))
        p0.b7 <- pnorm(q = 0, mean = mu.sk[2], sd = sqrt(Sigma.sk[2,2]))
      } else {
        p0.b0 <- pnorm(q = 0, mean = mu.sk, sd = sqrt(sigma2.sk))
        p0.b7 <- p0.b0
      }
      out[3] <- out[3] + pi.rhos[count] * (1-p0.b0)
      out[4] <- out[4] + pi.rhos[count] * p0.b0
      out[5] <- out[5] + pi.rhos[count] * (1-p0.b7)
      out[6] <- out[6] + pi.rhos[count] * p0.b7
      #The rest#
      if (s == 1) {
        mu.sk0 <- sxy.g[1]/(sxx.g[1,1] + 1/tau2.k); mu.sk7 <- sxy.g[2]/(sxx.g[2,2] + 1/tau2.k)
        p0.b0 <- pnorm(q = 0, mean = mu.sk0, sd = 1/sqrt(sxx.g[1,1] + 1/tau2.k))
        p0.b7 <- pnorm(q = 0, mean = mu.sk7, sd = 1/sqrt(sxx.g[2,2] + 1/tau2.k))
        out[3] <- out[3] + pi.10[k] * (1-p0.b0)
        out[4] <- out[4] + pi.10[k] * p0.b0
        out[5] <- out[5] + pi.01[k] * (1-p0.b7)
        out[6] <- out[6] + pi.01[k] * p0.b7
        out[11] <- out[11] + pi.10[k] * mu.sk0; out[12] <- out[12] + pi.01[k] * mu.sk7
      }
      
      ##Part of P(b0 > 0 & b7 <= 0); P(b0 < 0 & b7 >= 0); P(b0 <= 0, b7 > 0); P(b0 >= 0, b7 < 0), indices 7:10##
      if (rho.s < 1) {
        Sigma.tmp <- Sigma.sk; Sigma.tmp[1,2] <- -Sigma.tmp[1,2]; Sigma.tmp[2,1] <- -Sigma.tmp[2,1]
        mu.tmp <- mu.sk; mu.tmp[2] <- -mu.tmp[2]
        p.g0l0 <- pmvnorm(lower=c(0, 0), upper=c(Inf, Inf), mean=mu.tmp, sigma=Sigma.tmp, algorithm=TVPACK())
        p.l0g0 <- pmvnorm(lower=c(-Inf,-Inf), upper=c(0, 0), mean=mu.tmp, sigma=Sigma.tmp, algorithm=TVPACK())
        out[7] <- out[7] + pi.rhos[count] * p.g0l0
        out[8] <- out[8] + pi.rhos[count] * p.l0g0
        out[9] <- out[9] + pi.rhos[count] * p.l0g0
        out[10] <- out[10] + pi.rhos[count] * p.g0l0
      }
      #The rest#
      if (s == 1) {
        out[7] <- out[7] + pi.10[k] * (1-p0.b0)
        out[8] <- out[8] + pi.10[k] * p0.b0
        out[9] <- out[9] + pi.01[k] * (1-p0.b7)
        out[10] <- out[10] + pi.01[k] * p0.b7
      }
      
      count <- count + 1
    }
  }
  return(out)
}



###########Rarely used functions###########

##Estimate posterior probabilities for only a fraction of the sites

PostProbs.frac <- function(SXX, SXY, Post, Rho, var.vec, rho.vec, ind) {
  names.cols <- c("P(b0, b7 > 0)", "P(b0, b7 < 0)", "P(b0 > 0)", "P(b0 < 0)", "P(b7 > 0)", "P(b7 < 0)", "P(b0 > 0 & b7 <= 0)", "P(b0 < 0 & b7 >= 0)", "P(b0 <= 0, b7 > 0)", "P(b0 >= 0, b7 < 0)", "P(gamma > 0 | Data)", "P(gamma < 0 | Data)", "E(b_0 | Data)", "E(b_7 | Data)")
  
  p <- length(SXX)
  k.var <- length(var.vec)
  k.rho  <- length(rho.vec)
  K <- 1 + (2+k.rho)*k.var
  
  cl <- makeCluster(max(detectCores() - 1, 1))
  clusterExport(cl=cl, varlist=c("SXX", "SXY", "Post", "Rho", "k.var", "k.rho", "K", "var.vec", "rho.vec"), envir=environment())
  clusterEvalQ(cl=cl, expr={library(mvtnorm)})
  PPS <- matrix(NA, nrow=p, ncol=length(names.cols))
  PPS[ind,] <- t(parSapply(cl=cl, X=ind, FUN=Est.pps))
  stopCluster(cl=cl)
  
  colnames(PPS) <- names.cols
  return(PPS)
}
