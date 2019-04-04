require(foreach)

##Get empirical Bayesian estimates for pi, the prior probability each CpG lies in a given category
#lambda is the penalty term assigned to pi_{0,0}. The default is lambda = 10 (see ASH paper)
#The tolerance is the infinity norm < tol
EmBayesPi <- function(SXX, SXY, var.vec=c(0.05^2, 0.1^2, 0.2^2, 0.25^2), rho.vec=c(0, 1/3, 2/3, 1), lambda=10, pi0=NULL, tol=1e-5, max.iter=1e4, L=NULL) {
  p <- length(SXX)
  k.var <- length(var.vec); k.rho <- length(rho.vec)
  K <- 1 + k.var * (k.rho + 2)
  
  ##Create L, a K x p matrix##
  #L_{kg} = P(y_g | g is in group k)
  if (is.null(L)) {
    cl <- makeCluster(max(detectCores() - 1, 1))
    clusterExport(cl, c("SXX", "SXY", "var.vec", "rho.vec", "K", "k.var", "k.rho"), envir=environment())
    log.L <- rbind(rep(0,p), parSapply(cl, 1:p, Est.logL))
    stopCluster(cl)
    L <- exp(sweep(x = log.L, MARGIN = 2, STATS = apply(log.L, 2, max), FUN = "-"))   #The maximum entry in each column is 1
  }
    
  ll.vec <- rep(NA, max.iter)
  if (is.null(pi0)) {
    pi0 <- rep(1/K, K)
  }
  
  ##Run one iteration of turbo EM##
  #out.em.turbo <- turboem(par = pi0, fixptfn = update.em.pi, objfn=ll.em.pi, method="squarem", L=L, lambda=lambda, control.run=list(tol=tol, convtype="parameter", convfn.user=my.conv, maxiter=15e2, maxtime=1e16))
  #pi0 <- as.vector(out.em.turbo$pars)
  
  for (i in 1:max.iter) {
    pi0L <- L * pi0; sum.pi0L <- apply(pi0L, 2, sum)
    ll.vec[i] <- (sum(log( sum.pi0L )) + (lambda-1)*log(pi0[1]))/p
    #Plot ll#
    if (round(i/1e3) == (i/1e3)) {
      plot(which(!is.na(ll.vec)), ll.vec[!is.na(ll.vec)], xlab="Iteration", ylab="log-likelihood", type="l")
    }
    
    Probs.0 <- sweep(x = pi0L, MARGIN = 2, STATS = sum.pi0L, FUN = "/")
    Counts.0 <- apply(Probs.0, 1, sum); Counts.0[1] <- Counts.0[1] + lambda - 1
    pi1 <- Counts.0 / sum(Counts.0)
    if (max(abs(pi0-pi1)) < tol) {
      names.pi <- rep(NA, K)
      names.pi[1] <- "b_0=0, b_7=0"
      for (k in 1:k.var) {
        names.pi[k.var+1+k] <- paste0("b_0=0, tau=", as.character(round(sqrt(var.vec[k]), digits = 2)))
        names.pi[k+1] <- paste0("tau=", as.character(round(sqrt(var.vec[k]), digits = 2)), ", b_7=0")
        for (r in 1:k.rho) {
          names.pi[k.var*(2+r-1)+1+k] <- paste0("rho=", as.character(round(rho.vec[r], digits = 2)), ", tau=", as.character(round(sqrt(var.vec[k]), digits = 2)))
        }
      }
      names(pi0) <- names.pi
      return(list(prior.pi=pi0, Post.probs=Probs.0, out=1, ll=ll.vec[!is.na(ll.vec)], n.iter=i))
    }
    pi0 <- pi1
  }
  
  names.pi <- rep(NA, K)
  names.pi[1] <- "b_0=0, b_7=0"
  for (k in 1:k.var) {
    names.pi[k.var+1+k] <- paste0("b_0=0, tau=", as.character(round(sqrt(var.vec[k]), digits = 2)))
    names.pi[k+1] <- paste0("tau=", as.character(round(sqrt(var.vec[k]), digits = 2)), ", b_7=0")
    for (r in 1:k.rho) {
      names.pi[k.var*(2+r-1)+1+k] <- paste0("rho=", as.character(round(rho.vec[r], digits = 2)), ", tau=", as.character(round(sqrt(var.vec[k]), digits = 2)))
    }
  }
  names(pi1) <- names.pi
  return(list(prior.pi=pi1, Post.probs=Probs.0, out=0, ll=ll.vec, n.iter=i))
}


####Functions to parallelize####
#L is K x p. The columns of L go ((0,0), (var.vec,0), (0, var.vec), (var.vec, rho.vec[1]), ..., (var.vec, rho.vec[k.rho]))
#It is assume the last entry of rho.vec is 1
Est.logL <- function(g) {
  sxx <- SXX[[g]]
  sxy <- SXY[[g]]
  
  out <- rep(0, K-1)
  out[1:k.var] <- -1/2 * log(1 + var.vec * sxx[1,1]) + 1/2 * sxy[1]^2 / (1/var.vec + sxx[1,1])
  out[(k.var+1):(2*k.var)] <- -1/2 * log(1 + var.vec * sxx[2,2]) + 1/2 * sxy[2]^2 / (1/var.vec + sxx[2,2])
  out[(K-k.var):(K-1)] <- -1/2 * log(1 + var.vec*sum(sxx)) + 1/2 * sum(sxy)^2 / (1/var.vec + sum(sxx))
  
  count <- 2*k.var + 1
  for (k.r in 1:(k.rho-1)) {
    W <- diag(2); W[1,2] <- rho.vec[k.r]; W[2,1] <- W[1,2]
    Omega <- 1/(1-rho.vec[k.r]^2) * W; Omega[1,2] <- -Omega[1,2]; Omega[2,1] <- -Omega[2,1]
    for (k.v in 1:k.var) {
      v.k <- var.vec[k.v]
      tmp.v <- diag(2) + v.k * W %*% sxx
      out[count] <- -1/2 * log(tmp.v[1,1]*tmp.v[2,2] - tmp.v[1,2] * tmp.v[2,1]) + 1/2 * sum( sxy * solve(1/v.k * Omega + sxx, sxy) )
      count <- count + 1
    }
  }
  return(out)
}