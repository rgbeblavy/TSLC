# library(MASS); library(parallel); library(foreach); library(doParallel); library(ggplot2)

### This is a copy of the functions from the CMMB R package.
### These will be needed for running the functions in tslc.R, but are not new for this R package.
### See the following link: https://github.com/mbsohn/cmmb
get.init.lambda <- function(n.sample, n.feature, tol){
  f <- function(x, p) {(qnorm(1-x/p))^4 + 2*((qnorm(1-x/p))^2) - x;}
  k <- uniroot(f, lower=0, upper=n.feature-1, tol=tol, p=n.feature)$root
  lambda_0 <- sqrt(2/n.sample)*qnorm(1-k/n.feature)
  return(lambda_0)
}

comp_lasso_probit <- function(y, x, l_constraint, lam, mean_value, tol, max.iter){
  n <- nrow(x); p <- ncol(x)
  q <- 2*y-1; bet0 <- rep(0, p)
  gramC <- crossprod(l_constraint); diagC <- diag(gramC); dediagC <- gramC - diag(diagC)
  iter <- 0; ksi <- 0; ksi0 <- 1
  while (sum(abs(ksi-ksi0))>tol & iter<max.iter){
    ksi0 <- ksi; iter2 <- 0
    mu <- 1; bet <- rep(1,p)/p
    while (sum(abs(bet-bet0))>tol & iter2<max.iter){
      bet0 <- bet
      xbeta <- x %*% bet0
      phi <- dnorm(q*xbeta); PHI <- pnorm(q*xbeta)
      lmbd <- q*phi/PHI; LMBD <- diag(as.numeric(lmbd*(xbeta+lmbd)))
      u <- xbeta + solve(LMBD) %*% lmbd
      gramX <- crossprod(x, LMBD%*%x); diagX <- diag(gramX); dediagX <- gramX - diag(diagX)
      covXY <- crossprod(x, LMBD%*%u)
      term0 <- (covXY-dediagX%*%bet)/n - mu*(t(l_constraint)%*%ksi + dediagC%*%bet)
      term2 <- diagX/n + mu*diagC
      for(j in 1:p){
        term1 <- sign(term0[j])*max(0, abs(term0[j])-lam)
        bet[j] <- term1/term2[j]
        dif <- bet[j] - bet0[j]
        term0 <- term0 - dediagX[,j]*dif/n - dediagC[,j]*dif*mu
      }
      if(mean_value==TRUE){
        stp <- bet-bet0; t <- 1
        xbeta <- x %*% bet
        phi <- dnorm(q*xbeta); PHI <- pnorm(q*xbeta)
        lmbd <- q*phi/PHI
        ll_cur <- sum(log(PHI))
        eps <- crossprod(lmbd, x%*%stp); ll_new <- ll_cur + abs(eps)
        iter3 <- 1
        while ((ll_new+0.3*t*eps>ll_cur) & iter3<max.iter){
          bet_new <- bet0 + t*stp
          xbeta_new <- x %*% bet_new
          PHI_new <- pnorm(q*xbeta_new)
          ll_new <- sum(log(PHI_new))
          t <- 0.9*t
          iter3 <- iter3 + 1
        }
        bet <- bet_new
      }
      iter2 <- iter2 + 1
    }
    dif2 <- l_constraint%*%bet
    ksi <- ksi + dif2
    term0 <- term0 - mu*t(l_constraint)%*%dif2
    iter <- iter + 1
  }
  return(list(beta=bet, u=u, LMBD=LMBD))
}

comp_slr_probit <- function(y, X, contr, lam0, flag_mv, tol, max.iter){
  n <- nrow(X); p <- ncol(X); contr2 <- t(as.matrix(contr))
  if(is.null(lam0)) lam0 <- get.init.lambda(n, p, tol)
  sigma <- 1; sigma_2 <- 1; sigma_s <- 0.5; iter <- 1
  while (abs(sigma-sigma_s)>(tol*10) & iter<(max.iter/20)){
    iter <- iter + 1
    sigma <- (sigma_s + sigma_2)/2
    lam <- sigma*lam0
    rslt <- comp_lasso_probit(y, X, contr2, lam, flag_mv, tol, max.iter)
    bet2 <- rslt$beta
    u <- rslt$u
    LMBD <- rslt$LMBD
    s <- sum(abs(bet2)>tol); s <- min(s, n-1)
    sigma_s <- sqrt(t(u-X%*%bet2) %*% LMBD %*% (u-X%*%bet2))/sqrt(n-s-1)
    sigma_2 <- sigma
  }
  if(iter==max.iter) print("Not converge!")
  sigma <- sigma_s
  bet <- bet2
  return(list(beta=bet, u=u, LMBD=LMBD, sigma=sigma, lambda0=lam0))
}

comp.debias.probit <- function(y, M, x, p.adj.method="BH", lam0=NULL, flag_mv=FALSE, tol=1e-3, max.iter=2000){
  if(!is.vector(y)) y <- as.vector(y)
  n <- length(y); k <- ncol(M)
  if(is.null(x)){
    Z <- cbind(1, log(M))
    n.x <- 0
  } else{
    Z <- cbind(1, log(M), x)
    if(!is.matrix(x)) x <- as.matrix(x)
    n.x <- ncol(x)
  }
  n.vrs <- 1 + k + n.x
  contr <- c(0, rep(1/sqrt(k), k), rep(0, n.x))
  Z.til <- Z %*% (diag(n.vrs)-tcrossprod(contr))
  est.param <- comp_slr_probit(y, Z.til, contr, lam0, flag_mv, tol, max.iter)
  gam0 <- est.param$lambda0/3
  # zbeta <- Z.til %*% est.param$beta; q <- 2*y-1
  Sig <- crossprod(Z.til, est.param$LMBD%*%Z.til)/n
  Sig2 <- Sig - diag(diag(Sig))
  Q <- diag(n.vrs) - tcrossprod(contr)
  M.til <- matrix(0, n.vrs, n.vrs)
  for(i in 1:n.vrs){
    gam <- gam0/2
    while(gam<0.5){
      gam <- gam*2
      mi <- rep(1,n.vrs)
      mi0 <- rep(0,n.vrs)
      iter <- 1
      while(sum(abs(mi-mi0))>tol & iter<100){
        mi0 <- mi
        for(j in 1:n.vrs){
          v <- -Sig2[j,]%*%mi+Q[j,i]
          mi[j] <- sign(v)*max(0, abs(v)-gam)/Sig[j,j]
        }
        iter <- iter + 1
      }
      if(iter<100) break
    }
    M.til[i,] <- mi
  }
  M.til <- Q%*%M.til
  debias.B <- est.param$beta + M.til%*%t(Z.til)%*%(est.param$LMBD%*%(est.param$u-Z.til%*%est.param$beta))/n
  cov.debias.B <- as.numeric(est.param$sigma^2)*M.til%*%Sig%*%t(M.til)/n
  sd.debias.B <- sqrt(diag(cov.debias.B))
  rslt.df <- data.frame(Est=debias.B, SD=sd.debias.B, pval=2*pnorm(abs(debias.B)/sd.debias.B, lower.tail=FALSE))
  rslt.df$adjp <- p.adjust(rslt.df$pval, method=p.adj.method)
  if(n.x > 0){
    rownames(rslt.df) <- c("Intercept", paste0("M", 1:k), paste0("X", 1:n.x))
  } else{
    rownames(rslt.df) <- c("Intercept", paste0("M", 1:k))
  }
  return(rslt.df)
}

lasso_probit <- function(y, x, lam, mean_value, tol, max.iter){
  n <- nrow(x); p <- ncol(x)
  q <- 2*y-1; bet0 <- rep(0, p)
  iter1 <- 0
  mu <- 1; bet <- rep(1,p)/p
  while (sum(abs(bet-bet0))>tol & iter1<max.iter){
    bet0 <- bet
    xbeta <- x %*% bet0
    phi <- dnorm(q*xbeta); PHI <- pnorm(q*xbeta)
    lmbd <- q*phi/PHI; LMBD <- diag(as.numeric(lmbd*(xbeta+lmbd)))
    #u <- xbeta + solve(LMBD) %*% lmbd
    u <- xbeta + ginv(LMBD) %*% lmbd
    gramX <- crossprod(x, LMBD%*%x); diagX <- diag(gramX); dediagX <- gramX - diag(diagX)
    covXY <- crossprod(x, LMBD%*%u)
    term0 <- (covXY-dediagX%*%bet)/n
    term2 <- diagX/n
    for(j in 1:p){
      term1 <- sign(term0[j])*max(0, abs(term0[j])-lam)
      bet[j] <- term1/term2[j]
      dif <- bet[j] - bet0[j]
      term0 <- term0 - dediagX[,j]*dif/n
    }
    if(mean_value==TRUE){
      stp <- bet-bet0; t <- 1
      xbeta <- x %*% bet
      phi <- dnorm(q*xbeta); PHI <- pnorm(q*xbeta)
      lmbd <- q*phi/PHI
      ll_cur <- sum(log(PHI))
      eps <- crossprod(lmbd, x%*%stp); ll_new <- ll_cur + abs(eps)
      iter2 <- 1
      while ((ll_new+0.3*t*eps>ll_cur) & iter2<max.iter){
        bet_new <- bet0 + t*stp
        xbeta_new <- x %*% bet_new
        PHI_new <- pnorm(q*xbeta_new)
        ll_new <- sum(log(PHI_new))
        t <- 0.9*t
        iter2 <- iter2 + 1
      }
      bet <- bet_new
    }
    iter1 <- iter1 + 1
  }
  return(list(beta=bet, u=u, LMBD=LMBD))
}

slr_probit <- function(y, X, lam0, flag_mv, tol, max.iter){
  n <- nrow(X); p <- ncol(X)
  if(is.null(lam0)) lam0 <- get.init.lambda(n, p, tol)
  sigma <- 1; sigma_2 <- 1; sigma_s <- 0.5; iter <- 1
  while (abs(sigma-sigma_s)>(tol*10) & iter<(max.iter/20)){
    iter <- iter + 1
    sigma <- (sigma_s + sigma_2)/2
    lam <- sigma*lam0
    rslt <- lasso_probit(y, X, lam, flag_mv, tol, max.iter)
    bet2 <- rslt$beta
    u <- rslt$u
    LMBD <- rslt$LMBD
    s <- sum(abs(bet2)>tol); s <- min(s, n-1)
    sigma_s <- sqrt(t(u-X%*%bet2) %*% LMBD %*% (u-X%*%bet2))/sqrt(n-s-1)
    sigma_2 <- sigma
  }
  if(iter==max.iter) print("Not converge!")
  sigma <- sigma_s
  bet <- bet2
  return(list(beta=bet, u=u, LMBD=LMBD, sigma=sigma, lambda0=lam0))
}

debias.probit <- function(y, M, x, p.adj.method="BH", lam0=NULL, flag_mv=FALSE, tol=1e-3, max.iter=2000){
  if(!is.vector(y)) y <- as.vector(y)
  n <- length(y); k <- ncol(M)
  if(is.null(x)){
    Z <- cbind(1, log(M))
    n.x <- 0
  } else{
    Z <- cbind(1, log(M), x)
    if(!is.matrix(x)) x <- as.matrix(x)
    n.x <- ncol(x)
  }
  n.vrs <- 1 + k + n.x
  est.param <- slr_probit(y, Z, lam0, flag_mv, tol, max.iter)
  gam0 <- est.param$lambda0
  Sig <- crossprod(Z, est.param$LMBD%*%Z)/n
  Sig2 <- Sig - diag(diag(Sig))
  Q <- diag(n.vrs)
  M.hat <- matrix(0, n.vrs, n.vrs)
  for(i in 1:n.vrs){
    gam <- gam0/2
    stp.rho <- max(abs(Sig2))/diag(Sig)
    stp.crit <- stp.rho/(1+stp.rho)
    while(gam<stp.crit[i]){
      gam <- gam*2
      mi <- rep(1,n.vrs)
      mi0 <- rep(0,n.vrs)
      iter <- 1
      while(sum(abs(mi-mi0))>tol & iter<max.iter){
        mi0 <- mi
        for(j in 1:n.vrs){
          v <- -Sig2[j,]%*%mi+Q[j,i]
          mi[j] <- sign(v)*max(0, abs(v)-gam)/Sig[j,j]
        }
        iter <- iter + 1
      }
      if(iter<max.iter) break
    }
    M.hat[i,] <- mi
  }
  debias.B <- est.param$beta + M.hat%*%t(Z)%*%(est.param$LMBD%*%(est.param$u-Z%*%est.param$beta))/n
  cov.debias.B <- as.numeric(est.param$sigma^2)*M.hat%*%Sig%*%t(M.hat)/n
  sd.debias.B <- sqrt(diag(cov.debias.B))
  rslt.df <- data.frame(Est=debias.B, SD=sd.debias.B, pval=2*pnorm(abs(debias.B)/sd.debias.B, lower.tail=FALSE))
  rslt.df$adjp <- p.adjust(rslt.df$pval, method=p.adj.method)
  if(n.x > 0){
    rownames(rslt.df) <- c("Intercept", paste0("M", 1:k), paste0("X", 1:n.x))
  } else{
    rownames(rslt.df) <- c("Intercept", paste0("M", 1:k))
  }
  return(rslt.df)
}
