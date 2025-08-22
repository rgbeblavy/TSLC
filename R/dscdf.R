library(MASS); library(parallel); library(foreach); library(doParallel); library(ggplot2)

get.init.lambda <- function(n.sample, n.feature, tol){
  f <- function(x, p) {(qnorm(1-x/p))^4 + 2*((qnorm(1-x/p))^2) - x;}
  k <- uniroot(f, lower=0, upper=n.feature-1, tol=tol, p=n.feature)$root
  lambda_0 <- sqrt(2/n.sample)*qnorm(1-k/n.feature)
  return(lambda_0)
}

comp_lasso_probit <- function(y, x, l_constraint, lam, tol, max.iter, mean_value){
  n <- nrow(x); p <- ncol(x)
  q <- 2*y-1; bet0 <- rep(0, p)
  gramC <- crossprod(l_constraint); diagC <- diag(gramC); dediagC <- gramC - diag(diagC)
  mu <- 1; iter <- 0; ksi <- 0; ksi0 <- 1
  while (sum(abs(ksi-ksi0))>tol & iter<max.iter){
    ksi0 <- ksi; iter2 <- 0
    bet <- rep(1,p)/p
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

comp_slr_probit <- function(y, X, contr, lam0, tol, max.iter, mean_value){
  n <- nrow(X); p <- ncol(X); contr2 <- t(as.matrix(contr))
  if(is.null(lam0)) lam0 <- get.init.lambda(n, p, tol)
  sigma <- 1; sigma_2 <- 1; sigma_s <- 0.5; iter <- 1
  while (abs(sigma-sigma_s)>(tol*10) & iter<(max.iter/20)){
    iter <- iter + 1
    sigma <- (sigma_s + sigma_2)/2
    lam <- sigma*lam0
    rslt <- comp_lasso_probit(y, X, contr2, lam, tol, max.iter, mean_value)
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

comp.debias.probit <- function(y, M, x=NULL, p.adj.method="BH", lam0=NULL, tol=1e-3, max.iter=2000, mean_value=FALSE){
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
  est.param <- comp_slr_probit(y, Z.til, contr, lam0, tol, max.iter, mean_value)
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
    rownames(rslt.df) <- c("Intercept", paste0("F", 1:k), paste0("X", 1:n.x))
  } else{
    rownames(rslt.df) <- c("Intercept", paste0("F", 1:k))
  }
  return(rslt.df)
}

######################################
comp_lasso <- function(y, x, l_constraint, lam, tol, max.iter){
  n <- nrow(x); p <- ncol(x); k <- dim(l_constraint)[1]
  gramC <- crossprod(l_constraint); gramX <- crossprod(x)
  diagC <- diag(gramC); diagX <- diag(gramX)
  dediagC <- gramC - diag(diagC); dediagX <- gramX - diag(diagX)
  covXY <- crossprod(x, y)
  mu <- 1; bet <- rep(1,p)/p; bet0 <- rep(0,p); iter <- 0
  ksi <- rep(0,k); ksi0 <- rep(1,k)
  term0 <- (covXY-dediagX%*%bet)/n - mu*(t(l_constraint)%*%ksi + dediagC%*%bet)
  term2 <- diagX/n + mu*diagC
  while (sum(abs(ksi-ksi0))>tol && iter<max.iter){
    ksi0 <- ksi; iter2 <- 0; bet0 <- bet0 + 1
    while (sum(abs(bet-bet0))>tol && iter2<1000){
      bet0 <- bet
      for(j in 1:p){
        term1 <- sign(term0[j])*max(0, abs(term0[j])-lam)
        bet[j] <- term1/term2[j]
        dif <- bet[j] - bet0[j]
        term0 <- term0 - dediagX[,j]*dif/n - dediagC[,j]*dif*mu
      }
      iter2 <- iter2 + 1
    }
    dif2 <- l_constraint%*%bet
    ksi <- ksi + dif2
    term0 <- term0 - mu*t(l_constraint)%*%dif2
    iter <- iter + 1
    return(bet)
  }
}

### Scaled lasso for compositional data
comp_slr <- function(y, X, contr, lam0, tol, max.iter){
  n <- nrow(X); p <- ncol(X)
  mn.y <- mean(y); mn.X <- colMeans(X)
  cent.X <- (X-tcrossprod(rep(1, n), mn.X))
  contr2 <- t(as.matrix(contr))
  cent.y <- y-mn.y
  if(is.null(lam0)) lam0 <- get.init.lambda(n, p, tol)
  sigma <- 1; sigma_2 <- 1; sigma_s <- 0.5; iter <- 1
  while (abs(sigma-sigma_s)>(tol*10) & iter<(max.iter/20)){
    iter <- iter + 1
    sigma <- (sigma_s + sigma_2)/2
    lam <- sigma*lam0
    bet2 <- comp_lasso(cent.y, cent.X, contr2, lam, tol, max.iter)
    s <- sum(abs(bet2)>0.001)
    s <- min(s, n-1)
    sigma_s <- base::norm(cent.y-cent.X%*%bet2, type="2")/sqrt(n-s-1)
    sigma_2 <- sigma
  }
  if(iter>max.iter) print("Not converge!")
  sigma <- sigma_s
  bet <- bet2
  intercp <- mn.y - mn.X%*%bet
  return(list(beta=bet, intercept=intercp, sigma=sigma, lambda0=lam0))
}

comp.debias <- function(y, M, x=NULL, p.adj.method="BH", lam0=NULL, tol=1e-3, max.iter=2000){
  if(!is.vector(y)) y <- as.vector(y)
  n <- length(y); k <- ncol(M)
  if(is.null(x)){
    Z <- log(M)
    n.x <- 0
  } else{
    Z <- cbind(log(M), x)
    if(!is.matrix(x)) x <- as.matrix(x)
    n.x <- ncol(x)
  }
  n.vrs <- k + n.x
  contr <- c(rep(1/sqrt(k), k), rep(0, n.x))
  Z.til <- Z %*% (diag(n.vrs)-tcrossprod(contr))
  est.param <- comp_slr(y, Z.til, contr, lam0, tol, max.iter)
  cent.y <- scale(y, scale=FALSE)
  cent.Z.til <- scale(Z.til, scale=FALSE)
  gam0 <- est.param$lambda0/3; Sig <- crossprod(cent.Z.til)/n
  Sig2 <- Sig - diag(diag(Sig))
  Q <- diag(n.vrs) - tcrossprod(contr)
  M.til <- matrix(0, n.vrs, n.vrs)
  for(i in 1:n.vrs){
    gam <- gam0/2
    while(gam<0.5){
      gam <- gam*2
      mi <- rep(1,n.vrs)
      mi0 <- rep(0,n.vrs)
      iter <- 1;
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
    M.til[i,] <- mi
  }
  M.til <- Q%*%M.til
  debias.B <- est.param$beta + M.til%*%t(cent.Z.til)%*%(cent.y-cent.Z.til%*%est.param$beta-rep(est.param$intercept,n))/n
  cov.debias.B <- as.numeric(est.param$sigma^2)*M.til%*%Sig%*%t(M.til)/n
  sd.debias.B <- sqrt(diag(cov.debias.B))
  rslt.df <- data.frame(Est=debias.B, SD=sd.debias.B, pval=2*pnorm(abs(debias.B)/sd.debias.B, lower.tail=FALSE))
  rslt.df$adjp <- p.adjust(rslt.df$pval, method=p.adj.method)
  if(n.x > 0){
    rownames(rslt.df) <- c(paste0("T", 1:k), paste0("X", 1:n.x))
  } else{
    rownames(rslt.df) <- c(paste0("T", 1:k))
  }
  return(rslt.df)
}







