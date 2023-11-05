if(!require(Matrix)) install.packages("Matrix")
if(!require(splines)) install.packages("splines")


shrinkage1 <- function(y, kappa) 
{
  #beta <- sapply(y, function(i) max(0, i - kappa) - max(0, -i - kappa)) they are equivalent 
  beta <- sapply(y, function(i) sign(i)*max(0, abs(i) - kappa)) 
  beta <- matrix(beta, ncol = 1)
}

shrinkage2 <- function(y, kappa)
{
  y.norm <- norm(y, type="2")
  beta <- (y/y.norm)*max(0, y.norm-kappa)
  beta
}

b_update <- function(U, y, b.old, z.old, u.old, intercept.old, inverseL, lambda2, rho, group) 
{
  p <- length(unique(group))
  for (l in 1:p)
  {
    l.index <- which(names(group)==paste0("x",l))
    
    #hat(U_{l})
    Ul <- U[, l.index] 
    Uhat <- rbind(Ul, sqrt(rho)*t(inverseL))
    
    #hat(r_{-l})
    U_l <- U[, -l.index]
    r_l <-  y -  intercept.old - U_l %*% b.old[-l.index,]
    rhat <- rbind(r_l, sqrt(rho)*matrix((z.old[l.index,] + u.old[l.index,]/rho), ncol=1))
    
    #nul is proximal parameter for Ul
    nul <- 5*max(eigen(t(Ul)%*%Ul + rho*diag(length(l.index)))$values)
    
    #shrinkage function
    bl <- b.old[l.index,]
    b.old[l.index,] <-  shrinkage2(bl - t(Uhat) %*% (Uhat %*% bl - rhat)/nul, lambda2[l]/nul)
  }
  b.old
}

z_update <- function(b.new, u.old, inverseL, lambda1, rho, group)
{
  z.new <- Db.new <- matrix(0,1,nrow=length(group))
  p <-  length(unique(group))
  for (l in 1:p){
    l.index <- which(names(group)==paste0("x",l))
    Db.new[l.index,] <- inverseL %*% b.new[l.index,]
    z.new[l.index,] <- shrinkage1(Db.new[l.index,] + u.old[l.index,]/rho, lambda1[l]/rho) 
  }
  return(list(z=z.new, Db=Db.new))
}

admm <- function(U, y, intercept, inverseL, group, lambda1, lambda2, rho, maxit=2000, tol = 0.0001) 
{
  #initialization
  b.new <- matrix(0, nrow=length(group))
  z.new <- matrix(0, nrow=length(group))
  u.new <- matrix(0, nrow=length(group))
  intercept.new <- 0
  
  for (nit in 1:maxit) {
    
    b.old <- b.new
    z.old <- z.new
    u.old <- u.new
    intercept.old <- intercept.new
    
    #b-update
    b.new <- b_update(U, y, b.old, z.old, u.old, intercept.old, inverseL, lambda2, rho, group)
    
    #z-update 
    z <- z_update(b.new, u.old, inverseL, lambda1, rho, group)
    z.new <- z$z
    Db.new <- z$Db
    
    if (intercept) {
      p <- length(unique(group))
      yhat <- 0
      for (l in 1:p)
      {
        l.index <- which(names(group)==paste0("x",l))
        Ul <- U[, l.index] 
        yhat <- Ul %*% solve(inverseL) %*% z.new[l.index, ]
      }
      intercept.new <- mean(y-yhat)
    } else {
      intercept.new <- 0
    }
    
    #u-update
    u.new <- u.old + rho*(Db.new - z.new)
    
    #stop criterion
    epsilon <- norm(b.new-b.old, type = "2")/max(1, norm(b.new, type = "2"))
    
    if (epsilon <= tol) break
    if (nit == maxit) print("The maximum iteration is reached.")
  }
  output <- list(beta=b.new, Dbeta=z.new, intercept=intercept.new)
  return(output)
}


pre.smooth <- function(Xt, tps=NULL, nbasis){
  
  if(is.matrix(Xt)) { # there is only one functional predictor, convert it into list
    X <- list()
    X[[1]] <- Xt
  } 
  if (!is.list(Xt)){
    stop("Xt must be a list of matrices")
  }
  
  X <- Xt
  p <- length(X)
  
  if(is.null(tps))  tps <- seq(0, 1, length.out = ncol(X[[1]]))
  
  rng <- c(tps[1], tps[length(tps)])
  zspl <- create.bspline.basis(rng, norder=5, nbasis=nbasis) #Z(t)
  
  X.sm <- list()
  for (i in 1:p)  X.sm[[i]] <- Data2fd(argvals = tps, t(X[[i]]),zspl,nderiv=2,lambda=1e-9)
  for (i in 1:p) X[[i]] <- t(eval.fd(tps, X.sm[[i]]))
  
  return(X)
}

#obtain adaptive weights from initial estimatora
#lambda is a vector of proposed tunning parameter
initial.estimator <- function(U.ad, y, Omega.diag, lambdas) {
  
  n <- length(y)
  idx <- sample(1:n,n) # shuffle data
  Kfold <- 5
  s <- ceiling(n/Kfold)
  err.cv <- matrix(0,1,Kfold)
  
  mse.lam <- rep(0, length(lambdas))
  for (i in 1:length(lambdas)){
    mse <- rep(0, Kfold)
    for (k in 1:Kfold)
    {
      test.index <- idx[(s*(k-1)+1):min(s*k,n)]
      train.index <- idx[-((s*(k-1)+1):min(s*k,n))]  
      U.train <- U.ad[train.index,]
      U.test <- U.ad[-train.index,]
      b.ad <- solve.default(t(U.train) %*% U.train + lambdas[i]*Omega.diag, tol=1*10^(-32)) %*% t(U.train) %*% y[train.index]
      mse[k] <- mean((y[-train.index] - U.test %*% b.ad)^2)
    }
    mse.lam[i] <- mean(mse)
  }
  opt.lam <- lambdas[which.min(mse.lam)]
  print(paste("The optimal lambda for initial estimator is", opt.lam))
  b.ad <- solve.default(t(U.ad) %*% U.ad + opt.lam*Omega.diag, tol=1*10^(-32)) %*% t(U.ad) %*% y
  b.ad
}

#main algorithm 
FadDoS <- function(Xt, y, intercept=TRUE, tps=NULL, nbasis, phi, lambda1, lambda2, 
                   adaptive=FALSE, lambdas=c(1e-3,1e-4,1e-5,1e-6), standardize=FALSE, maxit=5000, tol=0.0005) 
{ 
  # there is only one functional predictor, convert it into list
  if(is.matrix(Xt)) { 
    X <- list()
    X[[1]] <- Xt
  } 
  
  if (!is.list(Xt)){
    stop("Xt must be a list of matrices")
  }
  
  X <- Xt
  p <- length(X)
  n <- length(y)
  
  if(is.null(tps))  tps <- seq(0, 1, length.out = ncol(X[[1]]))
  
  # basis spline set up 
  rng <- c(tps[1], tps[length(tps)])
  bspl <- create.bspline.basis(rangeval = rng, norder=5, nbasis=nbasis)
  Bt <- eval.basis(tps, bspl)
  
  # differences of each element which is a suppport grid of the corresponding observed function
  delt <- tps[2] - tps[1] 
  #Psi <-  inprod(bspl, bspl)
  Psi <- delt * t(Bt) %*% Bt  #Psi matrix: approximate norm of bsplines assuming dense design
  dBt <- matrix(NA,nrow(Bt),ncol(Bt))
  
  if (!require(splines)) install.packages("splines")
  
  for (k in 1:ncol(Bt)) # computation of second derivatives of B(t)
  {
    iS <- interpSpline(tps, Bt[,k])
    dBt[,k] <- predict(iS, tps, deriv = 2)$y
  }
  #Omega <- inprod(dBt, dBt)
  Omega <- delt * t(dBt) %*% dBt #Omega matrix: approximate norm of 2nd deriv of bsplines assuming dense design
  
  # define group
  group <- rep(1:p, each=nbasis)
  names(group) <- rep(paste0("x", 1:p), each=nbasis)
  
  if (adaptive == TRUE) {
    
    U.ad <- NULL
    for (i in 1:p){
      Ui <-  delt * X[[i]] %*% Bt #U=X(t)B(t)
      U.ad <- cbind(U.ad, Ui)
    }
    
    if (intercept) U.ad <- cbind(1, U.ad)
    
    if (intercept) {
      Omega.list <- lapply(1:(p+1), function(x) Omega)
      Omega.list[[1]] <- 1
      Omega.diag <- bdiag(Omega.list)
    } else {
      Omega.list <- lapply(1:p, function(x) Omega)
      Omega.diag <- bdiag(Omega.list)
    }
    
    if (standardize) U.ad <- scale(U.ad)
    b.ad <- initial.estimator(U.ad, y, Omega.diag, lambdas=lambdas)
    
    #initial estimator
    w <- v <- rep(0, p)
    for (i in 1:p) {
      beta.hat.ad <-  Bt %*% b.ad[(1+(i-1)*nbasis):(i*nbasis)]
      w[i] <- 1/norm(beta.hat.ad, type="1")
      v[i] <- 1/norm(beta.hat.ad, type="2")
    } 
    lambda1 <- lambda1*w
    lambda2 <- lambda2*v
  } else {
    lambda1 <- rep(lambda1, p)
    lambda2 <- rep(lambda2, p)
  }
  
  result <- FadDoS.alg(X,y,intercept,Bt,group,delt,phi,lambda1,lambda2,Psi,Omega,standardize,maxit=maxit,tol=tol)
  output <- list(intercept=result$intercept, beta.hat=result$beta.hat, yhat=result$yhat,
                 phi=phi, lambda1=lambda1, lambda2=lambda2)
  
  class(output) <- "FadDoS"
  return(output)
}

# phi/lambda1/lambda2 should be vectors of user-specified tunning parameters
cv.FadDoS <- function(Xt, y, intercept=TRUE, tps=NULL, nbasis, phi, lambda1, lambda2, 
                      adaptive=FALSE, lambdas=c(1e-3,1e-4,1e-5,1e-6), K=5, standardize=FALSE, maxit=5000, tol=0.0005)
{
  
  if(is.matrix(Xt)) {  # there is only one functional predictor, convert it into list
    X <- list()
    X[[1]] <- Xt
  } 
  
  if (!is.list(Xt)){
    stop("Xt must be a list of matrices")
  }
  
  X <- Xt
  p <- length(X)
  n <- length(y)
  
  if(is.null(tps))  tps <- seq(0, 1, length.out = ncol(X[[1]]))
  
  # basis spline set up 
  rng <- c(tps[1], tps[length(tps)])
  bspl <- create.bspline.basis(rangeval = rng, norder=5, nbasis=nbasis)
  Bt <- eval.basis(tps, bspl)
  
  # differences of each element which is a suppport grid of the corresponding observed function
  delt <- tps[2] - tps[1] 
  #Psi <-  inprod(bspl, bspl)
  Psi <- delt * t(Bt) %*% Bt  #Psi matrix: approximate norm of bsplines assuming dense design
  dBt <- matrix(NA,nrow(Bt),ncol(Bt))
  
  if (!require(splines)) install.packages("splines")
  
  for (k in 1:ncol(Bt)) # computation of second derivatives of B(t)
  {
    iS <- interpSpline(tps, Bt[,k])
    dBt[,k] <- predict(iS, tps, deriv = 2)$y
  }
  #Omega <- inprod(dBt, dBt)
  Omega <- delt * t(dBt) %*% dBt #Omega matrix: approximate norm of 2nd deriv of bsplines assuming dense design
  
  # define group
  group <- rep(1:p, each=nbasis)
  names(group) <- rep(paste0("x", 1:p), each=nbasis)
  
  #initial estimator
  w <- v <- rep(0, p)
  
  if (adaptive == TRUE) {
    U.ad <- NULL
    for (i in 1:p){
      Ui <-  delt * X[[i]] %*% Bt #U=X(t)B(t)
      U.ad <- cbind(U.ad, Ui)
    }
    
    Omega.list <- lapply(1:p, function(x) Omega)
    Omega.diag <- bdiag(Omega.list)
    
    if (standardize) U.ad <- scale(U.ad)
    b.ad <- initial.estimator(U.ad, y, Omega.diag, lambdas=lambdas)
    
    #initial estimator
    w <- v <- rep(0, p)
    for (i in 1:p) {
      beta.hat.ad <-  Bt %*% b.ad[(1+(i-1)*nbasis):(i*nbasis)]
      w[i] <- 1/norm(beta.hat.ad, type="1")
      v[i] <- 1/norm(beta.hat.ad, type="2")
    } 
  } else {
    w <- rep(1, p)
    v <- rep(1, p)
  }
  
  err.phi <- list()
  
  for (nphi in 1:length(phi)){
    
    err.tune <- matrix(rep(0,length(lambda1)*length(lambda2)), nrow=length(lambda1))
    for (nlam1 in 1:length(lambda1)){
      for (nlam2 in 1:length(lambda2)){
        
        lambda1.w <- w*lambda1[nlam1]
        lambda2.v <- v*lambda2[nlam2]
        
        idx <- sample(1:n,n) # shuffle data
        Kfold <- K
        s <- ceiling(n/Kfold)
        
        err.cv <- matrix(0,1,Kfold)
        for(k in 1:Kfold)
        {
          test.index <- idx[(s*(k-1)+1):min(s*k,n)]
          train.index <- idx[-((s*(k-1)+1):min(s*k,n))]
          newX.train <- lapply(X, function(l) l[train.index,])
          fit.cv <- FadDoS.alg(newX.train, y[train.index], intercept,Bt,group,delt,phi=phi[nphi],
                               lambda1=lambda1.w, lambda2=lambda2.v,Psi,Omega,standardize,maxit=maxit,tol=tol)
          newX.test <- lapply(X, function(l) l[test.index,])
          yhat <- predict.FadDoS(fit.cv, newXt=newX.test, tps=NULL)
          rss <- sum((yhat-y[test.index])^2)
          err.cv[k] <- rss
        }
        err <- sum(err.cv)
        err.tune[nlam1,nlam2] <- err 
      }
    }
    err.phi[[nphi]] <- err.tune
  }
  
  phi.index <- which.min(unlist(lapply(err.phi, function(l) min(l))))
  phi.opt <- phi[phi.index]
  
  h <- which(err.phi[[phi.index]]==min(err.phi[[phi.index]]), arr.ind = T)
  lambda1.opt <- lambda1[h[1]]
  lambda2.opt <- lambda2[h[2]]
  
  if (adaptive == TRUE){
    lam1 <- lambda1.opt*w
    lam2 <- lambda2.opt*v
  } else {
    lam1 <- rep(lambda1.opt, p)
    lam2 <- rep(lambda2.opt, p)
  }
  
  result <- FadDoS.alg(X,y,intercept,Bt,group,delt,phi=phi.opt,lambda1=lam1,lambda2=lam2,Psi,Omega,
                       standardize,maxit=maxit,tol=tol)
  output <- list(intercept=result$intercept, beta.hat=result$beta.hat, yhat=result$yhat, score=err.phi,
                 phi=phi.opt, lambda1=lambda1.opt, lambda2=lambda2.opt)
  
  return(output)
}


#lambda1/lambda2 are length p vectors of tuning parameters 
FadDoS.alg <- function(X,y,intercept,Bt,group,delt,phi,lambda1,lambda2,Psi,Omega,standardize=FALSE,maxit,tol)
{
  p <- length(X)
  
  K <- (Psi + phi * Omega)/(delt^2) ## K matrix
  inverseL <- backsolve(chol(K), x = diag(ncol(K)))  ## inverse of cholesky of K
  
  U <- NULL
  for (i in 1:p){
    Ui <- (1/delt) * delt * X[[i]] %*% Bt %*% inverseL #U=X(t)B(t)inv(L^T)
    U <- cbind(U, Ui)
  }
  
  if (standardize) U <- scale(U)
  
  fit <- admm(U=U, y=y, intercept=intercept, inverseL=inverseL, group=group, 
              lambda1=lambda1, lambda2=lambda2, rho=10, maxit=maxit, tol=tol) 
  intercept <- fit$intercept
  coefficient <- fit$Dbeta #coefficient
  names(coefficient) <- names(group)
  
  beta.hat <- list()
  for (i in 1:p) beta.hat[[i]] <-  (1/delt) * Bt %*% coefficient[which(names(group)==paste0("x",i))]
  
  Jmat <- list()
  for (i in 1:p) Jmat[[i]] <- delt * X[[i]] %*% beta.hat[[i]] 
  yhat <- Reduce("+", Jmat, accumulate=F)
  
  output <- list(intercept=intercept, beta.hat=beta.hat, yhat=yhat)
  class(output) <- "FadDoS"
  return(output)
  
}

predict.FadDoS <- function(object, newXt, tps=NULL)
{
  if (!is.list(newXt))  stop("newXt must be a list of matrices")
  if(is.null(tps))  tps <- seq(0, 1, length.out = ncol(newXt[[1]]))
  p <- length(newXt)
  
  # differences of each element which is a suppport grid of the corresponding observed function
  delt <- tps[2] - tps[1] 
  Jmat <- list()
  for (i in 1:p){
    Jmat[[i]] <- delt * newXt[[i]] %*% object$beta.hat[[i]] 
  }
  
  y <- Reduce("+", Jmat, accumulate=F) + object$intercept
  y
}



