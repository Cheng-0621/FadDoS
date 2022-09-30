#generate simulation data

generate.test.data <- function(nruns=10, N, beta1, beta2, beta3, p=10)
{
  N <- N
  intercept <- 1
  rangeval <- c(0,1)
  time <- seq(0, 1, length.out = 1000)
  dat <- data.generator(beta1,beta2,beta3,N,nruns,intercept,rangeval,time,p=p)
  return(dat)
}

data.generator <- function(beta1, beta2, beta3,N,nSimu,intercept,rangeval,time,norder=4,nknots=50, p=10)
{
  if (!require(fda)) install.packages("fda")
  nbasis=nknots+norder-2
  data.basis <- create.bspline.basis(rangeval=rangeval,norder=norder,nbasis=nbasis) #B(t)
  
  #N=sample size for dataset
  
  G3 <- G2 <- G1 <- matrix(0,nbasis,1)
  for(j in 1:nbasis) G1[j] <- inner.prod(beta1,data.basis,j) #computation of B(t)*\beta(t)
  for(j in 1:nbasis) G2[j] <- inner.prod(beta2,data.basis,j)
  for(j in 1:nbasis) G3[j] <- inner.prod(beta3,data.basis,j)
  
  dat <- list()
  
  for(i in 1:nSimu)
  {
    cMat <- list()
    for (j in 1:p){
      cMat[[j]] <- matrix(rnorm(N*nbasis),N,nbasis) #a_{ij} 
    }
    
    cMat2 <- cMat[[1]]
    cMat3 <- cMat[[2]]
    yTru <- intercept + cMat2 %*% G2 + cMat3 %*% G3  #X(t)\beta(t)=\sum(a_{ij}B(t)*\beta(t))
    ySig <- sd(yTru)
    y <-  yTru + matrix(rnorm(N),N,1) * ySig / 2
    
    #construct datasets fd object
    Xt <- list()
    for (j in 1:p){
      Xt[[j]] <- cMat[[j]] %*% t(eval.basis(time, data.basis)) #X(t) matrix
    }
    names(Xt) <- paste0("x", 1:p)
    dat[[i]] <- list(x=Xt,y=y,intercept=intercept)
  }
  dat
}

#true coefficient functions
beta0.func <- function(t){
  y <- matrix(0,1,length(t))
  y
}

beta1.func <- function(t){
  t1 <- t[t<=1/3]
  t3 <- t[t>=2/3]
  y <- matrix(0,1,length(t))
  y[t<=1/3] <- 2*sin(3*pi*t1)
  y[t>=2/3] <- -2*sin(3*pi*t3)
  y
}

beta2.func <- function(t){
  y <- 1.5*t^2 + 2*sin(3*pi*t) 
  y
}

#computation of basis function
inner.prod <- function(f,basis,j)
{
  rng <- getbasisrange(basis)
  knots <- c(rng[1],basis$params,rng[2])
  nbasis <- basis$nbasis
  norder <- basis$nbasis - length(knots) + 2
  
  a <- rng[1]
  if(j-norder > 0) a <- knots[j-norder+1]
  
  b <- rng[2]
  if (j <= nbasis-norder) b <- knots[j+1]
  
  bfun <- function(t)
  {
    mat <- eval.basis(t,basis)
    z <- t(mat[,j])
  }
  
  y <- integrate(function(t) {f(t)*bfun(t)},a,b)
  y$value
}



