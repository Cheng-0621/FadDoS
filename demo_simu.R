source("generate.data.R")
source("admm.R")

################################
####Generate Simulation Data####
################################

dat <- generate.test.data(nruns = 1, N=1200, beta1 = beta0.func, beta2 = beta1.func, beta3 = beta2.func)
train.index <- sample(1:1250, 200)
Xt <- lapply(dat[[1]]$x, function(i) i[train.index,])
y  <- dat[[1]]$y[train.index,] 


################################
###Fit Multivariate FLR Model###
################################

#pre-specified hyperparameters
phi <- 5e-5
lambda1 <- 1000
lambda2 <- 5
time <- seq(0, 1, length.out = 1000)

result <- FadDoS(Xt=Xt, y=y, intercept=T, tps=time, nbasis=30, phi=phi, lambda1 = lambda1, lambda2 = lambda2,
                 adaptive=TRUE, maxit=5000, tol=0.0005, lambdas = c(1e-3,1e-4,1e-5,1e-6))

#coefficient functions estimates
betahat <- result$beta.hat

#draw all coefficient functions estimates
par(mfrow=c(1,3), mar=c(5, 5.5, 1.5, 1.2))
plot(x=time, y=beta1.func(time), type="l", lty=1, col=grey(0.9), lwd=2,
     xlab="t", ylab=expression(widehat(beta)(t)), cex.lab=1.5, cex.axis=1.2, cex.axis=1, ylim=c(-3,3))
polygon(x=c(time,1), y=c(beta1.func(time),0), col = grey(0.95), border = NA)
lines(x=time, y=betahat[[1]], col="#BC3C29FF", lwd=3, lty=2) 

plot(x=time, y=beta2.func(time), type="l", lty=1, col=grey(0.9), lwd=2,
     xlab="t", ylab=expression(widehat(beta)(t)), cex.lab=1.5, cex.axis=1.2, cex.axis=1, ylim=c(-3,3))
polygon(x=c(time,1), y=c(beta2.func(time),0), col = grey(0.95), border = NA)
lines(x=time, y=betahat[[2]], col="#BC3C29FF", lwd=3, lty=2) 

plot(x=time, y=beta0.func(time), type="l", lty=1, col=grey(0.9), lwd=2,
     xlab="t", ylab=expression(widehat(beta)(t)), cex.lab=1.5, cex.axis=1.2, cex.axis=1, ylim=c(-3,3))
polygon(x=c(time,1), y=c(beta0.func(time),0), col = grey(0.95), border = NA)
for (i in 3:10) lines(x=time, y=betahat[[i]], col="#BC3C29FF", lwd=3, lty=2) 


##########################
#####Cross Validation#####
##########################

phi <- seq(3e-5, 7-5, length.out=5)
lambda1 <- seq(1000, 1200, length.out=5)
lambda2 <- seq(3, 5, length.out=5)

cv.result <- cv.FadDoS(Xt=Xt, y=y, intercept=T, nbasis=30, tps=time, phi=phi, lambda1 = lambda1, lambda2 = lambda2,
                       adaptive = TRUE, K = 5, maxit = 5000, tol=0.0005)

est.int <- cv.result$intercept
est.int

sse <-cv.result$score #sse matrix of CV: row for lambda1, col for lambda2, list index for phi
sse

#coefficient functions estimates
cv.betahat <- cv.result$beta.hat

#optimal parameters
cv.result$phi
cv.result$lambda1
cv.result$lambda2

#draw all coefficient functions estimates
par(mfrow=c(1,3), mar=c(5, 5.5, 1.5, 1.2))
plot(x=time, y=beta1.func(time), type="l", lty=1, col=grey(0.9), lwd=2,
     xlab="t", ylab=expression(widehat(beta)(t)), cex.lab=1.5, cex.axis=1.2, cex.axis=1, ylim=c(-3,3))
polygon(x=c(time,1), y=c(beta1.func(time),0), col = grey(0.95), border = NA)
lines(x=time, y=cv.betahat[[1]], col="#0072B5FF", lwd=3, lty=2) 

plot(x=time, y=beta2.func(time), type="l", lty=1, col=grey(0.9), lwd=2,
     xlab="t", ylab=expression(widehat(beta)(t)), cex.lab=1.5, cex.axis=1.2, cex.axis=1, ylim=c(-3,3))
polygon(x=c(time,1), y=c(beta2.func(time),0), col = grey(0.95), border = NA)
lines(x=time, y=cv.betahat[[2]], col="#0072B5FF", lwd=3, lty=2) 

plot(x=time, y=beta0.func(time), type="l", lty=1, col=grey(0.9), lwd=2,
     xlab="t", ylab=expression(widehat(beta)(t)), cex.lab=1.5, cex.axis=1.2, cex.axis=1, ylim=c(-3,3))
polygon(x=c(time,1), y=c(beta0.func(time),0), col = grey(0.95), border = NA)
for (i in 3:10) lines(x=time, y=cv.betahat[[i]], col="#0072B5FF", lwd=3, lty=2) 

