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

