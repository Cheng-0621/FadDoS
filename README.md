# Functional Adaptive Double-Sparsity Estimator for Multivariate Functional Linear Regression Models (FadDoS)

## Description
Sensor devices have been increasingly used in engineering and health studies recently,
and the captured multi-dimensional signals can be assessed under multivariate functional
linear regression (FLR) model. Kinect sensor is one of the advanced motion sensing devices providing human activity tracking of twenty-five joints, which can be used to recognize activities,
characterize mobility and evaluate physical health. Motivated by high-dimensional Kinect
sensor signals, we propose a novel functional adaptive double-sparsity (FadDoS) estimator,
which can achieve global sparsity via functional variable selection and local sparsity via sparse
coefficient estimation simultaneously. The combination of global and local sparsity is termed as
double-sparsity. The proposed estimator achieve double-sparsity through functional generalization
of sparse group lasso and adaptive penalization. By solving an optimization problem that combines
regularization and smoothing splines in a single objective function, we can provide the estimates
of coefficient functions with nice double-sparsity and smoothness control.


## Main Functions
* `FadDoS`: Fits multivariate FLR models with known tuning parameters and obtains coefficient function estimates with double-sparsity property.

```
FadDoS(Xt, y, intercept=TRUE, tps=NULL, nbasis, phi, lambda1, lambda2, adaptive=FALSE, lambdas, maxit, tol)
```

* `cv.FadDoS`: Does k-fold cross-validation for FadDoS, produces coefficient function estimates by optimal tuning parameter. 

```
cv.FadDoS(Xt, y, intercept=TRUE, tps=NULL, nbasis, phi, lambda1, lambda2,  adaptive=FALSE, lambdas, K, maxit, tol)
```

###Arguments
* Xt: 