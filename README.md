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
The two main functions of the FadDoS estimator are:

* `FadDoS`: Fits multivariate FLR models with known tuning parameters and obtains coefficient function estimates with double-sparsity property.

```
FadDoS(Xt, y, intercept=TRUE, tps=NULL, nbasis, phi, lambda1, lambda2, adaptive=FALSE, lambdas, maxit, tol)
```

* `cv.FadDoS`: Does k-fold cross-validation (CV) for FadDoS, produces coefficient function estimates by optimal tuning parameter. 

```
cv.FadDoS(Xt, y, intercept=TRUE, tps=NULL, nbasis, phi, lambda1, lambda2,  adaptive=FALSE, lambdas, K, maxit, tol)
```
### Arguments
* `Xt`: List of multiple functional predictors. 
* `y`: Vector of scalar responses. 
* `intercept`: Intercepts of multivariate FLR model. The default is `TRUE`.
* `tps`: Domain of functional predictors. If `NULL` (default), the domain will be $[0,1]$. 
* `nbasis`: The number of B-spline basis specified by users. 
* `phi`: Smoothness parameter. In `cv.FadDoS`, it is a vector of tuning parameters. 
* `lambda1`: Local sparsity parameter. In `cv.FadDoS`, it is a vector of tuning parameters. 
* `lambda2`: Global sparsity parameter. In `cv.FadDoS`, it is a vector of tuning parameters. 
* `adaptive`: Adaptive penalization is used. If `FALSE`, it is FDoS estimator. 
* `lambdas`: Tuning parameters for initial estimators for computing adaptive weights.
* `K`: The value of `K` used for the K-fold CV.
* `maxit`: Maximum iteration to stop the algorithm. 
* `tol`: Tolerance to stop the algorithm. 

### Values
* `intercept`: The estimated intercept. 
* `beta.hat`:  The list of estimated double-sparsity coefficient functions.
* `yhat`: The fitted values of response.
* `phi`: Optimal smoothness parameter. 
* `lambda1`: Optimal local sparsity parameter. 
* `lambda2`: Optimal global sparsity parameter
* `score`:  The matrix of SSE over all possible parameters for K-fold CV only in `cv.FadDoS`, 

## Examples

Suppose that $Y_{i}$ be the scalar response and $X_{ij}(t)$ be the $j$th functional covariate for subject $i$ observed at time $t$ in domain $\mathcal{T}$, the multivariate FLR model is as follows:
$$Y_{i}  =\mu + \sum_{j=1}^{10}\int_{0}^{1}X_{ij}(t)\beta_{j}(t)dt + \epsilon_{i}, \ \ i=1,\dots,n,$$
 where $\epsilon_{i} \sim N(0,1)$. Three different types of coefficient functions are considered, each representing a unique condition: 
 
(i) $\beta_{1}(t)$ has a zero subregion

$$
\beta_{1}(t) = \begin{cases}
    2\sin{(3\pi t)} & 0 \leq t \leq 1/3 \\
    0 & \text{if} \ 1/3 < t < 2/3 \\
    -2\sin{(3\pi t)} & 2/3 \leq t \leq 1;
\end{cases}
$$

(ii) $\beta_{2}(t)$ has no zero subregion but two crossings at zero, such that $\beta_{2}(t) = 1.5t^{2} + 2\sin{(3\pi t)}$; 

(iii) $\beta_{j}(t)=0$ for $j=3,\dots,10$, indicating that the functional covariate has no contribution to the response throughout the entire time domain. 

We first generate simulation data and randomly select 200 samples as our training set.
```
dat <- generate.test.data(nruns = 1, N=1200, beta1 = beta0.func, beta2 = beta1.func, beta3 = beta2.func)
train.index <- sample(1:1200, 200)
Xt <- lapply(dat[[1]]$x, function(i) i[train.index,])
y  <- dat[[1]]$y[train.index,] 
```


With pre-specified hyper-parameters, we are able to fit multivariate FLR model. The call the FadDoS estimator by using `FadDoS`.

```
phi <- 5e-5 #smoothness parameter
lambda1 <- 1000 #local sparsity parameter
lambda2 <- 5 #global sparsity parameter

result <- FadDoS(Xt=Xt, y=y, intercept=T, nbasis=30, phi=phi, lambda1 = lambda1, lambda2 = lambda2, adaptive=TRUE, maxit=5000, tol=0.0005, lambdas = c(1e-3,1e-4,1e-5,1e-6))
```
Generally, we would like to find optimal parameters using K-fold cross-validation for multivariate FLR models. We can use 	`cv.FadDoS` as follows. 

```
phi <- seq(3e-5, 7-5, length.out=5)
lambda1 <- seq(1000, 1200, length.out=5)
lambda2 <- seq(3, 5, length.out=5)

cv.result <- cv.FadDoS(Xt=Xt, y=y, intercept=T, nbasis=30, tps=time, phi=phi, lambda1 = lambda1, lambda2 = lambda2, adaptive = TRUE, K = 5, maxit = 5000, tol=0.0005)
```
## Files 
* `admm.R`: The main algorithm of our proposed FadDoS estimator. 
* `generate.data.R`: The procedures to generate training and testing data in simulation studies. 
* `demo.simu.R`: A demo script for a simulation study of the FadDoS estimator.
* `TUG.demo.rds`: A example dataset of Time Up \& Go test of the elderly. 
* `TUG.demo.R`:  A 3D animation of the example dataset. 

## Authors
* Cheng Cao, Jiguo Cao, Hailiang Wang, Kwok-Leung Tsui, Xinyue Li




