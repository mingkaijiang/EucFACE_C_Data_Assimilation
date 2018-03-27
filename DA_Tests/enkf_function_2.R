enkf_function <- function(y, mod, GGfunction, FFfunction, size=50,
                       logLik=FALSE, simplify=FALSE) {
    
    #This function implements the Ensemble Kalman Filter (EnKF) as described
    #by Gillijns and colleagues in their 2006 paper "What Is the Ensemble Kalman
    #Filter and How Well Does it Work?".
    
    #Note that this EnKF assumes (time invariant) additive noise.
    
    #Last update EnKF function: 2016/06/06
    
    #Arguments:
    #-y is the data. The data can be a vector (one-dimensional measurements /
    # univariate time series) or a matrix (multidimensional measurements / 
    # multivariate time series). In the matrix case, each of the columns contains
    # the measurements on a single dimension (or, similarly, a single [univariate]
    # time series).
    # Missing values (NA's) are allowed.
    #-mod is a list with the components m0 (initial state estimates),
    # C0 (error covariances of the initial state estimates), V (measurement noise),
    # and W (process noise).
    #-GGfunction is a function with arguments x and k. The GGfunction specifies the
    # state transition function. The transition function is usually denoted as
    # f(x[k], u[k]), where x are the states estimates and u the control input at
    # time step k.
    # See details below for how to handle control input in the GGfunction.
    #-FFfunction is a function with arguments x and k. The FFfunction specifies the
    # observation/measurement function. The measurement function is usually denoted
    # as h(x[k]).
    #-size is a number specifying the ensemble size.
    #-logLik: if TRUE, then return the negative loglikelihood of the specified model.
    #-simplify: if TRUE, then do not include the data (=y) in the output.
    
    #Details:
    #It is possible to include control input in the GGfunction.
    #As an example, let us assume that there are two states (x1 and x2) and that
    #we have 70 time steps at which we want to obtain state estimates. At each of
    #these 70 time steps there will be (external) control input.
    #
    #We may specify an arbitrary GGfunction with control input as follows:
    #GGfunction <- function (x, k){
    #x1 <- x[1]; x2 <- x[2]
    #c(2*x1 + 1*x2 - 4,
    #  3*x2 - 1*x2 - 2)}
    #where, at each time step, the control input for the first state equation
    #is -4, and the control input for the second state equation -2.
    #
    #Alternatively, we may store the control input for the two state equations
    #at each of the 70 time steps in a 70x2 matrix.
    #Subsequently, we specify the GGfunction as:
    #GGfunction <- function (x, k){
    #x1 <- x[1]; x2 <- x[2]
    #c(2*x1 + 1*x2 - u[k,1],
    #  3*x2 - 1*x2 - u[k,2])}
    #where we have stored the control input at each of the 70 time steps
    #in a matrix called u.
    #
    #Note that the latter specification of the GGfunction is useful when the control
    #input is known to vary over time.
    
    mod1 <- mod
    y <- as.matrix(y)
    ym <- nrow(y)
    yAttr <- attributes(y)
    p <- length(mod$m0)
    
    m <- cbind(mod$m0, matrix(0, ncol = ncol(y), nrow = length(mod$m0)))
    a <- matrix(0, ncol = ncol(y), nrow = length(mod$m0))
    f <- matrix(0, ncol = ncol(y), nrow = nrow(y))
    C <- vector(1 + ncol(y), mode = "list")
    ll <- 0
    
    #extract variances of a posteriori state estimates at t=0
    C[[1]] <- diag(mod$C0)
    
    #generate ensemble at t=0
    xa <- rmvnorm(n=size, mean=m[,1], sigma=as.matrix(mod$C0))
    
    for (i in seq(length = ncol(y))) {
        
        if (!any(whereNA <- is.na(y[i, ]))) {
            
            ##time update (forecast step)
            wks <- rmvnorm(n=size, sigma=as.matrix(mod$W))
            #a priori state estimate
            xf <- t(matrix(apply(xa, 1, GGfunction, k=i), nrow=p)) + wks
            a[i, ] <- apply(xf, 2, mean)
            #covariance of a priori state estimate:
            #1/(size-1) * tcrossprod(t(xf)-a[i, ])
            
            ##measurement update
            vks <- rmvnorm(n=size, sigma=as.matrix(mod$V))
            #predicted measurement
            yp <- t(matrix(apply(xf, 1, FFfunction, k=i), nrow=ym))
            f[i, ] <- apply(yp, 2, mean)
            #covariance of predicted measurement
            Qy <- 1/(size-1) * tcrossprod(t(yp)-f[i, ]) + mod$V
            
            #cross covariance between a priori state estimate and predicted measurement
            Qxy <- 1/(size-1) * tcrossprod(t(xf)-a[i, ], t(yp)-f[i, ])
            
            ##a posteriori estimates (analysis step)
            
            #Kalman gain
            Kk <- crossprod(t(Qxy), solve(Qy, tol=1e-30))
            #a posteriori state estimate
            yk <- y[i, ] + t(vks)
            xa <- t(t(xf) + crossprod(t(Kk), as.matrix(yk-t(yp))))
            m[i + 1, ] <- apply(xa, 2, mean)
            #extract a posteriori error variance
            C[[i + 1]] <- diag(1/(size-1) * tcrossprod(t(xa)-m[i + 1, ]))
            
            #compute log-likelihood
            if (logLik) {
                e <- as.matrix(y[i, ] - f[i,])
                ll <- ll + ym*log(2*pi) + sum(log(eigen(Qy)$values)) + 
                    crossprod(e, tcrossprod(solve(Qy, tol=1e-30), t(e)))}
        }
        
        else {
            if (all(whereNA)) {
                
                ##time update (forecast step)
                wks <- rmvnorm(n=size, sigma=as.matrix(mod$W))
                #a priori state estimate
                xf <- t(matrix(apply(xa, 1, GGfunction, k=i), nrow=p)) + wks
                a[i, ] <- apply(xf, 2, mean)
                #covariance of a priori state estimate
                Qx <- 1/(size-1) * tcrossprod(t(xf)-a[i, ])
                
                ##measurement update
                vks <- rmvnorm(n=size, sigma=as.matrix(mod$V))
                #predicted measurement
                yp <- t(matrix(apply(xf, 1, FFfunction, k=i), nrow=ym))
                f[i, ] <- apply(yp, 2, mean)
                
                ##a posteriori estimates
                
                xa <- xf
                
                #a posteriori state estimate
                m[i + 1, ] <- a[i, ]
                #extract a posteriori error variance
                C[[i + 1]] <- diag(Qx)
            }
            
            else {
                good <- !whereNA
                
                ##time update (forecast step)
                wks <- rmvnorm(n=size, sigma=as.matrix(mod$W))
                #a priori state estimate
                xf <- t(matrix(apply(xa, 1, GGfunction, k=i), nrow=p)) + wks
                a[i, ] <- apply(xf, 2, mean)
                #covariance of a priori state estimate:
                #1/(size-1) * tcrossprod(t(xf)-a[i, ])
                
                ##measurement update
                vks <- rmvnorm(n=size, sigma=as.matrix(mod$V))
                #predicted measurement
                yp <- t(matrix(apply(xf, 1, FFfunction, k=i), nrow=ym))
                f[i, ] <- apply(yp, 2, mean)
                #covariance of predicted measurement
                Qy <- 1/(size-1) * tcrossprod(t(yp[, good])-f[i, good]) + mod$V[good, good]
                
                #cross covariance between a priori state estimate and predicted measurement
                Qxy <- 1/(size-1) * tcrossprod(t(xf)-a[i, ], t(yp[, good])-f[i, good])
                
                ##a posteriori estimates (analysis step)
                
                #Kalman gain
                Kk <- crossprod(t(Qxy), solve(Qy, tol=1e-30))
                #a posteriori state estimate
                yk <- y[i, ] + t(vks)
                xa <- t(t(xf) + crossprod(t(Kk), as.matrix(yk[good, ]-t(yp)[good, ])))
                m[i + 1, ] <- apply(xa, 2, mean)
                #extract a posteriori error variance
                C[[i + 1]] <- diag(1/(size-1) * tcrossprod(t(xa)-m[i + 1, ]))
                
                #compute log-likelihood
                if (logLik) {
                    e <- as.matrix(y[i, good] - f[i, good])
                    ll <- ll + sum(good)*log(2*pi) + sum(log(eigen(Qy)$values)) + 
                        crossprod(e, tcrossprod(solve(Qy, tol=1e-30), t(e)))}
            }
        }
    }
    ans <- list(m = m, C = C, a = a, f = f)
    
    attributes(ans$f) <- yAttr
    
    if (logLik)
        ans <- c(ans, logLik = 0.5*ll)
    
    if (simplify) 
        ans <- c(mod = list(mod1), size = list(size),
                 GGfunction = list(GGfunction), FFfunction = list(FFfunction),
                 xa = list(xa), ans)
    else {
        attributes(y) <- yAttr
        ans <- c(y = list(y), mod = list(mod1), size = list(size), 
                 GGfunction = list(GGfunction), FFfunction = list(FFfunction),
                 xa = list(xa), ans)
    }
    return(ans)
}