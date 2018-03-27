EnKFforecast <- function (filterData, nAhead=1) {
    
    #This function predicts future system states and observations for
    #the ensemble Kalman Filter (EnKF).
    
    #Last update UKFforecast function: 2016/06/06
    
    #Arguments:
    #-filterdata is an object as returned by EnKF.
    #-nAhead is the number of steps ahead for which to predict system states
    # and observations.
    
    mod <- filterData
    
    GGfunction <- mod$GGfunction
    FFfunction <- mod$FFfunction
    size <- mod$size
    xa <- mod$xa
    
    modFuture <- mod$mod
    nobs <- nrow(as.matrix(mod$a))
    ym <- ncol(mod$f)
    lastObsIndex <- NROW(mod$m)
    modFuture$C0 <- mod$C[[lastObsIndex]]
    modFuture$m0 <- mod$m[lastObsIndex,]
    
    mod <- modFuture
    p <- length(mod$m0)
    a <- rbind(mod$m0, matrix(0, nAhead, p))
    R <- vector("list", nAhead + 1)
    R[[1]] <- mod$C0
    f <- matrix(0, nAhead, ym)
    Q <- vector("list", nAhead)
    
    for (it in 1:nAhead) {
        
        ##future states
        
        ##time update (forecast step)
        wks <- rmvnorm(n=size, sigma=as.matrix(mod$W))
        #a priori state estimate
        xf <- t(matrix(apply(xa, 1, GGfunction, k=nobs + it), nrow=p)) + wks
        a[it + 1, ] <- apply(xf, 2, mean)
        #variance of a priori state estimate
        R[[it + 1]] <- diag(1/(size-1) * tcrossprod(t(xf)-a[it, ]))
        
        ##future observations
        
        ##measurement update
        vks <- rmvnorm(n=size, sigma=as.matrix(mod$V))
        #predicted measurement
        yp <- t(matrix(apply(xf, 1, FFfunction, k=nobs + it), nrow=ym))
        f[it, ] <- apply(yp, 2, mean)
        #variance of predicted measurement
        Q[[it]] <- diag(1/(size-1) * tcrossprod(t(yp)-f[it, ]) + mod$V)
        
        xa <- xf
    }
    a <- a[-1, , drop = FALSE]
    R <- R[-1]
    ans <- list(a = a, R = R, f = f, Q = Q)
    
    return(ans)
}