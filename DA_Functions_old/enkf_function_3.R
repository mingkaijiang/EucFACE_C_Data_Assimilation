enkf_function_3 <- function(y, mod, GGfunction, FFfunction, size=50,
                       logLik=FALSE, simplify=FALSE) {
    
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
        
        if (!any(whereNA <- is.na(y[,i]))) {
            
            ##time update (forecast step)
            wks <- rmvnorm(n=size, sigma=as.matrix(mod$W))
            
            #a priori state estimate
            xf <- t(matrix(apply(xa, 1, GGfunction, k=i), nrow=p)) + wks
            a[,i] <- apply(xf, 2, mean)
            
            #covariance of a priori state estimate:
            #1/(size-1) * tcrossprod(t(xf)-a[i, ])
            
            ##measurement update
            vks <- rmvnorm(n=size, sigma=as.matrix(mod$V))
            
            #predicted measurement
            yp <- t(matrix(apply(xf, 1, FFfunction, k=i), nrow=ym))
            f[,i] <- apply(yp, 2, mean)
            
            #covariance of predicted measurement
            Qy <- 1/(size-1) * tcrossprod(t(yp)-f[,i]) + mod$V
            
            #cross covariance between a priori state estimate and predicted measurement
            Qxy <- 1/(size-1) * tcrossprod(t(xf)-a[,i], t(yp)-f[,i])
            
            ##a posteriori estimates (analysis step)
            
            #Kalman gain
            Kk <- crossprod(t(Qxy), solve(Qy, tol=1e-30))
            #a posteriori state estimate
            yk <- y[,i] + t(vks)
            xa <- t(t(xf) + crossprod(t(Kk), as.matrix(yk-t(yp))))
            m[,i + 1] <- apply(xa, 2, mean)
            
            #extract a posteriori error variance
            C[[i + 1]] <- diag(1/(size-1) * tcrossprod(t(xa)-m[,i + 1]))
            
            #compute log-likelihood
            if (logLik) {
                e <- as.matrix(y[,i] - f[,i])
                ll <- ll + ym*log(2*pi) + sum(log(eigen(Qy)$values)) + 
                    crossprod(e, tcrossprod(solve(Qy, tol=1e-30), t(e)))}
        }
        
        else {
            if (all(whereNA)) {
                
                ##time update (forecast step)
                wks <- rmvnorm(n=size, sigma=as.matrix(mod$W))
                #a priori state estimate
                xf <- t(matrix(apply(xa, 1, GGfunction, k=i), nrow=p)) + wks
                a[,i] <- apply(xf, 2, mean)
                
                #covariance of a priori state estimate
                Qx <- 1/(size-1) * tcrossprod(t(xf)-a[,i])
                
                ##measurement update
                vks <- rmvnorm(n=size, sigma=as.matrix(mod$V))
                
                #predicted measurement
                yp <- t(matrix(apply(xf, 1, FFfunction, k=i), nrow=ym))
                f[,i] <- apply(yp, 2, mean)
                
                ##a posteriori estimates
                xa <- xf
                
                #a posteriori state estimate
                m[,i + 1] <- a[,i]
                
                #extract a posteriori error variance
                C[[i + 1]] <- diag(Qx)
            }
            
            else {
                good <- !whereNA
                
                ##time update (forecast step)
                wks <- rmvnorm(n=size, sigma=as.matrix(mod$W))
                
                #a priori state estimate
                xf <- t(matrix(apply(xa, 1, GGfunction, k=i), nrow=p)) + wks
                a[,i] <- apply(xf, 2, mean)
                
                #covariance of a priori state estimate:
                #1/(size-1) * tcrossprod(t(xf)-a[i, ])
                
                ##measurement update
                vks <- rmvnorm(n=size, sigma=as.matrix(mod$V))
                
                #predicted measurement
                yp <- t(matrix(apply(xf, 1, FFfunction, k=i), nrow=ym))
                f[,i] <- apply(yp, 2, mean)
                
                #covariance of predicted measurement
                Qy <- 1/(size-1) * tcrossprod(t(yp[good,])-f[good,i]) + mod$V[good, good]
                
                #cross covariance between a priori state estimate and predicted measurement
                Qxy <- 1/(size-1) * tcrossprod(t(xf)-a[,i], t(yp[good,])-f[good, i])
                
                ##a posteriori estimates (analysis step)
                
                #Kalman gain
                Kk <- crossprod(t(Qxy), solve(Qy, tol=1e-30))
                
                #a posteriori state estimate
                yk <- y[,i] + t(vks)
                xa <- t(t(xf) + crossprod(t(Kk), as.matrix(yk[,good]-t(yp)[,good])))
                m[,i + 1] <- apply(xa, 2, mean)
                
                #extract a posteriori error variance
                C[[i + 1]] <- diag(1/(size-1) * tcrossprod(t(xa)-m[,i + 1]))
                
                #compute log-likelihood
                if (logLik) {
                    e <- as.matrix(y[good, i] - f[good, i])
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