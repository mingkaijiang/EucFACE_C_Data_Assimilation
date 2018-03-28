analysis_3 <- function(A, s, p, obs, i,
                       err_var, err_type, 
                       err_var_obs, err_var_obs, 
                       ens_var, q,
                       ens_var_obs, q_obs) {
    
    ### The standard analysis eqn: A = A + Pe H^T(H Pe H^T + Re)^-1 (D - H A)    (eq 48)
    ### can be reformed using D' = D - HA, Pe = A'(A')^T, Re = YY^T
    ### (where Y symbolises Gamma) such that it
    ### becomes A = A + A' A'^T H^T(HA' A'^T H^T + YY^T)^-1 D'
    
    # Minimum of nrobs and nrens
    nrmin = min(s$nrobs+1, s$nrens)
    sig = seq(0, nrmin)
    D_mean = seq(0, s$nrobs)
    E_mean = seq(0, s$nrobs)
    S_mean = seq(0, s$nrobs)
    HA_mean = seq(0, s$nrobs)
    #A_mean = seq(0, s$nrobs)
    A_mean = seq(0, s$nrens)

    ## identity matrix
    I = diag(1, s$nrens,s$nrens)
    U = matrix(0, s$nrobs,nrmin)
    
    # measurement matrix
    D = matrix(0, s$nrobs,s$nrens)
    S = matrix(0, s$nrobs,nrmin)
    
    # measurement ensemble of perturbation
    E = matrix(0, s$nrobs,s$nrens)
    
    # measurement operator that maps the model state to the measurement
    H = matrix(0, s$nrobs,s$ndims)
    
    HA = matrix(0, s$nrobs,s$nrens)
    ES = matrix(0, s$nrobs,s$nrens)
    X1 = matrix(0, nrmin,s$nrens)
    X2 = matrix(0, nrmin,s$nrens)
    X3 = matrix(0, s$nrobs,s$nrens)
    X4 = matrix(0, s$nrens,s$nrens)
    Reps = matrix(0, s$ndims,s$nrobs)
    A_tmp = matrix(0, s$ndims,s$nrens)
    # N1 <- matrix(1/s$nrens, s$nrens, s$nrens)
    
    ## ensemble perturbation matrix
    A_dash = matrix(0, s$ndims, s$nrens)
    
    ## calculate the ensemble mean at each timestep
    A_mean <- apply(A, 1, mean)
    
    ## ensemble perturbation matrix, i.e. eq 51
    A_dash <- A - A_mean

    ## ensemble covariance matrix
    C_A <- 1 / (s$nrens - 1) * tcrossprod(A_dash)

    ## Measurement matrix
    for (j in 1:s$ndims) {
            D[j,] <- obs[j,i] + rnorm(s$nrens, 0.0, err_var_obs[j,i])  # need to update error variance
    }
    
    ## Measurement ensemble of perturbations
    for (j in 1:s$ndims) {
        E[j,] <- D[j,] - obs[j,i] 
    }
    
    ## Measurement error covariance
    C_D <- 1 / (s$nrens - 1) * tcrossprod(E)
    
    
    ## Measurement operator
    H <- 
        
    ## D' = D - HA
    D_dash <- D - H %*% A
    
    ## S = HA'
    S <- H %*% A_dash
            
    ## C matrix = SS^T + (N-1)C_D    
    C <- tcrossprod(S) + (s$nrens - 1) * C_D
            
    ## X matrix = I + S^T * C^-1 * D'    
    X <- I + t(S) %*% C^-1 %*% D_dash     
    
    ## analysis matrix
    A_A <- A %*% X
    
        
    ## predicted measurement
    A_pred <- t(matrix(apply(A, 1, FFfunction, k=i), nrow=nrow(B)))
    B_mean[,i] <- apply(A_pred, 1, mean)
    
    ## covariance of predicted measurement    ## eq 13
    Qy <- tcrossprod(t(A_pred)-B_mean[,i]) / nrmin
    
    ## cross covariance between a priori state estimate and predicted measurement   ## eq 14
    Qxy <- 1/(s$nrens-1) * tcrossprod(t(A)-A_mean[,i], t(A_pred)-B_mean[,i])
    
    ## a posteriori estimates (analysis step)
    # Kalman gain
    Kk <- crossprod(t(Qxy), solve(Qy, tol=1e-30))
    
    # a posteriori state estimate
    yk <- B[,i] # + t(vks)
    xa <- t(t(A) + crossprod(t(Kk), as.matrix(yk-t(A_pred))))
    A[,i + 1] <- apply(xa, 1, mean)
    
    # extract a posteriori error variance
    ens_var <- diag(1/(s$nrens-1) * tcrossprod(t(xa)-A[,i + 1]))
    
    
    sig_sum = 0.0
    sig_sum1 = 0.0
    

    # Compute HA' + E
    ES = S + E
    
    # Compute SVD of HA'+E -> U and sig, using Eispack
    U <- svd(ES)
    
    
    for (i in c(1:nrmin)) {
        sig[i] <- sig[i] ** 2
    }
         

    sigsum = sum(sig)
    sigsum1 = 0.0
    nrsigma = 0
    for (i in c(1:nrmin)) {
        if (sigsum1/sigsum < 0.999) {
            nrsigma = nrsigma + 1
            sigsum1 = sigsum1 + sig[i]
            sig[i] = 1.0/sig[i]
        } else {
            sig[i:nrmin] = 0.0
        }
    }

    # Compute X1
    for (j in c(1:s$nrobs)) {
        for (k in c(1:nrmin)) {
            X1[k, j] <- sig[k] * U$v[j, k]
        }
    }

    X2 <- X1 %*% D

    # Compute X3 = U * X2
    X3 <- U %*% X2
    
    # Compute final analysis
    if(2 * ndim * nrobs > nrens * (nrobs+ndim)) {
        X4 <- S %*% X3
        for (i in 1:s$nrens) {
            X4[i,i] = X4[i,i] + 1.0
        }
        
        iblkmax = min(s$ndim, 200)
        multa(A, X4, ndim, nrens, iblkmax)
    } else {
        Reps <- A %*% S
        A <- Reps %*% X3
    }
    
    
    return(list(A=A, ens_var=ens_var, q=q,
                B=obs, ens_var_obs=ens_var_obs, q_obs=q_obs))
}
    


multa <- function(A, X, ndim, nrens, iblkmax) {
    for (ia in seq(1, ndim, iblkmax)) {
        ib = min(ia + iblkmax - 1, ndim)
        v(1:(ib - ia + 1), 1:nrens) = A(ia:ib, 1:nrens)
        A[ia, 1] <- v[1,1] %*% X[1,1]
    }

}