analysis_3 <- function(A, s, obs, i) {
    
    ### The standard analysis eqn: A = A + Pe H^T(H Pe H^T + Re)^-1 (D - H A)    (eq 48)
    ### can be reformed using D' = D - HA, Pe = A'(A')^T, Re = YY^T
    ### (where Y symbolises Gamma) such that it
    ### becomes A = A + A' A'^T H^T(HA' A'^T H^T + YY^T)^-1 D'
    
    B_tmp = matrix(0, s$ndims, s$nrobs)
    
    
    
    
    
    
    sig_sum = 0.0
    sig_sum1 = 0.0
    
    # Minimum of nrobs and nrens
    nrmin = min(s$nrobs+1, s$nrens)
    
    I = matrix(0, s$nrens,s$nrens)
    U = matrix(0, s$nrobs,nrmin)
    D = matrix(0, s$nrobs,s$nrens)
    S = matrix(0, s$nrobs,nrmin)
    E = matrix(0, s$nrobs,nrmin)
    H = matrix(0, s$nrobs,s$ndims)
    HA = matrix(0, s$nrobs,s$nrens)
    ES = matrix(0, s$nrobs,s$nrens)
    X1 = matrix(0, nrmin,s$nrens)
    X2 = matrix(0, nrmin,s$nrens)
    X3 = matrix(0, s$nrobs,s$nrens)
    X4 = matrix(0, s$nrens,s$nrens)
    Reps = matrix(0, s$ndims,s$nrobs)
    A_tmp = matrix(0, s$ndims,s$nrens)
    A_dash = matrix(0, s$ndims,s$nrens)
    
    sig = seq(0, nrmin)
    D_mean = seq(0, s$nrobs)
    E_mean = seq(0, s$nrobs)
    S_mean = seq(0, s$nrobs)
    HA_mean = seq(0, s$nrobs)
    A_mean = seq(0, s$nrobs)

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
}
    


multa <- function(A, X, ndim, nrens, iblkmax) {
    for (ia in seq(1, ndim, iblkmax)) {
        ib = min(ia + iblkmax - 1, ndim)
        v(1:(ib - ia + 1), 1:nrens) = A(ia:ib, 1:nrens)
        A[ia, 1] <- v[1,1] %*% X[1,1]
    }

}