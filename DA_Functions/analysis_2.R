analysis <- function(A, D, E, S, s) {
    
    ### The standard analysis eqn: A = A + Pe H^T(H Pe H^T + Re)^-1 (D - H A)    (eq 48)
    ### can be reformed using D' = D - HA, Pe = A'(A')^T, Re = YY^T
    ### (where Y symbolises Gamma) such that it
    ### becomes A = A + A' A'^T H^T(HA' A'^T H^T + YY^T)^-1 D'
    
    # holding arrays
    sig <- rep(0, nrmin)
    D_mean <- rep(0, s$nrobs)
    E_mean <- rep(0, s$nrobs)
    S_mean <- rep(0, s$nrobs)
    HA_mean <- rep(0, s$nrobs)
    A_mean <- rep(0, s$nrobs)
    
    # significance
    sigsum <- sum(sig)
    sigsum1 <- 0.0
    nrsigma <- 0
    
    # Minimum of nrobs - number of observations, and nrens - number of ensembles
    nrmin <- min(s$nrobs+1, s$nrens)
    
    # Matrix holding innovation
    D <- matrix(0, s$nrobs,s$nrens)
    
    # Matrix holding HA'
    # S <- matrix(0, s$nrobs,nrmin)  
    S <- matrix(0, s$nrobs,s$nrens)
    
    # Matrix holding observation perturbations
    # E <- matrix(0, s$nrobs,nrmin)
    E <- matrix(0, s$nrobs,s$nrens)
    
    # matrix holding HA' + E
    # ES <- matrix(0, s$nrobs,nrmin)
    ES <- matrix(0, s$nrobs,s$nrens)
    
    # Store SVD of HA' + E
    U <- matrix(0, s$nrobs,nrmin)
    
    # identity matrix
    I <- diag(1, s$nrens,s$nrens)
    
    H <- matrix(0, s$nrobs,s$ndims)
    HA <- matrix(0, s$nrobs,s$nrens)
    
    X1 <- matrix(0, nrmin,s$nrens)
    
    # X2 = X1 * D
    X2 <- matrix(0, nrmin,s$nrens)
    
    # X3 = U * X2
    X3 <- matrix(0, s$nrobs,s$nrens)
    
    # X4 = (HA')^T * X3
    X4 <- matrix(0, s$nrens,s$nrens)
    Reps <- matrix(0, s$ndims,s$nrobs)
    A_tmp <- matrix(0, s$ndims,s$nrens)
    A_dash <- matrix(0, s$ndims,s$nrens)

    # Compute HA' + E
    # But what is S and E? Where are the data?
    ES <- S + E
    
    # Compute SVD of HA'+E -> U and sig, using Eispack
    test <- svd(ES)
    U <- test$u         # the original U is a 10 by 11 matrix
    D <- diag(test$d)   # the original D is a 10 by 200 matrix

    # convert to eigenvalues
    for (j in c(1:nrmin)) {
        sig[j] <- sig[j] ** 2
    }
         
    # compute number of significant eigenvalues
    for (j in c(1:nrmin)) {
        if (sigsum1/sigsum < 0.999) {
            nrsigma <- nrsigma + 1
            sigsum1 <- sigsum1 + sig[j]
            sig[j] <- 1.0/sig[j]
        } else {
            sig[j:nrmin] <- 0.0
        }
    }

    # Compute X1
    for (j in c(1:s$nrobs)) {                   # j = 10
        for (k in c(1:nrmin)) {                 # k = 11
            X1[k, j] <- sig[k] * U$v[j, k]      # X1 = 11, 200, sig = 11, U$v = 10, 11
        }
    }

    # Compute X2 = X1 * D
    X2 <- X1 %*% D

    # Compute X3 = U * X2
    X3 <- U %*% X2
    
    # Compute final analysis
    if(2 * s$ndims * s$nrobs > s$nrens * (s$nrobs+s$ndims)) {
        # case nrobs = large
        # compute X4 = (HA')^T * X3
        X4 <- S %*% X3
        
        # compute X5 = X4 + I (stored in X4)
        for (i in 1:s$nrens) {
            X4[i,i] <- X4[i,i] + 1.0
        }
        
        # Compute A = A * X5
        iblkmax <- min(s$ndims, 200)
        A <- multa(A, X4, s$ndims, s$nrens, iblkmax)
        
    } else {
        # Case with nrobs = small
        # Compute representers Reps = A' * S^T
        Reps <- A %*% S
        
        # Compute A = A + Reps * X3
        A <- Reps %*% X3
    }
}
    


multa <- function(A, X, ndim, nrens, iblkmax) {
    for (ia in seq(1, ndim, iblkmax)) {
        ib <- min(ia + iblkmax - 1, ndim)
        v(1:(ib - ia + 1), 1:nrens) <- A(ia:ib, 1:nrens)
        #dgemm(’n’, ’n’, ib - ia + 1, nrens, nrens, 1.0, v(1,1), iblkmax, X(1,1), nrens,0.0, A(ia, 1), ndim)
        A[ia, 1] <- v[1,1] %*% X[1,1]
    }

}