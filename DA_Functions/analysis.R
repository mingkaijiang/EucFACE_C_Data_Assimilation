analysis <- function(A, s, obs, i) {
    
    ### The standard analysis eqn: A = A + Pe H^T(H Pe H^T + Re)^-1 (D - H A)    (eq 48)
    ### can be reformed using D' = D - HA, Pe = A'(A')^T, Re = YY^T
    ### (where Y symbolises Gamma) such that it
    ### becomes A = A + A' A'^T H^T(HA' A'^T H^T + YY^T)^-1 D'
    
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
    
    # # A = U*Σ*VT for real routines
    # dgesvd (’S’,           # jobu: the first min(m, n) columns of U (the left singular vectors) are returned in the array u;
    #          ’N’,          # if jobvt = 'N', no rows of VT/VH (no right singular vectors) are computed.
    #          nrobs,        # The number of rows of the matrix A (m≥ 0).
    #          nrens,        # The number of columns in A (n≥ 0).
    #          ES,           # a(lda,*) is an array containing the m-by-n matrix A.
    #          nrobs,        # lda: The leading dimension of the array a.
    #          sig,          # s: Array, size at least max(1, min(m,n)). Contains the singular values of A sorted so that s(i) ≥ s(i+1).
    #          U,            # u(ldu,*); the second dimension of u must be at least max(1, min(m, n)) if jobu = 'S'.
    #          nrobs,        # The leading dimensions of the output arrays u and vt, respectively.
    #          V,            # The dimension of the array work.
    #          nrens,        # If jobvt = 'N'or 'O', vt is not referenced.
    #          work,         # work is a workspace array, its dimension max(1, lwork)
    #          lwork,        # The dimension of the array work.
    #          ierr)         # If info = 0, the execution is successful.
    
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

    # Compute X2 = X1 * D
    # multiplying matrices
    # dgemm(’n’,                 # Character indicating that the matrices A and B should not be transposed or conjugate transposed before multiplication.
    #        ’n’, 
    #        nrmin,              # M: A: M rows by K columns
    #        nrens,              # N: B: K rows by N columns
    #        nrobs,              # K: C: M rows by N columns
    #        1.0,                # Real value used to scale the product of matrices A and B.
    #        X1,                 # Array used to store matrix A.
    #        nrmin,              # Leading dimension of array A, or the number of elements between successive columns (for column major storage) in memory. In the case of this exercise the leading dimension is the same as the number of rows.
    #        D,                  # Array used to store matrix B.
    #        nrobs,              # Leading dimension of array B, 
    #        0.0,                # Real value used to scale matrix C.
    #        X2,                 # Array used to store matrix C
    #        nrmin)              # Leading dimension of array C
    X2 <- X1 %*% D

    # Compute X3 = U * X2
    # dgemm(’n’, ’n’, nrobs, nrens, nrmin, 1.0, U, nrobs, X2, nrmin, 0.0, X3, nrobs)
    X3 <- U %*% X2
    
    # Compute final analysis
    if(2 * ndim * nrobs > nrens * (nrobs+ndim)) {
        # dgemm(’t’, ’n’, nrens, nrens, nrobs, 1:0, S, nrobs, X3, nrobs, 0.0, X4, nrens)
        X4 <- S %*% X3
        for (i in 1:s$nrens) {
            X4[i,i] = X4[i,i] + 1.0
        }
        
        iblkmax = min(s$ndim, 200)
        multa(A, X4, ndim, nrens, iblkmax)
    } else {
        #dgemm(’n’, ’t’, ndim, nrobs, nrens, 1.0, A, ndim, S, nrobs, 0.0, Reps, ndim)
        Reps <- A %*% S
        #dgemm(’n’, ’n’, ndim, nrens, nrobs, 1.0, Reps, ndim, X3, nrobs, 1.0, A, ndim)
        A <- Reps %*% X3
    }
}
    


multa <- function(A, X, ndim, nrens, iblkmax) {
    for (ia in seq(1, ndim, iblkmax)) {
        ib = min(ia + iblkmax - 1, ndim)
        v(1:(ib - ia + 1), 1:nrens) = A(ia:ib, 1:nrens)
        #dgemm(’n’, ’n’, ib - ia + 1, nrens, nrens, 1.0, v(1,1), iblkmax, X(1,1), nrens,0.0, A(ia, 1), ndim)
        A[ia, 1] <- v[1,1] %*% X[1,1]
    }

}