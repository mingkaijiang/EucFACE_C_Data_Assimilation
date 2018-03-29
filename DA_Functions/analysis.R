analysis <- function(A, s, obs, i,
                       err_var, err_type, 
                       err_var_obs, err_type_obs, 
                       ens_var, q) {
    
    ### The standard analysis eqn: A = A + Pe H^T(H Pe H^T + Re)^-1 (D - H A)    (eq 48)
    ### can be reformed using D' = D - HA, Pe = A'(A')^T, Re = YY^T
    ### (where Y symbolises Gamma) such that it
    ### becomes A = A + A' A'^T H^T(HA' A'^T H^T + YY^T)^-1 D'
    
    ## Define matrix and vector
    nrmin <- min(nrobs, s$nrens)         # Minimum of nrobs and nrens
    nrsigma <- 0
    sigsum <- 0
    sigsum1 <- 0
    
    I <- diag(1, s$nrens,s$nrens)        # identity matrix
    U <- matrix(0, nrobs,nrmin)          # local variable
    D <- matrix(0, nrobs,s$nrens)        # matrix holding innovations
    S <- matrix(0, nrobs,s$nrens)        # matrix holding product of HA'
    E <- matrix(0, nrobs,s$nrens)        # Matrix holding observation uncertainty
    H <- matrix(0, nrobs,s$ndims)        # Matrix holding observation operator
    HA <- matrix(0, nrobs,s$nrens)      
    ES <- matrix(0, nrobs,s$nrens)       # Matrix holding product of HA' + E
    X1 <- matrix(0, nrmin,nrobs)         # local variable
    X2 <- matrix(0, nrmin,s$nrens)       # local variable
    X3 <- matrix(0, nrobs,s$nrens)       # local variable
    X4 <- matrix(0, s$nrens,s$nrens)     # local variable
    A_dash <- matrix(0, s$ndims,s$nrens)
    Reps <- matrix(0, s$ndims,nrobs) 
    VT <- matrix(0, s$nrens, s$nrens)
    
    sig <- rep(0, nrmin)
    D_mean <- rep(0, nrobs)
    E_mean <- rep(0, nrobs)
    HA_mean <- rep(0, nrobs)
    A_mean <- rep(0, s$ndims)
    
    ## Generate the H matrix: gsl_matrix_set(H, j, k, (o + j)->obsop[k]);
    for (j in 1:nrobs) {
        for (k in 1:s$ndims) {
            H[j,k] <- obsop[k]   
        }
    }

    ## Compute HA
    HA <- H %*% A
    
    ## Calculate observation uncertainty (defined as gamma in Eversen 2003)
    for (j in 1:nrobs) {
        if (err_type_obs[8] == 0) {
            E[j,] <- rnorm(s$nrens, mean=0, sd=err_var_obs[8,i]) 
        } else {
            #E[j,] <- rnorm(s$nrens, mean=obs[8,j], sd=abs(obs[8,j] * err_var_obs[8,i])) 
            E[j,] <- rnorm(s$nrens, mean=0, sd=abs(obs[8,j] * err_var_obs[8,i])) 
        }
    }
    E_mean <- rowMeans(E, na.rm=T)
    
    
    ## Compute the innovations, D', where D' = D - HA. eqn 53 evenson 2003 
    for (j in 1:nrobs) {
        for (k in 1:s$nrens) {
            ## Add the observation uncertainty to observation (eqn 48 Evenson 2003) before taking HA away from it */
            D[j,k] <- E[j,k] + obs[8,j] - HA[j,k]
        }
    }
    D_mean <- rowMeans(D, na.rm=T)
    
    
    ### Calculate HA' by taking mean_HA from HA. 
    ### H must be linear!  
    ### Store in matrix S as Evenson does. 
    ### Point 5 section 5.2 page 356 Evenson 2003. 

    ## First figure out mean value of HA 
    HA_mean <- rowMeans(HA, na.rm=T)
    
    ## HA' = HA - mean_HA 
    for (j in 1:nrobs) {
        for (k in 1:s$nrens) {
            S[j,k] <- HA[j,k] - HA_mean[j]
        }
    }
    
    ## Compute HA' + E (evenson page 356) , Ha' stored in S...so ES = S + E 
    ES <- S + E
    
    ### Compute SVD of HA'+E store in U, eqn 59 and pg 357 evenson example 2003
    ### LAPACK bits for SVD
    #lwork = 2 * max(3 * s$nrens + nrobs, 5 * s$nrens)
    #
    ## SVD: ES = U * SIGMA * V^T
    ## ES: m * n,     i.e. nrobs * nrens
    ## SIGMA: m * n, all zero except diagnomal elements
    ## U: m * m
    ## V^T: (V transposed) is an n * n orthogonal matrix
    ## call dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
    # dgesvd_("S",        # the first min(m,n) columns of U (the left singular vectors) are returned in the array u
    #         "N",        # no rows of V^T/V^H (no right singular vectors) are computed
    #         nrobs,      # number of rows of the matrix ES
    #         s$nrens,    # number of columns in ES 
    #         es,         # array containing the m * n matrix ES, the 2nd dimension of a must be at least max (1,n)
    #         nrobs,      # the leading dimension of the array es
    #         sig,        # Contains the singular values of ES sorted so that sig (i) >= sig(i+1)
    #         U,          # m by m output
    #         nrobs,      # the leading dimension of the output array u 
    #         VT,          # should be n by n, why we have n by m here?
    #         inrens,     # the leading dimension of the output array v
    #         work_,      # workspace array, its dimension max(1, lwork)
    #         lwork,      # the dimension of the array work
    #         ierr)       # output
    # 
    out <- svd(ES, nu = min(nrobs, s$nrens), nv = min(nrobs, s$nrens))
    U <- out$u   
    VT <- out$v
    sig <- out$d

    ## Convert to eigenvalues and work out sigsum - pg 357 evenson example 2003
    for (j in 1:nrmin) {
        sig[j] <- sig[j]^2
        sigsum <- sigsum + sig[j]
    }
    
    ## Compute number of significant eigenvalues - pg 357 evenson example 2003
        for (j in 1:nrmin) {
            if ((sigsum1 / sigsum) < 0.999) {
                nrsigma <- nrsigma + 1
                sigsum1 <- sigsum1 + sig[j]
                sig[j] <- 1.0 / sig[j]
            } else {
                for (k in j:nrmin) {
                    sig[k] <- 0.0 
                }
            }
        }
    
    ## Compute X1 = sig * U - pg 357 evenson example 2003
        for (j in 1:nrobs) {
            for (k in 1:nrmin) {
                X1[k,j] <- sig[k] *  U[j,k] 
            }
        }
    
    ## Compute X2 = X1 * D - pg 357 evenson example 2003 
    X2 <- X1 %*% D

    ## Compute X3 = U * X2 - pg 357 evenson example 2003
    X3 <- U %*% X2

    ## Compute final analysis - pg 357 evenson example 2003
    ## X4 = (HA')^T * X3
    ## X5 = X4 +I
    ## A = A + A' * X5
    if ((2 * s$ndims * nrobs) > (s$nrens * (nrobs + s$ndims))) {
        ## compute X4 = (HA')^T * X3 - note S matrix is HA' 
        X4 <- t(S) %*% X3

        ## Compute X5 = X4 + I (store in X4) 
        for (k in 1:s$nrens) {
            X4[k,k] <- X4[k,k] + 1.0;
        }
        X4 <- X4 + I
        
        ## Compute A = A * X5 (note X5 stored in X4 -> see Evenson) */
        A <- A %*% X4
        
    } else {
        ## Compute representers Reps = A' * S^T
        ## Calculate the mean of the ensemble required to work out A' 
        A_mean <- rowMeans(A, na.rm=T)
        
        ## Figure out A' 
        for (j in 1:s$ndims) {
            for (k in 1:s$nrens) {
                A_dash[j,k] <- A[j,k] - A_mean[j]
            }
        }
        
        ## Compute representers Reps = A' * S^T 
        Reps <- A_dash %*% t(S)
        
        ## Compute A = A + Reps * X3 
        A <- A + Reps %*% X3
    }
    
    return(A)
}
    
