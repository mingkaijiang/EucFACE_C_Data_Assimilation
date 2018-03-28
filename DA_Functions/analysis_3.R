analysis_3 <- function(A, s, p, obs, i,
                       err_var, err_type, 
                       err_var_obs, err_type_obs, 
                       ens_var, q,
                       ens_var_obs, q_obs) {
    
    ### The standard analysis eqn: A = A + Pe H^T(H Pe H^T + Re)^-1 (D - H A)    (eq 48)
    ### can be reformed using D' = D - HA, Pe = A'(A')^T, Re = YY^T
    ### (where Y symbolises Gamma) such that it
    ### becomes A = A + A' A'^T H^T(HA' A'^T H^T + YY^T)^-1 D'
    
    ## Define matrix and vector
    nrmin = min(nrobs, s$nrens)         # Minimum of nrobs and nrens
    nrsigma = 0
    sigsum = 0
    sigsum1 = 0
    
    I = diag(1, s$nrens,s$nrens)        # identity matrix
    U = matrix(0, nrobs,nrmin)          # local variable
    D = matrix(0, nrobs,s$nrens)        # matrix holding innovations
    S = matrix(0, nrobs,s$nrens)        # matrix holding product of HA'
    E = matrix(0, nrobs,s$nrens)        # Matrix holding observation uncertainty
    H = matrix(0, nrobs,s$ndims)        # Matrix holding observation operator
    HA = matrix(0, nrobs,s$nrens)      
    ES = matrix(0, nrobs,s$nrens)       # Matrix holding product of HA' + E
    X1 = matrix(0, nrmin,nrobs)         # local variable
    X2 = matrix(0, nrmin,s$nrens)       # local variable
    X3 = matrix(0, nrobs,s$nrens)       # local variable
    X4 = matrix(0, s$nrens,s$nrens)     # local variable
    A_tmp = matrix(0, s$ndims,s$nrens)
    A_dash = matrix(0, s$ndims,s$nrens)
    Reps = matrix(0, s$ndims,nrobs)    
    
    sig = rep(0, nrmin)
    D_mean = rep(0, nrobs)
    E_mean = rep(0, nrobs)
    S_mean = rep(0, nrobs)
    HA_mean = rep(0, nrobs)
    A_mean = rep(0, s$ndims)
    
    ## Generate the H matrix
    for (j in 1:nrobs) {
        for (k in 1:ndims) {
            H[j,k] <- 
        }
    }

    
    # calculate A_mean: the ensemble mean that stores in each column of A_mean
    A_mean <- A %*% N1
    
    # # alternative, result is the same for A_dash
    #A_mean <- rep(0, s$ndims)  
    #for (k in 1:s$ndims) {
    #    A_mean[k] <- mean(A[k,])
    #}
    
    # Ensemble perturbation matrix
    A_dash <- A - A_mean
    
    # Ensemble covariance matrix
    PE <- tcrossprod(A_dash) / (s$ndims - 1)
    
    # D: ensemble of observation perturbation
    # a matrix of nrobs * nrens
    D <- obs
    
#    # New A_mean
#    A_mean = seq(0, s$nrobs)
#    
#    ## calculate A_mean: why nrobs?
#    #A_mean <- 
#
#    ## assign measurement operator matrix
#    #for (j in 1:s$nrobs) {               # need to create number of observations
#        for (k in 1:s$ndims) {
#            H[,k] <- obs[k,i]             # in the original script, this is not a obs value, i.e. obsop[x]
#        }
#    #}
#    
#    ## assign observation error matrix (nrobs * nrobs matrix)
#    ## question: why nrobs not ndims?
#    ##           why squared?
#    ##           right now err_var_obs is a ndims by ndays matrix
#    for (j in 1:s$nrobs) {
#        if (err_type_obs[j] == 0) {
#            R[j,j] <- (err_var_obs[j,i]) ^ 2
#        } else {
#            R[j,j] <- (err_var_obs[j,i] * abs(A_mean[j])) ^ 2
#        } 
#    }
#    
#    ## set up observation vector
#    for (j in 1:s$nrobs) {
#        d[j] <- obs[j,i]       # again, the original obs matrix is a ndims by ndays matrix, 
#    }                          # where are the additional observations?
#                               # which dims to use? 
#
#    ## (H.Pf.H^T + R)          
#                               # mod_err, matrix, this should be a ndims * ndims matrix
#    tmp <- mod_err %*% t(H)    # ndims * nrobs
#    tmp1 <- H %*% tmp          # nrobs * nrobs
#    tmp1 <- tmp1 + R           # nrobs * nrobs
#    
#    ## compuete Pf.H^T
#    xx <- mod_err %*% t(H)    # ndims * nrobs
#    
#    ## Pf.H^T * (H.Pf.H^T + R)^-1   Kalman gain matrix: ndims * nrobs
#    K <- xx %*% tmp1^-1        # ndims * nrobs
#    
#    ### model_analysis = model_forecast + K * (obs - model_forecast)
#    ## H * model_sv i.e. H.Psi_f
#    xxx <- H %*% A             # a vector of nrobs
#    
#    ## d - H.Psi_f. 
#    d <- d - xxx               # a vector of nrobs
#    
#    ## Psi_f + K * (d - H.Psi_f)
#    A <- A + K %*% d # a vector of ndims
#    
#    ### analysis of error covariance: Pa = Pf - K.H.Pf
#    ## H.Pf
#    crap <- H %*% mod_err      # nrobs * ndims
#    
#    ## K * H.Pf
#    crap2 <- H %*% crap        # ndims * ndims
#    
#    ## Pf - K.H.Pf
#    mod_err <- mod_err - crap2 # ndims * ndims
    
    return(list(A=A, ens_var=mod_err))
}
    
