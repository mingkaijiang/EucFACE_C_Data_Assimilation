analysis_3 <- function(A, s, p, obs, i,
                       err_var, err_type, 
                       err_var_obs, err_type_obs, 
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
    A_mean = seq(0, s$nrobs)

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
    
    # measurement error matrix
    R = matrix(0, s$nrobs, s$nrobs)
    
    HA = matrix(0, s$nrobs,s$nrens)
    ES = matrix(0, s$nrobs,s$nrens)
    X1 = matrix(0, nrmin,s$nrens)
    X2 = matrix(0, nrmin,s$nrens)
    X3 = matrix(0, s$nrobs,s$nrens)
    X4 = matrix(0, s$nrens,s$nrens)
    A_tmp = matrix(0, s$ndims,s$nrens)
    
    
    
    
    
    ## calculate A_mean: why nrobs?
    A_mean <- 

    ## assign measurement operator matrix
    #for (j in 1:s$nrobs) {               # need to create number of observations
        for (k in 1:s$ndims) {
            H[,k] <- obs[k,i]             # in the original script, this is not a obs value, i.e. obsop[x]
        }
    #}
    
    ## assign observation error matrix (nrobs * nrobs matrix)
    ## question: why nrobs not ndims?
    ##           why squared?
    ##           right now err_var_obs is a ndims by ndays matrix
    for (j in 1:s$nrobs) {
        if (err_type_obs[j] == 0) {
            R[j,j] <- (err_var_obs[j,i]) ^ 2
        } else {
            R[j,j] <- (err_var_obs[j,i] * abs(A_mean[j])) ^ 2
        } 
    }
    
    ## set up observation vector
    for (j in 1:s$nrobs) {
        d[j] <- obs[j,i]       # again, the original obs matrix is a ndims by ndays matrix, 
    }                          # where are the additional observations?
                               # which dims to use? 

    ## (H.Pf.H^T + R)          
                               # mod_err, matrix, this should be a ndims * ndims matrix
    tmp <- mod_err %*% t(H)    # ndims * nrobs
    tmp1 <- H %*% tmp          # nrobs * nrobs
    tmp1 <- tmp1 + R           # nrobs * nrobs
    
    ## compuete Pf.H^T
    xx <- mod_err %*% t(H)m    # ndims * nrobs
    
    ## Pf.H^T * (H.Pf.H^T + R)^-1   Kalman gain matrix: ndims * nrobs
    K <- xx %*% tmp1^-1        # ndims * nrobs
    
    ### model_analysis = model_forecast + K * (obs - model_forecast)
    ## H * model_sv i.e. H.Psi_f
    xxx <- H %*% A             # a vector of nrobs
    
    ## d - H.Psi_f. 
    d <- d - xxx               # a vector of nrobs
    
    ## Psi_f + K * (d - H.Psi_f)
    A <- A + K %*% d # a vector of ndims
    
    ### analysis of error covariance: Pa = Pf - K.H.Pf
    ## H.Pf
    crap <- H %*% mod_err      # nrobs * ndims
    
    ## K * H.Pf
    crap2 <- H %*% crap        # ndims * ndims
    
    ## Pf - K.H.Pf
    mod_err <- mod_err - crap2 # ndims * ndims
    
    return(list(A=A, ens_var=mod_err))
}
    
