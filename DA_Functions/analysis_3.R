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
    
    ## Generate the H matrix: gsl_matrix_set(H, j, k, (o + j)->obsop[k]);
    for (j in 1:nrobs) {
        for (k in 1:ndims) {
            H[j,k] <- obs[k,j]   # This should be a pointer operator function
        }
    }

    ## Compute HA
    HA <- H %*% A
    
    ## Calculate observation uncertainty (defined as gamma in Eversen 2003)
    for (j in 1:nrobs) {
        mean.value <- 0
        sum.value <- 0
        for (k in 1:s$nrens) {
            if ((o + i)->err_type_obs[j] == 0) {
                tmp = gsl_ran_gaussian(rand_num_gen, (o + i)->stdev)
            } else {
                tmp = gsl_ran_gaussian(rand_num_gen, (o + i)->stdev * fabs((o + i)->val))
            }
            sum.value = sum.value + tmp
            E[j,k] <- tmp
        }
        mean.value = sum.value / s$nrens
        E_mean[j] <- mean.value
    }
    
    ## Compute the innovations, D', where D' = D - HA. eqn 53 evenson 2003 
    for (j in 1:nrobs) {
        mean.value = 0
        sum.value  = 0
        for (k in 1:s$nrens) {
            ## Add the observation uncertainty to observation (eqn 48 Evenson 2003) before taking HA away from it */
            tmp = E[j,k] + (o + j)->val - HA[j,k]
            sum.value = sum.value + tmp
            D[j,k] <- tmp
        }
        mean.value = sum.value / s$nrens
        D_mean[j] <- mean.value
    }
    
    ### Calculate HA' by taking mean_HA from HA. 
    ### H must be linear!  
    ### Store in matrix S as Evenson does. 
    ### Point 5 section 5.2 page 356 Evenson 2003. 

    ## First figure out mean value of HA 
    for (j in 1:nrobs) {
        mean.value = 0
        sum.value  = 0
        for (k in 1:s$nrens) {
            sum.value = sum.value + HA[j,k] 
        }
        mean.value = sum.value / s$nrens
        HA_mean[j] <- mean.value
    }
    
    ## HA' = HA - mean_HA 
    for (j in 1:nrobs) {
        for (k in 1:s$nrens) {
            tmp = HA[j,k] - HA_mean[j]
            S[j,k] <- tmp
        }
    }
    
    ## Compute HA' + E (evenson page 356) , Ha' stored in S...so ES = S + E 
    ES <- S + E
    
    
    ### Compute SVD of HA'+E store in U, eqn 59 and pg 357 evenson example 2003
    ## LAPACK bits for SVD
    lwork = 2 * max(3 * s$nrens + nrobs, 5 * s$nrens)
    
    
    
    
    
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
    
