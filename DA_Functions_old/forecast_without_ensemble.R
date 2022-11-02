forecast_without_ensemble <- function(s, p, met, i, A, err_var, err_type, ens_var, q) {
    
    A_tmp = matrix(0, s$ndims, ndays)
    A_mean = rep(0, s$ndims)
    
    # generate model prediction
    lai = pmax(0.1, A[s$POS_CF,i]  / p$sla)
    gpp = acm(met, p, lai, i)
    
    A_tmp[s$POS_GPP,i] <- gpp
    A_tmp[s$POS_RA,i] <- gpp * p$t2
    A_tmp[s$POS_AF,i] <- gpp * p$t3 * (1. - p$t2)
    A_tmp[s$POS_AR,i] <- gpp * p$t4 * (1. - p$t2)
    A_tmp[s$POS_AW,i] <- gpp * (1. - p$t3 - p$t4) * (1. - p$t2)
    A_tmp[s$POS_LF,i] <- A[s$POS_CF,i] * p$t5
    A_tmp[s$POS_LW,i] <- A[s$POS_CW,i] * p$t6
    A_tmp[s$POS_LR,i] <- A[s$POS_CR,i] * p$t7
    A_tmp[s$POS_RH1,i] <- exp(0.0693 * met$temp[i]) * A[s$POS_CL,i] * p$t8
    A_tmp[s$POS_RH2,i] <- exp(0.0693 * met$temp[i]) * A[s$POS_CS,i] * p$t9
    A_tmp[s$POS_D,i] <- exp(0.0693 * met$temp[i]) * A[s$POS_CL,i] * p$t1
    A_tmp[s$POS_CF,i] <- A[s$POS_CF,i] + A[s$POS_AF,i] - A[s$POS_LF,i]
    A_tmp[s$POS_CW,i] <- A[s$POS_CW,i] + A[s$POS_AW,i] - A[s$POS_LW,i]
    A_tmp[s$POS_CR,i] <- A[s$POS_CR,i] + A[s$POS_AR,i] - A[s$POS_LR,i]
    A_tmp[s$POS_CL,i] <- A[s$POS_CL,i] + A[s$POS_LF,i] + A[s$POS_LR,i] - A[s$POS_RH1,i] - A[s$POS_D,i]
    A_tmp[s$POS_CS,i] <- A[s$POS_CS,i] + A[s$POS_D,i] + A[s$POS_LW,i] - A[s$POS_RH2,i]
    
    A[,i] <- A_tmp[,i]

    # Calculate ensemble average
    # for (j in 1:s$ndims) {
    #     A_mean[j] <- sum(A[j,]) / s$nrens
    # }

    # add the noise model into the ensemble
    for (j in 1:s$ndims) {
        A[j,i] <- A[j,i] + sqrt(p$delta_t) * p$rho * sqrt(ens_var[j]) * q[j,i]
        #for (k in 1:s$nrens) {
        #    A[j,k] <- A[j,k] + sqrt(p$delta_t) * p$rho * sqrt(ens_var[j]) * q[j,k]
        #}
    }
        
    # simulate the time evolution of model errors
    for (j in 1:s$ndims) {
        q_previous_time_step <- q[j, i]
        q[j,i] <- p$alpha * q_previous_time_step + sqrt(1.0 - p$alpha^2) * rnorm(1, 0.0, 1.0)
        #for (k in 1:s$nrens) {
        #    q_previous_time_step <- q[j,k]
        #    q[j,k] <- p$alpha * q_previous_time_step + sqrt(1.0 - p$alpha^2) * rnorm(1, 0.0, 1.0)
        #}
    }

    # calculate the error variance
    #ens_var <- assign_model_errors(s, ens_var, err_var, err_type, A_mean)
    ens_var <- assign_model_errors(s, ens_var, err_var, err_type, A[,i])
    
    return(list(A=A, ens_var=ens_var, q=q))
}
