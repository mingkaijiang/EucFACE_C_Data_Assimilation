forecast <- function(s, p, met, i, A, err_var, err_type, ens_var, q) {
    
    A_tmp = matrix(0, s$ndims, s$nrens)
    A_mean = seq(0, s$ndims)
    
    # generate model prediction
    lai = max(0.1, A[s$POS_CF,]  / p$sla)
    gpp = acm(met, p, lai, i)
    
    A_tmp[s$POS_GPP,] <- gpp
    A_tmp[s$POS_RA,] <- gpp * p$t2
    A_tmp[s$POS_AF,] <- gpp * p$t3 * (1. - p$t2)
    A_tmp[s$POS_AR,] <- gpp * p$t4 * (1. - p$t2)
    A_tmp[s$POS_AW,] <- gpp * (1. - p$t3 - p$t4) * (1. - p$t2)
    A_tmp[s$POS_LF,] <- A[s$POS_CF,] * p$t5
    A_tmp[s$POS_LW,] <- A[s$POS_CW,] * p$t6
    A_tmp[s$POS_LR,] <- A[s$POS_CR,] * p$t7
    A_tmp[s$POS_RH1,] <- exp(0.0693 * met$temp[i]) * A[s$POS_CL,] * p$t8
    A_tmp[s$POS_RH2,] <- exp(0.0693 * met$temp[i]) * A[s$POS_CS,] * p$t9
    A_tmp[s$POS_D,] <- exp(0.0693 * met$temp[i]) * A[s$POS_CL,] * p$t1
    A_tmp[s$POS_CF,] <- A[s$POS_CF,] + A[s$POS_AF,] - A[s$POS_LF,]
    A_tmp[s$POS_CW,] <- A[s$POS_CW,] + A[s$POS_AW,] - A[s$POS_LW,]
    A_tmp[s$POS_CR,] <- A[s$POS_CR,] + A[s$POS_AR,] - A[s$POS_LR,]
    A_tmp[s$POS_CL,] <- A[s$POS_CL,] + A[s$POS_LF,] + A[s$POS_LR,] - A[s$POS_RH1,] - A[s$POS_D,]
    A_tmp[s$POS_CS,] <- A[s$POS_CS,] + A[s$POS_D,] + A[s$POS_LW,] - A[s$POS_RH2,]
    
    A <- duplicate(A_tmp, shallow = FALSE)
    
    # Calculate ensemble average
    for (j in 1:s$ndims) {
        A_mean[j] <- sum(A[j,]) / s$nrens
    }

    # add the noise model into the ensemble
    for (j in 1:s$ndims) {
        for (k in 1:s$nrens) {
            A[j,k] <- A[j,k] + sqrt(p$delta_t) * p$rho * sqrt(ens_var[j]) * q[j,k]
        }
    }
        
    
    # simulate the time evolution of model errors
    for (j in 1:s$ndims) {
        for (k in 1:s$nrens) {
            q_previous_time_step <- q[j,k]
            q[j,k] <- p$alpha * q_previous_time_step + sqrt(1.0 - p$alpha^2) * rnorm(1, 0.0, 1.0)
        }
    }

    # calculate the error variance
    ens_var <- assign_model_errors(s, ens_var, err_var, err_type, A_mean)
    
    return(list(A=A, ens_var=ens_var, q=q))
}
