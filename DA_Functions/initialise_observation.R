initialise_observation <- function(p, s, A) {
    
    # right now all these observations are fake values
    # will need to update them based on EucFACE data
    for (k in seq(100, 1000, 100)) {
        A[s$POS_RA,k] <- 1.0 + rnorm(1, 0.0, 0.1 * 1.0)
        A[s$POS_AF,k] = 0.3 + rnorm(1, 0.0, 0.1 * 0.3)
        A[s$POS_AW,k] = 0.3 + rnorm(1, 0.0, 0.1 * 0.3)
        A[s$POS_AR,k] = 0.3 + rnorm(1, 0.0, 0.1 * 0.3)
        A[s$POS_LF,k] = 0.3 + rnorm(1, 0.0, 0.1 * 0.3)
        A[s$POS_LW,k] = 0.3 + rnorm(1, 0.0, 0.1 * 0.3)
        A[s$POS_LR,k] = 0.3 + rnorm(1, 0.0, 0.1 * 0.3)
        A[s$POS_RH1,k] = 0.3 + rnorm(1, 0.0, 0.1 * 0.3)
        A[s$POS_RH2,k] = 0.3 + rnorm(1, 0.0, 0.1 * 0.3)
        A[s$POS_D,k] = 0.3 + rnorm(1, 0.0, 0.1 * 0.3)
        A[s$POS_GPP,k] = 1.0 + rnorm(1, 0.0, 0.1 * 1.0)
        A[s$POS_CF,k] = p$cf0 + rnorm(1, 0.0, 0.1 * p$cf0)
        A[s$POS_CW,k] = p$cw0 + rnorm(1, 0.0, 0.1 * p$cw0)
        A[s$POS_CR,k] = p$cr0 + rnorm(1, 0.0, 0.1 * p$cr0)
        A[s$POS_CL,k] = p$cl0 + rnorm(1, 0.0, 0.1 * p$cl0)
        A[s$POS_CS,k] = p$cs0 + rnorm(1, 0.0, 0.1 * p$cs0)
    }
    
    
    return(A)
}
    
