initialise_observation <- function(p, s, A) {
    
    # right now all these observations are fake values
    # will need to update them based on EucFACE data
    A[s$POS_RA,] = 1.0 + rnorm(ndays, 0.0, 0.1 * 1.0)
    A[s$POS_AF,] = 0.3 + rnorm(ndays, 0.0, 0.1 * 0.3)
    A[s$POS_AW,] = 0.3 + rnorm(ndays, 0.0, 0.1 * 0.3)
    A[s$POS_AR,] = 0.3 + rnorm(ndays, 0.0, 0.1 * 0.3)
    A[s$POS_LF,] = 0.3 + rnorm(ndays, 0.0, 0.1 * 0.3)
    A[s$POS_LW,] = 0.3 + rnorm(ndays, 0.0, 0.1 * 0.3)
    A[s$POS_LR,] = 0.3 + rnorm(ndays, 0.0, 0.1 * 0.3)
    A[s$POS_RH1,] = 0.3 + rnorm(ndays, 0.0, 0.1 * 0.3)
    A[s$POS_RH2,] = 0.3 + rnorm(ndays, 0.0, 0.1 * 0.3)
    A[s$POS_D,] = 0.3 + rnorm(ndays, 0.0, 0.1 * 0.3)
    A[s$POS_GPP,] = 1.0 + rnorm(ndays, 0.0, 0.1 * 1.0)
    A[s$POS_CF,] = p$cf0 + rnorm(ndays, 0.0, 0.1 * p$cf0)
    A[s$POS_CW,] = p$cw0 + rnorm(ndays, 0.0, 0.1 * p$cw0)
    A[s$POS_CR,] = p$cr0 + rnorm(ndays, 0.0, 0.1 * p$cr0)
    A[s$POS_CL,] = p$cl0 + rnorm(ndays, 0.0, 0.1 * p$cl0)
    A[s$POS_CS,] = p$cs0 + rnorm(ndays, 0.0, 0.1 * p$cs0)
    
    return(A)
}
    
