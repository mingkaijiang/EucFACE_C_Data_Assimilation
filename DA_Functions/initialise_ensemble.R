initialise_ensemble <- function(p, s, A) {
    
    ### Why initialize with the following equations? Are they mostly random numbers, and
    ### since they are evolving over time, the initial condition doesn't matter that much?
    A[s$POS_RA,] = 1.0 + rnorm(s$nrens, 0.0, 0.1 * 1.0)
    A[s$POS_AF,] = 0.3 + rnorm(s$nrens, 0.0, 0.1 * 0.3)
    A[s$POS_AW,] = 0.3 + rnorm(s$nrens, 0.0, 0.1 * 0.3)
    A[s$POS_AR,] = 0.3 + rnorm(s$nrens, 0.0, 0.1 * 0.3)
    A[s$POS_LF,] = 0.3 + rnorm(s$nrens, 0.0, 0.1 * 0.3)
    A[s$POS_LW,] = 0.3 + rnorm(s$nrens, 0.0, 0.1 * 0.3)
    A[s$POS_LR,] = 0.3 + rnorm(s$nrens, 0.0, 0.1 * 0.3)
    A[s$POS_RH1,] = 0.3 + rnorm(s$nrens, 0.0, 0.1 * 0.3)
    A[s$POS_RH2,] = 0.3 + rnorm(s$nrens, 0.0, 0.1 * 0.3)
    A[s$POS_D,] = 0.3 + rnorm(s$nrens, 0.0, 0.1 * 0.3)
    A[s$POS_GPP,] = 1.0 + rnorm(s$nrens, 0.0, 0.1 * 1.0)
    A[s$POS_CF,] = p$cf0 + rnorm(s$nrens, 0.0, 0.1 * p$cf0)
    A[s$POS_CW,] = p$cw0 + rnorm(s$nrens, 0.0, 0.1 * p$cw0)
    A[s$POS_CR,] = p$cr0 + rnorm(s$nrens, 0.0, 0.1 * p$cr0)
    A[s$POS_CL,] = p$cl0 + rnorm(s$nrens, 0.0, 0.1 * p$cl0)
    A[s$POS_CS,] = p$cs0 + rnorm(s$nrens, 0.0, 0.1 * p$cs0)
    
    return(A)
}
    
