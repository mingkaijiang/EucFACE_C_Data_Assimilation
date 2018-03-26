initialise_observation <- function(p, s, B) {
    
    # right now all these observations are fake values
    # will need to update them based on EucFACE data
    for (k in seq(100, 1000, 100)) {
        B[s$POS_RA,k] <- 1.0 + rnorm(1, 0.0, 0.1 * 1.0)
        B[s$POS_AF,k] = 0.3 + rnorm(1, 0.0, 0.1 * 0.3)
        B[s$POS_AW,k] = 0.3 + rnorm(1, 0.0, 0.1 * 0.3)
        B[s$POS_AR,k] = 0.3 + rnorm(1, 0.0, 0.1 * 0.3)
        B[s$POS_LF,k] = 0.3 + rnorm(1, 0.0, 0.1 * 0.3)
        B[s$POS_LW,k] = 0.3 + rnorm(1, 0.0, 0.1 * 0.3)
        B[s$POS_LR,k] = 0.3 + rnorm(1, 0.0, 0.1 * 0.3)
        B[s$POS_RH1,k] = 0.3 + rnorm(1, 0.0, 0.1 * 0.3)
        B[s$POS_RH2,k] = 0.3 + rnorm(1, 0.0, 0.1 * 0.3)
        B[s$POS_D,k] = 0.3 + rnorm(1, 0.0, 0.1 * 0.3)
        B[s$POS_GPP,k] = 1.0 + rnorm(1, 0.0, 0.1 * 1.0)
        B[s$POS_CF,k] = p$cf0 + rnorm(1, 0.0, 0.1 * p$cf0)
        B[s$POS_CW,k] = p$cw0 + rnorm(1, 0.0, 0.1 * p$cw0)
        B[s$POS_CR,k] = p$cr0 + rnorm(1, 0.0, 0.1 * p$cr0)
        B[s$POS_CL,k] = p$cl0 + rnorm(1, 0.0, 0.1 * p$cl0)
        B[s$POS_CS,k] = p$cs0 + rnorm(1, 0.0, 0.1 * p$cs0)
    }
    
    
    return(B)
}
    
