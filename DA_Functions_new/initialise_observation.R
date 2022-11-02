initialise_observation <- function(p, s, B) {
    
    # right now all these observations are fake values
    # will need to update them based on EucFACE data
    #for (k in seq(100, 1000, 100)) {
        B[s$POS_RA,] <- 1.0 + rnorm(1, 0.0, 0.1 * 1.0)
        B[s$POS_AF,] = 0.3 + rnorm(1, 0.0, 0.1 * 0.3)
        B[s$POS_AW,] = 0.3 + rnorm(1, 0.0, 0.1 * 0.3)
        B[s$POS_AR,] = 0.3 + rnorm(1, 0.0, 0.1 * 0.3)
        B[s$POS_LF,] = 0.3 + rnorm(1, 0.0, 0.1 * 0.3)
        B[s$POS_LW,] = 0.3 + rnorm(1, 0.0, 0.1 * 0.3)
        B[s$POS_LR,] = 0.3 + rnorm(1, 0.0, 0.1 * 0.3)
        B[s$POS_RH1,] = 0.3 + rnorm(1, 0.0, 0.1 * 0.3)
        B[s$POS_RH2,] = 0.3 + rnorm(1, 0.0, 0.1 * 0.3)
        B[s$POS_D,] = 0.3 + rnorm(1, 0.0, 0.1 * 0.3)
        B[s$POS_GPP,] = 1.0 + rnorm(1, 0.0, 0.1 * 1.0)
        B[s$POS_CF,] = p$cf0 + rnorm(1, 0.0, 0.1 * p$cf0)
        B[s$POS_CW,] = p$cw0 + rnorm(1, 0.0, 0.1 * p$cw0)
        B[s$POS_CR,] = p$cr0 + rnorm(1, 0.0, 0.1 * p$cr0)
        B[s$POS_CL,] = p$cl0 + rnorm(1, 0.0, 0.1 * p$cl0)
        B[s$POS_CS,] = p$cs0 + rnorm(1, 0.0, 0.1 * p$cs0)
    #}
    
    
    return(B)
}
    
