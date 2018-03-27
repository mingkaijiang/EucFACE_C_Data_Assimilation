initialise_obs_error_variance <- function(s, err_var) {
    
    # right now all the error variance are fake values
    # will need to update them based on EucFACE data
#    for (k in seq(100, 1000, 100)) {
        err_var[s$POS_RA,] = 0.05  # %
        err_var[s$POS_AF,] = 0.05  # %
        err_var[s$POS_AW,] = 0.05  # %
        err_var[s$POS_AR,] = 0.05  # %
        err_var[s$POS_LF,] = 0.1  # g C m-2 d-1
        err_var[s$POS_LW,] = 0.1  # g C m-2 d-1
        err_var[s$POS_LR,] = 0.1  # g C m-2 d-1
        err_var[s$POS_CF,] = 0.05  # %
        err_var[s$POS_CW,] = 0.05  # %
        err_var[s$POS_CR,] = 0.05  # %
        err_var[s$POS_RH1,] = 0.05 # %
        err_var[s$POS_RH2,] = 0.05 # %
        err_var[s$POS_D,] = 0.05   # %
        err_var[s$POS_CL,] = 0.05  # %
        err_var[s$POS_CS,] = 0.05  # %
        err_var[s$POS_GPP,] = 0.05 # %
#    }

    return(err_var)
}
