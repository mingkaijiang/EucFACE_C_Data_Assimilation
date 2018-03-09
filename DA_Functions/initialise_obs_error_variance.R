initialise_obs_error_variance <- function(s, err_var) {
    
    # right now all the error variance are fake values
    # will need to update them based on EucFACE data
    for (k in seq(100, 1000, 100)) {
        err_var[s$POS_RA,k] = 0.05  # %
        err_var[s$POS_AF,k] = 0.05  # %
        err_var[s$POS_AW,k] = 0.05  # %
        err_var[s$POS_AR,k] = 0.05  # %
        err_var[s$POS_LF,k] = 0.1  # g C m-2 d-1
        err_var[s$POS_LW,k] = 0.1  # g C m-2 d-1
        err_var[s$POS_LR,k] = 0.1  # g C m-2 d-1
        err_var[s$POS_CF,k] = 0.05  # %
        err_var[s$POS_CW,k] = 0.05  # %
        err_var[s$POS_CR,k] = 0.05  # %
        err_var[s$POS_RH1,k] = 0.05 # %
        err_var[s$POS_RH2,k] = 0.05 # %
        err_var[s$POS_D,k] = 0.05   # %
        err_var[s$POS_CL,k] = 0.05  # %
        err_var[s$POS_CS,k] = 0.05  # %
        err_var[s$POS_GPP,k] = 0.05 # %
    }

    return(err_var)
}
