initialise_error_variance <- function(s, err_var) {
    
    # default error variances for the state vector elements
    err_var[s$POS_RA] = 0.2  # %
    err_var[s$POS_AF] = 0.2  # %
    err_var[s$POS_AW] = 0.2  # %
    err_var[s$POS_AR] = 0.2  # %
    err_var[s$POS_LF] = 0.5  # g C m-2 d-1
    err_var[s$POS_LW] = 0.5  # g C m-2 d-1
    err_var[s$POS_LR] = 0.5  # g C m-2 d-1
    err_var[s$POS_CF] = 0.2  # %
    err_var[s$POS_CW] = 0.2  # %
    err_var[s$POS_CR] = 0.2  # %
    err_var[s$POS_RH1] = 0.2 # %
    err_var[s$POS_RH2] = 0.2 # %
    err_var[s$POS_D] = 0.2   # %
    err_var[s$POS_CL] = 0.2  # %
    err_var[s$POS_CS] = 0.2  # %
    err_var[s$POS_GPP] = 0.2 # %

    return(err_var)
}
