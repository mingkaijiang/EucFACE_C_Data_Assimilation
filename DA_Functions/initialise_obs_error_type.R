initialise_obs_error_type <- function(s, err_type) {
    ## Define types of error 
    ## % as 1
    ## absolute values as 0
    
    err_type[s$POS_RA] = 1
    err_type[s$POS_AF] = 1
    err_type[s$POS_AW] = 1
    err_type[s$POS_AR] = 1
    err_type[s$POS_LF] = 1
    err_type[s$POS_LW] = 1
    err_type[s$POS_LR] = 1
    err_type[s$POS_CF] = 1
    err_type[s$POS_CW] = 1
    err_type[s$POS_CR] = 1
    err_type[s$POS_RH1] = 1
    err_type[s$POS_RH2] = 1
    err_type[s$POS_D] = 1
    err_type[s$POS_CL] = 1
    err_type[s$POS_CS] = 1
    err_type[s$POS_GPP] = 1
    
    return(err_type)
}
