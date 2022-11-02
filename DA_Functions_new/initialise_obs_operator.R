initialise_obs_operator <- function(obs) {
    
    #obsop <- rep(0, s$ndims)
    #
    #### This follows example of William
    #obsop[s$POS_RA] = 0.0
    #obsop[s$POS_AF] = 0.0 
    #obsop[s$POS_AW] = 0.0 
    #obsop[s$POS_AR] = 0.0 
    #obsop[s$POS_LF] = 0.0 
    #obsop[s$POS_LW] = 0.0 
    #obsop[s$POS_LR] = 0.0 
    #obsop[s$POS_CF] = 1.0 
    #obsop[s$POS_CW] = 0.0 
    #obsop[s$POS_CR] = 0.0 
    #obsop[s$POS_RH1] = 0.0
    #obsop[s$POS_RH2] = 0.0
    #obsop[s$POS_D] = 0.0 
    #obsop[s$POS_CL] = 0.0 
    #obsop[s$POS_CS] = 0.0 
    #obsop[s$POS_GPP] = 0.0 
    
    
    obsop <- matrix(0, nrow(obs), s$ndims)
    obsop = ifelse(!is.na(obs), 1, 0)
    
    return(obsop)
}
