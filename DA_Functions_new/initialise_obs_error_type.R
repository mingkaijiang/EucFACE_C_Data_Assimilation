initialise_obs_error_type <- function(s, err_type, obs) {
    ## Define types of error 
    ## % as 1
    ## absolute values as 0
    
    #err_type[s$POS_RA] = 1
    #err_type[s$POS_AF] = 1
    #err_type[s$POS_AW] = 1
    #err_type[s$POS_AR] = 1
    #err_type[s$POS_LF] = 1
    #err_type[s$POS_LW] = 1
    #err_type[s$POS_LR] = 1
    #err_type[s$POS_CF] = 1
    #err_type[s$POS_CW] = 1
    #err_type[s$POS_CR] = 1
    #err_type[s$POS_RH1] = 1
    #err_type[s$POS_RH2] = 1
    #err_type[s$POS_D] = 1
    #err_type[s$POS_CL] = 1
    #err_type[s$POS_CS] = 1
    #err_type[s$POS_GPP] = 1
    
    
    err_type[s$POS_RA,] = ifelse(!is.na(obs[,s$POS_RA]), 1, NA)   # %
    err_type[s$POS_AF,] = ifelse(!is.na(obs[,s$POS_AF]), 1, NA)  # %           # this is NA in Williams
    err_type[s$POS_AW,] = ifelse(!is.na(obs[,s$POS_AW]), 1, NA)  # %           # this is NA in Williams
    err_type[s$POS_AR,] = ifelse(!is.na(obs[,s$POS_AR]), 1, NA)  # %           # this is NA in Williams
    err_type[s$POS_LF,] = ifelse(!is.na(obs[,s$POS_LF]), 1, NA)   # %            
    err_type[s$POS_LW,] = ifelse(!is.na(obs[,s$POS_LW]), 1, NA)   # %             
    err_type[s$POS_LR,] = ifelse(!is.na(obs[,s$POS_LR]), 1, NA)   # %           # this is NA in Williams 
    err_type[s$POS_CF,] = ifelse(!is.na(obs[,s$POS_CF]), 1, NA)   # %
    err_type[s$POS_CW,] = ifelse(!is.na(obs[,s$POS_CW]), 1, NA)   # %
    err_type[s$POS_CR,] = ifelse(!is.na(obs[,s$POS_CR]), 1, NA)   # %
    err_type[s$POS_RH1,] = ifelse(!is.na(obs[,s$POS_RH1]), 1, NA) # %           # this is NA in Williams 
    err_type[s$POS_RH2,] = ifelse(!is.na(obs[,s$POS_RH2]), 1, NA) # %           # this is NA in Williams 
    err_type[s$POS_D,] = ifelse(!is.na(obs[,s$POS_D]), 1, NA)   # %           # this is NA in Williams 
    err_type[s$POS_CL,] = ifelse(!is.na(obs[,s$POS_CL]), 1, NA)    # %
    err_type[s$POS_CS,] = ifelse(!is.na(obs[,s$POS_CS]), 1, NA)    # %
    err_type[s$POS_GPP,] = ifelse(!is.na(obs[,s$POS_GPP]), 1, NA)   # %
    
    return(err_type)
}
