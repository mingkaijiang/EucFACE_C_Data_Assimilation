initialise_obs_error_variance <- function(s, err_var_obs, obs) {
    

        #err_var_obs[s$POS_RA,] = 0.5   # %
        #err_var_obs[s$POS_AF,] = 0.05  # %           # this is NA in Williams
        #err_var_obs[s$POS_AW,] = 0.05  # %           # this is NA in Williams
        #err_var_obs[s$POS_AR,] = 0.05  # %           # this is NA in Williams
        #err_var_obs[s$POS_LF,] = 0.2   # %            
        #err_var_obs[s$POS_LW,] = 0.2   # %             
        #err_var_obs[s$POS_LR,] = 0.2   # %           # this is NA in Williams 
        #err_var_obs[s$POS_CF,] = 0.1   # %
        #err_var_obs[s$POS_CW,] = 0.1   # %
        #err_var_obs[s$POS_CR,] = 0.3   # %
        #err_var_obs[s$POS_RH1,] = 0.05 # %           # this is NA in Williams 
        #err_var_obs[s$POS_RH2,] = 0.05 # %           # this is NA in Williams 
        #err_var_obs[s$POS_D,] = 0.05   # %           # this is NA in Williams 
        #err_var_obs[s$POS_CL,] = 0.3   # %
        #err_var_obs[s$POS_CS,] = 0.3   # %
        #err_var_obs[s$POS_GPP,] = 0.3  # %
        
        
        err_var_obs[s$POS_RA,] = ifelse(!is.na(obs[,s$POS_RA]), 0.5, NA)   # %
        err_var_obs[s$POS_AF,] = ifelse(!is.na(obs[,s$POS_AF]), 0.05, NA)  # %           # this is NA in Williams
        err_var_obs[s$POS_AW,] = ifelse(!is.na(obs[,s$POS_AW]), 0.05, NA)  # %           # this is NA in Williams
        err_var_obs[s$POS_AR,] = ifelse(!is.na(obs[,s$POS_AR]), 0.05, NA)  # %           # this is NA in Williams
        err_var_obs[s$POS_LF,] = ifelse(!is.na(obs[,s$POS_LF]), 0.2, NA)   # %            
        err_var_obs[s$POS_LW,] = ifelse(!is.na(obs[,s$POS_LW]), 0.2, NA)   # %             
        err_var_obs[s$POS_LR,] = ifelse(!is.na(obs[,s$POS_LR]), 0.2, NA)   # %           # this is NA in Williams 
        err_var_obs[s$POS_CF,] = ifelse(!is.na(obs[,s$POS_CF]), 0.1, NA)   # %
        err_var_obs[s$POS_CW,] = ifelse(!is.na(obs[,s$POS_CW]), 0.1, NA)   # %
        err_var_obs[s$POS_CR,] = ifelse(!is.na(obs[,s$POS_CR]), 0.3, NA)   # %
        err_var_obs[s$POS_RH1,] = ifelse(!is.na(obs[,s$POS_RH1]), 0.05, NA) # %           # this is NA in Williams 
        err_var_obs[s$POS_RH2,] = ifelse(!is.na(obs[,s$POS_RH2]), 0.05, NA) # %           # this is NA in Williams 
        err_var_obs[s$POS_D,] = ifelse(!is.na(obs[,s$POS_D]), 0.05, NA)   # %           # this is NA in Williams 
        err_var_obs[s$POS_CL,] = ifelse(!is.na(obs[,s$POS_CL]), 0.3, NA)    # %
        err_var_obs[s$POS_CS,] = ifelse(!is.na(obs[,s$POS_CS]), 0.3, NA)    # %
        err_var_obs[s$POS_GPP,] = ifelse(!is.na(obs[,s$POS_GPP]), 0.3, NA)   # %


    return(err_var_obs)
}
