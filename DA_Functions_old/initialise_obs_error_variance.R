initialise_obs_error_variance <- function(s, err_var) {
    
#    for (k in seq(100, 1000, 100)) {
        err_var[s$POS_RA,] = 0.5   # %
        err_var[s$POS_AF,] = 0.05  # %           # this is NA in Williams
        err_var[s$POS_AW,] = 0.05  # %           # this is NA in Williams
        err_var[s$POS_AR,] = 0.05  # %           # this is NA in Williams
        err_var[s$POS_LF,] = 0.2   # %            
        err_var[s$POS_LW,] = 0.2   # %             
        err_var[s$POS_LR,] = 0.2   # %           # this is NA in Williams 
        err_var[s$POS_CF,] = 0.1   # %
        err_var[s$POS_CW,] = 0.1   # %
        err_var[s$POS_CR,] = 0.3   # %
        err_var[s$POS_RH1,] = 0.05 # %           # this is NA in Williams 
        err_var[s$POS_RH2,] = 0.05 # %           # this is NA in Williams 
        err_var[s$POS_D,] = 0.05   # %           # this is NA in Williams 
        err_var[s$POS_CL,] = 0.3   # %
        err_var[s$POS_CS,] = 0.3   # %
        err_var[s$POS_GPP,] = 0.3  # %
#    }

    return(err_var)
}
