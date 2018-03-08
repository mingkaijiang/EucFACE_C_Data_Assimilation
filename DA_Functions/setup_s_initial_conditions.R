setup_s_initial_conditions <- function() {
    
    s <- c()
    
    s$nrobs = 0 # Number of observations
    s$ndims = 16  # Dimension of model state
    s$nrens = 200 # Number of ensemble members
    s$max_params = 15
    s$seed = 0
    
    s$POS_RA = 0
    s$POS_AF = 1
    s$POS_AW = 2
    s$POS_AR = 3
    s$POS_LF = 4
    s$POS_LW = 5
    s$POS_LR = 6
    s$POS_CF = 7
    s$POS_CW = 8
    s$POS_CR = 9
    s$POS_RH1 = 10
    s$POS_RH2 = 11
    s$POS_D = 12
    s$POS_CL = 13
    s$POS_CS = 14
    s$POS_GPP = 15
    
    return(s)
    
}