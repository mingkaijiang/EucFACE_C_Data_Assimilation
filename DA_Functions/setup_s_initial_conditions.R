setup_s_initial_conditions <- function() {
    
    s <- c()
    
    ### Number of observations
    s$nrobs = 0   
    
    ### Dimension of model state
    s$ndims = 16  
    
    ### Number of ensemble members
    s$nrens = 200 
    
    s$max_params = 15
    s$seed = 0
    
    ### Position of each variable in the data matrix
    ## Autotrophic respiration
    s$POS_RA = 0
    
    ## Allocation coefficients
    s$POS_AF = 1
    s$POS_AW = 2
    s$POS_AR = 3
    
    ## Litter production rates
    s$POS_LF = 4
    s$POS_LW = 5
    s$POS_LR = 6
    
    ## Plant C mass
    s$POS_CF = 7
    s$POS_CW = 8
    s$POS_CR = 9
    
    ## Rh of fresh litter
    s$POS_RH1 = 10
    
    ## Rh of SOM/WD
    s$POS_RH2 = 11
    
    ## Dcomposition rate of litter
    s$POS_D = 12
    
    ##  Litter mass
    s$POS_CL = 13
    
    ## Soil C mass
    s$POS_CS = 14
    
    ## GPP
    s$POS_GPP = 15
    
    return(s)
    
}