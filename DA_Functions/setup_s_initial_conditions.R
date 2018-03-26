setup_s_initial_conditions <- function() {
    
    s <- c()
    
    ### Number of observations
    s$nrobs = 10   
    
    ### Dimension of model state
    s$ndims = 16  
    
    ### Number of ensemble members
    s$nrens = 200 
    
    s$max_params = 15
    s$seed = 0
    
    ### Position of each variable in the data matrix
    ## Autotrophic respiration
    s$POS_RA = 1
    
    ## Allocation coefficients
    s$POS_AF = 2
    s$POS_AW = 3
    s$POS_AR = 4
    
    ## Litter production rates
    s$POS_LF = 5
    s$POS_LW = 6
    s$POS_LR = 7
    
    ## Plant C mass
    s$POS_CF = 8
    s$POS_CW = 9
    s$POS_CR = 10
    
    ## Rh of fresh litter
    s$POS_RH1 = 11
    
    ## Rh of SOM/WD
    s$POS_RH2 = 12
    
    ## Dcomposition rate of litter
    s$POS_D = 13
    
    ##  Litter mass
    s$POS_CL = 14
    
    ## Soil C mass
    s$POS_CS = 15
    
    ## GPP
    s$POS_GPP = 16
    
    return(s)
    
}