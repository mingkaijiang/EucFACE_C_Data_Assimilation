setup_p_initial_conditions <- function() {
    
    p <- c()
    
    ### dalec model values
    ## decomposition rate constant
    p$t1 = 4.41E-06
    
    ## autotrophic respiration as a fraction of GPP
    p$t2 = 0.473267
    
    ## Fraction of NPP allocated to foliage
    p$t3 = 0.314951
    
    ## fraction of NPP allocated to fineroot
    p$t4 = 0.434401
    
    ## Turnover rater of foliage
    p$t5 = 0.00266518
    
    ## Turnover rate of woody matter
    p$t6 = 2.06E-06
    
    ## Turnover rate of fineroot
    p$t7 = 2.48E-03
    
    ## Mineralization rate of fresh litter
    p$t8 = 2.28E-02
    
    ## Mineralization rate of soil organic matter and woody debris
    p$t9 = 2.65E-06
    
    ### Initial conditions
    ## Foliage C pool
    p$cf0 = 57.7049
    
    ## Wood (stems and coarse roots)
    p$cw0 = 769.863
    
    ## Fineroot
    p$cr0 = 101.955
    
    ## Fresh foliar and fineroot litter
    p$cl0 = 40.4494
    
    ## Soil organic matter plus woody debris
    p$cs0 = 9896.7
    
    # acm parameterisation
    p$a0 = 2.155;
    p$a1 = 0.0142;
    p$a2 = 217.9;
    p$a3 = 0.980;
    p$a4 = 0.155;
    p$a5 = 2.653;
    p$a6 = 4.309;
    p$a7 = 0.060;
    p$a8 = 1.062;
    p$a9 = 0.0006;
    
    # location - oregon
    p$lat = 44.4
    p$sla = 111.
    
    # timestep
    p$delta_t = 1.0
    
    # specified time decorrelation length s.
    p$tau = 1.0
    
    # The factor a should be related to the time step used, eqn 32
    # (i.e. this is zero)
    p$alpha = 1.0 - (p$delta_t / p$tau)
    
    # setup factor to ensure variance growth over time becomes independent of
    # alpha and delta_timestep (as long as the dynamical model is linear).
    p$rho = setup_stochastic_model_error(p)
 
    
    return(p)
    
}