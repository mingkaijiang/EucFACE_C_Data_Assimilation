initialise_nrobs <- function(obs) {
    
    nrobs <- rep(0, dim(obs)[1])
    
    ### This function is hard-wiring number of observations
    ### in the future could be based on real observation data
    #nrobs[s$POS_RA] = sum(!is.na(obs[s$POS_RA,]))
    #nrobs[s$POS_AF] = sum(!is.na(obs[s$POS_AF,]))
    #nrobs[s$POS_AW] = sum(!is.na(obs[s$POS_AW,]))
    #nrobs[s$POS_AR] = sum(!is.na(obs[s$POS_AR,]))
    #nrobs[s$POS_LF] = sum(!is.na(obs[s$POS_LF,]))
    #nrobs[s$POS_LW] = sum(!is.na(obs[s$POS_LW,]))
    #nrobs[s$POS_LR] = sum(!is.na(obs[s$POS_LR,]))
    #nrobs[s$POS_CF] = sum(!is.na(obs[s$POS_CF,]))
    #nrobs[s$POS_CW] = sum(!is.na(obs[s$POS_CW,]))
    #nrobs[s$POS_CR] = sum(!is.na(obs[s$POS_CR,]))
    #nrobs[s$POS_RH1] = sum(!is.na(obs[s$POS_RH1,]))
    #nrobs[s$POS_RH2] = sum(!is.na(obs[s$POS_RH2,]))
    #nrobs[s$POS_D] = sum(!is.na(obs[s$POS_D,]))
    #nrobs[s$POS_CL] = sum(!is.na(obs[s$POS_CL,]))
    #nrobs[s$POS_CS] = sum(!is.na(obs[s$POS_CS,]))
    #nrobs[s$POS_GPP] = sum(!is.na(obs[s$POS_GPP,]))
    
    for (i in 1:dim(obs)[1]) {
        nrobs[i] <- sum(!is.na(obs[i,]))
    }

    
    
    return(nrobs)
}
