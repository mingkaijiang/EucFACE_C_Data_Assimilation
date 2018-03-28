initialise_nrobs <- function(nrobs) {
    
    nrobs[s$POS_RA] = 4
    nrobs[s$POS_AF] = 4      # This is empty in Williams
    nrobs[s$POS_AW] = 4      # This is empty in Williams
    nrobs[s$POS_AR] = 4      # This is empty in Williams
    nrobs[s$POS_LF] = 18
    nrobs[s$POS_LW] = 18
    nrobs[s$POS_LR] = 1      # This is empty in Williams
    nrobs[s$POS_CF] = 4
    nrobs[s$POS_CW] = 3
    nrobs[s$POS_CR] = 2
    nrobs[s$POS_RH1] = 1     # This is empty in Williams
    nrobs[s$POS_RH2] = 1     # This is empty in Williams
    nrobs[s$POS_D] = 1       # This is empty in Williams
    nrobs[s$POS_CL] = 1
    nrobs[s$POS_CS] = 1
    nrobs[s$POS_GPP] = 1096
}
