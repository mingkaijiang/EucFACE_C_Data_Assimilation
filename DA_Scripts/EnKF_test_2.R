
####---- Clear the console ----####
rm(list=ls(all=TRUE))

####----  source all the input variables and functions ----####
source("include/prepare.R")

####----  Set up the model stuffs ----####
### read in met data
met = read.csv("Martin_Python/data/dalec_drivers.OREGON.no_obs.csv")
ndays <- nrow(met)

### initialise parameters and structures
p <- setup_p_initial_conditions()
s <- setup_s_initial_conditions()

### TRUE STATES
### Generate the true states T(x,t).
A <- matrix(0, ncol=ndays, nrow=s$ndims)                         # specify matrix for states
A <- initialise_true_state_matrix(p, s, A)


### Setup error holding arrays
err_var <- rep(0, s$ndims)
err_type <- rep(0, s$ndims)
ens_var <- rep(0, s$ndims)

### initialize error stuffs
err_var <- initialise_error_variance(s, err_var)
err_type <- initialise_error_type(s, err_type)


### initialize noise matrix 
q <- matrix(rnorm(s$ndims*s$nrens, 0.0, 1.0), s$ndims, s$nrens)


####----  Set up the observation stuffs ----####
### Create the observational matrix, for each state and day
B <- matrix(NA, s$ndims, ndays)

### initialize measurement matrix
B <- initialise_observation(p, s, B)

### initialize measurement error variance and type
err_var_obs <- matrix(0, s$ndims, ndays)
err_type_obs <- rep(0, s$ndims)
ens_var_obs <- matrix(0, s$ndims, ndays)

err_var_obs <- initialise_obs_error_variance(s, err_var_obs)
err_type_obs <- initialise_obs_error_type(s, err_type_obs)

### intialize measurement noise matrix
q_obs <- matrix(rnorm(s$ndims*ndays, 0.0, 1.0), s$ndims, ndays)

####-------- State transition matrix ---------####
A_tmp <- matrix(0, ncol=s$ndims, nrow=s$ndims)

#Entries of the interior nodes in the transition matrix
for(i in 2:(s$ndims-1)) {
    A_tmp[i, (i-1):(i+1)] <- c(1, -2, 1)
}

#Simulate true states
for (i in 2:ndays) {
    #states
    A[, i] <- A[,i-1] + (A_tmp%*%A[,i-1]) + )
}


#####----  set up storage df to store the simulation output, with uncertainties ----####
#ensembleDF <- matrix(0, nrow=ndays, ncol=(1+s$ndims*2))
#ensembleDF <- as.data.frame(ensembleDF)
#colnames(ensembleDF) <- c("Days", "RA", "AF", "AW", "AR", "LF", "LW", "LR",
#                          "CF", "CW", "CR", "RH1", "RH2", "D", "CL", "CS", "GPP", 
#                          'RA_STDEV', "AF_STDEV", "AW_STDEV", "AR_STDEV", 
#                          "LF_STDEV", "LW_STDEV", "LR_STDEV", "CF_STDEV", 
#                          "CW_STDEV", "CR_STDEV", "RH1_STDEV", "RH2_STDEV", 
#                          "D_STDEV", "CL_STDEV", "CS_STDEV", "GPP_STDEV")
#ensembleDF$Days <- c(1:ndays)
#
#####---- Run the model ----####
#for (i in 1:ndays) {
#    
#    ## Forecast model
#    out <- forecast(s, p, met, i, A, err_var,
#                    err_type, ens_var, q)
#    
#    ## split the list of the forecast model
#    A <- out$A
#    ens_var <- out$ens_var
#    q <- out$q
#    
#    # Save output
#    ensembleDF[i, 2:(s$ndims*2+1)] <- dump_output(s, A)
#    
#}
#
#### this is the predicted values
#A <- t(as.matrix(ensembleDF[,c("RA", "AF", "AW", "AR", "LF", "LW", "LR",
#                               "CF", "CW", "CR", "RH1", "RH2", "D", "CL", "CS", "GPP")]))


####----  Plotting ----####
### prepare obs matrix to df
obsDF <- matrix(NA, nrow=ndays, ncol=(1+s$ndims*2))
obsDF <- as.data.frame(obsDF)
colnames(obsDF) <- c("Days", "RA", "AF", "AW", "AR", "LF", "LW", "LR",
                          "CF", "CW", "CR", "RH1", "RH2", "D", "CL", "CS", "GPP", 
                          'RA_STDEV', "AF_STDEV", "AW_STDEV", "AR_STDEV", 
                          "LF_STDEV", "LW_STDEV", "LR_STDEV", "CF_STDEV", 
                          "CW_STDEV", "CR_STDEV", "RH1_STDEV", "RH2_STDEV", 
                          "D_STDEV", "CL_STDEV", "CS_STDEV", "GPP_STDEV")
obsDF$Days <- c(1:ndays)
obsDF[,2:(s$ndims+1)] <- as.data.frame(t(B))
obsDF[,(s$ndims+1):ncol(obsDF)] <- as.data.frame(t(B)*t(err_var_obs))


### plotting    
ggplot() +
    geom_ribbon(data=ensembleDF, aes(x = Days, ymin=CF-CF_STDEV, 
                  ymax=CF+CF_STDEV), fill="grey", alpha=1) +
    geom_line(data=ensembleDF, aes(y = CF, x=Days), color = "black") +
    geom_point(data=obsDF, aes(y = CF, x=Days), color="red") +
    geom_errorbar(data=obsDF, aes(ymin=CF-CF_STDEV, ymax=CF+CF_STDEV, x=Days), 
                  width=0.1, position=position_dodge(0.05), color="brown")
