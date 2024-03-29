#### This is a simple model for assimilating the C cycle reported in Williams et al. 2005 GCB
#### Note:
#### 1. Need to add EnKF function (either developing our own or using existing R function)
#### 2. Need to modify structure to incorporate the EucFACE data (mean and uncertainty = currently error is estimated)
#### 3. Need to replace the climate forcing file with EucFACE data
#### More advanced:
#### 4. Incorporate appropriate function to calculate GPP here within the script
#### 5. Then we can do water-carbon coupling.
#### 6. What about machine learning algorithm?


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

### Setup matrix holding ensemble members
A <- matrix(0, s$ndims, s$nrens)
err_var <- rep(0, s$ndims)
err_type <- rep(0, s$ndims)
ens_var <- rep(0, s$ndims)

### initialize noise matrix 
q <- matrix(rnorm(s$ndims*s$nrens, 0.0, 1.0), s$ndims, s$nrens)

### initialize A
A <- initialise_ensemble(p, s, A)

### initialize error stuffs
err_var <- initialise_error_variance(s, err_var)
err_type <- initialise_error_type(s, err_type)

####----  Set up the observation stuffs ----####
#### Note: I think obsDF should switch rows and columns, so that variable names are headers, but maybe not.
####       For now stick with the current code. 
obsDF <- read.csv("observation/obs_cf_1.csv", header=F)
nrobs <- ncol(obsDF)-1
obs <- as.matrix(obsDF[,2:ncol(obsDF)])
obsop <- initialise_obs_operator()
nrobs <- initialise_nrobs(obs)

#### initialize measurement error variance and type
err_var_obs <- matrix(0, s$ndims, ndays)
err_type_obs <- rep(0, s$ndims)

err_var_obs <- initialise_obs_error_variance(s, err_var_obs)
err_type_obs <- initialise_obs_error_type(s, err_type_obs)


####----  set up storage df to store the simulation output, with uncertainties ----####
ensembleDF <- matrix(0, nrow=ndays, ncol=(1+s$ndims*2))
ensembleDF <- as.data.frame(ensembleDF)
colnames(ensembleDF) <- c("Days", "RA", "AF", "AW", "AR", "LF", "LW", "LR",
                          "CF", "CW", "CR", "RH1", "RH2", "D", "CL", "CS", "GPP", 
                          'RA_STDEV', "AF_STDEV", "AW_STDEV", "AR_STDEV", 
                          "LF_STDEV", "LW_STDEV", "LR_STDEV", "CF_STDEV", 
                          "CW_STDEV", "CR_STDEV", "RH1_STDEV", "RH2_STDEV", 
                          "D_STDEV", "CL_STDEV", "CS_STDEV", "GPP_STDEV")
ensembleDF$Days <- c(1:ndays)

####---- Run the model ----####
for (i in 1:ndays) {
    
    ## Forecast model
    out <- forecast(s, p, met, i, A, err_var,
                    err_type, ens_var, q)
    
    ## split the list of the forecast model
    A <- out$A                 # model prediction ensemble member
    ens_var <- out$ens_var     # model variance
    q <- out$q                 # model error
    
    A <- analysis(A, s, obs, i, 
                  err_var, err_type, 
                  err_var_obs, err_type_obs, 
                  ens_var, q)
    
    ## Save output
    ensembleDF[i, 2:(s$ndims*2+1)] <- dump_output(s, A)
    
}


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
obsDF$CF <- mean(obs[s$POS_CF,], na.rm=T)
obsDF$CF_STDEV <- sd(obs[s$POS_CF,], na.rm=T)


### plotting    
ggplot() +
    geom_ribbon(data=ensembleDF, aes(x = Days, ymin=CF-CF_STDEV, 
                  ymax=CF+CF_STDEV), fill="grey", alpha=1) +
    geom_line(data=ensembleDF, aes(y = CF, x=Days), color = "black") +
    geom_point(data=obsDF, aes(y = CF, x=Days), color="red")  +
    geom_ribbon(data=obsDF, aes(ymin=CF-CF_STDEV, ymax=CF+CF_STDEV, x=Days), 
                  fill="brown", alpha=1/10)
