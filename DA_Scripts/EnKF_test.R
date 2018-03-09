
### Clear the console
rm(list=ls(all=TRUE))

### source all the input variables and functions
source("include/prepare.R")

### prepare observation data matrix
obs <- c()

### read in met data
met = read.csv("Martin_Python/data/dalec_drivers.OREGON.no_obs.csv")

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

### set up storage df to store the simulation output, with uncertainties
ndays <- nrow(met)
ensembleDF <- matrix(0, nrow=ndays, ncol=(1+s$ndims*2))
ensembleDF <- as.data.frame(ensembleDF)
colnames(ensembleDF) <- c("Days", "RA", "AF", "AW", "AR", "LF", "LW", "LR",
                          "RH1", "RH2", "D", "GPP", "CF", "CW", "CR", "CL", "CS",
                          'RA_STDEV', "AF_STDEV", "AW_STDEV", "AR_STDEV", 
                          "LF_STDEV", "LW_STDEV", "LR_STDEV", "RH1_STDEV", 
                          "RH2_STDEV", "D_STDEV", "GPP_STDEV", "CF_STDEV", 
                          "CW_STDEV", "CR_STDEV", "CL_STDEV", "CS_STDEV")
ensembleDF$Days <- c(1:ndays)

### Run the model
for (i in 1:ndays) {
    
    ## Forecast model
    out <- forecast(s, p, met, i, A, err_var,
                    err_type, ens_var, q)
    
    ## break the list of the forecast model
    A <- out$A
    ens_var <- out$ens_var
    q <- out$q
    
    ## Recalcualte model forecast where observations are avaliable
    # if (s$nrobs > 0) {
    #     analysis(A, c, obs)
    # }
    
    # Save output
    ensembleDF[i, 2:(s$ndims*2+1)] <- dump_output(s, A)
    
}

### Plotting
ggplot(ensembleDF) +
    geom_ribbon(aes(x = Days, ymin=CF-CF_STDEV, 
                  ymax=CF+CF_STDEV, fill="st.dev"), alpha=1) +
    geom_line(aes(y = CF, x=Days, color = "black")) +
    scale_colour_manual("", values = "blue") +
    scale_fill_manual("", values = "grey")
