#### This is a simple model for assimilating the C cycle reported in Williams et al. 2005 GCB
#### Note: revising the script to make sure it is correctly written
#### Previous functions stored in: DA_Functions_old folder
#### Previous run script is named as: run_program_old.R



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
obsDF <- read.csv("observation/obs_old.csv", header=T)
#nrobs <- ncol(obsDF)-1
obs <- as.matrix(obsDF[,2:ncol(obsDF)])
obsop <- initialise_obs_operator(obs)
nrobs <- initialise_nrobs(obs)

#### initialize measurement error variance and type
err_var_obs <- matrix(0, s$ndims, ndays)
err_type_obs <- matrix(0, s$ndims, ndays)

err_var_obs <- initialise_obs_error_variance(s, err_var_obs, obs)
err_type_obs <- initialise_obs_error_type(s, err_type_obs, obs)


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
    
    A <- analysis(A, s, obs, i, nrobs, obsop,
                  err_var, err_type, 
                  err_var_obs, err_type_obs, 
                  ens_var, q)
    
    ## Save output
    ensembleDF[i, 2:(s$ndims*2+1)] <- dump_output(s, A)
    
}

write.csv(ensembleDF, "output/fitted_output.csv", row.names=F)



####----  Plotting ----####
### prepare obs matrix to df
pobsDF <- matrix(NA, nrow=ndays, ncol=(1+s$ndims*2))
pobsDF <- as.data.frame(pobsDF)
colnames(pobsDF) <- c("Days", "RA", "AF", "AW", "AR", "LF", "LW", "LR",
                      "CF", "CW", "CR", "RH1", "RH2", "D", "CL", "CS", "GPP", 
                      'RA_STDEV', "AF_STDEV", "AW_STDEV", "AR_STDEV", 
                      "LF_STDEV", "LW_STDEV", "LR_STDEV", "CF_STDEV", 
                      "CW_STDEV", "CR_STDEV", "RH1_STDEV", "RH2_STDEV", 
                      "D_STDEV", "CL_STDEV", "CS_STDEV", "GPP_STDEV")
pobsDF$Days <- c(1:ndays)
pobsDF[,2:17] <- obs[,1:16]


### write fitted and observed output together
subDF <- pobsDF[,1:17]
colnames(subDF) <- c("Days", "obsRA", "obsAF", "obsAW", "obsAR", "obsLF", "obsLW", "obsLR",
                     "obsCF", "obsCW", "obsCR", "obsRH1", "obsRH2", "obsD", "obsCL", "obsCS", "obsGPP")

outDF <- merge(ensembleDF, subDF, by=c("Days"))

write.csv(outDF, "output/fitted_output_with_obs.csv", row.names=F)


### plotting    
p1 <- ggplot() +
    geom_ribbon(data=ensembleDF, aes(x = Days, ymin=CF-CF_STDEV, 
                  ymax=CF+CF_STDEV), fill="grey", alpha=1) +
    geom_line(data=ensembleDF, aes(y = CF, x=Days), color = "black") +
    geom_point(data=pobsDF, aes(y = CF, x=Days), color="red"); p1  #+
    #geom_ribbon(data=pobsDF, aes(ymin=CF-CF_STDEV, ymax=CF+CF_STDEV, x=Days), 
    #              fill="brown", alpha=1/10)


p2 <- ggplot() +
    geom_ribbon(data=ensembleDF, aes(x = Days, ymin=CW-CW_STDEV, 
                                     ymax=CW+CW_STDEV), fill="grey", alpha=1) +
    geom_line(data=ensembleDF, aes(y = CW, x=Days), color = "black") +
    geom_point(data=pobsDF, aes(y = CW, x=Days), color="red") ; p2 #+
#geom_ribbon(data=pobsDF, aes(ymin=CF-CF_STDEV, ymax=CF+CF_STDEV, x=Days), 
#              fill="brown", alpha=1/10)


p3 <- ggplot() +
    geom_ribbon(data=ensembleDF, aes(x = Days, ymin=CR-CR_STDEV, 
                                     ymax=CR+CR_STDEV), fill="grey", alpha=1) +
    geom_line(data=ensembleDF, aes(y = CR, x=Days), color = "black") +
    geom_point(data=pobsDF, aes(y = CR, x=Days), color="red"); p3



p4 <- ggplot() +
    geom_ribbon(data=ensembleDF, aes(x = Days, ymin=GPP-GPP_STDEV, 
                                     ymax=GPP+GPP_STDEV), fill="grey", alpha=1) +
    geom_line(data=ensembleDF, aes(y = GPP, x=Days), color = "black") +
    geom_point(data=pobsDF, aes(y = GPP, x=Days), color="red"); p4


p5 <- ggplot() +
    geom_ribbon(data=ensembleDF, aes(x = Days, ymin=RA-RA_STDEV, 
                                     ymax=RA+RA_STDEV), fill="grey", alpha=1) +
    geom_line(data=ensembleDF, aes(y = RA, x=Days), color = "black") +
    geom_point(data=pobsDF, aes(y = RA, x=Days), color="red"); p5


p6 <- ggplot() +
    geom_ribbon(data=ensembleDF, aes(x = Days, ymin=AF-AF_STDEV, 
                                     ymax=AF+AF_STDEV), fill="grey", alpha=1) +
    geom_line(data=ensembleDF, aes(y = AF, x=Days), color = "black") +
    geom_point(data=pobsDF, aes(y = AF, x=Days), color="red"); p6



p7 <- ggplot() +
    geom_ribbon(data=ensembleDF, aes(x = Days, ymin=AW-AW_STDEV, 
                                     ymax=AW+AW_STDEV), fill="grey", alpha=1) +
    geom_line(data=ensembleDF, aes(y = AW, x=Days), color = "black") +
    geom_point(data=pobsDF, aes(y = AW, x=Days), color="red"); p7



p8 <- ggplot() +
    geom_ribbon(data=ensembleDF, aes(x = Days, ymin=AR-AR_STDEV, 
                                     ymax=AR+AR_STDEV), fill="grey", alpha=1) +
    geom_line(data=ensembleDF, aes(y = AR, x=Days), color = "black") +
    geom_point(data=pobsDF, aes(y = AR, x=Days), color="red"); p8



p9 <- ggplot() +
    geom_ribbon(data=ensembleDF, aes(x = Days, ymin=LF-LF_STDEV, 
                                     ymax=LF+LF_STDEV), fill="grey", alpha=1) +
    geom_line(data=ensembleDF, aes(y = LF, x=Days), color = "black") +
    geom_point(data=pobsDF, aes(y = LF, x=Days), color="red"); p9


p10 <- ggplot() +
    geom_ribbon(data=ensembleDF, aes(x = Days, ymin=LW-LW_STDEV, 
                                     ymax=LW+LW_STDEV), fill="grey", alpha=1) +
    geom_line(data=ensembleDF, aes(y = LW, x=Days), color = "black") +
    geom_point(data=pobsDF, aes(y = LW, x=Days), color="red"); p10


p11 <- ggplot() +
    geom_ribbon(data=ensembleDF, aes(x = Days, ymin=LR-LR_STDEV, 
                                     ymax=LR+LR_STDEV), fill="grey", alpha=1) +
    geom_line(data=ensembleDF, aes(y = LR, x=Days), color = "black") +
    geom_point(data=pobsDF, aes(y = LR, x=Days), color="red"); p10



p12 <- ggplot() +
    geom_ribbon(data=ensembleDF, aes(x = Days, ymin=RH1-RH1_STDEV, 
                                     ymax=RH1+RH1_STDEV), fill="grey", alpha=1) +
    geom_line(data=ensembleDF, aes(y = RH1, x=Days), color = "black") +
    geom_point(data=pobsDF, aes(y = RH1, x=Days), color="red"); p12



p13 <- ggplot() +
    geom_ribbon(data=ensembleDF, aes(x = Days, ymin=RH2-RH2_STDEV, 
                                     ymax=RH2+RH2_STDEV), fill="grey", alpha=1) +
    geom_line(data=ensembleDF, aes(y = RH2, x=Days), color = "black") +
    geom_point(data=pobsDF, aes(y = RH2, x=Days), color="red"); p13


p14 <- ggplot() +
    geom_ribbon(data=ensembleDF, aes(x = Days, ymin=D-D_STDEV, 
                                     ymax=D+D_STDEV), fill="grey", alpha=1) +
    geom_line(data=ensembleDF, aes(y = D, x=Days), color = "black") +
    geom_point(data=pobsDF, aes(y = D, x=Days), color="red"); p14



p15 <- ggplot() +
    geom_ribbon(data=ensembleDF, aes(x = Days, ymin=CL-CL_STDEV, 
                                     ymax=CL+CL_STDEV), fill="grey", alpha=1) +
    geom_line(data=ensembleDF, aes(y = CL, x=Days), color = "black") +
    geom_point(data=pobsDF, aes(y = CL, x=Days), color="red"); p15


p16 <- ggplot() +
    geom_ribbon(data=ensembleDF, aes(x = Days, ymin=CS-CS_STDEV, 
                                     ymax=CS+CS_STDEV), fill="grey", alpha=1) +
    geom_line(data=ensembleDF, aes(y = CS, x=Days), color = "black")  +
    geom_point(data=pobsDF, aes(y = CS, x=Days), color="red"); p16


