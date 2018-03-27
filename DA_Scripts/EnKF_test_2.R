
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
q <- matrix(rnorm(s$ndims*ndays, 0.0, 1.0), s$ndims, ndays)


####-------- State transition matrix ---------####
A_tmp <- matrix(0, ncol=s$ndims, nrow=s$ndims)

#Entries of the interior nodes in the transition matrix
for(i in 2:(s$ndims-1)) {
    A_tmp[i, (i-1):(i+1)] <- c(1, -2, 1)
}

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
    out <- forecast_without_ensemble(s, p, met, i, A, err_var,
                                     err_type, ens_var, q)
    
    ## split the list of the forecast model
    A[,i] <- out$A[,i]
    ens_var <- out$ens_var
    q[,i] <- out$q[,i]
    
    # Save output
    ensembleDF[i, 2:(s$ndims*2+1)] <- dump_output(s, A)
    
}

### this is the predicted values
A <- t(as.matrix(ensembleDF[,c("RA", "AF", "AW", "AR", "LF", "LW", "LR",
                               "CF", "CW", "CR", "RH1", "RH2", "D", "CL", "CS", "GPP")]))

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

### Delete some data
B <- B[c(1,3,4,5,6,7,8,9,10,12,13,14,15,16),]
np <- nrow(B)
err_var_obs <- err_var_obs[c(1,3,4,5,6,7,8,9,10,12,13,14,15,16),]
err_type_obs <- err_type_obs[c(1,3,4,5,6,7,8,9,10,12,13,14,15,16)]
q_obs <- q_obs[c(1,3,4,5,6,7,8,9,10,12,13,14,15,16),]


##FORECASTING WITH THE ENSEMBLE KALMAN FILTER (EnKF)
ex1 <- list(m0=A[,1], # initial state estimates
            #error covariances of the initial state estimates:
            #this matrix reflects the uncertainty in our initial state estimates
            #you may change the values in this matrix and see how they influence
            #the behavior of the Kalman filter
            C0=diag(err_var),
            #measurement noise
            V=diag(q_obs[,1]),
            #process noise
            W=diag(q[,1]))


#Specify the state transition function:
#WARNING: always use arguments x and k when specifying the GGfunction
GGfunction <- function (A, k){
    A + (A_tmp%*%A)}

#Specify the observation/measurement function:
#WARNING: always use arguments x and k when specifying the FFfunction
FFfunction <- function (A, k){
    A[c(1,3,4,5,6,7,8,9,10,12,13,14,15,16)]}


##Compute the filtered (a posteriori) state estimates with the EnKF
#and employ 10 ensemble members in the EnKF 
enkf1 <- enkf_function_2(y=B, mod=ex1, size=10,
                       GGfunction=GGfunction, FFfunction=FFfunction)

#As a comparison, increase the size of the ensemble to 20
enkf2 <- enkf_function_2(y=B, mod=ex1, size=20,
                       GGfunction=GGfunction, FFfunction=FFfunction)

#And, finally, an EnKF with 100 ensemble members
enkf3 <- enkf_function_2(y=B, mod=ex1, size=100,
                       GGfunction=GGfunction, FFfunction=FFfunction)



#Plot the filtered state estimates at t=5
plot(1:ndays, A["RA",], type="l", col=c(gray(level=.5)),
     ylim=range(c(A["RA",], enkf1$m["RA",], enkf2$m["RA",], enkf3$m["RA",])),
     xlab="Node, i", ylab="Temperature (i)", main="t=5")
lines(1:s$ndims, enkf1$m[,501], lty=2, col="blue", lwd=1)
lines(1:s$ndims, enkf2$m[,501], lty=2, col="red", lwd=1)
lines(1:s$ndims, enkf3$m[,501], lty=2, col="darkgreen", lwd=1)
legend("topright", lty=c(1, 2, 2, 2),
       col=c(gray(level=.5), "blue", "red", "darkgreen"),
       legend=c("true state", "EnKf with 10 members",
                "EnKf with 20 members", "EnKf with 100 members"),
       bty="n", y.intersp=1.2, cex=.7)
