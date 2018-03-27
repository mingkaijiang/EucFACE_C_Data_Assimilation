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

##TRUE STATES
# Generate the true states T(x,t).
A <- matrix(0, ncol=ndays, nrow=s$ndims)                         # specify matrix for states
A[,1] <- c(200, rep(100, s$ndims-2), 150) + rnorm(s$ndims, 0, 1)   # initial conditions (at t=0)

#Specify noises (variances)
wk <- 0.001                                              # process noise
vk <- 1                                                  # measurement noise

#State transition matrix
A_tmp <- matrix(0, ncol=s$ndims, nrow=s$ndims)

#Entries of the interior nodes in the transition matrix
for(i in 2:(s$ndims-1)) {
    A_tmp[i, (i-1):(i+1)] <- c(1, -2, 1)
}

#Simulate true states
for (i in 2:ndays) {
    #states
    A[, i] <- A[,i-1] + (A_tmp%*%A[,i-1]) + rnorm(s$ndims, 0, sqrt(wk))
}

#A <- t(as.matrix(ensembleDF[,c("RA", "AF", "AW", "AR", "LF", "LW", "LR",
#                                "CF", "CW", "CR", "RH1", "RH2", "D", "CL", "CS", "GPP")]))

##MEASUREMENTS
np <- length(c(1,3,5,7,9)) #number of nodes
dataEx1 <- A[c(1,3,5,7,9),]

#Add measurement noise to the selected states
for (i in 1:ndays) {
    dataEx1[,i] <- dataEx1[,i] + rnorm(np, 0, sqrt(vk))
}


##FORECASTING WITH THE ENSEMBLE KALMAN FILTER (EnKF)
ex1 <- list(m0=A[,1], # initial state estimates
            #error covariances of the initial state estimates:
            #this matrix reflects the uncertainty in our initial state estimates
            #you may change the values in this matrix and see how they influence
            #the behavior of the Kalman filter
            C0=diag(rep(0.1, s$ndims)),
            #measurement noise
            V=diag(rep(vk, np)),
            #process noise
            W=diag(rep(wk, s$ndims)))

#Specify the state transition function:
#WARNING: always use arguments x and k when specifying the GGfunction
GGfunction <- function (A, k){
    A + (A_tmp%*%A)}

#Specify the observation/measurement function:
#WARNING: always use arguments x and k when specifying the FFfunction
FFfunction <- function (A, k){
    A[c(1,3,5,7,9)]}


##Compute the filtered (a posteriori) state estimates with the EnKF
#and employ 10 ensemble members in the EnKF 
enkf1 <- enkf_function(y=dataEx1, mod=ex1, size=10,
              GGfunction=GGfunction, FFfunction=FFfunction)

#As a comparison, increase the size of the ensemble to 20
enkf2 <- enkf_function(y=dataEx1, mod=ex1, size=20,
              GGfunction=GGfunction, FFfunction=FFfunction)

#And, finally, an EnKF with 100 ensemble members
enkf3 <- enkf_function(y=dataEx1, mod=ex1, size=100,
              GGfunction=GGfunction, FFfunction=FFfunction)



#Plot the filtered state estimates at t=5
plot(1:s$ndims, A[,500], type="l", col=c(gray(level=.5)),
     ylim=range(c(A[,500], enkf1$m[,501], enkf2$m[,501], enkf3$m[,501])),
     xlab="Node, i", ylab="Temperature (i)", main="t=5")
lines(1:s$ndims, enkf1$m[,501], lty=2, col="blue", lwd=1)
lines(1:s$ndims, enkf2$m[,501], lty=2, col="red", lwd=1)
lines(1:s$ndims, enkf3$m[,501], lty=2, col="darkgreen", lwd=1)
legend("topright", lty=c(1, 2, 2, 2),
       col=c(gray(level=.5), "blue", "red", "darkgreen"),
       legend=c("true state", "EnKf with 10 members",
                "EnKf with 20 members", "EnKf with 100 members"),
       bty="n", y.intersp=1.2, cex=.7)

#Error plot
e1 <- sapply(1:(nrow(A)-10), function (i) mean((A[,i]-enkf1$m[,i+1])^2))
e2 <- sapply(1:(nrow(A)-10), function (i) mean((A[,i]-enkf2$m[,i+1])^2))
e3 <- sapply(1:(nrow(A)-10), function (i) mean((A[,i]-enkf3$m[,i+1])^2))
rangeError <- range(cbind(e1, e2, e3))
plot(1:(nrow(A)-10), e1, type="l", lty=1, col="blue", log="y",
     ylim=rangeError, ylab="Mean Squared Error", xlab="Time index, k")
lines(1:(nrow(A)-10), e2, lty=1, col="red")
lines(1:(nrow(A)-10), e3, lty=1, col="darkgreen")
legend("topright", lty=c(2, 2, 2),
       col=c("blue", "red", "darkgreen"),
       legend=c("EnKf with 10 members",
                "EnKf with 20 members", "EnKf with 100 members"),
       bty="n", y.intersp=1.2, cex=.7)


