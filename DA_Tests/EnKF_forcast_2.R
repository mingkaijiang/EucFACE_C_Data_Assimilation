####----  source all the input variables and functions ----####
source("include/prepare.R")
library(mvtnorm)
source("DA_Tests/enkf_function_2.R")

####----  Set up the model stuffs ----####
### read in met data
met = read.csv("Martin_Python/data/dalec_drivers.OREGON.no_obs.csv")
ndays <- nrow(met)

### initialise parameters and structures
p <- setup_p_initial_conditions()
s <- setup_s_initial_conditions()

##TRUE STATES

# Generate the true states T(x,t).
dt <- 1                                                  # time step size used in the explicit method
tt <- ndays                                              # upper bound of time window (that is, max number of time steps)
st <- seq(1, tt, by=dt)                                  # lower time bounds of the integration intervals
ns <- length(st)                                         # number of integrations
nc <- s$ndims                                            # number of nodes
A <- matrix(0, ncol=ns, nrow=nc)                         # specify matrix for states
dx <- 1                                                  # spacing between nodes
A[,1] <- c(200, rep(100, nc-2), 150) + rnorm(nc, 0, 1)   # initial conditions (at t=0)

#Specify noises (variances)
wk <- 0.001                                              # process noise
vk <- 1                                                  # measurement noise

#State transition matrix
A_tmp <- matrix(0, ncol=nc, nrow=nc)

#Entries of the interior nodes in the transition matrix
for(i in 2:(nc-1)) {
    A_tmp[i, (i-1):(i+1)] <- c(1, -2, 1)
}

#Control input matrix
B <- matrix(0, ncol=2, nrow=nc)
B[1, 1] <- B[16, 2] <- 1                                 

#Heat source
timeHS <- 400 #time steps during which the heat sources act on the bar
hs <- c(0, c(rep(1, timeHS), rep(0, length(st)-timeHS)))

#Simulate true states
for (i in 2:ns) {
    #states
    A[, i] <- A[, i-1] + (A_tmp%*%A[, i-1])*dt +
        B%*%c(2*hs[i], .5*hs[i]) + rnorm(nc, 0, sqrt(wk))
}

A <- t(as.matrix(ensembleDF[,c("RA", "AF", "AW", "AR", "LF", "LW", "LR",
                                "CF", "CW", "CR", "RH1", "RH2", "D", "CL", "CS", "GPP")]))

##MEASUREMENTS

#We will create measurements for the first 510 time steps.
#More specifically, we assume that we take measurements from t=.01 to t=5.1.
#Again, let k be the time step index and dt the size of the time steps.
#In that way, we may compute the time (t) as k*dt. Thus, for k=1 and dt=1/100
#the corresponding time is 1*(1/100)=.01

#Extract the temperature states for t=.01 to t=5.1.
#Note that x[1, ] represents the temperature states at t=0. Thus, for t=.01
#to t=5.1 we need to select x[2:511, ]
A <- A[,2:511]
#Select temperatures at nodes 10, 20, ..., 90
np <- length(c(1:16)) #number of nodes
dataEx1 <- A[1:16,]

#Add measurement noise to the selected states
for (i in 1:510) {
    dataEx1[,i] <- dataEx1[,i] + rnorm(np, 0, sqrt(vk))
}


#Of these 510 steps, we will use the first 500 for fitting the EnKF
#and the last 10 steps for forecasting.
dataPred1 <- dataEx1[,501:510]
dataEx1 <- dataEx1[,1:500]

##FORECASTING WITH THE ENSEMBLE KALMAN FILTER (EnKF)
ex1 <- list(m0=A[,1], # initial state estimates
            #error covariances of the initial state estimates:
            #this matrix reflects the uncertainty in our initial state estimates
            #you may change the values in this matrix and see how they influence
            #the behavior of the Kalman filter
            C0=diag(rep(1e+1, nc)),
            #measurement noise
            V=diag(rep(vk, nc)),
            #process noise
            W=diag(rep(wk, nc)))

#Specify the state transition function:
#WARNING: always use arguments x and k when specifying the GGfunction
GGfunction <- function (A, k){
    A + (A_tmp%*%A)*dt + B%*%c(2*hs[k+1], .5*hs[k+1])}

#Specify the observation/measurement function:
#WARNING: always use arguments x and k when specifying the FFfunction
FFfunction <- function (A, k){
    A[1:16]}


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
plot(1:nc, x[500,], type="l", col=c(gray(level=.5)),
     ylim=range(c(x[500,], enkf1$m[501,], enkf2$m[501,], enkf3$m[501,])),
     xlab="Node, i", ylab="Temperature (i)", main="t=5")
lines(1:nc, enkf1$m[501,], lty=2, col="blue", lwd=1)
lines(1:nc, enkf2$m[501,], lty=2, col="red", lwd=1)
lines(1:nc, enkf3$m[501,], lty=2, col="darkgreen", lwd=1)
legend("topright", lty=c(1, 2, 2, 2),
       col=c(gray(level=.5), "blue", "red", "darkgreen"),
       legend=c("true state", "EnKf with 10 members",
                "EnKf with 20 members", "EnKf with 100 members"),
       bty="n", y.intersp=1.2, cex=.7)

#Error plot
e1 <- sapply(1:(nrow(x)-10), function (i) mean((x[i,]-enkf1$m[i+1,])^2))
e2 <- sapply(1:(nrow(x)-10), function (i) mean((x[i,]-enkf2$m[i+1,])^2))
e3 <- sapply(1:(nrow(x)-10), function (i) mean((x[i,]-enkf3$m[i+1,])^2))
rangeError <- range(cbind(e1, e2, e3))
plot(1:(nrow(x)-10), e1, type="l", lty=1, col="blue", log="y",
     ylim=rangeError, ylab="Mean Squared Error", xlab="Time index, k")
lines(1:(nrow(x)-10), e2, lty=1, col="red")
lines(1:(nrow(x)-10), e3, lty=1, col="darkgreen")
legend("topright", lty=c(2, 2, 2),
       col=c("blue", "red", "darkgreen"),
       legend=c("EnKf with 10 members",
                "EnKf with 20 members", "EnKf with 100 members"),
       bty="n", y.intersp=1.2, cex=.7)



#Compute forecasts for 10 time steps ahead
#(If this forecasting function generates an error, then decrease nAhead
#from 10 to, for instance, 5)
fc1 <- EnKFforecast(filterData=enkf1, nAhead=10)
fc2 <- EnKFforecast(filterData=enkf2, nAhead=10)
fc3 <- EnKFforecast(filterData=enkf3, nAhead=10)

##FUTURE SYSTEM STATES
#95% confidence intervals for forecasts

#EnKF with 10 members: future system states 10 steps ahead
seFoA <- sqrt(fc1$R[[10]])
ciFoA <- fc1$a[10, ] + qnorm(.05/2)*seFoA%o%c(1, -1)

#Plot future system states including confidence intervals
plot(1:nc, x[510, ], type="l", lty=2, col="red",
     xlab="Node, i", ylab="Temperature (i)",
     main="EnKF with 10 members")
lines(1:nc, fc1$a[10, ], col="blue", lwd=1)
lines(1:nc, ciFoA[,1], col="blue", lty=2, lwd=1)
lines(1:nc, ciFoA[,2], col="blue", lty=2, lwd=1)
legend("topright", lty=c(2, 1), 
       col=c("red", "blue"),
       legend=c("true state", "future state"),
       bty="n", y.intersp=1.2, cex=.7)



#EnKF with 20 members: future system states 10 steps ahead
seFoA <- sqrt(fc2$R[[10]])
ciFoA <- fc2$a[10, ] + qnorm(.05/2)*seFoA%o%c(1, -1)

#Plot future system states including confidence intervals
plot(1:nc, x[510, ], type="l", lty=2, col="red",
     xlab="Node, i", ylab="Temperature (i)",
     main="EnKF with 20 members")
lines(1:nc, fc2$a[10, ], col="blue", lwd=1)
lines(1:nc, ciFoA[,1], col="blue", lty=2, lwd=1)
lines(1:nc, ciFoA[,2], col="blue", lty=2, lwd=1)
legend("topright", lty=c(2, 1), 
       col=c("red", "blue"),
       legend=c("true state", "future state"),
       bty="n", y.intersp=1.2, cex=.7)



#EnKF with 100 members: future system states 10 steps ahead
seFoA <- sqrt(fc3$R[[10]])
ciFoA <- fc3$a[10, ] + qnorm(.05/2)*seFoA%o%c(1, -1)

#Plot future system states including confidence intervals
plot(1:nc, x[510, ], type="l", lty=2, col="red",
     xlab="Node, i", ylab="Temperature (i)",
     main="EnKF with 100 members")
lines(1:nc, fc3$a[10, ], col="blue", lwd=1)
lines(1:nc, ciFoA[,1], col="blue", lty=2, lwd=1)
lines(1:nc, ciFoA[,2], col="blue", lty=2, lwd=1)
legend("topright", lty=c(2, 1), 
       col=c("red", "blue"),
       legend=c("true state", "future state"),
       bty="n", y.intersp=1.2, cex=.7)



##FUTURE OBSERVATIONS
#95% confidence intervals for forecasts

#EnKF with 10 members: future observations
seFoX <- sqrt(fc1$Q[[10]])
ciFoX <- fc1$f[10, ] + qnorm(.05/2)*seFoX%o%c(1, -1)

#Plot prediction intervals of future observations
plot(1:nc, x[510, ], type="l", lty=2, col=gray(level=.7),
     xlab="Node, i", ylab="Temperature (i)",
     main="EnKF with 10 members")
points(seq(10,90,10), dataPred1[10,], col="red", cex=.5)
arrows(seq(10,90,10), ciFoX[,1], seq(10,90,10) ,ciFoX[,2],
       code=3, length=0.1, angle=90, col=gray(level=.4))
legend("topright", lty=c(2, NA, 1), pch=c(NA, 1, NA),
       col=c(gray(level=.7), "red", gray(level=.7)),
       legend=c("true state", "measurements", "prediction interval"),
       bty="n", y.intersp=1.2, cex=.7)


#EnKF with 20 members: future observations
seFoX <- sqrt(fc2$Q[[10]])
ciFoX <- fc2$f[10, ] + qnorm(.05/2)*seFoX%o%c(1, -1)

#Plot prediction intervals of future observations
plot(1:nc, x[510, ], type="l", lty=2, col=gray(level=.7),
     xlab="Node, i", ylab="Temperature (i)",
     main="EnKF with 20 members")
points(seq(10,90,10), dataPred1[10,], col="red", cex=.5)
arrows(seq(10,90,10), ciFoX[,1], seq(10,90,10) ,ciFoX[,2],
       code=3, length=0.1, angle=90, col=gray(level=.4))
legend("topright", lty=c(2, NA, 1), pch=c(NA, 1, NA),
       col=c(gray(level=.7), "red", gray(level=.7)),
       legend=c("true state", "measurements", "prediction interval"),
       bty="n", y.intersp=1.2, cex=.7)



#EnKF with 100 members: future observations
seFoX <- sqrt(fc3$Q[[10]])
ciFoX <- fc3$f[10, ] + qnorm(.05/2)*seFoX%o%c(1, -1)

#Plot prediction intervals of future observations
plot(1:nc, x[510, ], type="l", lty=2, col=gray(level=.7),
     xlab="Node, i", ylab="Temperature (i)",
     main="EnKF with 100 members")
points(seq(10,90,10), dataPred1[10,], col="red", cex=.5)
arrows(seq(10,90,10), ciFoX[,1], seq(10,90,10) ,ciFoX[,2],
       code=3, length=0.1, angle=90, col=gray(level=.4))
legend("topright", lty=c(2, NA, 1), pch=c(NA, 1, NA),
       col=c(gray(level=.7), "red", gray(level=.7)),
       legend=c("true state", "measurements", "prediction interval"),
       bty="n", y.intersp=1.2, cex=.7)





#Future system states and observations for Node 10 separately.

#EnKF with 100 members: future system states for node 10
seFoA <- unlist(lapply(fc3$R, function (i) sqrt(i[10])))
ciFoA <- fc3$a[1:10, 10] + qnorm(.05/2)*seFoA%o%c(1, -1)

#Plot future system states including 95% confidence intervals for Node 10
plot(1:10, x[501:510, 10], type="l", lty=2, col="red",
     ylim=range(ciFoA),
     xlab="Step ahead, k", ylab="Temperature (k)",
     main="EnKF with 100 members")
lines(1:10, fc3$a[1:10, 10], col="blue", lwd=1)
lines(1:10, ciFoA[,1], col="blue", lty=2, lwd=1)
lines(1:10, ciFoA[,2], col="blue", lty=2, lwd=1)
legend("topright", lty=c(2, 1, 2), 
       col=c("red", "blue", "blue"),
       legend=c("true state", "future state", "confidence interval"),
       bty="n", y.intersp=1.2, cex=.7)


#EnKF with 100 members: future observations for node 10
seFoX <- unlist(lapply(fc3$Q, function (i) sqrt(i[1])))
ciFoX <- fc3$f[1:10, 1] + qnorm(.05/2)*seFoX%o%c(1, -1)

#Plot 95% prediction intervals for future observations of Node 10
plot(1:10, x[501:510, 10], type="l", lty=2, col=gray(level=.7),
     ylim=range(ciFoX),
     xlab="Step ahead, k", ylab="Temperature (k)",
     main="EnKF with 100 members")
points(1:10, dataPred1[1:10, 1], col="red", cex=.5)
lines(1:10, fc3$f[1:10, 1], col="blue", lwd=1)
lines(1:10, ciFoX[,1], col="blue", lty=2, lwd=1)
lines(1:10, ciFoX[,2], col="blue", lty=2, lwd=1)
legend("topleft", lty=c(2, NA, 1, 2), pch=c(NA, 1, NA, NA),
       col=c(gray(level=.7), "red", "blue", "blue"),
       legend=c("true state", "measurement", "predicted", "prediction interval"),
       bty="n", y.intersp=1.2, cex=.7)

