
### Clear the console
rm(list=ls(all=TRUE))

### source all the input variables and functions
source("include/prepare.R")

### prepare observation data matrix
obs <- c()

### read in met data
met = read.csv("Martin_Python/data/dalec_drivers.OREGON.no_obs.csv")

### initialise structures
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

### Run predictive model
ndays <- nrow(met)
ensembleDF <- data.frame(ndays, NA, NA)
colnames(ensembleDF) <- c("Days", "ensemble_member_avg", "ensemble_member_stdev_error")
ensembleDF$Days <- 1:ndays

for (i in 1:nrow(met)) {
    out <- forecast(s, p, met, i, A, err_var,
                    err_type, ens_var, q)
    
    A <- out$A
    ens_var <- out$ens_var
    q <- out$q
    
    # Recalcualte model forecast where observations are avaliable
    # if (s$nrobs > 0) {
    #     analysis(A, c, obs)
    # }
    
    # Print to screen outputs
    ensembleDF[i, 2:3] <- dump_output(s, A)
    
}

### Plotting
ggplot(ensembleDF) +
    geom_ribbon(aes(x = Days, ymin=ensemble_member_avg-ensemble_member_stdev_error, 
                  ymax=ensemble_member_avg+ensemble_member_stdev_error, fill="st.dev"), alpha=1) +
    geom_line(aes(y = ensemble_member_avg, x=Days, color = "black")) +
    scale_colour_manual("", values = "blue") +
    scale_fill_manual("", values = "grey")
