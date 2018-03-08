dump_output <- function(s, A) {
    
    x = sum(A[s$POS_CF,])
    x2 = sum(A[s$POS_CF,]^2)
    
    ensemble_member_avg = x / s$nrens
    ensemble_member_stdev_error = sqrt((x2 - (x^2) / s$nrens) / s$nrens)
    
    print.out <- c(ensemble_member_avg, ensemble_member_stdev_error)
    
    # print(print.out)
    
    return(print.out)
}
