dump_output <- function(s, A) {
    
    x = sum(A[c.POS_CF,])
    x2 = sum(A[c.POS_CF,]^2)
    
    ensemble_member_avg = x / s$nrens
    ensemble_member_stdev_error = sqrt((x2 - (x^2) / s$nrens) / s$nrens)
    
    print(ensemble_member_avg, ensemble_member_stdev_error)
}
