dump_output <- function(s, A) {
    
    out <- matrix(0, ncol=s$ndims*2, nrow=1)
    out <- as.data.frame(out)

    for (j in 1:s$ndims) {
        x = sum(A[j,])
        x2 = sum(A[j,]^2)
        
        out[, j] = x / s$nrens
        out[, s$ndims+j] = sqrt((x2 - (x^2) / s$nrens) / s$nrens)
    }

    return(out)
}
