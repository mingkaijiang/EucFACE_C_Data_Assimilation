acm <- function(met, p, lai, i) {
    
    trange <- 0.5 * (met$maxt[i] - met$mint[i])
    gs <- abs(met$psid[i])^p$a9 / (p$a5 * met$rtot[i] + trange)
    pp <- lai * met$nit[i] / gs * p$a0 * exp(p$a7 * met$maxt[i])
    qq <- p$a2 - p$a3
    ci <- 0.5 * (met$ca[i] + qq - pp + sqrt((met$ca[i] + qq - pp)**2.0 - 4.0 * (met$ca[i] * qq - pp * p$a2)))
    e0 <- p$a6 * lai^2 / (lai**2 + p$a8)
    dec <- -23.4 * cos((360.0 * (met$doy[i]+ 10.0) / 365.0) * pi / 180.0) * pi / 180.0
    m <- tan(p$lat * pi / 180.0) * tan(dec)
    
    if (m >= 1.0) {
        dayl <- 24.0
    } else if (m <= -1.0) {
        dayl <- 0.0
    } else {
        dayl <- 24.0 * acos(-m) / pi
    }
    
    cps = e0 * met$rad[i] * gs * (met$ca[i] - ci) / (e0 * met$rad[i] + gs * (met$ca[i] - ci))
    gpp = cps * (p$a1 * dayl + p$a4)
    
    return(gpp)
}