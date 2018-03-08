analysis <- function(A, s, obs) {
    
    ### The standard analysis eqn: A = A + Pe H^T(H Pe H^T + Re)^-1 (D - H A) can
    ### be reformed using D' = D - HA, Pe = A'(A')^T, Re = YY^T
    ### (where Y symbolises Gamma) such that it
    ### becomes A = A + A' A'^T H^T(HA' A'^T H^T + YY^T)^-1 D'
    
    sig_sum = 0.0
    sig_sum1 = 0.0
    
    # Minimum of nrobs and nrens
    nrmin = min(s$nrobs+1, s$nrens)
    
    I = matrix(0, s$nrens,s$nrens)
    U = matrix(0, s$nrobs,nrmin)
    D = matrix(0, s$nrobs,s$nrens)
    S = matrix(0, s$nrobs,nrmin)
    E = matrix(0, s$nrobs,nrmin)
    H = matrix(0, s$nrobs,s$ndims)
    HA = matrix(0, s$nrobs,s$nrens)
    ES = matrix(0, s$nrobs,s$nrens)
    X1 = matrix(0, nrmin,s$nrens)
    X2 = matrix(0, nrmin,s$nrens)
    X3 = matrix(0, s$nrobs,s$nrens)
    X4 = matrix(0, s$nrens,s$nrens)
    Reps = matrix(0, s$ndims,s$nrobs)
    A_tmp = matrix(0, s$ndims,s$nrens)
    A_dash = matrix(0, s$ndims,s$nrens)
    
    sig = seq(0, nrmin)
    D_mean = seq(0, s$nrobs)
    E_mean = seq(0, s$nrobs)
    S_mean = seq(0, s$nrobs)
    HA_mean = seq(0, s$nrobs)
    A_mean = seq(0, s$nrobs)
}
    
