setup_stochastic_model_error <- function(p) {
### Set up stochastic model error according to eqn 42 - Evenson, 2003. This
### ensures that the variance growth over time becomes independent of alpha and
### delta_t (as long as the dynamical model is linear).
    
    # number of timesteps per time unit
    n = 1.0
    
    num = (1.0 - p$alpha)^2
    den = n - 2.0 * p$alpha * n * p$alpha^2 + (2.0 * p$alpha)^(n + 1.0)
    
    return(sqrt(1.0 / p$delta_t * num / den))
}
