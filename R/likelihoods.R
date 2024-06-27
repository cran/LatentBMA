# ---
# --- LOG LIKELIHOODS
# ---

ll_poisson  = function(y, z, Ni = NULL){dpois(y, exp(z), log=T)}
ll_binomial = function(y, z, Ni = NULL){dbinom(y, Ni, 1/(1+exp(-z)), log=T)}
