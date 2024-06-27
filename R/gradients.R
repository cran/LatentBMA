# ---
# --- LOG LIKELIHOOD GRADIENTS
# ---

gr_poisson  = function(y, z, Ni = NULL){y-exp(z)}
gr_binomial = function(y, z, Ni = NULL){(y-(Ni-y)*exp(z))/(1+exp(z))}
