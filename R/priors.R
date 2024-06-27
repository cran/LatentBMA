# ---
# --- MODEL SPACE PRIOR
# ---

BetaBinomial    = function(a, b, P, P.j){

  lgamma(a + b) - (lgamma(a) + lgamma(b)) + lgamma(a + P.j) + lgamma(b + P - P.j) - lgamma(a + b + P)

}

# ---
# --- g-PRIORS
# ---

g_hyperg3       = function(g, N=NULL, b=NULL, c=NULL){-(3/2) * log(1+g)}
g_hyperg4       = function(g, N=NULL, b=NULL, c=NULL){-(4/2) * log(1+g)}
g_hypergn3      = function(g, N=NULL, b=NULL, c=NULL){-(3/2) * log(1+g/N)}
g_hypergn4      = function(g, N=NULL, b=NULL, c=NULL){-(4/2) * log(1+g/N)}
g_horseshoe     = function(g, N=NULL, b=0.5,  c=0.5){lgamma(b+c) - (lgamma(b) + lgamma(c)) +(b-1)*log(g)-(b+c)*log(g+1)}
g_zellnersiow   = function(g, N=NULL, b=NULL, c=NULL){ (0.5 * log(N/2) - lgamma(0.5)) + (- 0.5 - 1) * log(g) - ((N/2) / g) }

