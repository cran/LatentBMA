# ---
# --- MODEL SPACE PRIOR
# ---

BetaBinomial    = function(a, b, P, P.j){

  lgamma(a + b) - (lgamma(a) + lgamma(b)) + lgamma(a + P.j) + lgamma(b + P - P.j) - lgamma(a + b + P)

}

# ---
# --- g-PRIORS
# ---

g_zellnersiow    = function(g, N=NULL, P=NULL, b=NULL, c=NULL){ (0.5 * log(N/2) - lgamma(0.5)) + (- 0.5 - 1) * log(g) - ((N/2) / g) }
g_hyperg2.5      = function(g, N=NULL, P=NULL, b=1, c=2.5/2-1){lgamma(b + c) - (lgamma(b) + lgamma(c)) + (b-1) * log(g) - (b+c) * log(g+1)}
g_hyperg3        = function(g, N=NULL, P=NULL, b=1, c=3/2-1){lgamma(b + c) - (lgamma(b) + lgamma(c)) + (b-1) * log(g) - (b+c) * log(g+1)}
g_hyperg3.5      = function(g, N=NULL, P=NULL, b=1, c=3.5/2-1){lgamma(b + c) - (lgamma(b) + lgamma(c)) + (b-1) * log(g) - (b+c) * log(g+1)}
g_hyperg4        = function(g, N=NULL, P=NULL, b=1, c=4/2-1){lgamma(b + c) - (lgamma(b) + lgamma(c)) + (b-1) * log(g) - (b+c) * log(g+1)}
g_horseshoe      = function(g, N=NULL, P=NULL, b=0.5,  c=0.5){lgamma(b+c) - (lgamma(b) + lgamma(c)) +(b-1)*log(g)-(b+c)*log(g+1)}
g_benchmarkN0.01 = function(g, N=NULL, P=NULL, b=0.01,  c=0.01){lgamma(b*N+c) - (lgamma(b*N) + lgamma(c)) +(b*N-1)*log(g)-(b*N+c)*log(g+1)}
g_benchmarkN0.1  = function(g, N=NULL, P=NULL, b=0.1,  c=0.1){lgamma(b*N+c) - (lgamma(b*N) + lgamma(c)) +(b*N-1)*log(g)-(b*N+c)*log(g+1)}
g_benchmarkN1    = function(g, N=NULL, P=NULL, b=1,  c=1){lgamma(b*N+c) - (lgamma(b*N) + lgamma(c)) +(b*N-1)*log(g)-(b*N+c)*log(g+1)}
g_benchmarkP0.01 = function(g, N=NULL, P=NULL, b=0.01,  c=0.01){lgamma(b*(P^2)+c) - (lgamma(b*(P^2)) + lgamma(c)) +(b*(P^2)-1)*log(g)-(b*(P^2)+c)*log(g+1)}
g_benchmarkP0.1  = function(g, N=NULL, P=NULL, b=0.1,  c=0.1){lgamma(b*(P^2)+c) - (lgamma(b*(P^2)) + lgamma(c)) +(b*(P^2)-1)*log(g)-(b*(P^2)+c)*log(g+1)}
g_benchmarkP1    = function(g, N=NULL, P=NULL, b=1,  c=1){lgamma(b*(P^2)+c) - (lgamma(b*(P^2)) + lgamma(c)) +(b*(P^2)-1)*log(g)-(b*(P^2)+c)*log(g+1)}
g_hypergn2.5     = function(g, N=NULL, P=NULL, b=NULL, c=NULL){log(2.5-2)-log(2*N)- (2.5/2) * log(1+g/N)}
g_hypergn3       = function(g, N=NULL, P=NULL, b=NULL, c=NULL){log(3-2)-log(2*N)- (3/2) * log(1+g/N)}
g_hypergn3.5     = function(g, N=NULL, P=NULL, b=NULL, c=NULL){log(3.5-2)-log(2*N)- (3.5/2) * log(1+g/N)}
g_hypergn4       = function(g, N=NULL, P=NULL, b=NULL, c=NULL){log(4-2)-log(2*N)- (4/2) * log(1+g/N)}
