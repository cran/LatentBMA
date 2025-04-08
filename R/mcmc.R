gammaProposal <- function(incl.curr) {
  
  P      <- length(incl.curr)      # Total number of predictors
  p_curr <- sum(incl.curr)         # Number of currently included predictors
  incl.prop <- incl.curr
  
  # Determine possible moves based on the current model:
  # Only addition if p_curr==0, only deletion if p_curr==P,
  # else each of add, swap, and delete has probability 1/3.
  if(p_curr == 0) {
    move <- 1  # Only addition is possible
  } else if(p_curr == P) {
    move <- 3  # Only deletion is possible
  } else {
    move <- sample(c(1, 2, 3), 1)  # All moves are possible with equal probability 1/3.
  }
  
  if(move == 1) {  # Addition move
    # Determine the probability for an addition move:
    # π_A(M_k) = 1 if model is empty; 1/3 otherwise.
    pi_A <- ifelse(p_curr == 0, 1, 1/3)
    
    # Choose a predictor (not in the model) uniformly.
    add_candidates <- which(incl.curr == 0)
    selector <- sample(add_candidates, 1)
    incl.prop[selector] <- 1
    p_prop <- sum(incl.prop)  # Proposed model size (p_curr + 1)
    
    # In the reverse move (deletion from the proposed model):
    # π_D(M*) = 1 if the proposed model is full; 1/3 otherwise.
    pi_D_rev <- ifelse(p_prop == P, 1, 1/3)
    
    # Compute log probabilities.
    # Forward: q(M*|M_k) = π_A(M_k) * [1/(P - p_curr)]
    # Reverse: q(M_k|M*) = π_D(M*) * [1/(p_prop)]
    Curr2Prop <- log(pi_A) - log(P - p_curr)
    Prop2Curr <- log(pi_D_rev) - log(p_prop)
    
  } else if(move == 2) {  # Swap move
    # Swap is only available if 0 < p_curr < P.
    # π_S(M_k) is 1/3 but cancels out in the correction factor.
    
    # Select one predictor to remove uniformly from the model.
    del_candidates <- which(incl.curr == 1)
    selector_del <- sample(del_candidates, 1)
    # And one predictor to add uniformly from those not in the model.
    add_candidates <- which(incl.curr == 0)
    selector_add <- sample(add_candidates, 1)
    
    incl.prop[selector_del] <- 0
    incl.prop[selector_add] <- 1
    # For swap moves, the forward and reverse moves are symmetric.
    Curr2Prop <- 0
    Prop2Curr <- 0
    
  } else if(move == 3) {  # Deletion move
    # Determine the probability for a deletion move:
    # π_D(M_k) = 1 if the model is full; 1/3 otherwise.
    pi_D <- ifelse(p_curr == P, 1, 1/3)
    
    # Choose a predictor to delete uniformly.
    del_candidates <- which(incl.curr == 1)
    selector <- sample(del_candidates, 1)
    incl.prop[selector] <- 0
    p_prop <- sum(incl.prop)  # Proposed model size (p_curr - 1)
    
    # In the reverse move (addition to the proposed model):
    # π_A(M*) = 1 if the proposed model is empty; 1/3 otherwise.
    pi_A_rev <- ifelse(p_prop == 0, 1, 1/3)
    
    # Compute log probabilities.
    # Forward: q(M*|M_k) = π_D(M_k) * [1/(p_curr)]
    # Reverse: q(M_k|M*) = π_A(M*) * [1/(P - p_prop)]
    Curr2Prop <- log(pi_D) - log(p_curr)
    Prop2Curr <- log(pi_A_rev) - log(P - p_prop)
  }
  
  return(list(incl.prop = incl.prop, 
              Prop2Curr = Prop2Curr, 
              Curr2Prop = Curr2Prop, 
              move = move))
}

# log marginal likelihood and OLS solutions
ML_OLS    = function(y, X, XX, g, TSS){

  X.star  = cbind(1, X)
  C       = chol(XX)
  XXinv   = chol2inv(C)
  betahat = XXinv %*% crossprod(X.star, y)

  SSR     = sum((y - X.star %*% betahat)^2)
  R2      = 1 - (SSR/TSS)
  P.j     = ncol(X)
  N       = nrow(X)
  LML     = 0.5*(N - 1 - P.j) * log(1+g) - 0.5 * (N-1) * log((TSS+g*SSR))

  return(list(LML = LML, beta = betahat, TSS=TSS, SSR=SSR, XXinv=XXinv, R2=R2, C=C))


}


# correction term for barker proposal - taken from giacomo zanella's github
CorrTermBarker = function(curr, prop, GradCurr, GradProp, zeros){

  beta1 = c(-GradProp * (curr - prop))
  beta2 = c(-GradCurr * (prop - curr))

  return((# compute acceptance with log_sum_exp trick for numerical stability
    -(pmax(beta1,zeros)+log1p(exp(-abs(beta1))))+
      (pmax(beta2,zeros)+log1p(exp(-abs(beta2))))
  ))
}


# main MCMC function
ULLGM_BMA_MCMC     = function(X,
                              y,
                              Ni             = NULL,
                              nsave          = 10000,
                              nburn          = 10000,
                              m              = NULL,
                              LogLikelihood  = NULL,
                              LogLikGradient = NULL,
                              gprior         = NULL,
                              model          = NULL,
                              verbose        = TRUE){

  # demean explanatory variables
  X           = apply(X, 2, function(x) x - mean(x))

  # gibbs parameters
  ntot        = nsave + nburn

  # dimensions
  P           = max(0, ncol(X))             # number of explanatory variables
  N           = length(y)                   # total number of observations

  # beta binomial model size prior: P/2 by default
  m           = ifelse(is.null(m), P/2, m)

  # if g prior not fixed, trigger updating of g (initialize at g=N)
  if(is.function(gprior)){randomg=TRUE; g = N} else{randomg=FALSE}


  # fixed g-priors
  if(!randomg){

    if(gprior=="RIC"){    g    = P^2}
    if(gprior=="UIP"){    g    = N}
    if(gprior=="BRIC"){   g    = max(P^2, N)}
    if(gprior=="SQRT-N"){ g    = sqrt(N)}

  }


  # starting values
  z                             = rnorm(N)
  gamma                         = sample(c(0,1), P, replace=T)
  X.curr                        = X[, gamma==1,drop=F]
  beta                          = rnorm(sum(gamma) + 1)
  sig2                          = 0.1
  theta                         = numeric(length=P+1)
  theta[c(T, gamma==1)]         = beta
  theta[c(F, gamma==0)]         = 0
  fit.mean                      = cbind(1, X.curr) %*% beta
  mz  = mean(z)
  TSS = crossprod(z-mz)
  
  # for adaptive barker proposal
  tau                           = matrix(0.1, N)
  accept                        = matrix(0, N)

  # for adaptive sampling of g
  tau.g                         = 0.1

  # storage
  gamma.store                   = array(NA, c(nsave, P))
  theta.store                   = array(NA, c(nsave, P + 1))
  sig2.store                    = array(NA, c(nsave, 1))
  g.store                       = array(NA, c(nsave, 1))
  ML.store                      = array(NA, c(nsave, 1))

  # some useful precomputations
  ones                          = rep(1, N)
  zeros                         = rep(0, N)
  XX                            = crossprod(cbind(1, X))

  # progress bar
  if(verbose){

    pb <- progress_bar$new(
      format = " sampling... [:bar] :percent eta: :eta",
      total = ntot, clear = FALSE, width= 60)

  }

  # take time
  starttime = Sys.time()

  for(irep in 1:ntot){

   

    # ---
    # --- between model step
    # ---

    # compute marginal likelihood and ols quantities at current model
    # note this needs to be recomputed in each iteration when g is random!
    curr.ix           = c(1, 1 + which(gamma==1))
    comp.curr         = ML_OLS(y = z, X = X.curr, XX = XX[curr.ix, curr.ix], g = g, TSS = TSS)
    ML                = comp.curr$LML
    beta.ols          = comp.curr$beta
    TSS               = comp.curr$TSS
    SSR               = comp.curr$SSR
    XXinv             = comp.curr$XXinv
    R2                = comp.curr$R2
    C                 = comp.curr$C

    # get proposal model (add-swap-delete move)
    GP                = gammaProposal(gamma)
    gamma.prop        = GP$incl.prop
    X.prop            = X[, gamma.prop==1, drop=F]

    # Check if X'X is invertible (full rank condition)
    XtX = XX[c(T, gamma.prop==1), c(T, gamma.prop==1), drop=F]
    
    rk  = qr(XtX)$rank
    #tolerance = 1e-50
    while (rk != ncol(XtX)) {
      
      # X'X is not invertible, propose a new gamma
      GP                = gammaProposal(gamma)
      gamma.prop        = GP$incl.prop
      X.prop            = X[, gamma.prop==1, drop=F]
      XtX               = XX[c(T, gamma.prop==1), c(T, gamma.prop==1), drop=F]
      rk                = qr(XtX)$rank
      
    }

    # compute relevant posterior quantities and marginal likelihood for proposed model
    prop.ix           = c(1, 1+which(gamma.prop==1))
    comp.prop         = ML_OLS(y = z, X = X.prop, XX = XX[prop.ix, prop.ix], g = g, TSS=TSS)
    ML.prop           = comp.prop$LML
    beta.ols.prop     = comp.prop$beta
    TSS.prop          = comp.prop$TSS
    SSR.prop          = comp.prop$SSR
    XXinv.prop        = comp.prop$XXinv
    R2.prop           = comp.prop$R2
    C.prop            = comp.prop$C

    # compute loglikelihood difference
    LogLikDiff        = ML.prop - ML

    # compute correction term from ADS proposal
    CorrTerm          = GP$Prop2Curr - GP$Curr2Prop

    # compute log prior difference
    LogPriorDiff      = BetaBinomial(a = 1, b=(P-m)/m, P = P, P.j = sum(gamma.prop)) - BetaBinomial(a = 1, b = (P-m)/m, P = P, P.j = sum(gamma))

    # compute acceptance probability
    alpha             = min(1, exp(LogLikDiff + LogPriorDiff + CorrTerm))

    # accept or reject
    if(runif(1) < alpha){

      gamma     = gamma.prop
      X.curr    = X.prop
      ML        = ML.prop
      beta.ols  = beta.ols.prop
      TSS       = TSS.prop
      SSR       = SSR.prop
      XXinv     = XXinv.prop
      R2        = R2.prop
      C         = C.prop

    }



    # ---
    # --- within model step
    # ---

      # if prior on g is specified, sample g
      if(randomg){

        # proposal
        g.prop        = exp(rnorm(1, log(g), sqrt(tau.g)))
        lprior.curr   = gprior(g, N = N, P = P)
        lprior.prop   = gprior(g.prop, N=N, P = P)

        # log likelihood (difference of marginal likelihood of current model for current versus proposed g )
        P.j           = sum(gamma)
        ll.curr       = 0.5*(N - 1 - P.j) * log(1+g) - 0.5 * (N-1) * log((1+g*(1-R2)) * TSS)
        ll.prop       = 0.5*(N - 1 - P.j) * log(1+g.prop) - 0.5 * (N-1) * log((1+g.prop*(1-R2)) * TSS)

        # acceptance probability (including jacobian due to log-scale proposal)
        alpha.g            = min(1, exp(lprior.prop - lprior.curr + ll.prop - ll.curr + log(1/g) - log(1/g.prop)))

        # accept or reject
        if(runif(1) < alpha.g){g = g.prop}

        # adapt scale parameter
        l.tau.g.2          = log(tau.g^2) + (irep ^ (-0.6)) * (alpha.g - 0.234)
        tau.g              = sqrt(exp(l.tau.g.2))

      }

      # recompute shrinkage factor
      delta       = g / (1+g)

      # sample error variance
      cn     = (delta * SSR + (1 - delta) * TSS)
      sig2   = 1 / rgamma(1, 0.5 * (N - 1), 0.5 * cn)

      # sample intercept
      beta     = rnorm(1, mz, sqrt(sig2/N))
      theta[1] = beta[1]

      # sample coefficients
      if(sum(gamma)>0){

        beta                       = c(beta, mnormt::rmnorm(1, delta * beta.ols[-1], sig2 * delta * XXinv[-1,-1]))
        theta[c(F, gamma == 1)]    = beta[-1]

      }

      theta[c(F, gamma  == 0)]   = 0
      fit.mean                  = cbind(1, X.curr) %*% beta
      
      
      # ---
      # --- update latent z
      # ---
      
      # compute gradient of log posterior at current value
      Zmu        = z - fit.mean
      eZ         = exp(z)
      GradCurr   = LogLikGradient(y, z, Ni) - Zmu/sig2
      
      # barker proposal
      zi         = tau * rnorm(N)
      b          = 2 * (runif(N) < (1 / (1 + exp(-GradCurr * zi)))) - 1
      z.prop     = z + zi * b
      
      # compute gradient of log posterior at proposal value
      Zpmu        = z.prop - fit.mean
      eZp         = exp(z.prop)
      GradProp    = LogLikGradient(y, z.prop, Ni) - Zpmu/sig2
      
      # barker proposal correction term
      CorrTerm        = CorrTermBarker(z, z.prop, GradCurr, GradProp, zeros)
      
      # log prior difference
      sqrtsig2        = sqrt(sig2)
      LogPriDiff      = 0.5 * ( (Zmu / sqrtsig2)^2 - (Zpmu / sqrtsig2)^2 )
      
      # log likelihood difference
      LogLikDiff      = LogLikelihood(y, z.prop, Ni) - LogLikelihood(y, z, Ni)
      
      # acceptance probability
      zeta            = pmin(ones, exp(LogLikDiff + LogPriDiff + CorrTerm))
      
      # accept or reject
      acc             = ifelse(runif(N) < zeta, 1, 0)
      z[acc == 1]     = z.prop[acc == 1]
      
      # global scale adaption
      l.tau2          = log(tau^2) + (irep ^ (-0.6)) * (zeta - 0.57)
      tau             = sqrt(exp(l.tau2))
      
      # update some important quantities
      mz  = mean(z)
      TSS = crossprod(z-mz)


    # storage
    if(irep > nburn){

      theta.store[irep - nburn,]      = theta
      gamma.store[irep - nburn,]      = gamma
      sig2.store[irep - nburn,]       = sig2
      g.store[irep-nburn,]            = g
      ML.store[irep-nburn,]           = ML

    }

    if(verbose){

      pb$tick()

    }
  }


  # take time
  endtime = Sys.time()

  # put results in list
  results = list(y = y,
                 X = X,
                 model = model,
                 nsave = nsave,
                 nburn = nburn,
                 beta = theta.store,
                 gamma = gamma.store,
                 sig2 = sig2.store,
                 g = g.store,
                 ML = ML.store,
                 time = endtime - starttime)

  return(results)

}

