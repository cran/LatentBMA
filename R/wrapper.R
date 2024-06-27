#' @name ULLGM_BMA
#'
#' @title Bayesian Model Averaging for Poisson Log-Normal and Binomial Logistic-Normal Regression Models
#'
#' @description \code{ULLGM_BMA} estimates Bayesian regression models using either a Poisson log-normal (PLN) or binomial logistic-normal (BiL) regression framework. It accounts for model uncertainty via Bayesian model averaging.
#'
#' @usage
#' ULLGM_BMA(X,
#'           y,
#'           model    = "PLN",
#'           gprior   = "BRIC",
#'           nsave    = 10000,
#'           nburn    = 2000,
#'           Ni       = NULL,
#'           m        = NULL,
#'           verbose  = TRUE)
#'
#' @param X A n x p design matrix where n is the number of observations and p is the number of explanatory variables.
#' @param y A n x 1 response vector. For PLN and BiL models, this is a count response.
#' @param model Indicates the model to be estimated. Options are \code{"PLN"} for the Poisson log-normal model and \code{"BiL"} for the binomial logistic-normal model. Default is \code{"PLN"}.
#' @param gprior Specifies the g-prior to be used. Options under fixed g are \code{"BRIC"} (g = max(n, p^2)), \code{"UIP"} (g = n), \code{"RIC"} (g = p^2), \code{"SQRT-N"} (g = sqrt(n)). Options under random g are the hyper-g and hyper-g/n priors of Liang et al. (2008) (\code{"hyper-g(a=3)"}, \code{"hyper-g(a=4)"}, \code{"hyper-g/n(a=3)"}, \code{"hyper-g/n(a=4)"}), the prior of Zellner & Siow (1980) (\code{"zellnersiow"}), and a Beta(0.5, 0.5) prior on g/(1+g) (\code{"horseshoe"}). Default is \code{"BRIC"}.
#' @param nsave The number of saved posterior samples. Defaults to 10,000.
#' @param nburn The number of initial burn-in samples. Defaults to 2,000.
#' @param Ni A vector containing the number of trials for each observation when estimating a binomial logistic-normal model. Required if \code{model = "BiL"}.
#' @param m The prior expected model size as per the beta-binomial prior of Ley and Steel (2009). Defaults to \code{p/2}, representing a uniform prior on model size.
#' @param verbose Logical indicator of whether progress should be printed during estimation. Default is \code{TRUE}.
#'
#' @return A list containing the inputs and selected posterior simulation outputs, such as posterior chains for the coefficients and inclusion vectors.
#' @note All explanatory variables in \code{X} are automatically demeaned within the function. All models do automatically include an intercept term.
#' @references
#' Liang, F., Paulo, R., Molina, G., Clyde, M. A., & Berger, J. O. (2008). Mixtures of g priors for Bayesian variable selection. Journal of the American Statistical Association, 103(481), 410-423.
#'
#' Zellner, A., & Siow, A. (1980). Posterior odds ratios for selected regression hypotheses. Trabajos de estadística y de investigación operativa, 31, 585-603.
#'
#' Ley, E., & Steel, M. F. J. (2009). On the effect of prior assumptions in Bayesian model averaging with applications to growth regression. Journal of Applied Econometrics, 24(4), 651-674.
#'
#' @examples
#' # Load package
#' library(LatentBMA)
#'
#' # Example 1: Estimate a PLN model under a BRIC prior with m = p/2 using simulated data
#' # Note: Use more samples for actual analysis
#' # Note: nsave = 250 and nburn = 250 are for demonstration purposes
#' X <- matrix(rnorm(100*20), 100, 20)
#' z <- 2 + X %*% c(0.5, -0.5, rep(0, 18)) + rnorm(100, 0, sqrt(0.25))
#' y <- rpois(100, exp(z))
#' results_pln <- ULLGM_BMA(X = X, y = y, model = "PLN", nsave = 250, nburn = 250)
#'
#' # Example 2: Estimate a BiL model under a Zellner-Siow prior with m = 5 using simulated data
#' # Note: Use more samples for actual analysis
#' # Note: nsave = 250 and nburn = 250 are for demonstration purposes
#' X  <- matrix(rnorm(100*20), 100, 20)
#' Ni <- rep(50, 100)
#' z  <- 2 + X %*% c(0.5, -0.5, rep(0, 18)) + rnorm(100, 0, sqrt(0.25))
#' y  <- rbinom(100, Ni, 1 / (1 + exp(-z)))
#' results_bil <- ULLGM_BMA(X = X, y = y, Ni = Ni, model = "BiL", nsave = 250, nburn = 250,
#'                          m = 5, gprior = "zellnersiow")
#'
#' @import progress mnormt stats
#' @export
ULLGM_BMA     = function(X,
                         y,
                         model          = "PLN",
                         gprior         = "BRIC",
                         nsave          = 10000,
                         nburn          = 2000,
                         Ni             = NULL,
                         m              = NULL,
                         verbose        = TRUE){


# ---
# --- BASIC CHECKS
# ---

X = as.matrix(X)
y = as.matrix(y)

if(is.null(Ni) & model == "BiL"){stop("Number of trials Ni not specified.")}
if(any(Ni == 1)){warning("Single number of trials (Ni=1) leads to non-identifiability of model.")}
if(qr(X)$rank<ncol(X)){warning("Provided design matrix is not of full rank.")}
if(any(colSums(X)==nrow(X))){warning("It seems an intercept term has been included in X. Remove the intercept to avoid multicollinearity.")}

# ---
# --- SET UP PRIORS, LIKELIHOODS, GRADIENTS
# ---

if(gprior=="zellnersiow"){    g_prior = g_zellnersiow}
if(gprior=="horseshoe"){      g_prior = g_horseshoe}
if(gprior=="hyper-g(a=3)"){   g_prior = g_hyperg3}
if(gprior=="hyper-g(a=4)"){   g_prior = g_hyperg4}
if(gprior=="hyper-g/n(a=3)"){ g_prior = g_hypergn3}
if(gprior=="hyper-g/n(a=4)"){ g_prior = g_hypergn4}
if(gprior=="RIC"){            g_prior = "RIC"}
if(gprior=="BRIC"){           g_prior = "BRIC"}
if(gprior=="UIP"){            g_prior = "UIP"}
if(gprior=="SQRT-N"){         g_prior = "SQRT-N"}

if(model == "PLN"){

  LogLikelihood  = ll_poisson
  LogLikGradient = gr_poisson

}

if(model == "BiL"){

  LogLikelihood  = ll_binomial
  LogLikGradient = gr_binomial

}


# ---
# --- RUN MCMC
# ---

mcmc_results = ULLGM_BMA_MCMC(X,
                              y,
                              Ni             = Ni,
                              nsave          = nsave,
                              nburn          = nburn,
                              m              = m,
                              LogLikelihood  = LogLikelihood,
                              LogLikGradient = LogLikGradient,
                              gprior         = g_prior,
                              model          = model,
                              verbose        = verbose)


# ---
# --- OUTPUT
# ---

return(mcmc_results)



}
