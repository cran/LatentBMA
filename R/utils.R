#' @name summaryBMA
#'
#' @title Summary Tables for \code{ULLGM_BMA} Estimation Results
#'
#' @description \code{summaryBMA} produces a table with estimated posterior means, standard deviations, and posterior inclusion probabilities (PIPs) for the results of a \code{ULLGM_BMA} estimation.
#'
#' @usage
#' summaryBMA(x,
#'            variable_names = NULL,
#'            digits         = 3,
#'            sort           = FALSE,
#'            type           = "pandoc")
#'
#' @param x The output object of \code{ULLGM_BMA}.
#' @param variable_names A character vector specifying the names of the columns of X.
#' @param digits Number of digits to round the table to. Defaults to 3.
#' @param sort Logical, indicating whether the table should be sorted by PIPs. Default is \code{FALSE}.
#' @param type A character string indicating the format of the table. Options are \code{'pandoc'} (default), \code{'latex'}, or \code{'html'}.
#'
#' @return Returns a 'knitr::kable' object containing the summary table.
#'
#' @author Gregor Zens
#'
#' @examples
#' # Load package
#' library(LatentBMA)
#'
#' # Example: Estimate a PLN model under a BRIC prior with m = p/2 using simulated data
#' # Note: Use more samples for actual analysis
#' # Note: nsave = 250 and nburn = 250 are for demonstration purposes
#' X <- matrix(rnorm(100*20), 100, 20)
#' z <- 2 + X %*% c(0.5, -0.5, rep(0, 18)) + rnorm(100, 0, sqrt(0.25))
#' y <- rpois(100, exp(z))
#' results_pln <- ULLGM_BMA(X = X, y = y, model = "PLN", nsave = 250, nburn = 250)
#' summaryBMA(results_pln)
#'
#' @import knitr
#' @export
summaryBMA = function(x,
                      variable_names = NULL,
                      digits = 3,
                      sort   = FALSE,
                      type   = "pandoc"
                      ){

  # specify variable names
  if(is.null(variable_names)){variable_names=colnames(x$X)}
  if(is.null(variable_names)){variable_names=paste0("x", 1:ncol(x$X))}

  # summarize results
  table.df              = data.frame(name = c("Intercept", variable_names, "sigma^2", "g", "Model Size"), beta.mean = c(colMeans(x$beta), mean(x$sig2), mean(x$g), mean(rowSums(x$gamma))), beta.sd = c(apply(x$beta,2,sd), sd(x$sig2), sd(x$g), sd(rowSums(x$gamma))), pip = c(1, colMeans(x$gamma), 1, 1, 1))
  table.df$beta.mean    = format(round(table.df$beta.mean, digits), nsmall = digits)
  table.df$beta.sd      = format(round(table.df$beta.sd, digits), nsmall = digits)
  table.df$pip          = format(round(table.df$pip, digits), nsmall = digits)
  table.df$pip[c(1, length(table.df$pip)-2, length(table.df$pip)-1, length(table.df$pip))] = "-"

if(sort){
table.df = table.df[order(table.df$pip, decreasing=T),]
}

colnames(table.df) = c("Variable", "Posterior Mean", "Posterior SD", "PIP")

return(
knitr::kable(table.df, type, booktabs=T, row.names = F, linesep = "", align="r")
)

}


#' @name topModels
#'
#' @title Extract Top Models from \code{ULLGM_BMA} Estimation Results
#'
#' @description \code{topModels} produces a table of the top n models from a \code{ULLGM_BMA} object, sorted by posterior model probabilities.
#'
#' @usage
#' topModels(x,
#'           variable_names = NULL,
#'           type           = "pandoc",
#'           digits         = 3,
#'           n              = 5)
#'
#' @param x The output object of \code{ULLGM_BMA}.
#' @param variable_names A character vector specifying the names of the columns of X.
#' @param digits Number of digits to round the table to. Defaults to 3.
#' @param n Number of top models to be returned. Defaults to 5.
#' @param type A character string indicating the format of the table. Options are \code{'pandoc'} (default), \code{'latex'}, or \code{'html'}.
#'
#' @return Returns a 'knitr::kable' object containing the table of top models.
#'
#' @author Gregor Zens
#'
#' @examples
#' # Load package
#' library(LatentBMA)
#'
#' # Example: Estimate a PLN model under a BRIC prior with m = p/2 using simulated data
#' # Note: Use more samples for actual analysis
#' # Note: nsave = 250 and nburn = 250 are for demonstration purposes
#' X <- matrix(rnorm(100*20), 100, 20)
#' z <- 2 + X %*% c(0.5, -0.5, rep(0, 18)) + rnorm(100, 0, sqrt(0.25))
#' y <- rpois(100, exp(z))
#' results_pln <- ULLGM_BMA(X = X, y = y, model = "PLN", nsave = 250, nburn = 250)
#' # Top 5 models
#' topModels(results_pln)
#'
#' @import knitr reshape2
#' @export
topModels = function(x,
                     variable_names = NULL,
                     type = "pandoc",
                     digits = 3,
                     n = 5){

  # specify variable names
  if(is.null(variable_names)){variable_names=colnames(x$X)}
  if(is.null(variable_names)){variable_names=paste0("x", 1:ncol(x$X))}

  # Step 1: Convert each row to a character string
  model_strings <- apply(x$gamma, 1, function(row) paste(row, collapse = ""))

  # Step 2: Count the occurrences of each unique model
  model_counts <- table(model_strings)

  # Step 3: Sort and find the top 5 models
  top_models <- sort(model_counts, decreasing = TRUE)[1:n]

  # Step 4: Calculate the proportion of appearances
  total_models <- length(model_strings)
  top_model_proportions <- top_models / total_models

  # Print the top 5 models and their proportions
  plot.df = data.frame(t(apply(do.call("cbind", lapply(1:length(top_models), function(i) unlist(strsplit(names(top_models[i]), "")))),2,as.numeric)))

  colnames(plot.df) = variable_names
  plot.df = plot.df[,colSums(plot.df)>0]
  plot.df$PMP = as.numeric(top_model_proportions)

  plot.df$model = paste0("Model #", 1:n)

  plot.df[plot.df == 1] = "x"
  plot.df[plot.df == 0] = ""
  plot.df$PMP = round(as.numeric(plot.df$PMP),digits)

  knitr::kable(t(plot.df[,c(ncol(plot.df), 1:(ncol(plot.df)-1))]), booktabs=T, type, linesep="")



}

#' @name plotPIP
#'
#' @title Visualization of Posterior Inclusion Probabilities
#'
#' @description \code{plotPIP} produces a visualization of the posterior inclusion probabilities (PIPs) extracted from \code{ULLGM_BMA} results.
#'
#' @usage
#' plotPIP(x,
#'         variable_names = NULL,
#'         sort           = TRUE)
#'
#' @param x The output object of \code{ULLGM_BMA}.
#' @param variable_names A character vector specifying the names of the columns of X.
#' @param sort Logical, indicating whether the plot should be sorted by PIP. Defaults to \code{TRUE}.
#'
#' @return Returns a 'ggplot2::ggplot' object.
#'
#' @author Gregor Zens
#'
#' @examples
#' # Load package
#' library(LatentBMA)
#'
#' # Example: Estimate a PLN model under a BRIC prior with m = p/2 using simulated data
#' # Note: Use more samples for actual analysis
#' # Note: nsave = 250 and nburn = 250 are for demonstration purposes
#' X <- matrix(rnorm(100*20), 100, 20)
#' z <- 2 + X %*% c(0.5, -0.5, rep(0, 18)) + rnorm(100, 0, sqrt(0.25))
#' y <- rpois(100, exp(z))
#' results_pln <- ULLGM_BMA(X = X, y = y, model = "PLN", nsave = 250, nburn = 250)
#' plotPIP(results_pln)
#'
#' @import ggplot2 reshape2
#' @export
plotPIP = function(x,
                   variable_names = NULL,
                   sort = TRUE){


  # specify variable names
  if(is.null(variable_names)){variable_names=colnames(x$X)}
  if(is.null(variable_names)){variable_names=paste0("x", 1:ncol(x$X))}

  # plot PIPs
  plot.df = data.frame(names = variable_names, pips = colMeans(x$gamma))
  plot.df = reshape2::melt(plot.df, id.vars="names")
  if(sort){plot.df$names = factor(plot.df$names, levels=(plot.df$names)[order(plot.df$value)])}

  return(
  ggplot(plot.df, aes(names, value)) +
    geom_hline(yintercept=0.5, col="grey20", lty=3) +
    geom_point() +
    geom_linerange(aes(ymin=0, ymax=value)) +
    coord_flip() +
    theme_bw() +
    xlab("") +
    ylab("PIP Estimate")
  )

}


#' @name plotBeta
#'
#' @title Visualization of Posterior Means of Coefficients
#'
#' @description \code{plotBeta} produces a visualization of the estimated posterior means of the coefficients, extracted from \code{ULLGM_BMA} results.
#'
#' @usage
#' plotBeta(x,
#'          variable_names = NULL,
#'          sort           = TRUE)
#'
#' @param x The output object of \code{ULLGM_BMA}.
#' @param variable_names A character vector specifying the names of the columns of X.
#' @param sort Logical, indicating whether the plot should be sorted by posterior mean. Defaults to \code{TRUE}.
#'
#' @return Returns a 'ggplot2::ggplot' object.
#'
#' @author Gregor Zens
#'
#' @examples
#' # Load package
#' library(LatentBMA)
#'
#' # Example: Estimate a PLN model under a BRIC prior with m = p/2 using simulated data
#' # Note: Use more samples for actual analysis
#' # Note: nsave = 250 and nburn = 250 are for demonstration purposes
#' X <- matrix(rnorm(100*20), 100, 20)
#' z <- 2 + X %*% c(0.5, -0.5, rep(0, 18)) + rnorm(100, 0, sqrt(0.25))
#' y <- rpois(100, exp(z))
#' results_pln <- ULLGM_BMA(X = X, y = y, model = "PLN", nsave = 250, nburn = 250)
#' plotBeta(results_pln)
#'
#' @import ggplot2 reshape2
#' @export
plotBeta = function(x,
                    variable_names = NULL,
                    sort = TRUE){


  # specify variable names
  if(is.null(variable_names)){variable_names=colnames(x$X)}
  if(is.null(variable_names)){variable_names=paste0("x", 1:ncol(x$X))}

  # plot PIPs
  plot.df = data.frame(names = variable_names, pips = colMeans(x$beta[,-1]))
  plot.df = reshape2::melt(plot.df, id.vars="names")
  if(sort){plot.df$names = factor(plot.df$names, levels=(plot.df$names)[order(plot.df$value)])}

  return(
    ggplot(plot.df, aes(names, value)) +
      geom_hline(yintercept=0, col="grey20", lty=3) +
      geom_point() +
      geom_linerange(aes(ymin=0, ymax=value)) +
      coord_flip() +
      theme_bw() +
      xlab("") +
      ylab("Beta Posterior Mean")
  )

}

#' @name plotModelSize
#'
#' @title Visualization of Model Size Posterior Distribution
#'
#' @description \code{plotModelSize} produces a visualization of the posterior distribution of model size, extracted from \code{ULLGM_BMA} results.
#'
#' @usage
#' plotModelSize(x)
#'
#' @param x The output object of \code{ULLGM_BMA}.
#'
#' @return Returns a 'ggplot2::ggplot' object visualizing the posterior distribution of model size.
#'
#' @author Gregor Zens
#'
#' @examples
#' # Load package
#' library(LatentBMA)
#'
#' # Example: Estimate a PLN model under a BRIC prior with m = p/2 using simulated data
#' # Note: Use more samples for actual analysis
#' # Note: nsave = 250 and nburn = 250 are for demonstration purposes
#' X <- matrix(rnorm(100*20), 100, 20)
#' z <- 2 + X %*% c(0.5, -0.5, rep(0, 18)) + rnorm(100, 0, sqrt(0.25))
#' y <- rpois(100, exp(z))
#' results_pln <- ULLGM_BMA(X = X, y = y, model = "PLN", nsave = 250, nburn = 250)
#' plotModelSize(results_pln)
#'
#' @import ggplot2 reshape2
#' @export
plotModelSize = function(x){

  plot.df = data.frame(S       = 1:ncol(x$X),
                       ms      = tabulate((rowSums(x$gamma)), ncol(x$X) ) / nrow((x$gamma)))

  plot.df = reshape2::melt(plot.df, id.vars="S")

  return(
  ggplot(plot.df, aes(x = S, y = value)) +
    geom_linerange(aes(ymin=0,  ymax=value)) +
    geom_point() +
    theme_bw() +
    xlab("Model Size") +
    ylab("Posterior Distribution")
)


}


#' @name tracePlot
#'
#' @title Traceplots for Selected Parameters
#'
#' @description \code{tracePlot} produces traceplots for selected parameters, extracted from \code{ULLGM_BMA} results.
#'
#' @usage
#' tracePlot(x, parameter = "beta", index = 1)
#'
#' @param x The output object of \code{ULLGM_BMA}.
#' @param parameter Specifies which parameter should be considered for the traceplot. Options are \code{"beta"} for coefficients, \code{"alpha"} for the intercept (the default), \code{"modelsize"} for model size, and \code{"sigma2"} for the error variance.
#' @param index If \code{parameter = "beta"}, specifies which coefficient should be shown. Defaults to 1, corresponding to the covariate in the first column of X.
#'
#' @return Returns a 'ggplot2::ggplot' object.
#'
#' @author Gregor Zens
#'
#' @export
tracePlot = function(x, parameter = "beta", index = 1){

  if(parameter == "alpha"){

    plot.df = data.frame(draws = 1:nrow(x$beta), trace = x$beta[,1])
    name = expression(alpha)

  }

  if(parameter == "sigma2"){

    plot.df = data.frame(draws = 1:nrow(x$beta), trace = x$sig2)
    name = expression(sigma^2)

  }

  if(parameter == "beta"){

    plot.df = data.frame(draws = 1:nrow(x$beta), trace = x$beta[,index+1])
    name = expression(beta[index])

  }

  if(parameter == "modelsize"){

    plot.df = data.frame(draws = 1:nrow(x$beta), trace = rowSums(x$gamma))
    name = "Model Size"

  }



  return(
    ggplot(plot.df, aes(x = draws, y = trace)) +
      geom_line() +
      theme_bw() +
      xlab("Iteration") +
      ylab(name)
  )


}



# to avoid CRAN note
globalVariables(c("S", "value", "draws"))
