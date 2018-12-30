#' Apply the Pearson correlation test
#'
#' @param formula   an object of class "formula"
#' @param data      a data frame object
#'
#' @return An itest object.
#'
#' @export
tmod_pearson_correlation_test <- function(formula, data)
{
  # grab data and reparse for logistic regression
  mf <- stats::model.frame(formula, data = data)
  y <- stats::model.response(mf)
  response_d <- attr(attr(mf, "terms"), "response")
  x <- mf[[-response_d]]

  # check numeric variables
  if (!is.numeric(x))
  {
    stop("Non-numeric predictor variable provided.")
  }
  if (!is.numeric(y))
  {
    stop("Non-numeric response variable provided.")
  }

  # run the correlation text
  newdf <- data.frame(y = y, x = x, stringsAsFactors=FALSE)
  mod <- stats::cor.test(~ x + y, data = newdf)

  # produce output
  df <- as.numeric(mod$parameter['df'])
  out <- list(name = "Pearson's product-moment correlation test",
              statistic_name = sprintf("t(%d)", df),
              statistic_value = as.numeric(mod$statistic),
              null = "True correlation is zero",
              alternative = "True correlation is non-zero",
              pvalue = mod$p.value,
              parameter_name = "(Pearson) correlation coefficient",
              pointest = as.numeric(mod$estimate),
              cint = c(mod$conf.int[1], mod$conf.int[2])
              )

  return(structure(out, class = "itest"))
}

#' Apply the Spearman rho correlation test
#'
#' @param formula   an object of class "formula"
#' @param data      a data frame object
#'
#' @return An itest object.
#'
#' @export
tmod_spearman_correlation_test <- function(formula, data)
{
  # grab data and reparse for logistic regression
  mf <- stats::model.frame(formula, data = data)
  y <- stats::model.response(mf)
  response_d <- attr(attr(mf, "terms"), "response")
  x <- mf[[-response_d]]

  # check numeric variables
  if (!is.numeric(x))
  {
    stop("Non-numeric predictor variable provided.")
  }
  if (!is.numeric(y))
  {
    stop("Non-numeric response variable provided.")
  }

  # run the correlation text
  newdf <- data.frame(y = y, x = x, stringsAsFactors=FALSE)
  mod <- stats::cor.test(~ x + y, data = newdf, method = "spearman")

  # produce output
  df <- as.numeric(mod$parameter['df'])
  out <- list(name = "Spearman's rank correlation rho test",
              statistic_name = "S",
              statistic_value = as.numeric(mod$statistic),
              null = "True rho is zero",
              alternative = "True rho is non-zero",
              pvalue = mod$p.value,
              parameter_name = "Spearman's correlation (rho)",
              pointest = as.numeric(mod$estimate)
              )

  return(structure(out, class = "itest"))
}

#' Apply the Kendall tau correlation test
#'
#' @param formula   an object of class "formula"
#' @param data      a data frame object
#'
#' @return An itest object.
#'
#' @export
tmod_kendall_correlation_test <- function(formula, data)
{
  # grab data and reparse for logistic regression
  mf <- stats::model.frame(formula, data = data)
  y <- stats::model.response(mf)
  response_d <- attr(attr(mf, "terms"), "response")
  x <- mf[[-response_d]]

  # check numeric variables
  if (!is.numeric(x))
  {
    stop("Non-numeric predictor variable provided.")
  }
  if (!is.numeric(y))
  {
    stop("Non-numeric response variable provided.")
  }

  # run the correlation text
  newdf <- data.frame(y = y, x = x, stringsAsFactors=FALSE)
  mod <- stats::cor.test(~ x + y, data = newdf, method = "kendall")

  # produce output
  df <- as.numeric(mod$parameter['df'])
  out <- list(name = "Kendall's rank correlation tau test",
              statistic_name = "z",
              statistic_value = as.numeric(mod$statistic),
              null = "True tau is zero",
              alternative = "True tau is non-zero",
              pvalue = mod$p.value,
              parameter_name = "Kendall's rank correlation (tau)",
              pointest = as.numeric(mod$estimate)
              )

  return(structure(out, class = "itest"))
}
