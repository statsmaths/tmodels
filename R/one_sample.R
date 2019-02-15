#' Apply One Sample T-Test
#'
#' @param formula   an object of class "formula"
#' @param data      a data frame object
#' @param h0        mean to test as the null-hypothesis
#'
#' @return An itest object.
#'
#' @export
tmod_t_test_one <- function(formula, data, h0 = 0)
{
  # grab data and reparse for logistic regression
  mf <- stats::model.frame(formula, data = data)
  y <- stats::model.response(mf)

  # check numeric response
  if (!is.numeric(y))
  {
    stop("Non-numeric response variable provided.")
  }

  # run the t-test
  mod <- stats::t.test(y, mu = h0)

  # produce output
  digits <- getOption("digits")
  df <- format(signif(as.numeric(mod$parameter['df']), max(1L, digits - 2L)))
  out <- list(name = "One Sample t-test",
              statistic_name = sprintf("t(%s)", df),
              statistic_value = as.numeric(mod$statistic),
              null = sprintf("Mean is equal to %f", h0),
              alternative = sprintf("Mean is not equal to %f", h0),
              pvalue = mod$p.value,
              parameter_name = "Average value",
              pointest = as.numeric(mod$estimate),
              cint = c(mod$conf.int[1], mod$conf.int[2])
              )

  return(structure(out, class = "itest"))
}


#' Apply the one-sample Wilcoxon Test
#'
#' @param formula   an object of class "formula"
#' @param data      a data frame object
#' @param h0        mean to test as the null-hypothesis
#'
#' @return An itest object.
#'
#' @export
tmod_wilcoxon_test <- function(formula, data, h0 = 0)
{
  # grab data and reparse for logistic regression
  mf <- stats::model.frame(formula, data = data)
  y <- stats::model.response(mf)

  # check numeric response
  if (!is.numeric(y))
  {
    stop("Non-numeric response variable provided.")
  }

  # run the t-test
  mod <- stats::wilcox.test(y, mu = h0)

  # produce output
  out <- list(name = "Wilcoxon One-Sample Test",
              statistic_name = "V",
              statistic_value = as.numeric(mod$statistic),
              null = sprintf("Location shift is equal to %f", h0),
              alternative = sprintf("Location shift is not equal to %f", h0),
              pvalue = mod$p.value
              )

  return(structure(out, class = "itest"))
}

