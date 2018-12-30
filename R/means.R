#' Apply a two-sample t-test
#'
#' @param formula   an object of class "formula"
#' @param data      a data frame object
#'
#' @return An itest object.
#'
#' @export
tmod_t_test <- function(formula, data)
{
  # grab data and reparse for logistic regression
  mf <- stats::model.frame(formula, data = data)
  y <- stats::model.response(mf)
  response_d <- attr(attr(mf, "terms"), "response")
  x <- factor(mf[[-response_d]])
  xlevels <- levels(x)

  # check numeric response
  if (!is.numeric(y))
  {
    stop("Non-numeric response variable provided.")
  }

  # check number of groups in predictor variable
  if (nlevels(x) != 2) stop("The predictor variable has ", nlevels(x),
                            " unique values. This test requires exactly 2.")

  # run the t-test
  newdf <- data.frame(y = y, x = x, stringsAsFactors=FALSE)
  mod <- stats::t.test(y ~ x, data = newdf)

  # produce output
  digits <- getOption("digits")
  df <- format(signif(as.numeric(mod$parameter['df']), max(1L, digits - 2L)))
  out <- list(name = "Two Sample t-test",
              statistic_name = sprintf("t(%s)", df),
              statistic_value = as.numeric(mod$statistic),
              null = "Difference in true means is zero",
              alternative = "Difference in true means is non-zero",
              pvalue = mod$p.value,
              parameter_name = sprintf("Difference in means (%s - %s)",
                                       xlevels[1], xlevels[2]),
              pointest = as.numeric(diff(rev(mod$estimate))),
              cint = c(mod$conf.int[1], mod$conf.int[2])
              )

  return(structure(out, class = "itest"))
}

#' Apply the two-sample Mann-Whitney Test
#'
#' @param formula   an object of class "formula"
#' @param data      a data frame object
#'
#' @return An itest object.
#'
#' @export
tmod_mann_whitney_test <- function(formula, data)
{
  # grab data and reparse for logistic regression
  mf <- stats::model.frame(formula, data = data)
  y <- stats::model.response(mf)
  response_d <- attr(attr(mf, "terms"), "response")
  x <- factor(mf[[-response_d]])
  xlevels <- levels(x)

  # check numeric response
  if (!is.numeric(y))
  {
    stop("Non-numeric response variable provided.")
  }

  # check number of groups in predictor variable
  if (nlevels(x) != 2) stop("The predictor variable has ", nlevels(x),
                            " unique values. This test requires exactly 2.")

  # run the t-test
  newdf <- data.frame(y = y, x = x, stringsAsFactors=FALSE)
  mod <- stats::wilcox.test(y ~ x, data = newdf)

  # produce output
  out <- list(name = "Mann-Whitney rank sum test with continuity correction",
              statistic_name = "W",
              statistic_value = as.numeric(mod$statistic),
              null = "True location shift is equal to zero",
              alternative = "True location shift is not equal to zero",
              pvalue = mod$p.value
              )

  return(structure(out, class = "itest"))
}

#' Apply the one-way Analysis of Variance test
#'
#' @param formula   an object of class "formula"
#' @param data      a data frame object
#'
#' @return An itest object.
#'
#' @export
tmod_one_way_anova_test <- function(formula, data)
{
  # grab data and reparse for logistic regression
  mf <- stats::model.frame(formula, data = data)
  y <- stats::model.response(mf)
  response_d <- attr(attr(mf, "terms"), "response")
  x <- factor(mf[[-response_d]])
  xlevels <- levels(x)

  # check numeric response
  if (!is.numeric(y))
  {
    stop("Non-numeric response variable provided.")
  }

  # run the t-test
  newdf <- data.frame(y = y, x = x, stringsAsFactors=FALSE)
  mod <- stats::anova(stats::lm(y ~ x, data = newdf))

  # produce output
  out <- list(name = "One-way Analysis of Variance (ANOVA)",
              statistic_name = sprintf("F(%d, %d)", mod$Df[1], mod$Df[2]),
              statistic_value = as.numeric(mod[['F value']][1]),
              null = "True means are the same in each group.",
              alternative = "True means are the same in each group.",
              pvalue = mod[['Pr(>F)']][1]
              )

  return(structure(out, class = "itest"))
}


#' Apply the Kruskal-Wallis test
#'
#' @param formula   an object of class "formula"
#' @param data      a data frame object
#'
#' @return An itest object.
#'
#' @export
tmod_kruskal_wallis_test <- function(formula, data)
{
  # grab data and reparse for logistic regression
  mf <- stats::model.frame(formula, data = data)
  y <- stats::model.response(mf)
  response_d <- attr(attr(mf, "terms"), "response")
  x <- factor(mf[[-response_d]])
  xlevels <- levels(x)

  # check numeric response
  if (!is.numeric(y))
  {
    stop("Non-numeric response variable provided.")
  }

  # run the t-test
  newdf <- data.frame(y = y, x = x, stringsAsFactors=FALSE)
  mod <- stats::kruskal.test(y ~ x, data = newdf)

  # produce output
  digits <- getOption("digits")
  df <- format(signif(as.numeric(mod$parameter['df']), max(1L, digits - 2L)))
  out <- list(name = "Kruskal-Wallis rank sum test",
              statistic_name = sprintf("chi-squared(%s)", df),
              statistic_value = as.numeric(mod$statistic),
              null = "Location parameters are the same in each group.",
              alternative = "Location parameters are not the same in each group.",
              pvalue = mod$p.value
              )

  return(structure(out, class = "itest"))
}



