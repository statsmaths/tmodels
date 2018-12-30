#' Apply Odds Ratio Test for Two or More Groups
#'
#' @param formula   an object of class "formula"
#' @param data      a data frame object
#' @param effect    optional name of the 'effect' category in y
#' @param treatment optional name of the 'treatment' category in x
#'
#' @return An itest object.
#'
#' @export
tmod_odds_ratio <- function(formula, data, effect = NULL, treatment = NULL)
{
  # grab data and reparse for logistic regression
  mf <- stats::model.frame(formula, data = data)
  y <- factor(stats::model.response(mf))
  response_d <- attr(attr(mf, "terms"), "response")
  g <- factor(mf[[-response_d]])

  # fix levels (if needed) and store the results
  if (!is.null(effect))
  {
    if (!(effect %in% levels(y)))
    {
      stop("An effect variable was selected that does not exist in the response variable.")
    }
    y <- stats::relevel(y, effect)
    y <- factor(y, levels = rev(levels(y)))
  }
  if (!is.null(treatment))
  {
    if (!(treatment %in% levels(g)))
    {
      stop("A treatment variable was selected that does not exist in the response variable.")
    }
    y <- stats::relevel(y, treatment)
  }
  glevels <- levels(g)
  ylevels <- levels(y)

  # create new data frame
  yi <- (as.integer(y) - 1L)
  newdf <- data.frame(yi = yi, g = g, stringsAsFactors=FALSE)

  # check for 2xk contingency table
  if (nlevels(y) != 2) stop("The response variable has ", nlevels(y),
                            " unique values. This test requires exactly 2.")

  k <- nlevels(g)
  if (k <= 1)
  {
    stop("Odds ratio test requires 2 or more groups.")
  }

  # select correct test
  if (k == 2)
  {
    out <- .odds_ratio_2x2(newdf, glevels, ylevels)
  } else {
    out <- .odds_ratio_2xk(newdf, glevels, ylevels)
  }

  return(out)
}

.odds_ratio_2x2 <- function(newdf, glevels, ylevels)
{
  # run logistic regression
  model <- stats::glm(yi ~ g, data = newdf, family = stats::binomial())

  # extract output components
  pval <- summary(model)$coefficients[-1,4]
  or <- as.numeric(exp(stats::coef(model))[-1])
  se <- summary(model)$coefficients[-1,2]
  Z <- summary(model)$coefficients[-1,3]
  ci <- exp(stats::confint(model, trace=FALSE))[-1, ]
  attr(ci, "conf.level") <- 0.95

  out <- list(name = "Odds Ratio Test (2 groups)",
              statistic_name = "Z",
              statistic_value = Z,
              null = "odds ratio is equal to 1",
              alternative = "odds ratio is not equal to 1",
              pvalue = pval,
              parameter_name = sprintf("odds ratio (%s / %s) of %s",
                                       glevels[1], glevels[2], ylevels[2]),
              pointest = or,
              se = se,
              cint = ci
              )

  return(structure(out, class = "itest"))
}

.odds_ratio_2xk <- function(newdf, glevels, ylevels)
{
  # run logistic regression
  model <- stats::glm(yi ~ g, data = newdf, family = stats::binomial())

  # extract output components
  amod <- stats::anova(model, test = "Rao")
  pval <- rev(amod[["Pr(>Chi)"]])[1]
  rao <- rev(amod[["Rao"]])[1]
  df <- rev(amod[["Df"]])[1]

  out <- list(name = "Odds Ratio Test (3 or more groups)",
              statistic_name = sprintf("chi-squared(%d)", df),
              statistic_value = rao,
              null = "all pairs of odds ratios are equal to 1",
              alternative = "one or more odds ratios are not equal to 1",
              pvalue = pval
              )

  return(structure(out, class = "itest"))
}

#' Apply Fisher's Exact Test for Count Data
#'
#' @param formula   an object of class "formula"
#' @param data      a data frame object
#'
#' @return An itest object.
#'
#' @export
tmod_fisher_test <- function(formula, data)
{
  # grab data and reparse for logistic regression
  mf <- stats::model.frame(formula, data = data)
  y <- factor(stats::model.response(mf))
  response_d <- attr(attr(mf, "terms"), "response")
  g <- factor(mf[[-response_d]])

  # create new data frame
  yi <- (as.integer(y) - 1L)
  newdf <- data.frame(yi = yi, g = g, stringsAsFactors=FALSE)

  # run fisher test
  mod <- stats::fisher.test(newdf$g, newdf$yi)

  # produce output
  out <- list(name = "Fisher's Exact Test for Count Data",
              null = "Group variable and response categories are independent",
              alternative = "Group variable and response categories are dependent",
              pvalue = mod$p.value
              )

  return(structure(out, class = "itest"))
}

#' Apply Pearson's Chi-squared test with Yates' continuity correction
#'
#' @param formula   an object of class "formula"
#' @param data      a data frame object
#'
#' @return An itest object.
#'
#' @export
tmod_chi_squared_test <- function(formula, data)
{
  # grab data and reparse for logistic regression
  mf <- stats::model.frame(formula, data = data)
  y <- factor(stats::model.response(mf))
  response_d <- attr(attr(mf, "terms"), "response")
  g <- factor(mf[[-response_d]])

  # create new data frame
  yi <- (as.integer(y) - 1L)
  newdf <- data.frame(yi = yi, g = g, stringsAsFactors=FALSE)

  # run fisher test
  mod <- stats::chisq.test(newdf$g, newdf$yi)

  # produce output
  df <- as.numeric(mod$parameter['df'])
  out <- list(name = "Pearson's Chi-squared test with Yates' continuity correction",
              null = "Group variable and response categories are independent",
              alternative = "Group variable and response categories are dependent",
              pvalue = mod$p.value,
              statistic_name = sprintf("chi-squared(%d)", df),
              statistic_value = mod$statistic
              )

  return(structure(out, class = "itest"))
}

