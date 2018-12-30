#' @export
print.itest <- function(x, digits = getOption("digits"), ...)
{

  cat("\n")
  cat(strwrap(x$name), sep = "\n")
  cat("\n")

  cat(sprintf("\tH0: %s\n", x$null))
  cat(sprintf("\tHA: %s\n\n", x$alternative))

  if (!is.null(x$statistic_value))
  {
    value <- format(signif(x$statistic_value, max(1L, digits - 2L)))
    cat(sprintf("\tTest statistic: %s = %s\n", x$statistic_name, value))
  }

  if (!is.null(x$pvalue))
  {
    pval <- format.pval(x$pvalue, digits = max(1L, digits - 3L))
    cat(sprintf("\tP-value: %s\n", pval))
  }

  if (!is.null(x$pointest))
  {
    point <- format(signif(x$pointest, max(1L, digits - 2L)))
    cat(sprintf("\n\tParameter: %s\n", x$parameter_name))
    cat(sprintf("\tPoint estimate: %s\n", point))
  }

  if (!is.null(x$se))
  {
    sef <- format(signif(x$se, max(1L, digits - 2L)))
    cat(sprintf("\tStandard Error: %s\n", sef))
  }

  if (!is.null(x$cint))
  {
    ci <- format(signif(x$cint, max(1L, digits - 2L)))
    cat(sprintf("\tConfidence interval: [%s, %s]\n", ci[1], ci[2]))
  }

  cat("\n")
}

#' Produce a contingency table from a data and formula object.
#'
#' @param formula   an object of class "formula"
#' @param data      a data frame object
#'
#' @return A table object with named rows and columns
#'
#' @export
tmod_contingency <- function(formula, data)
{
  # grab data and reparse for logistic regression
  mf <- stats::model.frame(formula, data = data)
  y <- factor(stats::model.response(mf))
  response_d <- attr(attr(mf, "terms"), "response")
  g <- factor(mf[[-response_d]])

  table(g, y, dnn = c("Predictor", "Response"))
}


#' Compute group means.
#'
#' @param formula   an object of class "formula"
#' @param data      a data frame object
#'
#' @return A named vector of means.
#'
#' @export
tmod_mean_by_group <- function(formula, data)
{
  # grab data and reparse for logistic regression
  mf <- stats::model.frame(formula, data = data)
  y <- as.numeric(stats::model.response(mf))
  response_d <- attr(attr(mf, "terms"), "response")
  g <- factor(mf[[-response_d]])

  tab <- tapply(y, g, mean)
  names(tab) <- sprintf("Mean of %s", names(tab))
  tab
}

#' Produce sample proportions.
#'
#' @param formula   an object of class "formula"
#' @param data      a data frame object
#'
#' @return A matrix of proportions. Column totals should (approximately) sum to one.
#'
#' @export
tmod_prop_by_group <- function(formula, data)
{
  # grab data and reparse for logistic regression
  mf <- stats::model.frame(formula, data = data)
  y <- factor(stats::model.response(mf))
  response_d <- attr(attr(mf, "terms"), "response")
  g <- factor(mf[[-response_d]])

  mat <- matrix(NA_real_, nrow = nlevels(y), ncol = nlevels(g))
  for (i in seq_len(nlevels(y)))
  {
    mat[i, ] <- tapply(y == levels(y)[i], g, mean)
  }

  colnames(mat) <- levels(g)
  rownames(mat) <- levels(y)
  mat
}



