#' Apply a linear regression model
#'
#' @param formula   an object of class "formula"
#' @param data      a data frame object
#'
#' @return An itest object.
#'
#' @export
tmod_linear_regression <- function(formula, data)
{

  # construct objects to produce regression model
  mf <- stats::model.frame(formula, data = data)
  X <- stats::model.matrix(mf, data=data)
  y <- stats::model.response(mf)
  id <- attr(X, "assign")

  # build regression and anova models
  model1 <- .tmod_lm_fit(X, y, mf)
  model2 <- .tmod_lm_fit(X[,id != 1,drop=FALSE], y, mf)
  amod <- stats::anova(model2, model1)

  # details of the independent variable
  iv_type <- as.character(attr(model1$terms, "dataClasses")[2])
  nlevels <- sum(id == 1)
  iv_var <- paste(attr(model1$terms, "term.labels")[1], collapse="; ")
  nu_vars <- paste(attr(model1$terms, "term.labels")[-1], collapse="; ")
  v_var <- names(attr(model1$terms, "dataClasses")[1])

  # produce outputs
  if (iv_type == "numeric")
  {
    out <- as.numeric(stats::coef(summary(model1))[2,])
    digits <- getOption("digits")
    st <- "change in %s for unit change in %s\n\t            controlling for -- %s"
    out <- list(name = "Linear regression; T-Test",
                statistic_name = "t",
                statistic_value = out[3],
                null = "Change in conditional mean is zero",
                alternative = "Change in conditional mean is non-zero",
                pvalue = out[4],
                parameter_name = sprintf(st, v_var, iv_var, nu_vars),
                pointest = out[1],
                cint = as.numeric(stats::confint(model1)[2,])
                )
  } else if (iv_type %in% c("character", "factor") & nlevels == 1) {
    out <- as.numeric(stats::coef(summary(model1))[2,])
    levels <- model1$xlevels[[iv_var]]
    digits <- getOption("digits")
    st <- "mean(%s|%s) - mean(%s|%s)\n\t           after controlling for -- %s"
    out <- list(name = "Linear regression; T-Test",
                statistic_name = "t",
                statistic_value = out[3],
                null = "Difference in conditional mean is zero",
                alternative = "Difference in conditional mean is non-zero",
                pvalue = out[4],
                parameter_name = sprintf(st, v_var, levels[2], v_var, levels[1], nu_vars),
                pointest = out[1],
                cint = as.numeric(stats::confint(model1)[2,])
                )
  } else {
    digits <- getOption("digits")
    st <- "mean(%s|%s) - mean(%s|%s)\n\t           after controlling for -- %s"
    out <- list(name = "Linear regression; F-Test",
                statistic_name = sprintf("F(1, %d)", abs(amod$Df[2])),
                statistic_value = as.numeric(amod$F[2]),
                null = "Difference in conditional mean is zero",
                alternative = "Difference in true means is non-zero",
                pvalue = amod[['Pr(>F)']][2]
                )
  }

  return(structure(out, class = "itest"))

}

.tmod_lm_fit <- function(X, y, mf)
{
  mt <- attr(mf, "terms")
  z <- stats::lm.fit(X, y)
  class(z) <- "lm"
  z$na.action <- attr(mf, "na.action")
  z$offset <- NULL
  z$contrasts <- attr(X, "contrasts")
  z$xlevels <- stats::.getXlevels(mt, mf)
  z$call <- NULL
  z$terms <- mt
  return(z)
}

#' Apply a logistic regression model
#'
#' @param formula   an object of class "formula"
#' @param data      a data frame object
#' @param effect    optional name of the 'effect' category in y
#'
#' @return An itest object.
#'
#' @export
tmod_logistic_regression <- function(formula, data, effect=NULL)
{

  # construct objects to produce regression model
  mf <- stats::model.frame(formula, data = data)
  X <- stats::model.matrix(mf, data = data)
  y <- stats::model.response(mf)
  id <- attr(X, "assign")

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

  # turn y into a categorical variables and make sure just two levels
  yo <- factor(y)
  ylevels <- levels(yo)
  y <- (as.integer(yo) - 1L)
  if (nlevels(yo) != 2) stop("The response variable has ", nlevels(yo),
                            " unique values. This test requires exactly 2.")

  # build regression and anova models
  model1 <- .tmod_glm_fit(X, y, mf)
  model2 <- .tmod_glm_fit(X[,id != 1,drop=FALSE], y, mf)
  amod <- stats::anova(model2, model1, test = "Rao")

  # details of the independent variable
  iv_type <- as.character(attr(model1$terms, "dataClasses")[2])
  nlevels <- sum(id == 1)
  iv_var <- paste(attr(model1$terms, "term.labels")[1], collapse="; ")
  nu_vars <- paste(attr(model1$terms, "term.labels")[-1], collapse="; ")
  v_var <- names(attr(model1$terms, "dataClasses")[1])

  if (iv_type == "numeric")
  {
    out <- as.numeric(stats::coef(summary(model1))[2,])
    digits <- getOption("digits")
    st <- "change in OR(%s=%s) for unit change in %s\n\t           controlling for -- %s"
    out <- list(name = "Logistic regression; Z-Test",
                statistic_name = "Z",
                statistic_value = out[3],
                null = "Change in conditional log odds is zero",
                alternative = "Change in conditional log odds is non-zero",
                pvalue = out[4],
                parameter_name = sprintf(st, v_var, ylevels[1], iv_var, nu_vars),
                pointest = out[1],
                cint = c(out[1] - 1.96 *  out[2],  out[1] + 1.96 *  out[2])
                )
  } else if (iv_type %in% c("character", "factor") & nlevels == 1) {
    out <- as.numeric(stats::coef(summary(model1))[2,])
    levels <- model1$xlevels[[iv_var]]
    digits <- getOption("digits")
    st <- "OR(%s=%s|%s) - OR(%s=%s|%s)\n\t           after controlling for -- %s"
    out <- list(name = "Logistic regression; Z-Test",
                statistic_name = "Z",
                statistic_value = out[3],
                null = "Difference in conditional log odds is zero",
                alternative = "Difference in conditional log odds is non-zero",
                pvalue = out[4],
                parameter_name = sprintf(st, v_var, ylevels[1], levels[2],
                                         v_var, ylevels[1], levels[1], nu_vars),
                pointest = out[1],
                cint = c(out[1] - 1.96 *  out[2],  out[1] + 1.96 *  out[2])
                )
  } else {
    digits <- getOption("digits")
    st <- "OR(%s=%s|%s) - OR(%s=%s|%s)\n\t           after controlling for -- %s"
    out <- list(name = "Logistic regression; F-Test",
                statistic_name = "Rao",
                statistic_value = as.numeric(amod$Rao[2]),
                null = "Difference in conditional log odds is zero",
                alternative = "Difference in true log odds is non-zero",
                pvalue = amod[['Pr(>Chi)']][2]
                )
  }



  return(structure(out, class = "itest"))

}

.tmod_glm_fit <- function(X, y, mf)
{
  mt <- attr(mf, "terms")
  z <- stats::lm.fit(X, y)
  class(z) <- "lm"
  z$na.action <- attr(mf, "na.action")
  z$offset <- NULL
  z$contrasts <- attr(X, "contrasts")
  z$xlevels <- stats::.getXlevels(mt, mf)
  z$call <- NULL
  z$terms <- mt
  return(z)
}

#' Get a linear regression model table
#'
#' @param formula   an object of class "formula"
#' @param data      a data frame object
#'
#' @return An itest object.
#'
#' @export
tmod_lin_reg_table <- function(formula, data)
{
  object <- stats::lm(formula, data)
  ci_vals <- stats::confint(object, level = 0.95)
  co <- stats::coef(summary(object))

  co <- cbind(co[,1], ci_vals, co[,4])
  colnames(co) <- c("Estimate", "CI Low", "CI High", "Pr(>|t|)")
  signif(co, 3)
}

