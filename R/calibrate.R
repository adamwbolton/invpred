calibrate <- function(object, ...) {
  UseMethod("calibrate")
}



calibrate.lm <- function(object, y0, interval = c("inversion", "Wald", "none"),
                         level = 0.95, mean.response = FALSE,
                         adjust = c("none", "Bonferroni", "Scheffe"), k, ...) {

  # Check model formula for correctness
  xname <- all.vars(stats::formula(object)[[3L]])
  yname <- all.vars(stats::formula(object)[[2L]])
  if (length(xname) != 1L) {
    stop("Only one independent variable allowed.")
  }
  if (length(yname) != 1L) {
    stop("Only one dependent variable allowed.")
  }

  # Check for intercept using terms object from model fit. Alternatively, this
  # can also be checked by testing if the first column name in model.matrix is
  # equal to "(Intercept)".
  if (!attr(stats::terms(object), "intercept")) {
    stop(paste(deparse(substitute(object)), "must contain an intercept."))
  }

  # Extract x values and y values from model frame
  mf <- stats::model.frame(object)
  if (ncol(mf) != 2) {
    stop("calibrate only works for the simple linear regression model.")
  }
  x <- stats::model.matrix(object)[, 2]
  y <- stats::model.response(mf)

  # Eta - mean response or mean of observed respone values
  eta <- mean(y0)  # mean of new observations
  m <- length(y0)  # number of new observations
  if (mean.response && m > 1) stop("Only one mean response value allowed.")

  # Fit a simple linear regression model and compute necessary components
  b <- unname(object$coefficients)
  n <- length(r <- object$residuals)  # sample size and residuals
  DF <- (DF1 <- n - 2) + (DF2 <- m - 1)  # degrees of freedom
  var1 <- sum(r ^ 2) / object$df.residual  # stage 1 variance estimate
  var2 <- if (m == 1) 0 else stats::var(y0)  # stage 2 variance estimate
  var.pooled <- (DF1 * var1 + DF2 * var2) / DF  # pooled estimate of variance
  sigma.pooled <- sqrt(var.pooled)  # sqrt of pooled variance estimate
  ssx <- sum((x - mean(x))^2)  # sum-of-squares for x, Sxx
  x0.mle <- (eta - b[1L])/b[2L]  # MLE of x0

  # Return point estimate only
  interval <- match.arg(interval)
  if (interval == "none") return(x0.mle)

  # Adjustment for simultaneous intervals
  adjust <- match.arg(adjust)  # FIXME: Does simultaneous work for m > 1?
  crit <- if (m != 1 || adjust == "none") stats::qt((1 + level)/2, n+m-3) else {
    switch(adjust,
           "Bonferroni" = stats::qt((level + 2*k - 1) / (2*k), n+m-3),
           "Scheffe"    = sqrt(k * stats::qf(level, k, n+m-3)))
  }

  # Inversion interval --------------------------------------------------------
  if (interval == "inversion") {

    c1 <- b[2L]^2 - (sigma.pooled^2 * crit^2)/ssx
    c2 <- if (mean.response) {
      c1/n + (eta - mean(y))^2/ssx
    } else {
      c1*(1/m + 1/n) + (eta - mean(y))^2/ssx
    }
    c3 <- b[2L] * (eta - mean(y))
    c4 <- crit * sigma.pooled

    # FIXME: catch errors and throw an appropriate warning
    if (c1 < 0 && c2 <= 0) {

      warning("The calibration line is not well determined.", call. = FALSE)
      lwr <- -Inf
      upr <- Inf

    } else {

      lwr <- mean(x) + (c3 - c4*sqrt(c2))/c1
      upr <- mean(x) + (c3 + c4*sqrt(c2))/c1
      if (c1 < 0 && c2 > 0) {
        stop(paste("The calibration line is not well determined. The resulting \nconfidence region is the union of two semi-infinite intervals:\n(", -Inf, ",",
                   round(upr, 4), ") U (", round(lwr, 4), ",", Inf, ")"),
             call. = FALSE)
      }

    }
    res <- list("estimate" = x0.mle,
                "lower"    = lwr,
                "upper"    = upr,
                "interval" = interval)

  }

  # Wald interval
  if (interval == "Wald") {

    # Compute standard error for Wald interval
    se <- if (mean.response) {
      abs((sigma.pooled/b[2]))*sqrt((1/n + (x0.mle - mean(x))^2/ssx))
    } else {
      abs((sigma.pooled/b[2]))*sqrt((1/m + 1/n + (x0.mle - mean(x))^2/ssx))
    }

    # Store results in a list
    res <- list("estimate" = x0.mle,
                "lower"    = x0.mle - crit * se,
                "upper"    = x0.mle + crit * se,
                "se"       = se,
                "interval" = interval)

  }

  # Assign class label and return results
  class(res) <- "invest"
  res

}

