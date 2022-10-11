#' @keywords internal
computeInversionInterval <- function(object, ...) {
  UseMethod("computeInversionInterval")
}


#' @keywords internal
computeInversionInterval.lm <- function(object, multi, x0.name, var1, var2,var.pooled, m,
                                        eta, crit, x0.est, mean.response, interval,
                                        newdata, lower, upper, extendInt, tol,
                                        maxiter) {

  # Pivotal quantity for linear model (pg. 95)
  rootfun <- function(x) {
    nd <- if (multi) {
      cbind(newdata, makeData(x, x0.name))  # append newdata
    } else {
      makeData(x, x0.name)
    }
    pred <- stats::predict(object, newdata = nd, se.fit = TRUE)
    denom <- if (mean.response) {
      pred$se.fit ^ 2
    } else if (interval %in% c("welch", "bf")) {
      var2/m + pred$se.fit^2
    } else if (interval == "t"){
      rat <- var.pooled/var1
      var.pooled/m + rat*pred$se.fit^2
    }
    (eta - pred$fit) ^ 2 / denom - crit ^ 2
  }

  # Compute lower and upper confidence limits
  lwr <- try(stats::uniroot(rootfun, interval = c(lower, x0.est),
                            extendInt = extendInt, tol = tol,
                            maxiter = maxiter)$root, silent = TRUE)

  upr <- try(stats::uniroot(rootfun, interval = c(x0.est, upper),
                            extendInt = extendInt, tol = tol,
                            maxiter = maxiter)$root, silent = TRUE)

  # Provide (informative) error message if confidence limits not found
  if (inherits(lwr, "try-error")) {
    stop(paste("Lower confidence limit not found in the search interval (",
               lower, ", ", upper,
               "). ", "Try tweaking the values of lower and upper. ",
               "Use plotFit for guidance.", sep = ""),
         call. = FALSE)
  }
  if (inherits(upr, "try-error")) {
    stop(paste("Upper confidence limit not found in the search interval (",
               lower, ", ", upper,
               "). ", "Try tweaking the values of lower and upper. ",
               "Use plotFit for guidance.", sep = ""),
         call. = FALSE)
  }

  # Return list of results
  list("estimate" = x0.est, "lower" = lwr, "upper" = upr,
       "interval" = "inversion")

}
