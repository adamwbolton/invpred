#' @keywords internal
computeInverseEstimate <- function(object, ...) {
  UseMethod("computeInverseEstimate")
}


#' @keywords internal
computeInverseEstimate.lm <- function(object, multi, x0.name, newdata, eta,
                                      lower, upper, extendInt, tol, maxiter) {

  # The function f(x0; betas) - y0
  rootfun <- function(x, object, multi, newdata, x0.name, eta) {
    nd <- if (multi) {
      cbind(newdata, makeData(x, x0.name))  # append newdata
    } else {
      makeData(x, x0.name)
    }
    stats::predict(object, newdata = nd, type = "response") - eta
  }

  # Try to solve the equation f(x0; betas) - y0 = 0 for x0
  res <- try(
    stats::uniroot(rootfun, interval = c(lower, upper), object = object,
                   multi = multi, newdata = newdata, x0.name = x0.name,
                   eta = eta, extendInt = extendInt, tol = tol,
                   maxiter = maxiter)$root,
    silent = FALSE
  )

  # Provide (informative) error message if point estimate is not found
  if (inherits(res, "try-error")) {
    stop(paste("Point estimate not found in the search interval (", lower,
               ", ", upper, "). ",
               "Try tweaking the values of lower and upper.", sep = ""),
         call. = FALSE)
  } else {
    res
  }

}
