#' Inverse prediction for linear models
#'
#' Provides point estimate and prediction intervals for unknown predictor value(s)
#' that corresponds to its observed proxy measurement or specified value of the mean response.
#'
#' @param object An object that inherits from class \code{\link[stats]{lm}}
#'
#' @param y0 The value of the observed response(s).  \code{y0} should be on
#' the scale of the response variable (e.g., a number between 0 and 1 for
#' binomial families).
#'
#' @param interval The type of interval required.
#'
#' @param level A numeric scalar between 0 and 1 giving the confidence level
#'   for the interval to be calculated.
#'
#' @param mean.response Logical indicating whether confidence intervals should
#' correspond to an individual response (\code{FALSE}) or a mean response
#' (\code{TRUE}). For \code{\link[stats]{glm}} objects, this is always
#' \code{TRUE}.
#'
#' @param x0.name For multiple linear regression, a character string giving the
#' name of the predictor variable of interest.
#'
#' @param newdata For multiple linear regression, a \code{data.frame} giving the
#' values of interest for all other predictor variables (i.e., those other than
#' \code{x0.name}).
#'
#' @param data An optional data frame. This is required if \code{object$data} is
#' \code{NULL}.
#'
#' @param nsim Positive integer specifying the number of bootstrap simulations;
#' the bootstrap B (or R).
#'
#' @param boot.type Character string specifying the type of bootstrap to use
#' when \code{interval = "percentile"}. Options are \code{"parametric"} and
#' \code{"nonparametric"}.
#'
#' @param seed Optional argument to \code{\link{set.seed}}.
#'
#' @param progress Logical indicating whether to display a text-based progress
#' bar during the bootstrap simulation.
#'
#' @param lower The lower endpoint of the interval to be searched.
#'
#' @param upper The upper endpoint of the interval to be searched.
#'
#' @param extendInt Character string specifying if the interval
#' \code{c(lower, upper)} should be extended or directly produce an error when
#' the inverse of the prediction function does not have differing signs at the
#' endpoints. The default, \code{"no"}, keeps the search interval and hence
#' produces an error. Can be abbreviated. See the documentation for the
#' \code{base} R function \code{uniroot} for details.
#'
#' @param q1 Optional lower cutoff to be used in forming confidence intervals.
#' Only used when \code{object} inherits from class \code{\link[nlme]{lme}}.
#' Defaults to \code{stats::qnorm((1+level)/2)}.
#'
#' @param q2 Optional upper cutoff to be used in forming confidence intervals.
#' Only used when \code{object} inherits from class \code{\link[nlme]{lme}}.
#' Defaults to \code{stats::qnorm((1-level)/2)}.
#'
#' @param tol The desired accuracy passed on to \code{\link[stats]{uniroot}}.
#' Recommend a minimum of \code{1e-10}.
#'
#' @param maxiter The maximum number of iterations passed on to \code{uniroot}.
#'
#' @param ... Additional optional arguments. At present, no optional arguments
#' are used.
#'
#' @return Returns an object of class \code{"invest"} or, if
#' \code{interval = "percentile"}, of class \code{c("invest", "bootCal")}. The
#' generic function \code{{plot}} can be used to plot the output
#' of the bootstrap simulation when \code{interval = "percentile"}.
#'
#'   An object of class \code{"invest"} containing the following components:
#'   \itemize{
#'     \item \code{estimate} The estimate of x0.
#'     \item \code{lwr} The lower confidence limit for x0.
#'     \item \code{upr} The upper confidence limit for x0.
#'     \item \code{se} An estimate of the standard error (Wald and percentile
#'                     intervals only).
#'     \item \code{bias} The bootstrap estimate of bias (percentile interval
#'                       only).
#'     \item \code{bootreps} Vector of bootstrap replicates (percentile
#'                           interval only).
#'     \item \code{nsim} The number of bootstrap replicates (percentile
#'                       interval only).
#'     \item \code{interval} The method used for calculating \code{lower} and
#'           \code{upper} (only used by \code{{print}} method).
#'   }
#' @rdname invpred
#'
#' @export
invpred <- function(object, y0, ...) {
  UseMethod("invpred")
}


#' @rdname invpred
#'
#' @export
invpred.lm <- function(object, y0,
                      interval = c("welch", "t", "bf", "wald"),
                      level = 0.95, mean.response = FALSE,
                      x0.name, newdata, data, lower, upper,
                      extendInt = "no", tol = .Machine$double.eps^0.25,
                      maxiter = 1000, k) {


  # Extract data, variable names, etc.
  .data  <- if (!missing(data)) {
    data
  } else {
    eval(object$call$data, envir = parent.frame())
  }
  yname <- all.vars(stats::formula(object)[[2]])

  # Predictor variable(s)
  multi <- FALSE
  xnames <- intersect(all.vars(stats::formula(object)[[3]]), colnames(.data))
  if (length(xnames) != 1) {
    multi <- TRUE
    if (missing(x0.name)) {
      stop("'x0.name' is missing, please select a valid predictor variable")
    }
    xnames <- xnames[xnames != x0.name]
    if (missing(newdata)) {
      # FIXME: Should the user be warned here about the default behavior?
      newdata <- as.data.frame(lapply(.data[, xnames], FUN = function(x) {
        if (is.numeric(x)) {
          stats::median(x, na.rm = TRUE)
        } else {
          names(which.max(table(x, useNA = "always")))
        }
      }))
    }
    if (!is.data.frame(newdata)) {
      stop("'newdata' must be a data frame")
    }
    if (nrow(newdata) != 1) {
      stop("'newdata' must have a single row")
    }
    if (ncol(newdata) != length(xnames) || !all(xnames %in% names(newdata))) {
      stop(paste0("'newdata' must contain a column for each predictor variable",
                  " used by ", deparse(substitute(object)),
                  " (except ", x0.name, ")"))
    }
  } else {
    x0.name <- xnames
  }

  # End-points for 'uniroot'
  if (missing(lower)) lower <- min(.data[, x0.name])  # lower limit default
  if (missing(upper)) upper <- max(.data[, x0.name])  # upper limit default

  # Define constants
  m <- length(y0)  # number of unknowns
  if(mean.response && m > 1) {
    stop("Only one mean response value allowed.")
  }
  eta <- mean(y0)  # mean unknown
  n <- length(stats::resid(object))  # in case of missing values
  p <- length(stats::coef(object))   # number of regression coefficients
  df1 <- n - p  # stage 1 degrees of freedom
  df2 <- m - 1  # stage 2 degrees of freedom
  var1 <- stats::sigma(object)^2  # stage 1 variance estimate
  var2 <- if (m == 1) {
    0
  } else {
    stats::var(y0)  # stage 2 variance estimate
  }


  sd1 <- sqrt(var1)
  sd2 <- sqrt(var2)

  var.pooled <- (df1 * var1 + df2 * var2) / (df1 + df2)  # pooled estimate

  # Calculate point estimate by inverting fitted model
  x0.est <- computeInverseEstimate(object, multi = multi, x0.name = x0.name,
                                   newdata = newdata, eta = eta, lower = lower,
                                   upper = upper, extendInt = extendInt,
                                   tol = tol, maxiter = maxiter)

  # Match arguments
  interval <- match.arg(interval)

  # Return point estimate only
  if (interval == "none") {
    return(stats::setNames(x0.est, x0.name))
  }


  crit <- if (interval == "t") {
    stats::qt((level + 1) / 2, df1 + df2)
  } else if (m > 1) {
    if (interval == "welch") {
      df <- ws.dof(df1, df2, sd1, sd2)
      stats::qt((level + 1) / 2, df)
    } else if (interval == "bf") {
      asht::qbf(p = (level+1)/2, s1 = sd1, s2 =sd2, n1 = df1 + 1, n2 = df2 + 1)
    }
  } else {
    stop("Must have m > 1 to estimate sd2.")
  }

  # Calculate confidence limits
  res <-
    # Inversion confidence/prediction interval
    computeInversionInterval(object, multi = multi, x0.name = x0.name,
                             var1 = var1, var2 = var2, var.pooled=var.pooled, interval,
                             m = m, eta = eta, crit = crit, x0.est = x0.est,
                             newdata = newdata, mean.response = mean.response,
                             lower = lower, upper = upper,
                             extendInt = extendInt, tol = tol,
                             maxiter = maxiter)


  # Assign class label(s) and return result
  class(res) <- "invest"
  res

}
