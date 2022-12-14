#' @keywords internal
#' @export
makeData <- function(x, label) {
  stats::setNames(data.frame(x), label)
}


#' @keywords internal
makeX <- function(object, newdata) {
  # Fixed effects model matrix
  stats::model.matrix(eval(object$call$fixed)[-2], data = newdata)
}


#' @keywords internal
makeZ <- function(object, newdata) {
  # Random effects model matrix
  Q <- object$dims$Q  # number of grouping levels
  mCall <- object$call  # list containing image of the nlme call
  fixed <- eval(eval(mCall$fixed)[-2])  # fixed effects formula
  reSt <- object$modelStruct$reStruct  # random effects structure
  mfArgs <- list(formula = nlme::asOneFormula(stats::formula(reSt), fixed),
                 data = newdata, na.action = stats::na.fail,
                 drop.unused.levels = TRUE)
  dataMix <- do.call("model.frame", mfArgs)
  stats::model.matrix(reSt, dataMix)
}


#' @keywords internal
sigma.lme <- function(object, ...) {
  object$sigma  # estimated standard deviation of the within-group error
}


#' @keywords internal
varY <- function(object, newdata) {

  # FIXME: What if object$call$correlation is not NULL?

  # Unconditional response variance: Var[Y] = Var[X*beta + Z*alpha + error]
  Zmat <- makeZ(object, newdata)  # random effects design matrix
  Gmat <- nlme::getVarCov(object)  # random effects variance-covariance matrix
  var.y <- Zmat %*% Gmat %*% t(Zmat) + stats::sigma(object) ^ 2  # ZGZ' + (sigma^2)I
  if (is.matrix(var.y)) {
    unname(diag(var.y))
  } else {
    var.y
  }
}



ws.dof <- function(df1, df2, sd1,sd2) {
  k1 <- df1 + 1
  k2 <- df2 + 1
  ((sd1^2/k1 + sd2^2/k2)^2/((sd1^2/k1)^2/(df1) + (sd2^2/k2)^2/(df2)))
}
