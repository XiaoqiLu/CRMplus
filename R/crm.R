#' Compute the Log Probability of Dose-Toxicity Data, Given Toxicity Rates
#'
#' @param dose.tox counts of dose-toxicity, an integer matrix where nrow = # of dose levels, and ncol = # of toxicity levels
#' @param tox.rate toxicity rates computed from model parameter, a matrix of the same dimension as dose.tox, rowSums(tox.rate) should be ones
#'
#' @return
#' @export
#'
#' @examples
LogLik <- function(dose.tox, tox.rate)
{
  loglik <- sum(dose.tox * log(tox.rate))
  return(loglik)
}

#' Compute Toxicity Rates for Power Working Model, Given Parameters and Skeleton
#'
#' @param b parameters of power working model, a real vector of length = (# of toxicity levels - 1)
#' @param skeleton mapped skeleton for dose levels, a vector of length = # of dose levels, should be between 0 and 1
#'
#' @return
#' @export
#'
#' @examples
PowerWorkingModel <- function(b, skeleton)
{
  cum.tox.rate <- cbind(1, outer(skeleton, cumsum(exp(b)), "^"), 0)
  m <- ncol(cum.tox.rate)
  tox.rate <- cum.tox.rate[, 1 : (m - 1)] - cum.tox.rate[, 2 : m]
  return(tox.rate)
}

#' Compute Toxicity Rates for Working Model
#'
#' @param param parameters of working model
#' @param skeleton mapped skeleton for dose levels
#' @param model working model, character, default = "power"
#' @param ... other (fixed) parameters of working model
#'
#' @return
#' @export
#'
#' @examples
WorkingModel <- function(param, skeleton, model = "power", ...)
{
  switch(model,
         power = PowerWorkingModel(param, skeleton),
         paste(model, "model currently NOT supported"))
}

# test --------------------------------------------------------------------

dose.tox <- matrix(c(1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1), 4, 3, byrow = TRUE)
tox.rate <- WorkingModel(c(-1, 1), c(0.1, 0.5, 0.6, 0.9), model = "power")
print(tox.rate)

print(LogLik(dose.tox, tox.rate))
