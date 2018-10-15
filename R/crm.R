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
  cum.tox.rate <- outer(skeleton, cumsum(exp(b)), "^")
  tox.rate <- cbind(1, cum.tox.rate) - cbind(cum.tox.rate, 0)
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
tox.rate <- WorkingModel(c(-1, 1), c(0.1, 0.5, 0.6, 0.9), model = "logistic")
print(tox.rate)

print(LogLik(dose.tox, tox.rate))
