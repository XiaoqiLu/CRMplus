#' Compute the Log Probability of Dose-Toxicity Data, Given Toxicity Rates
#'
#' @param dose.tox counts of dose-toxicity, an integer matrix where nrow = # of dose levels, and ncol = (# of toxicity levels + 1)
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
#' @param b parameters of power working model, a real vector of length = # of toxicity levels
#' @param skeleton mapped skeleton for dose levels, a vector of length = # of dose levels, should be between 0 and 1
#'
#' @return
#' @export
#'
#' @examples
WorkingModelPower <- function(b, skeleton)
{
  cum.tox.rate <- outer(skeleton, cumsum(exp(b)), "^")
  tox.rate <- cbind(1, cum.tox.rate) - cbind(cum.tox.rate, 0)
  return(list(tox.rate = tox.rate, cum.tox.rate = cum.tox.rate))
}

#' Compute Toxicity Rates for Working Model
#'
#' @param param parameters of working model, a real vector of length = # of toxicity levels
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
         power = WorkingModelPower(param, skeleton),
         paste(model, "model currently NOT supported"))
}

#' Compute Maimum Tolerated Dose (MTD) from Working Model
#'
#' @inheritParams WorkingModel
#' @param target target cumulative toxicity rates, a vector of length = # of toxicity levels
#'
#' @return
#' @export
#'
#' @examples
MTD <- function(param, skeleton, model = "power", target, ...)
{
  cum.tox.rate <- WorkingModel(param, skeleton, model, ...)$cum.tox.rate
  # mtd <- min(colSums(cum.tox.rate <= matrix(target, length(skeleton), length(target), byrow = TRUE)))
  mtd <- min(apply(abs(cum.tox.rate - matrix(target, length(skeleton), length(target), byrow = TRUE)), 2, which.min))
  return(mtd)
}

#' Fit Working Model by Maximum Likelihood Estimate (MLE)
#'
#' @inheritParams LogLik
#' @inheritParams WorkingModel
#'
#' @return
#' @export
#'
#' @examples
FitWorkingModelMLE <- function(dose.tox, skeleton, model = "power", ...)
{
  param <- rep(0, ncol(dose.tox) - 1)

  Fun <- function(x)
  {
    tox.rate <- WorkingModel(x, skeleton, model, ...)$tox.rate
    loglik <- LogLik(dose.tox, tox.rate)
    return(- loglik)
  }

  if (length(param) == 1)
  {
    param <- optim(param, Fun, method = "Brent", lower = -10, upper = 10)$par
  } else
  {
    param <- optim(param, Fun)$par
  }

  return(param)
}

#' Fit Working Model
#'
#' @inheritParams FitWorkingModelMLE
#' @param method fitting method, default = "mle"
#'
#' @return
#' @export
#'
#' @examples
FitWorkingModel <- function(dose.tox, skeleton, model = "power", method = "mle", ...)
{
  switch(method,
         mle = FitWorkingModelMLE(dose.tox, skeleton, model, ...),
         paste(method, "method currently NOT supported"))
}

#' Convert Trial Sequence to Dose-Toxicity Data (Counts)
#'
#' @param dose dose level sequence
#' @param tox toxicity sequence
#' @param n.dose number of possible dose levels
#' @param n.tox number of possible toxicity levels (0 not included)
#'
#' @return
#' @export
#'
#' @examples
Seq2DoseTox <- function(dose, tox, n.dose, n.tox)
{
  dose.tox <- matrix(0, n.dose, n.tox + 1)
  for (i in seq_along(dose))
  {
    dose.tox[dose[i], tox[i] + 1] <- dose.tox[dose[i], tox[i] + 1] + 1
  }
  return(dose.tox)
}

#' Dose Assignment/Recommendation Given MTD Estiamte and Highest Assigned Dose, Avoid Skipping Dose
#'
#' @param mtd MTD estimate
#' @param max.dose highest assigned dose so far
#'
#' @return
#' @export
#'
#' @examples
Assign <- function(mtd, max.dose)
{
  return(min(max(mtd, 1), max.dose + 1))
}

#' Continual Reassessment Method (CRM)
#'
#' @inheritParams FitWorkingModel
#' @inheritParams MTD
#' @inheritParams Assign
#'
#' @return
#' @export
#'
#' @examples
CRM <- function(dose.tox, skeleton, model = "power", method = "mle", target, max.dose = NULL, ...)
{
  if (is.null(max.dose))
  {
    dose.hist <- rowSums(dose.tox)
    if (sum(dose.hist) == 0)
    {
      max.dose <- 0
    } else
    {
      max.dose <- max(seq_along(dose.hist)[dose.hist > 0])
    }
  }

  param <- FitWorkingModel(dose.tox, skeleton, model, method, ...)
  mtd <- MTD(param, skeleton, model, target, ...)
  dose.next <- Assign(mtd, max.dose)

  return(list(param = param, mtd = mtd, dose.next = dose.next))
}


#' Parallel CRM
#'
#' @param skeleton mapped skeleton for dose levels, for parallel CRM it's a matrix of (# of dose levels)-by-(# of toxicity types)
#' @inheritParams CRM
#'
#' @return
#' @export
#'
#' @examples
ParallelCRM <- function(dose.tox, skeleton, model = "power", method = "mle", target, max.dose = NULL, ...)
{
  if (is.null(max.dose))
  {
    dose.hist <- rowSums(dose.tox)
    if (sum(dose.hist) == 0)
    {
      max.dose <- 0
    } else
    {
      max.dose <- max(seq_along(dose.hist)[dose.hist > 0])
    }
  }

  mtd <- nrow(dose.tox)
  param <- NULL
  for (k in 1 : length(target))
  {
    res <- CRM(dose.tox[, c(1, k + 1)], skeleton[, k], model, method, target[k], max.dose, ...)
    param <- cbind(param, res$param)
    mtd <- min(mtd, res$mtd)
  }
  dose.next <- Assign(mtd, max.dose)

  return(list(param = param, mtd = mtd, dose.next = dose.next))
}

#' CRM Simulation for One Subject
#'
#' @inheritParams Seq2DoseTox
#' @param n.trial # of trials in a sequence
#' @param Generator generator function to simulate trial data, takes dose as input and returns a random toxicity
#' @inheritParams CRM
#' @param crm CRM type, default = "regular" (use CRM function), alternative is "parallel" (use ParallelCRM function)
#'
#' @return
#' @export
#'
#' @examples
Simulation <- function(n.dose, n.tox, n.trial, Generator, skeleton, model = "power", method = "mle", target, crm = "regular", ...)
{
  dose.tox <- matrix(0, n.dose, n.tox + 1)
  max.dose <- 0
  dose.next <- 1

  for (i.trial in 1 : n.trial)
  {
    max.dose <- max(max.dose, dose.next)

    tox.next <- Generator(dose.next)
    dose.tox[dose.next, tox.next + 1] <- dose.tox[dose.next, tox.next + 1] + 1

    res <- switch(crm,
                  regular = CRM(dose.tox, skeleton, model, method, target, max.dose, ...),
                  parallel = ParallelCRM(dose.tox, skeleton, model, method, target, max.dose, ...),
                  paste(crm, "currently NOT supported"))

    mtd <- res$mtd
    dose.next <- res$dose.next
  }

  return(list(dose.tox = dose.tox, mtd = mtd))
}

# test --------------------------------------------------------------------

# param <- c(-1, 1)
# skeleton <- c(0.1, 0.5, 0.6, 0.9)
# dose <- c(1, 1, 2, 2)
# tox <- c(0, 0, 0, 0)
# dose.tox <- Seq2DoseTox(dose, tox, 4, 2)
# tox.rate <- WorkingModel(param, skeleton, model = "power")$tox.rate
# print(tox.rate)
#
# print(LogLik(dose.tox, tox.rate))
# mtd <- MTD(param, skeleton, target = c(0.5, 0.25))
# print(mtd)
#
# param.fit <- FitWorkingModelMLE(dose.tox, skeleton)
# print(param.fit)
#
# print(WorkingModel(param.fit, skeleton)$tox.rate)
#
# print(CRM(dose.tox[, c(1, 2)], skeleton, target = 0.5))
# print(ParallelCRM(dose.tox, cbind(skeleton, skeleton), target = c(0.5, 0.25)))
