#' Compute (Cumulative) Toxicity Rates for Working Model
#'
#' @param param parameters of working model, a real vector of length = # of tocicity levels (L)
#' @param skeleton mapped skeleton for dose levels, a vector of length = # of dose levels (K), should be between 0 and 1
#' @param working.model working model, character, default = "power" (currently only supports "power")
#'
#' @return
#' @export
#'
#' @examples
ToxRate <- function(param, skeleton, working.model = "power")
{
  switch(working.model,
         power = ToxRatePower(param, skeleton),
         paste(working.model, "model currently NOT supported"))
}

#' Compute Toxicity Rates for Power Working Model, Given Parameters and Skeleton
#'
#' @param b parameters of power working model, a real vector of length = # of toxicity levels (L)
#' @inheritParams ToxRate
#'
#' @return
#' @export
#'
#' @examples
ToxRatePower <- function(b, skeleton)
{
  tox.rate <- outer(skeleton, cumsum(exp(b)), "^")
  return(tox.rate)
}

#' Compute the Log Probability of Dose-Toxicity Data, Given Toxicity Rates
#'
#' @param dose.tox counts of dose-toxicity, an integer matrix where nrow = # of dose levels (K), and ncol = # of toxicity levels (L) + 1.
#' @param tox.rate (cumulative) toxicity rates computed from model parameter, a matrix where nrow = # of dose levels (K), and ncol = # of toxicity levels (L).
#'
#' @return
#' @export
#'
#' @examples
LogProb <- function(dose.tox, tox.rate)
{

  log.prob <- sum((dose.tox + sqrt(.Machine$double.eps)) * log(cbind(1, tox.rate) - cbind(tox.rate, 0)))
  return(log.prob)
}

#' Compute the Log Likelihood of Parameter, Given Dose-Toxicity Data and Skeleton, a wrapper for LogProb()
#'
#' @inheritParams LogProb
#' @inheritParams ToxRate
#'
#' @return
#' @export
#'
#' @examples
LogLik <- function(param, dose.tox, skeleton, working.model = "power")
{
  tox.rate <- ToxRate(param, skeleton, working.model)
  log.lik <- LogProb(dose.tox, tox.rate)

  return(log.lik)
}

#' Compute Maximum Tolerated Dose (MTD) from Toxicity Rates and Target
#'
#' @param target targets for (cumulative) toxicity rates, a vector of length = # of toxicity levels (L)
#' @inheritParams LogProb
#'
#' @return
#' @export
#'
#' @examples
MTD <- function(tox.rate, target)
{
  mtd <- min(apply(abs(t(tox.rate) - target), 1, which.min))
  return(mtd)
}

#' Fit Working Model by Maximum Likelihood Estimate (MLE)
#'
#' @param param.ini initial value of paramameter, default = rep(0, L)
#' @inheritParams LogLik
#'
#' @return
#' @export
#'
#' @examples
FitModel <- function(param.ini = NULL, dose.tox, skeleton, working.model = "power")
{
  if (is.null(param.ini)) param.ini <- rep(0, ncol(dose.tox) - 1)

  if (length(param.ini) == 1)
  {
    param <- optim(param.ini, LogLik, dose.tox = dose.tox, skeleton = skeleton, working.model = working.model,
                   method = "Brent", lower = -5, upper = 5, control = list(fnscale = -1))$par
  } else
  {
    param <- optim(param.ini, LogLik, dose.tox = dose.tox, skeleton = skeleton, working.model = working.model,
                   method = "Nelder-Mead", control = list(fnscale = -1))$par
  }

  return(param)
}

#' Convert Trial Sequence to Dose-Toxicity Data (Counts)
#'
#' @param dose dose level sequence
#' @param tox toxicity sequence
#' @param K # of possible dose levels
#' @param L # of possible toxicity levels (0 not included)
#'
#' @return
#' @export
#'
#' @examples
Seq2DoseTox <- function(dose, tox, K, L)
{
  dose.tox <- table(factor(dose, 1 : K), factor(tox, 0 : L))
  return(dose.tox)
}

#'  Dose Assignment/Recommendation Given MTD Estiamte and Trial Sequence, with Restrictions
#'
#' @param mtd MTD
#' @param cohort.size default = 1
#' @param allow.dose.skipping default = FALSE
#' @param allow.escalation.on.toxicity default = FALSE
#' @inheritParams Seq2DoseTox
#'
#' @return
#' @export
#'
#' @examples
Assign <- function(mtd, dose, tox, cohort.size = 1, allow.dose.skipping = FALSE, allow.escalation.on.toxicity = FALSE)
{
  assign <- mtd

  if (!allow.dose.skipping)
  {
    max.dose <- max(dose)
    if (assign > max.dose + 1) assign <- max.dose + 1
  }

  if (!allow.escalation.on.toxicity)
  {
    n <- length(dose)
    if (any(tox[(n - cohort.size + 1) : n] > 0)) assign <- min(assign, dose[n])
  }

  return(assign)
}

#' Continual Reassessment Method (CRM) with Multiple Toxicity Constraints
#'
#' @inheritParams Assign
#' @inheritParams MTD
#' @inheritParams FitModel
#'
#' @return
#' @export
#'
#' @examples
CRM <- function(param.ini = NULL, dose, tox, skeleton, target, working.model = "power",
                    cohort.size = 1, allow.dose.skipping = FALSE, allow.escalation.on.toxicity = FALSE)
{
  K <- length(skeleton)
  L <- length(target)

  dose.tox <- Seq2DoseTox(dose, tox, K, L)
  param <- FitModel(param.ini, dose.tox, skeleton, working.model)
  tox.rate <- ToxRate(param, skeleton, working.model)
  mtd <- MTD(tox.rate, target)
  assign <- Assign(mtd, dose, tox, cohort.size, allow.dose.skipping, allow.escalation.on.toxicity)

  return(list(param = param, tox.rate = tox.rate, mtd = mtd, assign = assign))
}

#' Single-Trial Simulation for CRM()
#'
#' @param n # of subjects (sample size) in each trial, should be a multiple of cohort.size
#' @param true.tox.rate true (cumulative) toxicity rates, a real matrix where nrow = # of dose levels (K), and ncol = # of toxicity levels (L).
#' @param seed seed for RNG
#' @inheritParams CRM
#'
#' @return
#' @export
#'
#' @examples
SimOneCRM <- function(n, true.tox.rate, skeleton, target, working.model = "power",
                      cohort.size = 1, allow.dose.skipping = FALSE, allow.escalation.on.toxicity = FALSE, seed = NULL)
{
  if (!is.null(seed)) set.seed(seed)

  K <- length(skeleton)
  L <- length(target)

  param <- matrix(NA, n, L)
  tox.rate <- array(NA, c(n, K, L))
  mtd <- rep(NA, n)

  dose <- NULL
  tox <- NULL

  param.prev <- NULL
  dose.next <- 1
  for (i.cohort in 1 : (n / cohort.size))
  {
    # generate data
    dose <- c(dose, rep(dose.next, cohort.size))
    for (j in 1 : cohort.size)
    {
      tox.next <- sum(runif(1) < true.tox.rate[dose.next, ])
      tox <- c(tox, tox.next)
    }

    # crm
    crm <- CRM(param.prev, dose, tox, skeleton, target, working.model, cohort.size, allow.dose.skipping, allow.escalation.on.toxicity)
    param[i.cohort * cohort.size, ] <- param.prev <- crm$param
    tox.rate[i.cohort * cohort.size, , ] <- crm$tox.rate
    mtd[i.cohort * cohort.size] <- crm$mtd
    dose.next <- crm$assign
  }

  return(list(dose = dose, tox = tox, param = param, tox.rate = tox.rate, mtd = mtd))
}

#' Multi-Trial Wrapper for SimOneCRM()
#'
#' @param n.sim # of simulations/trials
#' @inheritParams SimOneCRM
#'
#' @return
#' @export
#'
#' @examples
SimCRM <- function(n.sim, n, true.tox.rate, skeleton, target, working.model = "power",
                   cohort.size = 1, allow.dose.skipping = FALSE, allow.escalation.on.toxicity = FALSE, seed = NULL)
{
  if (!is.null(seed)) set.seed(seed)

  dose <- matrix(NA, n.sim, n)
  tox <- matrix(NA, n.sim, n)
  mtd <- rep(NA, n.sim)

  for (i.sim in 1 : n.sim)
  {
    sim <- SimOneCRM(n, true.tox.rate, skeleton, target, working.model, cohort.size, allow.dose.skipping, allow.escalation.on.toxicity, seed)
    dose[i.sim, ] <- sim$dose
    tox[i.sim, ] <- sim$tox
    mtd[i.sim] <- sim$mtd[n]
  }

  true.mtd <- MTD(true.tox.rate, target)
  pca <- colMeans(dose == true.mtd)
  pcs <- mean(mtd == true.mtd)

  return(list(dose = dose, tox = tox, mtd = mtd, true.mtd = true.mtd, pca = pca, pcs = pcs))
}
