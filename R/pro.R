#' Generate Patient-Reported Outcome (PRO), Given Toxicity Rates from Clinician, Patient, and Clinician/Patient
#'
#' @param prob.c toxicity rate from clinician
#' @param prob.p toxicity rate from patient
#' @param prob.cop toxicity rate from clinician or patient
#'
#' @return
#' @export
#'
#' @examples
RTox <- function(prob.c, prob.p, prob.cop)
{
  p00 <- 1 - prob.cop # none: 1
  p01 <- prob.cop - prob.c # patient only: 2
  p10 <- prob.cop - prob.p # clinician only: 3
  p11 <- prob.c + prob.p - prob.cop # both: 4

  ind <- sample(1 : 4, 1, prob = c(p00, p01, p10, p11))
  tox.c <- as.numeric(ind > 2) # 3 and 4
  tox.p <- as.numeric(ind %% 2 == 0) # 2 and 4

  return(list(tox.c = tox.c, tox.p = tox.p))
}

#' PRO-CRM marginal decoding (and, of course, marginal modeling)
#'
#' @param tox.c toxicity outcome from clinitian
#' @param tox.p toxicity outcome from patient
#' @inheritParams CRM
#'
#' @return
#' @export
#'
#' @examples
PROCRMmar <- function(param.ini = NULL, dose, tox.c, tox.p, skeleton, target)
{
  crm1 <- CRM(param.ini[1], dose, tox.p, skeleton[, 1], target[1])
  crm2 <- CRM(param.ini[2], dose, tox.c, skeleton[, 2], target[2])

  param <- c(crm1$param, crm2$param)
  tox.rate <- cbind(crm1$tox.rate, crm2$tox.rate)
  mtd <- min(crm1$mtd, crm2$mtd)
  assign <- min(crm1$assign, crm2$assign)

  return(list(param = param, tox.rate = tox.rate, mtd = mtd, assign = assign))
}

#' PRO-CRM joint decoding marginal modeling
#'
#' @param tox.c toxicity outcome from clinitian
#' @param tox.p toxicity outcome from patient
#' @inheritParams CRM
#'
#' @return
#' @export
#'
#' @examples
PROCRMjntmar <- function(param.ini = NULL, dose, tox.c, tox.p, skeleton, target)
{
  crm1 <- CRM(param.ini[1], dose, pmax(tox.c, tox.p), skeleton[, 1], target[1])
  crm2 <- CRM(param.ini[2], dose, tox.c, skeleton[, 2], target[2])

  param <- c(crm1$param, crm2$param)
  tox.rate <- cbind(crm1$tox.rate, crm2$tox.rate)
  mtd <- min(crm1$mtd, crm2$mtd)
  assign <- min(crm1$assign, crm2$assign)

  return(list(param = param, tox.rate = tox.rate, mtd = mtd, assign = assign))
}

#' PRO-CRM joint decoding joint modeling
#'
#' @param tox.c toxicity outcome from clinitian
#' @param tox.p toxicity outcome from patient
#' @inheritParams CRM
#'
#' @return
#' @export
#'
#' @examples
PROCRMjntjnt <- function(param.ini = NULL, dose, tox.c, tox.p, skeleton, target)
{
  tox <- tox.c + pmax(tox.c, tox.p)
  crm <- CRM(param.ini, dose, tox, skeleton, target)

  return(crm)
}

#' PRO-CRM wrapper
#'
#' @param method one of "mar", "jntmar", "jntjnt"
#' @inheritParams PROCRMmar
#'
#' @return
#' @export
#'
#' @examples
PROCRM <- function(param.ini = NULL, dose, tox.c, tox.p, skeleton, target, method)
{
  switch(method,
         mar = PROCRMmar(param.ini, dose, tox.c, tox.p, skeleton, target),
         jntmar = PROCRMjntmar(param.ini, dose, tox.c, tox.p, skeleton, target),
         jntjnt = PROCRMjntjnt(param.ini, dose, tox.c, tox.p, skeleton, target))
}

#' Single-Trial Simulation for PROCRM()
#'
#' @param n # of subjects (sample size) in each trial, should be a multiple of cohort.size
#' @param scenario.table scenario table (toxicity rates) for C, P, and CoP
#' @param seed seed for RNG
#' @inheritParams PROCRM
#'
#' @return
#' @export
#'
#' @examples
SimOnePROCRM <- function(n, scenario.table, skeleton, target, method, seed = NULL)
{
  if (!is.null(seed)) set.seed(seed)

  K <- nrow(scenario.table)

  param <- matrix(NA, n, 2)
  tox.rate <- array(NA, c(n, K, 2))
  mtd <- rep(NA, n)

  dose <- NULL
  tox.c <- NULL
  tox.p <- NULL

  param.prev <- NULL
  dose.next <- 1
  for (i in 1 : n)
  {
    # generate data
    dose <- c(dose, dose.next)
    tox.next <- RTox(scenario.table[dose.next, 1], scenario.table[dose.next, 2], scenario.table[dose.next, 3])
    tox.c <- c(tox.c, tox.next$tox.c)
    tox.p <- c(tox.p, tox.next$tox.p)

    # pro-crm
    crm <- PROCRM(param.prev, dose, tox.c, tox.p, skeleton, target, method)
    param[i, ] <- param.prev <- crm$param
    tox.rate[i, , ] <- crm$tox.rate
    mtd[i] <- crm$mtd
    dose.next <- crm$assign
  }

  return(list(dose = dose, tox.c = tox.c, tox.p = tox.p, param = param, tox.rate = tox.rate, mtd = mtd))
}

#' Multi-Trial Wrapper for SimOnePROCRM()
#'
#' @param n.sim # of simulations/trials
#' @inheritParams SimOnePROCRM
#'
#' @return
#' @export
#'
#' @examples
SimPROCRM <- function(n.sim, n, scenario.table, skeleton, target, method, seed = NULL)
{
  if (!is.null(seed)) set.seed(seed)

  dose <- matrix(NA, n.sim, n)
  tox.c <- matrix(NA, n.sim, n)
  tox.p <- matrix(NA, n.sim, n)
  mtd <- rep(NA, n.sim)

  for (i.sim in 1 : n.sim)
  {
    sim <- SimOnePROCRM(n, scenario.table, skeleton, target, method, seed)
    dose[i.sim, ] <- sim$dose
    tox.c[i.sim, ] <- sim$tox.c
    tox.p[i.sim, ] <- sim$tox.p
    mtd[i.sim] <- sim$mtd[n]
  }

  return(list(dose = dose, tox.c = tox.c, tox.p = tox.p, mtd = mtd))
}
