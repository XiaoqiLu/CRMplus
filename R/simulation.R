library(dfcrm)

set.seed(1024)

calibration <- read.csv("scenario/calibration.csv")
n.sim <- 100

Skeleton <- function(n, target, calibration)
{
  delta <- calibration$Delta[(calibration$N == n) & (calibration$Target == target)]
  return(getprior(delta, target, 3, 5))
}

SimWrapper <- function(n.sim, n.trial, Generator, skeleton, target, crm)
{
  mtd.sim <- NULL
  first.tox.sim <- NULL
  for (i.sim in 1 : n.sim)
  {
    sim <- Simulation(n.dose = 5, n.tox = length(target), n.trial = n.trial, Generator = Generator, skeleton = skeleton, model = "power", method = "mle", target = target, crm = crm)
    mtd.sim <- c(mtd.sim, sim$mtd)
    first.tox.sim <- rbind(first.tox.sim, sim$first.tox)
  }
  cat("\n(%)", format(1 : 5, width = 8), "\n")
  cat("MTD:", format(table(factor(mtd.sim, levels = 1 : 5)) / n.sim * 100, nsmall = 2, width = 8), "\n")
  cat("\n(%)", format(c(1 : n.trial, Inf), width = 3), "\n")
  for (tox in 0 : length(target))
  {
    cat(paste0("N", tox, ":"), format(round(table(factor(first.tox.sim[, tox + 1], levels = c(1 : n.trial, Inf))) / n.sim * 100), width = 3), "\n")
  }
  return(list(mtd = mtd.sim, first.tox = first.tox.sim))
}

for (i.scn in 1 : 7)
{
  scn <- read.csv(paste0("scenario/scenario_0", i.scn, ".csv"))
  cat("\n=================== scenario", i.scn, "==========================\n")

  cat("Clinician:           ", format(scn$Clinician, nsmall = 2, width = 6), "\n")
  cat("Patient:             ", format(scn$Patient, nsmall = 2, width = 6), "\n")
  cat("Clinician OR Patient:", format(scn$Clinician.Patient, nsmall = 2, width = 6), "\n")

  GenSimple <- function(dose) return(ToxGenerator(dose, "simple", scn))
  GenMarginal <- function(dose) return(ToxGenerator(dose, "marginal", scn))
  GenJoint <- function(dose) return(ToxGenerator(dose, "joint", scn))

  for (n.trial in c(21, 40, 18))
  {
    cat("\n----------- n.trial", n.trial, "-----------\n")
    skeleton25 <- Skeleton(n.trial, 0.25, calibration)
    skeleton35 <- Skeleton(n.trial, 0.35, calibration)
    skeleton50 <- Skeleton(n.trial, 0.50, calibration)

    # regular: simple encoding
    cat("\n[regular] regular - simple encoding\n")
    sim <- SimWrapper(n.sim = n.sim, n.trial = n.trial, Generator = GenSimple, skeleton = skeleton25, target = 0.25, crm = "regular")

    # parallel: marginal encoding
    cat("\n[marginal] parallel - marginal encoding\n")
    sim <- SimWrapper(n.sim = n.sim, n.trial = n.trial, Generator = GenMarginal, skeleton = cbind(skeleton35, skeleton25), target = c(0.35, 0.25), crm = "parallel")

    # regular: joint encoding
    cat("\n[joint] regular - joint encoding\n")
    sim <- SimWrapper(n.sim = n.sim, n.trial = n.trial, Generator = GenJoint, skeleton = skeleton25, target = c(0.50, 0.25), crm = "regular")

    # parallel: joint encoding
    cat("\n[naive] parallel - joint encoding\n")
    sim <- SimWrapper(n.sim = n.sim, n.trial = n.trial, Generator = GenJoint, skeleton = cbind(skeleton50, skeleton25), target = c(0.50, 0.25), crm = "parallel")

  }
}
