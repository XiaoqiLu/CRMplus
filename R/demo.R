source("R/crm.R")
source("R/pro.R")

seed <- 1

n.sim <- 10000

target.c <- 0.25
target.p <- 0.35
target.cop <- 0.50

calibration.table <- read.csv("scenario/calibration.csv", header = TRUE)

Skeleton <- function(n, target)
{
  require(dfcrm)
  hw <- calibration.table$Delta[(calibration.table$N == n) & (calibration.table$Target == target)]
  return(getprior(hw, target, 3, 5))
}



for (n in c(18, 40))
{
  cat("\n================= n", n, "==================\n")

  for (i.scenario in 1 : 7)
  {
    cat("\n-------------------- scenario", i.scenario, "--------------------\n")

    scenario.table <- as.matrix(read.csv(paste0("scenario/scenario_0", i.scenario, ".csv"), header = TRUE))

    cat("Clinician:           ", format(scenario.table[, 1], nsmall = 2, width = 6), "\n")
    cat("Patient:             ", format(scenario.table[, 2], nsmall = 2, width = 6), "\n")
    cat("Clinician OR Patient:", format(scenario.table[, 3], nsmall = 2, width = 6), "\n")

    # non-pro
    method <- "nopro"
    target <- target.c
    skeleton <- Skeleton(n, target)
    file.name <- paste0("result/n_", n, "_scenario_", i.scenario, "_method_", method, "RData")
    if (!file.exists(file.name))
    {
      sim <- SimCRM(n.sim, n, scenario.table[, 1, drop = FALSE], skeleton, target)
      save(sim, file = file.name)
    } else
    {
      load(file.name)
    }
    cat("No - PRO\n")
    print(table(factor(sim$mtd, 1 : 5)) / n.sim * 100)

    # marginal
    method <- "mar"
    target <- c(target.p, target.c)
    skeleton <- cbind(Skeleton(n, target[1]), Skeleton(n, target[2]))
    file.name <- paste0("result/n_", n, "_scenario_", i.scenario, "_method_", method, "RData")
    if (!file.exists(file.name))
    {
      sim <- SimPROCRM(n.sim, n, scenario.table, skeleton, target, method)
      save(sim, file = file.name)
    } else
    {
      load(file.name)
    }
    cat("Marginal\n")
    print(table(factor(sim$mtd, 1 : 5)) / n.sim * 100)

    # joint - joint
    method <- "jntjnt"
    target <- c(target.cop, target.c)
    skeleton <- Skeleton(n, target[2])
    file.name <- paste0("result/n_", n, "_scenario_", i.scenario, "_method_", method, "RData")
    if (!file.exists(file.name))
    {
      sim <- SimPROCRM(n.sim, n, scenario.table, skeleton, target, method)
      save(sim, file = file.name)
    } else
    {
      load(file.name)
    }
    cat("Joint\n")
    print(table(factor(sim$mtd, 1 : 5)) / n.sim * 100)

    # joint - marginal
    method <- "jntmar"
    target <- c(target.cop, target.c)
    skeleton <- cbind(Skeleton(n, target[1]), Skeleton(n, target[2]))
    file.name <- paste0("result/n_", n, "_scenario_", i.scenario, "_method_", method, "RData")
    if (!file.exists(file.name))
    {
      sim <- SimPROCRM(n.sim, n, scenario.table, skeleton, target, method)
      save(sim, file = file.name)
    } else
    {
      load(file.name)
    }
    cat("Joint/Marginal\n")
    print(table(factor(sim$mtd, 1 : 5)) / n.sim * 100)
  }
}

par(mfrow = c(1, 3))
n <- 18
i.scenario <- 3
seed <- 38

# i.scenario <- 5
# seed <- 3
# seed <- 5

cat("\n================= n", n, "==================\n")

cat("\n-------------------- scenario", i.scenario, "--------------------\n")

scenario.table <- as.matrix(read.csv(paste0("scenario/scenario_0", i.scenario, ".csv"), header = TRUE))

cat("Clinician:           ", format(scenario.table[, 1], nsmall = 2, width = 6), "\n")
cat("Patient:             ", format(scenario.table[, 2], nsmall = 2, width = 6), "\n")
cat("Clinician OR Patient:", format(scenario.table[, 3], nsmall = 2, width = 6), "\n")

# marginal
method <- "mar"
target <- c(target.p, target.c)
skeleton <- cbind(Skeleton(n, target[1]), Skeleton(n, target[2]))
sim <- SimOnePROCRM(n, scenario.table, skeleton, target, method, seed)
dose <- sim$dose
yc <- sim$tox.c
yp <- sim$tox.p
pt <- rep(1, n)
pt[(yc == 1) & (yp == 1)] <- 8
pt[(yc == 0) & (yp == 1)] <- 3
pt[(yc == 1) & (yp == 0)] <- 4
plot(1 : n, dose, type = "p", pch = pt, main = "Marginal", cex = 2, ylim = c(1, 5), xlab = "Patient", ylab = "Dose Level")
legend(x = "topright", legend = c("None", "Patient Only", "Clinician Only", "Both"), pch = c(1, 3, 4, 8))
cat(method, "\n")
print(round(t(sim$tox.rate[, , 1]) * 100, 1))
print(round(t(sim$tox.rate[, , 2]) * 100, 1))
# print(round(t(sim$param), 4))

# joint - joint
method <- "jntjnt"
target <- c(target.cop, target.c)
skeleton <- Skeleton(n, target[2])
sim <- SimOnePROCRM(n, scenario.table, skeleton, target, method, seed)
dose <- sim$dose
yc <- sim$tox.c
yp <- sim$tox.p
pt <- rep(1, n)
pt[(yc == 1) & (yp == 1)] <- 8
pt[(yc == 0) & (yp == 1)] <- 3
pt[(yc == 1) & (yp == 0)] <- 4
plot(1 : n, dose, type = "p", pch = pt, main = "Joint", cex = 2, ylim = c(1, 5), xlab = "Patient", ylab = "Dose Level")
legend(x = "topright", legend = c("None", "Patient Only", "Clinician Only", "Both"), pch = c(1, 3, 4, 8))
cat(method, "\n")
print(round(t(sim$tox.rate[, , 1]) * 100, 1))
print(round(t(sim$tox.rate[, , 2]) * 100, 1))
# print(round(t(sim$param), 4))

# joint - marginal
method <- "jntmar"
target <- c(target.cop, target.c)
skeleton <- cbind(Skeleton(n, target[1]), Skeleton(n, target[2]))
sim <- SimOnePROCRM(n, scenario.table, skeleton, target, method, seed)
dose <- sim$dose
yc <- sim$tox.c
yp <- sim$tox.p
pt <- rep(1, n)
pt[(yc == 1) & (yp == 1)] <- 8
pt[(yc == 0) & (yp == 1)] <- 3
pt[(yc == 1) & (yp == 0)] <- 4
plot(1 : n, dose, type = "p", pch = pt, main = "Joint/Marginal", cex = 2, ylim = c(1, 5), xlab = "Patient", ylab = "Dose Level")
legend(x = "topright", legend = c("None", "Patient Only", "Clinician Only", "Both"), pch = c(1, 3, 4, 8))
cat(method, "\n")
print(round(t(sim$tox.rate[, , 1]) * 100, 1))
print(round(t(sim$tox.rate[, , 2]) * 100, 1))
# print(round(t(sim$param), 4))
