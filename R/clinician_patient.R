#' Generate a code Representing Toxicity Reported by Clinician/Patient: 1-none; 2-patient only; 3-clinician only; 4-both
#'
#' @param dose dose level
#' @param scenario scenario table
#'
#' @return
#' @export
#'
#' @examples
Rand <- function(dose, scenario)
{
  prob.c <- scenario$Clinician[dose]
  prob.p <- scenario$Patient[dose]
  prob.cp <- scenario$Clinician.Patient[dose]
  p00 <- 1 - prob.cp
  p01 <- prob.cp - prob.c
  p10 <- prob.cp - prob.p
  p11 <- max(1 - p00 - p01 - p10, 0)
  # if (p00 < 0) cat("p00", p00, "\n")
  # if (p01 < 0) cat("p01", p01, "\n")
  # if (p10 < 0) cat("p10", p10, "\n")
  # if (p11 < 0) cat("p11", p11, "\n")

  x <- which(rmultinom(1, 1, prob = c(p00, p01, p10, p11)) == 1)

  return(x)
}

#' Simple Encoding
#'
#' @param x toxicity code: 1-4
#'
#' @return
#' @export
#'
#' @examples
EncodingSimple <- function(x)
{
  if (x <= 2)
  {
    return(0)
  } else
  {
    return(1)
  }
}

#' Marginal Encoding
#'
#' @param x toxicity code: 1-4
#'
#' @return
#' @export
#'
#' @examples
EncodingMarginal <- function(x)
{
  if (x == 1)
  {
    return(0)
  } else if (x == 2)
  {
    return(1)
  } else if (x == 3)
  {
    return(2)
  } else
  {
    return(c(1, 2))
  }
}

#' Joint Encoding
#'
#' @param x toxicity code: 1-4
#'
#' @return
#' @export
#'
#' @examples
EncodingJoint <- function(x)
{
  if (x == 1)
  {
    return(0)
  } else if (x == 2)
  {
    return (1)
  } else
  {
    return(2)
  }
}

#' Encoding (Simple/Marginal/Joint)
#'
#' @param x toxicity code: 1-4
#' @param type encoding type: simple/marginal/joint
#'
#' @return
#' @export
#'
#' @examples
Encoding <- function(x, type)
{
  switch(type,
         simple = EncodingSimple(x),
         marginal = EncodingMarginal(x),
         joint = EncodingJoint(x),
         paste(type, "endoging currently NOT supported"))
}

#' Toxicity Generator
#'
#' @inheritParams Rand
#' @inheritParams Encoding
#' @param scn scenario table, characterize the joint probability, usually read from file
#'
#' @return
#' @export
#'
#' @examples
ToxGenerator <- function(dose, type, scn)
{
  x <- Rand(dose, scn)
  return(Encoding(x, type))
}


# test --------------------------------------------------------------------

# scenario <- read.csv("scenario/scenario_01.csv")
# Generator <- function(dose)
# {
#   x <- Rand(dose, scenario)
#   return(EncodingJoint(x))
# }
#
# skeleton <- c(0.02, 0.09, 0.25, 0.44, 0.62)
# target <- c(0.50, 0.25)
# n.sim <- 1000
# mtd.sim <- rep(0, n.sim)
# pct <- 10
# for (i.sim in 1 : n.sim)
# {
#   sim <- Simulation(n.dose = 5, n.tox = 2, n.trial = 21, Generator, skeleton, target = target)
#   mtd.sim[i.sim] <- sim$mtd
#   if (i.sim / n.sim >= pct / 100)
#   {
#     cat(pct, "% finished\n")
#     pct <- pct + 10
#   }
# }


