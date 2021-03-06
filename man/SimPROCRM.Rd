% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pro.R
\name{SimPROCRM}
\alias{SimPROCRM}
\title{Multi-Trial Wrapper for SimOnePROCRM()}
\usage{
SimPROCRM(n.sim, n, scenario.table, skeleton, target, method,
  seed = NULL)
}
\arguments{
\item{n.sim}{# of simulations/trials}

\item{n}{# of subjects (sample size) in each trial, should be a multiple of cohort.size}

\item{scenario.table}{scenario table (toxicity rates) for C, P, and CoP}

\item{skeleton}{mapped skeleton for dose levels, a vector of length = # of dose levels (K), should be between 0 and 1}

\item{target}{targets for (cumulative) toxicity rates, a vector of length = # of toxicity levels (L)}

\item{method}{one of "mar", "jntmar", "jntjnt"}

\item{seed}{seed for RNG}
}
\description{
Multi-Trial Wrapper for SimOnePROCRM()
}
