\docType{package}
\name{clustermodelr}
\alias{clustermodelr}
\alias{clustermodelr-package}
\title{Model clustered, correlated data}
\description{
  clustermodelr provides a consistent, simple interface to
  model correlated data using a number of different methds:
  \describe{ \item{GEE:}{Generalized Estimating Equations
  with all correlation structures available from
  \code{geepack}. \code{\link{geer}}} \item{mixed-effect
  model:}{Mixed effect model in \link[lme4]{lme4} syntax
  \code{\link{mixed_modelr}}} \item{Liptak:}{Calculates the
  p-value for each entry in the cluster then combines the
  p-values adjusting for correlation.
  \code{\link{stouffer_liptakr}}} \item{bumping:}{something
  like bump-hunting but takes a putative "bump" and
  repeatedly compares coefficients of estimated covariates
  to the observed to assign significance.
  \code{\link{bumpingr}}} \item{SKAT:}{SKAT already accepts
  a matrix to test a null model. This just provides an
  interface that matches the rest of the functions in this
  package \code{\link{skatr}}} }
}
\details{
  Each of these functions will accept a formula like:

  \code{methylation ~ disease + age}

  (with a random intercept for mixed_modelr) where
  \code{methylation} need not be methylation values, but is
  assumed to be a matrix of correlated values.

  For each of these functions, the \strong{return value}
  will be a vector of:

  \code{c(covariate, p, coef.estimate)}

  where the covariate is taken as the first element on the
  RHS of the formula so \emph{disease} in the formula
  above.
}

