\name{mclust.lm.X}
\alias{mclust.lm.X}
\title{dispatch to one of the implemented cluster methods. potentially reading
the covariates from a file and parallelizing.}
\usage{
  mclust.lm.X(covs, meth, formula, X, gee.corstr = NULL,
    ..., mc.cores = 4)
}
\arguments{
  \item{covs}{covariate data.frame containing the terms in
  formula except "methylation" which is added
  automatically}

  \item{meth}{a list of matrices of correlated data or a
  single methylation matrix}

  \item{formula}{an R formula containing "methylation"}

  \item{X}{a matrix with columns matching those in meth.
  n_probes X n_samples.  Each row is tested by modifying
  \code{formula} so that it becomes the independent
  variable in the model and tested against methylation.}

  \item{gee.corstr}{if specified, the the corstr arg to
  geeglm.}

  \item{mc.cores}{the number of processors to use if meths
  is a list of matrices to test.}

  \item{...}{arguments sent to \code{\link{clust.lm}}}
}
\description{
  This method implements a sorted of methyl-eQTL with a
  formula specified as:
}
\details{
  \code{methylation ~ disease + age}

  each row in \code{X} is inserted into the model and
  tested so the model would be:

  \code{methylation ~ X[irow,] + disease + age}

  and the reported coefficent and p-value are from the
  X[irow,] covariate. This allows one to test a number of
  expression probes against a (number of) cluster of
  correlated methylation probes. Though we could also use
  this to test, for example a set of methylation probes
  against every OTU in a microbiome study. In this way, we
  could find DMRs related to the microbiome.

  See \code{\link{clust.lm}}
}

