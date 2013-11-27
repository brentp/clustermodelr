#' Perform the stouffer-liptak correction on a set of correlated pvalues
#'
#' Accepts a list of p-values and a correlation matrix and returns a
#' single p-value that combines the p-values accounting for their
#' correlation. Note that if \code{sigma} is  not positive-definitive,
#' then off diagonal elements will be multiplied by 0.9999.
#'
#' @param pvalues a vector of pvalues
#' @param sigma a matrix of shape \code{nrow=ncol=length(pvalues)}
#' @return combined p-value
#' @export
stouffer_liptak.combine = function(pvalues, sigma){
    stopifnot(length(pvalues) == nrow(sigma))
    qvalues = qnorm(pvalues, mean=0, sd=1, lower.tail=TRUE)
    C = try(chol(sigma), silent=TRUE)
    if(inherits(C, "try-error")){
        sigma[row(sigma) != col(sigma)] = sigma[row(sigma) != col(sigma)] * 0.999
        C = chol(sigma, pivot=TRUE)
    }
    Cm1 = solve(C) # C^-1
    qvalues = Cm1 %*% qvalues # Qstar = C^-1 * Q
    Cp = sum(qvalues) / sqrt(length(qvalues))
    pnorm(Cp, mean=0, sd=1, lower.tail=TRUE)
}

#' Combine p-values by the z-score correction.
#'
#' Accepts a list of p-values and a correlation matrix and returns a
#' single p-value that combines the p-values accounting for their
#' correlation. Note that if \code{sigma} is  not positive-definitive,
#' then off diagonal elements will be multiplied by 0.9999.
#'
#' @param pvalues a vector of pvalues
#' @param sigma a matrix of shape \code{nrow=ncol=length(pvalues)}
#' @return combined p-value
#' @export
zscore.combine = function(pvalues, sigma) {
    # from biseq paper
    n_probes = length(pvalues)
    z = mean(qnorm(pvalues, lower.tail=TRUE))
    sz = 1/n_probes * sqrt(n_probes + 2 * sum(sigma[lower.tri(sigma)]))
    pnorm(z/sz, lower.tail=TRUE)
}
