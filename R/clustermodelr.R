#' Model clustered, correlated data
#'
#' clustermodelr provides a consistent, simple interface to model correlated
#' data using a number of different methds:
#' \describe{ 
#' \item{GEE:}{Generalized Estimating Equations with all correlation structures
#'            available from \code{geepack}. \code{\link{geer}}}
#' \item{mixed-effect model:}{Mixed effect model in \link[lme4]{lme4} syntax
#'                           \code{\link{mixed_modelr}}}
#' \item{Liptak:}{Calculates the p-value for each entry in the cluster then
#'                combines the p-values adjusting for correlation.
#'                \code{\link{stouffer_liptakr}}}
#' \item{bumping:}{something like bump-hunting but takes a putative "bump" and
#'                repeatedly compares coefficients of estimated covariates to the observed
#'                to assign significance. \code{\link{bumpingr}}}
#' \item{SKAT:}{SKAT already accepts a matrix to test a null model. This just
#'             provides an interface that matches the rest of the functions in
#'             this package \code{\link{skatr}}}
#' }
#' 
#' @details
#' Each of these functions will accept a formula like:
#' 
#' \code{methylation ~ disease + age}
#'
#' (with a random intercept for mixed_modelr)
#' where \code{methylation} need not be methylation values, but is assumed to be
#' a matrix of correlated values.
#' 
#' For each of these functions, the \strong{return value} will be a vector of:
#' 
#' \code{c(covariate, p, coef.estimate)} 
#' 
#' where the covariate is taken as the first
#' element on the RHS of the formula so \emph{disease} in the formula above.
#'
#' @docType package
#' @name clustermodelr



suppressPackageStartupMessages(library(limma, quietly=TRUE))


#' Run lm on a single site
#' 
#' Unlike most of the function in this package, this function is used on a
#' single site.
#' @param covs covariate data.frame containing the terms in formula
#'        except "methylation" which is added automatically
#' @param methylation a single column matrix or a vector the same length
#'        as \code{nrow(covs)}
#' @param formula an R formula containing "methylation"
#' @return \code{list(covariate, p, coef)} where p and coef are for the coefficient
#'         of the first term on the RHS of the model.
#' @export
lmr = function(covs, methylation, formula){
    covs$methylation = methylation
    s = summary(lm(formula, covs))$coefficients
    covariate = rownames(s)[2]
    row = s[2,]
    list(covariate=covariate, p=row[['Pr(>|t|)']], coef=row[['Estimate']])
}

#' Perform the stouffer-liptak correction on a set of correlated pvalues
#'
#' Accepts a list of p-values and a correlation matrix and returns a
#' single p-value that combines the p-values accounting for their
#' correlation
#'
#' @param pvalues a vector of pvalues
#' @param sigma a matrix of shape \code{nrow=ncol=length(pvalues)}
#' @return combined p-value
#' @export
stouffer_liptak = function(pvalues, sigma){
    qvalues = qnorm(pvalues, mean=0, sd=1, lower.tail=TRUE)
    C = chol(sigma)
    Cm1 = solve(C) # C^-1
    qvalues = Cm1 %*% qvalues # Qstar = C^-1 * Q
    Cp = sum(qvalues) / sqrt(length(qvalues))
    pnorm(Cp, mean=0, sd=1, lower.tail=TRUE)
}


#' Run lm on each column in a cluster and combine p-values with the 
#' Stouffer-Liptak method
#' 
#' @param covs covariate data.frame containing the terms in formula
#'        except "methylation" which is added automatically
#' @param meth a matrix of correlated data.
#' @param formula an R formula containing "methylation"
#' @param cor.method either "spearman" or "pearson"
#' @return \code{list(covariate, p, coef)} where p and coef are for the coefficient
#'         of the first term on the RHS of the model.
#' @export
stouffer_liptakr = function(covs, meth, formula, cor.method="spearman"){
    # TODO: what if missing data in covariates.
    # set up another method that runs each in lm() and pulls the coefficents
    # and pvalues
    library(limma)
    covs$methylation = 1 # 
    sigma = cor(meth, method=cor.method)
    meth = t(meth)

    mod = model.matrix(formula, covs)
    covariate = colnames(mod)[1 + as.integer(colnames(mod)[1] == "(Intercept)")]

    fit = eBayes(lmFit(meth, mod))
    beta.orig = coefficients(fit)[,covariate]
    pvals = topTable(fit, coef=covariate)[,"P.Value"]
    beta.ave = sum(beta.orig) / length(beta.orig)

    p = stouffer_liptak(pvals, sigma)
    return(list(covariate=covariate, p=p, coef=beta.ave))
}

# for bumping
permute.residuals = function(mat, mod, mod0, iterations=100, p_samples=1, mc.cores=10){
    stopifnot(nrow(mod) == ncol(mat))

    reduced_lm = lmFit(mat, mod0)
    reduced_residuals = residuals(reduced_lm, mat)
    reduced_fitted = fitted(reduced_lm)

    fit = lmFit(mat, mod)

    coef.name = setdiff(colnames(mod), colnames(mod0))
    beta.orig = coefficients(fit)[,coef.name]

    rm(reduced_lm, fit); gc()
    nc = ncol(reduced_residuals)

    beta.list = mclapply(1:iterations, function(ix){
        mat_sim = reduced_fitted + reduced_residuals[,sample(1:nc)]
        ifit = lmFit(mat_sim, mod)
        icoef = coefficients(ifit)[,coef.name]
        w = ifit$sigma
        # get names as integer positions:
        names(icoef) = 1:length(icoef)
        names(w) = 1:length(w)
        sum.lowess(icoef, w)
    }, mc.cores=mc.cores)
        
    beta.sum = rep(0, n=iterations)
    for(i in 1:iterations){
        beta.sum[i] = beta.list[[i]]
    }
    beta.sum
}

sum.lowess = function(icoefs, weights, span=0.2){
    if(length(icoefs) < 3){ return(sum(icoefs)) }
    res = try(limma::loessFit(icoefs, as.integer(names(icoefs)),
                              span=span, weights=weights), silent=TRUE)
    if(class(res) == "try-error") return(sum(icoefs))
    return(sum(res$fitted))
}


#' Run a local bump-hunting algorithm
#' 
#' This performs a similar task to the Bump-Hunting algorithm, but here the
#' \code{meth} argument is a putative bump. The residuals of the null model
#' are shuffled added back to the null model and the beta's of that simulated
#' data are repeatedly stored and finally compare to the observed coefficient.
#' Due to the shufflings, this is much slower than the other functions in this
#' package.
#' 
#' @param covs covariate data.frame containing the terms in formula
#'        except "methylation" which is added automatically
#' @param meth a matrix of correlated data.
#' @param formula an R formula containing "methylation"
#' @param n_sims this is currently used as the minimum number of shuffled data
#'        sets to compare to. If the p-value is low, it will do more shufflings
#' @return \code{list(covariate, p, coef)} where p and coef are for the coefficient
#'         of the first term on the RHS of the model.
#' @export
bumpingr = function(covs, meth, formula, n_sims=100){
    suppressPackageStartupMessages(library('parallel', quietly=TRUE))
    covs$methylation = 1 # for formula => model.matrix

    mod = model.matrix(formula, covs)
    covariate = colnames(mod)[1 + as.integer(colnames(mod)[1] == "(Intercept)")]
    mod0 = mod[,!colnames(mod) == covariate, drop=FALSE]

    sim_beta_sums = permute.residuals(meth, mod, mod0, iterations=n_sims)
    stopifnot(length(sim_beta_sums) == n_sims)

    fit = lmFit(meth, mod)
    w = fit$sigma

    icoef = coefficients(fit)[,covariate]
    # get names as integer positions:
    names(icoef) = 1:length(icoef)
    names(w) = 1:length(w)
    beta_sum = sum.lowess(icoef, w)

    raw_beta_sum = sum(coefficients(fit)[,covariate])
    ngt = sum(abs(sim_beta_sums) >= abs(beta_sum))
    # progressive monte-carlo: only do lots of sims when it has a low p-value.
    if(ngt < 4 & n_sims == 100) return(bumpingr(covs, meth, formula, 2000))
    if(ngt < 10 & n_sims == 2000) return(bumpingr(covs, meth, formula, 5000))
    if(ngt < 10 & n_sims == 5000) return(bumpingr(covs, meth, formula, 15000))
    pval = (1 + ngt) / (1 + n_sims)
    return(list(covariate=covariate, p=pval, coef=raw_beta_sum / nrow(meth)))
}

#' Use Generalized Estimating Equations to assign significance to a cluster
#' of data.
#' 
#' @param covs covariate data.frame containing the terms in formula
#'        except "methylation" which is added automatically
#' @param formula an R formula containing "methylation"
#' @param idvar idvar sent to \code{geepack::geeglm}
#' @param corstr the corstr sent to \code{geepack::geeglm}
#' @return \code{list(covariate, p, coef)} where p and coef are for the coefficient
#'         of the first term on the RHS of the model.
#' @export
geer = function(covs, formula, idvar="CpG", corstr="ex"){
    # assume it's already sorted by CpG, then by id.
    if(idvar != "CpG" && corstr == "ar"){
        covs = covs[order(covs[,idvar], covs$CpG),]
    }
    suppressPackageStartupMessages(library('geepack', quietly=TRUE))
    stopifnot(!is.null(idvar))

    # NOTE, both of these are required as geeglm needs to deparse and
    # R CRAN has to make sure clustervar exists.
    clustervar = covs$clustervar = covs[,idvar]
    # can't do logistc with idvar of id, gives bad results for some reason
    s = summary(geeglm(formula, id=clustervar, data=covs, corstr=corstr))
    mm = model.matrix(formula, covs)
    covariate = colnames(mm)[1 + as.integer(colnames(mm)[1] == "(Intercept)")]
    row = s$coefficients[covariate,]
    return(list(covariate=covariate, p=row[['Pr(>|W|)']], coef=row[['Estimate']]))
}

#geer(read.csv('tt.csv'), methylation ~ disease, "id", "ex")

#' Use mixed-effects model in lme4 syntax to associate a covariate with a
#' cluster of data.
#' 
#' An example model would look like:
#'    methylation ~ disease + age + gender + (1|CpG) + (1|id)
#' To determine the associate of disease and methylation and allowing for
#' random intercepts by CpG site and by sample id.
#' 
#' @param covs covariate data.frame containing the terms in formula
#'        except "methylation" which is added automatically
#' @param formula an R formula containing "methylation"
#' @return \code{list(covariate, p, coef)} where p and coef are for the coefficient
#'         of the first term on the RHS of the model.
#' @export
mixed_modelr = function(covs, formula){
    suppressPackageStartupMessages(library('lme4', quietly=TRUE))
    suppressPackageStartupMessages(library('multcomp', quietly=TRUE))
    # automatically do logit regression.
    m = lmer(formula, covs)
    # take the first column unless it is intercept
    covariate = names(fixef(m))[1 + as.integer(names(fixef(m))[1] == "(Intercept)")]
    s = summary(glht(m, paste(covariate, "0", sep=" == ")))
    return(list(covariate=covariate, p=s$test$pvalues[[1]], coef=s$test$coefficients[[1]]))
}

#' Use SKAT to associate a covariate with a cluster of data.
#' 
#' An example model would look like:
#'    disease ~ 1
#' And the result would be testing if adding the correlated matrix of
#' data (with all the assumptions of SKAT) improves that null model.
#' 
#' @param covs covariate data.frame containing the terms in formula
#'        except "methylation" which is added automatically
#' @param meth a matrix of correlated data.
#' @param formula an R formula containing "methylation"
#' @param r.corr list of weights between kernel and rare variant test.
#' @return \code{list(covariate, p, coef)} where p and coef are for the coefficient
#'         of the first term on the RHS of the model.
#' @export
skatr = function(covs, meth, formula, r.corr=c(0.00, 0.015, 0.06, 0.15)){
    suppressPackageStartupMessages(library('SKAT', quietly=TRUE))
    covariate = all.vars(formula)[1]

    capture.output(obj <- SKAT_Null_Model(formula, out_type="D", data=covs))
    #sk <- SKAT(meth, obj, is_check_genotype=FALSE, method="davies", r.corr=0.6, kernel="linear")
    sk <- SKAT(meth, obj, is_check_genotype=FALSE, method="optimal.adj", kernel="linear",
            r.corr=r.corr)
    #sk <- SKAT(meth, obj, is_check_genotype=TRUE, method="optimal.adj", kernel="linear.weighted", weights.beta=c(1, 10))
    #sk <- SKAT(meth, obj, is_check_genotype=TRUE, method="optimal.adj", kernel="linear")
    #sink()
    return(list(covariate=covariate, p=sk$p.value, coef=NaN))
}

#' convert data to long format so that \code{covs} is replicated once for
#' each column in \code{meth}
#'
#' @param covs data.frame of covariates
#' @param meth matrix of methylation with same number of rows as \code{covs}
#' @export
expand.covs = function(covs, meth){
    if(!"id" %in% colnames(covs)) covs$id = 1:nrow(covs)
    n_samples = nrow(covs)
    meth = as.matrix(meth)
    stopifnot(nrow(meth) == n_samples)
    # e.g. meth is 68 patients * 4 CpGs
    #      covs is 68 patients * 5 covariates
    # need to replicated covs 4 times (1 per CpG)
    covs = covs[rep(1:nrow(covs), ncol(meth)),, drop=FALSE]
    cpgs = 1:ncol(meth)
    dim(meth) = NULL
    covs$methylation = meth
    covs$CpG = rep(cpgs, each=n_samples) # 1 1 1, 2 2 2, etc since CpG's are grouped.
    covs
}

#' dispatch to one of the implemented cluster methods
#' 
#' For every method except mixed_model, one or more of the arguments
#' must be specified. To run a linear model, simply send the formula
#' in lme4 syntax
#' 
#' @param covs covariate data.frame containing the terms in formula
#'        except "methylation" which is added automatically
#' @param meth a matrix of correlated data.
#' @param formula an R formula containing "methylation"
#' @param gee.corstr if specified, the the corstr arg to geeglm.
#'        gee.idvar must also be specified.
#' @param gee.idvar if specified, the cluster variable to geeglm
#' @param bumping if true then the bumping algorithm is used.
#' @param liptak if true then run the model on each probe in \code{meth}
#'        and perform the stouffer-liptak correction on the p-values
#' @param skat use the SKAT method to test associated. In this case, the
#'        model will look like: \code{disease ~ 1} and it will be tested
#'        against the methylation matrix
#' @export
clust.lm = function(covs, meth, formula,
                    gee.corstr=NULL, gee.idvar=NULL,
                    bumping=FALSE, liptak=FALSE, skat=FALSE){

    formula = as.formula(formula)

    if(ncol(meth) == 1 || is.vector(meth)){
        # just got one column, so we force it to use a linear model
        # remove random effects terms:
        lhs = grep("|", attr(terms(formula), "term.labels"), fixed=TRUE, value=TRUE, invert=TRUE)
        lhs = paste(lhs, collapse=" + ")
        formula = as.formula(paste("methylation", lhs, sep=" ~ "))
        return(lmr(covs, meth, formula))
    }

    # we assume there is one extra column for each CpG
    rownames(meth) = rownames(covs)

    if(bumping){ # wide
        return(bumpingr(covs, t(meth), formula))
    }
    if(skat){ # wide
        return(skatr(covs, meth, formula))
    }
    if(liptak){ # wide
        return(stouffer_liptakr(covs, meth, formula))
    }

    ###########################################
    # GEE and mixed models require long format.
    ###########################################
    covs = expand.covs(covs, meth) # TODO: make this send just the nrow, ncol
    is.mixed.model = any(grepl("|", attr(terms(formula), 'term.labels'), fixed=TRUE))
    # mixed-model
    if (is.null(gee.corstr)){
        stopifnot(is.mixed.model)
        return(mixed_modelr(covs, formula))
    # GEE
    } else if (!is.null(gee.corstr)){
        stopifnot(!is.null(gee.idvar))
        return(geer(covs, formula, idvar=gee.idvar, corstr=gee.corstr))
    # limma
    } else {
        # TODO this goes in the matrix section above and uses
        # duplicateCorrelation
        stop()
    }

}

#' used to communicate quickly from python
#' @export
#' @param bin.file file with binary data
read.bin = function(bin.file){
    conn = file(bin.file, 'rb')
    n_sites = readBin(conn, what=integer(), size=8, n=1)
    l = list()
    for(i in 1:n_sites){
        mdims = readBin(conn, what=integer(), size=8, n=2)
        nrow = mdims[1]
        ncol = mdims[2]
        dat = readBin(conn, what=numeric(), size=8, n=nrow * ncol)
        l[[i]] = matrix(dat, nrow=nrow, ncol=ncol, byrow=TRUE)
    }
    close(conn)
    l
}

#' dispatch to one of the implemented cluster methods. potentially reading
#' the covariates from a file and parallelizing.
#' 
#' See \code{\link{clust.lm}}
#' 
#' @param covs covariate data.frame containing the terms in formula
#'        except "methylation" which is added automatically
#' @param meths a list of matrices of correlated data.
#' @param formula an R formula containing "methylation"
#' @param gee.corstr if specified, the the corstr arg to geeglm.
#' @param mc.cores the number of processors to use if meths is a list of
#'        matrices to test.
#' @param ... arguments sent to \code{\link{clust.lm}}
#' @export
mclust.lm = function(covs, meths, formula, gee.corstr=NULL, ..., mc.cores=4){
    if(is.character(covs)) covs = read.csv(covs)

    # its a single entry, not list of matrices that we can parallelize
    if(is.matrix(meths) || is.data.frame(meths)){
        res = (clust.lm(covs, meths, formula, gee.corstr=gee.corstr, ...))
        return(data.frame(res))
    }

    suppressPackageStartupMessages(library('data.table', quietly=TRUE))
    suppressPackageStartupMessages(library('parallel', quietly=TRUE))

    cluster_ids = 1:length(meths)
    results = mclapply(cluster_ids, function(cs){
        res = try(clust.lm(covs, meths[[cs]], formula, gee.corstr=gee.corstr, ...))
        if(!inherits(res, "try-error")){
            res$cluster_id = cs
            return(res)
        }
        return(list(covariate=NA, p=NaN, coef=NaN, cluster_id=cs))
    }, mc.cores=mc.cores)
    results = rbindlist(results)
    rownames(results) = cluster_ids
    results
}


if(FALSE){
    covs = read.delim("clustercorr/tests/example-covariates.txt")
    covs$id = 1:nrow(covs)
    meth = read.csv('clustercorr/tests/example-meth.csv', row.names=1)

    #  check with only a single value
    meth = cbind(meth[,1])
    print(ncol(meth))

    print(mclust.lm(covs, meth, methylation ~ disease + (1|id) + (1|CpG)))

    print('liptak')
    print(mclust.lm(covs, meth, methylation ~ disease, liptak=TRUE))

    print(mclust.lm(covs, meth, methylation ~ disease + (1|id)))
    #print(clust.lm(covs, methylation ~ gene.E, gee.idvar="id", gee.corstr="ex"))
    print(mclust.lm(covs, meth, methylation ~ disease, gee.idvar="id", gee.corstr="ex"))
    print(mclust.lm(covs, meth, methylation ~ disease, gee.idvar="id", gee.corstr="ar"))
    print('bumping')
    print(mclust.lm(covs, meth, methylation ~ disease, bumping=TRUE))
    print(clust.lm(covs, as.matrix(meth), disease ~ 1, skat=TRUE))
}

#' read a matrix of numeric values with the first column as the row.names
#'
#' @param fname the file name of the Xpression dataset to read.
#' @export
readX = function(fname){
    X = as.matrix(read.delim(gzfile(fname), row.names=1))
    rownames(X) = gsub("-|:| ", ".", as.character(rownames(X)), perl=TRUE)
    X
} 

#' dispatch to one of the implemented cluster methods. potentially reading
#' the covariates from a file and parallelizing.
#'
#' This method implements a sorted of methyl-eQTL with a formula specified
#' as:
#' 
#'    \code{methylation ~ disease + age}
#'
#' each row in \code{X} is inserted into the model and tested so the model
#' would be:
#'
#'    \code{methylation ~ X[irow,] + disease + age}
#'
#' and the reported coefficent and p-value are from the X[irow,] covariate.
#' This allows one to test a number of expression probes against a (number
#' of) cluster of correlated methylation probes. Though we could also use
#' this to test, for example a set of methylation probes against every OTU
#' in a microbiome study. In this way, we could find DMRs related to the
#' microbiome.
#' 
#' See \code{\link{clust.lm}}
#' 
#' @param covs covariate data.frame containing the terms in formula
#'        except "methylation" which is added automatically
#' @param meth a list of matrices of correlated data or a single methylation
#'        matrix
#' @param formula an R formula containing "methylation"
#' @param X a matrix with columns matching those in meth. n_probes X n_samples.
#'        Each row is tested by modifying \code{formula} so that it becomes the
#'        independent variable in the model and tested against methylation.
#' @param gee.corstr if specified, the the corstr arg to geeglm.
#' @param mc.cores the number of processors to use if meths is a list of
#'        matrices to test.
#' @param ... arguments sent to \code{\link{clust.lm}}
#' @export
mclust.lm.X = function(covs, meth, formula, X, gee.corstr=NULL, ..., mc.cores=4){
    library(parallel)
    library(data.table)
    formula = as.formula(formula)
    if(is.character(covs)) covs = read.csv(covs)

    # if calling repeatedly, should be subsets of the expression matrix that are close to
    # (presumably) the methylation matrix being tested.
    if(is.character(X)){
        X = readX(X)
    }

    mc.cores = min(mc.cores, ncol(X))

    rnames = rownames(X)

    stopifnot(nrow(covs) %% ncol(X) == 0)
    n_each = nrow(covs) / ncol(X)

    # need this for when X_locs is not specified since we never readi
    # in the array in python

    # get a + b + c from y ~ a + b + x
    rhs = as.character(formula)[length(as.character(formula))]
    lhs = as.character(formula)[2]
    irows = 1:nrow(X)
    stopifnot(n_each >= 1)

    results = mclapply(irows, function(irow){
        X.row = rep(t(X[irow,]), n_each)
        covs2 = covs # make a copy so we dont end up with huge covs
        # add the expression column to the dataframe.
        covs2[,rnames[irow]] = X.row
        sformula = sprintf("%s ~ %s + %s", lhs, rnames[irow], rhs)
        # call with 1 core since we're already parallel here.
        res = mclust.lm(covs2, meth, as.formula(sformula),
                           gee.corstr=gee.corstr, ..., mc.cores=1)
        res$X = rnames[irow]
        res$model = sformula
        res
    }, mc.cores=mc.cores)
    rbindlist(results)
}

cprint = function(...) write(..., stdout())

#' generate correlated data
#'
#' @param rho numeric correlation value between 0 and 1
#' @param n_samples generate data for this many samples
#' @param n_sites generate data for this many sites (CpGs)
#' @param mean sent to \code{rnorm}
#' @param sd sent to \code{rnorm}
#' @return mat n_samples * n_sites matrix where \code{cor(mat[,1], mat[,2])} is
#'         on average equal to \code{rho}
#' @export
gen.correlated = function(rho, n_samples=100, n_sites=4, mean=0, sd=1){
    X = matrix(rnorm(n_samples * n_sites, mean=mean, sd=sd), nrow=n_samples)
    sigma = diag(n_sites)
    sigma <- rho ^ abs(row(sigma)-col(sigma))
    X %*% chol(sigma)
}


test_X = function(){
    covs = read.delim("clustercorr/tests/example-covariates.txt")
    covs$id = 1:nrow(covs)
    meth = read.csv('clustercorr/tests/example-meth.csv', row.names=1)
    #covs = covs[covs$cluster_set == 1,]
    X = read.delim(gzfile('clustercorr/tests/example-expression.txt.gz'), row.names=1)

    cprint("\nmixed-effects model")
    formula = methylation ~ disease + (1|id) + (1|CpG)
    df = mclust.lm.X(covs, meth, formula, X, testing=TRUE)
    print(head(df[order(as.numeric(df$p)),], n=5))

    cprint("\nGEE")
    formula = methylation ~ disease #+ (1|id) + (1|CpG)
    df = mclust.lm.X(covs, meth, formula, X, testing=TRUE, gee.corstr="ar", gee.idvar="id")
    print(head(df[order(as.numeric(df$p)),], n=5))

    cprint("\nbumping")
    formula = methylation ~ disease #+ (1|id) + (1|CpG)
    df = mclust.lm.X(covs, meth, formula, X, testing=TRUE, bumping=TRUE)
    print(head(df[order(as.numeric(df$p)),], n=5))

    cprint("\nliptak")
    formula = methylation ~ disease #+ (1|id) + (1|CpG)
    dfl = mclust.lm.X(covs, meth, formula, X, testing=TRUE, liptak=TRUE)
    print(head(dfl[order(dfl$p),], n=5))
    print(dfl[dfl$covariate == "A_33_P3403576",])

    # show that we get the same result (about with the linear model)
    # pvalue is  2.85844757130782e-06 for the clustered approach and
    # 7.88e-07 for looking at a single probe with a linear model in
    # the region. coefficients vary by ~ 0.001.
    probe = "A_33_P3403576"
    covs$X = t(X[probe,])
    ests = c()
    for(cname in colnames(meth)){
        covs$methylation = meth[,cname]
        cprint(paste0("\n", probe))
        s = summary(lm(methylation ~ X + disease, covs))$coefficients
        print(s)
        ests = c(ests, s['X', 'Estimate'])
    }
    print(mean(ests))

}
#test_X()

