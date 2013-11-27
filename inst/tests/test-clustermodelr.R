library(testthat)

covs = read.delim(system.file("extdata", "example-covariates.txt",
                              package="clustermodelr"))
covs$id = 1:nrow(covs)

meth = as.matrix(read.csv(system.file("extdata", "example-meth.csv",
                            package="clustermodelr"), row.names=1))

test_that("can run single column", {
    meth1 = cbind(meth[,1])
    formula = methylation ~ disease
    expect_that(length(lmr(covs, meth1, formula)), equals(3))
    expect_true("p" %in% names(lmr(covs, meth1, formula)))
    expect_true("coef" %in% names(lmr(covs, meth1, formula)))
})

test_that("can run models", {
    formula = methylation ~ disease
    expect_that(length(clust.lm(covs, meth, formula, combine="liptak")), equals(3))
    expect_that(length(clust.lm(covs, meth, formula, combine="z-score")), equals(3))
    expect_that(length(bumpingr(covs, meth, formula, n_sims=3)), equals(3))
    expect_that(length(clust.lm(covs, meth, formula, gee.idvar="id", gee.corstr="ar")), equals(3))
    expect_that(length(clust.lm(covs, meth, disease ~ 1, skat=TRUE)), equals(3))

    formula = methylation ~ disease + (1|id) + (1|CpG)
    expect_that(length(clust.lm(covs, meth, formula)), equals(3))

})


test_that("can run models on sparse data", {
    formula = methylation ~ case
    cases = gen.correlated(0.23, 20, 4, mean=0.01, sd=0.045)
    controls = gen.correlated(0.23, 20, 4, mean=0.0, sd=0.045)
    meth = rbind(cases, controls)
    colnames(meth) = paste0("probe_", 1:4)
    meth[33, 2] = NaN
    meth[11, 4] = NaN

    covs = data.frame(case=c(rep(1, 20), rep(0, 20)))
    rownames(meth) = rownames(covs) = paste0("sample_", 1:40)
    res = clust.lm(covs, meth, formula, combine="liptak")
    expect_true(!is.na(res$p))
    res = combiner(covs, meth, formula)
    expect_true(!is.na(res$p))
    bres = bumpingr(covs, meth, formula)
    expect_that(bres$coef, equals(res$coef))

    res = combiner(covs, meth, formula, combine.fn=zscore.combine)
    expect_true(!is.na(res$p))

})


cprint = function(...) write(..., stdout())

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
    df = mclust.lm.X(covs, meth, formula, X, testing=TRUE, gee.corstr="ar", gee.clustervar="id")
    print(head(df[order(as.numeric(df$p)),], n=5))

    cprint("\nbumping")
    formula = methylation ~ disease #+ (1|id) + (1|CpG)
    df = mclust.lm.X(covs, meth, formula, X, testing=TRUE, bumping=TRUE)
    print(head(df[order(as.numeric(df$p)),], n=5))

    cprint("\nliptak")
    formula = methylation ~ disease #+ (1|id) + (1|CpG)
    dfl = mclust.lm.X(covs, meth, formula, X, testing=TRUE, combine="liptak")
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

