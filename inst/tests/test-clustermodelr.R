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
})

test_that("can run models", {
    formula = methylation ~ disease
    expect_that(length(clust.lm(covs, meth, formula, liptak=TRUE)), equals(3))
    expect_that(length(clust.lm(covs, meth, formula, bumping=TRUE)), equals(3))
    expect_that(length(clust.lm(covs, meth, formula, gee.idvar="id", gee.corstr="ar")), equals(3))
    expect_that(length(clust.lm(covs, meth, disease ~ 1, skat=TRUE)), equals(3))
    formula = methylation ~ disease + (1|id) + (1|CpG)
    expect_that(length(clust.lm(covs, meth, formula)), equals(3))

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

