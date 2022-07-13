library(dearseq)
context("Gene Set Analysis with dgsa_seq wrapper")

test_that("Gene sets with no genes measure trigger warnings", {
  rm(list = ls())
  n <- 200
  r <- 12
  t <- matrix(rep(1:3), 4, ncol=1, nrow=r)
  sigma <- 0.4
  b0 <- 1
  b1 <- 0 #under the null
  y.tilde <- b0 + b1*t + rnorm(r, sd = sigma)
  y <- t(matrix(rnorm(n*r, sd = sqrt(sigma*abs(y.tilde))), ncol=n, nrow=r) +
           matrix(rep(y.tilde, n), ncol=n, nrow=r))
  x <- matrix(1, ncol=1, nrow=r)
  gs <- list(nrow(y) + 1:10)
  expect_warning(dgsa_seq(exprmat = y, covariates = x, variables2test = t,
                          genesets = gs, cov_variables2test_eff = matrix(1),
                          sample_group = rep(1:4, each=3),
                          which_test = "permutation",
                          which_weights = "none", preprocessed = TRUE))
  expect_warning(dgsa_seq(exprmat = y, covariates = x, variables2test = t,
                          genesets = gs, cov_variables2test_eff = matrix(1),
                          sample_group = rep(1:4, each=3),
                          which_test = "asymptotic",
                          which_weights = "none", preprocessed = TRUE))
})

test_that("NA runs", {
  rm(list = ls())
  n <- 1000 #number of genes
  nr=5 #number of measurements per subject (grouped data)
  ni=50 #number of subjects
  r <- nr*ni #number of measurements
  t <- matrix(rep(1:nr), ni, ncol=1, nrow=r) # the variable to be tested
  sigma <- 0.5
  b0 <- 1
  
  #under the null:
  b1 <- 0
  
  #create the matrix of gene expression
  y.tilde <- b0 + b1*t + rnorm(r, sd = sigma)
  y <- t(matrix(rnorm(n*r, sd = sqrt(sigma*abs(y.tilde))), ncol=n, nrow=r) +
           matrix(rep(y.tilde, n), ncol=n, nrow=r))
  
  # add NA
  y[10,8] = NA
  
  #no covariates
  x <- matrix(1, ncol=1, nrow=r)
  
  #run test
  #asymptotic test with preprocessed grouped data
  
  expect_warning(res <- dear_seq(exprmat=y, covariates=x, variables2test=t,
                                 sample_group=rep(1:ni, each=nr),
                                 which_test = "asymptotic",
                                 na.rm_dearseq = TRUE,
                                 which_weights='loclin', preprocessed=TRUE))
  
})

test_that("Returned as many p-values as there are genesets", {
  rm(list = ls())
  n <- 200
  r <- 12
  t <- matrix(rep(1:3), 4, ncol=1, nrow=r)
  sigma <- 0.4
  b0 <- 1
  b1 <- 0 #under the null
  y.tilde <- b0 + b1*t + rnorm(r, sd = sigma)
  y <- t(matrix(rnorm(n*r, sd = sqrt(sigma*abs(y.tilde))), ncol=n, nrow=r) +
           matrix(rep(y.tilde, n), ncol=n, nrow=r))
  x <- matrix(1, ncol=1, nrow=r)
  res <- dgsa_seq(exprmat = y, covariates = x, variables2test = t,
                  genesets = list(1:10, 11:20),
                  cov_variables2test_eff = matrix(1),
                  sample_group = rep(1:4, each=3), which_test = "asymptotic",
                  which_weights = "none", preprocessed = TRUE)
  expect_length(res$pvals$rawPval, n = length(res$genesets))
  expect_length(res$pvals$adjPval, n = length(res$genesets))
})


