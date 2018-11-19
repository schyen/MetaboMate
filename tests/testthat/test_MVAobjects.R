context("Panoptikum: read bruker, dummy var, MVA objects")

#### read bruker
path=system.file("extdata/", package = "MetaboMate")
readBruker(path)

test_that('read-in Bruker files', {
  expect_equal(length(which(is.na(X))), 0)
  expect_equal(length(which(is.na(ppm))), 0)
  expect_equal(nrow(X), 30)
  expect_equal(nrow(meta), 30)
  expect_equal(length(ppm), ncol(X))
})

#### dummy var
Y=factor(rep(LETTERS[1:3], each=5))
out=MetaboMate:::create_dummy_Y(Y)[[1]]

test_that('output dummy var from vector', {
  expect_equal(ncol(out), 3)
  expect_equal(length(unique(out[,1])), 2)
  expect_equal(length(unique(out[,2])), 2)
  expect_equal(length(unique(out[,3])), 2)
})

#### S4
s <- new("PCA_MetaboMate")
test_that('Number of slots in PCA object', {
  expect_equal(length(slotNames(s)), 5)
})

s <- new("OPLS_MetaboMate")
test_that('Number of slots in OPLS object', {
  expect_equal(length(slotNames(s)), 22)
})
