library(testthat)
library(MetaboMate)
path=system.file("extdata/", package = "MetaboMate")
print(path)
readBruker(path)


qcs=spec.quality(X, ppm, ppm.noise=c(14,14.2), plot=F)
qcs_e1=evaluate_promise(spec.quality(X, ppm[-1], ppm.noise=c(14,14.1), plot=F))



test_check("MetaboMate")
