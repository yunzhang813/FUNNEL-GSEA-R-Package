test_wMWUTest <- function() {
  data(sampleData)
  tst <- wMWUTest(index, stat, weight=weight, correlation=0.1)
  checkEqualsNumeric(as.vector(tst), c(0.97293933,0.02849134), tolerance=1.0e-4)
}
