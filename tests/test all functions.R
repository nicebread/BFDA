load("../finalSims/BFDA.0.5.RData")
load("../finalSims/BFDA.0.RData")

BFDA.analyze(BFDA.0.5, boundary=10)
BFDA.analyze(BFDA.0, boundary=10)