load("../finalSims/BFDA.0.one.RData")
load("../finalSims/BFDA.0.two.RData")

BFDA.analyze(BFDA.0.two, boundary=6)
BFDA.analyze(BFDA.0.one, boundary=6)

BFDA.analyze(BFDA.0.one, boundary=6, n.max=2000)