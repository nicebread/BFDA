load("../finalSims/BPA.0.one.RData")
load("../finalSims/BPA.0.two.RData")

BPA.analysis(BPA.0.two, boundary=6)
BPA.analysis(BPA.0.one, boundary=6)

BPA.analysis(BPA.0.one, boundary=6, n.max=2000)