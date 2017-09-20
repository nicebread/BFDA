expected.ES = 0.2

type="correlation"
n.min=20
n.max=350
stepsize=10
alternative="directional"
boundary=Inf
B=5000
design="sequential"
verbose=TRUE
cores=NA

BFDA <- function(expected.ES, type="t.between", n.min=20, n.max=200, stepsize=10, alternative="directional", boundary=Inf, B=5000, design="sequential", verbose=TRUE, cores=NA)

	cat(paste0("
Computing sequential BFDA for an expected ES of ", expected.ES, " for a ", type, " design with ", alternative, " hypothesis.
Sample sizes are computed from n.min = ", n.min, " to n.max = ", n.max, " with a stepsize of ", stepsize, ".
Number of Monte Carlo simulations: ", B, "
	"))

	if (is.na(cores)) {
		cores <- detectCores()-1
		cat(paste0("Using ", cores, " of ", detectCores(), " cores for simulations."))
	}

	sim.H1 <- BFDA.sim(expected.ES=expected.ES, type=type, n.min=n.min, n.max=n.max, alternative=alternative, boundary=boundary, B=B, verbose=verbose, cores=cores)

	sim.H0 <- BFDA.sim(expected.ES=0, type=type, n.min=n.min, n.max=n.max, alternative=alternative, boundary=boundary, B=B, verbose=verbose, cores=cores)
	
	BFDA.analyze(sim.H1, design="sequential", n.min=20, n.max=300, boundary=10)
	
	invisible(return()list(sim.H1, sim.H0))
}