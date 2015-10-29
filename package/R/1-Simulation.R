# source("1-Simulation.R", echo=TRUE)

#' Simulate ... TODO
#' @export
#' @import BayesFactor
#' @import dplyr
#' @import doParallel
#'
#' @param n.min Minimum n before optional stopping is started
#' @param n.max Maximum n - if that is reached without hitting a boundary, the run is aborted
#' @param boundary The Bayes factor (resp. its reciprocal) where the run is stopped as well. For a fixed-n design, set it to Inf
#' @param B Number of bootstrap samples (should be dividable by getDoParWorkers())
#' @param design "fixed.n" or "sequential". If design=="fixed.n", \code{n.min} and \code{boundary} are irrelevant, and all samples are drawn with n=n.max.
#' @param d The assumed true effect size. This can be a single number (this leads to a fixed assumed effect size, as in a classical power analysis) or a vector of numbers (e.g., \code{rnorm(100000, 0.5, 0.1)}). If it is a vector, the sampler draws a new effect size at each step. Hence, the provided distribution represents the uncertainty about the true effect size.
#' @param stepsize The number with which participants are added to the sample. If NA, the sample is increased +1 until it's 100, and +10 from that size on.
#' @param verbose Show output about progress?
#' @param cores number of parallel processes. If cores==1, no parallel framework is used.
#' @param alternative One of c("directional", "undirected") for directed (one-sided) or undirected (two-sided) hypothesis tests.
#' @param ETA Compute an estimate of the full simulation time? This adds some overhead to the simulation, so turn off for actual simulations.
#' @param ... Further parameters passed to the BF function
#'
#' @examples
#' \dontrun{
#' sim <- BPA.sim.ttest(d=0.5, n.min=20, n.max=300, boundary=Inf, 
#'				stepsize=1, design="sequential", B=1000, verbose=TRUE, cores=2)
#' save(sim, file="sim0.5.RData")
#' BPA.analysis(sim)
#' BPA.analysis(sim, boundary=6)
#' plot(sim, boundary=6)
#' plot(sim, boundary=6, n.max=80)
#'}
BPA.sim.ttest <- function(d, n.min=10, n.max=500, design="sequential", boundary=Inf, B=1000, stepsize=NA, alternative=c("directional", "undirected"), verbose=TRUE, cores=1, ETA=TRUE, ...) {

	alternative <- match.arg(alternative, c("directional", "undirected"))
	if (d==0 & alternative == "directional") {
		warning("Effect size is zero --> alternative is set to undirected!")
		alternative <- "undirected"
	}
	design <- match.arg(design, c("sequential", "fixed.n"))
	
	# Estimate the expected time for simulation
	if (ETA == TRUE) {
		print("Estimating duration of full simulation:")
		start.ETA <- Sys.time()
		BPA.sim.ttest(d=d, n.min=n.min, n.max=n.max, design=design, boundary=boundary, stepsize=stepsize, verbose=FALSE, cores=cores, alternative=alternative, ETA=FALSE, B=max(5, cores), ...)
		end.ETA <- Sys.time()
		print("A rough estimate of the necessary simulation time: ")
		# 1.3 is an empirically derived correction factor, probably due to parallel overhead?
		print((end.ETA-start.ETA)*(B/5)*1.3)
		print("Now starting main simulation.")
	}

	# register CPU cores for parallel processing
	registerDoParallel(cores=cores)	

	# define sample sizes that are simulated
	if (design=="fixed.n") {
		ns <- n.max
	} else if (design == "sequential"){
		if (is.na(stepsize)) {
			ns <- seq(n.min, min(n.max, 100), by=1)
			if (n.max > 100) ns <- c(ns, seq(105, min(n.max, 500), by=5))
			if (n.max > 500) ns <- c(ns, seq(510, n.max, by=10))
		} else {
			ns <- seq(n.min, n.max, by=stepsize)
		}
	}


	## ======================================================================
	## THE SIMULATION
	## ======================================================================

	start <- Sys.time()
	if (verbose==TRUE) print(paste0("Simulation started at ", start))
	flush.console()
	
	sim <- foreach(batch=1:getDoParWorkers(), .combine=rbind) %dopar% {

		max_b <- round(B/getDoParWorkers())
		res.counter <- 1

		# res saves the statistics at each step
		res <- matrix(NA, nrow=length(ns)*max_b, ncol=8, dimnames=list(NULL, c("id", "d", "boundary", "n", "logBF", "d.emp", "t.value", "p.value")))

		# fixed true ES? Simulate only one population
		if (length(d) == 1) {
			pop <- get_population(1000000, d)
			d1 <- d
		}

		for (b in 1:max_b) {
			# random true ES? Draw a new population at each step
			if (length(d) > 1) {
				d1 <- d[sample.int(length(d), 1)]
				pop <- get_population(1000000, d1)
			}

			if (verbose==TRUE)
				print(paste0(Sys.time(), ": batch = ", batch, "; d = ", round(d1, 2), "; Rep = ", b, "/", round(B/getDoParWorkers())))

			maxsamp <- pop[sample(nrow(pop), n.max), ]

			# res0 keeps the accumulating sample variables from this specific run
			res0 <- matrix(NA, nrow=length(ns), ncol=ncol(res), dimnames=dimnames(res))

			# increase sample size up to n.max (or only use n.max if design=="fixed.n")
			for (n in ns) {
				samp <- maxsamp[1:n, ]
				N <- nrow(samp)
				t0 <- t.test(samp[, 1], samp[, 2], var.equal=TRUE)

				if (alternative=="directional") {
					nullInterval <- c(0, Inf)
				} else {
					nullInterval <- NULL
				}
				
				# suppress the "t is large; approximation invoked" message
				suppressMessages({
					logBF0 <- BayesFactor::ttest.tstat(t0$statistic, N, N, nullInterval=nullInterval, ...)$bf
				})
				res0[which(ns == n), ] <- c(
					id	= batch*10^(floor(log(max_b, base=10))+2) + b,		# id is a unique id for each trajectory
					d	= d1,
					boundary = boundary,
					n	= N,
					logBF	= logBF0,
					d.emp	= 2*t0$statistic / sqrt(2*N-2),
					t.value	= t0$statistic,
					p.value	= t0$p.value)

				# if boundary is hit: stop sampling in this trajectory
				if (abs(logBF0) >= log(boundary)) {break;}
			} 	# of n

			res[res.counter:(res.counter+nrow(res0)-1), ] <- res0
			res.counter <- res.counter + nrow(res0)
		}  # of b's

		return(res)
	} # of %dopar%

	# reduce columns
	sim <- data.frame(sim)
	sim <- sim %>% filter(!is.na(d))

	if (verbose==TRUE) {
		end <- Sys.time()
		print(paste0("Simulation finished at ", end))
		cat("Duration: "); print(end - start)

	}
	
	res <- list(
		settings=list(
			n.min	= n.min, 
			n.max	= n.max, 
			design	= design,
			boundary= boundary,
			alternative = alternative,
			extra = list(...)
		),
		sim=sim
	)
	class(res) <- "BPA"
	return(res)
}



#' @export
#' @method print BPA
print.BPA <- function(x, ...) {

d.VAR <- var(x$sim$d)
if (d.VAR > 0) {
	QU <- paste0("; 25%/75% quantile = [", round(quantile(x$sim$d, prob=.25), 2), "; ", round(quantile(x$sim$d, prob=.75), 2), "]")
} else {QU <- ""}

cat(paste0("
Bayesian Power Analysis
--------------------------------	
Number of simulations: ", length(unique(x$sim$id)), "
Stopping boundary (evidential threshold): ", x$settings$boundary, "
Minimum n: ", x$settings$n.min, "
Maximum n: ", x$settings$n.max, "
Design: ", x$settings$design, "
Simulated true effect size (", ifelse(d.VAR > 0, "random", "fixed"), "): d=", round(mean(x$sim$d), 1), QU
))
}




# sim <- BPA.sim.ttest(d=0.5, n.min=20, n.max=300, boundary=Inf, stepsize=1, design="sequential", B=1000, verbose=TRUE, cores=2)
# save(sim, file="sim.0.5.RData")
#
# sim <- BPA.ttest.2(d=rnorm(100000, 0.5, 0.1), n.min=20, n.max=1000, boundary=Inf, stepsize=NA, design="sequential", B=10000, verbose=TRUE, cores=2)
# save(sim, file="sim.prior.0.5.RData")
#
# sim <- BPA.ttest.2(d=0, n.min=20, n.max=1000, boundary=Inf, stepsize=NA, design="sequential", B=10000, verbose=TRUE, cores=2)
# save(sim, file="sim.0.RData")

# sim <- BPA.ttest.2(d=0.5, n.min=20, n.max=400, boundary=Inf, stepsize=2, design="sequential", B=5000, verbose=TRUE, cores=2)
# save(sim, file="sim.0.5.step2RData")


# ---------------------------------------------------------------------
#

#sim <- BPA.sim.ttest(d=0, n.min=20, n.max=500, boundary=Inf, stepsize=NA, design="sequential", B=1000, verbose=TRUE, cores=2)
#save(sim, file="../../finalSims/sim.0b.RData")




# ---------------------------------------------------------------------
# Some plausibility checks

# At each cross-sectional cut, the success rate (via p values) should equal the calculated power
# This is only valid in the fixed-ES case!
# library(pwr)
#
# power.comparison <- data.frame()
# for (n in seq(min(sim$n), max(sim$n), by=10)) {
# 	print(n)
# 	power.analytically <- pwr.t.test(d=sim$d[1], n=n)$power
# 	power.sim <- sum(sim$p.value[sim$n==n]<.05)/length(sim$p.value[sim$n==n])
# 	power.comparison <- rbind(power.comparison, data.frame(power.analytically, power.sim))
# }
#
# plot(power.comparison, pch=21)
# abline(a=0, b=1, lty="dotted")
