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
#' @param rs
#' @param cores number of parallel processes. If cores==1, no parallel framework is used.
#'
#' @examples
#' \dontrun{
#' sim <- BPA.sim.ttest.2(d=0.5, n.min=20, n.max=300, boundary=Inf, stepsize=1, design="sequential", B=1000, verbose=TRUE, cores=2)
#' save(sim, file="sim0.5.RData")
#' BPA.analysis(sim)
#' BPA.analysis(sim, boundary=6)
#' plot(sim, boundary=6)
#' plot(sim, boundary=6, n.max=80)
#'}
BPA.sim.ttest.2 <- function(d, n.min=10, n.max=10000, design="sequential", boundary=10, B=1000, rs=1, stepsize=NA, verbose=TRUE, cores=1) {

	# register CPU cores for parallel processing
	registerDoParallel(cores=cores)

	design <- match.arg(design, c("sequential", "fixed.n"))

	# define sample sizes that are simulated
	if (design=="fixed.n") {
		ns <- n.max
	} else if (design == "sequential"){
		if (is.na(stepsize)) {
			# TODO: if n.max > 100 --> ns <- c(ns, seq())
			ns <- c(seq(n.min, 95, by=5), seq(100, n.max, by=10))
		} else {
			ns <- seq(n.min, n.max, by=stepsize)
		}
	}


	## ======================================================================
	## THE SIMULATION
	## ======================================================================

	start <- Sys.time()
	if (verbose==TRUE) print(paste0("Simulation started at ", start))

	sim <- foreach(batch=1:getDoParWorkers(), .combine=rbind) %dopar% {

		max_b <- round(B/getDoParWorkers())
		res.counter <- 1

		# res saves the statistics at each step
		res <- matrix(NA, nrow=length(ns)*length(rs)*max_b, ncol=10, dimnames=list(NULL, c("d", "r", "boundary", "batch", "rep", "n", "logBF", "d.emp", "t.value", "p.value")))

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

			for (r in rs) {

				# res0 keeps the accumulating sample variables from this specific run
				res0 <- matrix(NA, nrow=length(ns), ncol=ncol(res), dimnames=dimnames(res))

				# increase sample size up to n.max (or only use n.max if design=="fixed.n")
				for (n in ns) {
					samp <- maxsamp[1:n, ]
					N <- nrow(samp)
					t0 <- t.test(samp[, 1], samp[, 2], var.equal=TRUE)

					logBF0 <- BayesFactor::ttest.tstat(t0$statistic, N, N, rscale=r)$bf
					res0[which(ns == n), ] <- c(
						d	= d1,
						r 	= r,
						boundary = boundary,
						batch = batch,
						rep	= b,
						n	= N,
						logBF	= logBF0,
						d.emp	= 2*t0$statistic / sqrt(2*N-2),
						t.value	= t0$statistic,
						p.value	= t0$p.value)

					# if boundary is hit: stop sampling in this trajectory
					if (abs(logBF0) >= log(boundary)) {break;}
				} 	# of n
			} # of r's

			res[res.counter:(res.counter+nrow(res0)-1), ] <- res0
			res.counter <- res.counter + nrow(res0)
		}  # of b's

		return(res)
	} # of %dopar%

	# reduce columns
	sim <- data.frame(sim)
	sim <- sim %>% filter(!is.na(d)) %>% mutate(id=factor(batch):factor(rep)) %>% select(-batch, -rep)

	if (verbose==TRUE) {
		end <- Sys.time()
		print(paste0("Simulation finished at ", end))
		print(paste0("Duration: ", end - start))

	}
	
	res <- list(
		settings=list(
			n.min	= n.min, 
			n.max	= n.max, 
			design	= design,
			boundary= boundary
		),
		sim=sim
	)
	class(res) <- "BPA"
	return(res)
}



#' @export
#' @method print BPA
print.BPA <- function(x, ...) {
	cat(paste0("
Bayesian Power Analysis
--------------------------------	
Number of simulations: ", length(unique(x$sim$id)), "
Stopping boundary (evidential threshold): ", x$settings$boundary, "
Minimum n: ", x$settings$n.min, "
Maximum n: ", x$settings$n.max, "
Design: ", x$settings$design, "
"	
	))
}




# sim <- BPA.sim.ttest.2(d=0.5, n.min=20, n.max=300, boundary=Inf, stepsize=1, design="sequential", B=1000, verbose=TRUE, cores=2)
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
