#' Simulate ... TODO
#' @export
#' @import BayesFactor
#' @import dplyr
#' @import doParallel
#' @importFrom doRNG %dorng%
#'
#' @param n.min Minimum n before optional stopping is started
#' @param n.max Maximum n - if that is reached without hitting a boundary, the run is aborted
#' @param boundary The Bayes factor where a sequential run is stopped. For a fixed-n design, \code{boundary} is automatically set to Inf. Provide either two values for lower and upper boundary (e.g., c(1/3, 6). If only one value is provided, it automatically uses its reciprocal for the other boundary.
#' @param B Number of bootstrap samples; should be dividable by the numbers of \code{cores} (see also \code{getDoParWorkers()})
#' @param design "fixed.n" or "sequential". If design=="fixed.n", \code{n.min} and \code{boundary} are irrelevant, and all samples are drawn with n=n.max.
#' @param expected.ES The assumed true effect size. This can be a single number (this leads to a fixed assumed effect size, as in a classical power analysis) or a vector of numbers (e.g., \code{rnorm(100000, 0.5, 0.1)}). If it is a vector, the sampler draws a new effect size from this vector at each step. In this case, the provided distribution represents the uncertainty about the true effect size.
#' @param type Currently three designs are implemented: c("t.between", "t.paired", "correlation")
#' @param prior Define the prior distribution of your alternative hypothesis. Argument takes a list as an input: The first element of the list should be a character string defining the type of the distribution ("t", "Cauchy", or "normal"). The second element of the list is another list containing the parameters of the distribution. For the t-distribution, this is prior.location, prior.df, and prior.scale; for the Cauchy-distribution, it is prior.location and prior.scale; for the normal distribution, it is prior.mean and prior.variance. The default for t-tests is \code{list("Cauchy", list(prior.location = 0, prior.scale = sqrt(2)/2))}
#' @param stepsize The number with which participants are added to the sample. If NA, the sample is increased +1 until it's 100, and +10 from that size on.
#' @param verbose Show output about progress?
#' @param cores number of parallel processes. If cores==1, no parallel framework is used.
#' @param alternative One of c("directional", "undirected") for directed (one-sided) or undirected (two-sided) hypothesis tests in data analysis. Hence, this refers to the directionality of the analysis prior
#' @param ETA Compute an estimate of the full simulation time? This adds some overhead to the simulation, so turn off for actual simulations. NOT IMPLEMENTED YET
#' @param options.sample Further parameters passed to the data generating function (depending on the \code{type} of design). NOT IMPLEMENTED YET
#' @param seed The seed that is passed to the \code{dorng} function (which ensures reproducibility with parallel processing). If this parameter is set to \code{NULL}, a new seed is chosen at each run.
#' @param ... Further parameters passed to the BF.test function. Most importantly, the scale parameter \code{rscale} can be passed to adjust the width of the Cauchy analysis prior.
#'
#' @examples
#' \dontrun{
#' sim <- BFDA.sim(expected.ES=0.5, type="t.between", prior=list("Cauchy", list(prior.location=0, prior.scale=1)),
#'                 n.min=20, n.max=300, boundary=Inf, stepsize=1, design="sequential", B=1000, verbose=TRUE, cores=2)
#' save(sim, file="sim0.5.RData")
#' BFDA.analyze(sim)
#' BFDA.analyze(sim, boundary=6)
#' plot(sim, boundary=6)
#' plot(sim, boundary=6, n.max=80)
#'}


BFDA.sim <- function(expected.ES, type=c("t.between", "t.paired", "correlation"), prior = NULL, n.min=10, n.max=500, design=c("sequential", "fixed.n"), boundary=Inf, B=1000, stepsize=NA, alternative=c("two.sided", "greater", "less"), verbose=TRUE, cores=1, ETA=FALSE, options.sample=list(), seed=1234, ...) {
	
	# link to test specific functions
	# get() can reference a function by its (string) name
	sample.function <- get(paste0("sample.", type))
	select.function <- get(paste0("select.", type))
	BF.test.function <- get(paste0("BF.test.", type))
	freq.test.function <- get(paste0("freq.test.", type))
	prior.check.function <- get(paste0("prior.check.", type))

	alternative <- match.arg(alternative, c("two.sided", "greater", "less"))

	design <- match.arg(design, c("sequential", "fixed.n"))
	type <- match.arg(type, c("t.between", "t.paired", "correlation"))
	prior <- prior.check.function(prior)
	
	# # Estimate the expected time for simulation
	# if (ETA == TRUE) {
	# 	print("Estimating duration of full simulation:")
	# 	start.ETA <- Sys.time()
	# 	testRuns <- max(5, cores)
	# 	BFDA.sim(expected.ES=expected.ES, type=type, n.min=n.min, n.max=n.max, design=design, boundary=boundary, stepsize=stepsize, verbose=FALSE, cores=cores, alternative=alternative, ETA=FALSE, B=testRuns, ...)
	# 	end.ETA <- Sys.time()
	# 	print("A rough estimate of the necessary simulation time: ")
	# 	# 1.3 is an empirically derived correction factor, probably due to parallel overhead?
	# 	print((end.ETA-start.ETA)*(B/testRuns)*1.3)
	# 	print("Now starting main simulation.")
	# }

	# register CPU cores for parallel processing
	registerDoParallel(cores=cores)	

	# define sample sizes that are simulated
	if (design=="fixed.n") {
		ns <- n.max
		if (!is.infinite(boundary)) {
			warning("For fixed-n designs, boundary is automatically set to Inf")
			boundary <- Inf
		}
	} else if (design == "sequential"){
		if (is.na(stepsize)) {
			ns <- seq(n.min, min(n.max, 100), by=1)
			if (n.max > 100) ns <- c(ns, seq(105, min(n.max, 500), by=5))
			if (n.max > 500) ns <- c(ns, seq(510, n.max, by=10))
		} else {
			ns <- seq(n.min, n.max, by=stepsize)
		}
	}
	
	if (length(boundary) == 1) boundary <- c(1/boundary, boundary)
	boundary <- sort(boundary)
	logBoundary <- log(boundary)


	## ======================================================================
	## THE SIMULATION
	## ======================================================================

	start <- Sys.time()
	if (verbose==TRUE) print(paste0("Simulation started at ", start))
	flush.console()
	
	sim <- foreach(batch=1:getDoParWorkers(), .combine=rbind, .options.RNG=seed) %dorng% {

		max_b <- round(B/getDoParWorkers())
		res.counter <- 1

		# res saves the statistics at each step
		res <- matrix(NA, nrow=length(ns)*max_b, ncol=7, dimnames=list(NULL, c("id", "true.ES", "n", "logBF", "emp.ES", "statistic", "p.value")))

		# run max_b iterations in each parallel worker
		for (b in 1:max_b) {
			# Draw a new maximum sample at each step
			# If expected.ES has more than 1 value: draw a random value
			expected.ES.1 <- expected.ES[sample.int(length(expected.ES), 1)]
			maxsamp <- sample.function(n.max, expected.ES.1, options.sample)

			if (verbose==TRUE)
				print(paste0(Sys.time(), ": batch = ", batch, "; true.ES = ", round(expected.ES.1, 2), "; Rep = ", b, "/", round(B/getDoParWorkers())))			

			# res0 keeps the accumulating sample variables from this specific run
			res0 <- matrix(NA, nrow=length(ns), ncol=ncol(res), dimnames=dimnames(res))

			# increase sample size up to n.max (or only use n.max if design=="fixed.n")
			for (n in ns) {
				samp <- select.function(maxsamp, n)
				
				# do the frequentist test
				freq.test <- freq.test.function(samp, alternative)

				# do the BF test; supply freq.test to access t.value for faster computation
				logBF <- BF.test.function(samp, alternative, freq.test, prior, ...)
					
				res0[which(ns == n), ] <- c(
					id		= batch*10^(floor(log(max_b, base=10))+2) + b,		# id is a unique id for each trajectory
					true.ES	= expected.ES.1,
					n		= n,
					logBF	= logBF,
					emp.ES	= freq.test$emp.ES,
					statistic = freq.test$statistic,
					p.value	= freq.test$p.value)

				# if boundary is hit: stop sampling in this trajectory
				if (logBF<=logBoundary[1] | logBF >= logBoundary[2]) {break;}
			} 	# of n

			res[res.counter:(res.counter+nrow(res0)-1), ] <- res0
			res.counter <- res.counter + nrow(res0)
		}  # of b's

		return(res)		
	} # of %dorng%

	# reduce columns
	sim <- data.frame(sim)
	sim <- sim %>% filter(!is.na(true.ES))

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
			prior = prior,
			boundary= boundary,
			alternative = alternative,
			type=type,
			extra = list(...),
			packageVersion = packageVersion("BFDA")
		),
		sim=sim
	)
	class(res) <- "BFDA"
	return(res)
}



#' @export
#' @method print BFDA
print.BFDA <- function(x, ...) {

	versionCheck(x)

	ES.VAR <- var(x$sim$true.ES)
	if (nrow(x$sim) > 1 && ES.VAR > 0) {
		QU <- paste0("; 25%/75% quantile = [", round(quantile(x$sim$true.ES, prob=.25), 2), "; ", round(quantile(x$sim$true.ES, prob=.75), 2), "]")
	} else {QU <- ""}
	
	ES.type <- switch(x$settings$type,
		"t.between" = "Cohen's d",
		"t.paired" = "Cohen's d",
		"correlation" = "correlation",
		{paste0("ERROR: Test type ", x$settings$type, " not recognized.")}	#default
	)
	

	cat(paste0("
Bayes Factor Design Analysis
--------------------------------	
Number of simulations: ", length(unique(x$sim$id)), "
Stopping boundaries (evidential thresholds): ", round(x$settings$boundary[1], 2), "; ", round(x$settings$boundary[2], 2), "
Directionality of H1 analysis prior: ", x$settings$alternative, "
Minimum n: ", x$settings$n.min, "
Maximum n: ", x$settings$n.max, "
Design: ", x$settings$design, "
Simulated true effect size (", ifelse(nrow(x$sim) > 1 && ES.VAR > 0, "random", "fixed"), "): ", ES.type, " = ", round(mean(x$sim$true.ES), 3), QU,
"\n"
	))
}

