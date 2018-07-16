# make sure that the id is really unique!
#' Analyze a BFDA.sim object
#'
#' @param BFDA The result object from a \code{\link{BFDA.sim}} function call
#' @param design The function can analyze the simulations with regard to a fixed-n (\code{design="fixed"}) or in a sequential style (\code{design="sequential"})
#' @param n.min What is the minimum n that is sampled before optional stopping is started? Defaults to the smallest n in the BFDA object (only relevant for \code{design="sequential"}).
#' @param n.max What is the maximum n that is sampled in sequential sampling? Defaults to the largest n in the BFDA object (only relevant for \code{design="sequential"}).
#' @param n The fixed-n sample size that should be evaluated (only relevant for \code{design="fixed"})
#' @param boundary At which BF boundary should sequential trajectories stop? Either a single number (then the reciprocal is taken as the other boundary), or a vector of two numbers for lower and upper boundary (only relevant for \code{design="sequential"}).
#' @param verbose Print information about analysis?
#' @param alpha For a frequentist analysis in the fixed-n case: Use this alpha level.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' BFDA.analyze(sim, design="sequential", boundary=6, n.min=20, n.max=80)
#' }
BFDA.analyze <- function(BFDA, design=c("sequential", "fixed"), boundary, n.min=NULL, n.max=NULL, n=NULL, verbose=TRUE, alpha=.05) {

	# some preliminary checks
	versionCheck(BFDA)	
	design <- match.arg(design, c("sequential", "fixed"))
	if (design=="fixed" & is.null(n)) stop("If design='fixed' you have to specify n (not n.min nor n.max).")
	if (design=="fixed" & all(is.null(boundary))) stop("If design='fixed' you have to specify a BF value that is used as a decision criterion (i.e., set parameter 'boundary' to some value).")
		
	if (length(boundary) == 1) boundary <- sort(c(boundary, 1/boundary))
		
	if (design=="sequential") {
		return(do.call("BFDA.analyze.sequential", as.list(environment())))
	}
	if (design=="fixed") {
		return(do.call("BFDA.analyze.fixed", as.list(environment())))
	}
}



BFDA.analyze.sequential <- function(BFDA, boundary, n=NULL, n.min=NULL, n.max=NULL, design="sequential", verbose=TRUE, alpha=.05) {
	
	# set variables
	sim <- BFDA$sim
	if (is.null(n.max) || is.na(n.max) || is.infinite(n.max)) n.max <- max(sim$n)
	if (is.null(n.min) || is.na(n.min)) n.min <- min(sim$n)
	if (all(is.null(boundary))) boundary <- BFDA$settings$boundary
		
	# reduce simulation to relevant data
	sim <- sim %>% filter(n >= n.min, n <= n.max)
	
	if (length(boundary) == 1) boundary <- c(1/boundary, boundary)
	boundary <- sort(boundary)
	logBoundary <- log(boundary)
		
	if (boundary[2] > max(BFDA$settings$boundary)) warning(paste0("Error: The selected boundary (", boundary[2], ") for analysis is larger than the largest stopping boundary (", max(BFDA$settings$boundary), ") in the simulation stage. Cannot produce a meaningful analysis."))
				
		
	if (n.max > max(sim$n)) warning(paste0("Error: The selected n.max (", n.max, ") for analysis is larger than the largest n (", max(sim$n), ") in the simulation stage. Cannot produce a meaningful analysis."))


	# For the densities: Data frames of stopping times / stopping BFs
	n.max.hit <- sim %>% 
		group_by(id) %>% 
		filter(n == n.max, max(logBF) <= logBoundary[2] & min(logBF) >= logBoundary[1]) %>% 
		mutate(hitCondition = "n.max")
	
	# reduce to *first* break of a boundary
	boundary.hit <- sim %>% group_by(id) %>%
		filter(logBF>=logBoundary[2] | logBF<=logBoundary[1]) %>%
		filter(row_number()==1) %>% ungroup()	%>%
		mutate(hitCondition = "boundary")

	endpoint <- bind_rows(n.max.hit, boundary.hit)
	
	# compute counts of three outcome possibilities
	all.traj.n <- length(unique(sim$id))
	boundary.traj.n <- length(unique(boundary.hit$id))
	n.max.traj.n <- length(unique(n.max.hit$id))

	boundary.upper.traj.n <- length(unique(boundary.hit$id[boundary.hit$logBF>0]))
	boundary.lower.traj.n <- length(unique(boundary.hit$id[boundary.hit$logBF<0]))
	
	# sanity checks: all outcomes should sum to the overall number of trajectories
	if (all.traj.n != boundary.traj.n + n.max.traj.n | all.traj.n != n.max.traj.n + boundary.upper.traj.n + boundary.lower.traj.n) warning("outcomes do not sum up to 100%!")
		
	n.max.hit.frac <- n.max.traj.n/all.traj.n
	boundary.hit.frac <- boundary.traj.n/all.traj.n
	upper.hit.frac <- boundary.upper.traj.n/all.traj.n
	lower.hit.frac <- boundary.lower.traj.n/all.traj.n
	
	# ---------------------------------------------------------------------
	#  compute densities

	ns.upper <- boundary.hit$n[boundary.hit$logBF>0]
	if (length(ns.upper) >= 2) {
		d.top <- density(ns.upper, from=min(sim$n), to=max(ns.upper))
	} else {d.top <- NULL}
	
	ns.lower <- boundary.hit$n[boundary.hit$logBF<0]
	if (length(ns.lower) >= 2) {
		d.bottom <- density(ns.lower, from=min(sim$n), to=max(ns.lower))
	} else {d.bottom <- NULL}
	
	logBF.right <- n.max.hit$logBF
	if (length(logBF.right) >= 2) {
		d.right <- density(logBF.right, from=min(logBF.right), to=max(logBF.right))
	} else {d.right <- NULL}
	
	
	if (var(sim$n) == 0) {
		p.value <- sum(sim$p.value < alpha)/all.traj.n*100
	} else {p.value <- NA}
	# ---------------------------------------------------------------------
	# Output

	res <- list(
		settings = BFDA$settings,
		n.max=n.max,
		d.top=d.top,
		d.bottom = d.bottom,
		d.right = d.right,
		n.max.hit.frac = n.max.hit.frac,
		boundary.hit.frac = boundary.hit.frac,
		upper.hit.frac = upper.hit.frac,
		lower.hit.frac = lower.hit.frac,
		logBF.right = logBF.right,
		upper.hit.ids = unique(boundary.hit$id[boundary.hit$logBF>0]),
		lower.hit.ids = unique(boundary.hit$id[boundary.hit$logBF<0]),
		n.max.hit.ids = unique(n.max.hit$id),
		all.traj.n = all.traj.n, 
		boundary.traj.n = boundary.traj.n, 
		n.max.traj.n = n.max.traj.n,
		n.max.hit.logBF = n.max.hit$logBF,
		endpoint.n = endpoint$n,
		endpoint = endpoint,
		alpha = alpha,
		p.value = p.value,
		ASN = ceiling(mean(endpoint$n)),
		n.max.hit.H1 = sum(n.max.hit$logBF > log(3))/all.traj.n,
		n.max.hit.inconclusive = sum(n.max.hit$logBF < log(3) & n.max.hit$logBF > log(1/3))/all.traj.n,
		n.max.hit.H0 = sum(n.max.hit$logBF < log(1/3))/all.traj.n
	)
	
	class(res) <- "BFDAanalysisSequential"
	return(res)
}






#' Print a BFDA analysis
#' @export
#' @method print BFDAanalysisSequential
#' @param x A BFDA-analysis object (which is return from \code{BFDA.analyze})
#' @param digits Number of digits in display
#' @param ... (not used)
print.BFDAanalysisSequential <- function(x, ..., digits=1) {
with(x, {
	
	print(data.frame(
		outcome = c(paste0("Studies terminating at n.max (n=", n.max, ")"), "Studies terminating at a boundary", "--> Terminating at H1 boundary", "--> Terminating at H0 boundary"),
		percentage = paste0(c(round(n.max.hit.frac*100, digits), round(boundary.hit.frac*100, digits), round(upper.hit.frac*100, digits), round(lower.hit.frac*100, digits)), "%")))
	
	# If some studies stopped at n.max: report of categories of resulting Bayes factor
	
	# TODO: "(BF > 3)": insert actual boundary from parameter!
	if (n.max.traj.n > 0) {
		cat(paste0("\nOf ", round(n.max.traj.n/all.traj.n*100, digits), "% of studies terminating at n.max (n=", n.max, "):\n",
			round(sum(n.max.hit.logBF > log(3))/all.traj.n*100, digits), "% showed evidence for H1 (BF > 3)\n", 
			round(sum(n.max.hit.logBF < log(3) & n.max.hit.logBF > log(1/3))/all.traj.n*100, digits), "% were inconclusive (3 > BF > 1/3)\n", 
			round(sum(n.max.hit.logBF < log(1/3))/all.traj.n*100, digits), "% showed evidence for H0 (BF < 1/3)\n"
		))
	}
	
	# If some studies stopped at boundary: report ASN across all trials
	if (boundary.traj.n > 0) {	
		cat(paste0("\nAverage sample number (ASN) at stopping point (both boundary hits and n.max): n = ", ASN))
		cat("\n\nSample number quantiles (50/80/90/95%) at stopping point:\n")
		print(ceiling(quantile(endpoint.n, prob=c(0.5, 0.80, 0.90, 0.95))))
	}
})	
}





BFDA.analyze.fixed <- function(BFDA, n=NULL, n.min=NA, n.max=NA, design="fixed", boundary=NA, verbose=TRUE, alpha=.05) {
	# set variables
	sim <- BFDA$sim
	n.eval <- n	# rename, otherwise dplyr might get confused
		
	# reduce simulation to relevant data
	n.max.hit <- sim %>% filter(n == n.eval)

	if (n.eval > max(n.max.hit$n)) error(paste0("Error: The selected n (", n.eval, ") for analysis is larger than the largest n (", max(n.max.hit$n), ") in the simulation stage. Cannot produce a meaningful analysis."))

	# ---------------------------------------------------------------------
	# Output

	res <- list(
		settings = BFDA$settings,
		boundary = boundary,
		n.max.hit.logBF = n.max.hit$logBF
	)
	
	class(res) <- "BFDAanalysisFixed"
	return(res)
}






#' Print a BFDA analysis
#' @export
#' @method print BFDAanalysisFixed
#' @param x A BFDA-analysis object (which is return from \code{BFDA.analyze})
#' @param digits Number of digits in display
#' @param ... (not used)
print.BFDAanalysisFixed <- function(x, ..., digits=1) {
	with(x, {
		cat(paste0(
			round(sum(n.max.hit.logBF > log(max(boundary)))*100/length(n.max.hit.logBF), digits), "% showed evidence for H1 (BF > ", round(max(boundary), 4), ")\n", 
			round(sum(n.max.hit.logBF < log(max(boundary)) & n.max.hit.logBF > log(min(boundary)))*100/length(n.max.hit.logBF), digits), "% were inconclusive (", round(min(boundary), 4), " < BF < ", round(max(boundary), 4), ")\n", 
			round(sum(n.max.hit.logBF < log(min(boundary)))*100/length(n.max.hit.logBF), digits), "% showed evidence for H0 (BF < ", round(min(boundary), 4), ")\n"
			))
			
})}
