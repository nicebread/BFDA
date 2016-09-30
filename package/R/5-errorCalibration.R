# TODO: This function is a stub! DO NOT USE IT!

errorCalibration <- function(BFDA.H1, BFDA.H0, alpha=.05) {
	sim <- BFDA$sim
	if (is.na(n.max)) n.max <- max(sim$n)
	if (is.na(n.min)) n.min <- min(sim$n)
	if (all(is.na(boundary))) boundary <- max(sim$boundary)
		
	# reduce simulation to relevant data
	sim <- sim %>% filter(n >= n.min, n <= n.max)
	
	if (length(boundary) == 1) boundary <- sort(c(boundary, 1/boundary))
	logBoundary <- log(boundary)
		
	if (boundary[2] > max(sim$boundary)) warning(paste0("Error: The selected boundary (", boundary[2], ") for analysis is larger than the smallest stopping boundary (", min(sim$boundary), ") in the simulation stage. Cannot produce a meaningful analysis."))
				
		
	if (n.max > max(sim$n)) warning(paste0("Error: The selected n.max (", n.max, ") for analysis is larger than the largest n (", max(sim$n), ") in the simulation stage. Cannot produce a meaningful analysis."))

	# For the densities: Data frames of stopping times / stopping BFs
	n.max.hit <- sim %>% group_by(id) %>% filter(n == n.max, max(logBF) <= logBoundary[2] & min(logBF) >= logBoundary[1])

	# reduce to *first* break of a boundary
	boundary.hit <- sim %>% group_by(id) %>%
		filter(logBF>=logBoundary[2] | logBF<=logBoundary[1]) %>%
		filter(row_number()==1) %>% ungroup()	

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
		alpha = alpha,
		p.value = p.value,
		ASN = ceiling(mean(endpoint$n)),
		n.max.hit.H1 = sum(n.max.hit$logBF > log(3))/all.traj.n,
		n.max.hit.inconclusive = sum(n.max.hit$logBF < log(3) & n.max.hit$logBF > log(1/3))/all.traj.n,
		n.max.hit.H0 = sum(n.max.hit$logBF < log(1/3))/all.traj.n
	)
	
	class(res) <- "BFDAanalysis"
	return(res)
}