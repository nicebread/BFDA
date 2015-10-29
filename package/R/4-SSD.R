#' SSD = sample size determination
#' @details
#' Finds the optimal sample size to achieve desired level of power.

#' @export
#' @importFrom reshape2 melt
# TODO: Add link to proper function name
#' @param BPA A BPA simulation object, resulting from the \code{BPA.sim} function
#' @param criterion The critical boundary that defines the "success" of a study. This can be a single value, for example 1. In this case, all studies with BF < 1 are counted as "negative outcomes", and BF > 1 as positive outcomes (i.e., it quantifies the probability of getting a BF with the correct direction).
#' @param power The desired rate of true-positive results (in case of reality=H1)
#' @param alpha The desired rate of false-positive results (in case of reality=H0)
SSD <- function(BPA, criterion=c(1/3, 3), power=0.90, alpha=.025) {
	
	# just in case: order criterion
	criterion <- sort(criterion)
	logCriterion <- log(criterion)
	
	# determine the mode: H1 or H0-mode?
	analysis.mode <- ifelse(all(BPA$sim$d == 0), "H0", "H1")
	
	# sanity check: only valid, if same numbers of replications at each n.
	var.ns <- BPA$sim %>% group_by(n) %>% summarize(ns=n()) %>% summarise(var(.$ns))
	if (var.ns>0) {stop("Not numbers of replications at each n - consider setting boundary=Inf in your BPA.sim function.")}
	
	# categorize trajectories according to critical boundary
	if (length(logCriterion) == 1 && logCriterion == 1) {
		categories <- BPA$sim %>% group_by(n) %>% summarize(
			BF.positive = sum(logBF > 1)/n(),
			BF.negative = sum(logBF < 1)/n()
		)
	}
	if (length(logCriterion) == 1) {
		categories <- BPA$sim %>% group_by(n) %>% summarize(
			positive.results = sum(logBF > logCriterion)/n(),
			negative.results = sum(logBF <= logCriterion)/n()
		)
	}
	if (length(logCriterion) == 2) {
		categories <- BPA$sim %>% group_by(n) %>% summarize(
			positive.results = sum(logBF > logCriterion[2])/n(),
			inconclusive.results = sum(logBF <= logCriterion[2] & logBF > logCriterion[1])/n(),
			negative.results = sum(logBF <= logCriterion[1])/n()
		)
	}
	
	# find critical n where power [i.e., prob(positive result | assumed effect)] is desired power
	
	if (analysis.mode == "H1") {
		n.crit <- categories$n[which(categories$positive.results >= power)[1]]
	} else {
		n.crit <- categories$n[which(categories$positive.results <= alpha)[1]]
	}
	positive.results <- categories$positive.results[which(categories$n == n.crit)[1]]
	inconclusive.results <- categories$inconclusive.results[which(categories$n == n.crit)[1]]
	negative.results <- categories$negative.results[which(categories$n == n.crit)[1]]
	
	# ---------------------------------------------------------------------
	# The plot
	
	if (length(logCriterion) == 1) {
		colors <- c("green4", "red3")
		color.labels <- c(
			paste0("Posiive outcome with BF > ", round(criterion, 2)),
			paste0("Negative outcome with BF < ", round(criterion, 2))
		)
	}
	if (length(logCriterion) == 2) {
		colors <- c("green4", "orange1", "red3")
		color.labels <- c(
			paste0("Positive with BF > ", round(criterion[2], 2)),
			paste0("Inconclusive with ", round(criterion[1], 2), " < BF < ", round(criterion[2], 2)),
			paste0("Negative with BF < ", round(criterion[1], 2))
		)
	}
	
	# If delta==0, show output for H0 studies, else: output for H1 studies
	if (analysis.mode == "H0") {
		# output for null effect
		cat(paste0("A ", round(positive.results*100, 1), "% long-term rate of Type-I errors is achieved at n = ", n.crit, "\n",
		"This setting implies long-term rates of:\n", 
		"   ", 
		ifelse(length(logCriterion) == 2, 
			paste0(round(inconclusive.results*100, 1), "% inconclusive results and\n"),
			""),
		"   ", round(negative.results*100, 1), "% true-negative results."
		))
	} else {
		# output for true effect
		cat(paste0(round(positive.results*100, 1), "% power achieved at n = ", n.crit, "\n",
		"This setting implies long-term rates of:\n", 
		ifelse(length(logCriterion) == 2, 
			paste0(round(inconclusive.results*100, 1), "% inconclusive results and\n"),
			""),
		"   ", round(negative.results*100, 1), "% false-negative results."
		))
	}	
	
	
	# ---------------------------------------------------------------------
	# The plot
	
	# Compute some variables to allow stacked rects
	categories$p_i <- categories$positive.results + categories$inconclusive.results
	categories$n_next <- lead(categories$n)
	cat2 <- reshape2::melt(categories, id.vars="n")	
	cat2$n_next <- lead(cat2$n)
	
	p1 <- ggplot(categories, aes()) + geom_rect(aes(xmin=n, xmax=n_next, ymin=0, ymax=positive.results), fill=colors[1])
	if (length(logCriterion) == 2) {
		p1 <- p1 + geom_rect(aes(xmin=n, xmax=n_next, ymin=positive.results, ymax=p_i), fill=colors[2])
	}
	p1 <- p1 + geom_rect(aes(xmin=n, xmax=n_next, ymin=p_i, ymax=1), fill=colors[length(colors)])
	
	
	p1 <- p1 + theme_bw() + scale_y_continuous(labels=scales::percent) + xlab("n") + ylab("% of studies") + scale_fill_manual("Study outcome", values=colors, labels=color.labels)
	p1 <- p1 + annotate("segment", x=n.crit, xend=n.crit, y=0, yend=positive.results)
	p1 <- p1 + annotate("segment", x=min(cat2$n), xend=n.crit, y=positive.results, yend=positive.results, linetype="dashed")
	p1 <- p1 + annotate("point", x=n.crit, y=positive.results, size=3)
	
	if (analysis.mode == "H1") {
		p1 <- p1 + annotate("text", x=n.crit, y=positive.results, label=paste0(round(power*100), "% power achieved at n = ", n.crit), angle=90, size=3, hjust=1.1, vjust=1.5, fontface=2)
	} else {
		p1 <- p1 + annotate("text", x=n.crit, y=positive.results, label=paste0(round(alpha*100, 1), "% Type-I error rate achieved at n = ", n.crit), size=3, hjust=-0.1, vjust=-0.2, fontface=2)
	}
	
	return(p1)
}

