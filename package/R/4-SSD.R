#' SSD = sample size determination
#' @details
#' Finds the optimal fixed sample size to achieve desired level of power. (Note: this function does *not* assume optional stopping, but checks all simulated fixed sample sizes).

#' @export
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot

#' @param BFDA A BFDA simulation object, resulting from the \code{BFDA.sim} function
#' @param boundary The critical boundary that defines the "success" of a study. This can be a single value, for example 1. In this case, all studies with BF < 1 are counted as "negative outcomes", and BF > 1 as positive outcomes (i.e., it quantifies the probability of getting a BF with the correct direction).
#' @param power The desired rate of true-positive results (in case of reality=H1)
#' @param alpha The desired rate of false-positive results (in case of reality=H0)
#' @param plot Display plot? (TRUE/FALSE)
SSD <- function(BFDA, boundary=c(1/3, 3), power=0.90, alpha=.025, plot=TRUE) {

	# just in case: order boundary
	boundary <- sort(boundary)
	logBoundary <- log(boundary)
	
	# determine the mode: H1 or H0-mode?
	analysis.mode <- ifelse(all(BFDA$sim$true.ES == 0), "H0", "H1")
	
	# sanity check: only valid, if same numbers of replications at each n.
	var.ns <- BFDA$sim %>% group_by(n) %>% summarize(ns=n()) %>% summarise(var(.$ns))
	if (var.ns>0) {stop("Not the same number of replications at each n - please set `boundary=Inf` in your BFDA.sim function.")}
	
	# retrieve n.max from simulation object	
	n.max <- max(BFDA$sim$n)
	
	# categorize trajectories according to critical boundary
	if (length(logBoundary) == 1 && logBoundary == 1) {
		categories <- BFDA$sim %>% group_by(n) %>% summarize(
			BF.positive = sum(logBF > 1)/n(),
			BF.negative = sum(logBF < 1)/n()
		)
	} else if (length(logBoundary) == 1) {
		categories <- BFDA$sim %>% group_by(n) %>% summarize(
			positive.results = sum(logBF > logBoundary)/n(),
			negative.results = sum(logBF <= logBoundary)/n()
		)
	} else if (length(logBoundary) == 2) {
		categories <- BFDA$sim %>% group_by(n) %>% summarize(
			positive.results = sum(logBF > logBoundary[2])/n(),
			inconclusive.results = sum(logBF <= logBoundary[2] & logBF > logBoundary[1])/n(),
			negative.results = sum(logBF <= logBoundary[1])/n()
		)
	}
	
	# find critical n where power [i.e., prob(positive result | assumed effect)] is desired power
	
	if (analysis.mode == "H1") {
		n.crit <- categories$n[which(categories$positive.results >= power)[1]]
	} else if (analysis.mode == "H0") {
		n.crit <- categories$n[which(categories$positive.results <= alpha)[1]]
	}	
	
	
	if (is.na(n.crit)) {  # desired power has NOT been achieved
		positive.results <- categories$positive.results[nrow(categories)]
		inconclusive.results <- categories$inconclusive.results[nrow(categories)]
		negative.results <- categories$negative.results[nrow(categories)]
	} else {  # desired power has been achieved at a specific n
		positive.results <- categories$positive.results[which(categories$n == n.crit)[1]]
		inconclusive.results <- categories$inconclusive.results[which(categories$n == n.crit)[1]]
		negative.results <- categories$negative.results[which(categories$n == n.crit)[1]]
	}
	
	# ---------------------------------------------------------------------
	# The verbal output
	
	cat("Sample size determination for a fixed-n design:\n---------------------------------------\n\n")
	
	# If delta==0, show output for H0 studies, else: output for H1 studies
		if (analysis.mode == "H0") {
			# output for null effect
			
			if (is.na(n.crit)) {
				res <- paste0("Even at the maximal simulated n = ", n.max,", the achieved false positive error rate (", round(positive.results*100, 1), "%) did not reach the desired level of ", round(alpha*100, 1), "%\n", "\nThe maximal n has long-term rates of:\n")
			} else {
				res <- paste0("A >= ", round(alpha*100, 1), "% (actual: ", round(positive.results*100, 1), "%) long-term rate of Type-I errors is achieved at n = ", n.crit, "\nThis setting implies long-term rates of:\n")
			}
						
			res <- paste0(res, ifelse(length(logBoundary) == 2, 
				paste0(round(inconclusive.results*100, 1), "% inconclusive results and\n"), ""),
				"   ", round(negative.results*100, 1), "% true-negative results."
			)
			
			cat(res)
			
			
		} else if (analysis.mode == "H1") {
			# output for true effect
			if (is.na(n.crit)) {
				res <- paste0("Even at the maximal simulated n = ", n.max,", the achieved power (", round(positive.results*100, 1), "%) did not reach the desired level of ", round(power*100, 1), "%\n", "\nThe maximal n has long-term rates of:\n")
			} else {
				res <- paste0("A >= ", round(power*100, 1), "% (actual: ", round(positive.results*100, 1), "%) power achieved at n = ", n.crit, "\nThis setting implies long-term rates of:\n")
			}
						
			res <- paste0(res, ifelse(length(logBoundary) == 2, 
				paste0(round(inconclusive.results*100, 1), "% inconclusive results and\n"), ""),
				"   ", round(negative.results*100, 1), "% false-negative results."
			)
			
			cat(res)
		}	
			
		
	# ---------------------------------------------------------------------
	# The plot
	
	if (plot==TRUE) {
		if (length(logBoundary) == 1) {
			colors <- scales::alpha(c("green4", "red3"), alpha=0.5)
			color.labels <- c(
				paste0("Positive outcome with BF > ", round(boundary, 2)),
				paste0("Negative outcome with BF < ", round(boundary, 2))
			)
		}
		if (length(logBoundary) == 2) {
			colors <- scales::alpha(c("green4", "orange1", "red3"), alpha=0.5)
			color.labels <- c(
				paste0("Positive with BF > ", round(boundary[2], 2)),
				paste0("Inconclusive with ", round(boundary[1], 2), " < BF < ", round(boundary[2], 2)),
				paste0("Negative with BF < ", round(boundary[1], 2))
			)
		}		
	
		# crop right side of plot where everything is green/ or red
		if (analysis.mode == "H1") {
			n.display.max <- which(categories$positive.results==1)[1]
			if (!is.na(n.display.max)) categories <- categories[1:n.display.max, ]
		}
		if (analysis.mode == "H0") {
			n.display.max <- which(categories$negative.results==1)[1]
			if (!is.na(n.display.max)) categories <- categories[1:n.display.max, ]
		}
		if (is.na(n.crit)) n.display.max <- n.max
	
		# Compute some variables to allow stacked rects
		if (length(logBoundary) == 1) {
			categories$p_i <- categories$positive.results
		} else {
			categories$p_i <- categories$positive.results + categories$inconclusive.results
		}
		categories$n_next <- lead(categories$n)
		cat2 <- reshape2::melt(categories, id.vars="n")	
		cat2$n_next <- lead(cat2$n)
	
		p1 <- ggplot(categories, aes()) + geom_rect(aes(xmin=n, xmax=n_next, ymin=0, ymax=positive.results), fill=colors[1])
		if (length(logBoundary) == 2) {
			p1 <- p1 + geom_rect(aes(xmin=n, xmax=n_next, ymin=positive.results, ymax=p_i), fill=colors[2], color=NA)
		}
		p1 <- p1 + geom_rect(aes(xmin=n, xmax=n_next, ymin=p_i, ymax=1), fill=colors[length(colors)])
	
	
		p1 <- p1 + theme_bw() + scale_y_continuous(labels=scales::percent) + xlab("n") + ylab("% of studies") + scale_fill_manual("Study outcome", values=colors, labels=color.labels)
		p1 <- p1 + annotate("segment", x=n.crit, xend=n.crit, y=0, yend=positive.results)
		p1 <- p1 + annotate("segment", x=min(cat2$n), xend=n.crit, y=positive.results, yend=positive.results, linetype="dashed")
		p1 <- p1 + annotate("point", x=n.crit, y=positive.results, size=3)
	
		if (analysis.mode == "H1") {
			p1 <- p1 + annotate("text", x=n.crit, y=positive.results, label=paste0(round(power*100, 1), "% power achieved at n = ", n.crit), angle=90, size=3, hjust=1.1, vjust=1.5, fontface=2)
		} else {
			p1 <- p1 + annotate("text", x=n.crit, y=positive.results, label=paste0(round(alpha*100, 1), "% Type-I error rate achieved at n = ", n.crit), size=3, hjust=-0.1, vjust=-0.2, fontface=2)
		}
	
		return(p1)
	# of if plot==TRUE 
	} else {
		return(list(
			n.crit = n.crit,
			positive.results = positive.results,			
			inconclusive.results = inconclusive.results,			
			negative.results = negative.results
		))
	}
}

