# source("7-Leapfrog.R", echo=TRUE)

devtools::load_all("package")
library(truncdist)
library(doRNG)
library(tidyverse)


#' @param overall.N Maximum available N *per group*

leapfrog <- function(standardTreatment.ES, newTreatment.ES, overall.N, n.min, n.max, boundary, verbose=FALSE, ...) {
	
	N.used <- 0
	hops <- data.frame()
	
	if (length(boundary) == 1) boundary <- c(1/boundary, boundary)
	boundary <- sort(boundary)
	logBoundary <- log(boundary)

	armCounter <- 0
	for (treatment in newTreatment.ES) {
	
		N.left <- overall.N - N.used
		armCounter <- armCounter+1
		
		if (verbose==TRUE) {
			cat(paste0("Testing new treatment with true ES = ", round(treatment, 2), ". ", N.left, " participants left.\n"))
		}
	
		compare <- BFDA.sim(expected.ES = treatment - standardTreatment.ES, type="t.between", n.min=n.min, n.max=min(n.max, N.left), design="sequential", boundary=boundary, B=1, alternative="directional", verbose=FALSE, cores=1, ...)
	
		# look at final BF (either stopped at boundary or at n.max)
		outcome <- "ERROR"
		decision <- "ERROR"
		decisionQuality <- "ERROR"
		if (compare$sim$logBF[nrow(compare$sim)] <= logBoundary[1]) {
			outcome <- "not better"
			decision <- "stay with old standard treatment"
			decisionQuality <- ifelse(treatment > standardTreatment.ES, "false negative", "true negative")
		} else if (compare$sim$logBF[nrow(compare$sim)] >= logBoundary[2]) {
			outcome <- "better"
			decision <- "switch to new standard treatment"
			decisionQuality <- ifelse(treatment > standardTreatment.ES, "true positive", "false positive")
		} else if (compare$sim$n[nrow(compare$sim)] >= min(n.max, N.left)) {
			outcome <- "n.max reached, inconclusive"
			decision <- "stay with old standard treatment"
			decisionQuality <- ifelse(treatment > standardTreatment.ES, "false negative", "true negative")
		}
	
	
		if (outcome == "better") {
			if (verbose==TRUE) {
				cat(paste0("Outcome = ", outcome, ". Switching to new standard arm with ES = ", round(treatment, 2), ", which is a ", ifelse(treatment > standardTreatment.ES, "true positive", "false positive"), ".\n\n"))
			}
		
			standardTreatment.ES <- treatment
		}
	
		N.used <- N.used + compare$sim$n[nrow(compare$sim)]
	
		hops <- rbind(hops, data.frame(
			arm = armCounter,
			tested.ES = treatment,
			n.used = compare$sim$n[nrow(compare$sim)],
			N = N.used,
			logBF = compare$sim$logBF[nrow(compare$sim)],
			BF = exp(compare$sim$logBF[nrow(compare$sim)]),
			emp.ES = compare$sim$emp.ES[nrow(compare$sim)],
			outcome = outcome,
			decision = decision,
			standard.ES = standardTreatment.ES,
			decisionQuality = decisionQuality			
		))
		
		if ((N.used >= overall.N) | ((overall.N - N.used) < n.min)) {
			if (verbose==TRUE) {
				cat(paste0("All available participants used. Final effect size is ", round(standardTreatment.ES, 2), "\n"))
			}
			break;
		}
	}

	return(hops)
}




## ======================================================================
## Test a design space
## ======================================================================

params <- expand.grid(
	upper.boundary = c(3, 6, 10),
	lower.boundary = c(1/3, 1/6, 1/10),
	n.min = c(10, 20, 40),
	n.max = c(100, 150, 200)
)


set.seed(1234)
standardTreatment.ES <- 0.4
B <- 1000

#
# newTreatment.ES <- rtrunc(100, spec="norm", a=0.12, b = 0.72, mean = 0.42, sd = 0.15)
#
# x <- seq(0, 1, by=0.01)
# plot(x, dtrunc(x, spec="norm", a=0.12, b = 0.72, mean = 0.42, sd = 0.1), type="l", xlab="Cohen's d", ylab="")

registerDoParallel(cores=20)	

fullRes <- list()

# loop through design space
for (p in 1:nrow(params)) {
	print(paste0(p, "/", nrow(params)))
	flush.console()
	
	# repeat B times, split across CPUs
	sim <- foreach(batch=1:getDoParWorkers(), .combine=bind_rows) %dopar% {
		
		max_b <- round(B/getDoParWorkers())
		res0 <- list()
		
		for (b in 1:max_b) {
			LF <- leapfrog(standardTreatment.ES = 0.4, newTreatment.ES = rtrunc(100, spec="norm", a=0.12, b = 0.72, mean = 0.42, sd = 0.15), overall.N = 500, n.min = params[p, "n.min"], n.max = params[p, "n.max"], boundary = c(params[p, "lower.boundary"], params[p, "upper.boundary"]), verbose=FALSE, stepsize=5, rscale=0.4)
	
			res0[[b]] <- data.frame(
				batch = batch, 
				iteration = b,
				params[p, ],
				LF
			)		
		}
		
		res <- bind_rows(res0)
		return(res)
	} # of %dorng%
	
	fullRes[[p]] <- sim
	save(fullRes, file="fullRes.RData")
	
	print(mean(fullRes[[p]]$standard.ES))
	flush.console()
}

