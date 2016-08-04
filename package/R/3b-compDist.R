#' Plot a BFDA object for fixed n: Two densities for H1 and H0
#' @export
#' @importFrom scales alpha

# TODO: Add link to proper function name
#' @param BFDA.H1 A BFDA simulation object for H1, resulting from the \code{BFDA.sim} function
#' @param BFDA.H0 A BFDA simulation object for H0, resulting from the \code{BFDA.sim} function
#' @param n Fixed n at which BF distribution is evaluated
#' @param boundary Critical BF boundaries for H0 and H1
#' @param xlim limits on xaxis
#' @param noSplit If TRUE, the par(mfcol=...) command is skipped
#' @param cex Zoom factor for axes labels
#' @param cex.axis Zoom factor for tick labels
#' @param bw Black/white color scheme? (default:FALSE)

compDist <- function(BFDA.H1, BFDA.H0, n, boundary=c(1/6, 6), xlim=NA, noSplit=FALSE, cex=1.2, cex.axis=1, bw=FALSE) {

	boundary <- log(boundary)
	xlim <- log(xlim)

	N <- n	# rename, otherwise dplyr breaks ...
	logBF.H0 <-  BFDA.H0$sim %>% filter(n==N) %>% .$logBF
	logBF.H1 <-  BFDA.H1$sim %>% filter(n==N) %>% .$logBF
	allBF <- c(logBF.H0, logBF.H1)
	MAX <- quantile(allBF, prob=.95)
	
	ylim <- c(0, 1)	
	if (all(is.na(xlim))) xlim <- c(min(allBF), MAX)

	dens.H0 <- density(logBF.H0, from=xlim[1], to=xlim[2])
	dens.H1 <- density(logBF.H1, from=xlim[1], to=xlim[2])
	
	maxDens <- max(c(dens.H0$y, dens.H1$y))
	
	# normalize densities to max=1
	dens.H0$y <- dens.H0$y/maxDens
	dens.H1$y <- dens.H1$y/maxDens

	# split screen for two plots - omit for special layouts!
	if (noSplit==FALSE) par(mfrow=c(2, 1), mar=c(4, 2.3, 1.8, 1))

	# ---------------------------------------------------------------------
	# H1

	plot(NA, xlim=xlim, ylim=ylim, xlab="", ylab="", bty="n", axes=FALSE, main=expression(Under~H[1]))
	mtext("Density", side=2, line=1, cex=cex)

	# axes
	# automatic labeling of y-axis
	xaxis.at <- c(-log(30), -log(10), -log(3), log(1), log(3), log(10), log(30))
	xaxis.labels <- c("1/30", "1/10", "1/3", "1", "3", "10", "30")
	i <- 2
	repeat {
		if (xlim[2] >= log(10^i)) {
			xaxis.at <- c(xaxis.at, log(10^i))
			xaxis.labels <- c(xaxis.labels, as.character(10^i))
			i <- i+1
		} else {break;}
	}

	# set scale ticks
	axis(1, at = xaxis.at[inside(xaxis.at, c(-Inf, xlim[2]))],  labels=xaxis.labels[inside(xaxis.at, c(-Inf, xlim[2]))], las=1, cex.axis=cex.axis)	
	
	abline(v=boundary, lty="dashed")
	abline(v=log(1), lty="dotted")

	# Get the axis ranges, draw y-axis with arrow
	u <- par("usr")
	points(u[1], u[4], pch=17, xpd = TRUE)
	lines(c(u[1], u[1]), c(u[3], u[4]), xpd = TRUE)

	lines(dens.H1)

	if (bw==FALSE) {
		colors <- c("green4", "orange1", "red3")
	} else {
		colors <- c("grey80", "grey50", "grey10")
	}
	
	drawpoly(dens.H1, -Inf, boundary[1], col=scales::alpha(colors[3], 0.4))	
	if (length(boundary) == 2) {
		drawpoly(dens.H1, boundary[1], boundary[2], col=scales::alpha(colors[2], 0.4))
		drawpoly(dens.H1, boundary[2], Inf, col=scales::alpha(colors[1], 0.4))
	} else {
		drawpoly(dens.H1, boundary[1], Inf, col=scales::alpha(colors[1], 0.4))
	}
	

	# ---------------------------------------------------------------------
	# Under H0

	plot(NA, xlim=xlim, ylim=ylim, xlab="", ylab="", bty="n", axes=FALSE, main=expression(Under~H[0]))
	mtext("Density", side=2, line=1, cex=cex)
	mtext(expression(Bayes~factor~(BF[10])), side=1, line=2.7, cex=cex)

	# axes
	# automatic labeling of y-axis
	xaxis.at <- c(-log(30), -log(10), -log(3), log(1), log(3), log(10), log(30))
	xaxis.labels <- c("1/30", "1/10", "1/3", "1", "3", "10", "30")
	i <- 2
	repeat {
		if (xlim[2] >= log(10^i)) {
			xaxis.at <- c(xaxis.at, log(10^i))
			xaxis.labels <- c(xaxis.labels, as.character(10^i))
			i <- i+1
		} else {break;}
	}

	# Get the axis ranges, draw y-axis with arrow
	u <- par("usr")
	points(u[1], u[4], pch=17, xpd = TRUE)
	lines(c(u[1], u[1]), c(u[3], u[4]), xpd = TRUE)


	# set scale ticks
	axis(1, at = xaxis.at[inside(xaxis.at, c(-Inf, xlim[2]))],  labels=xaxis.labels[inside(xaxis.at, c(-Inf, xlim[2]))], las=1, cex.axis=cex.axis)
	abline(v=boundary, lty="dashed")
	abline(v=log(1), lty="dotted")


	# The densities
	lines(dens.H0)
	drawpoly(dens.H0, -Inf, boundary[1], col=scales::alpha(colors[1], 0.4))	
	if (length(boundary) == 2) {
		drawpoly(dens.H0, boundary[1], boundary[2], col=scales::alpha(colors[2], 0.4))
		drawpoly(dens.H0, boundary[2], Inf, col=scales::alpha(colors[3], 0.4))
	} else {
		drawpoly(dens.H0, boundary[1], Inf, col=scales::alpha(colors[3], 0.4))
	}
}

