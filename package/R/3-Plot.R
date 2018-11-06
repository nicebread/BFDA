#' Plot a BFDA object
#' @aliases plot.BFDA
#' @export
#' @importFrom scales alpha
#' @importFrom TeachingDemos shadowtext
#' @importFrom sfsmisc integrate.xy

# TODO: Add link to proper function name
#' @param BFDA A BFDA simulation object, resulting from the \code{\link{BFDA.sim}} function
#' @param forH1 If TRUE, use the BF_10 labels, otherwise the BF_01 labels
#' @param dens.right.offset How much should the right density (at n.max) be moved to the right?
#' @param dens.amplification Amplification factor for all densities.
#' @param n.min Minimum sample size
#' @param n.max Maximum sample size
#' @param boundary At which BF boundary should sequential trajectories stop? Either a single number (then the reciprocal is taken as the other boundary), or a vector of two numbers for lower and upper boundary
#' @param n.trajectories Number of demo trajectories that are plotted
#' @param xlim Define limits on x axis
#' @param ylim Define limits on y axis
#' @param xextension How much reltive space should be on the rights side of the plot? You can adjust this to change the distance of the BF category labels from the plot.
#' @param cat how many BF categories should be colored in the right density? Either 3 or 6.
#' @param yaxis.at Positions of the ticks on the y-axis. If NA, they are defined automatically.
#' @param yaxis.labels Labels of the ticks on the y-axis. If NA, they are defined automatically.
#' @param bw Black/white density on the right?
#' @param traj.selection Should a fixed set of trajectories be shown ("fixed"), or a selection that reflects the proportional of each stopping category (upper/lower/n.max hits; use "proportional").
#' @param n.max.label.position "fixed": Always centered at BF=1. "dynamic": Centered on the the peak of the right density.
#' @param cex.labels Zoom factor for axes text labels
#' @param cex.annotations Zoom factor for text annotations
#' @param ... Additional parameters passed through to the base plot function

#forH1 = TRUE; boundary=6; n.trajectories=60; n.max=500; dens.amplification=NA; dens.right.amplification=NA; plotratio=NA; cat=3; dens.right.offset=2; xlim=NA; ylim=NA; load("finalSims/sim.0.5.RData")

plotBFDA <- function(BFDA, boundary=10, n.trajectories=60, n.min=NA, n.max=NA, dens.amplification=1, cat=3, bw=FALSE, dens.right.offset=2, xlim=NA, ylim=NA, xextension=1.5, traj.selection="proportional", yaxis.at=NA, yaxis.labels=NA, forH1=TRUE, n.max.label.position="dynamic", cex.labels=1, cex.annotations=0.85, ...) {
	sim <- BFDA$sim
	traj.selection <- match.arg(traj.selection, c("proportional", "fixed"))
	n.max.label.position <- match.arg(n.max.label.position, c("dynamic", "fixed"))
	
	if (length(boundary) == 1) boundary <- sort(c(boundary, 1/boundary))
	logBoundary <- log(boundary)
	
		
	plotratio <- NA
	# ---------------------------------------------------------------------
	# Prepare data for plot


	if (is.na(n.max)) {
		n.max.compute <- n.max <- max(sim$n) # n.max.compute stores the maximum n for computation, n.max the display n.max
		if (!all(is.na(xlim))) {
			n.max <- min(xlim[2], max(sim$n))
		}
	} else {
		n.max.compute <- n.max
		if (!all(is.na(xlim))) {
			n.max <- min(xlim[2], n.max)
		}		
	}
	if (is.na(n.min)) {n.min <- min(sim$n)}
	
		
	indices <- BFDA.analyze(BFDA, boundary=boundary, design="sequential", n.min=n.min, n.max=n.max.compute, verbose=FALSE)
	
	# reduce data frame to actually relevant data
	sim <- sim %>% filter(n >= n.min & n <= n.max)

	# select some trajectories for displaying, relative for each stopping outcome
	
	set.seed(0xBEEF)	# for reproducibility: always sample the same trajectories
	if (traj.selection=="proportional") {
		seltraj <- 
		c(	as.character(sample(indices$upper.hit.ids, size=round(indices$upper.hit.frac*n.trajectories))),
		  	as.character(sample(indices$lower.hit.ids, size=round(indices$lower.hit.frac*n.trajectories))),
			as.character(sample(indices$n.max.hit.ids, size=round(indices$n.max.hit.frac*n.trajectories))))
		demo1 <- sim %>% filter(id %in% seltraj)
	}
	if (traj.selection=="fixed") {
		demo1 <- sim[sim$id %in% sample(unique(sim$id), n.trajectories), ]
	}

	# delete all points after first boundary hit
	unNA <- function(x) {x[is.na(x)]=Inf; return(x)}
	demo2 <- demo1 %>% group_by(id) %>% 
		mutate(firstbreak=unNA(which(logBF >= logBoundary[2] | logBF <= logBoundary[1])[1])) %>% 
		filter(row_number()<=firstbreak) %>% ungroup()

	# "Land" all final points on the boundary (looks better in plot)
	demo2$logBF[demo2$logBF > logBoundary[2]] <- logBoundary[2]
	demo2$logBF[demo2$logBF < logBoundary[1]] <- logBoundary[1]

	final_point_boundary <- demo2 %>% group_by(id) %>% filter(n == max(n), logBF==logBoundary[1] | logBF==logBoundary[2])
	final_point_n.max <- demo2 %>% filter(n == n.max, logBF<logBoundary[2] & logBF > logBoundary[1])
	
	
	# xlim: 1.5 times the boundary to leave space for densities and category labels
	if (all(is.na(ylim))) ylim <- c(max(min(sim$logBF), logBoundary[1])*1.7, min(max(sim$logBF), logBoundary[2])*1.7)*1.1		
	if (all(is.na(xlim))) {
		xlim <- c(min(sim$n), max(sim$n))
	}
	
		
	# stop();
	# indices$d.top <- density(runif(10000, 20, 40))
	# indices$d.bottom <- density(runif(10000, 20, 30))
	# indices$d.right <- density(runif(10000, log(1), log(3)))
	# indices$upper.hit.frac <- 0.6
	# indices$lower.hit.frac <- 0.2
	# indices$n.max.hit.frac <- 0.2
	
	## Compute densities & normalize densities, so that the areas sum to one
	# Compute upper density. Weigh by frequency, so that all three densities sum up to 1	
	#stop()
	if (!is.null(indices$d.top)) {
		d.top.area <- integrate.xy(indices$d.top$x, indices$d.top$y)
		dens.top <- as.data.frame(indices$d.top[c("x", "y")])
		dens.top <- dens.top[dens.top$y > 0.0001, ]
		# normalize density to area 1
		dens.top$y <- dens.top$y/d.top.area		
		dens.top$y <- dens.top$y*indices$upper.hit.frac	# weigh density
		dens.top <- dens.top[dens.top$x < n.max, ]
		#sum(diff(dens.top$x)*dens.top$y)
	} else {dens.top <- NULL}

	# Compute lower density
	if (!is.null(indices$d.bottom)) {
		d.bottom.area <- integrate.xy(indices$d.bottom$x, indices$d.bottom$y)
		dens.bottom <- as.data.frame(indices$d.bottom[c("x", "y")])
		dens.bottom <- dens.bottom[dens.bottom$y > 0.0001, ]
		# normalize density to area 1
		dens.bottom$y <- dens.bottom$y/d.bottom.area
		dens.bottom$y <- dens.bottom$y*indices$lower.hit.frac # weigh density
		dens.bottom <- dens.bottom[dens.bottom$x < n.max, ]
		#sum(diff(dens.bottom$x)*dens.bottom$y)
	} else {dens.bottom <- NULL}

	# Compute right density: Distribution of logBF	
	if (!is.null(indices$d.right)) {
		d.right.area <- integrate.xy(indices$d.right$x, indices$d.right$y)
		dens.right <- as.data.frame(indices$d.right[c("x", "y")])
		dens.right <- dens.right[dens.right$y > 0.0001, ]
		# normalize density to area 1
		dens.right$y <- dens.right$y/d.right.area
		dens.right$y <- dens.right$y*indices$n.max.hit.frac	# weigh density
		#sum(diff(dens.right$x)*dens.right$y)
	} else {dens.right <- NULL}
	
	# 
	#stop()
	
	
	# automatic labeling of y-axis
	if (all(is.na(yaxis.at))) {
		yaxis.at <- c(-log(30), -log(10), -log(3), log(1), log(3), log(10), log(30))
		yaxis.labels <- c("1/30", "1/10", "1/3", "1", "3", "10", "30")
		i <- 2
		repeat {
			if (ylim[2] >= log(10^i)) {
				yaxis.at <- c(yaxis.at, log(10^i))
				yaxis.labels <- c(yaxis.labels, as.character(10^i))
				i <- i+1
			} else {break;}
		}
	}
		
	stopifnot(length(yaxis.labels) == length(yaxis.at))
	
	

	## ======================================================================
	## The plot
	## ======================================================================

	# positions for the category labels
	labels.y <- c(4, 2.86, 1.7, .55, -.55, -1.7, -2.85, -4)
	
	par(mar=c(5, 5, 4, 6), xpd=TRUE)
	plot(NA, xlim=c(xlim[1], xlim[2]*xextension), ylim=ylim, xlab="", ylab="", bty="n", axes=FALSE, ...)
	
	# nice labels
	mtext("Sample size", side=1, line=2.5, cex=cex.labels)						# xlab
	mtext(expression(Bayes~factor~(BF[10])), side=2, line=3, cex=cex.labels)	# ylab
	
	# axes
	# set scale ticks
	ticks <- round(axTicks(1, c(par("xaxp")[1], xlim[2], par("xaxp")[3])))
	axis(1, at = ticks)
	axis(2, at = yaxis.at[inside(yaxis.at, ylim)],  labels=yaxis.labels[inside(yaxis.at, ylim)], las=2)	
	
	# draw the demo trajectories
	for (i in unique(demo2$id)) {
		lines(demo2$n[demo2$id==i], demo2$logBF[demo2$id==i], col=scales::alpha("grey30", 0.5))
	}
	
	# BF-boundaries
	if (all(is.finite(boundary))) {
		lines(x=c(xlim[1], n.max), y=c(logBoundary[1], logBoundary[1]), lty="solid", lwd=cex.labels)
		lines(x=c(xlim[1], n.max), y=c(logBoundary[2], logBoundary[2]), lty="solid", lwd=cex.labels)
	}
	
	## Annotation: horizontal lines at BF categories
	for (y in c(c(-log(c(30, 10, 3)), 0, log(c(3, 10, 30))))) {
		if (inside(y, ylim))
			lines(x=c(xlim[1], xlim[2]*xextension), y=rep(y, 2), lty="dotted", col="grey20")
	}	
			
	
	# ---------------------------------------------------------------------
	#  Draw densities
	
	
	# ---------------------------------------------------------------------
	#  Compute amplification factors: 
	
	# get aspect ratio so that densities are correctly scaled
	# TODO: Unnecessary?
	w <- par("pin")[1]/diff(par("usr")[1:2])
	h <- par("pin")[2]/diff(par("usr")[3:4])
	asp <- w/h
	
	
	if (!is.null(c(dens.bottom$y, dens.top$y))) {
		dens.amp <- 30*dens.amplification
	} else {
		dens.amp <- 300*dens.amplification
	}


	# Upper density
	if (!is.null(dens.top)) {
		poly.top <- rbind(data.frame(x=dens.top$x[1], y=0), dens.top, data.frame(x=max(dens.top$x), y=0))
		polygon(poly.top$x, (poly.top$y)*dens.amp + logBoundary[2], col="grey80")
		lines(dens.top$x, (dens.top$y)*dens.amp + logBoundary[2], col="grey30")
	}

	# Lower density
	if (!is.null(dens.bottom)) {
		poly.bottom <- rbind(data.frame(x=dens.bottom$x[1], y=0), dens.bottom, data.frame(x=max(dens.bottom$x), y=0))
		polygon(poly.bottom$x, -(poly.bottom$y)*dens.amp + logBoundary[1], col="grey80")
		lines(dens.bottom$x, -(dens.bottom$y)*dens.amp + logBoundary[1], col="grey30")
	}


	## Right side
	if (!is.null(dens.right) & indices$n.max.hit.frac > 0.005) {
		
		# Percentages of categories
		if (cat==3) {
			stoppingBF.perc <- rev(round(table(cut(indices$logBF.right, breaks=c(Inf, log(3), -log(3), -Inf)))/indices$all.traj.n, 3)*100)
			bound <- c(max(dens.right$x), log(3), -log(3), min(dens.right$x))
		}
		if (cat==6) {
			stoppingBF.perc <- rev(round(table(cut(indices$logBF.right, breaks=c(Inf, c(log(c(30, 10, 3)), 0, -log(c(3, 10, 30)), -Inf))))/indices$all.traj.n, 3)*100)
			bound <- c(max(dens.right$x), log(10), log(3), 0, -log(3), -log(10), min(dens.right$x))
		}

		if (bw==TRUE & cat==6) 
			colors <- c("grey80", "grey65", "grey50", "grey35", "grey20", "grey10")
		if (bw==TRUE & cat==3) 
			colors <- c("grey80", "grey50", "grey10")
		if (bw==FALSE & cat==6) 
			colors <- c("green4", "olivedrab3", "orange1", "darkorange2", "orangered2", "red3")
		if (bw==FALSE & cat==3) 
			colors <- c("green4", "orange1", "red3")
		

		# color polygons of the distribution
		
		for (i in 1:cat) {
			poly <- dens.right %>% filter(x < bound[i], x > bound[i+1])
			
			if (is.finite(poly[1, 1])) {
				poly <- rbind(data.frame(x=min(poly$x), y=0), poly, data.frame(x=max(poly$x), y=0))
				
				# flip x and y axis: the density (y) is plotted on the plot's x axis
				polygon(x=poly$y*dens.amp + n.max + dens.right.offset*(strheight("X")/asp), y=poly$x, col=colors[i], border=NA) 
				
				if (stoppingBF.perc[i] >= 1) {
					# find the point in the density which splits the area into two equal halves				
					median.area.pos <- which((cumsum(poly$y) > sum(poly$y)/2) == TRUE)[1]
					X <- poly$y[median.area.pos]
					Y <- poly$x[median.area.pos]				
					
					# shadowtext adds a white contour around the characters
					TeachingDemos::shadowtext(x=n.max + dens.right.offset*(strheight("X")/asp) + X*dens.amp, y=Y, labels=paste0(stoppingBF.perc[i], "%"), bg="white", col="black", adj=-0.3, cex=cex.annotations, r=0.1)
				}
			}
		} # of i

		# draw the right density
		if (!is.null(dens.right)) {
			lines(dens.right$y*dens.amp + n.max + dens.right.offset*(strheight("X")/asp), dens.right$x, col="grey30")
		}	
	
		# n.max boundary
		if (!is.na(n.max)) {lines(x=rep(n.max, 2), y=c(logBoundary[1], logBoundary[2]), lty="solid", lwd=cex.labels)}
	
		# add final points at n.max hit
		points(final_point_n.max$n, final_point_n.max$logBF, pch=16, cex=cex.annotations)
	} # of if (indices$n.max.hit.frac > 0)


	# add final points at boundary hit
	points(final_point_boundary$n, final_point_boundary$logBF, pch=16, cex=cex.annotations)
	

	if (!is.null(dens.right) & indices$n.max.hit.frac > 0.005) {
		xmax <- n.max + dens.right.offset*(strheight("X")/asp) + max(dens.right$y)*dens.amp + strwidth("XXXXXX")
	} else {
		xmax <- n.max + strheight("X")/asp
	}

	par(xpd=NA)	# allow drawing outside the plot
	if (inside(labels.y[8], ylim))
		text(xmax, labels.y[8], bquote(Very~strong~H[.(ifelse(forH1==TRUE,0,1))]), adj=c(0, 0.5), cex=cex.annotations)
	if (inside(labels.y[7], ylim))
		text(xmax, labels.y[7], bquote(Strong~H[.(ifelse(forH1==TRUE,0,1))]), adj=c(0, 0.5), cex=cex.annotations)
	if (inside(labels.y[6], ylim))
		text(xmax, labels.y[6], bquote(Moderate~H[.(ifelse(forH1==TRUE,0,1))]), adj=c(0, 0.5), cex=cex.annotations)
	if (inside(labels.y[5], ylim))
		text(xmax, labels.y[5], bquote(Anecdotal~H[.(ifelse(forH1==TRUE,0,1))]), adj=c(0, 0.5), cex=cex.annotations)
	if (inside(labels.y[1], ylim))
		text(xmax, labels.y[1], bquote(Very~strong~H[.(ifelse(forH1==TRUE,1,0))]), adj=c(0, 0.5), cex=cex.annotations)
	if (inside(labels.y[2], ylim))
		text(xmax, labels.y[2], bquote(Strong~H[.(ifelse(forH1==TRUE,1,0))]), adj=c(0, 0.5), cex=cex.annotations)
	if (inside(labels.y[3], ylim))
		text(xmax, labels.y[3], bquote(Moderate~H[.(ifelse(forH1==TRUE,1,0))]), adj=c(0, 0.5), cex=cex.annotations)
	if (inside(labels.y[4], ylim))
		text(xmax, labels.y[4], bquote(Anecdotal~H[.(ifelse(forH1==TRUE,1,0))]), adj=c(0, 0.5), cex=cex.annotations)
	par(xpd=TRUE)
	
	
	# Write labels of stopping percentages
	TeachingDemos::shadowtext(col="black", bg="white", x=xlim[1], y=logBoundary[2], label=bquote(paste(.(round(indices$upper.hit.frac*100)), "% stopped at ", H[1], " boundary")), cex=cex.annotations, adj=c(-0.1, -0.3))
	TeachingDemos::shadowtext(col="black", bg="white", x=xlim[1], y=logBoundary[1], label=bquote(paste(.(round(indices$lower.hit.frac*100)), "% stopped at ", H[0], " boundary")), cex=cex.annotations, adj=c(-0.1, 1.3))
	
	if (indices$n.max.hit.frac > 0.005 & all(is.finite(boundary))) {		
		if (n.max.label.position == "dynamic") {
			# find maximum of dens.right; center the label on this point
			Y <- dens.right$x[which.max(dens.right$y)]
		} else {
			Y <- log(1)
		}
		
		TeachingDemos::shadowtext(col="black", bg="white", x=n.max, y=Y, label=bquote(paste(.(round(indices$n.max.hit.frac*100)), "% stopped at ", n[max], " boundary")), cex=cex.annotations, srt=90, adj=c(0.5, 1.3))
	}
}


#' @export
plot.BFDA <- function(x, ...) {
	suppressWarnings(plotBFDA(x, ...))
}

