# TODO: add.density.at = c(20, 40, 60); NA = default = n.max

#' Plot a BPA object
#' @aliases plot.BPA
#' @export
#' @importFrom scales alpha
#' @importFrom TeachingDemos shadowtext

#' @param BPA A BPA simulation object
#' @param forH1 If TRUE, use the BF_10 labels, otherwise the BF_01 labels
#' @param dens.right.offset How much should the right density be moved to the right?
#' @param dens.amplification Amplification factor for all densities.
#' @param n.max Maximum sample size
#' @param boundary The BF threshold where to stop (resp. its reciprocal)
#' @param n.trajectories Number of demo trajectories
#' @param xlim Define limits on x axis
#' @param ylim Define limits on y axis
#' @param xextension How much reltive space should be on the rights side of the plot? You can adjust this to change the distance of the BF category labels from the plot.
#' @param cat how many BF categories should be colored in the right density? Either 3 or 6.
#' @param yaxis.at Positions of the ticks on the y-axis. If NA, they are defined automatically.
#' @param yaxis.labels Labels of the ticks on the y-axis. If NA, they are defined automatically.
#' @param bw Black/white density on the right?
#' @param traj.selection Should a fixed set of trajectories be shown ("fixed"), or a selection that reflects the proportional of each stopping category (upper/loewr/n.max hits; use "proportional").

#forH1 = TRUE; boundary=6; n.trajectories=60; n.max=500; dens.amplification=NA; dens.right.amplification=NA; plotratio=NA; cat=3; dens.right.offset=2; xlim=NA; ylim=NA; load("finalSims/sim.0.5.RData")



plotBPA <- function(BPA, boundary=30, n.trajectories=60, n.max=NA, dens.amplification=1, cat=3, bw=FALSE, dens.right.offset=5, xlim=NA, ylim=NA, xextension=1.8, traj.selection="proportional", yaxis.at=NA, yaxis.labels=NA, forH1=TRUE) {
	sim <- BPA$sim
	traj.selection <- match.arg(traj.selection, c("proportional", "fixed"))
		
	plotratio <- NA
	# ---------------------------------------------------------------------
	# Prepare data for plot

	logBoundary <- log(boundary)
	if (is.na(n.max)) {n.max.actual <- max(sim$n)} else {n.max.actual <- n.max}
		
	indices <- BPA.analysis(BPA, boundary=boundary, n.max=n.max.actual, verbose=FALSE)
	
	# reduce data frame to actually relevant data
	sim <- sim %>% filter(n <= n.max.actual)

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
		mutate(firstbreak=unNA(which(abs(logBF)>=logBoundary)[1])) %>% 
		filter(row_number()<=firstbreak) %>% ungroup()

	# "Land" all final points on the boundary (looks better in plot)
	demo2$logBF[demo2$logBF > logBoundary] <- logBoundary
	demo2$logBF[demo2$logBF < -logBoundary] <- -logBoundary

	final_point_boundary <- demo2 %>% group_by(id) %>% filter(n == max(n), abs(logBF)==logBoundary)
	final_point_n.max <- demo2 %>% filter(n == n.max.actual, logBF<logBoundary)
	
	
	## Compute densities & normalize densities, so that the areas sum to one
	# Compute upper density. Weigh by frequency, so that all three densities sum up to 1	
	if (!is.null(indices$d.top)) {
		d.top.area <- sum(diff(indices$d.top$x)*indices$d.top$y)
		dens.top <- as.data.frame(indices$d.top[c("x", "y")])
		dens.top <- dens.top[dens.top$y > 0.0001, ]
		# normalize density to area 1
		dens.top$y <- dens.top$y/d.top.area
		# weigh density
		dens.top$y <- dens.top$y*indices$upper.hit.frac	
		#sum(diff(dens.top$x)*dens.top$y)
	} else {dens.top <- NULL}

	# Compute lower density
	if (!is.null(indices$d.bottom)) {
		d.bottom.area <- sum(diff(indices$d.bottom$x)*indices$d.bottom$y)
		dens.bottom <- as.data.frame(indices$d.bottom[c("x", "y")])
		dens.bottom <- dens.bottom[dens.bottom$y > 0.0001, ]
		# normalize density to area 1
		dens.bottom$y <- dens.bottom$y/d.bottom.area
		dens.bottom$y <- dens.bottom$y*indices$lower.hit.frac # weigh density
		#sum(diff(dens.bottom$x)*dens.bottom$y)
	} else {dens.bottom <- NULL}

	# Compute right density: Distribution of logBF	
	if (!is.null(indices$d.right)) {
		d.right.area <- sum(diff(indices$d.right$x)*indices$d.right$y)
		dens.right <- as.data.frame(indices$d.right[c("x", "y")])
		dens.right <- dens.right[dens.right$y > 0.0001, ]
		# normalize density to area 1
		dens.right$y <- dens.right$y/d.right.area
		dens.right$y <- dens.right$y*indices$n.max.hit.frac	# weigh density
		#sum(diff(dens.right$x)*dens.right$y)
	} else {dens.right <- NULL}
	
	# 
	
	# ---------------------------------------------------------------------
	#  Compute amplification factors: 
	
	# xlim: 1.5 times the boundary to leave space for densities and category labels
	if (all(is.na(ylim))) ylim <- c(max(min(sim$logBF), -logBoundary)*1.7, min(max(sim$logBF), logBoundary)*1.7)*1.1		
	if (all(is.na(xlim))) {
		xlim <- c(min(sim$n), max(sim$n))
	}

	
	# automatic labeling of y-axis
	if (is.na(yaxis.at)) {
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
	plot(NA, xlim=c(xlim[1], xlim[2]*xextension), ylim=ylim, xlab="", ylab="", bty="n", axes=FALSE)
	
	# nice labels
	mtext("Sample size", side=1, line=2.5, cex=1.2)						# xlab
	mtext(expression(Bayes~factor~(BF[10])), side=2, line=3, cex=1.2)	# ylab
	
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
	if (is.finite(boundary)) {
		lines(x=c(xlim[1], n.max.actual), y=c(-logBoundary, -logBoundary), lty="solid", lwd=1.2)
		lines(x=c(xlim[1], n.max.actual), y=c(+logBoundary, +logBoundary), lty="solid", lwd=1.2)
	}
	
	## Annotation: horizontal lines at BF categories
	for (y in c(c(-log(c(30, 10, 3)), 0, log(c(3, 10, 30))))) {
		if (inside(y, ylim))
			lines(x=c(xlim[1], xlim[2]*xextension), y=rep(y, 2), lty="dotted", col="grey20")
	}	
			
	
	# ---------------------------------------------------------------------
	#  Draw densities
	
	
	# get aspect ratio so that densities are correctly scaled
	w <- par("pin")[1]/diff(par("usr")[1:2])
	h <- par("pin")[2]/diff(par("usr")[3:4])
	asp <- w/h
	
	
	# The highest density point on top/bottom should be maximally as high as the boundary/2
	# AND the rightmost density peak should be maximally 1/2 the range of n
	highest.dens <- max(c(ifelse(is.null(dens.bottom), NA, dens.bottom$y), ifelse(is.null(dens.top), NA, dens.top$y), ifelse(is.null(dens.right), NA, dens.right$y)), na.rm=TRUE)
	
	# compute density amplification; multiply with general factor
	if (!is.na(highest.dens)) {
		dens.amp.boundary <-  ((logBoundary/2)/highest.dens)*dens.amplification
		dens.amp.n.max <- (dens.amp.boundary/asp)*dens.amplification
	}
	
	#stop();


	# Upper density
	if (!is.null(dens.top)) {
		poly.top <- rbind(data.frame(x=dens.top$x[1], y=0), dens.top, data.frame(x=max(dens.top$x), y=0))
		polygon(poly.top$x, poly.top$y*dens.amp.boundary + logBoundary, col="grey80")
		lines(dens.top$x, dens.top$y*dens.amp.boundary + logBoundary, col="grey30")
	}

	# Lower density
	if (!is.null(dens.bottom)) {
		poly.bottom <- rbind(data.frame(x=dens.bottom$x[1], y=0), dens.bottom, data.frame(x=max(dens.bottom$x), y=0))
		polygon(poly.bottom$x, -poly.bottom$y*dens.amp.boundary - logBoundary, col="grey80")
		lines(dens.bottom$x, -dens.bottom$y*dens.amp.boundary - logBoundary, col="grey30")
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
			colors <- c("grey60", "grey80", "grey60", "grey80", "grey60", "grey80")
		if (bw==TRUE & cat==3) 
			colors <- c("grey60", "grey80", "grey60")
		if (bw==FALSE & cat==6) 
			colors <- c("green4", "olivedrab3", "orange1", "darkorange2", "orangered2", "red3")
		if (bw==FALSE & cat==3) 
			colors <- c("green4", "orange1", "red3")
		

		# color polygons of the distribution
		
		for (i in 1:cat) {
			poly <- dens.right %>% filter(x < bound[i], x > bound[i+1])
			
			if (is.finite(poly[1, 1])) {
				poly <- rbind(data.frame(x=min(poly$x), y=0), poly, data.frame(x=max(poly$x), y=0))
				polygon(poly$y*dens.amp.n.max + n.max.actual + dens.right.offset, poly$x, col=colors[i], border=NA) 
				
				if (stoppingBF.perc[i] >= 1) {
					# find the point in the density which splits the area into two equal halves				
					median.area.pos <- which((cumsum(poly$y) > sum(poly$y)/2) == TRUE)[1]
					X <- poly$y[median.area.pos]
					Y <- poly$x[median.area.pos]				
					
					# shadowtext adds a white contour around the characters
					TeachingDemos::shadowtext(x=n.max.actual + dens.right.offset*1.3 + X*dens.amp.n.max, y=Y, labels=paste0(stoppingBF.perc[i], "%"), bg="white", col="black", adj=-0.3, cex=0.8, r=0.1)
				}
			}
		} # of i

		# draw the right density
		if (!is.null(dens.right)) {
			lines(dens.right$y*dens.amp.n.max + n.max.actual + dens.right.offset, dens.right$x, col="grey30")
		}
	
	
		# n.max boundary
		if (!is.na(n.max)) {lines(x=rep(n.max, 2), y=c(-logBoundary, logBoundary), lty="solid", lwd=1.2)}
	
		# add final points at n.max hit
		points(final_point_n.max$n, final_point_n.max$logBF, pch=16, cex=.8)
	} # of if (indices$n.max.hit.frac > 0)


	# add final points at boundary hit
	points(final_point_boundary$n, final_point_boundary$logBF, pch=16, cex=.8)
	

	par(xpd=NA)	# allow drawing outside the plot
	if (inside(labels.y[8], ylim))
		text(xlim[2]*xextension, labels.y[8], bquote(Very~strong~H[.(ifelse(forH1==TRUE,0,1))]), adj=c(1, 0.5), cex=.8)
	if (inside(labels.y[7], ylim))
		text(xlim[2]*xextension, labels.y[7], bquote(Strong~H[.(ifelse(forH1==TRUE,0,1))]), adj=c(1, 0.5), cex=.8)
	if (inside(labels.y[6], ylim))
		text(xlim[2]*xextension, labels.y[6], bquote(Moderate~H[.(ifelse(forH1==TRUE,0,1))]), adj=c(1, 0.5), cex=.8)
	if (inside(labels.y[5], ylim))
		text(xlim[2]*xextension, labels.y[5], bquote(Anecdotal~H[.(ifelse(forH1==TRUE,0,1))]), adj=c(1, 0.5), cex=.8)
	if (inside(labels.y[1], ylim))
		text(xlim[2]*xextension, labels.y[1], bquote(Very~strong~H[.(ifelse(forH1==TRUE,1,0))]), adj=c(1, 0.5), cex=.8)
	if (inside(labels.y[2], ylim))
		text(xlim[2]*xextension, labels.y[2], bquote(Strong~H[.(ifelse(forH1==TRUE,1,0))]), adj=c(1, 0.5), cex=.8)
	if (inside(labels.y[3], ylim))
		text(xlim[2]*xextension, labels.y[3], bquote(Moderate~H[.(ifelse(forH1==TRUE,1,0))]), adj=c(1, 0.5), cex=.8)
	if (inside(labels.y[4], ylim))
		text(xlim[2]*xextension, labels.y[4], bquote(Anecdotal~H[.(ifelse(forH1==TRUE,1,0))]), adj=c(1, 0.5), cex=.8)
	par(xpd=TRUE)
	
	
	# Write labels of stopping percentages
	TeachingDemos::shadowtext(col="black", bg="white", x=xlim[1], y=logBoundary, label=bquote(paste(.(round(indices$upper.hit.frac*100)), "% stopped at ", H[1], " boundary")), cex=.8, adj=c(-0.1, -0.3))
	TeachingDemos::shadowtext(col="black", bg="white", x=xlim[1], y=-logBoundary, label=bquote(paste(.(round(indices$lower.hit.frac*100)), "% stopped at ", H[0], " boundary")), cex=.8, adj=c(-0.1, 1.3))
	
	if (indices$n.max.hit.frac > 0.005 & is.finite(boundary)) {		
		# find maximum of dens.right; center the label on this point
		Y <- dens.right$x[which.max(dens.right$y)]
		TeachingDemos::shadowtext(col="black", bg="white", x=n.max.actual, y=Y, label=bquote(paste(.(round(indices$n.max.hit.frac*100)), "% stopped at ", n[max], " boundary")), cex=.8, srt=90, adj=c(0.5, 1.3))
	}
}


#' @export
plot.BPA <- function(x, ...) {
	plotBPA(x, ...)
}

