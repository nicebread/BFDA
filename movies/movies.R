makeMovie <- function(BFDA, n.start, n.stop, by, fname="d0.5", fps=15, ...) {
	# create a new directory for the pictures
	if (!dir.exists(fname)) dir.create(fname)

	# create the picture sequence
	picName <- paste(fname, "/", fname, "_%03d.jpg", sep="")
	jpeg(picName, width=800, height=500, quality=95)
	for (n in seq(n.start, n.stop, by=by)) {
		print(n)
		plot(BFDA, n.max=n, ...)
	}
	graphics.off()

	# delete any existing movie file
	unlink(paste(fname, "/", fname,".avi",sep=""))

	# point system to R's working directory
	system(paste("cd ", gsub(" ", "\\ ", getwd(), fixed=TRUE)))

	# show & execute the command line expression for ffmpeg to glue the pictures together
	print(paste(paste0("ffmpeg -r ", fps, " -i ", fname, "/", fname, "_%03d.jpg -qscale 0 -r 25 ",  paste0(fname, "/", fname, ".avi"))))
	system(paste(paste0("ffmpeg -r ", fps, " -i ", fname, "/", fname, "_%03d.jpg -qscale 0 -r 25 ", paste0(fname, "/", fname, ".avi"))))
}

load("../finalSims/BFDA.0.5.RData")

makeMovie(sim, n.start=30, n.stop=300, by=1, fname="d0.5", fps=10, boundary=Inf, xlim=c(10, 300), ylim=c(log(1/40), 25), xextension=1.6, yaxis.at=c(-log(30), -log(10), -log(3), log(1), log(3), log(10), log(30), log(100), log(1000), log(10000)), yaxis.labels=c("1/30", "1/10", "1/3", "1", "3", "10", "30", "100", "1000", "10000"))


makeMovie(sim, n.start=30, n.stop=250, by=1, fname="d0.5.bound10", fps=10, boundary=10, xlim=c(10, 270), ylim=c(log(1/40), log(40)), traj.selection="fixed", dens.amplification=1.7, xextension=1.2, n.max.label.position="fixed", dens.right.offset=2)

# GIF animation (less pictures)
makeMovie(BFDA.0.5, n.start=30, n.stop=170, by=3, fname="GIF1", fps=10, boundary=10, xlim=c(10, 270), ylim=c(log(1/40), log(40)), traj.selection="fixed", dens.amplification=1.7, xextension=1.2, n.max.label.position="fixed", dens.right.offset=2, cex.labels=1.1, cex.annotations=0.95)