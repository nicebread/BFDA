makeMovie <- function(BPA, n.start, n.stop, by, fname="d0.5", fps=15, ...) {
	# create a new directory for the pictures
	dir.create(fname)

	# create the picture sequence
	picName <- paste(fname, "/", fname, "_%03d.jpg", sep="")
	jpeg(picName, width=800, height=500, quality=95)
	for (n in seq(n.start, n.stop, by=by)) {
		print(n)
		plot(BPA, n.max=n, ...)
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

load("../finalSims/sim0.5_step1.RData")
makeMovie(sim, n.start=30, n.stop=300, by=1, fname="d0.5", fps=10, boundary=Inf, xlim=c(10, 300), ylim=c(log(1/40), 25), xextension=1.6)

makeMovie(sim, n.start=30, n.stop=250, by=1, fname="d0.5.bound10", fps=10, boundary=10, xlim=c(10, 300), ylim=c(log(1/40), log(40)), traj.selection="fixed")
