# does a value lie in between two others?
inside <- function(x, R) Vectorize({x >= R[1] & x <= R[2]})

# get population with a specific standardized mean differences
get_population <- function(n, d, SD=1) {
	x <- rnorm(n, 0, SD)
	y <- rnorm(n, 0, SD) - d
	return(cbind(x, y))
}


drawpoly <- function(dens, from, to, ...) {
	poly <- data.frame(x=dens$x, y=dens$y)
	poly <- poly %>% filter(x>from & x<to)
	poly <- rbind(c(x=min(poly$x), y=0), poly, c(x=max(poly$x), 0))
	polygon(poly$x, poly$y, ...)
}


# simple wrapper: formats a number in f.2 format
f2 <- function(x, digits=2, prepoint=0, skipZero=FALSE) {
	
	if (skipZero == TRUE) {zero <- "."} else {zero <- "0."}
	
	if (length(dim(x)) == 2) {
		apply(x, 2, function(x2) {gsub("0.", zero, sprintf(paste("%",prepoint,".",digits,"f",sep=""), x2) , fixed=TRUE)})
	} else {
		gsub("0.", zero, sprintf(paste("%",prepoint,".",digits,"f",sep=""), x) , fixed=TRUE)
	}
}