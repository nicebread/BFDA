# does a value lie in between two others?
inside <- function(x, R) Vectorize({x >= R[1] & x <= R[2]})

# get population with a specific standardized mean differences
get_population <- function(n, d, SD=1) {
	x <- rnorm(n, 0, SD)
	y <- rnorm(n, 0, SD) - d
	return(cbind(x, y))
}
