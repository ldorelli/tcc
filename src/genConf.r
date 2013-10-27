genConfNorm <- function (file, pop,  t0, np) {
	cat (paste(t0, np, sep=' '), file=file, sep="\n")
	for (i in 1:pop) {
		x = 0
		x = rnorm(1, 1.5, 1)
		x = x 
		cat (paste(rnorm (1, 0, pi/6), x, sep=' '), file = file, append = TRUE, sep="\n")
	}
}

genConfUnif <- function (file, pop, t0, np) {
	cat (paste(t0, np, sep=' '), file=file, sep="\n")
	for (i in 1:pop) {
		x = 0
		while (x == 0) {
			x = runif(1, -1, 1)
		}
		cat (paste(runif (1, 0, 1), x, sep=' '), file = file, append = TRUE, sep="\n")
	}
}
