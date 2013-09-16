genConfNorm <- function (file, pop, t0, sigma, np, step) {
	cat (paste(t0, sigma, np, step, sep=' '), file=file, sep="\n")
	for (i in 1:pop) {
		x = 0
		while (x == 0) {
			x = rnorm(1, 0, 1)
		}
		cat (paste(rnorm (1, pi, pi/3), x, sep=' '), file = file, append = TRUE, sep="\n")
	}
}
