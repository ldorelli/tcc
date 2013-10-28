pdf ("fases.pdf")

for (i in 1:100)
{
	str = paste (i-1, ".theta", sep='')
	tab = read.table (str)
	x = runif(1, 0, 1)
	y = runif(1, 0, 1)
	z = runif(1, 0, 1)
	if (i == 1) {
        plot (tab[,1], tab[,2], ylim=c(-3.14,3.14), col=rgb(x,y,z,0.4), cex=.25, xlab="t", ylab="Fase")
    }
	else {
        points (tab[,1], tab[,2], ylim=c(-3.14,3.14), col=rgb(x,y,z,0.4), cex=.25, xlab="t", ylab="Fase")
    }
}

dev.off ()
