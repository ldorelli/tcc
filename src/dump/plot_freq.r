pdf ("frequencias.pdf")

for (i in 1:100)
{
	str = paste (i-1, ".freq", sep='')
	tab = read.table (str)
	x = runif(1, 0,1)
	y = runif(1, 0,1)
	z = runif(1, 0,1)
	if (i == 1) {
        plot (tab[,1], tab[,2], col=rgb(x,y,z, 0.5), ylim=c(0, 3.14), cex=.25, xlab="t", ylab="Frequência")
    }
	else {
        points (tab[,1], tab[,2], col=rgb(x, y, z, 0.5), cex=.25, ylim=c(0, 3.14), xlab="t", ylab="Frequência") 
    }
}

dev.off ()
