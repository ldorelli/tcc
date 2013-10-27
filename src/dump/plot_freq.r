pdf ("frequencias.pdf")

for (i in 1:100)
{
	str = paste (i-1, ".freq", sep='')
	tab = read.table (str)
	if (i == 1) {
        plot (tab[,1], tab[,2], col=i, cex=.25, xlab="t", ylab="fase")
    }
	else {
        points (tab[,1], tab[,2], col=i, cex=.25, xlab="t", ylab="fase")
    }
}

dev.off ()
