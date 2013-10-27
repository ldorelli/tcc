V = commandArgs()
print(V)
x = as.numeric(V[3])
step = as.numeric(V[5])
start = as.numeric(V[4])
times = as.numeric(V[6])
name = paste (x, "_grap.pdf", sep='')
pdf (name)

symbols = c(0, 1, 18, 17, 15, 23)

print (x)

leg = NULL
for (i in 1:times)
{
    idx = start + step*(i-1);
    leg = c (leg, start+step*(i-1))
	str = paste (x, "_", idx, sep='')
	tab = read.table (str)
    print (start+step*(i-1))
    print (symbols[i])
	if (i == 1) {
        plot (tab[,1], tab[,2], col='black', ylim=c(0,1), xlab='t', ylab='r', pch = symbols[i], type="b")
    }
	else {
        points (tab[,1], tab[,2], col='black', ylim=c(0,1), xlab='t',ylab='r', pch=symbols[i], t="b")
    }
}
legend ("bottomright", legend=leg, pch=symbols)

dev.off ()
