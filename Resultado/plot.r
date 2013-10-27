x = commandArgs()[3]

name = paste (x, "_grap.pdf", sep='')
pdf (name)

symbols = c(0, 1, 18, 17, 15, 23)

leg = NULL
for (i in 1:10)
{
    idx = 2 + 4*i;
    if (idx > 26) break;
    leg = c (leg, 2+4*i)
	str = paste (x, "_", idx, sep='')
	tab = read.table (str)
    print (2+4*i)
    print (symbols[i])
    print (";-;")
	if (i == 1) {
        plot (tab[,1], tab[,2], col='black', xlab='t', ylab='r', pch = symbols[i], type="b")
    }
	else {
        points (tab[,1], tab[,2], col='black', xlab='t',ylab='r', pch=symbols[i], t="b")
    }
}
print(symbols);
print(leg)
legend ("bottomright", legend=leg, pch=symbols)

dev.off ()
