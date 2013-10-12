import subprocess

def roda (pop, tipo, grau, output):
	tipo = str(tipo)
	params = ""
	if tipo == 1:
		params = str(grau/2)
	elif tipo == 2:
		params = str(grau/(pop-1))
	else:
		params = str(grau/2) + ",0.25"
	cmd = "./genGraph.exe " + str(pop) 
	cmd += " " + tipo + " " + params + " " + output
	subprocess.call(cmd, shell=True)
	sig = 0.0
	while sig <= 1:
		cmd = "./kakaroto.exe " + str(sig) + " 0.02"
		f = open ("../Resultado/" + tipo + "_" + str(grau), 'a')
		f.write(str(sig)+" ")
		f.flush()
		subprocess.call(cmd, shell=True, stdout=f)
		sig += 0.01

tipo = 1
while tipo != 4:
	grau = 6
	while grau != 30:
		roda (100, tipo, grau, "../networks/waw.el")
		grau += 2
	tipo += 1
