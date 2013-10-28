#!/usr/bin/python

import subprocess
import os
import sys

def roda (pop, tipo, grau, output, directory, args, sigmax, sigstep):
	params = ""
	sigmax = float(sigmax)
	sigstep = float(sigstep)
	if tipo == 1:
		params = str(grau/2.)
	elif tipo == 2:
		params = str(float(grau)/(pop-1))
	else:
		params = str(grau/2.) + ",0.25"
	cmd = "./genGraph.exe " + str(pop) 
	cmd += " " + str(tipo) + " " + params + " " + output
	subprocess.call(cmd, shell=True)
	
	cmd = "./genplot.sh ../networks/waw"
	subprocess.call(cmd, shell=True)
	cmd = "cp ../networks/waw.pdf ../Resultado/" + directory + "/" + str(tipo)  + "_" + str(grau) +  ".pdf";
	print cmd
	subprocess.call(cmd, shell=True)
	
	sig = 0.0
	x = 0
	while sig <= sigmax:
		cmd = "./kakaroto.exe " + str(sig) + " " + args
		print "Rodando " + cmd 
		f = open ("../Resultado/" + directory + "/" + str(tipo) + "_" + str(grau), 'a')
		f.write(str(sig)+" ")
		f.flush()
		subprocess.call(cmd, shell=True, stdout=f)
		sig += sigstep
		if x%10 == 0:	
			print x
			os.chdir("./dump")
			cmd = "R < plot_freq.r --no-save"
			subprocess.call(cmd, shell=True)
			cmd = "cp frequencias.pdf ../../Resultado/" + directory + "/frequencias" + str(tipo) + "_" + str(grau) + "_" + str(sig) + ".pdf";
			subprocess.call(cmd, shell=True)
			cmd = "R < plot_fase.r --no-save"
			subprocess.call(cmd, shell=True)
			cmd = "cp fases.pdf ../../Resultado/" + directory + "/fases" + str(tipo) + "_" + str(grau) + "_" + str(sig) +  ".pdf";
			subprocess.call(cmd, shell=True)
			os.chdir("..")
		x = x+1

print sys.argv
tipo = 1
subprocess.call ("rm ../Resultado/" + sys.argv[1] + "/*", shell=True);	
subprocess.call ("cp ../Resultado/plot.r ../Resultado/" + sys.argv[1] + "/", shell=True);	
while tipo != 4:
	grau = 6
	while grau <= 18:
		roda (1000, tipo, grau, "../networks/waw.el", sys.argv[1], sys.argv[2],
				sys.argv[3], sys.argv[4])
		grau += 6
	print "================================================================================="
	tipo += 1
