import subprocess
import sys

if len(sys.argv) < 1:
	print 'Usage: rodar.py + <1.filename>'
else:
	subprocess.call("rm ../Resultado/" + sys.argv[1], shell=True)

	sig = 0.0
	while sig <= 1:
		cmd = "./kakaroto.exe " + str(sig) + " 0.02"
		if len(sys.argv) > 2:
			cmd += " " + sys.argv[2]
		print cmd
		f = open("../Resultado/" + sys.argv[1], 'a')
		f.write (str(sig)+" ")
		f.flush()
		subprocess.call(cmd, shell=True, stdout=f)
		sig += 0.01
	print '\a'
