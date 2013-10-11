import subprocess
import sys

subprocess.call("rm ../Resultado/res.r", shell=True)

sig = 0.0
while sig <= 1:
	cmd = "./kakaroto.exe " + str(sig) + " 0.02"
	if len(sys.argv) > 1:
		cmd += " " + sys.argv[1]
	print cmd
	f = open("../Resultado/res.r", 'a')
	f.write (str(sig)+" ")
	f.flush()
	subprocess.call(cmd, shell=True, stdout=f)
	sig += 0.01
print '\a'
