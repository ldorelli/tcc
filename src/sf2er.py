import subprocess
import sys

if len(sys.argv) < 1:
	print 'Usage: sf2er.py + <1.dirname>'
else:
	subprocess.call("rm ../Resultado/" + sys.argv[1] + "/*", shell=True)
	
	for i in range(0, 101):
		alpha = i*0.01
		subprocess.call("./genGraph.exe 100 5 3," + str(alpha) + " \"../networks/sf2er.el\"", shell=True)
		sig = 0.0
		while sig <= 1:
			cmd = "./kakaroto.exe " + str(sig) + " 0.02 0 "
			if len(sys.argv) > 2:
				cmd += " " + sys.argv[2]
			print cmd
			f = open("../Resultado/" + sys.argv[1]+"/sf2er"+str(alpha), 'a')
			f.write (str(sig)+" ")
			f.flush()
			subprocess.call(cmd, shell=True, stdout=f)
			sig += 0.05
		f.close()
	print '\a'
