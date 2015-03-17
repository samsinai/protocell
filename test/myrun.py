from numpy import *
import os

for q in list(arange(0.0,1.001,0.1))+list(arange(0.91,0.991,0.01)):
        dirname = 'aE0.9_d00.0404_q%.3f_v1e-19' % (q)
        os.system('mkdir ' + dirname)
        os.system("cp run.pbs IT0.pcs " + dirname)
        os.chdir(dirname)
        f=open("input.pci", "w")
        f.write("""it 0
        max_it 100000000
        crep 10000000
        erep 1000
        seed 0
        log 0
        timebrep 1

        avgrep 10
        avgstart 50000000

        revive 4
        modelnr 2

	intcomp 1

	aA 0
	aC 0
	aD 0
	aE 0.9
	aF 0
        k_s 1
        q %.3f
        k_phi 0.04
        k_theta 0
        k_z 0.01
        z0 1.0
        k_d0 0.0404
        k_d1 0.0
        m 1000000000
        vol 1e-19
""" % (q))
        f.close()
        os.system("qsub run.pbs")
        os.chdir("..")
        
