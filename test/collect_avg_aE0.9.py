from glob import glob

files = glob('aE0.9_d00.0404_q?.*_v1e-19/averages.pco')

data = {}
for filename in files:
    f = open(filename)
    l = f.readlines()
    f.close()

    s = filename.split('/')
    q = s[0]
    q = q[16:]
    q = float(q[:q.find('_')])

    nAs = float(l[0].split('\t')[1])
    nBs = float(l[1].split('\t')[1])
    nCs = float(l[2].split('\t')[1])
    nDs = float(l[3].split('\t')[1])
    nEs = float(l[4].split('\t')[1])
    nFs = float(l[5].split('\t')[1])
    nZs = float(l[6].split('\t')[1])
    ncells = float(l[7].split('\t')[1])
    ntot = nAs+nBs+nCs+nDs+nEs+nFs
    fracA = nAs/ntot
    fracB = nBs/ntot
    fracC = nCs/ntot
    fracD = nDs/ntot
    fracE = nEs/ntot
    fracF = nFs/ntot

    filename2 = filename.replace('averages','output')
    f = open(filename2)
    nrrelives = int(f.read().split()[-1])
    f.close()

    data[q] = [q,fracE, fracB, ntot, ncells, nrrelives, nZs]


keys = data.keys()
keys.sort()
f = open("overview_d00.0404_aE0.9.txt","w")
for k in keys:
    f.write("%f\t%f\t%f\t%f\t%f\t%d\t%f\n" % tuple(data[k]))
f.close()

