import sys
import os
from numpy import *

if len(sys.argv)!=2:
    print("""\n Converts SIESTA EIG into python numpy npy \n\nUsage: \n\t> python3 %s Systemlabel.EIG   \n\n\te.g. python3  %s h2o-relaxation.EIG\n"""%(sys.argv[0], sys.argv[0]))
    exit()

fn = sys.argv[1]
d = open(fn, 'r')

EF = float(d.readline())
print(EF)
norb = int(d.readline().split()[0])
print(norb)

eigs = array([float(f) for f in d.readline().split()[1:]])
i=1
for l in d.readlines():
    teig = array([float(f) for f in l.split()])
    eigs = concatenate((eigs, teig), 0)

eigsep = zeros_like(eigs)

eigsep[0] = nan
eigsep[1:] = [eigs[j]-eigs[j-1] for j in range(1, eigsep.shape[0])] 

newfn = os.path.join(os.path.dirname(fn), "eig")
data = zeros((eigs.shape[0], 2))
data[:,0] = eigs
data[:,1] = eigsep

occ = eigs[eigs<EF]
print(eigsep[occ.shape[0]])
savez(newfn, E=eigs, dE=eigsep)
#savetxt(newfn, data)
savetxt(newfn, eigs-EF)
