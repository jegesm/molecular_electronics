from numpy import *
from scipy.io import FortranFile

import sys
import os

def rd(i=1):
    return fromfile(f, dtype=int32, count=i)


if len(sys.argv)!=2:
    print("""\n Converts SIESTA WFSX into python numpy npz \n\nUsage: \n\t> python3 %s Systemlabel.WFSX   \n\n\te.g. python3  %s h2o-relaxation.WFSX\n"""%(sys.argv[0], sys.argv[0]))
    exit()

fn = sys.argv[1]
with open(fn,'rb') as f:
    rd(1)
    nk = fromfile(f, dtype=int32, count=1)[0];
    gamma = fromfile(f, dtype=int32, count=1)[0]; rd(2)
    nspin = fromfile(f, dtype=int32, count=1)[0]; rd(2)
    nuotot = fromfile(f, dtype=int32, count=1)[0]; rd(1)
    #print('#nk', nk, gamma, nspin, nuotot)

    iaorb = zeros(nuotot, dtype=int32)
    labelfis = zeros(nuotot, dtype=str)
    iphorb = zeros(nuotot, dtype=int32)
    cnfigio = zeros(nuotot, dtype=int32)
    symfio = zeros(nuotot, dtype=str)

    thress = 0.00000001
    rd(1)
    for i in range(nuotot):
        iaorb[i] = fromfile(f, dtype=int32, count=1)[0]
        labelfis[i] = f.read(20)
        iphorb[i] = fromfile(f, dtype=int32, count=1)[0]
        cnfigio[i] = fromfile(f, dtype=int32, count=1)[0]
        symfio[i] = f.read(20)
    rd(1)

    rd(1)
    for iik in range(nk):
       for iispin in range(nspin):
            ik = fromfile(f, dtype=int32, count=1)[0]
            k = fromfile(f, dtype=float64, count=3)
            fromfile(f, dtype=float64, count=1)
            rd(2)

            #assert ik==iik+1, 'error in index of k-point'
            ispin = fromfile(f, dtype=int32, count=1)[0]
            rd(2)
            #assert ispin==iispin+1, 'error in index of spin'
            nwflist = fromfile(f, dtype=int32, count=1)[0]
            #print("#ik ", ik, k, "ispin ", ispin, "nfw ", nwflist)
            rd(1)


    # Loop over wavefunctions
            psi = zeros((nuotot, nwflist))
            energy = zeros(nwflist)
            for iw in range(nwflist):
                rd(1)
                indwf = fromfile(f, dtype=int32, count=1)[0]
                rd(2)
                energy[iw] = fromfile(f, dtype=float64, count=1)[0]
                #print("#energy ",  energy[iw], "indwf ", indwf)
                #print("Atom  Species Orb-global  Orb-in-atom Orb-type      Re(psi)   Im(psi)")
                rd(2)
                psi[:, iw] = fromfile(f, dtype=float32, count=nuotot)
                #print(dot(psi,psi))

                #print(psi[abs(psi*psi)>thress][:5])
                #print(psi[:5])
                rd(1)


newfn = os.path.join(os.path.dirname(fn), "wfs")
savez(newfn, E=energy, Psi=psi, iaorb=iaorb, iphorb=iphorb, labelfis=labelfis, cnfigio=cnfigio, symfio=symfio)


