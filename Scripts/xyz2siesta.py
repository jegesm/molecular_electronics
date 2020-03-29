from string import Template
from numpy import *
from shutil import copyfile
from scipy import spatial
import sys, os

UTILS_DIR = os.getenv("UTILS_DIR")
sys.path.append(UTILS_DIR)
from siesta_utils import *

import ase
from ase.io import read


if len(sys.argv) != 7:
    print("\n Converts a coordinate file (xyz, pdb, mol2 ...) into SIESTA input \n\n USAGE %s Molecule.mol2 inputTemplate franemum psdir  numFrameoutputDir" % (
    sys.argv[0]))
    sys.exit()

fxyz=sys.argv[1]
finp = sys.argv[2]
nm = int(sys.argv[3])
potdir = sys.argv[4]
fdir = sys.argv[5]

try:
    os.mkdir(fdir)
except:
    pass

PT = {}
for i,m in enumerate(Masses):
    PT[Symbols[i]] = [i, m]



mol = read(fxyz, nm)
pos = mol.get_positions()
mins = pos.min(axis=0)
maxs = pos.max(axis=0)
cell = get_cell(pos)


sym = mol.get_chemical_symbols()
species = get_species(sym, PT)
species_dict = {}
species_str = ""
for l in species:
    species_dict[l[1]]=l
    species_str+="{0} {1} {2}\n".format(l[0], l[2], l[1])

#ux = mol.get_cell()[0,0]
#ux = os.getenv("UX")
uz = float(sys.argv[6])
siesta_temp = Template(open(finp, 'r').read())
with open(os.path.join(fdir, "torun.fdf"), "w") as f:
    f.write(siesta_temp.substitute(NAT=len(pos), BASISSIZE="dzp", CELL=cell, NSP=len(set(sym)),
     SPECIES=species_str, UZ=uz))

coords=""
with open(os.path.join(fdir, "coords"), "w") as f:

    for i,s in enumerate(sym):
        c = pos[i]
        coords += "{0} {1} {2} {3}\n".format(c[0],c[1],c[2],species_dict[s][0])
    f.write(coords)


for s in list(set(sym)):
    print(s+".psf")
    copyfile(os.path.join(potdir, s+".psf"), os.path.join(fdir, s+".psf"))
