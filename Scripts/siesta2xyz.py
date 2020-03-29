from string import Template
import sys, os
from numpy import *

UTILS_DIR = os.getenv("UTILS_DIR")
sys.path.append(UTILS_DIR)
from siesta_utils import *

ASE_DIR = os.getenv("ASE_DIR")
sys.path.append(ASE_DIR)
from ase_ext.io import read
import ase_ext as ase


if len(sys.argv)!=3:
    print("""\n Converts SIESTA STRUCT_OUT into xyz and xsf file with ase (cell information included)\n\nUsage: \n\t> python3 %s struct_out outfile  \n\n\te.g. python3  %s siesta.STRUCT_OUT siesta.out\n"""%(sys.argv[0], sys.argv[0]))
    exit()

fstr = sys.argv[1]
fout = sys.argv[2]

coo, sym, cell = read_struct(fstr, fout)
strsyms = ""
for s in sym:
    strsyms += s

mol = ase.Atoms(strsyms, coo)
mol.set_cell(cell)
mol.write("temp.xyz")
mol.write("temp.xsf")
