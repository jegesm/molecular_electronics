from string import Template
from numpy import *
from scipy import spatial
import sys, os
PYYAEHMOP_DIR = os.getenv("PYYAEHMOP_DIR")
YAEHMOP_DIR = os.getenv("YAEHMOP_DIR")
sys.path.append(PYYAEHMOP_DIR)
from Yaehmop import *
ASE_DIR = os.getenv("ASE_DIR")
sys.path.append(ASE_DIR)
from ase_ext.io import read
import ase_ext as ase

if len(sys.argv)!=3:
    print("""\n Converts a coordinate file (xyz, pdb, mol2 ...) into YAEHMOP input \n\nUsage: \n\t> python3%s xyz xyzmodelid  \n\n\te.g. python3  %s molecule.xyz 0 \n"""%(sys.argv[0], sys.argv[0]))
    exit()


fn = sys.argv[1]
id = sys.argv[2]
mol = read(fn, id)
ymol = Yaehmop(mol=mol, yaehmop_binary=os.path.join(YAEHMOP_DIR,"bind"), param_file=os.path.join(YAEHMOP_DIR,"eht_param.dat"))

ymol.create_input(just_matrices=False)
ymol.run()
