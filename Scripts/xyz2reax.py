from string import Template
from numpy import *
from scipy import spatial
import sys, os
ASE_DIR = os.getenv("ASE_DIR")
sys.path.append(ASE_DIR)
from ase_ext.io import read


reax_data_temp = Template("""# Position data for reax force field

 $NAT atoms

 $NATT atom types

 $XLO $XHI xlo xhi
 $YLO $YHI ylo yhi
 $ZLO $ZHI zlo zhi

 Masses

 $MASSES

 Atoms

 $COORDS
 """)

def get_masses(syms):
    strm = ""
    chemical_symbols = ""
    for i, s in enumerate(syms):
        chemical_symbols += s.symbol + " "
        strm += "\t%d\t%f\t#%s\n"%(i+1, s.get('mass'), s.symbol)
    return strm, chemical_symbols


def get_coords(mol, indatoms):
    strc = ""
    for a in mol:
        symindex = where(indatoms==a.number)[0]

        strc += "%d %d 0 %f %f %f\n"%(a.index+1, symindex+1, a.x, a.y, a.z)
    return strc


if len(sys.argv)!=5:
    print("""\n Converts a coordinate file (xyz, pdb, mol2 ...) into lammps reax input \n\nUsage: \n\t> python3%s xyz reaxtemplateFile xyzmodelid outputdir \n\n\te.g. python3  %s molecule.xyz in.reax.temp 0 Reax-relax\n"""%(sys.argv[0], sys.argv[0]))
    exit()

finxyz = sys.argv[1]
finp = sys.argv[2]
nm = int(sys.argv[3])
fdir = sys.argv[4]

mol = read(finxyz, nm)

pos = mol.get_positions()
mins = pos.min(axis=0)
maxs = pos.max(axis=0)
#cell = get_cell(pos)

#maxs = (maxs - mins)/2

symlist = set(mol.get_chemical_symbols())
symlist = sorted(symlist)
syms = ase.Atoms(symlist)
indatoms = syms.get_atomic_numbers()
masses, chemical_symbols = get_masses(syms)
with open(os.path.join(fdir, "data"), "w") as f:
    f.write(reax_data_temp.substitute(NAT=mol.get_number_of_atoms(),\
                 NATT=len(set(mol.get_atomic_numbers())),\
                 XLO=mins[0]-5, XHI=maxs[0]+5,\
                 YLO=mins[1]-5, YHI=maxs[1]+5,\
                 ZLO=mins[2]-5, ZHI=maxs[2]+5,\
                 MASSES=masses,\
                 COORDS=get_coords(mol, indatoms)))

reax_temp = Template(open(finp, 'r').read())
with open(os.path.join(fdir, "in.reax"), "w") as f:
    f.write(reax_temp.substitute(SYM=chemical_symbols))

