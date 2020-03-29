import sys, os
from numpy import *
from shutil import copyfile
from string import Template

from ase.io import read
from ase import build, Atoms
from ase.constraints import FixAtoms

UTILS_DIR = os.getenv("UTILS_DIR")
sys.path.append(UTILS_DIR)
from siesta_utils import *

ASE_DIR = os.getenv("ASE_DIR")
sys.path.append(ASE_DIR)
from ase_ext.io import write


def attach_leads(mol):
    leadl = build.fcc100('Au', [3, 3, 6])
    leadr = build.fcc100('Au', [3, 3, 6])
    #num_leadatoms = leadl.get_number_of_atoms()

    dAuN = dAuC = 1.6
    Lcl = leadl[49]
    Lcr = leadr[4]
    Mcl = mol[0]
    Mcr = mol[-4]
    xs = Lcl.x - (Mcl.x + Mcr.x) / 2
    ys = Lcl.y - (Mcl.y + Mcr.y) / 2

    mol.translate([xs, ys, Lcl.z - Mcl.z + dAuN])
    leadr.translate([0.4, -0.4, Mcr.z - Lcr.z + dAuC])

    return leadl + mol + leadr


def extra_with_leads():
    dd = {'WriteCoorInitial': 'T',
'DM.MixingWeight': 0.07,
'WriteMDhistory': ' T',
'MaxSCFIterations': 4,
'DM.NumberPulay': 4,

'SaveHS': 'T',

'MD.TypeOfRun':         'cg',
'MD.NumCGsteps':        0,
#'MD.MaxForceTol         0.1555740000E-01 Ry/Bohr
'WriteCoorXmol':        'T',
'WriteForces':          'T',
'AllocReportLevel':     2,
'WriteMDhistory':       'T',
'WriteMDXmol':          'T',
}
    return dd


def leadmollead(fname="", finp="", potdir=""):

    if len(sys.argv)!=4:
        print("""\n Wtih given molecule coordinates places it into a molecular junction with 2 leads and generates input files for SIESTA  
        \n\nUsage: \n\t> python3 %s N X  \n\n\te.g. python3  %s 18 2\n"""%(sys.argv[0], sys.argv[0]))
        exit()

    if len(sys.argv) == 4:
        fname = sys.argv[1]
        finp = sys.argv[2]
        potdir = sys.argv[3]
    elif fname=="" or  finp=="" or potdir=="":
        print("""\n Wtih given molecule coordinates places it into a molecular junction with 2 leads and generates input files for SIESTA  
                \n\nUsage: \n\t> python3 %s N X  \n\n\te.g. python3  %s 18 2\n""" % (sys.argv[0], sys.argv[0]))
        exit()

    pullmol = read(fname, format='xyz')
    pullmol = attach_leads(pullmol)
    uz = pullmol.get_positions()[:, 2].max() - pullmol.get_positions()[:, 2].min() + 2.04
    pullmol.set_cell([[25, 0, 0], [0, 25, 0], [0, 0, uz]])

    fdir = os.path.join(os.path.dirname(fname), 'siesta')
    try:
        os.mkdir(fdir)
    except:
        pass

    PT = {}
    for i,m in enumerate(Masses):
        PT[Symbols[i]] = [i, m]

    pos = pullmol.get_positions()
    mins = pos.min(axis=0)
    maxs = pos.max(axis=0)
    cell = pullmol.get_cell()


    sym = pullmol.get_chemical_symbols()
    species = get_species(sym, PT)
    species_dict = {}
    species_str = ""
    for l in species:
        species_dict[l[1]]=l
        species_str+="{0} {1} {2}\n".format(l[0], l[2], l[1])

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



if __name__ == "__main__":
    leadmollead()