import sys, os
from numpy import *

from ase.io import read, write
from ase import build, Atoms
from ase.constraints import FixAtoms

UTILS_DIR = os.getenv("UTILS_DIR")
sys.path.append(UTILS_DIR)
from Stretch_mol_with_Siesta import Pull_Mol


ASE_DIR = os.getenv("ASE_DIR")
sys.path.append(ASE_DIR)
#from ase_ext.io import write
#import ase_ext as ase


class Pull_chain_with_leads(Pull_Mol):

    def __init__(self, mol, label='temp'):
        super().__init__(mol)
        self.num_leadatoms = 0
        self.label = label

    def attach_leads(self):
        leadl = build.fcc100('Au', [3, 3, 6])
        leadr = build.fcc100('Au', [3, 3, 6])
        self.num_leadatoms = leadl.get_number_of_atoms()

        dAuC = 1.6
        Lcl = leadl[49]
        Lcr = leadr[4]
        Mcl = self.mol[0]
        Mcr = self.mol[-1]
        xs = Lcl.x - (Mcl.x + Mcr.x) / 2
        ys = Lcl.y - (Mcl.y + Mcr.y) / 2

        self.mol.translate([xs, ys, Lcl.z - Mcl.z + dAuC])
        leadr.translate([0.4, -0.4, Mcr.z - Lcr.z + dAuC])

        self.mol = leadl + self.mol + leadr

    def update_mol(self):
        new_mol = read('.%s/%s.xyz'%(self.label, self.label))
        new_pos = new_mol.get_positions()
        self.mol.set_positions(new_pos)


def extra_for_leads():
    dd = {'WriteCoorInitial': 'T',
'DM.MixingWeight': 0.03,
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

def extra_for_opt():
    dd = {'WriteCoorInitial': 'T',
'DM.MixingWeight': 0.1,
'WriteMDhistory': ' T',
'MaxSCFIterations': 100,
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


def create_chain(numl):
    dC = 1.1
    pos = zeros((numl,3))
    pos[:,2] = linspace(0.0, dC*(numl-1), numl)
    return Atoms("C"*numl, positions=pos)

if len(sys.argv)<3:
    print("""\n Generates and N long carbon chain it into a molecular junction with 2 leads\n\nUsage: \n\t> python3 %s N  \n\n\te.g. python3  %s 18 \n"""%(sys.argv[0], sys.argv[0]))
    exit()

numl = int(sys.argv[1])
numstep = int(sys.argv[2])


# Generate the molecule
xyz = create_chain(numl)
#ujxyz = Atoms(xyz.get_chemical_symbols(), positions=xyz.get_positions())

# molecule is generated

# Start the stretching
pullmol = Pull_Mol(xyz)
for step in range(numstep):

    c = FixAtoms(indices=[0, -1])
    pullmol.label = 'psnap-%02d'%step
    pullmol.mol.set_constraint(c)
    pullmol.add_extra_siesta_arguments(extra_for_opt())
    pullmol.optimize()
    #pullmol.mol.write("psnap-%02d.xyz" % step)

    whole_setup = Pull_chain_with_leads(pullmol.mol, label="step-%02d"%step)
    whole_setup.attach_leads()
    print(step)
    uz = whole_setup.mol.get_positions()[:, 2].max() - whole_setup.mol.get_positions()[:, 2].min() + 2.04
    whole_setup.mol.set_cell([[25, 0, 0], [0, 25, 0], [0, 0, uz]])
    whole_setup.add_extra_siesta_arguments(extra_for_leads())
    whole_setup.iter_md()
    pullmol.mol.set_constraint()
    pullmol.stretch_one()
