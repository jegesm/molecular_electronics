import sys, os
from numpy import *

from ase.io import read
from ase import build, Atoms
from ase.constraints import FixAtoms

UTILS_DIR = os.getenv("UTILS_DIR")
sys.path.append(UTILS_DIR)
from Glycine import gen_randglycine, gen_mol_from_fasta
from Stretch_Mol import Pull_Mol_DFTB


ASE_DIR = os.getenv("ASE_DIR")
sys.path.append(ASE_DIR)
from ase_ext.io import write


def calculator_NVT():
    return  Dftb(
    label='h2o',
    atoms=test,
    run_manyDftb_steps=True,
    Hamiltonian_MaxAngularMomentum_='',
    Hamiltonian_MaxAngularMomentum_O='"p"',
    Hamiltonian_MaxAngularMomentum_H='"s"',
    Driver_='VelocityVerlet',
    Driver_MDRestartFrequency=5,
    Driver_Velocities_='',
    Driver_Velocities_empty='<<+ "velocities.txt"',
    Driver_Steps=500,
    Driver_KeepStationary='Yes',
    Driver_TimeStep=8.26,
    Driver_Thermostat_='Berendsen',
    Driver_Thermostat_Temperature=0.00339845142,  # 800 deg Celcius
    # Driver_Thermostat_Temperature=0.0, # 0 deg Kelvin
    Driver_Thermostat_CouplingStrength=0.01)


if len(sys.argv)==1 or len(sys.argv)>4:
    print("""\n Generates and N long polyglycine of which has X random amino acid replacements. Stretches it and optimizes the geometry (with DFTB) 
    \n\nUsage: \n\t> python3 %s N X  \n\n\te.g. python3  %s 18 2\n"""%(sys.argv[0], sys.argv[0]))
    exit()


def relax_mol(numl=9, numr=0, numstep=70):
    # Generate the molecule
    if len(sys.argv) == 3 and type(sys.argv[1]) == type('s'):
        xyz, mol2, fasta = gen_mol_from_fasta(sys.argv[1])
        numstep = int(sys.argv[2])
    elif len(sys.argv) == 4:
        numl = int(sys.argv[1])
        numr = int(sys.argv[2])
        numstep = int(sys.argv[3])
        xyz, mol2, fasta = gen_randglycine(numl, numr)
    else:
        xyz, mol2, fasta = gen_randglycine(numl, numr)

    try:
        os.mkdir(fasta)
    except:
        pass

    os.chdir(fasta)
    # Only modified ase (ase_ext) can print in mol2 and use bonded information for initial optimization with Dreiding
    # Need to add this to next to the standard ase library
    mol2.write("%s.mol2"%fasta)

    xyz.pop(4)
    xyz.rotate([0, -pi/180*90, 0])
    #ujxyz = Atoms(xyz.get_chemical_symbols(), positions=xyz.get_positions())

    # molecule is generated


    # Start the stretching
    #pullmol = Pull_Mol_DFTB(ujxyz)
    pullmol = Pull_Mol_DFTB(xyz.get_chemical_symbols(), positions=xyz.get_positions())
    for step in range(numstep):

        c = FixAtoms(indices=[0, -4])
        pullmol.label = 'opt-%s-%03d'%(fasta, step)
        pullmol.set_constraint(c)
        pullmol.optimize()
        pullmol.set_constraint()
        pullmol.stretch_one()


if __name__ == "__main__":
    relax_mol()