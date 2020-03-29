import sys, os
from numpy import *

from ase.io import read
from ase import build, Atoms
from ase.calculators.dftb import Dftb

ASE_DIR = os.getenv("ASE_DIR")
sys.path.append(ASE_DIR)
from ase_ext.io import write

from ase.calculators.siesta import Siesta
from ase.units import Ry, eV


def create_chain(numl):
    dC = 1.44
    pos = zeros((numl, 3))
    pos[:, 2] = linspace(0.0, dC * numl, numl)
    return Atoms("C" * numl, positions=pos)

class Pull_Mol():

    def __init__(self, mol, label="temp"):
        self.mol = mol.copy()
        self.label = 'temp'
        self.extra_fdf = {}

    def add_extra_siesta_arguments(self, filen):
        if type(filen) == 'str':
            Data = open(filen, 'r').readlines()
            for d in Data:
                if len(d.split())>0:
                    key = d.split()[0]
                    value = d[len(key):-1]
                    self.extra_fdf[key] = value
        else:
            self.extra_fdf = filen


    def iter_md(self):
        calc = Siesta(label = self.label,
                   xc='CA',
                   mesh_cutoff=400 * eV,
                   energy_shift=0.03 * eV,
                   basis_set='SZ',
                   fdf_arguments=self.extra_fdf,
                   )


        tempdir = "."+self.label
        try:
            os.mkdir(tempdir)
        except:
            pass

        os.chdir(tempdir)
        self.mol.set_calculator(calc)
        E = self.mol.get_total_energy()
        os.chdir('../')

    def optimize(self, calculator='siesta'):
        from ase.optimize import QuasiNewton

        if calculator=='siesta':
            calc = Siesta(label = self.label,
                   xc='CA',
                   mesh_cutoff=400 * eV,
                   energy_shift=0.03 * eV,
                   basis_set='SZ',
                   fdf_arguments=self.extra_fdf,
                   )
        elif calculator=='dftb':
            extras = {}
            for S in set(self.mol.get_chemical_symbols()):
                key = 'Hamiltonian_MaxAngularMomentum_'+S
                if S=='H':
                    value = '"s"'
                else:
                    value = '"p"'
                extras[key]=value
            calc = Dftb(label=self.label,
                 atoms=self.mol, **extras)

        tempdir = "."+self.label
        try:
            os.mkdir(tempdir)
        except:
            pass

        os.chdir(tempdir)
        self.mol.set_calculator(calc)
        dyn = QuasiNewton(self.mol, trajectory=self.label+'.traj')
        dyn.run(fmax=0.5)
        os.chdir('../')

    def update_mol(self):
        self.mol = read('.%s/%s.xyz'%(self.label, self.label))

    def stretch_one(self, ds=1.01):
        pos = self.mol.get_positions()
        pos[:, 2] = pos[:, 2] * ds
        self.mol.set_positions(pos)




"""

An example script

xyz = Pull_Mol(create_chain(6))
xyz.mol.set_cell([[10, 0, 0], [0, 10, 0], [0, 0, 20]])

for r in range(5):
    xyz.iter_md()
    xyz.stretch_one()
    xyz.mol.write("snap-%02d.xyz"%r)
"""
