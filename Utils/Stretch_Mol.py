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


class Pull_Mol(Atoms):
    label = "temp"

    def set_label(self, label):
        self.label = label

    def stretch_one(self, ds=1.01):
        pos = self.get_positions()
        pos[:, 2] = pos[:, 2] * ds
        self.set_positions(pos)

    def optimize(self):
        from ase.optimize import QuasiNewton
        self.add_calculator()
        tempdir = "."+self.label
        try:
            os.mkdir(tempdir)
        except:
            pass

        os.chdir(tempdir)
        dyn = QuasiNewton(self, trajectory=self.label+'.traj')
        dyn.run(fmax=0.5)
        os.chdir('../')

    def iter_md(self):
        tempdir = "."+self.label
        self.add_calculator()
        try:
            os.mkdir(tempdir)
        except:
            pass

        os.chdir(tempdir)
        E = self.mol.get_total_energy()
        os.chdir('../')

    def add_calculator(self):
        return None

class Pull_Mol_Siesta(Pull_Mol):
    extra_fdf = {}

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


    def add_calculator(self):
        self.set_calculator(Siesta(label=self.label,
                      xc='CA',
                      mesh_cutoff=400 * eV,
                      energy_shift=0.03 * eV,
                      basis_set='SZ',
                      fdf_arguments=self.extra_fdf,
                      ))



# class Pull_Mol_DFTB2(Pull_Mol):
#
#     def set_calculator(self):
#         extras = {}
#         for S in set(self.mol.get_chemical_symbols()):
#             key = 'Hamiltonian_MaxAngularMomentum_'+S
#             if S=='H':
#                 value = '"s"'
#             else:
#                 value = '"p"'
#             extras[key]=value
#
#         self.calculator = Dftb(label=self.label, atoms=self.mol, **extras)


class Pull_Mol_DFTB(Pull_Mol):

    def add_calculator(self):
        extras = {}
        for S in set(self.get_chemical_symbols()):
            key = 'Hamiltonian_MaxAngularMomentum_'+S
            if S=='H':
                value = '"s"'
            else:
                value = '"p"'
            extras[key]=value

        self.set_calculator(Dftb(label=self.label, atoms=self, **extras))



"""

An example script

xyz = Pull_Mol(create_chain(6))
xyz.mol.set_cell([[10, 0, 0], [0, 10, 0], [0, 0, 20]])

for r in range(5):
    xyz.iter_md()
    xyz.stretch_one()
    xyz.mol.write("snap-%02d.xyz"%r)
"""
