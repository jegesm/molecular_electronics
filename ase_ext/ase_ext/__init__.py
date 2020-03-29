# Copyright 2008, 2009 CAMd
# (see accompanying license files for details).

"""Atomic Simulation Environment."""


from ase_ext.atom import Atom
from ase_ext.atoms import Atoms

_deprecate_things_from_ase_module = True

# Some day in the future, we will uncomment this line:
#__all__ = ['Atoms', 'Atom']  

from ase_ext.units import *
from ase_ext.io import read, write
from ase_ext.io.trajectory import PickleTrajectory
from ase_ext.dft import STM, monkhorst_pack, DOS
from ase_ext.optimize.mdmin import MDMin
from ase_ext.optimize.lbfgs import HessLBFGS
from ase_ext.optimize.fire import FIRE
from ase_ext.optimize.lbfgs import LBFGS, LBFGSLineSearch
from ase_ext.optimize.bfgs import BFGS
from ase_ext.optimize import QuasiNewton
from ase_ext.md.verlet import VelocityVerlet
from ase_ext.md.langevin import Langevin
from ase_ext.constraints import *
from ase_ext.calculators.lj import LennardJones
from ase_ext.calculators.emt import EMT
from ase_ext.calculators.siesta import Siesta
from ase_ext.calculators.dacapo import Dacapo
from ase_ext.calculators.vasp import Vasp
from ase_ext.calculators.aims import Aims, AimsCube
from ase_ext.calculators.turbomole import Turbomole
from ase_ext.calculators.dftb import Dftb
from ase_ext.neb import NEB, SingleCalculatorNEB
from ase_ext.visualize import view
from ase_ext.data import chemical_symbols, atomic_numbers, atomic_names, \
     atomic_masses, covalent_radii, reference_states
from ase_ext.data.molecules import molecule

from math import sqrt, pi
import numpy as np



if _deprecate_things_from_ase_module:
    import types
    from ase_ext.utils import deprecate as dep

    _locals = locals()

    for name, obj in list(_locals.items()):
        if name.startswith('_') or name in ['Atoms', 'Atom']:
            continue

        if isinstance(obj, float):
            if name == 'pi':
                pi = dep.DeprecatedFloat(pi, 'pi', 'math')
            else:
                _locals[name] = dep.DeprecatedFloat(obj, name, 'ase.units')
        elif isinstance(obj, types.ModuleType):
            pass  # how about np?  XXX
        elif hasattr(obj, '__module__'):
            module = obj.__module__
            if module.startswith('ase.optimize'):
                module = 'ase.optimize'
            elif module.startswith('ase.md'):
                module = 'ase.md'
            elif name == 'PickleTrajectory':
                module = 'ase.io'
            _locals[name] = dep.Deprecate(obj, name, module)
        else:
            pass  # how about atomic_numbers, covalent_radii, ... ? XXX

    np = dep.DeprecatedNumpyImport()

