"""Interfaces to different ASE compatible force-calculators."""

import numpy as np

from ase_ext import _deprecate_things_from_ase_module
from ase_ext.calculators.lj import LennardJones
from ase_ext.calculators.emt import EMT
from ase_ext.calculators.siesta import Siesta
from ase_ext.calculators.dacapo import Dacapo
from ase_ext.calculators.vasp import Vasp
from ase_ext.calculators.aims import Aims, AimsCube
from ase_ext.calculators.turbomole import Turbomole
from ase_ext.calculators.exciting import Exciting
from ase_ext.calculators.dftb import Dftb
from ase_ext.calculators.singlepoint import SinglePointCalculator
from ase_ext.calculators.test import numeric_force, numeric_forces, TestPotential

if _deprecate_things_from_ase_module:
    from ase_ext.utils.deprecate import Deprecate
    _locals = locals()
    for name in ['LennardJones', 'EMT', 'Siesta', 'Dacapo', 'Vasp',
                 'Aims', 'AimsCube', 'Turbomole', 'Exciting', 'Dftb',
                 'SinglePointCalculator', 'numeric_force', 'numeric_forces',
                 'TestPotential']:
        obj = _locals[name]
        _locals[name] = Deprecate(obj, name, obj.__module__, 'ase.calculators')
