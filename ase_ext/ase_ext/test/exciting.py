import os
from ase_ext import Atoms
from ase_ext.io import read, write
from ase_ext.calculators import Exciting
from ase_ext.units import Bohr, Hartree
from ase_ext.test import NotAvailable

try:
    import lxml
except ImportError:
    raise NotAvailable('This test need lxml module.')

a = Atoms('N3O',
          [(0, 0, 0), (1, 0, 0), (0, 0, 1), (0.5, 0.5, 0.5)],
          pbc=True)

raise NotAvailable('Problem with lxml module.')

write('geo.exi', a)
b = read('geo.exi')

print(a)
print(a.get_positions())
print(b)
print(b.get_positions())

calculator = Exciting(dir='excitingtestfiles',
                      kpts=(4, 4, 3),
                      maxscl=3,
                      #bin='/fshome/chm/git/exciting/bin/excitingser'
                      )
