from ase_ext import Atoms
from ase_ext.calculators.emt import EMT

h2 = Atoms('H2', positions=[(0, 0, 0), (0, 0, 1.1)],
           calculator=EMT())
print(h2.calc.get_numeric_forces(h2))
print(h2.get_forces())

