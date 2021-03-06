from math import pi, cos, sin, sqrt, acos
import numpy as np

from ase_ext.atom import Atom
from ase_ext.atoms import Atoms
from ase_ext.parallel import paropen
from ase_ext.units import Bohr


def read_I_info(fileobj, index=-1):
    if isinstance(fileobj, str):
        fileobj = open(fileobj)

    lines = fileobj.readlines()

    del lines[0]

    finished = False
    s = Atoms()
    while not finished:
        w = lines.pop(0).split()
        if w[0].startswith('"'):
            position = Bohr * np.array([float(w[3]), float(w[4]), float(w[5])])
            s.append(Atom(w[0].replace('"',''), position))
        else:
            finished = True

    return s
