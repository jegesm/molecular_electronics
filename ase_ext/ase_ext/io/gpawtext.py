from ase_ext.atoms import Atom, Atoms
from ase_ext.calculators.singlepoint import SinglePointCalculator


def read_gpaw_text(fileobj, index=-1):
    if isinstance(fileobj, str):
        fileobj = open(fileobj)

    def index_startswith(lines, string):
        for i, line in enumerate(lines):
            if line.startswith(string):
                return i
        raise ValueError

    lines = fileobj.readlines()
    images = []
    while True:
        try:
            i = lines.index('Unit Cell:\n')
        except ValueError:
            pass
        else:
            cell = []
            pbc = []
            for line in lines[i + 3:i + 6]:
                words = line.split()
                if len(words) == 5:  # old format
                    cell.append(float(words[2]))
                    pbc.append(words[1] == 'yes')
                else:                # new format with GUC
                    cell.append([float(word) for word in words[3:6]])
                    pbc.append(words[2] == 'yes')
            
        try:
            i = lines.index('Positions:\n')
        except ValueError:
            break

        atoms = Atoms(cell=cell, pbc=pbc)
        for line in lines[i + 1:]:
            words = line.split()
            if len(words) != 5:
                break
            n, symbol, x, y, z = words
            symbol = symbol.split('.')[0]
            atoms.append(Atom(symbol, [float(x), float(y), float(z)]))
        lines = lines[i + 5:]
        try:
            i = lines.index('-------------------------\n')
        except ValueError:
            e = None
        else:
            line = lines[i + 9]
            assert line.startswith('Zero Kelvin:')
            e = float(line.split()[-1])
        try:
            ii = index_startswith(lines, 'Total Charge:')
        except ValueError:
            q = None
        else:
            q = float(lines[ii].split()[2])
        try:
            ii = lines.index('Forces in eV/Ang:\n')
        except ValueError:
            f = None
        else:
            f = []
            for i in range(ii + 1, ii + 1 + len(atoms)):
                try:
                    x, y, z = lines[i].split()[-3:]
                    f.append((float(x), float(y), float(z)))
                except (ValueError, IndexError) as m:
                    raise IOError('Malformed GPAW log file: %s' % m)

        if len(images) > 0 and e is None:
            break

        if e is not None or f is not None:
            atoms.set_calculator(SinglePointCalculator(e, f, None, None, atoms)) ### Fixme magmoms
        if q is not None:
            n = len(atoms)
            atoms.set_charges([q / n] * n)

        images.append(atoms)
        lines = lines[i:]

    if len(images) == 0:
        raise IOError('Corrupted GPAW-text file!')
    
    return images[index]
