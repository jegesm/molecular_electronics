from math import pi, cos, sin, sqrt, acos

from ase_ext.atoms import Atoms
from ase_ext.parallel import paropen


def read_xyz(fileobj, index=-1):
    if isinstance(fileobj, str):
        fileobj = open(fileobj)

    lines = fileobj.readlines()
    L1 = lines[0].split()
    if len(L1) == 1:
#        del lines[:2]
        natoms = int(L1[0])
    else:
        natoms = len(lines)
    energy=0
    images = []
    while len(lines) >= natoms:
        positions = []
        symbols = []
        forces = []
        charges = []
        L2 = lines[1].split()
        for line in lines[2:natoms+2]:
            linesplit = line.split()
            if len(linesplit)==5:
                symbol, x, y, z, c = linesplit[:5]
                charges.append(c)
            else:
                symbol, x, y, z = linesplit[:4]
                
            symbols.append(symbol)
            positions.append([float(x), float(y), float(z)])
            if len(linesplit) > 9:
                fx,fy,fz = linesplit[7:10]
                forces.append([float(fx), float(fy), float(fz)])
            
        cell=[]
        if forces==[]:
         forces=None
        if len(lines[1].split())>0:
            if len(lines[1].split("\""))>1:
             L = lines[1].split("\"")[1].split()
            else:
             L = lines[1].split()
             if len(L)>9:
              energy=float(L[9])
            try:
             cell = ((float(L[0]),float(L[1]),float(L[2])),(float(L[3]),float(L[4]),float(L[5])),(float(L[6]),float(L[7]),float(L[8])))
            except ValueError:
             cell = ((0,0,0),(0,0,0),(0,0,0))
            images.append(Atoms(symbols=symbols, positions=positions,forces=forces, energy=energy,cell=cell))
        else:
            images.append(Atoms(symbols=symbols, positions=positions,forces=forces, energy=energy)) 
        if len(charges)>0:
            images[-1].set_charges(charges)
        del lines[:natoms + 2]
        if ((len(lines)>0) and (len(lines[0].split())>0)):
            natoms = int(lines[0].split()[0])
    #return images[index]
    return images

def write_xyz(fileobj, images, symbols=[]):
    if isinstance(fileobj, str):
        fileobj = open(fileobj, 'w')

    if not isinstance(images, (list, tuple)):
        images = [images]

    if symbols==[]:
     symbols = images[0].get_chemical_symbols()
    else: 
     natoms = len(symbols)
     
    for atoms in images:
        fileobj.write('%d\n' % atoms.get_number_of_atoms())
        Strcell=""
        for c in atoms.get_cell()[:]:
            Strcell+=('%f %f %f ' % (c[0],c[1],c[2]))
        fileobj.write(Strcell+"\n")

#        for s, (x, y, z), c in zip(symbols, atoms.get_positions(), atoms.get_charges()):
#            fileobj.write('%-2s %22.15f %22.15f %22.15f %22.15f\n' % (s, x, y, z, c))
        for s, (x, y, z) in zip(symbols, atoms.get_positions()):
            fileobj.write('%-2s %22.15f %22.15f %22.15f \n' % (s, x, y, z))
