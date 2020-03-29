from math import pi, cos, sin, sqrt, acos

from ase_ext.atoms import Atoms
from ase_ext.parallel import paropen


def read_xyz(fileobj, index=-1):
    if isinstance(fileobj, str):
        fileobj = open(fileobj)

    Info={}

    lines = fileobj.readlines()
    L1 = lines[0].split()
    if len(L1) == 1:
#        del lines[:2]
        natoms = int(L1[0])
    else:
        natoms = len(lines)
    images = []
    while len(lines) >= natoms:
        positions = []
        symbols = []
        forces = []
        L2 = lines[1].split()
        for line in lines[2:natoms+2]:
            linesplit = line.split()
            symbol, x, y, z = linesplit[:4]
            symbols.append(symbol)
            positions.append([float(x), float(y), float(z)])
            if len(linesplit) > 9:
                fx,fy,fz = linesplit[7:10]
                forces.append([float(fx), float(fy), float(fz)])
            
        cell=[]
        LC=0
        Info['energy']=None
        if forces==[]:
         forces=None
        if len(lines[1].split())>0:
            L = lines[1].split()
            if L[0].split("=")[0]=="Lattice":
             LC=0
            else:
             Info['energy']=float(L[0].split("=")[1])
             if L[1].split("=")[0]=="Lattice":
              LC=1
             else:
              Info['time']=float(L[1].split("=")[1])
              if L[2].split("=")[0]=="Lattice":
               LC=2
              else:
               Info['i_step']=int(L[2].split("=")[1])
               if L[3].split("=")[0]=="Lattice":
                LC=3
            cell = ((float(L[LC].split("\"")[1]),float(L[LC+1]),float(L[LC+2])),(float(L[LC+3]),float(L[LC+4]),float(L[LC+5])),(float(L[LC+6]),float(L[LC+7]),float(L[LC+8].split("\"")[0])))

            images.append(Atoms(symbols=symbols, positions=positions,forces=forces, energy=Info['energy'],cell=cell,Info=Info))
        else:
            images.append(Atoms(symbols=symbols, positions=positions,forces=forces)) 
        del lines[:natoms + 2]
        if len(lines)>0:
            natoms = int(lines[0].split()[0])
    #return images[index]
    return images

def write_xyz(fileobj, images):
    if isinstance(fileobj, str):
        fileobj = paropen(fileobj, 'w')

    if not isinstance(images, (list, tuple)):
        images = [images]

    symbols = images[0].get_chemical_symbols()
    natoms = len(symbols)
    for atoms in images:
        fileobj.write('%d\n' % natoms)
        Strcell="Lattice=\""
        for c in atoms.get_cell()[:]:
            Strcell+=('%f %f %f ' % (c[0],c[1],c[2]))
        fileobj.write(Strcell[:-1]+"\" Properties=species:S:1:pos:R:3\n")

        for s, (x, y, z) in zip(symbols, atoms.get_positions()):
            fileobj.write('%-2s %22.15f %22.15f %22.15f\n' % (s, x, y, z))
