from math import pi, cos, sin, sqrt, acos
import re
from numpy import array
from ase_ext.atoms import Atoms
from ase_ext.parallel import paropen


def read_mol2(fileobj, index=-1):
    subnumber = re.compile('[0-9]*')
    if isinstance(fileobj, str):
        fileobj = open(fileobj)
    
    lines = fileobj.readlines()
    positions = []
    symbols = []
    ffsymbols=[]
    charges = []
    bonds = []
#    for L in lines:
    while len(lines)>0:
     Line=lines[0].split()
     del lines[0]
     if len(Line)>0:
      if Line[0]=="@<TRIPOS>MOLECULE":     
       del lines[0]
       Line=lines[0].split()
       del lines[0]
       numa = int(Line[0])
       numb = int(Line[1])
      
      if Line[0]=="@<TRIPOS>ATOM":
       for i in range(numa):
        Line=lines[0].split()
        del lines[0]
        
        symbol = Line[1]
        symbol=subnumber.sub('', symbol)
        if len(symbol) > 1 and (symbol[1].isupper() or symbol[1].isdigit()):
         symbol = symbol[0]
        ffsymbol = Line[5]
                                           
        x, y, z = Line[2:5]
        c = Line[8]
        charges.append(float(c))
        symbols.append(symbol)
        ffsymbols.append(ffsymbol)
        positions.append([float(x), float(y), float(z)])
            
      if Line[0]=="@<TRIPOS>BOND":
       for i in range(numb):
        Line=lines[0].split()
        del lines[0] 
         
        b1, b2 = Line[1:3]
        if Line[3]=="ar":
         bord=1.5
        elif Line[3]=="am":
         bord=1
        else:
         bord=float(Line[3])
        bonds.append( [min(int(b1)-1,int(b2)-1) ,max(int(b1)-1,int(b2)-1) ,bord])

    At=Atoms(symbols=symbols, positions=positions) 
    At.set_bonds(bonds)
    At.set_ffsymbols(ffsymbols)
    At.set_charges(charges)
    return At

def write_mol2_charges(fileobj, images):
    if isinstance(fileobj, str):
        fileobj = paropen(fileobj, 'w')

    if not isinstance(images, (list, tuple)):
        images = [images]

    symbols = images[0].get_chemical_symbols()
    natoms = len(symbols)
    for atoms in images:
        fileobj.write('%d\n' % natoms)
        Strcell=""
        for c in atoms.get_cell()[:]:
            Strcell+=('%f %f %f ' % (c[0],c[1],c[2]))
        fileobj.write(Strcell+"\n")

        for s, (x, y, z), c in zip(symbols, atoms.get_positions(), atoms.get_charges()):
            fileobj.write('%-2s %22.15f %22.15f %22.15f %22.15f\n' % (s, x, y, z, c))


def write_mol2(fileobj, images, ffsymbols=[]):
    if isinstance(fileobj, str):
        fileobj = open(fileobj, 'w')

    if not isinstance(images, (list, tuple)):
        images = [images]

    symbols = images[0].get_chemical_symbols()
    if ffsymbols==[]:
      ffsymbols = symbols
      
    natoms = len(symbols)
    for atoms in images:
        header = """@<TRIPOS>MOLECULE
*****
 %d %d 0 0 0
SMALL
GASTEIGER

""" % (atoms.get_number_of_atoms(), len(atoms.get_bonds()))

        fileobj.write(header)
        Strcell=""
        for c in atoms.get_cell()[:]:
            Strcell+=('%f %f %f ' % (c[0],c[1],c[2]))
        #fileobj.write(Strcell+"\n")

        stratoms = "@<TRIPOS>ATOM\n"
        id=1
        for s, (x, y, z), fs, c in zip(symbols, atoms.get_positions(), ffsymbols, atoms.get_charges()):
            stratoms += '%-5d %5s %12.7f %12.7f %12.7f %5s 1 UNK   %12.7f\n' % (id, s, x, y, z, fs, c)
            id+=1
        fileobj.write(stratoms)
        
        strbonds = "@<TRIPOS>BOND\n"
        for iB, B in enumerate(atoms.get_bonds()):
            strbonds += '%-5d %5d %5d %5d \n' % (iB+1, B[0]+1, B[1]+1, int(B[2]))
        fileobj.write(strbonds)
        