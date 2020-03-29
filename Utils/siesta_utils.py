from string import Template
from numpy import *
from scipy import spatial
import sys

siesta_temp = Template("""
WriteCoorInitial T
LatticeConstant 1.0 Ang
DM.MixingWeight 0.1
WriteMDhistory T
Number_of_species $NSP
SystemLabel Dens
PAO.BasisSize $BASISSIZE
SolutionMethod diagon
AtomicCoordinatesFormat Ang
MaxSCFIterations 120
DM.NumberPulay 4
MeshCutoff 500.000000 eV
NumberOfAtoms $NAT

PAO.EnergyShift         0.002 eV

%block LatticeVectors
$CELL
%endblock LatticeVectors

%block Chemical_Species_label
$SPECIES%endblock Chemical_Species_label

%block AtomicCoordinatesAndAtomicSpecies < coords
#%endblock AtomicCoordinatesAndAtomicSpecies

#COOP.write T
#SaveHS T
#WriteKpoints T
#WriteKbands T

#UseSaveData T
MD.UseSaveXV T
MD.UseSaveCG T

SaveRho T

MD.TypeOfRun         cg
MD.NumCGsteps        0
MD.MaxForceTol         0.1555740000E-01 Ry/Bohr
WriteCoorXmol        T
WriteForces          T
AllocReportLevel     2
WriteMDhistory       T
WriteMDXmol          T

xc.functional           LDA
xc.authors              CA
""")

Masses=[1.01,  4.00,  6.94,  9.01, 10.81, 12.01, \
    14.01, 16.00, 19.00, 20.18, 22.99, 24.31,  \
    26.98, 28.09, 30.97, 32.07, 35.45, 39.95, \
    39.10, 40.08, 44.96, 47.88, 50.94, 52.00, \
    54.94, 55.85, 58.93, 58.69, 63.55, 65.39, \
    69.72, 72.61, 74.92, 78.96, 79.90, 83.80, \
    85.47, 87.62, 88.91, 91.22, 92.91, 95.94, \
    98.91,101.07,102.91,106.42,107.87,112.41, \
    114.82,118.71,121.75,127.60,126.90,131.29, \
   132.91,137.33,138.91,140.12,140.91,144.24, \
   146.92,150.36,151.97,157.25,158.93,162.50, \
   164.93,167.26,168.93,173.04,174.97,178.49, \
   180.95,183.85,186.21,190.20,192.22,195.08, \
   196.97,200.59,204.38,207.20,208.98,208.98, \
   209.99,222.02,223.02,226.03,227.03,232.04, \
   231.04,238.03,237.05,244.06];

Symbols = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",\
 "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", \
 "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",\
 "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc",\
 "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",\
 "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",\
 "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", \
 "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", \
 "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es",\
 "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt"]

def read_symbols(fout):
    data = open(fout, 'r').readlines()
    nums = 0
    for d in data:
        if d.rfind("Number_of_species") > -1:
            nums = int(d.split()[1])
            break

    numl = 0
    for id, d in enumerate(data):
        if d.rfind("Chemical_Species_label") > -1:
            numl = id+1
            break

    syms = {}
    for d in data[numl:numl+nums]:
        syms[d.split()[0]] = d.split()[2]

    return syms

def read_struct(fstr, fout):
    dstr = open(fstr, 'r').readlines()
    syms = read_symbols(fout)
    cstr = str(dstr[0][:-1])+str(dstr[1][:-1])+str(dstr[2][:-1])
    cell = fromstring(cstr, dtype=float, sep=' ').reshape((3,3))
    nat = int(dstr[3].split()[0])
    dstr = dstr[4:]

    coo = []
    sym = []
    for ia, a in enumerate(dstr):
        ll = a.split()
        sym.append(syms[ll[0]])
        coo.append([float(ll[2])*cell[0,0], float(ll[3])*cell[1,1], float(ll[4])*cell[2,2]])
        if ia==nat:
            break

    coo = array(coo)
    return coo, sym, cell

def read_xyz(fxyz):
    dxyz=open(fxyz,'r').readlines()

    nat=int(dxyz[0])
    others=dxyz[1]
    dxyz=dxyz[2:]

    coo = []
    sym = []
    for i in range(nat):
        ll = dxyz[i].split()
        sym.append(ll[0])
        coo.append([float(ll[1]), float(ll[2]), float(ll[3])])

    coo = array(coo)
    return coo, sym, others

def get_cell_from_xyz(data):
    cell = array([d for d in data.split()])
    cell.reshape((3,3))
    return cell

def get_cell(coo):
    mins = coo.min(axis=0)
    maxs = coo.max(axis=0)
    c=maxs-mins+15
    return "{0} 0 0\n0 {1} 0\n0 0 {2}".format(c[0], c[1], c[2])

def get_species(sym, PT):
    tt=list(set(sym))
    tt.sort()
    species=[]
    for i,s in enumerate(tt):
        species.append([i+1,s,PT[s][0]+1])
    return species

