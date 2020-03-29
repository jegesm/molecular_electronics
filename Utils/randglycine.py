import os, sys
import subprocess
from numpy import random
import shlex
from string import Template

ASE_DIR = os.getenv("ASE_DIR")
DREI_DIR = os.getenv("DREI_DIR")
sys.path.append(ASE_DIR)
sys.path.append(DREI_DIR)
from ase_ext import io
import Dreiding_forcefield as DFF

fn_relaxed_xyz = "relaxed.xyz"
lammpsdrei_template = """units           real
atom_style      full
boundary        p p p
lattice         sc 1.0

pair_style      lj/cut/coul/long 6.0 6.0
bond_style      harmonic
angle_style     cosine/squared
dihedral_style  harmonic
improper_style  umbrella
kspace_style    pppm  5e-05

read_data       data
timestep 0.05
thermo 500

fix    1 all nvt temp 200. 200. 300.0
minimize 1.0e-4 1.0e-6 100 1000
run 1

write_dump all xyz  %s modify element $SYM
""" % (fn_relaxed_xyz)


def gen_randglycine(length, num_randunit):
    letters = ["P", "A", "V", "L", "I", "M", "C", "F", "Y", "W", "H", "K", "R", "Q", "N", "E", "D", "S", "T"]
    base_str = "Q" + length * "G"
    newstr = base_str
    randlist = []
    randletter = []

    while len(list(set(randlist))) < num_randunit:
        rp = random.randint(2, len(base_str)-2)
        rl = random.randint(0, len(letters))
        randlist.append(rp)
        randletter.append(letters[rl])
        randlist = list(set(randlist))
    randletter = randletter[:len(randlist)]

    for R, L in zip(randlist, randletter):
         newstr = newstr[:R] + L + newstr[R+1:]

    return gen_mol_from_fasta(newstr)

def gen_mol_from_fasta(fasta_str):

    temp_dir = os.getenv('TEMP_DIR')
    if temp_dir:
        temp_dir += "/"+fasta_str
    else:
        temp_dir = os.path.join("/tmp", fasta_str)

    try:
        os.mkdir(temp_dir)
    except:
        pass

    fdfn = os.path.join(temp_dir, fasta_str)
    with open(fdfn+".fasta", 'w') as f:
        f.write(fasta_str+"\n")

    mol2out = fdfn+".mol2"
    mol2cout = fdfn+"-cut.mol2"
    babel = os.getenv("BABEL_EXE")
    subprocess.call(shlex.split('%s -ifasta %s.fasta --join -omol2 %s' % (babel, fdfn, mol2out)))
    init_mol2 = io.read(mol2out, format="mol2")
    nat = init_mol2.get_number_of_atoms()
    bonds = init_mol2.get_bonds()
    newbonds=[]
    ie = 18   #this is the atom number of Q
    iv=nat-2    #this is the atom number of the last G

    for b in bonds:
        if (b[0] > ie-1 and b[0] < iv) and (b[1] > ie-1 and b[1] < iv):
            newbonds.append([b[0]-ie, b[1]-ie, b[2]])

    newmol = init_mol2[ie:iv]
    newmol.set_bonds(newbonds)
    newmol.set_ffsymbols(init_mol2.get_ffsymbols()[ie:iv])
    newmol.write(mol2cout, format="mol2", symbols=newmol.get_ffsymbols())

    drei_param = os.path.join(DREI_DIR, "DREI_PARAM")
    drei_outdir = os.path.join(temp_dir, "drei")
  
    try:
        os.mkdir(drei_outdir)
    except:
        pass


    DFF.create_input(mol2cout, lammpsdrei_template, drei_param, drei_outdir)

    drei_comm = os.getenv("LAMMPS_EXE")
    #subprocess.call(shlex.split('bash -c "ls %s; cd %s; %s < in.drei |tee out"'%(drei_comm, drei_outdir, drei_comm)))
    subprocess.call(shlex.split('bash -c "cd %s; %s < in.drei > out"' % (drei_outdir, drei_comm)))
    #subprocess.call(shlex.split('bash -c "cp %s ./"' % os.path.join(drei_outdir, relaxed_xyz)))

    relaxed_xyz = io.read(os.path.join(drei_outdir, fn_relaxed_xyz), 0)[0]
    return relaxed_xyz, init_mol2, fasta_str[1:-1]
