#!/usr/bin/python

import sys, os
ASE_DIR = os.getenv("ASE_DIR")
sys.path.append(ASE_DIR)
from ase_ext.io.xsf import read_xsf
from ase_ext.io.cube import write_cube

if len(sys.argv)!=2:
    print("""\n Converts an XSF file into CUBE \n\nUsage: \n\t> python3 %s XSFfile   \n\n\te.g. python3  %s h2o-relaxation.XSF\n"""%(sys.argv[0], sys.argv[0]))
    exit()

fn=sys.argv[1]
xsf = read_xsf(fn, read_data=True)

print(xsf[1])
#write_cube(sys.stdout,xsf[1],xsf[0])
write_cube(open(fn[:-3]+"cube", 'w'), data=xsf[0], atoms=xsf[1], comment="ldos")