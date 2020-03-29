from string import Template
from numpy import *
from scipy import spatial
import sys, os
import ase
from ase.io import read

if len(sys.argv)!=2:
    print("""\n Generates polyglycine chain according to the ccordinates given in Fujita et al.  https://doi.org/10.1063/1.1673987 \n An 18 unit long chain gives a translational periodic piece\n\nUsage: \n\t> python3 %s number_of_units  \n\n\te.g. python3  %s 18\n"""%(sys.argv[0], sys.argv[0]))
    exit()


sstr="""
N	0.53887	1.49842	0.62146
H	0.37153	1.63083	-0.35551
C	1.48281	0.70221	1.08774
O	1.72237	0.51072	2.28922
C	2.29804	0.00000	0.00000
H	2.78556	0.68274	-0.54423
H	2.95737	-0.61988	0.42547"""