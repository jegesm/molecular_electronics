#!/usr/bin/python3
import sys
import subprocess as sp
import argparse
import signal
import sys

#Let underlying signal_handler tend to the interruption of the program
def signal_handler(signal, frame):
    return

def createinput( arg1, arg2 ):
    process = sp.Popen(["./input_pdb",arg1,arg2])
    (output, err) = process.communicate()
    exit_code = process.wait()

def huckel(arg):
    process = sp.Popen(["./bind", arg])
    (output, err) = process.communicate()
    exit_code = process.wait()

def configinfo():
    print('The configuration file contains lines with printing options of the following syntax:')
    print('[printing keyword]=[true|false]')
    print('Example: overlap=true')

def printlist():
    print('Printing options to be included in the configuration file:')
    print('distance matrix             Print the distance matrix')
    print('overlap population          Print the Mulliken overlap population matrix')
    print('reduced overlap population  Print the Mulliken reduced overlap population matrix')
    print('charge matrix               Print the charge matrix')
    print('wave functions              Print the wave functions for the molecule')
    print('net charges                 Print the net charges on the atoms, as determined using Mulliken population analysis')
    print('overlap                     Print the overlap matrix')
    print('hamil                       Print the hamiltonian matrix')
    print('electrostatic               Print the electrostatic contribution to the total energy')
    print('levels                      Toggles the printing of energy levels at each k point in an extended calculation')
    print('fermi                       Print the Fermi energy (useful in combination with the Walsh option)')
    print('orbital energy              Allows the energy of a particular orbital to be printed')
    print('orbital coeff               Allows the coefficient of a particular atomic orbital in a given molecular orbital to be printed')
    print('orbital mapping             Generates the scheme used to number the individual atomic orbitals in a calculation')
    print('dump overlap                Creates binary output for the overlap matrix')
    print('dump hamil                  Creates binary output for the Hamiltonian matrix')
    print('dump sparse                 Creates sparse formatted binary output for the Hamiltonian and overlap matrices')
    print('dump matrix market          Creates sparse Matrix Market formmated binary output for the Hamiltonian and overlap matrices')

def main(argv):
    signal.signal(signal.SIGINT, signal_handler)
    parser = argparse.ArgumentParser(description='Welcome to the Huckel Code application.')
    parser.add_argument('-i','--input', nargs = 2)
    parser.add_argument('-c','--calculate')
    parser.add_argument('-f','--configinfo', action='store_true')
    parser.add_argument('-p','--printlist', action='store_true')
    if len(sys.argv[1:])==0:
        parser.print_help()
        parser.exit()
    else:
        try:
            arguments = parser.parse_args()
        except:
            parser.print_help()
            parser.exit()
        if(arguments.input):
            createinput(arguments.input[0],arguments.input[1])
        elif(arguments.calculate):
            huckel(arguments.calculate)
        elif(arguments.printlist):
            printlist()
        elif(arguments.configinfo):
            configinfo()

if __name__ == "__main__":
    main(sys.argv[1:])
