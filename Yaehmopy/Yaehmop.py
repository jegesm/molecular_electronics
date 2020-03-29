from numpy import *
import os
import ase
from ase import io


class Yaehmop():
    """

    :param mol:  the file that contains the coordinates of the system of an ASE.Atoms object
    :param yaehmop_binary: alternative path to Yaehmop binary, otherwise reads YAEHMOP_DIR
    :param syslabel: label for existing files or set it for the newone
    :param sparse: defines what to calculate or look for
    :param symbols: if the input coordinate file has weird symbols that might not be recognised then you can give alternatives
    """
    yaehmop_binary = ""
    sparse = False
    syslabel = ""
    #yahlabel = ""
    mol = []
    symbols = []
    nat = 0
    
    @property
    def label_(self):
        return "%s.yah"%self.syslabel
        
    @property
    def label_ham_(self):
        if self.sparse:
            return self.label_ + ".SPARSE.HAM"
        else:
            return self.label_ + ".HAM"
    @property
    def label_ov_(self):
        if self.sparse:
            return self.label_ + ".SPARSE.OV"
        else:
            return self.label_ + ".OV"
        
    def __init__(self, mol="", yaehmop_binary="", syslabel="test", sparse=False, symbols=""):

        self.yaehmop_binary = yaehmop_binary
        if not os.path.isfile(self.yaehmop_binary):
            self.yaehmop_binary = os.path.join(os.getenv("YAEHMOP_DIR"), 'bind')
            if not self.yaehmop_binary:
                print(" WARNING: Yaehmop binary %s is missing"%self.yaehmop_binary)

        self.sparse = sparse
        self.syslabel = syslabel
        
        if symbols:
            self.symbols = symbols

        if mol:
            self.mol = mol
        #else:
        #    print("The xyz coordinates should have a shape of (natoms, 3)")
            
            self.nat = len(self.mol)
        #self.yahlabel = "%s.yah"%self.syslabel
        
    def load_mol(self, molname):
        """

        :param molname: the file that contains the coordinates of the system of an ASE.Atoms object
        :return: nothing
        """
        if type(molname)==ase.atoms.Atoms:
            self.mol = molname
        elif molname[-4:]=="mol2":
            self.mol = read_mol2(molname)
            self.syslabel = molname.split(".")[0]
        else :
            self.mol = io.read(molname)
            self.syslabel = os.path.basename(molname).split(".")[0]
        
        self.nat = len(self.mol)
        
        
        
    def load_ham(self, ov=False, size_of_int=uint32):
        """

        :param ov: which matrix to read: overlap if True, hamiltonian if False
        :param size_of_int: for older yaehmop codes, if you receive an empty matrix set this to uint64
        :return: the matrix stored in syslabel.HAM or syslabel.OV
        """
        if not ov:
            fp = open(self.label_ham_, 'rb')
            if not os.path.isfile(self.label_ham_):
                raise
        else:
            fp = open(self.label_ov_, 'rb')
        [dummy, n] = fromfile(fp, dtype=size_of_int, count=2)
        hamiltonian = fromfile(fp, dtype=double, count=n*n)
        hamiltonian = hamiltonian.reshape(n,n)
        hamiltonian = hamiltonian + hamiltonian.T - diag(diag(hamiltonian))
        fp.close()
        return hamiltonian

    def load_sparseham(self, ov=False):
        """

        :param ov: which matrix to read: overlap if True, hamiltonian if False
        :return: the sparse matrix stored in syslabel.SPARSE.HAM or syslabel.SPARSE.OV
        """
        from scipy import sparse
        if not self.sparse:
            print("SPARSE was not set!!!")
            return 
        if not ov:
            T = loadtxt(self.label_ham_,skiprows=1)
        else:
            T = loadtxt(self.label_ov_,skiprows=1)
        
        dim = int(max(T[:,0]))+1
        Hs = sparse.coo_matrix( (T[:,2],(T[:,0],T[:,1])), shape=(dim,dim) )
        dim = Hs.shape[0]
        dd = sparse.diags([1],[0],(dim,dim))
        Hdiag=Hs.multiply(dd)
        sHs=(Hs+Hs.transpose()-Hdiag).tocoo()
        return sHs

    def create_input(self, just_matrices=True):
        """

        :param just_matrices: the default is to not calculate the eigenvalues, because that could take up a lot of time
        :return: nothing, just creates the ipnut file for Yaehmop
        """
        
        Data="%s\nMolecular\nGeometry\n%d\n"%(self.syslabel, self.nat)
        if type(self.mol)==type(ase.Atoms()):
            for A in self.mol:
                Data+="{0}\t{1}\t{2}\t{3}\n".format(A.symbol, A.x, A.y, A.z)
        elif self.symbols:
            for i in range(self.nat):
                Data+="%s\t%f\t%f\t%f\n"(self.symbols[i]. self.mol[0], self.mol[1], self.mol[2])
        else:
             print("ERROR: you have to define the symbols too!!!")

        if just_matrices:
           Data+="Charge\n0\nJust Matrices\nPrint\nOrbital Mapping\nWave Functions\n\nEnd_print\n"
        else:
           Data+="Charge\n0\nPrint\nOrbital Mapping\nWave Functions\n\nEnd_print\n"
        if self.sparse:
            Data += "Dump Sparse\n"
        else:
            Data += "Dump Hamil\nDump Overlap\n"

        #Data+="Charge\n0\nJust Matrices\nPrint\nOrbital Mapping\n\nEnd_print\nDump Sparse"
        #Data+="Charge\n0\nPrint\nOrbital Mapping\n\nEnd_print\nDump Hamil\nDump Overlap"

        Fyah = open(self.label_, 'w')
        Fyah.write(Data)
        Fyah.close()

    def run(self):
        """

        :return: writes the log text to stdout
        """
        import subprocess
        #xyz=io.read(fname) 
        #!rm $yahname""*
        #self.create_input()        
        print(subprocess.getoutput("%s %s"%(self.yaehmop_binary, self.label_)))

    def getBlockFromOut(self, blockname):
        """

        :param blockname: which block to recover from the output
        :return: the given block as numpy array
        """
        Data = open(self.label_+".out",'r').readlines()
        nlfirst=-1
        nllast=-1
        for i in range(len(Data)):
            if nlfirst>-1:
                if Data[i].rfind(blockname)>-1:
                    nllast=i
                    break
            if Data[i].rfind(blockname)>-1:
                nlfirst=i+1

        return Data[nlfirst:nllast]

    def GetIorbs(fnout,mol):
        Data = getDataFromTextFile(fnout,"; Orbital Mapping","END")

        #Residue names
        #mol = Py.readfile("mol2",fnmol2).next()
        #mol = ase.io.read(fnmol2)

        residues=[]
        for a in mol:
            residues.append(a.OBAtom.GetResidue().GetName())

        #Coordinates
        coords=[A.coords for A in mol.atoms]
        try:
            iorbs = [ [int(L.split()[2][:-1]),L.split()[1][:-1]] for L in Data]
            iorbs = [ [L[0],L[1],residues[L[0]-1],coords[L[0]-1]] for L in iorbs ]
        except ValueError:
            print(L)
            exit()

        return iorbs

    def IorbstoRes(mol,iorbs,Seq=[]):
        if Seq==[]:
            if type(mol)=="str":
                tmol = Py.readfile("mol2",mol).next()
                tmol = ase.io.read(mol)
                mol=tmol
            residues=[]
            for a in mol:
                residues.append(a.OBAtom.GetResidue().GetName())

            sortresidues = SortRes(residues)
            residues=array(residues)
            #Uresidues=unique(residues)
            Seq=list(set(residues))


        Resorbs={}
        #for Res in Uresidues:
        for Res in Seq:
            tempri=[]
            for i in range(len(iorbs)):
                if iorbs[i][2]==Res:
                    tempri.append(i)            
            Resorbs[Res]=array(tempri)

        return Resorbs

    def SeqConnectivity(mol,plot=False):
        residues=[]
        for a in mol:
            residues.append(a.OBAtom.GetResidue().GetName())

        residues = list(set(residues))
        Sresidues=range(len(residues))

        Rdist=2
        Cmat = {}
        for Res in residues:
            Cmat[Res]=[]

        for a in mol:
            resa = a.OBAtom.GetResidue().GetName()
            ac = a.coords
            for b in mol:
                resb = b.OBAtom.GetResidue().GetName()
                bc = b.coords
                RR = ((ac[0]-bc[0])**2+(ac[1]-bc[1])**2+(ac[2]-bc[2])**2)
                if  RR < Rdist and (a.atomicnum*b.atomicnum>1):
                    Cmat[resa].append(resb)

        for K in Cmat.keys():
            Cmat[K] = list(set(Cmat[K]))
            Cmat[K].remove(K)

        for K in Cmat.keys():
            if len(Cmat[K])>2:
                print("WARNING: Amino-acid connected to too many!!!")
                print( Cmat)
                for K in Cmat.keys():
                    print(K,list(set(Cmat[K])))


        #ENDINGS
        First=-1
        for K in Cmat.keys():
            if len(Cmat[K])==1:
                First=K

        #print Cmat
        #GET SEQUENCE
        Seq=[First]
        #while len(Cmat[First])!=0:
        while len(Seq)!=len(residues):
            Next = Cmat[First][0]
            Cmat[Next].remove(First)
            Seq.append(Next)
            First=Next

        if plot:
            print("%14s"%mol.title),
            for S in Seq[:-1]:
                print(S+" --"),
            print(Seq[-1])

        #REVERSE THE SEQUENCE SO IT STARTS WITH THE FIRST LETTER
        if int(Seq[0][3])!=1:
            Seq.reverse()

        return Seq

    def GetEig(self):
        return loadtxt(self.getBlockFromOut("OccupationNumbers"))

    def GetFermi(self):
        eigs = self.GetEig()
        imin = where(eigs[:,2]<0.001)[0][0]
        imax = where(eigs[:,2]>0)[0][-1]
        return (eigs[imax, 1] + eigs[imin, 1]) / 2

    def GetDos(self, width=0.1, emin=None, emax=None, npts=401):
        """Electronic Density Of States object.


        :param width: float
            Width of guassian smearing.  Use width=0.0 for linear tetrahedron
            interpolation.

        :param emin: float
        :param emax: float
        :param npts: int
            Number of points.
        :return: DOS
        """

        npts = npts
        #w_k = calc.get_k_point_weights()
        e_skn = self.GetEig()[:,1] - self.GetFermi()

        if emin is None:
            emin = e_skn.min() - 5 * width
        if emax is None:
            emax = e_skn.max() + 5 * width

        energies = linspace(emin, emax, npts)

        def delta(energy):
            """Return a delta-function centered at 'energy'."""
            x = -((energies - energy) / width) ** 2
            return exp(x) / (sqrt(pi) * width)


        """Get array of DOS values.

        The *spin* argument can be 0 or 1 (spin up or down) - if not
        specified, the total DOS is returned.
        """


        dos = zeros((npts,2))
        dos[:,0] = energies
        for e in e_skn:
            #for e in e_n:
                dos[:,1] += delta(e)
        return dos

