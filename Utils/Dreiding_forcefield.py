import sys, os

ASE_DIR = os.getenv("ASE_DIR")
sys.path.append(ASE_DIR)
from ase_ext import io
from string import Template

from numpy import zeros, sin, sqrt, dot, linalg, ones
PI = 3.14159265

def DistancePer(vec, cell, icell=[]):
    
    ss = dot(icell, vec)
    ss[0]=ss[0]-round(ss[0])
    ss[1]=ss[1]-round(ss[1])
    ss[2]=ss[2]-round(ss[2])
    vec = dot(cell, ss)
    return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2])
    
#WITH LINKCELL ALGORTIHM
def GetBonds(pos, cell,  bondpc, symbols,  Param):
    ncl=12;
    rsq=0.0;

    #Bond Tolerance
    bndfac=1.0 + bondpc/100.0;

    #Min-Max
    maxx, maxy, maxz = pos.max(axis=0)
    minx, miny, minz = pos.min(axis=0)
    
    #GET RADISUSES
    Radius=[]
    for s in symbols:
        Radius.append(Param[s]['Radius'])
        
    #Introducing linkcells
    cellsize = 8.0;
    lcellx = int((maxx-minx+0.0001)/cellsize+1)
    lcelly = int((maxy-miny+0.0001)/cellsize+1)
    lcellz = int((maxz-minz+0.0001)/cellsize+1)
        
    icell=linalg.inv(cell)
    
    if(lcellx<3):  
        lcellx=3;
    if(lcelly<3):
        lcelly=3;
    if(lcellz<3):
        lcellz=3;
    linkcellsize = lcellx*lcelly*lcellz;
    
    linkcell1=[]
    
    numa=len(pos)
    il=0;
    for iz in range(lcellz):
        for iy in range(lcelly):
            for ix in range(lcellx):
                tlinkcell=[]
                for m in range(numa):
                    if((pos[m][0]>=ix*cellsize+minx) & 
                            (pos[m][0]<(ix+1)*cellsize+minx) &
                            (pos[m][1]>=iy*cellsize+miny) & 
                            (pos[m][1]<(iy+1)*cellsize+miny) & 
                            (pos[m][2]>=iz*cellsize+minz) &
                            (pos[m][2]<(iz+1)*cellsize+minz)):
                        
                            tlinkcell.append(m);
                            
                linkcell1.append(tlinkcell)
                il+=1
                
    tlinkc=-ones((numa, ncl), dtype=int)
    tBnd=[]
    tlist=zeros(numa, dtype=int)
    
    # determine bonds in molecule
    nBnd=0;
    for il in range (linkcellsize):
        #WE NEED TO GATHER ALL ATOMS IN THE 8 NEIGHBOURING CELLS
        neighbors = []
        for isy in range(-1, 2):
            for isx in range(-1, 2):
                for isz in range(-1, 2):
                    celln=(il+isz+isx*lcellx+isy*lcelly*lcellx+linkcellsize)%linkcellsize;
                    neighbors.extend(linkcell1[celln]);
        
        neighbors.sort()
        neighbors=list(set(neighbors))
        
        for m in range(len(linkcell1[il])-1, -1, -1):
            for n in range(len(neighbors)):
                matm=linkcell1[il][m];
                natm=neighbors[n];
                if(matm>natm):
                    ff=bndfac*(Radius[matm]+Radius[natm])
                    xdis=pos[matm][0]-pos[natm][0]
                    ydis=pos[matm][1]-pos[natm][1]
                    zdis=pos[matm][2]-pos[natm][2]
                    rsq = DistancePer((xdis, ydis, zdis), cell, icell);
                    
                    if(rsq<0.1):
                       print( "How can be 2 atoms in the same place? Atom ", matm," ",  natm,  "\n")
                    if(rsq<=ff):
                        ttBnd=[]
                        ttBnd.append(min(matm,natm));
                        ttBnd.append(max(matm,natm));
                        ttBnd.append(1)
                        tBnd.append(ttBnd);                            
                        kkk=max(tlist[matm],tlist[natm]);
                        if(kkk==ncl):
                            print( "Warning!!! HOW CAN 1 ATOM HAVE MORE THAN 12 BONDS?  ", matm, " " ,natm,  "\n")
                        
                        tlinkc[matm][tlist[matm]]=natm;
                        tlinkc[natm][tlist[natm]]=matm;
                        tlist[matm]+=1;
                        tlist[natm]+=1;
                        
    #WE NEED TO SORT LINKC
    for mm in range(len(tlinkc)):
        for mi in range(len(tlinkc[mm])):
            for mj in range(mi+1, len(tlinkc[mm])):
                if((tlinkc[mm][mi]>tlinkc[mm][mj]) & (tlinkc[mm][mj]>-1)):
                    temp = tlinkc[mm][mi]
                    tlinkc[mm][mi]=tlinkc[mm][mj];
                    tlinkc[mm][mi]=temp;
                
    #AND THEN CONVERT -1 TO 0  
    for mm in range(len(tlinkc)):
        for mmm in range(len(tlinkc[mm])):
            if(tlinkc[mm][mmm]==-1):
                tlinkc[mm][mmm]=0;
                        
    iBnd=tBnd;
    nBnd=len(iBnd);
#    for B in iBnd:
#        print( B)
    return iBnd, tlist
    
def GetLinkMat(bonds):
 linkc={}
 for B in bonds:
  if B[0] in linkc.keys():
   linkc[B[0]].append(B[1])
  else:
   linkc[B[0]]=[]
   linkc[B[0]].append(B[1])
  
  if B[1] in linkc.keys():
   linkc[B[1]].append(B[0])
  else:
   linkc[B[1]]=[]
   linkc[B[1]].append(B[0])
   
 return linkc
  
def GetAngles(linkc):
 Ang=[]
 for K in linkc.keys():
  for N1 in linkc[K]:
   for N2 in linkc[K]:  
    if N1!=N2:
     Ang.append( (min(N1,N2),K,max(N1,N2)) )

 Ang=set(Ang)
 return list(Ang)

def GetDihedrals(bonds,linkc):
 Dih=[]
 for B in bonds:
  for N1 in linkc[B[0]]:
   for N2 in linkc[B[1]]:  
    if (N1!=N2) and (N1!=B[1]) and (N2!=B[0]):
     if (N1<N2):
      Dih.append((N1,B[0],B[1],N2,B[2]))
     else:
      Dih.append((N2,B[1],B[0],N1,B[2]))

 Dih=set(Dih)
 return list(Dih)
 
def GetInversions(linkc,ffsymbols):
 Inv=[]   
 for K in linkc.keys():
  if len(ffsymbols[K].split("_"))>1:
   if (len(linkc[K]) == 3) & (ffsymbols[K].split("_")[1] in ["2","R"]):
     Inv.append((K, linkc[K][0], linkc[K][1], linkc[K][2]))
     
 Inv=set(Inv)
 return list(Inv)
                     
   
def ReadParameters(fname):
    Param = {}
    try:
     WParam = open(fname, 'r').readlines()
    except TypeError:
        print( "DREI_PARAM ERROR!\nPlease copy here the file that contains the parameters!")
        return

    for Line in WParam:
        L = Line.split()
        if (len(L) > 7) & (L[0][0] != "#"):
            #Symbol    At_mass    Charge    Radius    Angle    Epsilon     Sigma    Zeta
            Tparam = {'Symbol': L[0], 
                'Mass': float(L[1]), 
                'Charge': float(L[2]), 
                'Radius': float(L[3]), 
                'Angle': float(L[4]), 
                'Epsilon': float(L[5]), 
                'Sigma': float(L[6]), 
                'Zeta': float(L[7]) }
            Param[L[0]]=Tparam

    return Param

def  Dreidingbndass(bonds, ffsymbols, Param, type="harm"):
#// *********************************************************************
#// routine for assigning bond parameters in the dreiding force field
#// *********************************************************************
    Bonds = {}
    for Bpair in bonds:
        T1 = ffsymbols[Bpair[0]]
        T2 = ffsymbols[Bpair[1]]
        FFtype = T1+"-"+T2
     
        if FFtype in Bonds.keys():
            Bonds[FFtype]['id'].append(Bpair)
        else:
         
            P1 = Param[T1]
            P2 = Param[T2]
    
            params=[0,0,0]
            name = ""
            if(type=="morse"):
                name = "mors"
                params[0]=70.0* Bpair[2]
                params[1]=P1['Radius']+P2['Radius']-0.01
                params[2]=sqrt(5.0)
            elif(type=="harm"):
                name="harm"
                params[0]=700.0* Bpair[2]
                params[1]=P1['Radius']+P2['Radius']-0.01
                params[2]=0.0
         

            Bonds[FFtype]={'name':name,'params':params,'id':[Bpair], 'lab':-1}
    
    lab=1
    for B in  Bonds.keys():
     Bonds[B]['lab']=lab
     lab+=1
    return Bonds


def Dreidingangass(angles,ffsymbols, Param):
#// *********************************************************************
#// routine for assigning valence angle parameters in the dreiding force field
#// *********************************************************************
    Angles = {}
    for Apair in angles:
     T1 = ffsymbols[Apair[0]]
     T2 = ffsymbols[Apair[1]]
     T3 = ffsymbols[Apair[2]]
     FFtype = T1+"-"+T2+"-"+T3
                   
     P1 = Param[T1]
     P2 = Param[T2]
     P3 = Param[T3]
     th0=P2['Angle']
     params=[0,0,0]
     name = ""
        
     if (abs(th0-180.0) < 0.000001):
      name="cos "
      params[0]=100.0
      params[1]=0.0
      params[2]=1.0
     else:
      name="hcos"
      params[0]=100.0/pow(sin(PI*th0/180.0),2)
      params[1]=th0
      params[2]=0.0

     if FFtype in Angles.keys():
      Angles[FFtype]['id'].append(Apair)
     else:
      Angles[FFtype]={'name':name,'params':params,'id':[Apair], 'lab':-1}
    
    lab=1
    for A in  Angles.keys():
     Angles[A]['lab']=lab
     lab+=1
    return Angles      
    

def Dreidingdihass(dihedrals,ffsymbols, linkc, Param):
    #//*********************************************************************
    #// routine for assigning dihedral parameters in the dreiding force field
    #// *********************************************************************
    Dih = {}
    for Dquart in dihedrals:
     T1 = ffsymbols[Dquart[0]]
     T2 = ffsymbols[Dquart[1]]
     T3 = ffsymbols[Dquart[2]]
     T4 = ffsymbols[Dquart[3]]
     bord=Dquart[4]
     FFtype = T1+"-"+T2+"-"+T3+"-"+T4
     
     if len(T1.split("_"))==0:
      ausI="M"
     else:
      ausI = T1.split("_")[1]
      busI = T1.split("_")[0]
     if len(T2.split("_"))==0:
      ausJ="M"
     else:
      ausJ = T2.split("_")[1]
      busJ = T2.split("_")[0]
     if len(T3.split("_"))==0:
      ausK="M"
     else:
      ausK = T3.split("_")[1]
      busK = T3.split("_")[0]
     
     degen=(len(linkc[Dquart[1]])-1)*(len(linkc[Dquart[2]])-1)
     name=""
     bent="NEM"
     #A case
     if (bord==1) & (ausJ=="3") & (ausK=="3"):
      V=2; N=3; F=180;
      bent="A"
     #B case
     if (bord==1) & (( (ausJ in ["2","R"]) & (ausK in ["3"])) | 
        ( (ausK in ["2","R"]) & (ausJ in ["3"]))):
      V=1; N=6; F=0;
      bent="B"
     #C case
     if (bord==2) & (ausJ=="2") & (ausK=="2"):
      V=45; N=2; F=180;
      bent="C"
     #D case
     if (bord==1.5) & (ausJ=="R") & (ausK=="R"):
      V=25; N=2; F=180;
      bent="D"
     #E case
     if (bord==1) & (ausJ in ["2","R"]) & (ausK in ["2","R"]):
      V=5; N=2; F=180;
      bent="E"
     #F case
     if (bord==1) & (ausJ=="R") & (ausK=="R"):
      V=10; N=2; F=180;
      bent="F"
     #G case
     if (ausJ in ["1","M"]) | (ausK in ["1","M"]) :
      V=0.0;
      bent="G"
     #H case
     if (bord==1) & (ausJ=="3") & (ausK=="3") & (busJ in ["O","S","Se","Te"]) & (busK in ["O","S","Se","Te"] ):
      V=2; N=2; F=90;
      bent="H"
     #I case
     if (((ausJ=="3") & (busJ in ["O","S","Se","Te"]) & (ausK in ["2","R"]) & (busK not in ["O","S","Se","Te"])) |
        ((ausK=="3") & (busK in ["O","S","Se","Te"]) & (ausJ in ["2","R"]) & (busJ not in ["O","S","Se","Te"]))):
      V=2; N=2; F=180
      bent="I"
     #J case
     if (bord==1) & (ausI not in ["2","R"]) & (( (ausJ in ["2","R"]) & (ausK in ["3"])) | 
        ((ausK in ["2","R"]) & (ausJ in ["3"]))):
      V=2; N=3; F=180;     
      bent="J"
         
#     print( bent,FFtype,degen,bord)
     try:
      params=[V/2./degen,N,N*F+180]     
     except UnboundLocalError:
         print( Dquart, ffsymbols[Dquart[0]], ffsymbols[Dquart[1]], ffsymbols[Dquart[2]], ffsymbols[Dquart[3]])
     if FFtype in Dih.keys():
      Dih[FFtype]['id'].append(Dquart)
     else:
      Dih[FFtype]={'name':name,'params':params,'id':[Dquart], 'lab':-1}

    lab=1
    for D in  Dih.keys():
     Dih[D]['lab']=lab
     lab+=1                      
    return Dih
    
def Dreidinginvass(inversions,ffsymbols, Param):
    #// *********************************************************************
    #// routine for assigning inversion parameters in the dreiding force field
    #// *********************************************************************
    Inv = {}
    for Iquart in inversions:  
     T1 = ffsymbols[Iquart[0]]
     T2 = ffsymbols[Iquart[1]]
     T3 = ffsymbols[Iquart[2]]
     T4 = ffsymbols[Iquart[3]]
     FFtype = T1+"-"+T2+"-"+T3+"-"+T4
     
     name = ""
     params=[40.0]                              
     if FFtype in Inv.keys():
      Inv[FFtype]['id'].append(Iquart)
     else:
      Inv[FFtype]={'name':name,'params':params,'id':[Iquart], 'lab':-1}

    lab=1
    for I in  Inv.keys():
     Inv[I]['lab']=lab
     lab+=1
    return Inv
    
def Syb2Dre(ffsymbols):
 sym={"Mn": "Mn", "Fe": "Fe", "Du": "Du", "C.3": "C_3", "C.2": "C_2", "C.2": "C_2", "C.1": "C_1", "C.2": "C_2", "C.2": "C_2", "C.3": "C_3", "N.4": "N_3", "N.3": "N_3","N.2": "N_2",
    "N.pl3": "N_2", "N.2": "N_2","N.1": "N_1","N.4": "N_2","N.pl3": "N_2","N.pl3": "N_2","O.3": "O_3","O.3": "O_2","O.2": "O_2","O.co2": "O_3","S.3": "S_3","S.3": "S_3","S.3": "S_3",
    "S.2": "S_3","S.3": "S_3","S.o2": "S_3","S.o": "S_3","S.3": "S_3","B": "B_2","B": "B_2","B": "B_3","P.3": "P_3","P.3": "P_3","P.3": "P_3","P.3": "P_3","P.3": "P_3","H": "H_","H": "H_",
    "H": "H_","H": "H_","H": "H_","H": "H_","F": "F_","Cl": "Cl","Br": "Br","I": "I_","Ge": "Ge3","Se": "Se3","Te": "Te3","O.3": "O_3","N.2": "N_2","N.2": "N_2","N": "N_2","C.ar": "C_R",
    "N.ar": "N_R", "N.am": "N_2","Na": "Na","Ca": "Ca","Li": "Li","Al": "Al3","Si": "Si3","O.co2": "O_2","C.cat": "C_2","O.2": "O_2","Pd": "Du","Pt": "Du","Ni": "Ni","Cu": "Cu","N.4": "N_3",
    "N.pl3N": "N_2","O.co2O": "O_3","C.catC": "C_2","O.": "O_3","C.": "C_3","V": "Du","Zn": "Zn","C.": "C","N.": "N","O.": "O","S.": "S","R.ar": "C","R.": "C","Nr.": "C","O.co2": "O_3","S.o2": "S_3",
    "He": "He","Be": "Be","Ne": "Ne","Mg": "Mg","Sc": "Sc","Ti": "Ti","Cr": "Cr","Co": "Co","Ga": "Ga","As": "As","Kr": "Kr","Rb": "Rb","Sr": "Sr","Y": "Y","Zr": "Zr","Nb": "Nb","Mo": "Mo",
    "Tc": "Tc","Ru": "Ru","Rh": "Rh","Ag": "Ag","Cd": "Cd","In": "In","Sb": "Sb","Xe": "Xe","Cs": "Cs","Ba": "Ba","La": "La","Ce": "Ce","Pr": "Pr","Nd": "Nd","Pm": "Pm","Sm": "Sm",
    "Eu": "Eu","Gd": "Gd","Tb": "Tb","Dy": "Dy","Ho": "Ho","Er": "Er","Tm": "Tm","Yb": "Yb","Lu": "Lu","Hf": "Hf","Ta": "Ta","Re": "Re","Os": "Os","Ir": "Ir","Au": "Au","Hg": "Hg",
    "Tl": "Tl", "Bi": "Bi","Po": "Po","At": "At","Rn": "Rn","Fr": "Fr","Ra": "Ra","Ac": "Ac","Th": "Th","Pa": "Pa","U": "U","Np": "Np","Pu": "Pu","Am": "Am","Cm": "Cm","Bk": "Bk",
    "Cf": "Cf","Es": "Es","Fm": "Fm","Md": "Md","No": "No","Lr": "Lr","Rf": "Rf","Db": "Db","Sg": "Sg","Bh": "Bh","Hs": "Hs","Mt": "Mt","Ar": "Ar","B": "B_2","B": "B_3","N.": "N",
    "O.2": "O_R","O.": "O_1","O": "H_HB","O": "H_HB","O.": "O_3","C.3": "C_3"}
    
 newsym=[];
 for S in ffsymbols:
    if  S in sym.keys():
        newsym.append(sym[S])
    else:
        newsym.append(S)
  
 return newsym
 
 
def CreateLammps(mol, Paramname):
    import numpy as np
    Param = ReadParameters(Paramname)
    blist=[] #this tells you how many bonds an atom has
#    print(  "YTTY ", len(mol.get_bonds()))
    if len(mol.get_bonds())==0:
        Lbond, blist = GetBonds(mol.get_positions(), mol.cell, 15, mol.get_chemical_symbols(),  Param)
    else:
        Lbond = mol.get_bonds()
    
    if len(mol.get_ffsymbols())==[]:
        mol.set_ffsymbols(SpecifyAtoms(mol.get_chemical_symbols(), blist))
        
    mol.set_ffsymbols(Syb2Dre(mol.get_ffsymbols())) 
    ffsymbols=mol.get_ffsymbols()
    linkc=GetLinkMat(Lbond)
    Langle = GetAngles(linkc)
    Ldih = GetDihedrals(Lbond, linkc)
    Linv = GetInversions(linkc,ffsymbols)
    
    Bonds = Dreidingbndass(Lbond,ffsymbols,Param,type="harm")
    Angles = Dreidingangass(Langle,ffsymbols,Param)
    Dih = Dreidingdihass(Ldih,ffsymbols, linkc, Param)
    Inv=Dreidinginvass(Linv,ffsymbols,Param)
    
    atmtypes={}
    lab=1
    for S in set(ffsymbols):
        atmtypes[S]={'label':lab,'mass':Param[S]['Mass']}
        lab+=1

    # MAKE SURE COORDINATES FIT IN TE CELL
    print( " WARNING: cell is recalculated ")
    pos=mol.get_positions()
    mins=np.min(pos,axis=0)
    maxs=np.max(pos,axis=0)

    Str= "Created by the holy Pajton script\n"
    Str+= "\t%d atoms\n\t%d bonds\n\t%d angles\n\t%d dihedrals\n\t%d impropers\n" % (mol.get_number_of_atoms(), len(Lbond),len(Langle),len(Ldih),len(Linv))
    Str+= "\t%d atom types\n\t%d bond types\n\t%d angle types\n\t%d dihedral types\n\t%d improper types\n" % (len(atmtypes),len(Bonds),len(Angles),len(Dih),len(Inv))
    Str+= "\n\t%f %f xlo xhi\n" % (mins[0]-6.0,maxs[0]+6.0)
    Str+= "\t%f %f ylo yhi\n" % (mins[1]-6.0,maxs[1]+6.0)
    Str+= "\t%f %f zlo zhi\n" % (mins[2]-6.0,maxs[2]+6.0)
    Str+= "\n Masses\n\n"

    chemical_symbols = ""
    #print("Insert into .in file: dump_modify "),
    ksl = ["X"]*len(atmtypes)
    for M in atmtypes.keys():
        ksl[atmtypes[M]['label']-1] = M
    #for i in range(len(atmtypes.keys())):
        #print(atmtypes.keys())
        #M = atmtypes[str(i+1)]
    for ie, E in enumerate(ksl):
        #Str+= " %d %f #%s\n" % (atmtypes[M]['label'],atmtypes[M]['mass'],M)
        Str += " %d %f #%s\n" % (ie+1, atmtypes[E]['mass'], E)
        symbol = E.split("_")[0]
        #print(symbol),
        chemical_symbols += symbol + " "
    
    if len(Bonds)>0: 
        Str+= "\n Bond Coeffs\n\n"
        for B in  Bonds.keys():
            BT=Bonds[B]['params']
            Str+= " %d %f %f # %s\n" % (Bonds[B]['lab'],2*BT[0],BT[1],B)
    
    if len(Angles)>0: 
        Str+= "\n Angle Coeffs\n\n"
        for B in  Angles.keys():
            BT=Angles[B]['params']
            Str+= " %d %f %f # %s\n" % (Angles[B]['lab'],2*BT[0],BT[1],B) 
 
    if len(Dih)>0: 
        Str+= "\n Dihedral Coeffs\n\n"
        for B in  Dih.keys():
            BT=Dih[B]['params']
            #Str+= " %d %f %d %d # %s\n" % (Dih[B]['lab'],BT[0],np.cos(BT[1]/180.*np.pi),BT[2],B)
            Str+= " %d %f %d %d # %s\n" % (Dih[B]['lab'],BT[0],round(np.cos(BT[2]/180.*np.pi)),BT[1],B)
    
    if len(Inv)>0:
        Str+= "\n Improper Coeffs\n\n"
        for B in  Inv.keys():
            BT=Inv[B]['params']
            Str+= " %d %f 0 # %s\n" % (Inv[B]['lab'],BT[0],B)
    
    Str+= "\n Atoms\n\n"
    id=1
    ch=mol.get_charges()
    pos=mol.get_positions()
    for i in range(mol.get_number_of_atoms()):
        Str+= "%d 1 %d %f %f %f %f 0 0 0\n" % (id,atmtypes[ffsymbols[i]]['label'],ch[i],pos[i][0], pos[i][1], pos[i][2])
        id+=1
    
    if len(Lbond)>0: 
        Str+= "\n Bonds\n\n"
        id=1
        for B in  Bonds:
            for bid in Bonds[B]['id']:
                Str+= " %d %d %d %d\n" % (id, Bonds[B]['lab'],bid[0]+1, bid[1]+1)
                id+=1
    
    if len(Langle)>0: 
        Str+= "\n Angles\n\n"
        id=1
        id=1
        for B in  Angles:
            for aid in Angles[B]['id']:
                Str+= " %d %d %d %d %d\n" % (id, Angles[B]['lab'],aid[0]+1, aid[1]+1, aid[2]+1)
                id+=1

    if len(Ldih)>0: 
        Str+= "\n Dihedrals\n\n"
        id=1
        for B in  Dih:
            for did in Dih[B]['id']:
                Str+= " %d %d %d %d %d %d\n" % (id, Dih[B]['lab'],did[0]+1, did[1]+1, did[2]+1, did[3]+1)
                id+=1
    
    if len(Linv)>0:
        Str+= "\n Impropers\n\n"
        id=1
        for B in  Inv:
            for iid in Inv[B]['id']:
                Str+= " %d %d %d %d %d %d\n" % (id, Inv[B]['lab'],iid[0]+1, iid[1]+1, iid[2]+1, iid[3]+1)
                id+=1

    Str+= "\n Pair Coeffs\n\n"
    for K in atmtypes.keys():
        Str+= "%d %s %s\n" % (atmtypes[K]['label'],Param[K]['Epsilon'],float(Param[K]['Sigma'])/(2**(1/6.)))
    
    #Fdata=open("data", 'w')
    #Fdata.write(Str)
    #Fdata.close()

    return Str, chemical_symbols
    
#THIS IS ONLY NEEDED IF YOU CAN'T MAKE OPENBABEL OR OTHER SOFTWARES TO RECOGNIZE YOUR ATOM TYPES (LIKE A BIG FLAKE OF GRAPHENE)
def SpecifyAtoms(symbols, blist):
    ffsymbols=[]
    for si in range(len(symbols)):
        sym=symbols[si]
        if sym=="H":
            if blist[si]==1:
                sym=sym+"_"
        
        elif sym=="C":
            if blist[si]==1:
                sym=sym+"_1";
            elif blist[si]==2:
                sym=sym+"_2";
            elif blist[si]==3:
                sym=sym+"_R"; #OR SP2
            elif blist[si]==4:
                sym=sym+"_3";
            else:
                print( si, blist[si], sym)
            
                
        elif sym=="O":
            if blist[si]==1:
                sym=sym+"_2";
            if blist[si]==2:
                sym=sym+"_2";
            if blist[si]==3:
                sym=sym+"_2";
            if blist[si]==4:
                sym=sym+"_2";
            
        elif sym=="N":
            if blist[si]==1: sym=sym+"_1";
            if blist[si]==2: sym=sym+"_2";
            if blist[si]==3: sym=sym+"_2";
            if blist[si]==4: sym=sym+"_3";
            
        elif sym=="P":
            if blist[si]==2: sym=sym+"_2";
            if blist[si]==3: sym=sym+"_3";
            
        elif sym=="S":
            if blist[si]==1: sym=sym+"_3";
            if blist[si]==2: sym=sym+"_3";
            if blist[si]==3: sym=sym+"_3";
            if blist[si]==4: sym=sym+"_3";
            
        elif sym=="N":
            if blist[si]==2: sym=sym+"_1";
            if blist[si]==3: sym=sym+"_2";
            if blist[si]==4: sym=sym+"_3";
            
        elif sym=="Al":
            sym=sym+"3";
                

        ffsymbols.append(sym)
    
    return ffsymbols



def create_input(mol2_input, drei_input_temp, Paramfile, outdir):
    mol = io.read(mol2_input, format="mol2")
    data, chemical_symbols = CreateLammps(mol, Paramfile)

    with open(os.path.join(outdir, "data"), "w") as f:
        f.write(data)

    if type(drei_input_temp) == type("str"):
        lammps_temp = Template(drei_input_temp)
    else:
        lammps_temp = Template(open(drei_input_temp, 'r').read())

    with open(os.path.join(outdir, "in.drei"), "w") as f:
        f.write(lammps_temp.substitute(SYM=chemical_symbols))


def main():
    if len(sys.argv)!=5:
        print( "ERROR: Not enough arguments\n\t\t\t USAGE %s Molecule.mol2 inp_tempFile Dreiding_paramfile outputDir" %(sys.argv[0]))
        sys.exit()

    mol2_input = sys.argv[1]
    drei_input = sys.argv[2]
    Paramfile = sys.argv[3]
    outdir = sys.argv[4]

    create_input(mol2_input, drei_input, Paramfile, outdir)
    
    
if __name__ == "__main__":
    main()
