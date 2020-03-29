#!/usr/bin/python3.5
import sys
#import time
import getopt


def reindex(ii, n):

    if ii<0:
      ii+=n
    else:
      ii-=1

    return ii

def settype(var):
    try:
     var=int(var)
    except ValueError:
     try:
      var=float(var)
     except ValueError:
      print("ERROR with argument " + var)
      exit
      
    return var

def Tokenize(args,vtype,*Sep):
    if args:
     NewList = []
     if len(Sep)>0:
       FirstSep=Sep[0]
     if len(Sep)>1:
      SecondSep=Sep[1]
    
     List = str(args).split(FirstSep)
     if vtype=='str':
      return List
     
     if vtype=='range':
      List[0]=float(List[0])
      List[1]=float(List[1])
      List[2]=int(List[2])        
      return List
      
     if vtype=='vec':
      List[0]=float(List[0])
      List[1]=float(List[1])
      List[2]=float(List[2])        
      return List
    
     if vtype=='matrix':
      for L in List:
       NewList.append([float(r) for r in L.split(SecondSep)])
      return NewList

     for L in List:
       if len(Sep)>1:
         BE =L.split(SecondSep)
       if len(BE)==1:
         NewList.append(settype(BE[0]))
       if len(BE)==2:
         NewList.append(settype(BE[0]))
         NewList.append(settype(BE[1]))
     
     return NewList
      
   
#Facilitate interpreting input
def GetOpt(opts,param,default=None):
 Result=False
 for lab, arg in enumerate(opts):
  if arg[0] == param:   
   if arg[1]=="":
    Result=True
   else:
    Result = arg[1]
   del opts[lab]   
 return Result
                         
   
def GetInput():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "htdn:f:o:vl", \
        ["help", "output=","leads=","natoms=","erange=","scatter=","scissor=",\
        "anderson=","vgate=","efshift=","bias=","trange=","ham=","emout=","leadout=",\
        "bfield=","mode=","eftemp=","thermosweep=","restartthermo=","thermospinmeth=",\
        "thermoefrange=","thermousesavetrm=", "thermoef=", "thermotrange=","thermotemp=",\
        "efrange=","path_em=","path_leads=","hamiltonianprovider=","bias_accuracy=",\
        "verbose=","brshift=","biaspotvec=","gatepotvec=", "krange=","workingdir=", "extraleads="]) 
#    except getopt.GetoptError, err:
    except getopt.GetoptError:

    # print help information and exit:
        print(str(err)) # will print something like "option -a not recognized"
        print("type \"%s -h --help\" for further help" % (sys.argv[0]))
        sys.exit(2)
    if len(opts)==0:
     print("No options were provided!\nType \"%s -h --help\" for further help" % (sys.argv[0]))
     sys.exit(2)
       
    Cmds={'Mode':{'help':" --mode      (Explained in Scope in the manual) \n" + \
           "\t(0) Band structure calculation \n" + \
           "\t(1) Only sweep in energy (2) 1D, 2D conductance histograms" + \
           "\n\t(3) Zero voltage thermal properties (4) Finite Voltage TR\n\n"   
           , 'value':None,'vtype':'int'},
          'ERange':{'help':"emin emax numsteps", 'value':None,'vtype':'range'},
          'TRange':{'help':"tmin tmax numsteps ! ONLY IN MODE=3 !",'value':None,'vtype':'range'},
          'kRange':{'help':"kmin kmax numsteps ! ONLY IN MODE=0 !",'value':None,'vtype':'range'},
          'BField':{'help':"Bx By Bz", 'value':None,'vtype':'vec'},
          'ThermoEFRange':{'help':"efmin efmax numsteps ! ONLY IN MODE=3 !", 'value':None,'vtype':'range'},
          'ThermoTRange':{'help':"tmin tmax numsteps ! ONLY IN MODE=3 !", 'value':None,'vtype':'range'},
          'ThermoSweep':{'help':"", 'value':None,'vtype':'int'},
          'ThermoEF':{'help':"", 'value':None,'vtype':'int'},
          'ThermoTemp':{'help':"", 'value':None,'vtype':'int'},
          'ThermoUseSaveTrm':{'help':"", 'value':None,'vtype':'int'},
          'ThermoSpinMeth':{'help':"", 'value':None,'vtype':'int'},
          'Bias':{'help':"mnV maxV numsteps  last_atom_in_left_electrode last_atom_in_right_electrode ! ONLY IN MODE=4", 'value':None,'vtype':'matrix'},
          'HamiltonianProvider':{'help':"HamiltonianProvider dft or tbm", 'value':None,'vtype':'str'},
          'Path_EM':{'help':"Output file of dft calculation", 'value':None,'vtype':'str'},
          'Path_Leads':{'help':"Output file of dft calculation for leads", 'value':None,'vtype':'str'},
          'scissors':{'help':"0/1 numelectrons Ehomo Elumo screening_charge_distance", 'value':None,'vtype':'cust'},
          'leads':{'help':"first atom in the first layer, last atom in last layer, number of PL, terminating PL, ideal lead -1/0/1, rand -1/0/1, Bands 0/1", 'value':None,'vtype':'cust'},
          'Natoms':{'help':"Number of atoms in the system", 'value':None,'vtype':'int'},
          'Vgate':{'help':"", 'value':None,'vtype':'float'},
          'WorkingDir':{'help':"", 'value':None,'vtype':'str'},
          'ExtraLeads':{'help':"Use external leads?", 'value':None,'vtype':'str'},
          'Bias_accuracy':{'help':"", 'value':None,'vtype':'int'},
          'Verbose':{'help':"detail of output", 'value':0,'vtype':'int'},
          'Brshift':{'help':"Which atoms are effected by the bias shift", 'value':0,'vtype':'cust'},
          'BiasPotVec':{'help':"Value of bias shift on atoms", 'value':0,'vtype':'cust'},
          'GatePotVec':{'help':"Value of gating on atoms", 'value':0,'vtype':'cust'}
          }
   
    ""
    Help="\nHELP TO THE INPUTCREATOR FOR GOLLUM\n\n"
    for cmd in Cmds:
     Help+="\n\t%s\n\n"%(Cmds[cmd]['help'])
     optname=cmd.lower()
     if Cmds[cmd]['vtype']=='int' or Cmds[cmd]['vtype']=='float':
       Cmds[cmd]['value']=settype(GetOpt(opts,"--"+optname))
     elif Cmds[cmd]['vtype']=='str' :
       Cmds[cmd]['value']=Tokenize(GetOpt(opts,"--"+optname),'str',",")
     elif Cmds[cmd]['vtype']=='range' :
       Cmds[cmd]['value']=Tokenize(GetOpt(opts,"--"+optname),'range',",",":")
     elif Cmds[cmd]['vtype']=='vec' :
       Cmds[cmd]['value']=Tokenize(GetOpt(opts,"--"+optname),'vec',",",":")
     elif Cmds[cmd]['vtype']=='matrix' :
       Cmds[cmd]['value']=Tokenize(GetOpt(opts,"--"+optname),'matrix',";",",")
     
    if GetOpt(opts,"-h") or GetOpt(opts,"--help"):
        print(Help)  
    else:
        Pleads=[]
        NOP=GetOpt(opts,"--leads")
        if NOP:
            leads=NOP.split()
            for L in leads:
                T=Tokenize(L,'',",",":")
                Pleads.append(T)

        Cmds['leads']['value'] = Pleads
     
        Cmds['Brshift']['value'] = []
        NOP=GetOpt(opts,"--brshift")
        if NOP:
            brshift=[]
            for L in NOP.split():  
                T=Tokenize(L,'',",",":")
                brshift.append(T)
            Cmds['Brshift']['value'] = brshift
      
        Cmds['BiasPotVec']['value']=[]
        NOP=GetOpt(opts,"--biaspotvec")
        if NOP: 
            biasvec=[]
            for L in NOP.split():  
                T=Tokenize(L,'',",",":")
                biasvec.append(T)
            Cmds['BiasPotVec']['value'] = biasvec
     
        Cmds['GatePotVec']['value']=[]
        NOP=GetOpt(opts,"--gatepotvec")
        if NOP: 
            gatevec=[]
            for L in NOP.split():  
                T=Tokenize(L,'',",",":")
                gatevec.append(T)
     
            Cmds['GatePotVec']['value' ] = gatevec
    return Cmds   



#leadp num of PL, terminating PL, var -1/0/1, var -1/0/1, print bands 0/1
def GenLeadp(leads):

 strL=""
 numl=0
 for L in leads:
  numl+=1
  strL+="%d %d %d\n" % (L[2],L[3],L[4])
  
 return numl,strL


#/storage/users/manrique/developing/newtrans/bin/Debug2/trcodecb -Fg Extended_Molecule -emin -0.01 -emax 0.01 -n 2 -l 1:4 2 -l 7:11 

#Create input for fortran 
def CreateFGInp(Cmds):

    natoms=Cmds['Natoms']['value']

    keys=list(Cmds.keys())
    fginp=""

    if Cmds['Path_Leads']['value']:
        fginp+="# name: Path_Leads\n# type: string\n# rows: %d\n" %len(Cmds['Path_Leads']['value'])
        for i in range(len(Cmds['Path_Leads']['value'])):
            fginp+="%s\n"%Cmds['Path_Leads']['value'][i]
    keys.remove('Path_Leads')

    for cmd in keys:
        if Cmds[cmd]['value']!=None and cmd!='leads' and cmd!='Natoms':
  
            if Cmds[cmd]['vtype']=='int' or Cmds[cmd]['vtype']=='float':
                fginp+="# name: %s\n# type: scalar\n%s \n" %(cmd,str(Cmds[cmd]['value']))
            elif Cmds[cmd]['vtype']=='str':
                fginp+="# name: %s\n# type: string\n" %cmd
                fginp+="# rows: 1\n%s\n"%Cmds[cmd]['value'][0]
            elif Cmds[cmd]['vtype']=='vec' or Cmds[cmd]['vtype']=='range':
                fginp+="""# name: %s
# type: matrix
# rows: 1
# columns: %d
"""%(cmd,len(Cmds[cmd]['value']))
                for i in range(len(Cmds[cmd]['value'])):
                    fginp+="%s "%str(Cmds[cmd]['value'][i])   
                fginp+="\n"
            elif Cmds[cmd]['vtype']=='matrix':
                fginp+="""# name: %s
# type: matrix
# rows: %d
# columns: %d
"""%(cmd,len(Cmds[cmd]['value']),len(Cmds[cmd]['value'][0]))
                for i in range(len(Cmds[cmd]['value'])):
                    for j in range(len(Cmds[cmd]['value'][i])):
                        fginp+="%s "%str(Cmds[cmd]['value'][i][j])   
                    fginp+="\n"


    numl,strL=GenLeadp(Cmds['leads']['value'])
    fginp+="""# name: leadp
# type: matrix
# rows: %d
# columns: 3
%s""" %(numl,strL)

    fginp+="""# name: atom
# type: matrix
# rows: %d
# columns: 6
""" % (natoms)

    atoms = Gen_atom(natoms,Cmds['leads']['value'])
    natoms = len(atoms)
    
    if Cmds['Brshift']['value']!=[]:
        for il in range(numl):
            B = Cmds['Brshift']['value'][il]
            for i in range(reindex(B[0],natoms),reindex(B[1],natoms)+1):
                atoms[i][2] = il+1

    if Cmds['BiasPotVec']['value']!=[]:
        for i in range(natoms):
            atoms[i][3] = Cmds['BiasPotVec']['value'][0][i]

    for i in range(natoms):
        fginp+="%d %d %d %d %4.2f %f\n" % (i+1,atoms[i][0],atoms[i][1],atoms[i][2],atoms[i][3],0)
 
    F=open("input","w")
    F.write(fginp)
    F.close()


def CreateFGInp_octave(Cmds):
    keys=list(Cmds.keys())
    oct_script=""
    
    
    oct_script+="inputP.HamiltonianProvider={%s};\n"%(Cmds['HamiltonianProvider']['value'])
    keys.remove('HamiltonianProvider')
    
    oct_script+="inputP.Path_EM={%s};\n"%(Cmds['Path_EM']['value'])
    keys.remove('Path_EM')
    
    oct_script+="inputP.Path_Leads=["
    if Cmds['Path_Leads']['value']:
        for il in range(len(Cmds['Path_Leads']['value'])):
            numL,  lstr = Cmds['Path_Leads']['value'][il].split()
            oct_script+="{'%s'};"%(Cmds['Path_Leads']['value'][il])
        oct_script+="];\n"
    keys.remove('Path_Leads')
    
    for cmd in keys:
        if Cmds[cmd]['value']!=[] and Cmds[cmd]['value']!=None:
            if Cmds[cmd]['vtype']=='matrix':
                oct_script+="inputP.%s=["%cmd
                for V in Cmds[cmd]['value']:
                    oct_script+="%s;"%(V)
                oct_script+='];\n' 
            else:
                oct_script+="inputP.%s=%s;\n"%(cmd,Cmds[cmd]['value'])
      

    oct_script+="inputP.leadp=zeros(%d,%d);\n"%(len(Cmds['leads']['value']), 3)
    n=0
    for L in Cmds['leads']['value']:
        n+=1
        oct_script+="inputP.leadp(%d,:)=[%4d %4d %4d];\n"% (n, L[2],L[3],L[4])
    
    oct_script+="inputP.atom=zeros(%d, %d);\n"%(Cmds['Natoms']['value'], 6)
    atoms = Gen_atom(Cmds['Natoms']['value'], Cmds['leads']['value'])
    natoms = len(atoms)
    for i in range(natoms):
        oct_script+="inputP.atom(%d,:)=[ %4d %4d %4d %4d %4d %4d];\n"% (i+1, i+1,atoms[i][0],atoms[i][1],0,0,0)
    
    #oct_script+="save('input.txt')\n"
    #oct_script+="save('input.mat','-mat')\n"
     
    F=open("gollum_input.m","w")
    F.write(oct_script)
    F.close()


def Gen_atom(natoms,leads):
    #create atom vector
    atoms = [[0,0,0,0,0] for i in range(1,natoms+1)]  
    
    cutmin=10000
    cutmax=0
    #Insert leads
    labL=1
    for Lead in leads:
        #  FA=Lead[0]-1
        FA = reindex(Lead[0], natoms)
        #  LA=Lead[1]-1
        LA = reindex(Lead[1], natoms)
        NA=(abs(LA-FA)+1)%float(Lead[2])
        cutmin=min(cutmin, FA)
        cutmin=min(cutmin, LA)
        cutmax=max(cutmax, FA)
        cutmax=max(cutmax, LA)
        
        if NA>0.00001:
            print( "ERROR with Lead layers: NA = ",NA)
            print( "The %d lead has uneven number of atoms" % (labL))
            sys.exit()
        else:
            PL=int(Lead[2])
            NA=int((abs(LA-FA)+1)/PL)

#  if FA<0:
#   FA+=natoms+1
#  if LA<0:
#   LA+=natoms+1
  
        lab=1
        if LA<FA:
           ld=-1
        else:
           ld=1

        for labpl in range(PL):
            for i in range(NA):
                atoms[FA+ld*(i+labpl*NA)][:2]=[labL,labpl+1]
            lab+=1

        labL+=1
    
    return atoms[cutmin:cutmax+1]
    #return atoms
    
 

if __name__ == '__main__':
    Cmds=GetInput()
    CreateFGInp(Cmds)
    CreateFGInp_octave(Cmds)
