#!/usr/bin/python

import math
import cmath
import numpy as np
import os,sys
import shutil
from numpy import *
import scipy.linalg  as scl
import copy
import sys
import time



def reindex(ii, jj, n):
    if ii<0:
        ii = ii+n
    else:
        ii = ii-1
    if jj<0:
        jj = jj+n
    else:
        jj = jj-1
            
    if(ii==jj) :
        print ("Empty index range!")

    if(ii>jj):
        ii, jj = jj, ii 
    
    return ii, jj
 
        
def CalcBand(fname,kp, Emin, Emax):
    Ham = ReadOctave(fname)
    ik=kp
    e=Emin
    
    iorb = np.array(Ham['iorb'], dtype='int16')
    iorb[:, 0]=iorb[:, 0]-1   #octave,matlab conversion
    
    #Reading in one
    dim=len(iorb)
    ef=float(Ham['FermiE'])

    #sh.n = iorb.max(0)[0]+1 
    n = iorb.max(0)[0]+1 
    
    #atom_orbitals = [len(np.nonzero(iorb[:, 0] == I)[0]) for I in range(sh.n)]
    H, S =     ReadHS(fname, dim, ik)
    kp=Ham['kpoints_EM'][ik][:2]
    #print kp
    #print "Energy ", e+ef
    #Hk=H+H*np.exp(1.0j*kp)+H.H*np.exp(-1.0j*kp)
    #Sk=S+S*np.exp(1.0j*kp)+S.H*np.exp(-1.0j*kp)
    Hk=H
    Sk=S
    K=Hk
    Ev,vl=np.linalg.eig(K)
    Fev=open("sav-"+str(ik), 'w')
    for E in Ev:
        Fev.write('%6.3f %6.3f %-8.4f\n' % (kp[0], kp[1],E-ef))
    Fev.close()
    
    Etry=-1000
    vk=0
    for i in range(len(Ev)):
        if Etry-Ev[i]>0:
            vk=vl[i]
    
    return Dos(Hk, Sk, vk, ef)
    
    

def Dos(H, S, vk, Ef):

    Ginv = H-(Ef-1)*S
    G = np.linalg.inv(Ginv)
    dos = -imag(np.diag(G))
    return dos
    
    

#----------------------------------------------------------
#UTILITIES
#----------------------------------------------------------

def eikr(H, S, kp, xij):
    
    kxij=kp*xij
    eikr = np.exp(1j*kxij)
    Hk=H*eikr
    Sk=S*eikr
    return Hk, Sk


def vkk2(fname, ikp):
    Ham = ReadOctave(fname)
   
    iorb = np.array(Ham['iorb'], dtype='int16')
    iorb[:, 0]=iorb[:, 0]-1   #octave,matlab conversion
    
    #Reading in one
    dim=len(iorb)
    ef=float(Ham['FermiE'])

    natom = iorb.max(0)[0]+1 
    niorb = iorb.shape[0]
    
    nkpnt=Ham['kpoints_EM'].shape[0]
    #Evs=np.zeros((nkpnt, niorb))
    #Evecs=np.zeros((nkpnt, niorb, niorb))
    H, S =     ReadHS(fname, dim, ikp)
    print(Ham['kpoints_EM'][ikp])
    #Hk, Sk = eikr(H, S, Ham['kpoints_EM'][ikp][1:], Ham['iorb'][:, 5:8])
    tev,tevec = scl.eig(H,S)
    #print(np.sum((H-H.H )))
        
    Evs = tev

    Evecs = np.absolute(tevec)**2        
    
    return  Evs, Evecs



def vkk(fname):
    if fname[-3:]=='mat':
        Ham = ReadMatlab(fname)
        ind=array(Ham['HSM'][:,0],dtype=int32)
    else:
        Ham = ReadOctave(fname)
    
    
    
    iorb = np.array(Ham['iorb'], dtype='int16')
    iorb[:, 0]=iorb[:, 0]-1   #octave,matlab conversion
    
    #Reading in one
    dim=len(iorb)
    ef=float(Ham['FermiE'])

    natom = iorb.max(0)[0]+1 
    niorb = iorb.shape[0]
    
    nkpnt=Ham['kpoints_EM'].shape[0]
    Evs=np.zeros((nkpnt, niorb))
    Evecs=np.zeros((nkpnt, niorb, niorb))
    Evecs2=np.zeros((nkpnt, niorb, niorb))
    for i in range(Ham['kpoints_EM'].shape[0]):
        if fname[-3:]=='mat':
            H=zeros((dim, dim), dtype=complex)
            S=zeros((dim, dim), dtype=complex)
            kHam = Ham['HSM'][where(ind[:]==i+1)]
            ii = array(kHam[:, 1]-1, dtype=int)
            jj = array(kHam[:, 2]-1, dtype=int)
            H[ii, jj]= kHam[:, 5] + 1.0j*kHam[:,6]
            S[ii, jj] = kHam[:,  3] + 1.0j*kHam[:,4]
        else:
            H, S =     ReadHS(fname, dim, i)
            
        tev,tevec = scl.eig(H,S)
        #print(np.sum((H-H.H )))
        
        Evs[i] = tev

        Evecs[i] = np.absolute(tevec)**2        
        tevec2=tevec
        for j in range(tevec.shape[0]):
            tevec2[:,j]=tevec[:,j]/diag(S)
        Evecs2[i] = tevec2
    
    return  Evs, Evecs, Evecs2

# LOAD EVERYTHING AT ONCE
def LoadData(dirs="./",  threshold=0.5,  atomorbs=[]):
    allSavok={};allBands={};allIorb={};allEf={};allEvs={}
    
    for D in os.listdir(dirs):
        if os.path.isfile(base+D+"/Eigs.npy"):
            print(D)
            Evs,SiestaBands,Evecs,iorb,Ef = Readdata(base+D)
    
            allIorb[D] = iorb
            allEf[D]=Ef
            allEvs[D]=Evs
            usav=zeros(Evs.shape)

            for i in range(Evs.shape[0]):
                usav[i] = SiestaBands[i::Evs.shape[0],1] 
            allSavok[D]=usav
    
            #MELYIK PALYAKAT NEZZUK
            #BITEI: AHOL A RENDSZAM NAGYOBB MINT 6
            Orbs={}
            Orbs['BiTeI']=where(allIorb[D][:,8]>6)
            #Orbs['BN']=where(iorb[:,9]<12)
            Orbs['graphene']=where(allIorb[D][:,8]==6)
            #Orbs['bngr']=np.where(iorb[:,0]>3)[0]

            #Sort out the bands
            #MENNYIRE LEGYEN LOKALIZALT A HFV AZ ADOTT ATOMOKON?
            threshold = 0.6
    
            Bands={}
            lab=0
            for key in Orbs.keys():
                Bands[key]={}
                for ti in range(Evs.shape[1]):
                    Bands[key][ti]=[]

                ind=np.ones(Evs.shape[1], dtype=bool)

                #LOOP k points
                for j in range(Evs.shape[0]):                                
                    #print(Orbs[key][0])
                    ind = where(sum(Evecs[j,Orbs[key][0],:],axis=0)>threshold)
    
                    #LOOP evs
                    for ti in ind[0]:            
                        Bands[key][ti].append([j/Evs.shape[0],Evs[j,ti]])        
        
                allBands[D]=Bands
                
                lab+=1
                
    return  allEvs, allSavok, allBands, allIorb,  allEf

def Readdata(bdir):
    Ham = ReadOctave(bdir+'/Extended_Molecule')
    iorb = array(Ham['iorb'], dtype='int16')
    Ef = float(Ham['FermiE'])
    Evs = load(bdir+"/Eigs.npy")
    Evecs = load(bdir+"/Evecs.npy")
    Sav = loadtxt(bdir+"/savok")
    del Ham
    
    return Evs,Sav,Evecs,iorb,Ef

def gsdfg():    
    Etry=-1000
    vk=0 
    Fvk=open("sajt-"+str(ik), 'w')
    inull=0
    for i in range(len(Ev)):
        if Etry-Ev[i]>0:
            inull=i
            
        Etry=Ev[i]
    
    for ie in range(inull-3, inull+3):
        for v in vl:
            Fvk.write('%6d %6.3f '%(ie, np.abs(v)))
        
        Fvk.write('\n')
        
    Fvk.close()
        
    Fev.close()
    
    
    

    



#////////////
def ReadHS(fname, dim,  kp):
    File = open(fname,"r")
    T = {}
    Line1 = File.readline().rstrip().split( )
    S = np.matrix(np.zeros((dim, dim)), dtype='complex64')
    H = np.matrix(np.zeros((dim, dim)), dtype='complex64')
    while len(Line1)>0:
        if (Line1[0]=="#"):
            if ((Line1[2]=="HSM") or (Line1[2]=="HS")):
                rows = File.readline().rstrip().split()
                rows = int(File.readline().rstrip().split()[2])
                cols = int(File.readline().rstrip().split()[2])
                break
        Line1 = File.readline().rstrip().split()
    
    L = File.readline().rstrip().split()
    try:
        while len(L)>0:                
            if int(L[0])==kp+1:
#                for i in range(rows):
                while int(L[0])==kp+1:
                    
                    ii = int(L[1])-1
                    jj = int(L[2])-1
                    H[ii,jj] = float(L[5])+1.0j*float(L[6]);
                    S[ii,jj] = float(L[3])+1.0j*float(L[4]); 
                    L = File.readline().rstrip().split()
                break
            L = File.readline().rstrip().split()
    except IndexError:
        print (L)
        
#    print H[0, 0]
    File.close()
    #return copy.deepcopy(H), copy.deepcopy(S)
    return H, S
    
def ReadMatlab(fname):
    from scipy.io import loadmat
    
    extmol = loadmat(fname)
    
    return extmol
    
def ReadOctave(fname):
    File = open(fname,"r")
    T = {}
    Line1 = File.readline().rstrip().split( )
    while len(Line1)>0:
        if Line1[0]=="#":
            if (Line1[2]!="HSM"):
                Line2 = File.readline().rstrip().split()
                if Line2[2]=="scalar":
                    T[Line1[2]] = File.readline().rstrip().split()[0]
                if Line2[2]=="matrix":
                    rows = int(File.readline().rstrip().split()[2])
                    cols = int(File.readline().rstrip().split()[2])
                    Aux = []
                    for i in range(rows):
                        Line = File.readline().rstrip().split()
                        Aux.append(np.mat(Line, dtype='float64'))
                        
    
                    T[Line1[2]] = np.reshape(Aux, (rows, cols))
                    
        Line1 = File.readline().rstrip().split()
    File.close()

    return T

def ReadHHH(fname):
    Data = open(fname,"r").readlines()
    dim=int(Data[-1].split()[0])
    S = np.matrix(np.zeros((dim, dim)), dtype='complex64')
    H = np.matrix(np.zeros((dim, dim)), dtype='complex64')
    for L in Data:
        D=L.split()
        i=int(D[0])-1
        j=int(D[1])-1
        S[i, j]=float(D[2])+1.0j*float(D[3])
        H[i, j]=float(D[4])+1.0j*float(D[5])
    return H, S

def ReadHaux(fname):
    Data = open(fname,"r").readlines()
    dim=max(int(Data[-1].split()[0])*2, int(Data[-1].split()[1])*2)
    print(dim)
    S = np.matrix(np.zeros((dim, dim)), dtype='complex64')
    H = np.matrix(np.zeros((dim, dim)), dtype='complex64')
    for L in Data:
        D=L.split()
        i=(int(D[0])-1)*2
        j=(int(D[1])-1)*2
        S[i, j]=float(D[2])+1.0j*float(D[3])
        S[i+1, j+1]=float(D[4])+1.0j*float(D[5])
        for si in range(2):
            for sj in range(2):
                H[i+si, j+sj]=float(D[6+4*si+2*sj])+1.0j*float(D[7+4*si+2*sj])
    return H, S


def Get_Orb_data(dirname,El): 
    Orb_dat={}
    for E in El:
        Orb_dat[E]={}
        fn=E+".ion"
        orbdat=open(dirname+"/"+fn,'r').readlines()
    for LN in range(len(orbdat)):
        if orbdat[LN].find("#orbital l, n, z")>-1:
            line=orbdat[LN].split()
            lnzp=line[0]+line[1]+line[2]+line[3]
            OD=[]  
            LN+=2
        while orbdat[LN].find("#")==-1:
            OD.append((float(orbdat[LN].split()[0]),float(orbdat[LN].split()[1])))
            LN+=1
        LN-=1
    #FOR NOW
#    if (line[2]=="1" and line[0]=="0"):
        OD=np.array(OD)
        OD.resize((len(OD),2))
    #NOW WE CREATE A FUNCTION, WHICH WILL INTERPOLATE VALUE TO THE GIVEN DISTANCE
        Orb_dat[E][lnzp]=interp1d(OD[:,0],OD[:,1],bounds_error=False,fill_value=0.0,copy=False)

    return Orb_dat


if __name__== "__main__":
    fname=sys.argv[1]
    kp=int(sys.argv[2])
    Emin, Emax=float(sys.argv[3]), float(sys.argv[4])
    vkk(fname,kp, Emin, Emax)
