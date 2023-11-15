#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 12:22:44 2023

@author: Kevin Jiguet-Covex
"""


import sys, os.path
import crystallography as cr
from functools import partial
import concurrent.futures
from os.path import join
import ctypes
import numpy as np
import time
import platform

OS = platform.system()
LibraryName =""

if OS == "Linux":
    LibraryName = 'c_UMDprocess.so'
elif OS == "Windows":
    LibraryName = 'c_UMDprocess.dll'
elif OS == "Darwin":
    LibraryName = 'c_UMDprocess.dylib'




current_path=os.path.abspath(__file__)
path_split=current_path.split('/')
path_red=path_split[1:-1]
path_new=''
for u in path_red:
    path_new+='/'+u



read_lib = ctypes.cdll.LoadLibrary(join(path_new, LibraryName))
read_lib.defineBonds.argtypes = [ctypes.POINTER(ctypes.c_char),ctypes.c_int,ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_int),ctypes.c_int,ctypes.c_int,ctypes.c_int]
read_lib.defineBonds.restype = ctypes.POINTER(ctypes.c_int)
read_lib.read_umd_values.argtypes = [ctypes.POINTER(ctypes.c_char),ctypes.c_int,ctypes.c_int]
read_lib.read_umd_values.restype = ctypes.POINTER(ctypes.c_double)
read_lib.free_memory.argtypes = [ctypes.POINTER(ctypes.c_double)]
read_lib.free_memory.restype = None


def read_snapshot_values_C(octettop,octetbot,step,natom,File,X,mode,maxStep):#Extracts the coordinates from a snapshot 
                                                                #whose data is contained between the octets "octettop" and "octetbot" in the umd file (thus the prep_read_coord function)
                                                                #The parameter X is the first index of the relevant 3 coordinate on each line in the umd file
    if(100*step//maxStep!=100*(step-1)//maxStep):
        sys.stdout.write("\rReading coordinates. Progress : "+str(100*step//maxStep)+"%\t")
        sys.stdout.flush()
    ff=open(File,"r")
    ff.seek(octettop,0)
    snapshot=ff.read(octetbot-octettop)
    ff.close()
    snapPointer = ctypes.c_char_p(snapshot.encode('utf-8'))
    CList=read_lib.read_umd_values(snapPointer,natom,X)
    SnapshotValues=[]
    if X==-1:
        k=12
    else :
        k=3
    #Converts the C data into a python list
    if mode=="lists":#Returns a list of 3-elements list (x,y,z) for each atom
        for i in range(natom):
            List=[]
            for j in range(k):
                List.append(round(CList[k*i+j],5))
            SnapshotValues.append(np.array(List))
    elif mode=="line":#Returns a unique list of values for all the atoms, atom after atom
        for i in range(k*natom):
            SnapshotValues.append(round(CList[i],5))
    elif mode=="chunks":#Returns a unique list of values for all the atoms, axis after axis
        for j in range(k):
            for i in range(natom):
                SnapshotValues.append(round(CList[k*i+j],5))
    
    read_lib.free_memory(CList)
  #  if step==300 :
#        print(SnapshotValues[41:45])
 #       sys.exit()
        
    return SnapshotValues

def read_snapshotbonds_C(octettop,octetbot,CentIndexes,AdjIndexes,File,natom):#Converts the bonding information of a given snapshot into a list of atoms 
                                                                               #the data is contained between the octets "octettop" and "octetbot" in the bonding file
    ff=open(File,"r")
    ff.seek(octettop,0)
    snapshot=ff.read(octetbot-octettop)
    ff.close()
    
    snapPointer = ctypes.c_char_p(snapshot.encode('utf-8'))
    CIp = (ctypes.c_int * len(CentIndexes))(*CentIndexes)
    AIp = (ctypes.c_int * len(AdjIndexes))(*AdjIndexes)

    
    BList = read_lib.defineBonds(snapPointer,len(snapshot),CIp,AIp,int(len(CentIndexes)/2),int(len(AdjIndexes)/2),natom)
    SnapshotBondingTable = []
    SnapshotBondIndexes = []
   
    #Converts the C data into two python lists
    for i in range(2,BList[0]+1):
        SnapshotBondingTable.append(BList[i])#This one contains the atoms
    for i in range(BList[0]+2,BList[0]+BList[1]+2):
        SnapshotBondIndexes.append(BList[i])#This one contains the information about which is bound to which in the other list
    return SnapshotBondingTable,SnapshotBondIndexes

def Octets_from_File(File,key,nlines):#Returns the indexes of octets pertaining to each snapshot in the file.
    ff=open(File,"r")                 #"key" is the pattern that makrs the separation between the successive snapshots in the file
    OctIndexes = []                   #"nlines" is the number of lines (= atoms) in each snapshot

    print("Preparing the file ",File," to be read...")
    while True :
        line = ff.readline()
        if not line :
            break
        line=line.strip().split()
        if len(line)>0 and line[0]==key:
            OctIndexes.append(ff.tell())
            for _ in range(nlines):
                ff.readline()
    OctIndexes.append(ff.tell())
    ff.close()
    return OctIndexes

def Crystallization(File,readlength = ""):#builds the Crystal
    MyCrystal = cr.Lattice()
    ff=open(File,"r")
    while True:
        line = ff.readline()
        if not line: break
            #print(line,len(line))
        if len(line) > 1:
            line=line.strip()
            entry=line.split()
            if entry[0] == 'natom':
                MyCrystal.natom = int(entry[1])
                MyCrystal.typat = [0 for _ in range(MyCrystal.natom)]
            elif entry[0] == 'ntypat':
                MyCrystal.ntypat = int(entry[1])
                MyCrystal.types = [0 for _ in range(MyCrystal.ntypat)]
                MyCrystal.elements = ['X' for _ in range(MyCrystal.ntypat)]
                MyCrystal.masses = [1.0 for _ in range(MyCrystal.ntypat)]
                MyCrystal.zelec = [1.0 for _ in range(MyCrystal.ntypat)]
            elif entry[0] == 'types':
                for ii in range(MyCrystal.ntypat):
                    MyCrystal.types[ii]=int(entry[ii+1])
            elif entry[0] == 'elements':
                for ii in range(MyCrystal.ntypat):
                    MyCrystal.elements[ii]=entry[ii+1]
            elif entry[0] == 'masses':
                for ii in range(MyCrystal.ntypat):
                    MyCrystal.masses[ii]=float(entry[ii+1])
            elif entry[0] == 'Zelectrons':
                for ii in range(MyCrystal.ntypat):
                    MyCrystal.zelec[ii]=float(entry[ii+1])
            elif entry[0] == 'typat':
                for ii in range(MyCrystal.natom):
                    MyCrystal.typat[ii] = int(entry[ii+1])
            elif entry[0] == 'timestep':
                TimeStep = float(entry[1])
            elif entry[0] == 'time':
                TimeInit = float(entry[1])
            elif entry[0] == 'acell':
                MyCrystal.acell=[float(entry[1]),float(entry[2]),float(entry[3])]
            elif entry[0] == 'rprim_a':
                MyCrystal.rprim[0]=[float(entry[1]),float(entry[2]),float(entry[3])]
            elif entry[0] == 'rprim_b':
                MyCrystal.rprim[1]=[float(entry[1]),float(entry[2]),float(entry[3])]
            elif entry[0] == 'rprim_c':
                MyCrystal.rprim[2]=[float(entry[1]),float(entry[2]),float(entry[3])]
            elif entry[0] == 'rprimd_a':
                MyCrystal.rprimd[0]=[float(entry[1]),float(entry[2]),float(entry[3])]
            elif entry[0] == 'rprimd_b':
                MyCrystal.rprimd[1]=[float(entry[1]),float(entry[2]),float(entry[3])]
            elif entry[0] == 'rprimd_c':
                MyCrystal.rprimd[2]=[float(entry[1]),float(entry[2]),float(entry[3])]
            if entry[0] == 'atoms:':
                break
    ff.close()
    
    if readlength == "WithLength" :
    
        ff = open(File,"rb")
    
        ff.seek(0,2)-1
        octet = ff.tell()
    
        char = ff.read(1)
    
        while char in [b'0',b'1',b'2',b'3',b'4',b'5',b'6',b'7',b'8',b'9',b'\t',b'.',b'\n',b" ",b"-",b"e",b"+",b'']:
            octet -= 1
            ff.seek(octet)
            char = ff.read(1)

        timeBool = [0,0,0,0,0]
        timeName = [b't',b'i',b'm',b'e',b' ']
        index = 4
    
        while 0 in timeBool :
            ff.seek(octet)
            char = ff.read(1)
            octet -= 1
            if char in timeName and index==timeName.index(char):
            
                timeBool[index]=1
                index -=1

            else :
                index = 4
                timeBool=[0,0,0,0,0]
    
        octet += 6
        time = 0
        decimal = 0
        ff.seek(octet)
        char = ff.read(1)
        while char in [b'0',b'1',b'2',b'3',b'4',b'5',b'6',b'7',b'8',b'9']:
            time*=10
            time+=int(char)
            octet +=1 
            ff.seek(octet)
            char = ff.read(1)

        if char == b'.':
            exp=0
            octet +=1 
            ff.seek(octet)
            char = ff.read(1)
            while char in [b'0',b'1',b'2',b'3',b'4',b'5',b'6',b'7',b'8',b'9']:
                decimal*=10
                decimal+=int(char)
                exp+=1
                octet +=1 
                ff.seek(octet)
                char = ff.read(1)
            decimal/=10**exp

        time += decimal
    
        ff.close()
    
        length = int((time-TimeInit)/TimeStep)+1
    
        return MyCrystal, TimeStep, length

    else :

        return MyCrystal, TimeStep

def read_values(UMDfile,key,mode="line",Nsteps=1,firststep = 0,laststep = None,cutoff="all",nCores = None):
    
    t=time.time()
    MyCrystal,TimeStep = Crystallization(UMDfile)
    OctIndexesRaw = Octets_from_File(UMDfile, 'atoms:', MyCrystal.natom)
    
    if cutoff=="all":
        cutoff = len(OctIndexesRaw)-1
    
    OctIndexes=OctIndexesRaw[:cutoff+1]

    if laststep == None :
        OctIndexes = OctIndexes[firststep:]
    else :
        OctIndexes = OctIndexes[firststep:laststep+1]
    
    OctTop = [OctIndexes[i*Nsteps] for i in range(int((len(OctIndexes)-1)/Nsteps))]#Indicates the beginning of each relevant snapshot
    OctBot = [OctIndexes[i*Nsteps+1] for i in range(int((len(OctIndexes)-1)/Nsteps))]#Indicates it's end

    if key == "xred":
        X=0
        print('Extracting the reduced coordinates from the file '+UMDfile+'... ')
 
    elif key == "xcart":
        X=3
        print('Extracting the cartesian coordinates from the file '+UMDfile+'... ')

    elif key == "absxcart":
        X=6
        print('Extracting the absolute coordinates from the file '+UMDfile+'... ')

    elif key == "velocity":
        X=9
        print('Extracting the velocities from the file '+UMDfile+'... ')
    elif key == "everything":
        X=-1

    else :
        print("Parameter < ",key," > not recognized")
        sys.exit()

    if mode!="line" and mode !="lists" and mode!="chunks":
        print("Parameter < ",mode," > not recognized")
        sys.exit()
    
        
    
    readvalues = partial(read_snapshot_values_C,natom=MyCrystal.natom,File=UMDfile,X=X,mode=mode,maxStep=len(OctTop)-1)
    if nCores != None :
        with concurrent.futures.ProcessPoolExecutor(max_workers = nCores) as executor :
            SnapshotsValuesList = list(executor.map(readvalues,OctTop,OctBot,[step for step in range(len(OctTop))]))        
    else:
        with concurrent.futures.ProcessPoolExecutor() as executor :
            SnapshotsValuesList = list(executor.map(readvalues,OctTop,OctBot,[step for step in range(len(OctTop))]))        


    print("Values extracted in ",(time.time()-t)," s")
    return MyCrystal, SnapshotsValuesList, TimeStep, len(OctIndexes)-1


def data_type(SnapshotsValuesList,natom,mode="line",datatype="list"):

    if mode=="atoms":
        SnapshotsValuesX = [[[x[3*i],x[3*i+1],x[3*i+2]]  for i in range(natom)] for x in SnapshotsValuesList]
    else :
        SnapshotsValuesX = SnapshotsValuesList

    if datatype=="array":
            SnapshotsValues = [np.array(x) for x in SnapshotsValuesX]
    else :
            SnapshotsValues = SnapshotsValuesX
    
    return SnapshotsValues

def read_bonds(BondFile,centEls,adjEls,Nsteps=1,mode="Indexes",nCores = None):
    
    MyCrystal,TimeStep = Crystallization(BondFile)

    if centEls == ['all']:
        centEls = MyCrystal.elements
    if adjEls == ['all']:
        adjEls = MyCrystal.elements

    CentElWitness = [False for _ in range(len(centEls))]
    AdjElWitness = [False for _ in range(len(adjEls))]

        
    AdjacentElements = []
    CentralElements = []


    for i in range(len(MyCrystal.elements)):
        el = MyCrystal.elements[i]
        trueEl = ""
        for ch in el :
            if not ch.isdigit():
                trueEl+=ch
        if trueEl in centEls :
            CentralElements.append(el)
            CentElWitness[centEls.index(trueEl)]=True

        elif el in centEls :
            CentralElements.append(el)
            CentElWitness[centEls.index(el)]=True

        if trueEl in adjEls:
            AdjacentElements.append(el)
            AdjElWitness[adjEls.index(trueEl)]=True
        elif el in adjEls:
            AdjacentElements.append(el)
            AdjElWitness[adjEls.index(el)]=True
            
    CentIndexes = []
    AdjIndexes = []
    
    for el in CentralElements :
        indCent = MyCrystal.elements.index(el)
        CentMin = sum(MyCrystal.types[:indCent])
        CentMax = CentMin + MyCrystal.types[indCent] -1
        CentIndexes+=[CentMin,CentMax]
        

    for el in AdjacentElements :
        indAdj = MyCrystal.elements.index(el)
        AdjMin = sum(MyCrystal.types[:indAdj])
        AdjMax = AdjMin + MyCrystal.types[indAdj] -1
        AdjIndexes+=[AdjMin,AdjMax]
        
   
#    if centEls == ["all"]:
#        CentIndexes = [0,MyCrystal.natom-1]
#        CentElWitness = [True]
#    if adjEls == ["all"]:
#        AdjIndexes = [0,MyCrystal.natom-1]
#        AdjElWitness = [True]
        
    noAt = []
    for wit in range(len(CentElWitness)) :
        if CentElWitness[wit] == False :
            noAt.append(centEls[wit])

    for wit in range(len(AdjElWitness)) :
        if AdjElWitness[wit] == False :
            noAt.append(adjEls[wit])

    if noAt != []:
        if len(noAt) ==1 :
            print("ERROR : element ",noAt[0]," not present in simulation")
        else :
            print("ERROR : elements ",noAt," not presents in simulation")
        return -1,noAt,0,[0,0],0
        

#    print(CentMin,CentMax,AdjMin,AdjMax,indCent,indAdj)

    OctIndexes = Octets_from_File(BondFile,"step",MyCrystal.natom)

    OctTop = [OctIndexes[i*Nsteps] for i in range(int((len(OctIndexes)-1)/Nsteps))]#Indicates the beginning of each relevant snapshot
    OctBot = [OctIndexes[i*Nsteps+1] for i in range(int((len(OctIndexes)-1)/Nsteps))]#Indicates it's end

    print("Extracting bonds from file ",BondFile)
    
    bondsRed = partial(read_snapshotbonds_C,CentIndexes=CentIndexes,AdjIndexes=AdjIndexes,File=BondFile,natom=MyCrystal.natom)

#    Data=[]
#    for i in range(len(OctTop)):
#        Data.append(bondsRed(OctTop[i],OctBot[i]))            
    if nCores != None :
        with concurrent.futures.ProcessPoolExecutor(max_workers = nCores) as executor :
            Data = list(executor.map(bondsRed,OctTop,OctBot))
    else:
        with concurrent.futures.ProcessPoolExecutor() as executor :
            Data = list(executor.map(bondsRed,OctTop,OctBot))
            
    Bonds = [D[0] for D in Data]
    BondsIndexes = [D[1] for D in Data]
    if mode == "Indexes":
        return CentIndexes, AdjIndexes, MyCrystal, [Bonds, BondsIndexes], TimeStep*Nsteps

    elif mode == "Dictionary":
        BondDicos = []
        
        for i in range(len(Bonds)) :
            Dico = {}
            SnapBonds = Bonds[i]
#            print(Bonds[i])
            SnapBondsIndexes = BondsIndexes[i]
#            print(SnapBondsIndexes)
#            sys.exit()
            for j in range(len(CentIndexes)//2):
                CentMin = CentIndexes[2*j]
                CentMax = CentIndexes[2*j+1]
                for atom in range(CentMin,CentMax+1):
                    bondsList = SnapBonds[SnapBondsIndexes[atom]:SnapBondsIndexes[atom+1]]
                    if len(bondsList)>1:
                        Dico[atom]=bondsList
            BondDicos.append(Dico)
        
        return CentIndexes, AdjIndexes, MyCrystal, BondDicos, TimeStep*Nsteps

def read_stresses_4visc(UMDfile):
    ff=open(UMDfile,"r")
    AllSnapshots=[]
    
    line = ff.readline().strip().split()
    natom=int(line[1])
    TimeStep=0
    
    while True : 
        SnapshotCrystal = cr.Lattice()
        flagheader=1
        while flagheader :
            l = ff.readline()
            if not l :
                l=ff.readline()
                if not l :#If we encounter two empty lines in a row, we're at the end of the file
                    print('len of allsnapshots is',len(AllSnapshots))
                    return AllSnapshots,TimeStep
            line = l.strip().split()
            if len(line)>0 :
                if line[0]=='timestep':
                    TimeStep = float(line[1])
                if line[0]=='atoms:':
                    flagheader=0
                elif line[0]=='rprimd_a':
                    SnapshotCrystal.rprimd[0]=[float(line[ii]) for ii in range(1,4)]
                elif line[0]=='rprimd_b':
                    SnapshotCrystal.rprimd[1]=[float(line[ii]) for ii in range(1,4)]
                elif line[0]=='rprimd_c':
                    SnapshotCrystal.rprimd[2]=[float(line[ii]) for ii in range(1,4)]
                elif line[0]=='StressTensor':
                    SnapshotCrystal.stress=[float(line[i]) for i in range(1,7)]
                elif line[0]=='Temperature':
                    SnapshotCrystal.temperature=float(line[1])
 
        SnapshotCrystal.cellvolume = SnapshotCrystal.makevolume()
        AllSnapshots.append(SnapshotCrystal)
        
        for _ in range(natom):
            next(ff)

def headerumd():
    print ('\n * * * * * ')
    print ('    UMD package for analyzing Molecular Dynamics simulations.')
    print ('distributed under GNU-GNU General Public License 3')
    print ('authors: Razvan Caracas, Kevin Jiguet-Covex, Tim Bögels, Anne Davis, Xi Zhu, Emma Stoutenbourg')
    print ('please cite as: ')
    print ('    cite as: Caracas, R. et al. Analyzing Melts and Fluids from Ab Initio Molecular Dynamics Simulations with the UMD Package. J. Vis. Exp. (175), e61534,doi:10.3791/61534 (2021)' )
    print (' ')
    
def print_header(FileName,MyCrystal):
    newfile = FileName + '.umd.dat'
    nf = open(newfile,'a')
    string = 'natom ' + str(MyCrystal.natom) + '\n'
    nf.write(string)
    string = 'ntypat ' + str(MyCrystal.ntypat) + '\n'
    nf.write(string)
    string = 'no.of.electrons ' + str(MyCrystal.noelectrons) + '\n'
    nf.write(string)
    string = 'types '
    for ii in range(MyCrystal.ntypat):
        string = string + str(MyCrystal.types[ii]) + ' '
    string = string + '\n'
    nf.write(string)
    string = 'elements '
    for ii in range(MyCrystal.ntypat):
        string = string + str(MyCrystal.elements[ii]) + ' '
    string = string + '\n'
    nf.write(string)
    string = 'masses '
    for ii in range(MyCrystal.ntypat):
        string = string + str(MyCrystal.masses[ii]) + ' '
    string = string + '\n'
    nf.write(string)
    string = 'Zelectrons '
    for ii in range(MyCrystal.ntypat):
        string = string + str(MyCrystal.zelec[ii]) + ' '
    string = string + '\n'
    nf.write(string)
    string = 'typat '
    for ii in range(MyCrystal.natom):
        string = string + str(MyCrystal.typat[ii]) + ' '
    string = string + '\n\n'
    nf.write(string)
    string = 'lambda_ThermoInt' + ' ' +  str(MyCrystal.lambda_ThermoInt)
    string = string + '\n\n'
    nf.write(string)
    nf.close()

def print_snapshots(FileName,MyCrystal,TimeStep,CurrentTime,diffcoords):
    newfile = FileName + '.umd.dat'
    nf = open(newfile,'a')
    string = 'timestep ' + str(TimeStep) + ' fs\n' + 'time ' + str(CurrentTime) + ' fs\n'
    nf.write(string)
    string = 'InternalEnergy ' + str(round(MyCrystal.internalenergy,6)) + ' eV\n'
    nf.write(string)
    string = 'ElectronicEntropy ' + str(round(MyCrystal.electronicentropy,6)) + ' eV\n'
    nf.write(string)
    string = 'KineticEnergyIons ' + str(round(MyCrystal.kineticenergy,6)) + ' eV\n'
    nf.write(string)
    string = 'EnergyWithDrift ' + str(round(MyCrystal.energywithdrift,6)) + ' eV\n'
    nf.write(string)
    string = 'Enthalpy ' + str(round(MyCrystal.enthalpy,6)) + ' eV\n'
    nf.write(string)
    string = 'Magnetization ' + str(round(MyCrystal.magnetization,6)) + ' Magneton-Bohr\n'
    nf.write(string)
    string = 'Temperature ' + str(round(MyCrystal.temperature,1)) + ' K\n'
    nf.write(string)
    string = 'Pressure ' + str(round(MyCrystal.pressure,4)) + ' GPa\n'
    nf.write(string)
    string = 'Density ' + str(round(MyCrystal.density,3)) + ' g.cm-3\n'
    nf.write(string)
    string = 'StressTensor '
    for ii in range(6):
        string = string + str(round(MyCrystal.stress[ii],4)) + ' '
    string = string + ' GPa\n'
    nf.write(string)
#    string = 'acell ' + str(round(MyCrystal.acell[0],3)) + ' ' + str(round(MyCrystal.acell[1],3) + ' ' + str(round(MyCrystal.acell[2],3)) + ' A\n'
    string = 'acell ' + str(MyCrystal.acell[0]) + ' ' + str(MyCrystal.acell[1]) + ' ' + str(MyCrystal.acell[2]) + ' A\n'
    nf.write(string)
    string = 'rprim_a ' + str(MyCrystal.rprimd[0][0]/MyCrystal.acell[0]) + '  ' +str(MyCrystal.rprimd[0][1]/MyCrystal.acell[0]) + '  ' +str(MyCrystal.rprimd[0][2]/MyCrystal.acell[0]) + '\n'
    nf.write(string)
    string = 'rprim_b ' + str(MyCrystal.rprimd[1][0]/MyCrystal.acell[1]) + '  ' +str(MyCrystal.rprimd[1][1]/MyCrystal.acell[1]) + '  ' +str(MyCrystal.rprimd[1][2]/MyCrystal.acell[1]) + '\n'
    nf.write(string)
    string = 'rprim_c ' + str(MyCrystal.rprimd[2][0]/MyCrystal.acell[2]) + '  ' +str(MyCrystal.rprimd[2][1]/MyCrystal.acell[2]) + '  ' +str(MyCrystal.rprimd[2][2]/MyCrystal.acell[2]) + '\n'
    nf.write(string)
    string = 'rprimd_a ' + str(MyCrystal.rprimd[0][0]) + '  ' +str(MyCrystal.rprimd[0][1]) + '  ' +str(MyCrystal.rprimd[0][2]) + ' A\n'
    nf.write(string)
    string = 'rprimd_b ' + str(MyCrystal.rprimd[1][0]) + '  ' +str(MyCrystal.rprimd[1][1]) + '  ' +str(MyCrystal.rprimd[1][2]) + ' A\n'
    nf.write(string)
    string = 'rprimd_c ' + str(MyCrystal.rprimd[2][0]) + '  ' +str(MyCrystal.rprimd[2][1]) + '  ' +str(MyCrystal.rprimd[2][2]) + ' A\n'
    nf.write(string)
    string = 'atoms: reduced*3 cartesian*3(A) abs.diff.*3(A) velocity*3(A/fs) force*3(eV/A) charge(no.elec) magnetization(magneton-Bohr) \n'
    nf.write(string)
    for iatom in range(MyCrystal.natom):
        string = str(round(MyCrystal.atoms[iatom].xred[0],5)) + ' ' + str(round(MyCrystal.atoms[iatom].xred[1],5)) + ' ' + str(round(MyCrystal.atoms[iatom].xred[2],5)) + ' '
        string = string + str(round(MyCrystal.atoms[iatom].xcart[0],5)) + ' ' + str(round(MyCrystal.atoms[iatom].xcart[1],5)) + ' ' + str(round(MyCrystal.atoms[iatom].xcart[2],5)) + ' '
        string = string + str(round(diffcoords[iatom][0],5)) + ' ' + str(round(diffcoords[iatom][1],5)) + ' ' + str(round(diffcoords[iatom][2],5)) + ' '
        string = string + str(round(MyCrystal.atoms[iatom].vels[0],5)) + ' ' + str(round(MyCrystal.atoms[iatom].vels[1],5)) + ' ' + str(round(MyCrystal.atoms[iatom].vels[2],5)) + ' '
        string = string + str(round(MyCrystal.atoms[iatom].forces[0],5)) + ' ' + str(round(MyCrystal.atoms[iatom].forces[1],5)) + ' ' + str(round(MyCrystal.atoms[iatom].forces[2],5)) + ' '
        string = string + str(MyCrystal.atoms[iatom].charge) + ' ' + str(MyCrystal.atoms[iatom].magnet) + '\n'
        nf.write(string)
    string='\n'
    nf.write(string)
    nf.close()
    return(CurrentTime,TimeStep)

def read_bigheader_umd(umdfile, short=0):
    MyCrystal = cr.Lattice()
    with open(umdfile,'r') as ff:
        while True:
            line = ff.readline()
            if not line: break
            #print(line,len(line))
            if len(line) > 1:
                line=line.strip()
                entry=line.split()
                if entry[0] == 'natom':
                    MyCrystal.natom = int(entry[1])
                    MyCrystal.typat = [0 for _ in range(MyCrystal.natom)]
                if entry[0] == 'ntypat':
                    MyCrystal.ntypat = int(entry[1])
                    MyCrystal.types = [0 for _ in range(MyCrystal.ntypat)]
                    MyCrystal.elements = ['X' for _ in range(MyCrystal.ntypat)]
                    MyCrystal.masses = [1.0 for _ in range(MyCrystal.ntypat)]
                    MyCrystal.zelec = [1.0 for _ in range(MyCrystal.ntypat)]
                if entry[0] == 'types':
                    for ii in range(MyCrystal.ntypat):
                        MyCrystal.types[ii]=int(entry[ii+1])
                if entry[0] == 'elements':
                    for ii in range(MyCrystal.ntypat):
                        MyCrystal.elements[ii]=entry[ii+1]
                if entry[0] == 'masses':
                    for ii in range(MyCrystal.ntypat):
                        MyCrystal.masses[ii]=float(entry[ii+1])
                if entry[0] == 'Zelectrons':
                    for ii in range(MyCrystal.ntypat):
                        MyCrystal.zelec[ii]=float(entry[ii+1])
                if entry[0] == 'typat':
                    for ii in range(MyCrystal.natom):
                        MyCrystal.typat[ii] = int(entry[ii+1])
                if entry[0] == 'timestep':
                    TimeStep = float(entry[1])
                    break
    return(MyCrystal,TimeStep)

def copyhead(Target,File,Nsteps):
    
    ff = open(File,"r")
    fa = open(Target,"w")
    
    while True : 
        line = ff.readline()
        if not line :
            ff.close()
            fa.close()
            break
        l = line.strip().split()
        
        if len(l)>0 :
            if l[0] == "timestep":
                line = l[0] +" "+ str(float(l[1])*Nsteps) + " "+ l[2] +"\n"
            if l[0] == "atoms:":
                ff.close()
                fa.close()
                break                
        fa.write(line)
