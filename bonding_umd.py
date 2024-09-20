#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 10:55:47 2023

@author: Kevin Jiguet-Covex
"""

import sys,getopt,os.path
import math
import time
from functools import partial
import concurrent.futures
import numpy as np
import ctypes
from os.path import join
import umd_processes_fast as umdpf
import platform

OS = platform.system()
LibraryName =""

if OS == "Linux":
    LibraryName = 'c_bonds_full.so'
elif OS == "Windows":
    LibraryName = 'c_bonds_full.dll'
elif OS == "Darwin":
    LibraryName = 'c_bonds_full.dylib'
    
current_path=os.path.abspath(__file__)
path_split=current_path.split('/')
path_red=path_split[1:-1]
path_new=''
for u in path_red:
    path_new+='/'+u


#Definition of the library and the functions
fullbond_lib = ctypes.cdll.LoadLibrary(join(path_new, LibraryName))


fullbond_lib.compute_Bonds_full.argtypes = [ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_int),ctypes.c_int,ctypes.c_int,ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.c_int,ctypes.POINTER(ctypes.c_double),ctypes.c_int,ctypes.c_int]
fullbond_lib.compute_Bonds_full.restype = ctypes.POINTER(ctypes.POINTER(ctypes.c_int))

fullbond_lib.free_memory.argtypes = [ctypes.POINTER(ctypes.POINTER(ctypes.c_int))]
fullbond_lib.free_memory.restypes = None

fullbond_lib.free_memory_int.argtypes = [ctypes.POINTER(ctypes.c_int)]
fullbond_lib.free_memory.restypes = None


def read_inputfile(InputFile,MyCrystal):#Creates a matrix from the .dat file containing the bond length for each pair of atom types
    BondTable = [[0.0 for _ in range(MyCrystal.ntypat)] for _ in range(MyCrystal.ntypat)]
    with open(InputFile) as ff:
        for line in ff:
            line=line.strip()
            entry=line.split()
            if (len(entry)==3):
                for ii in range(MyCrystal.ntypat):
                    if MyCrystal.elements[ii]==entry[0]:
                        for jj in range(MyCrystal.ntypat):
                            if MyCrystal.elements[jj]==entry[1]:
                                BondTable[ii][jj]=float(entry[2])*float(entry[2])
                                BondTable[jj][ii]=float(entry[2])*float(entry[2])
    return BondTable

def WriteBonding_C_full(MySnapshotL,step,maxSteps,MyCrystal,BondTable,timestep,natom,numCells,acell,specifics,ortho=True):#Calculates bonds in a given snapshot
    
    if(100*step//maxSteps != 100*(step-1)//maxSteps or True):
        sys.stdout.write("\rDetermining bonds between atoms. Progress : "+str(100*step//maxSteps)+"%") 
        sys.stdout.flush()

    #Preparing data to be used by the C function, pointers creation
    MSp = np.array(list(MySnapshotL)).ctypes.data_as(ctypes.POINTER(ctypes.c_double))    #C pointers
    BTp = np.array(BondTable).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    CrystalTypes = []
    for ty in range(MyCrystal.ntypat):
        CrystalTypes+=[ty for _ in range(MyCrystal.types[ty])]
        
    CTp = (ctypes.c_int * natom)(*CrystalTypes)
        
    if specifics == None :#No restriction of the spatial interval
        specifics = [0]
        specification = 0
        SPp = (ctypes.c_double * 1)(*specifics)
    else:
        specification = 1
        SPp = (ctypes.c_double * 6)(*specifics)
    
    if ortho :#Orthogonal cell
        RVp= np.array([0]).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        Cellp = (ctypes.c_double *3)(*acell)
        diag = 1
    else:
        Cellp = (ctypes.c_double * 9)(*acell[0])
        RVp = (ctypes.c_double * 9)(*acell[1])
        diag = 0

    Bonds = fullbond_lib.compute_Bonds_full(MSp,BTp,CTp,natom,MyCrystal.ntypat,Cellp,RVp,numCells,SPp,specification,diag)
               
    BondsList=[[at] for at in range(natom)]


    #Converting bond profile from C to python list
    #The first number of each sub-list is the length of said sub-list  
    for at in range(natom):
        L = Bonds[at][0]
        for ineighbor in range(1,L+1):
            BondsList[Bonds[at][ineighbor]].append(at)
        fullbond_lib.free_memory_int(Bonds[at])#Freeing the memory of each sub-list
    
    fullbond_lib.free_memory(Bonds)#Freeing the memory of the list
    return BondsList        

def main(argv):
    umdpf.headerumd()
    UMDname=''
    Nsteps = 1
    InputFile = ''
    header = ''
    numCells=None
    maxlength=None
    specifics = None
    nCores = None
    try:
        opts, arg = getopt.getopt(argv,"hf:s:l:i:n:p:k:",["fUMDfile","sSampling_Frequency","lMaxLength","iInputFile","nNumCells","pSpecifics","knCores"])
    except getopt.GetoptError:
        print ('bonding_umd.py -f <UMD_filename> -s <Sampling_Frequency> -l <MaxLength> -i <InputFile> -n <NumCells> -p <Specifics> -k <nCores>')
        sys.exit(2)
    if opts == [] :
        opts = [('-h','')]
    for opt, arg in opts:
        if opt == '-h':
            print ('Computation of the bonding map for each snapshot of a umd file')
            print ('bonding_umd.py -f <UMD_filename> -s <Sampling_Frequency> -l <MaxLength> -i <InputFile> -n <NumCells> -p <Specifics> -k <nCores>')
            print ('default values: -f output.umd.dat -s 1 -l None -p None -k None')
            print ('The input file contains the bond lengths for the different atom pairs. \n The option -l overwrites all the values of this file.')
            print ('-l : unique bonding length. Use in the case of not having an input file -i. ')
            print ('-n : states the number of sub-cells the script will work with. The default vomatically computed for optimal performances.')
            print ('-p : restricts the computation of the bonds in a specific place of the simulation. Takes a list of floats under the form -p [xmin,ymin,zmin,xmax,ymax,zmax], with the list contains the maximal and minimal coordinates values for the atoms to be considered. Default None (no restriction).')
            print ('-k : defines the number of cores that will be used in the parallel calculation. By default, this number is automatically calculated by Python.')
            sys.exit()
        elif opt in ("-l","--lmaxlength"):
            if(arg!=""):
                maxlength = float(arg)    
                header = header + "-l=" + str(maxlength)
        elif opt in ("-f", "--fUMDfile"):
            UMDname = str(arg)
            header = header + 'FILE: -f=' + UMDname
        elif opt in ("-s","--sNsteps"):
            Nsteps = int(arg)
            header = header + ' -s=' + arg
            print('Will sample the trajectories every ',Nsteps,' steps')
        elif opt in ("-i", "--iInputFile"):
            InputFile = str(arg)
            header = header + ' -i=' + InputFile
            print ('Bonding cutoffs to be read from file ',InputFile)
        elif opt in ("-n", "--nNumCells"):
            numCells = int(arg)
            header = header + ' -n=' + arg
        elif opt in ("-p","--pSpecifics"):
            specifics = eval(arg)
        elif opt in ("-k","--kCores"):
            nCores = int(arg)

    if not (os.path.isfile(UMDname)):
        print ('the UMD files ',UMDname,' does not exist')            
        sys.exit()

    sys.stdout.write(str(os.cpu_count())+" cores are available. ")
    if nCores != None :
        sys.stdout.write("We will use "+str(nCores)+" of them.\n")        
    else :
        sys.stdout.write("The number that will be used is automatically determined.\n")


    start=time.time()

    MyCrystal,AllSnapshotsL,TimeStep,lensnap = umdpf.read_values(UMDname,"xcart","line",Nsteps,nCores = nCores)#Creating the list of the xcart coordinates from the UMD file
    rprimd=MyCrystal.rprimd
    acell=MyCrystal.acell

    ortho = False    
    if [rprimd[0][1],rprimd[0][2],rprimd[1][0],rprimd[1][2],rprimd[2][0],rprimd[2][1]] == [0,0,0,0,0,0]:#If the rprimd matrix is diagonal, the cell is otrhogonal
        ortho = True

    if maxlength==None and len(InputFile)>0 :
        BondTable = read_inputfile(InputFile,MyCrystal)#Converts the input file into a matrix
    elif maxlength!=None :#If a unique bonding length is specified
        BondTable = [[maxlength**2 for _ in range(MyCrystal.ntypat)] for _ in range(MyCrystal.ntypat)]#Creates a matrix filled with -l value


    M=math.sqrt(max([max(Bondlengths) for Bondlengths in BondTable]))#maximal length of any bond

    if numCells == None :#If the number of subcells isn't specified, it is automatically calculated. The size of a subcell is smaller than the biggest bond length.
        numCells = int(min(acell)/M)
        print("Number of cells automatically fixed to "+str(numCells))
        
    if acell[0]/numCells<M or acell[1]/numCells<M or acell[2]/numCells<M:
        print('WARNING : one or more dimension(s) of the sub-cells smaller than the greatest bond length.')            

    natom=MyCrystal.natom#Creating the output file
    if UMDname[-8:] == ".umd.dat":
        FileAll=UMDname[:-8]+'.bonding.dat'
    else :
        FileAll=UMDname+'.bonding.dat'
            
    umdpf.copyhead(FileAll,UMDname,Nsteps)     #Copies the header of the UMD file into the output file
#    ortho = True
    if ortho == False :
        V0 = np.array(MyCrystal.rprimd[0])
        V1 = np.array(MyCrystal.rprimd[1])
        V2 = np.array(MyCrystal.rprimd[2])
    
        #Vectors of the reciprocal space
        RecVecs0 = np.cross(V1,V2)/np.dot(np.cross(V1,V2),V0)
        RecVecs1 = np.cross(V2,V0)/np.dot(np.cross(V2,V0),V1)
        RecVecs2 = np.cross(V0,V1)/np.dot(np.cross(V0,V1),V2)
        
        RecVecs = list(RecVecs0)+list(RecVecs1)+list(RecVecs2)
        CellVecs = list(V0)+list(V1)+list(V2)
        
        acell = [CellVecs,RecVecs]

    WriteBondingRed=partial(WriteBonding_C_full,maxSteps=len(AllSnapshotsL)-1,MyCrystal=MyCrystal,BondTable=BondTable,timestep=TimeStep,natom=natom,numCells=numCells,acell=acell,specifics = specifics,ortho=ortho)
        
    
#    Lines = []
#    for i in range(len(AllSnapshotsL)):
#        Lines.append(WriteBondingRed(AllSnapshotsL[i],i))
#    print(Lines[0])
#    sys.exit()

    if nCores != None :    
        with concurrent.futures.ProcessPoolExecutor(max_workers = nCores) as executor :#Calculation of the bonds
            Lines=list(executor.map(WriteBondingRed,AllSnapshotsL,[step for step in range(len(AllSnapshotsL))]))#Calculates the bond profile for each snapshot
    else :
        with concurrent.futures.ProcessPoolExecutor() as executor :#Calculation of the bonds
            Lines=list(executor.map(WriteBondingRed,AllSnapshotsL,[step for step in range(len(AllSnapshotsL))]))#Calculates the bond profile for each snapshot

    print("\n")
    step=0
    print("Writing...")#Writing the bonds
    fa=open(FileAll,'a')
    
    if specifics != None :
        line = "Space interval ["+str(specifics[0])+","+str(specifics[3])+"] "+"["+str(specifics[1])+","+str(specifics[4])+"] "+"["+str(specifics[2])+","+str(specifics[5])+"]\n"
        fa.write(line)
    for lines in Lines :                
        header = 'time '+str(step*TimeStep*Nsteps)+' fs\nstep '+str(step)+'\n'
        fa.write(header)
        for line in lines :
            l=''
            for atom in line :
                l+=str(atom)+'\t'
            l+='\n'
            fa.write(l)
                
        fa.write('end\n')
        step+=1
                            
    fa.close()
    
    
    print ('Bonds written in <',FileAll,'> file')

    end = time.time()
    print("runtime:",end-start," s")    

if __name__ == "__main__":
    main(sys.argv[1:])
