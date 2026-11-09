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


fullbond_lib.compute_Bonds_full.argtypes = [ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_int),ctypes.c_int,ctypes.c_int,ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_double),ctypes.c_int]
fullbond_lib.compute_Bonds_full.restype = ctypes.POINTER(ctypes.POINTER(ctypes.c_int))

fullbond_lib.free_memory.argtypes = [ctypes.POINTER(ctypes.POINTER(ctypes.c_int))]
fullbond_lib.free_memory.restypes = None

fullbond_lib.free_memory_int.argtypes = [ctypes.POINTER(ctypes.c_int)]
fullbond_lib.free_memory.restypes = None

fullbond_lib.bonds_text.argtypes = fullbond_lib.compute_Bonds_full.argtypes+[ctypes.POINTER(ctypes.c_long)]
fullbond_lib.bonds_text.restype = ctypes.c_void_p
fullbond_lib.free_text.argtypes = [ctypes.c_void_p]
fullbond_lib.free_text.restype = None


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

def bonding_batch(task,UMDname,natom,BondTable,CrystalTypes,ntypat,numCells,specifics,first,Nsteps,TimeStep):#The text of the bonding file for
    #the snapshots start to stop-1 of the selection (task = (start, stop, byte offset of the snapshot start in the UMD file)), read here from the UMD file
    start, stop, offset = task
    X, cells = umdpf.read_atoms_at(UMDname,offset,natom,stop-start,Nsteps,3,3,cells=True)
    BTp = np.ascontiguousarray(BondTable,dtype=float).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    CTp = (ctypes.c_int * natom)(*CrystalTypes)
    SPp = (ctypes.c_double * 6)(*specifics) if specifics is not None else (ctypes.c_double * 1)(0)
    NCp = (ctypes.c_int * 3)(*numCells)
    length = ctypes.c_long(0)
    parts = []
    for k in range(stop-start):
        step = start+k
        Cellp = (ctypes.c_double * 9)(*cells[k].flatten())#lattice vectors and reciprocal vectors (s_k = x.rvec_k), any cell
        RVp = (ctypes.c_double * 9)(*np.linalg.inv(cells[k]).T.flatten())
        x = np.ascontiguousarray(X[k])
        text = fullbond_lib.bonds_text(x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),BTp,CTp,natom,ntypat,Cellp,RVp,NCp,SPp,int(specifics is not None),ctypes.byref(length))
        if not text:
            raise MemoryError('bonds of the snapshot '+str(step))
        parts += ['time '+str(first*TimeStep+step*TimeStep*Nsteps)+' fs\nstep '+str(step)+'\n',ctypes.string_at(text,length.value).decode(),'end\n']
        fullbond_lib.free_text(text)
    return ''.join(parts)

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
    first = 0
    try:
        opts, arg = getopt.getopt(argv,"hf:s:l:i:c:n:p:k:",["fUMDfile","sSampling_Frequency","lMaxLength","iInitialStep","cInputFile","nNumCells","pSpecifics","knCores"])
    except getopt.GetoptError:
        print ('bonding_umd.py -f <UMD_filename> -c <InputFile> -l <MaxLength> -s <Sampling_Frequency> -i <InitialStep> -n <NumCells> -p <Specifics> -k <nCores>')
        sys.exit(2)
    if opts == [] :
        opts = [('-h','')]
    for opt, arg in opts:
        if opt == '-h':
            print ('Computation of the bonding map for each snapshot of a umd file')
            print ('bonding_umd.py -f <UMD_filename> -c <InputFile> -l <MaxLength> -s <Sampling_Frequency> -i <InitialStep> -n <NumCells> -p <Specifics> -k <nCores>')
            print ('default values: -f output.umd.dat -s 1 -i 0 -l None -p None -k None')
            print ('-f all : every *.umd.dat file of the current folder, one after the other.')
            print ('-c : input file with the bond lengths of the pairs of elements, one line A B cutoff, e.g. <base>.maximum_bonds.inp of analyze_gofr.py')
            print ('     (formerly -i, still accepted : -i followed by a file name). The option -l overwrites all the values of this file.')
            print ('-l : unique bonding length. Use in the case of not having an input file -c. ')
            print ('-s : one snapshot out of s is used. Default 1.')
            print ('-i : InitialStep. number of initial snapshots discarded. Default 0. The bonding file then starts with the snapshot i of the UMD')
            print ('     file (line first_snapshot i of its header, used by speciation_and_angles.py for the angles and by msd_clusters_umd.py).')
            print ('-n : number of sub-cells along each direction. Default : along each of a, b, c, the largest number whose sub-cells are wider than the longest bond.')
            print ('Any cell, orthogonal or not, changing or not (NPT) : each pair is bonded once, in its nearest image (reduced coordinates), when it is shorter than its cutoff.')
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
        elif opt in ("-i", "--iInitialStep") and arg.strip().isdigit():
            first = int(arg)
            header = header + ' -i=' + arg
        elif opt in ("-c", "--cInputFile", "-i", "--iInitialStep"):#-i <file> : the former name of -c
            InputFile = str(arg)
            header = header + ' -c=' + InputFile
            print ('Bonding cutoffs to be read from file ',InputFile)
        elif opt in ("-n", "--nNumCells"):
            numCells = int(arg)
            header = header + ' -n=' + arg
        elif opt in ("-p","--pSpecifics"):
            specifics = eval(arg)
        elif opt in ("-k","--kCores"):
            nCores = int(arg)

    if not (os.path.isfile(UMDname)):
        umdpf.missing_umd(UMDname)
        sys.exit()

    sys.stdout.write(str(os.cpu_count())+" cores are available. ")
    if nCores != None :
        sys.stdout.write("We will use "+str(nCores)+" of them.\n")        
    else :
        sys.stdout.write("The number that will be used is automatically determined.\n")


    start=time.time()

    MyCrystal,TimeStep = umdpf.Crystallization(UMDname)
    natom=MyCrystal.natom
    offsets, cells = umdpf.index_snapshots(UMDname,natom)#byte offset and lattice vectors (rows) of each snapshot (NPT)
    offsets, cells = offsets[first::Nsteps], cells[first::Nsteps]#the snapshots used
    if len(offsets) == 0:
        print('No snapshot after the snapshot',first)
        sys.exit()
    if np.isnan(cells).any():
        print('ERROR : snapshots of',UMDname,'without rprimd or acell lines')
        sys.exit()
    widths = (np.abs(np.linalg.det(cells))[:,None]/np.linalg.norm(np.cross(cells[:,[1,2,0]],cells[:,[2,0,1]]),axis=2)).min(axis=0)#smallest distances between opposite faces
    if maxlength==None and len(InputFile)>0 :
        BondTable = read_inputfile(InputFile,MyCrystal)#Converts the input file into a matrix
    elif maxlength!=None :#If a unique bonding length is specified
        BondTable = [[maxlength**2 for _ in range(MyCrystal.ntypat)] for _ in range(MyCrystal.ntypat)]#Creates a matrix filled with -l value
    else :
        print('Give the bond lengths with -c or -l')
        sys.exit(2)

    M=math.sqrt(max([max(Bondlengths) for Bondlengths in BondTable]))#maximal length of any bond

    if numCells == None :#If the number of subcells isn't specified, it is the largest one whose sub-cells are wider than the longest bond
        numCells = [max(1,min(64,int(w/M))) if M > 0 else 1 for w in widths]
        print("Number of sub-cells along a, b, c automatically fixed to "+' '.join(str(n) for n in numCells))
    else :
        numCells = [numCells]*3

    if any(w/n < M for w,n in zip(widths,numCells)):
        print('WARNING : one or more width(s) of the sub-cells smaller than the greatest bond length : bonds will be missed.')
    if 2*M > widths.min():
        print('WARNING : the greatest bond length is larger than half the width of the cell : only the nearest image of each pair is bonded.')

    if UMDname[-8:] == ".umd.dat":#the output file
        FileAll=UMDname[:-8]+'.bonding.dat'
    else :
        FileAll=UMDname+'.bonding.dat'

    umdpf.copyhead(FileAll,UMDname,Nsteps,first)     #Copies the header of the UMD file into the output file
    CrystalTypes = [ty for ty in range(MyCrystal.ntypat) for _ in range(MyCrystal.types[ty])]
    #Batches of snapshots, sized from the bytes of a snapshot in the file : each worker reads its snapshots from the UMD file and returns the
    #text of their bonds, written in order as the batches come (the memory holds a few batches, whatever the size of the run)
    nsnap = len(offsets)
    snapbytes = (os.path.getsize(UMDname)-offsets[0])/(nsnap*Nsteps)
    ranges = umdpf.batch_ranges(nsnap,snapbytes*Nsteps,nCores)
    worker = partial(bonding_batch,UMDname=UMDname,natom=natom,BondTable=BondTable,CrystalTypes=CrystalTypes,ntypat=MyCrystal.ntypat,
                     numCells=numCells,specifics=specifics,first=first,Nsteps=Nsteps,TimeStep=TimeStep)
    print('Determining the bonds of',nsnap,'snapshots in',len(ranges),'batches of up to',ranges[0][1]-ranges[0][0],'snapshots')
    with open(FileAll,'a') as fa:
        if specifics != None :
            fa.write("Space interval ["+str(specifics[0])+","+str(specifics[3])+"] "+"["+str(specifics[1])+","+str(specifics[4])+"] "+"["+str(specifics[2])+","+str(specifics[5])+"]\n")
        for (a,b),text in zip(ranges,umdpf.ordered_map(worker,[(a,b,offsets[a]) for a,b in ranges],nCores)):
            fa.write(text)
            sys.stdout.write("\rDetermining bonds between atoms. Progress : "+str(100*b//nsnap)+"%")
            sys.stdout.flush()
    print("\n")
    print ('Bonds written in <',FileAll,'> file')

    end = time.time()
    print("runtime:",end-start," s")    

if __name__ == "__main__":
    if not umdpf.run_on_all(sys.argv[1:],main):#-f all : every UMD file of the folder
        main(sys.argv[1:])
