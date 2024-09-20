#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 09:55:21 2023

@author: Kevin Jiguet-Covex
"""
import math
import platform
import sys,getopt,os.path
import umd_processes_fast as umdpf
import time
from functools import partial
import concurrent.futures
from os.path import join
import ctypes
import numpy as np
import matplotlib.pyplot as plt


OS = platform.system()
LibraryName =""

if OS == "Linux":
    LibraryName = 'c_gofr_new_alt.so'
elif OS == "Windows":
    LibraryName = 'c_gofr_new_alt.dll'
elif OS == "Darwin":
    LibraryName = 'c_gofr_new_alt.dylib'

current_path=os.path.abspath(__file__)
path_split=current_path.split('/')
path_red=path_split[1:-1]
path_new=''
for u in path_red:
    path_new+='/'+u


gofr_lib = ctypes.cdll.LoadLibrary(join(path_new, LibraryName))

gofr_lib.compute_gofr.argtypes = (ctypes.POINTER(ctypes.c_double),ctypes.c_int,ctypes.POINTER(ctypes.c_int),ctypes.c_double,ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.c_int,ctypes.c_double,ctypes.c_int)
gofr_lib.compute_gofr.restype = ctypes.POINTER(ctypes.c_int)


gofr_lib.free_memory.argtypes = [ctypes.POINTER(ctypes.c_int)]
gofr_lib.free_memory.restype = None

def compute_gofr_snapshot(MySnapshot,step,maxStep,Types,ntypes,discrete,acell,ndivx,maxlength,ortho):#Computes the distance matrix for each snapshot
    
    if(100*step//maxStep != 100*(step-1)//maxStep):
        sys.stdout.write('\rCalculating gofr. Progress : '+str(100*step//maxStep)+"%")
        sys.stdout.flush()
        
    if ortho :#If the unit cell is orthogonal
        RVp = np.array([0]).ctypes.data_as(ctypes.POINTER(ctypes.c_double))#Will not be used but is required as an argument
        acp = (ctypes.c_double *3)(*acell)#Lengths of each vector
        diag = 1
    else:
        acp = (ctypes.c_double * 9)(*acell[0])#Direct vectors
        RVp = (ctypes.c_double * 9)(*acell[1])#Reciprocal vectors
        diag = 0

    #Creation of the pointers
    Tp = (ctypes.c_int * (ntypes+2))(*Types)
    MSp = (ctypes.c_double * (3*Types[-1]))(*MySnapshot)

    gofr = np.zeros((ntypes*ntypes*(ndivx+1)),dtype=int)#The numpy matrix is directly given to the c function and retrieved afterwards

    gofr_C = gofr_lib.compute_gofr(MSp,ntypes,Tp,discrete,acp,RVp,ndivx+1,maxlength,diag)
    
    for i in range(ntypes*ntypes*(ndivx+1)):
                   gofr[i]=gofr_C[i]
    
    gofr_lib.free_memory(gofr_C)
    
    return gofr

def main(argv):
    start=time.time()
    umdpf.headerumd()
    nCores = None
    umdfile = ''
    Nsteps = 1
    discrete = 0.01            #delta_r  = width of bins in histogram
    InitialStep = 0           #in case we want to skip additional timesteps
    try:
        opts, arg = getopt.getopt(argv,"hf:s:d:i:k:", ["fumdfile", "Sampling_Frequency", "ddiscrete", "iInitialStep","knCores"])
    except getopt.GetoptError:
        print ('gofrs_umd.py -f <umdfile>  -s <No.steps>  -d <discretization.interval> -i <InitialStep> -k <nCores>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('gofrs_umd.py program to compute pair-distribution fuctions g(r) and print all the results in one file. Usage:')
            print ('gofrs_umd.py -f <umdfile> -s <Sampling_Frequency> -d <discretization.interval> -i <InitialStep> -b <UseGPU>')
            print ('umdfile = name of the umd file. Default = output.umd.dat')
            print ('Sampling_Frequency = frequency of sampling the trajectory. Default = 1')
            print ('discretization.interval = for plotting the g(r). Default = 0.01 Angstroms')
            print ('InitialStep = number of initial steps of the trajectory to be removed. Default = 0')
            print ("-k : number of cores to be used during the parallelization processes. By default, it's left to the machine.")
            sys.exit()

        elif opt in ("-f", "--fumdfile"):
            umdfile = str(arg)
            #print('umdfile = ',umdfile)
        elif opt in ("-s", "--sNsteps"):
            Nsteps = int(arg)
        elif opt in ("-d", "--ddiscrete"):
            discrete = float(arg)
            #print('the distance delta_r we use to compute the number of atoms around the central one is ',discrete, 'Angstroms')
        elif opt in ("-i", "--iInitialStep"):
            InitialStep = int(arg)
            #print('I will skip  ',initial,' timesteps. ')
        elif opt in("-k-","--knCores"):
            nCores = int(arg)

    MyCrystal,SnapshotsXCart,TimeStep,length = umdpf.read_values(umdfile,"xcart",Nsteps = Nsteps,firststep = InitialStep,nCores =nCores)

    rprimd=MyCrystal.rprimd
    acell=MyCrystal.acell

    ortho = False    
    if [rprimd[0][1],rprimd[0][2],rprimd[1][0],rprimd[1][2],rprimd[2][0],rprimd[2][1]] == [0,0,0,0,0,0]:
        ortho = True

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
        maxlength = min(np.linalg.norm(V0),np.linalg.norm(V2),np.linalg.norm(V2))


    else :
        acell = MyCrystal.acell
        maxlength = min(MyCrystal.acell[:2])



    Types=[sum(MyCrystal.types[:i]) for i in range(MyCrystal.ntypat)]
    Types.append(MyCrystal.natom)
    Types.append(MyCrystal.natom)
    ndivx = int(maxlength / (2 * discrete))
        
    ntypes = len(Types)-2
    compute_gofr_red = partial(compute_gofr_snapshot,maxStep=length-1,Types=Types,ntypes=ntypes,discrete=discrete,acell=acell,ndivx=ndivx,maxlength=maxlength,ortho=ortho)


    if nCores != None :
        with concurrent.futures.ProcessPoolExecutor(max_workers=nCores) as executor :
            Gofrs = list(executor.map(compute_gofr_red,SnapshotsXCart,[step for step in range(len(SnapshotsXCart))]))
    else:
        with concurrent.futures.ProcessPoolExecutor() as executor :
            Gofrs = list(executor.map(compute_gofr_red,SnapshotsXCart,[step for step in range(len(SnapshotsXCart))]))
        
    print("\n")

    Gofr = sum(Gofrs)
    G = np.zeros((ntypes,ntypes,ndivx+1),dtype=int)
    #We fill a 3D matrix with the values of the array     
    for i in range(ntypes):
        for j in range(i,ntypes):
            for k in range(ndivx+1):
                G[i][j][k]=Gofr[k*(ntypes*ntypes)+i*ntypes+j]
                G[j][i][k]=Gofr[k*(ntypes*ntypes)+i*ntypes+j]
        
    gofrname = umdfile[:-8]+".gofr.dat"
    fa = open(gofrname,'w')
        
    header = "dist\t"
    for iel in range(ntypes):
        for jel in range(ntypes):
            header+=MyCrystal.elements[iel]+"-"+MyCrystal.elements[jel]+'\tInt('+MyCrystal.elements[iel]+"-"+MyCrystal.elements[jel]+')\t'
    fa.write(header[:-1]+"\n")
        
    Int = np.zeros((ntypes,ntypes))
    volcell = maxlength**3#Surrounding volume that we took in account (for each atom) for the calculation of the gofr
        
    for kk in range(ndivx+1):
        smallr = kk*discrete
        bigr = smallr + discrete
        volshell = 4/3*math.pi*(bigr**3-smallr**3)#Spatial normalization of the number of atoms at a given radius       
        string = str(round(kk*discrete+discrete/2,5))+"\t"
        for iel in range(ntypes):
            for jel in range(ntypes):
                temp = G[iel][jel][kk]*volcell/(volshell*MyCrystal.types[iel]*MyCrystal.types[jel]*length)
                if temp==0:
                    temp = int(0)
                Int[iel][jel]+=temp/volcell*4*math.pi*MyCrystal.types[jel]*discrete*smallr**2
                if Int[iel][jel]==0:
                    st = str(int(0))
                else :
                    st = str(Int[iel][jel])
                string+=str(temp)+'\t'+st+"\t"
        fa.write(string+"\n")
                    
    fa.close()
        
    print("gofr file successfully created under the name "+gofrname)
    print("time :",time.time()-start)

if __name__ == "__main__":
   main(sys.argv[1:])
