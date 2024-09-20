#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### Author : Kevin Jiguet-Covex



import platform
import sys,getopt,os.path
import umd_processes_fast as umdpf
import time
import numpy as np
import concurrent.futures
from functools import partial
import ctypes
from os.path import join
import math


current_path=os.path.abspath(__file__)#For all this to work, the file c_msd.so must be in the same directory than this script
path_split=current_path.split('/')
path_red=path_split[1:-1]
path_new=''
for u in path_red:
    path_new+='/'+u

OS = platform.system()
LibraryName =""

if OS == "Linux":
    LibraryName = 'c_msd.so'
elif OS == "Windows":
    LibraryName = 'c_msd.dll'
elif OS == "Darwin":
    LibraryName = 'c_msd.dylib'

msd_lib = ctypes.cdll.LoadLibrary(join(path_new, LibraryName))
msd_lib.compute_msd_tilted.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.POINTER(ctypes.c_double)]
msd_lib.compute_msd_tilted.restype = ctypes.POINTER(ctypes.c_double)

msd_lib.compute_msd_tilted_multi.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.POINTER(ctypes.c_double),ctypes.c_int]
msd_lib.compute_msd_tilted_multi.restype = ctypes.POINTER(ctypes.c_double)


msd_lib.free_memory.argtypes = [ctypes.POINTER(ctypes.c_double)]
msd_lib.free_memory.restype = None

msd_lib.compute_msd.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int]
msd_lib.compute_msd.restype = ctypes.POINTER(ctypes.c_double)



def msdAtom_C(n,PosAr,numsteps,hh,vv,ballistic,natom):#To calculate the MSD without projection
    if((100*n)//natom!=(100*(n-1))//natom):
        sys.stdout.write("\rCalculating MSD. Progress : "+str((100*n)//natom)+"%")
        sys.stdout.flush()
    nitmax=(numsteps//2-ballistic)
    Posp = (PosAr.flatten()).ctypes.data_as(ctypes.POINTER(ctypes.c_double))#Creation of the pointer for the C function
    msd=[]
    msdP=msd_lib.compute_msd(Posp,hh,vv,ballistic,nitmax)#Calculation of the msd
    for i in range(int(nitmax/vv)):
        msd.append(msdP[i])
    msd_lib.free_memory(msdP)#liberation of memory
    return msd

def msdAtom_CAxes(n,PosAr,numsteps,hh,vv,ballistic,Axes,natom):#To calculate the MSD projected along 3 axes
    if((100*n)//natom!=(100*(n-1))//natom):
        sys.stdout.write("\rCalculating MSD. Progress : "+str((100*n)//natom)+"%")
        sys.stdout.flush()
    nitmax=(numsteps//2-ballistic)
    #Creation of the pointers
    Axesp = (9 * ctypes.c_double)(*Axes)
    Posp = (PosAr.flatten()).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    msd=msd_lib.compute_msd_tilted(Posp,hh,vv,ballistic,nitmax,Axesp)#Calculation of the msd
    msdA1,msdA2,msdA3 = [], [], []
    for i in range(int(nitmax/vv)):
        msdA1.append(msd[i])
        msdA2.append(msd[i+int(nitmax/vv)])
        msdA3.append(msd[i+2*int(nitmax/vv)])
    msd_lib.free_memory(msd)#liberation of memory
    return msdA1,msdA2,msdA3

def msdAtom_CAxes_multi(n,PosAr,numsteps,hh,vv,ballistic,Axes,nAxes,natom):#To calculate the msd along an arbitrary number of axes

    if((100*n)//natom!=(100*(n-1))//natom):
        sys.stdout.write("\rCalculating MSD. Progress : "+str((100*n)//natom)+"%")
        sys.stdout.flush()
    nitmax=(numsteps//2-ballistic)
    #Creation of the pointers
    Axesp = (3*nAxes * ctypes.c_double)(*Axes)
    Posp = (PosAr.flatten()).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    msd=msd_lib.compute_msd_tilted_multi(Posp,hh,vv,ballistic,nitmax,Axesp,int(len(Axes)/3))#Calculation of the msd
    msds = [[] for _ in range(nAxes)]
    for i in range(int(nitmax/vv)):
        for j in range(int(len(Axes)/3)):
            msds[j].append(msd[i+j*int(nitmax/vv)])
    msd_lib.free_memory(msd)#liberation of memory
    return msds


def main(argv):
    hh = 1
    vv = 1
    ballistic = 0
    TimeStep = 1
    umdpf.headerumd()
    start=time.time()
    umdfile=''
    Mode = "elements"
    Auto = False
#    Axes=[1,0,0,0,1,0,0,0,1]
    Axes = None
    nCores=None
    try:
        opts, arg = getopt.getopt(argv,"hf:z:v:b:m:a:k:",["fumdfile","zHorizontalJump","vVerticalJump","bBallistic","mMode","aAxes","knCores"])
    except getopt.GetoptError:
        print ('msd_umd.py -f <UMD_filename> -z <HorizontalJump> -v <VerticalJump> -b <Ballistic> -m <Mode> -a <Axes>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('Program to compute the Mean Square Displacement. Usage: ')
            print ('msd_umd.py -f <UMD_filename> -z <HorizontalJump> -v <VerticalJump> -b <Ballistic> -m <Mode> -a <Axes>')
            print ('UMD_filename = input file with the trajectory, in UMD format.')
            print (' As the MSD is measured with a sliding window of various size up to half the trajectory\'s length, ')
            print (' we can accelerate the calculation using a selected reduced sampling ')
            print ('-z : HorizontalJump. discretization for the start of the sampling window.')
            print ('-v :VerticalJump. discretization for the length of the sampling window.')
            print ('-b : Ballistic. estimation of the ballistic part of the trajectory. Default is 0. Typical values of 100 are sufficient.')
            print ("-m : mode. 'atoms' or 'elements' (default) ; to either print the msd of each individual atom (+ the mean msd of each element) or only the msd of the elements.")
            print ('-a : Axes. a list of float that describes the axes along which the msd should be calculated (they can be more than 3). If this argument is absent, the default will be the traditional cartesian coordinates, and the execution will be slightly faster. Write < -a Auto > to use the axes as they are defined by the initial crystal lattice (from the UMD file).')
            print ('Example of use : -a [1,0,0,0,1,0,0,1,1] will calculate the msd along the axes [1,0,0], [0,1,0] and [0,1,1].')
            sys.exit()
        elif opt in ("-f", "--fumdfile"):
            umdfile = str(arg)
            print ('I will use the ',umdfile,' for input','\n')
        elif opt in ("-z", "--zHorizontalJump"):
            hh = int(arg)
            if hh<=0 :
                print('ERROR : horizontal jump -z has to be strictly positive.')
                sys.exit()
        elif opt in ("-v", "--vVerticalJump"):
            vv = int(arg)
            if vv<=0 :
                print('ERROR : vertical jump -v has to be strictly positive.')                
                sys.exit()
        elif opt in ("-b", "--bBallistic"):
            ballistic = int(arg)
        elif opt in ("-m","--mMode"):
            Mode=str(arg)
        elif opt in ("-a","--aAxes"):
            if arg == "Auto":
                Auto =True
                nAxes = 3
            else :
                Axes = eval(arg)
            
                nAxes = int(len(Axes)/3)
                for i in range(nAxes):
                    li = math.sqrt(Axes[3*i]**2+Axes[3*i+1]**2+Axes[3*i+2]**2)
                    Axes[3*i]/=li
                    Axes[3*i+1]/=li
                    Axes[3*i+2]/=li
        elif opt in ("-k","--nCores"):
            nCores = int(arg)

            
    if (os.path.isfile(umdfile)):
        
        MyCrystal,SnapshotsValues,TimeStep,numsteps = umdpf.read_values(umdfile,"absxcart","chunks")

        if Auto :#Creation of the axes as the lattice vectors
            Axes = MyCrystal.rprim[0] + MyCrystal.rprim[1] + MyCrystal.rprim[2]    

        
        ArrayAbsXcart = np.array(SnapshotsValues).flatten()

        ListAbsXcart = [ArrayAbsXcart[at::MyCrystal.natom] for at in range(MyCrystal.natom)]#Each sub-list pertains to an atom and contains the list of x,y,z coordinates for each snapshot, for this atom

    
        print ('Number of atoms of each type is ',MyCrystal.types)
        Instants = []
        #Separation of the algorithm's path in accordance with the parameters who have been provided by the user
        
        if Axes !=None :
            msdAtomRed=partial(msdAtom_CAxes_multi,numsteps=numsteps,hh=hh,vv=vv,ballistic=ballistic,Axes=Axes,nAxes=nAxes,natom=MyCrystal.natom-1)
        else :
            msdAtomRed=partial(msdAtom_C,numsteps=numsteps,hh=hh,vv=vv,ballistic=ballistic,natom=MyCrystal.natom-1)

        if nCores!=None:
            with concurrent.futures.ProcessPoolExecutor(max_workers=nCores) as executor:
                msdArray=list(executor.map(msdAtomRed,[iatom for iatom in range(MyCrystal.natom)],[ListAbsXcart[iatom] for iatom in range(MyCrystal.natom)]))
        else :
            with concurrent.futures.ProcessPoolExecutor() as executor:                  
                msdArray=list(executor.map(msdAtomRed,[iatom for iatom in range(MyCrystal.natom)],[ListAbsXcart[iatom] for iatom in range(MyCrystal.natom)]))#MSD at each time, atom by atom


        msdfile = umdfile[:-8] + '.msd'
                              
        print("\n")
        headerstring='MSD : -f '+umdfile+' -m '+ Mode +'-a '+str(Axes)+'\n'
        if(Mode=="elements"):
            #print("enter elements")
            if Axes != None :
                msdfile+='.axes'
                MSD=np.array([[np.zeros(len(msdArray[0][0])) for _ in range(nAxes+1)] for _ in range(MyCrystal.ntypat)])
                for iatom in range(len(msdArray)) :
                    for kk in range(nAxes):
                        MSD[MyCrystal.typat[iatom]][1+kk]+=msdArray[iatom][kk]#Sum of MSDs atomic type by atomic type, axis by axis
                        MSD[MyCrystal.typat[iatom]][0]+=msdArray[iatom][kk]#Sum of MSDs of all axes
            else :
                MSD=np.array([[np.zeros(len(msdArray[0]))] for _ in range(MyCrystal.ntypat)])#If there's no axes to take in account, we just sum the MSD type by type
                for iatom in range(len(msdArray)) :
                    MSD[MyCrystal.typat[iatom]][0]+=msdArray[iatom]
                    
            msdfile+='.elements.dat'
            f = open(msdfile,'w')
            f.write(headerstring)
            string='time (fs)\t'
            if Axes!=None :
                for itypat in range(MyCrystal.ntypat):
                    string +=  MyCrystal.elements[itypat]+'\t'
                    for i in range(nAxes):
                        string +=  MyCrystal.elements[itypat] +"("+ str(i)+")"+'\t'
            else :
                for itypat in range(MyCrystal.ntypat):
                    string += MyCrystal.elements[itypat] + '\t'
            string = string + '\n'
            f.write(string)
            weight=int((int(numsteps/2)-ballistic)/hh)-1#number of time steps in MSD
            for ii in range(len(MSD[0][0])):
                instant = (float(ii)*TimeStep*vv)+ballistic#Time of occurence reconstructed
                Instants.append(instant)
                string = str(instant)
                for jj in range(MyCrystal.ntypat):#We then write the list of (normalized) MSD(s) for each element at this given time
                    string+='\t' + str(MSD[jj][0][ii]/(float(MyCrystal.types[jj])*float(weight)))
                    if Axes != None :
                        for kk in range(1,nAxes+1):
                            string += '\t' + str(MSD[jj][kk][ii]/(float(MyCrystal.types[jj])*float(weight)))
                string = string + '\n'
                f.write(string)
            print ('MSDs for elements printed in file ',msdfile)
           
        elif(Mode=="atoms"):#Same thing as above, but we print the MSD of individual atoms.
            #print("enter atoms")
            if Axes != None :
                msdfile+=".axes"
                MSD=np.array([[np.zeros(len(msdArray[0][0])) for _ in range(nAxes+1)] for _ in range(MyCrystal.natom)])
                for iatom in range(len(msdArray)) :
                    for kk in range(nAxes):
                        MSD[iatom][1+kk]=msdArray[iatom][kk]
                        MSD[iatom][0]+=msdArray[iatom][kk]
            else :
                MSD=np.array([[np.zeros(len(msdArray[0]))] for _ in range(MyCrystal.natom)])
                for iatom in range(len(msdArray)) :
                    MSD[iatom][0] = msdArray[iatom]
            msdfile+='.atoms.dat'
            f = open(msdfile,'w')
            f.write(headerstring)
            string='time_(fs)\t'
            lastEl=''
            for iatom in range(MyCrystal.natom):
                if lastEl!=MyCrystal.elements[MyCrystal.typat[iatom]]:
                    indexat = 0
                    lastEl=MyCrystal.elements[MyCrystal.typat[iatom]]
                label = lastEl+str(indexat)
                if Axes != None :
                    string +=  label + '\t' + label + '(1)\t' + label + '(2)\t'+ label + '(3)\t'
                else :
                    string += label +'\t'
                
                indexat+=1
                
            string = string + '\n'
            f.write(string)

            weight=int((int(numsteps/2)-ballistic)/hh)-1

            for ii in range(len(MSD[0][0])):     
                instant = (float(ii)*TimeStep*vv)+ballistic
                Instants.append(instant)
                string = str(instant)
                for jj in range(MyCrystal.natom):
                    string+='\t' + str(MSD[jj][0][ii]/(float(weight)))
                    if Axes != None :
                        for kk in range(1,nAxes+1):
                            string += '\t' + str(MSD[jj][kk][ii]/(float(weight)))
                string +='\n'
                f.write(string)
            print ('MSDs for atoms printed in file ',msdfile)

        
        return [msdfile,MSD,Instants,MyCrystal.elements,Axes]
         
    else:
        print ('umd file ',umdfile,'does not exist')
        sys.exit()
    
    end=time.time()
    
    print("runtime :",end-start)
 

if __name__ == "__main__":
    main(sys.argv[1:])
