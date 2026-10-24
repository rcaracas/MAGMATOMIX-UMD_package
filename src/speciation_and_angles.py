#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 09:44:43 2023
"""
###
##AUTHORS: KEVIN JIGUET
##Angles in any cell (NPT), bond lifetimes, times with sampling : Razvan Caracas
###

import platform
import sys,getopt,os.path
import umd_processes_fast as umdpf
import time
from functools import partial
import concurrent.futures
from os.path import join
import ctypes
import numpy as np


OS = platform.system()
if OS == "Linux":
    LibraryName = 'c_clusters.so'
elif OS == "Windows":
    LibraryName = 'c_clusters.dll'
else:
    LibraryName = 'c_clusters.dylib'

clust_lib = ctypes.cdll.LoadLibrary(join(os.path.dirname(os.path.abspath(__file__)),LibraryName))
clust_lib.fullclustering.argtypes = [ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_int),ctypes.c_int,ctypes.c_int,ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_int),ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int]
clust_lib.fullclustering.restype = ctypes.POINTER(ctypes.c_int)

clust_lib.free_memory_int.argtypes = [ctypes.POINTER(ctypes.c_int)]
clust_lib.free_memory_int.restype = None


def in_set(Indexes,natom):#Boolean mask of the atoms designated by Indexes = [first,last, first,last, ...]
    mask = np.zeros(natom,dtype=bool)
    for i in range(len(Indexes)//2):
        mask[Indexes[2*i]:Indexes[2*i+1]+1] = True
    return mask


def bond_angles(X,cell,center,neighbours):#Angles (degrees) neighbour i - center - neighbour j, for i<j in the order of the list ; nearest images in any cell
    d = (X[neighbours]-X[center])@np.linalg.inv(cell)
    d -= np.rint(d)
    v = d@cell
    i,j = np.triu_indices(len(neighbours),1)
    cos = np.sum(v[i]*v[j],axis=1)/(np.linalg.norm(v[i],axis=1)*np.linalg.norm(v[j],axis=1))
    return np.degrees(np.arccos(np.clip(cos,-1.0,1.0)))


def bond_lives(Bonds,BondsIndexes,MyCrystal,CentIndexes,AdjIndexes):#Lifetimes, in snapshots, of the bonds central-adjacent, by pair of elements
    natom = MyCrystal.natom
    central = in_set(CentIndexes,natom)
    adjacent = in_set(AdjIndexes,natom)
    lives = {}
    def record(bond,n):
        name = MyCrystal.elements[MyCrystal.typat[bond[0]]]+"-"+MyCrystal.elements[MyCrystal.typat[bond[1]]]
        lives.setdefault(name,[]).append(n)
    current = {}#bond -> number of consecutive snapshots up to the current one
    for SnapBonds,SnapIndexes in zip(Bonds,BondsIndexes):
        bonds = {}
        for at in np.nonzero(central)[0]:
            for ne in SnapBonds[SnapIndexes[at]:SnapIndexes[at+1]]:
                if adjacent[ne] and not (central[ne] and adjacent[at] and ne < at):#each bond once, also when both atoms are central and adjacent
                    bonds[(at,ne)] = current.get((at,ne),0)+1
        for bond,n in current.items():
            if bond not in bonds:
                record(bond,n)
        current = bonds
    for bond,n in current.items():#bonds alive at the end
        record(bond,n)
    return lives


def write_bond_lives(FileName,lives,TimeStep,BondFile,Centrals,Adjacents):
    with open(FileName,"w") as fb:
        fb.write("Bonds lifetime from file "+BondFile+"\tCentral elements : "+", ".join(Centrals)+"  Coordinating elements : "+", ".join(Adjacents)+"\n\n")
        for name in lives :
            counts = np.bincount(lives[name])
            ntot = counts.sum()
            fb.write("Bond : "+name+"\n"+"lifetime (fs)\tquantity\tpercentage\n")
            for n in range(1,len(counts)):
                fb.write(str(n*TimeStep)+"\t"+str(counts[n])+"\t"+str(counts[n]/ntot*100)+"\n")
            fb.write("\n")
    print("Bonds lifetimes written in file ",FileName)


def write_angles(FileName,header,Bonds,BondsIndexes,Snapshots,cells,MyCrystal,CentIndexes,AdjIndexes,TimeStep,Nsteps):#All the angles adjacent - central - adjacent of every snapshot
    adjacent = in_set(AdjIndexes,MyCrystal.natom)
    centrals = np.nonzero(in_set(CentIndexes,MyCrystal.natom))[0]
    TotSum = 0.0
    TotN = 0
    with open(FileName,'w') as fa:
        fa.write(header)
        for step in range(len(Bonds)) :
            SnapBonds,SnapIndexes = Bonds[step],BondsIndexes[step]
            X = np.array(Snapshots[step],dtype=float).reshape(-1,3)
            fa.write('\nstep : '+str(step*Nsteps)+'\ntime :'+str(step*TimeStep)+" fs \n")
            line = 'Angles : \n'
            SnapshotAngles = []
            for at in centrals :
                neighbours = [ne for ne in SnapBonds[SnapIndexes[at]:SnapIndexes[at+1]] if adjacent[ne]]
                if len(neighbours) > 1 :
                    angles = bond_angles(X,cells[step],at,neighbours)
                    i,j = np.triu_indices(len(neighbours),1)
                    for k in range(len(angles)):
                        line += str(float(angles[k]))+"\t("+str(neighbours[i[k]])+"-"+str(at)+"-"+str(neighbours[j[k]])+") \n"
                    SnapshotAngles += list(angles)
            fa.write(line)
            if SnapshotAngles :
                fa.write("mean : "+str(float(np.mean(SnapshotAngles)))+"\tstdev : "+str(float(np.std(SnapshotAngles)))+"\n")
                TotSum += sum(SnapshotAngles)
                TotN += len(SnapshotAngles)
            fa.write("end of step\n")
        fa.write("\n\n Mean of all angles : "+str(TotSum/TotN if TotN else 0.0))
    print("Angles written in file ",FileName)


def analysis_subtab(Clusters,Step):#Creates a dictionnary whose keys are the individual clusters, and the values are a list of list of their *consecutive* steps of existence.
    population={}
    for Snapshot in Clusters :
        for cluster in Snapshot :
            if str(cluster) in population.keys():
                if population[str(cluster)][-1][-1]==Step-1:
                    population[str(cluster)][-1].append(Step)
                else :
                    population[str(cluster)].append([Step])
            else :
                population[str(cluster)]=[[Step]]
        Step+=1
    return population


def clustering(SnapshotBonds,SnapshotBondIndexes,SnapshotXCart,cell,step,maxSteps,natom,nAts,CentIndexes,OutIndexes,r,AngleCalc):#Create the species from the Bond file

    if(maxSteps > 0 and 100*step//maxSteps != 100*(step-1)//maxSteps):
        sys.stdout.write("\rBuilding species. Progress : "+str(100*step//maxSteps)+"%")
        sys.stdout.flush()

    MaxNrBoundAtoms=SnapshotBondIndexes[-1]#maximum number of bound atoms for any atoms
    #Preparing data for the C script

    SBp = (ctypes.c_int * len(SnapshotBonds))(*SnapshotBonds)
    BIp = (ctypes.c_int * (len(SnapshotBondIndexes)-1))(*SnapshotBondIndexes[:-1])
    CIp = (ctypes.c_int * len(CentIndexes))(*CentIndexes)
    OIp = (ctypes.c_int * len(OutIndexes))(*OutIndexes)

    Np = clust_lib.fullclustering(SBp,BIp,natom,nAts,CIp,OIp,int(len(CentIndexes)/2),int(len(OutIndexes)/2),MaxNrBoundAtoms,r)#Computes the clusters
    Clusters = [[]]
    length = Np[0]

    #Converting C data into a python list of clusters, separated by -1 ; with r = 1 the central atom comes first
    for i in range(1, length + 1):
        atom=Np[i]
        if atom ==-1:
            Clusters.append([])
        else:
            Clusters[-1].append(atom)
    clust_lib.free_memory_int(Np)#Freeing memory

    if Clusters ==[[]]:
        return [],{}

    Angles = {}#Angles within each cluster, keyed by the sorted cluster
    if r==1 and AngleCalc:
        X = np.array(SnapshotXCart,dtype=float).reshape(-1,3)
        for cluster in Clusters :
            if len(cluster) > 2 :
                Angles[str(sorted(cluster))] = [float(a) for a in bond_angles(X,cell,cluster[0],cluster[1:])]

    return [sorted(c) for c in Clusters],Angles


def is_in(at,Indexes):
    for i in range(int(len(Indexes)/2)):
        if (at<=Indexes[2*i+1] and at>=Indexes[2*i]):
            return True
    return False


def main(argv):
    umdpf.headerumd()
    BondFile=''
    UMDFile=''
    Nsteps = 1
    Centrals = []
    Adjacents = []
    minlife = 1
    rings = 1
    header = ''
    start=time.time()
    flagangles = 0
    flagpopul = 1
    flagbondlife = 0
    nCores = None
    usage = 'speciation_and_angles.py -b <bond_filename> -f <umd_filename> -s <Sampling_Frequency> -c <Centrals> -a <Adjacents> -m <MinLife> -r <Rings> -p <Population> -t <Angles> -l <Bondslife> -k <nCores>'
    try:
        opts, arg = getopt.getopt(argv,"hb:f:s:c:a:m:r:t:l:p:k:",["bBondFile","fUMDFile","sSampling_Frequency", "cCentral","aAdjacent","mMinlife","rRings","tAngles","lBondslife","pPopulation","knCores"])
    except getopt.GetoptError:
        print (usage)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('speciation_and_angles.py program to identify speciation, calculate angles, and evaluate the lifetime of bonds')
            print (usage)
            print ('default values: -s 1 -m 1 -r 1 -p 1 -t 0 -l 0')
            print ('-b : interatomic bonding file, computed with bonding_umd.py.')
            print ('-c and -a : central and adjacent elements respectively, separated by commas. If one is \'all\', every atom will be taken in account for this role.')
            print ('-r : rings = 1 default, all adjacent atoms bind to central atoms ; rings = 0, polymerization, all adjacent atoms bind to central AND other adjacent atoms ; rings = x>0, all adjacent atoms bind to central then to other adjacent atoms to form a xth-coordination polyhedra')
            print ('-s : one snapshot of the bonding file out of s is used. Default 1.')
            print ('-m : minimal duration of existence for a chemical species to be taken into account (fs) ; default 1')
            print ('-p : 1 (default) or 0. Chemical species, written in <base>.r<r>.popul.dat (each cluster with its beginning, end and lifetime), <base>.r<r>.stat.dat')
            print ('     (total time and proportion of each species) and <base>.r<r>.step.dat (the species of every snapshot).')
            print ('-t : 0 (default) or 1. The angles adjacent - central - adjacent of every central atom in every snapshot, written in <base>.angles.dat ;')
            print ('     with -r 1 and -p 1, the mean angles of each cluster are also added to the popul and stat files. Needs the umd file (-f).')
            print ('     Nearest images in any cell, with the cell of each snapshot (NPT).')
            print ('-l : 0 (default) or 1. Lifetimes of the bonds central - adjacent, by pair of elements, written in <base>.bondslives.dat.')
            print ('-f : umd file, needed for the angles. The output files are named after it if given, otherwise after the bonding file.')
            print ('-k : number of cores. Default is all.')
            sys.exit()
        elif opt in ("-b", "--bBondFile"):
            BondFile = str(arg)
            header = header + 'FILE: -f=' + BondFile
        elif opt in ("-s","--sNsteps"):
            Nsteps = int(arg)
            header = header + ' -s=' + arg
            print('Will sample the MD trajectory every ',Nsteps,' steps')
        elif opt in ("-m","--mMinlife"):
            minlife = float(arg)
        elif opt in ("-f","--fUMDFile"):
            UMDFile = str(arg)
        elif opt in ("-c","--Central"):
            header = header + ' -c=' + arg
            Centrals = arg.split(",")
        elif opt in ("-a","--Adjacent"):
            header = header + ' -a=' + arg
            Adjacents = arg.split(",")
        elif opt in("-t","--tAngles"):
            flagangles = int(arg)
        elif opt in("-p","--pPopulation"):
            flagpopul = int(arg)
        elif opt in("-l","--lBondslife"):
            flagbondlife = int(arg)
        elif opt in ("-r","--rRings"):
            rings = int(arg)
            if rings == 0:
                print ('Calculation of polymerized coordination polyhedra')
            elif rings > 0:
                print ('Calculation of order '+str(arg)+" coordination sphere")
            else :
                print ('Undefined calculation')
                print ('-r should be a positive integer')
                sys.exit()
        elif opt in("-k","--kCores"):
            nCores = int(arg)

    header = header + ' -r=' + str(rings)

    if not (os.path.isfile(BondFile)):
        print ('the bond file ',BondFile,' does not exist')
        sys.exit()

    if Centrals == [] or Adjacents == []:
        print ('give the central (-c) and adjacent (-a) elements')
        sys.exit()

    if flagangles==1 and (not (os.path.isfile(UMDFile))):
        print ('the UMD file ',UMDFile,' does not exist')
        sys.exit()

    sys.stdout.write(str(os.cpu_count())+" cores are available. ")
    if nCores != None :
        sys.stdout.write("We will use "+str(nCores)+" of them.\n")
    else :
        sys.stdout.write("The number that will be used is automatically determined.\n")

    #The output files are named after the UMD file if given (-f is optional), otherwise after the bonding file
    if UMDFile != '':
        BaseName = UMDFile[:-8] if UMDFile.endswith('.umd.dat') else os.path.splitext(UMDFile)[0]
    else:
        BaseName = BondFile[:-12] if BondFile.endswith('.bonding.dat') else os.path.splitext(BondFile)[0]

    #TimeStep : time between two of the snapshots used (the timestep of the bonding file times Nsteps)
    CentIndexes,AdjIndexes,MyCrystal,[Bonds,BondsIndexes],TimeStep = umdpf.read_bonds(BondFile,Centrals,Adjacents,Nsteps,nCores=nCores)

    if CentIndexes == -1 :#For the GUI to display an error message whenever an element is missing
        return AdjIndexes

    if flagbondlife==1 :
        print("Calculating the lifetimes of the bonds...")
        lives = bond_lives(Bonds,BondsIndexes,MyCrystal,CentIndexes,AdjIndexes)
        write_bond_lives(BaseName+".bondslives.dat",lives,TimeStep,BondFile,Centrals,Adjacents)

    if flagangles==1 :#positions and cells of the UMD snapshots matching the bonding snapshots used
        MyCrystalUMD,TimeStepUMD = umdpf.Crystallization(UMDFile)
        TimeRatio = max(1,int(round(TimeStep/TimeStepUMD)))
        MyCrystalUMD,SnapshotsXCart,TimeStepUMD,length = umdpf.read_values(UMDFile,"xcart","line",TimeRatio,nCores=nCores)
        nsnap = min(len(Bonds),len(SnapshotsXCart))
        if nsnap < len(Bonds):
            print('WARNING : the UMD file has fewer snapshots than the bonding file : only the first ',nsnap,' are used.')
            Bonds,BondsIndexes = Bonds[:nsnap],BondsIndexes[:nsnap]
        SnapshotsXCart = SnapshotsXCart[:nsnap]
        cells = umdpf.read_cells(UMDFile,nsnap,TimeRatio)
        write_angles(BaseName+".angles.dat",header+"\n",Bonds,BondsIndexes,SnapshotsXCart,cells,MyCrystal,CentIndexes,AdjIndexes,TimeStep,Nsteps)
    else :
        SnapshotsXCart = [[] for _ in range(len(Bonds))]#Blank lists which will not be used but are needed as arguments
        cells = [None for _ in range(len(Bonds))]

    if flagpopul==0 :
        print("total duration : ",time.time()-start, " s")
        return True#For the GUI

    AllElements = []
    nAts = 0
    if Centrals != ["all"]:
        for el in Centrals :
            if el not in AllElements :
                AllElements.append(el)
                nAts += MyCrystal.types[MyCrystal.elements.index(el)]
    else :
        AllElements = MyCrystal.elements
        nAts = MyCrystal.natom

    if Adjacents != ["all"]:
        for el in Adjacents :
            if el not in AllElements :
                AllElements.append(el)
                nAts += MyCrystal.types[MyCrystal.elements.index(el)]
    else :
        AllElements = MyCrystal.elements
        nAts = MyCrystal.natom

    clusteringRed=partial(clustering,maxSteps=len(Bonds)-1,natom=MyCrystal.natom,nAts=nAts,CentIndexes=CentIndexes,OutIndexes=AdjIndexes,r=rings,AngleCalc=flagangles)

    with concurrent.futures.ProcessPoolExecutor(max_workers = nCores) as executor :
        Data=list(executor.map(clusteringRed,Bonds,BondsIndexes,SnapshotsXCart,cells,[step for step in range(len(Bonds))])) #Computes the clusters of atoms for each snapshot separately

    clusters,Angles = map(list,zip(*Data))

    population=analysis_subtab(clusters,0)
    FileAll = BaseName +'.r' + str(rings) + '.popul.dat'
    FileStat = BaseName + '.r' + str(rings) + '.stat.dat'
    FileStep = BaseName + '.r' + str(rings) + '.step.dat'
    header+="\n"

    print("\nWriting...")#writing the clusters in a file

    fs = open(FileStep,'w')
    fs.write(header)
    headstring = "Clusters\tNumber of atoms\tComposition\n"
    for step in range(len(clusters)) :
        st = "step\t"+str(step*Nsteps)+"\ntime\t"+str(step*TimeStep)+"\n"
        fs.write(st)
        fs.write(headstring)
        Clusts = clusters[step]
        vapor = [at for at in range(MyCrystal.natom)]#The lonely atoms are considered to be vapor, and the others will be gradually pruned from this list

        for clust in Clusts :
            index=[0 for _ in range(MyCrystal.ntypat)]
            name=''
            for at in clust:
                index[MyCrystal.typat[at]]+=1

            for elem in range(len(index)) :
                if index[elem]!=0:
                    name+=MyCrystal.elements[elem]+'_'+str(index[elem])
            newstring = name +'\t'+ str(len(clust)) +'\t'+ str(clust)+'\n'
            fs.write(newstring)
            for at in clust :
                vapor[at]=-1
        for at in vapor :#The remaining not -1 atoms are considered to be vapors
            if at !=-1 and (is_in(at,CentIndexes) or is_in(at,AdjIndexes)) :
                fs.write(MyCrystal.elements[MyCrystal.typat[at]]+'_1\t1\t['+str(at)+']\n')

    fs.close()

    dicoStats={}#Clusters by name (chemical species) ; values as sub-dictionaries whose keys are 'lifetime' (combined life time of the molecules of this specie) and '#atoms' (number of atoms of this specie)
    dicoTimes={}#Clusters (individual atoms) by step of appearance. values as list whose items are under the form : [individual atoms,name of the chemical specie,last step, life time]
    dicoMeanAngle={}#Mean angle, cluster by cluster. The keys are the names of clusters (individual atoms), the valuee are under the form [sum of (mean of angles for each instance of the cluser),number of instances]
    total=0
    for key in population :
        index=[0 for _ in range(MyCrystal.ntypat)]
        atomsCluster=key.strip('[]').split(', ')
        name=''
        for at in atomsCluster:
            index[MyCrystal.typat[int(at)]]+=1

        for elem in range(len(index)) :
            if index[elem]!=0:
                name+=MyCrystal.elements[elem]+'_'+str(index[elem])

        for life in population[key]:
            lifetime = TimeStep*(life[-1]-life[0]+1)#each snapshot used stands for TimeStep
            if lifetime>minlife :
                if name in dicoStats:
                    dicoStats[name]['lifetime']+=lifetime
                else :
                    dicoStats[name]={'lifetime':lifetime, '#atoms': len(atomsCluster)}
                total+=lifetime

                if life[0] in dicoTimes :
                    dicoTimes[life[0]].append([key,name,life[-1],lifetime])
                else:
                    dicoTimes[life[0]]=[[key,name,life[-1],lifetime]]

    if rings==1 and flagangles:
        newstring="Formula\tBegin (step)\tEnd(step)\tlifetime (fs)\tcomposition\tAngles\n"
    else :
        newstring="Formula\tBegin (step)\tEnd(step)\tlifetime (fs)\tcomposition\n"

    fa = open(FileAll,'w')
    fa.write(header)
    fa.write(newstring)
    for ii in range(len(Bonds)):
        if ii in dicoTimes:
            for data in dicoTimes[ii]:
                if data[3]>minlife :
                    newstring=data[1]+"\t"+str(ii*Nsteps)+"\t"+str(data[2]*Nsteps)+"\t"+str(data[3])+"\t"+data[0]

                    if rings==1 and flagangles :
                        if data[0] in Angles[ii].keys():
                            angles=np.array([np.array(Angles[jj][data[0]]) for jj in range(ii,data[2]+1)])
                            Meanangles = np.sum(angles,0)/len(angles)
                            totmean = np.sum(Meanangles,0)/len(Meanangles)

                            if data[1] in dicoMeanAngle:
                                dicoMeanAngle[data[1]][0]+=totmean
                                dicoMeanAngle[data[1]][1]+=1
                            else:
                                dicoMeanAngle[data[1]]=[totmean,1]
                        else:
                            Meanangles=np.array([])

                        for alpha in Meanangles:
                            newstring+="\t"+str(alpha)

                    newstring+="\n"
                    fa.write(newstring)
    fa.close()


    if(rings==1 and flagangles):
        newstring="Cluster\tTime(fs)\tPercent\tNumber of atoms\tMean Angles\n"
    else:
        newstring="Cluster\tTime(fs)\tPercent\tNumber of atoms\n"

    fa=open(FileStat,'w')
    fa.write(header)
    fa.write(newstring)
    for clust in dicoStats :
        newstring=clust+"  \t"+str(dicoStats[clust]['lifetime'])+"\t"+str(float(dicoStats[clust]['lifetime'])/float(total))+"\t"+str(dicoStats[clust]['#atoms'])
        if rings == 1 and clust in dicoMeanAngle:
            newstring+='\t'+str(dicoMeanAngle[clust][0]/dicoMeanAngle[clust][1])
        newstring+='\n'
        fa.write(newstring)
    fa.close()

    print ('Population written in file :  ',FileAll)
    print ('Statistics written in file :  ',FileStat)
    print ('Species for each snapshot written in file :',FileStep)

    end=time.time()
    print("total duration : ",end-start, " s")
    return True

if __name__ == "__main__":
    main(sys.argv[1:])
