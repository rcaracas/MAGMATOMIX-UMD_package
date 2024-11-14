#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 09:44:43 2023
"""
#!/usr/bin/env python3
###
##AUTHORS: KEVIN JIGUET
###

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


current_path=os.path.abspath(__file__)#For all this to work, the file < c_clusters.so > must be in the same directory than gofr_umd
path_split=current_path.split('/')
path_red=path_split[1:-1]
path_new=''
for u in path_red:
    path_new+='/'+u

OS = platform.system()
LibraryName =""


if OS == "Linux":
    LibraryName = 'c_clusters.so'

elif OS == "Windows":
    LibraryName = 'c_clusters.dll'

elif OS == "Darwin":
    LibraryName = 'c_clusters.dylib'
    


clust_lib = ctypes.cdll.LoadLibrary(join(path_new, LibraryName))
clust_lib.fullclustering.argtypes = [ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_int),ctypes.c_int,ctypes.c_int,ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_int),ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int]
clust_lib.fullclustering.restype = ctypes.POINTER(ctypes.c_int)

clust_lib.angles.argtypes = [ctypes.POINTER(ctypes.c_int),ctypes.c_int,ctypes.c_int,ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_int),ctypes.c_int,ctypes.POINTER(ctypes.c_int),ctypes.c_int,ctypes.c_double,ctypes.c_double,ctypes.c_double]
clust_lib.angles.restype = ctypes.POINTER(ctypes.c_double)

clust_lib.free_memory_int.argtypes = [ctypes.POINTER(ctypes.c_int)]
clust_lib.free_memory_int.restype = None

clust_lib.free_memory_double.argtypes = [ctypes.POINTER(ctypes.c_double)]
clust_lib.free_memory_double.restype = None

def analysis_BondLives(BondFile,Centrals,Adjacents,Nsteps,nCores):#Creates a dictionary whose entries are the atoms making a given angle, and the values are lists of initial and final times for the lives of this angle
    
    CentIndexes,AdjIndexes,MyCrystal,Bonds,TimeStep = umdpf.read_bonds(BondFile,Centrals,Adjacents,Nsteps,mode = "Dictionary",nCores=nCores)

    DicoBonds = {}
    for centInd in CentIndexes[::2] :
        for adjInd in AdjIndexes[::2] :
            centEl = MyCrystal.elements[MyCrystal.typat[centInd]]
            adjEl = MyCrystal.elements[MyCrystal.typat[adjInd]]
            key = centEl+"-"+adjEl
            DicoBonds[key] = []

    DicoIndBonds = {}
    for step in range(len(Bonds)) :
        SnapshotBonds = Bonds[step]
        for el in range(len(CentIndexes)//2):
            for centAt in range(CentIndexes[2*el],CentIndexes[2*el+1]+1):
                if centAt in SnapshotBonds :
                    bonds = SnapshotBonds[centAt]
                    for adjAt in bonds : 
                        bond = [centAt,adjAt]
                        if str(bond) in DicoIndBonds.keys():
                            if DicoIndBonds[str(bond)][-1][-1] == step -1:
                                DicoIndBonds[str(bond)][-1][-1] +=1
                            else :
                                DicoIndBonds[str(bond)].append([step,step])
                        else :
                            DicoIndBonds[str(bond)]=[[step,step]]
 
    for key in DicoIndBonds :
        bond = eval(key)
        name = MyCrystal.elements[MyCrystal.typat[bond[0]]]+"-"+MyCrystal.elements[MyCrystal.typat[bond[1]]]
        listLives = []
        times = DicoIndBonds[key]
        for life in times :
            listLives.append((life[1]-life[0]+1))
        DicoBonds[name]+=listLives

    return DicoBonds,TimeStep
    

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
    
def compute_angles(Bonds,MyCrystal,AllSnapshots,CentMin,CentMax,OutMin,OutMax,fa,TimeStep,Nsteps):#Calculates and writes the angles in a file.
    acell=MyCrystal.acell
    TotMean=0
    TotNangles=0
    for step in range(len(Bonds)) :
        SnapshotAngles=[]
        SnapshotBonds=Bonds[step]
        MySnapshot=AllSnapshots[step]
        fa.write('\nstep : '+str(step*Nsteps)+'\ntime :'+str(step*TimeStep*Nsteps)+" fs \n")
        line='Angles : \n'
        Nangles=0
        for internAt in SnapshotBonds :
            if internAt<=CentMax and internAt>=CentMin:
                for iexternAt in SnapshotBonds[internAt]:
                    if iexternAt<=OutMax and iexternAt>=OutMin:
                        for jexternAt in SnapshotBonds[internAt]:
                            if jexternAt<=OutMax and jexternAt>=OutMin:
                                if iexternAt<jexternAt :    
                                    Nangles+=1
                        
                                    xEi,yEi,zEi=MySnapshot[iexternAt][0],MySnapshot[iexternAt][1],MySnapshot[iexternAt][2]
                                    xEj,yEj,zEj=MySnapshot[jexternAt][0],MySnapshot[jexternAt][1],MySnapshot[jexternAt][2]
                                    xI,yI,zI=MySnapshot[internAt][0],MySnapshot[internAt][1],MySnapshot[internAt][2]
                        
                                    dx = xEi-xEj
                                    dy = yEi-yEj
                                    dz = zEi-zEj
                        
                                    valx = min(dx**2, (acell[0]-dx)**2, (acell[0]+dx)**2)
                                    valy = min(dy**2, (acell[1]-dy)**2, (acell[1]+dy)**2)
                                    valz = min(dz**2, (acell[2]-dz)**2, (acell[2]+dz)**2)
                                    distEiEj = valx + valy + valz               
                        
                                    dx = xI-xEj
                                    dy = yI-yEj
                                    dz = zI-zEj                        
                        
                                    valx = min(dx**2, (acell[0]-dx)**2, (acell[0]+dx)**2)
                                    valy = min(dy**2, (acell[1]-dy)**2, (acell[1]+dy)**2)
                                    valz = min(dz**2, (acell[2]-dz)**2, (acell[2]+dz)**2)
                                    distIEj = valx + valy + valz               
                        
                                    dx = xEi-xI
                                    dy = yEi-yI
                                    dz = zEi-zI                        
                        
                                    valx = min(dx**2, (acell[0]-dx)**2, (acell[0]+dx)**2)
                                    valy = min(dy**2, (acell[1]-dy)**2, (acell[1]+dy)**2)
                                    valz = min(dz**2, (acell[2]-dz)**2, (acell[2]+dz)**2)
                                    distIEi = valx + valy + valz               
        
                                    angle=math.acos((distIEi+distIEj-distEiEj)/(2*math.sqrt(distIEi*distIEj)))/math.pi*180
                                    
                                    line+=str(angle)+"\t("+str(iexternAt)+"-"+str(internAt)+"-"+str(jexternAt)+") \n"
                        
                                    SnapshotAngles.append(angle)
        fa.write(line)

        TotNangles+=Nangles
        if Nangles!=0:
            mean=sum(SnapshotAngles)/Nangles                  
            msd=0
            TotMean+=sum(SnapshotAngles)
            for angle in SnapshotAngles:
                msd+=(angle-mean)**2
            msd/=mean
            msd=math.sqrt(msd)
        
            fa.write("mean : "+str(mean)+"\tmsd : "+str(msd)+"\n")
        fa.write("end of step\n")

    if(TotNangles!=0):
        TotMean/=TotNangles

    fa.write("\n\n Mean of all angles : "+str(TotMean))
    fa.close()
        
def clustering(SnapshotBonds,SnapshotBondIndexes,SnapshotXCart,step,maxSteps,natom,nAts,CentIndexes,OutIndexes,acell,r,AngleCalc):#Create the species from the Bond file
    
    if(100*step//maxSteps != 100*(step-1)//maxSteps):
        sys.stdout.write("\rBuilding species. Progress : "+str(100*step//maxSteps)+"%")        
        sys.stdout.flush()
        
    M=SnapshotBondIndexes[-1]#maximum number of bound atoms for any atoms
    #Preparing data for the C script        

    SBp = (ctypes.c_int * len(SnapshotBonds))(*SnapshotBonds)
    BIp = (ctypes.c_int * (len(SnapshotBondIndexes)-1))(*SnapshotBondIndexes[:-1])
    CIp = (ctypes.c_int * len(CentIndexes))(*CentIndexes)
    OIp = (ctypes.c_int * len(OutIndexes))(*OutIndexes)
    
    Np = clust_lib.fullclustering(SBp,BIp,natom,nAts,CIp,OIp,int(len(CentIndexes)/2),int(len(OutIndexes)/2),M,r)#Computes the clusters
    Clusters = [[]]
    Angles = {}
    length = Np[0]
    #Converting C data into a python list of clusters
    for i in range(1,length+1):
        atom=Np[i]
        if atom ==-1 :
            Clusters.append([])
        else :   
            Clusters[-1].append(atom)

    
    if r==1 and AngleCalc:#Calculating the angles within each cluster, as a dictionary whose entries are clusters and values are the angles
        SXp = (ctypes.c_double * len(SnapshotXCart))(*SnapshotXCart)
        Ap = clust_lib.angles(Np,len(Clusters),M,SXp,CIp,int(len(CentIndexes)/2),OIp,int(len(OutIndexes)/2),acell[0],acell[1],acell[2])
        index=0
        flagnew=1
        for i in range(1,int(Ap[0])):
            if(Ap[i]==-1):
                index=index+1
                flagnew=1
            elif flagnew==1:
                flagnew=0
                Angles[str(Clusters[index])]=[Ap[i]]
            else :
                Angles[str(Clusters[index])].append(Ap[i])

        clust_lib.free_memory_double(Ap)#Freeing memory    

    clust_lib.free_memory_int(Np)#Freeing memory

    if Clusters ==[[]]:
        return [],{}
    
    return Clusters,Angles
    

def is_in(at,Indexes):
    for i in range(int(len(Indexes)/2)):
        if (at<=Indexes[2*i+1] and at>=Indexes[2*i]):
            return True
    return False

def main(argv):
    umdpf.headerumd()
    BondFile='bonding.umd.dat'
    UMDFile=''
    Nsteps = 1
    Centrals = []
    Adjacents = []
    minlife = 1
    rings = 1
    header = ''
    start=time.time()
    t = 0
    p = 1
    b = 0
    nCores = None
    try:
        opts, arg = getopt.getopt(argv,"hb:f:s:c:a:m:r:t:l:p:k:",["bBondFile","fUMDFile","sSampling_Frequency", "cCentral","aAdjacent","mMinlife","rRings","pPopulation","tAngles","lBondslife","pPopulation","knCores"])
    except getopt.GetoptError:
        print ('speciation_and_angles.py -b <bond_filename> -f <umd_filename> -s <Sampling_Frequency> -c <Cations> -a <Anions> -m <MinLife> -r <Rings> -p <Population> -t <Angles> -l <Bondslife> -k <nCores>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('speciation_and_angles.py program to identify speciation, calculate angles, and evaluate the lifetime of bonds')
            print ('speciation_and_angles.py -b <bond_filename> -f <umd_filename> -s <Sampling_Frequency> -c <Cations> -a <Anions> -m <MinLife> -r <Rings> -p <Population> -t <Angles> -l <Bondslife> -k <nCores>')
            print ('default values: -f bonding.umd.dat -s 1 -m 5 -r 1')
            print ('the bond file contains the bonds relations for each snapshot, which is computed with bonding_umd.py.')
            print ('-c and -a : central and adjacent elements respectively. If one is \'all\', every atom will be taken in account for this role.')
            print ('-r : rings = 1 default, all adjacent atoms bind to central atoms ; rings = 0, polymerization, all adjacent atoms bind to central AND other adjacent atoms ; rings = x>0, all adjacent atoms bind to central then to other adjacent atoms to form a xth-coordination polyhedra')
            print ('-m : minimal duration of existence for a chemical species to be taken into account (fs) ; default 1')
            print ('-t : nothing (default) or 1. In the latter case, and only if -r = 1,the angles of the molecules will be computed, their summit being the central atom. If this option is activated, the user needs to provide the corresponding umd file.')
            print ('-l : nothing (default) or 1. In the latter case, the life time of each type of bonds will be evaluated and printed in a separate file.')
            print ('-f : umd file used to calculate the angles (optional).')
            print ('-b : interatomic bonding file.')
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
        elif opt in ("-f","fUMDFile"):
            UMDFile = str(arg)
        elif opt in ("-c","--Central"):
            header = header + ' -c=' + arg
            Centrals = arg.split(",")
            #print ('Cation list is: ',Cations)
        elif opt in ("-a","--Adjacent"):
            header = header + ' -a=' + arg
            Adjacents = arg.split(",")
            #print ('Anion list is: ',Anions)
        elif opt in("-t","--tAngles"):
            t = int(arg)
        elif opt in("-p","--pPopulation"):
            p = int(arg)
        elif opt in("-l","--lBondslife"):
            b = int(arg)
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

    if t==1 and (not (os.path.isfile(UMDFile))):
        print ('the UMD file ',UMDFile,' does not exist')            
        sys.exit()
    
    sys.stdout.write(str(os.cpu_count())+" cores are available. ")
    if nCores != None :
        sys.stdout.write("We will use "+str(nCores)+" of them.\n")        
    else :
        sys.stdout.write("The number that will be used is automatically determined.\n")

            
    if b==1 :
        print("Calculating the lifetimes of the bonds...")
        DicoBonds,TimeStep = analysis_BondLives(BondFile,Centrals,Adjacents,Nsteps,nCores)
        DicoLives = {}
        for name in DicoBonds :
            times = DicoBonds[name]
            if times !=[]:
                times.sort()
                timebins=[i+1 for i in range(times[-1])]
                timenums=[0 for i in range(times[-1])]
                for life in times :
                    timenums[life-1]+=1    
                DicoLives[name] = [timebins,timenums]
            else:
                DicoLives[name] = [[],[]]
        FileName = BondFile[:-4]+".bondslives.dat"
        fb = open(FileName,"w")
        centstring = ""
        for el in Centrals[:-1] :
            centstring += el+", "
        centstring+=Centrals[-1]
        outstring = ""
        for el in Adjacents[:-1] :
            outstring += el+", "
        outstring+=Adjacents[-1]
        header = "Bonds lifetime from file "+BondFile+"\tCentral elements : "+centstring+"  Coordinating elements : "+outstring+"\n\n"
        fb.write(header)
        for name in DicoLives :
            [timebins,timenums]=DicoLives[name]
            Bondheader = "Bond : "+name+"\n"+"lifetime (fs)\tquantity\tpercentage\n"
            fb.write(Bondheader)
            ntot = sum(timenums)
            for timebin in timebins :
                line = ""
                line+=str((timebin-1)*TimeStep+1)+"\t"+str(timenums[timebin-1])+"\t"+str(timenums[timebin-1]/ntot*100)+"\n"
                fb.write(line)
            fb.write("\n")
        print("Bonds lifetimes written in file ",FileName)
        fb.close()
 
    if t==0 and p==0 :
        return True#For the GUI
    if t==1 and p==0 :
        if len(Centrals)==1 and len(Adjacents)==1:
            CentIndexes,AdjIndexes,MyCrystal,Bonds,TimeStep = umdpf.read_bonds(BondFile,Centrals,Adjacents,Nsteps,"Dictionary")
            if CentIndexes == -1 :#For the GUI to display an error message whenever an element is missing
                return AdjIndexes
            
            MyCrystalUMD,TimeStepUMD = umdpf.Crystallization(UMDFile)
            TimeRatio = int(TimeStep/TimeStepUMD)#Used to match the snapshots of UMD and bonding file, in case there has been a sampling of the snapshots in between

            (MyCrystalUMD,AllSnapshots,TimeStepUMD,length)=umdpf.read_values(UMDFile,"xcart","lists",Nsteps*TimeRatio)
            
            CentMin = CentIndexes[0]
            CentMax = CentIndexes[1]
            AdjMin = AdjIndexes[0]
            AdjMax = AdjIndexes[1]
        
            File=str(UMDFile)+".angles.dat"
        
            fa=open(File,'w')
            fa.write(header)
            compute_angles(Bonds,MyCrystal,AllSnapshots,CentMin,CentMax,AdjMin,AdjMax,fa,TimeStep,Nsteps)    
        
            print("Angles computed successfully. File created under the name ",File)
            print("total duration : ",(time.time()-start))
            return True
        
        else:
            print("Too many elements to compute the angles only. Regular identification of the chemical species will be proceeded.")
 
    CentIndexes,AdjIndexes,MyCrystal,[Bonds,BondsIndexes],TimeStep = umdpf.read_bonds(BondFile,Centrals,Adjacents,Nsteps,nCores=nCores)

    if CentIndexes == -1 :#For the GUI to display an error message whenever an element is missing
        return AdjIndexes

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
        
        
    if (t==1) and (rings==1):
        MyCrystalUMD,TimeStepUMD = umdpf.Crystallization(UMDFile)
        TimeRatio = int(TimeStep/TimeStepUMD)#Used to match the snapshots of UMD and bonding file, in case there has been a sampling of the snapshots in between
        MyCrystalUMD,SnapshotsXCart,TimeStepUMD,length = umdpf.read_values(UMDFile,"xcart","line",Nsteps*TimeRatio,nCores=nCores)
    else :
        SnapshotsXCart = [[] for _ in range(len(Bonds))]#Blank list which will not be used but is needed as an argument
    
    clusteringRed=partial(clustering,maxSteps=len(Bonds)-1,natom=MyCrystal.natom,nAts=nAts,CentIndexes=CentIndexes,OutIndexes=AdjIndexes,r=rings,acell=MyCrystal.acell,AngleCalc=t)

#    DataII = []
#    for i in range(len(Bonds)):
#        DataII.append(clusteringRed(Bonds[i],BondsIndexes[i],SnapshotsXCart[i],i))
#    print(DataII[0])
#    sys.exit()

    #Separating the algoritm's path accordingly to the user's will
    if nCores != None :
        with concurrent.futures.ProcessPoolExecutor(max_workers = nCores) as executor :
            Data=list(executor.map(clusteringRed,Bonds,BondsIndexes,SnapshotsXCart, [step for step in range(len(Bonds))])) #Computes the clusters of atoms for each snapshot separately
    else : 
        with concurrent.futures.ProcessPoolExecutor() as executor :
            Data=list(executor.map(clusteringRed,Bonds,BondsIndexes,SnapshotsXCart, [step for step in range(len(Bonds))])) #Computes the clusters of atoms for each snapshot separately
            
    clusters,Angles = map(list,zip(*Data))
        
    population=analysis_subtab(clusters,0)
        #Creating the output files        
    FileAll = UMDFile[:-8] +'.r' + str(rings) + '.popul.dat'
    FileStat = UMDFile[:-8] + '.r' + str(rings) + '.stat.dat'
    FileStep = UMDFile[:-8] + '.r' + str(rings) + '.step.dat'
    header+="\n"            
    
    print("\nWriting...")#writing the clusters in a file
    
    fs = open(FileStep,'w')
    fs.write(header)        
    headstring = "Clusters\tNumber of atoms\tComposition\n"
    for step in range(len(clusters)) :
        st = "step\t"+str(step*Nsteps)+"\ntime\t"+str(step*Nsteps*TimeStep)+"\n"
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
            
    dicoNames={}#Clusters by name (chemical species) ; values as list whose items are under the form : [individual atoms,initial step, last step, life time]
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
            if (Nsteps*TimeStep*(life[-1]-life[0]+1))>minlife :
                
                if name in dicoNames: 
                    dicoNames[name].append([key,life[0],life[-1],(life[-1]-life[0]+1)*TimeStep])
                    dicoStats[name]['lifetime']+=Nsteps*TimeStep*(life[-1]-life[0]+1)
                    total+=Nsteps*TimeStep*(life[-1]-life[0]+1)
                else :
                    dicoNames[name]=[[key,life[0],life[-1],(life[-1]-life[0]+1)*TimeStep]]
                    dicoStats[name]={'lifetime':Nsteps*TimeStep*(life[-1]-life[0]+1), '#atoms': len(atomsCluster)}
                    total+=Nsteps*TimeStep*(life[-1]-life[0]+1)

                if life[0] in dicoTimes :
                    dicoTimes[life[0]].append([key,name,life[-1],(life[-1]-life[0]+1)*Nsteps*TimeStep])
                else:
                    dicoTimes[life[0]]=[[key,name,life[-1],(life[-1]-life[0]+1)*Nsteps*TimeStep]]
    
    if rings==1 and t:
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
                    
                    if rings==1 and t :
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

        
    if(rings==1 and t):
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

if __name__ == "__main__":
    main(sys.argv[1:])
