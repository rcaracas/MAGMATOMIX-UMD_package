#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  2 14:54:47 2023

@author: Kevin Jiguet-Covex
"""

import sys,getopt,os.path
import crystallography as cr
import umd_process as umdp
import time
import numpy as np
import concurrent.futures
from functools import partial
import ctypes
from os.path import join
from scipy.fftpack import dct, fftfreq
import math



def correlation(TimeMatrix,timestep):
  
    # TimeMatrix should be in matrix format
    # entry1 entry2 entry3 ..
    # 0 1 2 3 4 5 6 
    # 1 1 2 3 4 5 6
    #.....
#    print("TimeMatrix=",TimeMatrix)
    
    nostep = len(TimeMatrix)
    noentries = len(TimeMatrix[1])
    maxtau = int(nostep / 2)
    autocorrelation = np.empty((maxtau,noentries))
    fft_correlation = np.empty((maxtau,noentries)) 
    temp = 1.0/np.arange(nostep,nostep-maxtau,-1)
    normalization = np.diag(temp)        
    
    for ientry in range(noentries): 

        temp1 = np.correlate(TimeMatrix[:,ientry],TimeMatrix[:,ientry],mode='full')[len(TimeMatrix[:,ientry])-1:] #same as that of [len(TimeMatrix[ientry])-1:]
        #although it's fast, it's not normalized
        autocorrelation[:,ientry] = np.matmul(normalization,temp1[0:maxtau]) 
        fft_correlation[:,ientry] = dct(autocorrelation[:,ientry],1)/2.0 * 2.0 * timestep #this gives exactly the same answer as above and hope you know what is time shifting in discrete fourier transform(that is the reason to have np.abs)

#    print("autoc=",autocorrelation[0])
    # we calculate frequency
    freq = fftfreq(2*maxtau,d=timestep)[:maxtau] 
#    print(fft_correlation)       
    return autocorrelation,fft_correlation,freq                	

def main(argv):
    start=time.time()
    maxsize=0
    minsize=0
    umdfile=''
    popfile=''
    minlife=0
    centralatom=''
    outeratom=''
    temperature=5000
    try:
        opts, args = getopt.getopt(argv,"hf:p:c:o:t:s:S:T:",["fumdfile","pPopulation","cCentralatom","oOuteratom","tMintime","sMinsize","SMaxsize","TTemperature"])
    except getopt.GetoptError:
        print ('msd_umd.py -f <umdfilename> -p <populfilename> -c <centralatom> -o <outeratom> -t <mintime> -s <minsize> -S <maxsize> -T <temperature>')
        sys.exit(2)
    print("opts=",opts)

    for o in opts:
        opt=o[0]
        arg=o[1]
        if opt in ('-h', "--help"):
            print ('vibr_clusters_umd.py program to compute the atomic velocity self-correlation of atomic clusters')
            print ('and extract relevant properties: vibrational spectrum and diffusion coefficient')
            print ('vibr_clusters_umd.py -f <umdfilename> -p <populfilename> -c <centralatom> -o <outeratom> -t <mintime> -s <minsize> -S <maxsize> -T <temperature>')
            print ('the program needs an UMD file and a population file')
            sys.exit()
        elif opt in ("-f", "--fumdfile"):
            umdfile = str(arg)
            print ('I will use the ',umdfile,' for input','\n')
        elif opt in ("-p", "--pPopulation"):
            popfile = str(arg)
            print ('I will use the ',popfile,' to determine the population','\n')
        elif opt in ("-c", "--cCentralatom"):
            centralatom = str(arg)
        elif opt in ("-o", "--oOuteratom"):
            outeratom = str(arg)
        elif opt in ("-t", "--tMintime"):
            minlife=float(arg)
        elif opt in ("-S", "--SMaxsize"):
            maxsize=float(arg)
        elif opt in ("-s", "--sMinsize"):
            minsize=float(arg)
        elif opt in ("-t","--tTemperature"):
            temperature = float(arg)

    

            
            
    Boltzmann=8.6173303*10**(-5) #eV
    Avogadro=6.022140857*(10**(23)) # mol-1
    au_angstrom_squre = 1.0 / Avogadro * (10**(-3)) * (10**10) / (1.602176634 * 10**-19) # unit is eV
    kB_T = temperature * Boltzmann

    MyCrystal = cr.Lattice()
        
        
    Centralatoms = []
    Outeratoms = []
    AllAtoms=[]
    ff=open(umdfile,'r')
    while True:
        line = ff.readline()
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
        if entry[0] == "rprimd_c":
            break
    
    
    for iatom in range(MyCrystal.natom):
        if MyCrystal.elements[MyCrystal.typat[iatom]]==centralatom:
            Centralatoms.append(iatom)    
            AllAtoms.append(iatom)#contains the list with the index of the central atoms from the 0 ... natom
        if MyCrystal.elements[MyCrystal.typat[iatom]]==outeratom:
            Outeratoms.append(iatom) 
            AllAtoms.append(iatom)#contains the list with the index of the coordinating atoms from the 0 ... natom
        
    print("searching for clusters containing only the following atoms : ", AllAtoms)
    
    dicoVel={x:[] for x in (AllAtoms)}
    dicoPos={x:[] for x in (AllAtoms)}

    timelec=time.time()

    NumSteps=0
    while True : 
        line = ff.readline()
        if not line: break
        #print(line,len(line))
        if len(line) > 1:
            line=line.strip()
            entry=line.split()
            if entry[0]=='atoms:':
                NumSteps+=1
                for iatom in range(MyCrystal.natom):
                    line = ff.readline()
                    if iatom in dicoVel.keys():
                        entry = line.strip().split()
                        dicoVel[iatom].append([float(entry[9]),float(entry[10]),float(entry[11])])#Velocities
                        dicoPos[iatom].append([float(entry[6]),float(entry[7]),float(entry[8])])#Absolute coordinate
    ff.close()
    timelecb=time.time()
    dicoClusters={}
    Names=[]
    dicoClusterTypes={t:{a:'' for a in AllAtoms} for t in range(NumSteps)}
        
    fp = open(popfile,'r')   
    fp.readline()
    fp.readline()
    while True : 
        line=fp.readline()
        if not line : break
        if len(line)>1 : 
            entry=line.strip().split()
            if float(entry[3])>=minlife and minsize<=len(entry[4:])<=maxsize:
                Cluster="["+(line.split('[')[-1]).strip("\n")
                if Cluster in dicoClusters.keys():
                    dicoClusters[Cluster]['Times'].append([int(entry[1]),int(entry[2])])
                    for t in range(int(entry[1]),int(entry[2])+1):
                        for at in eval(Cluster):
                            dicoClusterTypes[t][at] = dicoClusters[Cluster]['Molecule']
                else :
                    clusterType = [0 for _ in range(MyCrystal.ntypat)]
                    speciesName = ''
                    for atom in eval(Cluster):
                        clusterType[MyCrystal.typat[atom]]+=1                    
                    for indice in range(len(clusterType)) :
                        if clusterType[indice] > 0 :
                            speciesName += MyCrystal.elements[indice]+'_'+str(clusterType[indice])+" "
                    if speciesName not in Names : 
                        Names.append(speciesName)
                    dicoClusters[Cluster]={'Times':[[int(entry[1]),int(entry[2])]],'Molecule':speciesName}

                    for t in range(int(entry[1]),int(entry[2])+1):
                        for at in eval(Cluster):
                            dicoClusterTypes[t][at] = dicoClusters[Cluster]['Molecule']

    #sys.exit()

    timelecf=time.time()
    dicoTimeMatrixes = {cluster:[] for cluster in Names}
    
    if(MyCrystal.elements.index(centralatom)>MyCrystal.elements.index(outeratom)):
        centralIndex=-1
    else :
        centralIndex=0
        
    timemat=time.time()
    for cluster in dicoClusters.keys():

        Atoms = eval(cluster)
        central = Atoms.pop(centralIndex)
        Atoms = [central]+Atoms #make sure the first atom is the central atom
        
        for life in dicoClusters[cluster]['Times']:
#            print(life,len(Atoms))
            timematrix=np.empty((life[1]-life[0]+1,3*len(Atoms)))
            timematrixRad=np.empty((life[1]-life[0]+1,len(Atoms)-1))
            timematrixAng=np.empty((life[1]-life[0]+1,2*len(Atoms)-2))
 
            for t in range(life[1]-life[0]+1):
                for iatom in range(len(Atoms)) :
                    for ii in range(3):
                        timematrix[t][3*iatom+ii] = float(dicoVel[Atoms[iatom]][t][ii])
                for iatom in range(0,len(Atoms)-1) :#Calculating the radial and angular components of the velocities of the coordinating atoms
                    for ii in range(3):
                        relativeVel = np.array([dicoVel[Atoms[iatom+1]][t][0]-dicoVel[Atoms[0]][t][0],dicoVel[Atoms[iatom+1]][t][1]-dicoVel[Atoms[0]][t][1],dicoVel[Atoms[iatom+1]][t][2]-dicoVel[Atoms[0]][t][2]])
                        relativePos = np.array([dicoPos[Atoms[iatom+1]][t][0]-dicoPos[Atoms[0]][t][0],dicoPos[Atoms[iatom+1]][t][1]-dicoPos[Atoms[0]][t][1],dicoPos[Atoms[iatom+1]][t][2]-dicoPos[Atoms[0]][t][2]])
                        relativeDir = np.array(relativePos/np.linalg.norm(relativePos))
                        
                        
                        
                    
                        radVel = np.matmul(relativeVel,relativeDir)#norm of the radial component
                        angVel = relativeVel - radVel                    #norm of the angular component
                        
                        angVelTheta = np.matmul(np.cross(np.array([1,0,0]),relativeDir),angVel)
                        angVelPhi = np.matmul(angVel - angVelTheta,np.cross(relativeDir,np.cross(np.array([1,0,0]),relativeDir)))
                        
                        timematrixAng[t][2*iatom] =  angVelTheta
                        timematrixAng[t][2*iatom+1] = angVelPhi
                        timematrixRad[t][iatom] =  radVel
                        
            dicoTimeMatrixes[dicoClusters[cluster]['Molecule']].append([timematrix,timematrixAng,timematrixRad])
#    print("mats=",dicoTimeMatrixes)
    timematf=time.time()
#    sys.exit()    
    dicoFreqMatrixes = {cluster:[] for cluster in Names}
    dicoFreqAngMatrixes = {cluster:[] for cluster in Names}
    dicoFreqRadMatrixes = {cluster:[] for cluster in Names}
    taut=time.time()
    for name in dicoTimeMatrixes.keys():
        for TimeMatrix in dicoTimeMatrixes[name]:
            # correlation for different atom in each cluster, direction
            autocorrelation,fft_correlation,freq = correlation(TimeMatrix[0],TimeStep)
            dicoFreqMatrixes[name].append([fft_correlation,freq,autocorrelation])
            autocorrelation,fft_correlation,freq = correlation(TimeMatrix[1],TimeStep)
            dicoFreqAngMatrixes[name].append([fft_correlation,freq,autocorrelation])
            autocorrelation,fft_correlation,freq = correlation(TimeMatrix[2],TimeStep)
            dicoFreqRadMatrixes[name].append([fft_correlation,freq,autocorrelation])
    tautf=time.time()
#    print("freqMat=",dicoFreqMatrixes)
    
    dicoFreqClust = {name:None for name in Names}#Will contain the vibrational spectrum
    dicoFreqAngClust = {name:None for name in Names}
    dicoFreqRadClust = {name:None for name in Names}
#    print("keys=",dicoFreqMatrixes.keys())
    timemerge=time.time()
    for clusterType in dicoFreqMatrixes.keys():
#        print("cluster",clusterType)
        ListFreqMat = [x[0] for x in dicoFreqMatrixes[clusterType]]
        ListFreqAngMat = [x[0] for x in dicoFreqAngMatrixes[clusterType]]
        ListFreqRadMat = [x[0] for x in dicoFreqRadMatrixes[clusterType]]
        ListCorrMat = [x[2] for x in dicoFreqMatrixes[clusterType]]


        FreqMat = ListFreqMat[0].copy()
        FreqAngMat = ListFreqAngMat[0].copy()
        FreqRadMat = ListFreqRadMat[0].copy()
        CorrMat =dicoFreqMatrixes[clusterType][0][2].copy()        
        frequencies = dicoFreqMatrixes[clusterType][0][1].copy()
        index=0

        for mat in ListFreqMat[1:]:
            index+=1
            if len(mat)>len(FreqMat) :
                FreqMat = mat.copy()
                FreqAngMat = ListFreqAngMat[index].copy()
                FreqRadMat = ListFreqRadMat[index].copy()
                frequencies = dicoFreqMatrixes[clusterType][index][1].copy()
                CorrMat = dicoFreqMatrixes[clusterType][index][2]
                
        Coeffs=[1 for _ in range (len(FreqMat))]#For further normalization
        CoeffsCorr=[1 for _ in range (len(FreqMat))]
        freqstep=float(1/len(FreqMat))    
 #       print("fmat=",FreqMat,cluster)
        for mat in range(len(ListFreqMat)):
            N=len(ListFreqMat[mat])
            for k in range(N):
                f=float(k)/N
                freqbin = int(f/freqstep + 0.5) 
     #           print(mat[k])
                FreqMat[freqbin]+=ListFreqMat[mat][k]
                FreqAngMat[freqbin]+=ListFreqAngMat[mat][k]
                FreqRadMat[freqbin]+=ListFreqRadMat[mat][k]
                CorrMat[k]+=ListFreqMat[mat][k]
                
                Coeffs[freqbin]+=1
                CoeffsCorr[k]+=1
            
        for k in range(len(FreqMat)):#Normalization
            FreqMat[k]/=Coeffs[k]
            FreqAngMat[k]/=Coeffs[k]
            FreqRadMat[k]/=Coeffs[k]
            CorrMat[k]/=CoeffsCorr[k]
        
        print("CoeffsCorr=",CoeffsCorr)
        
        (T,n)=np.shape(FreqMat)
        average_CorrMat = np.zeros((T,int(n/3)))
        average_FreqMat = np.zeros((T,int(n/3))) #Averaged over x,y,z               
        average_FreqAngMat = np.zeros((T,int(n/3)-1))
        FreqRadMat*=MyCrystal.masses[MyCrystal.elements.index(outeratom)]* 2 * au_angstrom_squre / kB_T
        for ii in range(T):
            average_FreqMat[ii][0]=(FreqMat[ii][0]+FreqMat[ii][1]+FreqMat[ii][2])*MyCrystal.masses[MyCrystal.elements.index(centralatom)]* 2 * au_angstrom_squre / kB_T  #normalize to 3N -3
            average_CorrMat[ii][0]=(CorrMat[ii][0]+CorrMat[ii][1]+CorrMat[ii][2])*MyCrystal.masses[MyCrystal.elements.index(centralatom)]
            
            for at in range(1,int(n/3)):
                average_FreqMat[ii][at]=(FreqMat[ii][3*at]+FreqMat[ii][3*at+1]+FreqMat[ii][3*at+2])*MyCrystal.masses[MyCrystal.elements.index(outeratom)]* 2 * au_angstrom_squre / kB_T  #normalize to 3N -3
                average_FreqAngMat[ii][at-1] = (FreqAngMat[ii][2*(at-1)] + FreqAngMat[ii][2*(at-1)+1])*MyCrystal.masses[MyCrystal.elements.index(outeratom)]* 2 * au_angstrom_squre / kB_T
                average_CorrMat[ii][at]=(CorrMat[ii][3*at]+CorrMat[ii][3*at+1]+CorrMat[ii][3*at+2])*MyCrystal.masses[MyCrystal.elements.index(outeratom)]                
                                
        
        diffusion_coefficient=np.zeros(int(n/3))
        diffusion_coefficient[0]=average_FreqMat[0][0] / 12.0 * kB_T *(1.602176634 * 10**-19) /(1.0 / Avogadro * (10**(-3)) * MyCrystal.masses[MyCrystal.elements.index(centralatom)]) * 10**(-15)
        for at in range(1,int(n/3)):
            diffusion_coefficient[at]=average_FreqMat[0][at] / 12.0 * kB_T *(1.602176634 * 10**-19) /(1.0 / Avogadro * (10**(-3)) * MyCrystal.masses[MyCrystal.elements.index(outeratom)]) * 10**(-15)
        
        temp=average_CorrMat[0]
        for ii in range(len(average_CorrMat[0])):
            average_CorrMat[:,ii]/=temp[ii]
        
#        print("acm=",average_CorrMat[:,0])
#        print("tcm=",(CorrMat[:,0]))#+CorrMat[:,1]+CorrMat[:,2])*MyCrystal.masses[MyCrystal.elements.index(outeratom)])
        dicoFreqClust[clusterType]=[diffusion_coefficient,average_FreqMat,frequencies,average_FreqAngMat,FreqRadMat,average_CorrMat]
        print(Coeffs,clusterType)
    timemergef=time.time()

#    print("dicofreq=",dicoFreqClust)    
    
    
    
    specfilename = umdfile[:-8] + '.vibr_clust.dat'
    nf = open(specfilename,'w')
    vaffilename = umdfile[:-8] + '.vels_clust.scf.dat'
    nv = open(vaffilename,'w')
    
    header = "vibrational spectrum of "+centralatom+" (central atom) and "+outeratom+" (coordinating atoms)"
    headerVel = "average correlation of atoms per species : "+centralatom+" (central atom) and "+ outeratom+" (coordinating atoms)"
    
    nf.write(header+'\n')
    nv.write(headerVel+'\n')
    
    for species in dicoFreqClust.keys() :
        l=species.split()
        katoms = 0         
        headerstring = '\n'+l[0]+l[1] + ' : Frequency(cm^-1)\t'
        headerstringVel = '\n'+l[0]+l[1] + ': Vel_AutoCorr_fct(fs)\t'
        for el in l :
            print("el=",el)
            print(el.split("_"))
            element=el.split("_")[0]
            for ii in range(int(el.split("_")[1])):    
                katoms+=1
                headerstring = headerstring + 'VDos_' + element +"-"+ str(ii) + '\t'
                headerstringVel = headerstringVel + element +"-"+ str(ii) + '\t'
            if element==outeratom :
                headerstring = headerstring + 'Total_DOS_'+element+'\t'
#                headerstringVel = headerstringVel + 'Total_Vel_'+element+'\t'
        for ii in range(int(l[1].split("_")[1])):    
            headerstring = headerstring + 'VDosAng_' + l[1].split("_")[0] +"-"+ str(ii) + '\t'
        
        for ii in range(int(l[1].split("_")[1])):    
            headerstring = headerstring + 'VDosRad_' + l[1].split("_")[0] +"-"+ str(ii) + '\t'
        
        
        nf.write(headerstring+'\n'+'\n')
        nv.write(headerstringVel+'\n'+'\n')
        
        average_FreqMat = dicoFreqClust[species][1]        
        Ang_fft_correlation = dicoFreqClust[species][3]        
        Rad_fft_correlation = dicoFreqClust[species][4]        
        average_correlation = dicoFreqClust[species][5]

        freq = dicoFreqClust[species][2]
        print("freq=",freq)
        for ii in range(len(average_FreqMat)):
            string=str(freq[ii]*33356.41)
            stringVel = str(ii*TimeStep)
            for jj in range(katoms):
                string=string+'\t'+str(average_FreqMat[ii][jj])
                stringVel = stringVel +'\t' + str(average_correlation[ii][jj])
            totalDosCoord=0
            for jj in range(1,katoms):
                totalDosCoord+=average_FreqMat[ii][jj]
        
            string+='\t'+str(totalDosCoord)
            for jj in range(katoms-1):
                string=string+'\t'+str(Ang_fft_correlation[ii][jj])
            for jj in range(katoms-1):
                string=string+'\t'+str(Rad_fft_correlation[ii][jj])
            
            #            string = string + '\t' + str(total_fft_correlation[ii]) + '\n'
            nf.write(string+'\n')
            nv.write(stringVel+'\n')
    nf.close()
    nv.close()
#    print("total time : ",time.time()-start," time calc :",tautf-taut,"time merge :", timemergef-timemerge," timemat :",timematf-timemat," timelec :", timelecf-timelec, "timelecf =",timelecb-timelec)
     
        


if __name__ == "__main__":
    main(sys.argv[1:])
