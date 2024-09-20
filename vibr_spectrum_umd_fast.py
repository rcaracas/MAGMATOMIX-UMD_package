#!/usr/bin/env python3
"""
Created on Tue Jun  6 11:00:43 2023

@author: Zhi Li, Razvan Caracas, Kevin Jiguet-Covex
"""

import platform
import sys, getopt, os.path
import numpy as np
from scipy.fftpack import dct, fftfreq
import umd_processes_fast as umdpf
import concurrent.futures
from functools import partial
import time
import ctypes
from os.path import join


current_path=os.path.abspath(__file__)#For all this to work, the file c_autocorrelation_vib.so must be in the same directory than this script
path_split=current_path.split('/')
path_red=path_split[1:-1]
path_new=''
for u in path_red:
    path_new+='/'+u

OS = platform.system()
LibraryName =""

if OS == "Linux":
    LibraryName = 'c_autocorrelation_vib.so'
elif OS == "Windows":
    LibraryName = 'c_autocorrelation_vib.dll'
elif OS == "Darwin":
    LibraryName = 'c_autocorrelation_vib.dylib'


autc_lib = ctypes.cdll.LoadLibrary(join(path_new, LibraryName))
autc_lib.compute_autocorrelation.argtypes = [ctypes.POINTER(ctypes.c_double),ctypes.c_int,ctypes.c_int]
autc_lib.compute_autocorrelation.restype = ctypes.POINTER(ctypes.c_double)


def correlation_par_C(posList,timestep,temperature):

#     pos = posList.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    pos =(ctypes.c_double * len(posList))(*list(posList))
    nostep = len(posList)
    maxtau=int(nostep/2)
    autocorrelationP = autc_lib.compute_autocorrelation(pos,nostep,maxtau)
    autocorrelation = [autocorrelationP[i] for i in range(maxtau)]
    fft_correlation = dct(autocorrelation,1)/2.0 * 2.0 * timestep #this gives exactly the same answer as above and hope you know what is time shifting in discrete fourier transform(that is the reason to have np.abs)

    # we calculate frequency
    return autocorrelation,fft_correlation

def correlation_par(posList,normalization,timestep,temperature):


    nostep = len(posList)
    maxtau=int(nostep/2)
    temp = 1.0/np.arange(nostep,nostep-maxtau,-1)
    normalization = np.diag(temp)
    temp1 = np.correlate(posList,posList,mode='full')[len(posList)-1:]
    autocorrelation = np.matmul(normalization,temp1[0:maxtau])     
    fft_correlation = dct(autocorrelation,1)/2.0 * 2.0 * timestep #this gives exactly the same answer as above and hope you know what is time shifting in discrete fourier transform(that is the reason to have np.abs)

    # we calculate frequency
#    print("autocpar=",autocorrelation)
    return autocorrelation,fft_correlation

def correlation(TimeMatrix,timestep,temperature):
  
    # TimeMatrix should be in matrix format
    # entry1 entry2 entry3 ..
    # 0 1 2 3 4 5 6 
    # 1 1 2 3 4 5 6
    #.....
    
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

    # we calculate frequency
    freq = fftfreq(2*maxtau,d=timestep)[:maxtau] 
    return autocorrelation,fft_correlation,freq                	

def main(argv):
    start_time=time.time()
    umdpf.headerumd()
    umdfile = 'output.umd.dat'
    temperature = 5000
#    start_time  = time.time()
    #---------------------------------------------
    #           Warning
    #maxtau == max correlation time, should be 
    #           less than 1/2*total_simulation_steps        
    #---------------------------------------------
    try:
        opts, arg = getopt.getopt(argv,"hf:t:",["help","fumdfile"])
    except getopt.GetoptError:
        print ('Usage: vel_correlation_umd.py -f <UMDfilename> -t <temperature>')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', "--help"):
            print ('vel_correlation_umd.py program to compute the atomic velocity self-correlation')
            print ('  and extract rlevant properties: vibrational spectrum and diffusion coefficient')
            print ('vel_correlation_umd.py -f <umdfilename> -t <temperature>')
            print (' the program needs an UMD file prepared and the temperature')
            print (' defaults are: output.umd.dat and 5000K')
            sys.exit()
        elif opt in ("-f", "--fumdfile"):
            umdfile = arg
            print ('I will use the ',umdfile,' for input')
        elif opt in ("-t"):
            temperature = float(arg)
            print('the temperature is',temperature)
    if not (os.path.isfile(umdfile)):
        print ('umd file ',umdfile,'does not exist')
        sys.exit()


    Boltzmann=8.6173303*10**(-5) #eV
    Avogadro=6.022140857*(10**(23)) # mol-1
    au_angstrom_squre = 1.0 / Avogadro * (10**(-3)) * (10**10) / (1.602176634 * 10**-19) # unit is eV
    kB_T = temperature * Boltzmann

    # read umdfile
    (MyCrystal,AllSnapsVels,TimeStep,nostep)=umdpf.read_values(umdfile,"velocity")
    # make TimeMatrix file
    AllSnapshots=AllSnapsVels#[:10000]    
#    nostep=10000
    TimeList = [[np.array([AllSnapshots[i][3*at] for i in range(nostep)]) , np.array([AllSnapshots[i][3*at+1] for i in range(nostep)]), np.array([AllSnapshots[i][3*at+2] for i in range(nostep)])] for at in range(MyCrystal.natom)] 
    
    # correlation for different atom, direction


    maxtau=int(nostep/2)
    freq = fftfreq(2*maxtau,d=TimeStep)[:maxtau] 
    print("Computing correlations")
    corred_C = partial(correlation_par_C,timestep=TimeStep,temperature=temperature)
    with concurrent.futures.ProcessPoolExecutor() as executor :
        Listautocorr_C = list(executor.map(corred_C,[TimeList[i//3][i%3] for i in range(3*MyCrystal.natom)]))
    autocorrelation = np.transpose(np.array([x[0] for x in Listautocorr_C]))
    fft_correlation = np.transpose(np.array([x[1] for x in Listautocorr_C]))
    
    print("Averaging")
    # average over species, and mass weighted
    average_correlation = np.zeros((len(autocorrelation),MyCrystal.ntypat))
    total_average_correlation = np.zeros(len(autocorrelation))

    for ii in range(len(autocorrelation)):
        icounter = -1
        for jj in range(MyCrystal.natom):
            for zz in range(3):
                icounter = icounter + 1 
                average_correlation[ii][MyCrystal.typat[jj]] = average_correlation[ii][MyCrystal.typat[jj]] + autocorrelation[ii][icounter] * MyCrystal.masses[MyCrystal.typat[jj]]   
    total_average_correlation =  np.sum(average_correlation,axis=1) 
#    print("ac=",average_correlation)
    temp = total_average_correlation[0]
 #   print("temp=",total_average_correlation)
    total_average_correlation[:] = total_average_correlation[:] / temp
    temp = np.copy(average_correlation[0]) #normalization
    for ii in range(MyCrystal.ntypat):
        average_correlation[:,ii] = average_correlation[:,ii] / temp[ii]

    # average over species, and mass weighted
    diffusion_coefficient = np.zeros(MyCrystal.ntypat)
    average_fft_correlation = np.zeros((len(autocorrelation),MyCrystal.ntypat))
    total_fft_correlation = np.zeros(len(autocorrelation))
    for ii in range(len(fft_correlation)):
        icounter = -1
        for jj in range(MyCrystal.natom):
            for zz in range(3):
                icounter = icounter + 1
                average_fft_correlation[ii][MyCrystal.typat[jj]] = average_fft_correlation[ii][MyCrystal.typat[jj]] + fft_correlation[ii][icounter] * MyCrystal.masses[MyCrystal.typat[jj]]

    for ii in range(MyCrystal.ntypat):
        average_fft_correlation[:,ii] = average_fft_correlation[:,ii] * 2 * au_angstrom_squre / kB_T  #normalize to 3N -3
        diffusion_coefficient[ii] = average_fft_correlation[0][ii] / 12.0 / MyCrystal.types[ii] * kB_T *(1.602176634 * 10**-19) /(1.0 / Avogadro * (10**(-3)) * MyCrystal.masses[ii]) * 10**(-15)
    total_fft_correlation = np.sum(average_fft_correlation,axis=1)
    print('sum of average_fft_correlation:','theory gives a value (3N-3 =',3*MyCrystal.natom-3,') ,numerical calcualtion gives',np.trapz(np.sum(average_fft_correlation,axis=1),freq))
    print ('\n Diffusion:')
    print('For elements:',MyCrystal.elements)
    print('kB_T=',kB_T)
    print('corr=',average_fft_correlation[0])
    print('The diffusion coefficients are',diffusion_coefficient, 'in m^2/s')
    #---------------------------------------------
    vaffilename = umdfile[:-8] + '.vels.scf.dat'
    nf = open(vaffilename,'w')
    headerstring = 'Vel_AutoCorr_Fct(fs)\t'
    for itype in range(MyCrystal.ntypat):
        headerstring = headerstring + MyCrystal.elements[itype] + '\t'
    headerstring = headerstring + 'Average\n'
    nf.write(headerstring)
    for ii in range(len(average_correlation)):
        string=str(ii*TimeStep)
        for jj in range(MyCrystal.ntypat):
            string=string + '\t'+ str(average_correlation[ii][jj])
        string = string + '\t' + str(total_average_correlation[ii]) + '\n'
        nf.write(string)
    nf.close()

    specfilename = umdfile[:-8] + '.vibr.dat'
    nf = open(specfilename,'w')
    headerstring = 'Frequency(cm^-1)\t'
    for itype in range(MyCrystal.ntypat):
        headerstring = headerstring + 'VDos_' + MyCrystal.elements[itype] + '\t'
    headerstring = headerstring + 'Total_DOS\n'
    nf.write(headerstring)
    for ii in range(len(fft_correlation)):
        string=str(freq[ii]*33356.41)
        for jj in range(MyCrystal.ntypat):
            string=string+'\t'+str(average_fft_correlation[ii][jj])
        string = string + '\t' + str(total_fft_correlation[ii]) + '\n'
        nf.write(string)
    nf.close()

    end_time = time.time()
    #print('it takes ',(end_time-start_time)/60, 'mins to finish this job fast!')

if __name__ == "__main__":
   main(sys.argv[1:])
