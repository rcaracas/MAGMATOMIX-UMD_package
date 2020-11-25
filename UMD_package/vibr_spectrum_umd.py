#!/usr/bin/env python3

###
##AUTHORS: Zhi Li, Razvan Caracas
###

import sys, getopt, os.path
import numpy as np
from scipy.fftpack import dct, fftfreq
from . import umd_process as umd

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

    # recenter data ( deduct average)
    #TimeMatrix = TimeMatrix - TimeMatrix.mean(axis=0, keepdims=True) 
    #normalization factor
    temp = 1.0/np.arange(nostep,nostep-maxtau,-1)
    normalization = np.diag(temp)
    for ientry in range(noentries): 

        temp1 = np.correlate(TimeMatrix[:,ientry],TimeMatrix[:,ientry],mode='full')[len(TimeMatrix[:,ientry])-1:] #same as that of [len(TimeMatrix[ientry])-1:]
        #although it's fast, it's not normalized
        autocorrelation[:,ientry] = np.matmul(normalization,temp1[0:maxtau]) 
        #fft_correlation[:,ientry] = np.abs(fft(np.correlate(TimeMatrix[:,ientry],TimeMatrix[:,ientry],mode='full')))[:len(TimeMatrix[ientry])] #attention!!! that correlation is a even function, thats why I choose temp[:-1] instead of temp
        fft_correlation[:,ientry] = dct(autocorrelation[:,ientry],1)/2.0 * 2.0 * timestep #this gives exactly the same answer as above and hope you know what is time shifting in discrete fourier transform(that is the reason to have np.abs)

    # we calculate frequency
    freq = fftfreq(2*maxtau,d=timestep)[:maxtau] 
           
    return autocorrelation,fft_correlation,freq                	

def main(argv):
    umd.headerumd()
    umdfile = 'output.umd.dat'
    temperature = 5000
#    start_time  = time.time()
    #---------------------------------------------
    #           Warning
    #maxtau == max correlation time, should be 
    #           less than 1/2*total_simulation_steps        
    #---------------------------------------------
    try:
        opts, arg = getopt.getopt(argv,"hf:t:",["help","fumdfile="])
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
    (MyCrystal,AllSnapshots,TimeStep)=umd.readumd(umdfile,0)
    # make TimeMatrix file
    TimeMatrix = np.empty((len(AllSnapshots),3*MyCrystal.natom))
    for ii in range(len(AllSnapshots)):
        icounter = -1
        for jj in range(MyCrystal.natom):
            for zz in range(3):
                icounter = icounter + 1
                TimeMatrix[ii][icounter]=AllSnapshots[ii].atoms[jj].vels[zz]

    # correlation for different atom, direction
    autocorrelation,fft_correlation,freq = correlation(TimeMatrix,TimeStep,temperature)

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
    temp = total_average_correlation[0]
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
    print('The diffusion coefficients are',diffusion_coefficient, 'in m^2/s')
    #---------------------------------------------
    vaffilename = umdfile[:-8] + '.vels.scf.dat'
    nf = open(vaffilename,'w')
    headerstring = 'Vel_AutoCorr_Fct(fs)\t'
    for itype in range(MyCrystal.ntypat):
        headerstring = headerstring + MyCrystal.elements[itype] + '\t'
    headerstring = headerstring + '\n'
    nf.write(headerstring)
#    nf.write('velocity autocorrelation function(correlation_time(fs) species1 species2 ... total)\n')
    for ii in range(len(average_correlation)):
        string=str(ii*TimeStep)
        for jj in range(MyCrystal.ntypat):
            string=string + '    '+ str(average_correlation[ii][jj])
        string = string + '    ' + str(total_average_correlation[ii]) + '\n'
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

#    end_time = time.time()
#    print('it takes ',(end_time-start_time)/60, 'mins to finish this job!')

if __name__ == "__main__":
   main(sys.argv[1:])
