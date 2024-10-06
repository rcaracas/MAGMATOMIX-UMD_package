#!/usr/bin/env python3
# -*- coding: utf-8 -*-
### Author : Razvan Caracas
import sys, getopt, os, time, subprocess
import umd_processes_fast as umdpf
import crystallography as cr
import numpy as np
import matplotlib.pyplot as plt


def _compute_autocorrelation(array, window_size):
    # Helper function to compute the autocorrelation of a 1D array within the given window size.
    n = len(array)
    # Compute the number of non-NaN entries for each window length
    counts = np.correlate(~np.isnan(array[:n - window_size]), np.ones(window_size), mode='valid')
    array = array - np.nanmean(array)
    autocorr = np.correlate(array, array, mode='full')
    autocorr = autocorr[n - 1:]
    autocorr[:window_size] /= counts
    return autocorr[:window_size]


def perform_fft(autocorr_values):
    # Perform Fast Fourier Transform (FFT) on autocorrelation values.
    fft_values = np.fft.fft(autocorr_values, axis=-1)
    return fft_values


def plotaurocorrelation(autocorr_values,fft_values):
    # Plot autocorrelation
    plt.figure(figsize=(8, 6))
    plt.subplot(2, 1, 1)
    for idx in range(autocorr_values.shape[0]):
        plt.plot(autocorr_values[idx], label=f'Array {idx + 1}')
    plt.title('Autocorrelation')
    plt.xlabel('Lag')
    plt.ylabel('Autocorrelation Value')
    plt.legend()

    # Plot FFT of autocorrelation
    plt.subplot(2, 1, 2)
    for idx in range(fft_values.shape[0]):
        plt.plot(np.abs(fft_values[idx]), label=f'Array {idx + 1}')
    plt.title('FFT of Autocorrelation')
    plt.xlabel('Frequency')
    plt.ylabel('FFT Value (Magnitude)')
    plt.legend()

    plt.tight_layout()
    plt.show()


def prefactor(option,volume,temperature):
    print('computing prefactors')
    kb = 1.3806503e-23  # m2kg/s2K
    if option == 0:
        return  volume/(kb*temperature*1e15)   #(volume / (kb * T * (1e15)))


def main():
    UMDfile = 'OUTCAR.umd.dat'
    firststep = 0
    window_size = 1000
    property = 1  # 1=vibrations from velocities, 2=vibrations from positions
    flaggraph = 0
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hf:i:t:w:r:g")
    except getopt.GetoptError:
        print('spectra_umd.py -f <umdfile> -i <InitialStep> -t <Temperature> -w <WindowSize> -r <Property> -g <Graphical>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('spectra_umd.py program to compute the auto-correlation function of various time-dependent values')
            print('spectra_umd.py -f <umdfile> -i <InitialStep> -t <Temperature> -w <WindowSize> -r <Property> -g <Graphical>')
            print('default values: -f OUTCAR.umd.dat -i 0 -w 1000 -t 300 -r 0 -g 0')
            print('-r options:')
            print('0 : viscosity  (Default)')
            print('1 : vibrational spectrum from atomic velocities')
            print('2 : vibrational spectrum from atomic positions')
            print('-g options:')
            print('0 : no graphical output (default) ')
            print('1 : graphical output ')
            sys.exit()
        elif opt in ("-f"):
            UMDfile = str(arg)
        elif opt in ("-i"):
            firststep = int(arg)
        elif opt in ("-w"):
            window_size = int(arg)
        elif opt in ("-r"):
            property = int(arg)
        elif opt in ("-g"):
            flaggraph = int(arg)
            if flaggraph != 1:
                flaggraph = 0
    if os.path.isfile(UMDfile):  # Changed from UMDname to UMDfile
        print('The first ', firststep, 'timesteps will be discarded')
        MyCrystal = cr.Lattice()
        AllSnapshots = [cr.Lattice]
        (MyCrystal, AllSnapshots, TimeStep, length) = umdpf.read_values(UMDfile, "cellvolume", "line", 1, firststep,laststep=None, cutoff="all", nCores=None)
        volume = np.mean(AllSnapshots[0].cellvolume)/(1e30)
        (MyCrystal, AllSnapshots, TimeStep, length) = umdpf.read_values(UMDfile, "temperature", "line", 1, firststep,laststep=None, cutoff="all", nCores=None)
        temperature = np.mean(AllSnapshots[0].temperature)/(1e30)

        if property == 0:
            stresses = []
            (MyCrystal, AllSnapshots, TimeStep, length) = umdpf.read_values(UMDfile, "velocity", "line", 1, firststep,
                                                                            laststep=None, cutoff="all", nCores=None)
            for isnap in len(AllSnapshots):
                stresses.append(AllSnapshots[isnap].stress[3],AllSnapshots[isnap].stress[4],AllSnapshots[isnap].stress[5])
            autocorr_values = _compute_autocorrelation(stresses, window_size)
            fft_values = perform_fft(autocorr_values)
            prefactor(property,volume,temperature)
        elif property == 1:  # vibrational spectrum from atomic velocities
            velsx = velsy = velsz = []
            (MyCrystal, AllSnapshots, TimeStep, length) = umdpf.read_values(UMDfile, "velocity", "line", 1, firststep,laststep=None, cutoff="all", nCores=None)
            for isnap in len(AllSnapshots):
                velsx.append(AllSnapshots[isnap].atoms[:].vels[0])
                velsy.append(AllSnapshots[isnap].atoms[:].vels[1])
                velsz.append(AllSnapshots[isnap].atoms[:].vels[2])
            autocorr_valuesx = _compute_autocorrelation(velsx, window_size)
            fft_valuesx = perform_fft(autocorr_valuesx)
            autocorr_valuesy = _compute_autocorrelation(velsy, window_size)
            fft_valuesy = perform_fft(autocorr_valuesy)
            autocorr_valuesz = _compute_autocorrelation(velsz, window_size)
            fft_valuesz = perform_fft(autocorr_valuesz)
            prefactor(property,volume,temperature)
        elif property == 2:  # vibrational spectrum from atomic positions
            atomsx = atomsy = atomsz = []
            (MyCrystal, AllSnapshots, TimeStep, length) = umdpf.read_values(UMDfile, "xcart", "line", 1, firststep,
                                                                            laststep=None, cutoff="all", nCores=None)
            for isnap in len(AllSnapshots):
                atomsx.append(AllSnapshots[isnap].atoms[:].xcart[0])
                atomsy.append(AllSnapshots[isnap].atoms[:].xcart[1])
                atomsz.append(AllSnapshots[isnap].atoms[:].xcart[2])
            autocorr_valuesx = _compute_autocorrelation(atomsx, window_size)
            fft_valuesx = perform_fft(autocorr_valuesx)
            autocorr_valuesy = _compute_autocorrelation(atomsy, window_size)
            fft_valuesy = perform_fft(autocorr_valuesy)
            autocorr_valuesz = _compute_autocorrelation(atomsz, window_size)
            fft_valuesz = perform_fft(autocorr_valuesz)
            prefactor(property,volume,temperature)

        else:
            print("Only options 0, 1 and 2 are valid for -r ")
            sys.exit()
    else:
        print('the umdfile ', UMDfile, ' does not exist')
        sys.exit()

    # Perform autocorrelation on the extracted data
    # autocorr_values = _compute_autocorrelation(extracted_data,window_size)
    # Perform FFT on autocorrelation values
    # fft_values = perform_fft(autocorr_values)



if __name__ == "__main__":
    main()
