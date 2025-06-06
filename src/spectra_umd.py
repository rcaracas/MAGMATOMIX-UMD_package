#!/usr/bin/env python3
# -*- coding: utf-8 -*-
### Author : Razvan Caracas, with help from chatGPT
import sys, getopt, os
from os import times

from PIL.ImagePalette import random
from numpy.ma.core import array
from scipy.signal import freqs

import umd_processes_fast as umdpf
import crystallography as cr
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from scipy.ndimage import uniform_filter1d


# Placeholder function for autocorrelation and FFT
def _compute_autocorrelation(arr, max_window_size=None):
    """
    Compute the self-correlation function of a 1D array.

    Parameters:
    arr (np.ndarray): Input 1D array with shape (nsteps,).
    max_window_size (int, optional): Maximum window size for the autocorrelation length.

    Returns:
    np.ndarray: An array containing the self-correlation function for the nsteps dimension.
    """
    nsteps = len(arr)
    correlation = np.zeros(nsteps)
    #print('In compute_autocorrelation with arr of size', nsteps,' and window size', max_window_size)
    # Loop over each time lag
    for lag in range(nsteps):
        if max_window_size is not None and lag >= max_window_size:
            break
        valid_range = nsteps - lag
        correlation[lag] = np.mean(arr[:nsteps - lag] * arr[lag:nsteps])
    if max_window_size is not None:
        correlation = correlation[:max_window_size]
    #print('size of correlation is ', len(correlation))
    return correlation

def compute_vibrational_spectrum(combined_autocorr, dt):
    """
    Compute the vibrational spectrum from pre-computed velocity autocorrelation data.

    Parameters:
    combined_autocorr : ndarray
        Array of shape (N_atomic_types, N_steps) representing the pre-computed velocity autocorrelation functions for atomic types.
    dt : float
        Time step of the simulation in picoseconds (ps).

    Returns:
    result : ndarray
        Array with positive wavenumbers in the first column and the spectrum of each atomic type in subsequent columns.
    """
    # Number of atomic types and time steps
    print('  *** IN compute_vibrational_spectrum ')
    combined_autocorr=np.array(combined_autocorr)
    N_atomic_types, N_steps = combined_autocorr.shape
    print('no atomic tpyes is ', N_atomic_types)
    #print('no of steps is ', N_steps)

    # Define the frequency domain
    frequencies = np.fft.fftfreq(N_steps, d=dt)  # Frequencies in THz if dt is in fs
    #print('some frequencies are ', frequencies[120:129])
    #print('no. of frequencies is ', len(frequencies))

    # Convert frequencies from THz to cm^-1
    c_in_thz_cm = 0.00029979  # Speed of light in THz·cm
    wavenumbers = frequencies / c_in_thz_cm

    # Extract positive frequencies
    positive_indices = frequencies >= 0
    #print('some of the positive indices are ', positive_indices[120:129])
    positive_wavenumbers = wavenumbers[positive_indices]
    #print('some of the positive wavenumbers are ', positive_wavenumbers[120:129])
    #print('numer of positive wavenumbers is ', len(positive_wavenumbers))

    # Perform FFT on each autocorrelation function and extract the real part
    spectra = positive_wavenumbers.reshape(1,-1)
    #print('size of spectra is ', spectra.shape)
    for i in range(N_atomic_types):
        vacf = combined_autocorr[i]
        fft_vacf = np.fft.fft(vacf)
        real_fft_vacf = np.real(fft_vacf)
        positive_real_fft_vacf = real_fft_vacf[positive_indices]
        #spectrum = positive_real_fft_vacf / np.max(positive_real_fft_vacf)
        # Optional: smooth the spectrum to reduce noise
        #smoothed_spectrum = uniform_filter1d(spectrum, size=5)
        spectra = np.vstack([spectra, positive_real_fft_vacf.reshape(1, -1)])

    print('size of result is ', spectra.shape)

    print('  *** OUT OF compute_vibrational_spectrum ')
    return spectra

def plotaurocorrelation(autocorr_byatoms, fft_values):
    autocorr_byatoms = np.array(autocorr_byatoms)
    fft_values = np.array(fft_values)
    print('new len of fft is ', len(fft_values[0]))
    # Plot autocorrelation
    plt.figure(figsize=(8, 6))
    plt.subplot(2, 1, 1)
    for idx in range(autocorr_byatoms.shape[0]):
        plt.plot(autocorr_byatoms[idx], label=f'Atom type {idx + 1}')
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

def prefactor(option, volume, temperature, timestep):
    #print('computing prefactors, for option ', option)
    kb = 1.3806503  # m2kg/s2K
    if option == 0:
        print('Prefactor for viscosity is ', volume * timestep / (kb * temperature) * 1.0E-04)
        return 1.0E-04 * volume * timestep / (kb * temperature)
    elif option in [10, 11, 12, 13]:
        print('computing prefactors, for option ', option)
        Boltzmann = 8.6173303 * 10 ** (-5)  # eV
        Avogadro = 6.022140857 * (10 ** (23))  # mol-1
        au_angstrom_square = 1.0 / Avogadro * (10 ** (-3)) * (10 ** 10) / (1.602176634 * 10 ** -19)  # unit is eV
        print(' au_angstrom_square is ', au_angstrom_square)
        return (2 * au_angstrom_square / (temperature * Boltzmann))

def process_velocity_data(AllSnapshots, MyCrystal, property, window_size, pref):
    print('In process velocity: size of AllSnaphots', len(AllSnapshots), ' by ', len(AllSnapshots[0]), ' where the number of atoms is ',MyCrystal.natom)
    velsx = velsy = velsz = []
    autocorr_byatomsx = autocorr_byatomsy = autocorr_byatomsz = []
    fft_valuesx = fft_valuesy = fft_valuesz = []
    for iatom in range(MyCrystal.natom):
        velsx = []
        velsy = []
        velsz = []
        for isnap in range(len(AllSnapshots)):
            velsx.append(AllSnapshots[isnap][iatom * 3])
            velsy.append(AllSnapshots[isnap][iatom * 3 + 1])
            velsz.append(AllSnapshots[isnap][iatom * 3 + 2])
        velsx = np.array(velsx)
        velsy = np.array(velsy)
        velsz = np.array(velsz)

#        print('in process_velocity_data, size of velsx ', len(velsx))
        autocorr_byatomsx.append( _compute_autocorrelation(velsx, window_size) * MyCrystal.masses[MyCrystal.typat[iatom]] )
        autocorr_byatomsy.append( _compute_autocorrelation(velsy, window_size) * MyCrystal.masses[MyCrystal.typat[iatom]] )
        autocorr_byatomsz.append( _compute_autocorrelation(velsz, window_size) * MyCrystal.masses[MyCrystal.typat[iatom]] )
#        print('in process_velocity_data, size of autocorr_byatomsx ', len(autocorr_byatomsx),' by ', len(autocorr_byatomsx[0]))
    print('end of process velocity, with sizes:', len(velsx), ' by ', len(velsy), ' by ', len(velsz),' and ',len(autocorr_byatomsx),' by ', len(autocorr_byatomsy),' by ', len(autocorr_byatomsz))
    print('\n')
    return autocorr_byatomsx, autocorr_byatomsy, autocorr_byatomsz

def process_atomic_type_data(autocorr_byatomsx, autocorr_byatomsy, autocorr_byatomsz, MyCrystal):
    autocorr_bytypex = np.zeros((len(autocorr_byatomsx[0]),MyCrystal.ntypat))
    autocorr_bytypey = np.zeros((len(autocorr_byatomsx[0]),MyCrystal.ntypat))
    autocorr_bytypez = np.zeros((len(autocorr_byatomsx[0]),MyCrystal.ntypat))
    print('\nin process atomic_type_data, with natoms', MyCrystal.natom)
    print('in process atomic_type_data, with ntypat', MyCrystal.ntypat)
    print('in process atomic_type_data, with typat', MyCrystal.typat)

    for itype in range(MyCrystal.ntypat):
        type_indices = [i for i in range(MyCrystal.natom) if MyCrystal.typat[i] == itype]
        autocorr_bytypex[:][itype] = np.mean([autocorr_byatomsx[i] for i in type_indices], axis=0)
        autocorr_bytypey[:][itype] = np.mean([autocorr_byatomsy[i] for i in type_indices], axis=0)
        autocorr_bytypez[:][itype] = np.mean([autocorr_byatomsz[i] for i in type_indices], axis=0)
        print("finished computing autocorrelation for atomic type ", itype)
    print("finished with process_atomic_type_data")
    return autocorr_bytypex, autocorr_bytypey, autocorr_bytypez

def write_to_files(autocorr_data, fft_data, TimeStep, filename):
    nsteps = len(autocorr_data[0])
    time_values = np.arange(nsteps) * TimeStep
    print(' dimensions for prinitng of autocorr_data are ', len(autocorr_data), ' by ', len(autocorr_data[0]))
    print(' dimensions for prinitng of fft_data are ', len(fft_data), ' by ', len(fft_data[0]))
    with open(filename, 'w') as f:
        f.write("Time\tAutocorr\tFFT\n")
        for i in range(len(fft_data[0])):
            autocorr_values = "\t".join([str(autocorr_data[j][i]) for j in range(len(autocorr_data))])
            fft_values = "\t".join([str(np.abs(fft_data[j][i])) for j in range(len(fft_data))])
            f.write(f"{time_values[i]}\t{autocorr_values}\t{fft_values}\n")

def main():
    """
    Main function to execute the UMD spectrum analysis tool.

    :return: None
    """
    UMDfile = 'OUTCAR.umd.dat'
    firststep = 0
    window_size = 10000
    property = 0  # 0=viscosity; 1x=vibrations from velocities, 2x=vibrations from positions
    flaggraph = 0
    mean_volume = 1.0
    mean_temperature = 300.0
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
            print('10 : combined atomic velocities analysis')
            print('11 : separated by atomic types and directions (velocities)')
            print('12 : separated by atoms, mixed directions (velocities)')
            print('13 : separated by both atoms and directions (velocities)')
            print('-g options:')
            print('0 : no graphical output (default) ')
            print('1 : graphical output ')
            sys.exit()
        elif opt in ("-f"):
            UMDfile = str(arg)
        elif opt in ("-i"):
            firststep = int(arg)
        elif opt in ("-t"):
            mean_temperature = float(arg)
        elif opt in ("-w"):
            window_size = int(arg)
        elif opt in ("-r"):
            property = int(arg)
        elif opt in ("-g"):
            print("graphical output",arg)
            flaggraph = int(arg)
            if flaggraph != 1:
                flaggraph = 0
    if os.path.isfile(UMDfile):  # Changed from UMDname to UMDfile
        if firststep > 0:
            print('The first ', firststep, 'timesteps will be discarded')
        MyCrystal = cr.Lattice()
        AllSnapshots = []
        pref = 1.0
        if property == 0:  # viscosity part
            stresses = []
            (MyCrystal, AllSnapshots, TimeStep) = umdpf.read_stresses_4visc(UMDfile)
            print('after reading the stresses')
            print('Total length is ', len(AllSnapshots))
            temperatures = [ snapshot.temperature for snapshot in AllSnapshots]
            mean_temperature = np.mean(temperatures)
            print('Mean Temperature: ', mean_temperature, ' K')
            volumes = [snapshot.cellvolume for snapshot in AllSnapshots]
            mean_volume = np.mean(volumes)
            print('Cell volume is ', mean_volume, ' in angstroms^3')

            # Calculate the prefactor based on mean volume and mean temperature
            pref = prefactor(property, mean_volume, mean_temperature, TimeStep)

            # Extract stress components from snapshots
            stresses_xy = np.array([snapshot.stress[3] for snapshot in AllSnapshots])
            stresses_yz = np.array([snapshot.stress[4] for snapshot in AllSnapshots])
            stresses_zx = np.array([snapshot.stress[5] for snapshot in AllSnapshots])

            #print("shears as ",stresses_xy[1:10])

            # Compute autocorrelations
            window_size = min(window_size,int(len(stresses_xy)/2))
            autocorr_xy = _compute_autocorrelation(stresses_xy, window_size)
            autocorr_yz = _compute_autocorrelation(stresses_yz, window_size)
            autocorr_zx = _compute_autocorrelation(stresses_zx, window_size)

            # Integrate the autocorrelations and multiply by prefactor
            intvaluesxy = np.array([simps(autocorr_xy[:i], dx=1) for i in range(1, len(autocorr_xy) + 1)]) * pref
            intvaluesyz = np.array([simps(autocorr_yz[:i], dx=1) for i in range(1, len(autocorr_yz) + 1)]) * pref
            intvalueszx = np.array([simps(autocorr_zx[:i], dx=1) for i in range(1, len(autocorr_zx) + 1)]) * pref

            # Calculate the mean of the integral values for each step
            visc_mean = np.array([(intvaluesxy[i] + intvaluesyz[i] + intvalueszx[i]) / 3 for i in range(int(len(intvaluesxy)))])
            viscosity_integrals = [intvaluesxy, intvaluesyz, intvalueszx, visc_mean]
            # Store autocorrelations for plotting and writing
            autocorr_stresses = [autocorr_xy, autocorr_yz, autocorr_zx]

            # Plot autocorrelation and viscosity integrals
            plotaurocorrelation(autocorr_stresses, viscosity_integrals)

            # Write results to file
            viscfile = UMDfile[:-7] + 'visc.dat'
            write_to_files(autocorr_stresses, viscosity_integrals, TimeStep, viscfile)
        elif property in [10, 11]:  # vibrational spectrum from atomic velocities
            fftforprinting = []
            (MyCrystal, AllSnapshots, TimeStep, length) = umdpf.read_values(UMDfile, "velocity", "line", 1, firststep, laststep=None, cutoff="all", nCores=None)
            print("Spectra from velocities: size of AllSnapshots", len(AllSnapshots), " by ", len(AllSnapshots[0]), " where the number of atoms is ", MyCrystal.natom)

            autocorr_byatomsx, autocorr_byatomsy, autocorr_byatomsz = process_velocity_data(AllSnapshots, MyCrystal, property, window_size, pref)
            autocorr_byatomsx = np.array(autocorr_byatomsx)
            autocorr_byatomsy = np.array(autocorr_byatomsy)
            autocorr_byatomsz = np.array(autocorr_byatomsz)
            autocorr_bytypex, autocorr_bytypey, autocorr_bytypez = process_atomic_type_data(autocorr_byatomsx,autocorr_byatomsy,autocorr_byatomsz,MyCrystal)
            print('size of autocorr_bytypex is ', len(autocorr_bytypex), ' by ', len(autocorr_bytypex[0]))
            print('size of autocorr_bytypey is ', len(autocorr_bytypez), ' by ', len(autocorr_bytypey[0]))
            print('size of autocorr_bytypez is ', len(autocorr_bytypex), ' by ', len(autocorr_bytypez[0]))

            pref = prefactor(property, 1.0, mean_temperature, TimeStep)
            print(' prefactor for velocities is ', pref)

            if property == 10:  # separated by atomic types
                combined_autocorr = [np.sum([autocorr_bytypex[i], autocorr_bytypey[i], autocorr_bytypez[i]], axis=0) for i in range(len(autocorr_bytypex))]
                #print('size of combined_autocorr is ', len(combined_autocorr), ' by ', len(combined_autocorr[0]))
                spectrum = compute_vibrational_spectrum(combined_autocorr,TimeStep)
                #print("sizes for printing are ",correlforprinting.shape, ' and ', fftforprinting.shape)
                plotaurocorrelation(combined_autocorr, spectrum)
                #plotaurocorrelation(correlforprinting,fftforprinting)
                spectrumfile = UMDfile[:-7] + 'vibr_spec.dat'
                write_to_files(combined_autocorr, spectrum, TimeStep, spectrumfile)
                #write_to_files(correlforprinting, spectrum, TimeStep, spectr umfile)
            elif property == 11:  # separated by atomic types and directions
                autocorr_bytypex, autocorr_bytypey, autocorr_bytypez = process_atomic_type_data(autocorr_byatomsx, autocorr_byatomsy, autocorr_byatomsz, MyCrystal)
                autocorr_bytypex = np.array(autocorr_bytypex)
                autocorr_bytypey = np.array(autocorr_bytypey)
                autocorr_bytypez = np.array(autocorr_bytypez)
                combined_autocorr = [[1 / 3.0 * np.sum(autocorr_bytypex[i]), 1 / 3.0 * np.sum(autocorr_bytypey[i]), 1 / 3.0 * np.sum(autocorr_bytypez[i])] for i in range(len(autocorr_bytypex))]
                fft_combined = [perform_fft(ac,pref) for ac in combined_autocorr]
                plotaurocorrelation(combined_autocorr, fft_combined)

        elif property == 2:  # vibrational spectrum from atomic positions
            print('Computing the vibrational spectra form atomic positions is in construction.')
        else:
            print("Only options 0, 10, 11, 12, and 13 are valid for -r ")
            sys.exit()
        if flaggraph == 1:
            plotaurocorrelation(autocorr_byatoms, fft_values)
    else:
        print('the umdfile ', UMDfile, ' does not exist')
        sys.exit()
if __name__ == "__main__":
    main()
