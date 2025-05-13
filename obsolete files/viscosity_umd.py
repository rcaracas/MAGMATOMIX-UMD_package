#!/usr/bin/env python3
# -*- coding: utf-8 -*-


### Author : Razvan Caracas, with help from chatGPT
import sys, getopt, os
import umd_processes_fast as umdpf
import crystallography as cr
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson



# Placeholder function for autocorrelation and int
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
    # Loop over each time lag
    for lag in range(nsteps):
        if max_window_size is not None and lag >= max_window_size:
            break
        correlation[lag] = np.mean((arr[:nsteps - lag] * arr[lag:nsteps]))
    if max_window_size is not None:
        correlation = correlation[:max_window_size]

    return correlation


def plotaurocorrelation(UMDfile,autocorr_bydir, int_values,zero_autocorr):
    """
    

    Parameters
    ----------
    UMDfile : str
    name of the umd.dat file used
    autocorr_bydir : list
        autocorrelation by direction (XY, YZ, ZX)
    int_values : list
        list of the integral of the autocorrelation
    zero_autocorr : int
        first time an autocorrelation reach zero

    Returns
    -------
    None. Plot the graphic

    """
    autocorr_bydir = np.array(autocorr_bydir)
    int_values = np.array(int_values)
    print('new len of int is ', len(int_values[0]))
    # Plot autocorrelation
    plt.figure(figsize=(8, 6))
    plt.subplot(2, 1, 1)
    plt.plot(autocorr_bydir[0,:], label='Stress XY')
    plt.plot(autocorr_bydir[1,:], label='Stress YZ')
    plt.plot(autocorr_bydir[2,:], label='Stress ZX')
    plt.plot([zero_autocorr],[0], 'or')
    plt.title('Autocorrelation')
    plt.xlabel('Lag')
    plt.ylabel('Autocorrelation Value')
    plt.legend()
    # Plot int of autocorrelation
    plt.subplot(2, 1, 2)
    plt.plot(np.abs(int_values[0,:]), label='Stress XY')
    plt.plot(np.abs(int_values[1,:]), label='Stress YZ')
    plt.plot(np.abs(int_values[2,:]), label='Stress ZX')
    plt.plot(np.abs(int_values[3,:]), 'r',label='Average')
    plt.plot([zero_autocorr],np.abs(int_values[3,zero_autocorr]), 'or')
    plt.title('Viscosity')
    plt.xlabel('$\tau$ (fs)')
    plt.ylabel('Viscosity (Pa.s)')
    plt.legend()
    plt.tight_layout()
    plt.savefig(UMDfile+'visco.pdf')
    plt.show()

def prefactor(volume, temperature, timestep):
    """
    Compute the frefactor of the Green-Kubo equation. 
    Please note that the time step is required to correct integral of the GK equation as the code do the integration over the step and not over time !

    Parameters
    ----------
    volume : float
        volume of the cell in angstrom
    temperature : float
        average temperature of the simulation in Kelvin
    timestep : float
        timestep use in the simulation.

    Returns
    -------
    float
        prefactor needed to compute viscosity, unit are in SI unit to give the viscosity in Pa.s

    """
    #print('computing prefactors, for option ', option)
    kb = 1.3806503  # m2kg/s2K
    print('Prefactor for viscosity is ', volume * timestep / (kb * temperature) * 1.0E-04)
    return 1.0E-04 * volume * timestep / (kb * temperature)


def write_to_files_viscosity(autocorr_data, int_data, TimeStep, filename):
    nsteps = len(autocorr_data[0])
    print(len(autocorr_data[0]))
    print(np.shape(int_data))
    time_values = np.arange(nsteps) * TimeStep
    with open(filename, 'w') as f:
        f.write("Time\tAutocorr\tint\n")
        for i in range(nsteps):
            autocorr_values = "\t".join([str(autocorr_data[j][i]) for j in range(len(autocorr_data))])
            int_values = "\t".join([str(np.abs(int_data[j][i])) for j in range(len(int_data))])
            f.write(f"{time_values[i]}\t{autocorr_values}\t{int_values}\n")
            

def main():
    """
    Main function to execute the UMD viscosity tool.

    :return: None
    """
    UMDfile = 'OUTCAR.umd.dat'
    firststep = 0
    window_size = 1000
    #property = 0  # 0=viscosity; 1x=vibrations from velocities, 2x=vibrations from positions
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
            print('default values: -f OUTCAR.umd.dat -i 0 -w 1000 -t 300 -g 0')
            print('-g options:')
            print('0 : graphical output (default) ')
            print('1 : no graphical output ')
            sys.exit()
        elif opt in ("-f"):
            UMDfile = str(arg)
        elif opt in ("-i"):
            firststep = int(arg)
        elif opt in ("-t"):
            mean_temperature = float(arg)
        elif opt in ("-w"):
            window_size = int(arg)
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
        pref = prefactor(mean_volume, mean_temperature, TimeStep)

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
        i,j,k=0,0,0
        try:
            while autocorr_xy[i]>0 and i<window_size:
                i+=1
            
            while autocorr_xy[j]>0 and j<window_size:
                j+=1
            while autocorr_xy[k]>0 and k<window_size:
                k+=1
            zero_autocorr=max(i,j,k)
            print(zero_autocorr)
            # Integrate the autocorrelations and multiply by prefactor
            intvaluesxy = np.array([simps(autocorr_xy[:i], dx=1) for i in range(1, len(autocorr_xy) + 1)]) * pref
            intvaluesyz = np.array([simps(autocorr_yz[:i], dx=1) for i in range(1, len(autocorr_yz) + 1)]) * pref
            intvalueszx = np.array([simps(autocorr_zx[:i], dx=1) for i in range(1, len(autocorr_zx) + 1)]) * pref

            # Calculate the mean of the integral values for each step
            visc_mean = np.array([(intvaluesxy[i] + intvaluesyz[i] + intvalueszx[i]) / 3 for i in range(int(len(intvaluesxy)))])
            viscosity_integrals = [intvaluesxy, intvaluesyz, intvalueszx, visc_mean]
            
            #Print the viscosity average and std 
            print('Viscosity along XY=',np.mean(intvaluesxy[zero_autocorr::]),' Pa.s')
            print('Viscosity along YZ =',np.mean(intvaluesyz[zero_autocorr::]),' Pa.s')
            print('Viscosity along ZX=',np.mean(intvalueszx[zero_autocorr::]),' Pa.s')
            print('Viscosity =',np.mean(visc_mean[zero_autocorr::]),' +/- ', max([np.mean(visc_mean[zero_autocorr::])-np.mean(intvaluesxy[zero_autocorr::]),np.mean(visc_mean[zero_autocorr::])-np.mean(intvalueszx[zero_autocorr::]),np.mean(visc_mean[zero_autocorr::])-np.mean(intvaluesyz[zero_autocorr::])]),' Pa.s')
            # Store autocorrelations for plotting and writing
            max_diff=max([np.mean(visc_mean[zero_autocorr::])-np.mean(intvaluesxy[zero_autocorr::]),np.mean(visc_mean[zero_autocorr::])-np.mean(intvalueszx[zero_autocorr::]),np.mean(visc_mean[zero_autocorr::])-np.mean(intvaluesyz[zero_autocorr::])])*100/np.mean(visc_mean[zero_autocorr::])
            print('max difference between viscosty along specific direction and average visocosity is :',max_diff,' %')
            if max_diff>20:
                print('WARNING !!! convergence seems to be bad ! You MUST check the graph to (-g 1)')
            autocorr_stresses = [autocorr_xy, autocorr_yz, autocorr_zx]

            # Plot autocorrelation and viscosity integrals
            if flaggraph==0:
                plotaurocorrelation(UMDfile,autocorr_stresses, viscosity_integrals,zero_autocorr)

            # Write results to file
            viscfile = UMDfile[:-7] + 'visc.dat'
            write_to_files_viscosity(autocorr_stresses, viscosity_integrals, TimeStep, viscfile)
                    
        except:
            print("ERORR !!at least one of the auto correlation never converge !")

        
        
    else:
        print('the umdfile ', UMDfile, ' does not exist')
        sys.exit()
        
if __name__ == "__main__":
    main()
