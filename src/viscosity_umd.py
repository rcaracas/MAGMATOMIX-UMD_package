#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Dec  8 14:39:40 2022
@author: Emma Stoutenburg
"""

import matplotlib.pyplot as plt

import sys, getopt, os.path

import numpy as np

import umd_processes_fast as umdpf

from scipy.integrate import simps

import pandas as pd


kb = 1.3806503e-23  #m2kg/s2K

# manual autocorrelation function written by Zhi Li which does not non-dimensionalize the stresses!:

# This function returns the autocorrelation function of array <tab> over <length> time steps

# with the method of origin shifts and after the step <origin>

def autocorrel(tab,firststep,originstep,length):

    tablength=np.shape(tab)[0]

    # print(f'length of umd is {tablength}')

    numborigin=0

    tabautocorr=np.zeros(length)

    tabautocorrsigma=np.zeros(length)

    for orig in np.arange(firststep,tablength-length,originstep):

        numborigin = numborigin+1

        tabautocorr=tabautocorr+(tab[orig]-np.mean(tab[firststep:-1]))*(tab[orig:orig+length]-np.mean(tab[firststep:-1]))

    tabautocorr=tabautocorr/float(numborigin)

    return tabautocorr

def BuildStressTensor(AllSnapshots,firststep):

    TotTimeStep=len(AllSnapshots)

    StressTensorXY_GPa=np.zeros(TotTimeStep)

    StressTensorXZ_GPa=np.zeros(TotTimeStep)

    StressTensorYZ_GPa=np.zeros(TotTimeStep)

    for ii in range(firststep,TotTimeStep):

        StressTensorXY_GPa[ii] = AllSnapshots[ii].stress[3]

        StressTensorYZ_GPa[ii] = AllSnapshots[ii].stress[4]

        StressTensorXZ_GPa[ii] = AllSnapshots[ii].stress[5]

    StressTensorXY = StressTensorXY_GPa*(1e9)

    StressTensorXZ = StressTensorXZ_GPa*(1e9)

    StressTensorYZ = StressTensorYZ_GPa*(1e9)

    return StressTensorXY,StressTensorXZ,StressTensorYZ

def autocor_int_eta(AllSnapshots, step_initial, T, volume, UMDname, lagn):

    pp = UMDname[:-8] + '.visc.dat'

    StressTensorXY,StressTensorXZ,StressTensorYZ = BuildStressTensor(AllSnapshots, step_initial)

    intvaluesxy = []

    intvaluesxz = []

    intvaluesyz = []

    setlagn = lagn  # should this be varied? is it appropriate to use for "length" in the og fx?

    tau = range(0,setlagn)

    corxy = autocorrel(StressTensorXY, firststep=step_initial, originstep=1, length=setlagn)

    corxz = autocorrel(StressTensorXZ, firststep=step_initial, originstep=1, length=setlagn)

    coryz = autocorrel(StressTensorYZ, firststep=step_initial, originstep=1, length=setlagn)

    for i in range(1,setlagn+1):

        # Use Simpson's rule to integrate the autocorrelation functions

        areaxy = simps(corxy[0:i], dx=1)

        areaxz = simps(corxz[0:i], dx=1)

        areayz = simps(coryz[0:i], dx=1)

        intvaluesxy.append(areaxy)

        intvaluesxz.append(areaxz)

        intvaluesyz.append(areayz)

        #print ('area: ',area)

    coravg = []

    intvalavg = []

    for i in range(0, len(intvaluesxy)):

        avgcor = (corxy[i] + corxz[i] + coryz[i]) / 3

        avgint = (intvaluesxy[i] + intvaluesxz[i] + intvaluesyz[i]) / 3

        coravg.append(avgcor)

        intvalavg.append(avgint)

    # for viscosity calc: V in m3, kb in J/K, T in K, P in Pa and dividing by 1e15 for fs -> s to get Pa*s

    eta_xy = (volume / (kb * T * (1e15))) * np.asarray(intvaluesxy)

    eta_xz = (volume / (kb * T * (1e15))) * np.asarray(intvaluesxz)

    eta_yz = (volume / (kb * T * (1e15))) * np.asarray(intvaluesyz)

    eta_avg = (volume / (kb * T * (1e15))) * np.asarray(intvalavg)

    viscdf = pd.DataFrame({'tau': tau, 'autocorxy_Pa**2':corxy, 'autocorxz_Pa**2':corxz, 'autocoryz_Pa**2': coryz,

                           'intautocorxy_Pa*s':eta_xy, 'intautocorxz_Pa*s':eta_xz,

                           'intautocoryz_Pa*s':eta_yz, 'viscxy-xz-yzavg_Pa*s': eta_avg})

    viscdf.to_csv(pp, index=None, sep=' ')

    return None

def main(argv):

    umdpf.headerumd()

    UMDname = ''

    firststep = 0

    windowmax=2000

    try:

        opts, arg = getopt.getopt(argv, "hf:i:n:")

    except getopt.GetoptError:

        print(

            'viscosity_umd.py -f <umdfile> -i <Initial TimeStep> -n <lag for autocorrelation>')

        sys.exit(2)

    for opt, arg in opts:

        if opt == '-h':

            print(

                'viscosity_umd_firstneg.py to compute the viscosity of the fluid from the self-correlation of stresses\n'

                'assumes timesteps are in fs and stresses are in GPa to returns viscosity in Pa*s')

            print('viscosity_umd.py -f <umdfile> -i <InitialStep> -n <lag for autocorrelation>')

            print('umdfile = input file with the trajectory, in UMD format.')

            print('initialstep = Initial TimeStep from umd file. Default = 0')

            print('lag for autocorrelation = default 2000 fs')

            sys.exit()

        elif opt in ("-f"):

            UMDname = str(arg)

        elif opt in ("-i"):

            firststep = int(arg)

        elif opt in ("-n"):

            windowmax = int(arg)

    if (os.path.isfile(UMDname)):

        print('The first ', firststep, 'timesteps will be discarded')

        AllSnapshots = []

        (AllSnapshots, TimeStep) = umdpf.read_stresses_4visc(UMDname)

        print('Len(allsnapshots),timestep : ', len(AllSnapshots), TimeStep)

        volume = (AllSnapshots[0].cellvolume)/(1e30)  # m3

        print('Volume from the umd file = ', volume*(1e30), 'Å3')

        T = 0.0

        for ii in range(firststep, len(AllSnapshots)):

            T = T + AllSnapshots[ii].temperature

        #            print('current T at step ',ii,' is:', AllSnapshots[ii].T)

        T = T / (len(AllSnapshots) - firststep)

        print('Average T from the umd file = ', T, 'K')

        autocor_int_eta(AllSnapshots, firststep, T, volume, UMDname, windowmax)

    else:

        print('the umdfile ', UMDname, ' does not exist')

        sys.exit()

if __name__ == "__main__":

    args = sys.argv[1:]

    main(args)
