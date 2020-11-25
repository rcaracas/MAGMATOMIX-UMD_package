#!/usr/bin/env python3

# This script is to compute autocorrelations of the stress tensor
# to extract the viscosity of a LIQUID

import sys, getopt, os.path
import numpy as np
from scipy.optimize import curve_fit
from . import umd_process as umdp

HaeV = 27.21138505
bohrAng = 0.52917721092
HaB3GPa = 2.9421912e4
AtomMass= 1./5.485799090e-4
kBHa    = 1./315775.13
fsau    = 1./0.02418884326505
GPaPa   = 1.e9
aus     = 0.02418884326505e-15
eVTHz   = 1.602176565e-19/6.626070040e-34/1e12
eVHz    = 1.602176565e-19/6.626070040e-34
eV_ovcm = 1.602176565e-19/6.626070040e-34/29979245800.
eVJoule = 1.602176565e-19
hbar    = 6.626070040e-34/(2*np.pi)
epsil0  = 8.854187817e-12
c = 299792458.

# This function returns the autocorrelation function of array <tab> over <length> time steps 
# with the method of origin shifts and after the step <origin>  
def autocorrel(tab,firststep,originstep,length):
    tablength=np.shape(tab)[0]
    numborigin=0
    tabautocorr=np.zeros(length)
    tabautocorrsigma=np.zeros(length)
    for orig in np.arange(firststep,tablength-length,originstep):
        numborigin = numborigin+1
        #tabautocorr=tabautocorr+\
        #             (tab[orig]-np.mean(tab[orig:orig+length]))*\
        #             (tab[orig:orig+length]-np.mean(tab[orig:orig+length]))
        #tabautocorrsigma=tabautocorrsigma+(\
        #             (tab[orig]-np.mean(tab[orig:orig+length]))*\
        #             (tab[orig:orig+length]-np.mean(tab[orig:orig+length])))**2
        tabautocorr=tabautocorr+\
                     (tab[orig]-np.mean(tab[firststep:-1]))*\
                     (tab[orig:orig+length]-np.mean(tab[firststep:-1]))
        tabautocorrsigma=tabautocorrsigma+(\
                     (tab[orig]-np.mean(tab[firststep:-1]))*\
                     (tab[orig:orig+length]-np.mean(tab[firststep:-1])))**2
    tabautocorr=tabautocorr/float(numborigin)
    tabautocorrsigma=tabautocorrsigma/float(numborigin)
    tabautocorrsigma=np.sqrt(np.abs(tabautocorrsigma-tabautocorr**2))
    normalization=np.mean(tab[np.arange(firststep,tablength-length,originstep)]**2)
    return tabautocorr, tabautocorrsigma, numborigin,normalization

# Definition of the fitting functions and different formulas to be used in the viscosity determination
def fitfunc(t,b1,b2,tau1,tau2,w):
    return b1*np.exp(-t/tau1)+b2*np.exp(-t/tau2)*np.sin(w*t)
        
def intfitfunc(t,b1,b2,tau1,tau2,w):
    return -b1*tau1*np.exp(-t/tau1)+b2*np.exp(-t/tau2)/(1./tau2**2+w**2)*(-w*np.cos(w*t)-np.sin(w*t)/tau2)


def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

def viscosity(b1,b2,tau1,tau2,w):
    return b1*tau1+abs(b2)*tau2/(1+w**2*tau2**2)

def errviscosity(b1,b2,tau1,tau2,w,cov_mat):
    vec=[]
    vec=np.append(vec,[tau1])
    vec=np.append(vec,[tau2/(1+w**2*tau2**2)])
    vec=np.append(vec,[b1])
    vec=np.append(vec,[abs(b2)*(1-w**2*tau2**2)/(1+w**2*tau2**2)**2])
    vec=np.append(vec,[-2*w*b2*tau2**3/(1+w**2*tau2**2)**2])
    error=0.0
    for i in np.arange(5):
        for j in np.arange(5):
            error=error+(vec[i]*vec[j]*cov_mat[i,j])
    error = np.sqrt(error)
    return error

#This function builds the complete stress tensor including the kinetic
def BuildStressTensor(AllSnapshots,firststep):
    TotTimeStep=len(AllSnapshots)
    StressTensorXY=np.zeros(TotTimeStep)
    StressTensorXZ=np.zeros(TotTimeStep)
    StressTensorYZ=np.zeros(TotTimeStep)
    StressTensorMean=np.zeros(TotTimeStep)
    for ii in range(firststep,TotTimeStep):
        StressTensorXY[ii] = AllSnapshots[ii].stress[3]
        StressTensorYZ[ii] = AllSnapshots[ii].stress[4]
        StressTensorXZ[ii] = AllSnapshots[ii].stress[5]
        StressTensorMean[ii]=(StressTensorXY[ii]+StressTensorXZ[ii]+StressTensorYZ[ii])/3.
    return StressTensorXY,StressTensorXZ,StressTensorYZ,StressTensorMean

# This function computes the autocorrelation of the stresstensor and extract the viscosity
# by mean of a fit
def Stress2Viscosity(StressTensor,firststep,originshift,length,TimeStep,temperature,volume):
    autocorr, autocorrsigma, numborigin, normalization=\
             autocorrel(StressTensor,firststep,originshift,length)

    # First guess of the parameters
    a0=autocorr[0]
    a0dot=(autocorr[1]-autocorr[0])/TimeStep
    ind=find_nearest(autocorr,a0*np.exp(-1.))
    tau1=ind*TimeStep
    tau2=tau1*2
    b1=tau2/(tau2-tau1)*(a0+a0dot*tau2)
    b2=a0-b1

    print('In stress2Viscosity, values for volume,bohrAng,kBHa,temperature,HaB3GPa,GPaPa,aus = ',volume,bohrAng,kBHa,temperature,HaB3GPa,GPaPa,aus)
    print('where the fit parameters are: a0,a0dot,b1,b2,tau1,tau2,b1/(b2*tau1) = ',a0,a0dot,b1,b2,tau1,tau2,b1/(b2*tau1))
    # Fit and extraction of the values in Pa.s
    fit=curve_fit(fitfunc,np.arange(length)*TimeStep,autocorr,\
                  p0=(b1,b2,tau1,tau2,b1/(b2*tau1)),sigma=autocorrsigma,absolute_sigma=True, maxfev = 50000)
    eta=viscosity(fit[0][0],fit[0][1],fit[0][2],fit[0][3],fit[0][4])*\
                 (volume/bohrAng**3)/(kBHa*temperature)/HaB3GPa*GPaPa*aus
    etaerror=errviscosity(fit[0][0],fit[0][1],fit[0][2],fit[0][3],fit[0][4],fit[1])*\
                 (volume/bohrAng**3)/(kBHa*temperature)/HaB3GPa*GPaPa*aus

    return autocorr, autocorrsigma, fit, eta, etaerror                 
                
def ViscosityAnalysis(AllSnapshots,TimeStep,firststep,originshift,length,temperature,UMDname):
    volume=AllSnapshots[0].cellvolume
    StressTensorXY,StressTensorXZ,StressTensorYZ,StressTensorMean=BuildStressTensor(AllSnapshots,firststep)
    
    autocorrXY, autocorrsigmaXY, fitXY, etaXY, etaerrorXY =\
            Stress2Viscosity(StressTensorXY,firststep,originshift,length,TimeStep,temperature,volume)
    autocorrXZ, autocorrsigmaXZ, fitXZ, etaXZ, etaerrorXZ =\
            Stress2Viscosity(StressTensorXZ,firststep,originshift,length,TimeStep,temperature,volume)
    autocorrYZ, autocorrsigmaYZ, fitYZ, etaYZ, etaerrorYZ =\
            Stress2Viscosity(StressTensorYZ,firststep,originshift,length,TimeStep,temperature,volume)
    autocorr, autocorrsigma, fit, eta, etaerror =\
            Stress2Viscosity(StressTensorMean,firststep,originshift,length,TimeStep,temperature,volume)
    
        
    print('Extracted viscosity from fit:')
    print('Viscosity from tensor XY : ', etaXY, ' +- ', etaerrorXY,' Pa.s')
    print('Viscosity from tensor XZ : ', etaXZ, ' +- ', etaerrorXZ,' Pa.s')
    print('Viscosity from tensor YZ : ', etaYZ, ' +- ', etaerrorYZ,' Pa.s')
    print('Viscosity from average   : ', eta, ' +- ', etaerror,' Pa.s !!!! NOT TO USE')
        
#    matplotlib.rcParams.update({'font.size': 20})
#    mpl.rcParams['xtick.major.size'] = 12
#    mpl.rcParams['xtick.major.width'] = 2
#    mpl.rcParams['xtick.minor.size'] = 6
#    mpl.rcParams['xtick.minor.width'] = 2
#    mpl.rcParams['ytick.major.size'] = 12
#    mpl.rcParams['ytick.major.width'] = 2
#    mpl.rcParams['ytick.minor.size'] = 6
#    mpl.rcParams['ytick.minor.width'] = 2
#    mpl.rcParams['axes.linewidth'] = 2
#    mpl.rcParams['legend.fontsize'] = 18
#    mpl.rcParams['figure.subplot.bottom'] = 0.10    # the bottom of the subplots of the figure
#    mpl.rcParams['figure.subplot.top'] = 0.9    # the bottom of the subplots of the figure
#    mpl.rcParams['figure.subplot.right'] = 0.95    # the bottom of the subplots of the figure
#    mpl.rcParams['figure.subplot.left'] = 0.15    # the bottom of the subplots of the figure
#    mpl.rcParams['figure.subplot.wspace'] = 0.25
#    mpl.rcParams['figure.subplot.hspace'] = 0.25
#    mpl.rcParams['legend.frameon'] = False
#    mpl.rcParams['legend.handlelength'] = 2.5
#
#    fig=plt.subplot(2,1,1)
#
#    plt.plot(np.arange(length)*TimeStep,autocorrXY, color='blue',label='XY') 
#    plt.fill_between(np.arange(length)*TimeStep,autocorrXY-autocorrsigmaXY, autocorrXY+autocorrsigmaXY ,color='blue', alpha=0.15)
#    plt.plot(np.arange(length)*TimeStep,fitfunc(np.arange(length)*TimeStep,fitXY[0][0],fitXY[0][1],fitXY[0][2],fitXY[0][3],fitXY[0][4]), color='blue', linestyle='--', linewidth=2)
#    
#    plt.plot(np.arange(length)*TimeStep,autocorrXZ, color='green',label='XZ') 
#    plt.fill_between(np.arange(length)*TimeStep,autocorrXZ-autocorrsigmaXZ, autocorrXZ+autocorrsigmaXZ, color='green', alpha=0.15)
#    plt.plot(np.arange(length)*TimeStep,fitfunc(np.arange(length)*TimeStep,fitXZ[0][0],fitXZ[0][1],fitXZ[0][2],fitXZ[0][3],fitXZ[0][4]), color='green', linestyle='--', linewidth=2)
#    
#    plt.plot(np.arange(length)*TimeStep,autocorrYZ, color='purple',label='YZ') 
#    plt.fill_between(np.arange(length)*TimeStep,autocorrYZ-autocorrsigmaYZ, autocorrYZ+autocorrsigmaYZ, color='purple', alpha=0.15)
#    plt.plot(np.arange(length)*TimeStep,fitfunc(np.arange(length)*TimeStep,fitYZ[0][0],fitYZ[0][1],fitYZ[0][2],fitYZ[0][3],fitYZ[0][4]), color='purple', linestyle='--', linewidth=2)
#    
#    plt.plot(np.arange(length)*TimeStep,autocorr, color='red',label='Average') 
#    plt.fill_between(np.arange(length)*TimeStep,autocorr-autocorrsigma,autocorr+autocorrsigma,color='red',alpha=0.15)
#    plt.plot(np.arange(length)*TimeStep,fitfunc(np.arange(length)*TimeStep,fit[0][0],fit[0][1],fit[0][2],fit[0][3],fit[0][4]), color='red', linestyle='--', linewidth=2)
#    
#    fig.set_xlabel("Time [fs]")
#    fig.set_ylabel('Autocorrelation')
#
#
#
#    fig=plt.subplot(2,1,2)
#
#    plt.plot(np.arange(length)*TimeStep,np.cumsum(autocorrXY)*TimeStep, color='blue',label='XY') 
#    plt.plot(np.arange(length)*TimeStep, \
#         intfitfunc(np.arange(length)*TimeStep,fitXY[0][0],fitXY[0][1],fitXY[0][2],fitXY[0][3],fitXY[0][4]) \
#         - intfitfunc(0,fitXY[0][0],fitXY[0][1],fitXY[0][2],fitXY[0][3],fitXY[0][4]), \
#         color='blue', linestyle='--', linewidth=2)
#    
#    plt.plot(np.arange(length)*TimeStep,np.cumsum(autocorrXZ)*TimeStep, color='green',label='XZ') 
#    plt.plot(np.arange(length)*TimeStep, \
#         intfitfunc(np.arange(length)*TimeStep,fitXZ[0][0],fitXZ[0][1],fitXZ[0][2],fitXZ[0][3],fitXZ[0][4]) \
#         - intfitfunc(0,fitXZ[0][0],fitXZ[0][1],fitXZ[0][2],fitXZ[0][3],fitXZ[0][4]), \
#         color='green', linestyle='--', linewidth=2)
#    
#    plt.plot(np.arange(length)*TimeStep,np.cumsum(autocorrYZ)*TimeStep, color='purple',label='YZ') 
#    plt.plot(np.arange(length)*TimeStep, \
#         intfitfunc(np.arange(length)*TimeStep,fitYZ[0][0],fitYZ[0][1],fitYZ[0][2],fitYZ[0][3],fitYZ[0][4]) \
#         - intfitfunc(0,fitYZ[0][0],fitYZ[0][1],fitYZ[0][2],fitYZ[0][3],fitYZ[0][4]), \
#         color='purple', linestyle='--', linewidth=2)
#    
#            #plt.plot(np.arange(length)*TimeStep,np.cumsum(autocorr)*TimeStep, color='red',label='Average')
#            #plt.fill_between(np.arange(length)*TimeStep,np.cumsum(autocorr-autocorrsigma), np.cumsum(autocorr+autocorrsigma), color='red', alpha=0.15)
#            #plt.plot(np.arange(length)*TimeStep,intfitfunc(np.arange(length)*TimeStep,fit[0][0],fit[0][1],fit[0][2],fit[0][3],fit[0][4]), color='red', linestyle='--', linewidth=2)
#    
#    fig.set_xlabel("Time [fs]")
#    fig.set_ylabel('Integrated Autocorrelation')
#
#    plotname = UMDname[:-7] + 'pdf'
#    #print ('\n\nthe name of the umd file ',UMDname,' from which the pdf file is',plotname)
#    plt.savefig(plotname, bbox_inches='tight')

    intxy = np.cumsum(autocorrXY)*TimeStep
    intyz = np.cumsum(autocorrYZ)*TimeStep
    intxz = np.cumsum(autocorrXZ)*TimeStep
    Tottime = np.arange(length)*TimeStep
    dataname = UMDname[:-7] + 'visc.dat'
    header = 'Time\tautocorrXY\tautocorrYZ\tautocorrXZ\tInt_of_AutocorrXY\tInt_of_AutocorrYZ)\tInt_of_autocorrXZ\tAverage_Ints\n'
    fp = open(dataname,'w')
    fp.write(header)
    for itime in range(length):
        wline = str(itime) + '\t' + str(autocorrXY[itime]) + '\t' + str(autocorrYZ[itime]) + '\t' + str(autocorrXZ[itime]) + '\t' + str(intxy[itime]) + '\t' + str(intyz[itime]) + '\t' + str(intxz[itime]) + '\t' + str( (intxy[itime]+intyz[itime]+intxz[itime])/3 ) + '\n'
        fp.write(wline)
    fp.close()

    return


        
def main():
    argv = sys.argv[1:]
    umdp.headerumd()
    iterstep = 1
    firststep = 0
    UMDname = ''
    originshift=1
    length=1000000
    temperature = 1000
    try:
        opts, arg = getopt.getopt(argv,"hf:t:i:s:o:l:")
    except getopt.GetoptError:
        print ('viscosity_umd.py -f <umdfile> -i <Initial TimeStep>  -s <Sampling_Frequency> -o <frequency of origin shift> -l <correlation timelength> -t <temperature>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('viscosity_umd.py to compute the viscosity of the fluid from a fit on the self-correlation of stresses')
            print (' viscosity_umd.py -f <umdfile> -i <InitialStep>  -s <Sampling_Frequency> -o <frequency_of_origin_shift> -l <max_correlation_timelength>')
            print ('umdfile = input file with the trajectory, in UMD format.')
            print ('Initial TimeStep from umd file. Default = 0')
            print ('Sampling_Frequency = frequency of sampling the trajectory. Default = 1')
            print ('Frequency_of_origin_shift. Default = 1')
            print ('-l = maximum timelength for the correlation. Default = 1,000,000>')
            print ('-t - temperature. Default = 1000K')
            sys.exit()
        elif opt in ("-f"):
            UMDname = str(arg)
        elif opt in ("-t"):
            temperature = int(arg)
        elif opt in ("-i"):
            firststep = int(arg)
        elif opt in ("-n"):
            iterstep = int(arg)
        elif opt in ("-s"):
            originshift = int(arg)
        elif opt in ("-l"):
            length = int(arg)
    if (os.path.isfile(UMDname)):
        print('The first ',firststep, 'timesteps will be discarded')
        print('The umd.dat file is read every ',iterstep,' timesteps')
#        AllSnapshots = [cr.Lattice]
        AllSnapshots = []
        (AllSnapshots,TimeStep)=umdp.read_stresses_4visc(UMDname)
        
        print('Len(allsnapshots),timestep : ',len(AllSnapshots),TimeStep)
        volume=AllSnapshots[0].cellvolume
        temperature = 0.0
        for ii in range(firststep,len(AllSnapshots)):
            temperature = temperature + AllSnapshots[ii].temperature
#            print('current temperature at step ',ii,' is:', AllSnapshots[ii].temperature)
        temperature = temperature/(len(AllSnapshots)-firststep)
        print('Average temperature from the umd file = ',temperature)
        
        if (length>len(AllSnapshots)):
            print('The length requested (', length,')  is too long for the whole duration of the trajectory(', len(AllSnapshots),')')
            print('Will impose the length of the trajectory')
            length = int(len(AllSnapshots)/2)
        ViscosityAnalysis(AllSnapshots,TimeStep,firststep,originshift,length,temperature,UMDname)
    else:
        print ('the umdfile ',umdfile,' does not exist')
        sys.exit()

if __name__ == "__main__":
   main()
