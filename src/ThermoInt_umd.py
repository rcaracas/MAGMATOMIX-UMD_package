#!/usr/bin/env python3
"""
Authors: Adrien Saurety, Emma Stoutenbourg

This code provide an automatic treatment for thermodynamic integration in VASP.

The user runs a simulation for each lambda value needed, uses the VaspParser-ML.py from UMD package on each OUTCAR, move all the files umd.dat in a common directory 
and runs this code.

"""

####librairies
import numpy as np 
import sys
import getopt
import os
import subprocess
import matplotlib.pyplot as plt
import umd_process as umdp
import math
from scipy.integrate import trapz
import gpflow
import tensorflow as tf
####fundamuntal constants, anychange here may cause the collapse of the universe. Be careful.
kb=1.3806e-23  #boltzman constant eV/K
h=6.626e-34 #planck constant J.s
avogadro=6.0221E23 #avogadro number
ev2j=1.6e-19  #Joule to eV conversion

#### functions

def is_number(s):
    """
    Check if the input s is or not a number.
    
    Parameters
    ----------
    s: input to test

    Returns
    -------
    Bool, True: s is a number
          False: s not a number

    """
    try:
        
        float(s)
        return True
    except ValueError:
        return False

def grep_atolist(FileName): 
    """
    Extract the atomic species in the simulation.
    
    Parameters
    ----------
    FileName : .umd.dat file
        File obtains with parsing OUTCAR files from VASP with VaspParser-ML.py from UMD package.

    Returns
    -------
    type_list : list
        List of atomic mass of each atomic type.
    numberofatoms : list
        List of the number of atoms of each atomic type.

    """
    Pattern='masses'
    patterns=subprocess.check_output(['grep',Pattern,FileName]).decode("utf-8")
    greps=patterns.split('\n')
    elems=greps[0].split()
    mass_list=[]
    for i in elems[1::]:
        mass_list+=[float(i)]
    Pattern='types'
    patterns=subprocess.check_output(['grep',Pattern,FileName]).decode("utf-8")
    greps=patterns.split('\n')
    elems=greps[0].split()
    type_list=[]
    for i in elems[1::]:
        type_list+=[int(i)]

    return type_list, mass_list

def grep_acell(FileName): 
    """
    Extract the volume of the simulated cell.
    
    Parameters
    ----------
    FileName : .umd.dat file
        File obtains with parsing OUTCAR files from VASP with VaspParser-ML.py from UMD package.

    Returns
    -------
    vol : float
         Volume of the cell.
         
    """
    Pattern='acell'
    patterns=subprocess.check_output(['grep',Pattern,FileName]).decode("utf-8")
    greps=patterns.split('\n')
    elems=greps[0].split()
    a=float(elems[1])
    b=float(elems[2])
    c=float(elems[3])
    vol=a*b*c*1E-30
    return vol

def grep_lambda(FileName):
    """
    Extract the lambda -aka coupling constant- used in each simulation.
    
    Parameters
    ----------
    FileName : .umd.dat file
        File obtains with parsing OUTCAR files from VASP with VaspParser-ML.py from UMD package.

    Returns
    -------
    lamb: float
        The coupling constant lambda of the simulation.
        
        
    """
    Pattern='lambda_ThermoInt'
    patterns=subprocess.check_output(['grep',Pattern,FileName]).decode("utf-8")
    greps=patterns.split('\n')
    elems=greps[0].split()
    lamb=float(elems[1])
    return lamb

def grep_pattern(FileName, Pattern, SkipSteps): 
    """
    Extract one property of the system along the simulation.
    
    Parameters
    ----------
    FileName : .umd.dat file
        File obtains with parsing OUTCAR files from VASP with VaspParser-ML.py from UMD package.
    SkipSteps : int
        Number of step to remove for thermalization

    Returns
    -------
    average: float
        average value of the property
    stdev: float
        standard deviation of the property
    """
    data = []
    average = 0 
    stdev = 0 
    anchor=Pattern.split()[0]
    patterns=subprocess.check_output(['grep',Pattern,FileName]).decode("utf-8")
    greps=patterns.split('\n')
    for isteps in range(SkipSteps+1,len(greps)):
        elems=greps[isteps].split()
        for ii in range(len(elems)):
            if elems[ii] == anchor:
                for jj in range(ii+1,len(elems)):
                    if is_number(elems[jj]):
                        data.append(float(elems[jj]))
                        break
    average = sum(data)/len(data)
    stdev = np.std(data)
    return average, stdev

def getmolmass(FileName_list):
    """
    Extract the molar mass and the number of atom of the simulation
    
    Parameters
    ----------
    FileName : .umd.dat file
        File obtains with parsing OUTCAR files from VASP with VaspParser-ML.py from UMD package.
   
    Returns
    -------
    M: float
        molar mass of the simu
    numberofatoms: int
        number of atom in the simulation
    
    """
    tl, ml= grep_atolist(FileName_list[0])
    M=0
    numberofatoms = 0
    for i in range(len(tl)):
        M+=tl[i]*ml[i]
    for i in range(len(tl)):
        numberofatoms+=tl[i]
    return M, numberofatoms

def ordering(lambda_list,av_toten_list,std_toten_list):
    '''
    Order the list of average and std of electronic energy by increasing lambda.

    Parameters
    ----------
    lambda_list : list
        list of the lambda of the simulations
    av_toten_list : list
        list of the average electronic energy 
    std_toten_list : list
        list of the std electronic energy 

    Returns
    -------
    llf : TYPE
       list of the lambda of the simulations ordered
    tlf : TYPE
        list of the average electronic energy ordered
    slf : TYPE
        list of the std electronic energy ordered

    '''
    llf=sorted(lambda_list)
    tlf=[]
    slf=[]
    for i in range(len(llf)):
        lamb=llf[i]
        for j in range(len(llf)):
            if lamb==lambda_list[j]:
                tlf+=[av_toten_list[j]]
                slf+=[std_toten_list[j]]
    return llf, tlf, slf


def ig_free_energy(FileName_list, SkipSteps, type_list, mass_list,T):
    """
    Extract the lambda -aka coupling constant- used in each simulation.
    
    Parameters
    ----------
    FileName : .umd.dat file
        File obtains with parsing OUTCAR files from VASP with VaspParser-ML.py from UMD package.
    SkipSteps : int
        Number of step to remove for thermalization
    type_list : list
        list of the number of atoms of each atomic type
    mass_list : list
        list of the mass of each atomic type
    T : float
        average temperature of the simulation

    Returns
    -------
    f0 the free energy of the ideal gas corresponding to our system.
        
    """
    vol=grep_acell(FileName_list[0])
    f0=0
    for i in range(len(type_list)):
        N=type_list[i] #number of atom of first type
        m=(mass_list[i]/1000)*(1/avogadro)
        BWL=h/(np.sqrt(2*3.1416*m*kb*T)) #Broglie Wave Length
        f0+=((-1)*kb*T*(N*math.log(vol/(BWL**3))-math.log(math.factorial(N))))*6.24e18
    print("Ideal gas free energy", f0," eV")
    return f0

def getET(FileName_list, SkipSteps):
    """
    Extract the internal energy and temperature of the ab initio simulation (i.e.  the one with lambda=1)
    
    Parameters
    ----------
    FileName : .umd.dat file
        File obtains with parsing OUTCAR files from VASP with VaspParser-ML.py from UMD package.
    SkipSteps : int
        Number of step to remove for thermalization

    Returns
    -------
    E : float
        internal energy
    T : float
        temperature
    """
    for f in FileName_list:
        if grep_lambda(f)==1.0:
            E,stdE=grep_pattern(f,'InternalEnergy',SkipSteps)
            T,stdT=grep_pattern(f,'Temperature',SkipSteps)
    return E,T

def integral_varchange_gp(FileName_list, SkipSteps):
    '''
    Integrate the value of toten of all the umd file

    Parameters
    ----------
    FileName_list : list of .umd.dat files
        File obtains with parsing OUTCAR files from VASP with VaspParser-ML.py from UMD package.
    SkipSteps : int
        Number of step to remove for thermalization

    Returns
    -------
    inte : float
        value of the integrale
    var : 
        variance of the integral (used to compute the errorbar on the free energy)

    '''
    lambda_list=[]
    av_toten_list=[]
    std_toten_list=[]
    inte=0
    k=0.8
    for f in FileName_list:
        lambda_list+=[grep_lambda(f)]
        toten=grep_pattern(f,'TotalElectronicEnergy',SkipSteps)
        av_toten_list+=[toten[0]]
        std_toten_list+=[toten[1]]
    lambda_list, av_toten_list,std_toten_list=ordering(lambda_list, av_toten_list, std_toten_list)
    plt.figure()
    plt.errorbar(lambda_list,av_toten_list,yerr=std_toten_list,marker='o',linestyle='')
    plt.xlabel('$\lambda$')
    plt.ylabel('internal energy (eV)')
    plt.savefig('TI_value.pdf')
    try:
        x_list=[-1]+[2*(i**(1-k))-1 for i in lambda_list]
        av_along_x_list=[0.0]+[i*(j**k)*(1/(2*(1-k))) for i,j in zip(av_toten_list,lambda_list)]
        plt.figure()
        X=np.c_[x_list]
        Y=np.c_[av_along_x_list]
        plt.plot(X[0], Y[0], "kx", mew=2,)
        plt.plot(X[1], Y[1], "kx", mew=2,)
        model = gpflow.models.GPR(
            (X, Y),
            gpflow.kernels.Matern32(variance=30000),
            mean_function=gpflow.functions.Polynomial(3),
            )
        opt = gpflow.optimizers.Scipy()
        opt.minimize(model.training_loss, model.trainable_variables)
        nb_pts=1000
        Xplot = np.linspace(-1.0, 1.0, nb_pts)[:, None]

        f_mean, f_var = model.predict_f(Xplot, full_cov=False)
        y_mean, y_var = model.predict_y(Xplot)

        f_lower = f_mean - 1.96 *np.sqrt(f_var)
        f_upper = f_mean + 1.96* np.sqrt(f_var)

        plt.plot(X[2::], Y[2::], "ro", mew=2, label="data")
        plt.plot(Xplot, f_mean, "-", color="C0", label="mean")
        plt.plot(Xplot, f_lower, "--", color="C0", label="f 95% confidence")
        plt.plot(Xplot, f_upper, "--", color="C0")
        plt.fill_between(
            Xplot[:, 0], f_lower[:, 0], f_upper[:, 0], color="C0", alpha=0.1
        )
        plt.legend()
        plt.savefig('TI_cv.pdf')
        fm=tf.Variable(f_mean).numpy()
        fm=np.reshape(fm,nb_pts)
        x=tf.Variable(Xplot).numpy()
        x=np.reshape(x,nb_pts)
        inte=trapz(fm,x)
        fmu=tf.Variable(f_upper).numpy()
        fmu=np.reshape(fmu,nb_pts)
        var=trapz(fmu,x)-inte
        return inte, var
    except:
        print("ERROR!! you need at least one simulation with lambda below 0.1 !!!")

def free_energy(FileName_list, SkipSteps,T):
    '''
   Integrate the value of toten of all the umd file

   Parameters
   ----------
   FileName_list : list of .umd.dat files
       File obtains with parsing OUTCAR files from VASP with VaspParser-ML.py from UMD package.
   SkipSteps : int
       Number of step to remove for thermalization

   Returns
   -------
   f : list of [inte, var]
   with 
       inte : float
           value of the integrale
       var : float
           variance of the integral (used to compute the errorbar on the free energy)

    '''
    type_list, mass_list=grep_atolist(FileName_list[0])
    f0=ig_free_energy(FileName_list, SkipSteps,type_list, mass_list,T)
    inte,var=integral_varchange_gp(FileName_list, SkipSteps)
    f=[f0+inte,f0+inte+var]
    print("Free energy of coupling AIMD to IG", inte,"eV")
    return f

def entropy_calc(F,E,T):
    '''
    compute the entroy 

    Parameters
    ----------
    F : float
        Free energy
    E : float
        internal energy
    T : Tfloat
        temperature

    Returns
    -------
    S : float
        entropy

    '''
    S=np.zeros(len(F))
    for k in range(len(F)):
        S[k]=(E-F[k])/T
    return S

def atounit2natunitmol(f,numberofatoms):
    '''
    change the unit from atomic unit to natural unit

    Parameters
    ----------
    f : float
        free energy
    numberofatoms : list
        number of atom of each type

    Returns
    -------
    list
        free energy and free energy error in natural unit (kJ/atom)
    '''
    F=f[0]
    Ferror=f[1]-f[0]
    F=F*ev2j*avogadro/numberofatoms/1000
    Ferror=Ferror*ev2j*avogadro/numberofatoms/1000
    return [F, Ferror]

def atounit2natunitkg(f,M):
    '''
    change the unit from atomic unit to SI unit (J/kg)

    Parameters
    ----------
    f : float
        free energy
    M: float
        molar mass

    Returns
    -------
    list
        free energy and free energy error in SI unit (kJ/kg)

    '''
    F=f[0]
    Ferror=f[1]-f[0]
    F=F*ev2j*avogadro/(M/1000)
    Ferror=Ferror*ev2j*avogadro/(M/1000)
    return [F, Ferror]

def main(argv):
    FileName_list = []
    FileName_list += [each for each in os.listdir() if each.endswith('.umd.dat')]
    SkipSteps= 1
    unit=False
    gui=True
    umdp.headerumd()
    
    try:
        opts, arg = getopt.getopt(argv,"h:f:i:u:g",["fFileName,iInitialStep,uUnit,gGUI"])
    except getopt.GetoptError:
        print('')
        sys.exit()
    for opt, arg in opts:
        if opt == '-h':
            print('')
            sys.exit()
        if opt in ("-f", "--fFileName"):
            FileName_list= arg
        if opt in ("-i", "--iInitialStep"):
            SkipSteps= arg
        if opt in ("-u", "--uUnit"):
            try:
                arg=bool(arg)
                unit=arg
            except ValueError:
                print('unvalid unit, autoswitch to default (SI unit)')
        if opt in ("-g", "--gGUI"):
            try:
                arg=bool(arg)
                gui=arg
            except ValueError:
                print('unvalid command, autoswitch to default (plot not shown)')
    

    M, numberofatoms=getmolmass(FileName_list)
    E,T=getET(FileName_list,500)
    print('internal energy ab initio system = ',E,' eV')
    print('average temperature in the full ab initio system = ',T,' K')
    F=free_energy(FileName_list, SkipSteps,T)
    S=entropy_calc(F,E,T)
    if unit:
        print('Free Energy =',atounit2natunitmol(F,M, numberofatoms)[0],'+/-',abs(atounit2natunitmol(F,M,numberofatoms)[1]),'kJ/mol')
        print('Entropy =',atounit2natunitmol(S,M, numberofatoms)[0],'+/-',abs(atounit2natunitmol(S,M,numberofatoms)[1]),'kJ/K/mol')
    else:
        print('Free Energy =',atounit2natunitkg(F,M)[0],'+/-',abs(atounit2natunitkg(F,M)[1]),'J/kg')
        print('Entropy =',atounit2natunitkg(S,M)[0],'+/-',abs(atounit2natunitkg(S,M)[1]),'J/K/kg')
    if gui:
        plt.show()
        


   

if __name__ == "__main__":
   main(sys.argv[1:])


