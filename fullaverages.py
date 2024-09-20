#!/usr/bin/env python3
###
##AUTHORS: ANAIS KOBSCH, RAZVAN CARACAS
###


import glob
import sys
import getopt
import subprocess
import numpy as np
import crystallography as cr
import umd_processes_fast as umdpf



def selection_value(list_sigma, list_err):
    """ Selection of the value of the error on the mean that converged, or taking the highest value """
    level_convergence = 6 #number of block transformation we want the sigme to be constant to say it converged
    stable = 0
    breakvalue = 0
    finalbreakvalue = 0
    #we loop throught all the block transformations
    for i in range(len(list_sigma)-level_convergence): 
        #and for each block transfor, we test the convergence with the next one
        for ii in range(i,len(list_sigma)-1):
            #we loop through all the values in the error range around the sigma value[i]
            for j in np.arange(list_sigma[i]-list_err[i],list_sigma[i]+list_err[i]+list_err[i]/10,list_err[i]/10):
                #if this value is in the range of the next value (i+1), then we add 1 to the counter of stability
                if j <= list_sigma[ii+1]+list_err[ii+1] and j >= list_sigma[ii+1]-list_err[ii+1]:
                    stable += 1
                    breakvalue = 1
                    break
            #we test if we exited the last loop throught a break or not
            #if we did not encounter a break, then the stable value is reset to 0 (we did not start to converge) and we move to the next i
            if breakvalue == 0:
                stable = 0
                break
            else:
                breakvalue = 0
                #if the sigma value converged, then we record the value and break the loop on block transformations
                if stable >= level_convergence:
                    final_sigma = '{:1.1e}'.format(list_sigma[i+1]) 
                    final_err = '{:1.1e}'.format(list_err[i+1])
                    finalbreakvalue = 1
                    #print('encounter finalbreak')
                    break
        if finalbreakvalue == 1:
            break
    #we test if we exited the loop with a break or not
    #if we did not encounter a break, then the final value taken is the max of the list and we indicate a > sign in the string
    if finalbreakvalue == 0:
        final_sigma = '>' + '{:1.1e}'.format(max(list_sigma))
        final_err = '>' + '{:1.1e}'.format(max(list_err))
    return final_sigma, final_err

def blocking_transformation(nsteps,data):
    """ Transformation of a data set into a new one, of half length by averaging data two-by-two"""
    new_data = []
    new_nsteps = 0
    for i in range(0,nsteps//2):
        new_data.append((data[2*i] + data[2*i+1])/2 )
    if nsteps%2 == 1:  #=0 if nsteps is even, =1 if nsteps is odd --> we add the last data point
        new_data.append(data[-1])
    new_nsteps = len(new_data)
    return new_nsteps, new_data
        
def blocking_method(data,variance):
    """Compute the sequence of average on the mean and its error using blocking method"""
    nsteps = len(data)
    #initialization of the sequences
    list_sigma = [(variance/(nsteps - 1))**(1/2)]
    list_err  = [1/( (2*(nsteps-1))**(1/2) )]
    #loop while the size of the new data set is above 2
    while nsteps > 2:
        #transformation of the data and nsteps using blocking method
        nsteps, data = blocking_transformation(nsteps,data)              
        #calculation of c0/(n-1) and addition to the sequence
        average = sum(data)/nsteps
        variance = 0
        for ii in range(nsteps):
            variance = variance + (data[ii]-average)**2
        variance = variance/nsteps
        list_sigma.append((variance/(nsteps-1))**(1/2))
        list_err.append(1/( (2*(nsteps-1))**(1/2) ))
    return list_sigma, list_err

def convert_to_float(string):
    if string[0] == '>':
        value = float(string[1:])
    else:
        value = float(string[:])
    return value 

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def grep_pattern(FileName, Pattern, SkipSteps):
    """Extraction of values with grep and creation of corresponding list of data"""
    data = []
    average = 0 
    stdev = 0 
    variance = 0
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
    for ii in range(len(data)):
        variance = variance + (data[ii]-average)**2
    variance = variance/len(data)
    stdev = np.sqrt(variance)
    return len(data), data, average, variance, stdev


def headerfile(firstfile,units):
    """ creation of the newfile with correct header """
    with open(firstfile, 'r') as ff:
        while True:
            line = ff.readline()
            if len(line) > 1:
                line=line.strip()
                entry=line.split()
                if entry[0] == 'natom':
                    natom = int(entry[1])
                if entry[0] == 'types':
                    types = entry[1:]
                if entry[0] == 'elements':
                    elements = entry[1:]
                    break
    if units == 0:
        firstline = 'file\tnsteps\tDensity(g/cm3)\tstdev_rho\terr_rho\tPressure(GPa)\tstdev_P\terr_P\tTemperature(K)\tstdev_T\terr_T\tInternEnergy(eV/atm)\tstdev_E\terr_E\tCvm(Nkb)\tstdev_Cvm\n'
    else:
        firstline = 'file\tnsteps\tDensity(g/cm3)\tstdev_rho\terr_rho\tPressure(GPa)\tstdev_P\terr_P\tTemperature(K)\tstdev_T\terr_T\tTemperature(eV)\tstdev_T\terr_T\tInternEnergy(eV/atm)\tstdev_E\terr_E\tInternEnergy(eV/unitcell)\tstdev_E\terr_E\tInternEnergy(J/g)\tstdev_E\terr_E\tCvm(Nkb)\tstdev_Cvm\tCvm(J/K/mol)\tstdev_Cvm\tCvm(J/K/g)\tstdev_Cvm\tCv(J/K/unitcell)\tstdev_Cv\n'
    # creation of the header
    newfilename = 'fullthermo.dat'
    print('The file ',newfilename,' is created')
    f = open(newfilename,'w')
    f.write("elements\t"+"\t".join(x for x in elements)+ "\n")
    f.write("number\t"+"\t".join(x for x in types)+ "\n")
    f.write(firstline)
    return f, natom      #I return the newly created files f along with natom


def main(argv):
    """     ********* Main program *********     """
    SkipSteps=0
    units = 0
    #parameters
    Na=6.022e23            #Avogadro constant
    eVtoJ = 1.6e-19        #1eV = 1.6e-19 J
    KtoeV = 0.000086173324 #1K = 8.62E-5 eV
    kb = 1.38064852e-23    #boltzmann constant
    umdpf.headerumd()
    try:
        options, arg = getopt.getopt(argv,"hi:u:",["iInitialStep","units"])
    except getopt.GetoptError:
        print('fullaverages.py -i <InitialStep> -u <units version: 0 = minimal (default), 1 = full units>')
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('fullaverages.py -i <InitialStep> -u <units version: 0 = minimal (default), 1 = full units>')
            print('Program to extract and average thermodynamic data from all the umd files in the current folder')
            print('As ouptut it produces a fullthermo.dat file. On each line there are specified in order: ')
            print('    - the umd filename, the number of snapshots')
            print('    - 3 columns per variable (rho,P,T,E)')
            print('        - average')
            print('        - stdev of data to the mean')
            print('        - statistical error on the mean')
            print('    - 2 columns for the Cv(Nkb) with the calculated value and statistical error')
            print('WARNING: please note that the Cv formulation is valdid *only* for NVT simulations')
            print(' ')
            sys.exit()
        elif opt in ("-i", "--iInitialStep"):
            SkipSteps = int(arg)
        elif opt in ("-u", "--units"):
            units = int(arg)
    #Initialization
    files = sorted(glob.glob('*.umd.dat'))                             #list every umd files in alphabetic order
    f, natom = headerfile(files[0], units)             #I create the first newfile for averages    
    #Calculation for each file and writing of the newfile
    for file in files:
        results = [] #initialization of the final result array, used to print each line of result in the new file
        is_E = 0 #initialization of the indicator of energy (in order to compute the Cv)
        is_KE = 0
        print('Averaging in file',file)
        results.append(file)
        #********** Extraction of the Data for the Density, Pressure, Temperature and Energy
        # In each case we try to:
        #     1) extract the data from the umd, create the array and compute the average and stendard deviation
        #     2) compute the sequence of values obtained from (c0/(n-1))**(1/2) (see blocking method in Flyvbjerg1989)
        #     3) select the values to write in the file
        #
        #**** Density
        try:
            nsteps_rho, Density, Average_rho, Variance_rho, stdev_rho  = grep_pattern(file,"Density", SkipSteps)
            list_sigma_rho, list_err_rho = blocking_method(Density,Variance_rho)
            final_sigma_rho, final_err_rho = selection_value(list_sigma_rho, list_err_rho)
        except subprocess.CalledProcessError:
            print ('    WARNING !!! Density is not defined in the UMD file')
            list_sigma_rho, list_err_rho = (np.zeros(3),np.zeros(3))
            Average_rho,Variance_rho,stdev_rho,final_sigma_rho = 0,0,0,0
        #**** Pressure
        try:
            nsteps_P, Pressure, Average_P, Variance_P, stdev_P  = grep_pattern(file,"Pressure", SkipSteps)
            list_sigma_P, list_err_P = blocking_method(Pressure,Variance_P)
            final_sigma_P, final_err_P = selection_value(list_sigma_P, list_err_P)
        except subprocess.CalledProcessError:
            print ('    WARNING !!! Pressure is not defined in the UMD file')
            nsteps_P, Pressure, Average_P, Variance_P, stdev_P  = (0.0,0.0,0.0,0.0,0.0)
            list_sigma_P, list_err_P = (np.zeros(3),np.zeros(3))
            Average_P, Variance_P, stdev_P,final_sigma_P = 0,0,0,0
        #**** Temperature
        try:
            nsteps_T, Temperature, Average_T, Variance_T, stdev_T  = grep_pattern(file,"Temperature", SkipSteps)
            list_sigma_T, list_err_T = blocking_method(Temperature,Variance_T)
            final_sigma_T, final_err_T = selection_value(list_sigma_T, list_err_T)
        except subprocess.CalledProcessError:
            print ('    WARNING !!! Temperature is not defined in the UMD file')
            list_sigma_T, list_err_T = (np.zeros(3),np.zeros(3))
            Average_T, Variance_T, stdev_T,final_sigma_T = 0,0,0,0
            
        #**** Energies
        try:
            #this energy does not include the kinetic energy of ions, it is the energy without entropy of VASP
            nsteps_EWE, EnergyWithoutEntropy, Average_EWE, Variance_EWE, stdev_EWE  = grep_pattern(file,"InternalEnergy", SkipSteps)
            list_sigma_EWE, list_err_EWE = blocking_method(EnergyWithoutEntropy,Variance_EWE)
            final_sigma_EWE, final_err_EWE = selection_value(list_sigma_EWE, list_err_EWE)
            is_E = 1
        except subprocess.CalledProcessError:
            print ('    WARNING !!! InternalEnergy is not defined in the UMD file')
            nsteps_EWE, EnergyWithoutEntropy, Average_EWE, Variance_EWE, stdev_EWE  = (0.0,0.0,0.0,0.0,0.0)
            list_sigma_EWE, list_err_EWE = (np.zeros(3),np.zeros(3))
        try:
            nsteps_K, KineticEnergy, Average_K, Variance_K, stdev_K  = grep_pattern(file,"KineticEnergyIons", SkipSteps)
            is_KE = 1
        except subprocess.CalledProcessError:
            print ('    WARNING !!! KineticEnergy is not defined in the UMD file')
        #the internal energy as used in the Hugoniot calculation is the internal energy which includes the kinetic energy of ions
        if is_E == 1 and is_KE == 1:
            Energy = [KineticEnergy[ii] + EnergyWithoutEntropy[ii] for ii in range(len(KineticEnergy))]
            nsteps_E = len(Energy)
            Average_E = sum(Energy)/nsteps_E
            Variance_E = 0
            for ii in range(nsteps_E):
                Variance_E = Variance_E + (Energy[ii]-Average_E)**2
            Variance_E = Variance_E/nsteps_E
            stdev_E = np.sqrt(Variance_E)
            list_sigma_E, list_err_E = blocking_method(Energy,Variance_E)
            final_sigma_E, final_err_E = selection_value(list_sigma_E, list_err_E)
        else:
            print ('    WARNING !!! InternalEnergy does not include the kinetic energy of ions!!')
            nsteps_E, Energy, Average_E, Variance_E, stdev_E  = nsteps_EWE, EnergyWithoutEntropy, Average_EWE, Variance_EWE, stdev_EWE 
            list_sigma_E, list_err_E = list_sigma_EWE, list_err_EWE
            final_sigma_E, final_err_E = selection_value(list_sigma_E, list_err_E)
        #**** For the calculation of Cv and conversions we need the mass of the cell and natoms
        #--> we read the header of the umd file using umd_process store the elements information in the class MyCrystal
        (MyCrystal,TimeStep) = umdpf.read_bigheader_umd(file)
        natom = MyCrystal.natom
        #***** Calculation of Cv and statistical error using the bootstrap method
        #*** The calculation of Cv depends on the data we have (E, KE)
        # If we have E but not KE, then we take the ideal gas kinetic part to compute the Cv
        if is_E == 1 and is_KE == 0: 
            print ('    For information: Cv is computed using the ideal gas kinetic part')
            Cv_bootstrap = []
            for n in range(1000):
                #print('*****Bootstrap N°',n)
                Energy_rand = np.random.choice(Energy,nsteps_E,replace = True)
                variance_Erand = 0
                variance_Erand = np.mean(Energy_rand**2) - np.mean(Energy_rand)**2
                Cv = variance_Erand * (eVtoJ)**2 / (kb * Average_T**2)   + 3/2 * natom * kb #in J/unitcell/K
                Cv_bootstrap.append(Cv)
            Cv_bootstrap = np.asarray(Cv_bootstrap)    
            Cv = np.mean(Cv_bootstrap) #an estimator of the mean of Cv, in J/unitcell/K
            Cv_stdev = np.sqrt(np.mean(Cv_bootstrap**2) - np.mean(Cv_bootstrap)**2 )  #the error on the previous mean, in J/unitcell/K
        # If we have E and KE, then we take both variances to compute the Cv
        elif is_E == 1 and is_KE == 1: 
            print ('    For information: Cv is computed using both internal and kinetic energy variances')
            Cv_bootstrap = []
            for n in range(1000):
                #print('*****Bootstrap N°',n)
                Energy_rand = np.random.choice(Energy,nsteps_E,replace = True)
                variance_Erand = 0
                variance_Erand = np.mean(Energy_rand**2) - np.mean(Energy_rand)**2
                Cv = variance_Erand  * (eVtoJ)**2 / (kb * Average_T**2)   #in J/unitcell/K
                Cv_bootstrap.append(Cv)
            Cv_bootstrap = np.asarray(Cv_bootstrap)    
            Cv = np.mean(Cv_bootstrap) #an estimator of the mean of Cv, in J/unitcell/K
            Cv_stdev = np.sqrt(np.mean(Cv_bootstrap**2) - np.mean(Cv_bootstrap)**2 )  #the error on the previous mean, in J/unitcell/K
        #if we do not have the energy, then we do not compute the Cv
        else:
            print ('    WARNING !!! Cv is not computed')
            Cv, Cv_stdev = (0,0)
        Cvm = Cv * Na / natom   #in J/K/mol
        Cvm_stdev = Cv_stdev * Na / natom  #in J/K/mol
        Cvm_Nkb = Cv / (kb * natom)         #in Nakb
        Cvm_Nkb_stdev = Cv_stdev / (kb * natom) #in Nakb

        #****** Write result line
        #**** For the conversions we need the mass of the cell 
        #--> we read the header of the umd file using umd_process store the elements information in the class MyCrystal
        (MyCrystal,TimeStep) = umdpf.read_bigheader_umd(file)
        mass = 0
        for itypat in range(MyCrystal.ntypat):
            (atomicname,atomicsymbol, atomicnumber,MyCrystal.masses[itypat])=cr.Elements2rest(MyCrystal.elements[itypat])
            mass = mass + MyCrystal.masses[itypat] * MyCrystal.types[itypat]
        mass = mass / Na
        #convert eV/unitcell to eV/atom or to J/g
        if final_sigma_E[0] == '>':
            final_sigma_E_conv = '>'+'{:1.1e}'.format(float(final_sigma_E[1:])/natom)
            final_sigma_E_mass = '>'+'{:1.1e}'.format(float(final_sigma_E[1:])*eVtoJ / mass)
        else:
            final_sigma_E_conv = '{:1.1e}'.format(float(final_sigma_E[:])/natom)    
            final_sigma_E_mass = '{:1.1e}'.format(float(final_sigma_E[:])*eVtoJ / mass)    
        if units == 0:           
            results.extend([str(nsteps_P),'{:1.2f}'.format(Average_rho), '{:1.1e}'.format(stdev_rho), str(final_sigma_rho),'{:1.2e}'.format(Average_P),'{:1.1e}'.format(stdev_P),str(final_sigma_P), '{:1.0f}'.format(Average_T), '{:1.0f}'.format(stdev_T), str(final_sigma_T), '{:1.2e}'.format(Average_E/natom), '{:1.1e}'.format(stdev_E/natom), str(final_sigma_E_conv), '{:1.2f}'.format(Cvm_Nkb), '{:1.1e}'.format(Cvm_Nkb_stdev)])
        else: #additional conversion of units
            #conversion of eV/unitcell to J/g
            Average_E_mass = Average_E *eVtoJ / mass
            stdev_E_mass = stdev_E *eVtoJ  / mass
            #conversion of J/K to J/K/g
            Cvm_mass = Cv / mass
            Cvm_mass_stdev = Cv_stdev / mass
            #conversion K to eV
            if final_sigma_T[0] == '>':
                final_sigma_T_conv = '>'+'{:1.1e}'.format(float(final_sigma_T[1:])*KtoeV)
            else:
                final_sigma_T_conv = '{:1.1e}'.format(float(final_sigma_T[:])*KtoeV)
            results.extend([str(nsteps_P),'{:1.2f}'.format(Average_rho), '{:1.1e}'.format(stdev_rho), str(final_sigma_rho),'{:1.2e}'.format(Average_P),'{:1.1e}'.format(stdev_P),str(final_sigma_P), '{:1.0f}'.format(Average_T), '{:1.0f}'.format(stdev_T), str(final_sigma_T), '{:1.2e}'.format(Average_T*KtoeV), '{:1.1e}'.format(stdev_T*KtoeV), str(final_sigma_T_conv), '{:1.2e}'.format(Average_E/natom), '{:1.1e}'.format(stdev_E/natom), str(final_sigma_E_conv),   '{:1.3e}'.format(Average_E), '{:1.1e}'.format(stdev_E), str(final_sigma_E),   '{:1.3e}'.format(Average_E_mass), '{:1.1e}'.format(stdev_E_mass), str(final_sigma_E_mass),   '{:1.2f}'.format(Cvm_Nkb), '{:1.1e}'.format(Cvm_Nkb_stdev),   '{:1.2e}'.format(Cvm), '{:1.1e}'.format(Cvm_stdev),  '{:1.2e}'.format(Cvm_mass), '{:1.1e}'.format(Cvm_mass_stdev), '{:1.2e}'.format(Cv), '{:1.1e}'.format(Cv_stdev)])
        f.write("\t".join(x for x in results)+ "\n")
    f.close()
    

#   ********* Execution *********   
if __name__ == "__main__":
    main(sys.argv[1:])

