###
#AUTHORS: ANAIS KOBSCH
###



"""     ********* Importation of the packages and modules used here *********     """
import sys
import getopt
import glob
import os
import numpy as np
from scipy import stats


def headerfile(firstfile, dirpath, headfile):
    """creation of the newfile with correct header"""
    firstline = ['atom']  #beginning of the first line of the file gofr
    secondline = ['file']
    # creation of the header from the first line of the first file
    with open(firstfile, 'r') as f:
        line = f.readline()
    atoms = line.strip('time_(fs)').split() #list of atoms
    for atom in atoms:
        for i in range(0,5):
            firstline.append(atom)
        secondline.extend(['D(m2/s)','D_stdev','R_squared','slope','intercept'])
    newfilename = dirpath+'/diffusivities.txt' 
    print('The file ',newfilename,' is created')
    f = open(newfilename,'w')
    if headfile != '':
        with open(headfile,'r') as hf:
            for line in hf:
                f.write(line)
    f.write("\t".join(x for x in firstline)+ "\n")
    f.write("\t".join(x for x in secondline)+ "\n")
    return f, atoms      #I return the newly created files f along with the list of element couples


def calculation_diffusivities(t,data,steps):
    """     ********* Calculation of the self diffusivities (using linear regression) *********     """
    i = 0
    while t[i] < steps:
        i = i+1
    NewTemps = t[i:]
    NewData = data[i:]
    #least squared linear regression giving the correlation coefficient
    slope, intercept, r_value, p_value, std_err = stats.linregress(NewTemps, NewData) #intercept = y(x=0), std_err = error on the slope
    R_squared = r_value**2
    diffusivity = slope * 10**(-5) / 6.
    diff_err = std_err * 10**(-5) / 6.
    return diffusivity , diff_err , R_squared , slope , intercept



def main():
    argv = sys.argv[1:]
    """     ********* Main program *********     """
    SkipTemps = 0
    headfile = ''
    try:
        options,arg = getopt.getopt(argv,"hs:f:",["sSkipTemps","file"])
    except getopt.GetoptError:
        print('plot-msd.py -s <SkipTemps>(time in fs, default=0) -f <header file (default none)>')
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('analyze_msd.py program to compute self diffusivities from msd.dat file and write them in a .txt file per subfolder')
            print('analyze_msd.py -s <SkipTemps>(time in fs, default=0) -f <header file (default none)>')
            print('')
            print('headerfile is a txt file containing the element and number (useful to compute the densities)')
            print('')
            print('This code produces a diffusivities.txt file with:')
            print('     - a column with filename')
            print('     - 5 columns per atom type')
            print('         - Diffusivity (m^2/s)')
            print('         - stdev on diffusivity (m^2/s)')
            print('         - R^2  coefficient associated to the fit')     
            print('         - slope of the fit') 
            print('         - intercept (y(x=0))') 
            sys.exit()
        elif opt in ("-s","--sSkipTemps"):
            SkipTemps = int(arg)
        elif opt in ('-f','--file'):
            headfile = str(arg)
    for dirpath, dirnames, filenames in os.walk(os.curdir):
        files = sorted(glob.glob(dirpath+'/*.msd.dat')) #I list every gofr files in alphabetic order
        if files != []:
            f, atoms = headerfile(files[0], dirpath, headfile)                          #I create the first newfile for gofr and save the list of element couples 
            for file in files:
                results = [file]
                for atom in atoms:
                    data = np.loadtxt(file,skiprows=1,usecols=atoms.index(atom)+1,unpack=True)
                    Temps = np.loadtxt(file,skiprows=1,usecols=0,unpack=True)
                    D, D_stdev, R, slope, intercept = calculation_diffusivities(Temps,data,SkipTemps) #we compute the diffusivities and stdev using linear fit
                    results.extend([str(D), str(D_stdev), str(R), str(slope), str(intercept)])
                f.write("\t".join(x for x in results)+ "\n")                  #we write in the file the result line
            f.close



if __name__ == "__main__":
   main()
