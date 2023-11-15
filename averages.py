#!/usr/bin/env python3
"""
@author: Razvan Caracas, Znais Kobsch
"""

import numpy as np 
import sys
import getopt
import os
import subprocess
import matplotlib.pyplot as plt
import umd_processes_fast as umdpf


def is_number(s):
    try:
        
        float(s)
        return True
    except ValueError:
        return False

def grep_pattern(FileName, Pattern, SkipSteps): 
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
    #print('Averages over ',len(data),' steps for', Pattern, 'are: mean = ',average,' variance = ', variance, ' stdev = ', stdev)
    return data, average, stdev

def plot_function(name,ax,data,average,stdev):
    bleuspecial = '#00ccff'
    x=np.arange(1,len(data)+1,1)
    if name == "Pressure":
        unit = "GPa"
    elif name == "Energy":
        unit = "eV"
    elif name == "Temperature":
        unit = "K"
    else:
        unit = "??"
    plt.minorticks_on()
    ax.plot(x,data,'.', color=bleuspecial)
    ax.plot((0,x[len(x)-1]),(average, average), 'm-', linewidth=4)
    txtaverage = 'Average = ' + str(round(average)) + ' $\pm$ ' + str(round(stdev)).rstrip('0').rstrip('.') + ' ' + unit           #for temperature
    ax.text(0.02,0.01, txtaverage,transform=ax.transAxes, fontsize=12)
    ax.autoscale(enable=True,axis='both',tight=True)


def main(argv):
    FileName = ''
    SkipSteps=0
    umdpf.headerumd()
    try:
        opts, arg = getopt.getopt(argv,"hf:i:",["fFileName,iInitialStep"])
    except getopt.GetoptError:
        print('averages.py -f <UMDFile> -i <InitialStep>')
        sys.exit()
    for opt, arg in opts:
        if opt == '-h':
            print('averages.py -f <UMDFile> -i <InitialStep>')
            print('averages.py program to extract the average and the spread of the numerical values of Pressure, Temperature, and InternalEnergy')
            sys.exit()
        elif opt in ("-f", "--fFileName"):
            FileName = str(arg)
        elif opt in ("-i", "--iInitialStep"):
            SkipSteps = int(arg)
    #print(FileName)
    if (os.path.isfile(FileName)):
        #****Calculations and creation of arrays
        Pressure, Average_P, stdev_P  = grep_pattern(FileName,"Pressure", SkipSteps)
        Temperature, Average_T, stdev_T  = grep_pattern(FileName,"Temperature", SkipSteps)
        Energy, Average_E, stdev_E  = grep_pattern(FileName,"InternalEnergy", SkipSteps)
        arraytowrite = FileName + '\t#PstdTstdEstd\t' + str(len(Pressure)) + '\t' + str(Average_P) + '\t' + str(stdev_P) + '\t' + str(Average_T) + '\t' + str(stdev_T) + '\t' + str(Average_E) + '\t' + str(stdev_E)
        #print('Averages over no.steps, P, stdev_P, T, stdev_T, E, stdev_E',len(Pressure),Average_P,stdev_P,Average_T,stdev_T,Average_E,stdev_E)
        print(arraytowrite)
    else:
        print('No input file or file "',FileName,'" does not exist')
        sys.exit()


if __name__ == "__main__":
   main(sys.argv[1:])




