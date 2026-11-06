#!/usr/bin/python
"""
Created on Mon Aug  14 14:54:47 2023

@author: Kevin Jiguet-Covex
"""

import numpy
import sys
import getopt
import os
import subprocess
import umd_processes_fast as umdpf

class Ex(Exception):
    pass

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def main(argv):
    FileName = ''
    Pattern = ''
    data = []
    SkipSteps=0
    average = 0 
    stdev = 0 
    variance = 0 
    anchor=''
    MaxStep=None
    try:
        opts, arg = getopt.getopt(argv,"hf:p:s:m:",["fFileName","pPattern","sSkipSteps","mMaxStep"])
    except getopt.GetoptError:
        print('average.py -f <FileName> -p <Pattern> -s <SkipSteps> -m <MaxStep>')
    for opt, arg in opts:
        if opt == '-h':
            print('average.py program to extract and average numerical values')
            print('average.py -f <FileName> -p <Pattern> -s <SkipSteps> -m <MaxStep>')
            sys.exit()
        elif opt in ("-f", "--fFileName"):
            FileName = str(arg)
        elif opt in ("-p", "--pPattern"):
            Pattern = str(arg)
            anchor=Pattern.split()[0]
        elif opt in ("-s", "--sSkipSteps"):
            SkipSteps = int(arg)
        elif opt in ("-m", "--mMaxStep"):
            MaxStep = int(arg)


    if(anchor==''):
        raise Ex('No parameter given ; no mean calculated. Use -p to specify a parameter.')
    if os.path.isfile(FileName): 
        try :
            patterns=subprocess.check_output(['grep',Pattern,FileName])
        except :
            raise Ex('Parameter '+ Pattern+" not found in file.")
        (MyCrystal,TimeStep) = umdpf.Crystallization(FileName)
        patternsstr=patterns.decode()
        greps=patternsstr.split('\n')
        unit=None
        if MaxStep ==None :
            imax = len(greps)
        else:
            imax = min(MaxStep,len(greps))
        for isteps in range(SkipSteps+1,imax):
            elems=greps[isteps].split()
            for ii in range(len(elems)):
                if elems[ii] == anchor:
                    for jj in range(ii+1,len(elems)):
                        if is_number(elems[jj]):
                            data.append(float(elems[jj]))
                            break
                if unit==None:
                    unit=elems[-1]
        average = sum(data)/len(data)
        for ii in range(len(data)):
            variance = variance + (data[ii]-average)**2
        variance = variance/len(data)
        stdev = numpy.sqrt(variance)
        return data,average,variance,stdev,TimeStep,unit
        print('Averages over ',len(data),' ares: mean = ',average,' variance = ', variance, ' stdev = ', stdev)
    else:
        raise Ex('File not found')
#        sys.exit()


if __name__ == "__main__":
   main(sys.argv[1:])

