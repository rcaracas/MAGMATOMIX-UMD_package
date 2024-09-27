#!/usr/bin/env python3
"""
Created on Tue May  2 14:54:47 2022

@author: Tim Bögels
"""

import numpy
import sys
import getopt
import os
import subprocess


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
    try:
        opts, arg = getopt.getopt(argv,"hf:p:s:",["fFileName","pPattern","sSkipSteps"])
    except getopt.GetoptError:
        print('average.py -f <FileName> -p <Pattern> -s <SkipSteps> ')
        sys.exit()
    for opt, arg in opts:
        if opt == '-h':
            print ('average.py program to extract and average numerical values')
            print ('average.py -f <FileName> -p <Pattern> -s <SkipSteps>')
            sys.exit()
        elif opt in ("-f", "--fFileName"):
            FileName = str(arg)
        elif opt in ("-p", "--pPattern"):
            Pattern = str(arg)
            anchor=Pattern.split()[0]
        elif opt in ("-s", "--sSkipSteps"):
            SkipSteps = int(arg)
    if (anchor==''):
        print('No parameter given ; no mean calculated. Use -p to specify a parameter.')
        sys.exit()
    if (os.path.isfile(FileName)): 
        try :
            patterns=subprocess.check_output(['grep',Pattern,FileName])
        except subprocess.CalledProcessError :
            print('Parameter _',Pattern,'_ not found. You may try a different name.')
            sys.exit()
        patternsstr=patterns.decode()
        greps=patternsstr.split('\n')
        for isteps in range(SkipSteps+1,len(greps)):
            elems=greps[isteps].split()
            for ii in range(len(elems)):
                if elems[ii] == anchor:
                    for jj in range(ii+1,len(elems)):
                        if is_number(elems[jj]):
                            if(jj==len(elems)-1):
                                data.append(float(elems[jj]))
                                break
                            else :
                                StrTensor=elems[jj:]
                                Tensor=[]
                                for elem in StrTensor :
                                    if is_number(elem):
                                        Tensor.append(float(elem))
                                data.append(Tensor)
                                break
        if(type(data[0]==list)):
            average = [0.0 for x in range(len(data[0]))]
            stdev = [0.0 for x in range(len(data[0]))]
            #print(data)
            for ii in range(len(data[0])):
                average[ii]=sum([ll[ii]/len(data) for ll in data])
            variance = [0.0 for x in range(len(data[0]))]
            for jj in range(len(data[0])):
                    variance[jj] = sum([(data[zz][jj]-average[jj])**2/len(data) for zz in range(len(data))])
            stdev = numpy.sqrt(variance)
        else:
            average = sum(data)/len(data)
            for ii in range(len(data)):
                variance = variance + (data[ii]-average)**2
            variance = variance/len(data)
            stdev = numpy.sqrt(variance)
        print ('Average of ',Pattern,' over ',len(data),' ares: mean = ',average,' variance = ', variance, ' stdev = ', stdev)
    else:
        print ('No input file or file ',FileName,'does not exist')
        sys.exit()


if __name__ == "__main__":
   main(sys.argv[1:])

