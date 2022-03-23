#!/usr/bin/python
import numpy
import sys
import getopt
import re
import os
import subprocess


#while True :
#    try:
#        l = f.readline().split(None,10)
#        data = float(''.join(l[1:]))
#        print l[:1], data
#        total.append(data)
#    except ValueError:
#        average = sum(total)/len(total)
#        print 'Arithmetic average is: ', average
#        break

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
    index = 0
    average = 0 
    stdev = 0 
    variance = 0 
    try:
        opts, arg = getopt.getopt(argv,"hf:p:s:",["fFileName","pPattern","sSkipSteps"])
    except getopt.GetoptError:
        print 'average.py -f <FileName> -p <Pattern> -s <SkipSteps> '
        sys.exit(d)
    for opt, arg in opts:
        if opt == '-h':
            print 'average.py program to extract and average numerical values'
            print 'average.py -f <FileName> -p <Pattern> -s <SkipSteps>'
            sys.exit()
        elif opt in ("-f", "--fFileName"):
            FileName = str(arg)
        elif opt in ("-p", "--pPattern"):
            Pattern = str(arg)
            anchor=Pattern.split()[0]
#            print 'anchor is ', anchor
        elif opt in ("-s", "--sSkipSteps"):
            SkipSteps = int(arg)
    if (os.path.isfile(FileName)): 
#        print 'will read ',FileName,' file'
#        print 'will look for ',Pattern
        patterns=subprocess.check_output(['grep',Pattern,FileName])
#        print 'lenght of patterns is',len(patterns)
#        print 'patterns is',patterns
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
#        print 'Arithmetic average is: ', average
        for ii in range(len(data)):
            variance = variance + (data[ii]-average)**2
        variance = variance/len(data)
        stdev = numpy.sqrt(variance)
#        print 'Standard deviation is ',stdev
        print 'Averages over ',len(data),' ares: mean = ',average,' variance = ', variance, ' stdev = ', stdev
    else:
        print 'No input file or file ',FileName,'does not exist'
        sys.exit()


if __name__ == "__main__":
   main(sys.argv[1:])

