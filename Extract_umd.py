#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 15:13:21 2023

@author: k
"""
import sys,getopt,itertools
import umd_processes_fast as umdpf


def main(argv):
    umdfile=''
    parameter=None
    try:
        opts, arg = getopt.getopt(argv,"hf:p:",["fUMDfile","pParameter"])
    except getopt.GetoptError:
        print ('Extract_umd.py -f <UMD_filename> -p <pParameter>')
        sys.exit(2)
    if opts == [] :
        opts = [('-h','')]
    for opt, arg in opts:
        if opt == "-h":
          print("Extract_umd.py to extract the list of values taken by a parameter snapshot by snapshot")
          print("-f : name of the umd file")
          print("-p : name of the parameter (Temperature, Pressure, InternalEnergy...)")
        elif opt == "-f":
            umdfile=str(arg)
        elif opt == "-p":
            parameter = str(arg)

    MyCrystal,TimeStep = umdpf.Crystallization(umdfile)

    Var=[]
    Times=[]
    tunit=None
    unit=None
    ff=open(umdfile,"r")
    while True :
        line=ff.readline()
        if not line :
            break
        line=line.strip().split()
        if len(line)>0 and line[0]=="time":
            Times.append(float(line[1]))
            if tunit==None:
                tunit=line[2]
        elif len(line)>0 and line[0]==parameter :#If the line corresponds to the parameter
            Var.append(float(line[-2]))#Value of the parameter
            if unit==None:
                unit=line[-1]#Unity of the parameter
            for _ in range(MyCrystal.natom):
                ff.readline()
    return [Var,Times,unit,tunit,TimeStep]
        

if __name__ == "__main__" :
    main(sys.argv[1:])