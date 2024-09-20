#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 10:07:50 2023

@author: k
"""


### Author : Kevin Jiguet-Covex


import sys,getopt,os.path
import umd_processes_fast as umdpf
import time
import numpy as np
import concurrent.futures
from functools import partial
import ctypes
from os.path import join
import math
import msd_umd
import matplotlib.pyplot as plt


def mean_MSD(msd,instants):
    Mean=0
    for x in range(1,len(msd)):
        Mean+=msd[x]/instants[x]
    Mean/=(len(msd)-1)
    print(Mean)
 
    return Mean

def main(argv):    
    nData = 20
    nCores=None
    try:
        opts, arg = getopt.getopt(argv,"hf:n:a:k:",["fumdfile","nnData","aAxes","knCores"])
    except getopt.GetoptError:
        print("MSDdistrib.py -f <umdfile> -n <nData> -a <Axes> -k <nCores>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print ('MSDdistrib to calculate the msd in a system along a distribution of axes betxeen two liminar axes')
            print ('-f : UMD file')
            print ('-n : number of axes in the distribution (default 20)')
            print ('-a : Explicit expression of the liminar axes under the form [xa,xb,xc,ya,yb,yc]')
            print ('-k : number of cores used in parallelization')
            sys.exit(2)
        elif opt == "-f":
            File = str(arg)
        elif opt in "-n":
            nData = int(arg)
        elif opt in "-a":
            Axes = eval(arg)
        elif opt in "-k":
            nCores=int(arg)

    vectors = []
    Thetas = []
    A1 = np.array(Axes[:3])
    A2 = np.array(Axes[3:])
    
    mat = []
    
    for i in range(nData):#Creation of the axes
        theta = math.pi/nData*i
        v = math.cos(theta)*A1+math.sin(theta)*A2
        vectors.append(v)
        mat+=list(v)
        Thetas.append(theta/math.pi)
        
    args=["-f",File,"-z","10","-m","elements","-a",str(mat)]

    if nCores!=None:
        args+=["-k",str(nCores)]
        
    [msdfile,MSD,Instants,Elements,Axes] = msd_umd.main(args)#Computation (and file creation) of the msd

    MSDsEl = [[] for _ in range(len(Elements))]
    for el in range(len(Elements)):#Calculation of the mean MSD for each axe
        coeff = 0
        MSDsEl
        for theta in range(nData):
            for t in range(len(Instants)-1):
                m=MSD[el][1+theta][1+t]
                t=Instants[1+t]
                coeff += m/t
            coeff/=len(Instants)-1
            MSDsEl[el].append(coeff)
                                        
    return Thetas, MSDsEl , Elements, msdfile
        
if __name__ == "__main__":
    main(sys.argv[1:])
