#!/usr/bin/env python3
"""
Created on Tue Jun  6 11:00:43 2023

@author: Razvan Caracas
"""

import sys,getopt,os
import crystallography as cr
import umd_processes_fast as umdpf

def print_poscar(MyCrystal,AllSnapshots,UMDname,firststep,laststep,iterstep):
    print(len(AllSnapshots))
    print(len(AllSnapshots[0]))
    for istep in range(firststep,laststep,iterstep):
        print(istep)
        poscarfile = UMDname[:-7] + '_' + str(istep) + '_.POSCAR'
        ff = open(poscarfile,'w')
        string = UMDname + ' snapshot ' + str(istep) + '\n'
        ff.write(string)
        string = '  1.0 \n'
        ff.write(string)
        string = '   ' + str(MyCrystal.rprimd[0][0]) + '  '  + str(MyCrystal.rprimd[0][1]) + '  '  + str(MyCrystal.rprimd[0][2]) + '\n'
        ff.write(string)
        string = '   ' + str(MyCrystal.rprimd[1][0]) + '  '  + str(MyCrystal.rprimd[1][1]) + '  '  + str(MyCrystal.rprimd[1][2]) + '\n'
        ff.write(string)
        string = '   ' + str(MyCrystal.rprimd[2][0]) + '  '  + str(MyCrystal.rprimd[2][1]) + '  '  + str(MyCrystal.rprimd[2][2]) + '\n'
        ff.write(string)
        string = '  '
        for itype in range(MyCrystal.ntypat):
            string = string + MyCrystal.elements[itype] + '  '
        string = string + '\n'
        ff.write(string)
        string = '  '
        for itype in range(MyCrystal.ntypat):
            string = string + str(MyCrystal.types[itype]) + '  '
        string = string + '\n'
        ff.write(string)
        string = 'Direct\n'
        ff.write(string)
        for iatom in range(MyCrystal.natom):
            string = str(AllSnapshots[istep][12*iatom]) + ' ' + str(AllSnapshots[istep][12*iatom+1]) + ' ' + str(AllSnapshots[istep][12*iatom+2]) +'\n'
            ff.write(string)
        string = '\n'
        ff.write(string)
        for iatom in range(MyCrystal.natom):
            string = str(AllSnapshots[istep][12*iatom+9]) + ' ' + str(AllSnapshots[istep][12*iatom+10]) + ' ' + str(AllSnapshots[istep][12*iatom+11]) +'\n'
            ff.write(string)
        string = '\n'
        ff.write(string)
        ff.close()

def main(argv):
    iterstep = 1
    firststep = 0
    laststep = 10000000
    UMDname = 'output.umd.dat'
    umdpf.headerumd()
    try:
        opts, arg = getopt.getopt(argv,"hf:i:l:s:")
    except getopt.GetoptError:
        print ('umd2poscar.py -f <umdfile> -i <InitialStep> -l <LastStep> -s <Sampling_Frequency>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('umd2poscar.py program to extract POSCAR snapshots from the umd file')
            print ('umd2poscar.py -f <umdfile> -i <InitialStep> -l <LastStep> -s <Sampling_Frequency>')
            print (' default values: -f output.umd.dat -i 0 -l 10000000 -s 1')
            sys.exit()
        elif opt in ("-f"):
            UMDname = str(arg)
        elif opt in ("-i"):
            firststep = int(arg)
        elif opt in ("-l"):
            laststep = int(arg)
        elif opt in ("-s"):
            iterstep = int(arg)
    if (os.path.isfile(UMDname)):
        print('The first ',firststep, 'timesteps will be discarded')
        print('The POSCAR file contains every ',iterstep,' timesteps')
        MyCrystal = cr.Lattice()
        AllSnapshots = [cr.Lattice]
        (MyCrystal,AllSnapshots,TimeStep,length)=umdpf.read_values(UMDname,"everything",mode="line",Nsteps=iterstep)
        if laststep > len(AllSnapshots):
            laststep = len(AllSnapshots)
        print_poscar(MyCrystal,AllSnapshots,UMDname,firststep,laststep,iterstep)
    else:
        print ('the umdfile ',UMDname,' does not exist')
        sys.exit()

if __name__ == "__main__":
   main(sys.argv[1:])
