#!/usr/bin/env python3

###
##AUTHORS: RAZVAN CARACAS
###

import sys,getopt,os
import crystallography as cr
import umd_process as umdp

def print_poscar(MyCrystal,AllSnapshots,UMDname,firststep,laststep,iterstep):
    for istep in range(firststep,laststep,iterstep):
        poscarfile = UMDname[:-7] + '_' + str(istep) + '_.POSCAR'
        ff = open(poscarfile,'w')
        string = UMDname + ' snapshot ' + str(istep) + '\n'
        ff.write(string)
        string = '  1.0 \n'
        ff.write(string)
        string = '   ' + str(AllSnapshots[istep].rprimd[0][0]) + '  '  + str(AllSnapshots[istep].rprimd[0][1]) + '  '  + str(AllSnapshots[istep].rprimd[0][2]) + '\n'
        ff.write(string)
        string = '   ' + str(AllSnapshots[istep].rprimd[1][0]) + '  '  + str(AllSnapshots[istep].rprimd[1][1]) + '  '  + str(AllSnapshots[istep].rprimd[1][2]) + '\n'
        ff.write(string)
        string = '   ' + str(AllSnapshots[istep].rprimd[2][0]) + '  '  + str(AllSnapshots[istep].rprimd[2][1]) + '  '  + str(AllSnapshots[istep].rprimd[2][2]) + '\n'
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
            string = str(AllSnapshots[istep].atoms[iatom].xred[0]) + ' ' + str(AllSnapshots[istep].atoms[iatom].xred[1]) + ' ' + str(AllSnapshots[istep].atoms[iatom].xred[2]) +'\n'
            ff.write(string)
        string = '\n'
        ff.write(string)
        for iatom in range(MyCrystal.natom):
            string = str(AllSnapshots[istep].atoms[iatom].vels[0]) + ' ' + str(AllSnapshots[istep].atoms[iatom].vels[1]) + ' ' + str(AllSnapshots[istep].atoms[iatom].vels[2]) +'\n'
            ff.write(string)
        string = '\n'
        ff.write(string)
        ff.close()

def main(argv):
    iterstep = 1
    firststep = 0
    laststep = 10000000
    UMDname = 'output.umd.dat'
    umdp.headerumd()
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
        (MyCrystal,AllSnapshots,TimeStep)=umdp.readumd(UMDname)
        if laststep > len(AllSnapshots):
            laststep = len(AllSnapshots)
        print_poscar(MyCrystal,AllSnapshots,UMDname,firststep,laststep,iterstep)
    else:
        print ('the umdfile ',umdfile,' does not exist')
        sys.exit()

if __name__ == "__main__":
   main(sys.argv[1:])
