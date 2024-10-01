#!/usr/bin/env python3
"""
Created on Tue Jun  6 11:00:43 2023

@author: Razvan Caracas
"""

import sys,getopt,os
import crystallography as cr
import umd_processes_fast as umdpf

def print_xyz(outputfile,MyCrystal,AllSnapshots,UMDname):
    ff = open(outputfile,'w')
    print("number of snapshots is ",len(AllSnapshots))
    for istep in range(len(AllSnapshots)):
        string = str(MyCrystal.natom) + '\n' + UMDname +'\n'
        ff.write(string)
        for iatom in range(MyCrystal.natom):
            string = MyCrystal.elements[MyCrystal.typat[iatom]] + ' ' + str(AllSnapshots[istep][3*iatom]) + ' ' + str(AllSnapshots[istep][3*iatom+1]) + ' ' + str(AllSnapshots[istep][3*iatom+2]) +'\n'
            ff.write(string)
    ff.close()

def print_poscar(outputfile,MyCrystal,AllSnapshots,UMDname):
    for istep in range(len(AllSnapshots)):
        print(istep)
        poscarfile = outputfile[:-6] + '_' + str(istep) + '_.POSCAR'
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

def print_lammps(outputfile,MyCrystal,AllSnapshots,UMDname):
    print("Printing Lammps")
    string = UMDname[:-8]+'.lammps'
    # header=str(natom)+'\t atoms\n'
    # header+=str(int(typat[-1])+1)+'\t atom types\n'
    # header+="0.0\t"+Lx+"\txlo xhi\n"
    # header+="0.0\t"+Ly+"\tylo yhi\n"
    # header+="0.0\t"+Lz+"\tzlo zhi\n\n\n"
    # header+="Atoms\n\n"
    #
    # newstring+=line[3]+'\t'+line[4]+'\t'+line[5]+'\n'#Only the xcart (cartesian) ones

def print_deepmd(ooutputfile,MyCrystal,AllSnapshots,UMDname):
    print("Printing deepMD")


def main(argv):
    iterstep = 1
    firststep = 0
    lastsep = 1
    Timestep = 1.0
    ksnapshot = -1
    outoption = 1
    flagoutname = 0
    UMDname = 'output.umd.dat'
    umdpf.headerumd()
    try:
        opts, arg = getopt.getopt(argv,"hf:i:s:k:o:r:")
    except getopt.GetoptError:
        #print('input parameters ',arg,' have an error')
        print('Please use as: ')
        print ('umd2out.py -f <umdfile> -i <InitialStep> -s <Sampling_Frequency> -k <Particular_Snapshot> -o <outputfile> -r <output_option>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('umd2out.py program to extract various files from the umd file')
            print ('uumd2out.py -f <umdfile> -i <InitialStep> -s <Sampling_Frequency> -k <Particular_Snapshot> -r <output_option>')
            print (' default values: -f output.umd.dat -i 0 -s 1 -o UMDfile.xyz -r 1 -k -1')
            print('\n For output, the options are: ')
            print('   -r 1 = xyz file')
            print('   -r 2 = poscar file')
            print('   -r 3 = lammps file')
            print('   -r 4 = deepmd training file')
            sys.exit()
        elif opt in ("-f"):
            UMDname = str(arg)
        elif opt in ("-i"):
            firststep = int(arg)
        elif opt in ("-k"):
            ksnapshot = int(arg)
            firststep = ksnapshot
            laststep = ksnapshot
        elif opt in ("-s"):
            iterstep = int(arg)
        elif opt in ("-r"):
            outoption = int(arg)


    if (os.path.isfile(UMDname)):
        print('The first ',firststep, 'timesteps will be discarded')
        print('The POSCAR file contains every ',iterstep,' timesteps')
        MyCrystal = cr.Lattice()
        AllSnapshots = [cr.Lattice]

        if outoption == 1:
            if ksnapshot == -1:
                (MyCrystal, AllSnapshots, TimeStep, length) = umdpf.read_values(UMDname, "xcart", mode="line",Nsteps=iterstep)
                outputfile = UMDname[:-7] + 'xyz'
                print_xyz(outputfile, MyCrystal, AllSnapshots, UMDname)
            else:
                (MyCrystal, AllSnapshots, TimeStep, length) = umdpf.read_values(UMDname, "xcart", mode="line",Nsteps=1,firststep=ksnapshot,laststep=ksnapshot+1)
                outputfile = UMDname[:-7] + str(ksnapshot) + '.xyz'
                print_xyz(outputfile, MyCrystal, AllSnapshots, UMDname)
        elif outoption == 2:
            (MyCrystal, AllSnapshots, TimeStep, length) = umdpf.read_values(UMDname, "everything", mode="line",
                                                                            Nsteps=iterstep)
            outputfile = UMDname[:-7] + 'poscar'
            print("output file is ", outputfile)
            print_poscar(outputfile,MyCrystal,AllSnapshots,UMDname)
        elif outoption == 3:
            (MyCrystal, AllSnapshots, TimeStep, length) = umdpf.read_values(UMDname, "everything", mode="line",
                                                                            Nsteps=iterstep)
            outputfile = UMDname[:-7] + 'lammps'
            print_lammps(outputfile,MyCrystal,AllSnapshots,UMDname)
        elif outoption == 4:
            (MyCrystal, AllSnapshots, TimeStep, length) = umdpf.read_values(UMDname, "everything", mode="line",
                                                                            Nsteps=iterstep)
            outputfile = UMDname[:-7] + 'deepmd'
            print_deepmd(outputfile,MyCrystal,AllSnapshots,UMDname)

    else:
        print ('the umdfile ',UMDname,' does not exist')
        sys.exit()

if __name__ == "__main__":
   main(sys.argv[1:])
