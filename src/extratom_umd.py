#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### Author : Razvan Caracas

import sys, getopt, os, time, subprocess
import umd_processes_fast as umdpf
import crystallography as cr

def print_atom_line(FileName,OutFile,TimeStep,iatom):
    # read umdfile line by line and only print the acell, rprim and the iatom lines into the new umd file
    istep = 0
    fout = open(OutFile, 'a')
    with open(FileName, 'r') as infile:
        for line in infile:
            if line.startswith("acell"):
                string = 'timestep ' + str(TimeStep) + ' fs\n'
                fout.write(string)
                istep += 1
                string = 'time ' + str(istep * TimeStep) + ' fs\n'
                fout.write(string)
                fout.write(line)
                for _ in range(6):
                    line = next(infile)
                    fout.write(line)
            elif line.startswith("atoms:"):
                fout.write(line)
                # Move iatom lines forward
                for _ in range(iatom):
                    line = next(infile)  # Skip iatom lines
                atom_str = next(infile)  # The next line is the atom we want
                fout.write(atom_str)
                fout.write("\n")
    fout.close()


def main(argv):
    TimeStep = 1
    iatom = 0
    umdpf.headerumd()
    start=time.time()
    umdfile='my.umd.dat'
    try:
        opts, arg = getopt.getopt(argv,"hf:n:",["fumdfile","nAtom"])
    except getopt.GetoptError:
        print ('extratom_umd.py -f <UMD_filename> -n <Extracted Atom>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('Program to Extract one atom from the umd file. Usage: ')
            print ('extratom_umd.py -f <UMD_filename> -n <Extracted Atom> ')
            print ('UMD_filename = input file with the trajectory, in UMD format.')
            print ('-n : The index of the atom to be extracted. Default = 0th atom')
            sys.exit()
        elif opt in ("-f", "--fumdfile"):
            umdfile = str(arg)
            #print ('I will use the ',umdfile,' for input','\n')
            umdoutfile = f'{umdfile[:-8]}.{iatom}.umd.dat'
        elif opt in ("-n", "--nAtom"):
            iatom = int(arg)
            if isinstance(iatom, int):
                print("Extracting atom no. ", iatom)
            else:
                print("Error: IATOM must be a positive integer greater than 0")
                sys.exit()


    if (os.path.isfile(umdfile)):

        # read tghe header of the umd file
        MyCrystal = cr.Lattice()
        AllSnapshots = [cr.Lattice]
        (MyCrystal, TimeStep) = umdpf.read_bigheader_umd(umdfile)
        #print("natom is ",MyCrystal.natom)

        #create new file
        umdoutfile = f'{umdfile[:-8]}.{iatom}.umd.dat'
        #print atomic header of the new file
        umdpf.print_atomic_header(umdoutfile, MyCrystal, iatom)
        #pint the atom iotself
        print_atom_line(umdfile,umdoutfile,TimeStep,iatom)
        print("Atom ", iatom, " successfully extracted to file ", umdoutfile)

    else:
        print('the umdfile ', umdfile, ' does not exist')
        sys.exit()

    print("time :",time.time()-start)

if __name__ == "__main__":
   main(sys.argv[1:])
