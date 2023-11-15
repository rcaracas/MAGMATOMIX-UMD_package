#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 11:00:43 2023

@author: Kevin Jiguet-Covex
"""
import umd_processes_fast as umdpf
import sys, getopt


def main(argv):
    umdpf.headerumd()
    UMDname = ''
    Nsteps = 1
    try:
        opts, arg = getopt.getopt(argv,"hf:s:",['fUMDfile'])
    except getopt.GetoptError:
        print ('umd_to_lammps.py -f <umdfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('umd_to_lammps.py to convert a umd file into a LAMMPS-type file')
            print ('umd_to_lammps -f umdFile -s SamplingFrequency')
            sys.exit()
        elif opt in ("-f"):
            UMDname = str(arg)
        elif opt in ("-s"):
            Nsteps = int(arg)
            
    fa = open(UMDname[:-8]+'.lammps','w')#creating the lammps file
    ff = open(UMDname,'r')
    #We read the header to extract some information about the system, such as the number of atoms of each type
    
    while True :
        line = ff.readline().strip().split()
        if len(line)>0:
            if line[0]=="natom" :
                natom = int(line[1])
            elif line[0]=="typat":
                typat = [line[i] for i in range(1,natom+1)]                
            elif line[0]=="rprimd_a":
                Lx = line[1]
            elif line[0]=="rprimd_b":
                Ly = line[2]
            elif line[0]=="rprimd_c":
                Lz = line[3]
            if line[0]=="atoms:":
                break
    
    headervar = UMDname[:-8]+'.lammps'
    header=str(natom)+'\t atoms\n'
    header+=str(int(typat[-1])+1)+'\t atom types\n'
    header+="0.0\t"+Lx+"\txlo xhi\n"
    header+="0.0\t"+Ly+"\tylo yhi\n"
    header+="0.0\t"+Lz+"\tzlo zhi\n\n\n"
    header+="Atoms\n\n"
    
    atomindex=1
    snapshotindex = 0
    ff.close()
    ff=open(UMDname,'r')
    #We fill the file with the atoms coordinates
    while True :
        line = ff.readline()
        if not line : 
            break
        l = line.strip().split()
        if len(l)>0:
            if l[0]=='time':
                Time = l[1]
            elif l[0]=='atoms:':
                if(int(snapshotindex/Nsteps)==snapshotindex/Nsteps):                
                    print('Converting snapshot from time ',Time)
                    fa.write(headervar+'\t time : '+Time+'\n\n')
                    fa.write(header)
                    for atom in range(natom):
                        newstring = '     '+str(atomindex)+'  '+str(int(typat[atom])+1)+'\t'
                        line = ff.readline().strip().split()
                        newstring+=line[3]+'\t'+line[4]+'\t'+line[5]+'\n'#Only the xcart (cartesian) ones
                        fa.write(newstring)
                        atomindex+=1
                    fa.write('\n')
                    atomindex=1
                snapshotindex+=1
                    
    
    fa.close()
    ff.close()   
    print('LAMMPS file successfully created under the name <'+UMDname[:-8]+'.lammps'+'>')
        

if __name__ == "__main__":
   main(sys.argv[1:])
