#!/usr/bin/env python3

import numpy as np
import crystallography as cr


##
## UMD READING FUNCTIONS
##
## one function reads the entire umd
## then a series of functions can read several parts of the umd (like the main header, the info of each snapshot, the atomic positions, the atomic velocities, the forces, etc
##


def readumd(umdfile, short=0):
    niter = 0
    MyCrystal = cr.Lattice()
    AllSnapshots = []
    MySnapshot = cr.Lattice()
    with open(umdfile,'r') as ff:
        while True:
            line = ff.readline()
            if not line: break
            #print(line,len(line))
            if len(line) > 1:
                line=line.strip()
                entry=line.split()
                if entry[0] == 'natom':
                    MyCrystal.natom = int(entry[1])
                    MyCrystal.typat = [0 for _ in range(MyCrystal.natom)]
                if entry[0] == 'ntypat':
                    MyCrystal.ntypat = int(entry[1])
                    MyCrystal.types = [0 for _ in range(MyCrystal.ntypat)]
                    MyCrystal.elements = ['X' for _ in range(MyCrystal.ntypat)]
                    MyCrystal.masses = [1.0 for _ in range(MyCrystal.ntypat)]
                    MyCrystal.zelec = [1.0 for _ in range(MyCrystal.ntypat)]
                if entry[0] == 'types':
                    for ii in range(MyCrystal.ntypat):
                        MyCrystal.types[ii]=int(entry[ii+1])
                if entry[0] == 'elements':
                    for ii in range(MyCrystal.ntypat):
                        MyCrystal.elements[ii]=entry[ii+1]
                if entry[0] == 'masses':
                    for ii in range(MyCrystal.ntypat):
                        MyCrystal.masses[ii]=float(entry[ii+1])
                if entry[0] == 'Zelectrons':
                    for ii in range(MyCrystal.ntypat):
                        MyCrystal.zelec[ii]=float(entry[ii+1])
                if entry[0] == 'typat':
                    for ii in range(MyCrystal.natom):
                        MyCrystal.typat[ii] = int(entry[ii+1])
                if entry[0] == 'timestep':
                    TimeStep = float(entry[1])
                if entry[0] == 'acell':
                    for ii in range(3):
                        MySnapshot.acell[ii] = float(entry[ii+1])
                if entry[0] == 'acell':
                    for ii in range(3):
                        MySnapshot.acell[ii] = float(entry[ii+1])
                if entry[0] == 'rprim_a':
                    for ii in range(3):
                        MySnapshot.rprim[0][ii] = float(entry[ii+1])
                if entry[0] == 'rprim_b':
                    for ii in range(3):
                        MySnapshot.rprim[1][ii] = float(entry[ii+1])
                if entry[0] == 'rprim_c':
                    for ii in range(3):
                        MySnapshot.rprim[2][ii] = float(entry[ii+1])
                if entry[0] == 'rprimd_a':
                    for ii in range(3):
                        MySnapshot.rprimd[0][ii] = float(entry[ii+1])
                if entry[0] == 'rprimd_b':
                    for ii in range(3):
                        MySnapshot.rprimd[1][ii] = float(entry[ii+1])
                if entry[0] == 'rprimd_c':
                    for ii in range(3):
                        MySnapshot.rprimd[2][ii] = float(entry[ii+1])

                if entry[0] == 'atoms:':
                    #print('current iteration no.',niter)
                    MySnapshot.atoms = [cr.Atom() for _ in range(MyCrystal.natom)]
                    for iatom in range(MyCrystal.natom):
                        line = ff.readline()
                        line=line.strip()
                        entry=line.split()
                        for jj in range(3):
                            MySnapshot.atoms[iatom].xred[jj] = float(entry[jj])
                            MySnapshot.atoms[iatom].xcart[jj] = float(entry[jj+3])
                            MySnapshot.atoms[iatom].absxcart[jj] = float(entry[jj+6])
                            MySnapshot.atoms[iatom].vels[jj] = float(entry[jj+9])
                    #print(MySnapshot.atoms[iatom].xcart[0])
                    AllSnapshots.append(MySnapshot)
                    if short == 1:
                        break  #if we selected the short option then we will stop to read the umd file after the first iteration
                    MySnapshot = cr.Lattice()
                    niter += 1
    #!!! remove first element from AllSnapshots
    print('len of allsnapshots is',len(AllSnapshots))
    return(MyCrystal,AllSnapshots,TimeStep)


def read_bigheader_umd(umdfile, short=0):
    MyCrystal = cr.Lattice()
    with open(umdfile,'r') as ff:
        while True:
            line = ff.readline()
            if not line: break
            #print(line,len(line))
            if len(line) > 1:
                line=line.strip()
                entry=line.split()
                if entry[0] == 'natom':
                    MyCrystal.natom = int(entry[1])
                    MyCrystal.typat = [0 for _ in range(MyCrystal.natom)]
                if entry[0] == 'ntypat':
                    MyCrystal.ntypat = int(entry[1])
                    MyCrystal.types = [0 for _ in range(MyCrystal.ntypat)]
                    MyCrystal.elements = ['X' for _ in range(MyCrystal.ntypat)]
                    MyCrystal.masses = [1.0 for _ in range(MyCrystal.ntypat)]
                    MyCrystal.zelec = [1.0 for _ in range(MyCrystal.ntypat)]
                if entry[0] == 'types':
                    for ii in range(MyCrystal.ntypat):
                        MyCrystal.types[ii]=int(entry[ii+1])
                if entry[0] == 'elements':
                    for ii in range(MyCrystal.ntypat):
                        MyCrystal.elements[ii]=entry[ii+1]
                if entry[0] == 'masses':
                    for ii in range(MyCrystal.ntypat):
                        MyCrystal.masses[ii]=float(entry[ii+1])
                if entry[0] == 'Zelectrons':
                    for ii in range(MyCrystal.ntypat):
                        MyCrystal.zelec[ii]=float(entry[ii+1])
                if entry[0] == 'typat':
                    for ii in range(MyCrystal.natom):
                        MyCrystal.typat[ii] = int(entry[ii+1])
                if entry[0] == 'timestep':
                    TimeStep = float(entry[1])
                    break
    return(MyCrystal,TimeStep)


def read_header_snapshots(umdfile, short=0):#Ne semble pas être utilisée
    AllSnapshots = []
    MyCrystal=cr.Lattice()
    MySnapshot = cr.Lattice()
    with open(umdfile,'r') as ff:
        while True:
            line = ff.readline()
            if not line: break
            #print(line,len(line))
            if len(line) > 1:
                line=line.strip()
                entry=line.split()
                if entry[0] == 'natom':
                    MyCrystal.natom = int(entry[1])
                    MyCrystal.typat = [0 for _ in range(MyCrystal.natom)]
                if entry[0] == 'ntypat':
                    MyCrystal.ntypat = int(entry[1])
                    MyCrystal.types = [0 for _ in range(MyCrystal.ntypat)]
                    MyCrystal.elements = ['X' for _ in range(MyCrystal.ntypat)]
                    MyCrystal.masses = [1.0 for _ in range(MyCrystal.ntypat)]
                    MyCrystal.zelec = [1.0 for _ in range(MyCrystal.ntypat)]
                if entry[0] == 'types':
                    for ii in range(MyCrystal.ntypat):
                        MyCrystal.types[ii]=int(entry[ii+1])
                if entry[0] == 'elements':
                    for ii in range(MyCrystal.ntypat):
                        MyCrystal.elements[ii]=entry[ii+1]
                if entry[0] == 'masses':
                    for ii in range(MyCrystal.ntypat):
                        MyCrystal.masses[ii]=float(entry[ii+1])
                if entry[0] == 'Zelectrons':
                    for ii in range(MyCrystal.ntypat):
                        MyCrystal.zelec[ii]=float(entry[ii+1])
                if entry[0] == 'typat':
                    for ii in range(MyCrystal.natom):
                        MyCrystal.typat[ii] = int(entry[ii+1])
                if entry[0] == 'timestep':
                    TimeStep = float(entry[1])
                if entry[0] == 'acell':
                    for ii in range(3):
                        MySnapshot.acell[ii] = float(entry[ii+1])
                if entry[0] == 'acell':
                    for ii in range(3):
                        MySnapshot.acell[ii] = float(entry[ii+1])
                if entry[0] == 'rprim_a':
                    for ii in range(3):
                        MySnapshot.rprim[0][ii] = float(entry[ii+1])
                if entry[0] == 'rprim_b':
                    for ii in range(3):
                        MySnapshot.rprim[1][ii] = float(entry[ii+1])
                if entry[0] == 'rprim_c':
                    for ii in range(3):
                        MySnapshot.rprim[2][ii] = float(entry[ii+1])
                if entry[0] == 'rprimd_a':
                    for ii in range(3):
                        MySnapshot.rprimd[0][ii] = float(entry[ii+1])
                if entry[0] == 'rprimd_b':
                    for ii in range(3):
                        MySnapshot.rprimd[1][ii] = float(entry[ii+1])
                if entry[0] == 'rprimd_c':
                    for ii in range(3):
                        MySnapshot.rprimd[2][ii] = float(entry[ii+1])

                if entry[0] == 'atoms:':
                    AllSnapshots.append(MySnapshot)
    return(MyCrystal,AllSnapshots,TimeStep)


def read_absxcart(umdfile, short=0):
    niter = 0
    MyCrystal = cr.Lattice()
    AllSnapshots = []
    MySnapshot = cr.Lattice()
    with open(umdfile,'r') as ff:
        while True:
            line = ff.readline()
            if not line: break
            #print(line,len(line))
            if len(line) > 1:
                line=line.strip()
                entry=line.split()
                if entry[0] == 'natom':
                    MyCrystal.natom = int(entry[1])
                    MyCrystal.typat = [0 for _ in range(MyCrystal.natom)]
                if entry[0] == 'ntypat':
                    MyCrystal.ntypat = int(entry[1])
                    MyCrystal.types = [0 for _ in range(MyCrystal.ntypat)]
                    MyCrystal.elements = ['X' for _ in range(MyCrystal.ntypat)]
                    MyCrystal.masses = [1.0 for _ in range(MyCrystal.ntypat)]
                    MyCrystal.zelec = [1.0 for _ in range(MyCrystal.ntypat)]
                if entry[0] == 'types':
                    for ii in range(MyCrystal.ntypat):
                        MyCrystal.types[ii]=int(entry[ii+1])
                if entry[0] == 'elements':
                    for ii in range(MyCrystal.ntypat):
                        MyCrystal.elements[ii]=entry[ii+1]
                if entry[0] == 'masses':
                    for ii in range(MyCrystal.ntypat):
                        MyCrystal.masses[ii]=float(entry[ii+1])
                if entry[0] == 'Zelectrons':
                    for ii in range(MyCrystal.ntypat):
                        MyCrystal.zelec[ii]=float(entry[ii+1])
                if entry[0] == 'typat':
                    for ii in range(MyCrystal.natom):
                        MyCrystal.typat[ii] = int(entry[ii+1])
                if entry[0] == 'timestep':
                    TimeStep = float(entry[1])
                if entry[0] == 'acell':
                    for ii in range(3):
                        MySnapshot.acell[ii] = float(entry[ii+1])
                if entry[0] == 'acell':
                    for ii in range(3):
                        MySnapshot.acell[ii] = float(entry[ii+1])
                if entry[0] == 'rprim_a':
                    for ii in range(3):
                        MySnapshot.rprim[0][ii] = float(entry[ii+1])
                if entry[0] == 'rprim_b':
                    for ii in range(3):
                        MySnapshot.rprim[1][ii] = float(entry[ii+1])
                if entry[0] == 'rprim_c':
                    for ii in range(3):
                        MySnapshot.rprim[2][ii] = float(entry[ii+1])
                if entry[0] == 'rprimd_a':
                    for ii in range(3):
                        MySnapshot.rprimd[0][ii] = float(entry[ii+1])
                if entry[0] == 'rprimd_b':
                    for ii in range(3):
                        MySnapshot.rprimd[1][ii] = float(entry[ii+1])
                if entry[0] == 'rprimd_c':
                    for ii in range(3):
                        MySnapshot.rprimd[2][ii] = float(entry[ii+1])

                if entry[0] == 'atoms:':
                    #print('current iteration no.',niter)
                    MySnapshot.atoms = [cr.Atom() for _ in range(MyCrystal.natom)]
                    for iatom in range(MyCrystal.natom):
                        line = ff.readline()
                        line=line.strip()
                        entry=line.split()
                        for jj in range(3):
                            MySnapshot.atoms[iatom].absxcart[jj] = float(entry[jj+6])
                    #print(MySnapshot.atoms[iatom].xcart[0])
                    AllSnapshots.append(MySnapshot)
                    if short == 1:
                        break  #if we selected the short option then we will stop to read the umd file after the first iteration
                    MySnapshot = cr.Lattice()
                    niter += 1
    #!!! remove first element from AllSnapshots
    print('len of allsnapshots is',len(AllSnapshots))
    return(MyCrystal,AllSnapshots,TimeStep)


def read_xcart(umdfile, short=0):
    niter = 0
    MyCrystal = cr.Lattice()
    AllSnapshots = []
    MySnapshot = cr.Lattice()
    with open(umdfile,'r') as ff:
        while True:
            line = ff.readline()
            if not line: break
            #print(line,len(line))
            if len(line) > 1:
                line=line.strip()
                entry=line.split()
                if entry[0] == 'natom':
                    MyCrystal.natom = int(entry[1])
                    MyCrystal.typat = [0 for _ in range(MyCrystal.natom)]
                if entry[0] == 'ntypat':
                    MyCrystal.ntypat = int(entry[1])
                    MyCrystal.types = [0 for _ in range(MyCrystal.ntypat)]
                    MyCrystal.elements = ['X' for _ in range(MyCrystal.ntypat)]
                    MyCrystal.masses = [1.0 for _ in range(MyCrystal.ntypat)]
                    MyCrystal.zelec = [1.0 for _ in range(MyCrystal.ntypat)]
                if entry[0] == 'types':
                    for ii in range(MyCrystal.ntypat):
                        MyCrystal.types[ii]=int(entry[ii+1])
                if entry[0] == 'elements':
                    for ii in range(MyCrystal.ntypat):
                        MyCrystal.elements[ii]=entry[ii+1]
                if entry[0] == 'masses':
                    for ii in range(MyCrystal.ntypat):
                        MyCrystal.masses[ii]=float(entry[ii+1])
                if entry[0] == 'Zelectrons':
                    for ii in range(MyCrystal.ntypat):
                        MyCrystal.zelec[ii]=float(entry[ii+1])
                if entry[0] == 'typat':
                    for ii in range(MyCrystal.natom):
                        MyCrystal.typat[ii] = int(entry[ii+1])
                if entry[0] == 'timestep':
                    TimeStep = float(entry[1])
                if entry[0] == 'acell':
                    for ii in range(3):
                        MySnapshot.acell[ii] = float(entry[ii+1])
                if entry[0] == 'acell':
                    for ii in range(3):
                        MySnapshot.acell[ii] = float(entry[ii+1])
                if entry[0] == 'rprim_a':
                    for ii in range(3):
                        MySnapshot.rprim[0][ii] = float(entry[ii+1])
                if entry[0] == 'rprim_b':
                    for ii in range(3):
                        MySnapshot.rprim[1][ii] = float(entry[ii+1])
                if entry[0] == 'rprim_c':
                    for ii in range(3):
                        MySnapshot.rprim[2][ii] = float(entry[ii+1])
                if entry[0] == 'rprimd_a':
                    for ii in range(3):
                        MySnapshot.rprimd[0][ii] = float(entry[ii+1])
                if entry[0] == 'rprimd_b':
                    for ii in range(3):
                        MySnapshot.rprimd[1][ii] = float(entry[ii+1])
                if entry[0] == 'rprimd_c':
                    for ii in range(3):
                        MySnapshot.rprimd[2][ii] = float(entry[ii+1])

                if entry[0] == 'atoms:':
                    #print('current iteration no.',niter)
                    MySnapshot.atoms = [cr.Atom() for _ in range(MyCrystal.natom)]
                    for iatom in range(MyCrystal.natom):
                        line = ff.readline()
                        line=line.strip()
                        entry=line.split()
                        for jj in range(3):
                            MySnapshot.atoms[iatom].xcart[jj] = float(entry[jj+3])
                    #print(MySnapshot.atoms[iatom].xcart[0])
                    AllSnapshots.append(MySnapshot)
                    if short == 1:
                        break  #if we selected the short option then we will stop to read the umd file after the first iteration
                    MySnapshot = cr.Lattice()
                    niter += 1
    #!!! remove first element from AllSnapshots
    print('len of allsnapshots is',len(AllSnapshots))
    return(MyCrystal,AllSnapshots,TimeStep)


def read_stresses_4visc(umdfile, short=0):
    niter = 0
    AllSnapshots = []
    MySnapshot = cr.Lattice()
    with open(umdfile,'r') as ff:
        while True:
            line = ff.readline()
            if not line: break
            #print(line,len(line))
            if len(line) > 1:
                line=line.strip()
                entry=line.split()
                if entry[0] == 'timestep':
                    TimeStep = float(entry[1])
                if entry[0] == 'rprimd_a':
                    for ii in range(3):
                        MySnapshot.rprimd[0][ii] = float(entry[ii+1])
                if entry[0] == 'rprimd_b':
                    for ii in range(3):
                        MySnapshot.rprimd[1][ii] = float(entry[ii+1])
                if entry[0] == 'rprimd_c':
                    for ii in range(3):
                        MySnapshot.rprimd[2][ii] = float(entry[ii+1])
                if entry[0] == 'Temperature':
                    MySnapshot.temperature = float(entry[1])
                if entry[0] == 'StressTensor':
                    MySnapshot.stress = [0.0 for x in range(6)]
                    for ii in range(6):
                        MySnapshot.stress[ii] = float(entry[ii+1])

                if entry[0] == 'atoms:':
                    #print('current iteration no.',niter)
                    MySnapshot.cellvolume = MySnapshot.makevolume()
                    #print(MySnapshot.atoms[iatom].xcart[0])
                    AllSnapshots.append(MySnapshot)
                    if short == 1:
                        break  #if we selected the short option then we will stop to read the umd file after the first iteration
                    MySnapshot = cr.Lattice()
                    niter += 1
    #!!! remove first element from AllSnapshots
    print('len of allsnapshots is',len(AllSnapshots))
    return(AllSnapshots,TimeStep)



def print_header(FileName,MyCrystal):
    newfile = FileName + '.umd.dat'
    nf = open(newfile,'a')
    string = 'natom ' + str(MyCrystal.natom) + '\n'
    nf.write(string)
    string = 'ntypat ' + str(MyCrystal.ntypat) + '\n'
    nf.write(string)
    string = 'no.of.electrons ' + str(MyCrystal.noelectrons) + '\n'
    nf.write(string)
    string = 'types '
    for ii in range(MyCrystal.ntypat):
        string = string + str(MyCrystal.types[ii]) + ' '
    string = string + '\n'
    nf.write(string)
    string = 'elements '
    for ii in range(MyCrystal.ntypat):
        string = string + str(MyCrystal.elements[ii]) + ' '
    string = string + '\n'
    nf.write(string)
    string = 'masses '
    for ii in range(MyCrystal.ntypat):
        string = string + str(MyCrystal.masses[ii]) + ' '
    string = string + '\n'
    nf.write(string)
    string = 'Zelectrons '
    for ii in range(MyCrystal.ntypat):
        string = string + str(MyCrystal.zelec[ii]) + ' '
    string = string + '\n'
    nf.write(string)
    string = 'typat '
    for ii in range(MyCrystal.natom):
        string = string + str(MyCrystal.typat[ii]) + ' '
    string = string + '\n\n'
    nf.write(string)
    string = 'lambda_ThermoInt' + ' ' +  str(MyCrystal.lambda_ThermoInt)
    string = string + '\n\n'
    nf.write(string)
    nf.close()


def print_snapshots(FileName,MyCrystal,TimeStep,CurrentTime,diffcoords):
    newfile = FileName + '.umd.dat'
    nf = open(newfile,'a')
    string = 'timestep ' + str(TimeStep) + ' fs\n' + 'time ' + str(CurrentTime) + ' fs\n'
    nf.write(string)
    string = 'InternalEnergy ' + str(round(MyCrystal.internalenergy,6)) + ' eV\n'
    nf.write(string)
    string = 'ElectronicEntropy ' + str(round(MyCrystal.electronicentropy,6)) + ' eV\n'
    nf.write(string)
    string = 'KineticEnergyIons ' + str(round(MyCrystal.kineticenergy,6)) + ' eV\n'
    nf.write(string)
    string = 'EnergyWithDrift ' + str(round(MyCrystal.energywithdrift,6)) + ' eV\n'
    nf.write(string)
    string = 'Enthalpy ' + str(round(MyCrystal.enthalpy,6)) + ' eV\n'
    nf.write(string)
    string = 'Magnetization ' + str(round(MyCrystal.magnetization,6)) + ' Magneton-Bohr\n'
    nf.write(string)
    string = 'Temperature ' + str(round(MyCrystal.temperature,1)) + ' K\n'
    nf.write(string)
    string = 'Pressure ' + str(round(MyCrystal.pressure,4)) + ' GPa\n'
    nf.write(string)
    string = 'Density ' + str(round(MyCrystal.density,3)) + ' g.cm-3\n'
    nf.write(string)
    string = 'StressTensor '
    for ii in range(6):
        string = string + str(round(MyCrystal.stress[ii],4)) + ' '
    string = string + ' GPa\n'
    nf.write(string)
#    string = 'acell ' + str(round(MyCrystal.acell[0],3)) + ' ' + str(round(MyCrystal.acell[1],3) + ' ' + str(round(MyCrystal.acell[2],3)) + ' A\n'
    string = 'acell ' + str(MyCrystal.acell[0]) + ' ' + str(MyCrystal.acell[1]) + ' ' + str(MyCrystal.acell[2]) + ' A\n'
    nf.write(string)
    string = 'rprim_a ' + str(MyCrystal.rprimd[0][0]/MyCrystal.acell[0]) + '  ' +str(MyCrystal.rprimd[0][1]/MyCrystal.acell[0]) + '  ' +str(MyCrystal.rprimd[0][2]/MyCrystal.acell[0]) + '\n'
    nf.write(string)
    string = 'rprim_b ' + str(MyCrystal.rprimd[1][0]/MyCrystal.acell[1]) + '  ' +str(MyCrystal.rprimd[1][1]/MyCrystal.acell[1]) + '  ' +str(MyCrystal.rprimd[1][2]/MyCrystal.acell[1]) + '\n'
    nf.write(string)
    string = 'rprim_c ' + str(MyCrystal.rprimd[2][0]/MyCrystal.acell[2]) + '  ' +str(MyCrystal.rprimd[2][1]/MyCrystal.acell[2]) + '  ' +str(MyCrystal.rprimd[2][2]/MyCrystal.acell[2]) + '\n'
    nf.write(string)
    string = 'rprimd_a ' + str(MyCrystal.rprimd[0][0]) + '  ' +str(MyCrystal.rprimd[0][1]) + '  ' +str(MyCrystal.rprimd[0][2]) + ' A\n'
    nf.write(string)
    string = 'rprimd_b ' + str(MyCrystal.rprimd[1][0]) + '  ' +str(MyCrystal.rprimd[1][1]) + '  ' +str(MyCrystal.rprimd[1][2]) + ' A\n'
    nf.write(string)
    string = 'rprimd_c ' + str(MyCrystal.rprimd[2][0]) + '  ' +str(MyCrystal.rprimd[2][1]) + '  ' +str(MyCrystal.rprimd[2][2]) + ' A\n'
    nf.write(string)
    string = 'atoms: reduced*3 cartesian*3(A) abs.diff.*3(A) velocity*3(A/fs) force*3(eV/A) charge(no.elec) magnetization(magneton-Bohr) \n'
    nf.write(string)
    for iatom in range(MyCrystal.natom):
        string = str(round(MyCrystal.atoms[iatom].xred[0],5)) + ' ' + str(round(MyCrystal.atoms[iatom].xred[1],5)) + ' ' + str(round(MyCrystal.atoms[iatom].xred[2],5)) + ' '
        string = string + str(round(MyCrystal.atoms[iatom].xcart[0],5)) + ' ' + str(round(MyCrystal.atoms[iatom].xcart[1],5)) + ' ' + str(round(MyCrystal.atoms[iatom].xcart[2],5)) + ' '
        string = string + str(round(diffcoords[iatom][0],5)) + ' ' + str(round(diffcoords[iatom][1],5)) + ' ' + str(round(diffcoords[iatom][2],5)) + ' '
        string = string + str(round(MyCrystal.atoms[iatom].vels[0],5)) + ' ' + str(round(MyCrystal.atoms[iatom].vels[1],5)) + ' ' + str(round(MyCrystal.atoms[iatom].vels[2],5)) + ' '
        string = string + str(round(MyCrystal.atoms[iatom].forces[0],5)) + ' ' + str(round(MyCrystal.atoms[iatom].forces[1],5)) + ' ' + str(round(MyCrystal.atoms[iatom].forces[2],5)) + ' '
        string = string + str(MyCrystal.atoms[iatom].charge) + ' ' + str(MyCrystal.atoms[iatom].magnet) + '\n'
        nf.write(string)
    string='\n'
    nf.write(string)
    nf.close()
    return(CurrentTime,TimeStep)


def sort_umd(MyCrystal):
    iflag = -1
    #print ('There are ',MyCrystal.natom,'  atoms to order')
    notypatoms = 0
    AtomicOrdering = [-1 for _ in range(MyCrystal.natom)]
    for iatom in range(MyCrystal.natom):
        #print('Atom no. ', iatom,' with symbol ',MyCrystal.atoms[iatom].symbol,' and ordering ',AtomicOrdering[iatom])
        if AtomicOrdering[iatom] == -1:
            #print ('new atomic type')
            iflag = iflag + 1
            notypatoms = notypatoms + 1
            AtomicOrdering[iatom] = iflag
            MyCrystal.typat[iatom] = notypatoms
            #print ('Atomic ordering of iatom ',iatom,' is ',AtomicOrdering[iatom])
            for jatom in range(iatom+1,MyCrystal.natom):
                if MyCrystal.atoms[jatom].symbol == MyCrystal.atoms[iatom].symbol:
                    #print ('small atom no ',jatom,MyCrystal.atoms[jatom].symbol,' same atomic type as big atom ',iatom,MyCrystal.atoms[iatom].symbol)
                    MyCrystal.typat[jatom] = notypatoms
                    iflag = iflag + 1
                    AtomicOrdering[jatom] = iflag
                    #print ('Atomic ordering of jatom ',jatom,' is ',AtomicOrdering[jatom],MyCrystal.atoms[jatom].symbol)
    return(np.argsort(AtomicOrdering))


def headerumd():
    print ('\n * * * * * ')
    print ('    UMD package for analyzing Molecular Dynamics simulations.')
    print ('distributed under GNU-GNU General Public License 3')
    print ('please cite as: ')
    print ('    Razvan Caracas, Jean-Alexis Hernandez, Anais Kobsch, Zhi Li, Natalia Solomatova, and Francois Soubiran' )
    print ('    UMD: An open source package for analyzing Molecular Dynamics simulations' )
    print ('    Journal of Visualized Experiments, in press/media prep (2020)' )
    print (' ')



