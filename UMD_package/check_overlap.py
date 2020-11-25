#!/usr/bin/env python3
###
#AUTHORS: ANAIS KOBSCH
###

import sys,getopt,os.path,math, re
import numpy as np
from . import crystallography as cr
from . import umd_process as umdp


def print_header(umdfile,MyCrystal,radius):
    """Create the new file overlap.dat and write data in it """
    newfile = umdfile[:-8] + '.overlap.dat' #name of the new file
    nf = open(newfile,'w')
    string = 'natom ' + str(MyCrystal.natom) + '\n'
    nf.write(string)
    string = 'ntypat ' + str(MyCrystal.ntypat) + '\n'
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
    string = 'typat '
    for ii in range(MyCrystal.natom):
        string = string + str(MyCrystal.typat[ii]) + ' '
    string = string + '\n\n'
    nf.write(string)
    nf.write('Radius in Angstroms\n')
    for elem in MyCrystal.elements:
        nf.write( elem + '\t' + str(round(radius[elem],2)) + '\n' )
    nf.write('\n')
    nf.close()
    return(newfile)


def print_overlap(newfile,overlap,dmin,istep,MyCrystal,radius,level):
    """Write the iteration number in the new file, along with dmin and overlap fraction in case of overlap """
    nf = open(newfile,'a')
    string = 'Iteration '+str(istep)+'\t'
    string2 = 'OK'
    for ii in range(MyCrystal.ntypat):
        for jj in range(ii,MyCrystal.ntypat):
            Vtot = 4/3 * np.pi * radius[MyCrystal.elements[ii]]**3 + 4/3 * np.pi * radius[MyCrystal.elements[jj]]**3
            if overlap[ii][jj] > level/100*Vtot:
                string2 = 'NOT-ok'
                break
        else:
            continue  # only executed if the inner loop did NOT break
        break  # only executed if the inner loop DID break
    nf.write(string+string2+'\n')
    if string2 == 'NOT-ok':
        string = 'dmin\t' +  '\t'.join(str(MyCrystal.elements[i])  for i in range(len(MyCrystal.elements))) +'\n'
        nf.write(string)
        for i in range(len(dmin)):
            string =  str(MyCrystal.elements[i])+ '\t' + '\t'.join(str(round(dmin[i][j],2))  for j in range(len(dmin[i]))) + '\n'
            nf.write(string)
        nf.write('overlap_volume_%\n')
        maximum = 0
        for ii in range(MyCrystal.ntypat):
            for jj in range(ii,MyCrystal.ntypat):
                Vtot = 4/3 * np.pi * radius[MyCrystal.elements[ii]]**3 + 4/3 * np.pi * radius[MyCrystal.elements[jj]]**3
                string = str(MyCrystal.elements[ii])+'-'+str(MyCrystal.elements[jj]) + '\t' + str(round(overlap[ii][jj]/Vtot * 100,2)) + '\n'
                nf.write(string)
                if overlap[ii][jj]/Vtot * 100 > maximum:
                    maximum = overlap[ii][jj]/Vtot * 100
                    maxelems = str(MyCrystal.elements[ii])+'-'+str(MyCrystal.elements[jj])
        print('iter', istep,'for ',maxelems, 'max overlap of ', round(maximum,2), '%')
    
    


def computeoverlap(dmin, radius, MyCrystal):
    """For every couple of atoms, compute the fraction of overlap between atoms"""
    #Initialization of the matrix overlap with same dimension as dmin
    overlap = []
    for i in range(MyCrystal.ntypat):
        overlap.append([])
        for j in range(MyCrystal.ntypat):
            overlap[i].append(float('nan'))
    #compute distances and fill the matrix  overlap
    for iatom in range(MyCrystal.ntypat):
        for jatom in range(iatom,MyCrystal.ntypat):
            #distance overlap
            #overlap[iatom][jatom] = dmin[iatom][jatom] / (radius[MyCrystal.elements[iatom]] + radius[MyCrystal.elements[jatom]] )
            #volume overlap
            overlap[iatom][jatom] = (np.pi * (radius[MyCrystal.elements[iatom]] + radius[MyCrystal.elements[jatom]] - dmin[iatom][jatom])**2 * (dmin[iatom][jatom]**2 + 2*dmin[iatom][jatom] * radius[MyCrystal.elements[iatom]] + 2*dmin[iatom][jatom] * radius[MyCrystal.elements[jatom]] - 3*radius[MyCrystal.elements[iatom]]**2  - 3*radius[MyCrystal.elements[jatom]]**2 + 6*radius[MyCrystal.elements[iatom]]*radius[MyCrystal.elements[jatom]] ) ) / (12*dmin[iatom][jatom] ) 
#    print('overlap = ',overlap)
    return overlap


def computealldmin(MyCrystal,MySnapshot):
    """For every couple of atoms, compute the minimal distance between atoms"""
    #Initialization of the matrix dmin
    dmin=[]
    for i in range(MyCrystal.ntypat):
        dmin.append([])
        for j in range(MyCrystal.ntypat):
            dmin[i].append(999)
    #compute distances and fill the matrix dmin
    for iatom in range (MyCrystal.natom):
        for jatom in range(iatom,MyCrystal.natom):
            #compute distances along x, y and z between the two atoms
            dx = MySnapshot.atoms[jatom].xcart[0] - MySnapshot.atoms[iatom].xcart[0]
            dy = MySnapshot.atoms[jatom].xcart[1] - MySnapshot.atoms[iatom].xcart[1]
            dz = MySnapshot.atoms[jatom].xcart[2] - MySnapshot.atoms[iatom].xcart[2]
            #if these distances are too large, then we correct them by translating one of the two atoms of one unit cell
            #if we want to work on a non cubic cell, then these corrections need to be changed
            #if (dx < -maxlength/2): dx = dx + maxlength  #these two lines are replaced by the one below: if |dx|>acell/2 then int(...) gives ± 1 
            #if (dx >  maxlength/2): dx = dx - maxlength
            dx = dx-MyCrystal.acell[0]*int(dx/(0.5*MyCrystal.acell[0]))
            dy = dy-MyCrystal.acell[0]*int(dy/(0.5*MyCrystal.acell[0]))
            dz = dz-MyCrystal.acell[0]*int(dz/(0.5*MyCrystal.acell[0]))
            d1 = math.sqrt(dx**2 + dy**2 + dz**2) 
            #update the matrix containing the minimal distance for each type of atomic couple
            if d1 < dmin[MyCrystal.typat[iatom]][MyCrystal.typat[jatom]]:
                if d1 == 0: continue #we skip the case in which we consider one atom with itself
                else:
                    dmin[MyCrystal.typat[iatom]][MyCrystal.typat[jatom]] = d1
#    print('dmin = ',dmin)           
    return(dmin)


def read_radius(radiusfile):
    """Read the file radius.inp, extract and convert the radius of atomic spheres in angstroms"""
    AngtoBohr = 1.8897259886 #1A in Bohr
    radius = {}
    with open(radiusfile, 'r') as r:
        while True:
            line = r.readline()
            if not line: break
            entry = line.split()
            if re.match('[Bb]',entry[2]): #if unit is in Bohr then we convert to angstroms
                radius[entry[0]] = float(entry[1]) / AngtoBohr
            else:
                radius[entry[0]] = float(entry[1])
    #print('radius in angstrom are', radius)
    return(radius)
                


def read_umd(umdfile,radius,Nsteps,InitialStep,level):
    """Read umd file and store data into classes """
    umdp.headerumd()
    niter = 0
    istep = 0
    MyCrystal = cr.Lattice()
    with open(umdfile,'r') as ff:   #read the head and half of the 1st iter in order to get acell
        while True:
            line = ff.readline()
            if not line: break
            #print(line,len(line))
            if len(line) > 1:
                line=line.strip()
                entry=line.split()
                if entry[0] == 'natom':
                    MyCrystal.natom = int(entry[1])
                    #print('natom=',MyCrystal.natom)
                    MyCrystal.typat = [0 for _ in range(MyCrystal.natom)]
                if entry[0] == 'ntypat':
                    MyCrystal.ntypat = int(entry[1])
                    MyCrystal.types = [0 for _ in range(MyCrystal.ntypat)]
                    MyCrystal.elements = ['X' for _ in range(MyCrystal.ntypat)]
                if entry[0] == 'types':
                    for ii in range(MyCrystal.ntypat):
                        MyCrystal.types[ii]=int(entry[ii+1])
                if entry[0] == 'elements':
                    for ii in range(MyCrystal.ntypat):
                        MyCrystal.elements[ii]=entry[ii+1]
                if entry[0] == 'typat':
                    for ii in range(MyCrystal.natom):
                        MyCrystal.typat[ii] = int(entry[ii+1])
                if entry[0] == 'acell':
                    for i in range(3):
                        MyCrystal.acell[i] = float(entry[i+1])
    newfile = print_header(umdfile,MyCrystal,radius) #create the new file
    with open(umdfile,'r') as ff:   #read the entire file to get the atomic positions
        while True:
            line = ff.readline()
            if not line: break
            #print(line,len(line))
            if len(line) > 1:
                line=line.strip()
                entry=line.split()          
                if entry[0] == 'atoms:':
                    istep += 1
                    #print('current iteration no.',istep)
                    MySnapshot = cr.Lattice()
                    MySnapshot.atoms = [cr.Atom() for _ in range(MyCrystal.natom)]
                    for iatom in range(MyCrystal.natom):
                        line = ff.readline()
                        line=line.strip()
                        entry=line.split()
                        for jj in range(3):
                            MySnapshot.atoms[iatom].xcart[jj]=float(entry[jj+3])
                        #print(MySnapshot.atoms[iatom].xcart[0])
                    if istep > InitialStep:
                        if niter % Nsteps == 0: #if the remainder of niter/Nsteps is 0, then I compute the dmin
                            #print (' at step ',istep,' I am calling dmin ')
                            dmin = computealldmin(MyCrystal,MySnapshot)
                            overlap = computeoverlap(dmin, radius, MyCrystal)
                            print_overlap(newfile,overlap,dmin,istep,MyCrystal,radius,level)
                        niter += 1
    
        

def main(argv):
    umdfile='output.umd.dat'
    radiusfile = 'radius.inp' #file containing 3 columns: the element symbol, the value of RCORE and the unit (Angstrom or Bohr) for each element in the simu
    Nsteps = 1
    InitialStep = 0            #in case we want to skip additional timesteps
    level = 12      #default value for the overlap percentage
    try:
        opts, arg = getopt.getopt(argv,"hf:r:s:i:l:",["fumdfile","radiusfile","sSampling_Frequency","iInitialStep","level_overlap"])
    except getopt.GetoptError:
        print ('check_overlap.py -f <UMD_filename> -r <radius_filename> -s <Sampling_Frequency> -i <InitialStep> -l <level_of_overlap_%>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('check_overlap.py program to compute the minimum distance between atomic couples and print it along with the volume overlap % in case of overlap of the atomic spheres')
            print ('check_overlap.py -f <UMD_filename> -r <radius_filename>  -s <Sampling_Frequency> -i <InitialStep> -l <acceptable_level_of_overlap_%>')
            print ('Default values are: umdfile = output.umd.dat, radiusfile = radius.inp , Sampling_Frequency = 1, InitialStep = 0, level = 12 (%) ')
            print ('Radius.inp file must contain for each atom, on one line the element symbol, the value of the pseudopotential/PAW spheres, and the unit (Angstrom or Bohr)')
            sys.exit()
        elif opt in ("-f", "--fumdfile"):
            umdfile = str(arg)
            #print('umdfile = ',umdfile)
        elif opt in ("-r","--radiusfile"):
            radiusfile = str(arg)
            #print('radius file = ',radiusfile)
        elif opt in ("-s", "--sNsteps"):
            Nsteps = int(arg)
            #print('we compute the distances every ',Nsteps,' steps')
        elif opt in ("-i", "--iInitialStep"):
            InitialStep = int(arg)
                #print('I will skip  ',initial,' timesteps. ')
        elif opt in ("-l","--level_overlap"):
            level = float(arg)
    if (os.path.isfile(radiusfile)): 
        radius = read_radius(radiusfile)
    else:
        print ('the radius.inp file ',radiusfile,' does not exist')
        sys.exit()
    if (os.path.isfile(umdfile)):
        print('******** for the umd file',umdfile)
        read_umd(umdfile,radius,Nsteps,InitialStep,level)
    else:
        print ('the umd file ',umdfile,' does not exist')
        sys.exit()


if __name__ == "__main__":
   main(sys.argv[1:])
