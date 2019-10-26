#!/usr/bin/env python
# -*- coding: utf-8 -*-

###
#AUTHORS: RAZVAN CARACAS, ANAIS KOBSCH
###

import sys,getopt,numpy,os.path,math
import crystallography as cr
import umd_process as umd

def print_gofrs(umdfile,MyCrystal,ndivx,discrete,normalization,maxlength,gofr):
    """Normalize gofr, compute integral, create the new file gofrs.dat and write data in it """
    newfile = umdfile[:-8] + '.gofr.dat' #name of the new file
    nf = open(newfile,'w')
    string = 'dist'+' '
    for ii in range(MyCrystal.ntypat):
        for jj in range(MyCrystal.ntypat):
            string = string + str(MyCrystal.elements[ii])+'-'+MyCrystal.elements[jj]+' Int('+MyCrystal.elements[ii]+'-'+MyCrystal.elements[jj]+') '
    string = string + '\n'
    nf.write(string)
    intgofr = [[0.0 for jatom in range(MyCrystal.ntypat)] for iatom in range (MyCrystal.ntypat)]
    for kk in range(ndivx+1):
        #print ('step no. ',kk)
        smallr = float(kk)*discrete
        bigr = smallr + discrete
        volshell = 4.0 * numpy.pi * (bigr**3 - smallr**3) / 3.0
        volcell = maxlength**3 #in case of non cubic cell, we'll have to use MyCrystal.cellvolume, and so import the rprimd in the read_umd
        string = str(smallr+discrete/2)+' ' 
        for iatom in range(MyCrystal.ntypat):
            for jatom in range(MyCrystal.ntypat):
                gofr[iatom][jatom][0] = 0 #we delete this element since it correspond to the distance = 0 (iatom with itself) 
                gofr[iatom][jatom][kk] = float(gofr[iatom][jatom][kk]) * volcell / (volshell * MyCrystal.types[iatom] * MyCrystal.types[jatom] * normalization )
                intgofr[iatom][jatom] += gofr[iatom][jatom][kk]*smallr**2*discrete*4*numpy.pi * MyCrystal.types[jatom]/volcell #cumulative integral of gofr in spherical coordinates 
                string = string + str(gofr[iatom][jatom][kk]) + ' ' + str(intgofr[iatom][jatom]) + ' '
        string = string + '\n'
        nf.write(string)
    nf.close()


def computeallgofr(MyCrystal,MySnapshot,discrete,maxlength,gofr):
    """For every couple of atoms, compute the minimal distance betweens atoms and add 1 to the corresponding gofr[r,r+delta_r] slot """
    for iatom in range (MyCrystal.natom):
        for jatom in range(MyCrystal.natom):
            #if you change this line to for jatom in range(iatom+1,MyCrystl.natom) you speed up everything by two. 
            #just don't forget about the double counting
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
            if (d1 < maxlength/2): #if the distance between the two atoms is below half of the minimum distance of the unit cell, then we use it to compute the gofr
                gofr[MyCrystal.typat[iatom]][MyCrystal.typat[jatom]][int(d1/discrete)] += 1
    return(gofr)


def read_umd(umdfile,Nsteps,discrete,InitialStep):
    """Read umd file and store data into classes """
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
                    print('natom=',MyCrystal.natom)
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
                    maxlength = min(MyCrystal.acell[0],MyCrystal.acell[1],MyCrystal.acell[2])
    ndivx = int(maxlength / (2*discrete))
    print ('acell',MyCrystal.acell[0],'discrete',discrete,'ndivx ',ndivx)
    gofr = [[[0 for kk in range(ndivx+1)] for jatom in range(MyCrystal.ntypat)] for iatom in range(MyCrystal.ntypat )]
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
                        if niter % Nsteps == 0: #if the remainder of niter/Nsteps is 0, then I compute the gofr
                            #print (' at step ',istep,' I am calling gofr ')
                            gofr = computeallgofr(MyCrystal,MySnapshot,discrete,maxlength,gofr)
                        niter += 1
    normalization = niter / Nsteps  #number of time steps actually used in the calculation of gofr
    print_gofrs(umdfile,MyCrystal,ndivx,discrete,normalization,maxlength,gofr)
        

def main(argv):
    umd.headerumd()
    umdfile='output.umd.dat'
    Nsteps = 1
    discrete = 0.01            #delta_r  = width of bins in histogram
    InitialStep = 0            #in case we want to skip additional timesteps
    try:
        opts, arg = getopt.getopt(argv,"hf:s:d:i:",["fumdfile","Sampling_Frequency","ddiscrete","iInitialStep"])
    except getopt.GetoptError:
        print ('gofrs_umd.py -f <umdfile>  -s <No.steps>  -d <discretization.interval> -i <InitialStep>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('gofrs_umd.py program to compute pair-distribution fuctions g(r) and print all the results in one file. Usage:')
            print ('gofrs_umd.py -f <umdfile> -s <Sampling_Frequency> -d <discretization.interval> -i <InitialStep>')
            print ('umdfile = name of the umd file. Default = output.umd.dat')
            print ('Sampling_Frequency = frequency of sampling the trajectory. Default = 1')
            print ('discretization.interval = for plotting the g(r). Default = 0.01 Angstroms')
            print ('InitialStep = number of initial steps of the trajectory to be removed. Default = 0')
            sys.exit()
        elif opt in ("-f", "--fumdfile"):
            umdfile = str(arg)
            #print('umdfile = ',umdfile)
        elif opt in ("-s", "--sNsteps"):
            Nsteps = int(arg)
        elif opt in ("-d", "--ddiscrete"):
            discrete = float(arg)
            #print('the distance delta_r we use to compute the number of atoms around the central one is ',discrete, 'Angstroms')
        elif opt in ("-i", "--iInitialStep"):
            InitialStep = int(arg)
                #print('I will skip  ',initial,' timesteps. ')
    if (os.path.isfile(umdfile)): 
        read_umd(umdfile,Nsteps,discrete,InitialStep)
    else:
        print ('the umd file ',umdfile,' does not exist')
        sys.exit()


if __name__ == "__main__":
   main(sys.argv[1:])
