#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 13:28:39 2022

@author: timbogels
"""

import numpy as np
import sys,getopt,os
import crystallography as cr
import umd_process as up
import random


#https://stackabuse.com/calculating-euclidean-distance-with-numpy/
#np.linalg.norm() # for dinstance matrix

def shuffler(MyCrystal,AllSnapshots,UMDname,istep,nbox,close,min_distance):
    box_rib = 1/nbox
    allboxcoords = []
    box_id=[]
    displ_matrix=[]
    xval=yval=zval=0

    #Define the box coordinates
    for p in range(0,nbox):
        z = p * box_rib
        for q in range(0,nbox):
            y = q * box_rib
            for r in range(0,nbox):
                x = r * box_rib
                boxcoord = [x,y,z]
                allboxcoords.append(boxcoord)
    #print(allboxcoords)
    #print()

    #Make a list of the box ID that can be shuffled, e.g. index 0 becomes 15 means box 1 becomes box 16.
    for i in range(0,nbox**3):
        box_id.append(i)

    random.shuffle(box_id)
    #print(box_id)
    #print()

    #Build the displacement matrix that contains the displacement needed for all atoms in each box ID, e.g. the values on index 0 show what displacement all atoms in box 1 go through to go to their new location
    for i in range (0,nbox**3):
        new_id = box_id[i]
        x = ((allboxcoords[new_id])[0])-((allboxcoords[i])[0])
        y = ((allboxcoords[new_id])[1])-((allboxcoords[i])[1])
        z = ((allboxcoords[new_id])[2])-((allboxcoords[i])[2])
        displ_vector=[x,y,z]
        displ_matrix.append(displ_vector)

    #print(displ_matrix)

    for iatom in range(MyCrystal.natom):
        for p in range(0,nbox):
            if 0+p*box_rib <= AllSnapshots[istep].atoms[iatom].xred[2] < 0+(p+1)*box_rib:
                zval=p
        for q in range(0,nbox):
            if 0+q*box_rib <= AllSnapshots[istep].atoms[iatom].xred[1] < 0+(q+1)*box_rib:
                yval=q
        for r in range(0,nbox):
            if 0+r*box_rib <= AllSnapshots[istep].atoms[iatom].xred[0] < 0+(r+1)*box_rib:
                xval=r
        box_index = xval*nbox**0+yval*nbox**1+zval*nbox**2
        #print(xval,yval,zval,box_index,box_id[box_index])
        #print(AllSnapshots[istep].atoms[iatom].xred[0],AllSnapshots[istep].atoms[iatom].xred[1],AllSnapshots[istep].atoms[iatom].xred[2])
        AllSnapshots[istep].atoms[iatom].xred[2] += displ_matrix[box_index][2]
        AllSnapshots[istep].atoms[iatom].xred[1] += displ_matrix[box_index][1]
        AllSnapshots[istep].atoms[iatom].xred[0] += displ_matrix[box_index][0]
        #print(AllSnapshots[istep].atoms[iatom].xred[0],AllSnapshots[istep].atoms[iatom].xred[1],AllSnapshots[istep].atoms[iatom].xred[2])
    count = 0
    for i in range(MyCrystal.natom):
        for j in range(i,MyCrystal.natom):
            atom_1=np.array((AllSnapshots[istep].atoms[i].xred[0]*AllSnapshots[istep].rprimd[0][0],AllSnapshots[istep].atoms[i].xred[1]*AllSnapshots[istep].rprimd[0][0],AllSnapshots[istep].atoms[i].xred[2]*AllSnapshots[istep].rprimd[0][0]))
            atom_2=np.array((AllSnapshots[istep].atoms[j].xred[0]*AllSnapshots[istep].rprimd[0][0],AllSnapshots[istep].atoms[j].xred[1]*AllSnapshots[istep].rprimd[0][0],AllSnapshots[istep].atoms[j].xred[2]*AllSnapshots[istep].rprimd[0][0]))
            distance = np.linalg.norm(atom_1-atom_2)
            if distance < min_distance and i != j:
                print('Atoms are too close:', distance, 'atoms', i,j)
                count +=1
    if count >= 1:
        close = True
    else:
        close = False
                 
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
        string = '%.7s' % str(AllSnapshots[istep].atoms[iatom].xred[0]) + ' ' + '%.7s' % str(AllSnapshots[istep].atoms[iatom].xred[1]) + ' ' + '%.7s' % str(AllSnapshots[istep].atoms[iatom].xred[2]) +'\n'
        ff.write(string)
    string = '\n'
    ff.write(string)
    for iatom in range(MyCrystal.natom):
        string = str(AllSnapshots[istep].atoms[iatom].vels[0]) + ' ' + str(AllSnapshots[istep].atoms[iatom].vels[1]) + ' ' + str(AllSnapshots[istep].atoms[iatom].vels[2]) +'\n'
        ff.write(string)
    string = '\n'
    ff.write(string)
    ff.close()
    return(close)
    
def main(argv):
    istep = 0
    close = False
    nbox = 3
    lbox = 0
    min_distance = 0.5
    UMDname = 'output.umd.dat'
    up.headerumd()
    try:
        opts, arg = getopt.getopt(argv,"hf:i:b:l:d:")
    except getopt.GetoptError:
        print ('umd2poscar.py -f <umdfile> -i <Snapshot> -n <Box_gridsize> -l <Box_length> -d <Atom_distance_allowance>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('umd2poscar.py program to extract POSCAR snapshots from the umd file')
            print ('umd2poscar.py -f <umdfile> -i <Snapshot> -b <Box_gridsize> -l <Box_length> -d <Atom_distance_allowance>')
            print (' default values: -f output.umd.dat -i 0 -b 3 -l 0 (not used) -d 0.5')
            sys.exit()
        elif opt in ("-f"):
            UMDname = str(arg)
        elif opt in ("-i"):
            istep = int(arg)
        elif opt in ("-b"):
            nbox = int(arg)
        elif opt in ("-l"):
            lbox = float(arg)
        elif opt in ("-d"):
            min_distance = float(arg)
    if (os.path.isfile(UMDname)):
        print('The POSCAR file contains ',istep,' timestep, shuffled')
        MyCrystal = cr.Lattice()
        AllSnapshots = [cr.Lattice]
        (MyCrystal,AllSnapshots,TimeStep)=up.readumd(UMDname)
        if lbox > 0:
            nbox = int(AllSnapshots[istep].rprimd[0][0] // lbox)
            print('The nearest Box_gridsize for your given length = ',nbox)
        if nbox == 1:
            print('There is not much to shuffle with 1 box,exiting')
            sys.exit()
        print('number of boxes is ',nbox**3)
        close = shuffler(MyCrystal,AllSnapshots,UMDname,istep,nbox,close,min_distance)
        while True:
            if close == True:
                print("Some atoms are too close, retrying")
                close = shuffler(MyCrystal,AllSnapshots,UMDname,istep,nbox,close,min_distance)
            if close == False:
                print("Check done, writing POSCAR")
                exit()         
    else:
        print ('the umdfile ',UMDname,' does not exist')
        sys.exit()

if __name__ == "__main__":
   main(sys.argv[1:])

'''   
--------------------------------------
Code outline:
    
User gives an input UMD file, the desired snapshot number and has a choice for how to pick boxes to shuffle
Using the umd2poscar we get the poscar first
The user either specififies a number of boxes (x=y=z) so then the box length is lattice length / # boxes
Or the user specifies the box length, where it then modulus operator (%) is applied to get # boxes (rounded down), used boxes is given
We have a matrix containing all the x,y,z coordinates, we check each index of xred for each iatom to give it a box index, so we can see iatom p has box index q

'''
        
