#!/usr/bin/env python


###
##AUTHORS: RAZVAN CARACAS, ANAIS KOBSCH
###

import sys,getopt,numpy,os.path,math,itertools
import crystallography as cr
import umd_process as umd
import datetime
from subprocess import call



def msd(MyCrystal,AllSnapshots,TimeStep,hh,vv,ballistic,umdfile):
    #compute msd

#    print 'the ii loop will start at ',ballistic+hh,' go to ',niter,' with a step of ',hh
#    print 'the jj loop will start at ii +',ballistic+vv,' go to ',niter,' with a step of ',vv
    
    niter = len(AllSnapshots)
    msd = [[0.0 for x in range(niter-ballistic+1)] for x in range(MyCrystal.ntypat+1)]
    weight = [0 for x in range(niter)]
    counter = [0 for x in range(niter-ballistic)]
    ref = [0.0 for x in range(3)]
    currZ = 0
    for iatom in range(MyCrystal.natom):
        currZ = MyCrystal.typat[iatom]
        for ii in range(ballistic+hh,int(niter/2),hh):          #initial time t
            for kk in range(3): 
                ref[kk]=AllSnapshots[ii].atoms[iatom].absxcart[kk]     #coorindates of the reference atom at time t
        #print('reference is ',ref)
            for jj in range(ballistic,int(niter/2),vv):         #new time tao
                msd[currZ][jj] = msd[currZ][jj] + (AllSnapshots[ii+jj].atoms[iatom].absxcart[0]-ref[0])**2 +(AllSnapshots[ii+jj].atoms[iatom].absxcart[1]-ref[1])**2 +(AllSnapshots[ii+jj].atoms[iatom].absxcart[2]-ref[2])**2       #distance between t+tao - t
                #print(iatom,ii,ii+jj,AllSnapshots[ii+jj].atoms[iatom].absxcart,'distant to',AllSnapshots[ii].atoms[iatom].absxcart)
                weight[jj]=weight[jj]+1                           #normalization to the number of t times (larger tao's are counted fewer times) mutiplied by the number of atoms .
                                                                #when applyging the final normalization this number has to be divided by the total number of atoms to give only the number of t's:
                                                                # weight = sum_no.atoms sum_t's (0 ... tao)
                                                                # normalization is no. of taos
                                                                # which is this weight divided by no.atoms
    msdfile = umdfile[:-8] + '.msd.dat'
    f = open(msdfile,'w')
    string='time_(fs)\t'
    for itypat in range(MyCrystal.ntypat):
        string=string + MyCrystal.elements[itypat] + '\t'
#    print 'header for msd is \n',string
#    print string
    string = string + '\n'
    f.write(string)
#    diffcoeff = [[0.0 for x in range(niter/2)] for x in range(len(typat))]
    for ii in range(int(niter/2)):                              # all possible times, max is total_simulation_time / 2
        if (weight[ii]>0):                                      # further check to be sure that tao's have been counted
            string = str(float(ii)*TimeStep)
#            stringD = ''
            for jj in range(MyCrystal.ntypat):
#                print ' iter ',ii,'current element',Elements[jj-1],'msd',str(msd[jj][ii]),'weighted by no.atoms',str(msd[jj][ii]/float(typat[jj-1]))
                string = string + '\t' + str(msd[jj][ii]/(float(MyCrystal.types[jj])*float(weight[ii])/float(MyCrystal.natom)))
#                print(msd[jj][ii],MyCrystal.types[jj],weight[ii],MyCrystal.natom)
#                stringD = stringD + '\t' + str( (msd[jj][ii] / (float(typat[jj-1])*float(weight[ii])))/float(ii))
#                string = string + '  ' + Elements[jj-1]+'\t' + str(msd[jj][ii]) + '\t' + str(msd[jj][ii]/float(weight[ii]))+'\t' + str((msd[jj][ii]/float(weight[ii]))/float(typat[jj-1]))
#            print string
#                diffcoeff[jj-1][ii] = 1E-05*msd[jj][ii] / (6*float(typat[jj-1])*float(weight[ii]/natom)*float(ii))
            string = string + '\n'
            f.write(string)
#            print ii,' steps were counted ',weight[ii],' times'
#    for jj in range(len(typat)):
#        print numpy.mean(diffcoeff[jj][:]),'\t',numpy.std(diffcoeff[jj][:])
    print ('MSDs printed in file ',msdfile)

def main(argv):
    XYZfile='file.xyz'
    hh = 1
    vv = 1
    ww = 1
    natom = 0
    typat =[]
    Elements = []
    ballistic = 0
    TimeStep = 1
    umd.headerumd()
    try:
        opts, arg = getopt.getopt(argv,"hf:z:v:b:",["fumdfile","zHorizontalJump","vVerticalJump","bBallistic"])
    except getopt.GetoptError:
        print ('msd_umd.py -f <XYZ_filename> -z <HorizontalJump> -v <VerticalJump> -b <Ballistic>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('Program to compute the Mean Square Displacement. Usage: ')
            print ('msd_umd.py -f <UMD_filename> -z <HorizontalJump> -v <VerticalJump> -b <Ballistic>')
            print ('UMD_filename = input file with the trajectory, in UMD format.')
            print (' As the MSD is measured with a sliding window of various size up to half the trajectory\'s length, ')
            print (' we can accelerate the calculation using a selected reduced sampling ')
            print ('HorizontalJump = discretization for the start of the sampling window.')
            print ('VerticalJump = discretization for the length of the sampling window.')
            print ('Ballistic = estimation of the ballistic part of the trajectory. Default is 0. Typical values of 100 are sufficient.')
            sys.exit()
        elif opt in ("-f", "--fumdfile"):
            umdfile = str(arg)
            print ('I will use the ',umdfile,' for input')
        elif opt in ("-z", "--zHorizontalJump"):
            hh = int(arg)
        elif opt in ("-v", "--vVerticalJump"):
            vv = int(arg)
        elif opt in ("-b", "--bBallistic"):
            ballistic = int(arg)
    if (os.path.isfile(umdfile)):
#        initstruct(XYZfile)
        MyCrystal = cr.Lattice()
        AllSnapshots = [cr.Lattice]
        (MyCrystal,AllSnapshots,TimeStep)=umd.read_absxcart(umdfile)
#        (MyCrystal,AllSnapshots,TimeStep)=up.readumd(umdfile)
#        print 'Elements are ',elem
        print ('Number of atoms of each type is ',MyCrystal.types)
#        print 'Atomic types are ',znucl
        msd(MyCrystal,AllSnapshots,TimeStep,hh,vv,ballistic,umdfile)
    else:
        print ('umd file ',umdfile,'does not exist')
        sys.exit()
 

if __name__ == "__main__":
   main(sys.argv[1:])

