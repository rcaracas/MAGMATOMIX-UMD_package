#!/usr/bin/env python3

###
##AUTHORS: RAZVAN CARACAS, ANAIS KOBSCH, NATALIA SOLOMATOVA
###

import sys,getopt,os.path,math
import crystallography as cr
import umd_process as umdp

def read_outcar(FileName,InitialStep):
    #read poscar file
    print ('Reading outcar file ',FileName,' from initial step ',InitialStep)
    switch = False
    flagrprimd = -1
    iatom = -1    
    istep = 0
    MyCrystal = cr.Lattice()
    TimeStep = 1.0
    newfile = FileName + '.umd.dat'
    nf = open(newfile,'w')
    nf.close()
    atomictype = 0
    CurrentTime = 0.0
    TimeStep = 0.0
    flagscalee = -1
    with open(FileName,'r') as ff:
        while True:
            line = ff.readline()
            if not line: break
            line=line.strip()
            entry=line.split()
            if (len(entry)>0):
                if (entry[0]=='ions'):          #determine how many atoms
                    jatom = 0
                    MyCrystal.natom = 0
                    MyCrystal.ntypat = len(entry)-4
                    MyCrystal.types = [0 for _ in range(MyCrystal.ntypat)]
                    for ii in range(len(entry)-4):
                        MyCrystal.natom = MyCrystal.natom + int(entry[4+ii])
                        MyCrystal.types[ii] = int(entry[4+ii])
                    
                    MyCrystal.typat = [0 for _ in range (MyCrystal.natom)]
                    for ii in range(len(entry)-4):
                        atomictype += 1
                        for jj in range(int(entry[4+ii])):
                            MyCrystal.typat[jatom] = atomictype-1
                            jatom = jatom + 1
                    MyCrystal.atoms = [cr.Atom() for _ in range(MyCrystal.natom)]
                    oldpos = [cr.Atom() for _ in range(MyCrystal.natom)]
                    MyCrystal.elements = ['X' for _ in range(MyCrystal.ntypat)]
                    MyCrystal.masses = [0.0 for _ in range(MyCrystal.ntypat)]
                    MyCrystal.zelec = [0.0 for _ in range(MyCrystal.ntypat)]
                    MyCrystal.stress = [0.0 for _ in range(6)]
                    phase = [[0.0,0.0,0.0] for _ in range(MyCrystal.natom)]
                    diffcoords = [[0.0,0.0,0.0] for _ in range(MyCrystal.natom)]
                    MyCrystal.magnetization = 0.0
                    print ('Total number of atoms = ',len(diffcoords))
                    break
    ff.close()
    flagmass = 0
    flagtitle = 0
    atn = ''
    ats = ''
    atno = ''
    with open(FileName,'r') as ff:
        while True:
            line = ff.readline()
            if not line: break
            line=line.strip()
            entry=line.split()
            if (len(entry)>0):
                if (len(entry)>1):
                    if entry[1] == 'Iteration':
                        break
#                if (entry[0]=='VRHFIN'):
#                    flagtitle += 1
#                    atomicsymbol = entry[1]
#                    MyCrystal.elements[flagtitle-1] = atomicsymbol[1:-1]
#                    print ('element ',flagtitle,' is ',MyCrystal.elements[flagtitle-1])
#                    (atn,ats,atno,MyCrystal.masses[flagtitle-1])=cr.Elements2rest(MyCrystal.elements[flagtitle-1])
#                    print ('just chekcing: ',atn,ats,atno)
                if (entry[0]=='POTCAR:'):
                    if (flagtitle < MyCrystal.ntypat) :
                        MyCrystal.elements[flagtitle] = entry[2].split('_')[0]
#                        print ('element ',flagtitle,' is ',MyCrystal.elements[flagtitle])
                        (atn,ats,atno,MyCrystal.masses[flagtitle])=cr.Elements2rest(MyCrystal.elements[flagtitle])
#                        print ('just chekcing: ',atn,ats,atno)
                        print ('element ',flagtitle,' is ',MyCrystal.elements[flagtitle],' with atomic number ',atno,' and mass ',MyCrystal.masses[flagtitle])
                        flagtitle += 1
                if (entry[0]=='POMASS'):
                    if flagmass < MyCrystal.ntypat:
#                        print('flagmass is',flagmass,' with the line ',entry)
#                        for ii in range(MyCrystal.ntypat):
                        MyCrystal.masses[flagmass] = float(entry[2][:-1])
                        MyCrystal.zelec[flagmass] = float(entry[5])
                        flagmass += 1
#                if (entry[0]=='ZVAL'):
#                    for ii in range(MyCrystal.ntypat):
#                        MyCrystal.zelec[ii]=float(entry[ii+2])
#                if (entry[0]=='Mass'):
#                    flagmass=1
                if (entry[0]=='NELECT'):
                    MyCrystal.noelectrons = float(entry[2])
                if (entry[0]=='SCALEE'):
                    MyCrystal.lambda_ThermoInt = float(entry[2])



    ff.close()
    umdp.print_header(FileName,MyCrystal)
    with open(FileName,'r') as ff:
        while True:
            line = ff.readline()
            if not line: break
            line=line.strip()
            entry=line.split()
            if (len(entry) == 6):
#                print ('a line of 6 elements: ',entry)
                if (entry[4] == 'magnetization'):
#                    print ('the line of 6 elements: ',entry)
                    MyCrystal.magnetization = float(entry[5])
            if(switch==False):
                if (len(entry)>0 and entry[0]=='POTIM'):
                    TimeStep=float(entry[2])
                if (len(entry)>=2 and entry[1]=='aborting'):
                    switch = True
                else:
                    continue
            if (switch==True and len(entry)>0):
#                print 'entry[0] is ',entry[0]
                jatom = 0                     
                if (entry[0]=='Total+kin.'):            #reading stress tensor, corrected for the internal kinetic pressure (the term contains it itself0
                    for ii in range(6):
                        MyCrystal.stress[ii]=float(entry[ii+1])/10.0            #transforms from kbars into GPa
                    MyCrystal.pressure = (float(entry[1])+float(entry[2])+float(entry[3]))/30.0 #transforms also from kbars into GPa
                if entry[0] == 'energy':   
                    if entry[2] == 'entropy=': #reading Kohn-Sham energy, contains all the electronic energy without the electronic entropy 
#                        MyCrystal.internalenergy=float(entry[3])  #this leads only part of the internal energy, one should add the kinetic energy of ions
                                                                  #the variance of {this energy + kinetic energy of ions}  yields Cv
                        flagscalee += 1
#                        print('reading energy withtout entropy, energy, lambda,flagscalee = ',MyCrystal.lambda_ThermoInt,entry[3],flagscalee)
                        if MyCrystal.lambda_ThermoInt < 1.0:
                            if flagscalee == 0:
                                MyCrystal.internalenergy=float(entry[3])
#                            if flagscalee == 1:
                            else:
                                flagscalee = -1
                        else:
                            MyCrystal.internalenergy=float(entry[3])
                if (entry[0]=='%'):
                    if (entry[1]=='ion-electron'):   #the term T*Sel in the formula F = E - T*Sel, with F = Kohn-Sham energy and E = Kohn-Sham energy without the electronic entropy 
                        MyCrystal.electronicentropy=MyCrystal.internalenergy - float(entry[4])
                if (entry[0] == 'magnetization'):       #reading magnetization of individual atoms
                    #print('reading magnetization')
                    line = ff.readline()
                    line = ff.readline()
                    line = ff.readline()
                    for ii in range(MyCrystal.natom):
                        line = ff.readline()
                        line=line.strip()
                        entry=line.split()
                        MyCrystal.atoms[ii].magnet = float(entry[len(entry)-1])
                        #print ('for atom ',iatom,' magnetization is ',MyCrystal.atoms[ii].magnet)
#                if (len(entry) == 6):
#                    print ('a line of 6 elements: ',entry)
#                    if (entry[0] == 'number'):
#                        print ('the line of 6 elements: ',entry)
#                        MyCrystal.magnetization = float(entry[5])
                if line == 'total charge':           #reading the atomi charges
                    #print('reading magnetization')
                    line = ff.readline()
                    line = ff.readline()
                    line = ff.readline()
                    for ii in range(MyCrystal.natom):
                        line = ff.readline()
                        line=line.strip()
                        entry=line.split()
                        MyCrystal.atoms[ii].charge = float(entry[len(entry)-1])
                if (entry[0] == 'kinetic'):               #reading kinetic energy of ions
                    if (len(entry) > 2):
                        if (entry[2] == 'EKIN'):
                            MyCrystal.kineticenergy = float(entry[4])
                if (entry[0] == 'total'):               #reading energy
                                                        #KS energy + thermostat + ion-kinetic terms
                                                        #used to check the drift in energy over a simulation 
                    if (len(entry) > 2):
                        if (entry[2] == 'ETOTAL'):
                            MyCrystal.energywithdrift = float(entry[4])
                            #ETOTAL is the last value we want to take in the umd, so we can put the switch off and increase the counter of steps
                            #Once istep >= InitialStep, we can compute the velocities, diffcoord etc. and print the umd
                            switch = False
                            istep = istep + 1
                            if istep>=InitialStep:
                                #print (' treating step no. ', istep)
                                for jatom in range(MyCrystal.natom):
                                    for ii in range(3):
                                        jump = MyCrystal.atoms[jatom].xcart[ii] - oldpos[jatom].xcart[ii]
                                        if jump > MyCrystal.acell[ii]/2:
                                            #print ('positive jump',jump,jatom,MyCrystal.atoms[jatom].xcart[ii],oldpos[jatom].xcart[ii])
                                            phase[jatom][ii] = phase[jatom][ii] - (MyCrystal.rprimd[ii][0]+MyCrystal.rprimd[ii][1]+MyCrystal.rprimd[ii][2])
                                            diffcoords[jatom][ii] = MyCrystal.atoms[jatom].xcart[ii] + phase[jatom][ii]
                                            jump = jump - (MyCrystal.rprimd[ii][0]+MyCrystal.rprimd[ii][1]+MyCrystal.rprimd[ii][2])
                                        elif jump < -MyCrystal.acell[ii]/2:
                                            #print ('negatie jump',jump,jatom,MyCrystal.atoms[jatom].xcart[ii],oldpos[jatom].xcart[ii])
                                            phase[jatom][ii] = phase[jatom][ii] + (MyCrystal.rprimd[ii][0]+MyCrystal.rprimd[ii][1]+MyCrystal.rprimd[ii][2])
                                            jump = jump + (MyCrystal.rprimd[ii][0]+MyCrystal.rprimd[ii][1]+MyCrystal.rprimd[ii][2])
                                            diffcoords[jatom][ii] = MyCrystal.atoms[jatom].xcart[ii] + phase[jatom][ii]
                                        else:
                                            diffcoords[jatom][ii] = MyCrystal.atoms[jatom].xcart[ii] + phase[jatom][ii]
                                        MyCrystal.atoms[jatom].vels[ii] = jump/TimeStep
                                        for jj in range(3):
                                            MyCrystal.atoms[jatom].xred[ii] = MyCrystal.atoms[jatom].xred[ii] + MyCrystal.gprimd[ii][jj]*MyCrystal.atoms[jatom].xcart[ii]
                                            while MyCrystal.atoms[jatom].xred[ii] >=1.0:
                                                MyCrystal.atoms[jatom].xred[ii] = MyCrystal.atoms[jatom].xred[ii] - 1.0
                                (CurrentTime,TimeStep) = umdp.print_snapshots(FileName,MyCrystal,TimeStep,(istep-InitialStep)*TimeStep,diffcoords)
                if (entry[0]=='kin.'):          #reading the temperature
                    if len(entry) == 7:
                        MyCrystal.temperature = float(entry[5])
                    elif len(entry) == 6:
                        mixedtemp = entry[4]
                        MyCrystal.temperature = float(mixedtemp[-8:])
                    else:
                        print('defect on the temperature line',entry)
                if (flagrprimd>-1):             #reaeding the unit cell
#                    print 'entries are ',entry
                    MyCrystal.rprimd[flagrprimd][0]=float(entry[0].split(',')[0])
                    MyCrystal.rprimd[flagrprimd][1]=float(entry[1].split(',')[0])
                    MyCrystal.rprimd[flagrprimd][2]=float(entry[2].split(')')[0])
                    MyCrystal.acell[flagrprimd]=math.sqrt(MyCrystal.rprimd[flagrprimd][0]*MyCrystal.rprimd[flagrprimd][0] + MyCrystal.rprimd[flagrprimd][1]*MyCrystal.rprimd[flagrprimd][1] + MyCrystal.rprimd[flagrprimd][2]*MyCrystal.rprimd[flagrprimd][2])
                    flagrprimd = flagrprimd+1
                if (flagrprimd == 3):
                    MyCrystal.gprimd = MyCrystal.makegprimd()
                    #MyCrystal.cellvolume = MyCrystal.makevolume()
                    MyCrystal.density = MyCrystal.getdensity()
                    #print('gprimd  ',MyCrystal.gprimd)
                    #print(MyCrystal.cellvolume)
                    flagrprimd = -1
                if (entry[0]=='direct'):
                    flagrprimd = flagrprimd + 1
                if (iatom>-1):                  #reading the atomic positions
                    #print ('current iatom is ',iatom)
                    MyCrystal.atoms[iatom].xred = [0.0,0.0,0.0]
                    MyCrystal.atoms[iatom].vels = [0.0,0.0,0.0]
                    oldpos[iatom].xcart = MyCrystal.atoms[iatom].xcart
                    MyCrystal.atoms[iatom].xcart = [float(entry[0]),float(entry[1]),float(entry[2])]
                    MyCrystal.atoms[iatom].forces = [float(entry[3]),float(entry[4]),float(entry[5])]
                    iatom = iatom + 1
                    if (iatom==MyCrystal.natom):
                        iatom = -1
                if (entry[0]=='POSITION'):
                    line = ff.readline()
                    iatom = 0
    return(CurrentTime,TimeStep)



def main(argv):
    OUTCARname='OUTCAR'
    Nsteps = 1
    InitialStep = 0
    CurrentTime = 0.0
    TimeStep = 0.0
    string = ''
    umdp.headerumd()
    try:
        opts, arg = getopt.getopt(argv,"hf:i:",["fOUTCARfile","iInitialStep"])
    except getopt.GetoptError:
        print ('VaspParser.py -f <OUTCAR_filename> -i <InitialStep>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('VaspParser.py program to translate VASP outcar files into UMD ascii format file')
            print ('VaspParser.py -f <OUTCAR_filename> -i <InitialStep>')
            print ('   default values: -f OUTCAR -i 0')
            sys.exit()
        elif opt in ("-f", "--fOUTCARfile"):
            OUTCARname = str(arg)
        elif opt in ("-i", "--iInitialStep"):
            InitialStep = int(arg)
    if (os.path.isfile(OUTCARname)):
        (CurrentTime,TimeStep) = read_outcar(OUTCARname,InitialStep)
        string = 'Total simulation time ' + str(CurrentTime) + ' fs, done in ' + str(int(CurrentTime/TimeStep)) + ' steps of ' + str(TimeStep) + ' fs each.'
        print(string)
    else:
        print ('Error. The outcar files ',OUTCARname,' does not exist')
        sys.exit()

 

if __name__ == "__main__":
   main(sys.argv[1:])



