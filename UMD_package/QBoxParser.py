#!/usr/bin/env python3
###
##AUTHORS: RAZVAN CARACAS
###

import sys,getopt,os.path,math
import crystallography as cr

def print_umd(FileName,MyCrystal,TimeStep,CurrentTime):
    newfile = FileName + '.umd.dat'
    nf = open(newfile,'a')
    string = 'timestep ' + str(round(TimeStep,4)) + ' fs\n' + 'time ' + str(round(CurrentTime,4)) + ' fs\n'
    nf.write(string)
    string = 'InternalEnergy ' + str(round(MyCrystal.internalenergy,6)) + ' eV\n'
    nf.write(string)
    string = 'EnergyWithDrift ' + str(round(MyCrystal.energywithdrift,6)) + ' eV\n'
    nf.write(string)
    string = 'Enthalpy ' + str(round(MyCrystal.enthalpy,6)) + ' eV\n'
    nf.write(string)
    string = 'Temperature ' + str(round(MyCrystal.temperature,1)) + ' K\n'
    nf.write(string)
    string = 'Pressure ' + str(round(MyCrystal.pressure,3)) + ' GPa\n'
    nf.write(string)
    string = 'StressTensor '
    for ii in range(6):
        string = string + str(round(MyCrystal.stress[ii]/10,2)) + ' '
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
    string = 'atoms: reduced*3 cartesian*3(A) abs.diff.*3(A) velocity*3(A/fs) force*3(eV/A)  \n'
    nf.write(string)
    for iatom in range(MyCrystal.natom):
        string = str(round(MyCrystal.atoms[iatom].xred[0],5)) + ' ' + str(round(MyCrystal.atoms[iatom].xred[1],5)) + ' ' + str(round(MyCrystal.atoms[iatom].xred[2],5)) + ' '
        string = string + str(round(MyCrystal.atoms[iatom].xcart[0],5)) + ' ' + str(round(MyCrystal.atoms[iatom].xcart[1],5)) + ' ' + str(round(MyCrystal.atoms[iatom].xcart[2],5)) + ' '
        string = string + str(round(MyCrystal.atoms[iatom].absxcart[0],5)) + ' ' + str(round(MyCrystal.atoms[iatom].absxcart[1],5)) + ' ' + str(round(MyCrystal.atoms[iatom].absxcart[2],5)) + ' '
        string = string + str(round(MyCrystal.atoms[iatom].vels[0],5)) + ' ' + str(round(MyCrystal.atoms[iatom].vels[1],5)) + ' ' + str(round(MyCrystal.atoms[iatom].vels[2],5)) + ' '
        string = string + str(round(MyCrystal.atoms[iatom].forces[0],5)) + ' ' + str(round(MyCrystal.atoms[iatom].forces[1],5)) + ' ' + str(round(MyCrystal.atoms[iatom].forces[2],5)) + ' '
        string = string + str(MyCrystal.atoms[iatom].magnet) + '\n'
        nf.write(string)
    string='\n'
    nf.write(string)
    nf.close()


def print_header(FileName,MyCrystal):
    newfile = FileName + '.umd.dat'
    nf = open(newfile,'a')
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
    nf.close()


def read_qbox(FileName,InitialStep,SYSfile):
    #read poscar file
    print ('reading output file ',FileName)
    print('initializing structure from file ',SYSfile)
    print('discarding ',InitialStep,' initial steps')
    GPaAcube = 1 / ( (0.529177249**3) * 29421.033 * 27.21138602)
    BohrtoAng = 0.5291772
    switch = False
    flagrprimd = -1
    flagstart = 0
    iatom = -1    
    istep = 0
    titleflag = 0
    MyCrystal = cr.Lattice()
    TimeStep = 1.0
    newfile = FileName + '.umd.dat'
    atomname = ''
    nf = open(newfile,'w')
    nf.close()
    MyCrystal.elements = []
    MyCrystal.natom = 0
    MyCrystal.ntypat = 0
    MyCrystal.types = []
    MyCrystal.typat = []
    if SYSfile == '':
        structName = FileName
        indexatom = 1
    else:
        structName = SYSfile
        indexatom = 0
    print('reading structure from',structName)
    with open(structName,'r') as ff:
        while True:
            line = ff.readline()
            if not line: break
            #print('new line is ',line)
            line=line.strip()
            entry=line.split()
            if (len(entry)>1):
                #print('with entry[1] is ',entry[1])
                if (entry[indexatom]=='atom'):          #determine how many atoms
                    #print(' ')
                    #jatom = 0
                    MyCrystal.natom = MyCrystal.natom + 1
                    #print('total number of atoms so far is',MyCrystal.natom,' with current atom ',entry[2])
                    atomname = ''.join([kk for kk in entry[1] if not kk.isdigit()])
                    #print ('newly renamed atom is ',atomname)
                    flagexistentatom = 0
                    if len(MyCrystal.elements)>0:
                        #print('the list of elements is not empty, but contains already ',len(MyCrystal.elements),' elements :',MyCrystal.elements)
                        for jj in range(len(MyCrystal.elements)):
                            #print ('  atomname ',atomname,'to be compared with type no. ',jj,' which is ',MyCrystal.elements[jj])
                            if atomname == MyCrystal.elements[jj]:
                                flagexistentatom = 1
                                indexoftype = jj
                        if flagexistentatom == 0:
                                MyCrystal.elements.append(atomname)
                                #print('new version of mycrystal.elements', MyCrystal.elements)
                                MyCrystal.types.append(1)
                                MyCrystal.ntypat = MyCrystal.ntypat + 1
                                #print('    newly added atom. ')
                                #print('    at this point elements are: ',MyCrystal.elements)
                                #print('    current atom list with ntypat ',MyCrystal.ntypat,' and types ',MyCrystal.types)
                        else:
                            #print('    element already existent ',MyCrystal.elements[jj],' with index no. ',indexoftype)
                            MyCrystal.types[indexoftype] += 1
                    else:
                        #print ('adding first element to the empty list', MyCrystal.elements)
                        MyCrystal.elements.append(atomname)
                        #print('first version of atomic types', MyCrystal.elements)
                        MyCrystal.ntypat = MyCrystal.ntypat + 1
                        MyCrystal.types = [1]
                        #print('current ntypat and types ',MyCrystal.ntypat,MyCrystal.types)
                if (entry[1]=='</wavefunction>'):
                    break
                if (entry[1]=='End'):
                    break
    ff.close()
    MyCrystal.typat = [0 for i in range(MyCrystal.natom)]
    iatom = 0
    for itypat in range(MyCrystal.ntypat):
        for itypes in range(MyCrystal.types[itypat]):
            MyCrystal.typat[iatom] = itypat
            iatom += 1
    
    print ('MyCrystal has been initiliazed as: ')
    MyCrystal.atoms = [cr.Atom() for _ in range(MyCrystal.natom)]
    oldpos = [[0.0,0.0,0.0] for _ in range(MyCrystal.natom)]
    print ('natom = ',MyCrystal.natom)
    print ('ntypat = ',MyCrystal.ntypat)
    print ('types = ',MyCrystal.types)
    print ('elements = ',MyCrystal.elements)
    print_header(FileName,MyCrystal)

    with open(FileName,'r') as ff:
        while True:
            line = ff.readline()
            if not line: break
            line=line.strip()
            entry=line.split()
            if len(entry)==4:
                #print('line with four entries ',line)
                if entry[2] == 'dt':
                    TimeStep = float(entry[3]) * 0.02418884326505
                    #print('timestep is',TimeStep)
            if len(entry)==2:
                if entry[0] == '<iteration':
                    counter = entry[1].split("\"")
                    flagstart +=1
            if flagstart > InitialStep and len(entry)>0:
                if entry[0] == '<sigma_xx>':
                    MyCrystal.stress[0] = float(entry[1])
                if entry[0] == '<sigma_yy>':
                    MyCrystal.stress[1] = float(entry[1])
                if entry[0] == '<sigma_zz>':
                    MyCrystal.stress[2] = float(entry[1])
                    MyCrystal.pressure = (MyCrystal.stress[0]+MyCrystal.stress[1]+MyCrystal.stress[2])/3.0
                if entry[0] == '<sigma_yz>':
                    MyCrystal.stress[3] = float(entry[1])
                if entry[0] == '<sigma_xz>':
                    MyCrystal.stress[4] = float(entry[1])
                if entry[0] == '<sigma_xy>':
                    MyCrystal.stress[5] = float(entry[1])
                if entry[0] == '<temp_ion>':
                    MyCrystal.temperature = float(entry[1])
                if entry[0] == '<etotal>':
                    MyCrystal.internalenergy=float(entry[1])
                if entry[0] == 'a=\"':
                    MyCrystal.rprimd[0][0] = float(entry[1]) * BohrtoAng
                    MyCrystal.rprimd[0][1] = float(entry[2]) * BohrtoAng
                    MyCrystal.rprimd[0][2] = float(entry[3][:-1]) * BohrtoAng
                    MyCrystal.acell[0]=math.sqrt(MyCrystal.rprimd[0][0]**2 + MyCrystal.rprimd[0][1]**2 + MyCrystal.rprimd[0][2]**2)
                    MyCrystal.rprim[0][0] = MyCrystal.rprimd[0][0]/MyCrystal.acell[0]
                    MyCrystal.rprim[0][1] = MyCrystal.rprimd[0][1]/MyCrystal.acell[0]
                    MyCrystal.rprim[0][2] = MyCrystal.rprimd[0][2]/MyCrystal.acell[0]
                if entry[0] == 'b=\"':
                    MyCrystal.rprimd[1][0] = float(entry[1]) * BohrtoAng
                    MyCrystal.rprimd[1][1] = float(entry[2]) * BohrtoAng
                    MyCrystal.rprimd[1][2] = float(entry[3][:-1]) * BohrtoAng
                    MyCrystal.acell[1]=math.sqrt(MyCrystal.rprimd[1][0]**2 + MyCrystal.rprimd[1][1]**2 + MyCrystal.rprimd[1][2]**2)
                    MyCrystal.rprim[1][0] = MyCrystal.rprimd[1][0]/MyCrystal.acell[1]
                    MyCrystal.rprim[1][1] = MyCrystal.rprimd[1][1]/MyCrystal.acell[1]
                    MyCrystal.rprim[1][2] = MyCrystal.rprimd[1][2]/MyCrystal.acell[1]
                if entry[0] == 'c=\"':
                    MyCrystal.rprimd[2][0] = float(entry[1]) * BohrtoAng
                    MyCrystal.rprimd[2][1] = float(entry[2]) * BohrtoAng
                    MyCrystal.rprimd[2][2] = float(entry[3][:-1]) * BohrtoAng
                    MyCrystal.acell[2]=math.sqrt(MyCrystal.rprimd[2][0]**2 + MyCrystal.rprimd[2][1]**2 + MyCrystal.rprimd[2][2]**2)
                    MyCrystal.rprim[2][0] = MyCrystal.rprimd[2][0]/MyCrystal.acell[2]
                    MyCrystal.rprim[2][1] = MyCrystal.rprimd[2][1]/MyCrystal.acell[2]
                    MyCrystal.rprim[2][2] = MyCrystal.rprimd[2][2]/MyCrystal.acell[2]
                    #print('rprimd is',MyCrystal.rprimd)
                    MyCrystal.cellvolume = MyCrystal.makevolume()
                    #print('cellvolume is',MyCrystal.cellvolume)
                    MyCrystal.gprimd = MyCrystal.makegprimd()
                    MyCrystal.cellvolume = MyCrystal.makevolume()
                    #print('cellvolume is',MyCrystal.cellvolume)
                if entry[0] == '<atomset>':
                    jatom = -1
                if entry[0] == '<position>':
                    jatom += 1
                    MyCrystal.atoms[jatom].absxcart[0] = float(entry[1]) * BohrtoAng
                    MyCrystal.atoms[jatom].absxcart[1] = float(entry[2]) * BohrtoAng
                    MyCrystal.atoms[jatom].absxcart[2] = float(entry[3]) * BohrtoAng
                    #MyCrystal.atoms[jatom].xred = MyCrystal.absxcart2xred(MyCrystal.atoms[jatom].absxcart)
                    #MyCrystal.atoms[jatom].xcart = MyCrystal.absxred2xcart(MyCrystal.atoms[jatom].xred)
                if entry[0] == '<velocity>':
                    MyCrystal.atoms[jatom].vels[0] = float(entry[1])
                    MyCrystal.atoms[jatom].vels[1] = float(entry[2])
                    MyCrystal.atoms[jatom].vels[2] = float(entry[3])
                if entry[0] == '<force>':
                    MyCrystal.atoms[jatom].forces[0] = float(entry[1])
                    MyCrystal.atoms[jatom].forces[1] = float(entry[2])
                    MyCrystal.atoms[jatom].forces[2] = float(entry[3])
                if entry[0] == '</iteration>':
                    MyCrystal.enthalpy = MyCrystal.internalenergy + MyCrystal.pressure * MyCrystal.cellvolume * GPaAcube
                    MyCrystal.allabsxcart2xred()
                    MyCrystal.allxred2xcart()
                    print_umd(FileName,MyCrystal,TimeStep,(flagstart-InitialStep)*TimeStep)
    return(CurrentTime,TimeStep)






def main(argv):
    OUTCARname='Qbox.out'
    SYSfile = 'sys.struct'
    Nsteps = 1
    InitialStep = 0
    umd.headerumd()
    print ('Parser of QBox output files to UMD format.')
    try:
        opts, arg = getopt.getopt(argv,"hf:i:s:",["fOUTfile","iInitialStep","sSYSfile"])
    except getopt.GetoptError:
        print ('QBoxParser.py program to extract relevant information from the QBox output file and write it into a umd file')
        print ('QboxParser.py -f <OUT_filename> -i <InitialStep> -s <SYS_structfile>')
        print ('  default values: -f Qbox.out -i -s sys.struct')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('QBoxParser program to translate QBox outcar files into UMD ascii format file')
            print ('QBoxParser -f <OUT_filename> -i <InitialStep> -s <SYS_structfile>')
            print ('OUT_filename = QBox output file. Default = Qbox.outp ')
            print ('SYS_structfile = QBox-formatted structure file. Default = sys.struct ')
            print ('InitialStep = number of steps to be discarded. Default = 0')
            sys.exit()
        elif opt in ("-f", "--fOUTCARfile"):
            OUTCARname = str(arg)
        elif opt in ("-i", "--iInitialStep"):
            InitialStep = int(arg)
        elif opt in ("-s", "--sSYSfile"):
            SYSfile = str(arg)
    if (os.path.isfile(OUTCARname)):
        (CurrentTime,TimeStep) = read_qbox(OUTCARname,InitialStep,SYSfile)
        string = 'Total simulation time ' + str(CurrentTime) + ' fs, done in ' + str(int(CurrentTime/TimeStep)) + ' steps of ' + str(TimeStep) + ' fs each.'
        print(string)
    else:
        print ('the outcar files ',OUTCARname,' does not exist')
        sys.exit()

 

if __name__ == "__main__":
   main(sys.argv[1:])



