#!/usr/bin/env python3
###
##AUTHORS: RAZVAN CARACAS
###

import sys,getopt,os.path,math
import numpy as np
from scipy.spatial.transform import Rotation
import crystallography as cr
import umd_processes_fast as umdpf



def ReadMoleculesFile(MoleculesFile):
    Multi_read = []
    AllIns_read = []
    MyMolec_read = cr.Lattice()
    flagatom = -1
    iatom = -1
    with open(MoleculesFile,'r') as ff:   #read the head and half of the 1st iter in order to get acell
        while True:
            line = ff.readline()
            if not line: break
            line=line.strip()
            entry=line.split()
#            print('current line is',line)
            if len(entry)>1:
                if entry[0] == 'ntime':
                    Multi_read.append(int(entry[1]))
#                    print('current molecule will be inserted',int(entry[1]),' times')
                if flagatom > 0:
                    iatom += 1
 #                   print('current line is :',line,': for atom no. ',iatom)
                    MyMolec_read.atoms[iatom].symbol = str(entry[0])
                    for ii in range(3):
                        MyMolec_read.atoms[iatom].xcart[ii] = float(entry[ii+1])
                    if iatom == MyMolec_read.natom-1 :
                        iatom = -1
                        flagatom = -1
#                        print('current molecule has ',MyMolec_read.natom,' atoms')
#                        for xx in range(3):
#                            for jatom in range(MyMolec_read.natom):
#                                MyMolec_read.atoms[MyMolec_read.natom].xcart[xx] = MyMolec_read.atoms[MyMolec_read.natom].xcart[xx] + MyMolec_read.atoms[jatom].xcart[xx]
#                            MyMolec_read.atoms[MyMolec_read.natom].xcart[xx] = MyMolec_read.atoms[MyMolec_read.natom].xcart[xx]/float(MyMolec_read.natom)
                        AllIns_read.append(MyMolec_read)
                        MyMolec_read = cr.Lattice()
                if entry[0] == 'natom':
                    MyMolec_read.natom = int(entry[1])
#                    print('current molecule has ',MyMolec_read.natom,' atoms')
                    MyMolec_read.atoms = [cr.Atom() for _ in range(MyMolec_read.natom)]
                    iatom = -1
                    flagatom = 1
#    for ii in range(len(AllIns_read)):
#        print(' molecule no. ',ii,' has ',AllIns_read[ii].natom,' atoms')
    return(Multi_read,AllIns_read)

           
def BuildEmptyBox(UnitCell,TotalNoAtoms):
    #print('building the empty box')
    MyNewCrystal = cr.Lattice()
    MyNewCrystal.natom = TotalNoAtoms
    MyNewCrystal.typat = [-1 for _ in range(MyNewCrystal.natom)]
    MyNewCrystal.atoms = [cr.Atom() for _ in range(TotalNoAtoms)]
    MyNewCrystal.acell[0] = UnitCell
    MyNewCrystal.acell[1] = UnitCell
    MyNewCrystal.acell[2] = UnitCell
    NoInsertedAtoms = 0
    return(MyNewCrystal,NoInsertedAtoms)
    
    
def BuildUMDBox(MyCrystal,MyUMDStructure,TotalNoAtoms):
    #print('copying the box from the UMD file')
    MyNewCrystal = cr.Lattice()
#    MyNewCrystal = MyUMDStructure
    print('in BuildUMDbox: MyCrystal.natom is ',MyCrystal.natom,'  and TotalNoAtoms is ',TotalNoAtoms)
    MyNewCrystal.natom = MyCrystal.natom + TotalNoAtoms
    MyNewCrystal.typat = [-1 for _ in range(MyNewCrystal.natom)]
    MyNewCrystal.atoms = [cr.Atom() for _ in range(MyCrystal.natom)]
    #print (' initial number of atoms is ',MyCrystal.natom)
    for iatom in range(MyCrystal.natom):
        MyNewCrystal.atoms[iatom].symbol = MyCrystal.elements[MyCrystal.typat[iatom]]
        #print (' next atom is ',MyNewCrystal.atoms[iatom].symbol)
    for iatom in range(TotalNoAtoms):
        MyNewCrystal.atoms.append(cr.Atom())
    #print ('in all there will be former ',MyCrystal.natom,' + new ',TotalNoAtoms,' = ',MyNewCrystal.natom,' atoms')
    NoInsertedAtoms = MyCrystal.natom
    return(MyNewCrystal,NoInsertedAtoms)


def PositionMolecule(MultiMolecules,AllMolecules,MyNewCrystal,MyCrystal,TotalNoAtoms,NoInsertedAtoms,Rcutoff,CurrStructs,header):
    #places the new molecules in the former structure
    #MultiMolecules stores how many molecules of each type need to be inserted
    #AllMolecules stores the actual structure of each molecule
    #TryMolec tries the position and rotation of the new molecule before approving its position
    #setting up dimensions
    print('starting to PositionMolecules with:')
    print('MyCrystal.natom is ',MyCrystal.natom)
    print('AllMolecules size is ',len(MultiMolecules))
    AtomicOrdering = []
    TryMolec = cr.Lattice()
    string = ''
    notrials = 0
    #print('already inserted', NoInsertedAtoms, '  atoms')
    #print('Inserting ',TotalNoAtoms,' atoms')
    for imolectype in range(len(MultiMolecules)):
        TryMolec.natom = AllMolecules[imolectype].natom
        TryMolec.atoms = [cr.Atom() for _ in range(TryMolec.natom)]
        print ('molecule no.', imolectype,' with ',TryMolec.natom,' atoms')
        for jmolec in range(MultiMolecules[imolectype]):
            flagpos = 1
            notrials = 0
            while flagpos>0 and notrials<10:
                notrials = notrials + 1
                print ('this is trial no. ',notrials,' for molecule no. ',jmolec)
                print ('DEALING with new molecule')
                OrigMolecule = np.random.rand(3)                #random position of central atom
                OrigRotationAxis = np.random.rand(3)          #random rotation oaxis f molecule
                OrigRotationAxis = OrigRotationAxis / np.linalg.norm(OrigRotationAxis)
                OrigRotationAngle = np.random.rand() * 2 * math.pi      #random rotation angle of molecule
                rot = Rotation.from_rotvec(OrigRotationAngle * OrigRotationAxis)
            #place all atoms of molecule
                print('random part. XYZ: ',OrigMolecule,' axis: ',OrigRotationAxis,' angle: ',OrigRotationAngle)
                for iatom in range(TryMolec.natom):
                    flagpos = 0
                    TryMolec.atoms[iatom].symbol = AllMolecules[imolectype].atoms[iatom].symbol
                    TryMolec.atoms[iatom].xcart = AllMolecules[imolectype].atoms[iatom].xcart
                    #print('before rotation :    ',TryMolec.atoms[iatom].xcart)
                    TryMolec.atoms[iatom].xcart = rot.apply(AllMolecules[imolectype].atoms[iatom].xcart)
                    #print ('after rotation :    ',TryMolec.atoms[iatom].xcart)
                    for ii in range(3):
                        TryMolec.atoms[iatom].xcart[ii] = TryMolec.atoms[iatom].xcart[ii] + OrigMolecule[ii]*MyNewCrystal.acell[ii]
                        
                    #checking to previous NoInsertedAtoms atoms
                    #print ('current atom ',iatom ,' after rotation :    ',TryMolec.atoms[iatom].symbol,TryMolec.atoms[iatom].xcart)
                    #print ('checking with respect to ',NoInsertedAtoms,' previous atoms that are at positions:')
                    #for jatom in range(NoInsertedAtoms):
                        #print (MyNewCrystal.atoms[jatom].symbol,'  ',MyNewCrystal.atoms[jatom].xcart)
                        
                    for icheckatom in range(NoInsertedAtoms):
                        #print ('comparing with ',NoInsertedAtoms,' previsouly inserted atoms')
                        for ix in range(-1,2):
                            for iy in range(-1,2):
                                for iz in range(-1,2):
                                    atref = [MyNewCrystal.atoms[icheckatom].xcart[0] + MyNewCrystal.acell[0]*float(ix), MyNewCrystal.atoms[icheckatom].xcart[1] + MyNewCrystal.acell[1]*float(iy), MyNewCrystal.atoms[icheckatom].xcart[2] + MyNewCrystal.acell[2]*float(iz)]
                                    rdist = np.linalg.norm(TryMolec.atoms[iatom].xcart - atref)
                                    #print ('distance to atom ',icheckatom,' at pos: ',MyNewCrystal.atoms[icheckatom].xcart,' = ',rdist)
                                    if rdist < float(Rcutoff):
                                        #print ('atom ',iatom,' of current molecule is at a distance of ',rdist)
                                        flagpos = flagpos + 1
                                        break
                    #print ('flagpos is ',flagpos)
                    if flagpos >0 :
                        break
                        
            if flagpos == 0:
                print ('adding new molecule no ',jmolec,' of type ',imolectype,' after ',notrials,' trials')
                for iatom in range(TryMolec.natom):
                    #print ('adding new atom no. ',NoInsertedAtoms,TryMolec.atoms[iatom].symbol,' at ',TryMolec.atoms[iatom].xcart)
                    MyNewCrystal.atoms[NoInsertedAtoms] = TryMolec.atoms[iatom]
                    NoInsertedAtoms = NoInsertedAtoms + 1
                    #print ('newly added atom ',MyNewCrystal.atoms[iatom].symbol,' at:    ',MyNewCrystal.atoms[iatom].xcart)
                    #print ('updated total number of inserted atoms: ',NoInsertedAtoms)
                TryMolec.natom = AllMolecules[imolectype].natom
                TryMolec.atoms = [cr.Atom() for _ in range(TryMolec.natom)]
            else:
                print('I could not insert your molecule even after 1000 trials. Most likely there is not enough space in the simulation cell to accomodate all molecules within an exclusion radius of ',Rcutoff)
                print('  I suggest you increase the unit cell UnitCell or the exclusion radius Rcutoff.')
                exit()

    print('done with generating, we move to ordering')
    AtomicOrdering = umdpf.sort_umd(MyNewCrystal)
    #print ('Ordered atoms are',AtomicOrdering)
    #for iatom in range(MyNewCrystal.natom):
        #print ('Atom no. ',iatom,' with symbol ',MyNewCrystal.atoms[iatom].symbol,' of type ',MyNewCrystal.typat[iatom])
    #print ('Writing ',MyNewCrystal.natom,' atoms in the ',filename,' XYZ file')
    if notrials < 1000 and flagpos == 0:
        filename = 'struct-' + str(CurrStructs) + '.xyz'
        f = open(filename,'w')
        f.write(str(MyNewCrystal.natom))
        string = '\nWriting ' + str(MyNewCrystal.natom) + ' with unit cell ' + str(MyNewCrystal.acell[0]) + ' ' + str(MyNewCrystal.acell[1]) + ' ' + str(MyNewCrystal.acell[2]) + '\n'
        f.write(string)
        string = ''
        print ('Writing ',MyNewCrystal.natom,' atoms in the ',filename,' XYZ file with unit cell ',MyNewCrystal.acell[0],' ',MyNewCrystal.acell[1],' ',MyNewCrystal.acell[2])
        for iatom in range(MyNewCrystal.natom):
            string = string + MyNewCrystal.atoms[AtomicOrdering[iatom]].symbol + '    '
            for ii in range(3):
                string = string + str(MyNewCrystal.atoms[AtomicOrdering[iatom]].xcart[ii]/MyNewCrystal.acell[ii]) + '  '
            string = string + '\n'
        f.write(string)
        f.close()
#        filename = 'struct-' + str(CurrStructs) + '.POSCAR'
#        f = open(filename,'w')
#        string = header + '  No. ' + str(CurrStructs) + '\n'
#        f.write(string)
#        f.write('  1.0\n')
#        string = '    ' + MyNewCrystal.acell[0]) + '0.0  0.0 \n'
#        f.write(string)
#        string = '    0.0 + ' + MyNewCrystal.acell[1]) + ' 0.0 \n'
#        f.write(string)
#        string = '    0.0  0.0' + MyNewCrystal.acell[02) + ' \n'
#        f.write(string)
#        for ii in range(MyNewCrystal.

def main(argv):
    umdpf.headerumd()
    UMDname='output.umd.dat'
    MoleculesFile='molecules.dat'
    Ksteps = -1
    UnitCell = 20.0
    Rcutoff = 2.0
    header = ''
    MyUMDStructure = cr.Lattice()
    MyCrystal = cr.Lattice()
    AllSnapshots = [cr.Lattice]
    AllMolecules = [cr.Lattice]
    MultiMolecules = [cr.Lattice]
    try:
        opts, arg = getopt.getopt(argv,"hf:s:i:a:r:",["fUMDfile","sSampling_Frequency","iMoleculesFile","aCubicUnitCell","rrCutoff"])
    except getopt.GetoptError:
        print ('insert_umd.py -f <UMD_filename> -s <Sampling_Frequency> -i <File_with_list_of_Molecules -a <Cubic_unit_cell> -r <cutoff_Radius>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('insert.py program to insert molecular impurities in a structure')
            print ('insert_umd.py -f <UMD_filename> -s <Sampling_Frequency> -i <File_with_list_of_Molecules -a <Cubic_unit_cell> -r <cutoff_Radius>')
            print (' -f specifies the UMD file')
            print (' -s option gies the sampling frequency. Default: -1 ')
            print ('    -1 generates one new structure. inserts molecules in an empty box')
            print ('    0 inserts molecules in the last snapshot')
            print ('    Ksteps inserts molecules every Ksteps snapshots')
            print (' -i the file containing the number and type of molecules to be inserted. Default = molecules.dat')
            print (' -a the size of the empty cubic unit cell, Default 20 angstroms ')
            print (' -r cutoff radius to prevent overlap of atomic spheres. Default 2.0 angstroms.')
            sys.exit()
        elif opt in ("-f", "--fUMDfile"):
            UMDname = str(arg)
            print('The umd file is ',UMDname)
            header = header + 'FILE: -f=' + UMDname
        elif opt in ("-i", "--iMoleculesFile"):
            MoleculesFile = str(arg)
            header = header + ' -i=' + MoleculesFile
            print ('The list with molecules is read from file ',MoleculesFile)
        elif opt in ("-a","--aCubicUnitCell"):
            UnitCell = float(arg)
            header = header + ' -s=' + arg
            print('Will insert molecules in a cubic unit cell of ',UnitCell,' angstroms')
        elif opt in ("-r","--rrCutoff"):
            Rcutoff = float(arg)
            header = header + ' -s=' + arg
            print('The cutoff distance for insertion is ',Rcutoff,' angstroms')
        elif opt in ("-s","--sKsteps"):
            Ksteps = int(arg)
            header = header + ' -s=' + arg
     #print('Parameters: UMDfile MoleculesFile CellUnit Frequency',UMDname, MoleculesFile, UnitCell, Ksteps,Rcutoff)
    
    #checks and reads the molecules.dat file
    if (os.path.isfile(MoleculesFile)):
        (MultiMolecules,AllMolecules) = ReadMoleculesFile(MoleculesFile)
        TotalNoAtoms = 0
        for imolectype in range(len(MultiMolecules)):
            TotalNoAtoms = TotalNoAtoms + MultiMolecules[imolectype] * AllMolecules[imolectype].natom
    else:
        print ('the Molecules File ',MoleculesFile,' does not exist')
        sys.exit()

    #insert molecules
    
    #inserts molecules in an empty box
    if Ksteps == -1:
        print('There are ',MultiMolecules,' molecules to be inserted in a cube with side ',UnitCell,' angstroms')
        #for ii in range(len(MultiMolecules)):
            #print('There are ',MultiMolecules[ii],' of molecules no. ',ii,' having no. atoms ',AllMolecules[ii].natom)
            #for kk in range(AllMolecules[ii].natom):
                #print(AllMolecules[ii].atoms[kk].symbol)
        (MyUMDStructure,NoInsertedAtoms) = BuildEmptyBox(UnitCell,TotalNoAtoms)
        CurrStructs = 0
        PositionMolecule(MultiMolecules,AllMolecules,MyUMDStructure,MyCrystal,TotalNoAtoms,NoInsertedAtoms,Rcutoff,CurrStructs,header)
    
    #inserts molecules in the last full snapshot
    elif Ksteps == 0:    #inserts molecules in the  last snapshot of the UMD file
        if not os.path.isfile(UMDname):
            print ('the UMD files ',UMDname,' does not exist')
            sys.exit()
        else:
            (MyCrystal, AllSnapshots, TimeStep, length)=umdpf.read_values(UMDname, "everything", mode="line", Nsteps=1)
            print ('The length of the simulation is ',len(AllSnapshots),' snapshots')
            print ('I will insert molecules in the last snapshot of the ',UMDname,' structure with ',MyCrystal.natom,' atoms')
            MyUMDStructure = AllSnapshots[len(AllSnapshots)-1]
            (MyNewCrystal,NoInsertedAtoms) = BuildUMDBox(MyCrystal,MyUMDStructure,TotalNoAtoms)
            CurrStructs = len(AllSnapshots)
            PositionMolecule(MultiMolecules,AllMolecules,MyNewCrystal,MyCrystal,TotalNoAtoms,NoInsertedAtoms,Rcutoff,CurrStructs,header)
    else:               #inserts molecules in the UMD file
        if not os.path.isfile(UMDname):
            print ('the UMD files ',UMDname,' does not exist')
            sys.exit()
        else:
            print ('I will insert molecules in ',os.path.isfile(UMDname),' structure every ',Ksteps,' steps')
            (MyCrystal,AllSnapshots,TimeStep)=umdpf.readumd(UMDname)
            print ('The length of the simulation is ',len(AllSnapshots),' snapshots')
            firststep = 0
            laststep = len(AllSnapshots)
            for istep in range(firststep,laststep,Ksteps):
                MyUMDStructure = AllSnapshots[istep]
                (MyNewCrystal,NoInsertedAtoms) = BuildUMDBox(MyCrystal,MyUMDStructure,TotalNoAtoms)
                PositionMolecule(MultiMolecules,AllMolecules,MyNewCrystal,MyCrystal,TotalNoAtoms,NoInsertedAtoms,Rcutoff,istep,header)


    
if __name__ == "__main__":
    main(sys.argv[1:])
