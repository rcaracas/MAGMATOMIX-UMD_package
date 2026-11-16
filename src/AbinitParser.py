#!/usr/bin/env python3
###
#AUTHORS: RAZVAN CARACAS
###
#Converts an ABINIT molecular dynamics output file into a UMD file. Layout checked on ABINIT 6.4.1, ionmov 9, constant cell ; other layouts are refused or flagged.

import sys,getopt,os.path,re
import numpy as np
import crystallography as cr
import umd_processes_fast as umdpf
from vibrations_umd import KB,AMU_A2FS2

BOHR = 0.529177210903               #bohr in A
HARTREE = 27.211386245988           #Hartree in eV
AUTIME = 2.418884326505E-02         #atomic unit of time in fs
EVA3_GPA = 160.2176634              #eV/A^3 in GPa
HAB3_GPA = HARTREE/BOHR**3*EVA3_GPA #Hartree/bohr^3 in GPa
CHECKED_IONMOV = (9,)               #the layout of the MD output was checked on these values of ionmov (the positions, forces and velocities of the first block of a step belong together)


def read_input(FileName):#The input variables echoed by ABINIT (-outvars), the cell (R(1..3), bohr) and the zion of the pseudopotentials, up to the first MD step
    variables = {}
    cell = []
    zion = []
    inoutvars = False
    foundstep = False
    keyword = ''
    with open(FileName,'r') as ff:
        for line in ff:
            if line.startswith(' MOLDYN STEP NUMBER'):
                foundstep = True
                break
            if '-outvars: echo values of preprocessed' in line:
                inoutvars = True
                continue
            if inoutvars and (line.startswith('====') or 'chkinp' in line):
                inoutvars = False
            if inoutvars:
                entry = line.split()
                if len(entry) > 0 and entry[0] in ('P','-'):      #memory and parallelism flags in front of the variable
                    entry = entry[1:]
                if len(entry) > 0 and entry[0][0].isalpha():
                    keyword = entry[0]
                    variables[keyword] = entry[1:]
                elif keyword != '':
                    variables[keyword] += entry
            match = re.search(r'R\(\d\)=\s*(\S+)\s+(\S+)\s+(\S+)\s+G\(\d\)',line)
            if match and len(cell) < 3:
                cell.append([float(match.group(ii)) for ii in (1,2,3)])
            if 'znucl, zion, pspdat' in line:
                zion.append(float(line.split()[1]))
    return(variables,cell,zion,foundstep)


def read_block(ff,natom):#The natom lines of three numbers following a title line ; None when the file ends before
    try:
        block = [[float(x) for x in ff.readline().split()] for _ in range(natom)]
        return np.array(block,dtype=float).reshape(natom,3)
    except ValueError:
        return None


def read_abinit(FileName,InitialStep,SkipStep):
    print ('Reading ABINIT output file ',FileName,' from initial step ',InitialStep,' writing every ',SkipStep)
    (variables,cell,zion,foundstep) = read_input(FileName)
    if 'ndtset' in variables:
        print('Error. This output has several datasets (ndtset) : give a run of one dataset.')
        sys.exit(1)
    if not foundstep or 'natom' not in variables or 'dtion' not in variables or len(cell) < 3:
        print('Error. No "MOLDYN STEP NUMBER" line or no input variables in ',FileName,' : not a molecular dynamics run in the layout of ABINIT 6.4.1,')
        print('the only one checked so far. Newer versions of ABINIT print the MD steps differently : send an example, so that the parser can be checked on it.')
        sys.exit(1)
    if int(variables.get('optcell',['0'])[0]) != 0:
        print('Error. Variable cell (optcell ',variables['optcell'][0],') : not handled yet, no such run was available to check the output of the cell.')
        sys.exit(1)
    ionmov = int(variables['ionmov'][0])
    if ionmov not in CHECKED_IONMOV:
        print('WARNING : ionmov ',ionmov,' was not checked (only ionmov ',CHECKED_IONMOV,'). Check that the velocities belong to the positions of each snapshot.')

    MyCrystal = cr.Lattice()
    natom = int(variables['natom'][0])
    typat_in = np.array([int(x)-1 for x in variables['typat']])
    order = np.argsort(typat_in,kind='stable')          #the atoms of the UMD file are grouped by element ; ABINIT keeps the order of the input
    MyCrystal.natom = natom
    MyCrystal.ntypat = int(variables['ntypat'][0])
    MyCrystal.types = [int(np.sum(typat_in == ii)) for ii in range(MyCrystal.ntypat)]
    MyCrystal.typat = typat_in[order].tolist()
    MyCrystal.elements = [cr.Elements2rest(int(round(float(x))))[1] for x in variables['znucl'][:MyCrystal.ntypat]]
    MyCrystal.masses = [float(x) for x in variables['amu'][:MyCrystal.ntypat]]
    MyCrystal.zelec = zion[:MyCrystal.ntypat]
    MyCrystal.noelectrons = sum(zz*nn for zz,nn in zip(MyCrystal.zelec,MyCrystal.types)) - float(variables.get('charge',['0.0'])[0])
    MyCrystal.atoms = [cr.Atom() for _ in range(natom)]
    MyCrystal.rprimd = (np.array(cell)*BOHR).tolist()
    MyCrystal.acell = np.linalg.norm(MyCrystal.rprimd,axis=1).tolist()
    MyCrystal.gprimd = MyCrystal.makegprimd()
    mass = np.array(MyCrystal.masses)[MyCrystal.typat]  #ABINIT's masses, getdensity below sets those of the table
    volume = MyCrystal.makevolume()
    rprimd = np.array(MyCrystal.rprimd)
    gprimd = np.array(MyCrystal.gprimd)
    dtion = float(variables['dtion'][0])                #atomic units
    TimeStep = round(dtion*AUTIME,6)                    #fs
    print ('Total number of atoms = ',natom,' elements ',MyCrystal.elements,' numbers ',MyCrystal.types,' time step ',TimeStep,' fs')
    open(FileName + '.umd.dat','w').close()
    umdpf.print_header(FileName,MyCrystal,'ABINIT : the UMD file originated from an ABINIT ab initio molecular dynamics simulation')
    MyCrystal.density = MyCrystal.getdensity()

    xcart = xnext = forces = vels = sigma = oldx = diffcoords = None
    magnetization = 0.0
    hasmagn = False
    written = 0
    lastwritten = -1
    unconverged = 0
    nostress = 0
    sumT = 0.0
    nT = 0
    target = []
    entropy = None
    with open(FileName,'r') as ff:
        while True:
            line = ff.readline()
            if not line: break
            entry = line.split()
            if len(entry) == 0:
                continue
            if line.startswith(' MOLDYN STEP NUMBER'):
                xcart = xnext = forces = vels = sigma = None
            elif entry[0] == 'iter' and 'Etot(hartree)' in entry:
                hasmagn = entry[-1] == 'magn'
            elif entry[0] == 'ETOT' and hasmagn:
                magnetization = float(entry[-1])
            elif 'was not enough SCF cycles to converge' in line:
                unconverged += 1
            elif entry[0] == '-kT*entropy':
                entropy = -float(entry[-1])*HARTREE
            elif 'Cartesian components of stress tensor (hartree/bohr^3)' in line and xcart is None:
                #the stress comes before the first block of the step ; sigma(1 1) sigma(3 2) / sigma(2 2) sigma(3 1) / sigma(3 3) sigma(2 1)
                rows = [ff.readline().split('=') for _ in range(3)]
                if all(len(row) == 3 for row in rows):
                    sigma = [float(row[1].split()[0]) for row in rows] + [float(row[2]) for row in rows]
            elif 'Cartesian coordinates (xcart)' in line:
                #each step prints two blocks, before and after the positions are moved : the positions of the second are those of the next step
                if xcart is None:
                    xcart = read_block(ff,natom)
                elif xnext is None:
                    xnext = read_block(ff,natom)
            elif 'Cartesian forces (fcart)' in line and forces is None:
                forces = read_block(ff,natom)
            elif 'Cartesian velocities (vel)' in line and vels is None:
                #those of the first block are centered on its positions (they follow 1/2 (forward + backward differences)) ;
                #those of the second, and the EKIN printed at the end of the step, belong to the next step
                vels = read_block(ff,natom)
            elif 'At the end of Moldyn step' in line:
                istep = int(line.split('step')[1].split(',')[0])
                potential = float(line.split('POT.En.=')[1].split()[0])*HARTREE
                ff.readline()
                temperatureline = ff.readline().split('=')
                if len(temperatureline) == 2:
                    target.append(float(temperatureline[1].split()[0]))
                if istep < InitialStep or xcart is None or forces is None or vels is None:
                    continue
                #the unwrapping is done at every MD step, also those not written (-s SkipStep) ; nearest image in reduced coordinates, for any cell
                x = xcart[order]*BOHR
                if oldx is None:
                    diffcoords = x.copy()
                else:
                    dcart = x - oldx
                    dcart -= np.rint(dcart @ gprimd.T) @ rprimd
                    diffcoords = diffcoords + dcart
                oldx = x
                #velocities in atomic units = bohr / atomic time ; the first step has none (they are initialized after its first block): forward difference
                if vels.any():
                    v = vels[order]*BOHR/AUTIME
                elif xnext is not None:
                    v = (xnext[order]-xcart[order])*BOHR/(dtion*AUTIME)
                else:
                    v = vels[order]
                kinetic = 0.5*AMU_A2FS2*np.sum(mass[:,None]*v*v)
                sumT += 2.0*kinetic/(3.0*natom*KB)
                nT += 1
                if (istep - InitialStep) % SkipStep != 0:
                    continue
                xred = x @ gprimd.T
                xred -= np.floor(xred)
                wrapped = xred @ rprimd
                f = forces[order]*HARTREE/BOHR
                for iatom in range(natom):
                    atom = MyCrystal.atoms[iatom]
                    atom.xred = xred[iatom].tolist()
                    atom.xcart = wrapped[iatom].tolist()
                    atom.vels = v[iatom].tolist()
                    atom.forces = f[iatom].tolist()
                MyCrystal.toten = potential                         #POT.En. = the free energy Etotal
                MyCrystal.electronicentropy = 0.0                   #ABINIT prints -kT*entropy at the end of the run only
                MyCrystal.kineticenergy = kinetic
                MyCrystal.internalenergy = potential + kinetic
                MyCrystal.energywithdrift = MyCrystal.internalenergy
                #3N degrees of freedom, as ABINIT (TKA of the end of each step)
                MyCrystal.temperature = 2.0*kinetic/(3.0*natom*KB)
                MyCrystal.magnetization = magnetization
                if sigma is None:
                    nostress += 1
                    MyCrystal.stress = [0.0 for _ in range(6)]
                    MyCrystal.pressure = 0.0
                else:
                    #UMD = positive in compression with the kinetic part, as the Total+kin. of VASP ; ABINIT's sigma = positive in tension, without the kinetic part
                    kin = AMU_A2FS2*EVA3_GPA/volume*np.array([np.sum(mass*v[:,ii]*v[:,jj]) for ii,jj in ((0,0),(1,1),(2,2),(0,1),(1,2),(2,0))])
                    #sigma = xx yy zz yz zx xy
                    MyCrystal.stress = (-HAB3_GPA*np.array([sigma[0],sigma[1],sigma[2],sigma[5],sigma[3],sigma[4]]) + kin).tolist()
                    MyCrystal.pressure = sum(MyCrystal.stress[:3])/3.0
                MyCrystal.enthalpy = MyCrystal.internalenergy + MyCrystal.pressure*volume/EVA3_GPA
                umdpf.print_snapshots(FileName,MyCrystal,TimeStep*SkipStep,(istep-InitialStep)*TimeStep,diffcoords)
                written += 1
                lastwritten = istep
    if written == 0:
        print('WARNING : no snapshot written : the run has no complete MD step from step ',InitialStep)
        return(0.0,TimeStep)
    if entropy is not None:
        print('ABINIT prints the electronic entropy only at the end of the run (T S = ',round(entropy,3),' eV in the last step) : ElectronicEntropy is 0 in every snapshot, and InternalEnergy = free energy + kinetic energy')
    if unconverged > 0:
        print('WARNING : ',unconverged,' MD steps ended with their SCF cycle not converged')
    if nostress > 0:
        print('WARNING : ',nostress,' snapshots have no stress tensor in the output : StressTensor and Pressure 0 there')
    meanT = sumT/nT
    print('Mean kinetic temperature ',round(meanT,1),' K')
    if len(target) > 0 and (meanT < 0.9*min(target) or meanT > 1.1*max(target)):
        print('WARNING : the run asked for ',min(target),' to ',max(target),' K (mditemp, mdftemp) : the kinetic temperature is far from it, the thermostat does not hold the run there')
    return((lastwritten-InitialStep)*TimeStep,TimeStep)


def main(argv):
    FileName = ''
    InitialStep = 0
    SkipStep = 1
    umdpf.headerumd()
    usage = 'AbinitParser.py -f <ABINIT_output_filename> -i <InitialStep> -s <SkipStep>'
    try:
        opts, arg = getopt.getopt(argv,"hf:i:s:",["fABINITfile","iInitialStep","sSkipStep"])
    except getopt.GetoptError:
        print (usage)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('AbinitParser.py program to translate ABINIT molecular dynamics output files into UMD ascii format file, <file>.umd.dat')
            print('Checked on ABINIT 6.4.1, ionmov 9 (Langevin), constant cell.')
            print (usage)
            print ('   default values: -i 0 -s 1')
            sys.exit()
        elif opt in ("-f", "--fABINITfile"):
            FileName = str(arg)
        elif opt in ("-i", "--iInitialStep"):
            InitialStep = int(arg)
        elif opt in ("-s", "--iSkipStep"):
            SkipStep = int(arg)
    if (os.path.isfile(FileName)):
        (CurrentTime,TimeStep) = read_abinit(FileName,InitialStep,SkipStep)
        print('Total simulation time ' + str(CurrentTime) + ' fs, done in ' + str(int(round(CurrentTime/TimeStep))) + ' steps of ' + str(TimeStep) + ' fs each.')
    else:
        print ('Error. The ABINIT output file ',FileName,' does not exist : give it with -f <file>')
        sys.exit(1)


if __name__ == "__main__":
   main(sys.argv[1:])
