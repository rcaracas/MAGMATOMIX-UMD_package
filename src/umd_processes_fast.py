#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 12:22:44 2023

@author: Kevin Jiguet-Covex
"""


import sys, os.path, re
import crystallography as cr
from functools import partial
import concurrent.futures
import ctypes
import numpy as np
import time
import platform

LibraryName = {'Linux':'c_UMDprocess.so', 'Windows':'c_UMDprocess.dll', 'Darwin':'c_UMDprocess.dylib'}.get(platform.system(),'c_UMDprocess.so')
read_lib = ctypes.cdll.LoadLibrary(os.path.join(os.path.dirname(os.path.abspath(__file__)),LibraryName))
read_lib.defineBonds.argtypes = [ctypes.POINTER(ctypes.c_char),ctypes.c_int,ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_int),ctypes.c_int,ctypes.c_int,ctypes.c_int]
read_lib.defineBonds.restype = ctypes.POINTER(ctypes.c_int)
read_lib.free_memory.argtypes = [ctypes.c_void_p]
read_lib.free_memory.restype = None
read_lib.umd_count.argtypes = [ctypes.c_char_p,ctypes.c_long,ctypes.POINTER(ctypes.c_long)]
read_lib.umd_count.restype = ctypes.c_long
read_lib.umd_read.argtypes = [ctypes.c_char_p,ctypes.c_long,ctypes.c_long,ctypes.c_long,ctypes.c_long,ctypes.c_int,ctypes.c_int,ctypes.c_void_p,ctypes.c_void_p]
read_lib.umd_read.restype = ctypes.c_long

COLUMNS = {'xred':(0,3), 'xcart':(3,3), 'absxcart':(6,3), 'velocity':(9,3), 'everything':(0,12)}   #first column and number of columns of the atom lines


def read_snapshotbonds_C(octettop,octetbot,CentIndexes,AdjIndexes,File,natom):#Converts the bonding information of a given snapshot into a list of atoms 
                                                                               #the data is contained between the octets "octettop" and "octetbot" in the bonding file
    with open(File,"rb") as ff:
        ff.seek(octettop,0)
        snapshot = ff.read(octetbot-octettop)
    CIp = (ctypes.c_int * len(CentIndexes))(*CentIndexes)
    AIp = (ctypes.c_int * len(AdjIndexes))(*AdjIndexes)
    BList = read_lib.defineBonds(snapshot,len(snapshot),CIp,AIp,len(CentIndexes)//2,len(AdjIndexes)//2,natom)
    SnapshotBondingTable = BList[2:BList[0]+1]                          #the bound atoms
    SnapshotBondIndexes = BList[BList[0]+2:BList[0]+BList[1]+2]         #where the bound atoms of each atom begin in SnapshotBondingTable
    read_lib.free_memory(BList)
    return SnapshotBondingTable,SnapshotBondIndexes

def Octets_from_File(File,key,nlines):#Byte offsets of the snapshots of a file : after each line beginning with key, and the end of the file (nlines
    OctIndexes = []                   #lines, the atoms, follow each key line)
    print("Preparing the file ",File," to be read...")
    with open(File,"rb") as ff:
        key = key.encode()
        while True :
            line = ff.readline()
            if not line :
                break
            line = line.split()
            if len(line)>0 and line[0]==key:
                OctIndexes.append(ff.tell())
                for _ in range(nlines):
                    ff.readline()
        OctIndexes.append(ff.tell())
    return OctIndexes

def Crystallization(File):#The header of a UMD file (or of a bonding file) and the cell of its first snapshot : returns a Lattice (natom, ntypat, types,
    #elements, masses, zelec, typat, noelectrons, lambda_ThermoInt, origin, acell, rprim, rprimd) and the timestep (fs ; 0.0 when absent)
    MyCrystal = cr.Lattice()
    TimeStep = 0.0
    with open(File,"r") as ff:
        for line in ff:
            entry = line.split()
            if len(entry) == 0:
                continue
            if entry[0] == 'natom':
                MyCrystal.natom = int(entry[1])
                MyCrystal.typat = [0 for _ in range(MyCrystal.natom)]
            elif entry[0] == 'ntypat':
                MyCrystal.ntypat = int(entry[1])
                MyCrystal.types = [0 for _ in range(MyCrystal.ntypat)]
                MyCrystal.elements = ['X' for _ in range(MyCrystal.ntypat)]
                MyCrystal.masses = [1.0 for _ in range(MyCrystal.ntypat)]
                MyCrystal.zelec = [1.0 for _ in range(MyCrystal.ntypat)]
            elif entry[0] == 'no.of.electrons':
                MyCrystal.noelectrons = float(entry[1])
            elif entry[0] == 'types':
                MyCrystal.types = [int(x) for x in entry[1:MyCrystal.ntypat+1]]
            elif entry[0] == 'elements':
                MyCrystal.elements = entry[1:MyCrystal.ntypat+1]
            elif entry[0] == 'masses':
                MyCrystal.masses = [float(x) for x in entry[1:MyCrystal.ntypat+1]]
            elif entry[0] == 'Zelectrons':
                MyCrystal.zelec = [float(x) for x in entry[1:MyCrystal.ntypat+1]]
            elif entry[0] == 'typat':
                MyCrystal.typat = [int(x) for x in entry[1:MyCrystal.natom+1]]
            elif entry[0] == 'lambda_ThermoInt':
                MyCrystal.lambda_ThermoInt = float(entry[1])
            elif entry[0] == 'origin':
                MyCrystal.origin = ' '.join(entry[1:])
            elif entry[0] == 'first_snapshot':
                MyCrystal.firstsnapshot = int(entry[1])
            elif entry[0] == 'timestep':
                TimeStep = float(entry[1])
            elif entry[0] == 'acell':
                MyCrystal.acell = [float(x) for x in entry[1:4]]
            elif entry[0] in ('rprim_a','rprim_b','rprim_c'):
                MyCrystal.rprim['abc'.index(entry[0][-1])] = [float(x) for x in entry[1:4]]
            elif entry[0] in ('rprimd_a','rprimd_b','rprimd_c'):
                MyCrystal.rprimd['abc'.index(entry[0][-1])] = [float(x) for x in entry[1:4]]
            elif entry[0] in ('atoms:','step'):                         #first snapshot of a UMD file, of a bonding file
                break
    return MyCrystal, TimeStep

def count_snapshots(UMDfile,natom):#Number of complete snapshots of a UMD file (a last snapshot cut short, by a run still being written, is left out)
    total = ctypes.c_long(0)
    nsnap = read_lib.umd_count(os.fsencode(UMDfile),natom,ctypes.byref(total))
    if nsnap < 0:
        print('Cannot open',UMDfile)
        sys.exit()
    if total.value > nsnap:
        print('WARNING : the last snapshot of',UMDfile,'is incomplete ; it is left out')
    return nsnap

def read_atoms(UMDfile,natom,first,last,sampling,col0,ncol,cells=False):#Columns col0 to col0+ncol-1 of the atom lines of the snapshots first,
    #first+sampling, ... before last (from 0), in one pass in C : array (snapshot, atom, column), and with cells the lattice vectors (rows a, b, c) of
    #each of these snapshots (rprimd, else acell, else NaN)
    nsel = max(0,-(-(last-first)//sampling))
    values = np.empty((nsel,natom,ncol))
    cellvalues = np.empty((nsel,3,3)) if cells else None
    nread = read_lib.umd_read(os.fsencode(UMDfile),natom,first,last,sampling,col0,ncol,values.ctypes.data,cellvalues.ctypes.data if cells else None)
    if nread == -1:
        print('Cannot open',UMDfile)
        sys.exit()
    if nread < -1:
        print('ERROR : an atom line of the snapshot',first+sampling*(-2-nread),'of',UMDfile,'has fewer than',col0+ncol,'numbers, or a value that is not a number')
        sys.exit()
    if nread != nsel:
        print('ERROR :',nread,'snapshots read in',UMDfile,'instead of',nsel)
        sys.exit()
    return (values, cellvalues) if cells else values

def read_values(UMDfile,key,mode="line",Nsteps=1,firststep = 0,laststep = None,cutoff="all",nCores = None):#Positions or velocities of the snapshots
    #firststep, firststep+Nsteps, ... before laststep (default : the end ; cutoff : only the first cutoff snapshots of the file), of the key
    #xred, xcart, absxcart (unwrapped), velocity, or everything (the four, 12 columns) : returns the header (Crystallization), an array of the values of
    #each snapshot (mode line : atom after atom ; lists : (snapshot, atom, column) ; chunks : column after column), the timestep, and the number of
    #snapshots from firststep to laststep. One pass in C (nCores is not used, kept for the callers)
    t = time.time()
    if key not in COLUMNS:
        print("Parameter < ",key," > not recognized")
        sys.exit()
    if mode not in ("line","lists","chunks"):
        print("Parameter < ",mode," > not recognized")
        sys.exit()
    MyCrystal,TimeStep = Crystallization(UMDfile)
    nsnap = count_snapshots(UMDfile,MyCrystal.natom)
    if cutoff != "all":
        nsnap = min(nsnap,cutoff)
    first = min(max(0,firststep),nsnap)
    last = nsnap if laststep is None else min(max(laststep,first),nsnap)
    col0,ncol = COLUMNS[key]
    print('Extracting the',key,'values of the snapshots',first,'to',last-1,'every',Nsteps,'from the file',UMDfile)
    values = read_atoms(UMDfile,MyCrystal.natom,first,last,Nsteps,col0,ncol)
    if mode == "line":
        values = values.reshape(len(values),-1)
    elif mode == "chunks":
        values = values.transpose(0,2,1).reshape(len(values),-1)
    print(len(values),'snapshots read in',round(time.time()-t,2),'s')
    return MyCrystal, values, TimeStep, last-first

def read_cells(UMDfile,nsnap,Nsteps=1,firststep=0):#Lattice vectors (rows a, b, c) of the nsnap snapshots read by read_values with the same Nsteps and
    #firststep : the cell of each snapshot (NPT), from its rprimd lines, else its acell
    MyCrystal,TimeStep = Crystallization(UMDfile)
    values, cells = read_atoms(UMDfile,MyCrystal.natom,firststep,firststep+Nsteps*nsnap,Nsteps,0,0,cells=True)
    if np.isnan(cells).any():
        print('ERROR : snapshots of',UMDfile,'without rprimd or acell lines')
        sys.exit()
    return cells

def read_snapshots(UMDfile,natom,first=0,sampling=1,fields=None,atoms=True):#Line by line, any file size : yields (snapshot, fields, atoms) for every
    #snapshot ; fields = {name : list of the numbers of its line} (the names of fields, default all), atoms = array natom x columns (reduced*3,
    #cartesian*3, unwrapped*3, velocity*3, force*3, charge, magnetization) of every sampling-th snapshot from first on, None for the others
    with open(UMDfile,'r') as f:
        values = None
        isnap = -1
        for line in f:
            word = line.split(' ',1)[0].strip()
            if word == 'timestep':
                isnap += 1
                values = {}
            if values is None:
                continue
            if word == 'atoms:':
                block = [next(f) for _ in range(natom)]
                array = None
                if atoms and isnap >= first and (isnap-first) % sampling == 0:
                    array = np.array(' '.join(block).split(),dtype=float).reshape(natom,-1)
                yield isnap, values, array
                values = None
            elif word and (fields is None or word in fields):
                numbers = []
                for x in line.split()[1:]:
                    try:
                        numbers.append(float(x))
                    except ValueError:
                        break
                values[word] = numbers


def cell_of(fields):#Lattice vectors (rows a, b, c, A) from the rprimd lines of a snapshot read by read_snapshots, else from acell (orthogonal), else None
    if all(k in fields for k in ('rprimd_a','rprimd_b','rprimd_c')):
        return np.array([fields['rprimd_a'][:3],fields['rprimd_b'][:3],fields['rprimd_c'][:3]])
    if len(fields.get('acell',[])) >= 3:
        return np.diag(fields['acell'][:3])
    return None


def xyz_frame(symbols,xcart,cell=None,info='',columns=()):#One frame of an extended XYZ file (plain XYZ readers see symbol x y z) : the number of atoms,
    #the comment line (Lattice of the cell, Properties of the columns, pbc, then info, words key=value), one line per atom : symbol, position (A),
    #then the columns, a list of (name, array natom x k)
    props = 'species:S:1:pos:R:3' + ''.join(':'+name+':R:'+str(np.shape(c)[1]) for name,c in columns)
    comment = 'Properties=' + props
    if cell is not None:
        comment = 'Lattice="' + ' '.join('%.8f' % x for x in np.ravel(cell)) + '" ' + comment + ' pbc="T T T"'
    if info:
        comment += ' ' + info
    lines = [str(len(symbols)),comment]
    data = np.hstack([np.asarray(xcart,dtype=float)]+[np.asarray(c,dtype=float) for name,c in columns])
    for el,row in zip(symbols,data):
        lines.append(el+' '+' '.join('%.6f' % x for x in row))
    return '\n'.join(lines)+'\n'

def read_bonds(BondFile,centEls,adjEls,Nsteps=1,nCores = None):#Bonds of every Nsteps-th snapshot of a bonding file, between the central and adjacent elements
    
    MyCrystal,TimeStep = Crystallization(BondFile)

    if centEls == ['all']:
        centEls = MyCrystal.elements
    if adjEls == ['all']:
        adjEls = MyCrystal.elements

    CentElWitness = [False for _ in range(len(centEls))]
    AdjElWitness = [False for _ in range(len(adjEls))]

        
    AdjacentElements = []
    CentralElements = []


    for i in range(len(MyCrystal.elements)):
        el = MyCrystal.elements[i]
        trueEl = ""
        for ch in el :
            if not ch.isdigit():
                trueEl+=ch
        if trueEl in centEls :
            CentralElements.append(el)
            CentElWitness[centEls.index(trueEl)]=True

        elif el in centEls :
            CentralElements.append(el)
            CentElWitness[centEls.index(el)]=True

        if trueEl in adjEls:
            AdjacentElements.append(el)
            AdjElWitness[adjEls.index(trueEl)]=True
        elif el in adjEls:
            AdjacentElements.append(el)
            AdjElWitness[adjEls.index(el)]=True
            
    CentIndexes = []
    AdjIndexes = []
    
    for el in CentralElements :
        indCent = MyCrystal.elements.index(el)
        CentMin = sum(MyCrystal.types[:indCent])
        CentMax = CentMin + MyCrystal.types[indCent] -1
        CentIndexes+=[CentMin,CentMax]
        

    for el in AdjacentElements :
        indAdj = MyCrystal.elements.index(el)
        AdjMin = sum(MyCrystal.types[:indAdj])
        AdjMax = AdjMin + MyCrystal.types[indAdj] -1
        AdjIndexes+=[AdjMin,AdjMax]
        
   
#    if centEls == ["all"]:
#        CentIndexes = [0,MyCrystal.natom-1]
#        CentElWitness = [True]
#    if adjEls == ["all"]:
#        AdjIndexes = [0,MyCrystal.natom-1]
#        AdjElWitness = [True]
        
    noAt = []
    for wit in range(len(CentElWitness)) :
        if CentElWitness[wit] == False :
            noAt.append(centEls[wit])

    for wit in range(len(AdjElWitness)) :
        if AdjElWitness[wit] == False :
            noAt.append(adjEls[wit])

    if noAt != []:
        if len(noAt) ==1 :
            print("ERROR : element ",noAt[0]," not present in simulation")
        else :
            print("ERROR : elements ",noAt," not presents in simulation")
        return -1,noAt,0,[0,0],0
        

#    print(CentMin,CentMax,AdjMin,AdjMax,indCent,indAdj)

    OctIndexes = Octets_from_File(BondFile,"step",MyCrystal.natom)

    OctTop = OctIndexes[0:-1:Nsteps]                                   #beginning of the snapshots 0, Nsteps, 2 Nsteps, ...
    OctBot = OctIndexes[1::Nsteps]                                      #and their end

    print("Extracting bonds from file ",BondFile)
    
    bondsRed = partial(read_snapshotbonds_C,CentIndexes=CentIndexes,AdjIndexes=AdjIndexes,File=BondFile,natom=MyCrystal.natom)

#    Data=[]
#    for i in range(len(OctTop)):
#        Data.append(bondsRed(OctTop[i],OctBot[i]))            
    if nCores != None :
        with concurrent.futures.ProcessPoolExecutor(max_workers = nCores) as executor :
            Data = list(executor.map(bondsRed,OctTop,OctBot))
    else:
        with concurrent.futures.ProcessPoolExecutor() as executor :
            Data = list(executor.map(bondsRed,OctTop,OctBot))
            
    Bonds = [D[0] for D in Data]
    BondsIndexes = [D[1] for D in Data]
    return CentIndexes, AdjIndexes, MyCrystal, [Bonds, BondsIndexes], TimeStep*Nsteps#TimeStep*Nsteps : time between two snapshots used

def read_stresses_4visc(UMDfile):#The header (Crystallization), one Lattice per snapshot (stress, temperature, rprimd, cellvolume) and the timestep
    MyCrystal,TimeStep = Crystallization(UMDfile)
    AllSnapshots = []
    for isnap,fields,atoms in read_snapshots(UMDfile,MyCrystal.natom,fields={'timestep','StressTensor','Temperature','acell','rprimd_a','rprimd_b','rprimd_c'},atoms=False):
        SnapshotCrystal = cr.Lattice()
        TimeStep = fields['timestep'][0]
        if 'StressTensor' in fields:
            SnapshotCrystal.stress = fields['StressTensor'][0:6]
        if 'Temperature' in fields:
            SnapshotCrystal.temperature = fields['Temperature'][0]
        cell = cell_of(fields)
        if cell is not None:
            SnapshotCrystal.rprimd = cell.tolist()
        SnapshotCrystal.cellvolume = SnapshotCrystal.makevolume()
        AllSnapshots.append(SnapshotCrystal)
    return MyCrystal, AllSnapshots, TimeStep

def headerumd():
    print ('\n * * * * * ')
    print ('    UMD package for analyzing Molecular Dynamics simulations, version 2.9')
    print ('distributed under GNU General Public License 3 and CeCILL v2')
    print ('https://github.com/rcaracas/MAGMATOMIX-UMD_package')
    print ('please cite as: Caracas, R. et al. Analyzing Melts and Fluids from Ab Initio Molecular Dynamics Simulations')
    print ('    with the UMD Package. J. Vis. Exp. (175), e61534, doi:10.3791/61534 (2021)')
    print (' ')

def former_outputs(folder,pattern):#The files of folder whose whole name matches the regular expression pattern (the names of the outputs of a script)
    return [os.path.join(folder,f) for f in sorted(os.listdir(folder or os.curdir)) if re.fullmatch(pattern,f)]

#The outputs of the scripts of the package, used by umd_clean.py : regular expressions of the whole file names. RUN_OUTPUTS follow the
#base <base> of a UMD file <base>.umd.dat of the same folder ; FOLDER_OUTPUTS are written once per folder. The UMD files themselves (parsers,
#build_umd.py, --extract, repair_umd.py) are never outputs here. Keep this list in step with the scripts.
SELECTION = r'(?:\.(?:all|(?:[A-Z][a-z]?)+)-(?:all|(?:[A-Z][a-z]?)+))?'          #gofr_umd.py -c/-a : <base>.<central>-<adjacent>
RUN_OUTPUTS = [r'\.check\.dat', r'\.overlap\.dat',                                                        #check_umd.py
               SELECTION+r'\.gofr\.dat',                                                                  #gofr_umd.py
               SELECTION+r'\.(?:maximum_bonds|average_bonds|coordinations)\.inp', SELECTION+r'\.gofr\.[A-Z][a-z]?-[A-Z][a-z]?\.png',#analyze_gofr.py
               r'\.bonding\.dat',                                                                         #bonding_umd.py
               r'\.r\d+\.(?:popul|stat|step)\.dat', r'\.angles\.dat', r'\.bondslives\.dat',               #speciation_and_angles.py
               r'\.r\d+\.(?:liquid|vapor)\.(?:popul|stat)\.dat',                                           #liquid_and_vapor_speciation.py
               r'\.r\d+(?:\.(?:liquid|vapor))?(?:\.(?:liquid|vapor))?\.lifetimes\.pdf',                    #lifetimes.py
               r'\.r\d+(?:\.(?:liquid|vapor))?__[A-Za-z0-9_]+\.msd\.dat',                                   #msd_clusters_umd.py
               r'\.msd\.dat', r'\.msd\.blocks\.dat', r'\.diffusion\.dat', r'\.diffusion\.tensor\.dat',       #msd_umd.py
               r'\.stress\.scf\.dat', r'\.visc\.dat',                                                      #viscosity_umd.py
               r'\.vels\.scf\.dat', r'\.vibr\.dat', r'\.vibr\.axes\.dat', r'\.vibr\.err\.dat', r'\.vacf\.diffusion\.dat',#vibrations_umd.py
               r'\.vdos\.thermo\.dat', r'\.2pt\.dat',                                                      #vdos_thermo.py
               r'(?:\.\d+)?\.xyz', r'(?:\.\d+)?\.lammpstrj', r'\._\d+_\.POSCAR\.vasp', r'\._\d+_\.lammps\.data']#umd2out.py
RUN_FOLDERS = [r'(?:\.\d+)?\.deepmd']                                    #umd2out.py : the files it writes inside (DEEPMD_FILES)
DEEPMD_FILES = [r'type\.raw', r'type_map\.raw', r'set\.\d{3}/(?:box|coord|energy|force|virial)\.npy']
FOLDER_OUTPUTS = [r'gofrs_summary\.txt',                                   #analyze_gofr.py
                  r'diffusivities\.txt', r'diffusivities_clusters\.txt',  #analyze_msd.py, analyze_msd_clusters.py
                  r'fullthermo\.dat',                                      #averages.py -a
                  r'thermo_integration\.(?:dat|pdf)',                      #thermo_integration_umd.py
                  r'lifetimes\.merged(?:\.(?:liquid|vapor))?\.pdf',        #lifetimes.py -m (default name)
                  r'.*all\.-t.*\.-L.*\.-R\d+\.-C=.*\.dat',                 #stat-concentrate.py
                  r'umd_workflow\.log', r'umd_workflow\.done', r'umd_workflow\.summary\.dat']#umd_workflow.py

def umd_outputs(folder):#The outputs of the package in folder (paths relative to it), from the UMD files there : never a UMD file, never a file
    names = sorted(os.listdir(folder or os.curdir))                      #whose name is not that of an output
    bases = [n[:-8] for n in names if n.endswith('.umd.dat')]
    builds = [b for b in bases if re.fullmatch(r'.*\.build\d+',b)]     #build_umd.py writes <base>.build<k>.xyz with its UMD file : kept with it
    found = []
    for n in names:
        path = os.path.join(folder,n)
        if n.endswith('.umd.dat') or n in [b+'.xyz' for b in builds]:
            continue
        if os.path.isdir(path) and not os.path.islink(path):
            if any(re.fullmatch(re.escape(b)+pattern,n) for b in bases for pattern in RUN_FOLDERS):
                for sub, dirs, files in os.walk(path):
                    for f in sorted(files):
                        rel = os.path.relpath(os.path.join(sub,f),path).replace(os.sep,'/')
                        if any(re.fullmatch(pattern,rel) for pattern in DEEPMD_FILES):
                            found.append(os.path.join(n,rel))
        elif any(re.fullmatch(pattern,n) for pattern in FOLDER_OUTPUTS) or \
             any(re.fullmatch(re.escape(b)+pattern,n) for b in bases for pattern in RUN_OUTPUTS):
            found.append(n)
    return found

def remove_former_outputs(paths):#Rule of the package : a script run again first removes the outputs of its former run, so that none is left stale
    for path in paths:            #(e.g. after a failure, or when the new run writes fewer files) ; only these named files, never anything else
        if os.path.lexists(path) and not os.path.isdir(path):
            os.remove(path)
            print('former output removed : ',path)

def write_poscar(filename,title,cell,elements,counts,xred,vels=None):#VASP5 POSCAR : cell (rows a, b, c, A), the atoms grouped by element (counts[k] atoms
    #of elements[k], in this order), direct coordinates, and Cartesian velocities (A/fs) if given ; the only POSCAR writer of the package
    with open(filename,'w') as ff:
        ff.write(title+'\n')
        ff.write('  1.0\n')
        for vec in cell:
            ff.write('  {:.10f}  {:.10f}  {:.10f}\n'.format(*vec))
        ff.write('  '+'  '.join(elements)+'\n')
        ff.write('  '+'  '.join(str(n) for n in counts)+'\n')
        ff.write('Direct\n')
        for x in xred:
            ff.write('  {:.10f}  {:.10f}  {:.10f}\n'.format(*x))
        if vels is not None:
            ff.write('Cartesian\n')
            for v in vels:
                ff.write('  {:.10e}  {:.10e}  {:.10e}\n'.format(*v))


def read_poscar(filename):#VASP5 POSCAR/CONTCAR (element line), reduced or Cartesian coordinates, optional selective dynamics and Cartesian velocities ;
    #returns the title, the cell (rows a, b, c, A), the element of each atom, the reduced coordinates and the velocities (A/fs ; None if absent)
    with open(filename,'r') as ff:
        lines = ff.readlines()
    scale = float(lines[1].split()[0])
    if scale < 0.0:
        print('Negative (volume-based) scaling factors in POSCAR files are not supported.')
        sys.exit()
    cell = np.array([[float(x) for x in lines[2+ii].split()[0:3]] for ii in range(3)])*scale
    elements = lines[5].split()
    if all(x.lstrip('-').isdigit() for x in elements):
        print('The POSCAR file',filename,'has no element-symbols line (VASP4 format). Please use a VASP5-style POSCAR, such as the ones exported by VESTA.')
        sys.exit()
    counts = [int(x) for x in lines[6].split()]
    natom = sum(counts)
    lineno = 7
    if lines[lineno].strip()[0] in 'sS':                                #selective dynamics
        lineno += 1
    mode = lines[lineno].strip()[0]
    coords = np.array([[float(x) for x in lines[lineno+1+ii].split()[0:3]] for ii in range(natom)])
    xred = coords if mode in 'dD' else (coords*scale) @ np.linalg.inv(cell)
    vels = None
    rest = lines[lineno+1+natom:]
    if len(rest) > natom:                                               #velocity block : a mode line (blank, Cartesian or Direct), then one line per atom
        if rest[0].strip()[:1] in ('d','D'):
            print('WARNING : velocities in direct coordinates in',filename,'are not read')
        else:
            try:
                vels = np.array([[float(x) for x in rest[1+ii].split()[0:3]] for ii in range(natom)])
            except (ValueError,IndexError):
                vels = None
    symbols = [el for el,n in zip(elements,counts) for _ in range(n)]
    return lines[0].strip(), cell, symbols, xred, vels


def print_header(FileName,MyCrystal,origin=None):#Appends the header of a UMD file to FileName.umd.dat ; origin : the line 'origin <origin>' at the end of
    #the header, the code that produced the run (every reader skips the header lines it does not know ; not the first line, which old readers take as
    #natom, and no field name in it, as averages.py greps the fields)
    string = 'natom ' + str(MyCrystal.natom) + '\n'
    string += 'ntypat ' + str(MyCrystal.ntypat) + '\n'
    string += 'no.of.electrons ' + str(MyCrystal.noelectrons) + '\n'
    string += 'types ' + ''.join(str(MyCrystal.types[ii]) + ' ' for ii in range(MyCrystal.ntypat)) + '\n'
    string += 'elements ' + ''.join(str(MyCrystal.elements[ii]) + ' ' for ii in range(MyCrystal.ntypat)) + '\n'
    string += 'masses ' + ''.join(str(MyCrystal.masses[ii]) + ' ' for ii in range(MyCrystal.ntypat)) + '\n'
    string += 'Zelectrons ' + ''.join(str(MyCrystal.zelec[ii]) + ' ' for ii in range(MyCrystal.ntypat)) + '\n'
    string += 'typat ' + ''.join(str(MyCrystal.typat[ii]) + ' ' for ii in range(MyCrystal.natom)) + '\n\n'
    string += 'lambda_ThermoInt ' + str(MyCrystal.lambda_ThermoInt) + '\n'
    if origin:
        string += 'origin ' + origin + '\n'
    string += '\n'
    with open(FileName + '.umd.dat','a') as nf:
        nf.write(string)

def print_snapshots(FileName,MyCrystal,TimeStep,CurrentTime,diffcoords):#Appends a snapshot to FileName.umd.dat : the fields of MyCrystal, and per atom
    #xred, xcart, the unwrapped positions diffcoords, velocities, forces (str(round(x,5))), charge, magnetization
    lines = ['timestep ' + str(TimeStep) + ' fs\n' + 'time ' + str(CurrentTime) + ' fs\n',
             'InternalEnergy ' + str(round(MyCrystal.internalenergy,6)) + ' eV\n',
             'TotalElectronicEnergy ' + str(round(MyCrystal.toten, 6)) + ' eV\n',                  #needed for thermo int
             'ElectronicEntropy ' + str(round(MyCrystal.electronicentropy,6)) + ' eV\n',
             'KineticEnergyIons ' + str(round(MyCrystal.kineticenergy,6)) + ' eV\n',
             'EnergyWithDrift ' + str(round(MyCrystal.energywithdrift,6)) + ' eV\n',
             'Enthalpy ' + str(round(MyCrystal.enthalpy,6)) + ' eV\n',
             'Magnetization ' + str(round(MyCrystal.magnetization,6)) + ' Magneton-Bohr\n',
             'Temperature ' + str(round(MyCrystal.temperature,1)) + ' K\n',
             'Pressure ' + str(round(MyCrystal.pressure,4)) + ' GPa\n',
             'Density ' + str(round(MyCrystal.density,3)) + ' g.cm-3\n',
             'StressTensor ' + ''.join(str(round(MyCrystal.stress[ii],4)) + ' ' for ii in range(6)) + ' GPa\n',
             'acell ' + str(MyCrystal.acell[0]) + ' ' + str(MyCrystal.acell[1]) + ' ' + str(MyCrystal.acell[2]) + ' A\n']
    for ii,name in enumerate(('a','b','c')):
        lines.append('rprim_' + name + ' ' + '  '.join(str(MyCrystal.rprimd[ii][jj]/MyCrystal.acell[ii]) for jj in range(3)) + '\n')
    for ii,name in enumerate(('a','b','c')):
        lines.append('rprimd_' + name + ' ' + '  '.join(str(MyCrystal.rprimd[ii][jj]) for jj in range(3)) + ' A\n')
    lines.append('atoms: reduced*3 cartesian*3(A) abs.diff.*3(A) velocity*3(A/fs) force*3(eV/A) charge(no.elec) magnetization(magneton-Bohr) \n')
    for iatom in range(MyCrystal.natom):
        atom = MyCrystal.atoms[iatom]
        lines.append(' '.join([str(round(x,5)) for x in (atom.xred[0],atom.xred[1],atom.xred[2],atom.xcart[0],atom.xcart[1],atom.xcart[2],
                                                          diffcoords[iatom][0],diffcoords[iatom][1],diffcoords[iatom][2],atom.vels[0],atom.vels[1],atom.vels[2],
                                                          atom.forces[0],atom.forces[1],atom.forces[2])]+[str(atom.charge),str(atom.magnet)]) + '\n')
    lines.append('\n')
    with open(FileName + '.umd.dat','a') as nf:
        nf.write(''.join(lines))
    return(CurrentTime,TimeStep)

def copyhead(Target,File,Nsteps,first=0):#Writes to Target the lines of File before its first "atoms:" line (header and fields of the first snapshot), the
    #timestep multiplied by Nsteps, and when first > 0 a line first_snapshot (the snapshot of File where the snapshots of Target start)
    with open(File,"r") as ff, open(Target,"w") as fa:
        for line in ff:
            l = line.split()
            if len(l)>0 :
                if l[0] == "timestep":
                    line = l[0] +" "+ str(float(l[1])*Nsteps) + " "+ l[2] +"\n"
                if l[0] == "atoms:":
                    break
            fa.write(line)
        if first:
            fa.write('first_snapshot '+str(first)+' : the snapshots of this file start at this snapshot of the UMD file\n')
