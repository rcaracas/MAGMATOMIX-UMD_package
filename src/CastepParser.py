#!/usr/bin/env python3
###
#AUTHORS: RAZVAN CARACAS
###
#Converts the trajectory file of a CASTEP molecular dynamics run (<seed>.md) into a UMD file. The .md file holds, for every frame, the time, the
#energies (<-- E), the temperature (<-- T), the pressure (<-- P), the cell (<-- h), its velocity (<-- hv), the stress (<-- S) and per atom the
#position (<-- R), the velocity (<-- V) and the force (<-- F), all in atomic units (Hartree, bohr, hbar/Hartree) ; converted to A, fs, eV, GPa.
#The masses are not in the file : they are found from its kinetic energies (1/2 sum m v^2 of the velocities, exact) ; else -m, or the table of elements.

import sys, getopt, os.path, gzip
import numpy as np
import crystallography as cr
import umd_processes_fast as umdpf
from vibrations_umd import KB, AMU_A2FS2
from CP2KParser import element_of, write_snapshot

BOHR = 0.529177210903               #bohr in A
HARTREE = 27.211386245988           #Hartree in eV
AUTIME = 2.418884326505E-02         #atomic unit of time in fs
EVA3_GPA = 160.2176634              #eV/A^3 in GPa
HAB3_GPA = HARTREE/BOHR**3*EVA3_GPA #Hartree/bohr^3 in GPa
AMU_ME = 1822.888486209             #atomic mass unit in electron masses
TAGS = ('E','T','P','h','hv','S','R','V','F')
FITFRAMES = 200                     #frames used to find the masses from the kinetic energies
USAGE = 'CastepParser.py -f <seed.md> -m <masses> -k <species> -o <output> -i <InitialStep> -s <SkipStep>'


def open_md(name):#the .md file, also compressed (.md.gz)
    return gzip.open(name,'rt') if name.endswith('.gz') else open(name,'r')


def stop(*message):
    print(*message)
    sys.exit()


#****The .md file

def scan_md(name):#First pass : for every frame its time (atomic units), the line of the file where it begins and the number of lines of each kind
    #(<-- E, <-- R ...) ; the species and number of every atom (R lines of the first frame) ; for the first FITFRAMES frames the kinetic energy (Hartree)
    #and the sum of v^2 of each species (atomic units), from which the masses are found ; the kinds of lines that are not known
    frames = []
    labels = []
    fit = []
    unknown = {}
    inheader = False
    stray = None                                                        #a line that is not understood : an error, but for the last line of the file
    with open_md(name) as f:
        for n, line in enumerate(f,1):
            if stray is not None and line.strip():
                stop('Line',stray[0],'of',name,'is neither a time nor a line of values (<--) :',stray[1],'\nNot a .md file of CASTEP ?')
            pos = line.rfind('<--')
            if pos >= 0:
                if not frames:
                    stop('Line',n,'of',name,': a line of values (<--) before the time of a frame. Not a .md file of CASTEP ?')
                tag = line[pos+3:].strip()
                counts = frames[-1][2]
                counts[tag] = counts.get(tag,0)+1
                if tag not in TAGS:
                    unknown[tag] = unknown.get(tag,0)+1
                elif tag == 'R' and len(frames) == 1:
                    labels.append(tuple(line.split()[:2]))
                elif tag in ('E','V') and len(frames) <= FITFRAMES:
                    tokens = line[:pos].split()
                    try:
                        if tag == 'E':
                            fit[-1][0] = float(tokens[2])
                        else:
                            fit[-1][1][tokens[0]] = fit[-1][1].get(tokens[0],0.0)+sum(float(x)**2 for x in tokens[2:5])
                    except (ValueError,IndexError):
                        pass                                            #a line cut short : the frame is found incomplete below
                continue
            word = line.strip().lower()
            if not word:
                continue
            if word in ('begin header','end header'):                   #the header (a comment)
                inheader = word == 'begin header'
                continue
            if inheader:
                continue
            tokens = word.split()
            try:
                time = float(tokens[0]) if len(tokens) == 1 else None
            except ValueError:
                time = None
            if time is None:
                stray = (n,line.strip()[:80])
                continue
            frames.append([time,n,{}])
            if len(frames) <= FITFRAMES:
                fit.append([None,{}])
    if stray is not None:                                               #the file ends in the middle of a line (a run still being written) : its last frame is cut
        if not frames:
            stop('Line',stray[0],'of',name,'is neither a time nor a line of values (<--) :',stray[1],'\nNot a .md file of CASTEP ?')
        frames[-1][2]['cut'] = 1
    return frames, labels, fit, {tag:n for tag,n in unknown.items() if tag}


def complete(counts,natom):#a frame with all its lines : the energies, the cell, the positions, and velocities, forces and stress whole or absent
    return 'cut' not in counts and counts.get('E',0) == 1 and counts.get('h',0) == 3 and counts.get('R',0) == natom and counts.get('V',0) in (0,natom) and \
           counts.get('F',0) in (0,natom) and counts.get('S',0) in (0,3) and counts.get('T',0) <= 1 and counts.get('P',0) <= 1


def last_branch(times,tolerance):#A run continued from an earlier time appends to the same .md file : the frames after the restart point are then written
    #again. True for the frames to keep : those with no later frame at the same or at an earlier time
    keep = np.ones(len(times),dtype=bool)
    if len(times) > 1:
        later = np.minimum.accumulate(np.array(times)[::-1])[::-1]
        keep[:-1] = np.array(times)[:-1] < later[1:]-tolerance
    return keep


def read_md(name,keep):#Second pass : yields the frames to keep, as {tag : array} : E (potential, Hamiltonian, kinetic), T, P, h (rows a, b, c), S, and
    #R, V, F (natom x 3) with labels (species and number of every atom) ; T, P, S, V, F None when the frame has no such lines
    def build(rows,start):
        try:
            atoms = np.array(rows['R'])
            frame = {'labels':atoms[:,:2], 'R':atoms[:,2:5].astype(float), 'E':[float(x) for x in rows['E'][0][:3]],
                     'h':np.array(rows['h'])[:,:3].astype(float)}
            for tag in ('V','F'):
                frame[tag] = None
                if tag in rows:
                    values = np.array(rows[tag])
                    if not np.array_equal(values[:,:2],atoms[:,:2]):
                        stop('The atoms of the lines <--',tag,'are not those of the lines <-- R, in the frame beginning at line',start,'of',name)
                    frame[tag] = values[:,2:5].astype(float)
            frame['S'] = np.array(rows['S'])[:,:3].astype(float) if 'S' in rows else None
            for tag in ('T','P'):
                frame[tag] = float(rows[tag][0][0]) if tag in rows else None
            if frame['R'].shape[1] != 3 or len(frame['E']) != 3 or frame['h'].shape != (3,3):
                raise ValueError
        except (ValueError,IndexError):
            stop('The frame beginning at line',start,'of',name,'cannot be read : its lines do not have the expected numbers')
        return frame
    with open_md(name) as f:
        index = -1
        rows = None
        start = 0
        inheader = False
        for n, line in enumerate(f,1):
            pos = line.rfind('<--')
            if pos >= 0:
                if rows is not None:
                    rows.setdefault(line[pos+3:].strip(),[]).append(line[:pos].split())
                continue
            word = line.strip().lower()
            if not word:
                continue
            if word in ('begin header','end header'):
                inheader = word == 'begin header'
                continue
            if inheader:
                continue
            if rows is not None:
                yield build(rows,start)
            index += 1
            start = n
            rows = {} if index < len(keep) and keep[index] else None
        if rows is not None:
            yield build(rows,start)


def round_step(value):#The time step in fs. CASTEP converts it with its own constants, which differ from version to version by 1e-9 : a step within
    #2e-8 of a number of 6 digits is that number (1.0000000002 and 0.9999999989 fs are 1 fs) ; other steps (20 atomic units...) are kept
    rounded = float('%.6g' % value)
    return rounded if abs(rounded-value) <= 2.0E-8*abs(value) else value


def fit_masses(fit,species):#The masses of the species (amu) that give the kinetic energies of the file from its velocities, KE = 1/2 sum m v^2 : exact
    #in a .md file (1e-15), so a few frames give the masses of the run to 9 digits. {} when they cannot be found (fewer frames than species, no velocities)
    names = list(dict.fromkeys(species))
    rows = [(KE,v2) for KE,v2 in fit if KE is not None and KE > 0.0 and all(v2.get(s,0.0) > 0.0 for s in names)]
    if len(rows) < len(names):
        return {}
    A = 0.5*np.array([[v2[s] for s in names] for KE,v2 in rows])
    KE = np.array([row[0] for row in rows])
    fitted, residuals, rank, singular = np.linalg.lstsq(A,KE,rcond=None)
    if rank < len(names) or fitted.min() <= 0.0 or np.abs(A @ fitted/KE-1.0).max() > 1.0E-6:
        return {}
    return {s:float('%.8g' % (m/AMU_ME)) for s,m in zip(names,fitted)}


#****Main

def main(argv):
    umdpf.headerumd()
    mdname = ''
    outname = None
    masses_opt = {}
    kinds = {}
    InitialStep = 0
    SkipStep = 1
    try:
        opts, arg = getopt.getopt(argv,"hf:m:k:o:i:s:")
    except getopt.GetoptError:
        print(USAGE)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('CastepParser.py : converts the trajectory file of a CASTEP molecular dynamics run (<seed>.md) into a UMD file')
            print(USAGE)
            print('-f : the .md file of the run (also compressed, .md.gz). It holds every frame in atomic units : the time, and the lines')
            print('     <-- E : potential energy, Hamiltonian energy (the conserved quantity of the ensemble), kinetic energy ;')
            print('     <-- T : temperature ; <-- P : pressure ; <-- h : cell vectors a, b, c ; <-- hv : their velocities (not used) ;')
            print('     <-- S : stress tensor ; <-- R, <-- V, <-- F : position, velocity and force of every atom.')
            print('     InternalEnergy = potential + kinetic energy (as VaspParser.py), TotalElectronicEnergy = potential energy, KineticEnergyIons,')
            print('     EnergyWithDrift = Hamiltonian energy, Temperature = that of the file, ElectronicEntropy 0 (not in the file).')
            print('     StressTensor = - S : the S of CASTEP is negative in compression and holds the kinetic part, the UMD stress is positive in')
            print('     compression with the kinetic part. Pressure = that of the file (= trace / 3), Enthalpy = InternalEnergy + P V.')
            print('     A run made without calculate_stress has no P and S lines : Pressure, StressTensor and Enthalpy are then 0.')
            print('-m : masses, e.g. -m H:2.01410,O:15.9994. They are not in the .md file : by default they are those that give the kinetic')
            print('     energies of the file from its velocities (exact : the masses of the run, isotopes included) ; the table of elements when')
            print('     they cannot be found (fewer frames than species, no velocities).')
            print('-k : element of the species, when their names are not element symbols, e.g. -k D:H,Fe_up:Fe')
            print('     (default : the element symbol the name begins with : Fe:up -> Fe, H:D -> H, O1 -> O)')
            print('-o : name of the UMD file (default <.md file>.umd.dat, without .gz)')
            print('-i : first frame written (default 0) ; -s : every s-th frame (default 1)')
            print('The atoms are grouped by element, in the order of the file within an element. Unwrapped positions : those of CASTEP, which does')
            print('not wrap the atoms during a run, followed from frame to frame (also the frames not written) ; with a changing cell (NPT, NPH)')
            print('each move is the change of the reduced coordinates in the cell of the frame, so that the atoms do not move with the cell.')
            print('A run continued from an earlier time appends to the same file : of the frames written several times the last ones are kept.')
            print('Checked on NVT and NVHug runs in orthogonal cells of constant shape ; not yet on a variable or a non-orthogonal cell.')
            sys.exit()
        elif opt == '-f':
            mdname = str(arg)
        elif opt == '-m':
            try:
                masses_opt = {item.split(':')[0].strip():float(item.split(':')[1]) for item in arg.split(',')}
            except (ValueError,IndexError):
                stop('-m : masses as element:mass, e.g. -m H:2.01410,O:15.9994')
        elif opt == '-k':
            kinds = {item.rpartition(':')[0].strip():item.rpartition(':')[2].strip() for item in arg.split(',')}
            if any(cr.element_entry(el) is None for el in kinds.values()):
                stop('-k : species as species:element, e.g. -k D:H,Fe_up:Fe')
        elif opt == '-o':
            outname = str(arg)
        elif opt == '-i':
            InitialStep = int(arg)
        elif opt == '-s':
            SkipStep = int(arg)
    if not os.path.isfile(mdname):
        stop(('The .md file '+mdname+' does not exist') if mdname else 'No .md file : give it with -f (<seed>.md, the trajectory of a CASTEP molecular dynamics run)')
    if SkipStep < 1 or InitialStep < 0:
        stop('-s must be at least 1 and -i at least 0')
    if outname is None:
        outname = (mdname[:-3] if mdname.endswith('.gz') else mdname)+'.umd.dat'
    if not outname.endswith('.umd.dat'):
        stop('The name of the UMD file must end with .umd.dat')
    base = outname[:-8]

    #First pass : the frames and their times ; a last frame cut short ; the frames written again after a restart
    frames, labels, fit, unknown = scan_md(mdname)
    natom = len(labels)
    if natom == 0:
        stop('No frame with positions (<-- R) in',mdname,': not a .md file of CASTEP ?')
    if len(set(labels)) != natom:
        stop('Several lines <-- R for the same atom in the first frame of',mdname,'(several beads per atom, path-integral MD ?) : not handled')
    species = [label[0] for label in labels]
    fit = [row for row, frame in zip(fit,frames) if complete(frame[2],natom)]
    nframes = len(frames)
    if not complete(frames[-1][2],natom):
        print('The last frame of',mdname,'(line',str(frames[-1][1])+') is not complete : left out')
        frames = frames[:-1]
    for time, line, counts in frames:
        if not complete(counts,natom):
            stop('The frame beginning at line',line,'of',mdname,'is not complete, or has not the',natom,'atoms of the first frame (',
                 ' '.join(str(n)+' <-- '+tag for tag,n in counts.items() if tag != 'cut'),')')
    if not frames:
        stop('No complete frame in',mdname)
    rawtimes = np.array([frame[0] for frame in frames])
    moves = np.abs(np.diff(rawtimes))
    smallest = moves[moves > 0.0].min() if np.any(moves > 0.0) else 0.0
    keep = last_branch(rawtimes,0.5*smallest)
    if not keep.all():
        back = [k for k in range(1,len(rawtimes)) if rawtimes[k] <= rawtimes[k-1]+0.5*smallest]
        print('WARNING : the time goes back',len(back),'times in',mdname,'(first from %.6f to %.6f fs) : a run continued from an earlier time,' %
              (rawtimes[back[0]-1]*AUTIME,rawtimes[back[0]]*AUTIME))
        print('   appended to the same file. Of the frames written several times the last ones are kept :',int((~keep).sum()),'frames left out')
    times = rawtimes[keep]
    nkept = len(times)
    #Times in fs. With the frames at whole numbers of their smallest interval (a constant time step) they are counted in intervals, so that the
    #rounding of the unit of time by CASTEP does not show (999.999999 fs for the 1000th step of 1 fs) ; from 0 when the first time is one too
    interval = np.diff(times).min() if nkept > 1 else 0.0
    counts = (times-times[0])/interval if interval > 0.0 else np.zeros(nkept)
    regular = interval > 0.0 and np.abs(counts-np.rint(counts)).max() < 1.0E-4
    if regular:
        counts = np.rint(counts).astype(np.int64)
        dt = round_step(interval*AUTIME)
        first = times[0]/interval
        timefs = np.round((np.rint(first)*dt if abs(first-np.rint(first)) < 1.0E-4 else times[0]*AUTIME)+counts*dt,6)
    else:
        dt = None
        timefs = np.round(times*AUTIME,6)
        if nkept > 1:
            print('The times of the frames are not whole numbers of one interval (time step changed during the run ?) : each time is converted as it is')
    written = list(range(InitialStep,nkept,SkipStep))
    if not written:
        stop('No frame to write : the file has',nkept,'frames, and -i is',InitialStep)

    def between(k):#timestep of the k-th snapshot written : the time from the snapshot written before it (for the first : to the second)
        a, b = (written[k-1],written[k]) if k > 0 else (written[0],written[min(1,len(written)-1)])
        if a == b:
            return round(dt*SkipStep,6) if dt is not None else 1.0
        return round(float(counts[b]-counts[a])*dt,6) if regular else round(float(timefs[b]-timefs[a]),6)

    #Elements and masses
    symbols = [element_of(s,kinds) for s in species]
    strangers = sorted(set(s for s,el in zip(species,symbols) if el is None))
    if strangers:
        stop('Unknown elements',strangers,': give the element of these species with -k, e.g. -k',strangers[0]+':O')
    renamed = sorted(set((s,el) for s,el in zip(species,symbols) if s != el))
    if renamed:
        print('species read as elements :',' '.join(s+'='+el for s,el in renamed),'(-k to change)')
    elements = list(dict.fromkeys(symbols))                             #in the order of first appearance
    elindex = np.array([elements.index(el) for el in symbols])
    order = np.lexsort((np.arange(natom),elindex))                      #grouped by element, in the order of the file within an element
    if set(masses_opt)-set(elements):
        stop('-m : no atoms of',' '.join(sorted(set(masses_opt)-set(elements))),'in the run')
    fitted = fit_masses(fit,species)
    masses = []
    source = []
    for el in elements:
        ofspecies = list(dict.fromkeys(s for s,e in zip(species,symbols) if e == el))
        values = [fitted[s] for s in ofspecies if s in fitted]
        if el in masses_opt:
            masses.append(masses_opt[el])
            source.append('-m')
        elif len(values) == len(ofspecies) and max(values)-min(values) <= 1.0E-4*max(values):
            masses.append(float('%.8g' % (sum(values)/len(values))))
            source.append('the kinetic energies of the file')
        else:
            if len(values) == len(ofspecies):
                print('WARNING : the species of',el,'(',' '.join(ofspecies),') have different masses in the run (',' '.join('%.8g' % m for m in values),
                      ') : the UMD file has one mass per element, that of the table (-m to give another)')
            masses.append(cr.Elements2rest(el)[3])
            source.append('the table of elements')
    for el, m, origin in zip(elements,masses,source):
        table = cr.Elements2rest(el)[3]
        if origin == 'the kinetic energies of the file' and abs(m-table) > 1.0E-3*table:
            print('Note : the mass of',el,'in the run is',m,'; the table of elements has',table,'(an isotope ?)')

    MyCrystal = cr.Lattice()
    MyCrystal.natom = natom
    MyCrystal.ntypat = len(elements)
    MyCrystal.elements = elements
    MyCrystal.types = [int((elindex == k).sum()) for k in range(len(elements))]
    MyCrystal.typat = list(elindex[order])
    MyCrystal.masses = masses
    MyCrystal.zelec = [0.0 for _ in elements]
    MyCrystal.noelectrons = 0.0
    MyCrystal.atoms = [cr.Atom() for _ in range(natom)]
    atommass = np.array(masses)[elindex[order]]
    open(outname,'w').close()
    umdpf.print_header(base,MyCrystal,'CASTEP : the UMD file originated from a CASTEP ab initio molecular dynamics simulation')
    print('Trajectory',mdname,':',nframes,'frames ;',('%.9g fs between frames' % dt) if dt is not None else 'one frame' if nkept == 1 else 'irregular times')
    print(natom,'atoms :',' '.join(str(n)+' '+el for n,el in zip(MyCrystal.types,elements)),'; masses',' '.join(str(m) for m in masses),
          'from',' and '.join(dict.fromkeys(source)))

    #Second pass : the snapshots
    reflabels = previous = firstcell = None
    boxmove = 0.0
    changing = skewed = False
    nwritten = 0
    lacking = {'V':0, 'F':0, 'S':0, 'T':0, 'pressure':0}
    ratios = []                                                         #kinetic energy of the velocities with the masses of the header / that of the file
    for index, frame in enumerate(read_md(mdname,keep)):
        if reflabels is None:
            reflabels = frame['labels']
        elif not np.array_equal(frame['labels'],reflabels):
            stop('The atoms of the frame at',timefs[index],'fs are not those of the first frame')
        cell = frame['h']*BOHR
        volume = abs(np.linalg.det(cell))
        if volume < 1.0E-9:
            stop('A cell without volume in the frame at',timefs[index],'fs')
        if firstcell is None:
            firstcell = cell
            skewed = np.abs(cell-np.diag(np.diag(cell))).max() > 1.0E-4*np.abs(cell).max()
        changing = changing or not np.array_equal(cell,firstcell)
        #Unwrapped positions : CASTEP does not wrap the atoms, so su is the reduced coordinate of the file ; followed from frame to frame (d - rint(d)), which
        #changes nothing then and removes whole cell vectors if the positions were wrapped. Each move is the change of su in the cell of the frame, so that
        #the atoms do not move with a changing cell (NPT ; as LAMMPSParser2umd.py and CP2KParser.py) : boxmove removes su dh. 0 for a constant cell
        s = (frame['R'][order]*BOHR) @ np.linalg.inv(cell)
        xred = s - np.floor(s)
        xcart = xred @ cell
        if previous is None:
            su = s
        else:
            d = s - previous
            su = previoussu + d - np.rint(d)
            if not np.array_equal(cell,previouscell):
                boxmove = boxmove + previoussu @ (previouscell-cell)
        absx = su @ cell + boxmove
        previous, previoussu, previouscell = s, su, cell
        if index < InitialStep or (index-InitialStep) % SkipStep:
            continue
        for tag in ('V','F','S','T'):
            lacking[tag] += frame[tag] is None
        vels = frame['V'][order]*(BOHR/AUTIME) if frame['V'] is not None else np.zeros((natom,3))
        forces = frame['F'][order]*(HARTREE/BOHR) if frame['F'] is not None else np.zeros((natom,3))
        PE, conserved, KE = (x*HARTREE for x in frame['E'])
        if KE > 0.0 and frame['V'] is not None:
            ratios.append(0.5*(atommass*(vels**2).sum(axis=1)).sum()*AMU_A2FS2/KE)
        T = frame['T']*HARTREE/KB if frame['T'] is not None else 2.0*KE/(max(3*natom-3,1)*KB)
        if frame['S'] is not None:                                      #UMD : positive in compression, with the kinetic part ; xx yy zz xy yz zx
            S = -frame['S']*HAB3_GPA
            stress = [S[0,0],S[1,1],S[2,2],0.5*(S[0,1]+S[1,0]),0.5*(S[1,2]+S[2,1]),0.5*(S[2,0]+S[0,2])]
        else:
            stress = [0.0 for _ in range(6)]
        haspressure = frame['P'] is not None or frame['S'] is not None
        lacking['pressure'] += not haspressure
        P = frame['P']*HAB3_GPA if frame['P'] is not None else sum(stress[0:3])/3.0
        #The cell is written to 1e-7 A, the digits of the main output of CASTEP : its versions differ by 1e-9 in the bohr (9.9999999945 A for a cell of 10 A)
        record = {'time':float(timefs[index]), 'cell':np.round(cell,7), 'internalenergy':PE+KE, 'toten':PE, 'kineticenergy':KE, 'energywithdrift':conserved,
                  'enthalpy':PE+KE+P*volume/EVA3_GPA if haspressure else 0.0, 'temperature':T, 'pressure':P, 'stress':stress,
                  'values':[xred.tolist(), xcart.tolist(), absx.tolist(), vels.tolist(), forces.tolist(), [0.0 for _ in range(natom)]]}
        write_snapshot(base,MyCrystal,atommass,record,between(nwritten))
        nwritten += 1
    print(nkept,'frames read,',nwritten,'written in',outname)
    if unknown:
        print('Lines not used :',' '.join('<-- '+tag+' ('+str(n)+' lines)' for tag,n in unknown.items()))
    for tag, effect in (('V','No velocities (<-- V) : velocities 0'),('F','No forces (<-- F) : forces 0'),
                        ('S','No stress (<-- S ; run without calculate_stress) : StressTensor 0'),
                        ('pressure','No pressure (<-- P) : Pressure and Enthalpy 0'),
                        ('T','No temperature (<-- T) : that of the kinetic energy, with 3N-3 degrees of freedom')):
        if lacking[tag]:
            print(effect+('' if lacking[tag] == nwritten else ' in '+str(lacking[tag])+' of the '+str(nwritten)+' snapshots'))
    if changing:
        print('The cell changes during the run : each snapshot has its cell. Note : not yet checked on a real CASTEP run with a variable cell.')
    if skewed:
        print('Note : the cell vectors are not along x, y, z. The lines <-- h are read as the vectors a, b, c (as the rows of lattice_cart) ; not yet')
        print('   checked on a real CASTEP run : compare the cell of the UMD file with that of the run.')
    if ratios:
        ratio = float(np.median(ratios))
        print('Check : kinetic energy of the velocities with these masses / kinetic energy of the file = %.6f (median of %d snapshots)' % (ratio,len(ratios)))
        if abs(ratio-1.0) > 1.0E-3:
            print('WARNING : they differ : the masses used here are not those of the run.',('Masses that give the kinetic energies of the file : -m '+
                  ','.join('%s:%.8g' % (s,m) for s,m in fitted.items())) if fitted else 'No set of masses gives the kinetic energies of the file : check the file.')


if __name__ == "__main__":
    main(sys.argv[1:])
