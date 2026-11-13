#!/usr/bin/env python3
###
#AUTHORS: RAZVAN CARACAS
###
#Converts a CP2K molecular dynamics run into a UMD file. CP2K writes a run in several files : the positions (<project>-pos-1.xyz), the velocities
#(-vel-1.xyz), the forces (-frc-1.xyz), the energies (-1.ener), the stress tensor (-1.stress) and the cell (-1.cell), matched here by their MD step.
#The xyz files have no cell : it comes from -c (numbers, the cell file, an input or output file of CP2K), or from the Lattice of extended xyz files.
#Units of CP2K by default : positions A, velocities bohr/a.u. of time, forces Hartree/bohr, energies Hartree, stress bar ; converted to A, fs, eV, GPa.

import sys, getopt, os.path, re
import numpy as np
import crystallography as cr
import umd_processes_fast as umdpf
from vibrations_umd import KB, AMU_A2FS2

BOHR = 0.529177210903               #bohr in A
HARTREE = 27.211386245988           #Hartree in eV
AUTIME = 2.418884326505E-02         #atomic unit of time in fs
EVA3_GPA = 160.2176634              #eV/A^3 in GPa
BAR = 1.0E-4                        #bar in GPa
NUMBER = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?'
USAGE = 'CP2KParser.py -f <positions> -v <velocities> -r <forces> -e <energies> -a <stress> -c <cell> -m <masses> -k <kinds> -t <timestep> -o <output> -i <InitialStep> -s <SkipStep>'


def number(text):
    return float(text.replace('D','E').replace('d','e'))


def element_of(name,kinds={}):#Element of an atom name of CP2K : the symbol, or for a kind name (O1, Hw, OW, Si_a...) the symbol it begins with
    if name in kinds:
        return kinds[name]
    if cr.element_entry(name) is not None:
        return name
    found = re.match(r'[A-Z][a-z]?',name)
    for candidate in ((found.group(0),found.group(0)[0]) if found else ()):
        if cr.element_entry(candidate) is not None:
            return candidate
    return None


#****Cell

def cell_from_parameters(a,b,c,alpha=90.0,beta=90.0,gamma=90.0):#Rows a, b, c with a along x and b in the xy plane, as CP2K builds the cell from ABC and ALPHA_BETA_GAMMA
    if max(abs(alpha-90.0),abs(beta-90.0),abs(gamma-90.0)) < 1e-12:
        return np.diag([float(a),float(b),float(c)])
    ca, cb, cg, sg = np.cos(np.radians(alpha)), np.cos(np.radians(beta)), np.cos(np.radians(gamma)), np.sin(np.radians(gamma))
    cy = (ca-cb*cg)/sg
    return np.array([[a,0.0,0.0],[b*cg,b*sg,0.0],[c*cb,c*cy,c*np.sqrt(max(1.0-cb*cb-cy*cy,0.0))]])


def read_cell_file(name):#<project>-1.cell : '# Step Time [fs] Ax Ay Az Bx By Bz Cx Cy Cz [Angstrom] Volume' -> {step : cell}, rows a, b, c
    cells = {}
    with open(name,'r') as f:
        for line in f:
            if line.startswith('#'):
                if 'Angstrom' not in line:
                    return {}
                continue
            tokens = line.split()
            try:
                cells[int(tokens[0])] = np.array([float(x) for x in tokens[2:11]]).reshape(3,3)
            except (ValueError,IndexError):
                return {}
    return cells


def read_cp2k_input(name):#&CELL of a CP2K input or restart file (ABC, ALPHA_BETA_GAMMA, or the vectors A, B, C ; [angstrom] or [bohr]), and the MASS and
    #ELEMENT of its &KIND sections ; the &CELL_REF section inside &CELL is not the cell. Returns the cell (or None), {element : masses of its kinds, None
    #without MASS}, and a message when the cell cannot be used
    sections = []
    values = {}
    kinds = {}
    kind = None
    with open(name,'r') as f:
        for line in f:
            line = re.split(r'[!#]',line)[0].strip()
            if not line:
                continue
            tokens = line.split()
            word = tokens[0].upper()
            if word.startswith('&END'):
                if sections:
                    sections.pop()
                continue
            if word.startswith('&'):
                sections.append(word[1:])
                if word == '&KIND' and len(tokens) > 1:
                    kind = tokens[1]
                    kinds[kind] = {'element':element_of(kind), 'mass':None}
                continue
            if sections and sections[-1] == 'CELL':
                unit = 1.0
                rest = tokens[1:]
                if rest and rest[0].startswith('['):
                    unitname = rest[0].strip('[]').lower()
                    if unitname in ('bohr','au','a.u.'):
                        unit = BOHR
                    elif unitname not in ('angstrom','deg','degree'):
                        return None, {}, 'unit '+rest[0]+' of '+word+' in its &CELL section'
                    rest = rest[1:]
                if word in ('ABC','A','B','C','ALPHA_BETA_GAMMA'):
                    try:
                        values[word] = [number(x)*(1.0 if word == 'ALPHA_BETA_GAMMA' else unit) for x in rest[0:3]]
                    except ValueError:
                        return None, {}, 'the line '+line+' of its &CELL section'
                elif word == 'MULTIPLE_UNIT_CELL' and rest[0:3] != ['1','1','1']:
                    return None, {}, 'MULTIPLE_UNIT_CELL '+' '.join(rest)+' (the cell of the run is a multiple of its &CELL)'
                elif word in ('CELL_FILE_NAME','CELL_FILE_FORMAT'):
                    return None, {}, 'its cell is in another file ('+word+')'
            elif sections and sections[-1] == 'KIND' and kind is not None and len(tokens) > 1:
                if word == 'MASS':
                    try:
                        kinds[kind]['mass'] = number(tokens[1])
                    except ValueError:
                        pass
                elif word == 'ELEMENT':
                    kinds[kind]['element'] = element_of(tokens[1].capitalize())
    masses = {}
    for entry in kinds.values():
        masses.setdefault(entry['element'],set()).add(entry['mass'])
    if all(k in values and len(values[k]) == 3 for k in ('A','B','C')):
        return np.array([values['A'],values['B'],values['C']]), masses, None
    if 'ABC' in values and len(values['ABC']) == 3:
        return cell_from_parameters(*values['ABC'],*values.get('ALPHA_BETA_GAMMA',[90.0,90.0,90.0])), masses, None
    return None, masses, None


def read_cp2k_output(name):#The first cell of the main output file of CP2K : 'CELL| Vector a [angstrom]: x y z |a| = ...' (printed with the input, at
    #the medium print level), else 'MD_INI| Cell lengths [ang]' and 'MD_INI| Cell angles [deg]' (printed at the start of every MD run, also at the low
    #print level ; a along x, b in the xy plane), with the digits printed there
    cell = {}
    lengths = angles = None
    with open(name,'r') as f:
        for line in f:
            found = re.match(r'\s*CELL\|\s+Vector ([abc]) \[angstrom\]:\s+('+NUMBER+r')\s+('+NUMBER+r')\s+('+NUMBER+')',line)
            if found and found.group(1) not in cell:
                cell[found.group(1)] = [float(found.group(k)) for k in (2,3,4)]
                if len(cell) == 3:
                    return np.array([cell['a'],cell['b'],cell['c']])
            found = re.match(r'\s*MD_INI\|\s+Cell (lengths \[ang\]|angles \[deg\])\s+('+NUMBER+r')\s+('+NUMBER+r')\s+('+NUMBER+')',line)
            if found:
                values = [float(found.group(k)) for k in (2,3,4)]
                if found.group(1).startswith('lengths') and lengths is None:
                    lengths = values
                elif found.group(1).startswith('angles') and angles is None:
                    angles = values
                if lengths is not None and angles is not None:
                    return cell_from_parameters(*lengths,*angles)
    return cell_from_parameters(*lengths) if lengths is not None else None


def read_cell_option(text):#-c : numbers (a | a,b,c | a,b,c,alpha,beta,gamma | the 9 components of a, b, c ; A), or a file (cell file, input, restart or
    #main output of CP2K). Returns the constant cell or None, the cells per step or None, the masses of the &KIND sections, what it was read from
    if not os.path.isfile(text):
        try:
            values = [float(x) for x in text.replace(',',' ').split()]
        except ValueError:
            print('-c',text,': neither a file nor numbers')
            sys.exit()
        if len(values) == 1:
            return np.diag(values*3), None, {}, '-c (cubic)'
        if len(values) == 3:
            return np.diag(values), None, {}, '-c (a, b, c)'
        if len(values) == 6:
            return cell_from_parameters(*values), None, {}, '-c (a, b, c, alpha, beta, gamma ; a along x, b in the xy plane)'
        if len(values) == 9:
            return np.array(values).reshape(3,3), None, {}, '-c (vectors a, b, c)'
        print('-c takes 1 (cubic), 3 (a,b,c), 6 (a,b,c,alpha,beta,gamma) or 9 (the vectors a, b, c) numbers, or a file')
        sys.exit()
    cells = read_cell_file(text)
    if cells:
        return None, cells, {}, 'the cell file '+text+' (the cell of every step)'
    cell, masses, problem = read_cp2k_input(text)
    if problem:
        print('The cell cannot be taken from',text,':',problem,'; give it with -c as numbers')
        sys.exit()
    if cell is not None:
        return cell, None, masses, 'the &CELL section of '+text
    cell = read_cp2k_output(text)
    if cell is not None:
        return cell, None, {}, 'the main output '+text+' (CELL| or MD_INI| lines, with the digits printed there)'
    print('No cell found in',text,': not a cell file (# Step Time Ax ... Cz Volume), nor an input with a &CELL section, nor an output with CELL| Vector or MD_INI| Cell lines')
    sys.exit()


#****Tables by MD step : energies and stress

def read_ener(name):#<project>-1.ener : '# Step Nr. Time[fs] Kin.[a.u.] Temp[K] Pot.[a.u.] Cons Qty[a.u.] UsedTime[s]' -> {step : (time, kin, temp, pot, cons)}
    #a step written twice (run restarted from an earlier step) : the last one. Returns also the number of such steps
    rows = {}
    repeated = 0
    with open(name,'r') as f:
        for line in f:
            if line.startswith('#'):
                header = re.sub(r'\s+','',line)
                if not all(k in header for k in ('Time','Kin.[a.u.]','Temp[K]','Pot.[a.u.]','ConsQty[a.u.]')) or \
                   not header.index('Time') < header.index('Kin.') < header.index('Temp') < header.index('Pot.') < header.index('ConsQty'):
                    print('The columns of the energy file',name,'are not Step, Time[fs], Kin.[a.u.], Temp[K], Pot.[a.u.], Cons Qty[a.u.] :')
                    print('  ',line.strip())
                    sys.exit()
                continue
            tokens = line.split()
            try:
                step = int(tokens[0])
                values = tuple(float(x) for x in tokens[1:6])
            except (ValueError,IndexError):
                continue
            if len(values) == 5:
                repeated += step in rows
                rows[step] = values
    return rows, repeated


def read_stress(name):#<project>-1.stress : '# Step Time [fs] xx xy xz yx yy yz zx zy zz [bar]', the total pressure tensor of CP2K (virial + kinetic part,
    #positive in compression) -> {step : xx yy zz xy yz zx in GPa}, the order of the UMD file
    rows = {}
    with open(name,'r') as f:
        for line in f:
            if line.startswith('#'):
                units = set(re.findall(r'\[(\w+)\]',line)) - {'fs'}
                if units and units != {'bar'}:
                    print('The stress file',name,'is not in bar :',line.strip())
                    sys.exit()
                continue
            tokens = line.split()
            try:
                step = int(tokens[0])
                xx, xy, xz, yx, yy, yz, zx, zy, zz = (float(x)*BAR for x in tokens[2:11])
            except (ValueError,IndexError):
                continue
            rows[step] = [xx,yy,zz,0.5*(xy+yx),0.5*(yz+zy),0.5*(zx+xz)]
    return rows


#****Trajectory files

def read_title(line):#' i = N, time = t, E = e' (XMOL), or 'Lattice="..." Properties=... Step=N Time=t Energy=e' (extended xyz) : step, time, energy, cell
    cell = None
    found = re.search(r'Lattice\s*=\s*"([^"]*)"',line)
    if found:
        values = found.group(1).split()
        if len(values) == 9:
            cell = np.array([float(x) for x in values]).reshape(3,3)
        step = re.search(r'Step\s*=\s*(\d+)',line)
        time = re.search(r'Time\s*=\s*('+NUMBER+')',line)
        energy = re.search(r'Energy\s*=\s*('+NUMBER+')',line)
    else:
        step = re.search(r'\bi\s*=\s*(\d+)',line)
        time = re.search(r'\btime\s*=\s*('+NUMBER+')',line)
        energy = re.search(r'\bE\s*=\s*('+NUMBER+')',line)
    return (int(step.group(1)) if step else None, float(time.group(1)) if time else None, float(energy.group(1)) if energy else None, cell)


def number_of_atoms(name,line):
    try:
        return int(line.split()[0])
    except ValueError:
        print('The file',name,'is not an xyz file : a line with the number of atoms was expected, found')
        print('  ',line.strip()[:100])
        sys.exit()


def scan_xyz(name):#Quick pass over an xyz file : the step and the time of every frame
    steps = []
    times = []
    with open(name,'r') as f:
        while True:
            line = f.readline()
            if not line.strip():
                return steps, times
            natom = number_of_atoms(name,line)
            step, time, energy, cell = read_title(f.readline())
            for _ in range(natom):
                f.readline()
            steps.append(step)
            times.append(time)


def last_branch(steps):#CP2K appends to its files when a run is restarted : the frames after the restart point are then written again. True for the frames
    #to keep : those with no later frame at the same or at an earlier step
    keep = np.ones(len(steps),dtype=bool)
    if len(steps) > 1 and all(s is not None for s in steps):
        later = np.minimum.accumulate(np.array(steps)[::-1])[::-1]
        keep[:-1] = np.array(steps)[:-1] < later[1:]
    return keep


def read_xyz(name,keep=None):#Yields the frames of a CP2K xyz file (XMOL or extended xyz), but those not in keep : step, time, energy, cell, names of
    #the atoms, array natom x 3. A last frame cut short (run still being written) is left out
    with open(name,'r') as f:
        index = -1
        while True:
            line = f.readline()
            if not line.strip():
                return
            natom = number_of_atoms(name,line)
            title = f.readline()
            lines = [f.readline() for _ in range(natom)]
            index += 1
            if keep is not None and index < len(keep) and not keep[index]:
                continue
            try:
                rows = np.array(' '.join(lines).split()).reshape(natom,-1)
                values = rows[:,1:4].astype(float)
                if not lines[-1] or values.shape != (natom,3):
                    raise ValueError
            except ValueError:
                print('The last frame of',name,'is not complete : left out')
                return
            yield read_title(title) + (rows[:,0],values)


def frames_by_step(name,restarted):#The frames of a velocity or force file, asked by increasing step : get(step,index,names) -> array, or None without a frame
    #at this step (by position in the file when the titles have no step)
    reader = read_xyz(name,last_branch(scan_xyz(name)[0]) if restarted else None)
    state = {'frame':next(reader,None), 'index':0}
    def get(step,index,names):
        def position(frame):
            return (frame[0],step) if step is not None and frame[0] is not None else (state['index'],index)
        while state['frame'] is not None and position(state['frame'])[0] < position(state['frame'])[1]:
            state['frame'] = next(reader,None)
            state['index'] += 1
        frame = state['frame']
        if frame is None or position(frame)[0] != position(frame)[1]:
            return None
        if not np.array_equal(frame[4],names):
            print('The atoms of',name,'at step',step,'are not those of the positions file')
            sys.exit()
        return frame[5]
    return get


#****Main

def main(argv):
    umdpf.headerumd()
    posname = velname = frcname = enername = stressname = ''
    cellopt = None
    outname = None
    masses_opt = {}
    kinds = {}
    timestep = None
    InitialStep = 0
    SkipStep = 1
    try:
        opts, arg = getopt.getopt(argv,"hf:v:r:e:a:c:m:k:t:o:i:s:")
    except getopt.GetoptError:
        print(USAGE)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('CP2KParser.py : converts a CP2K molecular dynamics run into a UMD file')
            print(USAGE)
            print('-f : positions, <project>-pos-1.xyz (&MOTION &PRINT &TRAJECTORY ; format XMOL, the default, or EXTXYZ ; in A, the default unit).')
            print('     With this standard name the other files of the run are found in the same folder when not given : <project>-vel-1.xyz,')
            print('     <project>-frc-1.xyz, <project>-1.ener, <project>-1.stress, <project>-1.cell. Files with other names : -v -r -e -a -c.')
            print('     All the files are matched by the MD step of their frames and lines, so they can be written at different intervals.')
            print('-v : velocities, <project>-vel-1.xyz (&PRINT &VELOCITIES ; in bohr/a.u. of time, the default unit). Without them : velocities 0.')
            print('-r : forces, <project>-frc-1.xyz (&PRINT &FORCES ; in Hartree/bohr, the default unit). Without them : forces 0.')
            print('-e : energies, <project>-1.ener (columns Step, Time[fs], Kin.[a.u.], Temp[K], Pot.[a.u.], Cons Qty[a.u.]).')
            print('     InternalEnergy = Pot. + Kin. (potential + kinetic, as VaspParser.py), TotalElectronicEnergy = Pot., KineticEnergyIons = Kin.,')
            print('     EnergyWithDrift = Cons Qty (the conserved quantity, with the energy of the thermostat), Temperature = Temp.')
            print('     Without it : the potential energy from the E of the positions file, kinetic energy and temperature (3N-3) from the velocities.')
            print('-a : stress tensor, <project>-1.stress (&PRINT &STRESS, with STRESS_TENSOR ANALYTICAL in &FORCE_EVAL ; in bar) : the total pressure')
            print('     tensor of CP2K, kinetic part included, positive in compression, as the StressTensor of the UMD file. Pressure = its trace / 3,')
            print('     Enthalpy = InternalEnergy + P V. Without it : pressure, stress and enthalpy 0.')
            print('-c : cell. The xyz files of CP2K have none (but EXTXYZ : the Lattice of every frame). One of :')
            print('     numbers in A : a (cubic) | a,b,c | a,b,c,alpha,beta,gamma (a along x, b in the xy plane, as CP2K) | the 9 components of a, b, c ;')
            print('     <project>-1.cell (&PRINT &CELL), the cell of every step (NPT) ;')
            print('     an input or restart file of CP2K (its &CELL section : ABC, ALPHA_BETA_GAMMA, or A, B, C ; also the MASS of its &KIND sections) ;')
            print('     the main output file of CP2K (its lines CELL| Vector a [angstrom] ..., else MD_INI| Cell lengths [ang] and angles [deg] ;')
            print('     the cell at the start of that output : for an NPT run, give the cell file instead).')
            print('-m : masses of the run, when they are not the standard ones, e.g. -m H:2.01410 for deuterium (default : MASS of the &KIND')
            print('     sections of the input given with -c, else the table of elements). With velocities and energies the parser compares')
            print('     the kinetic energy of the velocities with that of the energy file, and warns when the masses or the units do not match.')
            print('-k : element of the atom kinds, when the xyz files have kind names (PRINT_ATOM_KIND), e.g. -k OW:O,HW:H')
            print('     (default : the element symbol the name begins with : O1 -> O, Hw -> H, OW -> O).')
            print('-t : time step of the MD in fs (default : from the times of the energy file, else of the positions file)')
            print('-o : name of the UMD file (default <positions file>.umd.dat)')
            print('-i : first frame written (default 0) ; -s : every s-th frame (default 1)')
            print('The atoms are grouped by element, in the order of the xyz file within an element. Unwrapped positions : those of CP2K, which does')
            print('not wrap the atoms during a run, followed from frame to frame (also the frames not written) ; with a changing cell (NPT) each move is')
            print('the change of the reduced coordinates in the cell of the frame, so that the atoms do not move with the cell.')
            print('A run restarted from an earlier step appends to the same files : of the steps written several times the last ones are kept.')
            sys.exit()
        elif opt == '-f':
            posname = str(arg)
        elif opt == '-v':
            velname = str(arg)
        elif opt == '-r':
            frcname = str(arg)
        elif opt == '-e':
            enername = str(arg)
        elif opt == '-a':
            stressname = str(arg)
        elif opt == '-c':
            cellopt = str(arg)
        elif opt == '-m':
            try:
                masses_opt = {item.split(':')[0].strip():float(item.split(':')[1]) for item in arg.split(',')}
            except (ValueError,IndexError):
                print('-m : masses as element:mass, e.g. -m H:2.01410,O:15.9994')
                sys.exit()
        elif opt == '-k':
            kinds = {item.partition(':')[0].strip():item.partition(':')[2].strip() for item in arg.split(',')}
            if any(cr.element_entry(el) is None for el in kinds.values()):
                print('-k : kinds as kind:element, e.g. -k OW:O,HW:H')
                sys.exit()
        elif opt == '-t':
            timestep = float(arg)
        elif opt == '-o':
            outname = str(arg)
        elif opt == '-i':
            InitialStep = int(arg)
        elif opt == '-s':
            SkipStep = int(arg)
    if not os.path.isfile(posname):
        print(('The positions file '+posname+' does not exist') if posname else 'No positions file : give it with -f (<project>-pos-1.xyz)')
        sys.exit()
    if SkipStep < 1 or InitialStep < 0:
        print('-s must be at least 1 and -i at least 0')
        sys.exit()
    folder, name = os.path.split(posname)
    standard = re.match(r'(.*)-pos-(\d+)\.xyz$',name)                   #<project>-pos-1.xyz : the other files of the run, in the same folder
    if standard:
        found = []
        for current, tail in ((velname,'-vel-%s.xyz'),(frcname,'-frc-%s.xyz'),(enername,'-%s.ener'),(stressname,'-%s.stress'),(cellopt,'-%s.cell')):
            candidate = os.path.join(folder,standard.group(1)+tail % standard.group(2))
            found.append(candidate if not current and os.path.isfile(candidate) else current)
        velname, frcname, enername, stressname, cellopt = found
    for filename in (velname,frcname,enername,stressname):
        if filename and not os.path.isfile(filename):
            print('The file',filename,'does not exist')
            sys.exit()
    if outname is None:
        outname = posname+'.umd.dat'
    if not outname.endswith('.umd.dat'):
        print('The name of the UMD file must end with .umd.dat')
        sys.exit()
    base = outname[:-8]
    for filename in (posname,velname,frcname):                          #<base>.xyz is the XYZ output of umd2out.py for <base>.umd.dat
        if filename and os.path.abspath(filename) == os.path.abspath(base+'.xyz'):
            print('The name',outname,'cannot be used : umd2out.py would write its XYZ trajectory,',base+'.xyz',', over this file of CP2K, and')
            print('umd_clean.py would take it as an output of the package. Choose another name with -o (default',posname+'.umd.dat)')
            sys.exit()
    constantcell, cells, inputmasses, cellsource = read_cell_option(cellopt) if cellopt else (None,None,{},None)
    ener, repeated = read_ener(enername) if enername else ({},0)
    stressrows = read_stress(stressname) if stressname else {}
    print('Positions',posname,'; velocities',velname or 'none','; forces',frcname or 'none')
    print('Energies',enername or 'none','; stress',stressname or 'none','; cell from',cellsource or 'the positions file (extended xyz)')
    for filename, table, what in ((enername,ener,'energy'),(stressname,stressrows,'stress')):
        if filename and not table:
            print('No line of numbers in the',what,'file',filename)
            sys.exit()

    #Steps and times of the frames ; the frames written again after a restart ; the time step of the MD
    steps, times = scan_xyz(posname)
    keep = last_branch(steps)
    restarted = not keep.all()
    if restarted:
        back = [k for k in range(1,len(steps)) if steps[k] <= steps[k-1]]
        print('WARNING : the steps go back',len(back),'times in the positions file (first from step',steps[back[0]-1],'to',str(steps[back[0]])+') : a run restarted from')
        print('   an earlier step, appended to the same files. Of the steps written several times the last ones are kept :',int((~keep).sum()),'frames left out')
    elif repeated:
        print('WARNING :',repeated,'steps appear several times in the energy file (run restarted from an earlier step ?) : the last ones are used')
    if timestep is None:
        marks = sorted((s,v[0]) for s,v in ener.items()) if len(ener) > 1 else [(s,t) for s,t,k in zip(steps,times,keep) if k and s is not None and t is not None]
        if len(marks) > 1 and marks[-1][0] > marks[0][0]:
            timestep = (marks[-1][1]-marks[0][1])/(marks[-1][0]-marks[0][0])
            print('Time step of the MD : %.9g fs (from the times of the %s)' % (timestep,'energy file' if len(ener) > 1 else 'positions file ; -t to give it'))
    getvels = frames_by_step(velname,restarted) if velname else None
    getforces = frames_by_step(frcname,restarted) if frcname else None

    def between(count,lastcount,time,lasttime):                         #time between two snapshots written : MD steps x time step, else the difference of their times
        return round((count-lastcount)*timestep if timestep is not None else time-lasttime,6)

    MyCrystal = None
    nframes = nwritten = nocell = 0
    missing = {'velocities':[], 'forces':[], 'energies':[], 'stress':[], 'cell':[]}
    previous = origin = None
    boxmove = 0.0
    kinetics = []                                                       #per frame written : the sum of v^2 of each element (A^2/fs^2), then Kin. of the energy file (eV)
    for step, time, energy, framecell, names, positions in read_xyz(posname,keep):
        natom = positions.shape[0]
        if MyCrystal is None:                                           #first frame : elements, masses, order of the atoms
            if framecell is None and constantcell is None and cells is None:
                print('The xyz files of CP2K have no cell : give it with -c (numbers, the cell file <project>-1.cell, or an input or output file of CP2K)')
                sys.exit()
            symbols = [element_of(n,kinds) for n in names]
            unknown = sorted(set(n for n,s in zip(names,symbols) if s is None))
            if unknown:
                print('Unknown elements',unknown,': give the element of these atom kinds with -k, e.g. -k',unknown[0]+':O')
                sys.exit()
            renamed = sorted(set((n,s) for n,s in zip(names,symbols) if n != s))
            if renamed:
                print('atom kinds read as elements :',' '.join(n+'='+s for n,s in renamed),'(-k to change)')
            elements = list(dict.fromkeys(symbols))                     #in the order of first appearance
            elindex = np.array([elements.index(s) for s in symbols])
            order = np.lexsort((np.arange(natom),elindex))              #grouped by element, in the order of the file within an element
            masses = []
            source = []
            for el in elements:
                ofinput = inputmasses.get(el,set())
                if el in masses_opt:
                    masses.append(masses_opt[el])
                    source.append('-m')
                elif len(ofinput) == 1 and None not in ofinput:
                    masses.append(next(iter(ofinput)))
                    source.append('the &KIND sections of the input')
                else:
                    if len(ofinput) > 1:
                        print('WARNING : the kinds of',el,'do not all have the same MASS in the input : the UMD file has one mass per element, that of the table (-m to give it)')
                    masses.append(cr.Elements2rest(el)[3])
                    source.append('the table of elements')
            if set(masses_opt) - set(elements):
                print('-m : no atoms of',' '.join(sorted(set(masses_opt)-set(elements))),'in the run')
                sys.exit()
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
            groupindex = elindex[order]
            atommass = np.array(masses)[groupindex]
            refnames = names
            open(outname,'w').close()
            umdpf.print_header(base,MyCrystal,'CP2K : the UMD file originated from a CP2K molecular dynamics simulation')
            print(natom,'atoms :',' '.join(str(n)+' '+el for n,el in zip(MyCrystal.types,elements)),'; masses',' '.join(str(m) for m in masses),
                  'from',' and '.join(dict.fromkeys(source)))
        if natom != MyCrystal.natom or not np.array_equal(names,refnames):
            print('The atoms of the frame at step',step,'are not those of the first frame')
            sys.exit()
        index = nframes
        nframes += 1
        if cells is not None:                                           #the cell of every step ; a frame without its line : the cell of the frame before
            cell = cells.get(step)
            if cell is None:
                missing['cell'].append(step)
                if previous is None:
                    nocell += 1
                    continue
                cell = previouscell
        else:
            cell = constantcell if constantcell is not None else framecell
        if cell is None or abs(np.linalg.det(cell)) < 1e-9:
            print('No cell, or a cell without volume, at step',step)
            sys.exit()
        #Unwrapped positions : CP2K does not wrap the atoms, so su is the reduced coordinate of the file ; followed from frame to frame (d - rint(d)), which
        #changes nothing then and removes whole cell vectors if the positions were wrapped. Each move is the change of su in the cell of the frame, so that the
        #atoms do not move with a changing cell (NPT ; as LAMMPSParser2umd.py) : CP2K writes x = su h(t), and boxmove removes su dh. 0 for a constant cell
        s = positions[order] @ np.linalg.inv(cell)
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
        row = ener.get(step)
        count = step if step is not None else index
        if row is not None:                                             #6 decimals, against 3 in the titles of the xyz files
            time = row[0]
        elif timestep is not None:                                      #from the first time, by MD steps : the 3 decimals of the titles give irregular times
            if origin is None:                                          #when the time step is not a round number
                origin = (count,time if time is not None else count*timestep)
            time = origin[1]+(count-origin[0])*timestep
        elif time is None:
            print('No time in the title lines of',posname,': give the energy file with -e, or the time step of the MD with -t')
            sys.exit()
        time = round(time,6)
        vels = getvels(step,index,names) if getvels else None
        forces = getforces(step,index,names) if getforces else None
        stress = stressrows.get(step)
        for key, value, given in (('velocities',vels,velname),('forces',forces,frcname),('energies',row,enername),('stress',stress,stressname)):
            if given and value is None:
                missing[key].append(step)
        vels = vels[order]*(BOHR/AUTIME) if vels is not None else np.zeros((natom,3))
        forces = forces[order]*(HARTREE/BOHR) if forces is not None else np.zeros((natom,3))
        v2 = np.bincount(groupindex,weights=(vels**2).sum(axis=1),minlength=len(elements))
        kinetic = 0.5*(np.array(masses)*v2).sum()*AMU_A2FS2
        if row is not None:
            KE, T, PE, conserved = row[1]*HARTREE, row[2], row[3]*HARTREE, row[4]*HARTREE
            if kinetic > 0.0:
                kinetics.append(list(v2)+[KE])
        else:
            KE, T, PE = kinetic, 2.0*kinetic/(max(3*natom-3,1)*KB), (energy*HARTREE if energy is not None else 0.0)
            conserved = PE + KE
        P = sum(stress[0:3])/3.0 if stress is not None else 0.0
        record = {'time':time, 'cell':cell, 'internalenergy':PE+KE, 'toten':PE, 'kineticenergy':KE, 'energywithdrift':conserved,
                  'enthalpy':PE+KE+P*abs(np.linalg.det(cell))/EVA3_GPA if stress is not None else 0.0, 'temperature':T, 'pressure':P,
                  'stress':stress if stress is not None else [0.0 for _ in range(6)],
                  'values':[xred.tolist(), xcart.tolist(), absx.tolist(), vels.tolist(), forces.tolist(), [0.0 for _ in range(natom)]]}
        if nwritten == 0:                                               #the first snapshot waits for the second : its timestep is the time between them
            first, firstcount = record, count
        else:
            if nwritten == 1:
                write_snapshot(base,MyCrystal,atommass,first,between(count,firstcount,time,first['time']))
            write_snapshot(base,MyCrystal,atommass,record,between(count,lastcount,time,lasttime))
        lasttime, lastcount = time, count
        nwritten += 1
    if MyCrystal is None:
        print('No frame in the positions file',posname)
        sys.exit()
    if nwritten == 1:                                                   #a single snapshot : a structure
        write_snapshot(base,MyCrystal,atommass,first,round(timestep*SkipStep,6) if timestep is not None else 1.0)
    print(nframes,'frames read,',nwritten,'written in',outname)
    for key, filename, effect in (('velocities',velname,'their velocities are 0'),('forces',frcname,'their forces are 0'),
                                  ('energies',enername,'their potential energy is the E of the positions file, their kinetic energy and temperature those of the velocities'),
                                  ('stress',stressname,'their pressure, stress and enthalpy are 0'),
                                  ('cell',cellopt,'they have the cell of the frame before ('+str(nocell)+' before the first cell left out)')):
        if missing[key]:
            print('WARNING :',len(missing[key]),'frames have no',('frame' if key in ('velocities','forces') else 'line'),'at their step in',filename,
                  '(steps',missing[key][0],'to',str(missing[key][-1])+') :',effect)
    if not velname:
        print('No velocities file : velocities 0 ; KineticEnergyIons and Temperature',('from the energy file' if enername else '0'))
    if not frcname:
        print('No forces file : forces 0')
    if not enername:
        print('No energy file : potential energy from the E of the positions file, EnergyWithDrift = InternalEnergy (no conserved quantity)')
    if not stressname:
        print('No stress file : Pressure, StressTensor and Enthalpy 0')
    check_kinetic(kinetics,elements,masses)


def check_kinetic(kinetics,elements,masses):#The kinetic energy of the velocities (with the masses used) against Kin. of the energy file : a ratio off 1 means
    #other masses in the run (deuterium, isotopes : -m), or velocities not in bohr/a.u. of time. The masses that fit the energy file are then printed
    if not kinetics:
        return
    data = np.array(kinetics)
    v2, KE = data[:,:-1], data[:,-1]
    ratio = np.median(0.5*(v2*np.array(masses)).sum(axis=1)*AMU_A2FS2/KE)
    print('Check : kinetic energy of the velocities / Kin. of the energy file = %.6f (median of %d frames)' % (ratio,len(KE)))
    if abs(ratio-1.0) <= 1.0E-3:
        return
    print('WARNING : they differ. Either the masses of the run are not those used here (deuterium ? isotopes ? give them with -m), or the')
    print('   velocities are not in bohr/a.u. of time (the default unit of CP2K), or the two files are not from the same run.')
    if len(KE) >= 10*len(elements):
        fitted, residuals, rank, singular = np.linalg.lstsq(0.5*v2*AMU_A2FS2,KE,rcond=None)
        misfit = np.abs(0.5*(v2*fitted).sum(axis=1)*AMU_A2FS2/KE-1.0).max()
        factors = fitted/np.array(masses)
        if rank < len(elements) or misfit > 1.0E-4:
            print('   No set of masses reproduces the kinetic energies of the energy file (misfit %.1e) : check the units and the files' % misfit)
        elif np.ptp(factors) < 1.0E-3*abs(factors.mean()):              #every mass off by the same factor : the unit of the velocities
            print('   The kinetic energies agree when the velocities are multiplied by %.6g : they are in another unit' % np.sqrt(factors.mean()))
            print('   (velocities in A/fs instead of bohr/a.u. of time : 0.0457103) ; write them in the default unit of CP2K')
        elif factors.min() > 0.3 and factors.max() < 4.0:
            print('   Masses that reproduce the kinetic energies of the energy file (within %.0e) : -m' % max(misfit,1e-9),','.join('%s:%.6g' % (el,m) for el,m in zip(elements,fitted)))
        else:
            print('   The masses that would reproduce the kinetic energies of the energy file are not masses of these elements (',
                  ' '.join('%s %.4g' % (el,m) for el,m in zip(elements,fitted)),') : check the units and the files')


def write_snapshot(base,MyCrystal,atommass,record,TimeStep):#One snapshot by umd_processes_fast.print_snapshots
    for key in ('internalenergy','toten','kineticenergy','energywithdrift','enthalpy','temperature','pressure','stress'):
        setattr(MyCrystal,key,record[key])
    cell = record['cell']
    MyCrystal.density = atommass.sum()/6.02214076E23/(abs(np.linalg.det(cell))*1.0E-24)
    MyCrystal.rprimd = cell.tolist()
    MyCrystal.acell = list(np.linalg.norm(cell,axis=1))
    values = record['values']
    for i,atom in enumerate(MyCrystal.atoms):
        atom.xred, atom.xcart, atom.absxcart, atom.vels, atom.forces, atom.charge = (v[i] for v in values)
        atom.magnet = 0.0
    umdpf.print_snapshots(base,MyCrystal,TimeStep,record['time'],values[2])


if __name__ == "__main__":
    main(sys.argv[1:])
