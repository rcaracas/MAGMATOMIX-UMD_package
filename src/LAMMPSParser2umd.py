#!/usr/bin/env python3
###
#AUTHORS: KEVIN JIGUET-COVEX, RAZVAN CARACAS
###
#Converts a LAMMPS trajectory (dump custom/atom text file, or dump xyz) and its log file into a UMD file. The positions, velocities, forces and
#charges come from the dump, the cell from the box of every frame (any triclinic box, NPT), the energies, temperature, pressure and stress
#tensor from the thermo output of the log, matched to the frames by their step. Units metal or real, converted to A, fs, eV, GPa.

import sys, getopt, os.path, re
import numpy as np
import crystallography as cr
import umd_processes_fast as umdpf
from vibrations_umd import KB, AMU_A2FS2

KCALMOL = 0.0433641043                  #1 kcal/mol in eV
UNITS = {'metal':{'time':1000.0, 'energy':1.0, 'velocity':1.0E-3, 'force':1.0, 'pressure':1.0E-4, 'timestep':0.001},     #ps, eV, A/ps, eV/A, bar
         'real':{'time':1.0, 'energy':KCALMOL, 'velocity':1.0, 'force':KCALMOL, 'pressure':1.01325E-4, 'timestep':1.0}}  #fs, kcal/mol, A/fs, atm
STRESS = ('Pxx','Pyy','Pzz','Pxy','Pyz','Pxz')     #UMD order xx yy zz xy yz zx


#****Log file

def numbers(tokens):
    try:
        return [float(x) for x in tokens]
    except ValueError:
        return None


def read_log(logname):#Units, time step, masses, elements, box and thermo output (step -> {keyword : value}) of a LAMMPS log file
    log = {'units':None, 'timestep':None, 'masses':{}, 'elements':None, 'box':None, 'thermo':{}, 'columns':set(), 'repeated':0}
    header = None
    multistep = None
    with open(logname,'r') as f:
        for line in f:
            tokens = line.split()
            if not tokens:
                multistep = None
                continue
            if tokens[0] != 'pair_coeff' and re.search(r'pair_coeff\s+\*\s+\*',line):#the info of pair_style deepmd ends without a newline :
                tokens = line[line.index('pair_coeff'):].split()                    #the pair_coeff command follows on the same line
            if tokens[0] == 'units' and len(tokens) > 1:
                log['units'] = tokens[1]
            elif tokens[0] == 'timestep' and len(tokens) > 1 and numbers(tokens[1:2]):
                log['timestep'] = float(tokens[1])
            elif tokens[0] == 'mass' and len(tokens) > 2 and tokens[1].isdigit() and numbers(tokens[2:3]):
                log['masses'][int(tokens[1])] = float(tokens[2])
            elif tokens[0] == 'dump_modify' and 'element' in tokens:
                log['elements'] = [x for x in tokens[tokens.index('element')+1:] if cr.element_entry(x) is not None]
            elif tokens[0] == 'pair_coeff' and tokens[1:3] == ['*','*'] and log['elements'] is None:
                symbols = []
                for x in reversed(tokens[3:]):              #the element names at the end (eam/alloy, meam, tersoff, deepmd, ...)
                    if cr.element_entry(x) is None or x.isdigit():
                        break
                    symbols.insert(0,x)
                if symbols:
                    log['pair_elements'] = symbols
            elif (tokens[0] in ('orthogonal','triclinic') and len(tokens) > 2 and tokens[1] == 'box') and log['box'] is None:
                values = [float(x) for x in re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?',line)]
                if len(values) >= 6:
                    lo, hi = np.array(values[0:3]), np.array(values[3:6])
                    xy, xz, yz = values[6:9] if tokens[0] == 'triclinic' and len(values) >= 9 else (0.0,0.0,0.0)
                    log['box'] = (np.array([[hi[0]-lo[0],0.0,0.0],[xy,hi[1]-lo[1],0.0],[xz,yz,hi[2]-lo[2]]]), lo)
            elif tokens[0] == 'Step' and numbers(tokens[1:]) is None:            #thermo_style one / custom : a header, then one row per output step
                header = tokens
                continue
            elif tokens[0].startswith('-----') and len(tokens) > 2 and tokens[1] == 'Step':         #thermo_style multi : a block per step
                multistep = int(tokens[2].rstrip(','))
                header = None
                add_thermo(log,multistep,{})
                continue
            elif tokens[0] == 'Loop':
                header = multistep = None
                continue
            if header is not None and len(tokens) == len(header):
                values = numbers(tokens)
                if values is not None:
                    add_thermo(log,int(values[0]),dict(zip(header[1:],values[1:])))
            elif multistep is not None:
                pairs = re.findall(r'(\S+)\s*=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)',line)
                log['thermo'][multistep].update({k:float(v) for k,v in pairs})
                log['columns'].update(k for k,v in pairs)
    if log['elements'] is None and log.get('pair_elements'):
        log['elements'] = log['pair_elements']
    return log


def add_thermo(log,step,values):#The first row of a run repeats the last step of the previous one ; a later run reusing step numbers overwrites them
    if step in log['thermo'] and values and any(abs(log['thermo'][step].get(k,v)-v) > 1e-8*max(1.0,abs(v)) for k,v in values.items()):
        log['repeated'] += 1
    log['thermo'][step] = dict(values)
    log['columns'].update(values)


def read_stress_file(stressname):#Lines 'step pxx pyy pzz pxy pxz pyz' (fix ave/time output, LAMMPS order and pressure units ; '#' comments)
    stress = {}
    with open(stressname,'r') as f:
        for line in f:
            values = numbers(line.split()[:7]) if not line.startswith('#') else None
            if values is not None and len(values) == 7:
                pxx, pyy, pzz, pxy, pxz, pyz = values[1:]
                stress[int(values[0])] = [pxx,pyy,pzz,pxy,pyz,pxz]
    return stress


#****Dump file

def read_box(f,line):#ITEM: BOX BOUNDS [xy xz yz | abc origin] ... : cell (rows a, b, c) and origin
    rows = [numbers(next(f).split()) for _ in range(3)]
    words = line.split()
    if 'abc' in words:                                                  #general triclinic : a, b, c and the origin
        return np.array([r[0:3] for r in rows]), np.array([r[3] for r in rows])
    if 'xy' in words:                                                   #restricted triclinic : bounds of the tilted box, and tilts
        (xlob,xhib,xy), (ylob,yhib,xz), (zlo,zhi,yz) = rows
        xlo, xhi = xlob-min(0.0,xy,xz,xy+xz), xhib-max(0.0,xy,xz,xy+xz)
        ylo, yhi = ylob-min(0.0,yz), yhib-max(0.0,yz)
        return np.array([[xhi-xlo,0.0,0.0],[xy,yhi-ylo,0.0],[xz,yz,zhi-zlo]]), np.array([xlo,ylo,zlo])
    (xlo,xhi), (ylo,yhi), (zlo,zhi) = [r[0:2] for r in rows]
    return np.diag([xhi-xlo,yhi-ylo,zhi-zlo]), np.array([xlo,ylo,zlo])


def read_dump(dumpname,logbox):#Yields the frames : step, time (or None), cell, origin, column names, array natom x columns of strings
    with open(dumpname,'r') as f:
        first = f.readline()
        f.seek(0)
        if first.startswith('ITEM:'):
            step = time = cell = origin = None
            for line in f:
                if line.startswith('ITEM: TIMESTEP'):
                    step = int(next(f).split()[0])
                    time = None
                elif line.startswith('ITEM: TIME'):
                    time = float(next(f).split()[0])
                elif line.startswith('ITEM: NUMBER OF ATOMS'):
                    natom = int(next(f).split()[0])
                elif line.startswith('ITEM: BOX BOUNDS'):
                    cell, origin = read_box(f,line)
                elif line.startswith('ITEM: ATOMS'):
                    columns = line.split()[2:]
                    rows = np.array(' '.join(next(f) for _ in range(natom)).split()).reshape(natom,len(columns))
                    yield step, time, cell, origin, columns, rows
        else:                                                           #dump xyz : number of atoms, 'Atoms. Timestep: N', then 'type x y z'
            if logbox is None:
                print('The dump xyz format has no box : it is read from the log file (orthogonal/triclinic box = ...), which has none')
                sys.exit()
            while True:
                line = f.readline()
                if not line.strip():
                    break
                natom = int(line.split()[0])
                comment = f.readline().split()
                step = int(comment[-1]) if comment and comment[-1].lstrip('-').isdigit() else None
                rows = np.array(' '.join(f.readline() for _ in range(natom)).split()).reshape(natom,-1)[:,0:4]
                yield step, None, logbox[0], logbox[1], ['type','x','y','z'], rows


def columns_of(columns,names):
    return [columns.index(n) for n in names] if all(n in columns for n in names) else None


#****Main

def main(argv):
    umdpf.headerumd()
    dumpname = ''
    logname = ''
    stressname = ''
    outname = None
    units = None
    timestep = None
    elements_opt = None
    InitialStep = 0
    SkipStep = 1
    try:
        opts, arg = getopt.getopt(argv,"hf:l:a:o:u:t:e:i:s:",["fFile","lLog","aAllpress"])
    except getopt.GetoptError:
        print('LAMMPSParser2umd.py -f <dump file> -l <log file> -a <stress file> -o <output> -u <units> -t <timestep> -e <elements> -i <InitialStep> -s <SkipStep>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('LAMMPSParser2umd.py : converts a LAMMPS trajectory and its log file into a UMD file')
            print('LAMMPSParser2umd.py -f <dump file> -l <log file> -a <stress file> -o <output> -u <units> -t <timestep> -e <elements> -i <InitialStep> -s <SkipStep>')
            print('-f : dump file, text : dump custom/atom (ITEM: lines ; any box, also NPT) or dump xyz (box from the log). Columns read : id, type,')
            print('     element, x y z or xs ys zs (wrapped), xu yu zu or xsu ysu zsu (unwrapped), ix iy iz (image flags), vx vy vz, fx fy fz, q.')
            print('     Unwrapped positions from xu/xsu, else from the image flags, else by following the atoms from frame to frame (as')
            print('     VaspParser.py ; the dump must then be frequent enough that no atom moves half a cell between frames).')
            print('     The atoms are sorted by id and grouped by element (by id within an element), the same order in every frame.')
            print('-l : log file (optional) : units, timestep, mass, dump_modify ... element, pair_coeff * * ... <elements>, box, and the thermo')
            print('     output (thermo_style one, custom or multi), matched to the frames by their step. Read : Temp, Press, PotEng, KinEng,')
            print('     TotEng, Econserve, Enthalpy, Pxx Pyy Pzz Pxy Pxz Pyz. InternalEnergy = TotEng (potential + kinetic, as VaspParser.py),')
            print('     TotalElectronicEnergy = PotEng, KineticEnergyIons = KinEng (else TotEng - PotEng), EnergyWithDrift = Econserve')
            print('     (TotEng + thermostat energy) or TotEng.')
            print('     Without them : kinetic energy and temperature (3N-3 degrees of freedom) from the velocities, the other fields 0.')
            print('-a : stress file (optional ; replaces Pxx... of the log) : lines step pxx pyy pzz pxy pxz pyz, e.g. the output of')
            print('     fix ave/time 1 1 N c_thermo_press[1] c_thermo_press[2] ... c_thermo_press[6] file <stress file>')
            print('-o : name of the UMD file (default <dump file>.umd.dat)')
            print('-u : units, metal or real (default : the units of the log, else metal)')
            print('-t : time step of the MD in fs (default : the timestep of the log, else the LAMMPS default of the units)')
            print('-e : element of each atom type, separated by commas (type 1, 2, ... ; default : the element column of the dump, else')
            print('     dump_modify element or pair_coeff * * of the log, else the element of the mass of each type)')
            print('-i : first frame written (default 0) ; -s : every s-th frame (default 1)')
            sys.exit()
        elif opt in ('-f','--fFile'):
            dumpname = str(arg)
        elif opt in ('-l','--lLog'):
            logname = str(arg)
        elif opt in ('-a','--aAllpress'):
            stressname = str(arg)
        elif opt == '-o':
            outname = str(arg)
        elif opt == '-u':
            units = str(arg)
        elif opt == '-t':
            timestep = float(arg)
        elif opt == '-e':
            elements_opt = [x.strip() for x in arg.split(',')]
        elif opt == '-i':
            InitialStep = int(arg)
        elif opt == '-s':
            SkipStep = int(arg)
    if not os.path.isfile(dumpname):
        print('The dump file',dumpname,'does not exist')
        sys.exit()
    for name in (logname,stressname):
        if name and not os.path.isfile(name):
            print('The file',name,'does not exist')
            sys.exit()
    if SkipStep < 1 or InitialStep < 0:
        print('-s must be at least 1 and -i at least 0')
        sys.exit()
    log = read_log(logname) if logname else {'units':None,'timestep':None,'masses':{},'elements':None,'box':None,'thermo':{},'columns':set(),'repeated':0}
    units = units or log['units'] or 'metal'
    if units not in UNITS:
        print('Units',units,'are not supported : metal or real')
        sys.exit()
    U = UNITS[units]
    dt = timestep if timestep is not None else (log['timestep'] if log['timestep'] is not None else U['timestep'])*U['time']     #fs
    stressfile = read_stress_file(stressname) if stressname else {}
    if outname is None:
        outname = dumpname+'.umd.dat'
    if not outname.endswith('.umd.dat'):
        print('The name of the UMD file must end with .umd.dat')
        sys.exit()
    base = outname[:-8]
    print('Dump file',dumpname,'; log file',logname if logname else 'none','; units',units,'; time step %g fs' % dt)
    if log['repeated']:
        print('WARNING :',log['repeated'],'steps of the thermo output appear in several runs with different values (reset_timestep ?) : the last ones are used')

    MyCrystal = None
    missing = set()
    nframes = nwritten = nothermo = 0
    previous = None
    for step, time, cell, origin, columns, rows in read_dump(dumpname,log['box']):
        natom = rows.shape[0]
        if 'id' in columns:
            byid = np.argsort(rows[:,columns.index('id')].astype(int),kind='stable')
            rows = rows[byid]
        if MyCrystal is None:                                           #first frame : elements, masses, order of the atoms
            symbols = None
            if 'type' in columns:
                types = rows[:,columns.index('type')]
                if not all(t.isdigit() for t in types):                 #dump xyz with dump_modify element : symbols in the type column
                    symbols = list(types)
                    types = None
                else:
                    types = types.astype(int)
            else:
                types = None
            if 'element' in columns and elements_opt is None:
                symbols = list(rows[:,columns.index('element')])
            elif types is not None:
                ntype = types.max()
                if elements_opt is not None:
                    names, source = elements_opt, '-e'
                elif log['elements'] is not None:
                    names, source = log['elements'], 'log file'
                elif len(log['masses']) >= ntype:
                    names, source = [], 'masses of the log file'
                    for t in range(1,ntype+1):
                        z = min(range(1,104),key=lambda z: abs(cr.Elements2rest(z)[3]-log['masses'][t]))
                        names.append(cr.Elements2rest(z)[1])
                else:
                    print('No element names : give them with -e (type 1, 2, ...), or add an element column to the dump (dump_modify element)')
                    sys.exit()
                if len(names) < ntype:
                    print('The',source,'gives',len(names),'elements for',ntype,'atom types : give them with -e')
                    sys.exit()
                symbols = [names[t-1] for t in types]
                print('elements from the',source,':',' '.join(str(t+1)+'='+n for t,n in enumerate(names[:ntype])))
            elif symbols is None:
                print('No type or element column in the dump')
                sys.exit()
            if any(cr.element_entry(s) is None for s in set(symbols)):
                print('Unknown elements',sorted(set(s for s in symbols if cr.element_entry(s) is None)),': give them with -e')
                sys.exit()
            elements = list(dict.fromkeys(symbols))                     #in the order of first appearance
            elindex = np.array([elements.index(s) for s in symbols])
            order = np.lexsort((np.arange(natom),elindex))              #grouped by element, by id within an element
            masses = []
            for el in elements:
                mass = cr.Elements2rest(el)[3]
                if types is not None and log['masses']:
                    found = sorted(set(log['masses'].get(int(t)) for t,s in zip(types,symbols) if s == el) - {None})
                    if len(found) == 1:
                        mass = found[0]
                masses.append(mass)
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
            atommass = np.array([masses[k] for k in elindex[order]])
            ref_symbols = [symbols[i] for i in order]
            open(outname,'w').close()
            umdpf.print_header(base,MyCrystal,'LAMMPS : the UMD file originated from a LAMMPS molecular dynamics simulation')
            colpos = {k:columns_of(columns,v) for k,v in (('x',('x','y','z')),('xs',('xs','ys','zs')),('xu',('xu','yu','zu')),('xsu',('xsu','ysu','zsu')),
                                                          ('i',('ix','iy','iz')),('v',('vx','vy','vz')),('f',('fx','fy','fz')))}
            if not any(colpos[k] for k in ('x','xs','xu','xsu')):
                print('No positions in the dump (x y z, xs ys zs, xu yu zu or xsu ysu zsu)')
                sys.exit()
            unwrap = ('from xu yu zu' if colpos['xu'] else 'from xsu ysu zsu') if colpos['xu'] or colpos['xsu'] else ('from the image flags' if colpos['i'] else 'followed from frame to frame')
            print(natom,'atoms :',' '.join(str(n)+' '+el for n,el in zip(MyCrystal.types,elements)),'; unwrapped positions :',unwrap)
        if natom != MyCrystal.natom:
            print('Frame at step',step,'has',natom,'atoms instead of',MyCrystal.natom)
            sys.exit()
        rows = rows[order]
        if 'element' in columns and elements_opt is None and list(rows[:,columns.index('element')]) != ref_symbols:
            print('The elements of the atoms change at step',step)
            sys.exit()
        inverse = np.linalg.inv(cell)
        def floats(key):
            return rows[:,colpos[key]].astype(float)
        if colpos['x']:
            s = (floats('x')-origin) @ inverse
        elif colpos['xs']:
            s = floats('xs')
        else:
            s = (floats('xu')-origin) @ inverse if colpos['xu'] else floats('xsu')
        xred = s - np.floor(s)
        xcart = xred @ cell
        if colpos['xu']:
            absx = floats('xu')-origin
        elif colpos['xsu']:
            absx = floats('xsu') @ cell
        elif colpos['i']:
            absx = (s + floats('i')) @ cell
        elif previous is None:
            absx = xcart.copy()
        else:                                                           #nearest move between frames, in reduced coordinates of the current cell
            d = xred - previous
            absx = previousabs + (d - np.rint(d)) @ cell
        previous, previousabs = xred, absx
        if time is None:
            time = (step if step is not None else nframes)*dt
        else:
            time *= U['time']
        time = round(time,6)                                            #fs ; no rounding noise of the unit conversion in the timestep
        nframes += 1
        if nframes-1 < InitialStep or (nframes-1-InitialStep) % SkipStep:
            continue
        vels = floats('v')*U['velocity'] if colpos['v'] else np.zeros((natom,3))
        forces = floats('f')*U['force'] if colpos['f'] else np.zeros((natom,3))
        charges = rows[:,columns.index('q')].astype(float) if 'q' in columns else np.zeros(natom)
        thermo = log['thermo'].get(step,{})
        if logname and not thermo:
            nothermo += 1
        def get(key,factor):
            if key in thermo:
                return thermo[key]*factor
            missing.add(key)
            return None
        kinetic = 0.5*(atommass[:,None]*vels**2).sum()*AMU_A2FS2
        KE = get('KinEng',U['energy'])
        PE = get('PotEng',U['energy'])
        TE = get('TotEng',U['energy'])
        if KE is None:
            KE = TE - PE if TE is not None and PE is not None else kinetic   #TotEng - PotEng = KinEng of LAMMPS, else the velocities
        if TE is None and PE is not None:
            TE = PE + KE
        if PE is None and TE is not None:
            PE = TE - KE
        T = get('Temp',1.0)
        if T is None:
            T = 2.0*kinetic/(max(3*natom-3,1)*KB)
        stress = stressfile.get(step)
        if stress is None:
            stress = [thermo[k] for k in STRESS] if all(k in thermo for k in STRESS) else None
            if stress is None:
                missing.add('Pxx...')
        P = get('Press',U['pressure'])
        if P is None and stress is not None:
            P = sum(stress[0:3])/3.0*U['pressure']
        econ = get('Econserve',U['energy'])
        enth = get('Enthalpy',U['energy'])
        record = {'time':time, 'cell':cell, 'internalenergy':TE if TE is not None else 0.0, 'toten':PE if PE is not None else 0.0,
                  'kineticenergy':KE, 'energywithdrift':econ if econ is not None else (TE if TE is not None else 0.0),
                  'enthalpy':enth if enth is not None else 0.0, 'temperature':T, 'pressure':P if P is not None else 0.0,
                  'stress':[x*U['pressure'] for x in stress] if stress is not None else [0.0 for _ in range(6)],
                  'values':[xred.tolist(), xcart.tolist(), absx.tolist(), vels.tolist(), forces.tolist(), charges.tolist()]}
        if nwritten == 0:                                               #the first snapshot waits for the second : its timestep is the time between them
            first = record
        else:
            if nwritten == 1:
                write_snapshot(base,MyCrystal,atommass,first,round(time-first['time'],6))
            write_snapshot(base,MyCrystal,atommass,record,round(time-lasttime,6))
        lasttime = time
        nwritten += 1
    if nwritten == 1:
        write_snapshot(base,MyCrystal,atommass,first,dt)
    if MyCrystal is None:
        print('No frame in the dump file',dumpname)
        sys.exit()
    print(nframes,'frames read,',nwritten,'written in',outname)
    if nothermo:
        print('WARNING :',nothermo,'of the',nwritten,'frames written have no thermo output at their step (thermo and dump every N steps differ ?) :')
        print('   their energies, pressure and stress are 0 and their temperature is that of the velocities')
    if logname:
        notfound = sorted(missing & {'KinEng','PotEng','TotEng','Temp','Press','Pxx...'})
        if notfound:
            print('Not found in the thermo output (or no thermo output at the step of some frames) :',', '.join(notfound))
            print('   KineticEnergyIons = TotEng - PotEng, else computed from the velocities as the Temperature when missing ; the other fields 0')
    else:
        print('No log file : KineticEnergyIons and Temperature from the velocities, energies, pressure and stress 0')


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
