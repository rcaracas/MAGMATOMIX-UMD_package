###
#AUTHORS: RAZVAN CARACAS
###
#Readers of the output files for the validation suite. They are written independently of the package (no import of src/), so that
#they also check the files the package writes. Each returns plain numbers or strings.

import os, re, math

SCALARS = ('InternalEnergy','TotalElectronicEnergy','ElectronicEntropy','KineticEnergyIons','EnergyWithDrift','Enthalpy','Temperature',
           'Pressure','Density')


def number(x):
    try:
        v = float(x)
    except (TypeError,ValueError):
        return x
    return v


def umd_digest(path):#A UMD file read line by line : snapshots, atoms, means of the scalar fields, stress, cell volume, and the atoms of
    #the first and last snapshots (mean unwrapped displacement squared, rms velocity and force)
    natom = None
    nsnap = 0
    sums = {}
    stress = [0.0]*6
    nstress = 0
    volume = 0.0
    nvol = 0
    rprimd = {}
    first = last = None
    lasttime = None
    with open(path) as f:
        for line in f:
            entry = line.split()
            if not entry:
                continue
            key = entry[0]
            if key == 'natom' and natom is None:
                natom = int(entry[1])
            elif key in SCALARS:
                sums[key] = sums.get(key,0.0)+float(entry[1])
            elif key == 'time':
                lasttime = float(entry[1])
            elif key == 'StressTensor':
                stress = [s+float(x) for s,x in zip(stress,entry[1:7])]
                nstress += 1
            elif key in ('rprimd_a','rprimd_b','rprimd_c'):
                rprimd[key] = [float(x) for x in entry[1:4]]
                if len(rprimd) == 3:
                    a, b, c = rprimd['rprimd_a'], rprimd['rprimd_b'], rprimd['rprimd_c']
                    volume += abs(a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0]))
                    nvol += 1
                    rprimd = {}
            elif key == 'atoms:':
                nsnap += 1
                atoms = [[float(x) for x in next(f).split()] for _ in range(natom)]
                if first is None:
                    first = atoms
                last = atoms
    out = {'natom':natom,'snapshots':nsnap}
    if lasttime is not None:
        out['last_time'] = lasttime
    for key in SCALARS:
        if key in sums and nsnap:
            out['mean_'+key] = sums[key]/nsnap
    if nstress:
        for lab, s in zip(('xx','yy','zz','xy','yz','zx'),stress):
            out['mean_stress_'+lab] = s/nstress
    if nvol:
        out['mean_volume'] = volume/nvol
    if first is not None:
        n = len(last)
        out['msd_first_last'] = sum((b[6+k]-a[6+k])**2 for a,b in zip(first,last) for k in range(3))/n
        out['rms_velocity_last'] = math.sqrt(sum(b[9+k]**2 for b in last for k in range(3))/n) if len(last[0]) >= 12 else 'none'
        out['rms_force_last'] = math.sqrt(sum(b[12+k]**2 for b in last for k in range(3))/n) if len(last[0]) >= 15 else 'none'
        out['mean_xcart_last'] = sum(b[3+k] for b in last for k in range(3))/(3*n)
    return out


def rows(path,sep=None):#The lines of a text file split into fields, without the comment lines (#) and the empty ones
    out = []
    with open(path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            out.append([x.strip() for x in (line.rstrip('\n').split(sep) if sep else line.split())])
    return out


def table(path,first,sep=None):#{first column : {column name : value}} of a table whose header line begins with first
    names = None
    out = {}
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            entry = [x.strip() for x in (line.rstrip('\n').split(sep) if sep else line.split())]
            if entry[0].lstrip('#').strip() == first or (entry[0] == '#' and len(entry) > 1 and entry[1] == first):
                names = [x for x in entry if x != '#'] if entry[0] == '#' else entry
                names[0] = names[0].lstrip('#')
                continue
            if names and not line.startswith('#'):
                out[entry[0]] = {n:number(v) for n,v in zip(names,entry)}
    return out


def pairs(path):#files 'A B value' (bond cutoffs, coordinations of analyze_gofr.py) : {'A-B' : value}
    out = {}
    for entry in rows(path):
        if len(entry) == 3:
            out[entry[0]+'-'+entry[1]] = number(entry[2])
    return out


def find(text,pattern,group=1):#a number in a text (the output of a script) ; a pattern with one group
    found = re.search(pattern,text)
    if not found:
        raise ValueError('not found : '+pattern)
    return number(found.group(group))


def findall(text,pattern):
    return [number(x) for x in re.findall(pattern,text)]


def files(folder,pattern):#the files of a folder whose names match a regular expression, sorted
    return sorted(f for f in os.listdir(folder) if re.fullmatch(pattern,f))


def column(path,name,sep=None):#a column of a table with a header line (the first line that is not a comment), as numbers
    lines = rows(path,sep)
    names = lines[0]
    k = names.index(name)
    return [number(r[k]) for r in lines[1:] if len(r) > k]
