#!/usr/bin/env python3
"""
@author: Razvan Caracas, Tim Bögels
Builds new structures - UMD files and VASP POSCAR files - from a snapshot of a UMD file, a POSCAR file or an empty box, by
operations applied in the order of the command line : supercell, shuffle of sub-cells, insertion of molecules, scaling of the cell ;
the atoms without velocities get Maxwell-Boltzmann velocities.
Merges build_supercell.py, box_shuffle.py and insert_umd.py.
"""

import sys, getopt, os
import numpy as np
from scipy.spatial.transform import Rotation
import crystallography as cr
import umd_processes_fast as umdpf


class Structure:#One configuration : cell (rows a, b, c, A), element of each atom, reduced coordinates, velocities (A/fs), whether each atom has a velocity
    #(False for inserted atoms, and for a POSCAR without velocities) ; masses of the elements from the UMD header when read from a UMD file, else from
    #crystallography.Elements2rest
    def __init__(self,cell,symbols,xred,vels,title,massdict=None,hasvel=None):
        self.cell = np.array(cell,dtype=float)
        self.symbols = list(symbols)
        self.xred = np.array(xred,dtype=float).reshape(-1,3) % 1.0
        self.vels = np.array(vels,dtype=float).reshape(-1,3)
        self.title = title
        self.massdict = dict(massdict or {})
        self.hasvel = np.ones(len(self.symbols),dtype=bool) if hasvel is None else np.array(hasvel,dtype=bool)

    def new(self,cell,symbols,xred,vels,hasvel):#another structure with the same title and masses
        return Structure(cell,symbols,xred,vels,self.title,self.massdict,hasvel)

    def mass(self,el):
        if el not in self.massdict:
            self.massdict[el] = cr.Elements2rest(el)[3]
        return self.massdict[el]

    def masses(self):
        return np.array([self.mass(el) for el in self.symbols])

    def widths(self):#distances between opposite faces of the cell
        return abs(np.linalg.det(self.cell))/np.linalg.norm(np.cross(self.cell[[1,2,0]],self.cell[[2,0,1]]),axis=1)


#****Input

def count_snapshots(UMDname):
    with open(UMDname,'r') as f:
        return sum(1 for line in f if line.startswith('atoms:'))


def read_umd(UMDname,first,sampling):#Snapshots first, first+sampling, ... (only first if sampling is 0) ; also the Zelectrons of the elements
    nsnap = count_snapshots(UMDname)
    if first < 0:
        first += nsnap
    if not 0 <= first < nsnap:
        print('the snapshot',first,'is not in',UMDname,'(',nsnap,'snapshots)')
        sys.exit()
    if sampling > 0:
        MyCrystal,Snapshots,TimeStep,length = umdpf.read_values(UMDname,"everything",mode="line",Nsteps=sampling,firststep=first)
        steps = [first+sampling*k for k in range(len(Snapshots))]
    else:
        MyCrystal,Snapshots,TimeStep,length = umdpf.read_values(UMDname,"everything",mode="line",Nsteps=1,firststep=first,laststep=first+1)
        steps = [first]
    cells = umdpf.read_cells(UMDname,len(Snapshots),max(sampling,1),first)
    symbols = [MyCrystal.elements[t] for t in MyCrystal.typat]
    structures = []
    for step,values,cell in zip(steps,Snapshots,cells):
        values = np.array(values).reshape(MyCrystal.natom,12)            #xred, xcart, absxcart, velocities
        structures.append((step,Structure(cell,symbols,values[:,0:3],values[:,9:12],UMDname+' snapshot '+str(step),dict(zip(MyCrystal.elements,MyCrystal.masses)))))
    return structures, dict(zip(MyCrystal.elements,MyCrystal.zelec))


def read_poscar(PoscarFile):#VASP5 POSCAR/CONTCAR (element line), reduced or Cartesian coordinates, optional selective dynamics and Cartesian velocities
    with open(PoscarFile,'r') as ff:
        lines = ff.readlines()
    scale = float(lines[1].split()[0])
    if scale < 0.0:
        print('Negative (volume-based) scaling factors in POSCAR files are not supported.')
        sys.exit()
    cell = np.array([[float(x) for x in lines[2+ii].split()[0:3]] for ii in range(3)])*scale
    elements = lines[5].split()
    if all(x.lstrip('-').isdigit() for x in elements):
        print('The POSCAR file',PoscarFile,'has no element-symbols line (VASP4 format). Please use a VASP5-style POSCAR, such as the ones exported by VESTA.')
        sys.exit()
    counts = [int(x) for x in lines[6].split()]
    natom = sum(counts)
    lineno = 7
    if lines[lineno].strip()[0] in 'sS':                                #selective dynamics
        lineno += 1
    mode = lines[lineno].strip()[0]
    coords = np.array([[float(x) for x in lines[lineno+1+ii].split()[0:3]] for ii in range(natom)])
    xred = coords if mode in 'dD' else (coords*scale) @ np.linalg.inv(cell)
    vels = np.zeros((natom,3))
    hasvel = False
    lineno += 1+natom
    rest = [line for line in lines[lineno:]]
    if len(rest) > natom:                                               #velocity block : a mode line (blank, Cartesian or Direct), then one line per atom
        if rest[0].strip()[:1] in ('d','D'):
            print('WARNING : velocities in direct coordinates in',PoscarFile,'are not read')
        else:
            try:
                vels = np.array([[float(x) for x in rest[1+ii].split()[0:3]] for ii in range(natom)])
                hasvel = True
            except (ValueError,IndexError):
                vels = np.zeros((natom,3))
    symbols = [el for el,n in zip(elements,counts) for _ in range(n)]
    print('Read',natom,'atoms from',PoscarFile,'with elements',elements,'and counts',counts,'' if hasvel else '; no velocities')
    return Structure(cell,symbols,xred,vels,PoscarFile,hasvel=np.full(natom,hasvel))


#****Operations

def supercell(st,n):#Replicates the cell n[0] x n[1] x n[2] times along a, b, c ; the velocities are replicated too
    shifts = np.array([[ix,iy,iz] for ix in range(n[0]) for iy in range(n[1]) for iz in range(n[2])])
    xred = ((st.xred[:,None,:]+shifts[None,:,:])/n).reshape(-1,3)
    symbols = [el for el in st.symbols for _ in range(len(shifts))]
    vels = np.repeat(st.vels,len(shifts),axis=0)
    print('supercell',' x '.join(str(k) for k in n),':',len(symbols),'atoms')
    return st.new(st.cell*np.array(n)[:,None],symbols,xred,vels,np.repeat(st.hasvel,len(shifts)))


def distances(xred,cell):#Distances (A) between all the atoms, nearest images (exact up to half the smallest width of the cell), inf on the diagonal
    d = xred[:,None,:]-xred[None,:,:]
    d -= np.rint(d)
    r = np.linalg.norm(d @ cell,axis=2)
    np.fill_diagonal(r,np.inf)
    return r


def shortest_by_pair(st):#Shortest distance between the atoms of each pair of elements ; inf if no pair
    elements = sorted(set(st.symbols),key=st.symbols.index)
    typat = np.array([elements.index(el) for el in st.symbols])
    r = distances(st.xred,st.cell)
    dmin = np.full((len(elements),len(elements)),np.inf)
    for it in range(len(elements)):
        for jt in range(len(elements)):
            sub = r[np.ix_(typat == it,typat == jt)]
            if sub.size:
                dmin[it,jt] = sub.min()
    return elements, typat, dmin


def push_apart(xred,cell,threshold,masses,rng,maxiter=1000):#Moves apart the pairs of atoms closer than their minimum distance (threshold, natom x natom),
    #each along the line joining them, the lighter one more (by m_other/(m_i+m_j)), until no pair is too close ; returns the positions and the iterations
    inv = np.linalg.inv(cell)
    share = masses[None,:]/(masses[:,None]+masses[None,:])                #part of the missing distance covered by atom i (row) of the pair (i,j)
    x = xred.copy()
    for it in range(maxiter):
        d = x[:,None,:]-x[None,:,:]
        d -= np.rint(d)
        v = d @ cell                                                    #from j to i, nearest image
        r = np.linalg.norm(v,axis=2)
        np.fill_diagonal(r,np.inf)
        missing = threshold-r
        bad = missing > 0.0
        if not bad.any():
            return x, it
        same = bad & (r < 1e-8)                                         #two atoms at the same place : any direction
        if same.any():
            v[same] = rng.standard_normal((same.sum(),3))
            r = np.where(same,np.linalg.norm(v,axis=2),r)
        push = np.zeros_like(r)
        push[bad] = (missing[bad]*1.001+1e-6)*share[bad]/r[bad]       #slightly beyond the minimum, to converge
        x = (x+(push[:,:,None]*v).sum(axis=1) @ inv) % 1.0
    return x, maxiter


def signature(symbols,xred):#The set of (element, position) of a structure, to compare structures
    x = np.round(np.round(xred,6) % 1.0,6) % 1.0
    return frozenset(zip(symbols,map(tuple,x)))


def shuffle(st,arg,factor,mindist,ntries,rng):#Cuts the cell into sub-cells along a, b, c (reduced coordinates : congruent for any cell shape) and permutes them
    #at random with their atoms and velocities ; of ntries permutations the one with the smallest overlaps is kept, then the atoms too close are pushed apart.
    #The permutations that only translate the structure by sub-cells (identity, rigid shifts, swaps of identical copies of a supercell) are rejected : they
    #have no overlaps, so they would always be selected
    widths = st.widths()
    if arg.upper().endswith('A'):                                       #size of the sub-cells in A
        nb = np.array([max(1,int(w//float(arg[:-1]))) for w in widths])
    else:
        nb = np.array([int(v) for v in arg.split(',')]*(3 if ',' not in arg else 1))
    print('shuffle : cell widths (A)',' '.join('%.3f' % w for w in widths),', sub-cells',' x '.join(str(n) for n in nb),
          ', sub-cell widths (A)',' '.join('%.3f' % w for w in widths/nb))
    if nb.prod() < 2:
        print('There is nothing to shuffle with 1 sub-cell, exiting')
        sys.exit()
    elements, typat, original = shortest_by_pair(st)
    if factor is not None:
        threshold = np.where(np.isfinite(original),factor*original,0.0)
        print('   minimum distances (A) :',factor,'x the shortest distances of the structure :')
        for it in range(len(elements)):
            print('     '+'  '.join('%s-%s %.3f' % (elements[it],elements[jt],threshold[it,jt]) for jt in range(it,len(elements)) if threshold[it,jt] > 0.0))
    else:
        threshold = np.full_like(original,mindist)
        print('   minimum distance (A) :',mindist,'; shortest distance before the shuffle :','%.3f' % original.min())
    pairthreshold = np.minimum(threshold[typat][:,typat],distances(st.xred,st.cell)-1e-6)   #no pair pushed beyond its distance before the shuffle :
    #the pairs inside a sub-cell keep their distances (e.g. O-H bonds shorter than -d), only the contacts made by the shuffle are pushed apart
    idx = np.floor(st.xred*nb).astype(int) % nb                         #sub-cell of each atom along a, b, c
    flat = idx[:,0]+nb[0]*(idx[:,1]+nb[1]*idx[:,2])
    shifts = np.array([[ia,ib,ic] for ia in range(nb[0]) for ib in range(nb[1]) for ic in range(nb[2])])/nb
    translated = {signature(st.symbols,st.xred+t) for t in shifts}      #the structure translated by sub-cells
    best = None
    rejected = 0
    for itry in range(ntries):
        new = rng.permutation(nb.prod())[flat]
        newidx = np.stack([new % nb[0],(new//nb[0]) % nb[1],new//(nb[0]*nb[1])],axis=1)
        x = (st.xred+(newidx-idx)/nb) % 1.0
        if signature(st.symbols,x) in translated:
            rejected += 1
            continue
        overlap = np.clip(pairthreshold-distances(x,st.cell),0.0,None).sum()/2.0
        if best is None or overlap < best[0]:
            best = (overlap,x)
    if best is None:
        print('Every permutation tried only translates the structure : the sub-cells are copies of the same cell (after --supercell ?). Use')
        print('smaller sub-cells, or --shuffle before --supercell.')
        sys.exit()
    if rejected:
        print('  ',rejected,'of',ntries,'permutations only translate the structure : rejected')
    overlap, x = best
    moved, iterations = push_apart(x,st.cell,pairthreshold,st.masses(),rng)
    d = moved-x
    d -= np.rint(d)
    disp = np.linalg.norm(d @ st.cell,axis=1)
    nmoved = (disp > 1e-6).sum()
    print('   overlaps of the best of',ntries,'permutations','%.2f' % overlap,'A ;',nmoved,'atoms pushed apart in',iterations,'iterations, displacement max',
          '%.2f' % disp.max(),'A, rms over the moved atoms','%.2f' % (np.sqrt((disp[disp > 1e-6]**2).mean()) if nmoved else 0.0),'A')
    if iterations >= 1000:
        print('WARNING : some atoms are still too close : smaller -c or -d')
    after = shortest_by_pair(st.new(st.cell,st.symbols,moved,st.vels,st.hasvel))[2]
    short = ['%s-%s %.2f (%.2f)' % (elements[it],elements[jt],after[it,jt],original[it,jt]) for it in range(len(elements)) for jt in range(it,len(elements))
             if after[it,jt] < 0.8*original[it,jt]]
    if short:
        print('WARNING : shortest distances (A) below 0.8 x those before the shuffle (in parentheses) :',', '.join(short),': larger -d or -c ?')
    return st.new(st.cell,st.symbols,moved,st.vels,st.hasvel)


def read_molecules(MoleculesFile):#Molecules to insert : blocks 'ntime n' (how many), 'natom m', then m lines 'element x y z' (A)
    molecules = []
    with open(MoleculesFile,'r') as ff:
        lines = [line.split() for line in ff if line.strip()]
    ii = 0
    while ii < len(lines):
        if lines[ii][0] == 'ntime':
            ntime = int(lines[ii][1])
            natom = int(lines[ii+1][1])
            atoms = lines[ii+2:ii+2+natom]
            molecules.append((ntime,[a[0] for a in atoms],np.array([[float(x) for x in a[1:4]] for a in atoms])))
            ii += 2+natom
        else:
            ii += 1
    return molecules


def insert(st,MoleculesFile,rcut,rng,maxtrials=10000):#Inserts the molecules at random positions and orientations (uniform), no atom closer than rcut to
    #any other atom (nearest images) ; the new atoms have zero velocities
    molecules = read_molecules(MoleculesFile)
    inv = np.linalg.inv(st.cell)
    xred, symbols, vels, hasvel = st.xred.copy(), list(st.symbols), st.vels.copy(), st.hasvel.copy()
    for ntime,msymbols,mxcart in molecules:
        masses = np.array([st.mass(el) for el in msymbols])
        centred = mxcart-(masses[:,None]*mxcart).sum(axis=0)/masses.sum()   #rotations about the centre of mass
        for imol in range(ntime):
            for trial in range(1,maxtrials+1):
                new = ((Rotation.random(random_state=rng).apply(centred)) @ inv + rng.random(3)) % 1.0
                if len(xred) == 0:
                    break
                d = new[:,None,:]-xred[None,:,:]
                d -= np.rint(d)
                if np.linalg.norm(d @ st.cell,axis=2).min() >= rcut:
                    break
            else:
                print('I could not insert molecule',imol+1,'of',''.join(msymbols),'after',maxtrials,'trials : there is not enough space in the cell for all the')
                print('molecules with the exclusion radius -r',rcut,'A : increase the cell or decrease -r.')
                sys.exit()
            xred = np.vstack([xred,new])
            symbols += msymbols
            vels = np.vstack([vels,np.zeros((len(msymbols),3))])
            hasvel = np.concatenate([hasvel,np.zeros(len(msymbols),dtype=bool)])
        print('insert :',ntime,'x',''.join(msymbols))
    return st.new(st.cell,symbols,xred,vels,hasvel)


def scale(st,arg,density):#Scales the lengths of the cell by f, or fa, fb, fc along a, b, c, or isotropically to a density (g/cm3) ; reduced coordinates
    #and velocities unchanged, so all the distances, bonds included, scale too
    masses = st.masses()
    if density is not None:
        rho = masses.sum()/6.02214076E23/(abs(np.linalg.det(st.cell))*1E-24)
        f = np.full(3,(rho/density)**(1.0/3.0))
    else:
        f = np.array([float(v) for v in arg.split(',')]*(3 if ',' not in arg else 1))
    cell = st.cell*f[:,None]
    new = st.new(cell,st.symbols,st.xred,st.vels,st.hasvel)
    rho = masses.sum()/6.02214076E23/(abs(np.linalg.det(cell))*1E-24)
    r = distances(new.xred,cell) if len(st.symbols) > 1 else np.array([[np.inf]])
    i,j = np.unravel_index(np.argmin(r),r.shape)
    print('scale : lengths x',' '.join('%.5f' % x for x in f),', volume','%.3f' % abs(np.linalg.det(cell)),'A^3, density','%.4f' % rho,'g/cm3',
          ', shortest distance','%.3f' % r[i,j],'A ('+st.symbols[i]+'-'+st.symbols[j]+')' if np.isfinite(r[i,j]) else '')
    return new


def velocities(st,temperature,redraw,rng):#Maxwell-Boltzmann velocities at temperature for the atoms without velocities (all with redraw) : no net
    #momentum of the drawn atoms, their kinetic temperature set to temperature (3(n-1) degrees of freedom) ; one atom : a plain draw
    draw = np.ones(len(st.symbols),dtype=bool) if redraw else ~st.hasvel
    n = draw.sum()
    if n == 0:
        return st
    masses = st.masses()[draw]
    v = rng.standard_normal((n,3))*np.sqrt(8.617333262E-05*temperature/(masses*1.0364269E-4))[:,None]     #sqrt(kB T/m) in A/fs
    if n > 1:
        v -= (masses[:,None]*v).sum(axis=0)/masses.sum()
        ekin = 0.5*(masses[:,None]*v**2).sum()*1.0364269E-4
        v *= np.sqrt(1.5*(n-1)*8.617333262E-05*temperature/ekin)
    vels = st.vels.copy()
    vels[draw] = v
    print('velocities : Maxwell-Boltzmann at',temperature,'K for',n,'atoms'+(' (all)' if redraw else ' without velocities'))
    return st.new(st.cell,st.symbols,st.xred,vels,np.ones(len(st.symbols),dtype=bool))


#****Output

def ordered(st):#Atoms grouped by element (VASP and UMD), in the order of first appearance ; the order inside each element is kept
    elements = sorted(set(st.symbols),key=st.symbols.index)
    order = sorted(range(len(st.symbols)),key=lambda i: elements.index(st.symbols[i]))
    counts = [st.symbols.count(el) for el in elements]
    return elements, counts, np.array(order,dtype=int)


def write_poscar(filename,st):#Direct coordinates, Cartesian velocities in A/fs
    elements, counts, order = ordered(st)
    with open(filename,'w') as ff:
        ff.write(st.title+'\n')
        ff.write('  1.0\n')
        for vec in st.cell:
            ff.write('  {:.10f}  {:.10f}  {:.10f}\n'.format(*vec))
        ff.write('  '+'  '.join(elements)+'\n')
        ff.write('  '+'  '.join(str(n) for n in counts)+'\n')
        ff.write('Direct\n')
        for x in st.xred[order]:
            ff.write('  {:.10f}  {:.10f}  {:.10f}\n'.format(*x))
        ff.write('Cartesian\n')
        for v in st.vels[order]:
            ff.write('  {:.10e}  {:.10e}  {:.10e}\n'.format(*v))


def write_xyz(filename,st):
    elements, counts, order = ordered(st)
    xcart = st.xred @ st.cell
    with open(filename,'w') as ff:
        ff.write(str(len(st.symbols))+'\n'+st.title+'\n')
        for i in order:
            ff.write(st.symbols[i]+'  {:.6f}  {:.6f}  {:.6f}\n'.format(*xcart[i]))


def write_umd(filename,st,zelec):#One snapshot ; positions, velocities, cell, density and kinetic temperature, the other fields 0
    elements, counts, order = ordered(st)
    MyCrystal = cr.Lattice()
    MyCrystal.natom = len(st.symbols)
    MyCrystal.ntypat = len(elements)
    MyCrystal.types = counts
    MyCrystal.elements = elements
    MyCrystal.masses = [st.mass(el) for el in elements]
    MyCrystal.zelec = [zelec.get(el,0.0) for el in elements]
    MyCrystal.noelectrons = sum(z*n for z,n in zip(MyCrystal.zelec,counts))
    MyCrystal.typat = [elements.index(st.symbols[i]) for i in order]
    MyCrystal.acell = list(np.linalg.norm(st.cell,axis=1))
    MyCrystal.rprimd = [list(vec) for vec in st.cell]
    masses = st.masses()
    MyCrystal.density = masses.sum()/6.02214076E23/(abs(np.linalg.det(st.cell))*1E-24)
    ekin = 0.5*(masses[:,None]*st.vels**2).sum()*1.0364269E-4         #amu A^2/fs^2 -> eV
    MyCrystal.kineticenergy = ekin
    MyCrystal.temperature = 2.0*ekin/(3.0*len(masses)*8.617333262E-05)
    xcart = st.xred @ st.cell
    MyCrystal.atoms = [cr.Atom() for _ in range(MyCrystal.natom)]
    for inew,i in enumerate(order):
        MyCrystal.atoms[inew].xred = list(st.xred[i])
        MyCrystal.atoms[inew].xcart = list(xcart[i])
        MyCrystal.atoms[inew].vels = list(st.vels[i])
    base = filename[:-8]
    open(filename,'w').close()                                          #print_header and print_snapshots append
    umdpf.print_header(base,MyCrystal)
    umdpf.print_snapshots(base,MyCrystal,1.0,0.0,xcart[order])


def main(argv):
    infile = ''
    first = -1
    sampling = 0
    box = None
    operations = []
    factor = None
    mindist = 1.0
    ntries = 100
    rcut = 2.0
    temperature = 1500.0
    redraw = False
    nconfig = 1
    seed = None
    outputs = ['umd','poscar']
    name = ''
    umdpf.headerumd()
    usage = ('build_umd.py -f <UMD or POSCAR file> -i <snapshot> -s <sampling> -a <empty box> --supercell <na,nb,nc> --shuffle <sub-cells> '
             '--insert <molecules file> --scale <factors> --density <g/cm3> -d <min distance> -c <factor> --tries <tries> -r <cutoff> -t <temperature> '
             '--velocities all -n <configurations> --seed <seed> -o <outputs> --name <name>')
    try:
        opts, arg = getopt.getopt(argv,"hf:i:s:a:c:d:t:r:n:o:",["supercell=","shuffle=","insert=","scale=","density=","tries=","velocities=","seed=","name="])
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('build_umd.py : builds new structures (UMD, VASP POSCAR or xyz files) from a snapshot of a UMD file, a POSCAR file or an empty box,')
            print('by the operations --supercell, --shuffle, --insert, --scale and --density, applied in the order of the command line ; each can be')
            print('repeated. The atoms without velocities then get Maxwell-Boltzmann velocities (-t). Usage :')
            print(usage)
            print('Input :')
            print('-f : a UMD file (*.umd.dat) or a VASP5 POSCAR/CONTCAR file (element line ; reduced or Cartesian coordinates ; Cartesian velocities read).')
            print('-i : with a UMD file, the snapshot, counted from 0 ; negative from the end. Default -1 (the last).')
            print('-s : with a UMD file, every s-th snapshot from -i on (one structure each). Default 0 : only -i.')
            print('-a : an empty box instead of -f : a cube of side a, or a,b,c for an orthorhombic box (A).')
            print('Operations :')
            print('--supercell na,nb,nc : repeats the cell na x nb x nc times along a, b, c (velocities repeated too).')
            print('--shuffle n, or na,nb,nc, or <size>A (e.g. 4A) : cuts the cell into sub-cells along a, b, c (n per direction, or the widths of the')
            print('     cell divided by size) and permutes them at random with their atoms and velocities : the local structure inside the sub-cells is')
            print('     kept, the global configuration is randomized. Every permutation brings atoms too close at the interfaces of the sub-cells : of')
            print('     --tries random permutations (default 100) the one with the smallest overlaps is kept, then the atoms closer than -d (A, default')
            print('     1.0) are pushed apart along the line joining them, the lighter one more ; no pair is pushed beyond its distance before the')
            print('     shuffle, so the bonds inside the sub-cells are kept (e.g. O-H shorter than -d) and only the contacts made by the shuffle are')
            print('     pushed. -c f instead of -d : minimum distance of each pair of elements = f x its shortest distance before the shuffle (large')
            print('     for rare elements).')
            print('     Keeping the least overlapping permutation favours those that move fewer interfaces (pyrolite 2x2x2 : half as many atoms')
            print('     pushed with --tries 100 as with 1) ; --tries 1 gives a plain random permutation.')
            print('     Permutations that only translate the structure by sub-cells (e.g. swaps of identical copies after --supercell) are rejected.')
            print('--insert <file> : inserts molecules at random positions and orientations, no atom closer than -r (A, default 2.0) to any other')
            print('     (nearest images). The file has blocks "ntime n", "natom m", then m lines "element x y z" (A) ; e.g. data/molecules.dat.')
            print('     Choose -r for the system : it applies to all the pairs of a new atom and another atom.')
            print('--scale f, or fa,fb,fc : multiplies the lengths of the cell by f (or along a, b, c) ; the reduced coordinates and the velocities are')
            print('     kept, so all the distances, bonds included, scale too. The new density and the shortest distance are printed.')
            print('--density rho : scales the cell isotropically to the density rho (g/cm3), as --scale.')
            print('Velocities :')
            print('-t : temperature (K) of the Maxwell-Boltzmann velocities given to the atoms without velocities (inserted atoms, atoms of a POSCAR')
            print('     without velocities or of an empty box), with no net momentum and exactly this kinetic temperature. Default 1500.')
            print('--velocities all : new velocities for all the atoms (those of the snapshot are replaced).')
            print('Output :')
            print('-n : number of structures built from each input structure (the random operations differ). Default 1.')
            print('--seed : seed of the random numbers, to reproduce a build. Default : random.')
            print('-o : output formats, separated by commas : umd, poscar, xyz. Default umd,poscar.')
            print('--name : base of the output names. Default : the input file without .umd.dat or extension, or "box".')
            print('Files <name>[.<snapshot>].build<k>.umd.dat, .POSCAR, .xyz (k = 1 .. n). The atoms are grouped by element, in the order of first')
            print('appearance. UMD : one snapshot, with the cell, positions, velocities, density and kinetic temperature (the other fields are 0).')
            print('POSCAR : direct coordinates and Cartesian velocities (A/fs).')
            sys.exit()
        elif opt == '-f':
            infile = arg
        elif opt == '-i':
            first = int(arg)
        elif opt == '-s':
            sampling = int(arg)
        elif opt == '-a':
            box = [float(v) for v in arg.split(',')]
            box = box*3 if len(box) == 1 else box
        elif opt in ('--supercell','--shuffle','--insert','--scale','--density'):
            operations.append((opt[2:],arg))
        elif opt == '-c':
            factor = float(arg)
        elif opt == '-d':
            mindist = float(arg)
        elif opt == '--tries':
            ntries = int(arg)
        elif opt == '-t':
            temperature = float(arg)
        elif opt == '--velocities':
            if arg != 'all':
                print('--velocities takes only "all"')
                sys.exit()
            redraw = True
        elif opt == '-r':
            rcut = float(arg)
        elif opt == '-n':
            nconfig = int(arg)
        elif opt == '--seed':
            seed = int(arg)
        elif opt == '-o':
            outputs = arg.split(',')
        elif opt == '--name':
            name = arg

    zelec = {}
    if box is not None:
        inputs = [(None,Structure(np.diag(box),[],np.zeros((0,3)),np.zeros((0,3)),'empty box'))]
        base = name or 'box'
    elif not os.path.isfile(infile):
        print('the input file',infile,'does not exist : give a UMD or POSCAR file with -f, or an empty box with -a')
        sys.exit()
    else:
        with open(infile,'r') as f:
            isumd = f.readline().split()[:1] == ['natom']
        if isumd:
            inputs, zelec = read_umd(infile,first,sampling)
            base = name or (infile[:-8] if infile.endswith('.umd.dat') else infile)
        else:
            inputs = [(None,read_poscar(infile))]
            base = name or os.path.splitext(infile)[0]
    for step,st in inputs:
        if st.vels.size and not st.vels.any():
            print('WARNING : all the velocities of',st.title,'are 0 (first snapshot of VaspParser.py ?) : take another snapshot or --velocities all')
    for fmt in outputs:
        if fmt not in ('umd','poscar','xyz'):
            print('Unknown output format',fmt,': use umd, poscar, xyz')
            sys.exit()

    rng = np.random.default_rng(seed)
    for step,st0 in inputs:
        for k in range(1,nconfig+1):
            st = st0
            print('\n'+st.title+', structure',k,':',len(st.symbols),'atoms')
            for op,arg in operations:
                if op == 'supercell':
                    n = [int(v) for v in arg.split(',')]
                    st = supercell(st,n*3 if len(n) == 1 else n)
                elif op == 'shuffle':
                    st = shuffle(st,arg,factor,mindist,ntries,rng)
                elif op == 'insert':
                    st = insert(st,arg,rcut,rng)
                elif op == 'scale':
                    st = scale(st,arg,None)
                elif op == 'density':
                    st = scale(st,None,float(arg))
            st = velocities(st,temperature,redraw,rng)
            st.title = st.title+(' '+' '.join('--'+op+' '+arg for op,arg in operations) if operations else '')
            filebase = base+('.'+str(step) if step is not None else '')+'.build'+str(k)
            for fmt in outputs:
                if fmt == 'umd':
                    filename = filebase+'.umd.dat'
                    write_umd(filename,st,zelec)
                elif fmt == 'poscar':
                    filename = filebase+'.POSCAR'
                    write_poscar(filename,st)
                else:
                    filename = filebase+'.xyz'
                    write_xyz(filename,st)
                print('   written in',filename)


if __name__ == "__main__":
    main(sys.argv[1:])
