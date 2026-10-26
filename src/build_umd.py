#!/usr/bin/env python3
"""
@author: Razvan Caracas, Tim Bögels
Builds new structures - UMD files and VASP POSCAR files - from a snapshot of a UMD file, a POSCAR file or an empty box, by
operations applied in the order of the command line : supercell, shuffle of sub-cells, insertion of molecules, scaling of the cell,
replacement of atoms, removal of molecules ; the atoms without velocities get Maxwell-Boltzmann velocities.
--extract writes instead the trajectory of a selection of atoms (formerly extratom_umd.py).
Merges build_supercell.py, box_shuffle.py and insert_umd.py.
"""

import sys, getopt, os, re
import numpy as np
from scipy.spatial.transform import Rotation
from scipy.sparse.csgraph import connected_components
import crystallography as cr
import umd_processes_fast as umdpf
from vibrations_umd import KB, AMU_A2FS2                          #Boltzmann constant (eV/K), 1 amu.(A/fs)^2 in eV


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

def count_snapshots(UMDname):#complete snapshots
    return umdpf.count_snapshots(UMDname,umdpf.Crystallization(UMDname)[0].natom)


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


def read_poscar(PoscarFile):#VASP5 POSCAR/CONTCAR, read by umd_processes_fast.read_poscar
    title, cell, symbols, xred, vels = umdpf.read_poscar(PoscarFile)
    elements = sorted(set(symbols),key=symbols.index)
    print('Read',len(symbols),'atoms from',PoscarFile,'with elements',elements,'and counts',[symbols.count(el) for el in elements],'' if vels is not None else '; no velocities')
    hasvel = np.full(len(symbols),vels is not None)
    return Structure(cell,symbols,xred,vels if vels is not None else np.zeros((len(symbols),3)),PoscarFile,hasvel=hasvel)


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


def replace(st,arg,rng):#Replaces atoms of element A by element B : all, n random ones, or a fraction (A:B, A:B:n, A:B:0.25) ; the velocities are scaled by
    #sqrt(m_A/m_B), which keeps the kinetic energy of each atom
    parts = arg.split(':')
    if len(parts) not in (2,3):
        print('--replace takes A:B, A:B:n or A:B:fraction, e.g. Fe:Mg:4')
        sys.exit()
    old, new = parts[0], parts[1]
    candidates = np.array([i for i,el in enumerate(st.symbols) if el == old],dtype=int)
    if len(parts) == 2:
        n = len(candidates)
    elif '.' in parts[2]:
        n = int(round(float(parts[2])*len(candidates)))
    else:
        n = int(parts[2])
    if n > len(candidates):
        print('--replace',arg,': only',len(candidates),'atoms of',old)
        sys.exit()
    chosen = rng.choice(candidates,n,replace=False) if n < len(candidates) else candidates
    symbols = list(st.symbols)
    vels = st.vels.copy()
    vels[chosen] *= np.sqrt(st.mass(old)/st.mass(new))
    for i in chosen:
        symbols[i] = new
    print('replace :',n,'of',len(candidates),old,'by',new,'(velocities x','%.4f' % np.sqrt(st.mass(old)/st.mass(new))+')')
    return st.new(st.cell,symbols,st.xred,vels,st.hasvel)


def read_cutoffs(bondfile,length):#Bond cutoffs : a function of two elements, from a file of lines "A B cutoff" (e.g. maximum_bonds.inp of analyze_gofr.py ;
    #pairs absent : no bond) or one length for all the pairs
    if bondfile:
        cutoffs = {}
        with open(bondfile,'r') as f:
            for line in f:
                entry = line.split()
                if len(entry) >= 3:
                    cutoffs[(entry[0],entry[1])] = cutoffs[(entry[1],entry[0])] = float(entry[2])
        return lambda a,b: cutoffs.get((a,b),0.0)
    if length:
        return lambda a,b: length
    print('--remove needs the bond lengths : -l <length> or -b <file of lines "A B cutoff">')
    sys.exit()


def composition(formula):#'H2O' or 'H_2O_1' -> {'H': 2, 'O': 1}
    comp = {}
    for el,n in re.findall(r'([A-Z][a-z]?)_?(\d*)',formula):
        comp[el] = comp.get(el,0)+(int(n) if n else 1)
    return comp


def remove(st,arg,cutoff,rng):#Removes whole molecules of a formula : all, or n chosen at random (H2O, H2O:3) ; a molecule = a group of atoms connected by bonds
    #(nearest images) with exactly this composition, so O alone means free O atoms
    formula, _, count = arg.partition(':')
    target = composition(formula)
    elements = sorted(set(st.symbols),key=st.symbols.index)
    limit = np.array([[cutoff(a,b) for b in elements] for a in elements])
    typat = np.array([elements.index(el) for el in st.symbols])
    bonds = distances(st.xred,st.cell) < limit[typat][:,typat]
    ncomp, label = connected_components(bonds,directed=False)
    molecules = [np.nonzero(label == c)[0] for c in range(ncomp)]
    found = [m for m in molecules if composition(''.join(st.symbols[i] for i in m)) == target]
    n = len(found) if not count else int(count)
    if n > len(found):
        print('--remove',arg,': only',len(found),'molecules',formula,'found with these bond lengths')
        sys.exit()
    chosen = [found[k] for k in rng.choice(len(found),n,replace=False)] if n < len(found) else found
    drop = np.concatenate(chosen) if chosen else np.zeros(0,dtype=int)
    keep = np.setdiff1d(np.arange(len(st.symbols)),drop)
    print('remove :',n,'of',len(found),'molecules',formula,'(',len(drop),'atoms ) ;',len(molecules),'molecules in all, the largest of',max(len(m) for m in molecules),'atoms')
    return st.new(st.cell,[st.symbols[i] for i in keep],st.xred[keep],st.vels[keep],st.hasvel[keep])


def velocities(st,temperature,redraw,rng):#Maxwell-Boltzmann velocities at temperature for the atoms without velocities (all with redraw) : no net
    #momentum of the drawn atoms, their kinetic temperature set to temperature (3(n-1) degrees of freedom) ; one atom : a plain draw
    draw = np.ones(len(st.symbols),dtype=bool) if redraw else ~st.hasvel
    n = draw.sum()
    if n == 0:
        return st
    masses = st.masses()[draw]
    v = rng.standard_normal((n,3))*np.sqrt(KB*temperature/(masses*AMU_A2FS2))[:,None]                  #sqrt(kB T/m) in A/fs
    if n > 1:
        v -= (masses[:,None]*v).sum(axis=0)/masses.sum()
        ekin = 0.5*(masses[:,None]*v**2).sum()*AMU_A2FS2
        v *= np.sqrt(1.5*(n-1)*KB*temperature/ekin)
    vels = st.vels.copy()
    vels[draw] = v
    print('velocities : Maxwell-Boltzmann at',temperature,'K for',n,'atoms'+(' (all)' if redraw else ' without velocities'))
    return st.new(st.cell,st.symbols,st.xred,vels,np.ones(len(st.symbols),dtype=bool))


#****Extraction of atoms : a trajectory

def selection(arg,elements,typat):#Atoms of a selection : indices and ranges counted from 0, and element symbols, separated by commas (e.g. 12 ; 0-9,25 ; O,H ; Si,12)
    chosen = set()
    for item in arg.split(','):
        if item in elements:
            chosen.update(i for i,t in enumerate(typat) if elements[t] == item)
        elif re.fullmatch(r'\d+(-\d+)?',item):
            lo, _, hi = item.partition('-')
            chosen.update(range(int(lo),int(hi or lo)+1))
        else:
            print('--extract : "'+item+'" is neither an element of the file',elements,'nor an atom index or range')
            sys.exit()
    if max(chosen) >= len(typat):
        print('--extract : atom',max(chosen),'is not in the file (',len(typat),'atoms, counted from 0)')
        sys.exit()
    return sorted(chosen)


def extract(UMDname,arg,first,sampling,outname):#Trajectory of the selected atoms, every sampling-th snapshot from first on, read and written line by line : the atom lines
    #are copied, and so are the lines of each snapshot (energies, T, P, cell : those of the whole system), with the timestep multiplied by sampling
    MyCrystal,TimeStep = umdpf.Crystallization(UMDname)
    atoms = selection(arg,MyCrystal.elements,MyCrystal.typat)
    keep = set(atoms)
    present = sorted(set(MyCrystal.typat[i] for i in atoms))           #elements of the selection, in the order of the file
    newtypat = [present.index(MyCrystal.typat[i]) for i in atoms]
    counts = [newtypat.count(k) for k in range(len(present))]
    header = {'natom': str(len(atoms)), 'ntypat': str(len(present)),
              'types': ' '.join(str(n) for n in counts)+' ', 'elements': ' '.join(MyCrystal.elements[t] for t in present)+' ',
              'masses': ' '.join(str(MyCrystal.masses[t]) for t in present)+' ', 'Zelectrons': ' '.join(str(MyCrystal.zelec[t]) for t in present)+' ',
              'no.of.electrons': str(sum(MyCrystal.zelec[t]*n for t,n in zip(present,counts))), 'typat': ' '.join(str(t) for t in newtypat)+' '}
    written = 0
    with open(UMDname,'r') as fin, open(outname,'w') as fout:
        isnap = -1
        block = []
        inheader = True
        for line in fin:
            word = line.split(' ',1)[0].strip()
            if inheader:
                if word == 'timestep':
                    inheader = False
                else:
                    fout.write(word+' '+header[word]+'\n' if word in header else line)
                    continue
            if word == 'timestep':
                isnap += 1
                use = isnap >= first and (isnap-first) % sampling == 0
                block = []
            if not use:
                continue
            if word == 'timestep':
                fout.write('timestep '+str(TimeStep*sampling)+' fs\n')
            elif word == 'atoms:':
                fout.write(line)
                for iatom in range(MyCrystal.natom):
                    atomline = next(fin)
                    if iatom in keep:
                        fout.write(atomline)
                fout.write('\n')
                written += 1
            elif line.strip():
                fout.write(line)
    print('extract :',len(atoms),'atoms (',', '.join('%d %s' % (n,MyCrystal.elements[t]) for t,n in zip(present,counts)),') of',MyCrystal.natom,', ',written,
          'snapshots (every',sampling,'from',first,') written in',outname)


#****Output

def ordered(st):#Atoms grouped by element (VASP and UMD), in the order of first appearance ; the order inside each element is kept
    elements = sorted(set(st.symbols),key=st.symbols.index)
    order = sorted(range(len(st.symbols)),key=lambda i: elements.index(st.symbols[i]))
    counts = [st.symbols.count(el) for el in elements]
    return elements, counts, np.array(order,dtype=int)


def write_poscar(filename,st):#Direct coordinates, Cartesian velocities in A/fs, by umd_processes_fast.write_poscar
    elements, counts, order = ordered(st)
    umdpf.write_poscar(filename,st.title,st.cell,elements,counts,st.xred[order],st.vels[order])


def write_xyz(filename,st):#Extended XYZ (with the cell), by umd_processes_fast.xyz_frame
    elements, counts, order = ordered(st)
    xcart = st.xred @ st.cell
    with open(filename,'w') as ff:
        ff.write(umdpf.xyz_frame([st.symbols[i] for i in order],xcart[order],st.cell,'title="'+st.title.replace('"',"'")+'"'))


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
    ekin = 0.5*(masses[:,None]*st.vels**2).sum()*AMU_A2FS2
    MyCrystal.kineticenergy = ekin
    MyCrystal.temperature = 2.0*ekin/(3.0*len(masses)*KB)
    xcart = st.xred @ st.cell
    MyCrystal.atoms = [cr.Atom() for _ in range(MyCrystal.natom)]
    for inew,i in enumerate(order):
        MyCrystal.atoms[inew].xred = list(st.xred[i])
        MyCrystal.atoms[inew].xcart = list(xcart[i])
        MyCrystal.atoms[inew].vels = list(st.vels[i])
    base = filename[:-8]
    open(filename,'w').close()                                          #print_header and print_snapshots append
    umdpf.print_header(base,MyCrystal,'build_umd.py : the UMD file is a structure built by build_umd.py, '+st.title)
    umdpf.print_snapshots(base,MyCrystal,1.0,0.0,xcart[order])


def main(argv):
    infile = ''
    first = -1
    sampling = 0
    box = None
    operations = []
    factor = None
    mindist = 1.5
    ntries = 100
    rcut = 2.0
    temperature = 1500.0
    redraw = False
    bondfile = ''
    length = 0.0
    extraction = None
    firstgiven = False
    nconfig = 1
    seed = None
    outputs = ['umd','poscar']
    name = ''
    umdpf.headerumd()
    usage = ('build_umd.py -f <UMD or POSCAR file> -i <snapshot> -s <sampling> -a <empty box> --supercell <na,nb,nc> --shuffle <sub-cells> '
             '--insert <molecules file> --scale <factors> --density <g/cm3> --replace <A:B:n> --remove <formula:n> -d <min distance> -c <factor> '
             '--tries <tries> -r <cutoff> -l <bond length> -b <bond file> -t <temperature> --velocities all -n <configurations> --seed <seed> -o <outputs> --name <name>'
             '\n       build_umd.py -f <UMD file> --extract <atoms> -i <first snapshot> -s <sampling> --name <name>')
    try:
        opts, arg = getopt.getopt(argv,"hf:i:s:a:c:d:t:r:l:b:n:o:",["supercell=","shuffle=","insert=","scale=","density=","replace=","remove=","extract=","tries=","velocities=","seed=","name="])
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('build_umd.py : builds new structures (UMD, VASP POSCAR or xyz files) from a snapshot of a UMD file, a POSCAR file or an empty box,')
            print('by the operations --supercell, --shuffle, --insert, --scale, --density, --replace and --remove, applied in the order of the command')
            print('line ; each can be repeated. The atoms without velocities then get Maxwell-Boltzmann velocities (-t). Usage :')
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
            print('     1.5) are pushed apart along the line joining them, the lighter one more ; no pair is pushed beyond its distance before the')
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
            print('--replace A:B, A:B:n or A:B:x : replaces all the atoms of element A by element B, or n of them, or a fraction x (with a decimal point),')
            print('     chosen at random ; their velocities are scaled by sqrt(m_A/m_B), which keeps their kinetic energy.')
            print('--remove formula or formula:n : removes all the molecules of this formula (e.g. H2O, CO2, O for free O atoms), or n of them chosen at')
            print('     random. A molecule is a group of atoms connected by bonds (nearest images) with exactly this composition ; the bonds are the pairs')
            print('     closer than -l (A, one length for all the pairs) or than the cutoffs of the file -b (lines "A B cutoff", e.g. the')
            print('     maximum_bonds.inp of analyze_gofr.py ; absent pairs : no bond).')
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
            print('Extraction (instead of the operations) :')
            print('--extract <atoms> : writes the trajectory of a selection of atoms in a new UMD file : atom indices and ranges counted from 0, and')
            print('     element symbols, separated by commas (e.g. 12 ; 0-9,25 ; O,H ; Si,12). Every -s-th snapshot (default 1) from -i (default 0).')
            print('     The atom lines are copied as they are ; so are the lines of each snapshot (energies, temperature, pressure, cell), which are')
            print('     those of the whole system, with the timestep multiplied by -s. Read and written line by line : no limit on the size of the')
            print('     file. Output <name>.<atoms>.umd.dat (commas replaced by _), or --name.')
            sys.exit()
        elif opt == '-f':
            infile = arg
        elif opt == '-i':
            first = int(arg)
            firstgiven = True
        elif opt == '-s':
            sampling = int(arg)
        elif opt == '-a':
            box = [float(v) for v in arg.split(',')]
            box = box*3 if len(box) == 1 else box
        elif opt in ('--supercell','--shuffle','--insert','--scale','--density','--replace','--remove'):
            operations.append((opt[2:],arg))
        elif opt == '--extract':
            extraction = arg
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
        elif opt == '-l':
            length = float(arg)
        elif opt == '-b':
            bondfile = arg
        elif opt == '-n':
            nconfig = int(arg)
        elif opt == '--seed':
            seed = int(arg)
        elif opt == '-o':
            outputs = arg.split(',')
        elif opt == '--name':
            name = arg

    if extraction is not None:
        if operations or box is not None:
            print('--extract writes a trajectory : it cannot be combined with the operations or with -a')
            sys.exit()
        if not os.path.isfile(infile):
            print('the UMD file',infile,'does not exist')
            sys.exit()
        nsnap = count_snapshots(infile)
        start = (first if first >= 0 else first+nsnap) if firstgiven else 0
        base = infile[:-8] if infile.endswith('.umd.dat') else infile
        extract(infile,extraction,start,max(sampling,1),(name or base+'.'+extraction.replace(',','_'))+'.umd.dat')
        return

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
                elif op == 'replace':
                    st = replace(st,arg,rng)
                elif op == 'remove':
                    st = remove(st,arg,read_cutoffs(bondfile,length),rng)
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
