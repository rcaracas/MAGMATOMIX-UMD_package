#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### Author : Razvan Caracas
### Exchange of atoms between a condensed phase (melt, droplet, slab) and the gas, through their interface. The condensed phase is defined by a
### smoothed density of its elements (all but the volatiles) at every snapshot, as the instantaneous surface of Willard and Chandler (J. Phys. Chem. B
### 114, 1954, 2010) : Gaussians of width sigma, the surface where the density is half that of the bulk. Every atom is in the condensed phase (L),
### the interface (I : within delta of the surface) or the gas (G), by its distance to the surface ; the paths through the interface (bounces,
### condensations, visits, evaporations) and the time spent in it ; the exchanges, their rates, the fluxes in both directions in time and the
### equilibrium between them. Moved from composition_umd.py (2026-11-16), where the phase came from a grid of sub-cells


import sys, getopt, os.path, re
import time
import heapq
from itertools import groupby
from functools import partial
import numpy as np
from scipy import ndimage
import umd_processes_fast as umdpf
from composition_umd import statistics, VOLATILES
from vibrations_umd import KB, AMU_A2FS2                          #Boltzmann constant (eV/K), 1 amu (A/fs)^2 in eV

EV_A3_GPA = 160.21766208                #1 eV/A^3 in GPa
NA = 6.02214076e23                      #Avogadro constant, 1/mol
#Van der Waals a of gases, L^2 bar/mol^2 (CRC Handbook of Chemistry and Physics), by composition ; 1 L^2 bar/mol^2 = 0.1 Pa m^6/mol^2. Used in the
#modified Van der Waals law of the gas, P = N k T/V - a N^2/(NA^2 V^2) (no excluded volume ; Davis and Caracas 2024, Commun. Earth Environ.
#5, 344, eq. 1), for a mixture with a N^2 = (sum_i N_i a_i^1/2)^2 ; other species : a = 0 (ideal), or -A
VDW_A = {'H2':0.2476,'He':0.0346,'Ne':0.2135,'Ar':1.355,'Kr':2.349,'Xe':4.250,'N2':1.370,'O2':1.382,'CO':1.505,'CO2':3.640,'H2O':5.536,'CH4':2.283,
         'NH3':4.225,'H2S':4.490,'SO2':6.803,'HCl':3.716,'Cl2':6.579,'F2':1.171,'NO':1.358,'N2O':3.832,'NO2':5.354,'HF':9.565,'HCN':11.29,'C2H6':5.562,
         'SiH4':4.377,'Hg':8.200}


GAS, INTERFACE, CONDENSED = 0, 1, 2
PATHS = ('G-I-G','G-I-L','L-I-L','L-I-G')           #through the interface : bounce from the gas, condensation, visit from the melt, evaporation


#----------------------------------------------------------------- liquid and vapor species (formerly liquid_and_vapor_speciation.py)
#The species of a popul file (any -r, typically r1) separated into liquid and vapor, from the clusters of an r0 step file of the same bonding
#file : an atom is in the vapor at a step when its r0 cluster has fewer than -l atoms ; a species is in the liquid when one of its atoms is.
#exchange_umd.py -P <popul file> -S <r0 step file> -l <length> ; split() is used by lifetimes.py -l

def read_step(stepfile,length):#Liquid atoms of each step of an r0 step file, and the steps
    liquid = {}
    with open(stepfile,'r') as f:
        lines = f.readlines()
    step = None
    for line in lines[1:]:
        entry = line.rstrip('\n').split('\t')
        if entry[0] == 'step':
            step = int(entry[1])
            liquid[step] = set()
        elif step is not None and len(entry) == 3 and entry[1].isdigit():
            if int(entry[1]) >= length:
                liquid[step].update(int(k) for k in entry[2].strip('[]').split(','))
    return liquid


def read_popul(populfile):#Header values and species of a popul file
    with open(populfile,'r') as f:
        lines = f.readlines()
    found = re.search(r'-s=(\d+)',lines[0])
    nsteps = int(found.group(1)) if found else 1                             #sampling of the speciation run, in steps of the bonding file
    found = re.search(r'timestep\s+([0-9.eE+-]+)\s*fs',lines[0])
    timestep = float(found.group(1)) if found else None                     #time between two steps of the bonding file
    species = []
    for line in lines[2:]:
        entry = line.rstrip('\n').split('\t')
        if len(entry) < 5: continue
        species.append([entry[0],int(entry[1]),int(entry[2]),float(entry[3]),entry[4]])
    return nsteps, timestep, species


def species_runs(species,phasesteps,nsteps,dt):#Species of one phase : each run of consecutive snapshots of a species in this phase ; phasesteps : for each species, its steps in this phase
    runs = []
    for (formula,begin,end,lifetime,composition),steps in zip(species,phasesteps):
        for _,run in groupby(enumerate(steps),lambda x: x[1]//nsteps - x[0]):#consecutive snapshots
            run = [s for _,s in run]
            runs.append([formula,run[0],run[-1],len(run)*nsteps*dt,composition])
    return runs


def split(populfile,stepfile,length):#Vapor and liquid species of a popul file (formula, begin, end, lifetime, composition), from an r0 step file
    liquid = read_step(stepfile,length)
    nsteps, timestep, species = read_popul(populfile)
    if timestep is None:#older popul files : the time step from the step file
        steps = sorted(liquid)
        with open(stepfile,'r') as f:
            times = [float(line.split('\t')[1]) for line in f if line.startswith('time\t')]
        timestep = (times[1]-times[0])/(steps[1]-steps[0])
    vaporsteps = []
    liquidsteps = []
    for formula,begin,end,lifetime,composition in species:
        atoms = [int(j) for j in composition.strip('[]').split(',')]
        vapor = []
        liq = []
        for step in range(begin,end+1,nsteps):
            if step not in liquid:
                print('ERROR : step',step,'of',populfile,'is not in',stepfile,': use the same sampling (-s) in speciation_and_angles.py')
                sys.exit(1)
            if any(at in liquid[step] for at in atoms):
                liq.append(step)
            else:
                vapor.append(step)
        vaporsteps.append(vapor)
        liquidsteps.append(liq)
    return species_runs(species,vaporsteps,nsteps,timestep), species_runs(species,liquidsteps,nsteps,timestep)


def write_phase(filename,runs):#popul and stat files of one phase
    with open(filename+'.popul.dat','w') as nf:
        nf.write('Formula\tBegin (step)\tEnd (step)\tLifetime (fs)\tComposition\n')
        stats = {}
        total = 0.0
        for formula,begin,end,life,composition in runs:
            nf.write(formula+'\t'+str(begin)+'\t'+str(end)+'\t'+str(life)+'\t'+composition+'\n')
            total += life
            if formula in stats:
                stats[formula][0] += life
            else:
                stats[formula] = [life,len(composition.split(','))]
    with open(filename+'.stat.dat','w') as nf:
        nf.write('Cluster\tTime (fs)\tPercent\tNumber of atoms\n')
        for formula in stats:
            nf.write(formula+'\t'+str(stats[formula][0])+'\t'+str(stats[formula][0]/total)+'\t'+str(stats[formula][1])+'\n')



#----------------------------------------------------------------- molecules from a bonding file

def read_molecules(bondfile,natom):#Molecules of every step of a bonding file of bonding_umd.py (the connected bonds) : the UMD snapshot of its first step
    #(first_snapshot, 0 when absent), its timestep (to set against that of the UMD file), and per step the molecule of each atom (labels 0..n-1) and
    #the number of atoms of each molecule
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    first, timestep = 0, None
    labels, sizes = [], []
    rows, cols, started = [], [], False
    def close():
        n, lab = connected_components(coo_matrix((np.ones(len(rows)),(rows,cols)),shape=(natom,natom)),directed=False)
        labels.append(lab.astype(np.int32))
        sizes.append(np.bincount(lab,minlength=n).astype(np.int32))
    with open(bondfile) as f:
        for line in f:
            w = line.split()
            if not w:
                continue
            if w[0] == 'step':
                if started:
                    close()
                rows, cols, started = [], [], True
            elif not started:
                if w[0] == 'timestep':
                    timestep = float(w[1])
                elif w[0] == 'first_snapshot':
                    first = int(w[1])
            elif len(w) > 1 and w[0].isdigit():                          #an atom and its bound atoms (a time line is skipped)
                a = int(w[0])
                for b in w[1:]:
                    rows.append(a)
                    cols.append(int(b))
    if started:
        close()
    return first, timestep, labels, sizes


def composition(formula):#Composition of a formula, UMD style (O_1H_2) or chemical (H2O) : sorted tuple of (element, number)
    parts = re.findall(r'([A-Z][a-z]?)_(\d+)',formula) if '_' in formula else [(el,n or '1') for el, n in re.findall(r'([A-Z][a-z]?)(\d*)',formula)]
    comp = {}
    for el, n in parts:
        comp[el] = comp.get(el,0)+int(n)
    return tuple(sorted(comp.items()))


def vdw_table(extra):#Van der Waals a (L^2 bar/mol^2) by composition : the table VDW_A, and -A (formula:a, separated by commas) over it
    table = {composition(f):a for f, a in VDW_A.items()}
    for item in extra:
        f, a = item.rsplit(':',1)
        table[composition(f.strip())] = float(a)
    return table


def formula_of(atoms,typat,elements):#Name of a molecule from its atoms : El_n in the order of the elements of the UMD file, as speciation_and_angles.py
    counts = np.bincount(typat[atoms],minlength=len(elements))
    return ''.join(el+'_'+str(c) for el, c in zip(elements,counts) if c)


def grid_shape(cell,spacing):#Points of the density grid along a, b, c : at most spacing apart
    return np.maximum(4,np.ceil(np.linalg.norm(cell,axis=1)/spacing).astype(np.int64))


def smoothed_density(red,cell,shape,sigma):#Number density of the atoms of reduced coordinates red (atoms/A^3) on a periodic grid of shape, smoothed by
    #Gaussians of width sigma (A) in Cartesian space (any cell) ; deposited by cloud in cell ; returns the density and its Cartesian gradient
    size = int(np.prod(shape))
    s = (red % 1.0)*shape
    i0 = np.floor(s).astype(np.int64)
    f = s-i0
    rho = np.zeros(size)
    for dx in (0,1):
        for dy in (0,1):
            for dz in (0,1):
                w = (f[:,0] if dx else 1-f[:,0])*(f[:,1] if dy else 1-f[:,1])*(f[:,2] if dz else 1-f[:,2])
                idx = (((i0[:,0]+dx) % shape[0])*shape[1]+(i0[:,1]+dy) % shape[1])*shape[2]+(i0[:,2]+dz) % shape[2]
                rho += np.bincount(idx,weights=w,minlength=size)
    rho = rho.reshape(shape)/(abs(np.linalg.det(cell))/size)
    rec = 2.0*np.pi*np.linalg.inv(cell).T                             #k = n0 rec0 + n1 rec1 + n2 rec2 (rows)
    n0 = np.fft.fftfreq(shape[0])*shape[0]
    n1 = np.fft.fftfreq(shape[1])*shape[1]
    n2 = np.fft.rfftfreq(shape[2])*shape[2]
    kx = [n0[:,None,None]*rec[0][c]+n1[None,:,None]*rec[1][c]+n2[None,None,:]*rec[2][c] for c in range(3)]
    F = np.fft.rfftn(rho)*np.exp(-0.5*sigma**2*(kx[0]**2+kx[1]**2+kx[2]**2))
    dens = np.fft.irfftn(F,s=tuple(shape))
    grad = [np.fft.irfftn(1j*kx[c]*F,s=tuple(shape)) for c in range(3)]
    return dens, grad


def self_density(cell,shape,sigma,trials=8,seed=7):#Density that an atom gives at its own position (Gaussian, cloud in cell, linear interpolation),
    #to remove from the density at a melt atom, which then counts only the others (about 1/((2 pi)^1.5 sigma^3))
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(trials):
        red = rng.random((1,3))
        dens, grad = smoothed_density(red,cell,shape,sigma)
        vals.append(ndimage.map_coordinates(dens,((red % 1.0)*shape).T,order=1,mode='grid-wrap')[0])
    return float(np.mean(vals))


def ball(cell,shape,radius):#Footprint of the grid points within radius (A, Cartesian) of a point, for any cell : offsets along a, b, c
    steps = cell/np.array(shape)[:,None]                                #the grid vectors (rows)
    widths = abs(np.linalg.det(cell))/np.array([np.linalg.norm(np.cross(cell[1],cell[2])),np.linalg.norm(np.cross(cell[2],cell[0])),
                                                np.linalg.norm(np.cross(cell[0],cell[1]))])
    n = np.ceil(radius/(widths/np.array(shape))).astype(int)
    o = np.stack(np.meshgrid(*[np.arange(-k,k+1) for k in n],indexing='ij'),axis=-1)
    return np.linalg.norm(o @ steps,axis=-1) <= radius


def zone_fields(dens,foot):#Least and largest density within the footprint of every grid point (periodic) : the surface rho = level passes within the
    #radius of a point when the least is below the level and the largest above
    return ndimage.minimum_filter(dens,footprint=foot,mode='wrap'), ndimage.maximum_filter(dens,footprint=foot,mode='wrap')


def at_atoms(dens,grad,red,shape):#Density and the norm of its gradient at the atoms
    coords = ((red % 1.0)*shape).T
    rho = ndimage.map_coordinates(dens,coords,order=1,mode='grid-wrap')
    g = np.sqrt(sum(ndimage.map_coordinates(gc,coords,order=1,mode='grid-wrap')**2 for gc in grad))
    return rho, g


def state_batch(task,umdfile,offsets,natom,sampling,cells,melt,shape,sigma,rhob,delta,selfval,masses,typat,nel,foot=None):#For snapshots a to b of the selection :
    #the state of every atom (gas 0, interface 1, condensed 2) : interface when the surface (smoothed density of the melt elements = rhob/2) passes
    #within delta of it (the density within delta both below and above the level), condensed when all the density within delta is above, gas when
    #all is below ; the area of the surface (coarea : the integral of |grad rho| delta(rho - rhob/2), the delta a Gaussian in density as wide as
    #the change of density over one grid spacing across the surface ; a sharp shell thinner than the grid missed flat surfaces), the volume of the gas (density below rhob/2), and per element and state the atoms and the sum of m v^2 (amu A^2/fs^2 ;
    #temperature of the gas). A linear estimate of the distance, (rho - rhob/2)/|grad rho|, flagged 6 % of the atoms deep in a slab as interface
    a, b = task
    X = umdpf.read_atoms_at(umdfile,offsets[a],natom,b-a,sampling,0,12)  #reduced, Cartesian, unwrapped, velocity
    states = np.empty((b-a,natom),dtype=np.int8)
    outside = np.empty((b-a,natom),dtype=bool)                         #below the half-density surface (the density of the others at the atom)
    mv2s = np.empty((b-a,natom),dtype=np.float32)                       #m v^2 of each atom (amu A^2/fs^2 ; temperature of the gas molecules)
    area, gasvol = np.empty(b-a), np.empty(b-a)
    zones = np.empty((b-a,2))                                           #volumes of the gas zone and of the condensed zone (the atoms of G and L live there)
    kin = np.zeros((b-a,nel,3,2))
    level = 0.5*rhob
    spacing = abs(np.linalg.det(cells[a]))**(1.0/3.0)/np.prod(shape)**(1.0/3.0)
    width = spacing*rhob/(sigma*np.sqrt(2.0*np.pi))                     #the density change over one grid spacing across a flat surface
    for t in range(b-a):
        cell = cells[a+t]
        red = X[t,:,0:3]
        dens, grad = smoothed_density(red[melt],cell,shape,sigma)
        lo, hi = zone_fields(dens,foot if foot is not None else ball(cell,shape,delta))
        coords = ((red % 1.0)*shape).T
        lo_at = ndimage.map_coordinates(lo,coords,order=1,mode='grid-wrap')
        hi_at = ndimage.map_coordinates(hi,coords,order=1,mode='grid-wrap')
        states[t] = np.where(lo_at >= level,CONDENSED,np.where(hi_at < level,GAS,INTERFACE))
        outside[t] = ndimage.map_coordinates(dens,coords,order=1,mode='grid-wrap')-selfval*melt < level
        voxel = abs(np.linalg.det(cell))/dens.size
        zones[t] = np.count_nonzero(hi < level)*voxel, np.count_nonzero(lo >= level)*voxel
        gnorm = np.sqrt(grad[0]**2+grad[1]**2+grad[2]**2)
        area[t] = (gnorm*np.exp(-0.5*((dens-level)/width)**2)).sum()*voxel/(width*np.sqrt(2.0*np.pi))   #coarea, a smooth delta
        gasvol[t] = np.count_nonzero(dens < level)*voxel
        mv2 = masses[typat]*(X[t,:,9:12]**2).sum(axis=1)
        mv2s[t] = mv2
        for st in (GAS,INTERFACE,CONDENSED):
            sel = states[t] == st
            kin[t,:,st,0] = np.bincount(typat[sel],minlength=nel)
            kin[t,:,st,1] = np.bincount(typat[sel],weights=mv2[sel],minlength=nel)
    return states, outside, mv2s, area, gasvol, zones, kin


def bulk_density(umdfile,offsets,natom,cells,melt,shape,sigma,delta,selfval,nsample=40):#Density of the bulk of the condensed phase, on a sample of
    #the snapshots : the mean smoothed density at the grid points deep inside it, all the density within 2 delta of them above half the first
    #estimate (the 90th percentile of the density at the melt atoms, the atom left out). At grid points, not at atoms : at an atom the density of
    #the others misses the correlation hole of a condensed phase (11 % low on a slab of known density). The choice of the points barely depends on
    #the density averaged (the least density within 2 delta, deep in the bulk, is far above half) ; choosing with a linear distance to the
    #surface, (rho - rho_half)/|grad rho|, favoured the points where the density fluctuates up (8 % high), and choosing with half of the atoms
    #and averaging the other half is biased low in a condensed phase (the two halves exclude each other locally)
    picks = np.unique(np.linspace(0,len(offsets)-1,min(nsample,len(offsets))).astype(int))
    firsts, fields = [], []
    for k in picks:
        red = umdpf.read_atoms_at(umdfile,offsets[k],natom,1,1,0,3)[0][melt]
        dens, grad = smoothed_density(red,cells[k],shape,sigma)
        firsts.append(at_atoms(dens,grad,red,shape)[0]-selfval)
        fields.append((k,dens))
    first = np.percentile(np.concatenate(firsts),90)
    vals = np.concatenate([dens[ndimage.minimum_filter(dens,footprint=ball(cells[k],shape,2.0*delta),mode='wrap') >= 0.5*first] for k, dens in fields])
    return (vals.mean(), len(vals)) if len(vals) >= 20 else (first, len(vals))


def state_runs(seq,shortest):#Stays of an atom in its states (series of 0, 1, 2 in time) : [[state, length], ...] ; an inner stay shorter than shortest
    #snapshots between two stays in the same state is merged with them (the shortest first, then the earliest), so that crossings back and forth are
    #not counted ; a short stay between two different states is kept (it is a passage). A heap over a linked list : N log N
    change = np.nonzero(np.diff(seq.astype(np.int8)))[0]+1
    starts = np.concatenate([[0],change])
    length = np.diff(np.concatenate([starts,[len(seq)]])).tolist()
    state = [int(seq[st]) for st in starts]
    n = len(length)
    if shortest <= 1 or n <= 2:
        return [[s,l] for s, l in zip(state,length)]
    prev = list(range(-1,n-1))
    nxt = list(range(1,n+1))
    nxt[-1] = -1
    alive = [True]*n
    def eligible(k):
        return prev[k] >= 0 and nxt[k] >= 0 and length[k] < shortest and state[prev[k]] == state[nxt[k]]
    heap = [(length[k],k) for k in range(n) if eligible(k)]
    heapq.heapify(heap)
    while heap:
        l, k = heapq.heappop(heap)
        if not alive[k] or l != length[k] or not eligible(k):
            continue
        p, q = prev[k], nxt[k]
        length[p] += length[k]+length[q]
        alive[k] = alive[q] = False
        nxt[p] = nxt[q]
        if nxt[q] >= 0:
            prev[nxt[q]] = p
        for r in (prev[p],p,nxt[p]):
            if r >= 0 and eligible(r):
                heapq.heappush(heap,(length[r],r))
    return [[state[k],length[k]] for k in range(n) if alive[k]]


def follow_atom(runs):#From the stays of an atom : its committed phase (condensed or gas ; a stay in the interface belongs to the phase before it,
    #or after it at the start), the snapshots of its evaporations and condensations, its committed stays, and its paths through the interface
    committed = []
    last = next((s for s, l in runs if s != INTERFACE),None)
    for s, l in runs:
        if s != INTERFACE:
            last = s
        committed.append(last)
    evap, cond, paths, stays = [], [], [], []
    pos = 0
    for i, (s, l) in enumerate(runs):
        if i > 0 and committed[i] is not None and committed[i-1] is not None and committed[i] != committed[i-1]:
            (evap if committed[i] == GAS else cond).append(pos)
        if s == INTERFACE and 0 < i < len(runs)-1:
            paths.append(('LG'[runs[i-1][0] == GAS]+'-I-'+'LG'[runs[i+1][0] == GAS],l))
        elif s != INTERFACE and i > 0 and runs[i-1][0] != INTERFACE:
            paths.append(('direct L-G' if s == GAS else 'direct G-L',0))
        pos += l
    #committed stays : consecutive runs of the same committed phase joined
    cur, length, first = None, 0, True
    for (s, l), c in zip(runs,committed):
        if c == cur:
            length += l
        else:
            if cur is not None:
                stays.append((cur,length,first))
                first = False
            cur, length = c, l
    if cur is not None:
        stays.append((cur,length,first))
    stays = [(c,l,f,i == len(stays)-1) for i, (c,l,f) in enumerate(stays)]   #phase, length, cut by the start, cut by the end
    return committed, evap, cond, paths, stays


def main(argv):
    umdpf.headerumd()
    start = time.time()
    umdfile = ''
    Nsteps = 1
    InitialStep = 0
    nCores = None
    volatiles = None
    transient = 50.0
    window = None
    sigma = 2.4
    spacing = 1.0
    rhob = None
    delta = None
    Tuser = None
    vdwextra = []
    bondfile = None
    mode = 'field'
    length = 12
    populfile = None
    stepfile = None
    usage = 'exchange_umd.py -f <umdfile> -g <volatile elements> -q <sigma> -p <grid spacing> -r <bulk density> -x <interface half-width> -t <transient> -w <window> -T <temperature> -A <vdW a> -b <bonding file> -m <field|clusters> -l <length> -i <InitialStep> -s <Sampling> -k <nCores>\n'+\
            '       exchange_umd.py -P <popul file> -S <r0 step file> -l <length>   (liquid and vapor species)'
    try:
        opts, arg = getopt.getopt(argv,"hf:g:q:p:r:x:t:w:i:s:k:T:b:m:l:P:S:A:",["fumdfile","gvolatiles","qsigma","pspacing","rbulk","xdelta","ttransient","wwindow","iInitialStep","sSampling","knCores","Temperature",
                                                         "bbonding","mmode","llength","Ppopul","Sstep","Avdw"])
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('Program to analyze the exchange of atoms between a condensed phase (melt, droplet, slab) and the gas, through their interface. Usage :')
            print(usage)
            print('-f : UMD file. -f all : every *.umd.dat file of the current folder, one after the other.')
            print('-g : volatile elements, separated by commas, or none. Default : those of '+','.join(VOLATILES)+' present in the run. The condensed')
            print('     phase is defined by the density of the other elements (the melt) : a gas bubble full of volatiles inside the melt is gas.')
            print('-q : width of the Gaussians of the smoothed density, A. Default 2.4 (Willard and Chandler 2010).')
            print('-p : spacing of the density grid, A. Default 1.0.')
            print('-r : density of the bulk of the condensed phase, atoms of the melt elements per A^3. Default : the mean smoothed density at the grid')
            print('     points deep inside it, measured on 40 snapshots. The surface is where the smoothed density is half of it.')
            print('-x : half-width of the interface, A : an atom is in the interface when the surface passes within x of it (the smoothed density')
            print('     within x both below and above half the bulk), in the condensed phase when all of it is above, in the gas when all of it is below.')
            print('     Default 1.28 sigma (3.1 A) : for a flat surface, where its density is 0.1 and 0.9 of the bulk. A distance, not thresholds on the')
            print('     density at the atom, which fluctuates by ~7 % in the bulk of a melt (with thresholds 0.1 and 0.9, 40 % of the atoms of a bulk melt')
            print('     would be in the interface).')
            print('-t : a stay shorter than t fs between two stays in the same state (crossings back and forth) is merged with them. Default 50 fs ;')
            print('     0 : every crossing.')
            print('-w : length of the time windows of the fluxes, fs. Default : a tenth of the run.')
            print('-T : temperature of the gas, K. Default : from the velocities of the atoms in the gas (kinetic) ; without velocities, the mean')
            print('     Temperature of the run.')
            print('-A : Van der Waals a of gas species, L^2 bar/mol^2, formula:a separated by commas (UMD formulas as in molecules.dat, O_1H_1, or')
            print('     chemical, OH), added to or replacing those of the CRC Handbook (H2, He, Ne, Ar, Kr, Xe, N2, O2, CO, CO2, H2O, CH4, NH3, H2S, SO2,')
            print('     HCl, Cl2, F2, NO, N2O, NO2, HF, HCN, C2H6, SiH4, Hg). Species without a count as ideal (listed in the summary).')
            print('-b : bonding file of bonding_umd.py (default <base>.bonding.dat when it exists) : the molecules (connected bonds) of every step ;')
            print('     the gas by molecule (composition, partial pressures) and the evaporations and condensations by molecule (molecules.dat). An atom')
            print('     in the gas zone whose molecule has atoms that are not (bonded to the condensed phase, even by one atom) is adsorbed : it counts')
            print('     in the interface for everything (time in each state, paths, evaporations, gas, temperature) : an atom evaporates when its')
            print('     molecule detaches.')
            print('-m : definition of the condensed phase : field (default : the smoothed density, with the interface) or clusters (with -b : an atom is')
            print('     condensed when its molecule has at least -l atoms, as the r0 clusters of speciation_and_angles.py ; no interface).')
            print('-l : smallest molecule of the condensed phase with -m clusters, and of the liquid with -P. Default 12.')
            print('-P, -S : liquid and vapor species (formerly liquid_and_vapor_speciation.py ; no UMD file) : the species of the popul file -P (any -r,')
            print('     typically r1) separated with the r0 step file -S of the same bonding file (default : <prefix>.r0.step.dat) : an atom is in the')
            print('     vapor at a step when its r0 cluster has fewer than -l atoms ; a species is in the liquid when one of its atoms is. Its steps must be')
            print('     in the step file (same sampling -s of speciation_and_angles.py). Outputs <popul base>.vapor.popul.dat, .vapor.stat.dat,')
            print('     .liquid.popul.dat, .liquid.stat.dat (popul and stat formats).')
            print('-i : InitialStep. number of initial snapshots discarded. Default is 0.')
            print('-s : Sampling. one snapshot out of s is used. Default is 1.')
            print('-k : number of cores. Default is all.')
            print('Any cell, orthogonal or not, NPT : the density of every snapshot in its own cell.')
            print('Outputs, in the folder <base>.exchange (with -m clusters <base>.exchange.clusters) :')
            print('  summary.dat : the interface (area, volume of the gas), the fraction of the time of each element in the condensed phase, the')
            print('     interface and the gas ; the paths through the interface (G-I-G : bounce of a gas atom, G-I-L : condensation, L-I-L : visit of a')
            print('     melt atom, L-I-G : evaporation) with the time spent in the interface, the sticking coefficient G-I-L/(G-I-L + G-I-G) and the')
            print('     escape probability L-I-G/(L-I-G + L-I-L) ; the evaporations and condensations (an atom changes phase when it goes from the')
            print('     condensed phase to the gas or back, through the interface or not), the rate constants, the stays, the fluxes in both directions')
            print('     over the run (per ps and per nm^2) and the time from which they are equal (equilibrium)')
            print('  states.dat : per snapshot, the area of the interface, the volume of the gas, the atoms of each element in each state')
            print('  paths.dat : distribution of the times spent in the interface, per path and element')
            print('  fluxes.dat : per time window and element, the fluxes of evaporation and condensation (events per ps, errors), the net flux, the')
            print('     atoms in the gas and the area of the interface')
            print('  gas.dat : per snapshot, the volumes of the gas and condensed zones, the atoms of each element in the gas, the partial pressures')
            print('     The summary also gives the gas : its temperature (per element too), composition, partial pressures of the atoms (ideal gas,')
            print('     n k T in the gas zone ; atoms counted as independent particles : for molecules an upper bound), the partition coefficients')
            print('     between the condensed phase and the gas (densities, and mole fractions), and the Hertz-Knudsen flux with which the gas hits the')
            print('     interface, n (k T / 2 pi m)^1/2, against the arrivals of gas atoms counted at the interface and the condensation coefficient')
            print('  pressure.dat : the gas by molecules (with a bonding file ; else every atom its own particle), two definitions : zone (all the atoms')
            print('     of a molecule in the gas zone ; its volume) and surface (all below the half-density surface ; the volume below it) ; atoms outside')
            print('     whose molecule is not all outside (bonded to the condensed phase) are adsorbed, not gas. Per sample : gas volume, molecules,')
            print('     adsorbed atoms, pressure ideal N k T/V and modified Van der Waals N k T/V - a N^2/(NA^2 V^2) (Davis and Caracas 2024, Commun. Earth')
            print('     Environ. 5, 344 ; mixture : a N^2 = (sum N_i a_i^1/2)^2), T the kinetic temperature of the atoms of the gas molecules ; means')
            print('     with errors over the run in the summary')
            print('  molecules.dat (with a bonding file) : per molecule, its number in the gas, partial pressure, evaporations and condensations, fluxes,')
            print('     its number in the gas below the half-density surface and its Van der Waals a')
            print('  stays.dat : distribution of the durations of the stays in the gas and in the condensed phase')
            print('  atoms.dat : per atom, the fractions of the time in the gas and the interface, its evaporations and condensations')
            sys.exit()
        elif opt in ("-f","--fumdfile"):
            umdfile = str(arg)
        elif opt in ("-g","--gvolatiles"):
            volatiles = [] if arg.strip().lower() == 'none' else [el.strip() for el in arg.split(',') if el.strip()]
        elif opt in ("-q","--qsigma"):
            sigma = float(arg)
        elif opt in ("-p","--pspacing"):
            spacing = float(arg)
        elif opt in ("-r","--rbulk"):
            rhob = float(arg)
        elif opt in ("-x","--xdelta"):
            delta = float(arg)
        elif opt in ("-t","--ttransient"):
            transient = float(arg)
        elif opt in ("-w","--wwindow"):
            window = float(arg)
        elif opt in ("-i","--iInitialStep"):
            InitialStep = int(arg)
        elif opt in ("-s","--sSampling"):
            Nsteps = int(arg)
        elif opt in ("-k","--knCores"):
            nCores = int(arg)
        elif opt in ("-T","--Temperature"):
            Tuser = float(arg)
        elif opt in ("-b","--bbonding"):
            bondfile = arg
        elif opt in ("-m","--mmode"):
            mode = arg.strip().lower()
            if mode not in ('field','clusters'):
                print('-m : field or clusters')
                sys.exit(2)
        elif opt in ("-l","--llength"):
            length = int(arg)
        elif opt in ("-P","--Ppopul"):
            populfile = arg
        elif opt in ("-A","--Avdw"):
            vdwextra = [x for x in arg.split(',') if x.strip()]
        elif opt in ("-S","--Sstep"):
            stepfile = arg
    if populfile is not None:#liquid and vapor species of a popul file (formerly liquid_and_vapor_speciation.py) : no UMD file needed
        if stepfile is None:
            stepfile = re.sub(r'\.r\d+\.popul\.dat$','.r0.step.dat',populfile)
        for name in (populfile,stepfile):
            if not os.path.isfile(name):
                print('the file',name,'does not exist : -P <popul file> -S <r0 step file of speciation_and_angles.py -r 0>')
                sys.exit(2)
        base = populfile[:-10] if populfile.endswith('.popul.dat') else populfile
        umdpf.remove_former_outputs([base+phase+ext for phase in ('.vapor','.liquid') for ext in ('.popul.dat','.stat.dat')])#of a former run
        vapor, liquid = split(populfile,stepfile,length)
        write_phase(base+'.vapor',vapor)
        write_phase(base+'.liquid',liquid)
        print('Vapor and liquid species written in',base+'.vapor.popul.dat,',base+'.vapor.stat.dat,',base+'.liquid.popul.dat,',base+'.liquid.stat.dat')
        return
    if delta is None:
        delta = 1.2816*sigma                                            #0.1 and 0.9 of the bulk for a flat surface : erfc(delta/(sigma sqrt 2))/2 = 0.1
    umdfile = umdpf.umd_input(umdfile)                     #<run>.umd.dat.gz : decompressed into a temporary file, the outputs named from <run>.umd.dat
    if not os.path.isfile(umdpf.real_path(umdfile)):
        umdpf.missing_umd(umdfile)
        sys.exit()
    MyCrystal, TimeStep = umdpf.Crystallization(umdfile)
    elements = list(MyCrystal.elements)
    nel = len(elements)
    natom = MyCrystal.natom
    typat = np.array(MyCrystal.typat,dtype=np.int64)
    N = np.bincount(typat,minlength=nel)
    masses = np.array(MyCrystal.masses,dtype=float)
    if volatiles is None:
        volatiles = [el for el in elements if el in VOLATILES]
    unknown = [el for el in volatiles if el not in elements]
    if unknown:
        print('the volatile elements',','.join(unknown),'are not in the UMD file, whose elements are',','.join(elements))
        sys.exit(2)
    isvol = np.array([el in volatiles for el in elements])
    if isvol.all():
        print('all the elements are volatile : the condensed phase is defined by the density of all the atoms')
        isvol[:] = False
    melt = ~isvol[typat]

    offsets, cells = umdpf.index_snapshots(umdfile,natom)
    offsets, cells = offsets[InitialStep::Nsteps], cells[InitialStep::Nsteps]
    nsnap = len(offsets)
    if nsnap < 8:
        print('Only',nsnap,'snapshots after the snapshot',InitialStep,'with sampling',Nsteps,': at least 8 are needed')
        sys.exit()
    if np.isnan(cells).any():
        print('ERROR : snapshots of',umdfile,'without rprimd or acell lines')
        sys.exit()
    if 6*nsnap*natom > 4*2**30:
        print('the state, side and kinetic energy of every atom at every snapshot need','%.1f' % (6*nsnap*natom/2**30),'GB (more than 4 GB) : use -s or -i for fewer snapshots')
        sys.exit()
    base = umdfile[:-8] if umdfile.endswith('.umd.dat') else os.path.splitext(umdfile)[0]
    outdir = base+'.exchange'+('.clusters' if mode == 'clusters' else '')     #the two definitions of the condensed phase side by side
    if os.path.exists(outdir) and not os.path.isdir(outdir):
        print(outdir,'is a file, not a folder : remove or rename it')
        sys.exit(2)
    os.makedirs(outdir,exist_ok=True)
    names = {kind:os.path.join(outdir,kind+'.dat') for kind in ('summary','states','gas','pressure','paths','stays','atoms','fluxes','molecules')}
    umdpf.remove_former_outputs(list(names.values()))
    dt = (TimeStep if TimeStep > 0 else 1.0)*Nsteps
    if TimeStep <= 0:
        print('WARNING : no timestep in the UMD file : times are given in snapshots')
    timesall = (InitialStep+np.arange(nsnap)*Nsteps)*(TimeStep if TimeStep > 0 else 1.0)
    shape = grid_shape(cells[0],spacing)
    selfval = self_density(cells[0],shape,sigma)
    if rhob is None:
        rhob, ndeep = bulk_density(umdfile,offsets,natom,cells,melt,shape,sigma,delta,selfval)
        print('bulk density of the condensed phase :','%.5f' % rhob,'atoms/A^3 (mean at',ndeep,'grid points deep inside, on','%d' % min(40,nsnap),'snapshots)')
    print(nsnap,'snapshots ; density grid','%dx%dx%d' % tuple(shape),'; sigma','%g' % sigma,'A ; interface half-width','%.2f' % delta,'A ; volatile elements :',
          ','.join(volatiles) if isvol.any() else 'none')

    #States of the atoms, interface, gas, kinetic energies : in batches by worker processes (as gofr_umd.py)
    nworkers = max(1,min(nCores if nCores is not None else (os.cpu_count() or 1),nsnap))
    snapbytes = (os.path.getsize(umdpf.real_path(umdfile))-offsets[0])/nsnap
    tasks = umdpf.batch_ranges(nsnap,snapbytes,nworkers)
    constant = np.allclose(cells,cells[0])
    worker = partial(state_batch,umdfile=umdfile,offsets=offsets,natom=natom,sampling=Nsteps,cells=cells,melt=melt,shape=shape,sigma=sigma,rhob=rhob,
                     delta=delta,selfval=selfval,masses=masses,typat=typat,nel=nel,foot=ball(cells[0],shape,delta) if constant else None)
    states = np.empty((nsnap,natom),dtype=np.int8)
    outside = np.empty((nsnap,natom),dtype=bool)
    mv2s = np.empty((nsnap,natom),dtype=np.float32)
    area, gasvol = np.empty(nsnap), np.empty(nsnap)
    zones = np.empty((nsnap,2))
    kin = np.empty((nsnap,nel,3,2))
    for (a,b), result in zip(tasks,umdpf.ordered_map(worker,tasks,nworkers)):
        states[a:b], outside[a:b], mv2s[a:b], area[a:b], gasvol[a:b], zones[a:b], kin[a:b] = result
    volumes = np.abs(np.linalg.det(cells))
    #Molecules of the bonding file, matched to the snapshots analyzed : bonding step k is the UMD snapshot first + k sb
    if bondfile is None and os.path.isfile(base+'.bonding.dat'):
        bondfile = base+'.bonding.dat'
    if mode == 'clusters' and (bondfile is None or not os.path.isfile(bondfile)):
        print('-m clusters needs a bonding file (-b, or <base>.bonding.dat)')
        sys.exit(2)
    stepof = None
    if bondfile is not None:
        if not os.path.isfile(bondfile):
            print('the bonding file',bondfile,'does not exist')
            sys.exit(2)
        firstb, tb, mlabels, msizes = read_molecules(bondfile,natom)
        sb = max(1,int(round((tb or TimeStep or 1.0)/(TimeStep or 1.0))))   #UMD snapshots between two steps of the bonding file
        usnap = InitialStep+np.arange(nsnap)*Nsteps                    #UMD snapshot of each snapshot analyzed
        kof = np.clip(np.rint((usnap-firstb)/sb).astype(int),0,len(mlabels)-1)
        stepof = kof                                                    #the nearest bonding step of each snapshot analyzed
        off = np.abs(firstb+kof*sb-usnap)
        if (off > max(sb,Nsteps)).any():
            print('WARNING : the bonding file',bondfile,'does not cover the snapshots',int(usnap[off > max(sb,Nsteps)][0]),'... : the nearest steps are used')
        print('molecules from',bondfile,':',len(mlabels),'steps, one every',sb,'snapshots of the UMD file from the snapshot',firstb)
        if mode == 'clusters':#the condensed phase : the atoms of molecules of at least -l atoms ; no interface
            for t in range(nsnap):
                k = kof[t]
                states[t] = np.where(msizes[k][mlabels[k]] >= length,CONDENSED,GAS)
            print('condensed phase : the atoms of molecules of at least',length,'atoms (-m clusters)')
    #Adsorbed atoms : in the gas zone, but their molecule (bonding step nearest to the snapshot) has atoms that are not : bonded to the condensed phase,
    #even by one atom, they are not gas. They count as interface (I) for everything that follows (paths, exchanges, gas, temperature) : an atom
    #evaporates when its molecule detaches. Not with -m clusters (there a small molecule is gas by definition)
    adsorbed = np.zeros((nsnap,nel))                                    #per snapshot, adsorbed atoms of each element (zone definition)
    if stepof is not None and mode != 'clusters':
        for t in range(nsnap):
            k = stepof[t]
            lab = mlabels[k]
            isgas = states[t] == GAS
            notgas = np.bincount(lab,weights=(~isgas).astype(float),minlength=len(msizes[k]))
            ads_ = isgas&(notgas[lab] > 0)
            if ads_.any():
                states[t][ads_] = INTERFACE
                adsorbed[t] = np.bincount(typat[ads_],minlength=nel)
        for t in range(nsnap):                                          #atoms and m v^2 per element and state, with the adsorbed atoms in I
            for st in (GAS,INTERFACE,CONDENSED):
                sel = states[t] == st
                kin[t,:,st,0] = np.bincount(typat[sel],minlength=nel)
                kin[t,:,st,1] = np.bincount(typat[sel],weights=mv2s[t][sel].astype(float),minlength=nel)
        print('adsorbed atoms (in the gas zone, bonded to atoms that are not : counted in the interface), mean per snapshot :',
              ' '.join(el+' %.3g' % adsorbed[:,e].mean() for e, el in enumerate(elements) if adsorbed[:,e].any()) or 'none')
    with open(names['states'],'w') as f:
        f.write('# '+umdfile+' : per snapshot, the area of the interface (the surface where the smoothed density of the melt elements is half that of\n')
        f.write('# the bulk, '+'%.5f' % rhob+' atoms/A^3), the volume of the gas (below half), and the atoms of each element in the condensed phase (L),\n')
        f.write('# the interface (I, within '+'%.2f' % delta+' A of the surface) and the gas (G)\n')
        f.write('time_(fs)\tinterface_(nm2)\tgas_volume_(A3)\tgas_fraction\t'+'\t'.join(st+'_'+el for el in elements for st in ('L','I','G'))+'\n')
        f.write(''.join('%.2f\t%.4f\t%.2f\t%.5f\t' % (timesall[t],area[t]/100.0,gasvol[t],gasvol[t]/volumes[t])
                        +'\t'.join('%d' % kin[t,e,st,0] for e in range(nel) for st in (CONDENSED,INTERFACE,GAS))+'\n' for t in range(nsnap)))

    #Each atom followed : stays (short ones between equal states merged), committed phase, exchanges, paths through the interface
    shortest = max(1,int(np.ceil(transient/dt-1e-9))) if transient > 0 else 1
    exch = {e:{'evap':0,'cond':0,'tcond':0,'tgas':0,'pgas':[],'pcond':[]} for e in range(nel)}
    paths = {e:{p:[] for p in PATHS+('direct L-G','direct G-L')} for e in range(nel)}
    evapat, condat = [[] for e in range(nel)], [[] for e in range(nel)]
    ingas = np.zeros((nsnap,nel))
    timeinstate = np.zeros((nel,3))
    atomrows = []
    molevap, molcond = {}, {}                                           #formula : events, each atom counting 1/(atoms of its molecule)
    formcache = {}
    def molecule(a,t):#formula and size of the molecule of atom a at the bonding step nearest to the snapshot t
        k = stepof[t]
        lab = mlabels[k][a]
        if (k,lab) not in formcache:
            atoms = np.nonzero(mlabels[k] == lab)[0]
            formcache[(k,lab)] = (formula_of(atoms,typat,elements),len(atoms))
        return formcache[(k,lab)]
    for a in range(natom):
        e = typat[a]
        runs = state_runs(states[:,a],shortest)
        committed, ev, co, pth, stays = follow_atom(runs)
        if stepof is not None:
            for t in ev:                                                #the molecule it reaches the gas with
                f_, n_ = molecule(a,t)
                molevap[f_] = molevap.get(f_,0.0)+1.0/n_
            for t in co:                                                #the molecule it left the gas with
                f_, n_ = molecule(a,max(t-1,0))
                molcond[f_] = molcond.get(f_,0.0)+1.0/n_
        x = exch[e]
        x['evap'] += len(ev)
        x['cond'] += len(co)
        evapat[e] += ev
        condat[e] += co
        pos = 0
        for (s, l), c in zip(runs,committed):
            timeinstate[e,s] += l
            if c == GAS:
                ingas[pos:pos+l,e] += 1
                x['tgas'] += l
            elif c == CONDENSED:
                x['tcond'] += l
            pos += l
        for c, l, cutstart, cutend in stays:
            if not cutstart and not cutend:
                x['pgas' if c == GAS else 'pcond'].append(l)
        for p, l in pth:
            paths[e][p].append(l)
        nI = sum(l for s, l in runs if s == INTERFACE)
        nG = sum(l for (s, l), c in zip(runs,committed) if c == GAS)
        atomrows.append((a,elements[e],nG/nsnap,nI/nsnap,len(ev),len(co),'LIG'[[CONDENSED,INTERFACE,GAS].index(runs[0][0])],'LIG'[[CONDENSED,INTERFACE,GAS].index(runs[-1][0])]))
    moving = [e for e in range(nel) if exch[e]['evap']+exch[e]['cond'] > 0]
    maxlen = max([max(x['pgas']+x['pcond']) for x in exch.values() if x['pgas'] or x['pcond']]+[0])
    with open(names['stays'],'w') as f:
        f.write('# '+umdfile+' : number of complete stays (not cut by the start or the end of the run) of the atoms of each element in the gas and in\n')
        f.write('# the condensed phase (a passage through the interface belongs to the phase before it), against their duration\n')
        f.write('duration_(fs)\t'+'\t'.join('gas_'+elements[e]+'\tcondensed_'+elements[e] for e in moving)+'\n')
        hists = [(np.bincount(exch[e]['pgas'],minlength=maxlen+1),np.bincount(exch[e]['pcond'],minlength=maxlen+1)) for e in moving]
        f.write(''.join('%.2f\t' % (k*dt)+'\t'.join(str(int(g[k]))+'\t'+str(int(c[k])) for g, c in hists)+'\n' for k in range(1,maxlen+1)
                        if any(g[k] or c[k] for g, c in hists)))
    with open(names['atoms'],'w') as f:
        f.write('# '+umdfile+' : per atom (index from 0, in the order of the UMD file) : fraction of the time in the gas (its committed phase) and in the\n')
        f.write('# interface, evaporations and condensations, its state at the start and at the end (L condensed, I interface, G gas)\n')
        f.write('atom\telement\tgas_fraction\tinterface_fraction\tevaporations\tcondensations\tfirst\tlast\n')
        f.write(''.join('%d\t%s\t%.4f\t%.4f\t%d\t%d\t%s\t%s\n' % r for r in atomrows))
    pathmax = max([max(v) for pe in paths.values() for k, v in pe.items() if k in PATHS and v]+[0])
    withpaths = [e for e in range(nel) if any(paths[e][p] for p in PATHS)]
    with open(names['paths'],'w') as f:
        f.write('# '+umdfile+' : number of passages through the interface against the time spent in it, per path and element : G-I-G a gas atom\n')
        f.write('# bouncing back, G-I-L a condensation, L-I-L a melt atom visiting the surface, L-I-G an evaporation\n')
        f.write('time_in_interface_(fs)\t'+'\t'.join(p+'_'+elements[e] for e in withpaths for p in PATHS)+'\n')
        hp = [np.bincount(paths[e][p],minlength=pathmax+1) for e in withpaths for p in PATHS]
        f.write(''.join('%.2f\t' % (k*dt)+'\t'.join(str(int(h[k])) for h in hp)+'\n' for k in range(1,pathmax+1) if any(h[k] for h in hp)))

    #Fluxes per time window : events per ps (Poisson errors), net = evaporations - condensations (the gas gains atoms when > 0) ; equilibrium : from
    #the first window after which every net flux is within 2 errors of 0
    nwin = max(1,int(round(nsnap*dt/window))) if window else min(10,nsnap)
    wins = np.array_split(np.arange(nsnap),nwin)
    fluxrows, equil = [], {}
    for w in wins:
        span = len(w)*dt/1000.0                                        #ps
        row = [timesall[w[0]],timesall[w[-1]]+dt,area[w].mean()/100.0]
        for e in moving:
            ne = np.count_nonzero((np.array(evapat[e]) >= w[0])&(np.array(evapat[e]) <= w[-1]))
            nc = np.count_nonzero((np.array(condat[e]) >= w[0])&(np.array(condat[e]) <= w[-1]))
            row += [ne/span,np.sqrt(ne)/span,nc/span,np.sqrt(nc)/span,(ne-nc)/span,np.sqrt(ne+nc)/span,ingas[w,e].mean()]
        fluxrows.append(row)
    fluxrows = np.array(fluxrows)
    for k, e in enumerate(moving):
        net, err = fluxrows[:,3+7*k+4], fluxrows[:,3+7*k+5]
        ok = np.abs(net) <= 2.0*np.maximum(err,1e-300)
        equil[e] = next((j for j in range(len(wins)) if ok[j:].all()),None)
    with open(names['fluxes'],'w') as f:
        f.write('# '+umdfile+' : fluxes between the condensed phase and the gas per time window : evaporation (condensed -> gas) and condensation\n')
        f.write('# (gas -> condensed), events per ps with Poisson errors ; net = evaporation - condensation (> 0 : the gas gains atoms ; 0 within its\n')
        f.write('# error : equilibrium) ; atoms in the gas (mean over the window) ; area of the interface (mean, nm^2)\n')
        f.write('window_start_(fs)\twindow_end_(fs)\tinterface_(nm2)\t'+'\t'.join('\t'.join(c+'_'+elements[e] for c in ('evaporation_(1/ps)','evaporation_err',
                'condensation_(1/ps)','condensation_err','net_(1/ps)','net_err','in_gas')) for e in moving)+'\n')
        f.write(''.join('\t'.join('%.6g' % v for v in r)+'\n' for r in fluxrows))

    #Gas : temperature (kinetic, of the atoms in the gas ; else -T or the Temperature of the run), partial pressures n k T in the gas zone, partition
    #coefficients between the condensed zone and the gas zone, Hertz-Knudsen flux on the interface
    NG = kin[:,:,GAS,0]                                                 #snapshot x element : atoms in the gas
    NL = kin[:,:,CONDENSED,0]
    #snapshots whose velocities are 0 or off by more than a factor 2 (kinetic temperature of all the atoms against the median : first snapshots of
    #the parsers, learning steps of ML runs, a faulty snapshot) are left out of the temperature : one of them dominates a sum of v^2
    Tall = kin[:,:,:,1].sum(axis=(1,2))*AMU_A2FS2/(3.0*KB*natom)
    med = np.median(Tall[Tall > 0]) if (Tall > 0).any() else 0.0
    goodv = (Tall > 0.5*med)&(Tall < 2.0*med)
    nbadv = int(nsnap-goodv.sum())
    mv2G = kin[goodv][:,:,GAS,1].sum(axis=0)
    if Tuser is not None:
        Tgas, Tsource = Tuser, 'given (-T)'
    elif med > 0 and NG[goodv].sum() > 0:
        Tgas = mv2G.sum()*AMU_A2FS2/(3.0*KB*NG[goodv].sum())
        Tsource = 'kinetic, of the atoms in the gas'+(' ; %d snapshots with velocities 0 or off by more than a factor 2 left out' % nbadv if nbadv else '')
    else:
        temps = []
        with open(umdpf.real_path(umdfile)) as fu:
            for line in fu:
                if line.startswith('Temperature'):
                    temps.append(float(line.split()[1]))
        temps = temps[InitialStep::Nsteps]
        Tgas, Tsource = (float(np.mean(temps)) if temps else np.nan), 'Temperature of the run (no velocities)'
    Tel = np.where(NG[goodv].sum(axis=0) > 0,mv2G*AMU_A2FS2/(3.0*KB*np.maximum(NG[goodv].sum(axis=0),1)),np.nan) if Tsource.startswith('kinetic') else np.full(nel,np.nan)
    VG, VL = zones[:,0], zones[:,1]
    with np.errstate(invalid='ignore',divide='ignore'):
        Pel = NG*KB*Tgas/VG[:,None]*EV_A3_GPA                           #GPa, snapshot x element
        cG = NG/VG[:,None]
        cL = NL/VL[:,None]
    with open(names['gas'],'w') as f:
        f.write('# '+umdfile+' : per snapshot, the volume of the gas zone (all the density within '+'%.2f' % delta+' A below half the bulk) and of the condensed\n')
        f.write('# zone (all above), the atoms of each element in the gas, and their partial pressures n k T (ideal gas, T = '+'%.1f' % Tgas+' K, '+Tsource+')\n')
        f.write('time_(fs)\tgas_zone_(A3)\tcondensed_zone_(A3)\tP_total_(GPa)\t'+'\t'.join('N_gas_'+el for el in elements)+'\t'+'\t'.join('P_'+el+'_(GPa)' for el in elements)+'\n')
        f.write(''.join('%.2f\t%.2f\t%.2f\t%.6g\t' % (timesall[t],VG[t],VL[t],np.nansum(Pel[t]))+'\t'.join('%d' % v for v in NG[t])+'\t'
                        +'\t'.join('%.6g' % v for v in Pel[t])+'\n' for t in range(nsnap)))

    #The gas by molecules, two definitions : zone (all the atoms of the molecule in the gas zone, G ; volume : the gas zone) and surface (all its atoms
    #below the half-density surface ; volume : below the surface). An atom outside whose molecule is not all outside is adsorbed : not gas. Molecules
    #from the bonding file at each of its steps (the nearest snapshot analyzed), else every atom its own molecule at every snapshot. Pressure : ideal
    #N k T/V and modified Van der Waals N k T/V - (sum N_i a_i^1/2)^2/(NA^2 V^2), T the kinetic temperature of the atoms of the gas molecules (pooled
    #over the run, snapshots with faulty velocities left out ; -T)
    DEFS = ('zone','surface')
    vdw = vdw_table(vdwextra)
    if stepof is not None:
        steps = sorted(set(stepof.tolist()))
        tof = {k:int(np.nonzero(stepof == k)[0][np.argmin(np.abs(InitialStep+np.nonzero(stepof == k)[0]*Nsteps-(firstb+k*sb)))]) for k in steps}
        samples = [(tof[k],mlabels[k],msizes[k],k) for k in steps]
    else:
        single = np.arange(natom)
        samples = [(t,single,np.ones(natom,dtype=np.int32),None) for t in range(nsnap)]
    ns = len(samples)
    formcache = {}
    def formula(k,lab,m):
        if k is None:
            return elements[typat[m]]+'_1'
        if (k,m) not in formcache:
            formcache[(k,m)] = formula_of(np.nonzero(lab == m)[0],typat,elements)
        return formcache[(k,m)]
    gascount = {d:{} for d in DEFS}                                     #formula : number in the gas at each sample
    Ngas = {d:np.zeros(ns) for d in DEFS}
    Sa = {d:np.zeros(ns) for d in DEFS}                                 #sum of N_i a_i^1/2 (L bar^1/2/mol)
    Nnoa = {d:np.zeros(ns) for d in DEFS}                               #gas molecules without a
    ads = {d:np.zeros((ns,nel)) for d in DEFS}                          #adsorbed atoms of each element
    mv2gas, ngasat = {d:0.0 for d in DEFS}, {d:0 for d in DEFS}
    Vs = {'zone':np.array([zones[t,0] for t, *_ in samples]),'surface':np.array([gasvol[t] for t, *_ in samples])}
    tsample = np.array([t for t, *_ in samples])
    for i, (t,lab,sizes,k) in enumerate(samples):
        for d, flag in (('zone',states[t] == GAS),('surface',outside[t])):
            nout = np.bincount(lab,weights=flag.astype(float),minlength=len(sizes))
            gasm = (nout == sizes)&(sizes > 0)
            ingas = gasm[lab]
            ads[d][i] = adsorbed[t] if d == 'zone' and stepof is not None and mode != 'clusters' else np.bincount(typat[flag&~ingas],minlength=nel)
            if goodv[t]:
                mv2gas[d] += float(mv2s[t][ingas].sum())
                ngasat[d] += int(ingas.sum())
            for m in np.nonzero(gasm)[0]:
                f_ = formula(k,lab,m)
                gascount[d].setdefault(f_,np.zeros(ns))[i] += 1
                a_ = vdw.get(composition(f_))
                if a_ is None:
                    Nnoa[d][i] += 1
                else:
                    Sa[d][i] += np.sqrt(a_)
            Ngas[d][i] = gasm.sum()
    Tdef, Pid, Pvdw = {}, {}, {}
    for d in DEFS:
        if Tuser is not None or not Tsource.startswith('kinetic') or ngasat[d] == 0:
            Tdef[d] = Tgas
        else:
            Tdef[d] = mv2gas[d]*AMU_A2FS2/(3.0*KB*ngasat[d])
        with np.errstate(invalid='ignore',divide='ignore'):
            Pid[d] = Ngas[d]*KB*Tdef[d]/Vs[d]*EV_A3_GPA
            Pvdw[d] = Pid[d]-0.1*Sa[d]**2/(NA**2*(Vs[d]*1e-30)**2)*1e-9
    with open(names['pressure'],'w') as f:
        f.write('# '+umdfile+' : the gas by molecules'+(' (connected bonds of '+bondfile+', at each of its steps)' if stepof is not None else ' (no bonding file : every atom')
                +('' if stepof is not None else ' its own particle, every snapshot)')+', two definitions : zone (all the atoms of the molecule in the gas zone,\n')
        f.write('# the density within '+'%.2f' % delta+' A all below half the bulk ; volume of that zone) and surface (all its atoms below the half-density surface ;\n')
        f.write('# volume below it) ; adsorbed : atoms outside whose molecule is not all outside ; pressures ideal N k T/V and modified Van der Waals\n')
        f.write('# N k T/V - (sum N_i a_i^1/2)^2/(NA^2 V^2) (T = '+', '.join(d+' '+'%.1f' % Tdef[d] for d in DEFS)+' K)\n')
        f.write('time_(fs)\t'+'\t'.join('\t'.join(c+'_'+d for c in ('V_gas_(A3)','molecules','adsorbed_atoms','P_ideal_(GPa)','P_vdW_(GPa)')) for d in DEFS)+'\n')
        f.write(''.join('%.2f\t' % timesall[tsample[i]]+'\t'.join('%.2f\t%d\t%d\t%.6g\t%.6g' % (Vs[d][i],Ngas[d][i],ads[d][i].sum(),Pid[d][i],Pvdw[d][i]) for d in DEFS)+'\n'
                        for i in range(ns)))

    #Molecules in the gas (zone definition, as before ; also the surface count) : partial pressures n k T in the gas zone, events
    molrows = []
    if stepof is not None:
        Tps = nsnap*dt/1000.0
        VGs = Vs['zone']
        allforms = sorted(set(gascount['zone'])|set(molevap)|set(molcond),key=lambda f_: (-(gascount['zone'].get(f_,np.zeros(1)).mean()),f_))
        totalgas = sum(c for c in gascount['zone'].values()) if gascount['zone'] else np.zeros(ns)
        for f_ in allforms:
            c = gascount['zone'].get(f_,np.zeros(ns))
            with np.errstate(invalid='ignore',divide='ignore'):
                P = c*KB*Tgas/VGs*EV_A3_GPA
            sP = statistics(P) if len(P) >= 8 and np.isfinite(P).all() else (np.nanmean(P),np.nan,np.nan,False,0.0)
            ne, nc = molevap.get(f_,0.0), molcond.get(f_,0.0)
            a_ = vdw.get(composition(f_))
            molrows.append((f_,c.mean(),c.sum()/totalgas.sum() if totalgas.sum() > 0 else np.nan,sP[0],sP[2],ne,nc,ne/Tps,nc/Tps,
                            (ne-nc)/np.sqrt(ne+nc) if ne+nc > 0 else 0.0,gascount['surface'].get(f_,np.zeros(ns)).mean(),a_ if a_ is not None else np.nan))
        with open(names['molecules'],'w') as f:
            f.write('# '+umdfile+' : molecules (connected bonds of '+bondfile+') : mean number in the gas (all their atoms in the gas), fraction of the gas\n')
            f.write('# molecules, partial pressure n k T in the gas zone (T = '+'%.1f' % Tgas+' K), evaporations and condensations (each atom event counted 1/(atoms of\n')
            f.write('# its molecule) : molecules), fluxes per ps, net = (evaporations - condensations)/sqrt(sum) ; mean number in the gas below the half-density\n')
            f.write('# surface ; Van der Waals a (L^2 bar/mol^2 ; nan : none, ideal)\n')
            f.write('molecule\tin_gas\tgas_fraction\tP_(GPa)\tP_err\tevaporations\tcondensations\tevaporation_(1/ps)\tcondensation_(1/ps)\tnet\tin_gas_surface\ta_vdW\n')
            f.write(''.join('%s\t%.5g\t%.5g\t%.5g\t%.2g\t%.4g\t%.4g\t%.4g\t%.4g\t%.1f\t%.5g\t%.4g\n' % r for r in molrows))

    #Summary
    g = lambda v, f: (f % v) if np.isfinite(v) else 'nan'
    sta = statistics(area/100.0)
    stg = statistics(gasvol/volumes)
    out = []
    out.append('# '+umdfile+' : '+str(nsnap)+' snapshots from '+str(InitialStep)+' every '+str(Nsteps)+' ; '+'%.2f' % dt+' fs between them')
    out.append('# condensed phase : smoothed density of the elements '+','.join(el for el, v in zip(elements,isvol) if not v)+' (Gaussians of '+'%g' % sigma+' A, grid '
               +'%dx%dx%d' % tuple(shape)+') ; surface at half the bulk density '+'%.5f' % rhob+' atoms/A^3 ; interface : within '+'%.2f' % delta+' A of it')
    out.append('# stays shorter than '+'%g' % transient+' fs between two stays in the same state merged ; errors of the means from their autocorrelation time')
    out.append('')
    out.append('Interface (means over the run ; the area from the gradient of the density at the surface)')
    out.append('area_(nm2)\t'+'%.5g' % sta[0]+'\t'+('>' if sta[3] else '')+'%.2g' % sta[2])
    out.append('gas_volume_fraction\t'+'%.5g' % stg[0]+'\t'+('>' if stg[3] else '')+'%.2g' % stg[2])
    out.append('bulk_density_(atoms/A3)\t'+'%.6g' % rhob)
    out.append('')
    out.append('Fraction of the time of the atoms of each element in the condensed phase (L), the interface (I) and the gas (G)')
    out.append('element\tcondensed\tinterface\tgas')
    for e, el in enumerate(elements):
        tot = timeinstate[e].sum()
        out.append(el+'\t'+'%.4f\t%.4f\t%.4f' % (timeinstate[e,CONDENSED]/tot,timeinstate[e,INTERFACE]/tot,timeinstate[e,GAS]/tot))
    out.append('')
    out.append('Paths through the interface (complete passages) : number, mean time in the interface (fs, error std/sqrt(n)) ; G-I-G a gas atom bouncing')
    out.append('back, G-I-L a condensation, L-I-L a melt atom visiting the surface, L-I-G an evaporation ; direct : a change of phase between two snapshots')
    out.append('without a stay in the interface (sampling too coarse) ; sticking = G-I-L/(G-I-L + G-I-G), the probability that a gas atom reaching the')
    out.append('interface condenses ; escape = L-I-G/(L-I-G + L-I-L), the probability that a melt atom reaching it evaporates (binomial errors)')
    out.append('element\t'+'\t'.join(p+'\t'+p+'_time_(fs)\t'+p+'_time_err' for p in PATHS)+'\tdirect_L-G\tdirect_G-L\tsticking\tsticking_err\tescape\tescape_err')
    for e, el in enumerate(elements):
        row = el
        for p in PATHS:
            v = np.array(paths[e][p])*dt
            row += '\t%d\t%s\t%s' % (len(v),g(v.mean(),'%.4g') if len(v) else 'nan',g(v.std()/np.sqrt(len(v)),'%.2g') if len(v) > 1 else 'nan')
        row += '\t%d\t%d' % (len(paths[e]['direct L-G']),len(paths[e]['direct G-L']))
        for a_, b_ in (('G-I-L','G-I-G'),('L-I-G','L-I-L')):
            na, nb = len(paths[e][a_]), len(paths[e][b_])
            pr = na/(na+nb) if na+nb > 0 else np.nan
            row += '\t'+g(pr,'%.4f')+'\t'+(g(np.sqrt(pr*(1-pr)/(na+nb)),'%.2g') if na+nb > 0 else 'nan')
        out.append(row)
    out.append('')
    out.append('Exchange between the gas and the condensed phase (an atom changes phase when it goes from one to the other, through the interface or')
    out.append('not ; stays shorter than '+'%g' % transient+' fs merged) : evaporations (condensed -> gas) and condensations, rate constants')
    out.append('k = events / time the atoms of the element spent in the phase they leave (per ps, Poisson errors), mean stays (complete ones, errors')
    out.append('std / sqrt(stays)), fraction of the time in the gas ; net = evaporations - condensations in their combined error (|net| > 3 : the')
    out.append('partition is still changing)')
    out.append('element\tevaporations\tcondensations\tk_evaporation_(1/ps)\terror\tk_condensation_(1/ps)\terror\tmean_stay_condensed_(fs)\terror'
               '\tmean_stay_gas_(fs)\terror\tgas_fraction\tnet')
    for e, el in enumerate(elements):
        x = exch[e]
        tc, tg = x['tcond']*dt/1000.0, x['tgas']*dt/1000.0              #ps
        ke = x['evap']/tc if tc > 0 else np.nan
        kc = x['cond']/tg if tg > 0 else np.nan
        pc, pg = np.array(x['pcond'])*dt, np.array(x['pgas'])*dt
        net = (x['evap']-x['cond'])/np.sqrt(x['evap']+x['cond']) if x['evap']+x['cond'] > 0 else 0.0
        out.append(el+'\t%d\t%d\t' % (x['evap'],x['cond'])+'%.4g\t%.2g\t%.4g\t%.2g\t' % (ke,np.sqrt(x['evap'])/tc if tc > 0 else np.nan,kc,np.sqrt(x['cond'])/tg if tg > 0 else np.nan)
                   +('%.4g\t%.2g\t' % (pc.mean(),pc.std()/np.sqrt(len(pc)) if len(pc) > 1 else np.nan) if len(pc) else 'nan\tnan\t')
                   +('%.4g\t%.2g\t' % (pg.mean(),pg.std()/np.sqrt(len(pg)) if len(pg) > 1 else np.nan) if len(pg) else 'nan\tnan\t')
                   +'%.4f\t%.1f' % (tg/(tc+tg) if tc+tg > 0 else np.nan,net))
    out.append('')
    out.append('Fluxes in time (fluxes.dat, '+str(len(wins))+' windows of '+'%.4g' % (nsnap*dt/len(wins))+' fs) : evaporation and condensation fluxes over the run (events per ps,')
    out.append('and per ps and nm^2 of interface, mean '+'%.4g' % sta[0]+' nm^2), and the time from which the two are equal : every window from then on has')
    out.append('a net flux within 2 errors of 0 (equilibrium) ; not reached : the last window still has a net flux')
    out.append('element\tevaporation_(1/ps)\tcondensation_(1/ps)\tevaporation_(1/ps/nm2)\tcondensation_(1/ps/nm2)\tequilibrium_from_(fs)')
    T = nsnap*dt/1000.0
    for e in moving:
        x = exch[e]
        when = equil[e]
        out.append(elements[e]+'\t%.4g\t%.4g\t%.4g\t%.4g\t' % (x['evap']/T,x['cond']/T,x['evap']/T/sta[0] if sta[0] > 0 else np.nan,
                   x['cond']/T/sta[0] if sta[0] > 0 else np.nan)+('%.2f' % fluxrows[when][0] if when is not None else 'not reached'))
    out.append('')
    stv = lambda x: statistics(x[np.isfinite(x)]) if np.isfinite(x).sum() >= 8 else (np.nan,np.nan,np.nan,False,0.0)
    Ptot = np.nansum(Pel,axis=1)
    stP = stv(Ptot)
    out.append('Gas : temperature '+g(Tgas,'%.1f')+' K ('+Tsource+') ; zones (means) : gas '+g(VG.mean(),'%.5g')+' A^3, condensed '+g(VL.mean(),'%.5g')+' A^3')
    out.append('(all the density within '+'%.2f' % delta+' A below, or above, half the bulk) ; total pressure of the atoms in the gas '+g(stP[0],'%.5g')+' +/- '
               +g(stP[2],'%.2g')+' GPa = '+g(stP[0]*1e4,'%.5g')+' bar (ideal gas n k T ; the atoms counted as independent particles : an upper bound when')
    out.append('they are bound in molecules). Per element : temperature (kinetic, equipartition check), atoms in the gas, mole fraction of the gas,')
    out.append('partial pressure, densities in the condensed and gas zones, partition coefficient D = c(condensed)/c(gas) and K = x(condensed)/x(gas) ;')
    out.append('Hertz-Knudsen : the flux with which the gas hits the interface, J = n (k T / 2 pi m)^1/2 (per ps and nm^2), against the arrivals of gas')
    out.append('atoms at the interface counted (G-I-G, G-I-L, direct G-L ; their ratio : 1 for a gas of free atoms) and the condensation coefficient')
    out.append('alpha = condensation flux / J, to compare with the sticking coefficient of the paths')
    out.append('element\tT_gas_(K)\tN_gas\tx_gas\tP_(GPa)\tP_err\tc_condensed_(1/A3)\tc_gas_(1/A3)\tD\tD_err\tK\tK_err\tJ_HK_(1/ps/nm2)\tarrivals_(1/ps/nm2)'
               '\tarrivals/J_HK\talpha_HK\tsticking')
    xG_all = NG.sum(axis=1)
    xL_all = NL.sum(axis=1)
    Tps, A = nsnap*dt/1000.0, sta[0]
    for e, el in enumerate(elements):
        sP = stv(Pel[:,e])
        sL, sG = stv(cL[:,e]), stv(cG[:,e])
        D = sL[0]/sG[0] if sG[0] > 0 else np.nan
        Derr = D*np.hypot(sL[2]/sL[0] if sL[0] > 0 else 0.0,sG[2]/sG[0] if sG[0] > 0 else 0.0) if np.isfinite(D) else np.nan
        with np.errstate(invalid='ignore',divide='ignore'):
            xl, xg = NL[:,e]/xL_all, NG[:,e]/xG_all
        sxl, sxg = stv(xl), stv(xg)
        K = sxl[0]/sxg[0] if sxg[0] > 0 else np.nan
        Kerr = K*np.hypot(sxl[2]/sxl[0] if sxl[0] > 0 else 0.0,sxg[2]/sxg[0] if sxg[0] > 0 else 0.0) if np.isfinite(K) else np.nan
        n = sG[0]
        vth = np.sqrt(KB*Tgas/(2.0*np.pi*masses[e]*AMU_A2FS2)) if np.isfinite(Tgas) else np.nan     #A/fs
        J = n*vth*1e5 if np.isfinite(n) else np.nan                    #atoms/A^2/fs -> per ps and nm^2
        arrivals = len(paths[e]['G-I-G'])+len(paths[e]['G-I-L'])+len(paths[e]['direct G-L'])
        arr = arrivals/Tps/A if A > 0 else np.nan
        cond = exch[e]['cond']/Tps/A if A > 0 else np.nan
        ng, nl = len(paths[e]['G-I-L']), len(paths[e]['G-I-G'])
        out.append(el+'\t'+g(Tel[e],'%.1f')+'\t'+'%.4g' % NG[:,e].mean()+'\t'+g(sxg[0],'%.4f')+'\t'+g(sP[0],'%.5g')+'\t'+g(sP[2],'%.2g')+'\t'+g(sL[0],'%.5g')+'\t'
                   +g(sG[0],'%.5g')+'\t'+g(D,'%.4g')+'\t'+g(Derr,'%.2g')+'\t'+g(K,'%.4g')+'\t'+g(Kerr,'%.2g')+'\t'+g(J,'%.4g')+'\t'+g(arr,'%.4g')+'\t'
                   +g(arr/J if J > 0 else np.nan,'%.3f')+'\t'+g(cond/J if J > 0 else np.nan,'%.4f')+'\t'+g(ng/(ng+nl) if ng+nl > 0 else np.nan,'%.4f'))
    out.append('')
    if molrows:
        out.append('Molecules (connected bonds of '+bondfile+') : mean number in the gas (all their atoms in the gas), fraction of the gas molecules, partial')
        out.append('pressure n k T in the gas zone (GPa, error from the autocorrelation time), evaporations and condensations of molecules (each atom event')
        out.append('counted 1/(atoms of its molecule)), fluxes per ps, net = (evaporations - condensations)/sqrt(sum) ; the 10 most abundant in the gas')
        out.append('molecule\tin_gas\tgas_fraction\tP_(GPa)\tP_err\tevaporations\tcondensations\tevaporation_(1/ps)\tcondensation_(1/ps)\tnet')
        out += ['%s\t%.5g\t%.5g\t%.5g\t%.2g\t%.4g\t%.4g\t%.4g\t%.4g\t%.1f' % r[:10] for r in molrows[:10]]
        out.append('total_molecules_in_gas\t'+'%.5g' % sum(r[1] for r in molrows)+'\tpressure_of_the_molecules_(GPa)\t'+'%.5g' % sum(r[3] for r in molrows if np.isfinite(r[3])))
        out.append('')
    out.append('Gas pressure by molecules ('+('connected bonds of '+bondfile if stepof is not None else 'no bonding file : every atom its own particle')+', '+str(ns)
               +' samples), modified Van der Waals law P = N k T/V - a N^2/(NA^2 V^2) (Davis and Caracas 2024, Commun. Earth Environ. 5, 344),')
    out.append('a N^2 = (sum_i N_i a_i^1/2)^2 with a_i of the CRC Handbook (-A for others ; without a : ideal) ; two definitions of the gas : zone (all the')
    out.append('atoms of a molecule in the gas zone ; its volume) and surface (all its atoms below the half-density surface ; the volume below it) ; adsorbed :')
    out.append('atoms outside whose molecule is not all outside (molecules bonded to the condensed phase, and with the surface definition also the outermost')
    out.append('atoms of the condensed phase itself), not gas ; T : kinetic, of the atoms of the gas molecules ; errors of the means from the autocorrelation time')
    out.append('definition\tV_gas_(A3)\tV_err\tmolecules\tmolecules_err\tadsorbed_atoms\tT_gas_(K)\tP_ideal_(GPa)\tP_ideal_err\tP_vdW_(GPa)\tP_vdW_err\twithout_a')
    for d in DEFS:
        sV, sN, sI, sW = stv(Vs[d]), stv(Ngas[d]), stv(Pid[d]), stv(Pvdw[d])
        out.append(d+'\t'+g(sV[0],'%.5g')+'\t'+g(sV[2],'%.2g')+'\t'+g(sN[0],'%.5g')+'\t'+g(sN[2],'%.2g')+'\t'+'%.4g' % ads[d].sum(axis=1).mean()+'\t'
                   +g(Tdef[d],'%.1f')+'\t'+g(sI[0],'%.5g')+'\t'+g(sI[2],'%.2g')+'\t'+g(sW[0],'%.5g')+'\t'+g(sW[2],'%.2g')+'\t'
                   +g(Nnoa[d].sum()/Ngas[d].sum() if Ngas[d].sum() > 0 else np.nan,'%.3f'))
    out.append('adsorbed atoms per element (means) : '+' ; '.join(d+' '+' '.join(el+' %.3g' % ads[d][:,e].mean() for e, el in enumerate(elements) if ads[d][:,e].any())
               for d in DEFS))
    noa = {}
    for d in DEFS:
        for f_, c in gascount[d].items():
            if composition(f_) not in vdw:
                noa[f_] = max(noa.get(f_,0.0),c.mean())
    if noa:
        out.append('gas species without a (ideal) : '+' '.join(f_+' %.3g' % n for f_, n in sorted(noa.items(),key=lambda x: (-x[1],x[0]))[:10]))
    out.append('')
    with open(names['summary'],'w') as f:
        f.write('\n'.join(out)+'\n')
    print()
    print('interface :','%.4g' % sta[0],'nm^2 ; gas','%.3f' % stg[0],'of the volume ; gas temperature','%.1f' % Tgas,'K ('+Tsource+') ; pressure of the gas atoms',
          '%.4g' % stP[0],'GPa')
    print('gas pressure (modified Van der Waals, molecules) :',' ; '.join(d+' '+'%.4g' % stv(Pvdw[d])[0]+' +/- '+'%.2g' % stv(Pvdw[d])[2]+' GPa' for d in DEFS))
    print('evaporations / condensations :','  '.join(elements[e]+' %d/%d' % (exch[e]['evap'],exch[e]['cond']) for e in moving) or 'none')
    print('Files',', '.join(os.path.basename(n) for n in names.values()),'written in the folder',outdir)
    print('runtime : ','%.1f' % (time.time()-start),' s')


if __name__ == "__main__":
    if not umdpf.run_on_all(sys.argv[1:],main):#-f all : every UMD file of the folder
        main(sys.argv[1:])
