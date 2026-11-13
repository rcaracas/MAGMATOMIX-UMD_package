#!/usr/bin/env python3
###
#AUTHORS: RAZVAN CARACAS
###
#Orientations of the molecular units of a run (H2O, NH3, CH4, OH, SO4, CO3...) : a central atom with its bonded neighbours, from the bonding file.
#Correlations between the orientations of different units (against their distance, first shell, Kirkwood factor, total polarization), self
#correlation of the orientation of a unit in time, map of the orientations in the frame of the cell and its entropy.
#The orientation of a unit is a set of unit vectors V : its dipole axis (the sum of its bond vectors : H2O, NH3, OH), or, for units without a dipole
#(CH4, SO4, CO3), its bonds. The correlation of order l between two sets A and B is the mean of P_l(a.b) over the vectors a of A and b of B :
#for dipoles, l = 1 is the cosine of their angle ; a tetrahedron has no l = 1 or 2 by symmetry, its first orders are 3 and 4.

#Entropy of the units by the two-phase thermodynamic model (2PT ; Lin, Blanco, Goddard, J. Chem. Phys. 119, 11792 (2003) ; for molecules, Lin, Maiti,
#Goddard, J. Phys. Chem. B 114, 8191 (2010)) : from the velocities, the densities of states of the translation of the units (centre of mass), of their
#rotation (angular velocities in the body frame, weighted by the moments of inertia) and of their internal motions ; translation and rotation each
#split into a gas-like part (hard spheres ; rigid rotor) and a solid-like part (quantum oscillators), with their own fluidicity. The volume of a unit
#comes from the user, e.g. from a Bader analysis (ACF.dat files of snapshots of the run), else from the cell volume per unit.

import sys, getopt, os.path, time, re, itertools
import numpy as np
from scipy.special import lpmv, factorial
from scipy.optimize import linear_sum_assignment
import umd_processes_fast as umdpf
import analyze_gofr
from vibrations_umd import acf_fft, write_columns, vdos, read_temperatures, CM1, KB, AMU_A2FS2
import vdos_thermo

LMAX = 4                #orders of the correlations : 1 and 2 for a dipole axis, 1 to 4 for the bonds (those the shape of the unit allows)
ACTIVE = 0.05           #an order is kept for the bonds when its self value (the mean of P_l between the bonds of a unit) exceeds this
POLAR = 0.2             #a unit has a dipole axis when |sum of its unit bond vectors| / number of bonds exceeds this (water 0.61, NH3 0.37, CH4 0)
MAXEL = 4.0E6           #numbers per array in the batches of the FFTs and of the pairs
GAP = 50.0              #fs : default of -g
NCOS = 40               #bins of the cosine between the dipoles of neighbours
PAIRSNAPS = 2000        #snapshots used for the pair correlations by default
BOHR = 0.529177210903   #A
MATCHTOL = 0.02         #A : largest distance between an atom of an ACF.dat file and the same atom in its snapshot of the UMD file
LINEAR = 0.05           #a unit is linear when its smallest principal moment of inertia / the largest is below this (OH : 0 ; water 0.33)
SYMTOL = 0.15           #a permutation of the neighbours is a rotation of the unit when it maps the unit onto itself within this fraction of a bond
MAXPERM = 5040          #most permutations of the neighbours tested for the symmetry number


def harmonics(v,l):#The 2l+1 real functions R_lm of the unit vectors v (..., 3) : their scalar product for two vectors a, b is P_l(a.b) (addition theorem)
    z = np.clip(v[...,2],-1.0,1.0)
    phi = np.arctan2(v[...,1],v[...,0])
    out = [lpmv(0,l,z)]
    for m in range(1,l+1):
        f = np.sqrt(2.0*factorial(l-m)/factorial(l+m))*lpmv(m,l,z)
        out += [f*np.cos(m*phi), f*np.sin(m*phi)]
    return np.stack(out,axis=-1)


def lives_of(valid,maxgap,minlen):#valid (snapshot, unit) : the central atom has the composition asked. A life of a unit : its valid snapshots, with the
    #gaps of at most maxgap snapshots (it goes on as the same molecule, left out during the gap) ; a longer gap ends it. Lives of fewer than minlen
    #snapshots are dropped. Returns the lives (unit, first, last) and the id of the life of each unit at each snapshot (-1 outside a life)
    lid = -np.ones(valid.shape,dtype=np.int32)
    lives = []
    for i in range(valid.shape[1]):
        edges = np.flatnonzero(np.diff(np.r_[0,valid[:,i].astype(np.int8),0]))
        if len(edges) == 0:
            continue
        starts, ends = edges[0::2], edges[1::2]-1
        split = (starts[1:]-ends[:-1]-1) > maxgap
        for s,e in zip(starts[np.r_[True,split]],ends[np.r_[split,True]]):
            if e-s+1 >= minlen:
                lid[s:e+1,i] = len(lives)
                lives.append((i,int(s),int(e)))
    return lives, lid


def masked_acf(batches,nlags,slices):#Sum over the batches (x (n, L, components), zero where masked ; w (n, L), 1 where present) of the products
    #x(t).x(t+lag), by group of components (slices), and of the numbers of pairs of times w(t) w(t+lag) : numerators (group, lag), denominator (lag)
    num = np.zeros((len(slices),nlags))
    den = np.zeros(nlags)
    for x,w in batches:
        L = x.shape[1]
        nl = min(nlags,L)
        F = np.fft.rfft(x,n=2*L,axis=1)
        acf = np.fft.irfft(F*F.conj(),n=2*L,axis=1)[:,:nl]
        for q,sl in enumerate(slices):
            num[q,:nl] += acf[:,:,sl].sum(axis=(0,2))
        W = np.fft.rfft(w,n=2*L,axis=1)
        den[:nl] += np.fft.irfft(W*W.conj(),n=2*L,axis=1)[:,:nl].sum(axis=0)
    return num, np.rint(den)


def by_length(segments):#(key, first, last) grouped by their length, in pieces of at most MAXEL numbers for 8 components
    groups = {}
    for key,s,e in segments:
        groups.setdefault(e-s+1,[]).append((key,s))
    for L,group in groups.items():
        step = max(1,int(MAXEL//(16*L)))
        for k in range(0,len(group),step):
            yield L, group[k:k+step]


def unit_batches(segments,Q,used):#Batches of masked_acf for the units : segments (unit, first, last) of their lives
    for L,group in by_length(segments):
        yield np.stack([Q[s:s+L,i] for i,s in group]), np.stack([used[s:s+L,i] for i,s in group]).astype(float)


def bond_batches(segments,rows,btime,bharm):#Batches of masked_acf for the tagged bonds : segments (bond group, first, last), rows[group] the entries of the
    #bond (central atom, neighbour, life) in the arrays btime (snapshot) and bharm (harmonics of its direction)
    for L,group in by_length(segments):
        x = np.zeros((len(group),L,bharm.shape[1]))
        w = np.zeros((len(group),L))
        for j,(g,s) in enumerate(group):
            t = btime[rows[g]]
            keep = (t >= s) & (t < s+L)
            x[j,t[keep]-s] = bharm[rows[g]][keep]
            w[j,t[keep]-s] = 1.0
        yield x, w


def clip_segments(segments,t0,t1):#The parts of the segments (key, first, last) inside the snapshots t0 .. t1-1, at least 2 snapshots long
    out = []
    for key,s,e in segments:
        s, e = max(s,t0), min(e,t1-1)
        if e > s:
            out.append((key,s,e))
    return out


def normalized(num,den):#C(lag) / C(0) ; nan where no pair of times
    with np.errstate(invalid='ignore',divide='ignore'):
        c = np.where(den > 0,num/np.where(den > 0,den,1.0),np.nan)
    return c/c[...,0:1]


def block_error(values):#Standard error of the mean over the blocks (axis 0) ; nan where a block has none
    values = np.asarray(values,dtype=float)
    return np.std(values,axis=0,ddof=1)/np.sqrt(values.shape[0]) if values.shape[0] > 1 else np.full(values.shape[1:],np.nan)


def decay_time(t,c):#Time at which c falls to 1/e (linear interpolation) and mean of c over the last fifth of the window ; None when it does not fall
    c = np.asarray(c,dtype=float)
    ok = np.isfinite(c)
    t, c = t[ok], c[ok]
    plateau = float(np.mean(c[-max(1,len(c)//5):])) if len(c) else np.nan
    below = np.flatnonzero(c < np.exp(-1.0))
    if len(below) == 0 or below[0] == 0:
        return None, plateau
    k = below[0]
    return float(t[k-1]+(c[k-1]-np.exp(-1.0))*(t[k]-t[k-1])/(c[k-1]-c[k])), plateau


def map_entropy(counts):#Entropy of a distribution on the sphere, in kB, relative to the uniform one (0 : free rotation ; -ln(number of bins) : one bin),
    #from the counts in bins of equal solid angle
    total = counts.sum()
    if total == 0:
        return np.nan
    p = counts[counts > 0]/total
    return float(-(p*np.log(p*counts.size)).sum())


#*** Volume of a unit : a number, Bader volumes per element, or ACF.dat files of the Bader program matched to snapshots of the run

def volume_option(arg):#-V : ('number', A^3), ('elements', {element: A^3}) or ('files', the ACF.dat files)
    try:
        return 'number', float(arg)
    except ValueError:
        pass
    parts = arg.split(',')
    if all(re.fullmatch(r'[A-Z][a-z]?:[0-9.eE+-]+',p) for p in parts):
        return 'elements', {p.split(':')[0]:float(p.split(':')[1]) for p in parts}
    files = []
    for part in parts:
        if os.path.isdir(part):#every file with ACF in its name, in the folder and below (e.g. one folder per snapshot)
            found = sorted(os.path.join(d,n) for d,_,names in os.walk(part) for n in names if 'acf' in n.lower())
        elif any(ch in part for ch in '*?['):
            import glob
            found = sorted(glob.glob(part))
        else:
            found = [part] if os.path.isfile(part) else []
        if not found:
            print('-V : no ACF file found for',part)
            sys.exit()
        files += found
    return 'files', files


def read_acf(name):#ACF.dat of the Bader program (Henkelman group) : positions, atomic volumes (A or bohr, as the density it read) and vacuum volume
    pos, vol = [], []
    vacuum = 0.0
    with open(name,'r') as f:
        for line in f:
            entry = line.split()
            if len(entry) == 7 and entry[0].isdigit():
                pos.append([float(x) for x in entry[1:4]])
                vol.append(float(entry[6]))
            elif line.strip().startswith('VACUUM VOLUME'):
                vacuum = float(entry[-1])
    return np.array(pos).reshape(-1,3), np.array(vol), vacuum


def snapshot_numbers(name):#The numbers in the path of a file, those of its name first (the last first), then those of its folders : its snapshot
    numbers = []
    for part in [os.path.basename(name)]+os.path.dirname(name).split(os.sep)[::-1]:
        numbers += [int(x) for x in re.findall(r'\d+',part)][::-1]
    return list(dict.fromkeys(numbers))


def nearest_atoms(sa,su,cell,shift):#For each reduced position of sa (+ shift) the nearest of su, nearest images : its index and the distance (A)
    j = np.empty(len(sa),dtype=np.int64)
    d = np.empty(len(sa))
    step = max(1,int(MAXEL//(3*len(su))))
    for a in range(0,len(sa),step):
        ds = sa[a:a+step,None,:]+shift-su[None,:,:]
        ds -= np.rint(ds)
        r = np.linalg.norm(ds @ cell,axis=2)
        j[a:a+step] = np.argmin(r,axis=1)
        d[a:a+step] = r[np.arange(len(r)),j[a:a+step]]
    return j, d


def match_atoms(pos,X,cell):#The atom of the snapshot (positions X, cell) at each position of pos (A) : the same atoms in any order, wrapped or not,
    #possibly shifted as a whole (origin of a cube file). Returns the atoms and the largest distance, or None
    inv = np.linalg.inv(cell)
    sa, su = pos @ inv, X @ inv
    shifts = [np.zeros(3)]
    probe = np.linspace(0,len(sa)-1,min(8,len(sa))).astype(int)
    if nearest_atoms(sa[probe],su,cell,shifts[0])[1].max() >= MATCHTOL:#not at the same origin : the shifts that put the first atom on an atom
        shifts = [s for s in su-sa[0] if nearest_atoms(sa[probe],su,cell,s)[1].max() < MATCHTOL][:1]
    for shift in shifts:
        j, d = nearest_atoms(sa,su,cell,shift)
        if d.max() < MATCHTOL and len(np.unique(j)) == len(su):
            return j, float(d.max())
    return None


def bader_volumes(files,UMDFile,natom):#Volume of every atom (A^3) in each ACF.dat file, matched to the snapshot of the UMD file named in its path
    #(counted from 0, as umd2out.py names its POSCAR files) by the positions of the atoms : [(file, snapshot, unit of the file, largest distance,
    #atom volumes in the order of the UMD file, (sum of the volumes + vacuum) / volume of the cell)]
    offsets, cells = umdpf.index_snapshots(UMDFile,natom)
    col0, ncol = umdpf.COLUMNS['xcart']
    out = []
    for name in files:
        pos, vol, vacuum = read_acf(name)
        if len(pos) != natom:
            print('ERROR : the ACF file',name,'has',len(pos),'atoms, the UMD file',natom)
            sys.exit()
        tried, found = [], None
        for k in snapshot_numbers(name):
            if k >= len(offsets):
                continue
            tried.append(k)
            X = umdpf.read_atoms_at(UMDFile,offsets[k],natom,1,1,col0,ncol)[0]
            for unit,scale in (('A',1.0),('bohr',BOHR)):
                m = match_atoms(pos*scale,X,cells[k])
                if m is not None:
                    found = (k,unit,scale,m)
                    break
            if found is not None:
                break
        if found is None:
            print('ERROR : the atoms of',name,'are at the positions of none of the snapshots',tried,'named in its path (snapshots of',UMDFile,
                  'counted from 0, as umd2out.py -r poscar names them ; positions in A or bohr)')
            sys.exit()
        k, unit, scale, (j,dev) = found
        atomvol = np.zeros(natom)
        atomvol[j] = vol*scale**3
        out.append((name,k,unit,dev,atomvol,(vol.sum()+vacuum)*scale**3/abs(np.linalg.det(cells[k]))))
    return out


#*** Motions of the units : translation, rotation in a body frame, internal motions

def kabsch(P,Q,w):#Proper rotations R (n,3,3) with Q_a ~ R P_a for the centred sets P, Q (n, k, 3) with weights w (k), and the weighted rms deviations
    H = np.einsum('a,nai,naj->nij',w,P,Q)
    U, S, Wt = np.linalg.svd(H)
    D = np.ones((len(P),3))
    D[:,2] = np.sign(np.linalg.det(U)*np.linalg.det(Wt))
    R = np.einsum('nji,nj,nkj->nik',Wt,D,U)
    dev = Q-np.einsum('nij,naj->nai',R,P)
    return R, np.sqrt(np.einsum('a,nai->n',w,dev**2)/w.sum())


class Motions:#Velocities of the units split into the translation of the centre of mass, the rotation and the internal motions, at every snapshot
    #(signal : sqrt(M) V_cm, sqrt(I_k) omega_k, sqrt(m_a) u_a, so that the sum of the squares is twice the kinetic energy). The body frame of a unit
    #is followed by the rotation that best maps its shape at the start of its life (in its principal frame) onto its current shape ; when its
    #neighbours change (a proton hops) the new ones take the places they are nearest to. Linear units : angular velocity in the frame of the cell
    def __init__(self,N,msite,sitetype,nsnap,samples):
        self.m = msite
        self.M = msite.sum()
        self.k = len(msite)
        self.sitetype = sitetype
        self.atoms = -np.ones((N,self.k-1),dtype=np.int64)
        self.life = -np.ones(N,dtype=np.int64)
        self.ref = np.zeros((N,self.k,3))
        self.R = np.tile(np.eye(3),(N,1,1))
        self.signal = np.zeros((nsnap,N,6+3*self.k),dtype=np.float32)
        self.linear = None
        self.moments = np.zeros(3)
        self.count = 0
        self.sampletimes = samples
        self.samples = []

    def add(self,t,units,lives,centre,nb,rel,V):#units (n) of the snapshot t, their lives, central atoms, neighbours (n, k-1 : by element, then by
        #index), bond vectors (n, k-1, 3) ; V the velocities of the snapshot
        n = len(units)
        new = self.life[units] != lives
        atoms = self.atoms[units]
        same = ~new & (np.sort(atoms,axis=1) == np.sort(nb,axis=1)).all(axis=1)
        idx = np.take_along_axis(np.argsort(nb,axis=1),np.argsort(np.argsort(atoms,axis=1),axis=1),axis=1)
        sites = np.where(same[:,None],idx,np.arange(self.k-1)[None,:])#where each site is in nb
        atoms = np.where(same[:,None],atoms,nb)
        P = np.concatenate([np.zeros((n,1,3)),np.take_along_axis(rel,sites[:,:,None],axis=1)],axis=1)
        rho = P-(self.m[None,:,None]*P).sum(axis=1,keepdims=True)/self.M
        for q in np.flatnonzero(~new & ~same):#the neighbours changed during the life : each new one takes the nearest place of the rotated shape
            pred = self.ref[units[q]] @ self.R[units[q]].T
            for e in np.unique(self.sitetype):
                s = np.flatnonzero(self.sitetype == e)
                row, col = linear_sum_assignment(np.linalg.norm(pred[1+s][:,None,:]-rho[q,1+s][None,:,:],axis=2))
                atoms[q,s[row]] = nb[q,s[col]]
                rho[q,1+s[row]] = rho[q,1+s[col]].copy()
        self.atoms[units] = atoms
        self.life[units] = lives
        vel = V[np.concatenate([centre[:,None],atoms],axis=1)]
        vcm = (self.m[None,:,None]*vel).sum(axis=1)/self.M
        dv = vel-vcm[:,None,:]
        L = (self.m[None,:,None]*np.cross(rho,dv)).sum(axis=1)
        I = (self.m[None,:,None,None]*((rho**2).sum(axis=2)[:,:,None,None]*np.eye(3)-rho[:,:,:,None]*rho[:,:,None,:])).sum(axis=1)
        eig, vec = np.linalg.eigh(I)
        if self.linear is None:
            self.linear = bool(np.median(eig[:,0]/eig[:,2]) < LINEAR)
        omega = np.einsum('nij,nj->ni',np.linalg.pinv(I,rcond=1.0E-6),L)
        u = dv-np.cross(omega[:,None,:],rho)
        if new.any():#shape of the unit at the start of a life, in its principal frame (right-handed)
            E = vec[new]
            E[:,:,2] *= np.sign(np.linalg.det(E))[:,None]
            self.ref[units[new]] = np.einsum('nai,nij->naj',rho[new],E)
        R, dev = kabsch(self.ref[units],rho,self.m)
        self.R[units] = R
        if self.linear:
            rot = np.sqrt(0.5*np.trace(I,axis1=1,axis2=2))[:,None]*omega
        else:
            Ib = np.einsum('nji,njk,nki->ni',R,I,R)
            rot = np.sqrt(np.maximum(Ib,0.0))*np.einsum('nji,nj->ni',R,omega)
        self.signal[t,units] = np.hstack([np.sqrt(self.M)*vcm,rot,(np.sqrt(self.m)[None,:,None]*u).reshape(n,-1)])
        self.moments += eig.sum(axis=0)
        self.count += n
        if t in self.sampletimes:
            self.samples.append(rho[:16])


def symmetry_number(samples,m,sitetype):#Number of permutations of the neighbours (by element) that a proper rotation maps onto the unit, in most of
    #the sampled units : the rotational symmetry number (H2O 2, NH3 3, CH4 12, OH 1, CO2 2). None when there are too many permutations
    rho = np.concatenate(samples)
    blocks = [np.flatnonzero(sitetype == e) for e in np.unique(sitetype)]
    if np.prod([factorial(len(b)) for b in blocks]) > MAXPERM:
        return None
    bond = np.linalg.norm(rho[:,1:]-rho[:,:1],axis=2).mean()
    count = 0
    for combo in itertools.product(*[itertools.permutations(b) for b in blocks]):
        perm = np.arange(len(m))
        for b,p in zip(blocks,combo):
            perm[1+b] = 1+np.array(p)
        if np.median(kabsch(rho,rho[:,perm],m)[1]) < SYMTOL*bond:
            count += 1
    return count


def rigid_rotor(moments,T,sigma,linear):#Entropy of the classical rigid rotor (kB) : moments of inertia (amu A^2, increasing) ; symmetry number
    theta = vdos_thermo.H_SI**2/(8.0*np.pi**2*np.asarray(moments)*vdos_thermo.AMU*1.0E-20*vdos_thermo.KB_SI)#rotational temperatures (K)
    if linear:
        return 1.0+np.log(T/(sigma*np.sqrt(theta[1]*theta[2])))
    return np.log(np.sqrt(np.pi)*np.exp(1.5)/sigma*np.sqrt(T**3/np.prod(theta)))


def two_pt(nu,g,M,volume,T,SR,dof):#Entropies (kB per unit) of translation, rotation, internal motions from their DoS g (3, frequencies ; states per
    #cm^-1 per unit) : unit of mass M (amu) and volume (A^3) at T ; SR the entropy of the rigid rotor, dof its degrees of freedom (3, or 2 if linear)
    w = vdos_thermo.weights(nu,T)['S']
    out = {}
    gas, f, delta, y, shs = vdos_thermo.two_phase(nu,g[0],1,M,volume,T)
    out.update({'f_trn':f,'Delta_trn':delta,'y':y,'S_hs':shs,'gas_trn':gas,
                'S_trn':vdos_thermo.integrate(nu,gas,np.full_like(nu,shs/3.0),shs/3.0)+vdos_thermo.integrate(nu,g[0]-gas,*w)})
    gas, f, delta, _, _ = vdos_thermo.two_phase(nu,g[1]*3.0/dof,1,M,volume,T)#the same model for the rotations, normalized to 3 modes
    gas *= dof/3.0
    out.update({'f_rot':f,'Delta_rot':delta,'gas_rot':gas,
                'S_rot':vdos_thermo.integrate(nu,gas,np.full_like(nu,SR/dof),SR/dof)+vdos_thermo.integrate(nu,g[1]-gas,*w)})
    out['S_int'] = vdos_thermo.integrate(nu,g[2],*w)
    out['S_total'] = out['S_trn']+out['S_rot']+out['S_int']
    return out


def main(argv):
    umdpf.headerumd()
    start = time.time()
    UMDFile = BondFile = central = ''
    adjacents = []
    numbers = None
    mode = 'auto'
    maxgap = GAP
    minlife = 0.0
    Nsteps = 1
    nblocks = 4
    window = None
    dr = 0.05
    rshell = None
    nmap = 30
    pairsampling = None
    temperature = None
    volumespec = None
    sigma = None
    usage = 'rotations_umd.py -f <UMD file> -b <bonding file> -c <central element> -a <adjacent elements> -x <numbers of adjacent atoms> -u <dipole|bonds> -g <gap> -m <min life> -s <sampling> -n <blocks> -w <window> -d <dr> -r <first shell> -p <map bins> -j <pair sampling> -V <volume of a unit> -t <temperature> -y <symmetry number>'
    try:
        opts, arg = getopt.getopt(argv,"hf:b:c:a:x:u:g:m:s:n:w:d:r:p:j:V:t:y:")
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('rotations_umd.py : orientations of the molecular units of a run, their correlations between units and in time')
            print(usage)
            print('-f : UMD file ; -f all : every *.umd.dat file of the folder, each with its bonding file')
            print('-b : bonding file of bonding_umd.py (default <base>.bonding.dat)')
            print('-c : central element of the units (one) ; -a : adjacent elements, separated by commas. A unit is a central atom with its bonded')
            print('     adjacent atoms : -c O -a H water, -c N -a H ammonia, -c C -a H methane, -c S -a O sulfate.')
            print('-x : number of adjacent atoms of each element of -a, e.g. -x 2 for H2O, -x 3 for H3O or NH3, -x 3,1 with -a H,O')
            print('     (default : the most frequent composition of the central atoms that have neighbours). A snapshot of a central atom counts')
            print('     only when it has exactly this composition : CH3 is not CH4.')
            print('-u : orientation of a unit. dipole : the axis from the central atom along the sum of its bond vectors (the dipole of H2O, NH3,')
            print('     OH by symmetry ; geometric, the UMD file has no electronic dipoles). bonds : its bond directions, for units without a')
            print('     dipole (CH4, NH4, SO4, CO3). Default : dipole when the units have one, else bonds.')
            print('     The correlation of order l between two units is the mean of P_l(a.b) over their vectors a, b : for dipoles l = 1 is the')
            print('     cosine of their angle and l = 2 is P2 ; for bonds the orders 1 to 4 that the shape of the unit allows (a tetrahedron : 3 and 4).')
            print('-g : longest gap, in fs, during which a unit that lost its composition (a bond across the cutoff, a proton leaving and coming')
            print('     back) goes on as the same molecule ; it is left out during the gap. A longer gap ends it. Default',GAP)
            print('-m : shortest life of a unit, in fs ; shorter lives are left out of everything. Default 0')
            print('-s : every s-th snapshot of the bonding file. Default 1')
            print('-n : number of blocks of the run for the errors. Default 4')
            print('-w : length of the time correlations, in snapshots. Default half a block, at most half the run')
            print('-d : bin of the distance between units, in A. Default 0.05 ; distances up to half the smallest width of the cell')
            print('-r : radius of the first shell of units, in A. Default : the first minimum of the g(r) of the central atoms of the units')
            print('-p : bins of the map of the orientations along cos(theta) (twice as many along phi). Default 30')
            print('-j : the pair correlations use every j-th snapshot. Default : about',PAIRSNAPS,'snapshots over the run')
            print('Entropy (2PT), when the UMD file has velocities : translation, rotation and internal motions of the units, from their densities of')
            print('     states ; translation and rotation split into gas-like (hard spheres, rigid rotor) and solid-like (quantum oscillators) parts.')
            print('-V : volume of a unit, for its density in the 2PT model, in one of three forms :')
            print('     a number, in A^3 : -V 25.3')
            print('     Bader volumes per element, in A^3, summed over the atoms of a unit : -V O:16.5,H:4.2')
            print('     ACF.dat files of the Bader program, one per snapshot analysed : files, patterns in quotes, or folders (their files with ACF in')
            print('        the name, also in sub-folders), separated by commas : -V bader/ or -V \'ACF_*.dat\'. The number of the snapshot is in the')
            print('        name of the file or of its folder (counted from 0, as umd2out.py -r poscar names its files : run.umd.dat._1200_.POSCAR.vasp')
            print('        -> ACF_1200.dat or 1200/ACF.dat) ; the atoms of the file are matched to those of the snapshot by their positions (any')
            print('        order, A or bohr, any origin), so the Bader analysis can come from any code. Volume of a unit = sum of the volumes of its')
            print('        atoms in the snapshot, averaged over the units and the files.')
            print('     Default : volume of the cell / number of units (right only when the units fill the cell, e.g. pure water)')
            print('-t : temperature, K. Default : the mean Temperature of the snapshots used')
            print('-y : rotational symmetry number of the unit (H2O 2, NH3 3, CH4 12). Default : found from the shapes of the units')
            print('Outputs, <base>.<unit> with the unit named from the central atom and its neighbours (OH2, NH3, CH4) :')
            print('  .rot.pair.dat  : against the distance r of the central atoms : g(r), number of pairs, correlations between the two units')
            print('                   (cos and P2 of the angle of their dipoles, and dipolar = 3 (u_i.r)(u_j.r) - u_i.u_j ; K1 .. K4 for bonds),')
            print('                   Kirkwood factor of the units within r')
            print('  .rot.shell.dat : distribution of the cosine of the angle between the dipoles of first neighbours')
            print('  .rot.time.dat  : in time : self correlation of the orientation of a unit (C1, C2 ; bonds : the orders allowed), of its bonds')
            print('                   followed atom by atom (bond_C1, bond_C2 : they also see the moves that map the unit onto itself),')
            print('                   and of the total orientation of all the units (total_C1 ... : collective)')
            print('  .rot.map.dat   : density of the orientations on the sphere, in the Cartesian frame, relative to a uniform one')
            print('  .rot.2pt.dat   : densities of states of translation, rotation, internal motions of a unit, and the gas-like parts (2PT)')
            print('  .rot.dat       : summary : units and lives, first shell, Kirkwood factor, polarization, correlation times, entropy of the map,')
            print('                   entropies of translation, rotation and internal motions (2PT), fluidicities, volume, moments of inertia')
            sys.exit()
        elif opt == '-f':
            UMDFile = str(arg)
        elif opt == '-b':
            BondFile = str(arg)
        elif opt == '-c':
            central = str(arg)
        elif opt == '-a':
            adjacents = arg.split(',')
        elif opt == '-x':
            numbers = [int(x) for x in arg.split(',')]
        elif opt == '-u':
            mode = str(arg)
        elif opt == '-g':
            maxgap = float(arg)
        elif opt == '-m':
            minlife = float(arg)
        elif opt == '-s':
            Nsteps = int(arg)
        elif opt == '-n':
            nblocks = int(arg)
        elif opt == '-w':
            window = int(arg)
        elif opt == '-d':
            dr = float(arg)
        elif opt == '-r':
            rshell = float(arg)
        elif opt == '-p':
            nmap = int(arg)
        elif opt == '-j':
            pairsampling = int(arg)
        elif opt == '-V':
            volumespec = str(arg)
        elif opt == '-t':
            temperature = float(arg)
        elif opt == '-y':
            sigma = int(arg)
    if not os.path.isfile(UMDFile):
        umdpf.missing_umd(UMDFile)
        sys.exit()
    base = UMDFile[:-8] if UMDFile.endswith('.umd.dat') else os.path.splitext(UMDFile)[0]
    if not BondFile:
        BondFile = base+'.bonding.dat'
    if not os.path.isfile(BondFile):
        print('the bonding file',BondFile,'does not exist : compute it with bonding_umd.py, or give it with -b')
        sys.exit()
    if not central or ',' in central or not adjacents:
        print('give the central element (-c, one element) and the adjacent elements (-a) of the units, e.g. -c O -a H')
        sys.exit()
    if mode not in ('auto','dipole','bonds') or nblocks < 1 or Nsteps < 1 or dr <= 0 or nmap < 2 or maxgap < 0:
        print('-u is dipole or bonds ; -n, -s at least 1 ; -d above 0 ; -p at least 2 ; -g at least 0')
        sys.exit()
    if numbers is not None and len(numbers) != len(adjacents):
        print('-x takes one number per element of -a :',len(adjacents))
        sys.exit()
    volume = volume_option(volumespec) if volumespec else None

    #The bonds central - adjacent of every snapshot used, and the composition of every central atom
    CentIndexes,AdjIndexes,MyCrystal,[Bonds,BondsIndexes],TimeStep = umdpf.read_bonds(BondFile,[central],adjacents,Nsteps)
    if CentIndexes == -1:#an element missing from the run : read_bonds printed the error
        sys.exit()
    natom = MyCrystal.natom
    typat = np.asarray(MyCrystal.typat)
    iscentral = typat == MyCrystal.elements.index(central)
    adjtype = -np.ones(natom,dtype=np.int64)
    for k,el in enumerate(adjacents):
        adjtype[typat == MyCrystal.elements.index(el)] = k
    centres = np.flatnonzero(iscentral)
    N = len(centres)
    cid = -np.ones(natom,dtype=np.int64)
    cid[centres] = np.arange(N)
    MyCrystalUMD,TimeStepUMD = umdpf.Crystallization(UMDFile)
    ratio = max(1,int(round(TimeStep/TimeStepUMD))) if TimeStepUMD > 0 else Nsteps
    MyCrystalUMD,X,TimeStepUMD,length = umdpf.read_values(UMDFile,"xcart","lists",ratio,MyCrystal.firstsnapshot)
    nsnap = min(len(Bonds),len(X))
    if nsnap < len(Bonds):
        print('WARNING : the UMD file has fewer snapshots than the bonding file : only the first',nsnap,'are used')
    if nsnap < 2*nblocks:
        print('Only',nsnap,'snapshots : too few for',nblocks,'blocks')
        sys.exit()
    cells = umdpf.read_cells(UMDFile,nsnap,ratio,MyCrystal.firstsnapshot)
    Vel = umdpf.read_values(UMDFile,"velocity","lists",ratio,MyCrystal.firstsnapshot)[1][:nsnap]
    v2 = (Vel**2).sum(axis=(1,2))                                       #snapshots with velocities : the first one of a parser may have none, or,
    velok = (v2 > 0.1*np.median(v2)) & (v2 < 10.0*np.median(v2))        #written by former parsers, the positions in their place
    if velok.any() and not velok.all():
        print('Velocities of',int((~velok).sum()),'snapshots zero or more than 10 times off the others (the first snapshot of a parser) : left out of the entropy')

    def bonds_of(t):#The bonds central atom -> adjacent atom of the snapshot t
        offsets = np.asarray(BondsIndexes[t][:natom+1])
        at = np.repeat(np.arange(natom),np.diff(offsets))
        ne = np.asarray(Bonds[t][offsets[0]:offsets[-1]],dtype=np.int64)
        keep = iscentral[at] & (adjtype[ne] >= 0) & (at != ne)
        return at[keep], ne[keep]

    comp = np.zeros((nsnap,N,len(adjacents)),dtype=np.int16)
    for t in range(nsnap):
        at, ne = bonds_of(t)
        np.add.at(comp[t],(cid[at],adjtype[ne]),1)
    kinds, counts = np.unique(comp.reshape(-1,len(adjacents)),axis=0,return_counts=True)
    order = np.argsort(-counts)
    def label(c):
        return central+''.join(el+(str(n) if n != 1 else '') for el,n in zip(adjacents,c) if n > 0)
    print('Compositions of the',N,central,'atoms over',nsnap,'snapshots :',
          ', '.join('%s %.2f %%' % (label(kinds[k]),100.0*counts[k]/counts.sum()) for k in order[:6]))
    if numbers is None:
        bound = [k for k in order if kinds[k].sum() > 0]
        if not bound:
            print('No',central,'atom is bonded to',', '.join(adjacents),'in',BondFile)
            sys.exit()
        numbers = [int(x) for x in kinds[bound[0]]]
    species = label(numbers)
    nbonds = sum(numbers)
    if nbonds == 0:
        print('-x : a unit has at least one adjacent atom')
        sys.exit()
    valid = (comp == np.array(numbers)).all(axis=2)
    lives, lid = lives_of(valid,int(maxgap//TimeStep),max(1,int(np.ceil(minlife/TimeStep))))
    used = valid & (lid >= 0)
    if not used.any():
        print('No unit',species,'in the run (with lives of at least',minlife,'fs)')
        sys.exit()
    spans = np.array([e-s+1 for i,s,e in lives])
    gaps = int(sum(e-s+1 for i,s,e in lives)-used.sum())
    print('Units',species,': %.2f per snapshot (%d at most), %d lives of %.1f fs on average (%d over the whole run) ; gaps up to %g fs bridged : %d snapshots of units left out'
          % (used.sum()/nsnap,used.sum(axis=1).max(),len(lives),spans.mean()*TimeStep,int((spans == nsnap).sum()),maxgap,gaps))

    names = [base+'.'+species+ext for ext in ('.rot.pair.dat','.rot.shell.dat','.rot.time.dat','.rot.map.dat','.rot.dat','.rot.2pt.dat')]
    umdpf.remove_former_outputs(names)

    #Orientation of every unit at every snapshot : the bond vectors (nearest images, cell of the snapshot) ; dipole or bonds
    inverses = np.linalg.inv(cells)
    def vectors_of(t):#bonds of the units of the snapshot t : unit of each bond, neighbour atom, bond vector
        at, ne = bonds_of(t)
        keep = used[t,cid[at]]
        at, ne = at[keep], ne[keep]
        d = (X[t][ne]-X[t][at]) @ inverses[t]
        d -= np.rint(d)
        return cid[at], ne, d @ cells[t]
    if mode == 'auto':
        polar = []
        for t in range(0,nsnap,max(1,nsnap//50)):
            c, ne, v = vectors_of(t)
            if len(c):
                b = v/np.linalg.norm(v,axis=1)[:,None]
                total = np.zeros((N,3))
                np.add.at(total,c,b)
                polar.append(np.linalg.norm(total[used[t]],axis=1).mean()/nbonds)
        mode = 'dipole' if np.mean(polar) > POLAR else 'bonds'
        print('|sum of the bond directions| / number of bonds = %.3f : orientation of the units by their %s' % (np.mean(polar),'dipole axis' if mode == 'dipole' else 'bonds (no dipole)'))
    orders = [1,2] if mode == 'dipole' else list(range(1,LMAX+1))
    slices = []
    ncomp = 0
    for l in orders:
        slices.append(slice(ncomp,ncomp+2*l+1))
        ncomp += 2*l+1
    bslices = [slice(0,3),slice(3,8)]#the tagged bonds : orders 1 and 2

    widths = [abs(np.linalg.det(c))/np.linalg.norm(np.cross(c[[1,2,0]],c[[2,0,1]]),axis=1) for c in cells]
    rmax = 0.5*np.min(widths)
    nbins = int(rmax/dr)
    if pairsampling is None:
        pairsampling = max(1,nsnap//PAIRSNAPS)
    blockof = np.minimum(np.arange(nsnap)*nblocks//nsnap,nblocks-1)
    edges = np.arange(nbins+1)*dr
    shellvolume = 4.0/3.0*np.pi*(edges[1:]**3-edges[:-1]**3)
    paircount = np.zeros((nblocks,nbins))
    ideal = np.zeros((nblocks,nbins))
    pairsum = np.zeros((nblocks,len(orders),nbins))
    dipolar = np.zeros((nblocks,nbins))
    coshist = np.zeros((nblocks,nbins,NCOS))
    linkcount = np.zeros((nblocks,nbins))                               #linked pairs : a bond of one unit points to the other
    linksum = np.zeros((nblocks,len(orders),nbins))
    linkcos = np.zeros((nblocks,nbins,NCOS))
    place = -np.ones(N,dtype=np.int64)
    unitcount = np.zeros(nblocks)
    mapcount = np.zeros((nblocks,nmap*2*nmap))
    maphalf = np.zeros((2,nmap*2*nmap))                                 #the map of one unit out of two (even and odd units) : sampling of its entropy
    Q = np.zeros((nsnap,N,ncomp))
    U = np.zeros((nsnap,N,3)) if mode == 'dipole' else None
    btime, bunit, batom, bharm = [], [], [], []
    motions = None
    if velok.any():#translation, rotation and internal motions of the units, for their entropy
        masses = np.asarray(MyCrystalUMD.masses,dtype=float)
        sitetype = np.repeat(np.arange(len(adjacents)),numbers)
        msite = np.r_[masses[MyCrystalUMD.elements.index(central)],masses[[MyCrystalUMD.elements.index(adjacents[e]) for e in sitetype]]]
        motions = Motions(N,msite,sitetype,nsnap,set(range(0,nsnap,max(1,nsnap//50))))
    else:
        print('The UMD file has no velocities : no entropy of the units')
    for t in range(nsnap):
        c, ne, v = vectors_of(t)
        if len(c) == 0:
            continue
        b = v/np.linalg.norm(v,axis=1)[:,None]
        btime.append(np.full(len(c),t))
        bunit.append(c)
        batom.append(ne)
        bharm.append(np.hstack([harmonics(b,1),harmonics(b,2)]))
        if motions is not None and velok[t]:
            srt = np.lexsort((ne,adjtype[ne],c))
            cu = c[srt][::nbonds]
            motions.add(t,cu,lid[t,cu],centres[cu],ne[srt].reshape(-1,nbonds),v[srt].reshape(-1,nbonds,3),Vel[t])
        units = np.flatnonzero(used[t])
        if mode == 'dipole':
            total = np.zeros((N,3))
            np.add.at(total,c,v)
            u = total[units]/np.linalg.norm(total[units],axis=1)[:,None]
            U[t,units] = u
            Q[t,units] = np.hstack([harmonics(u,l) for l in orders])
            shown = u
        else:
            np.add.at(Q[t],c,np.hstack([harmonics(b,l) for l in orders])/nbonds)
            shown = b
        ib = blockof[t]
        cosbin = np.minimum((0.5*(shown[:,2]+1.0)*nmap).astype(int),nmap-1)
        phibin = np.minimum(((np.arctan2(shown[:,1],shown[:,0])+np.pi)/(2.0*np.pi)*2*nmap).astype(int),2*nmap-1)
        mapcount[ib] += np.bincount(cosbin*2*nmap+phibin,minlength=nmap*2*nmap)
        if mode == 'dipole':
            for half in (0,1):
                maphalf[half] += np.bincount((cosbin*2*nmap+phibin)[units % 2 == half],minlength=nmap*2*nmap)
        if t % pairsampling or len(units) < 2:
            continue
        #Pairs of units of this snapshot, each once : distance of their central atoms (nearest image), correlations of their orientations
        n = len(units)
        s = X[t][centres[units]] @ inverses[t]
        q = Q[t,units]
        unitcount[ib] += n
        ideal[ib] += 0.5*n*(n-1)*shellvolume/abs(np.linalg.det(cells[t]))
        chunk = max(1,int(MAXEL//(3*n)))
        for a in range(0,n,chunk):
            ds = s[a:a+chunk,None,:]-s[None,:,:]
            ds -= np.rint(ds)
            dv = ds @ cells[t]
            r = np.linalg.norm(dv,axis=2)
            kbin = (r/dr).astype(int)
            ok = (np.arange(a,min(a+chunk,n))[:,None] < np.arange(n)[None,:]) & (kbin < nbins)
            kb = kbin[ok]
            paircount[ib] += np.bincount(kb,minlength=nbins)
            for k,sl in enumerate(slices):
                value = q[a:a+chunk,sl] @ q[:,sl].T
                pairsum[ib,k] += np.bincount(kb,weights=value[ok],minlength=nbins)
                if mode == 'dipole' and k == 0:
                    cosine = value[ok]
                    np.add.at(coshist[ib],(kb,np.clip((0.5*(cosine+1.0)*NCOS).astype(int),0,NCOS-1)),1)
                    rhat = dv[ok]/r[ok][:,None]
                    ii, jj = np.nonzero(ok)
                    along = (u[a+ii]*rhat).sum(axis=1)*(u[jj]*rhat).sum(axis=1)
                    dipolar[ib] += np.bincount(kb,weights=3.0*along-cosine,minlength=nbins)
        #Linked pairs : a bond of the unit i points to the unit j : j is the nearest other unit of its adjacent atom, which is closer to the central
        #atom of j than the central atom of i is (the hydrogen bond, donor -> acceptor, of H2O, NH3, OH). Against the distance of the central atoms
        place[units] = np.arange(n)
        donor = place[c]
        sa = X[t][ne] @ inverses[t]
        for a in range(0,len(donor),chunk):
            i = donor[a:a+chunk]
            ds = sa[a:a+chunk,None,:]-s[None,:,:]
            ds -= np.rint(ds)
            dist = np.linalg.norm(ds @ cells[t],axis=2)
            dist[np.arange(len(i)),i] = np.inf
            j = np.argmin(dist,axis=1)
            dc = s[i]-s[j]
            dc -= np.rint(dc)
            rij = np.linalg.norm(dc @ cells[t],axis=1)
            kbin = (rij/dr).astype(int)
            ok = (dist[np.arange(len(i)),j] < rij) & (kbin < nbins)
            i, j, kb = i[ok], j[ok], kbin[ok]
            linkcount[ib] += np.bincount(kb,minlength=nbins)
            for k,sl in enumerate(slices):
                value = (q[i,sl]*q[j,sl]).sum(axis=1)
                linksum[ib,k] += np.bincount(kb,weights=value,minlength=nbins)
                if mode == 'dipole' and k == 0:
                    np.add.at(linkcos[ib],(kb,np.clip((0.5*(value+1.0)*NCOS).astype(int),0,NCOS-1)),1)
    del X, Vel

    #Self values of the orders (shape of the unit) ; for the bonds only the orders that the shape allows
    selfvalue = np.array([(Q[...,sl]**2).sum(axis=2)[used].mean() for sl in slices])
    if mode == 'bonds':
        active = [k for k in range(len(orders)) if selfvalue[k] > ACTIVE]
        print('mean of P_l between the bonds of a unit :',', '.join('l=%d %.3f' % (orders[k],selfvalue[k]) for k in range(len(orders))),
              '; orders used :',', '.join(str(orders[k]) for k in active))
        if not active:
            print('No order up to',LMAX,'describes the orientation of these units')
            sys.exit()
    else:
        active = [0,1]
    tag = ['cos','P2'] if mode == 'dipole' else ['K'+str(orders[k]) for k in active]

    #*** Pairs of units against their distance
    r = edges[:-1]+0.5*dr
    count = paircount.sum(axis=0)
    with np.errstate(invalid='ignore',divide='ignore'):
        g = np.where(ideal.sum(axis=0) > 0,count/ideal.sum(axis=0),0.0)
        K = np.array([np.where(count > 0,pairsum[:,k].sum(axis=0)/np.where(count > 0,count,1)/selfvalue[k],np.nan) for k in active])
        Kb = np.array([[np.where(paircount[b] > 0,pairsum[b,k]/np.where(paircount[b] > 0,paircount[b],1)/selfvalue[k],np.nan) for k in active] for b in range(nblocks)])
        Kerr = block_error(Kb)
        dip = np.where(count > 0,dipolar.sum(axis=0)/np.where(count > 0,count,1),np.nan)
    kirkwood = 1.0+np.cumsum(2.0*pairsum[:,active[0]].sum(axis=0)/selfvalue[active[0]])/unitcount.sum() if unitcount.sum() else np.full(nbins,np.nan)
    coordination = np.cumsum(2.0*count)/unitcount.sum() if unitcount.sum() else np.zeros(nbins)
    nused = max(1,len(range(0,nsnap,pairsampling)))
    labels = ['r_(A)','g','pairs_per_snapshot','linked_per_snapshot']
    columns = [r,g,count/nused,linkcount.sum(axis=0)/nused]
    for k,name in enumerate(tag):
        labels += [name,name+'_err']
        columns += [K[k],Kerr[k]]
    if mode == 'dipole':
        labels.append('dipolar')
        columns.append(dip)
    labels.append('Kirkwood')
    columns.append(kirkwood)
    what = ('cos, P2 : mean of u_i.u_j and of P2(u_i.u_j) for the dipole axes u of two units at this distance ; dipolar : mean of 3 (u_i.r)(u_j.r) - u_i.u_j, r the unit vector between them'
            if mode == 'dipole' else 'K_l : mean of P_l between the bonds of two units at this distance, divided by its value inside a unit ('+', '.join('l=%d %.4f' % (orders[k],selfvalue[k]) for k in active)+')')
    comments = ['rotations_umd.py, '+UMDFile+', units '+species+' ('+('dipole axis' if mode == 'dipole' else 'bonds')+') ; every '+str(pairsampling)+' snapshots',
                'r : distance between the central atoms of two units ; g : their pair distribution ; '+what,
                'linked : a bond of one unit points to the other (its adjacent atom is nearest to that unit, and closer to its central atom than its own central atom is : the hydrogen bond)',
                'Kirkwood : 1 + sum over the other units within r of '+tag[0]+' ; errors from '+str(nblocks)+' blocks']
    write_columns(names[0],comments,labels,columns)

    #First shell of units
    shell = None
    note = ''
    given = rshell is not None
    if rshell is None:
        values, note = analyze_gofr.first_shell(r,g,coordination,0.03) if count.sum() else (None,'no pairs')
        if values is not None:
            rshell = values['xmin']
    if rshell is not None:
        inside = edges[1:] <= rshell+1e-9
        if paircount[:,inside].sum(axis=1).min() > 0:
            perblock = np.array([[pairsum[b,k][inside].sum()/paircount[b][inside].sum()/selfvalue[k] for k in active] for b in range(nblocks)])
            shell = {'r':rshell, 'neighbours':2.0*count[inside].sum()/unitcount.sum(),
                     'values':[pairsum[:,k][:,inside].sum()/count[inside].sum()/selfvalue[k] for k in active], 'errors':block_error(perblock)}
            if linkcount[:,inside].sum(axis=1).min() > 0:
                perblock = np.array([[linksum[b,k][inside].sum()/linkcount[b][inside].sum()/selfvalue[k] for k in active] for b in range(nblocks)])
                shell.update({'links':linkcount[:,inside].sum()/unitcount.sum(),'linkerrors':block_error(perblock),
                              'linkvalues':[linksum[:,k][:,inside].sum()/linkcount[:,inside].sum()/selfvalue[k] for k in active]})
            if mode == 'dipole':
                h = coshist[:,inside].sum(axis=(0,1))
                hb = coshist[:,inside].sum(axis=1)
                centers = -1.0+(np.arange(NCOS)+0.5)*2.0/NCOS
                with np.errstate(invalid='ignore',divide='ignore'):
                    density = h/h.sum()*NCOS/2.0
                    derr = block_error(hb/hb.sum(axis=1)[:,None]*NCOS/2.0)
                    hl = linkcos[:,inside].sum(axis=(0,1))
                    hlb = linkcos[:,inside].sum(axis=1)
                    ldensity = hl/hl.sum()*NCOS/2.0 if hl.sum() else np.full(NCOS,np.nan)
                    lerr = block_error(hlb/hlb.sum(axis=1)[:,None]*NCOS/2.0) if hlb.sum(axis=1).min() > 0 else np.full(NCOS,np.nan)
                write_columns(names[1],['rotations_umd.py, '+UMDFile+', units '+species+' : cosine of the angle between the dipole axes of two units closer than %.3f A (first shell, %.2f neighbours)' % (rshell,shell['neighbours']),
                                        'probability density of the cosine (0.5 for independent orientations) : for all the pairs of the shell, and for the linked ones (a bond of one unit points to the other : the hydrogen bond)',
                                        'errors from '+str(nblocks)+' blocks'],
                              ['cos','density','density_err','pairs','density_linked','density_linked_err','linked'],[centers,density,derr,h,ldensity,lerr,hl])

    #*** Time correlations : self (units), tagged bonds, collective (sum over the units)
    blocklen = nsnap//nblocks
    nlags = min(nsnap//2,window if window is not None else max(2,blocklen//2))
    tlag = np.arange(nlags)*TimeStep
    btime, bunit, batom, bharm = np.concatenate(btime), np.concatenate(bunit), np.concatenate(batom), np.concatenate(bharm)
    blife = lid[btime,bunit]
    key = np.lexsort((btime,batom,blife))
    btime, batom, blife, bharm = btime[key], batom[key], blife[key], bharm[key]
    first = np.flatnonzero(np.r_[True,(blife[1:] != blife[:-1]) | (batom[1:] != batom[:-1])])
    last = np.r_[first[1:],len(blife)]
    rows = [slice(a,b) for a,b in zip(first,last)]
    bondsegments = [(g,lives[blife[a]][1],lives[blife[a]][2]) for g,a in enumerate(first)]
    M = Q.sum(axis=1)                                                    #total orientation of the units of each snapshot (zero where a unit is not used)
    def correlations(t0,t1,nl):
        num, den = masked_acf(unit_batches(clip_segments(lives,t0,t1),Q,used),nl,slices)
        selfc = normalized(num,den)
        num, den = masked_acf(bond_batches(clip_segments(bondsegments,t0,t1),rows,btime,bharm),nl,bslices)
        bondc = normalized(num,den)
        acf = acf_fft(M[t0:t1],nl)
        total = np.array([acf[:,sl].sum(axis=1) for sl in slices])
        return selfc, bondc, total/total[:,0:1], total[:,0]
    selfc, bondc, totalc, total0 = correlations(0,nsnap,nlags)
    nlb = min(nlags,blocklen//2)
    blocks = [correlations(b*blocklen,(b+1)*blocklen,nlb) for b in range(nblocks)]
    def errors(index,shape):
        err = np.full(shape,np.nan)
        if nlb > 1:
            err[:,:nlb] = block_error([blk[index] for blk in blocks])
        return err
    selferr, bonderr, totalerr = errors(0,selfc.shape), errors(1,bondc.shape), errors(2,totalc.shape)
    names_l = ['C'+str(orders[k]) for k in active]
    labels = ['time_(fs)']
    columns = [tlag]
    for k,name in zip(active,names_l):
        labels += [name,name+'_err']
        columns += [selfc[k],selferr[k]]
    for k,name in enumerate(('bond_C1','bond_C2')):
        labels += [name,name+'_err']
        columns += [bondc[k],bonderr[k]]
    for k,name in zip(active,names_l):
        labels += ['total_'+name,'total_'+name+'_err']
        columns += [totalc[k],totalerr[k]]
    write_columns(names[2],['rotations_umd.py, '+UMDFile+', units '+species+' ('+('dipole axis' if mode == 'dipole' else 'bonds')+') ; gaps up to '+str(maxgap)+' fs bridged',
                            'C_l : mean of P_l(u(0).u(t)) for the orientation of a unit, over its lives and all the time origins, 1 at t = 0'+('' if mode == 'dipole' else ' (over the bonds of the unit : it does not see a move that maps the unit onto itself)'),
                            'bond_C1, bond_C2 : the same for one bond, followed by its neighbour atom ; total_C_l : for the sum of the orientations of all the units (collective)',
                            'errors from '+str(nblocks)+' blocks, up to their length'],labels,columns)

    #*** Entropy of the units (2PT)
    lines2pt, rows2pt = [], []
    T = temperature
    if motions is not None and motions.count and T is None:
        temps = read_temperatures(UMDFile)
        sel = MyCrystal.firstsnapshot+ratio*np.arange(nsnap)
        sel = sel[sel < len(temps)]
        T = float(temps[sel].mean()) if len(sel) else 0.0
        if T <= 0.0:
            print('No Temperature in the UMD file : give it with -t for the entropy')
            motions = None
    if motions is not None and motions.count:
        k = motions.k
        dof = 2 if motions.linear else 3
        nint = 3*k-3-dof
        Munit = motions.M
        #Volume of a unit
        if volume is None:
            vols = np.array([abs(np.linalg.det(cells[t]))/used[t].sum() for t in range(nsnap) if used[t].any()])
            vunit, verr = float(vols.mean()), np.nan
            source = 'volume of the cell / number of units'
            filled = used.sum()*(1+nbonds)/(nsnap*natom)
            lines2pt.append('   volume of a unit : %.3f A^3 = volume of the cell / number of units' % vunit)
            if filled < 0.95:
                lines2pt.append('   WARNING : the atoms of the units are only %.1f %% of the atoms : this volume counts the other atoms too ; give the' % (100*filled))
                lines2pt.append('   volume of the units with -V (Bader volumes)')
        elif volume[0] == 'number':
            vunit, verr, source = volume[1], np.nan, 'given'
            lines2pt.append('   volume of a unit : %.3f A^3 (given)' % vunit)
        elif volume[0] == 'elements':
            missing = [el for el in [central]+adjacents if el not in volume[1]]
            if missing:
                print('-V : no volume for',', '.join(missing))
                sys.exit()
            vunit = volume[1][central]+sum(n*volume[1][el] for el,n in zip(adjacents,numbers))
            verr, source = np.nan, 'Bader volumes of the elements '+','.join('%s:%g' % kv for kv in volume[1].items())
            lines2pt.append('   volume of a unit : %.3f A^3, from the volumes of its atoms : %s' % (vunit,', '.join('%s %g' % kv for kv in volume[1].items())))
        else:
            unitvols, filemeans, perelement = [], [], {el:[] for el in [central]+adjacents}
            lines2pt.append('   volume of a unit from the Bader volumes of its atoms (ACF.dat files) :')
            for name,snap,unit,dev,atomvol,closure in bader_volumes(volume[1],UMDFile,natom):
                tb = (snap-MyCrystal.firstsnapshot)/ratio
                t = int(round(tb))
                if t < 0 or t >= nsnap:
                    lines2pt.append('      %s : snapshot %d is outside the snapshots analysed : left out' % (name,snap))
                    continue
                at, ne = bonds_of(t)
                keep = valid[t,cid[at]]
                uv = atomvol[centres].copy()
                np.add.at(uv,cid[at[keep]],atomvol[ne[keep]])
                units = np.flatnonzero(valid[t])
                unitvols += list(uv[units])
                filemeans.append(uv[units].mean() if len(units) else np.nan)
                for el in perelement:
                    perelement[el] += list(atomvol[typat == MyCrystal.elements.index(el)])
                lines2pt.append('      %s : snapshot %d (positions in %s, within %.4f A)%s, %d units, %.3f A^3 per unit ; sum of the volumes / cell %.4f'
                                % (name,snap,unit,dev,'' if tb == t else ' (units of the nearest snapshot of the bonding file)',len(units),filemeans[-1],closure))
                if abs(closure-1.0) > 0.02:
                    lines2pt.append('      WARNING : the volumes of this file do not fill the cell : the Bader analysis is not of this cell ?')
            if not unitvols:
                print('-V : no unit in the snapshots of the ACF files')
                sys.exit()
            vunit = float(np.mean(unitvols))
            fm = np.array(filemeans)[np.isfinite(filemeans)]
            verr = float(np.std(fm,ddof=1)/np.sqrt(len(fm))) if len(fm) > 1 else np.nan
            source = 'Bader volumes of '+str(len(fm))+' snapshots'
            lines2pt.append('      atoms : '+', '.join('%s %.3f +/- %.3f A^3' % (el,np.mean(v),np.std(v)) for el,v in perelement.items() if v))
            lines2pt.append('   volume of a unit : %.3f A^3, spread %.3f between units, error %.3f from %d snapshots' % (vunit,np.std(unitvols),verr,len(fm)))
        #Rigid rotor
        moments = motions.moments/motions.count
        found = symmetry_number(motions.samples,motions.m,motions.sitetype) if motions.samples else None
        sym = sigma if sigma is not None else (found if found else 1)
        if sigma is None and not found:
            lines2pt.append('   WARNING : symmetry number not found (too many permutations of the neighbours) : 1 ; give it with -y')
        SR = rigid_rotor(moments,T,sym,motions.linear)
        #Densities of states from the correlations of the motions over the lives, all origins ; blocks for the errors
        used2 = used & velok[:,None]
        slices2 = [slice(0,3),slice(3,6),slice(6,6+3*k)]
        def spectra(t0,t1,nl):
            num, den = masked_acf(unit_batches(clip_segments(lives,t0,t1),motions.signal,used2),nl,slices2)
            c = np.where(den > 0,num/np.where(den > 0,den,1.0),0.0)
            freq, S = vdos(c,TimeStep,4*nl,'hann')
            return freq, S*2.0*AMU_A2FS2/(KB*T)/CM1, c[:,0]
        nu, g, c0 = spectra(0,nsnap,nlags)
        res = two_pt(nu,g,Munit,vunit,T,SR,dof)
        keys = ('S_trn','S_rot','S_int','S_total','f_trn','f_rot')
        if nlb > 1:
            blockres = [two_pt(*spectra(b*blocklen,(b+1)*blocklen,nlb)[:2],Munit,vunit,T,SR,dof) for b in range(nblocks)]
            err = {key:float(block_error([[r[key]] for r in blockres])[0]) for key in keys}
        else:
            err = {key:np.nan for key in keys}
        larger = two_pt(nu,g,Munit,1.1*vunit,T,SR,dof)
        tkin = [c0[0]*AMU_A2FS2/(3.0*KB),c0[1]*AMU_A2FS2/(dof*KB),c0[2]*AMU_A2FS2/(nint*KB) if nint > 0 else np.nan]
        integral = [vdos_thermo.integrate(nu,g[q],np.ones_like(nu),1.0) for q in range(3)]
        write_columns(names[5],['rotations_umd.py, '+UMDFile+', units '+species+' : densities of states of a unit (states per cm^-1), T = %.2f K, volume %.3f A^3 (%s)' % (T,vunit,source),
                                'translation of the centre of mass, rotation (angular velocities in the body frame x moments of inertia), internal motions ; window %d snapshots, Hann taper' % nlags,
                                'gas_trn, gas_rot : their gas-like parts in the 2PT model (hard spheres, rigid rotor) ; the rest is solid-like (quantum oscillators)'],
                      ['frequency_(cm^-1)','translation','gas_trn','rotation','gas_rot','internal'],[nu,g[0],res['gas_trn'],g[1],res['gas_rot'],g[2]])
        J = vdos_thermo.R
        lines2pt[:0] = ['Entropy of a unit '+species+' (2PT) at %.2f K : %s, symmetry number %d%s' % (T,'linear' if motions.linear else 'non-linear',sym,'' if sigma is None else ' (given)')]
        lines2pt += ['   kinetic temperatures : translation %.1f K, rotation %.1f K, internal %.1f K (modes in the DoS : %.3f, %.3f, %.3f of %d, %d, %d)'
                     % (*tkin,*integral,3,dof,nint),
                     '   moments of inertia %.4f %.4f %.4f amu A^2 ; rigid rotor : %.4f kB' % (*moments,SR),
                     '   translation : Delta %.4f, fluidicity %.4f +/- %.4f, packing fraction %.4f' % (res['Delta_trn'],res['f_trn'],err['f_trn'],res['y']),
                     '   rotation    : Delta %.4f, fluidicity %.4f +/- %.4f' % (res['Delta_rot'],res['f_rot'],err['f_rot'])]
        for key,name in (('S_trn','translation'),('S_rot','rotation'),('S_int','internal'),('S_total','total')):
            lines2pt.append('   S %-11s : %8.4f +/- %.4f kB = %8.3f +/- %.3f J/(mol K)' % (name,res[key],err[key],res[key]*J,err[key]*J))
            rows2pt.append(('S_'+name+'_(kB)',res[key],err[key]))
        lines2pt.append('   10 %% more volume : S translation %+.4f kB, rotation %+.4f kB' % (larger['S_trn']-res['S_trn'],larger['S_rot']-res['S_rot']))
        lines2pt.append('   rotation of the free rigid rotor %.4f kB : the units lose %.4f kB of rotational entropy' % (SR,SR-res['S_rot']))
        if nint > 0 and TimeStep > 1.0/(2.0*4000.0)*CM1:
            lines2pt.append('   WARNING : %.2f fs between snapshots : frequencies above %.0f cm^-1 are folded back ; use every snapshot of the run' % (TimeStep,CM1/(2.0*TimeStep)))
        rows2pt = [('temperature_(K)',T,np.nan),('volume_per_unit_(A3)',vunit,verr),('symmetry_number',sym,np.nan),
                   ('T_kin_translation_(K)',tkin[0],np.nan),('T_kin_rotation_(K)',tkin[1],np.nan),('T_kin_internal_(K)',tkin[2],np.nan),
                   ('I_1_(amu.A2)',moments[0],np.nan),('I_2_(amu.A2)',moments[1],np.nan),('I_3_(amu.A2)',moments[2],np.nan),('S_rigid_rotor_(kB)',SR,np.nan),
                   ('Delta_translation',res['Delta_trn'],np.nan),('fluidicity_translation',res['f_trn'],err['f_trn']),('packing_fraction',res['y'],np.nan),
                   ('Delta_rotation',res['Delta_rot'],np.nan),('fluidicity_rotation',res['f_rot'],err['f_rot'])]+rows2pt
        rows2pt += [('S_total_(J/mol/K)',res['S_total']*J,err['S_total']*J),('dS_translation_10%_volume_(kB)',larger['S_trn']-res['S_trn'],np.nan),
                    ('dS_rotation_10%_volume_(kB)',larger['S_rot']-res['S_rot'],np.nan)]

    #*** Kirkwood factors of the whole cell, polarization, map and its entropy
    collective = [float(total0[k]/((Q[...,slices[k]]**2).sum(axis=2).sum(axis=1).mean())) for k in active]
    collective_b = [[float((M[b*blocklen:(b+1)*blocklen,slices[k]]**2).sum(axis=1).mean()/(Q[b*blocklen:(b+1)*blocklen,:,slices[k]]**2).sum(axis=(1,2)).mean()) for k in active] for b in range(nblocks)]
    collective_err = block_error(collective_b)
    density = mapcount.sum(axis=0)
    cosc = -1.0+(np.arange(nmap)+0.5)*2.0/nmap
    phic = -180.0+(np.arange(2*nmap)+0.5)*180.0/nmap
    write_columns(names[3],['rotations_umd.py, '+UMDFile+', units '+species+' : directions of the '+('dipole axes' if mode == 'dipole' else 'bonds')+' in the Cartesian frame of the UMD file',
                            'cos_theta : component along z ; phi : angle from x in the xy plane (degrees) ; bins of equal solid angle ; density relative to a uniform distribution (1)'],
                  ['cos_theta','phi_(deg)','density','counts'],
                  [np.repeat(cosc,2*nmap),np.tile(phic,nmap),density/density.sum()*density.size,density])
    entropy = map_entropy(density)
    entropy_err = float(block_error([[map_entropy(mapcount[b])] for b in range(nblocks)])[0])
    entropy_half = float(np.nanmean([map_entropy(maphalf[0]),map_entropy(maphalf[1])])) if mode == 'dipole' and maphalf.sum(axis=1).min() > 0 else np.nan

    #*** Summary
    rowsout = [('units_per_snapshot',used.sum()/nsnap,np.nan),('lives',len(lives),np.nan),('mean_life_(fs)',spans.mean()*TimeStep,np.nan)]
    print('\nUnits',species,'of',UMDFile,':',('dipole axes' if mode == 'dipole' else 'bonds'))
    if shell is not None:
        print('First shell of units : up to %.3f A (%s), %.2f neighbours' % (shell['r'],'given' if given else 'first minimum of g(r)',shell['neighbours']))
        rowsout += [('first_shell_(A)',shell['r'],np.nan),('neighbours',shell['neighbours'],np.nan)]
        for name,value,err in zip(tag,shell['values'],shell['errors']):
            print('   %-8s between first neighbours : %8.4f +/- %.4f' % (name,value,err))
            rowsout.append((name+'_first_shell',value,err))
        if 'links' in shell:
            print('   linked neighbours (a bond of the unit points to them ; for H2O, NH3, OH : the hydrogen bond it gives) : %.3f per unit' % shell['links'])
            rowsout.append(('links_per_unit',shell['links'],np.nan))
            for name,value,err in zip(tag,shell['linkvalues'],shell['linkerrors']):
                print('   %-8s between linked neighbours : %8.4f +/- %.4f' % (name,value,err))
                rowsout.append((name+'_linked',value,err))
    else:
        print('First shell of units not found ('+(note or 'no pairs in it')+') : give its radius with -r')
    for k,name,value,err in zip(active,tag,collective,collective_err):
        print('   Kirkwood factor of the cell, order %d : %8.4f +/- %.4f   (|sum over the units|^2 / sum of the units ; 1 without correlations)' % (orders[k],value,err))
        rowsout.append(('Kirkwood_'+str(orders[k]),value,err))
    if mode == 'dipole':
        pol = U.sum(axis=1).mean(axis=0)/(used.sum()/nsnap)
        polb = np.array([U[b*blocklen:(b+1)*blocklen].sum(axis=1).mean(axis=0)/max(1e-300,used[b*blocklen:(b+1)*blocklen].sum()/blocklen) for b in range(nblocks)])
        perr = block_error(polb)
        print('   mean dipole axis of a unit (polarization) : %.4f %.4f %.4f  +/- %.4f %.4f %.4f ; length %.4f' % (*pol,*perr,np.linalg.norm(pol)))
        for name,value,err in zip('xyz',pol,perr):
            rowsout.append(('polarization_'+name,value,err))
    for k,name in zip(active,names_l):
        tau, plateau = decay_time(tlag,selfc[k])
        print('   %s : falls to 1/e %s ; %.4f at the end of the window (%.0f fs)' % (name,('in %.1f fs' % tau) if tau is not None else 'not within the window',plateau,tlag[-1]))
        rowsout += [('tau_'+name+'_(fs)',tau if tau is not None else np.nan,np.nan),('end_'+name,plateau,np.nan)]
    for k,name in enumerate(('bond_C1','bond_C2')):
        tau, plateau = decay_time(tlag,bondc[k])
        print('   %s : falls to 1/e %s ; %.4f at the end of the window' % (name,('in %.1f fs' % tau) if tau is not None else 'not within the window',plateau))
        rowsout += [('tau_'+name+'_(fs)',tau if tau is not None else np.nan,np.nan),('end_'+name,plateau,np.nan)]
    for k,name in zip(active,names_l):
        tau, plateau = decay_time(tlag,totalc[k])
        rowsout += [('tau_total_'+name+'_(fs)',tau if tau is not None else np.nan,np.nan),('end_total_'+name,plateau,np.nan)]
    if mode == 'dipole':
        print('   entropy of the map of the dipole axes : %.4f +/- %.4f kB per unit, relative to free rotation (0) ; %d bins of %.4f sr' % (entropy,entropy_err,density.size,4.0*np.pi/density.size))
        print('      (one unit at a time : the correlations between units, e.g. the ice rules, lower it further)')
        rowsout.append(('map_entropy_(kB)',entropy,entropy_err))
        if np.isfinite(entropy_half):#fewer orientations sampled give a lower entropy, about as 1/number : with twice the units, about entropy + (entropy - half)
            print('      with one unit out of two : %.4f ; extrapolated to many units : %.4f' % (entropy_half,2.0*entropy-entropy_half))
            rowsout += [('map_entropy_half_units_(kB)',entropy_half,np.nan),('map_entropy_extrapolated_(kB)',2.0*entropy-entropy_half,np.nan)]
            if entropy-entropy_half > 0.05:
                print('      WARNING : the map is not sampled enough (the units do not turn during the run, or they are too few for its',density.size,'bins) :')
                print('      its entropy is too low ; the extrapolated value is a rough correction. Fewer bins (-p) need fewer orientations.')
    for line in lines2pt:
        print(line)
    rowsout += rows2pt
    with open(names[4],'w') as f:
        f.write('# rotations_umd.py, '+UMDFile+', bonding file '+BondFile+', units '+species+' ('+('dipole axis' if mode == 'dipole' else 'bonds')+')\n')
        f.write('# '+str(nsnap)+' snapshots every '+str(TimeStep)+' fs ; gaps up to '+str(maxgap)+' fs bridged, lives of at least '+str(minlife)+' fs ; errors from '+str(nblocks)+' blocks\n')
        f.write('quantity\tvalue\terror\n')
        for name,value,err in rowsout:
            f.write(name+'\t'+str(value)+'\t'+str(err)+'\n')
    print('Files',', '.join(os.path.basename(n) for n in names if os.path.isfile(n)),'written')
    print('runtime : %.1f s' % (time.time()-start))


if __name__ == "__main__":
    if not umdpf.run_on_all(sys.argv[1:],main,paired={'-b':'.bonding.dat'}):
        main(sys.argv[1:])
