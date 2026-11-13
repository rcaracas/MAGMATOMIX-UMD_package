#!/usr/bin/env python3
###
#AUTHORS: RAZVAN CARACAS
###
#Synthetic UMD files with known answers, for the validation suite : synthetic.py <kind> [seed]. Writes <kind>.umd.dat (several files for
#ti) and <kind>.truth.json, the exact answers, derived analytically here, independently of the package. The UMD format is written here
#too, as documented (not with umd_processes_fast), so that the readers of the package are checked on a file they did not write.
#Kinds :
#   idealgas    uniform random atoms (Ar, Ne) in a triclinic cell deformed at every snapshot : g(r) = 1 ; Int(r) = rho 4/3 pi r^3
#   fcc         fcc Cu with small random displacements, in its cubic cell (fcc_cubic) and its 60 deg primitive cell (fcc_prim) :
#               12 first neighbours
#   water       rigid H2O molecules (O-H 0.9572 A, H-O-H 104.52 deg) at random positions and orientations in a triclinic cell, moving
#               rigidly : species H2O, angle 104.52
#   randomwalk  Gaussian random walks with known D : isotropic (Ar 1e-8, Kr 2e-9 m2/s) and anisotropic (Xe, principal values
#               2, 1, 0.5 x 1e-9 m2/s along rotated axes) : D, D tensor
#   ou          Ornstein-Uhlenbeck velocities (VACF = kT/m exp(-t/tau)), positions their integral : D = kT tau / m, kinetic T
#   oscillators harmonic oscillators of known frequencies (Si 500, O 1100 cm-1) with thermal amplitudes : vDOS peaks, Einstein ZPE,
#               Cv, S, E
#   stress      stress tensor as Ornstein-Uhlenbeck processes of known variance and time : viscosity = V s^2 tau / kT ; T and P as
#               AR(1) series of known mean and correlation : errors of the averages
#   ti          atoms in harmonic wells in a periodic cubic box, sampled exactly at the Gauss-Legendre lambdas (k = 0.8) and lambda = 1 :
#               <U>_lambda and the free energy exact
#   faults      a clean random walk (faults_clean) and the same with a duplicated snapshot and a missed crossing of a cell face
#               (faults_bad) : check_umd.py
#   rotors      rigid molecules turning by rotational Brownian motion (independent Gaussian rotation vectors at every step : the mean of
#               P_l(u(0).u(t)) is f_l^t with f_l known) : rotors (H2O in pairs, parallel head to tail at 3.0 A or antiparallel at 3.6 A, in a
#               triclinic cell changing at every snapshot), methane (CH4 in pairs, same orientation at 4.0 A or inverted at 4.6 A),
#               sixstate (H2O with the dipole jumping between 6 directions), gaps (one H2O losing an H for 10 and for 100 snapshots)
#   rotors2pt   rigid molecules on a lattice, the velocity of their centre of mass and their angular velocity in the body frame (each component
#               weighted by its moment of inertia) as Ornstein-Uhlenbeck processes : the correlations of translation and rotation are 3 kT
#               exp(-t/tau) (2 kT for a linear molecule, its angular velocity in the frame of the cell, as rotations_umd.py takes it) ; the truth is
#               the 2PT entropy of the correlations the run really has, written here independently (and of exp(-t/tau), within the noise). The centres stay on the lattice (only the velocities are used). H2O (asymmetric top, symmetry
#               number 2, with ACF.dat files of the Bader program in bohr, atoms shuffled and the origin shifted), CH4 (spherical top, 12), OH (linear, 1)

import sys, json, math
import numpy as np

KB = 8.617333262e-5                 #eV/K
AMU_A2FS2 = 103.642696562           #amu A^2/fs^2 in eV (1.66053906660e-27 kg * 1e-10 m^2 / 1e-30 s^2 / 1.602176634e-19 J)
MASS = {'H':1.008,'C':12.011,'O':15.999,'Si':28.085,'Cu':63.546,'Ar':39.948,'Ne':20.180,'Kr':83.798,'Xe':131.293}
NAVOGADRO = 6.02214076e23
EV = 1.602176634e-19
HPLANCK = 6.62607015e-34
CLIGHT = 2.99792458e10              #cm/s


def cell_from(a,b,c,alpha,beta,gamma):#rows a, b, c ; a along x, b in the xy plane
    al, be, ga = (math.radians(x) for x in (alpha,beta,gamma))
    va = [a,0.0,0.0]
    vb = [b*math.cos(ga),b*math.sin(ga),0.0]
    cx = c*math.cos(be)
    cy = c*(math.cos(al)-math.cos(be)*math.cos(ga))/math.sin(ga)
    vc = [cx,cy,math.sqrt(c*c-cx*cx-cy*cy)]
    return np.array([va,vb,vc])


class UMD:#writer of a UMD file, as documented : header, then per snapshot the fields and one line per atom
    def __init__(self,name,elements,counts,lam=1.0):
        self.f = open(name,'w')
        self.elements = elements
        self.counts = counts
        self.natom = sum(counts)
        self.masses = np.repeat([MASS[e] for e in elements],counts)
        typat = ' '.join(str(k) for k,n in enumerate(counts) for _ in range(n))
        self.f.write('natom %d\nntypat %d\nno.of.electrons 0.0\ntypes %s \nelements %s \nmasses %s \nZelectrons %s \ntypat %s \n\n'
                     'lambda_ThermoInt %.12g\norigin synthetic.py : a synthetic file of the validation suite\n\n'
                     % (self.natom,len(elements),' '.join(map(str,counts)),' '.join(elements),' '.join(str(MASS[e]) for e in elements),
                        ' '.join('0.0' for _ in elements),typat,lam))

    def snapshot(self,time,dt,cell,xcart,unwrapped,vels,forces=None,fields=None):
        fields = dict(fields or {})
        volume = abs(np.linalg.det(cell))
        ke = 0.5*np.sum(self.masses[:,None]*vels**2)*AMU_A2FS2
        fields.setdefault('KineticEnergyIons',ke)
        fields.setdefault('Temperature',2*ke/(3*(self.natom-1)*KB) if self.natom > 1 else 0.0)
        fields.setdefault('Density',np.sum(self.masses)/NAVOGADRO/(volume*1e-24))
        u = fields.get('TotalElectronicEnergy',0.0)
        fields.setdefault('InternalEnergy',u+fields['KineticEnergyIons'])
        fields.setdefault('EnergyWithDrift',fields['InternalEnergy'])
        stress = fields.pop('StressTensor',[0.0]*6)
        fields.setdefault('Pressure',sum(stress[:3])/3)
        xred = np.linalg.solve(cell.T,xcart.T).T
        xred -= np.floor(xred)
        xcart = xred @ cell
        if forces is None:
            forces = np.zeros_like(xcart)
        acell = np.linalg.norm(cell,axis=1)
        lines = ['timestep %.12g fs\ntime %.12g fs\n' % (dt,time)]
        for key, unit in (('InternalEnergy','eV'),('TotalElectronicEnergy','eV'),('ElectronicEntropy','eV'),('KineticEnergyIons','eV'),
                          ('EnergyWithDrift','eV'),('Enthalpy','eV'),('Magnetization','Magneton-Bohr'),('Temperature','K'),
                          ('Pressure','GPa'),('Density','g.cm-3')):
            lines.append('%s %.10g %s\n' % (key,fields.get(key,0.0),unit))
        lines.append('StressTensor '+' '.join('%.8g' % x for x in stress)+'  GPa\n')
        lines.append('acell %.12g %.12g %.12g A\n' % tuple(acell))
        for k, n in enumerate('abc'):
            lines.append('rprim_%s %s\n' % (n,'  '.join('%.12g' % x for x in cell[k]/acell[k])))
        for k, n in enumerate('abc'):
            lines.append('rprimd_%s %s A\n' % (n,'  '.join('%.12g' % x for x in cell[k])))
        lines.append('atoms: reduced*3 cartesian*3(A) abs.diff.*3(A) velocity*3(A/fs) force*3(eV/A) charge(no.elec) magnetization(magneton-Bohr) \n')
        for i in range(self.natom):
            lines.append(' '.join('%.9f' % x for x in (*xred[i],*xcart[i],*unwrapped[i])) +' '+' '.join('%.9g' % x for x in (*vels[i],*forces[i]))+' 0.0 0.0\n')
        lines.append('\n')
        self.f.write(''.join(lines))

    def close(self):
        self.f.close()


def random_rotation(rng):#uniform rotation matrix (QR of a Gaussian matrix)
    q, r = np.linalg.qr(rng.normal(size=(3,3)))
    q *= np.sign(np.diag(r))
    if np.linalg.det(q) < 0:
        q[:,0] *= -1
    return q


def write_truth(name,truth):
    with open(name+'.truth.json','w') as f:
        json.dump(truth,f,indent=1)


#------------------------------------------------------------------------------------------------------------------------------- kinds

def idealgas(rng):
    elements, counts = ['Ar','Ne'], [60,40]
    n = sum(counts)
    base = cell_from(15.0,16.0,17.0,80.0,95.0,105.0)
    umd = UMD('idealgas.umd.dat',elements,counts)
    vols = []
    nsnap = 400
    for s in range(nsnap):
        cell = base @ (np.eye(3)+rng.uniform(-0.05,0.05,(3,3)))     #a new deformation at every snapshot
        vols.append(abs(np.linalg.det(cell)))
        x = rng.random((n,3)) @ cell
        umd.snapshot(float(s),1.0,cell,x,x,np.zeros((n,3)))
    umd.close()
    invV = float(np.mean(1/np.array(vols)))
    #Int(A-B)(r) = <N_B/V> 4/3 pi r^3 for B != A, <(N_A-1)/V> for A-A ; g normalized with each volume, so the mean of 1/V
    write_truth('idealgas',{'g':1.0,'rho':{'Ar':counts[0]*invV,'Ne':counts[1]*invV},'counts':dict(zip(elements,counts)),'snapshots':nsnap})


def fcc(rng):
    a = 3.615
    sigma = 0.04
    nsnap = 60
    #cubic : 3x3x3 conventional cells (108 atoms) ; primitive : 4x4x4 of the 60 deg cell (64 atoms)
    basis = np.array([[0,0,0],[0.5,0.5,0],[0.5,0,0.5],[0,0.5,0.5]])
    cubic = np.array([[i,j,k] for i in range(3) for j in range(3) for k in range(3)])
    sites_c = ((cubic[:,None,:]+basis[None,:,:]).reshape(-1,3))*a
    cell_c = np.eye(3)*3*a
    prim = np.array([[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]])*a
    grid = np.array([[i,j,k] for i in range(4) for j in range(4) for k in range(4)])
    sites_p = grid @ prim
    cell_p = 4*prim
    for name, sites, cell in (('fcc_cubic',sites_c,cell_c),('fcc_prim',sites_p,cell_p)):
        n = len(sites)
        umd = UMD(name+'.umd.dat',['Cu'],[n])
        for s in range(nsnap):
            x = sites+rng.normal(0,sigma,(n,3))
            umd.snapshot(float(s),1.0,cell,x,x,np.zeros((n,3)))
        umd.close()
    write_truth('fcc',{'first_shell':12,'second_shell':6,'nn_distance':a/math.sqrt(2),'second_distance':a,'cutoff':(a/math.sqrt(2)+a)/2})


def water(rng):
    roh, angle = 0.9572, 104.52
    th = math.radians(angle/2)
    mol = np.array([[0,0,0],[roh*math.sin(th),roh*math.cos(th),0],[-roh*math.sin(th),roh*math.cos(th),0]])
    nmol = 20
    cell = cell_from(13.0,14.0,15.0,75.0,100.0,110.0)
    inv = np.linalg.inv(cell)
    centers = []
    while len(centers) < nmol:                                          #centres at least 3.2 A apart (nearest images)
        c = rng.random(3)
        if all(np.linalg.norm((d-np.rint(d)) @ cell) > 3.2 for d in (c-np.array(centers) if centers else np.zeros((0,3)))):
            centers.append(c)
    centers = np.array(centers) @ cell
    rots = [random_rotation(rng) for _ in range(nmol)]
    umd = UMD('water.umd.dat',['O','H'],[nmol,2*nmol])
    nsnap = 50
    for s in range(nsnap):
        shift = rng.normal(0,0.02,(nmol,3))*0                             #rigid molecules, small common motion
        drift = s*np.array([0.013,0.021,-0.017])                         #the whole system drifts across the faces
        pos = np.array([c+drift+(mol @ R.T) for c,R in zip(centers,rots)])
        O = pos[:,0]
        H = pos[:,1:].reshape(-1,3)
        x = np.vstack([O,H])
        umd.snapshot(float(s),1.0,cell,x,x,np.zeros((3*nmol,3)))
    umd.close()
    write_truth('water',{'angle':angle,'roh':roh,'molecules':nmol,'snapshots':nsnap})


def randomwalk(rng):
    elements, counts = ['Ar','Kr','Xe'], [40,40,40]
    n = sum(counts)
    dt = 1.0
    nsnap = 3000
    D = {'Ar':1e-8,'Kr':2e-9}                                           #m2/s ; 1e-8 m2/s = 1e-3 A^2/fs
    R = random_rotation(rng)
    principal = np.array([2e-9,1e-9,0.5e-9])
    tensor = R @ np.diag(principal) @ R.T                               #D tensor of Xe (columns of R : principal directions)
    cov = {'Ar':np.eye(3)*2*D['Ar']*1e5*dt,'Kr':np.eye(3)*2*D['Kr']*1e5*dt,'Xe':2*tensor*1e5*dt}   #A^2 per step
    cell = cell_from(18.0,18.0,18.0,90,90,90)
    x = rng.random((n,3)) @ cell
    chol = {el:np.linalg.cholesky(c) for el,c in cov.items()}
    kinds = np.repeat(elements,counts)
    umd = UMD('randomwalk.umd.dat',elements,counts)
    for s in range(nsnap):
        if s:
            step = np.empty((n,3))
            for el in elements:
                m = kinds == el
                step[m] = rng.normal(size=(m.sum(),3)) @ chol[el].T
            x = x+step
            v = step/dt
        else:
            v = np.zeros((n,3))
        umd.snapshot(s*dt,dt,cell,x,x,v)
    umd.close()
    write_truth('randomwalk',{'D':D,'Xe_tensor':tensor.tolist(),'Xe_principal':principal.tolist(),'Xe_directions':R.T.tolist(),
                              'Xe_isotropic':float(principal.mean())})


def ou(rng):
    elements, counts = ['Ar','Ne'], [24,24]
    n = sum(counts)
    T = 1000.0
    tau = {'Ar':40.0,'Ne':15.0}                                         #fs
    dt = 1.0
    nsnap = 8000
    masses = np.repeat([MASS[e] for e in elements],counts)
    taus = np.repeat([tau[e] for e in elements],counts)
    sv = np.sqrt(KB*T/(masses*AMU_A2FS2))                              #A/fs, per component
    phi = np.exp(-dt/taus)
    cell = cell_from(20.0,20.0,20.0,90,90,90)
    x = rng.random((n,3)) @ cell
    v = rng.normal(size=(n,3))*sv[:,None]
    umd = UMD('ou.umd.dat',elements,counts)
    for s in range(nsnap):
        if s:
            vnew = phi[:,None]*v+np.sqrt(1-phi**2)[:,None]*sv[:,None]*rng.normal(size=(n,3))
            x = x+0.5*(v+vnew)*dt                                       #trapezoid : positions follow the velocities
            v = vnew
        umd.snapshot(s*dt,dt,cell,x,x,v,fields={'Temperature':T})
    umd.close()
    #D = <v_x^2> tau (per component) ; in m2/s : A^2/fs = 1e-5 m2/s
    Dexact = {el:KB*T/(MASS[el]*AMU_A2FS2)*tau[el]*1e-5 for el in elements}
    write_truth('ou',{'D':Dexact,'T':T,'tau':tau})


def oscillators(rng):
    elements, counts = ['Si','O'], [20,20]
    n = sum(counts)
    T = 300.0
    nu = {'Si':500.0,'O':1100.0}                                        #cm-1
    dt = 0.5
    nsnap = 8000
    masses = np.repeat([MASS[e] for e in elements],counts)
    w = np.repeat([2*math.pi*CLIGHT*nu[e]*1e-15 for e in elements],counts)   #rad/fs
    sv = np.sqrt(KB*T/(masses*AMU_A2FS2))
    #classical thermal amplitudes : v0 Gaussian with kT/m, x0 Gaussian with kT/(m w^2) ; x = x0 cos + v0/w sin
    x0 = rng.normal(size=(n,3))*(sv/w)[:,None]
    v0 = rng.normal(size=(n,3))*sv[:,None]
    cell = cell_from(16.0,16.0,16.0,90,90,90)
    sites = rng.random((n,3)) @ cell
    umd = UMD('oscillators.umd.dat',elements,counts)
    for s in range(nsnap):
        t = s*dt
        c, si = np.cos(w*t)[:,None], np.sin(w*t)[:,None]
        u = x0*c+v0/w[:,None]*si
        v = -x0*w[:,None]*si+v0*c
        f = -(masses*AMU_A2FS2*w**2)[:,None]*u
        umd.snapshot(t,dt,cell,sites+u,sites+u,v,forces=f,fields={'Temperature':T})
    umd.close()
    truth = {'nu':nu,'T':T}
    for el in elements:#Einstein solid, 3 modes per atom, quantum ; per atom : ZPE, E (with ZPE) in meV, S, Cv in kB
        e = HPLANCK*CLIGHT*nu[el]/EV                                   #eV
        xx = e/(KB*T)
        truth[el] = {'ZPE':3*e/2*1000,'E':3*(e/2+e/math.expm1(xx))*1000,'Cv':3*xx*xx*math.exp(xx)/math.expm1(xx)**2,
                     'S':3*(xx/math.expm1(xx)-math.log(1-math.exp(-xx))),'F':3*(e/2+KB*T*math.log(1-math.exp(-xx)))*1000}
    write_truth('oscillators',truth)


def stress(rng):
    elements, counts = ['Ar'], [2]
    n = 2
    T = 2000.0
    tau = 25.0                                                          #fs
    s_off = 1.5                                                         #GPa, standard deviation of xy, yz, zx
    dt = 1.0
    nsnap = 40000
    cell = cell_from(15.0,15.0,15.0,90,90,90)
    V = abs(np.linalg.det(cell))
    phi = math.exp(-dt/tau)
    sd = np.array([math.sqrt(2)*s_off]*3+[s_off]*3)                    #diagonal variance 2 s^2 : <N_xx N_xx> = 4/3 <s_xy s_xy> (isotropic)
    p0 = 5.0
    sig = rng.normal(size=6)*sd
    #T and P fields : AR(1) series of known mean and correlation (averages.py) ; P = p0 + mean of the diagonal
    phiT, sT, T0 = 0.95, 30.0, T
    tser = T0+rng.normal()*sT
    x = rng.random((n,3)) @ cell
    zero = np.zeros((n,3))
    umd = UMD('stress.umd.dat',elements,counts)
    for s in range(nsnap):
        if s:
            sig = phi*sig+math.sqrt(1-phi*phi)*sd*rng.normal(size=6)
            tser = T0+phiT*(tser-T0)+math.sqrt(1-phiT**2)*sT*rng.normal()
        st = sig.copy()
        st[:3] += p0
        umd.snapshot(s*dt,dt,cell,x,x,zero,fields={'Temperature':tser,'StressTensor':list(st),'Pressure':float(st[:3].mean())})
    umd.close()
    #eta = V/(kB T) int <s_xy(0) s_xy(t)> dt = V s^2 tau / (kB T) (OU, continuous limit ; the discrete sum differs by (dt/tau)^2/12)
    eta = V*1e-30*(s_off*1e9)**2*tau*1e-15/(KB*EV*T)
    #AR(1) : error of the mean = sT/sqrt(n) sqrt((1+phi)/(1-phi))
    write_truth('stress',{'eta':eta,'T_mean':T0,'T_err':sT/math.sqrt(nsnap)*math.sqrt((1+phiT)/(1-phiT)),'P_mean':p0,
                          'P_err':s_off*math.sqrt(2/3)/math.sqrt(nsnap)*math.sqrt((1+phi)/(1-phi)),'tau':tau,'V':V,'T':T})


def ti(rng):
    from scipy.special import erf, roots_legendre
    elements, counts = ['Ar','Ne'], [4,4]
    n = 8
    T = 1500.0
    kT = KB*T
    L = 7.0
    kspring = 2.0                                                       #eV/A^2
    nsnap = 2000
    cell = np.eye(3)*L
    sites = (np.array([[0,0,0],[0.5,0.5,0],[0.5,0,0.5],[0,0.5,0.5],[0.5,0.5,0.5],[0,0,0.5],[0,0.5,0],[0.5,0,0]])+0.25) @ cell
    masses = np.repeat([MASS[e] for e in elements],counts)
    sv = np.sqrt(kT/(masses*AMU_A2FS2))
    xnodes, _ = roots_legendre(10)
    lambdas = [((xx+1)/2)**5 for xx in xnodes]+[1.0]                  #k = 0.8 : lambda = ((x+1)/2)^(1/(1-k))
    a = L/2
    def mean_d2(lam):#<d^2> of a coordinate with the density exp(-lam k d^2/2kT) on [-a, a]
        if lam*kspring*a*a/kT < 1e-8:
            return a*a/3
        s2 = kT/(lam*kspring)
        z = a/math.sqrt(s2)
        return s2*(1-2*z*math.exp(-z*z/2)/(math.sqrt(2*math.pi)*erf(z/math.sqrt(2))))
    def sample(lam,size):#truncated Gaussian on [-a, a] by inversion of the erf
        from scipy.special import erfinv
        s = math.sqrt(kT/(lam*kspring))
        e = erf(a/(s*math.sqrt(2)))
        return s*math.sqrt(2)*erfinv(e*(2*rng.random(size)-1))
    truth = {'U':{}}
    for lam in lambdas:
        name = 'ti_%.6e' % lam
        umd = UMD(name+'.umd.dat',elements,counts,lam=lam)
        for s in range(nsnap):
            d = sample(lam,(n,3)) if lam*kspring*a*a/kT > 1e-8 else rng.uniform(-a,a,(n,3))
            U = 0.5*kspring*np.sum(d*d)
            v = rng.normal(size=(n,3))*sv[:,None]
            v -= np.sum(masses[:,None]*v,axis=0)/np.sum(masses)
            ke = 0.5*np.sum(masses[:,None]*v**2)*AMU_A2FS2
            umd.snapshot(float(s),1.0,cell,sites+d,sites+d,v,forces=-kspring*d,
                         fields={'TotalElectronicEnergy':U,'KineticEnergyIons':ke,'InternalEnergy':U+ke,'Temperature':T})
        umd.close()
        truth['U']['%.6e' % lam] = 3*n*0.5*kspring*mean_d2(lam)
    #Free energy : ideal gas (3N, -kT sum [N ln(V/Lambda^3) - ln N!]) + Delta F = -3N kT ln(z(1)/L)
    hbar2pi = HPLANCK
    Fig = 0.0
    for el, N in zip(elements,counts):
        m = MASS[el]*1e-3/NAVOGADRO
        lam_th = HPLANCK/math.sqrt(2*math.pi*m*KB*EV*T)*1e10            #A
        Fig += -kT*(N*math.log(L**3/lam_th**3)-math.lgamma(N+1))
    z1 = math.sqrt(2*math.pi*kT/kspring)*erf(a*math.sqrt(kspring/(2*kT)))
    dF = -3*n*kT*math.log(z1/L)
    truth.update({'F_ideal':Fig,'dF':dF,'F':Fig+dF,'E':truth['U']['%.6e' % 1.0]+1.5*(n-1)*kT,'T':T,'natom':n})
    #E : <U>_1 + <KE> ; KE with 3N-3 degrees of freedom (no net momentum) ; S = (E - F)/T
    truth['S_kB_per_atom'] = (truth['E']-truth['F'])/(kT*n)
    write_truth('ti',truth)


def faults(rng):
    elements, counts = ['Ar'], [50]
    n = 50
    dt = 1.0
    cell = cell_from(12.0,12.0,12.0,90,90,90)
    x0 = rng.random((n,3)) @ cell
    steps = rng.normal(0,0.02,(300,n,3))
    for name in ('faults_clean','faults_bad'):
        umd = UMD(name+'.umd.dat',elements,counts)
        x = x0.copy()
        t = 0.0
        for s in range(300):
            if s:
                x = x+steps[s]
            un = x.copy()
            if name == 'faults_bad' and s >= 150:
                un[7] -= cell[0]                                         #a missed crossing : the unwrapped position jumps by a
            umd.snapshot(t,dt,cell,x,un,steps[s]/dt if s else np.zeros((n,3)))
            if name == 'faults_bad' and s == 100:                       #a duplicated snapshot (same time, same positions)
                umd.snapshot(t,dt,cell,x,un,steps[s]/dt)
            t += dt
        umd.close()
    write_truth('faults',{'clean':'no warning for 1, 2, 3','bad':'warnings for 1 (duplicated time), 2 or 3 (missed crossing)'})


def rotation_of(v):#rotation matrix of the rotation vector v (Rodrigues)
    th = np.linalg.norm(v)
    if th < 1e-14:
        return np.eye(3)
    k = v/th
    K = np.array([[0,-k[2],k[1]],[k[2],0,-k[0]],[-k[1],k[0],0]])
    return np.eye(3)+math.sin(th)*K+(1-math.cos(th))*(K @ K)


def rotor_factor(l,sigma):#Mean of P_l(u.Ru) for the rotations R of Gaussian rotation vectors (sigma per component) : the character of the rotation by the
    #angle w, sin((l+1/2)w)/sin(w/2)/(2l+1), over the Maxwell distribution of w. After t independent steps : this factor to the power t
    w = np.linspace(1e-9,14*sigma,400001)
    y = np.sqrt(2/np.pi)*w**2*np.exp(-w**2/(2*sigma**2))/sigma**3*np.sin((l+0.5)*w)/np.sin(0.5*w)/(2*l+1)
    return float(np.sum(0.5*(y[1:]+y[:-1])*np.diff(w)))


def rotors(rng):
    roh, th = 0.9572, math.radians(104.52/2)
    h2o = np.array([[0,0,0],[roh*math.sin(th),roh*math.cos(th),0],[-roh*math.sin(th),roh*math.cos(th),0]])   #dipole axis along y, H-H along x
    flip = np.diag([-1.0,-1.0,1.0])                                     #half turn about z : dipole along -y
    truth = {}

    #rotors : 15 pairs A (the two dipoles along the line from the first to the second molecule, 3.0 A apart : cos 1, dipolar 2, both H of the first
    #point to the second) and 9 pairs B (dipoles antiparallel, perpendicular to the line, molecular planes perpendicular to it, 3.6 A apart : cos -1,
    #dipolar 1, no H pointing to the other). Each pair turns as a rigid body ; the cell changes at every snapshot, the pairs drift across its faces
    nA, nB, nsnap, sigma = 15, 9, 2000, 0.14
    npair = nA+nB
    grid = np.array([[(i+0.5)/4,(j+0.5)/3,(k+0.5)/2] for i in range(4) for j in range(3) for k in range(2)])
    R = [random_rotation(rng) for _ in range(npair)]
    umd = UMD('rotors.umd.dat',['O','H'],[2*npair,4*npair])
    for s in range(nsnap):
        cell = cell_from(40.0,33.0,30.0,84.0,97.0,101.0)*(1+0.01*math.sin(s/40.0))
        drift = s*np.array([0.011,-0.017,0.013])
        O, H = [], []
        for p in range(npair):
            if s:
                R[p] = rotation_of(rng.normal(0,sigma,3)) @ R[p]
            centre = grid[p] @ cell + drift
            if p < nA:
                mols = [h2o+np.array([0,-1.5,0]), h2o+np.array([0,1.5,0])]
            else:
                mols = [h2o+np.array([0,0,-1.8]), h2o @ flip.T+np.array([0,0,1.8])]
            for m in mols:
                x = centre + m @ R[p].T
                O.append(x[0])
                H += [x[1],x[2]]
        x = np.array(O+H)
        umd.snapshot(float(s),1.0,cell,x,x,np.zeros_like(x))
    umd.close()
    n = 2*npair
    truth['rotors'] = {'factors':{str(l):rotor_factor(l,sigma) for l in (1,2)},'molecules':n,'cos_3.0':1.0,'P2_3.0':1.0,'dipolar_3.0':2.0,
                       'cos_3.6':-1.0,'P2_3.6':1.0,'dipolar_3.6':1.0,'neighbours':1.0,'cos_shell':(2*nA-2*nB)/n,'links_per_unit':2*nA/n,
                       'cos_linked':1.0,'kirkwood_3.3':1+2*nA/n,'kirkwood_3.9':1+(2*nA-2*nB)/n,'shell_radius':4.0}

    #methane : 10 pairs with the same orientation, 4.0 A apart (K3 = K4 = 1), 6 pairs with the second molecule inverted, 4.6 A apart (K3 = -1, K4 = 1)
    t = np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]])/math.sqrt(3)*1.09
    nsame, ninv, sigma = 10, 6, 0.07
    npair = nsame+ninv
    grid = np.array([[(i+0.5)/4,(j+0.5)/2,(k+0.5)/2] for i in range(4) for j in range(2) for k in range(2)])
    cell = cell_from(40.0,22.0,21.0,90.0,90.0,90.0)
    R = [random_rotation(rng) for _ in range(npair)]
    umd = UMD('methane.umd.dat',['C','H'],[2*npair,8*npair])
    for s in range(nsnap):
        C, H = [], []
        for p in range(npair):
            if s:
                R[p] = rotation_of(rng.normal(0,sigma,3)) @ R[p]
            d = 2.0 if p < nsame else 2.3
            for sign, shape in ((-1,t),(1,t if p < nsame else -t)):
                c = grid[p] @ cell + np.array([0,0,sign*d]) @ R[p].T
                C.append(c)
                H += list(c + shape @ R[p].T)
        x = np.array(C+H)
        umd.snapshot(float(s),1.0,cell,x,x,np.zeros_like(x))
    umd.close()
    truth['methane'] = {'factors':{str(l):rotor_factor(l,sigma) for l in (1,2,3,4)},'self':{'1':0.0,'2':0.0,'3':5.0/9.0,'4':7.0/27.0},
                        'K3_4.0':1.0,'K4_4.0':1.0,'K3_4.6':-1.0,'K4_4.6':1.0}

    #sixstate : 36 H2O on a lattice, the dipole along one of 6 directions (+-x, +-y, +-z of a tilted frame) ; at every step, with probability p, a jump to
    #one of the 5 other directions : C1 = C2 = (1 - 6p/5)^t ; the map has 6 bins : entropy ln 6 - ln(number of bins) ; independent molecules
    prob = 0.02
    frame = random_rotation(rng)
    states = []
    for axis in range(3):
        for sign in (1,-1):
            y = sign*frame[:,axis]
            xdir = frame[:,(axis+1) % 3]
            states.append(np.column_stack([xdir,y,np.cross(xdir,y)]))  #body x, y (dipole), z -> these columns
    sites = np.array([[5.0*i+2.5,5.0*j+2.5,5.0*k+2.5] for i in range(4) for j in range(3) for k in range(3)])
    cell = np.diag([20.0,15.0,15.0])
    state = rng.integers(0,6,len(sites))
    umd = UMD('sixstate.umd.dat',['O','H'],[len(sites),2*len(sites)])
    for s in range(nsnap):
        if s:
            jump = rng.random(len(sites)) < prob
            state = np.where(jump,(state+rng.integers(1,6,len(sites))) % 6,state)
        O, H = [], []
        for c, k in zip(sites,state):
            x = c + h2o @ states[k].T
            O.append(x[0])
            H += [x[1],x[2]]
        x = np.array(O+H)
        umd.snapshot(float(s),1.0,cell,x,x,np.zeros_like(x))
    umd.close()
    truth['sixstate'] = {'factor':1-6*prob/5,'entropy_plus_ln_bins':math.log(6.0),'molecules':len(sites)}

    #gaps : 8 H2O turning ; the first loses an H (moved 3 A away along its bond) during the snapshots 100-109 and 300-399
    sites = np.array([[6.0*i+3,6.0*j+3,6.0*k+3] for i in range(2) for j in range(2) for k in range(2)])
    cell = np.diag([12.0,12.0,12.0])
    R = [random_rotation(rng) for _ in sites]
    umd = UMD('gaps.umd.dat',['O','H'],[len(sites),2*len(sites)])
    for s in range(500):
        O, H = [], []
        for m, c in enumerate(sites):
            if s:
                R[m] = rotation_of(rng.normal(0,0.05,3)) @ R[m]
            x = c + h2o @ R[m].T
            if m == 0 and (100 <= s < 110 or 300 <= s < 400):
                x[1] = x[0]+(x[1]-x[0])*(1+3.0/roh)
            O.append(x[0])
            H += [x[1],x[2]]
        x = np.array(O+H)
        umd.snapshot(float(s),1.0,cell,x,x,np.zeros_like(x))
    umd.close()
    truth['gaps'] = {'lives':{'-g 50':9,'-g 5':10,'-g 200':8,'-g 5 -m 150':8},'left_out':{'-g 50':10,'-g 5':0,'-g 200':110}}   #-m 150 : the lives 0-99 and 400-499 dropped
    write_truth('rotors',truth)



CM1 = 33356.40952                   #cm^-1 per fs^-1
BOHR = 0.529177210903

def entropy_2pt(corr,dof,T,mass,volume,dt,SR=None):#2PT entropy (kB per molecule) of a motion of correlation corr (eV, lags 0 .. window-1) (Lin et al.
    #2003, 2010) : spectrum of the correlation with a one-sided Hann taper over the window, fluidicity, gas-like Lorentzian (hard spheres for the
    #translation, rigid rotor of entropy SR for the rotation), solid-like quantum oscillators. Returns f and S
    from scipy.optimize import brentq
    kT = KB*T
    window = len(corr)
    t = np.arange(window)*dt
    c = corr*0.5*(1+np.cos(np.pi*np.arange(window)/(window-1)))
    nu = np.linspace(0.0,CM1/(2*dt),20001)                                            #cm^-1
    S = dt*(c[0]+2*np.cos(2*np.pi*np.outer(nu/CM1,t[1:])) @ c[1:])                       #eV fs : int over -tmax .. tmax
    g = 2*S/kT/CM1*3.0/dof                                                              #states per cm^-1, normalized to 3
    s0 = g[0]/CLIGHT                                                                    #per Hz
    m = mass*1.66053906660e-27
    kTJ = kT*EV
    delta = 2*s0/9*math.sqrt(math.pi*kTJ/m)*(1/(volume*1e-30))**(1/3)*(6/math.pi)**(2/3)
    f = brentq(lambda f: 2*delta**-4.5*f**7.5-6*delta**-3*f**5-delta**-1.5*f**3.5+6*delta**-1.5*f**2.5+2*f-2,1e-12,1.0)
    gas = g[0]/(1+(math.pi*g[0]*nu/(6*f))**2)
    x = HPLANCK*CLIGHT*np.maximum(nu,1e-9)/kTJ
    ws = x/np.expm1(x)-np.log(-np.expm1(-x))
    if SR is None:
        y = f**2.5/delta**1.5
        z = (1+y+y*y-y**3)/(1-y)**3
        SR = 2.5+math.log((2*math.pi*m*kTJ/HPLANCK**2)**1.5*volume*1e-30/f*z)+y*(3*y-4)/(1-y)**2
        wgas = SR/3
    else:
        wgas = SR/dof
    dn = nu[1]-nu[0]
    integ = lambda h: float(np.sum(0.5*(h[1:]+h[:-1]))*dn)
    solid = (g-gas)*dof/3.0
    #the first interval, where ln(x) diverges : the mean of the weight over (0, nu1)
    x1 = x[1]
    ws0 = 2.0-math.log(x1)+x1**2/72.0
    first = 0.5*(solid[0]+solid[1])*ws0*dn
    rest = integ(solid[1:]*ws[1:])
    return f, integ(gas*dof/3.0)*wgas+first+rest


def rotors2pt(rng):
    T, dt, nsnap, window, spacing = 300.0, 1.0, 2000, 300, 12.0
    taut, taur = 30.0, 20.0
    roh, th = 0.9572, math.radians(104.52/2)
    tet = np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]])/math.sqrt(3)*1.09
    shapes = {'h2o':(['O','H'],np.array([[0,0,0],[roh*math.sin(th),roh*math.cos(th),0],[-roh*math.sin(th),roh*math.cos(th),0]]),2,34.0),
              'ch4':(['C','H'],np.vstack([[0,0,0],tet]),12,40.0),
              'oh':(['O','H'],np.array([[0,0,0],[0,0,0.97]]),1,25.0)}
    lattice = np.array([[i,j,k] for i in range(3) for j in range(3) for k in range(3)])*spacing+spacing/2
    cell = np.eye(3)*3*spacing
    truth = {}
    for name,(elements,shape,sigma,volume) in shapes.items():
        nm = len(lattice)
        masses = np.array([MASS[elements[0]]]+[MASS[elements[1]]]*(len(shape)-1))
        M = masses.sum()
        ref = shape-(masses[:,None]*shape).sum(axis=0)/M
        I = (masses[:,None,None]*((ref**2).sum(axis=1)[:,None,None]*np.eye(3)-ref[:,:,None]*ref[:,None,:])).sum(axis=0)
        mom, axes = np.linalg.eigh(I)
        ref = ref @ axes                                                                    #principal frame
        linear = mom[0] < 1e-8*mom[2]
        kT = KB*T/AMU_A2FS2                                                                 #amu A^2/fs^2
        at, ar = math.exp(-dt/taut), math.exp(-dt/taur)
        V = rng.normal(0,math.sqrt(kT/M),(nm,3))
        sd = np.sqrt(kT/np.where(mom > 1e-8*mom[2],mom,np.inf))                             #no rotation about the axis of a linear molecule
        W = rng.normal(0,1,(nm,3))*sd
        R = [random_rotation(rng) for _ in range(nm)]
        counts = [nm,nm*(len(shape)-1)]
        umd = UMD('rotors2pt_'+name+'.umd.dat',elements,counts)
        snaps = {}
        trn, rot = [], []                                                                    #sqrt(M) V and sqrt(I_k) w_k : their correlations
        for s in range(nsnap):
            if s:
                V = at*V+math.sqrt(1-at*at)*rng.normal(0,math.sqrt(kT/M),(nm,3))
                for q in range(nm):
                    R[q] = R[q] @ rotation_of(W[q]*dt)                                      #body angular velocity : dR/dt = R [w]x
                W = ar*W+math.sqrt(1-ar*ar)*rng.normal(0,1,(nm,3))*sd
            trn.append(np.sqrt(M)*V)
            rot.append(np.sqrt(mom[2])*np.array([R[q] @ W[q] for q in range(nm)]) if linear else np.sqrt(mom)*W)   #linear : in the frame of the cell
            pos = np.array([lattice[q]+ref @ R[q].T for q in range(nm)])                   #(molecule, atom, 3)
            vel = np.array([V[q]+np.cross(R[q] @ W[q],ref @ R[q].T) for q in range(nm)])
            x = np.vstack([pos[:,0],pos[:,1:].reshape(-1,3)])
            v = np.vstack([vel[:,0],vel[:,1:].reshape(-1,3)])
            umd.snapshot(s*dt,dt,cell,x,x,v,fields={'Temperature':T})
            if name == 'h2o' and s in (500,1500):
                snaps[s] = x
        umd.close()
        dof = 2 if linear else 3
        theta = HPLANCK**2/(8*math.pi**2*np.maximum(mom,1e-30)*1.66053906660e-47*KB*EV)
        if linear:
            SR = 1+math.log(T/(sigma*theta[2]))
        else:
            SR = math.log(math.sqrt(math.pi)*math.exp(1.5)/sigma*math.sqrt(T**3/np.prod(theta)))
        def corr(x):#mean over the molecules and all the time origins of x(0).x(t), in eV : the correlations the run really has
            x = np.array(x)
            return np.array([(x[:nsnap-j]*x[j:]).sum(axis=2).mean() for j in range(window)])*AMU_A2FS2
        ct, cr = corr(trn), corr(rot)
        ft, St = entropy_2pt(ct,3,T,M,volume,dt)
        fr, Sr = entropy_2pt(cr,dof,T,M,volume,dt,SR)
        fte, Ste = entropy_2pt(3*KB*T*np.exp(-np.arange(window)*dt/taut),3,T,M,volume,dt)  #expected from exp(-t/tau) : within the noise of the run
        fre, Sre = entropy_2pt(dof*KB*T*np.exp(-np.arange(window)*dt/taur),dof,T,M,volume,dt,SR)
        truth[name] = {'symmetry_number':sigma,'linear':bool(linear),'moments':[float(x) for x in mom],'S_rigid_rotor':SR,'volume':volume,
                       'f_trn':ft,'f_rot':fr,'S_trn':St,'S_rot':Sr,'T':T,'window':window,'T_trn':ct[0]/(3*KB),'T_rot':cr[0]/(dof*KB),
                       'expected':{'f_trn':fte,'f_rot':fre,'S_trn':Ste,'S_rot':Sre}}
        if name == 'h2o':#ACF.dat files : O 22, H 6 A^3 (a unit 34), the rest vacuum ; in bohr, atoms shuffled, origin shifted, wrapped
            for s,x in snaps.items():
                order = rng.permutation(len(x))
                xs = (x+np.array([1.7,-2.3,0.9])) @ np.linalg.inv(cell)
                xs = (xs-np.floor(xs)) @ cell/BOHR
                vol = np.where(np.arange(len(x)) < nm,22.0,6.0)/BOHR**3
                vacuum = abs(np.linalg.det(cell))/BOHR**3-vol.sum()
                folder = 'bader2pt'
                import os
                os.makedirs(folder,exist_ok=True)
                with open(os.path.join(folder,'h2o_%d.ACF.dat' % s),'w') as f:
                    f.write('    #         X           Y           Z       CHARGE      MIN DIST   ATOMIC VOL\n '+'-'*80+'\n')
                    for i,a in enumerate(order):
                        f.write('%5d %11.4f %11.4f %11.4f %11.4f %11.4f %11.4f\n' % (i+1,*xs[a],8.0 if a < nm else 0.5,0.5,vol[a]))
                    f.write(' '+'-'*80+'\n    VACUUM CHARGE:               0.0000\n    VACUUM VOLUME:          %11.4f\n    NUMBER OF ELECTRONS:     %.4f\n' % (vacuum,9.0*nm))
    write_truth('rotors2pt',truth)

KINDS = {'idealgas':idealgas,'fcc':fcc,'water':water,'randomwalk':randomwalk,'ou':ou,'oscillators':oscillators,'stress':stress,'ti':ti,
         'faults':faults,'rotors':rotors,'rotors2pt':rotors2pt}

if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] not in KINDS:
        print('synthetic.py <kind> [seed] ; kinds :',', '.join(KINDS))
        sys.exit(2)
    KINDS[sys.argv[1]](np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 12345))
