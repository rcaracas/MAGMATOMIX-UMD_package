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

import sys, json, math
import numpy as np

KB = 8.617333262e-5                 #eV/K
AMU_A2FS2 = 103.642696562           #amu A^2/fs^2 in eV (1.66053906660e-27 kg * 1e-10 m^2 / 1e-30 s^2 / 1.602176634e-19 J)
MASS = {'H':1.008,'O':15.999,'Si':28.085,'Cu':63.546,'Ar':39.948,'Ne':20.180,'Kr':83.798,'Xe':131.293}
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


KINDS = {'idealgas':idealgas,'fcc':fcc,'water':water,'randomwalk':randomwalk,'ou':ou,'oscillators':oscillators,'stress':stress,'ti':ti,
         'faults':faults}

if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] not in KINDS:
        print('synthetic.py <kind> [seed] ; kinds :',', '.join(KINDS))
        sys.exit(2)
    KINDS[sys.argv[1]](np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 12345))
