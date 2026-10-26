#!/usr/bin/env python3
###
#AUTHORS: ADRIEN SAURETY, EMMA STOUTENBOURG, RAZVAN CARACAS
###
#Free energy and entropy by thermodynamic integration from the ideal gas to the ab initio system (VASP SCALEE = lambda) :
#   F = F_ideal gas + integral_0^1 <U>_lambda dlambda,   S = (E - F)/T
#One UMD file per lambda (VaspParser.py : lambda_ThermoInt in the header, TotalElectronicEnergy = TOTEN/SCALEE = U, the unscaled DFT free
#energy). <U>_lambda diverges as lambda -> 0 (the atoms overlap in the ideal gas), so the integral is taken over x = 2 lambda^(1-k) - 1
#(Dorner, Sukurma, Dellago, Kresse, PRL 121, 195701, 2018) : integral_-1^1 <U> lambda^k / (2(1-k)) dx, whose integrand vanishes at x = -1 when
#<U> grows more slowly than lambda^-k (k = 0.8). Gauss-Legendre quadrature in x when the lambdas are its nodes (the recipe of the VASP wiki),
#else a cubic spline through (x=-1, 0) and the points.

import sys, getopt, os, glob, math
import numpy as np
from scipy.interpolate import CubicSpline
import umd_processes_fast as umdpf
from averages import statistics, kinetic_included

KB = 8.617333262E-05        #eV/K
KBJ = 1.380649E-23          #J/K
H = 6.62607015E-34          #J s
NA = 6.02214076E23
EV = 1.602176634E-19        #J
FIELDS = ('TotalElectronicEnergy','InternalEnergy','KineticEnergyIons','ElectronicEntropy','Temperature','rprimd_a','rprimd_b','rprimd_c','acell')


def read_lambda(UMDname):#lambda_ThermoInt of the header (1 when absent)
    with open(UMDname,'r') as f:
        for line in f:
            if line.startswith('lambda_ThermoInt'):
                return float(line.split()[1])
            if line.startswith('timestep'):
                break
    return 1.0


def read_run(UMDname,skip):#Header, lambda, and the series of the fields from snapshot skip on
    MyCrystal, TimeStep = umdpf.Crystallization(UMDname)
    series = {k:[] for k in FIELDS}
    cells = []
    for isnap, fields, atoms in umdpf.read_snapshots(UMDname,MyCrystal.natom,fields=FIELDS,atoms=False):
        if isnap < skip:
            continue
        for k in ('TotalElectronicEnergy','InternalEnergy','KineticEnergyIons','ElectronicEntropy','Temperature'):
            series[k].append(fields[k][0] if fields.get(k) else np.nan)
        cells.append(umdpf.cell_of(fields))
    run = {k:np.array(series[k]) for k in ('TotalElectronicEnergy','InternalEnergy','KineticEnergyIons','ElectronicEntropy','Temperature')}
    run['volume'] = np.array([abs(np.linalg.det(c)) if c is not None else np.nan for c in cells])
    run['name'] = UMDname
    run['lambda'] = read_lambda(UMDname)
    run['crystal'] = MyCrystal
    return run


def ideal_gas(crystal,volume,T):#Helmholtz free energy of the ideal gas of the atoms (eV per cell) : -kT sum_el [N ln(V/Lambda^3) - ln N!]
    f = 0.0
    for el,n,m in zip(crystal.elements,crystal.types,crystal.masses):
        mass = float(m)/1000.0/NA                                       #kg
        wavelength = H/math.sqrt(2.0*math.pi*mass*KBJ*T)                #m
        f += -KB*T*(n*math.log(volume*1E-30/wavelength**3) - math.lgamma(n+1))
    return f


def gauss_legendre(lambdas,k):#Mask of the runs whose lambdas are the n Gauss-Legendre nodes in x (all the runs, or all but lambda = 1), and the
    #weights (lambdas sorted) ; None if they are not
    for use in (lambdas < 1.0, np.ones(len(lambdas),bool)):
        if use.sum() >= 2:
            nodes, weights = np.polynomial.legendre.leggauss(int(use.sum()))
            if np.allclose(2.0*lambdas[use]**(1.0-k)-1.0,nodes,atol=2e-3):
                return use, weights
    return None


def quadrature_weights(x,method):#Weights w with integral_-1^1 g dx = sum w_i g_i for the points x (sorted, including -1) : spline or trapezoid
    n = len(x)
    if method == 'trapezoid':
        w = np.zeros(n)
        dx = np.diff(x)
        w[:-1] += dx/2.0
        w[1:] += dx/2.0
        if x[-1] < 1.0:                                                 #constant extrapolation to x = 1
            w[-1] += 1.0-x[-1]
        return w
    w = np.zeros(n)
    for i in range(n):                                                  #the spline is linear in the data : integrate the spline of each unit vector
        e = np.zeros(n)
        e[i] = 1.0
        w[i] = CubicSpline(x,e,bc_type='not-a-knot' if n > 3 else 'natural',extrapolate=True).integrate(-1.0,1.0)
    return w


def main(argv):
    umdpf.headerumd()
    files = None
    skip = 0
    k = 0.8
    Tuser = None
    visual = False
    try:
        opts, arg = getopt.getopt(argv,"hf:i:k:t:v")
    except getopt.GetoptError:
        print('thermo_integration_umd.py -f <files> -i <InitialStep> -k <exponent> -t <temperature> -v')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('thermo_integration_umd.py : free energy and entropy by thermodynamic integration from the ideal gas (VASP SCALEE = lambda)')
            print('thermo_integration_umd.py -f <files> -i <InitialStep> -k <exponent> -t <temperature> -v')
            print('-f : UMD files, one per lambda, separated by commas (default : every *.umd.dat of the folder) ; same cell, atoms and')
            print('     temperature ; VaspParser.py writes lambda_ThermoInt (SCALEE) in the header and TotalElectronicEnergy = TOTEN/SCALEE')
            print('     (U, the DFT free energy unscaled ; files parsed between January and October 2026 have 0 there : re-parse the OUTCAR)')
            print('-i : snapshots skipped at the start of each run (default 0)')
            print('-k : exponent of the change of variable x = 2 lambda^(1-k) - 1 (default 0.8, Dorner et al. 2018 ; the integrand must vanish')
            print('     at lambda -> 0, i.e. <U> grow more slowly than lambda^-k). Gauss-Legendre quadrature in x when the lambdas are its n')
            print('     nodes, lambda = ((x+1)/2)^(1/(1-k)) (plus lambda = 1 for E), else a cubic spline in x ; the trapezoid rule is printed as a check.')
            print('     Use at least 8 nodes : for atoms in harmonic wells (exact) the quadrature errors of 4, 6, 8, 10 Gauss-Legendre nodes are 0.06,')
            print('     0.03, 0.002, 0.001 eV/atom (the spline 0.14, 0.05, 0.008, 0.0015) ; the integrand changes fast where the atoms localize.')
            print('-t : temperature (K) of the ideal gas and of S : the thermostat temperature (TEBEG) is best (default : the mean Temperature')
            print('     of all the runs ; with few atoms it differs from TEBEG by its noise, ~0.5 % for 8 atoms and 4000 snapshots)')
            print('-v : plot of the integrand, saved in thermo_integration.pdf')
            print('Errors : error of each <U> from its integrated autocorrelation time (as averages.py -a), propagated through the')
            print('quadrature weights ; the quadrature error itself is not included (compare with the trapezoid value and with more lambdas).')
            print('Output : the table of the runs and the results, also in thermo_integration.dat')
            sys.exit()
        elif opt == '-f':
            files = [x.strip() for x in arg.split(',') if x.strip()]
        elif opt == '-i':
            skip = int(arg)
        elif opt == '-k':
            k = float(arg)
        elif opt == '-t':
            Tuser = float(arg)
        elif opt == '-v':
            visual = True
    if files is None:
        files = sorted(glob.glob('*.umd.dat'))
    missing = [f for f in files if not os.path.isfile(f)]
    if missing or not files:
        print('Files not found :',', '.join(missing) if missing else 'no *.umd.dat in the folder')
        sys.exit()
    if not 0.0 < k < 1.0:
        print('-k must be between 0 and 1')
        sys.exit()

    runs = sorted((read_run(f,skip) for f in files),key=lambda r: r['lambda'])
    lambdas = np.array([r['lambda'] for r in runs])
    if len(set(lambdas)) < len(lambdas):
        print('Several files have the same lambda :',', '.join(r['name']+' ('+str(r['lambda'])+')' for r in runs))
        sys.exit()
    ref = runs[0]['crystal']
    for r in runs:
        c = r['crystal']
        if list(c.elements) != list(ref.elements) or list(c.types) != list(ref.types):
            print('The files do not have the same atoms :',r['name'])
            sys.exit()
        if len(r['TotalElectronicEnergy']) < 10:
            print(r['name'],': only',len(r['TotalElectronicEnergy']),'snapshots after -i',skip)
            sys.exit()
        if np.all(r['TotalElectronicEnergy'] == 0.0) or np.all(np.isnan(r['TotalElectronicEnergy'])):
            print(r['name'],': TotalElectronicEnergy is 0 or missing : re-parse the OUTCAR with the current VaspParser.py')
            sys.exit()
    natom = ref.natom
    lines = []
    def out(s=''):
        print(s)
        lines.append(s)

    #<U>_lambda of each run
    out('Thermodynamic integration from the ideal gas : %d runs, %d atoms (%s), -i %d, k = %g' % (len(runs),natom,' '.join('%d %s' % (n,el) for el,n in zip(ref.elements,ref.types)),skip,k))
    out('%-40s %10s %10s %18s %12s %10s %8s %8s' % ('file','lambda','x','<U> (eV)','error','T (K)','V (A^3)','snaps'))
    U, Uerr, temps, Terrs, vols = [], [], [], [], []
    notes = []
    for r in runs:
        mean, std, err, lower, drift = statistics(r['TotalElectronicEnergy'])
        T, Tstd, Terr, Tlower, Tdrift = statistics(r['Temperature'][~np.isnan(r['Temperature'])])
        Terrs.append(Terr)
        V = np.nanmean(r['volume'])
        U.append(mean); Uerr.append(err); temps.append(T); vols.append(V)
        x = 2.0*r['lambda']**(1.0-k)-1.0
        out('%-40s %10.6f %10.6f %18.6f %12.6f%s %10.1f %8.2f %8d' % (os.path.basename(r['name'])[-40:],r['lambda'],x,mean,err,'>' if lower else ' ',T,V,len(r['TotalElectronicEnergy'])))
        if drift > 4.0:
            notes.append('WARNING : '+r['name']+' : <U> drifts (a quarter of the run deviates by %.1f errors) : equilibrate longer or use -i' % drift)
        if lower:
            notes.append('WARNING : '+r['name']+' : fewer than 50 independent samples, the error of <U> is a lower bound')
    U, Uerr, temps, vols = map(np.array,(U,Uerr,temps,vols))
    if np.ptp(vols) > 1e-6*vols.max():
        notes.append('WARNING : the volumes of the runs differ (%.4f to %.4f A^3) : the integration must be at constant volume' % (vols.min(),vols.max()))
    Terrs = np.array(Terrs)
    if np.any(np.abs(temps-temps.mean()) > 4.0*np.sqrt(Terrs**2+1e-12) + 0.002*temps.mean()):
        notes.append('WARNING : the mean temperatures of the runs differ (%.1f to %.1f K) : the integration must be at constant temperature' % (temps.min(),temps.max()))

    #integral over x
    x = 2.0*lambdas**(1.0-k)-1.0
    g = U*lambdas**k/(2.0*(1.0-k))
    gerr = Uerr*lambdas**k/(2.0*(1.0-k))
    gl = gauss_legendre(lambdas,k)
    xs = np.concatenate(([-1.0],x)) if x[0] > -1.0 else x
    gs = np.concatenate(([0.0],g)) if x[0] > -1.0 else g
    es = np.concatenate(([0.0],gerr)) if x[0] > -1.0 else gerr
    results = {}
    for method in ('spline','trapezoid'):
        w = quadrature_weights(xs,method)
        results[method] = (w @ gs, math.sqrt(w**2 @ es**2))
    if gl is not None:
        use, weights = gl
        results['Gauss-Legendre'] = (weights @ g[use], math.sqrt(weights**2 @ gerr[use]**2))
        method = 'Gauss-Legendre'
    else:
        method = 'spline'
    if lambdas.max() < 1.0 - 1e-9:
        notes.append('WARNING : no run at lambda = 1 : E (for S) is missing'+('' if gl is not None else ', and the integral is extrapolated to x = 1'))
    if x[0] > -0.5:
        notes.append('WARNING : the smallest lambda is %g (x = %.2f) : the integrand is assumed to vanish at x = -1 far from the points' % (lambdas[0],x[0]))
    dF, dFerr = results[method]

    #ideal gas, F, E, S at lambda = 1
    one = runs[-1]
    T = Tuser if Tuser is not None else temps.mean()
    V = vols[-1]
    F0 = ideal_gas(ref,V,T)
    F = F0 + dF
    out()
    out('integral of <U> over lambda : %.6f +- %.6f eV (%s ; %s)' % (dF,dFerr,method,' ; '.join('%s %.6f +- %.6f' % (m,v[0],v[1]) for m,v in results.items() if m != method)))
    out('ideal gas free energy at T = %.2f K, V = %.4f A^3 : %.6f eV' % (T,V,F0))
    out('F = %.6f +- %.6f eV per cell = %.6f +- %.6f eV/atom' % (F,dFerr,F/natom,dFerr/natom))
    mass = sum(float(m)*n for m,n in zip(ref.masses,ref.types))                       #g/mol per cell
    out('    %.4f +- %.4f kJ/mol (of atoms) ; %.4f +- %.4f kJ/g' % (F/natom*EV*NA/1000.0,dFerr/natom*EV*NA/1000.0,F*EV*NA/mass/1000.0,dFerr*EV*NA/mass/1000.0))
    if lambdas[-1] == 1.0:
        IE, K = one['InternalEnergy'], one['KineticEnergyIons']
        included, how = kinetic_included(IE,np.nan_to_num(K),np.nan_to_num(one['TotalElectronicEnergy']),np.nan_to_num(one['ElectronicEntropy']))
        E = IE if included else IE + np.nan_to_num(K)
        Emean, Estd, Eerr, Elower, Edrift = statistics(E)
        S = (Emean - F)/T
        Serr = math.sqrt(Eerr**2 + dFerr**2)/T
        out('E (lambda = 1, InternalEnergy %s the kinetic energy : %s) = %.6f +- %.6f eV' % ('with' if included else '+',how,Emean,Eerr))
        out('S = (E - F)/T = %.6e +- %.6e eV/K per cell = %.4f +- %.4f kB/atom' % (S,Serr,S/(natom*KB),Serr/(natom*KB)))
        out('    %.4f +- %.4f J/K/mol (of atoms) ; %.4f +- %.4f J/K/g' % (S/natom*EV*NA,Serr/natom*EV*NA,S*EV*NA/mass,Serr*EV*NA/mass))
    for note in notes:
        out(note)
    with open('thermo_integration.dat','w') as f:
        f.write('\n'.join(lines)+'\n')
    print('Written in thermo_integration.dat')

    if visual:
        import matplotlib.pyplot as plt
        fig, (a1,a2) = plt.subplots(1,2,figsize=(11,4.5))
        a1.errorbar(lambdas,U,yerr=Uerr,fmt='o')
        a1.set_xlabel('lambda'); a1.set_ylabel('<U> (eV)'); a1.set_xscale('log')
        xx = np.linspace(-1,1,400)
        a2.plot(xx,CubicSpline(xs,gs,bc_type='not-a-knot' if len(xs) > 3 else 'natural')(xx),'-',label='cubic spline')
        a2.errorbar(xs,gs,yerr=es,fmt='o',label='<U> lambda^k / 2(1-k)')
        a2.set_xlabel('x = 2 lambda^(1-k) - 1'); a2.set_ylabel('integrand (eV)'); a2.legend()
        fig.tight_layout()
        fig.savefig('thermo_integration.pdf')
        print('Plot in thermo_integration.pdf')
        plt.show()


if __name__ == "__main__":
    main(sys.argv[1:])
