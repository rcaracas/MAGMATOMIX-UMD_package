#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### Authors : Razvan Caracas
### Harmonic thermodynamics from the vibrational densities of states written by vibrations_umd.py :
### zero-point energy, energy, free energy, entropy and heat capacity, quantum and classical, and the quantum corrections to classical MD
### Entropy of liquids with the two-phase thermodynamic (2PT) model : Lin, Blanco, Goddard, J. Chem. Phys. 119, 11792 (2003)


import sys, getopt, os.path, re
import numpy as np
from scipy.optimize import brentq
import umd_processes_fast as umdpf


KB = 8.617333262E-05#Boltzmann constant in eV/K
HC = 1.239841984E-04#h c in eV.cm : energy of 1 cm^-1
R = 8.314462618#gas constant in J/(mol K)
KB_SI = 1.380649E-23#J/K
H_SI = 6.62607015E-34#J.s
AMU = 1.66053906660E-27#kg
C_CM = 2.99792458E10#speed of light in cm/s


def read_umd_header(umdfile):#Elements, number of atoms and masses (amu) of each element, from the header of the UMD file
    elements = None
    types = None
    masses = None
    with open(umdfile,'r') as f:
        for line in f:
            entry = line.split()
            if len(entry) > 0 and entry[0] == 'elements':
                elements = entry[1:]
            elif len(entry) > 0 and entry[0] == 'types':
                types = [int(x) for x in entry[1:]]
            elif len(entry) > 0 and entry[0] == 'masses':
                masses = [float(x) for x in entry[1:]]
            if elements is not None and types is not None and masses is not None:
                break
    return elements, np.array(types), np.array(masses)


def mean_volume(umdfile):#Mean volume of the cell (A^3) over the snapshots, from the rprimd lines
    volumes = []
    cell = []
    with open(umdfile,'r') as f:
        for line in f:
            if line.startswith('rprimd_'):
                entry = line.split()
                cell.append([float(entry[1]),float(entry[2]),float(entry[3])])
                if len(cell) == 3:
                    volumes.append(abs(np.linalg.det(np.array(cell))))
                    cell = []
    return np.mean(volumes)


def read_vdos(vibrfile,elements):#Frequencies (cm^-1), vDOS of each element and their errors (None without them), and the comment lines
    comments = []
    lines = []
    with open(vibrfile,'r') as f:
        for line in f:
            if line.startswith('#'):
                comments.append(line[1:].strip())
            elif line.strip() != '':
                lines.append(line.split())
    labels = lines[0]
    data = np.array(lines[1:],dtype=float)
    g = np.array([data[:,labels.index('VDos_'+el)] for el in elements])
    gerr = None
    if all('VDos_'+el+'_err' in labels for el in elements):
        gerr = np.array([data[:,labels.index('VDos_'+el+'_err')] for el in elements])
    return data[:,0], g, gerr, comments


def integrate(nu,g,f,f0):#int_0^numax g(nu) f(nu) dnu, trapezoidal ; f0 = int_0^nu1 f(nu) dnu / nu1 for the first interval, where f may diverge (log) at nu = 0
    dnu = nu[1]-nu[0]
    first = 0.5*(g[...,0]+g[...,1])*f0*dnu
    rest = dnu*(0.5*g[...,1]*f[1]+(g[...,2:-1]*f[2:-1]).sum(axis=-1)+0.5*g[...,-1]*f[-1])
    return first+rest


def weights(nu,T):#Functions f(nu) whose integrals with the vDOS give the thermodynamics at T, and f0, their mean over the first interval (0, nu1)
    #quantum and classical E, F (eV), S, Cv (kB), zero-point energy (eV), and the quantum corrections dE, dF, dS, dCv
    x = HC*nu/(KB*T)#x = h nu/(kB T)
    x[0] = 1.0#placeholder at nu = 0, replaced by the analytic first interval
    half = 0.5*x
    coth = 1.0/np.tanh(half)
    lnsinh = np.where(x < 50.0,np.log(2.0*np.sinh(np.minimum(half,25.0))),half)#ln(2 sinh(x/2)), = x/2 for large x
    emx = np.exp(-np.minimum(x,700.0))
    x1 = x[1]
    #functions of x and their mean over the first interval (0, nu1), from their expansions at small x
    q = {}
    q['zpe'] = (HC*0.5*nu, HC*0.25*nu[1])
    q['E'] = (KB*T*half*coth, KB*T*(1.0+x1**2/36.0))#(h nu/2) coth(h nu/2kT) -> kT
    q['F'] = (KB*T*lnsinh, KB*T*(np.log(x1)-1.0+x1**2/72.0))#kT ln(2 sinh(x/2)) -> kT ln x
    q['S'] = (half*coth-lnsinh, 2.0-np.log(x1)+x1**2/72.0)#in kB : x/2 coth(x/2) - ln(2 sinh(x/2)) -> 1 - ln x + x^2/24
    q['Cv'] = (x**2*emx/(1.0-emx)**2, 1.0-x1**2/36.0)#in kB : x^2 e^x/(e^x-1)^2 -> 1
    q['E_cl'] = (KB*T*np.ones_like(x), KB*T)
    q['F_cl'] = (KB*T*np.log(x), KB*T*(np.log(x1)-1.0))
    q['S_cl'] = (1.0-np.log(x), 2.0-np.log(x1))
    q['Cv_cl'] = (np.ones_like(x), 1.0)
    for key in ('E','F','S','Cv'):#quantum corrections to classical MD (Berens et al. 1983)
        q['d'+key] = (q[key][0]-q[key+'_cl'][0], q[key][1]-q[key+'_cl'][1])
    return q


def harmonic(nu,g,T):#Integrals of the vDOS g(nu) (states per cm^-1) with the functions of weights
    return {key:integrate(nu,g,f,f0) for key,(f,f0) in weights(nu,T).items()}


def harmonic_errors(nu,gerr,T,ncorr):#Errors of the integrals : the errors of the vDOS are taken independent on a grid of ncorr points (one FWHM)
    dnu = nu[1]-nu[0]
    err = {}
    for key,(f,f0) in weights(nu,T).items():
        f = f.copy()
        f[0] = f0
        err[key] = np.sqrt(((f*gerr*dnu)**2).sum(axis=-1)*ncorr)
    return err


def fluidicity(delta):#Fraction f of gas-like modes, root in (0,1) of the hard-sphere equation of 2PT (Lin et al. 2003, eq. 35)
    func = lambda f: 2.0*delta**-4.5*f**7.5-6.0*delta**-3*f**5-delta**-1.5*f**3.5+6.0*delta**-1.5*f**2.5+2.0*f-2.0
    if func(1.0) <= 0.0:
        return 1.0
    return brentq(func,1.0E-12,1.0,xtol=1.0E-14)


def two_phase(nu,g,n,mass,volume,T):#2PT decomposition of the vDOS g (states per cm^-1, normalized to 3n) of n atoms of mass (amu) in volume (A^3) at T
    #returns the gas-like vDOS, f, Delta, the packing fraction y and the entropy per gas-like atom of the hard-sphere fluid (kB)
    s0 = g[0]
    m = mass*AMU
    delta = 2.0*(s0/C_CM)/(9.0*n)*np.sqrt(np.pi*KB_SI*T/m)*(n/(volume*1.0E-30))**(1.0/3.0)*(6.0/np.pi)**(2.0/3.0)#normalized diffusivity
    if delta <= 0.0:
        return np.zeros_like(g), 0.0, 0.0, 0.0, 0.0
    f = fluidicity(delta)
    ggas = s0/(1.0+(np.pi*s0*nu/(6.0*f*n))**2)#hard-sphere gas : Lorentzian with the same g(0), integrating to 3 f n
    y = f**2.5/delta**1.5#packing fraction of the hard spheres
    z = (1.0+y+y**2-y**3)/(1.0-y)**3#Carnahan-Starling compressibility
    lam3 = (2.0*np.pi*m*KB_SI*T/H_SI**2)**1.5#1/Lambda^3, in m^-3
    shs = 2.5+np.log(lam3*volume*1.0E-30/(f*n)*z)+y*(3.0*y-4.0)/(1.0-y)**2#per gas-like atom, with the ideal mixing through V/(f n)
    return ggas, f, delta, y, shs


def main(argv):
    umdfile = ''
    temperatures = None
    usage = 'vdos_thermo.py -f <UMD_filename> -t <Temperatures>'
    try:
        opts, arg = getopt.getopt(argv,"hf:t:",["fumdfile","tTemperatures"])
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('Program to compute the harmonic thermodynamics from the vibrational densities of states (vDOS) of vibrations_umd.py. Usage: ')
            print(usage)
            print('It reads <base>.vibr.err.dat (or <base>.vibr.dat, without errors), written by vibrations_umd.py from <base>.umd.dat, and')
            print('the number of atoms of each element from the header of the UMD file. The vDOS of each element is normalized to 3 N_El modes.')
            print('For each element and for all the atoms, per atom :')
            print('  quantum harmonic : zero-point energy ZPE, energy E (with the ZPE), free energy F, entropy S, heat capacity Cv')
            print('  classical harmonic : E_cl = 3 kB T, F_cl, S_cl, Cv_cl = 3 kB')
            print('  quantum corrections dE = E-E_cl, dF = F-F_cl, dS = S-S_cl, dCv = Cv-Cv_cl : add them to the classical values of the MD')
            print('  (Berens et al., J. Chem. Phys. 79, 2375, 1983), e.g. to the energy E and the Cv of averages.py -a.')
            print('The frequencies are those of the MD at its temperature, so they include the anharmonicity at that temperature ; at other')
            print('temperatures -t the frequencies are kept (quasi-harmonic only for solids, and only near the MD temperature).')
            print('Liquids : the diffusive modes give a vDOS > 0 at zero frequency, for which the harmonic S and F are not defined : a warning')
            print('is printed. The quantum corrections and Cv remain valid, since the quantum and classical oscillators are the same at low')
            print('frequency. The entropy of liquids is given by the two-phase thermodynamic model (2PT, Lin, Blanco, Goddard, J. Chem. Phys.')
            print('119, 11792, 2003), at the temperature of the vDOS : the vDOS of each element is split into a gas-like part, a hard-sphere')
            print('fluid with the same vDOS at zero frequency (Lorentzian, 3 f N_El modes ; the fluidicity f comes from the normalized')
            print('diffusivity Delta), and a solid-like part, quantum harmonic oscillators. S = S_gas (Carnahan-Starling hard spheres) + S_solid.')
            print('Each element is a 2PT system of its own, with its partial density N_El/V, so that S_gas includes the ideal entropy of mixing')
            print('(an approximation for mixtures : the cross-correlations of the velocities are neglected). For a solid f = 0 and S = S_harmonic.')
            print('The gas-like Lorentzian decays as 1/nu^2, so that it exceeds the vDOS at high frequency, where the solid-like vDOS becomes')
            print('negative (column negative, ~2 % of the modes for liquid argon) : the quantum S_2PT is barely affected, but a classical S_2PT')
            print('and quantum corrections of the solid-like modes would be ; the quantum corrections are those of the whole vDOS.')
            print('Test : Lennard-Jones liquid argon at 87.3 K and 1.395 g/cm^3 gives S_2PT = 56.2 J/(mol K) ; experiment ~55.6.')
            print('-f : UMD file from which vibrations_umd.py computed the vDOS.')
            print('-f all : every *.umd.dat file of the current folder, one after the other.')
            print('-t : temperatures in K, separated by commas, e.g. -t 300,1000,1500. Default is the temperature of the vDOS.')
            print('Output : <base>.vdos.thermo.dat : the quantities at each temperature, per atom, with their statistical errors (from the')
            print('         error band of the vDOS) when <base>.vibr.err.dat exists.')
            print('         <base>.2pt.dat : the 2PT entropies of each element and of all the atoms, per atom, at the temperature of the vDOS.')
            sys.exit()
        elif opt in ("-f","--fumdfile"):
            umdfile = str(arg)
        elif opt in ("-t","--tTemperatures"):
            try:
                temperatures = [float(x) for x in arg.split(',')]
            except ValueError:
                print('ERROR : the temperatures -t have to be numbers separated by commas, e.g. -t 300,1000.')
                sys.exit()
            if min(temperatures) <= 0.0:
                print('ERROR : the temperatures -t have to be strictly positive.')
                sys.exit()

    if not os.path.isfile(umdfile):
        umdpf.missing_umd(umdfile)
        sys.exit()
    base = umdfile[:-8] if umdfile.endswith('.umd.dat') else os.path.splitext(umdfile)[0]
    vibrfile = base+'.vibr.err.dat'
    if not os.path.isfile(vibrfile):
        vibrfile = base+'.vibr.dat'
        if not os.path.isfile(vibrfile):
            print('ERROR : neither ',base+'.vibr.err.dat',' nor ',vibrfile,' exist : run vibrations_umd.py on ',umdfile,' first.')
            sys.exit()
        print('No error file : the quantities are given without errors. Run vibrations_umd.py again for <base>.vibr.err.dat.')
    elements, counts, masses = read_umd_header(umdfile)
    nu, g, gerr, comments = read_vdos(vibrfile,elements)
    Tvdos = None
    fwhm = None
    for comment in comments:
        found = re.search(r'T = ([0-9.]+) K',comment)
        if found:
            Tvdos = float(found.group(1))
        found = re.search(r'resolution \(FWHM\) ([0-9.]+) cm',comment)
        if found:
            fwhm = float(found.group(1))
    if Tvdos is None:#vibr.dat has no comments : the mean temperature of the UMD file, as used by vibrations_umd.py by default
        T = [float(line.split()[1]) for line in open(umdfile) if line.startswith('Temperature')]
        Tvdos = np.mean(T)
    if temperatures is None:
        temperatures = [Tvdos]

    #Normalization to 3 N_El modes : vibrations_umd.py normalizes to 3 N_El T_kinetic/T
    dnu = nu[1]-nu[0]
    norms = (g.sum(axis=-1)-0.5*(g[:,0]+g[:,-1]))*dnu
    scale = 3.0*counts/norms
    g = g*scale[:,None]
    if gerr is not None:
        gerr = gerr*scale[:,None]
    print('vDOS read from ',vibrfile,' (computed at ','%.1f' % Tvdos,' K)')
    print('Normalized to 3 N_El modes : '+', '.join(el+' x'+'%.4f' % scale[ii] for ii,el in enumerate(elements)))
    negative = -(np.minimum(g,0.0).sum(axis=-1))*dnu/(3.0*counts)
    if negative.max() > 0.01:
        print('WARNING : the vDOS has negative parts (up to ','%.1f' % (100*negative.max()),'% of the modes) : use the Hann taper in vibrations_umd.py.')
    g0 = g[:,0]/(3.0*counts)
    peak = g.max(axis=-1)/(3.0*counts)
    diffusive = [el for ii,el in enumerate(elements) if g0[ii] > 0.02*peak[ii]]
    if diffusive != []:
        print('WARNING : the vDOS of ',', '.join(diffusive),' is not zero at zero frequency : diffusive modes (liquid, or a mobile atom).')
        print('          The harmonic S and F of these elements, and of all the atoms, are not meaningful ; dE, dF, dS, dCv, E and Cv are.')
        print('          Their entropy is given by the 2PT model at the end.')

    rows = []
    for T in temperatures:
        res = harmonic(nu,g,T)
        err = None
        if gerr is not None:#the errors of the vDOS are correlated over the resolution
            err = harmonic_errors(nu,gerr,T,max(1.0,(fwhm if fwhm is not None else dnu)/dnu))
        quantities = ['zpe','E','F','S','Cv','dE','dF','dS','dCv']
        units = {'zpe':1000.0,'E':1000.0,'F':1000.0,'dE':1000.0,'dF':1000.0,'S':1.0,'Cv':1.0,'dS':1.0,'dCv':1.0}
        print('\nT = ','%.1f' % T,' K ; per atom : ZPE, E, F, dE, dF in meV ; S, Cv, dS, dCv in kB (x R = ','%.3f' % R,' for J/(mol K) per mole of atoms)')
        print('{:<6s}{:>6s}'.format('','N')+''.join('{:>10s}'.format(q.upper() if q == 'zpe' else q) for q in quantities))
        total = {key:res[key].sum() for key in res}
        for ii,el in enumerate(list(elements)+['all']):
            n = counts[ii] if ii < len(elements) else counts.sum()
            vals = [(res[q][ii] if ii < len(elements) else total[q])/n*units[q] for q in quantities]
            print('{:<6s}{:>6d}'.format(el,int(n))+''.join('{:>10.3f}'.format(v) for v in vals))
            row = [T,el,int(n)]+vals
            if err is not None:
                errs = [(err[q][ii] if ii < len(elements) else np.sqrt((err[q]**2).sum()))/n*units[q] for q in quantities]
                row += errs
            rows.append(row)
        if err is not None:
            relerr = max(np.sqrt((err[q]**2).sum())/abs(total[q]) for q in ['E','S','Cv'] if total[q] != 0.0)
            print('statistical errors (from the error band of the vDOS) : at most ','%.2g' % (100*relerr),'% of E, S, Cv for all the atoms ; see the output file')
        print('cell : E = ','%.4f' % total['E'],' eV, F = ','%.4f' % total['F'],' eV, S = ','%.3f' % total['S'],' kB, Cv = ','%.3f' % total['Cv'],
              ' kB ; quantum corrections dE = ','%.4f' % total['dE'],' eV, dF = ','%.4f' % total['dF'],' eV, dCv = ','%.3f' % total['dCv'],' kB')

    outfile = base+'.vdos.thermo.dat'
    with open(outfile,'w') as f:
        f.write('# Harmonic thermodynamics from the vDOS of '+vibrfile+', computed at '+'%.1f' % Tvdos+' K, normalized to 3 N_El modes\n')
        f.write('# per atom : ZPE, E (with ZPE), F, dE, dF in meV ; S, Cv, dS, dCv in kB ; d = quantum - classical harmonic, to add to classical MD values\n')
        if diffusive != []:
            f.write('# diffusive modes (vDOS > 0 at zero frequency) for '+', '.join(diffusive)+' : their S and F, and those of all, are not meaningful\n')
        labels = ['T(K)','element','N','ZPE','E','F','S','Cv','dE','dF','dS','dCv']
        if gerr is not None:
            labels += [q+'_err' for q in ['ZPE','E','F','S','Cv','dE','dF','dS','dCv']]
        f.write('\t'.join(labels)+'\n')
        for row in rows:
            f.write('\t'.join(str(x) for x in row)+'\n')
    print('\nThermodynamics printed in file ',outfile)

    #Two-phase thermodynamic model, at the temperature of the vDOS
    volume = mean_volume(umdfile)
    w = weights(nu,Tvdos)
    labels = ['N','Delta','f','y','S_gas','S_solid','S_2PT','negative']
    table = []
    for ii in range(len(elements)):
        ggas, f, delta, y, shs = two_phase(nu,g[ii],counts[ii],masses[ii],volume,Tvdos)
        gsolid = g[ii]-ggas
        sgas = integrate(nu,ggas,np.full_like(nu,shs/3.0),shs/3.0)#each of the 3 f N gas-like modes carries S_HS/3
        ssolid = integrate(nu,gsolid,*w['S'])
        negative = -np.minimum(gsolid,0.0).sum()*dnu/(3.0*counts[ii])#where the Lorentzian tail of the gas-like vDOS exceeds the vDOS
        table.append([counts[ii],delta,f,y,sgas,ssolid,sgas+ssolid,negative])
    table = np.array(table)
    print('\nTwo-phase thermodynamic model (2PT) at ','%.1f' % Tvdos,' K, volume ','%.2f' % volume,' A^3 ; per atom, S in kB')
    print('S_2PT = S_gas (hard spheres) + S_solid (quantum oscillators) ; negative : fraction of the modes where the solid-like vDOS < 0')
    print('{:<6s}{:>6s}'.format('','N')+''.join('{:>10s}'.format(lab) for lab in labels[1:]))
    rows2pt = []
    for ii,el in enumerate(list(elements)+['all']):
        if ii < len(elements):
            n = counts[ii]
            vals = [table[ii,1],table[ii,2],table[ii,3]]+list(table[ii,4:7]/n)+[table[ii,7]]
        else:#all the atoms : f is the fraction of gas-like modes, Delta and y have no meaning
            n = counts.sum()
            vals = [np.nan,(table[:,2]*counts).sum()/n,np.nan]+list(table[:,4:7].sum(axis=0)/n)+[(table[:,7]*counts).sum()/n]
        print('{:<6s}{:>6d}'.format(el,int(n))+''.join('{:>10.4f}'.format(v) if jj in (0,1,2,6) else '{:>10.3f}'.format(v) for jj,v in enumerate(vals)))
        rows2pt.append([el,int(n)]+vals)
    stot = table[:,6].sum()
    print('cell : S_2PT = ','%.3f' % stot,' kB = ','%.4f' % (stot*KB),' eV/K ; T S_2PT = ','%.4f' % (Tvdos*stot*KB),' eV')
    print('       F = E_MD + dE - T S_2PT, with E_MD the mean energy of the MD (e.g. from averages.py) and dE the quantum correction above.')
    print('       The solid-like vDOS is negative where the 1/nu^2 tail of the gas-like Lorentzian exceeds the vDOS : this barely changes')
    print('       S_2PT (the quantum entropy of high frequencies is small), but makes a classical S_2PT and quantum corrections of the')
    print('       solid-like modes unreliable : the quantum corrections dE, dS, dCv are those of the whole vDOS, in the table above.')
    outfile2pt = base+'.2pt.dat'
    with open(outfile2pt,'w') as f:
        f.write('# Two-phase thermodynamic model (Lin, Blanco, Goddard 2003) from the vDOS of '+vibrfile+' at '+'%.1f' % Tvdos+' K, volume '+'%.4f' % volume+' A^3\n')
        f.write('# per atom : S in kB ; each element is a 2PT system with its partial density ; for all : f = fraction of gas-like modes\n')
        f.write('# negative : fraction of the modes where the solid-like vDOS (vDOS - gas-like Lorentzian) is negative\n')
        f.write('\t'.join(['element']+labels)+'\n')
        for row in rows2pt:
            f.write('\t'.join(str(x) for x in row)+'\n')
    print('2PT entropies printed in file ',outfile2pt)


if __name__ == "__main__":
    if not umdpf.run_on_all(sys.argv[1:],main):#-f all : every UMD file of the folder
        main(sys.argv[1:])
