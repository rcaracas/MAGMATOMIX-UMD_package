#!/usr/bin/env python3
###
#AUTHORS: ANAIS KOBSCH, RAZVAN CARACAS
###
#Quality checks of a run from its UMD file : integrity of the file (times, unwrapping, jumps), integration (energy drift, time step, momentum,
#net forces), thermodynamic state (drifts, temperature, cell, stress), electronic structure (overlap of the core spheres, smearing) and charge
#balance of the composition. Each check ends with OK, WARNING or skipped. The diffusion, viscosity, vibrations, g(r) and speciation keep their own
#checks in their scripts. Formerly check_overlap.py (the overlap of the core spheres only).

import sys, getopt, os.path, re
import numpy as np
import umd_processes_fast as umdpf
import crystallography as cr
from averages import statistics
from vibrations_umd import KB, AMU_A2FS2

FIELDS = ('timestep','time','InternalEnergy','ElectronicEntropy','KineticEnergyIons','EnergyWithDrift','Temperature','Pressure','StressTensor',
          'rprimd_a','rprimd_b','rprimd_c')


class Report:#Lines printed and kept for the report file ; the warnings are listed at the end
    def __init__(self):
        self.lines = []
        self.warnings = []
        self.number = 0

    def write(self,text=''):
        print(text)
        self.lines.append(text)

    def check(self,title):
        self.number += 1
        self.current = str(self.number)+'. '+title
        self.write('\n'+self.current)

    def result(self,status,text=''):
        self.write('   -> '+status+(' : '+text if text else ''))
        if status == 'WARNING':
            self.warnings.append(self.current+' : '+text)


#****Reading

def read_radius(radiusfile):#Lines 'element radius unit' (unit A or Bohr) ; radii in A
    AngtoBohr = 1.8897259886
    radius = {}
    with open(radiusfile,'r') as r:
        for line in r:
            entry = line.split()
            if len(entry) < 3:
                continue
            radius[entry[0]] = float(entry[1])/AngtoBohr if re.match('[Bb]',entry[2]) else float(entry[1])
    return radius


#****Geometry

def lens_volume(d,r1,r2):#Volume of the intersection of two spheres of radii r1, r2 whose centres are d apart
    if d >= r1+r2:
        return 0.0
    if d <= abs(r1-r2):
        return 4.0/3.0*np.pi*min(r1,r2)**3
    return np.pi*(r1+r2-d)**2*(d**2+2*d*(r1+r2)-3*(r1-r2)**2)/(12*d)


def shortest_distances(xred,cell,masks):#Shortest distance between the atoms of each pair of elements, nearest images (exact up to half the smallest width of the cell)
    d = xred[:,None,:]-xred[None,:,:]
    d -= np.rint(d)
    r = np.linalg.norm(d @ cell,axis=2)
    np.fill_diagonal(r,np.inf)
    return np.array([[r[np.ix_(mi,mj)].min() if mi.any() and mj.any() else np.inf for mj in masks] for mi in masks])


def mean_with_error(x,unit,fmt='%.4g'):
    mean, std, err, lower, drift = statistics(np.asarray(x,dtype=float))
    return mean, err, drift, (fmt % mean)+' +/- '+('>' if lower else '')+('%.2g' % err)+' '+unit


#****Checks

def main(argv):
    umdfile = ''
    radiusfile = ''
    first = 0
    sampling = 1
    level = 12.0
    charges = {}
    umdpf.headerumd()
    usage = 'check_umd.py -f <UMD file> -i <InitialStep> -s <Sampling> -r <radius file> -l <level_of_overlap_%> -q <formal charges, e.g. Fe:3,C:-4>'
    try:
        opts, arg = getopt.getopt(argv,"hf:i:s:r:l:q:")
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('check_umd.py : quality checks of a run from its UMD file, each ending with OK, WARNING or skipped :')
            print('  file      1. times (constant spacing, no duplicated, missing or backward snapshots) and repeated configurations')
            print('            2. unwrapping : the unwrapped minus the wrapped positions are whole cell vectors (to 1e-3 in reduced coordinates)')
            print('            3. jumps : no atom moves more than 2 max(|v|) dt + 0.2 A between two snapshots read')
            print('  dynamics  4. energy drift : slope of EnergyWithDrift (conserved with the thermostat) < 1 meV/atom/ps, and its jumps')
            print('            5. time step : the fastest atoms move less than 0.1 A per time step')
            print('            6. momentum : centre-of-mass velocity and displacement, compared with the displacements of the atoms (MSD)')
            print('            7. net force : |sum of the forces| < 0.05 x the rms force of an atom')
            print('  state     8. equilibration : no quarter of the run deviating from the mean T, P, E by more than 4 errors (as averages.py -a)')
            print('            9. temperature : the kinetic temperature of the velocities agrees with the Temperature field (5 %) ; zero velocities')
            print('           10. cell (NPT) : drift of the volume and of the lattice lengths')
            print('           11. stress : mean off-diagonal stresses 0 and diagonal ones equal within 4 errors (hydrostatic : liquids)')
            print('  electrons 12. overlap of the core spheres (with -r) : below -l % of the volume of the two spheres')
            print('           13. smearing : electronic entropy term T S_el < 1 meV/atom')
            print('  charges  14. formal charge balance of the composition (crystallography.Elements2charge ; -q to change)')
            print(usage)
            print('-f : the UMD file.')
            print('-f all : every *.umd.dat file of the current folder, one after the other.')
            print('-i : snapshots skipped at the beginning, for the checks of the atoms (2, 3, 5-7, 9, 12) and of the state (4, 8-11, 13). Default 0.')
            print('-s : every s-th snapshot for the checks of the atoms. Default 1 (all : the jumps are checked between consecutive snapshots).')
            print('-r : radius file for check 12 : one line per element, "element radius unit" (Angstrom or Bohr), e.g. "Si 1.9 Bohr" (RCORE of')
            print('     a VASP POTCAR, or the pseudopotential / PAW radius of another code). Without it, check 12 is skipped.')
            print('-l : acceptable overlap of the core spheres, in % of the sum of their volumes. Default 12.')
            print('-q : formal charges replacing those of the element table, e.g. -q Fe:3 or -q C:-4,N:-3 ; fractions allowed (Fe:2.5).')
            print('Output : the report, also in <base>.check.dat ; with -r, <base>.overlap.dat (one line per snapshot : OK/NOT-ok, then for each pair')
            print('of elements the shortest distance (A) and the overlap (%)).')
            sys.exit()
        elif opt == '-f':
            umdfile = arg
        elif opt == '-i':
            first = int(arg)
        elif opt == '-s':
            sampling = int(arg)
        elif opt == '-r':
            radiusfile = arg
        elif opt == '-l':
            level = float(arg)
        elif opt == '-q':
            for item in arg.split(','):
                el, q = item.split(':')
                charges[el] = float(q)
    if not os.path.isfile(umdfile):
        umdpf.missing_umd(umdfile)
        sys.exit()
    if radiusfile and not os.path.isfile(radiusfile):
        print('the radius file',radiusfile,'does not exist')
        sys.exit()
    base = umdfile[:-8] if umdfile.endswith('.umd.dat') else os.path.splitext(umdfile)[0]
    umdpf.remove_former_outputs([base+'.check.dat',base+'.overlap.dat'])#of a former run

    MyCrystal,TimeStep = umdpf.Crystallization(umdfile)
    natom = MyCrystal.natom
    elements = MyCrystal.elements
    typat = np.array(MyCrystal.typat)
    masses = np.array(MyCrystal.masses)[typat]
    masks = [typat == it for it in range(MyCrystal.ntypat)]
    radius = read_radius(radiusfile) if radiusfile else None
    if radius is not None:
        missing = [el for el in elements if el not in radius]
        if missing:
            print('no radius for',', '.join(missing),'in',radiusfile)
            sys.exit()
        pairs = [(i,j) for i in range(MyCrystal.ntypat) for j in range(i,MyCrystal.ntypat)]
        labels = [elements[i]+'-'+elements[j] for i,j in pairs]
        volumes = [4.0/3.0*np.pi*(radius[elements[i]]**3+radius[elements[j]]**3) for i,j in pairs]
        ov = open(base+'.overlap.dat','w')
        ov.write('# Overlap of the core spheres in '+umdfile+', every '+str(sampling)+' snapshots from '+str(first)+' ; acceptable overlap '+str(level)+' %\n')
        ov.write('# radii (A) : '+' '.join(el+' '+'%.4f' % radius[el] for el in elements)+'\n')
        ov.write('# overlap = volume of the intersection of the two spheres at the shortest distance, in % of the sum of their volumes\n')
        ov.write('snapshot\tstatus\t'+'\t'.join(lab+'_dmin(A)\t'+lab+'_overlap(%)' for lab in labels)+'\n')

    #****One pass over the file
    series = {k: [] for k in FIELDS}                                   #fields of every snapshot
    cells = []                                                          #cell of every snapshot (None without rprimd lines)
    steps = []
    repeated = 0
    unwrap = []                                                         #per snapshot read : largest non-integer part of (unwrapped - wrapped) in reduced coordinates
    jumps = []                                                          #(snapshot, atom, displacement, bound)
    maxdisp = np.zeros(MyCrystal.ntypat)
    speeds = []                                                         #per snapshot read : largest speed of each element (A/fs)
    Pcom = []
    Rcom = []
    netforce = []
    rmsforce = []
    Tvel = []
    Tfield = []
    zero = []
    dmins = []
    percents = []
    read = []
    first_abs = previous = None
    for isnap, fields, atoms in umdpf.read_snapshots(umdfile,natom,first,sampling,FIELDS):
        steps.append(isnap)
        for k in FIELDS:
            v = fields.get(k)
            series[k].append(v if k in ('StressTensor','rprimd_a','rprimd_b','rprimd_c') else (v[0] if v else np.nan))
        cell = umdpf.cell_of(fields)
        cells.append(cell)
        if atoms is None:
            continue
        read.append(isnap)
        xred, xcart, absx, vel = atoms[:,0:3], atoms[:,3:6], atoms[:,6:9], atoms[:,9:12]
        speed = np.linalg.norm(vel,axis=1)
        speeds.append([speed[m].max() for m in masks])
        if previous is not None:
            pstep, pabs, pspeed, ptime, pcart = previous
            frac = ((absx-pabs)-(xcart-pcart)) @ np.linalg.inv(cell)    #the unwrapped positions move as the wrapped ones, plus whole cell vectors at the crossings
            unwrap.append(np.abs(frac-np.rint(frac)).max())
            if np.array_equal(absx,pabs):
                repeated += 1
            dt = fields['time'][0]-ptime if 'time' in fields else TimeStep*sampling
            disp = np.linalg.norm(absx-pabs,axis=1)
            bound = 2.0*np.maximum(speed,pspeed)*abs(dt)+0.2
            for ia in np.nonzero(disp > bound)[0]:
                jumps.append((isnap,ia,disp[ia],bound[ia]))
            maxdisp = np.maximum(maxdisp,[disp[m].max() for m in masks])
        else:
            first_abs = absx.copy()
        previous = (isnap,absx,speed,fields['time'][0] if 'time' in fields else isnap*TimeStep,xcart)
        Pcom.append((masses[:,None]*vel).sum(axis=0))
        Rcom.append((masses[:,None]*absx).sum(axis=0)/masses.sum())
        if atoms.shape[1] >= 15:
            forces = atoms[:,12:15]
            netforce.append(np.linalg.norm(forces.sum(axis=0)))
            rmsforce.append(np.sqrt((forces**2).sum(axis=1).mean()))
        ekin = 0.5*(masses[:,None]*vel**2).sum()*AMU_A2FS2
        Tvel.append(2.0*ekin/(3.0*(natom-1)*KB))
        Tfield.append(fields['Temperature'][0] if 'Temperature' in fields else np.nan)
        zero.append(not vel.any())
        if radius is not None:
            dm = shortest_distances(xred % 1.0,cell,masks)
            d = [dm[i,j] for i,j in pairs]
            pc = [100.0*lens_volume(dd,radius[elements[i]],radius[elements[j]])/v if np.isfinite(dd) else 0.0 for dd,(i,j),v in zip(d,pairs,volumes)]
            dmins.append(d)
            percents.append(pc)
            ov.write(str(isnap)+'\t'+('NOT-ok' if max(pc) > level else 'OK')+'\t'+'\t'.join('%.4f\t%.3f' % (dd,p) for dd,p in zip(d,pc))+'\n')
    last_abs = previous[1] if previous is not None else None
    if radius is not None:
        ov.close()
    nsnap = len(steps)
    run = np.array(steps) >= first                                      #snapshots of the checks of the state
    times = np.array(series['time'],dtype=float)
    natoms_read = len(read)
    if natoms_read == 0:
        print('no snapshot read from',first,'every',sampling)
        sys.exit()

    Tvel = np.array(Tvel)
    Tfield = np.array(Tfield)
    ratio = Tvel/np.where(Tfield > 0,Tfield,np.nan)
    faulty = np.array(zero) | (ratio > 2.0) | (ratio < 0.5)           #velocities 0, or off by more than 2 (e.g. positions written as velocities by old parsers)
    good = ~faulty & ~np.isnan(ratio)
    novel = all(zero)                                                   #no velocities in the file (e.g. a LAMMPS dump without vx vy vz)
    NOVEL = 'all the velocities are 0 : not in the file (e.g. a LAMMPS dump without vx vy vz)'

    rep = Report()
    rep.write('Checks of '+umdfile+' : '+str(natom)+' atoms ('+', '.join('%d %s' % (n,el) for el,n in zip(elements,MyCrystal.types))+'), '+str(nsnap)
              +' snapshots, timestep '+str(TimeStep)+' fs ; atoms checked in '+str(natoms_read)+' snapshots (every '+str(sampling)+' from '+str(first)+')')

    #1. Times
    rep.check('Times and snapshots')
    dts = np.diff(times)
    ts = np.array(series['timestep'],dtype=float)
    nback = int((dts < 0).sum())
    ndup = int((dts == 0).sum())
    ngap = int((dts > 1.5*ts[1:]).sum())
    nirreg = int((np.abs(dts-ts[1:]) > 1e-3*ts[1:]).sum())-ngap-ndup-nback
    varying = np.ptp(ts) > 1e-6*np.abs(ts).max()
    rep.write('   time from '+str(times[0])+' to '+str(times[-1])+' fs ; timestep field '+('from '+str(ts.min())+' to '+str(ts.max()) if varying else 'constant '+str(ts[0]))+' fs')
    rep.write('   time differences : '+str(nback)+' backward, '+str(ndup)+' zero (duplicated time), '+str(ngap)+' gaps (> 1.5 timestep), '+str(nirreg)+' other irregular ; '
              +str(repeated)+' snapshots repeating the positions of the previous one read')
    if nback or ndup or ngap or nirreg or repeated or varying:
        bad = np.nonzero((dts <= 0) | (np.abs(dts-ts[1:]) > 1e-3*ts[1:]))[0][:5]
        rep.result('WARNING',('irregular times (first at snapshots '+', '.join(str(k+1) for k in bad)+')' if len(bad) else 'the timestep field changes')
                   +' : badly concatenated restarts ? the MSD, VACF and viscosity assume a constant time step')
    else:
        rep.result('OK')

    #2. Unwrapping
    rep.check('Unwrapping of the positions')
    if not unwrap:
        rep.result('skipped','fewer than 2 snapshots read')
    else:
        nbad = int((np.array(unwrap) > 1e-3).sum())
        rep.write('   between two snapshots read, the change of the unwrapped positions minus that of the wrapped ones, in reduced coordinates of the cell,')
        rep.write('   should be whole numbers (crossings of the faces) : largest deviation %.2e ; pairs of snapshots above 1e-3 : %d' % (max(unwrap),nbad))
        rep.result('WARNING','the unwrapped positions do not follow the wrapped ones : re-parse the OUTCAR or use repair_umd.py') if nbad else rep.result('OK')

    #3. Jumps
    rep.check('Jumps of the atoms between snapshots')
    if natoms_read < 2:
        rep.result('skipped','fewer than 2 snapshots read')
    else:
        rep.write('   largest displacement between two snapshots read (A) : '+', '.join('%s %.3f' % (el,d) for el,d in zip(elements,maxdisp)))
        if novel:
            rep.result('skipped',NOVEL+' : no bound for the displacements')
        elif jumps:
            j = max(jumps,key=lambda x: x[2]/x[3])
            rep.result('WARNING',str(len(jumps))+' atom displacements above 2 max(|v|) dt + 0.2 A, the largest atom '+str(j[1])+' at snapshot '+str(j[0])
                       +' : %.2f A (bound %.2f A) : restarts or wrapping errors ?' % (j[2],j[3]))
        else:
            rep.result('OK')

    #4. Energy drift
    rep.check('Energy drift (EnergyWithDrift, conserved with the thermostat)')
    E = np.array(series['EnergyWithDrift'],dtype=float)[run]
    t = times[run]
    if np.all(np.isnan(E)) or len(E) < 3:
        rep.result('skipped','no EnergyWithDrift in the file')
    elif not np.any(E):
        rep.result('skipped','EnergyWithDrift 0 in every snapshot : not in the file (e.g. a LAMMPS dump without its log)')
    else:
        dE = np.diff(E)
        typical = np.median(np.abs(dE))
        big = np.nonzero(np.abs(dE) > 20.0*typical)[0] if typical > 0 else np.array([],dtype=int)
        smooth = dE.copy()
        smooth[big] = 0.0                                               #the drift is fitted with the jumps removed
        slope = np.polyfit(t,np.concatenate([[0.0],np.cumsum(smooth)]),1)[0]*1000.0/natom*1000.0      #eV/fs -> meV/atom/ps
        rep.write('   drift %.3f meV/atom/ps over %.3g ps (jumps removed), i.e. %.3f meV/atom over the run ; %d jumps above 20 x the median change between snapshots'
                  % (slope,(t[-1]-t[0])/1000.0,slope*(t[-1]-t[0])/1000.0,len(big))+(' (largest %.2f meV/atom)' % (np.abs(dE[big]).max()*1000.0/natom) if len(big) else ''))
        if abs(slope) > 1.0 or len(big):
            text = []
            if abs(slope) > 1.0:
                text.append('drift above 1 meV/atom/ps : time step too long, or electronic convergence (EDIFF) too loose')
            if len(big):
                text.append(str(len(big))+' jumps above 20 x the median change (snapshots '+', '.join(str(k) for k in np.array(steps)[run][big[:6]+1])+(' ...' if len(big) > 6 else '')
                            +') : restarts (thermostat reset), learning steps of an on-the-fly ML run, or electronic steps not converged')
            rep.result('WARNING',' ; '.join(text))
        else:
            rep.result('OK')

    #5. Time step
    rep.check('Time step')
    if novel:
        rep.result('skipped',NOVEL)
    else:
        speeds = np.array(speeds)[~faulty] if (~faulty).any() else np.array(speeds)
        fastest = int(np.argmax(np.median(speeds,axis=0)))
        vfast = np.median(speeds[:,fastest])
        dtmax = 0.1/vfast if vfast > 0 else np.inf
        rep.write('   fastest element %s : median of its largest speed %.4f A/fs (largest %.4f) ; time step for 0.1 A per step : %.2f fs ; timestep of the file %s fs'
                  % (elements[fastest],vfast,speeds[:,fastest].max(),dtmax,TimeStep))
        if TimeStep > dtmax:
            rep.result('WARNING','the fastest atoms move %.2f A per timestep of the file : too long a time step if the snapshots are consecutive MD steps' % (vfast*TimeStep))
        else:
            rep.result('OK','(if the snapshots are not consecutive MD steps, compare the MD time step with %.2f fs)' % dtmax)

    #6. Momentum
    rep.check('Momentum and centre of mass')
    Pcom = np.array(Pcom)
    vcom = (np.linalg.norm(Pcom,axis=1)/masses.sum())[~faulty] if (~faulty).any() else np.linalg.norm(Pcom,axis=1)/masses.sum()
    Rcom = np.array(Rcom)
    dcom = np.linalg.norm(Rcom-Rcom[0],axis=1)
    rep.write('   centre-of-mass velocity : mean %.2e, largest %.2e A/fs ; centre-of-mass displacement : %.3f A at the end, %.3f A at most' % (vcom.mean(),vcom.max(),dcom[-1],dcom.max()))
    if natoms_read >= 2:
        msd_end = [((last_abs[m]-first_abs[m])**2).sum(axis=1).mean() for m in masks]
        slow = int(np.argmin(msd_end))
        rep.write('   squared displacement of the centre of mass at the end %.4f A^2 ; mean squared displacement of the atoms over the run : %s'
                  % (dcom[-1]**2,', '.join('%s %.3f' % (el,x) for el,x in zip(elements,msd_end))+' A^2'))
        if dcom[-1]**2 > 0.05*msd_end[slow]:
            rep.result('WARNING','the centre of mass moved %.3f A^2, %.0f %% of the MSD of %s : the MSD and D include this drift (thermostat without momentum conservation, or initial velocities with a net momentum)'
                       % (dcom[-1]**2,100.0*dcom[-1]**2/msd_end[slow],elements[slow]))
        else:
            rep.result('OK')
    else:
        rep.result('skipped','fewer than 2 snapshots read')

    #7. Net force
    rep.check('Net force')
    if not netforce:
        rep.result('skipped','no forces in the file')
    elif not np.any(np.array(rmsforce) > 0):
        rep.result('skipped','all the forces are 0 : not in the file (e.g. a LAMMPS dump without fx fy fz)')
    else:
        netforce = np.array(netforce)
        rmsforce = np.array(rmsforce)
        fratio = netforce/np.sqrt(natom)/np.where(rmsforce > 0,rmsforce,np.nan)   #error per atom that the net force implies, if spread over the atoms, / rms force
        rep.write('   |sum of the forces| : mean %.4f, largest %.4f eV/A ; rms force of an atom %.3f eV/A ; |sum|/sqrt(N) / rms force : mean %.4f'
                  % (netforce.mean(),netforce.max(),np.nanmean(rmsforce),np.nanmean(fratio)))
        if np.nanmean(fratio) > 0.02:
            rep.result('WARNING','the net force implies force errors above 2 % of the rms force of an atom : energy cutoff or FFT grid too coarse ?')
        else:
            rep.result('OK')

    #8. Equilibration
    rep.check('Equilibration of T, P, E (blocks, as averages.py -a)')
    drifting = []
    for name in ('Temperature','Pressure','InternalEnergy'):
        x = np.array(series[name],dtype=float)[run]
        if np.all(np.isnan(x)):
            rep.write('   '+name+' : not in the file')
            continue
        if not np.any(x):
            rep.write('   '+name+' : 0 in every snapshot, not in the file')
            continue
        mean, err, drift, text = mean_with_error(x,{'Temperature':'K','Pressure':'GPa','InternalEnergy':'eV'}[name])
        rep.write('   %-15s %s ; largest deviation of a quarter %.1f errors' % (name,text,drift))
        if drift > 4.0:
            drifting.append(name)
    rep.result('WARNING',', '.join(drifting)+' : a quarter of the run deviates by more than 4 errors : not equilibrated ? skip the first snapshots with -i') if drifting else rep.result('OK')

    #9. Temperature
    rep.check('Temperature of the velocities')
    if novel:
        rep.result('skipped',NOVEL)
    else:
        nzero = int(np.sum(zero))
        nfaulty = int(faulty.sum())-nzero
        if good.any():
            r = np.median(ratio[good])
            rep.write('   kinetic temperature of the velocities (3N-3 degrees of freedom) / Temperature field : median %.4f (means %.1f K / %.1f K, faulty snapshots left out)'
                      % (r,Tvel[good].mean(),Tfield[good].mean()))
        else:
            r = 1.0
        rep.write('   snapshots with all velocities 0 : %d ; with velocities off by more than a factor 2 : %d%s'
                  % (nzero,nfaulty,(' (snapshots '+', '.join(str(read[k]) for k in np.nonzero(faulty & ~np.array(zero))[0][:5])+')') if nfaulty else ''))
        later = faulty.copy()
        if read[0] == 0 and faulty[0]:
            later[0] = False
            rep.write('   the first snapshot has '+('zero velocities' if zero[0] else 'velocities off by more than a factor 2')+' : expected from the parsers (no previous step ;'
                      +' older versions wrote the positions there) ; left out of checks 5, 6, 9')
        text = []
        if abs(r-1.0) > 0.05:
            text.append('the velocities and the Temperature field differ by more than 5 %')
        nlater = int(later.sum())
        if nlater:
            text.append(str(nlater)+' snapshots after the first with velocities 0 or off by more than a factor 2 (VaspParser.py before October 2026 wrote 0 on the learning steps of ML runs) :'
                        +' re-parse the OUTCAR or repair_umd.py ; left out of checks 5, 6, 9')
        if text:
            rep.result('WARNING',' ; '.join(text))
        else:
            rep.result('OK')

    #10. Cell
    rep.check('Cell (NPT)')
    if any(c is None for c in cells):
        rep.result('skipped','no rprimd lines in the file')
    else:
        cells = np.array(cells)[run]
        volume = np.abs(np.linalg.det(cells))
        if np.ptp(cells,axis=0).max() == 0:
            rep.result('skipped','constant cell (NVT), volume %.3f A^3' % volume[0])
        else:
            drifting = []
            for name,x in (('volume',volume),('a',np.linalg.norm(cells[:,0],axis=1)),('b',np.linalg.norm(cells[:,1],axis=1)),('c',np.linalg.norm(cells[:,2],axis=1))):
                mean, err, drift, text = mean_with_error(x,'A^3' if name == 'volume' else 'A')
                rep.write('   %-7s %s ; largest deviation of a quarter %.1f errors' % (name,text,drift))
                if drift > 4.0:
                    drifting.append(name)
            rep.result('WARNING',', '.join(drifting)+' drifting : the cell is not equilibrated') if drifting else rep.result('OK')

    #11. Stress
    rep.check('Stress tensor (hydrostatic in a liquid)')
    S = [s for s,keep in zip(series['StressTensor'],run) if keep and s and len(s) >= 6]
    if not S:
        rep.result('skipped','no StressTensor in the file')
    elif not np.any(np.array(S)[:,:6]):
        rep.result('skipped','StressTensor 0 in every snapshot : not in the file (e.g. a LAMMPS log without Pxx Pyy Pzz Pxy Pxz Pyz)')
    else:
        S = np.array(S)[:,:6]
        stats = [mean_with_error(S[:,k],'GPa') for k in range(6)]
        rep.write('   '+' ; '.join(lab+' '+st[3] for lab,st in zip(('xx','yy','zz','xy','yz','zx'),stats)))
        text = []
        off = [lab for lab,st in zip(('xy','yz','zx'),stats[3:]) if abs(st[0]) > 4.0*st[1]]
        if off:
            text.append('mean off-diagonal stresses '+', '.join(off)+' not 0 (> 4 errors)')
        diag = [a+'-'+b for (a,sa),(b,sb) in ((('xx',stats[0]),('yy',stats[1])),(('yy',stats[1]),('zz',stats[2])),(('zz',stats[2]),('xx',stats[0]))) if abs(sa[0]-sb[0]) > 4.0*np.hypot(sa[1],sb[1])]
        if diag:
            text.append('diagonal stresses '+', '.join(diag)+' differ (> 4 errors)')
        rep.result('WARNING',' ; '.join(text)+' : non-hydrostatic : normal for a strained crystal ; for a liquid, run too short or cell not relaxed') if text else rep.result('OK')

    #12. Overlap of the core spheres
    rep.check('Overlap of the core spheres')
    if radius is None:
        rep.result('skipped','no radius file (-r)')
    else:
        dmins = np.array(dmins)
        percents = np.array(percents)
        rep.write('   {:<10s}{:>12s}{:>12s}{:>14s}{:>12s}{:>16s}'.format('pair','r1+r2 (A)','dmin (A)','max overlap %','snapshot','snapshots > -l'))
        for p,(lab,(i,j)) in enumerate(zip(labels,pairs)):
            if not np.isfinite(dmins[:,p]).any():
                continue
            kmax = int(np.argmax(percents[:,p]))
            rep.write('   {:<10s}{:>12.3f}{:>12.3f}{:>14.2f}{:>12d}{:>16d}'.format(lab,radius[elements[i]]+radius[elements[j]],np.nanmin(np.where(np.isfinite(dmins[:,p]),dmins[:,p],np.nan)),
                      percents[kmax,p],read[kmax],int((percents[:,p] > level).sum())))
        nbad = int((percents.max(axis=1) > level).sum())
        rep.write('   per snapshot in '+base+'.overlap.dat')
        rep.result('WARNING',str(nbad)+' snapshots with an overlap above '+str(level)+' % : pseudopotentials not transferable at this pressure ?') if nbad else rep.result('OK')

    #13. Smearing
    rep.check('Smearing (electronic entropy)')
    Sel = np.array(series['ElectronicEntropy'],dtype=float)[run]
    if np.all(np.isnan(Sel)):
        rep.result('skipped','no ElectronicEntropy in the file')
    else:
        tsel = np.nanmean(np.abs(Sel))*1000.0/natom
        rep.write('   |T S_el| : mean %.4f meV/atom, largest %.4f meV/atom' % (tsel,np.nanmax(np.abs(Sel))*1000.0/natom))
        rep.result('WARNING','above 1 meV/atom : smearing (SIGMA) too wide, or a metallic system (then the free energy is the relevant energy)') if tsel > 1.0 else rep.result('OK')

    #14. Charge balance
    rep.check('Formal charge balance of the composition')
    q = {}
    for el in elements:
        q[el] = charges.get(el,cr.Elements2charge(el))
    unknown = [el for el in elements if q[el] is None]
    if unknown:
        rep.result('skipped','no formal charge for '+', '.join(unknown)+' : give it with -q')
    else:
        total = sum(q[el]*n for el,n in zip(elements,MyCrystal.types))
        rep.write('   '+' + '.join('%d %s(%+g)' % (n,el,q[el]) for el,n in zip(elements,MyCrystal.types))+' = %+g' % total
                  +('   (changed with -q : '+', '.join(el+' %+g' % charges[el] for el in charges)+')' if charges else ''))
        rep.result('WARNING','the formal charges do not balance (%+g e) : other oxidation states (e.g. Fe3+ : -q Fe:3 ; C, N, S in molecules), or a charged cell' % total) if abs(total) > 1e-6 else rep.result('OK')

    rep.write('\nSummary : '+str(len(rep.warnings))+' warning(s)')
    for w in rep.warnings:
        rep.write('   '+w)
    with open(base+'.check.dat','w') as f:
        f.write('\n'.join(rep.lines)+'\n')
    print('Report written in',base+'.check.dat'+(' ; overlaps in '+base+'.overlap.dat' if radius is not None else ''))


if __name__ == "__main__":
    if not umdpf.run_on_all(sys.argv[1:],main):#-f all : every UMD file of the folder
        main(sys.argv[1:])
