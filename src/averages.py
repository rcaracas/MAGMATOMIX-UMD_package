#!/usr/bin/env python3
"""
@author: Razvan Caracas, Anais Kobsch, Tim Bögels
Averages of the fields of UMD files : one quantity or P, T, E of one file, or the thermodynamics (density, P, T, E with their
errors, Cv) of all the UMD files of a folder (formerly fullaverages.py)
"""

import numpy as np
import sys
import glob
import getopt
import os
import subprocess
import umd_processes_fast as umdpf
import crystallography as cr
from vibrations_umd import mean_error

KB = 8.617333262E-05        #Boltzmann constant, eV/K
EV = 1.602176634E-19        #J
NA = 6.02214076E23          #Avogadro constant


#labels for the components of the multi-valued UMD fields; StressTensor follows the VASP order xx yy zz xy yz zx
ComponentLabels = {'StressTensor':['xx','yy','zz','xy','yz','zx'], 'acell':['a','b','c']}

def is_number(s):
    try:

        float(s)
        return True
    except ValueError:
        return False

def grep_pattern(FileName, Pattern, SkipSteps, required=True):
    """Extraction of all the numerical values following Pattern, skipping the first SkipSteps snapshots.
    Returns an array with one row per component (1 row for scalars, 6 for StressTensor, ...) and the unit ;
    (None, '') when the field is absent and not required."""
    data = []
    unit = ''
    anchor=Pattern.split()[0]
    try:
        patterns=subprocess.check_output(['grep',Pattern,FileName]).decode("utf-8")
    except subprocess.CalledProcessError:
        if not required:
            return None, ''
        print('Parameter ',Pattern,' not found in ',FileName)
        sys.exit()
    greps=patterns.split('\n')
    for isteps in range(SkipSteps,len(greps)):
        elems=greps[isteps].split()
        for ii in range(len(elems)):
            if elems[ii] == anchor:
                values = []         #all the consecutive numbers after the anchor, e.g. the 6 components of StressTensor
                for jj in range(ii+1,len(elems)):
                    if is_number(elems[jj]):
                        values.append(float(elems[jj]))
                    elif values:
                        unit = elems[jj]
                        break
                if values:
                    data.append(values)
                break
    if not data:
        if not required:
            return None, ''
        print('No numerical values found after ',anchor,' in ',FileName)
        sys.exit()
    ncomp = min(len(values) for values in data)
    data = np.array([values[:ncomp] for values in data]).T     #one row per component
    return data, unit

def statistics(x,nblocks=4):#Mean, stdev, error of the mean (integrated autocorrelation time), whether the error is only a lower bound
    #(fewer than 50 independent samples), and the largest deviation of the means of nblocks blocks from the run mean, in errors of a block mean
    if np.ptp(x) == 0.0:                                                #constant, e.g. the density of an NVT run
        return x[0], 0.0, 0.0, False, 0.0
    err, tau = mean_error(x,with_tau=True)
    blocklen = len(x)//nblocks
    blocks = np.array([x[kk*blocklen:(kk+1)*blocklen].mean() for kk in range(nblocks)])
    drift = np.abs(blocks-x.mean()).max()/(err*np.sqrt(nblocks)) if err > 0.0 else 0.0
    return x.mean(), x.std(), err, len(x)/(2.0*tau) < 50.0, drift


def kinetic_included(IE,K,TEE,ES=None):#Whether InternalEnergy contains KineticEnergyIons (VaspParser.py since 2024 ; not in older UMD files), and how it was found.
    #TotalElectronicEnergy is the free energy TOTEN (VaspParser.py ; 0 in the files of its versions of 2026 before the fix), so InternalEnergy =
    #TotalElectronicEnergy + ElectronicEntropy (+ KineticEnergyIons) ; otherwise (thermodynamic integration, SCALEE < 1) from the correlation
    if TEE is not None and np.any(TEE != 0.0):
        E0 = TEE + (ES if ES is not None else 0.0)
        if np.allclose(IE,E0+K,atol=1e-3):
            return True, 'InternalEnergy = TotalElectronicEnergy + ElectronicEntropy + KineticEnergyIons'
        if np.allclose(IE,E0,atol=1e-3):
            return False, 'InternalEnergy = TotalElectronicEnergy + ElectronicEntropy'
    ratio = np.cov(IE,K)[0,1]/np.var(K) if np.var(K) > 0.0 else 0.0  #1 if included, 0 if not : potential and kinetic energies are uncorrelated in NVT
    return bool(ratio > 0.5), 'cov(InternalEnergy,KineticEnergyIons)/var(KineticEnergyIons) = ' + '%.2f' % ratio


def thermodynamics(FileName,SkipSteps):#Density, P, T, E (eV per cell, ionic kinetic energy included) and Cv (eV/K per cell) of one UMD file : statistics, and notes
    MyCrystal,TimeStep = umdpf.read_bigheader_umd(FileName)
    natom = MyCrystal.natom
    mass = sum(cr.Elements2rest(el)[3]*n for el,n in zip(MyCrystal.elements,MyCrystal.types))/NA      #g per cell
    fields = {}
    for name in ('Density','Pressure','Temperature','InternalEnergy','KineticEnergyIons','TotalElectronicEnergy','ElectronicEntropy','rprimd_a','acell'):
        data, unit = grep_pattern(FileName,name,SkipSteps,required=False)
        fields[name] = data
    notes = []
    results = {}
    if fields['Density'] is not None:
        rho = fields['Density'][0]
    else:
        cells = umdpf.read_cells(FileName,len(fields['Temperature'][0]),1,SkipSteps)
        rho = mass/(np.abs(np.linalg.det(cells))*1E-24)
        notes.append('no Density in the file : computed from the cell')
    results['rho'] = statistics(rho)
    for name,key in (('Pressure','P'),('Temperature','T')):
        if fields[name] is None:
            print('Parameter ',name,' not found in ',FileName)
            sys.exit()
        results[key] = statistics(fields[name][0])
    IE = fields['InternalEnergy'][0]
    if fields['KineticEnergyIons'] is not None:
        K = fields['KineticEnergyIons'][0]
        TEE = fields['TotalElectronicEnergy'][0] if fields['TotalElectronicEnergy'] is not None else None
        ES = fields['ElectronicEntropy'][0] if fields['ElectronicEntropy'] is not None else None
        included, how = kinetic_included(IE,K,TEE,ES)
        U = IE-K if included else IE                                     #potential energy
        E = U+K
        notes.append('InternalEnergy ' + ('includes' if included else 'does not include') + ' the kinetic energy of the ions (' + how + ')')
        T = results['T'][0]
        ratio = np.var(K)/(1.5*natom*KB**2*T**2)
        notes.append('fluctuations of the kinetic energy / canonical value 3/2 N kB^2 T^2 = ' + '%.2f' % ratio)
        if abs(ratio-1.0) > 0.25:
            notes.append('WARNING : the kinetic energy does not fluctuate as in the canonical ensemble (thermostat ?) : Cv is not reliable')
    else:
        U = E = IE
        notes.append('no KineticEnergyIons : InternalEnergy taken as the energy, and as the potential energy for Cv')
    results['E'] = statistics(E)
    #Cv of the canonical ensemble : fluctuations of the potential energy + the kinetic part 3/2 N kB (exact for NVT ; it avoids the noise of
    #the kinetic fluctuations, and the potential and kinetic energies are independent)
    T = results['T'][0]
    dU2 = (U-U.mean())**2
    mean, std, err, lower, drift = statistics(dU2)
    results['Cv'] = (mean/(KB*T**2) + 1.5*natom*KB, err/(KB*T**2), lower, drift)
    cell = fields['rprimd_a'] if fields['rprimd_a'] is not None else fields['acell']
    if cell is not None and np.ptp(cell,axis=1).max() > 0.0:
        notes.append('WARNING : the cell changes (NPT) : Cv from the energy fluctuations is valid only for NVT runs')
    return len(IE), natom, mass, results, notes


def error_string(err,lower,scale=1.0):
    return ('>' if lower else '') + '{:1.1e}'.format(err*scale)


def write_thermo(files,SkipSteps,units):#Thermodynamics of the UMD files, in fullthermo.dat
    if units == 0:
        firstline = 'file\tnsteps\tDensity(g/cm3)\tstdev_rho\terr_rho\tPressure(GPa)\tstdev_P\terr_P\tTemperature(K)\tstdev_T\terr_T\tInternEnergy(eV/atm)\tstdev_E\terr_E\tCvm(Nkb)\tstdev_Cvm\n'
    else:
        firstline = 'file\tnsteps\tDensity(g/cm3)\tstdev_rho\terr_rho\tPressure(GPa)\tstdev_P\terr_P\tTemperature(K)\tstdev_T\terr_T\tTemperature(eV)\tstdev_T\terr_T\tInternEnergy(eV/atm)\tstdev_E\terr_E\tInternEnergy(eV/unitcell)\tstdev_E\terr_E\tInternEnergy(J/g)\tstdev_E\terr_E\tCvm(Nkb)\tstdev_Cvm\tCvm(J/K/mol)\tstdev_Cvm\tCvm(J/K/g)\tstdev_Cvm\tCv(J/K/unitcell)\tstdev_Cv\n'
    MyCrystal,TimeStep = umdpf.read_bigheader_umd(files[0])
    f = open('fullthermo.dat','w')
    f.write('elements\t' + '\t'.join(MyCrystal.elements) + '\n')
    f.write('number\t' + '\t'.join(str(n) for n in MyCrystal.types) + '\n')
    f.write(firstline)
    for FileName in files:
        print('\nAveraging in file',FileName)
        Crystal,TimeStep = umdpf.read_bigheader_umd(FileName)
        if Crystal.elements != MyCrystal.elements or Crystal.types != MyCrystal.types:
            print('WARNING : composition',Crystal.elements,Crystal.types,'differs from the header of fullthermo.dat (',files[0],')')
        nsteps, natom, mass, results, notes = thermodynamics(FileName,SkipSteps)
        for note in notes:
            print('   ',note)
        for key,name in (('rho','Density'),('P','Pressure'),('T','Temperature'),('E','Energy'),('Cv','Cv')):
            drift = results[key][-1]
            if drift > 4.0:
                print('    WARNING : a quarter of the run differs from the mean',name,'by','%.1f' % drift,'errors : drift (equilibration ? skip the first snapshots with -i)')
            if results[key][-2]:
                print('    WARNING : fewer than 50 independent samples of',name,': its error is a lower bound (">")')
        rho, P, T, E, Cv = results['rho'], results['P'], results['T'], results['E'], results['Cv']
        Cv_Nkb = 1.0/(KB*natom)                          #eV/K per cell -> N kB
        line = [FileName, str(nsteps), '{:1.4f}'.format(rho[0]), '{:1.1e}'.format(rho[1]), error_string(rho[2],rho[3]),
                '{:1.4e}'.format(P[0]), '{:1.1e}'.format(P[1]), error_string(P[2],P[3]),
                '{:1.1f}'.format(T[0]), '{:1.1f}'.format(T[1]), error_string(T[2],T[3])]
        if units == 0:
            line += ['{:1.5e}'.format(E[0]/natom), '{:1.1e}'.format(E[1]/natom), error_string(E[2],E[3],1.0/natom),
                     '{:1.3f}'.format(Cv[0]*Cv_Nkb), error_string(Cv[1],Cv[2],Cv_Nkb)]
        else:
            line += ['{:1.4e}'.format(T[0]*KB), '{:1.1e}'.format(T[1]*KB), error_string(T[2],T[3],KB),
                     '{:1.5e}'.format(E[0]/natom), '{:1.1e}'.format(E[1]/natom), error_string(E[2],E[3],1.0/natom),
                     '{:1.6e}'.format(E[0]), '{:1.1e}'.format(E[1]), error_string(E[2],E[3]),
                     '{:1.5e}'.format(E[0]*EV/mass), '{:1.1e}'.format(E[1]*EV/mass), error_string(E[2],E[3],EV/mass),
                     '{:1.3f}'.format(Cv[0]*Cv_Nkb), error_string(Cv[1],Cv[2],Cv_Nkb),
                     '{:1.3e}'.format(Cv[0]*EV*NA/natom), error_string(Cv[1],Cv[2],EV*NA/natom),
                     '{:1.3e}'.format(Cv[0]*EV/mass), error_string(Cv[1],Cv[2],EV/mass),
                     '{:1.3e}'.format(Cv[0]*EV), error_string(Cv[1],Cv[2],EV)]
        print('    density','{:1.4f}'.format(rho[0]),'+/-',error_string(rho[2],rho[3]),'g/cm3, P','{:1.4g}'.format(P[0]),'+/-',error_string(P[2],P[3]),
              'GPa, T','{:1.1f}'.format(T[0]),'+/-',error_string(T[2],T[3]),'K, E','{:1.6g}'.format(E[0]/natom),'+/-',error_string(E[2],E[3],1.0/natom),
              'eV/atom, Cv','{:1.3f}'.format(Cv[0]*Cv_Nkb),'+/-',error_string(Cv[1],Cv[2],Cv_Nkb),'N kB')
        f.write('\t'.join(line) + '\n')
    f.close()
    print('\nThermodynamics written in fullthermo.dat')


def plot_series(series):
    """Plots one series at a time, with Previous/Next buttons to move between them.
    Each series is (title, unit, data, average, stdev) and gets its own y-scale, since e.g. diagonal and
    off-diagonal stresses, or pressure and energy, can differ by orders of magnitude."""
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Button
    fig, ax = plt.subplots()
    line, = ax.plot(series[0][2])
    avgline = ax.axhline(series[0][3], color='m')
    ax.set_xlabel('step')
    state = {'iseries':0}
    def draw():
        iseries = state['iseries']
        title, unit, data, average, stdev = series[iseries]
        line.set_data(np.arange(len(data)), data)
        avgline.set_ydata([average,average])
        if len(series) > 1:
            title = title + '  (' + str(iseries+1) + '/' + str(len(series)) + ')'
        ax.set_title(title + ':  mean = ' + '%.6g' % average + ' +/- ' + '%.3g' % stdev + ' ' + unit)
        ax.set_ylabel(title.split()[0] + ' (' + unit + ')' if unit else title.split()[0])
        ax.relim()
        ax.autoscale(enable=True)      #re-enabled in case the toolbar zoom switched it off
        fig.canvas.draw_idle()
    draw()
    if len(series) > 1:
        fig.subplots_adjust(bottom=0.2)
        def step(delta):
            state['iseries'] = (state['iseries'] + delta) % len(series)
            draw()
        bprev = Button(fig.add_axes([0.59, 0.03, 0.15, 0.06]), 'Previous')
        bnext = Button(fig.add_axes([0.75, 0.03, 0.15, 0.06]), 'Next')
        bprev.on_clicked(lambda event: step(-1))
        bnext.on_clicked(lambda event: step(1))
        fig._umd_buttons = (bprev, bnext)       #keeps the buttons alive, otherwise they are garbage-collected and stop responding
    plt.show()


def main(argv):
    files = []
    Pattern = ''
    SkipSteps = 0
    LongOutput = False
    Visual = False
    Thermo = False
    units = 0
    umdpf.headerumd()
    usage = 'averages.py -f <UMDFiles> -p <Pattern> -i <SkipSteps> (or -s <SkipSteps>) -l -v -a -u <units>'
    try:
        opts, arg = getopt.getopt(argv,"hf:p:i:s:lvau:")
    except getopt.GetoptError:
        print(usage)
        sys.exit()
    for opt, arg in opts:
        if opt == '-h':
            print(usage)
            print('averages.py program to extract the average and the spread of the numerical values in UMD files')
            print('   -f <UMDFiles>    the UMD files, separated by commas. Default: every *.umd.dat of the current folder')
            print('   -p <Pattern>     the quantity to average, e.g. Pressure, Temperature, StressTensor, acell')
            print('                    default: Pressure, Temperature, and InternalEnergy, written on one line per file as')
            print('                    file  #PstdTstdEstd  nsteps  P  stdev_P  T  stdev_T  E  stdev_E')
            print('                    multi-valued quantities (StressTensor, acell, ...) are averaged component by component')
            print('   -i or -s <SkipSteps>   number of initial snapshots to skip (default 0 = use all the snapshots)')
            print('   -l               long output: mean, variance, stdev and error of the mean printed for every quantity and component')
            print('   -v               display the evolution of every quantity and component (default: no display)')
            print('   -a               all the thermodynamics (formerly fullaverages.py), written in fullthermo.dat, one line per file:')
            print('                    the number of snapshots, then mean, stdev and error of the mean of the density, P, T and E, and Cv with its error.')
            print('                    E = the internal energy with the kinetic energy of the ions (included in InternalEnergy by the current')
            print('                    VaspParser.py, added from KineticEnergyIons for older UMD files : detected for each file).')
            print('                    Cv = <dU^2>/(kB T^2) + 3/2 N kB, U = the potential energy : valid only for NVT runs.')
            print('                    Errors of the mean from the integrated autocorrelation time ; ">" : fewer than 50 independent samples,')
            print('                    the error is a lower bound. Warnings when a quarter of the run deviates by more than 4 errors (drift).')
            print('   -u <units>       with -a, 0 = minimal (default), 1 = also T in eV, E per cell and in J/g, Cv in J/K/mol, J/K/g, J/K per cell')
            sys.exit()
        elif opt == '-f':
            files = arg.split(',')
        elif opt == '-p':
            Pattern = str(arg)
        elif opt in ('-i', '-s'):
            SkipSteps = int(arg)
        elif opt == '-l':
            LongOutput = True
        elif opt == '-v':
            Visual = True
        elif opt == '-a':
            Thermo = True
        elif opt == '-u':
            units = int(arg)
    if not files:
        files = sorted(glob.glob('*.umd.dat'))
    if not files:
        print('No UMD file in the current folder : give them with -f')
        sys.exit()
    for FileName in files:
        if not os.path.isfile(FileName):
            print('No input file or file "',FileName,'" does not exist')
            sys.exit()
    if Thermo:
        write_thermo(files,SkipSteps,units)
        return
    for FileName in files:
        averages(FileName,Pattern,SkipSteps,LongOutput,Visual)


def averages(FileName,Pattern,SkipSteps,LongOutput,Visual):#One quantity, or P, T and E, of one file
    if Pattern:
        Patterns = [Pattern]
    else:
        Patterns = ['Pressure', 'Temperature', 'InternalEnergy']
    #****Calculations
    series = []         #(title, unit, data, average, stdev) for every quantity and component
    short = []          #(nsteps, means, stdevs) for every quantity
    for pattern in Patterns:
        data, unit = grep_pattern(FileName, pattern, SkipSteps)
        anchor = pattern.split()[0]
        means = data.mean(axis=1)
        variances = data.var(axis=1)
        stdevs = np.sqrt(variances)
        short.append((data.shape[1], means, stdevs))
        ncomp = len(data)
        labels = ComponentLabels.get(anchor,[str(ii+1) for ii in range(ncomp)])
        for ic in range(ncomp):
            title = anchor + ' ' + labels[ic] if ncomp > 1 else anchor
            series.append((title, unit, data[ic], means[ic], stdevs[ic]))
            if LongOutput:
                print('Averages of ',title,' over ',data.shape[1],' steps: mean = ',means[ic],' variance = ', variances[ic], ' stdev = ', stdevs[ic],
                      ' error of the mean = ', mean_error(data[ic]), unit)
    if not LongOutput:
        if Pattern:
            nsteps, means, stdevs = short[0]
            arraytowrite = FileName + '\t#' + Pattern.split()[0] + '\t' + str(nsteps)
            for ic in range(len(means)):
                arraytowrite = arraytowrite + '\t' + str(means[ic]) + '\t' + str(stdevs[ic])
        else:
            arraytowrite = FileName + '\t#PstdTstdEstd\t' + str(short[0][0])
            for nsteps, means, stdevs in short:
                arraytowrite = arraytowrite + '\t' + str(means[0]) + '\t' + str(stdevs[0])
        print(arraytowrite)
    if Visual:
        plot_series(series)


if __name__ == "__main__":
   main(sys.argv[1:])
