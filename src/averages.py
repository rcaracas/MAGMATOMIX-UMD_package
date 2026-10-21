#!/usr/bin/env python3
"""
@author: Razvan Caracas, Anais Kobsch, Tim Bögels
"""

import numpy as np
import sys
import getopt
import os
import subprocess
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import umd_processes_fast as umdpf


#labels for the components of the multi-valued UMD fields; StressTensor follows the VASP order xx yy zz xy yz zx
ComponentLabels = {'StressTensor':['xx','yy','zz','xy','yz','zx'], 'acell':['a','b','c']}

def is_number(s):
    try:

        float(s)
        return True
    except ValueError:
        return False

def grep_pattern(FileName, Pattern, SkipSteps):
    """Extraction of all the numerical values following Pattern, skipping the first SkipSteps snapshots.
    Returns an array with one row per component (1 row for scalars, 6 for StressTensor, ...) and the unit."""
    data = []
    unit = ''
    anchor=Pattern.split()[0]
    try:
        patterns=subprocess.check_output(['grep',Pattern,FileName]).decode("utf-8")
    except subprocess.CalledProcessError:
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
        print('No numerical values found after ',anchor,' in ',FileName)
        sys.exit()
    ncomp = min(len(values) for values in data)
    data = np.array([values[:ncomp] for values in data]).T     #one row per component
    return data, unit

def plot_series(series):
    """Plots one series at a time, with Previous/Next buttons to move between them.
    Each series is (title, unit, data, average, stdev) and gets its own y-scale, since e.g. diagonal and
    off-diagonal stresses, or pressure and energy, can differ by orders of magnitude."""
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
    FileName = ''
    Pattern = ''
    SkipSteps = 0
    LongOutput = False
    Visual = False
    umdpf.headerumd()
    usage = 'averages.py -f <UMDFile> -p <Pattern> -i <SkipSteps> (or -s <SkipSteps>) -l -v'
    try:
        opts, arg = getopt.getopt(argv,"hf:p:i:s:lv")
    except getopt.GetoptError:
        print(usage)
        sys.exit()
    for opt, arg in opts:
        if opt == '-h':
            print(usage)
            print('averages.py program to extract the average and the spread of the numerical values in a UMD file')
            print('   -f <UMDFile>     the UMD file')
            print('   -p <Pattern>     the quantity to average, e.g. Pressure, Temperature, StressTensor, acell')
            print('                    default: Pressure, Temperature, and InternalEnergy, written on one line as')
            print('                    file  #PstdTstdEstd  nsteps  P  stdev_P  T  stdev_T  E  stdev_E')
            print('                    multi-valued quantities (StressTensor, acell, ...) are averaged component by component')
            print('   -i or -s <SkipSteps>   number of initial snapshots to skip (default 0 = use all the snapshots)')
            print('   -l               long output: mean, variance, and stdev printed for every quantity and component')
            print('   -v               display the evolution of every quantity and component (default: no display)')
            sys.exit()
        elif opt == '-f':
            FileName = str(arg)
        elif opt == '-p':
            Pattern = str(arg)
        elif opt in ('-i', '-s'):
            SkipSteps = int(arg)
        elif opt == '-l':
            LongOutput = True
        elif opt == '-v':
            Visual = True
    if not os.path.isfile(FileName):
        print('No input file or file "',FileName,'" does not exist')
        sys.exit()
    if Pattern:
        Patterns = [Pattern]
    else:
        Patterns = ['Pressure', 'Temperature', 'InternalEnergy']
    #****Calculations
    series = []         #(title, unit, data, average, stdev) for every quantity and component
    short = []          #(nsteps, averages, stdevs) for every quantity
    for pattern in Patterns:
        data, unit = grep_pattern(FileName, pattern, SkipSteps)
        anchor = pattern.split()[0]
        averages = data.mean(axis=1)
        variances = data.var(axis=1)
        stdevs = np.sqrt(variances)
        short.append((data.shape[1], averages, stdevs))
        ncomp = len(data)
        labels = ComponentLabels.get(anchor,[str(ii+1) for ii in range(ncomp)])
        for ic in range(ncomp):
            title = anchor + ' ' + labels[ic] if ncomp > 1 else anchor
            series.append((title, unit, data[ic], averages[ic], stdevs[ic]))
            if LongOutput:
                print('Averages of ',title,' over ',data.shape[1],' steps: mean = ',averages[ic],' variance = ', variances[ic], ' stdev = ', stdevs[ic], unit)
    if not LongOutput:
        if Pattern:
            nsteps, averages, stdevs = short[0]
            arraytowrite = FileName + '\t#' + Pattern.split()[0] + '\t' + str(nsteps)
            for ic in range(len(averages)):
                arraytowrite = arraytowrite + '\t' + str(averages[ic]) + '\t' + str(stdevs[ic])
        else:
            arraytowrite = FileName + '\t#PstdTstdEstd\t' + str(short[0][0])
            for nsteps, averages, stdevs in short:
                arraytowrite = arraytowrite + '\t' + str(averages[0]) + '\t' + str(stdevs[0])
        print(arraytowrite)
    if Visual:
        plot_series(series)


if __name__ == "__main__":
   main(sys.argv[1:])
