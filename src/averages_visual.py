#!/usr/bin/env python3
"""
Created on Tue May  2 14:54:47 2022

@author: Tim Bögels
"""

import numpy
import sys
import getopt
import os
import subprocess
import matplotlib.pyplot as plt
from matplotlib.widgets import Button


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

#labels for the components of the multi-valued UMD fields; StressTensor follows the VASP order xx yy zz xy yz zx
ComponentLabels = {'StressTensor':['xx','yy','zz','xy','yz','zx'], 'acell':['a','b','c']}

def plot_components(anchor,data,unit,averages,stdevs):
    ncomp = len(data)
    labels = ComponentLabels.get(anchor,[str(ii+1) for ii in range(ncomp)])
    fig, ax = plt.subplots()
    line, = ax.plot(data[0])
    avgline = ax.axhline(averages[0], color='m')
    ax.set_xlabel('step')
    state = {'icomp':0}
    def draw():
        ic = state['icomp']
        line.set_ydata(data[ic])
        avgline.set_ydata([averages[ic],averages[ic]])
        title = anchor
        if ncomp > 1:
            title = title + ' ' + labels[ic] + '  (' + str(ic+1) + '/' + str(ncomp) + ')'
        ax.set_title(title + ':  mean = ' + '%.6g' % averages[ic] + ' +/- ' + '%.3g' % stdevs[ic] + ' ' + unit)
        ax.set_ylabel(anchor + ' (' + unit + ')' if unit else anchor)
        ax.relim()
        ax.autoscale(enable=True)      #re-enabled in case the toolbar zoom switched it off
        fig.canvas.draw_idle()
    draw()
    if ncomp > 1:       #one component at a time, since diagonal and off-diagonal stresses can differ by orders of magnitude
        fig.subplots_adjust(bottom=0.2)
        def step(delta):
            state['icomp'] = (state['icomp'] + delta) % ncomp
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
    data = []
    SkipSteps=0
    d=0
    anchor=''
    try:
        opts, arg = getopt.getopt(argv,"hf:p:i:",["fFileName","pPattern","iSkipSteps"])
    except getopt.GetoptError:
        print('average.py -f <FileName> -p <Pattern> -i <SkipSteps> ')
        sys.exit(d)
    for opt, arg in opts:
        if opt == '-h':
            print('average.py program to extract and average numerical values')
            print('average_visual.py -f <FileName> -p <Pattern> -i <SkipSteps>')
            sys.exit()
        elif opt in ("-f", "--fFileName"):
            FileName = str(arg)
        elif opt in ("-p", "--pPattern"):
            Pattern = str(arg)
            anchor=Pattern.split()[0]
        elif opt in ("-i", "--iSkipSteps"):
            SkipSteps = int(arg)
    if(anchor==''):
        print('No parameter given ; no mean calculated. Use -p to specify a parameter.')
        sys.exit()
    if (os.path.isfile(FileName)): 
        try :
            patterns=subprocess.check_output(['grep',Pattern,FileName])
        except subprocess.CalledProcessError :
            print('Parameter ',Pattern,' not found.')
            sys.exit()
        patternsstr=patterns.decode()
        greps=patternsstr.split('\n')
        unit = ''
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
            print('No numerical values found after ',anchor)
            sys.exit()
        ncomp = min(len(values) for values in data)
        data = numpy.array([values[:ncomp] for values in data]).T     #one row per component
        averages = data.mean(axis=1)
        variances = data.var(axis=1)
        stdevs = numpy.sqrt(variances)
        for ic in range(ncomp):
            if ncomp > 1:
                print('Component ',ic+1,' of ',anchor,': ',end='')
            print('Averages over ',data.shape[1],' steps: mean = ',averages[ic],' variance = ', variances[ic], ' stdev = ', stdevs[ic])
        plot_components(anchor,data,unit,averages,stdevs)
    else:
        print('No input file or file ',FileName,'does not exist')
        sys.exit()


if __name__ == "__main__":
   main(sys.argv[1:])

