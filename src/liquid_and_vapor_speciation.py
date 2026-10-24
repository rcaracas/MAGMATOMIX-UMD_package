#!/usr/bin/env python3
# coding: utf-8

#Separates the species of a popul file (any -r, typically r1) into liquid and vapor, from the clusters of an r0 step file :
#an atom is in the vapor at a step when its r0 cluster has fewer than -l atoms ; a species is in the liquid when one of its atoms is.

import sys
import getopt
import re
from itertools import groupby


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


def write_phase(filename,species,phasesteps,nsteps,dt):#popul and stat files of one phase ; phasesteps : for each species, its steps in this phase
    with open(filename+'.popul.dat','w') as nf:
        nf.write('Formula\tBegin (step)\tEnd (step)\tLifetime (fs)\tComposition\n')
        stats = {}
        total = 0.0
        for (formula,begin,end,lifetime,composition),steps in zip(species,phasesteps):
            for _,run in groupby(enumerate(steps),lambda x: x[1]//nsteps - x[0]):#consecutive snapshots
                run = [s for _,s in run]
                life = len(run)*nsteps*dt
                nf.write(formula+'\t'+str(run[0])+'\t'+str(run[-1])+'\t'+str(life)+'\t'+composition+'\n')
                total += life
                if formula in stats:
                    stats[formula][0] += life
                else:
                    stats[formula] = [life,len(composition.split(','))]
    with open(filename+'.stat.dat','w') as nf:
        nf.write('Cluster\tTime (fs)\tPercent\tNumber of atoms\n')
        for formula in stats:
            nf.write(formula+'\t'+str(stats[formula][0])+'\t'+str(stats[formula][0]/total)+'\t'+str(stats[formula][1])+'\n')


def main(argv):
    stepfile = ''
    populfile = ''
    length = 12
    usage = 'liquid_and_vapor_speciation.py -s <file.r0.step.dat> -p <file.r1.popul.dat> -l <length>'
    try:
        opts, args = getopt.getopt(argv,"hs:p:l:")
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('Separates the species of a popul file into liquid and vapor species, from the r0 clusters of a step file of the same bonding file.')
            print(usage)
            print('-s : step file of speciation_and_angles.py -r 0 : an atom is in the vapor at a step when its r0 cluster has fewer than -l atoms.')
            print('-p : popul file (typically -r 1) : a species is in the liquid at a step when one of its atoms is. Its steps must be in the step file')
            print('     (same sampling -s of speciation_and_angles.py, or a denser one for the step file).')
            print('-l : smallest r0 cluster of the liquid, in atoms. Default 12.')
            print('Output : <popul base>.vapor.popul.dat, .vapor.stat.dat, .liquid.popul.dat, .liquid.stat.dat (popul and stat formats).')
            sys.exit()
        elif opt == '-s':
            stepfile = str(arg)
        elif opt == '-p':
            populfile = str(arg)
        elif opt == '-l':
            length = int(arg)

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

    base = populfile[:-10] if populfile.endswith('.popul.dat') else populfile
    write_phase(base+'.vapor',species,vaporsteps,nsteps,timestep)
    write_phase(base+'.liquid',species,liquidsteps,nsteps,timestep)
    print('Vapor and liquid species written in',base+'.vapor.popul.dat,',base+'.vapor.stat.dat,',base+'.liquid.popul.dat,',base+'.liquid.stat.dat')


if __name__ == "__main__":
    main(sys.argv[1:])
