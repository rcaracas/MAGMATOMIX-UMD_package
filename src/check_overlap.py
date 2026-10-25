#!/usr/bin/env python3
###
#AUTHORS: ANAIS KOBSCH, RAZVAN CARACAS
###
#Overlap of the core spheres (pseudopotential / PAW radii) of the atoms during a run : for each snapshot and each pair of elements, the shortest
#distance and the volume of the overlap of the two spheres at that distance, in % of the sum of their volumes. A quick check of the transferability
#of the pseudopotentials at high pressure. The radii come from the user (radius file), so it works for any ab initio code.

import sys, getopt, os.path, re
import numpy as np
import umd_processes_fast as umdpf


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


def main(argv):
    umdfile = 'output.umd.dat'
    radiusfile = 'radius.inp'
    Nsteps = 1
    InitialStep = 0
    level = 12.0
    umdpf.headerumd()
    usage = 'check_overlap.py -f <UMD_filename> -r <radius_filename> -s <Sampling_Frequency> -i <InitialStep> -l <level_of_overlap_%>'
    try:
        opts, arg = getopt.getopt(argv,"hf:r:s:i:l:",["fumdfile","radiusfile","sSampling_Frequency","iInitialStep","level_overlap"])
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('check_overlap.py : overlap of the core spheres (pseudopotential / PAW radii) of the atoms during a run. For each snapshot and each')
            print('pair of elements : the shortest distance (nearest images, cell of each snapshot) and the volume of the overlap of the two spheres')
            print('at that distance, in % of the sum of their volumes. A quick check of the pseudopotentials at high pressure, for any ab initio code.')
            print(usage)
            print('-f : the UMD file. Default output.umd.dat.')
            print('-r : radius file : one line per element, "element radius unit" (unit Angstrom or Bohr), e.g. "Si 1.9 Bohr" (RCORE of a VASP')
            print('     POTCAR, or the pseudopotential / PAW radius of another code). Default radius.inp.')
            print('-s : every s-th snapshot. Default 1.')
            print('-i : number of initial snapshots skipped. Default 0.')
            print('-l : acceptable overlap, in % of the sum of the volumes of the two spheres. Default 12.')
            print('Output : <base>.overlap.dat : the radii, then one line per snapshot : OK or NOT-ok (an overlap above -l), and for each pair of')
            print('elements the shortest distance (A) and the overlap (%). A summary per pair is printed.')
            sys.exit()
        elif opt in ("-f","--fumdfile"):
            umdfile = str(arg)
        elif opt in ("-r","--radiusfile"):
            radiusfile = str(arg)
        elif opt in ("-s","--sNsteps"):
            Nsteps = int(arg)
        elif opt in ("-i","--iInitialStep"):
            InitialStep = int(arg)
        elif opt in ("-l","--level_overlap"):
            level = float(arg)
    if not os.path.isfile(radiusfile):
        print('the radius file',radiusfile,'does not exist')
        sys.exit()
    if not os.path.isfile(umdfile):
        print('the umd file',umdfile,'does not exist')
        sys.exit()
    radius = read_radius(radiusfile)

    MyCrystal,Snapshots,TimeStep,length = umdpf.read_values(umdfile,"xred","line",Nsteps=Nsteps,firststep=InitialStep)
    cells = umdpf.read_cells(umdfile,len(Snapshots),Nsteps,InitialStep)
    missing = [el for el in MyCrystal.elements if el not in radius]
    if missing:
        print('no radius for',', '.join(missing),'in',radiusfile)
        sys.exit()
    typat = np.array(MyCrystal.typat)
    masks = [typat == it for it in range(MyCrystal.ntypat)]
    elements = MyCrystal.elements
    pairs = [(i,j) for i in range(MyCrystal.ntypat) for j in range(i,MyCrystal.ntypat)]
    labels = [elements[i]+'-'+elements[j] for i,j in pairs]
    volumes = [4.0/3.0*np.pi*(radius[elements[i]]**3+radius[elements[j]]**3) for i,j in pairs]

    base = umdfile[:-8] if umdfile.endswith('.umd.dat') else umdfile
    newfile = base+'.overlap.dat'
    dmins = []
    percents = []
    with open(newfile,'w') as nf:
        nf.write('# Overlap of the core spheres in '+umdfile+', every '+str(Nsteps)+' snapshots from '+str(InitialStep)+' ; acceptable overlap '+str(level)+' %\n')
        nf.write('# radii (A) : '+' '.join(el+' '+'%.4f' % radius[el] for el in elements)+'\n')
        nf.write('# overlap = volume of the intersection of the two spheres at the shortest distance, in % of the sum of their volumes\n')
        nf.write('snapshot\tstatus\t'+'\t'.join(lab+'_dmin(A)\t'+lab+'_overlap(%)' for lab in labels)+'\n')
        for k,(values,cell) in enumerate(zip(Snapshots,cells)):
            xred = np.array(values).reshape(MyCrystal.natom,3)
            dmin = shortest_distances(xred,cell,masks)
            d = [dmin[i,j] for i,j in pairs]
            pc = [100.0*lens_volume(dd,radius[elements[i]],radius[elements[j]])/v if np.isfinite(dd) else 0.0 for dd,(i,j),v in zip(d,pairs,volumes)]
            dmins.append(d)
            percents.append(pc)
            status = 'NOT-ok' if max(pc) > level else 'OK'
            nf.write(str(InitialStep+k*Nsteps)+'\t'+status+'\t'+'\t'.join('%.4f\t%.3f' % (dd,p) for dd,p in zip(d,pc))+'\n')
    dmins = np.array(dmins)
    percents = np.array(percents)
    print('\nOverlap of the core spheres over',len(Snapshots),'snapshots (every',Nsteps,'from',InitialStep,') ; acceptable overlap',level,'%')
    print('{:<10s}{:>12s}{:>12s}{:>14s}{:>12s}{:>16s}'.format('pair','r1+r2 (A)','dmin (A)','max overlap %','snapshot','snapshots > -l'))
    for p,(lab,(i,j)) in enumerate(zip(labels,pairs)):
        if not np.isfinite(dmins[:,p]).any():
            continue
        kmax = int(np.argmax(percents[:,p]))
        print('{:<10s}{:>12.3f}{:>12.3f}{:>14.2f}{:>12d}{:>16d}'.format(lab,radius[elements[i]]+radius[elements[j]],np.nanmin(np.where(np.isfinite(dmins[:,p]),dmins[:,p],np.nan)),
              percents[kmax,p],InitialStep+kmax*Nsteps,int((percents[:,p] > level).sum())))
    print('Overlaps written in',newfile)


if __name__ == "__main__":
    main(sys.argv[1:])
