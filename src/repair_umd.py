#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### Authors : Razvan Caracas
### Repairs the velocities, unwrapped positions and reduced coordinates of a UMD file written by an older VaspParser.py,
### recomputing them from the wrapped cartesian positions exactly as the current VaspParser.py does


import sys, getopt, os.path
import time
import numpy as np
import crystallography as cr
import umd_processes_fast as umdpf


def reciprocal(rprimd):#gprimd of crystallography.Lattice, the inverse of the cell used for the reduced coordinates
    cell = cr.Lattice()
    cell.rprimd = rprimd
    return np.array(cell.makegprimd())


def reduced(g,x):#sum_jj g[ii][jj] x[jj], summed from 0.0 in the same order as in VaspParser.py so that the results are identical
    return np.column_stack([0.0+g[ii][0]*x[:,0]+g[ii][1]*x[:,1]+g[ii][2]*x[:,2] for ii in range(3)])


def fmt(x):#the format of umd_processes_fast.print_snapshots : str(round(x,5)) of a python float
    return str(round(float(x),5))


def repair(umdfile,outfile):
    #Velocities, unwrapped positions (abs.diff.) and reduced coordinates are recomputed from the cartesian positions xcart, as in VaspParser.py :
    #v(n) = (x(n) - x(n-1) + phase)/dt, where a jump across a cell face is an integer change of the reduced displacement ; the first snapshot gets v = 0
    #Everything else (energies, cells, forces, charges ...) is copied unchanged
    natom = None
    rprimd = [None,None,None]
    timestep = None
    times = []
    xold = None
    phase = None
    nsnap = 0
    nzero = 0
    ncross = 0
    nmoved = 0
    nchanged = 0
    skip = None
    fin = open(umdfile,'r')
    fout = open(outfile,'w')
    for line in fin:
        entry = line.split()
        if len(entry) > 0:
            if entry[0] == 'natom':
                natom = int(entry[1])
            elif entry[0] == 'timestep':
                timestep = float(entry[1])
            elif entry[0] == 'time':
                times.append(float(entry[1]))
            elif entry[0] in ('rprimd_a','rprimd_b','rprimd_c'):
                rprimd['abc'.index(entry[0][-1])] = [float(entry[1]),float(entry[2]),float(entry[3])]
        fout.write(line)
        if len(entry) == 0 or entry[0] != 'atoms:':
            continue
        tokens = [fin.readline().split() for _ in range(natom)]
        data = np.array([[float(x) for x in tok[:15]] for tok in tokens])
        x = data[:,3:6]
        g = reciprocal(rprimd)
        if xold is None:#first snapshot : no previous position, no phase, zero velocity
            phase = np.zeros((natom,3))
            vels = np.zeros((natom,3))
        else:
            dt = times[-1]-times[-2]
            if skip is None:
                skip = dt/timestep
            if abs(dt-timestep) < 1.0E-6*timestep:#one MD step apart : the MD step, as VaspParser.py
                dt = timestep
            dcart = x-xold
            addphase = -np.rint(reduced(g,dcart))
            r = np.array(rprimd)
            addphase_cart = np.column_stack([0.0+addphase[:,0]*r[0][ii]+addphase[:,1]*r[1][ii]+addphase[:,2]*r[2][ii] for ii in range(3)])
            phase = phase+addphase_cart
            vels = (x-xold+addphase_cart)/dt
            ncross += int(np.count_nonzero(addphase))
        xred = reduced(g,x)
        xred -= np.floor(xred)
        absx = x+phase
        if not data[:,9:12].any():
            nzero += 1
        elif np.abs(vels-data[:,9:12]).max() > 2.0E-5:#beyond the rounding to 5 decimals
            nchanged += 1
        if np.abs(absx-data[:,6:9]).max() > 1.0E-4:
            nmoved += 1
        for iatom in range(natom):
            fout.write(' '.join(fmt(v) for v in xred[iatom])+' '+' '.join(tokens[iatom][3:6])+' '+' '.join(fmt(v) for v in absx[iatom])+' '
                       +' '.join(fmt(v) for v in vels[iatom])+' '+' '.join(tokens[iatom][12:])+'\n')
        xold = x
        nsnap += 1
    fin.close()
    fout.close()
    return nsnap, nzero, nchanged, ncross, nmoved, skip


def main(argv):
    umdpf.headerumd()
    start = time.time()
    umdfile = ''
    outfile = None
    usage = 'repair_umd.py -f <UMD_filename> -o <output_filename>'
    try:
        opts, arg = getopt.getopt(argv,"hf:o:",["fumdfile","ooutfile"])
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('Program to repair a UMD file written by VaspParser.py before October 23, 2026, without its OUTCAR. Usage: ')
            print(usage)
            print('Those versions set the velocities to zero, and missed the crossings of the cell faces, on the learning steps of on-the-fly')
            print('ML runs (the OUTCAR prints the positions twice), tested the crossings in a way valid only for orthogonal cells, wrote the')
            print('positions instead of the velocities on the first snapshot, and computed wrong reduced coordinates in non-orthogonal cells.')
            print('The cartesian positions xcart are right (they are those of the OUTCAR) : the velocities, the unwrapped positions (abs.diff.)')
            print('and the reduced coordinates are recomputed from them, exactly as the current VaspParser.py does. Everything else is copied.')
            print('The first snapshot gets zero velocities (no previous position ; vibrations_umd.py discards it).')
            print('If the UMD file was written every s steps (-s s), the velocities can only be the mean velocities over the s steps between')
            print('snapshots : fine for the unwrapped positions (MSD), not for the VACF and the vDOS. The file has to be parsed again for them.')
            print('-f : UMD file to repair.')
            print('-f all : every *.umd.dat file of the current folder, one after the other (each into <base>.repaired.umd.dat ; -o cannot be given then).')
            print('-o : output file. Default is <base>.repaired.umd.dat.')
            sys.exit()
        elif opt in ("-f","--fumdfile"):
            umdfile = str(arg)
        elif opt in ("-o","--ooutfile"):
            outfile = str(arg)
    if not os.path.isfile(umdfile):
        umdpf.missing_umd(umdfile)
        sys.exit()
    if outfile is None:
        base = umdfile[:-8] if umdfile.endswith('.umd.dat') else os.path.splitext(umdfile)[0]
        outfile = base+'.repaired.umd.dat'
    if os.path.abspath(outfile) == os.path.abspath(umdfile):
        print('ERROR : the output file cannot be the UMD file itself.')
        sys.exit()

    print('Repairing ',umdfile)
    nsnap, nzero, nchanged, ncross, nmoved, skip = repair(umdfile,outfile)
    print(nsnap,' snapshots')
    print(nzero,' snapshots had all their velocities zero ; ',nchanged,' other snapshots have changed velocities')
    print(ncross,' crossings of a cell face ; ',nmoved,' snapshots with changed unwrapped positions')
    if nchanged > 0.1*(nsnap-nzero):
        print('WARNING : the velocities of most snapshots changed : they were not the differences of consecutive positions. If the file')
        print('          was written by the current VaspParser.py with -s > 1, its velocities were right and should not be replaced.')
    if skip is not None and abs(skip-1.0) > 1.0E-6:
        print('WARNING : the snapshots are ','%g' % skip,' MD steps apart : the velocities are mean velocities over these steps, valid for')
        print('          the unwrapped positions (MSD) but not for the VACF and the vDOS. Parse the OUTCAR again with VaspParser.py for them.')
        print('          The timestep written in the file is still the MD step, while the analysis scripts expect the time between snapshots.')
    print('Repaired UMD file written in ',outfile)
    print('runtime : ','%.1f' % (time.time()-start),' s')


if __name__ == "__main__":
    if not umdpf.run_on_all(sys.argv[1:],main,single=('-o',)):#-f all : every UMD file of the folder
        main(sys.argv[1:])
