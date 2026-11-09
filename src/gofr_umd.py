#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### Authors : Kevin Jiguet-Covex, Razvan Caracas
### One method for all the cells : nearest image in reduced coordinates, cell of every snapshot (NPT), processes around the C library : Razvan Caracas
### Pair distribution functions g(r) and running coordination numbers of all the pairs of elements, or of central and adjacent elements (-c, -a)


import sys, getopt, os.path
import time
import platform
import ctypes
from functools import partial
import numpy as np
import umd_processes_fast as umdpf


OS = platform.system()
if OS == "Linux":
    LibraryName = 'c_gofr.so'
elif OS == "Windows":
    LibraryName = 'c_gofr.dll'
else:
    LibraryName = 'c_gofr.dylib'
gofr_lib = ctypes.cdll.LoadLibrary(os.path.join(os.path.dirname(os.path.abspath(__file__)),LibraryName))
array_double = np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
gofr_lib.gofr_histogram.argtypes = (array_double,ctypes.c_int,ctypes.c_int,np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS'),ctypes.c_int,
                                    array_double,array_double,ctypes.c_double,ctypes.c_int,array_double,array_double)
gofr_lib.gofr_histogram.restype = None


def half_width(cells):#Half the smallest distance between opposite faces of each cell : the nearest image is exact up to it
    volumes = np.abs(np.linalg.det(cells))
    areas = np.array([np.linalg.norm(np.cross(cells[:,(k+1)%3],cells[:,(k+2)%3]),axis=1) for k in range(3)])
    return 0.5*(volumes/areas).min(axis=0)


def histograms(red,typat,ntypes,cells,weights,dr,nbins,counts=None,wcounts=None):#Pair counts, and counts weighted by the volume of their snapshot,
    #from the C library ; added to counts and wcounts when given (flat arrays ntypes*ntypes*nbins)
    if counts is None:
        counts = np.zeros(ntypes*ntypes*nbins)
        wcounts = np.zeros(ntypes*ntypes*nbins)
    gofr_lib.gofr_histogram(np.ascontiguousarray(red),len(red),red.shape[1],typat,ntypes,np.ascontiguousarray(cells),
                            np.ascontiguousarray(weights),dr,nbins,counts,wcounts)
    return counts.reshape(ntypes,ntypes,nbins), wcounts.reshape(ntypes,ntypes,nbins)


def block_histograms(block,umdfile,offsets,natom,Nsteps,keep,typat,ntypes,cells,invcells,volumes,dr,nbins,batch):#The histograms of a block of
    #consecutive snapshots (indexes of the selection), read from the UMD file batch snapshots at a time and added in order : the memory holds one batch
    counts = np.zeros(ntypes*ntypes*nbins)
    wcounts = np.zeros(ntypes*ntypes*nbins)
    for a in range(block[0],block[-1]+1,batch):
        b = min(a+batch,block[-1]+1)
        X = umdpf.read_atoms_at(umdfile,offsets[a],natom,b-a,Nsteps,3,3)[:,keep]
        red = np.einsum('sni,sij->snj',X,invcells[a:b])#x = s M, with the lattice vectors as rows of M
        histograms(red,typat,ntypes,cells[a:b],volumes[a:b],dr,nbins,counts,wcounts)
    return counts.reshape(ntypes,ntypes,nbins), wcounts.reshape(ntypes,ntypes,nbins)


def main(argv):
    umdpf.headerumd()
    start = time.time()
    umdfile = ''
    Nsteps = 1
    discrete = 0.01
    InitialStep = 0
    nCores = None
    central = 'all'
    adjacent = 'all'
    usage = 'gofr_umd.py -f <umdfile> -s <Sampling_Frequency> -d <discretization.interval> -i <InitialStep> -c <central> -a <adjacent> -k <nCores>'
    try:
        opts, arg = getopt.getopt(argv,"hf:s:d:i:c:a:k:",["fumdfile","sSampling_Frequency","ddiscrete","iInitialStep","ccentral","aadjacent","knCores"])
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('Program to compute the pair distribution functions g(r) and the running coordination numbers of all the pairs of elements. Usage:')
            print(usage)
            print('-f : UMD file.')
            print('-f all : every *.umd.dat file of the current folder, one after the other.')
            print('-s : Sampling_Frequency. one snapshot out of s is used. Default is 1.')
            print('-d : width of the bins of the g(r), in A. Default is 0.01.')
            print('-i : InitialStep. number of initial snapshots discarded. Default is 0.')
            print('-c : central elements, separated by commas, e.g. O or Si,Mg. Default is all.')
            print('-a : adjacent elements, separated by commas, e.g. O. Default is all.')
            print('     The pairs computed are central-adjacent and their reverse (same g(r), the coordination seen from the other atom), as in')
            print('     analyze_gofr.py : -c O gives O-X and X-O for every element X ; -c Si,Mg -a O gives Si-O, Mg-O, O-Si, O-Mg. All the other')
            print('     pairs are written with g(r) = 0 and Int = 0, so analyze_gofr.py gives them a bond length of 0, i.e. no bond.')
            print('     Only the atoms of the central and adjacent elements enter the computation, so a subset is faster than all the pairs.')
            print('-k : number of cores. Default is all.')
            print('Any cell, orthogonal or not, and cells changing at every snapshot (NPT) : the nearest image of each pair is taken in')
            print('reduced coordinates with the cell of its snapshot, which is exact up to half the smallest distance between opposite faces')
            print('of the cell, the largest distance computed. g(r) is normalized with the volume of each snapshot and tends to 1 at large r.')
            print('Output : <base>.gofr.dat : dist (center of the bin, A), then for each pair A-B its g(r) and Int(A-B), the mean number of')
            print('         B atoms around an A atom up to the outer edge of the bin.')
            print('         With -c or -a : <base>.<central>-<adjacent>.gofr.dat, the elements written one after the other, e.g. -c Si,Mg -a O :')
            print('         <base>.SiMg-O.gofr.dat, -c O : <base>.O-all.gofr.dat ; all the pairs, in the order of the elements of the UMD file.')
            print('         The file of a former run with the same selection is removed first ; those of other selections are kept.')
            sys.exit()
        elif opt in ("-f","--fumdfile"):
            umdfile = str(arg)
        elif opt in ("-s","--sSampling_Frequency"):
            Nsteps = int(arg)
        elif opt in ("-d","--ddiscrete"):
            discrete = float(arg)
        elif opt in ("-i","--iInitialStep"):
            InitialStep = int(arg)
        elif opt in ("-c","--ccentral"):
            central = arg
        elif opt in ("-a","--aadjacent"):
            adjacent = arg
        elif opt in ("-k","--knCores"):
            nCores = int(arg)
    if not os.path.isfile(umdfile):
        umdpf.missing_umd(umdfile)
        sys.exit()
    MyCrystal = umdpf.Crystallization(umdfile)[0]
    selected = {}
    for name, value in (('central',central),('adjacent',adjacent)):
        if value == 'all':
            selected[name] = list(MyCrystal.elements)
        else:
            selected[name] = [el.strip() for el in value.split(',') if el.strip()]
            unknown = [el for el in selected[name] if el not in MyCrystal.elements]
            if unknown or not selected[name]:
                print('the',name,'elements',value,'are not all in the UMD file, whose elements are',','.join(MyCrystal.elements))
                sys.exit(2)

    gofrname = umdfile[:-8] if umdfile.endswith('.umd.dat') else os.path.splitext(umdfile)[0]
    if central != 'all' or adjacent != 'all':
        gofrname += '.'+('all' if central == 'all' else ''.join(selected['central']))+'-'+('all' if adjacent == 'all' else ''.join(selected['adjacent']))
    gofrname += '.gofr.dat'
    umdpf.remove_former_outputs([gofrname])#the file of a former run with the same selection ; the gofr files of other selections are kept
    natom = MyCrystal.natom
    offsets, cells = umdpf.index_snapshots(umdfile,natom)#byte offset and cell of each snapshot (NPT)
    offsets, cells = offsets[InitialStep::Nsteps], cells[InitialStep::Nsteps]#the snapshots used : from InitialStep, every Nsteps
    nsnap = len(offsets)
    if nsnap == 0:
        print('No snapshot after the snapshot',InitialStep)
        sys.exit()
    if np.isnan(cells).any():
        print('ERROR : snapshots of',umdfile,'without rprimd or acell lines')
        sys.exit()
    #Only the atoms of the central and adjacent elements are kept ; all the atoms of an element are, so the normalization is unchanged
    elements = MyCrystal.elements
    used = [i for i in range(MyCrystal.ntypat) if elements[i] in selected['central'] or elements[i] in selected['adjacent']]
    typat_all = np.array(MyCrystal.typat)
    keep = np.isin(typat_all,used)
    volumes = np.abs(np.linalg.det(cells))
    rmax = half_width(cells).min()
    nbins = int(rmax/discrete)
    ntypes = len(used)
    typat = np.searchsorted(used,typat_all[keep]).astype(np.int32)#types renumbered 0..ntypes-1 in the order of the UMD file
    counts_el = np.bincount(typat,minlength=ntypes)
    #All the pairs of the UMD file are written, in its order ; those not central-adjacent (or the reverse) are 0, i.e. no bond for analyze_gofr.py
    pairlist = [(elements[I]+'-'+elements[J],
                 (used.index(I),used.index(J)) if (elements[I] in selected['central'] and elements[J] in selected['adjacent']) or
                                                 (elements[J] in selected['central'] and elements[I] in selected['adjacent']) else None)
                for I in range(MyCrystal.ntypat) for J in range(MyCrystal.ntypat)]
    print(nsnap,' snapshots ; volume ','%.3f' % volumes.mean(),' A^3'+('' if np.ptp(volumes) < 1.0E-6*volumes.mean() else ' (from %.3f to %.3f)' % (volumes.min(),volumes.max()))
          +' ; g(r) up to ','%.4f' % (nbins*discrete),' A, half the smallest distance between opposite faces of the cell')

    #The snapshots are cut into blocks computed in parallel processes, each reading its block from the UMD file in batches sized from the bytes
    #of a snapshot (the memory holds one batch per process). Not threads : the C reader (strtod) then spends its time on a lock of the system
    nworkers = max(1,min(nCores if nCores is not None else (os.cpu_count() or 1),nsnap))
    blocks = [block for block in np.array_split(np.arange(nsnap),nworkers) if len(block)]
    snapbytes = (os.path.getsize(umdfile)-offsets[0])/nsnap
    batch = max(1,int(umdpf.BATCH_BYTES//max(1,snapbytes)))
    worker = partial(block_histograms,umdfile=umdfile,offsets=offsets,natom=natom,Nsteps=Nsteps,keep=keep,typat=typat,ntypes=ntypes,cells=cells,
                     invcells=np.linalg.inv(cells),volumes=volumes,dr=discrete,nbins=nbins,batch=batch)
    results = list(umdpf.ordered_map(worker,blocks,nworkers))
    counts = sum(r[0] for r in results)
    wcounts = sum(r[1] for r in results)

    #g_ij(r) = < V n_ij(r) > / (N_i (N_j - delta_ij) 4/3 pi ((r+dr)^3 - r^3)), n_ij counting the ordered pairs ; Int_ij = cumulated n_ij / N_i
    edges = discrete*np.arange(nbins+1)
    shells = 4.0/3.0*np.pi*(edges[1:]**3-edges[:-1]**3)
    pairs = counts_el[:,None]*(counts_el[None,:]-np.eye(ntypes))
    with np.errstate(divide='ignore',invalid='ignore'):
        g = np.where(pairs[:,:,None] > 0,wcounts/(nsnap*pairs[:,:,None]*shells),0.0)
    coord = np.cumsum(counts,axis=-1)/(nsnap*counts_el[:,None,None])

    with open(gofrname,'w') as fa:
        fa.write('dist\t'+'\t'.join(pair+'\tInt('+pair+')' for pair, ij in pairlist)+'\n')
        for k in range(nbins):
            fa.write(str(round(edges[k]+0.5*discrete,5))+'\t'+'\t'.join(str(g[ij[0],ij[1],k])+'\t'+str(coord[ij[0],ij[1],k]) if ij is not None else '0.0\t0.0'
                                                                    for pair, ij in pairlist)+'\n')
    print('gofr file successfully created under the name ',gofrname)
    print('runtime : ','%.1f' % (time.time()-start),' s')


if __name__ == "__main__":
    if not umdpf.run_on_all(sys.argv[1:],main):#-f all : every UMD file of the folder
        main(sys.argv[1:])
