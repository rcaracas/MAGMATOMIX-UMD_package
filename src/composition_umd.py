#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### Author : Razvan Caracas
### Earlier attempts, never in the package : bubbles.py (2017 : melt and gas from the clusters of a popul file), densitytracker.py (Tim Bogels,
### 2020 : the mass in a grid of cubes over time)
### Compositional fluctuations : the atoms of each element counted in a grid of sub-cells at every snapshot ; distributions of the local
### composition and density against ideal mixing, fluctuations and correlations of the elements, voids (bubbles, or the gas around a droplet)
### and dense regions, and how long the heterogeneities live


import sys, getopt, os.path, re, heapq
import time
from functools import partial
import numpy as np
from scipy import ndimage
from scipy.stats import binom, hypergeom
import umd_processes_fast as umdpf
import averages


def statistics(x):#averages.statistics (mean, stdev, error, lower bound, drift), a series constant to rounding taken as constant (error 0) :
    #e.g. the variance of the count of an element with a single atom, binomial at every snapshot
    x = np.asarray(x,dtype=float)
    if np.ptp(x) <= 1.0E-9*max(abs(x).max(),1.0E-300):
        return x.mean(), 0.0, 0.0, False, 0.0
    return averages.statistics(x)


AMU_A3 = 1.66053906660                  #g/cm3 for 1 amu/A^3
VOLATILES = ['H','He','C','N','F','Ne','S','Cl','Ar','Kr','Xe']   #the volatile elements by default, when present in the run
COLORS = ['#2a78d6','#eb6834','#1baf7a','#eda100','#e87ba4','#008300','#6250d6','#e34948']#one per element, in the order of the UMD file (as lifetimes.py)
OTHER = '#a3a29c'               #the elements beyond the eighth, together
TEXT, MUTED, SURFACE, VOIDSHADE = '#0b0b0b', '#52514e', '#fcfcfb', '#e4e3df'


def widths(cell):#Distances between the opposite faces of a cell (rows a, b, c)
    volume = abs(np.linalg.det(cell))
    return np.array([volume/np.linalg.norm(np.cross(cell[1],cell[2])),volume/np.linalg.norm(np.cross(cell[2],cell[0])),
                     volume/np.linalg.norm(np.cross(cell[0],cell[1]))])


def count_batch(task,umdfile,offsets,natom,sampling,grid,typat,nel,keepidx=False):#Counts of the atoms of each element in each sub-cell, snapshots a to b of
    #the selection : array (snapshot, sub-cell, element). Sub-cells in reduced coordinates, so the same fraction of the cell for any shape and
    #at every snapshot (NPT) ; sub-cell index (ia nb + ib) nc + ic
    a, b = task
    X = umdpf.read_atoms_at(umdfile,offsets[a],natom,b-a,sampling,0,3)         #reduced coordinates
    s = X-np.floor(X)
    idx = np.minimum((s*grid).astype(np.int64),grid-1)
    ncell = int(np.prod(grid))
    cellidx = (idx[...,0]*grid[1]+idx[...,1])*grid[2]+idx[...,2]
    flat = (np.arange(b-a)[:,None]*ncell+cellidx)*nel+typat[None,:]
    counts = np.bincount(flat.ravel(),minlength=(b-a)*ncell*nel).reshape(b-a,ncell,nel).astype(np.int32)
    return (counts, cellidx.astype(np.int32)) if keepidx else counts       #with keepidx the sub-cell of every atom (exchange gas-melt)


def periodic_regions(mask):#Connected regions of the True sub-cells of a grid (face neighbours, across the periodic faces too) : the region of
    #every sub-cell (0 outside the mask, else 1 .. number of regions) and the number of regions
    lab, n = ndimage.label(mask)
    if n == 0:
        return lab, 0
    parent = np.arange(n+1)
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for ax in range(3):
        if mask.shape[ax] < 2:
            continue
        first, last = np.take(lab,0,axis=ax), np.take(lab,-1,axis=ax)
        both = (first > 0)&(last > 0)
        for i, j in set(zip(first[both].tolist(),last[both].tolist())):
            ri, rj = find(i), find(j)
            if ri != rj:
                parent[max(ri,rj)] = min(ri,rj)
    roots = np.array([find(i) for i in range(n+1)])
    used, renum = np.unique(roots[1:],return_inverse=True)
    region = np.zeros(n+1,dtype=np.int64)
    region[1:] = renum+1
    return region[lab], len(used)


def acf(x,nlag):#Autocorrelation in time (axis 0) of each column of x, averaged over the columns, unbiased (each lag over its own number of
    #origins), lags 0 to nlag-1 ; the mean is not removed here (deviations from a global mean, so that a heterogeneity that stays is seen)
    n = x.shape[0]
    size = 1<<int(np.ceil(np.log2(2*n)))
    f = np.fft.rfft(x,n=size,axis=0)
    c = np.fft.irfft((f*f.conj()).real,n=size,axis=0)[:nlag].mean(axis=1)
    return c/(n-np.arange(nlag))


def decay_time(c,dt):#Time at which a normalized correlation first falls below 1/e (linear interpolation), nan when it does not
    below = np.nonzero(c < 1.0/np.e)[0]
    if len(below) == 0 or below[0] == 0:
        return np.nan
    k = below[0]
    return dt*(k-1+(c[k-1]-1.0/np.e)/(c[k-1]-c[k]))


def unbiased(cf,n):#Normalized correlation of series whose own mean was removed (cf, lags from 0) : removing the mean of a finite series lowers
    #every lag by the variance of that mean, c0 2 tau / n ; it is added back (tau and the correction found together, a few iterations)
    if cf[0] <= 0:
        return np.full(len(cf),np.nan)
    V = 0.0
    for _ in range(5):
        rho = (cf+V)/(cf[0]+V)
        cut = np.argmax(rho < 0.0) or len(rho)
        tau = 0.5+rho[1:cut].sum()
        V = (cf[0]+V)*2.0*tau/n
    return (cf+V)/(cf[0]+V)


def correlations(series,dt,nblocks=4):#For each series (dictionary name : array time x sub-cells, deviations from the mean of the run), up to half the
    #run : the normalized correlation in time (averaged over the sub-cells) and its 1/e time ; the persistent fraction, the part of the variance
    #that stays in place over the run, (mean over the sub-cells of their mean^2 - its noise) / variance, the noise of the mean of a sub-cell being
    #its variance x 2 tau / T (tau : integrated correlation time of the fluctuations about its own mean) ; the 1/e time of these fluctuations, which
    #decay even when the heterogeneity stays ; the correlation at the end (T/4 to T/2), which tends to the persistent fraction. Errors : 1/e
    #times from nblocks blocks, the persistent fraction from its spread over the sub-cells
    nsnap = next(iter(series.values())).shape[0]
    nlag = max(2,nsnap//2)
    curves, times = {}, {}
    for name, x in series.items():
        c = acf(x,nlag)
        curves[name] = c/c[0] if c[0] > 0 else np.full(nlag,np.nan)
        mean = x.mean(axis=0)                                           #per sub-cell
        y = x-mean
        cf = acf(y,nlag)
        rho = unbiased(cf,nsnap)
        blocklen = nsnap//nblocks
        tb, tf = [], []
        for k in range(nblocks):
            xb = x[k*blocklen:(k+1)*blocklen]
            cb = acf(xb,max(2,blocklen//2))
            tb.append(decay_time(cb/cb[0],dt) if cb[0] > 0 else np.nan)
            yb = xb-xb.mean(axis=0)
            cfb = acf(yb,max(2,blocklen//2))
            tf.append(decay_time(unbiased(cfb,blocklen),dt) if cfb[0] > 0 else np.nan)
        tb, tf = np.array(tb), np.array(tf)
        total = (x**2).mean()
        if total > 0 and cf[0] > 0:
            cut = np.argmax(rho < 0.0) or len(rho)
            tau = 0.5+rho[1:cut].sum()                                  #in snapshots
            q = mean**2-y.var(axis=0)*2.0*tau/nsnap                     #per sub-cell : its static part, noise removed
            persistent, perr = q.mean()/total, q.std()/np.sqrt(len(q))/total
        else:
            persistent, perr = np.nan, np.nan
        end = curves[name][nlag//2:].mean()
        times[name] = (decay_time(curves[name],dt),np.std(tb)/np.sqrt(nblocks) if np.all(np.isfinite(tb)) else np.nan,
                       persistent,perr,decay_time(rho,dt),np.std(tf)/np.sqrt(nblocks) if np.all(np.isfinite(tf)) else np.nan,end)
    return curves, times


def parse_cells(arg,grid,rng):#Sub-cells of -c : ia:ib:ic (indexes from 0 along a, b, c), several separated by commas (1:1:1,5:2:2), r (one at
    #random), rN (N at random, e.g. r5), all (every sub-cell, in order)
    allcells = [(i,j,k) for i in range(grid[0]) for j in range(grid[1]) for k in range(grid[2])]
    if arg == 'all':
        return allcells
    found = re.fullmatch(r'r(\d*)',arg)
    if found:
        n = min(int(found.group(1) or 1),len(allcells))
        return [allcells[m] for m in rng.choice(len(allcells),n,replace=False)]
    cells = []
    for part in arg.split(','):
        try:
            cell = tuple(int(x) for x in part.split(':'))
        except ValueError:
            cell = ()
        if len(cell) != 3 or any(c < 0 or c >= g for c, g in zip(cell,grid)):
            print('-c : ia:ib:ic with 0 <= ia <',grid[0],', 0 <= ib <',grid[1],', 0 <= ic <',grid[2],', several separated by commas ;',
                  'or r (one at random), r5 (five at random), all')
            sys.exit(2)
        if cell not in cells:
            cells.append(cell)
    return cells


def plot_subcells(cells,grid,times,C,elements,isvol,rho_all,rho_melt,limit,title,tracked=None):#For each sub-cell (ia, ib, ic) : the number of atoms of the
    #tracked elements (default all ; -a)
    #against time, stacked (the top line is the total), and below the density of the sub-cell with the void threshold, the voids shaded ;
    #Previous and Next buttons when there are several sub-cells (Next on the last one closes the window, as the other scripts). Returns the
    #figure and the function that draws sub-cell k (for tests)
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Button
    nel = len(elements)
    tracked = list(range(nel)) if tracked is None else list(tracked)
    order = [e for e in tracked if not isvol[e]]+[e for e in tracked if isvol[e]]   #melt elements at the bottom, volatiles on top
    shown = order if len(order) <= len(COLORS) else order[:len(COLORS)-1]
    colors = {e:COLORS[e] for e in shown} if nel <= len(COLORS) else {e:COLORS[k] for k, e in enumerate(shown)}
    rest = [e for e in order if e not in shown]
    fig = plt.figure(figsize=(10,6.5),facecolor=SURFACE)
    ax1 = fig.add_axes([0.08,0.47,0.74,0.42],facecolor=SURFACE)
    ax2 = fig.add_axes([0.08,0.14,0.74,0.27],facecolor=SURFACE,sharex=ax1)
    state = {'k':0}
    def draw(k):
        ia, ib, ic = cells[k]
        idx = (ia*grid[1]+ib)*grid[2]+ic
        for ax in (ax1,ax2):
            ax.clear()
            ax.set_facecolor(SURFACE)
            for side in ('top','right'):
                ax.spines[side].set_visible(False)
            for side in ('left','bottom'):
                ax.spines[side].set_color(MUTED)
            ax.tick_params(colors=MUTED,labelsize=8,length=2)
            ax.grid(axis='y',color=VOIDSHADE,linewidth=0.6)
            ax.set_axisbelow(True)
        layers = [C[:,idx,e] for e in shown]+([C[:,idx,rest].sum(axis=1)] if rest else [])
        labels = [elements[e] for e in shown]+(['other ('+','.join(elements[e] for e in rest)+')'] if rest else [])
        cols = [colors[e] for e in shown]+([OTHER] if rest else [])
        ax1.stackplot(times,layers,colors=cols,labels=labels,edgecolor=SURFACE,linewidth=0.4,step='post')    #counts : steps
        ax1.set_ylabel('atoms in the sub-cell'+('' if len(tracked) == nel else ' ('+','.join(elements[e] for e in order)+')'),color=TEXT,fontsize=9)
        ax1.set_ylim(bottom=0)
        handles, names = ax1.get_legend_handles_labels()
        means = [C[:,idx,e].mean() for e in shown]+([C[:,idx,rest].sum(axis=1).mean()] if rest else [])
        ax1.legend(handles[::-1],[n+'  '+'%.2f' % m for n, m in zip(names,means)][::-1],loc='upper left',bbox_to_anchor=(1.01,1.0),frameon=False,
                   fontsize=8,labelcolor=TEXT,title='element, mean',title_fontsize=8)
        void = rho_melt[:,idx] < limit
        ax2.fill_between(times,0,1,where=void,transform=ax2.get_xaxis_transform(),color=VOIDSHADE,linewidth=0,label='void')
        ax2.plot(times,rho_all[:,idx],color=TEXT,linewidth=1.0,label='all atoms')
        ax2.plot(times,rho_melt[:,idx],color=MUTED,linewidth=1.0,linestyle=(0,(4,2)),label='melt elements')
        ax2.axhline(limit,color=MUTED,linewidth=0.8,linestyle=(0,(1,2)),label='void threshold')
        ax2.set_ylim(bottom=0)
        ax2.set_xlim(times[0],times[-1])
        ax2.set_ylabel('density (g/cm$^3$)',color=TEXT,fontsize=9)
        ax2.set_xlabel('time (fs)',color=TEXT,fontsize=9)
        ax2.legend(loc='upper left',bbox_to_anchor=(1.01,1.0),frameon=False,fontsize=8,labelcolor=TEXT)
        plt.setp(ax1.get_xticklabels(),visible=False)
        fig.suptitle(title+' : sub-cell %d:%d:%d of %dx%dx%d' % (ia,ib,ic,*grid)+('   (%d of %d)' % (k+1,len(cells)) if len(cells) > 1 else '')
                     +'\nvoid %.0f %% of the time' % (100*void.mean()),color=TEXT,fontsize=10,x=0.08,ha='left')
        fig.canvas.draw_idle()
    def step(delta):
        k = state['k']+delta
        if k >= len(cells):
            plt.close(fig)
            return
        state['k'] = max(0,k)
        draw(state['k'])
    draw(0)
    if len(cells) > 1:
        bprev = Button(fig.add_axes([0.62,0.02,0.1,0.05]),'Previous')
        bnext = Button(fig.add_axes([0.73,0.02,0.1,0.05]),'Next')
        bprev.on_clicked(lambda event: step(-1))
        bnext.on_clicked(lambda event: step(1))
        fig._umd_buttons = (bprev,bnext)       #kept alive, otherwise garbage-collected and unresponsive
    return fig, draw


def phase_runs(series,shortest):#Periods of an atom in a phase (series of booleans in time) : [[phase, length], ...], the inner periods shorter than
    #shortest snapshots merged into the surrounding phase (the shortest first, then the earliest), so that crossings back and forth at the surface
    #are not exchanges. A heap over a linked list : N log N in the number of periods (merging one at a time from a list was N^2 : hours for an atom
    #flickering at the surface of a long run)
    change = np.nonzero(np.diff(series.astype(np.int8)))[0]+1
    starts = np.concatenate([[0],change])
    length = np.diff(np.concatenate([starts,[len(series)]])).tolist()
    phase = [bool(series[st]) for st in starts]
    n = len(length)
    if shortest > 1 and n > 2:
        prev = list(range(-1,n-1))
        nxt = list(range(1,n+1))
        nxt[-1] = -1
        alive = [True]*n
        heap = [(length[k],k) for k in range(1,n-1) if length[k] < shortest]
        heapq.heapify(heap)
        while heap:
            l, k = heapq.heappop(heap)
            if not alive[k] or l != length[k] or prev[k] < 0 or nxt[k] < 0:
                continue                                                #merged since, or now at an end
            p, q = prev[k], nxt[k]
            length[p] += length[k]+length[q]
            alive[k] = alive[q] = False
            nxt[p] = nxt[q]
            if nxt[q] >= 0:
                prev[nxt[q]] = p
            if prev[p] >= 0 and nxt[p] >= 0 and length[p] < shortest:
                heapq.heappush(heap,(length[p],p))
        return [[phase[k],length[k]] for k in range(n) if alive[k]]
    return [[ph,l] for ph, l in zip(phase,length)]


def track_regions(labels,nlabels,C,subvol,dt):#Connected regions (labels : snapshot x sub-cell, 0 outside, 1..nlabels[t]) followed in time by the
    #sub-cells they share : the largest overlaps first, a region of the next snapshot continues at most one of this snapshot and the reverse ; not
    #continued : merged (it overlaps a region of the next snapshot) or vanished ; new regions : split (they overlap one of this snapshot) or new.
    #Returns one record per region : start, end (snapshots), birth, death, the region it split from or merged into, volume and atoms per element
    nsnap, ncell = labels.shape
    nel = C.shape[2]
    regions = []
    current = {}                                                        #label at t : index in regions
    def new(t,birth,parent=-1):
        regions.append({'start':t,'end':t,'birth':birth,'death':'end','parent':parent,'into':-1,'vol':0.0,'volmax':0.0,'atoms':np.zeros(nel),'n':0})
        return len(regions)-1
    def accumulate(t,mapping):
        lab = labels[t]
        size = np.bincount(lab,minlength=nlabels[t]+1)
        atoms = np.array([np.bincount(lab,weights=C[t,:,e],minlength=nlabels[t]+1) for e in range(nel)]).T
        for l, r in mapping.items():
            g = regions[r]
            g['end'] = t
            v = size[l]*subvol[t]
            g['vol'] += v
            g['volmax'] = max(g['volmax'],v)
            g['atoms'] += atoms[l]
            g['n'] += 1
    current = {l:new(0,'start') for l in range(1,nlabels[0]+1)}
    accumulate(0,current)
    for t in range(1,nsnap):
        a, b = labels[t-1], labels[t]
        both = (a > 0)&(b > 0)
        pairs, counts = np.unique(np.stack([a[both],b[both]]),axis=1,return_counts=True) if both.any() else (np.zeros((2,0),dtype=int),np.zeros(0,dtype=int))
        order = sorted(range(len(counts)),key=lambda k: (-counts[k],pairs[0][k],pairs[1][k]))
        nxt, used = {}, set()
        for k in order:
            la, lb = int(pairs[0][k]), int(pairs[1][k])
            if lb not in nxt and la not in used:
                nxt[lb] = current[la]
                used.add(la)
        best_prev = {}                                                  #for each region of t, the region of t-1 it overlaps most, and the reverse
        best_next = {}
        for k in order:
            la, lb = int(pairs[0][k]), int(pairs[1][k])
            best_prev.setdefault(lb,la)
            best_next.setdefault(la,lb)
        for lb in range(1,nlabels[t]+1):
            if lb not in nxt:
                nxt[lb] = new(t,'split',current[best_prev[lb]]) if lb in best_prev else new(t,'new')
        for la, r in current.items():
            if la not in used:
                if la in best_next:
                    regions[r]['death'] = 'merged'
                    regions[r]['into'] = nxt[best_next[la]]
                else:
                    regions[r]['death'] = 'vanished'
        current = nxt
        accumulate(t,current)
    return regions


def quarters(x):#Means of the four quarters of a series, each with its error (autocorrelation time within the quarter), and the largest separation
    #of two quarters in their combined errors : a change during the run (segregation in progress) even when abrupt, which inflates the error of the
    #whole run and hides itself in the drift of averages.statistics
    blocks = np.array_split(np.asarray(x,dtype=float),4)
    st = [statistics(b) if len(b) >= 4 else (b.mean(),0.0,0.0,True,0.0) for b in blocks]
    means = np.array([v[0] for v in st])
    errs = np.array([v[2] for v in st])
    lower = any(v[3] for v in st)                                      #an error of a quarter only a lower bound : the separation is not reliable
    sep = 0.0
    for i in range(4):
        for j in range(i+1,4):
            comb = np.hypot(errs[i],errs[j])
            if comb > 0:
                sep = max(sep,abs(means[i]-means[j])/comb)
            elif means[i] != means[j]:
                sep = np.inf
    return means, errs, sep, lower


def episodes(state):#Periods of consecutive snapshots of each sub-cell in a state (array snapshot x sub-cell of booleans) : (durations in snapshots
    #of the complete periods, sub-cell of each, number of periods cut by the start or the end of the run)
    nsnap, ncell = state.shape
    padded = np.zeros((ncell,nsnap+2),dtype=np.int8)
    padded[:,1:-1] = state.T
    change = np.diff(padded,axis=1)
    sc, st = np.nonzero(change == 1)                                    #sorted by sub-cell, then time : starts and ends alternate in each sub-cell
    ec, et = np.nonzero(change == -1)
    complete = (st > 0)&(et < nsnap)
    return (et-st)[complete], sc[complete], int((~complete).sum())


def parse_states(arg,elements):#-e El:f,... : (element index, factor) ; f > 1 enriched (x >= f x bulk), f < 1 depleted (x <= f x bulk)
    states = []
    for part in arg.split(','):
        el, _, f = part.partition(':')
        try:
            factor = float(f)
        except ValueError:
            factor = -1.0
        if el.strip() not in elements or factor <= 0.0 or factor == 1.0:
            print('-e : element:factor separated by commas, e.g. Ni:1.5 (enriched, x >= 1.5 x the mole fraction of the cell) or Ni:0.5 (depleted) ;',
                  'the elements of the UMD file are',','.join(elements))
            sys.exit(2)
        states.append((elements.index(el.strip()),factor))
    return states


def fmt(value,error=None,digits=4):
    if error is None or not np.isfinite(error):
        return ('%.'+str(digits)+'g') % value
    return ('%.'+str(digits)+'g') % value+' +/- '+'%.2g' % error


def main(argv):
    umdpf.headerumd()
    start = time.time()
    umdfile = ''
    Nsteps = 1
    InitialStep = 0
    nCores = None
    grid = None
    length = 5.0
    volatiles = None
    threshold = 0.2
    densitybin = 0.05
    view = False
    cellarg = None
    atomarg = None
    statearg = None
    transient = 50.0
    usage = 'composition_umd.py -f <umdfile> -n <grid> -l <sub-cell length> -g <volatile elements> -d <void threshold> -i <InitialStep> -s <Sampling> -b <density bin> -k <nCores> -c <sub-cells> -a <elements> -e <El:factor> -t <transient> -v'
    try:
        opts, arg = getopt.getopt(argv,"hf:n:l:g:d:i:s:b:k:c:a:e:t:v",["fumdfile","ngrid","llength","gvolatiles","dthreshold","iInitialStep","sSampling","bbin","knCores","csubcells","aatoms","estates","ttransient","vview"])
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('Program to analyze the compositional fluctuations of a run : the atoms of each element are counted in a grid of sub-cells at every')
            print('snapshot. Usage :')
            print(usage)
            print('-f : UMD file. -f all : every *.umd.dat file of the current folder, one after the other.')
            print('-n : number of sub-cells along a, b, c : one number for the three, or na,nb,nc. Default : from -l.')
            print('-l : length of a sub-cell, in A : the grid has round(width/l) sub-cells along each lattice vector, at least 1 (width = distance')
            print('     between the opposite faces of the cell). Default 5 A, about 10 atoms of a silicate melt per sub-cell.')
            print('-g : volatile (gas) elements, separated by commas, or none. Default : those of '+','.join(VOLATILES)+' present in the run.')
            print('     The voids are defined by the density of the other elements (the melt) : the volatiles fill them.')
            print('-d : a sub-cell is a void when the density of its melt elements is below d times the density of the melt, the density seen by')
            print('     the mass of the melt elements (mean over the sub-cells weighted by their mass). Default 0.2.')
            print('-i : InitialStep. number of initial snapshots discarded. Default is 0.')
            print('-s : Sampling. one snapshot out of s is used. Default is 1.')
            print('-b : width of the bins of the density distribution, in g/cm3. Default 0.05.')
            print('-k : number of cores. Default is all.')
            print('-c : sub-cells followed one by one : ia:ib:ic (indexes from 0 along a, b, c, e.g. 0:3:2), several separated by commas (1:1:1,5:2:2),')
            print('     r (one at random), r5 (five at random), all. The table subcells.dat has one line per sub-cell : its mean')
            print('     number of atoms of each element of -a with errors, mole fractions and their ratio to the whole cell (enrichment), fraction of the')
            print('     time it is a void, density, and how long its content is kept ; with -c all, the map of the composition. Each sub-cell also gets')
            print('     the file subcell_<ia>_<ib>_<ic>.dat : per snapshot the number of atoms of each element of -a, their mole')
            print('     fractions, the density (all the atoms, the melt elements) and whether it is a void, with its means in the header. The files of other')
            print('     sub-cells, from former runs, are kept. More than 125 sub-cells (5x5x5) : a warning, as many files are written.')
            print('-a : elements followed in the files of -c and in the plot, separated by commas, e.g. Si,O. Default : all. The other outputs always')
            print('     take all the atoms (the voids depend on all the melt elements).')
            print('-e : composition states of a sub-cell, followed in time with the void and dense states : El:f, several separated by commas ; f > 1 :')
            print('     enriched, x_El >= f times the mole fraction of El in the cell (Ni:1.5) ; f < 1 : depleted, x_El <= f times it (Ni:0.5).')
            print('-t : regions (bubbles, droplets) living less than t fs are fluctuations of their surface, counted apart. Default 50 fs.')
            print('     The exchange of atoms between the gas and the condensed phase, its fluxes and equilibrium : exchange_umd.py.')
            print('-v : plot of the sub-cells of -c (without -c : every sub-cell, in random order) : the number of atoms of each element of -a against time')
            print('     (stacked : the top is their total), and the density of the sub-cell with the void threshold, the voids shaded ; Previous and')
            print('     Next buttons.')
            print('Any cell, orthogonal or not, and cells changing at every snapshot (NPT) : the sub-cells are cut in reduced coordinates, so they all')
            print('have the same volume, V/n, at every snapshot.')
            print('Outputs : all in the folder <base>.comp.<grid>, <grid> = na x nb x nc, e.g. pyrolite.comp.6x6x6 (one folder per grid ; delete it to remove')
            print('the analysis) :')
            print('  summary.dat : for each element the variance of its concentration in a sub-cell, N_El - x_El N_all, divided by that of a random')
            print('     alloy given the atoms of each sub-cell (hypergeometric) : 1 random, > 1 segregation (demixing), < 1 chemical order ; density')
            print('     fluctuations do not enter it (the measure for alloys) ; and the variance of its count divided by that of ideal mixing')
            print('     (each atom anywhere in the cell independently : a binomial law), 1 for an ideal gas, below 1 for a dense liquid, above 1 when')
            print('     the element gathers ; the correlation of the counts of every pair of elements (0 for ideal mixing, > 0 when they go together) ;')
            print('     voids and dense regions ; the times over which the local composition is kept (1/e of the correlations) ; errors of the means')
            print('     from their autocorrelation time, as averages.py.')
            print('  dist.dat : distribution of the number of atoms of each element (and of all) in a sub-cell, with the binomial law')
            print('  density.dat : distribution of the local density (all the atoms, the melt elements, the volatiles), g/cm3')
            print('  voids.dat : per snapshot, the voids (fraction of the volume, number of regions connected by faces, across the')
            print('     periodic faces too, and the largest), the dense regions (droplets in an expanded run), the density of each, and the fraction of')
            print('     the atoms of each element in the voids and in the largest dense region (partition between gas and condensed phase : with small')
            print('     sub-cells a heavy atom alone in the gas makes its sub-cell dense, so the fraction in the largest dense region is the steadier one)')
            print('  fluctuations.dat : per snapshot, the variance ratios of the counts and of the concentrations of each element and the correlations')
            print('     between the elements : segregation or mixing in progress (a quench) shows as a change in time')
            print('  episodes.dat : distribution of the durations of the periods a sub-cell stays void, dense, or in each state of -e (complete periods ;')
            print('     those cut by the start or the end of the run are counted in the summary but left out) ; mean durations in summary.dat')
            print('  time.dat : correlations in time of the counts of each element in a sub-cell, of the voids and of the density')
            print('     of the melt elements, from the deviations from the mean of the run : a heterogeneity that stays in place keeps them above 0.')
            sys.exit()
        elif opt in ("-f","--fumdfile"):
            umdfile = str(arg)
        elif opt in ("-n","--ngrid"):
            grid = [int(x) for x in arg.split(',')]
            if len(grid) == 1:
                grid = grid*3
            if len(grid) != 3 or min(grid) < 1:
                print('-n : one number, or three separated by commas, each at least 1')
                sys.exit(2)
        elif opt in ("-l","--llength"):
            length = float(arg)
        elif opt in ("-v","--vview"):
            view = True
        elif opt in ("-c","--csubcells"):
            cellarg = arg.strip().lower()
        elif opt in ("-t","--ttransient"):
            transient = float(arg)
        elif opt in ("-e","--estates"):
            statearg = arg.strip()
        elif opt in ("-a","--aatoms"):
            atomarg = [el.strip() for el in arg.split(',') if el.strip()]
        elif opt in ("-g","--gvolatiles"):
            volatiles = [] if arg.strip().lower() == 'none' else [el.strip() for el in arg.split(',') if el.strip()]
        elif opt in ("-d","--dthreshold"):
            threshold = float(arg)
        elif opt in ("-i","--iInitialStep"):
            InitialStep = int(arg)
        elif opt in ("-s","--sSampling"):
            Nsteps = int(arg)
        elif opt in ("-b","--bbin"):
            densitybin = float(arg)
        elif opt in ("-k","--knCores"):
            nCores = int(arg)
    umdfile = umdpf.umd_input(umdfile)                     #<run>.umd.dat.gz : decompressed into a temporary file, the outputs named from <run>.umd.dat
    if not os.path.isfile(umdpf.real_path(umdfile)):
        umdpf.missing_umd(umdfile)
        sys.exit()
    MyCrystal, TimeStep = umdpf.Crystallization(umdfile)
    elements = list(MyCrystal.elements)
    nel = len(elements)
    natom = MyCrystal.natom
    typat = np.array(MyCrystal.typat,dtype=np.int64)
    N = np.bincount(typat,minlength=nel)
    masses = np.array(MyCrystal.masses,dtype=float)
    if volatiles is None:
        volatiles = [el for el in elements if el in VOLATILES]
    unknown = [el for el in volatiles if el not in elements]
    if unknown:
        print('the volatile elements',','.join(unknown),'are not in the UMD file, whose elements are',','.join(elements))
        sys.exit(2)
    isvol = np.array([el in volatiles for el in elements])
    if isvol.all():
        print('all the elements are volatile : the voids are defined by the density of all the atoms')
        isvol[:] = False

    offsets, cells = umdpf.index_snapshots(umdfile,natom)
    offsets, cells = offsets[InitialStep::Nsteps], cells[InitialStep::Nsteps]
    nsnap = len(offsets)
    if nsnap < 8:
        print('Only',nsnap,'snapshots after the snapshot',InitialStep,'with sampling',Nsteps,': at least 8 are needed')
        sys.exit()
    if np.isnan(cells).any():
        print('ERROR : snapshots of',umdfile,'without rprimd or acell lines')
        sys.exit()
    if grid is None:
        grid = [max(1,int(round(w/length))) for w in widths(cells[0])]
    grid = np.array(grid,dtype=np.int64)
    ncell = int(np.prod(grid))
    if ncell < 2:
        print('the grid has a single sub-cell : give a smaller -l or a larger -n')
        sys.exit(2)
    states = parse_states(statearg,elements) if statearg else []
    rng = np.random.default_rng()
    chosen = parse_cells(cellarg,grid,rng) if cellarg else []
    viewcells = (chosen or [parse_cells('all',grid,rng)[m] for m in rng.permutation(int(np.prod(grid)))]) if view else []   #a plot only with -v
    if atomarg is None:
        tracked = list(range(nel))
    else:
        unknown = [el for el in atomarg if el not in elements]
        if unknown or not atomarg:
            print('-a : the elements',','.join(unknown),'are not in the UMD file, whose elements are',','.join(elements))
            sys.exit(2)
        tracked = [e for e, el in enumerate(elements) if el in atomarg]   #in the order of the UMD file
    tag = '%dx%dx%d' % tuple(grid)
    base = umdfile[:-8] if umdfile.endswith('.umd.dat') else os.path.splitext(umdfile)[0]
    outdir = base+'.comp.'+tag                                         #all the outputs of this run and grid : removed at once
    if os.path.exists(outdir) and not os.path.isdir(outdir):
        print(outdir,'is a file, not a folder : remove or rename it')
        sys.exit(2)
    os.makedirs(outdir,exist_ok=True)
    names = {kind:os.path.join(outdir,(kind or 'summary')+'.dat') for kind in ('','dist','density','voids','time','fluctuations','episodes','regions')}
    umdpf.remove_former_outputs(list(names.values()))
    volumes = np.abs(np.linalg.det(cells))
    subvol = volumes/ncell                                          #A^3, every sub-cell of a snapshot
    subwidth = widths(cells[0])/grid
    dt = (TimeStep if TimeStep > 0 else 1.0)*Nsteps
    if TimeStep <= 0:
        print('WARNING : no timestep in the UMD file : times are given in snapshots')
    print(nsnap,'snapshots ; grid',tag,'=',ncell,'sub-cells of','%.3f' % subvol.mean(),'A^3, widths','%.2f %.2f %.2f' % tuple(subwidth),'A ;',
          '%.2f' % (natom/ncell),'atoms per sub-cell on average')
    print('volatile elements :',','.join(volatiles) if isvol.any() else 'none','; melt elements :',','.join(el for el, v in zip(elements,isvol) if not v))

    #Counts : snapshot x sub-cell x element, read in batches by worker processes (as gofr_umd.py, bonding_umd.py)
    nworkers = max(1,min(nCores if nCores is not None else (os.cpu_count() or 1),nsnap))
    snapbytes = (os.path.getsize(umdpf.real_path(umdfile))-offsets[0])/nsnap
    tasks = umdpf.batch_ranges(nsnap,snapbytes,nworkers)
    worker = partial(count_batch,umdfile=umdfile,offsets=offsets,natom=natom,sampling=Nsteps,grid=grid,typat=typat,nel=nel)
    C = np.empty((nsnap,ncell,nel),dtype=np.int32)
    for (a,b), result in zip(tasks,umdpf.ordered_map(worker,tasks,nworkers)):
        C[a:b] = result
    assert np.all(C.sum(axis=1) == N[None,:])

    p = 1.0/ncell
    mu = N*p
    labels = elements+['all']
    Ntot = np.append(N,natom)

    #(1) Distributions of the counts, with the binomial law of ideal mixing, and of the local density
    Call = np.concatenate([C,C.sum(axis=2,keepdims=True)],axis=2)
    nmax = int(Call.max())
    hist = np.array([np.bincount(Call[:,:,e].ravel(),minlength=nmax+1)[:nmax+1] for e in range(nel+1)],dtype=float)/(nsnap*ncell)
    ks = np.arange(nmax+1)
    with open(names['dist'],'w') as f:
        f.write('# '+umdfile+' : '+str(nsnap)+' snapshots, grid '+tag+' ; P_X = probability that a sub-cell holds N atoms of X ; binom_X = the same\n')
        f.write('# for ideal mixing (every atom anywhere in the cell, independently) : binomial law of N_X atoms with probability 1/'+str(ncell)+'\n')
        f.write('N\t'+'\t'.join('P_'+l+'\tbinom_'+l for l in labels)+'\n')
        for k in ks:
            f.write(str(k)+'\t'+'\t'.join('%.6g\t%.6g' % (hist[e,k],binom.pmf(k,Ntot[e],p)) for e in range(nel+1))+'\n')
    rho_all = (C*masses).sum(axis=2)/subvol[:,None]*AMU_A3           #g/cm3, snapshot x sub-cell
    rho_melt = (C*(masses*~isvol)).sum(axis=2)/subvol[:,None]*AMU_A3
    rho_vol = rho_all-rho_melt
    top = max(rho_all.max(),densitybin)
    edges = np.arange(0.0,top+densitybin,densitybin)
    with open(names['density'],'w') as f:
        f.write('# '+umdfile+' : '+str(nsnap)+' snapshots, grid '+tag+' ; probability density of the density of a sub-cell (per g/cm3) : all the atoms,\n')
        f.write('# the melt elements, the volatile elements ('+(','.join(volatiles) if isvol.any() else 'none')+')\n')
        f.write('density_(g/cm3)\tP_all\tP_melt\tP_volatiles\n')
        hs = [np.histogram(r.ravel(),bins=edges,density=True)[0] for r in (rho_all,rho_melt,rho_vol)]
        for k in range(len(edges)-1):
            f.write('%.4f\t%.6g\t%.6g\t%.6g\n' % (0.5*(edges[k]+edges[k+1]),hs[0][k],hs[1][k],hs[2][k]))

    #(2) Fluctuations : variance of the counts against the binomial one, correlations between the elements ; per snapshot, then the statistics
    #of these series (errors from their autocorrelation time)
    D = Call-np.append(mu,natom*p)[None,None,:]
    cov = np.einsum('tce,tcf->tef',D,D,optimize=True)/ncell              #snapshot x element x element (all included)
    binvar = Ntot*p*(1.0-p)
    var = [statistics(cov[:,e,e]) for e in range(nel+1)]
    covstats = {(e,f):statistics(cov[:,e,f]) for e in range(nel+1) for f in range(nel+1)}
    meanvar = np.array([v[0] for v in var])
    corr = np.array([[covstats[(e,f)][0]/np.sqrt(meanvar[e]*meanvar[f]) if meanvar[e]*meanvar[f] > 0 else np.nan for f in range(nel+1)] for e in range(nel+1)])
    correrr = np.array([[covstats[(e,f)][2]/np.sqrt(meanvar[e]*meanvar[f]) if meanvar[e]*meanvar[f] > 0 and e != f else 0.0 for f in range(nel+1)] for e in range(nel+1)])

    #Concentration fluctuations : d = N_X - x_X N_all in a sub-cell ; for a random alloy (the atoms of each sub-cell given, their elements drawn
    #from the whole cell without replacement) <d^2> = N_all x (1-x) (N - N_all)/(N - 1) exactly (hypergeometric) : the ratio is 1, > 1 segregation
    total = C.sum(axis=2)
    xbulk = N/natom
    hyp = total*(natom-total)/(natom-1.0)                                #snapshot x sub-cell
    ccratio = np.full((nsnap,nel),np.nan)
    for e in range(nel):
        if 0.0 < xbulk[e] < 1.0:
            d = C[:,:,e]-xbulk[e]*total
            ccratio[:,e] = (d**2).sum(axis=1)/(xbulk[e]*(1.0-xbulk[e])*hyp.sum(axis=1))
    ccstats = [statistics(ccratio[:,e]) if np.isfinite(ccratio[:,e]).all() else (np.nan,np.nan,np.nan,False,0.0) for e in range(nel)]
    ccquarters = [quarters(ccratio[:,e]) if np.isfinite(ccratio[:,e]).all() else (np.full(4,np.nan),np.full(4,np.nan),0.0,True) for e in range(nel)]
    countratio = np.array([cov[:,e,e]/binvar[e] for e in range(nel+1)]).T  #snapshot x (elements, all)
    timesall = (InitialStep+np.arange(nsnap)*Nsteps)*(TimeStep if TimeStep > 0 else 1.0)
    pairs = [(e,f) for e in range(nel) for f in range(e+1,nel)]
    with np.errstate(invalid='ignore',divide='ignore'):
        paircorr = np.array([cov[:,e,f]/np.sqrt(cov[:,e,e]*cov[:,f,f]) for e, f in pairs]).T if pairs else np.zeros((nsnap,0))
    with open(names['fluctuations'],'w') as f:
        f.write('# '+umdfile+' : grid '+tag+' ; per snapshot : count_X = variance of the number of atoms of X in a sub-cell / ideal mixing (binomial) ;\n')
        f.write('# conc_X = variance of N_X - x_X N_all / random alloy (hypergeometric, 1 ; > 1 segregation, < 1 chemical order) ; corr_X-Y = correlation\n')
        f.write('# of the counts of X and Y over the sub-cells. A segregation or a mixing in progress shows as a change in time.\n')
        f.write('time_(fs)\t'+'\t'.join('count_'+l for l in labels)+'\t'+'\t'.join('conc_'+el for el in elements)+''.join('\tcorr_'+elements[e]+'-'+elements[g] for e, g in pairs)+'\n')
        f.write(''.join('%.2f\t' % timesall[t]+'\t'.join('%.5g' % v for v in countratio[t])+'\t'+'\t'.join('%.5g' % v for v in ccratio[t])
                        +''.join('\t%.5g' % v for v in paircorr[t])+'\n' for t in range(nsnap)))

    #(4) Voids : sub-cells whose density of melt elements is below threshold x the density of the melt (seen by its mass) ; connected regions
    rho_ref = (rho_melt**2).sum()/rho_melt.sum()
    void = rho_melt < threshold*rho_ref
    shape = tuple(grid)
    rows = []
    condensed = np.zeros((nsnap,ncell),dtype=bool)                       #the sub-cells of the largest dense region (the most melt atoms)
    voidlabels = np.zeros((nsnap,ncell),dtype=np.int32)                 #the regions of every snapshot, followed in time below
    denselabels = np.zeros((nsnap,ncell),dtype=np.int32)
    nvoid, ndense = np.zeros(nsnap,dtype=np.int64), np.zeros(nsnap,dtype=np.int64)
    for t in range(nsnap):
        m = void[t].reshape(shape)
        vlab, nvr = periodic_regions(m)
        dlab, ndr = periodic_regions(~m)
        vs = np.bincount(vlab.ravel(),minlength=nvr+1)[1:]
        ds = np.bincount(dlab.ravel(),minlength=ndr+1)[1:]
        meltatoms = C[t][:,~isvol].sum(axis=1).astype(float)              #per sub-cell, in the order of the grid (ia nb + ib) nc + ic
        dmelt = np.bincount(dlab.ravel(),weights=meltatoms,minlength=ndr+1)[1:]
        nv = int(m.sum())
        massvoid = (rho_all[t][void[t]]).sum()
        massdense = (rho_all[t][~void[t]]).sum()
        invoid = C[t][void[t]].sum(axis=0)/N
        voidlabels[t], denselabels[t], nvoid[t], ndense[t] = vlab.ravel(), dlab.ravel(), nvr, ndr
        if ndr:
            condensed[t] = dlab.ravel() == 1+int(np.argmax(dmelt))
        inlargest = C[t][condensed[t]].sum(axis=0)/N if ndr else np.zeros(nel)
        rows.append([nv/ncell,nvr,(vs.max() if nvr else 0)*subvol[t],ndr,(ds.max() if ndr else 0)*subvol[t],
                     dmelt.max()/meltatoms.sum() if ndr and meltatoms.sum() > 0 else np.nan,
                     massdense/(ncell-nv) if nv < ncell else np.nan,massvoid/nv if nv > 0 else np.nan]+list(invoid)+list(inlargest))
    rows = np.array(rows)
    voidcols = ['void_fraction','void_regions','largest_void_(A3)','dense_regions','largest_dense_(A3)','largest_dense_melt_fraction',
                'dense_density_(g/cm3)','void_density_(g/cm3)']+['in_voids_'+el for el in elements]+['in_largest_dense_'+el for el in elements]
    with open(names['voids'],'w') as f:
        f.write('# '+umdfile+' : grid '+tag+' ; a void is a sub-cell whose density of melt elements is below '+str(threshold)+' x '+'%.4f' % rho_ref+' g/cm3 ;\n')
        f.write('# regions connected by faces, across the periodic faces too ; densities of all the atoms in the dense sub-cells and in the voids ;\n')
        f.write('# largest_dense_melt_fraction : fraction of the atoms of the melt elements in the largest dense region (one droplet holding the melt :\n')
        f.write('# close to 1 ; the other dense regions are then molecules of the gas) ; in_voids_X : fraction of the atoms of X in the voids ;\n')
        f.write('# in_largest_dense_X : in the largest dense region (that with the most melt atoms) : the condensed phase, the rest being the gas\n')
        f.write('time_(fs)\t'+'\t'.join(voidcols)+'\n')
        for t in range(nsnap):
            f.write('%.2f\t' % ((InitialStep+t*Nsteps)*(TimeStep if TimeStep > 0 else 1.0))+'\t'.join('%.6g' % x for x in rows[t])+'\n')

    #(2b) States of the sub-cells and how long they last : void, dense, and the composition states of -e (enriched or depleted in an element) ;
    #for these, the fraction of the time expected for a random alloy (hypergeometric, given the atoms of each sub-cell)
    statedefs = [('void',void),('dense',~void)]
    expected = {}
    with np.errstate(invalid='ignore',divide='ignore'):
        for e, factor in states:
            name = ('enriched_' if factor > 1.0 else 'depleted_')+elements[e]+'_%g' % factor
            limit = factor*xbulk[e]*total                                #in atoms
            if factor > 1.0:
                statedefs.append((name,(total > 0)&(C[:,:,e] >= limit-1e-9)))
                kmin = np.ceil(limit-1e-9)
                prob = np.where(total > 0,hypergeom.sf(kmin-1,natom,N[e],total),0.0)
            else:
                statedefs.append((name,(total > 0)&(C[:,:,e] <= limit+1e-9)))
                prob = np.where(total > 0,hypergeom.cdf(np.floor(limit+1e-9),natom,N[e],total),0.0)
            expected[name] = float(prob.mean())
    epi = {name:episodes(st) for name, st in statedefs}
    maxlen = max([int(d.max()) for d, c, n in epi.values() if len(d)]+[0])
    with open(names['episodes'],'w') as f:
        f.write('# '+umdfile+' : grid '+tag+' ; number of complete periods (not cut by the start or the end of the run) during which a sub-cell stays\n')
        f.write('# in a state, against their duration : void (density of the melt elements below the threshold), dense, and the states of -e\n')
        f.write('duration_(fs)\t'+'\t'.join(name for name, st in statedefs)+'\n')
        hist = [np.bincount(epi[name][0],minlength=maxlen+1) for name, st in statedefs]
        f.write(''.join('%.2f\t' % (k*dt)+'\t'.join(str(int(h[k])) for h in hist)+'\n' for k in range(1,maxlen+1)))

    #(4b) Lifetimes of the regions : the void regions (bubbles ; the gas around a droplet) and the dense regions (droplets ; the melt around bubbles)
    #followed in time ; those living less than -t are fluctuations of the surface, counted apart from the means
    tracked_regions = {'void':track_regions(voidlabels,nvoid,C,subvol,dt),'dense':track_regions(denselabels,ndense,C,subvol,dt)}
    shortlife = max(1,int(np.ceil(transient/dt-1e-9))) if transient > 0 else 1
    with open(names['regions'],'w') as f:
        f.write('# '+umdfile+' : grid '+tag+' ; the connected void and dense regions followed in time by the sub-cells they share (the largest overlaps first) :\n')
        f.write('# birth : start (present at the first snapshot), new, split (from the region parent) ; death : end (present at the last snapshot), vanished,\n')
        f.write('# merged (into the region into) ; complete : neither cut by the start nor by the end ; volume and atoms per element : means over its life\n')
        f.write('region\tkind\tbirth_(fs)\tdeath_(fs)\tlifetime_(fs)\tcomplete\tbirth\tdeath\tparent\tinto\tmean_volume_(A3)\tmax_volume_(A3)\t'
                +'\t'.join('N_'+el for el in elements)+'\n')
        rows_out, k = [], 0
        offset = {}
        for kind in ('void','dense'):
            offset[kind] = k
            k += len(tracked_regions[kind])
        for kind in ('void','dense'):
            for i, g in enumerate(tracked_regions[kind]):
                life = g['end']-g['start']+1
                complete = int(g['birth'] != 'start' and g['death'] != 'end')
                rows_out.append('%d\t%s\t%.2f\t%.2f\t%.2f\t%d\t%s\t%s\t%d\t%d\t%.2f\t%.2f\t' % (offset[kind]+i,kind,timesall[g['start']],timesall[g['end']],life*dt,complete,
                                g['birth'],g['death'],offset[kind]+g['parent'] if g['parent'] >= 0 else -1,offset[kind]+g['into'] if g['into'] >= 0 else -1,
                                g['vol']/g['n'],g['volmax'])+'\t'.join('%.3f' % (v/g['n']) for v in g['atoms'])+'\n')
        f.write(''.join(rows_out))

    #(5) Correlations in time of the local composition, from the deviations from the mean of the run
    series = {'count_'+l:D[:,:,e].astype(float) for e, l in enumerate(labels)}
    series['void'] = void.astype(float)-void.mean()
    series['melt_density'] = rho_melt-rho_melt.mean()
    curves, decay = correlations(series,dt)
    nlag = len(next(iter(curves.values())))
    with open(names['time'],'w') as f:
        f.write('# '+umdfile+' : grid '+tag+' ; correlations in time, normalized to 1 at t = 0, of the deviations from the mean of the run, averaged over\n')
        f.write('# the sub-cells : count_X of the number of atoms of X in a sub-cell, void of the void indicator, melt_density of the density of the melt elements\n')
        f.write('time_(fs)\t'+'\t'.join(series)+'\n')
        for k in range(nlag):
            f.write('%.2f\t' % (k*dt)+'\t'.join('%.6g' % curves[name][k] for name in series)+'\n')

    #Summary
    out = []
    out.append('# '+umdfile+' : '+str(nsnap)+' snapshots from '+str(InitialStep)+' every '+str(Nsteps)+' ; '+'%.2f' % dt+' fs between them')
    out.append('# grid '+tag+' : '+str(ncell)+' sub-cells of '+'%.3f' % subvol.mean()+' A^3, widths '+'%.2f %.2f %.2f' % tuple(subwidth)+' A ; cell '+'%.3f' % volumes.mean()+' A^3')
    out.append('# volatile elements : '+(','.join(volatiles) if isvol.any() else 'none')+' ; void : density of the melt elements below '+str(threshold)
               +' x '+'%.4f' % rho_ref+' g/cm3 (the density of the melt seen by its mass)')
    out.append('# errors : of the means over the run, from their autocorrelation time ; > : a lower bound (fewer than 50 independent samples)')
    out.append('')
    out.append('Fluctuations of the number of atoms in a sub-cell : variance / variance of ideal mixing (binomial) : 1 ideal gas, < 1 dense liquid, > 1 gathering')
    out.append('element\tatoms\tvariance\terror\tbinomial_variance\tratio\terror_ratio\tfraction_in_voids\tfraction_in_largest_dense')
    fracvoid = [statistics(rows[:,8+e]) for e in range(nel)]
    fraclarge = [statistics(rows[:,8+nel+e]) for e in range(nel)]
    for e, l in enumerate(labels):
        v = var[e]
        ratio, rerr = v[0]/binvar[e], v[2]/binvar[e]
        out.append(l+'\t'+str(Ntot[e])+'\t'+'%.5g' % v[0]+'\t'+('>' if v[3] else '')+'%.2g' % v[2]+'\t'+'%.5g' % binvar[e]+'\t'
                   +'%.4f' % ratio+'\t'+('>' if v[3] else '')+'%.2g' % rerr+'\t'+('%.4f' % fracvoid[e][0] if e < nel else '')
                   +'\t'+('%.4f' % fraclarge[e][0] if e < nel else ''))
    out.append('')
    out.append('Fluctuations of the concentration in a sub-cell : variance of N_X - x_X N_all / that of a random alloy (hypergeometric) : 1 random,')
    out.append('> 1 segregation (demixing), < 1 chemical order ; density fluctuations do not enter it. The means of the four quarters of the run, each')
    out.append('with its own error, and change = the largest separation of two quarters in their combined errors (> 4 : the fluctuations change')
    out.append('during the run : segregation or mixing in progress ; when it happened : the quarters, and fluctuations.dat)')
    out.append('element\tx\tconcentration_ratio\terror\tquarter_1\tquarter_2\tquarter_3\tquarter_4\tquarter_errors\tchange')
    for e, el in enumerate(elements):
        c, (qm, qe, sep, low) = ccstats[e], ccquarters[e]
        out.append(el+'\t'+'%.5f' % xbulk[e]+'\t'+('%.4f' % c[0] if np.isfinite(c[0]) else 'nan')+'\t'+('>' if c[3] else '')+'%.2g' % c[2]+'\t'
                   +'\t'.join('%.4f' % m for m in qm)+'\t'+' '.join(('>' if low else '')+'%.2g' % x_ for x_ in qe)+'\t'+'%.1f' % sep)
    countq = [quarters(countratio[:,e]) for e in range(nel)]
    changed = [e for e in range(nel) if ccquarters[e][2] > 4.0 or countq[e][2] > 4.0]
    reliable = [e for e in changed if N[e] >= 10 and not (ccquarters[e][2] > 4.0 and ccquarters[e][3]) and not (countq[e][2] > 4.0 and countq[e][3])]
    if reliable:
        out.append('WARNING : the fluctuations of '+','.join(elements[e] for e in reliable)+' change during the run (two quarters differ by more than 4')
        out.append('errors) : segregation or mixing in progress, or not equilibrated ; see the quarters and fluctuations.dat')
    uncertain = [e for e in changed if e not in reliable]
    if uncertain:
        out.append('Note : the quarters of '+','.join(elements[e] for e in uncertain)+' also differ by more than 4 errors, but these errors are uncertain (fewer')
        out.append('than 10 atoms, or fewer than 50 independent samples in a quarter : > before the errors) : not a reliable sign of a change')
    out.append('')
    out.append('Correlation of the counts of two elements in a sub-cell (0 ideal mixing, > 0 together, < 0 apart) ; errors in the second table')
    out.append('\t'+'\t'.join(labels))
    for e, l in enumerate(labels):
        out.append(l+'\t'+'\t'.join('%.4f' % corr[e,f] for f in range(nel+1)))
    out.append('errors\t'+'\t'.join(labels))
    for e, l in enumerate(labels):
        out.append(l+'\t'+'\t'.join('%.2g' % correrr[e,f] for f in range(nel+1)))
    out.append('')
    out.append('Voids and dense regions (means over the run)')
    for k, name in enumerate(voidcols[:8]):
        st = statistics(rows[:,k][np.isfinite(rows[:,k])]) if np.isfinite(rows[:,k]).sum() >= 8 else None
        out.append(name+'\t'+('%.5g' % st[0]+'\t'+('>' if st[3] else '')+'%.2g' % st[2] if st else 'nan'))
    out.append('melt_density_seen_by_its_mass_(g/cm3)\t'+'%.5g' % rho_ref)
    out.append('')
    out.append('States of the sub-cells and how long they last : complete periods (cut ones : by the start or the end of the run, left out of the')
    out.append('means) ; mean duration with its error (std / sqrt(periods)) ; fraction of the time the sub-cells are in the state, and for the composition')
    out.append('states the fraction expected for a random alloy (hypergeometric, given the atoms of each sub-cell). Distributions in episodes.dat')
    out.append('state\tperiods\tcut_periods\tmean_(fs)\terror\tmedian_(fs)\tlongest_(fs)\ttime_fraction\trandom_alloy_fraction')
    for name, st in statedefs:
        d, cc, ncut = epi[name]
        out.append(name+'\t'+str(len(d))+'\t'+str(ncut)+'\t'+('%.4g' % (d.mean()*dt) if len(d) else 'nan')+'\t'
                   +('%.2g' % (d.std()/np.sqrt(len(d))*dt) if len(d) > 1 else 'nan')+'\t'+('%.4g' % (np.median(d)*dt) if len(d) else 'nan')+'\t'
                   +('%.4g' % (d.max()*dt) if len(d) else 'nan')+'\t'+'%.4f' % st.mean()+'\t'+('%.4f' % expected[name] if name in expected else ''))
    out.append('')
    out.append('Regions followed in time (regions.dat) : the connected void regions (bubbles, or the gas around a droplet) and dense regions (droplets, or')
    out.append('the melt around bubbles), followed by the sub-cells they share ; those living less than '+'%g' % transient+' fs (-t) are fluctuations, counted apart ;')
    out.append('lifetime : of the complete regions (neither cut by the start nor by the end of the run) living at least -t, mean with error (std/sqrt(n)) ;')
    out.append('births (new, by a split) and deaths (vanished, merged) of the regions living at least -t ; volume and atoms : means over their lives')
    out.append('kind\tregions\tshort_lived\tcomplete\tmean_lifetime_(fs)\terror\tmedian_(fs)\tlongest_(fs)\tnew\tsplit\tvanished\tmerged\tmean_volume_(A3)\t'
               +'\t'.join('N_'+el for el in elements))
    for kind in ('void','dense'):
        regs = tracked_regions[kind]
        lives = np.array([g['end']-g['start']+1 for g in regs])
        lasting = [g for g, l in zip(regs,lives) if l >= shortlife]
        comp = np.array([g['end']-g['start']+1 for g in lasting if g['birth'] != 'start' and g['death'] != 'end'])*dt
        count = lambda key, val: sum(1 for g in lasting if g[key] == val)
        vol = np.mean([g['vol']/g['n'] for g in lasting]) if lasting else np.nan
        atoms = np.mean([g['atoms']/g['n'] for g in lasting],axis=0) if lasting else np.full(nel,np.nan)
        out.append(kind+'\t%d\t%d\t%d\t' % (len(lasting),len(regs)-len(lasting),len(comp))
                   +('%.4g\t%s\t%.4g\t%.4g' % (comp.mean(),'%.2g' % (comp.std()/np.sqrt(len(comp))) if len(comp) > 1 else 'nan',np.median(comp),comp.max()) if len(comp) else 'nan\tnan\tnan\tnan')
                   +'\t%d\t%d\t%d\t%d\t' % (count('birth','new'),count('birth','split'),count('death','vanished'),count('death','merged'))
                   +'%.4g\t' % vol+'\t'.join('%.3f' % v for v in atoms))
    out.append('')
    out.append('Times over which the local composition is kept : decay = 1/e time of the correlation of the deviations from the mean of the run (nan : not')
    out.append('reached in half the run, a heterogeneity that stays) ; persistent = the part of the variance that stays in place over the run (0 : every')
    out.append('sub-cell alike on average ; noise of the means removed) ; fluctuating_decay = 1/e time of the fluctuations of each sub-cell about its')
    out.append('own mean, which decay even when the heterogeneity stays ; correlation_at_end = correlation from T/4 to T/2, close to persistent when the')
    out.append('fluctuations have decayed. Errors : of the times from 4 blocks of the run, of the persistent fraction from its spread over the sub-cells')
    out.append('series\tdecay_(fs)\terror\tpersistent\terror\tfluctuating_decay_(fs)\terror\tcorrelation_at_end')
    g = lambda v, f: (f % v) if np.isfinite(v) else 'nan'
    for name in series:
        d = decay[name]
        out.append(name+'\t'+g(d[0],'%.4g')+'\t'+g(d[1],'%.2g')+'\t'+g(d[2],'%.4f')+'\t'+g(d[3],'%.2g')+'\t'+g(d[4],'%.4g')+'\t'+g(d[5],'%.2g')+'\t'+g(d[6],'%.4f'))
    with open(names[''],'w') as f:
        f.write('\n'.join(out)+'\n')

    print()
    print('variance of the counts / ideal mixing :','  '.join(l+' '+'%.3f' % (var[e][0]/binvar[e]) for e, l in enumerate(labels)))
    vf = statistics(rows[:,0])
    print('voids :','%.3f' % vf[0],'of the volume ;','%.3g' % np.nanmean(rows[:,1]),'void regions and','%.3g' % np.nanmean(rows[:,3]),
          'dense regions per snapshot, the largest with','%.3f' % np.nanmean(rows[:,5]),'of the melt atoms ; density of the dense regions',
          '%.3f' % np.nanmean(rows[:,6]),'g/cm3')
    print('fraction of the atoms in the largest dense region (the condensed phase ; 1 - it : the gas) :',
          '  '.join(el+' '+'%.3f' % fraclarge[e][0] for e, el in enumerate(elements)))
    if isvol.any():
        print('fraction of the volatile atoms in the voids :','  '.join(el+' '+'%.3f' % fracvoid[e][0] for e, el in enumerate(elements) if isvol[e]))
    print('persistent fractions (the heterogeneity that stays over the run) :','  '.join(name+' '+'%.3f' % decay[name][2] for name in series))
    print('1/e times of the fluctuations of the sub-cells (fs) :','  '.join(name+' '+'%.4g' % decay[name][4] for name in series))
    print('Files',', '.join(os.path.basename(n) for n in names.values()),'written in the folder',outdir)
    if chosen:#per sub-cell of -c : one row of the table subcells.dat, and the file of its time series (the files of other
        #sub-cells, from former runs, are kept, as umd2out.py -k)
        if len(chosen) > 125:
            print('WARNING :',len(chosen),'sub-cells : as many files subcell_<ia>_<ib>_<ic>.dat are written in',outdir,'; the table subcells.dat',
                  'has one line per sub-cell')
        times = (InitialStep+np.arange(nsnap)*Nsteps)*(TimeStep if TimeStep > 0 else 1.0)
        total = C.sum(axis=2)
        bulk = N/natom                                                  #mole fractions of the whole cell
        tablename = os.path.join(outdir,'subcells.dat')
        columns = ['subcell']+[c for e in tracked for c in ('N_'+elements[e],'N_'+elements[e]+'_err','x_'+elements[e],'x_'+elements[e]+'/bulk')] \
                  +['N_all','N_all_err','void_time_fraction','void_time_fraction_err','density_all_(g/cm3)','density_all_err','density_melt_(g/cm3)',
                    'density_melt_err','decay_N_all_(fs)']+[c for name, st in statedefs for c in (name+'_time_fraction',name+'_mean_period_(fs)')]
        rows = []
        for ia, ib, ic in chosen:
            idx = (ia*grid[1]+ib)*grid[2]+ic
            stN = {e:statistics(C[:,idx,e]) for e in tracked}
            stall = statistics(total[:,idx])
            stvoid = statistics(void[:,idx].astype(float))
            strho, stmelt = statistics(rho_all[:,idx]), statistics(rho_melt[:,idx])
            dev = (total[:,idx]-total[:,idx].mean()).astype(float)[:,None]    #its own fluctuations : how long its content is kept
            c = acf(dev,max(2,nsnap//2))
            tau = decay_time(c/c[0],dt) if c[0] > 0 else np.nan
            x = {e:stN[e][0]/stall[0] if stall[0] > 0 else np.nan for e in tracked}
            values = [v for e in tracked for v in (stN[e][0],stN[e][2],x[e],x[e]/bulk[e])]+[stall[0],stall[2],stvoid[0],stvoid[2],strho[0],strho[2],stmelt[0],stmelt[2],tau]
            for name, st in statedefs:
                d, cc, ncut = epi[name]
                mine = d[cc == idx]
                values += [st[:,idx].mean(),mine.mean()*dt if len(mine) else np.nan]
            rows.append(('%d:%d:%d' % (ia,ib,ic),values))
            cellname = os.path.join(outdir,'subcell_%d_%d_%d.dat' % (ia,ib,ic))
            with np.errstate(invalid='ignore',divide='ignore'):
                frac = np.where(total[:,idx,None] > 0,C[:,idx,tracked]/total[:,idx,None],np.nan)
            with open(cellname,'w') as f:
                f.write('# '+umdfile+' : sub-cell %d:%d:%d of the grid %s (reduced coordinates %s to %s along a, b, c), %.3f A^3 on average\n'
                        % (ia,ib,ic,tag,'%.4f %.4f %.4f' % tuple(np.array([ia,ib,ic])/grid),'%.4f %.4f %.4f' % tuple((np.array([ia,ib,ic])+1)/grid),subvol.mean()))
                f.write('# means over the run (errors from the autocorrelation time) : '+' ; '.join(elements[e]+' %.3f +/- %.2g (x %.4f, %.3f x the cell)' % (stN[e][0],stN[e][2],x[e],x[e]/bulk[e]) for e in tracked)
                        +' ; all %.3f +/- %.2g atoms\n' % (stall[0],stall[2]))
                f.write('# void %.3f of the time ; density %.4f +/- %.2g g/cm3 (melt elements %.4f) ; its number of atoms decays to 1/e of its fluctuations in %s fs\n'
                        % (stvoid[0],strho[0],strho[2],stmelt[0],'%.4g' % tau if np.isfinite(tau) else 'more than half the run'))
                f.write('# N_X : atoms of X in the sub-cell ; x_X = N_X / all the atoms of the sub-cell (nan when empty) ; density of all the atoms and of the\n')
                f.write('# melt elements (g/cm3) ; void : 1 when the density of the melt elements is below '+'%.4f' % (threshold*rho_ref)+' g/cm3\n')
                f.write('time_(fs)\t'+'\t'.join('N_'+elements[e] for e in tracked)+'\tN_all\t'+'\t'.join('x_'+elements[e] for e in tracked)
                        +'\tdensity_all_(g/cm3)\tdensity_melt_(g/cm3)\tvoid\n')
                lines = []
                for t in range(nsnap):
                    lines.append('%.2f\t' % times[t]+'\t'.join(str(int(C[t,idx,e])) for e in tracked)+'\t'+str(int(total[t,idx]))+'\t'
                                 +'\t'.join('%.4f' % x_ for x_ in frac[t])+'\t%.5f\t%.5f\t%d\n' % (rho_all[t,idx],rho_melt[t,idx],void[t,idx]))
                f.write(''.join(lines))
        with open(tablename,'w') as f:#rewritten with the sub-cells of this run
            f.write('# '+umdfile+' : '+str(nsnap)+' snapshots, grid '+tag+' ; one line per sub-cell of -c (ia:ib:ic, from 0 along a, b, c) : means over the run\n')
            f.write('# of the number of atoms N_X with its error (from the autocorrelation time ; a lower bound when the run holds fewer than 50 independent\n')
            f.write('# samples), the mole fraction x_X = <N_X>/<N_all> and its ratio to that of the whole cell (> 1 : enriched), the fraction of the time\n')
            f.write('# the sub-cell is a void, its density, and the time over which its number of atoms decays to 1/e of its fluctuations (nan : more than\n')
            f.write('# half the run). Whole cell : '+' '.join(elements[e]+' x %.4f' % bulk[e] for e in tracked)+' ; mean '+'%.4f' % (natom/ncell)+' atoms per sub-cell\n')
            f.write('\t'.join(columns)+'\n')
            for name, values in rows:
                f.write(name+'\t'+'\t'.join('%.6g' % v for v in values)+'\n')
        print(len(chosen),'sub-cell(s) : subcells.dat and',', '.join('subcell_%d_%d_%d.dat' % c for c in chosen[:3])+(' ...' if len(chosen) > 3 else ''),'in',outdir)
    print('runtime : ','%.1f' % (time.time()-start),' s')
    if viewcells:
        import matplotlib.pyplot as plt
        times = (InitialStep+np.arange(nsnap)*Nsteps)*(TimeStep if TimeStep > 0 else 1.0)
        fig, draw = plot_subcells(viewcells,grid,times,C,elements,isvol,rho_all,rho_melt,threshold*rho_ref,os.path.basename(base),tracked)
        plt.show()


if __name__ == "__main__":
    if not umdpf.run_on_all(sys.argv[1:],main):#-f all : every UMD file of the folder
        main(sys.argv[1:])
