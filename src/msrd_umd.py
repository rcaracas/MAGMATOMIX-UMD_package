#!/usr/bin/env python3
###
#AUTHORS: RAZVAN CARACAS
###
#Mean square relative displacements (MSRD) of the bonds of coordination polyhedra (SiO4, AlO4, MgO6, ...) : for every bond central - adjacent
#atom of a polyhedron, the change of its bond vector b = r_adjacent - r_central between a time origin t0 and t0 + t, projected on the bond direction
#at t0 (longitudinal, parallel) and on two directions perpendicular to it (transverse) :
#   MSRD_par(t) = < [(b(t0+t) - b(t0)) . e_par(t0)]^2 >,  MSRD_perp1(t), MSRD_perp2(t) the same along e_perp1(t0), e_perp2(t0)
#averaged over the bonds and all the time origins. The bond vector is a relative position : the translation of the polyhedron cancels in it.
#Transverse directions : e_perp1 in the plane central - adjacent - second central atom when the adjacent atom is bonded to another central atom
#(a bridging O : the Si-O-Si plane, e_perp1 its bending, e_perp2 out of the plane), else in the plane of the bond and of the nearest other adjacent
#atom of the polyhedron (the O-Si-O plane) ; e_perp2 = e_par x e_perp1.
#Two frames : lab, the directions of the cell, in which the transverse MSRD also holds the libration and the rotational diffusion of the polyhedron ;
#body, the frame of the polyhedron, its rigid rotation between t0 and t0 + t removed (the rotation that best maps, mass-weighted, its shape at
#the start of its life onto its current shape, Kabsch) : the internal deformation only. Their difference is the rotation of the polyhedron.
#A polyhedron is followed by its central atom as long as it has the same adjacent atoms (a bond that breaks or an exchanged ligand ends its life ;
#a gap of at most -g fs during which it does not have the composition asked is bridged, left out).
#Static values, as measured by EXAFS : the distribution of the bond length r (mean, variance = the Debye-Waller factor of the distance, third
#cumulant), and the mean square relative displacements about the mean positions of the life of each bond, parallel and perpendicular to the mean
#bond (sigma2_par, sigma2_perp) ; the distance R0 between the mean positions differs from <r> by about sigma2_perp / (2 R0).

import sys, getopt, os.path, time
import numpy as np
import umd_processes_fast as umdpf
from rotations_umd import kabsch, block_error, clip_segments
from vibrations_umd import write_columns

GAP = 50.0              #fs : default of -g
MAXEL = 2.0E6           #numbers per array in the batches of the FFTs
ENDPAIRS = 0.1          #MSRD/2 of the summary : at the last lag that has at least this fraction of the pairs of times of lag 0
COLLINEAR = 0.05        #the reference vector of e_perp1 is not used when its component perpendicular to the bond is below this fraction of its length
FRAMES = ('lab','body')
CLASSES = ('bridging','nonbridging')
DIRS = ('par','perp1','perp2')
PAIRS = ((0,0,1.0),(1,1,1.0),(2,2,1.0),(0,1,2.0),(0,2,2.0),(1,2,2.0))   #components a, b of x_a x_b, and their weight in the double sum


COLORS = ['#2a78d6','#eb6834','#1baf7a','#eda100']     #par, perp1, perp2, perp : the categorical palette of the package (lifetimes.py), in its order
TEXT, MUTED, SURFACE, GRID = '#0b0b0b', '#52514e', '#fcfcfb', '#e4e3df'


def read_outputs(timefile,sumfile):#The MSRD in time (header -> column) and the summary (quantity -> (value, error)) written by main
    lines = [l for l in open(timefile) if not l.startswith('#')]
    head = lines[0].rstrip('\n').split('\t')
    data = np.array([[float(x) for x in l.split('\t')] for l in lines[1:]]).reshape(-1,len(head))
    columns = {h:data[:,k] for k,h in enumerate(head)}
    summary = {}
    for l in open(sumfile):
        entry = l.rstrip('\n').split('\t')
        if len(entry) == 3 and not l.startswith('#') and entry[0] != 'quantity':
            summary[entry[0]] = (float(entry[1]),float(entry[2]))
    return columns, summary


def plot_msrd(timefile,sumfile,central,adjacents,title,interactive=True):#One page per adjacent element : rows lab and body frames, columns all bonds, bridging,
    #non-bridging ; the MSRD along the bond and across it with their error bands, 2 sigma2 (static) dotted, the static values in text boxes.
    #interactive : buttons, Done for one page, Previous and Next for several (Next reads Done on the last one and closes the window ; none in the
    #PDF of -p). Returns the figure and the function that draws page k (for tests and for the PDF)
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    col, summ = read_outputs(timefile,sumfile)
    t = col['time_(fs)']
    fig, axes = plt.subplots(2,3,figsize=(13,7.6),facecolor=SURFACE)
    fig.subplots_adjust(left=0.07,right=0.985,top=0.75,bottom=0.11 if interactive else 0.08,wspace=0.18,hspace=0.30)
    names = {'par':'parallel (along the bond)','perp1':'perp1 (in the %s-%%s-%s or %%s-%s-%%s plane)' % (central,central,central),
             'perp2':'perp2 (normal to that plane)','perp':'perp = perp1 + perp2'}
    def value(key):
        return summ.get(key,(np.nan,np.nan))
    def draw(k):
        el = adjacents[k]
        pre = '%s-%s' % (central,el)
        for r,frame in enumerate(('lab','body')):
            row = axes[r]
            for ax in row:
                ax.clear()
            for c,tag in enumerate(('all','bridging','nonbridging')):
                ax = axes[r,c]
                ax.set_facecolor(SURFACE)
                for side in ('top','right'):
                    ax.spines[side].set_visible(False)
                for side in ('left','bottom'):
                    ax.spines[side].set_color(MUTED)
                ax.tick_params(colors=MUTED,labelsize=8,length=2)
                ax.grid(axis='y',color=GRID,linewidth=0.6)
                ax.set_axisbelow(True)
                comps = ('par','perp') if c == 0 else ('par','perp1','perp2')
                stem = '%s_%s_' % (frame,el) if c == 0 else '%s_%s_%s_' % (frame,el,tag)
                frac = value('%s_%s_fraction' % (pre,tag))[0]
                if not np.isfinite(frac) or frac == 0:
                    ax.text(0.5,0.5,'no %s bonds' % ('' if c == 0 else tag.replace('nonbridging','non-bridging')),transform=ax.transAxes,ha='center',va='center',color=MUTED,fontsize=9)
                    ax.set_xticks([])
                    ax.set_yticks([])
                else:
                    top = 0.0
                    for comp in comps:
                        color = COLORS[('par','perp1','perp2','perp').index(comp)]
                        y, e = col.get(stem+comp), col.get(stem+comp+'_err')
                        if y is None:
                            continue
                        ok = np.isfinite(y)
                        ax.plot(t[ok],y[ok],color=color,linewidth=2.0,solid_capstyle='round')
                        top = max(top,np.nanmax(y[ok]) if ok.any() else 0.0)
                        if e is not None:
                            band = ok & np.isfinite(e)
                            ax.fill_between(t[band],(y-e)[band],(y+e)[band],color=color,alpha=0.18,linewidth=0)
                    s2par = value('%s_%s_%s_sigma2_par_(A2)' % (pre,tag,frame))
                    s2perp = value('%s_%s_%s_sigma2_perp_(A2)' % (pre,tag,frame))
                    R0 = value('%s_%s_%s_R0_(A)' % (pre,tag,frame))[0]
                    ax.axhline(2*s2par[0],color=COLORS[0],linewidth=1.2,linestyle=(0,(1,2)))
                    if c == 0:
                        ax.axhline(2*s2perp[0],color=COLORS[3],linewidth=1.2,linestyle=(0,(1,2)))
                        statics = [2*s2par[0],2*s2perp[0]]
                    else:#perp1 + perp2 reach 2 sigma2_perp : each, on average, sigma2_perp
                        ax.axhline(s2perp[0],color=MUTED,linewidth=1.2,linestyle=(0,(1,2)))
                        statics = [2*s2par[0],s2perp[0]]
                    ax.text(0.03,0.97,'$\\sigma^2_\\parallel$ = %.4f $\\pm$ %.4f $\\AA^2$\n$\\sigma^2_\\perp$ = %.4f $\\pm$ %.4f $\\AA^2$\n$R_0$ = %.4f $\\AA$'
                            % (s2par[0],s2par[1],s2perp[0],s2perp[1],R0),transform=ax.transAxes,ha='left',va='top',fontsize=8,color=TEXT,
                            bbox=dict(boxstyle='round,pad=0.4',facecolor='white',edgecolor=GRID,linewidth=0.8))
                    ax.set_xlim(0,t[-1])
                    ax.set_ylim(0,1.6*max([top]+[v for v in statics if np.isfinite(v)]))     #the upper part free for the text box
                if r == 0:
                    rm, s2r = value('%s_%s_r_mean_(A)' % (pre,tag)), value('%s_%s_sigma2_r_(A2)' % (pre,tag))
                    head = {'all':'all bonds','bridging':'bridging bonds','nonbridging':'non-bridging bonds'}[tag]
                    sub = '' if not np.isfinite(frac) or frac == 0 else '\n%.1f %% of the bonds\n$\\langle r \\rangle$ = %.4f $\\AA$, $\\sigma^2_r$ = %.4f $\\pm$ %.4f $\\AA^2$' % (100*frac,rm[0],s2r[0],s2r[1])
                    ax.set_title(head+sub,color=TEXT,fontsize=9,loc='left')
                if c == 0:
                    ax.set_ylabel(('lab frame' if frame == 'lab' else 'body frame (rotation removed)')+'\nMSRD ($\\AA^2$)',color=TEXT,fontsize=9)
                if r == 1:
                    ax.set_xlabel('time (fs)',color=TEXT,fontsize=9)
        handles = ([Line2D([],[],color=COLORS[q],linewidth=2.0) for q in (0,1,2,3)]
                   +[Line2D([],[],color=COLORS[0],linewidth=1.2,linestyle=(0,(1,2))),Line2D([],[],color=COLORS[3],linewidth=1.2,linestyle=(0,(1,2))),
                     Line2D([],[],color=MUTED,linewidth=1.2,linestyle=(0,(1,2)))])
        labels = [names['par'],names['perp1'] % (el,el,el),names['perp2'],names['perp'],'2 $\\sigma^2_\\parallel$ static : the plateau of parallel once decorrelated',
                  '2 $\\sigma^2_\\perp$ static : the plateau of perp','$\\sigma^2_\\perp$ static : the mean plateau of perp1 and perp2']
        for old in fig.legends:
            old.remove()
        fig.legend(handles,labels,loc='upper left',bbox_to_anchor=(0.07,0.925),ncol=3,frameon=False,fontsize=8.5,labelcolor=TEXT)
        lives, life = value('lives')[0], value('mean_life_(fs)')[0]
        fig.suptitle('%s : %s-%s bonds of the polyhedra, mean square relative displacements%s\n%d lives of %.0f fs on average ; static $\\sigma^2$ about the mean positions of each life ; bands : errors from blocks'
                     % (title,central,el,'   (%d of %d)' % (k+1,len(adjacents)) if len(adjacents) > 1 else '',lives,life),color=TEXT,fontsize=10,x=0.07,ha='left',y=0.985)
        fig.canvas.draw_idle()
    draw(0)
    if interactive:
        umdpf.page_buttons(fig,len(adjacents),draw,[0.76,0.015,0.1,0.04],[0.87,0.015,0.1,0.04])
    return fig, draw


def lives_with_ligands(valid,lig,maxgap,minlen):#A life of a polyhedron : its central atom with the same adjacent atoms (lig (snapshot, unit, k),
    #sorted), the gaps of at most maxgap snapshots without the composition bridged ; another set of adjacent atoms ends it. Lives of fewer than minlen
    #snapshots are dropped. Returns the lives (unit, first, last) and the id of the life of each unit at each snapshot (-1 outside a life)
    lid = -np.ones(valid.shape,dtype=np.int32)
    lives = []
    for i in range(valid.shape[1]):
        tv = np.flatnonzero(valid[:,i])
        if len(tv) == 0:
            continue
        change = (lig[tv[1:],i] != lig[tv[:-1],i]).any(axis=1) | (np.diff(tv)-1 > maxgap)
        cut = np.flatnonzero(change)
        for a,b in zip(np.r_[0,cut+1],np.r_[cut,len(tv)-1]):
            s, e = int(tv[a]), int(tv[b])
            if e-s+1 >= minlen:
                lid[s:e+1,i] = len(lives)
                lives.append((i,s,e))
    return lives, lid


def frame_vectors(b,v):#Unit vectors e_par, e_perp1, e_perp2 (..., 3 directions, 3) of the bonds b, e_perp1 from the component of v perpendicular to b ;
    #zero where b is zero (masked)
    r = np.linalg.norm(b,axis=-1,keepdims=True)
    ep = np.divide(b,r,out=np.zeros_like(b),where=r > 0)
    w = v-(v*ep).sum(axis=-1,keepdims=True)*ep
    n = np.linalg.norm(w,axis=-1,keepdims=True)
    e1 = np.divide(w,n,out=np.zeros_like(w),where=n > 0)
    return np.stack([ep,e1,np.cross(ep,e1)],axis=-2)


class Sums:#Sums over the bonds and the time origins of a set of lives : numerators of the MSRD (frame, element, class, direction, lag), numbers of
    #pairs of times (element, class, lag) ; static sums per (frame, element, class) and of the bond lengths per (element, class)
    def __init__(self,nel,nl):
        self.nl = nl
        self.num = np.zeros((2,nel,2,3,nl))
        self.den = np.zeros((nel,2,nl))
        self.par = np.zeros((2,nel,2))                      #sum of the squared deviations from the mean position of the life, along the mean bond
        self.dev = np.zeros((2,nel,2))                      #sum of the squared deviations
        self.R0 = np.zeros((2,nel,2))                       #sum of the distances between the mean positions
        self.n = np.zeros((nel,2))                          #number of (bond, snapshot)
        self.r = np.zeros((nel,2,3))                        #sums of r, r^2, r^3

    def add(self,x,e,m,cls,el,nel):#x (frame, n, k, L, 3) bond vectors, e (frame, n, k, L, 3, 3) their frames at each time, m (n, L) mask,
        #cls (n, k, L) 1 for bridging, el (n, k) element of each bond
        n, k, L = cls.shape
        nl = min(self.nl,L)
        mk = np.broadcast_to(m[:,None,:],(n,k,L))
        Fm = np.fft.rfft(mk,n=2*L,axis=2)
        weights = [[mk*((el[:,:,None] == q) & (cls == (1 if c == 0 else 0))) for c in range(2)] for q in range(nel)]
        for q in range(nel):
            for c in range(2):
                W = np.fft.rfft(weights[q][c],n=2*L,axis=2)
                self.den[q,c,:nl] += np.fft.irfft((W.conj()*Fm).sum(axis=(0,1)),n=2*L)[:nl]
        for f in range(2):
            xf = x[f]
            Fy = np.fft.rfft(xf,n=2*L,axis=2)                                          #(n, k, F, 3)
            FA = [np.fft.rfft(xf[...,a]*xf[...,b],n=2*L,axis=2) for a,b,_ in PAIRS]
            for d in range(3):
                ed = e[f][...,d,:]                                                      #(n, k, L, 3)
                p = (xf*ed).sum(axis=-1)                                                #projection of the bond at t0 on the direction
                for q in range(nel):
                    for c in range(2):
                        w = weights[q][c]
                        if not w.any():
                            continue
                        spec = np.fft.rfft(w*p*p,n=2*L,axis=2).conj()*Fm
                        for a in range(3):
                            spec -= 2.0*np.fft.rfft(w*p*ed[...,a],n=2*L,axis=2).conj()*Fy[...,a]
                        for (a,b,cab),F in zip(PAIRS,FA):
                            spec += cab*np.fft.rfft(w*ed[...,a]*ed[...,b],n=2*L,axis=2).conj()*F
                        self.num[f,q,c,d,:nl] += np.fft.irfft(spec.sum(axis=(0,1)),n=2*L)[:nl]
        #Static : deviations from the mean position of the life of each bond (in each frame) ; bond lengths
        cnt = mk.sum(axis=2)                                                            #(n, k)
        r = np.linalg.norm(x[0],axis=-1)                                                #(n, k, L), zero where masked
        for f in range(2):
            mean = x[f].sum(axis=2)/np.maximum(cnt,1)[...,None]                         #(n, k, 3)
            R0 = np.linalg.norm(mean,axis=-1)
            u = np.divide(mean,R0[...,None],out=np.zeros_like(mean),where=R0[...,None] > 0)
            dv = (x[f]-mean[:,:,None,:])*mk[...,None]
            dpar = (dv*u[:,:,None,:]).sum(axis=-1)**2
            dall = (dv**2).sum(axis=-1)
            for q in range(nel):
                for c in range(2):
                    w = weights[q][c]
                    self.par[f,q,c] += (w*dpar).sum()
                    self.dev[f,q,c] += (w*dall).sum()
                    self.R0[f,q,c] += (w*R0[:,:,None]).sum()
                    if f == 0:
                        self.n[q,c] += w.sum()
                        self.r[q,c] += [(w*r).sum(),(w*r**2).sum(),(w*r**3).sum()]


def batches(segments,B,V,C,Rot,used,lid):#Lives (unit, first, last) grouped by their length, in pieces of at most MAXEL numbers : x, e, m, cls of
    #Sums.add (the lab frame and the body frame), the units of the batch
    groups = {}
    for i,s,e in segments:
        groups.setdefault(e-s+1,[]).append((i,s))
    k = B.shape[2]
    for L,group in groups.items():
        step = max(1,int(MAXEL//(60*k*L)))
        for a in range(0,len(group),step):
            part = group[a:a+step]
            m = np.stack([used[s:s+L,i] for i,s in part]).astype(float)                             #(n, L)
            b = np.stack([B[s:s+L,i] for i,s in part]).transpose(0,2,1,3)*m[:,None,:,None]          #(n, k, L, 3)
            v = np.stack([V[s:s+L,i] for i,s in part]).transpose(0,2,1,3)
            cls = np.stack([C[s:s+L,i] for i,s in part]).transpose(0,2,1)*m[:,None,:].astype(np.int8)
            R = np.stack([Rot[s:s+L,i] for i,s in part])                                            #(n, L, 3, 3) : current = R ref
            e = frame_vectors(b,v)*m[:,None,:,None,None]                                            #(n, k, L, 3, 3)
            xb = np.einsum('nkti,ntij->nktj',b,R)                                                   #body frame : R^T b
            eb = np.einsum('nktdi,ntij->nktdj',e,R)
            yield np.stack([b,xb]), np.stack([e,eb]), m, cls, part


def main(argv):
    umdpf.headerumd()
    start = time.time()
    UMDFile = BondFile = central = ''
    adjacents = []
    numbers = None
    maxgap = GAP
    minlife = 0.0
    Nsteps = 1
    nblocks = 4
    window = None
    view = pdf = False
    usage = 'msrd_umd.py -f <UMD file> -b <bonding file> -c <central element> -a <adjacent elements> -x <numbers of adjacent atoms> -g <gap> -m <min life> -s <sampling> -n <blocks> -w <window> -v -p'
    try:
        opts, arg = getopt.getopt(argv,"hf:b:c:a:x:g:m:s:n:w:vp")
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('msrd_umd.py : mean square relative displacements of the bonds of coordination polyhedra (SiO4, AlO4, MgO6...), along each bond')
            print('     (parallel) and across it (two transverse directions), against time and static (EXAFS Debye-Waller factors)')
            print(usage)
            print('-f : UMD file ; -f all : every *.umd.dat file of the folder, each with its bonding file')
            print('-b : bonding file of bonding_umd.py (default <base>.bonding.dat)')
            print('-c : central element of the polyhedra (one) ; -a : adjacent elements, separated by commas : -c Si -a O for SiO4')
            print('-x : number of adjacent atoms of each element of -a, e.g. -x 4 for SiO4, -x 6 for MgO6 (default : the most frequent composition')
            print('     of the central atoms that have neighbours). A snapshot of a central atom counts only when it has exactly this composition.')
            print('-g : longest gap, in fs, during which a polyhedron that lost its composition goes on, left out during the gap, when it comes back')
            print('     with the same adjacent atoms ; a longer gap, or other adjacent atoms, end its life. Default',GAP)
            print('-m : shortest life of a polyhedron, in fs ; shorter lives are left out. Default 0')
            print('-s : every s-th snapshot of the bonding file. Default 1')
            print('-n : number of blocks of the run for the errors. Default 4')
            print('-w : length of the MSRD in time, in snapshots. Default half a block, at most half the run')
            print('-v : plot : MSRD against time, lab and body frames, all, bridging and non-bridging bonds, with their errors, 2 sigma2 (static)')
            print('     dotted and the static values in text boxes ; one page per adjacent element (Previous, Next ; Done on the last page)')
            print('-p : the same plot written to <base>.<unit>.msrd.pdf, one page per adjacent element')
            print('For each bond central - adjacent atom, b = r_adjacent - r_central (the translation of the polyhedron cancels in it) ; MSRD(t) =')
            print('     < [(b(t0+t) - b(t0)) . e(t0)]^2 > over the bonds and all the time origins t0 of the lives, e(t0) along the bond (par) or across it :')
            print('     perp1 in the plane of the bond and of a second central atom bonded to the adjacent atom (bridging : the Si-O-Si plane), else of')
            print('     the nearest other adjacent atom of the polyhedron (non-bridging : the O-Si-O plane) ; perp2 normal to that plane.')
            print('     lab : directions of the cell (the transverse MSRD includes the rotation of the polyhedron) ; body : the rigid rotation of the')
            print('     polyhedron between t0 and t0+t removed (Kabsch, mass-weighted) : its internal deformation only.')
            print('Outputs, <base>.<unit> with the unit named from the central atom and its neighbours (SiO4, MgO6) :')
            print('  .msrd.time.dat : MSRD(t) in A^2, lab and body frames, per adjacent element : all bonds (par, perp = perp1 + perp2), bridging and')
            print('                   non-bridging bonds (par, perp1, perp2) by the state of the bond at t0, with errors from the blocks ; pairs of times')
            print('  .msrd.pdf      : with -p, the plot of -v')
            print('  .msrd.dat      : summary : polyhedra and lives, fraction of bridging bonds, bond length <r>, its variance (sigma2_r, the Debye-')
            print('                   Waller factor of the distance) and third cumulant, distance between the mean positions R0, static sigma2_par and')
            print('                   sigma2_perp about the mean positions of each life (lab and body), MSRD/2 at the end of the window, or at the')
            print('                   last lag with at least %d %% of the pairs of times of lag 0 when the lives are shorter' % (100*ENDPAIRS))
            sys.exit()
        elif opt == '-f':
            UMDFile = str(arg)
        elif opt == '-b':
            BondFile = str(arg)
        elif opt == '-c':
            central = str(arg)
        elif opt == '-a':
            adjacents = arg.split(',')
        elif opt == '-x':
            numbers = [int(x) for x in arg.split(',')]
        elif opt == '-g':
            maxgap = float(arg)
        elif opt == '-m':
            minlife = float(arg)
        elif opt == '-s':
            Nsteps = int(arg)
        elif opt == '-n':
            nblocks = int(arg)
        elif opt == '-w':
            window = int(arg)
        elif opt == '-v':
            view = True
        elif opt == '-p':
            pdf = True
    UMDFile = umdpf.umd_input(UMDFile)                     #<run>.umd.dat.gz : decompressed into a temporary file, the outputs named from <run>.umd.dat
    if not os.path.isfile(umdpf.real_path(UMDFile)):
        umdpf.missing_umd(UMDFile)
        sys.exit()
    base = UMDFile[:-8] if UMDFile.endswith('.umd.dat') else os.path.splitext(UMDFile)[0]
    if not BondFile:
        BondFile = base+'.bonding.dat'
    if not os.path.isfile(BondFile):
        print('the bonding file',BondFile,'does not exist : compute it with bonding_umd.py, or give it with -b')
        sys.exit()
    if not central or ',' in central or not adjacents:
        print('give the central element (-c, one element) and the adjacent elements (-a) of the polyhedra, e.g. -c Si -a O')
        sys.exit()
    if nblocks < 1 or Nsteps < 1 or maxgap < 0 or (window is not None and window < 2):
        print('-n, -s at least 1 ; -g at least 0 ; -w at least 2')
        sys.exit()
    if numbers is not None and len(numbers) != len(adjacents):
        print('-x takes one number per element of -a :',len(adjacents))
        sys.exit()

    #The bonds central - adjacent of every snapshot used, and the composition of every central atom
    CentIndexes,AdjIndexes,MyCrystal,[Bonds,BondsIndexes],TimeStep = umdpf.read_bonds(BondFile,[central],adjacents,Nsteps)
    if CentIndexes == -1:#an element missing from the run : read_bonds printed the error
        sys.exit()
    natom = MyCrystal.natom
    typat = np.asarray(MyCrystal.typat)
    iscentral = typat == MyCrystal.elements.index(central)
    adjtype = -np.ones(natom,dtype=np.int64)
    for k,el in enumerate(adjacents):
        adjtype[typat == MyCrystal.elements.index(el)] = k
    centres = np.flatnonzero(iscentral)
    N = len(centres)
    cid = -np.ones(natom,dtype=np.int64)
    cid[centres] = np.arange(N)
    MyCrystalUMD,TimeStepUMD = umdpf.Crystallization(UMDFile)
    ratio = max(1,int(round(TimeStep/TimeStepUMD))) if TimeStepUMD > 0 else Nsteps
    MyCrystalUMD,X,TimeStepUMD,length = umdpf.read_values(UMDFile,"xcart","lists",ratio,MyCrystal.firstsnapshot)
    nsnap = min(len(Bonds),len(X))
    if nsnap < len(Bonds):
        print('WARNING : the UMD file has fewer snapshots than the bonding file : only the first',nsnap,'are used')
    if nsnap < 2*nblocks:
        print('Only',nsnap,'snapshots : too few for',nblocks,'blocks')
        sys.exit()
    cells = umdpf.read_cells(UMDFile,nsnap,ratio,MyCrystal.firstsnapshot)
    inverses = np.linalg.inv(cells)
    masses = np.asarray(MyCrystalUMD.masses,dtype=float)

    def bonds_of(t):#The bonds central atom -> adjacent atom of the snapshot t
        offsets = np.asarray(BondsIndexes[t][:natom+1])
        at = np.repeat(np.arange(natom),np.diff(offsets))
        ne = np.asarray(Bonds[t][offsets[0]:offsets[-1]],dtype=np.int64)
        keep = iscentral[at] & (adjtype[ne] >= 0) & (at != ne)
        return at[keep], ne[keep]

    comp = np.zeros((nsnap,N,len(adjacents)),dtype=np.int16)
    for t in range(nsnap):
        at, ne = bonds_of(t)
        np.add.at(comp[t],(cid[at],adjtype[ne]),1)
    kinds, counts = np.unique(comp.reshape(-1,len(adjacents)),axis=0,return_counts=True)
    order = np.argsort(-counts)
    def label(c):
        return central+''.join(el+(str(n) if n != 1 else '') for el,n in zip(adjacents,c) if n > 0)
    print('Compositions of the',N,central,'atoms over',nsnap,'snapshots :',
          ', '.join('%s %.2f %%' % (label(kinds[k]),100.0*counts[k]/counts.sum()) for k in order[:6]))
    if numbers is None:
        bound = [k for k in order if kinds[k].sum() > 0]
        if not bound:
            print('No',central,'atom is bonded to',', '.join(adjacents),'in',BondFile)
            sys.exit()
        numbers = [int(x) for x in kinds[bound[0]]]
    species = label(numbers)
    k = sum(numbers)
    if k == 0:
        print('-x : a polyhedron has at least one adjacent atom')
        sys.exit()
    valid = (comp == np.array(numbers)).all(axis=2)
    del comp
    lig = -np.ones((nsnap,N,k),dtype=np.int64)                         #adjacent atoms of each valid polyhedron, sorted (so grouped by element)
    for t in range(nsnap):
        at, ne = bonds_of(t)
        keep = valid[t,cid[at]]
        at, ne = at[keep], ne[keep]
        srt = np.lexsort((ne,cid[at]))
        lig[t,cid[at[srt][::k]]] = ne[srt].reshape(-1,k)
    lives, lid = lives_with_ligands(valid,lig,int(maxgap//TimeStep),max(1,int(np.ceil(minlife/TimeStep))))
    used = valid & (lid >= 0)
    if not used.any():
        print('No polyhedron',species,'in the run (with lives of at least',minlife,'fs)')
        sys.exit()
    spans = np.array([e-s+1 for i,s,e in lives])
    gaps = int(spans.sum()-used.sum())
    print('Polyhedra',species,': %.2f per snapshot (%d at most), %d lives of %.1f fs on average (%d over the whole run) ; gaps up to %g fs bridged : %d snapshots of polyhedra left out'
          % (used.sum()/nsnap,used.sum(axis=1).max(),len(lives),spans.mean()*TimeStep,int((spans == nsnap).sum()),maxgap,gaps))
    sitetype = adjtype[lig[used][0]]                                   #element of each site : the same for every polyhedron (atoms grouped by element)
    msite = np.r_[masses[MyCrystalUMD.elements.index(central)],masses[[MyCrystalUMD.elements.index(adjacents[q]) for q in sitetype]]]
    nel = len(adjacents)

    names = [base+'.'+species+ext for ext in ('.msrd.time.dat','.msrd.dat','.msrd.pdf')]
    umdpf.remove_former_outputs(names if pdf else names[:2])

    #Bond vectors (nearest images, cell of the snapshot), reference vectors of e_perp1, bridging or not, rotation of each polyhedron
    B = np.zeros((nsnap,N,k,3))
    V = np.zeros((nsnap,N,k,3))
    C = np.zeros((nsnap,N,k),dtype=np.int8)
    Rot = np.tile(np.eye(3),(nsnap,N,1,1))
    ref = np.zeros((N,k+1,3))
    life = -np.ones(N,dtype=np.int64)
    fallback = 0
    triclusters = 0
    for t in range(nsnap):
        units = np.flatnonzero(used[t])
        if len(units) == 0:
            continue
        at, ne = bonds_of(t)
        nc = np.bincount(ne,minlength=natom)                           #central atoms bonded to each adjacent atom
        L = lig[t,units]                                               #(n, k)
        cen = centres[units]
        d = (X[t][L]-X[t][cen][:,None,:]) @ inverses[t]
        d -= np.rint(d)
        b = d @ cells[t]                                               #(n, k, 3)
        B[t,units] = b
        #Non-bridging : the nearest other adjacent atom of the polyhedron
        if k > 1:
            dist = np.linalg.norm(b[:,:,None,:]-b[:,None,:,:],axis=3)
            dist[:,np.arange(k),np.arange(k)] = np.inf
            other = np.take_along_axis(b,np.argmin(dist,axis=2)[:,:,None],axis=1)
        else:
            other = np.broadcast_to(np.where(np.abs(b[...,:1]) < 0.9*np.linalg.norm(b,axis=2,keepdims=True),[1.0,0.0,0.0],[0.0,1.0,0.0]),b.shape).copy()
        v = other
        bridging = nc[L] >= 2
        C[t,units] = bridging
        if bridging.any():#the nearest other central atom bonded to the adjacent atom : the vector adjacent -> second central atom
            key = at*natom+ne
            srt = np.argsort(ne,kind='stable')
            ne_s, at_s = ne[srt], at[srt]
            ub, ok = np.nonzero(bridging)
            for u_,j_ in zip(ub,ok):
                o = L[u_,j_]
                a0, a1 = np.searchsorted(ne_s,o), np.searchsorted(ne_s,o,side='right')
                cands = at_s[a0:a1]
                cands = cands[cands != cen[u_]]
                if len(cands) > 1:
                    triclusters += 1
                dc = (X[t][cands]-X[t][o]) @ inverses[t]
                dc -= np.rint(dc)
                dc = dc @ cells[t]
                w = dc[np.argmin(np.linalg.norm(dc,axis=1))]
                bh = b[u_,j_]/np.linalg.norm(b[u_,j_])
                if np.linalg.norm(w-(w @ bh)*bh) > COLLINEAR*np.linalg.norm(w):
                    v[u_,j_] = w
                else:
                    fallback += 1
        V[t,units] = v
        #Body frame : the rotation that maps the shape of the polyhedron at the start of its life onto its current shape (mass-weighted)
        P = np.concatenate([np.zeros((len(units),1,3)),b],axis=1)
        rho = P-(msite[None,:,None]*P).sum(axis=1,keepdims=True)/msite.sum()
        new = life[units] != lid[t,units]
        ref[units[new]] = rho[new]
        life[units] = lid[t,units]
        Rot[t,units] = kabsch(ref[units],rho,msite)[0]
    del X
    if triclusters:
        print(triclusters,'bonds (bond x snapshot) to adjacent atoms bonded to three or more central atoms : e_perp1 toward the nearest other central atom')
    if fallback:
        print(fallback,'bonds (bond x snapshot) with a bridging angle within',COLLINEAR,'rad of 180 deg : e_perp1 from the nearest other adjacent atom instead')
    if k == 1:
        print('One adjacent atom per polyhedron : the non-bridging transverse directions are taken from the axes of the cell (perp1, perp2 have no')
        print('     meaning of their own there ; perp = perp1 + perp2 has)')

    #*** MSRD(t) and static values, over the run and over the blocks
    blocklen = nsnap//nblocks
    nlags = min(nsnap//2,window if window is not None else max(2,blocklen//2))
    nlb = min(nlags,blocklen//2)
    tlag = np.arange(nlags)*TimeStep
    def accumulate(t0,t1,nl):
        S = Sums(nel,nl)
        for x,e,m,cls,part in batches(clip_segments(lives,t0,t1),B,V,C,Rot,used,lid):
            el = np.broadcast_to(sitetype[None,:],(len(part),k))
            S.add(x,e,m,cls,el,nel)
        return S
    def msrd(S):#MSRD (frame, element, class 0 1 2 = all, bridging, non-bridging, direction, lag) and pairs (element, class, lag)
        num = np.concatenate([S.num.sum(axis=2,keepdims=True),S.num],axis=2)
        den = np.concatenate([S.den.sum(axis=1,keepdims=True),S.den],axis=1)
        with np.errstate(invalid='ignore',divide='ignore'):
            return np.where(den[None,:,:,None,:] > 0,num/np.where(den > 0,den,1.0)[None,:,:,None,:],np.nan), den
    def static(S):#per (element, class) : n, <r>, C2, C3 ; per (frame, element, class) : R0, sigma2_par, sigma2_perp
        n = np.concatenate([S.n.sum(axis=1,keepdims=True),S.n],axis=1)
        r = np.concatenate([S.r.sum(axis=1,keepdims=True),S.r],axis=1)
        par = np.concatenate([S.par.sum(axis=2,keepdims=True),S.par],axis=2)
        dev = np.concatenate([S.dev.sum(axis=2,keepdims=True),S.dev],axis=2)
        R0 = np.concatenate([S.R0.sum(axis=2,keepdims=True),S.R0],axis=2)
        with np.errstate(invalid='ignore',divide='ignore'):
            m1, m2, m3 = r[...,0]/n, r[...,1]/n, r[...,2]/n
            out = {'r':m1,'C2':m2-m1**2,'C3':m3-3*m1*m2+2*m1**3,'R0':R0/n,'s2par':par/n,'s2perp':(dev-par)/n}
            out['fraction'] = n/n[:,:1]
        return out
    S = accumulate(0,nsnap,nlags)
    M, pairs = msrd(S)
    st = static(S)
    if nlb > 1:
        blockS = [accumulate(b*blocklen,(b+1)*blocklen,nlb) for b in range(nblocks)]
        Merr = np.full(M.shape,np.nan)
        Merr[...,:nlb] = block_error([msrd(s)[0] for s in blockS])
        bst = [static(s) for s in blockS]
        sterr = {key:block_error([b[key] for b in bst]) for key in st}
    else:
        Merr = np.full(M.shape,np.nan)
        sterr = {key:np.full(np.shape(st[key]),np.nan) for key in st}
    print('Windows : %d snapshots (%.1f fs) for the run, %d for the errors from %d blocks' % (nlags,nlags*TimeStep,nlb,nblocks))

    #*** MSRD in time
    tags = ('all',)+CLASSES
    labels = ['time_(fs)']
    columns = [tlag]
    for f,frame in enumerate(FRAMES):
        for q,el in enumerate(adjacents):
            for c,tag in enumerate(tags):
                if c == 0:
                    comps = (('par',M[f,q,0,0],Merr[f,q,0,0]),('perp',M[f,q,0,1]+M[f,q,0,2],None))
                else:
                    comps = tuple((DIRS[d],M[f,q,c,d],Merr[f,q,c,d]) for d in range(3))
                for name,val,err in comps:
                    if err is None:#perp of all the bonds : its error from the blocks of perp1 + perp2
                        err = np.full(nlags,np.nan)
                        if nlb > 1:
                            err[:nlb] = block_error([(lambda m: m[f,q,0,1]+m[f,q,0,2])(msrd(s)[0]) for s in blockS])
                    head = '%s_%s_%s' % (frame,el,name) if c == 0 else '%s_%s_%s_%s' % (frame,el,tag,name)
                    labels += [head,head+'_err']
                    columns += [val,err]
    for q,el in enumerate(adjacents):
        for c,tag in enumerate(tags):
            labels.append('pairs_%s_%s' % (el,tag))
            columns.append(pairs[q,c])
    write_columns(names[0],['msrd_umd.py, '+UMDFile+', bonding file '+BondFile+', polyhedra '+species+' ; %d lives ; gaps up to %g fs bridged' % (len(lives),maxgap),
                            'MSRD(t) = < [(b(t0+t) - b(t0)) . e(t0)]^2 > in A^2, b = r_adjacent - r_central, over the bonds and all the time origins of the lives',
                            'par : e along the bond at t0 ; perp1 : across it, in the plane of the bond and of the second central atom of a bridging adjacent atom (Si-O-Si), else of the nearest other adjacent atom (O-Si-O) ; perp2 : normal to that plane ; perp = perp1 + perp2',
                            'lab : directions of the cell (rotation of the polyhedron included) ; body : the rigid rotation of the polyhedron between t0 and t0+t removed (Kabsch, mass-weighted)',
                            'bridging / nonbridging : state of the bond at t0 ; pairs : number of (origin, end) pairs of times over the run ; errors from %d blocks, up to their length' % nblocks],
                  labels,columns)

    #*** Summary
    rows = [('polyhedra_per_snapshot',used.sum()/nsnap,np.nan),('lives',len(lives),np.nan),('mean_life_(fs)',spans.mean()*TimeStep,np.nan),
            ('lives_over_the_whole_run',int((spans == nsnap).sum()),np.nan)]
    print('\nBonds of the polyhedra',species,'of',UMDFile,'(A, A^2 ; errors from %d blocks)' % nblocks)
    for q,el in enumerate(adjacents):
        for c,tag in enumerate(tags):
            pre = '%s-%s_%s' % (central,el,tag)
            if not np.isfinite(st['r'][q,c]):
                print('   %-22s : none' % pre)
                continue
            print('   %-22s : %5.1f %% of the bonds ; <r> %.5f +/- %.5f ; sigma2_r %.6f +/- %.6f ; C3 %.3e'
                  % (pre,100*st['fraction'][q,c],st['r'][q,c],sterr['r'][q,c],st['C2'][q,c],sterr['C2'][q,c],st['C3'][q,c]))
            rows += [(pre+'_fraction',st['fraction'][q,c],sterr['fraction'][q,c]),(pre+'_r_mean_(A)',st['r'][q,c],sterr['r'][q,c]),
                     (pre+'_sigma2_r_(A2)',st['C2'][q,c],sterr['C2'][q,c]),(pre+'_C3_(A3)',st['C3'][q,c],sterr['C3'][q,c])]
            kept = np.flatnonzero(pairs[q,c] >= ENDPAIRS*pairs[q,c,0])         #the last lag with at least ENDPAIRS of the pairs of lag 0
            le = int(kept[-1]) if len(kept) else 0
            rows.append((pre+'_end_lag_(fs)',tlag[le],np.nan))
            for f,frame in enumerate(FRAMES):
                R0, sp, sq = st['R0'][f,q,c], st['s2par'][f,q,c], st['s2perp'][f,q,c]
                end = M[f,q,c,:,le]/2.0
                print('      %-4s : R0 %.5f ; sigma2_par %.6f +/- %.6f, sigma2_perp %.6f +/- %.6f ; <r> - R0 %.6f, sigma2_perp/(2 R0) %.6f ; MSRD/2 at %.0f fs : par %.6f, perp %.6f'
                      % (frame,R0,sp,sterr['s2par'][f,q,c],sq,sterr['s2perp'][f,q,c],st['r'][q,c]-R0,sq/(2*R0),tlag[le],end[0],end[1]+end[2]))
                rows += [(pre+'_'+frame+'_R0_(A)',R0,sterr['R0'][f,q,c]),(pre+'_'+frame+'_sigma2_par_(A2)',sp,sterr['s2par'][f,q,c]),
                         (pre+'_'+frame+'_sigma2_perp_(A2)',sq,sterr['s2perp'][f,q,c]),(pre+'_'+frame+'_r_minus_R0_(A)',st['r'][q,c]-R0,np.nan),
                         (pre+'_'+frame+'_MSRD_par_end_half_(A2)',end[0],Merr[f,q,c,0,le]/2.0)]
                if c == 0:
                    rows.append((pre+'_'+frame+'_MSRD_perp_end_half_(A2)',end[1]+end[2],np.nan))
                else:
                    rows += [(pre+'_'+frame+'_MSRD_perp1_end_half_(A2)',end[1],Merr[f,q,c,1,le]/2.0),
                             (pre+'_'+frame+'_MSRD_perp2_end_half_(A2)',end[2],Merr[f,q,c,2,le]/2.0)]
    print('   sigma2 : about the mean positions of each life (EXAFS : sigma2_par the Debye-Waller factor of the bond, sigma2_perp its transverse part) ;')
    print('   lives short against the vibrations give them too low : compare with MSRD/2 at the last lag with at least %d %% of the pairs of lag 0' % (100*ENDPAIRS))
    print('   (the end of the window, or before when the lives are shorter), 2 sigma2 once decorrelated.')
    print('   lab sigma2_perp includes the libration (and, in a liquid, the rotation) of the polyhedron ; body removes its rigid rotation.')
    with open(names[1],'w') as f:
        f.write('# msrd_umd.py, '+UMDFile+', bonding file '+BondFile+', polyhedra '+species+'\n')
        f.write('# '+str(nsnap)+' snapshots every '+str(TimeStep)+' fs ; gaps up to '+str(maxgap)+' fs bridged, lives of at least '+str(minlife)+' fs ; window '+str(nlags)+' snapshots ; errors from '+str(nblocks)+' blocks\n')
        f.write('quantity\tvalue\terror\n')
        for name,value,err in rows:
            f.write(name+'\t'+str(value)+'\t'+str(err)+'\n')
    if pdf:
        import matplotlib
        if not view:
            matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        fig, draw = plot_msrd(names[0],names[1],central,adjacents,os.path.basename(base),interactive=False)
        with PdfPages(names[2]) as out:
            for q in range(nel):
                draw(q)
                out.savefig(fig,facecolor=SURFACE)
        plt.close(fig)
    print('Files',', '.join(os.path.basename(n) for n in (names if pdf else names[:2])),'written')
    print('runtime : %.1f s' % (time.time()-start))
    if view:
        import matplotlib.pyplot as plt
        fig, draw = plot_msrd(names[0],names[1],central,adjacents,os.path.basename(base))
        plt.show()


if __name__ == "__main__":
    if not umdpf.run_on_all(sys.argv[1:],main,paired={'-b':'.bonding.dat'}):
        main(sys.argv[1:])
