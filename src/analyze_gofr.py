#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### Authors : Natalia Solomatova, Anaïs Kobsch, Razvan Caracas
### Merge of analyze_gofr.py, analyze_gofr_semi_automatic.py and analyze_1gofr.py ; automatic detection of the first shell : Razvan Caracas
### First maximum, first minimum, coordination number and mean distance of the first shell of every pair of elements, from the gofr files of gofr_umd.py


import sys, getopt, os, glob
import numpy as np
from scipy.ndimage import gaussian_filter1d, median_filter
from scipy.signal import find_peaks


NOISE_PEAKS = 10.0      #a peak must be more prominent than this many times the noise of the smoothed g(r)
TOP_MAX = 0.25          #the maximum is fitted over the top quarter of the first peak
BOTTOM_MIN = 0.1        #the minimum is fitted over the bottom tenth of the dip after the first peak
SHALLOW = 0.5           #warning when g(min)/g(max) is above this : the end of the first shell is ill-defined


def read_gofr(gofrfile):#Pairs, distances, g(r) and running coordination numbers of a gofr file of gofr_umd.py
    with open(gofrfile,'r') as f:
        header = f.readline().split()
    data = np.loadtxt(gofrfile,skiprows=1,ndmin=2)
    return header[1::2], data[:,0], data[:,1::2], data[:,2::2]


def smoothed_noise(g,s):#Local noise of the smoothed g(r) : from the median of |second differences| of g over +/- 25 bins (robust to the
    d2 = np.abs(np.diff(g,2))                                    #curvature of sharp peaks), reduced by the gaussian of s bins
    d2 = np.concatenate(([d2[0]],d2,[d2[-1]]))
    sigma = 1.4826*median_filter(d2,size=51,mode='nearest')/np.sqrt(6.0)
    return sigma/np.sqrt(2.0*np.sqrt(np.pi)*max(s,0.5))


def refine_max(r,g,gs,ip,im):#Position and height of the first maximum : cubic fitted to g over the top quarter of the peak
    h = gs[ip]-gs[im]
    lo = ip
    while lo > 0 and gs[lo-1] >= gs[ip]-TOP_MAX*h:
        lo -= 1
    hi = ip
    while hi < im and gs[hi+1] >= gs[ip]-TOP_MAX*h:
        hi += 1
    lo, hi = min(lo,ip-2), max(hi,ip+2)
    if lo < 0 or hi >= len(r):
        return r[ip], gs[ip]
    p = np.polyfit(r[lo:hi+1]-r[ip],g[lo:hi+1],3)
    x = np.linspace(r[lo],r[hi],2001)
    y = np.polyval(p,x-r[ip])
    k = np.argmax(y)
    if k == 0 or k == len(x)-1:#no maximum inside the window
        return r[ip], gs[ip]
    return x[k], y[k]


def refine_min(r,g,gs,ip,im,iend):#Position of the first minimum : parabola fitted to g over the bottom tenth of the minimum, else the middle of this region
    h = gs[ip]-gs[im]
    lo = im
    while lo > ip and gs[lo-1] <= gs[im]+BOTTOM_MIN*h:
        lo -= 1
    hi = im
    while hi < iend and gs[hi+1] <= gs[im]+BOTTOM_MIN*h:
        hi += 1
    lo, hi = min(lo,im-2), max(hi,im+2)
    if lo <= ip or hi >= len(r):
        return r[im]
    p = np.polyfit(r[lo:hi+1]-r[im],g[lo:hi+1],2)
    if p[0] > 0:
        x = r[im]-0.5*p[1]/p[0]
        if r[lo] < x < r[hi]:
            return x
    return 0.5*(r[lo]+r[hi])#flat bottom (e.g. g = 0 between the shells of a molecule) : any cutoff inside gives the same bonds


def shell_values(r,g,gs,coord,xmax,ymax,xmin):#Coordination number and mean distance of the first shell ending at xmin
    dr = r[1]-r[0]
    inshell = r < xmin
    r2g = (r[inshell]**2*g[inshell]).sum()
    return dict(xmax=xmax,ymax=ymax,xmin=xmin,ymin=float(np.interp(xmin,r,gs)),
                coord=float(np.interp(xmin,r+0.5*dr,coord)),#Int(A-B) of gofr_umd.py runs up to the outer edge of each bin
                rmean=(r[inshell]**3*g[inshell]).sum()/r2g if r2g > 0 else np.nan)


def first_shell(r,g,coord,sigma):#First maximum and minimum of one g(r) ; returns the shell values and a note, or None and the reason
    dr = r[1]-r[0]
    if not np.any(g > 0):
        return None, 'no pairs'
    s = sigma/dr
    gs = gaussian_filter1d(g,s,mode='nearest') if s > 0 else g.copy()
    noise = smoothed_noise(g,s)
    peaks, props = find_peaks(gs,prominence=0)
    peaks = peaks[props['prominences'] > NOISE_PEAKS*noise[peaks]]
    first = peaks[gs[peaks] > 1.0]#a shell is denser than the average : g > 1 at its maximum
    if len(first) == 0:
        return None, 'no peak above the noise and above 1'
    ip = first[0]
    later = peaks[peaks > ip]
    iend = later[0] if len(later) > 0 else len(r)-1
    im = ip+int(np.argmin(gs[ip:iend+1]))
    if im >= len(r)-3:
        return None, 'no minimum up to '+'%.2f' % r[-1]+' A'
    xmax, ymax = refine_max(r,g,gs,ip,im)
    xmin = refine_min(r,g,gs,ip,im,iend)
    values = shell_values(r,g,gs,coord,xmax,ymax,xmin)
    note = ''
    if values['ymin'] > SHALLOW*values['ymax']:
        note = 'shallow minimum'
    return values, note


def curves(an,pair):#g(r), its smoothed curve and the running coordination number of a pair
    k = an['pairs'].index(pair)
    r, g = an['r'], an['G'][:,k]
    gs = gaussian_filter1d(g,an['sigma']/(r[1]-r[0]),mode='nearest') if an['sigma'] > 0 else g.copy()
    return r, g, gs, an['C'][:,k]


def set_shell(an,pair,values,note):#Sets the first shell of a pair and of its reverse : same g(r), so same positions, each its own coordination number
    a, b = pair.split('-')
    for p in (pair,b+'-'+a):
        if p not in an['results'] or (p != pair and an['results'][p]['note'] == 'not analyzed'):
            continue
        if values is None:
            an['results'][p] = dict(values=None,note=note)
        else:
            r, g, gs, coord = curves(an,p)
            an['results'][p] = dict(values=shell_values(r,g,gs,coord,values['xmax'],values['ymax'],values['xmin']),note=note)


def analyze_file(gofrfile,sigma,couples):#First shell of every pair of a gofr file
    pairs, r, G, C = read_gofr(gofrfile)
    an = dict(file=gofrfile,pairs=pairs,r=r,G=G,C=C,sigma=sigma,results={})
    for pair in pairs:
        a, b = pair.split('-')
        if couples and pair not in couples and b+'-'+a not in couples:
            an['results'][pair] = dict(values=None,note='not analyzed')
        elif pair not in an['results']:
            r, g, gs, coord = curves(an,pair)
            an['results'][pair] = {}
            set_shell(an,pair,*first_shell(r,g,coord,sigma))
    return an


def print_results(an):
    print('\n',an['file'])
    print('%-8s %8s %8s %8s %8s %8s %8s  %s' % ('pair','xmax(A)','g(xmax)','xmin(A)','g(xmin)','coord','rmean(A)','note'))
    for pair in an['pairs']:
        v, note = an['results'][pair]['values'], an['results'][pair]['note']
        if v is None:
            print('%-8s %8s %8s %8s %8s %8s %8s  %s' % (pair,'-','-','-','-','-','-',note))
        else:
            print('%-8s %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f  %s' % (pair,v['xmax'],v['ymax'],v['xmin'],v['ymin'],v['coord'],v['rmean'],note))


def write_inputs(an):#The three files read by bonding_umd.py -i : the first minima (bond cutoffs), the first maxima and the coordination numbers
    base = an['file'][:-9] if an['file'].endswith('.gofr.dat') else os.path.splitext(an['file'])[0]
    for suffix, key in (('.maximum_bonds.inp','xmin'),('.average_bonds.inp','xmax'),('.coordinations.inp','coord')):
        with open(base+suffix,'w') as f:
            for pair in an['pairs']:
                v = an['results'][pair]['values']
                f.write('\t'.join(pair.split('-'))+'\t'+('%.4f' % v[key] if v is not None else '0')+'\n')
    print('Files ',base+'.maximum_bonds.inp, .average_bonds.inp and .coordinations.inp written')


def write_summary(analyses,dirpath):#One line per gofr file of the folder ; the columns are the union of the pairs of all the files
    pairs = []
    for an in analyses:
        pairs += [p for p in an['pairs'] if p not in pairs]
    keys = ('xmax','ymax','xmin','coord','ymin','rmean')
    summary = os.path.join(dirpath,'gofrs_summary.txt')
    with open(summary,'w') as f:
        f.write('file\t'+'\t'.join(p for p in pairs for _ in keys)+'\n')
        f.write('file\t'+'\t'.join(k for _ in pairs for k in keys)+'\n')
        for an in analyses:
            line = [os.path.basename(an['file'])]
            for p in pairs:
                v = an['results'][p]['values'] if p in an['results'] else None
                line += ['%.4f' % v[k] for k in keys] if v is not None else ['nan']*len(keys)
            f.write('\t'.join(line)+'\n')
    print('File ',summary,' written')


def draw_pair(ax,an,pair):#g(r), its smoothed curve, the first maximum and minimum ; coordination numbers and mean distance in the title
    r, g, gs, coord = curves(an,pair)
    v, note = an['results'][pair]['values'], an['results'][pair]['note']
    ax.clear()
    ax.plot(r,g,'o',markersize=2,color='0.6',markerfacecolor='white',label='g(r)')
    if an['sigma'] > 0:
        ax.plot(r,gs,'-',color='k',linewidth=1,label='smoothed')
    rev = pair.split('-')[1]+'-'+pair.split('-')[0]
    title = os.path.basename(an['file'])+'   '+pair
    if v is not None:
        ax.plot(v['xmax'],v['ymax'],'o',color='tab:blue',markersize=7,label='maximum %.3f A' % v['xmax'])
        ax.plot(v['xmin'],v['ymin'],'o',color='tab:red',markersize=7,label='minimum %.3f A' % v['xmin'])
        ax.axvline(v['xmin'],color='tab:red',linestyle='--',linewidth=0.8)
        title += '\ncoord '+pair+' = %.3f' % v['coord']
        if rev != pair and rev in an['results'] and an['results'][rev]['values'] is not None:
            title += ',  '+rev+' = %.3f' % an['results'][rev]['values']['coord']
        title += ',  mean distance %.3f A' % v['rmean']
    else:
        title += '\nno first shell : '+note
    if note and v is not None:
        title += '  ('+note+')'
    ax.set_title(title,fontsize='medium')
    ax.set_xlabel('distance (A)')
    ax.set_ylabel('g(r)')
    ax.set_ylim(-0.02*max(g.max(),1.0),1.05*max(g.max(),1.0))
    ax.grid(alpha=0.3)
    ax.legend(loc='upper right',fontsize='small')


def unique_pairs(an):#A-B and B-A have the same g(r) : one of them is shown
    shown = []
    for pair in an['pairs']:
        a, b = pair.split('-')
        if b+'-'+a not in shown and an['results'][pair]['note'] != 'not analyzed':
            shown.append(pair)
    return shown


def clicked(an,pair,x,kind,snap):#Maximum or minimum moved to a click : with snap, the extremum of the smoothed g(r) within 0.1 A, refined
    r, g, gs, coord = curves(an,pair)                              #as the automatic ones ; without, the clicked distance
    i = int(np.argmin(np.abs(r-x)))
    if snap:
        near = np.where(np.abs(r-x) <= 0.1)[0]
        i = near[np.argmax(gs[near])] if kind == 'max' else near[np.argmin(gs[near])]
    v = an['results'][pair]['values']
    if kind == 'max':
        im = int(np.argmin(np.abs(r-v['xmin']))) if v is not None else len(r)-1
        xmax, ymax = refine_max(r,g,gs,i,max(im,i+1)) if snap else (x,float(np.interp(x,r,gs)))
        return dict(xmax=xmax,ymax=ymax,xmin=v['xmin'] if v is not None else np.nan)
    ip = int(np.argmin(np.abs(r-v['xmax']))) if v is not None else 0
    xmin = refine_min(r,g,gs,min(ip,i-1),i,len(r)-1) if snap else x
    return dict(xmax=v['xmax'] if v is not None else np.nan,ymax=v['ymax'] if v is not None else np.nan,xmin=xmin)


def check_plots(analyses):#One plot per pair, with Previous/Next ; left click moves the maximum, right click the minimum
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Button
    items = [(an,pair) for an in analyses for pair in unique_pairs(an)]
    if not items:
        return
    fig, ax = plt.subplots(figsize=(9,6))
    fig.subplots_adjust(bottom=0.2)
    state = {'i':0,'snap':True,'pending':{}}
    def draw():
        an, pair = items[state['i']]
        draw_pair(ax,an,pair)
        fig.suptitle('%d/%d    left click: maximum, right click: minimum' % (state['i']+1,len(items)),x=0.98,ha='right',fontsize='small')
        fig.canvas.draw_idle()
    def step(delta):
        state['i'] = (state['i']+delta) % len(items)
        state['pending'] = {}
        draw()
    def onclick(event):
        if event.inaxes is not ax or event.button not in (1,3) or (fig.canvas.toolbar is not None and fig.canvas.toolbar.mode != ''):
            return
        an, pair = items[state['i']]
        kind = 'max' if event.button == 1 else 'min'
        new = clicked(an,pair,event.xdata,kind,state['snap'])
        if an['results'][pair]['values'] is None:#no shell yet : both the maximum and the minimum are needed
            state['pending'][kind] = new
            if len(state['pending']) < 2:
                ax.axvline(event.xdata,color='tab:blue' if kind == 'max' else 'tab:red',linewidth=0.8)
                fig.canvas.draw_idle()
                return
            new = dict(xmax=state['pending']['max']['xmax'],ymax=state['pending']['max']['ymax'],xmin=state['pending']['min']['xmin'])
            state['pending'] = {}
        set_shell(an,pair,new,'set by hand')
        draw()
    def discard(event):
        an, pair = items[state['i']]
        set_shell(an,pair,None,'discarded')
        draw()
    def automatic(event):
        an, pair = items[state['i']]
        r, g, gs, coord = curves(an,pair)
        set_shell(an,pair,*first_shell(r,g,coord,an['sigma']))
        draw()
    def togglesnap(event):
        state['snap'] = not state['snap']
        bsnap.label.set_text('Snap: on' if state['snap'] else 'Snap: off')
        fig.canvas.draw_idle()
    buttons = []
    for label, pos, func in (('Discard',0.10,discard),('Automatic',0.26,automatic),('Snap: on',0.42,togglesnap),
                             ('Previous',0.62,lambda e: step(-1)),('Next',0.78,lambda e: step(1))):
        bt = Button(fig.add_axes([pos,0.04,0.14,0.07]),label)
        bt.on_clicked(func)
        buttons.append(bt)
    bsnap = buttons[2]
    fig._umd_buttons = buttons      #keeps the buttons alive, otherwise they are garbage-collected and stop responding
    fig.canvas.mpl_connect('button_press_event',onclick)
    draw()
    plt.show()


def save_plots(an):#One png per pair
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    base = an['file'][:-9] if an['file'].endswith('.gofr.dat') else os.path.splitext(an['file'])[0]
    for pair in unique_pairs(an):
        fig, ax = plt.subplots(figsize=(8,5))
        draw_pair(ax,an,pair)
        fig.tight_layout()
        fig.savefig(base+'.gofr.'+pair+'.png',dpi=150)
        plt.close(fig)
    print('Plots ',base+'.gofr.<pair>.png written')


def main(argv):
    gofrfile = ''
    couples = []
    sigma = 0.03
    view = False
    plots = False
    usage = 'analyze_gofr.py [-f <gofr file>] [-a <pairs>] [-g <smoothing>] [-v] [-p]'
    try:
        opts, arg = getopt.getopt(argv,"hf:a:g:vp",["fgofrfile","aatoms","gsmoothing"])
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('First maximum, first minimum, coordination number and mean distance of the first shell of every pair of elements,')
            print('from the gofr files of gofr_umd.py. Usage:')
            print(usage)
            print('-f : one gofr file. Default : every *.gofr.dat in the current folder and below it, with one summary per folder.')
            print('-a : pairs to analyze, e.g. Si-O,Mg-O (B-A is included with A-B). The other pairs get 0 in the .inp files, i.e. no bond.')
            print('     Default : all the pairs.')
            print('-g : sigma of the gaussian smoothing used to find the peaks, in A. Default 0.03. 0 : no smoothing. The positions are then')
            print('     fitted to the g(r) itself. Increase it for noisy g(r) (few atoms of an element, short runs).')
            print('-v : check each pair in a plot (Previous/Next). Left click moves the maximum (blue), right click the minimum (red);')
            print('     with Snap on they go to the extremum of the smoothed g(r) within 0.1 A of the click, fitted as the automatic ones,')
            print('     with Snap off they stay where clicked. For a pair without a shell, click both. Discard : no shell for this pair')
            print('     (0 in the .inp files). Automatic : back to the automatic values. A-B and B-A are changed together. Clicks are ignored')
            print('     while the zoom or pan of the toolbar is active. Close the window to write the files.')
            print('-p : save one plot per pair, <base>.gofr.<pair>.png.')
            print('Method : the g(r) is smoothed; its first peak more prominent than '+str(NOISE_PEAKS)+' times the local noise (from the scatter')
            print('of the g(r) around the smoothed curve) is the first shell, and its first minimum is the lowest point before the next such')
            print('peak. The maximum is a cubic fitted to the g(r) over the top quarter of the peak; the minimum a parabola over its bottom')
            print('tenth, or the middle of this region when it is flat (g = 0 between the shells of a molecule: any cutoff inside it gives')
            print('the same bonds). The coordination number is Int(A-B) at the minimum; the mean distance is the mean of r over the shell,')
            print('int r^3 g / int r^2 g up to the minimum. A minimum higher than '+str(SHALLOW)+' times the maximum is flagged: there the end of')
            print('the shell, hence the bonds, depend on the cutoff. Pairs without a peak above the noise, or without a minimum before the')
            print('largest distance of the g(r), get 0. Check the noisy pairs with -v.')
            print('Outputs, for each gofr file <base>.gofr.dat (A B value, the format read by bonding_umd.py -i):')
            print('   <base>.maximum_bonds.inp : first minima, the bond cutoffs for bonding_umd.py and the speciation.')
            print('   <base>.average_bonds.inp : first maxima (most probable bond lengths).')
            print('   <base>.coordinations.inp : coordination numbers.')
            print('and gofrs_summary.txt in each folder : one line per file, for each pair xmax ymax xmin coord ymin rmean (nan: no shell).')
            sys.exit()
        elif opt in ("-f","--fgofrfile"):
            gofrfile = str(arg)
        elif opt in ("-a","--aatoms"):
            couples = arg.split(',')
        elif opt in ("-g","--gsmoothing"):
            sigma = float(arg)
        elif opt == '-v':
            view = True
        elif opt == '-p':
            plots = True
    if gofrfile != '':
        if not os.path.isfile(gofrfile):
            print('gofr file ',gofrfile,' does not exist')
            sys.exit()
        folders = [(os.path.dirname(gofrfile) or os.curdir,[gofrfile])]
    else:
        folders = []
        for dirpath, dirnames, filenames in os.walk(os.curdir):
            files = sorted(glob.glob(os.path.join(dirpath,'*.gofr.dat')))
            if files:
                folders.append((dirpath,files))
        if not folders:
            print('no *.gofr.dat file in the current folder or below it')
            sys.exit()
    analyses = []
    for dirpath, files in folders:
        for f in files:
            an = analyze_file(f,sigma,couples)
            print_results(an)
            analyses.append((dirpath,an))
    if view:
        check_plots([an for _, an in analyses])
        for _, an in analyses:#the values after the checks
            print_results(an)
    for dirpath, files in folders:
        group = [an for d, an in analyses if d == dirpath]
        for an in group:
            write_inputs(an)
            if plots:
                save_plots(an)
        write_summary(group,dirpath)


if __name__ == "__main__":
    main(sys.argv[1:])
