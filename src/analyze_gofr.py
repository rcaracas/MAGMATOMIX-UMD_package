#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### Authors : Natalia Solomatova, Anaïs Kobsch, Razvan Caracas
### Merge of analyze_gofr.py, analyze_gofr_semi_automatic.py and analyze_1gofr.py ; automatic detection of the first shell : Razvan Caracas
### First maximum, first minimum, coordination number and mean distance of the first shell of every pair of elements, from the gofr files of gofr_umd.py


import sys, getopt, os, glob, re
import numpy as np
from scipy.ndimage import gaussian_filter1d, median_filter
from scipy.signal import find_peaks
import umd_processes_fast as umdpf


NOISE_PEAKS = 10.0      #a peak must be more prominent than this many times the noise of the smoothed g(r)
TOP_MAX = 0.25          #the maximum is fitted over the top quarter of the first peak
BOTTOM_MIN = 0.1        #the minimum is fitted over the bottom tenth of the dip after the first peak
SHALLOW = 0.5           #warning when g(min)/g(max) is above this : the end of the first shell is ill-defined
MAX_FIT = 1.0           #the cubics are fitted over at most this width (A), centered on the peak or the dip


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


def refine_max(r,g,gs,ip,im):#Position and height of the first maximum : cubic fitted to g over the top quarter of the peak, at most MAX_FIT wide
    half = int(round(0.5*MAX_FIT/(r[1]-r[0])))
    h = gs[ip]-gs[im]
    lo = ip
    while lo > 0 and gs[lo-1] >= gs[ip]-TOP_MAX*h:
        lo -= 1
    hi = ip
    while hi < im and gs[hi+1] >= gs[ip]-TOP_MAX*h:
        hi += 1
    lo, hi = max(min(lo,ip-2),ip-half), min(max(hi,ip+2),ip+half)
    if lo < 0 or hi >= len(r):
        return r[ip], gs[ip]
    p = np.polyfit(r[lo:hi+1]-r[ip],g[lo:hi+1],3)
    x = np.linspace(r[lo],r[hi],2001)
    y = np.polyval(p,x-r[ip])
    k = np.argmax(y)
    if k == 0 or k == len(x)-1:#no maximum inside the window
        return r[ip], gs[ip]
    return x[k], y[k]


def refine_min(r,g,gs,ip,im,iend):#Position of the first minimum : cubic fitted to g over the bottom tenth of the dip, at most MAX_FIT wide ;
    half = int(round(0.5*MAX_FIT/(r[1]-r[0])))                   #the lowest point of the smoothed g(r) when the cubic has no minimum inside
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
    if g[im] == 0:#no pairs at the minimum (g = 0 between the shells of a molecule) : the middle of the gap, any cutoff inside gives the same bonds
        zlo, zhi = im, im
        while zlo > ip and g[zlo-1] == 0:
            zlo -= 1
        while zhi < len(r)-1 and g[zhi+1] == 0:
            zhi += 1
        return 0.5*(r[zlo]+r[zhi])
    lo, hi = max(lo,im-half), min(hi,im+half)
    p = np.polyfit(r[lo:hi+1]-r[im],g[lo:hi+1],3)
    x = np.linspace(r[lo],r[hi],2001)
    y = np.polyval(p,x-r[im])
    k = np.argmin(y)
    if k == 0 or k == len(x)-1:#no minimum inside the window
        return r[im]
    return x[k]


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


def selected(a,b,couples,central,adjacent):#Pair A-B to be analyzed : listed in couples (or B-A), and A central and B adjacent or the reverse
    if couples and a+'-'+b not in couples and b+'-'+a not in couples:
        return False
    return (central == 'all' or a in central) and (adjacent == 'all' or b in adjacent) or \
           (central == 'all' or b in central) and (adjacent == 'all' or a in adjacent)


def analyze_file(gofrfile,sigma,couples,central='all',adjacent='all'):#First shell of every selected pair of a gofr file
    pairs, r, G, C = read_gofr(gofrfile)
    an = dict(file=gofrfile,pairs=pairs,r=r,G=G,C=C,sigma=sigma,results={})
    for pair in pairs:
        a, b = pair.split('-')
        if not selected(a,b,couples,central,adjacent):
            an['results'][pair] = dict(values=None,note='not analyzed')
        elif pair not in an['results']:
            r, g, gs, coord = curves(an,pair)
            an['results'][pair] = {}
            set_shell(an,pair,*first_shell(r,g,coord,sigma))
    return an


def print_results(an):#One line per pair A-B, with the coordination numbers of A-B and of B-A (same g(r), same shell) ;
    print('\n',an['file'])          #the pairs not analyzed or without pairs on one line
    print('%-8s %8s %8s %8s %8s %8s %8s %8s  %s' % ('pair A-B','xmax(A)','g(xmax)','xmin(A)','g(xmin)','CN A-B','CN B-A','rmean(A)','note'))
    others = []
    for pair in an['pairs']:
        a, b = pair.split('-')
        rev = b+'-'+a
        note = an['results'][pair]['note']
        if rev != pair and rev in an['pairs'][:an['pairs'].index(pair)]:#already on the line of B-A
            continue
        if note in ('not analyzed','no pairs'):
            others.append(pair)
            continue
        v = an['results'][pair]['values']
        vrev = an['results'][rev]['values'] if rev != pair and rev in an['results'] else None
        if v is None:
            print('%-8s %8s %8s %8s %8s %8s %8s %8s  %s' % (pair,'-','-','-','-','-','-','-',note))
        else:
            print('%-8s %8.3f %8.3f %8.3f %8.3f %8.3f %8s %8.3f  %s' % (pair,v['xmax'],v['ymax'],v['xmin'],v['ymin'],v['coord'],
                                                                    '%.3f' % vrev['coord'] if vrev is not None else '-',v['rmean'],note))
    if others:
        print('not analyzed or without pairs (0 in the .inp files) :',', '.join(others))


def output_names(gofrfile):#The outputs of analyze_gofr.py for a gofr file : the three .inp files and the plots of -p
    base = gofrfile[:-9] if gofrfile.endswith('.gofr.dat') else os.path.splitext(gofrfile)[0]
    folder, name = os.path.split(base)
    return [base+suffix for suffix in ('.maximum_bonds.inp','.average_bonds.inp','.coordinations.inp')] \
           +umdpf.former_outputs(folder,re.escape(name)+r'\.gofr\.[A-Z][a-z]?-[A-Z][a-z]?\.png')


def write_inputs(an):#The three files read by bonding_umd.py -c : the first minima (bond cutoffs), the first maxima and the coordination numbers ;
    base = an['file'][:-9] if an['file'].endswith('.gofr.dat') else os.path.splitext(an['file'])[0]#every pair of the elements of the file,
    elements = list(dict.fromkeys(el for pair in an['pairs'] for el in pair.split('-')))     #0 (no bond) for those without a shell or not
    for suffix, key in (('.maximum_bonds.inp','xmin'),('.average_bonds.inp','xmax'),('.coordinations.inp','coord')):#in the file
        with open(base+suffix,'w') as f:
            for pair in (a+'-'+b for a in elements for b in elements):
                v = an['results'][pair]['values'] if pair in an['results'] else None
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


def unique_pairs(an):#A-B and B-A have the same g(r) : one of them is shown ; not the pairs not analyzed or without pairs (g = 0 : not
    shown = []                                        #computed by gofr_umd.py -c/-a, or a single atom of an element with itself)
    for pair in an['pairs']:
        a, b = pair.split('-')
        if b+'-'+a not in shown and an['results'][pair]['note'] not in ('not analyzed','no pairs'):
            shown.append(pair)
    return shown


def fit_click(r,gs,g,x,kind,ip=None):#Manual input : the automatic fits started from the clicked peak or dip instead of the one found :
    i = int(np.argmin(np.abs(r-x)))                     #from the click, up the smoothed g(r) to the top of its peak (down to the
    sign = 1 if kind == 'max' else -1                   #bottom of its dip), then the maximum of a cubic over the top quarter of
    while True:                                         #the peak (the minimum of a cubic over the bottom tenth of the dip)
        best = max((k for k in (i-1,i,i+1) if 0 <= k < len(r) and (kind == 'max' or k > ip)),key=lambda k: sign*gs[k])
        if best == i:
            break
        i = best
    if kind == 'max':
        im = i+int(np.argmin(gs[i:])) if i < len(r)-1 else i
        return refine_max(r,g,gs,i,max(im,i+1))
    xmin = refine_min(r,g,gs,min(ip,i-1),i,len(r)-1)
    return xmin, float(np.interp(xmin,r,gs))


def check_plots(analyses):#One plot per pair, with Previous/Next ; two left clicks, the maximum then the minimum, in the mode of the buttons
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Button
    items = [(an,pair) for an in analyses for pair in unique_pairs(an)]
    if not items:
        return
    fig, ax = plt.subplots(figsize=(9,6))
    fig.subplots_adjust(bottom=0.2)
    state = {'i':0,'mode':'manual','clicks':[],'message':''}
    names = {'manual':'Manual input','all manual':'All manual'}
    def draw():
        an, pair = items[state['i']]
        draw_pair(ax,an,pair)
        for x, y in state['clicks']:#the maximum already clicked, waiting for the minimum
            ax.axvline(x,color='tab:blue',linewidth=0.8)
        todo = 'click the minimum' if state['clicks'] else 'click the maximum, then the minimum'
        if state['i'] == len(items)-1:
            todo += '    (last pair : Next closes and writes the files)'
        fig.suptitle('%d/%d    %s : %s%s' % (state['i']+1,len(items),names[state['mode']],todo,
                                            ('    '+state['message']) if state['message'] else ''),x=0.98,ha='right',fontsize='small')
        for mode, bt in (('manual',buttons[2]),('all manual',buttons[3])):
            bt.color = 'lightsteelblue' if state['mode'] == mode else '0.85'
            bt.ax.set_facecolor(bt.color)
        fig.canvas.draw_idle()
    def step(delta):#one pass over the pairs : Next on the last one closes the window (the files are then written), Previous on the first stays
        if state['i']+delta >= len(items):
            plt.close(fig)
            return
        state['i'] = max(0,state['i']+delta)
        state['clicks'], state['message'] = [], ''
        draw()
    def onclick(event):
        if event.inaxes is not ax or event.button != 1 or (fig.canvas.toolbar is not None and fig.canvas.toolbar.mode != ''):
            return
        an, pair = items[state['i']]
        r, g, gs, coord = curves(an,pair)
        state['message'] = ''
        kind = 'min' if state['clicks'] else 'max'
        if state['mode'] == 'manual':
            ip = int(np.argmin(np.abs(r-state['clicks'][0][0]))) if state['clicks'] else None
            x, y = fit_click(r,gs,g,event.xdata,kind,ip)
        else:
            x, y = event.xdata, (event.ydata if kind == 'max' else float(np.interp(event.xdata,r,gs)))
        if kind == 'max':
            state['clicks'] = [(x,y)]
            draw()
            return
        xmax, ymax = state['clicks'][0]
        state['clicks'] = []
        if x <= xmax:
            state['message'] = 'the minimum must be after the maximum : click both again'
        else:
            set_shell(an,pair,dict(xmax=xmax,ymax=ymax,xmin=x),names[state['mode']].lower())
        draw()
    def setmode(mode):
        state['mode'], state['clicks'], state['message'] = mode, [], ''
        draw()
    def discard(event):
        an, pair = items[state['i']]
        set_shell(an,pair,None,'discarded')
        state['clicks'], state['message'] = [], ''
        draw()
    def automatic(event):
        an, pair = items[state['i']]
        r, g, gs, coord = curves(an,pair)
        set_shell(an,pair,*first_shell(r,g,coord,an['sigma']))
        state['clicks'], state['message'] = [], ''
        draw()
    buttons = []
    for label, pos, func in (('Discard',0.04,discard),('Automatic',0.195,automatic),('Manual input',0.35,lambda e: setmode('manual')),
                             ('All manual',0.505,lambda e: setmode('all manual')),('Previous',0.69,lambda e: step(-1)),('Next',0.845,lambda e: step(1))):
        bt = Button(fig.add_axes([pos,0.04,0.135,0.07]),label)
        bt.on_clicked(func)
        buttons.append(bt)
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
    central = 'all'
    adjacent = 'all'
    sigma = 0.03
    view = False
    plots = False
    usage = 'analyze_gofr.py [-f <gofr file>] [-c <central>] [-a <adjacent>] [-g <smoothing>] [-v] [-p]'
    try:
        opts, arg = getopt.getopt(argv,"hf:c:a:g:vp",["fgofrfile","ccentral","aadjacent","gsmoothing"])
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('First maximum, first minimum, coordination number and mean distance of the first shell of every pair of elements,')
            print('from the gofr files of gofr_umd.py. Usage:')
            print(usage)
            print('-f : one gofr file. Default : every *.gofr.dat in the current folder and below it, with one summary per folder.')
            print('-c : central elements, separated by commas, e.g. Si,Mg. Default is all.')
            print('-a : adjacent elements, separated by commas, e.g. O. Default is all. As in gofr_umd.py and speciation_and_angles.py :')
            print('     the pairs central-adjacent are analyzed, and their reverse (same g(r), the coordination seen from the other atom),')
            print('     e.g. -c Si,Mg -a O : Si-O, Mg-O, O-Si, O-Mg ; -c O : O-X and X-O for every element X.')
            print('     -a with pairs, e.g. -a Si-O,Mg-O (the former use, still accepted) : these pairs and their reverse.')
            print('     The other pairs get 0 in the .inp files, i.e. no bond. Default : all the pairs.')
            print('-g : sigma of the gaussian smoothing used to find the peaks, in A. Default 0.03. 0 : no smoothing. The positions are then')
            print('     fitted to the g(r) itself. Increase it for noisy g(r) (few atoms of an element, short runs).')
            print('-v : check each pair in a plot. Two left clicks set the first shell : the first one the maximum (blue), the second one the')
            print('     minimum (red, the bond cutoff), in the mode of the buttons :')
            print('       Manual input (default) : the fits of the automatic analysis, on the peak and the dip clicked instead of those found :')
            print('                    from the first click up the smoothed g(r) to the top of its peak, the maximum of a 3rd-order polynomial')
            print('                    fitted over the top quarter of the peak ; from the second click down to the bottom of its dip, the minimum')
            print('                    of a 3rd-order polynomial over the bottom tenth of the dip (fits at most '+str(MAX_FIT)+' A wide).')
            print('       All manual : the clicked points are taken as they are.')
            print('     Discard : no shell for this pair (0 in the .inp files). Automatic : back to the automatic values. Previous/Next : the')
            print('     other pairs and files, one pass : each pair once (A-B and B-A, which have the same g(r), are changed together), and Next')
            print('     on the last pair closes the window and writes the files (so does closing the window). Clicks are ignored while the zoom')
            print('     or pan of the toolbar is active.')
            print('-p : save one plot per pair, <base>.gofr.<pair>.png.')
            print('Method : the g(r) is smoothed by a gaussian (-g); its first peak more prominent than '+str(NOISE_PEAKS)+' times the local noise (from the scatter')
            print('of the g(r) around the smoothed curve) is the first shell, and its first minimum is the lowest point before the next such')
            print('peak. The maximum is a cubic fitted to the g(r) over the top quarter of the peak; the minimum a cubic over the bottom')
            print('tenth of the dip (both fits at most '+str(MAX_FIT)+' A wide, centered on the peak or the dip; a cubic follows asymmetric dips),')
            print('or the lowest point of the smoothed g(r) when the cubic has no minimum inside, or the middle of the gap when g = 0 there')
            print('(between the shells of a molecule: any cutoff inside it gives the same bonds). The coordination number is Int(A-B) at')
            print('the minimum; the mean distance is the mean of r over the shell, int r^3 g / int r^2 g up to the minimum. A minimum')
            print('higher than '+str(SHALLOW)+' times the maximum is flagged: there the end of the shell, hence the bonds, depend on the cutoff.')
            print('Pairs without a peak above the noise, or without a minimum before the largest distance of the g(r), get 0. Check the')
            print('noisy pairs with -v.')
            print('Outputs, for each gofr file <base>.gofr.dat (A B value, the format read by bonding_umd.py -c):')
            print('   <base>.maximum_bonds.inp : first minima, the bond cutoffs for bonding_umd.py and the speciation.')
            print('   <base>.average_bonds.inp : first maxima (most probable bond lengths).')
            print('   <base>.coordinations.inp : coordination numbers.')
            print('and gofrs_summary.txt in each folder : one line per file, for each pair xmax ymax xmin coord ymin rmean (nan: no shell).')
            print('The .inp files list every pair of the elements of the gofr file, with 0 (no bond) for the pairs without a shell, not')
            print('analyzed, or with g(r) = 0 (not computed by gofr_umd.py -c/-a). The outputs of a former run on the same gofr file (the')
            print('three .inp files and the plots of -p) are removed first.')
            sys.exit()
        elif opt in ("-f","--fgofrfile"):
            gofrfile = str(arg)
        elif opt in ("-c","--ccentral"):
            central = 'all' if arg == 'all' else [el.strip() for el in arg.split(',') if el.strip()]
        elif opt in ("-a","--aadjacent") and '-' in arg:#pairs : the former use of -a
            couples = [pair.strip() for pair in arg.split(',') if pair.strip()]
        elif opt in ("-a","--aadjacent"):
            adjacent = 'all' if arg == 'all' else [el.strip() for el in arg.split(',') if el.strip()]
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
            umdpf.remove_former_outputs(output_names(f))#the outputs of a former run on this gofr file
            an = analyze_file(f,sigma,couples,central,adjacent)
            if all(an['results'][pair]['note'] == 'not analyzed' for pair in an['pairs']):
                print('\n',f,': no pair of the selection (elements of the file :',
                      ','.join(dict.fromkeys(el for pair in an['pairs'] for el in pair.split('-'))),')')
            analyses.append((dirpath,an))
    if view:
        check_plots([an for _, an in analyses])
    for _, an in analyses:#once, after the checks with -v
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
