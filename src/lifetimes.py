#!/usr/bin/env python3
# coding: utf-8
### Ana Anzulović (first version), Razvan Caracas
### Lifetimes of every individual species of popul files (speciation_and_angles.py -r 1), in a grid of panels :
### one column per central element, one row per number of adjacent atoms (MO_x), each species sorted by decreasing lifetime ;
### the single atoms in row 0, the free adjacent atoms in columns of their own ; optionally split into liquid and vapor

import sys, getopt, glob, os, re
import numpy as np
import liquid_and_vapor_speciation as lvs

COLORS = ['#2a78d6','#eb6834','#1baf7a','#eda100','#e87ba4','#008300','#6250d6','#e34948']#one per central element, in this fixed order
TEXT = '#0b0b0b'
MUTED = '#52514e'


def read_popul(populfile):#Central and adjacent elements of the header (None if absent), and the species : formula, lifetime (fs)
    centrals = adjacents = None
    species = []
    with open(populfile,'r') as f:
        lines = f.readlines()
    start = 0
    for k,line in enumerate(lines[:3]):
        found = re.search(r'-c=(\S+)',line)
        if found: centrals = found.group(1).split(',')
        found = re.search(r'-a=(\S+)',line)
        if found: adjacents = found.group(1).split(',')
        if line.startswith('Formula'): start = k+1
    for line in lines[start:]:
        entry = line.rstrip('\n').split('\t')
        if len(entry) >= 5:
            species.append((entry[0],float(entry[3])))
    return centrals, adjacents, species


def composition(formula):#'Mg_1O_4' -> {'Mg': 1, 'O': 4}
    return {el: int(n) for el,n in re.findall(r'([A-Za-z]+)_(\d+)',formula)}


def sort_species(species,centrals,minlife):#Lifetimes by (element, number of adjacent atoms), and the formulas : species with one central atom, and single atoms (x = 0)
    lifetimes = {}
    formulas = {}
    skipped = 0
    for formula,life in species:
        if life < minlife: continue
        comp = composition(formula)
        if sum(comp.values()) == 1:                                         #single atom, central or adjacent
            key = (list(comp)[0],0)
            lifetimes.setdefault(key,[]).append(life)
            formulas.setdefault(key,set()).add(formula)
            continue
        central = [el for el in comp if el in centrals]
        if len(central) != 1 or comp[central[0]] != 1:
            skipped += 1                                                    #several central atoms, or several atoms and no central one
            continue
        key = (central[0],sum(comp.values())-1)
        lifetimes.setdefault(key,[]).append(life)
        formulas.setdefault(key,set()).add(formula)
    return lifetimes, formulas, skipped


def nice(v):#Round upper limit of an axis
    e = 10**np.floor(np.log10(v))
    for m in (1,1.2,1.5,2,2.5,3,4,5,6,8,10):
        if m*e >= v*0.999: return m*e
    return 10*e


def label(formula,central):#'Mg_1O_4' -> MgO$_4$, the central element first
    comp = composition(formula)
    text = central
    for el,n in comp.items():
        if el != central:
            text += el+('$_{%d}$' % n if n > 1 else '')
    return text


def grid(lifetimes,centrals,adjacents,rows):#Columns (central elements present, then free adjacent atoms) and rows (numbers of adjacent atoms) of a figure
    columns = [el for el in centrals if any(key[0] == el for key in lifetimes)]
    free = sorted(key[0] for key in lifetimes if key[0] not in centrals)
    free.sort(key=lambda el: adjacents.index(el) if adjacents and el in adjacents else len(adjacents or []))
    columns += free
    if rows is None:
        xs = [key[1] for key in lifetimes if key[1] > 0]
        rows = list(range(min(xs),max(xs)+1)) if xs else []
    if any(key[1] == 0 for key in lifetimes) and 0 not in rows:
        rows = [0]+[x for x in rows if x > 0]                               #the single atoms are always shown
    return columns, rows


def draw(fig,title,lifetimes,formulas,centrals,adjacents,layout):#fig : a figure or a subfigure, cleared first ; layout : columns and rows
    fig.clf()
    columns, rows = layout
    if not columns or not rows:
        fig.text(0.5,0.5,'no species',ha='center',color=MUTED)
        fig.suptitle(title,fontsize=10,color=MUTED)
        return
    axes = fig.subplots(len(rows),len(columns),squeeze=False)
    w,h = fig.figure.get_size_inches() if hasattr(fig,'figure') else fig.get_size_inches()
    fig.subplots_adjust(left=1.0/w,right=1-0.6/w,bottom=0.75/h,top=1-0.8/h,wspace=0.45,hspace=0.45)#margins in inches
    adj = adjacents[0] if adjacents and len(adjacents) == 1 else 'X'
    order = centrals+[el for el in (adjacents or []) if el not in centrals]  #colour of each element, from the header
    for j,el in enumerate(columns):
        free = el not in centrals                                           #column of the free atoms of an adjacent element : row 0 only
        color = COLORS[(order.index(el) if el in order else len(order)+j) % len(COLORS)]
        for i,x in enumerate(rows):
            ax = axes[i,j]
            if free and x != 0:
                ax.axis('off')
                if j == len(columns)-1:
                    ax.text(1.12,0.5,str(x),transform=ax.transAxes,ha='left',va='center',fontsize=11,color=TEXT)
                continue
            for side in ('top','right'):
                ax.spines[side].set_visible(False)
            for side in ('left','bottom'):
                ax.spines[side].set_color(MUTED)
            ax.tick_params(colors=MUTED,labelsize=8,length=2)
            lives = np.sort(lifetimes.get((el,x),[]))[::-1]
            if len(lives):
                ax.stairs(lives,np.arange(len(lives)+1),fill=True,color=color,linewidth=0)
                top = nice(lives[0])
                ax.set_xlim(0,len(lives))
                ax.set_ylim(0,top)
                ax.set_xticks([0,len(lives)])
                ax.set_yticks([0,top])
                ax.set_yticklabels(['0','%g' % top])
                names = sorted(formulas[(el,x)])
                name = label(names[0],el) if len(names) == 1 else el+'X$_{%d}$' % x#X : several adjacent elements
                avg = lives.mean()
            else:
                ax.set_xlim(0,10)
                ax.set_ylim(0,10)
                ax.set_xticks([0,10])
                ax.set_yticks([0,10])
                name = el+(adj+'$_{%d}$' % x if x > 1 else adj if x == 1 else '')
                avg = 0.0
            ax.text(0.97,0.95,name+'\n'+'%.0f fs' % avg,transform=ax.transAxes,ha='right',va='top',fontsize=9,color=TEXT,
                    bbox=dict(boxstyle='square,pad=0.15',facecolor='white',edgecolor='none',alpha=0.8))
            if i == 0:
                ax.set_title(el if not free else 'free '+el,fontsize=12,color=TEXT)
            if j == len(columns)-1:
                ax.text(1.12,0.5,str(x),transform=ax.transAxes,ha='left',va='center',fontsize=11,color=TEXT)
            if i == len(rows)-1 and not free:
                ax.set_xlabel(el+(adj+'$_x$' if adj != 'X' else '')+' species',fontsize=10,color=TEXT)
    fig.text(0.3/w,0.5,'Lifetime of species (fs)',rotation=90,ha='center',va='center',fontsize=12,color=TEXT)
    fig.suptitle(title,fontsize=10,color=MUTED,y=1-0.15/h)


def add_figures(figures,title,outfile,sets,phases,centrals,adjacents,minlife,rows):#One figure per phase, all with the same columns and rows
    sorted_sets = [sort_species(sets[phase],centrals,minlife) for phase in phases]
    union = {key: True for lifetimes,formulas,skipped in sorted_sets for key in lifetimes}
    layout = grid(union,centrals,adjacents,rows)
    for phase,(lifetimes,formulas,skipped) in zip(phases,sorted_sets):
        if not phase:
            phaseout = outfile
        elif outfile.endswith('.lifetimes.pdf'):#<popul base>.liquid.lifetimes.pdf
            phaseout = outfile[:-14]+'.'+phase+'.lifetimes.pdf'
        else:#merged : <output>.liquid.pdf
            phaseout = os.path.splitext(outfile)[0]+'.'+phase+'.pdf'
        figures.append((title+(' : '+phase if phase else ''),phaseout,lifetimes,formulas,centrals,adjacents,layout))
        print(title+(' '+phase if phase else ''),':',sum(len(v) for v in lifetimes.values()),'species,',skipped,'with several central atoms left out')


def main(argv):
    files = []
    stepfiles = []
    centralsopt = None
    minlife = 0.0
    rows = None
    length = None
    merge = False
    view = False
    output = ''
    usage = 'lifetimes.py -f <popul files> -c <central elements> -t <min lifetime> -x <xmin,xmax> -l <length> -s <r0 step files> -m -o <output> -v'
    try:
        opts, args = getopt.getopt(argv,"hf:c:t:x:l:s:mo:v")
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('Lifetimes of every individual species of popul files of speciation_and_angles.py -r 1, in a grid of panels : one column per central')
            print('element, one row per number of adjacent atoms x (MO_x), the species sorted by decreasing lifetime, with their mean lifetime. Usage:')
            print(usage)
            print('-f : popul files, separated by commas. Default : every *.r1.popul.dat of the current folder.')
            print('-c : central elements (columns), separated by commas. Default : the -c of the popul header (needed for the vapor/liquid popul files')
            print('     of liquid_and_vapor_speciation.py, which have no header).')
            print('-t : species living less than t fs are left out. Default 0.')
            print('-x : rows from xmin to xmax adjacent atoms, e.g. -x 3,9. Default : all the x present. The single atoms (row 0) are always shown.')
            print('-l : splits the species into liquid and vapor, as liquid_and_vapor_speciation.py : an atom is in the vapor when its r0 cluster has')
            print('     fewer than l atoms, a species is in the liquid when one of its atoms is. One figure for the liquid and one for the vapor.')
            print('-s : with -l, the r0 step files (speciation_and_angles.py -r 0, same bonding file and -s), one per popul file, separated by commas.')
            print('     Default : <prefix>.r0.step.dat for each <prefix>.r<r>.popul.dat.')
            print('-m : merges the files into one figure, for independent runs at the same P and T. Default : one figure per file.')
            print('-o : name of the merged figure (with -m). Default lifetimes.merged.pdf.')
            print('-v : shows the figures, with Previous/Next buttons.')
            print('Output : <popul base>.lifetimes.pdf per file (or the merged one) ; with -l, <popul base>.liquid.lifetimes.pdf and .vapor.lifetimes.pdf')
            print('(merged : <output>.liquid.pdf and .vapor.pdf). The single atoms are in row 0 : the central ones in their column, the adjacent ones')
            print('(not bound to a central atom) in columns of their own (free O...). Species with several central atoms are left out. A panel with')
            print('several formulas (several adjacent elements, e.g. CO2 and CH2) is labelled CX2.')
            sys.exit()
        elif opt == '-f':
            files = arg.split(',')
        elif opt == '-c':
            centralsopt = arg.split(',')
        elif opt == '-t':
            minlife = float(arg)
        elif opt == '-x':
            xmin,xmax = [int(v) for v in arg.split(',')]
            rows = list(range(xmin,xmax+1))
        elif opt == '-l':
            length = int(arg)
        elif opt == '-s':
            stepfiles = arg.split(',')
        elif opt == '-m':
            merge = True
        elif opt == '-o':
            output = arg
        elif opt == '-v':
            view = True
    if not files:
        files = sorted(glob.glob('*.r1.popul.dat'))
    if not files:
        print('no popul file : give them with -f')
        sys.exit()
    if length is not None:
        if not stepfiles:
            stepfiles = [re.sub(r'\.r\d+\.popul\.dat$','.r0.step.dat',populfile) for populfile in files]
        if len(stepfiles) != len(files):
            print('give one r0 step file per popul file with -s')
            sys.exit()
        for stepfile in stepfiles:
            if not os.path.isfile(stepfile):
                print('the r0 step file',stepfile,'does not exist : give it with -s')
                sys.exit()
    phases = ['liquid','vapor'] if length is not None else ['']

    import matplotlib
    if not view:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    figures = []#(title, output file, lifetimes, formulas, centrals, adjacents)
    pooled = {phase: [] for phase in phases}
    allcentrals = []
    alladjacents = []
    for k,populfile in enumerate(files):
        centrals, adjacents, species = read_popul(populfile)
        centrals = centralsopt or centrals
        if not centrals:
            print('no central elements in the header of',populfile,': give them with -c')
            sys.exit()
        if length is not None:
            vapor, liquid = lvs.split(populfile,stepfiles[k],length)
            sets = {'liquid': [(sp[0],sp[3]) for sp in liquid], 'vapor': [(sp[0],sp[3]) for sp in vapor]}
        else:
            sets = {'': species}
        allcentrals += [el for el in centrals if el not in allcentrals]
        alladjacents += [el for el in (adjacents or []) if el not in alladjacents]
        base = populfile[:-10] if populfile.endswith('.popul.dat') else os.path.splitext(populfile)[0]
        if merge:
            for phase in phases:
                pooled[phase] += sets[phase]
            continue
        add_figures(figures,populfile,base+'.lifetimes.pdf',sets,phases,centrals,adjacents,minlife,rows)
    if merge:
        add_figures(figures,'merged : '+', '.join(files),output or 'lifetimes.merged.pdf',pooled,phases,allcentrals,alladjacents or None,minlife,rows)

    for title,outfile,lifetimes,formulas,centrals,adjacents,layout in figures:
        columns, figrows = layout
        fig = plt.figure(figsize=(1.9*max(1,len(columns))+1.6,1.35*max(1,len(figrows))+1.2))
        draw(fig,title,lifetimes,formulas,centrals,adjacents,layout)
        fig.savefig(outfile)
        plt.close(fig)
        print('figure written in',outfile)

    if view:
        from matplotlib.widgets import Button
        state = {'k': 0}
        fig = plt.figure(figsize=(12,9))
        panel = fig.subfigures(1,1)#redrawn at each step ; the buttons stay on fig
        def show():
            title,outfile,lifetimes,formulas,centrals,adjacents,layout = figures[state['k']]
            draw(panel,title+'  (%d/%d)' % (state['k']+1,len(figures)),lifetimes,formulas,centrals,adjacents,layout)
            fig.canvas.draw_idle()
        def step(d):
            state['k'] = (state['k']+d) % len(figures)
            show()
        bprev = Button(fig.add_axes([0.70,0.01,0.1,0.04]),'Previous')
        bnext = Button(fig.add_axes([0.81,0.01,0.1,0.04]),'Next')
        bprev.on_clicked(lambda event: step(-1))
        bnext.on_clicked(lambda event: step(1))
        fig._umd_buttons = (bprev,bnext)
        show()
        plt.show()


if __name__ == "__main__":
    main(sys.argv[1:])
