#!/usr/bin/env python3
"""Graphical abstract of MAGMATOMIX/UMD 3.0 (one 16:9 slide), in the style of the flowchart:

    python3 make_graphical_abstract.py <folder with the workflow outputs of pyrolite.a13.T3000K.outcar.umd.dat> [output.svg]

The small curves of the three panels are the g(r), MSD and vDOS of that run (e.g. the folder of the case
workflow_pyrolite of the validation suite). Default output: umd-graphical-abstract.svg next to this script.
"""
import os, sys, math
import numpy as np

W, H = 1920, 1080
PURPLE, MAGENTA, BLACK, GRAY, LIGHT = '#782F9A', '#ED127B', '#231F20', '#6d6a6b', '#f4eef7'
FONT = "Helvetica, 'Helvetica Neue', Arial, sans-serif"
BASE = 'pyrolite.a13.T3000K.outcar'
out = []

def esc(s):
    return s.replace('&','&amp;').replace('<','&lt;').replace('>','&gt;')

def text(x,y,lines,size,color=BLACK,anchor='start',weight='normal',spacing=1.25,mono=False):
    if isinstance(lines,str):
        lines = [lines]
    font = "Menlo, 'DejaVu Sans Mono', monospace" if mono else FONT
    tspans = ''.join('<tspan x="%g" dy="%g">%s</tspan>' % (x,0 if i == 0 else size*spacing,esc(l)) for i,l in enumerate(lines))
    out.append('<text x="%g" y="%g" font-family="%s" font-size="%g" fill="%s" text-anchor="%s" font-weight="%s">%s</text>'
               % (x,y,font,size,color,anchor,weight,tspans))

def rect(x,y,w,h,stroke=BLACK,fill='white',sw=3,r=14,dash=None):
    d = ' stroke-dasharray="%s"' % dash if dash else ''
    out.append('<rect x="%g" y="%g" width="%g" height="%g" rx="%g" fill="%s" stroke="%s" stroke-width="%g"%s/>' % (x,y,w,h,r,fill,stroke,sw,d))

def arrow(x0,y0,x1,y1,color=BLACK,sw=4):
    a = math.atan2(y1-y0,x1-x0)
    L, Wd = 22, 10
    bx, by = x1-0.8*L*math.cos(a), y1-0.8*L*math.sin(a)
    out.append('<path d="M %g %g L %g %g" stroke="%s" stroke-width="%g" fill="none"/>' % (x0,y0,bx,by,color,sw))
    p = [(x1,y1),(x1-L*math.cos(a)+Wd*math.sin(a),y1-L*math.sin(a)-Wd*math.cos(a)),(x1-L*math.cos(a)-Wd*math.sin(a),y1-L*math.sin(a)+Wd*math.cos(a))]
    out.append('<path d="M %s Z" fill="%s"/>' % (' L '.join('%g %g' % q for q in p),color))

def table(path):
    with open(path) as f:
        lines = [l for l in f if not l.startswith('#')]
    return lines[0].split(), np.loadtxt(lines[1:],ndmin=2)

def curve(x,y,box,xlim=None,ylim=None,color=PURPLE,sw=3.5,logx=False,logy=False,fill=False):
    """Draw y(x) inside box=(x0,y0,w,h); a thin axis line at the bottom and left."""
    bx, by, bw, bh = box
    x = np.asarray(x,float); y = np.asarray(y,float)
    if logx: x = np.log10(x)
    if logy: y = np.log10(y)
    xlim = xlim or (x.min(),x.max())
    ylim = ylim or (min(0,y.min()),y.max()*1.05)
    keep = (x >= xlim[0]) & (x <= xlim[1])
    x, y = x[keep], y[keep]
    X = bx + (x-xlim[0])/(xlim[1]-xlim[0])*bw
    Y = by + bh - (y-ylim[0])/(ylim[1]-ylim[0])*bh
    d = 'M ' + ' L '.join('%.1f %.1f' % q for q in zip(X,Y))
    if fill:
        out.append('<path d="%s L %.1f %.1f L %.1f %.1f Z" fill="%s" fill-opacity="0.12" stroke="none"/>' % (d,X[-1],by+bh,X[0],by+bh,color))
    out.append('<path d="%s" fill="none" stroke="%s" stroke-width="%g" stroke-linejoin="round"/>' % (d,color,sw))

def axes(box):
    bx, by, bw, bh = box
    out.append('<path d="M %g %g L %g %g L %g %g" fill="none" stroke="%s" stroke-width="2"/>' % (bx,by-6,bx,by+bh,bx+bw+6,by+bh,GRAY))

def main(folder,svgfile):
    p = lambda ext: os.path.join(folder,BASE+ext)

    # title
    text(80,108,'MAGMATOMIX/UMD 3.0',80,PURPLE,weight='bold')
    text(82,160,'from molecular-dynamics runs to properties, with their errors',34,BLACK)

    # inputs
    rect(80,250,300,370,BLACK,'white',3)
    text(230,307,['ab initio and','machine-learning','molecular dynamics'],27,BLACK,'middle',weight='bold',spacing=1.2)
    text(230,450,['VASP','ABINIT','LAMMPS','QBox'],29,PURPLE,'middle',spacing=1.3)
    arrow(380,435,462,435)

    # the UMD file (a page with a folded corner) and the checks
    x0, y0, w, h, f = 470, 250, 250, 300, 46
    out.append('<path d="M %g %g L %g %g L %g %g L %g %g L %g %g Z" fill="%s" stroke="%s" stroke-width="3" stroke-linejoin="round"/>'
               % (x0,y0,x0+w-f,y0,x0+w,y0+f,x0+w,y0+h,x0,y0+h,LIGHT,PURPLE))
    out.append('<path d="M %g %g L %g %g L %g %g" fill="none" stroke="%s" stroke-width="3" stroke-linejoin="round"/>' % (x0+w-f,y0,x0+w-f,y0+f,x0+w,y0+f,PURPLE))
    for i in range(6):
        yy = y0+150+i*20
        out.append('<path d="M %g %g L %g %g" stroke="%s" stroke-width="3" stroke-opacity="0.35"/>' % (x0+30,yy,x0+w-30-(i%3)*30,yy,PURPLE))
    text(x0+w/2,y0+85,'.umd.dat',34,PURPLE,'middle',weight='bold')
    text(x0+w/2,y0+122,'one format for all codes',20,BLACK,'middle')
    text(x0+w/2,y0+h+48,['check_umd.py','14 quality checks'],24,BLACK,'middle',spacing=1.3)

    # the three panels
    px, pw, gap, py, ph = 810, 330, 20, 235, 410
    arrow(735,435,px-10,435)
    # one command runs the three panels
    x1 = px+3*pw+2*gap
    out.append('<path d="M %g %g L %g %g L %g %g L %g %g" fill="none" stroke="%s" stroke-width="2.5"/>' % (px,py-22,px,py-38,x1,py-38,x1,py-22,PURPLE))
    out.append('<rect x="%g" y="%g" width="460" height="40" fill="white"/>' % ((px+x1)/2-230,py-58))
    text((px+x1)/2,py-30,'one command : umd_workflow.py',25,PURPLE,'middle',mono=True)
    panels = [('STRUCTURE',['g(r), bonds, coordination','speciation, angles','lifetimes of species','liquid / vapor']),
              ('TRANSPORT',['diffusion and its tensor','diffusion of species','viscosity','velocity autocorrelation']),
              ('THERMODYNAMICS',['averages, heat capacity','vibrational spectra','2PT entropy','free energy (integration)'])]
    boxes = []
    for i,(title,items) in enumerate(panels):
        x = px+i*(pw+gap)
        rect(x,py,pw,ph,BLACK,'white',3,dash='12 9')
        text(x+pw/2,py+50,title,28,BLACK,'middle',weight='bold')
        text(x+28,py+290,items,22,MAGENTA,spacing=1.35)
        boxes.append((x+40,py+80,pw-80,140))
    # g(r) of Si-O, Mg-O, O-O
    names, g = table(p('.gofr.dat'))
    axes(boxes[0])
    for pair,col in (('Si-O',PURPLE),('Mg-O',MAGENTA)):
        curve(g[:,0],g[:,names.index(pair)],boxes[0],xlim=(1.0,5.0),ylim=(0,10.5),color=col,sw=3)
    # MSD of O and Si, log-log
    names, m = table(p('.msd.dat'))
    axes(boxes[1])
    for el,col in (('O',PURPLE),('Si',MAGENTA)):
        curve(m[1:,0],m[1:,names.index(el)],boxes[1],xlim=(0,np.log10(m[-1,0])),ylim=(-4,1.3),color=col,sw=3,logx=True,logy=True)
    # total vDOS
    names, v = table(p('.vibr.err.dat'))
    axes(boxes[2])
    curve(v[:,0],v[:,names.index('Total_DOS')],boxes[2],xlim=(0,1400),color=PURPLE,sw=3,fill=True)
    for (bx,by,bw,bh),cap in zip(boxes,['g(r) of Si–O, Mg–O','MSD of O, Si (log–log)','vDOS']):
        text(bx+bw/2,by+bh+32,cap,18,GRAY,'middle')

    # what is new in 3.0
    new = [('any cell','NVT and NPT runs,','non-orthogonal cells'),
           ('± errors','on every statistical','quantity, calibrated'),
           ('validated','34 cases, synthetic runs','with exact answers'),
           ('6× faster','parallel workers,','memory bounded')]
    bw, bgap, by0 = (W-160-3*30)/4, 30, 700
    for i,(big,l1,l2) in enumerate(new):
        x = 80+i*(bw+bgap)
        rect(x,by0,bw,140,PURPLE,LIGHT,0,r=16)
        text(x+30,by0+58,big,40,PURPLE,weight='bold')
        text(x+30,by0+96,[l1,l2],22,BLACK,spacing=1.3)

    # footer
    # schools: the first one written, the rest of the line left empty for the next ones
    text(80,915,'Schools:',30,MAGENTA,weight='bold')
    text(215,915,'Bucharest, April 2027',30,BLACK)
    out.append('<path d="M 560 922 L %g 922" stroke="%s" stroke-width="2" stroke-dasharray="3 9" stroke-linecap="round"/>' % (W-80,GRAY))

    out.append('<path d="M 80 970 L %g 970" stroke="%s" stroke-width="2"/>' % (W-80,GRAY))
    text(80,1025,'github.com/rcaracas/MAGMATOMIX-UMD_package',26,PURPLE,mono=True)
    text(W-80,1025,'next, 3.x: towards 10⁹ atoms',28,MAGENTA,'end',weight='bold')

    svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">' % (W,H,W,H),
           '<rect width="%d" height="%d" fill="white"/>' % (W,H)] + out + ['</svg>']
    with open(svgfile,'w') as f:
        f.write('\n'.join(svg)+'\n')
    print('written',svgfile)

if __name__ == '__main__':
    if len(sys.argv) not in (2,3):
        sys.exit(__doc__)
    here = os.path.dirname(os.path.abspath(__file__))
    main(sys.argv[1],sys.argv[2] if len(sys.argv) == 3 else os.path.join(here,'umd-graphical-abstract.svg'))
