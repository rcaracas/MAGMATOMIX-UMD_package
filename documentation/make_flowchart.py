#!/usr/bin/env python3
#UMD flowchart as an SVG with editable text, in the style of Figure-umd-flowchart.ai (2021) : python3 make_flowchart.py umd-flowchart.svg
import sys, math

W, H = 1381, 1420
PURPLE, MAGENTA, BLACK, BROWN = '#782F9A', '#ED127B', '#231F20', '#6F4324'
FONT = "Helvetica, 'Helvetica Neue', Arial, sans-serif"
CMDFONT = "'Myriad Pro', 'Helvetica Neue', Helvetica, Arial, sans-serif"
out = []

def esc(s):
    return s.replace('&','&amp;').replace('<','&lt;').replace('>','&gt;')

def text(x,y,lines,size,color,anchor='start',weight='normal',font=FONT,rotate=0,spacing=1.15,cls=''):
    if isinstance(lines,str):
        lines = [lines]
    tr = ' transform="rotate(%g %g %g)"' % (rotate,x,y) if rotate else ''
    tspans = ''.join('<tspan x="%g" dy="%g">%s</tspan>' % (x,0 if i == 0 else size*spacing,esc(l)) for i,l in enumerate(lines))
    out.append('<text x="%g" y="%g" font-family="%s" font-size="%g" fill="%s" text-anchor="%s" font-weight="%s"%s>%s</text>'
               % (x,y,font,size,color,anchor,weight,tr,tspans))

def cmd(x,y,lines,anchor='start',rotate=0,size=17):     #a command of the package (purple)
    text(x,y,lines,size,PURPLE,anchor,font=CMDFONT,rotate=rotate)

def prop(x,y,lines,anchor='start',size=17.5):            #a property (magenta capitals)
    text(x,y,lines,size,MAGENTA,anchor,font=CMDFONT,spacing=1.12)

def box(x,y,w,h,label,size=21,sw=3.2):
    out.append('<rect x="%g" y="%g" width="%g" height="%g" fill="white" stroke="%s" stroke-width="%g"/>' % (x,y,w,h,BLACK,sw))
    if label:
        labels = label if isinstance(label,list) else [label]
        y0 = y + h/2 - (len(labels)-1)*size*0.55 + size*0.35
        text(x+w/2,y0,labels,size,BLACK,'middle',spacing=1.1)
    return (x,y,w,h)

def stack(x,y,w,h,labels,size=21):                          #boxes on top of each other (several files of one script)
    return [box(x,y+i*h,w,h,l,size) for i,l in enumerate(labels)]

def panel(x,y,w,h):
    out.append('<rect x="%g" y="%g" width="%g" height="%g" fill="none" stroke="%s" stroke-width="1" stroke-dasharray="11.8 11.8"/>' % (x,y,w,h,BLACK))

def line(points,head=True,sw=2):                            #black line through points, arrow head at the end
    d = 'M ' + ' L '.join('%g %g' % p for p in points)
    out.append('<path d="%s" fill="none" stroke="%s" stroke-width="%g"/>' % (d,BLACK,sw))
    if head:
        (x0,y0),(x1,y1) = points[-2],points[-1]
        a = math.atan2(y1-y0,x1-x0)
        L, Wd = 16, 6.5
        p = [(x1,y1),(x1-L*math.cos(a)+Wd*math.sin(a),y1-L*math.sin(a)-Wd*math.cos(a)),
             (x1-0.75*L*math.cos(a),y1-0.75*L*math.sin(a)),(x1-L*math.cos(a)-Wd*math.sin(a),y1-L*math.sin(a)+Wd*math.cos(a))]
        out.append('<path d="M %s Z" fill="%s"/>' % (' L '.join('%g %g' % q for q in p),BLACK))

def gray(x0,y0,x1,y1):                                      #the gray arrow from a file to a property
    a = math.atan2(y1-y0,x1-x0)
    L, Wd = 40, 17
    bx, by = x1-0.72*L*math.cos(a), y1-0.72*L*math.sin(a)
    out.append('<path d="M %g %g L %g %g" stroke="%s" stroke-width="3"/>' % (x0,y0,bx,by,BLACK))
    p = [(x1,y1),(x1-L*math.cos(a)+Wd*math.sin(a),y1-L*math.sin(a)-Wd*math.cos(a)),(bx,by),
         (x1-L*math.cos(a)-Wd*math.sin(a),y1-L*math.sin(a)+Wd*math.cos(a))]
    out.append('<path d="M %s Z" fill="url(#grayfill)" stroke="%s" stroke-width="1.5" stroke-linejoin="round"/>' % (' L '.join('%g %g' % q for q in p),BLACK))

def title(x,y,lines,anchor='start'):
    text(x,y,lines,33,BLACK,anchor,weight='bold',spacing=1.12)

#---------------------------------------------------------------- panels and titles
TOP = 345
DY = TOP-300                                                #the lower panels are drawn DY lower
panel(262,6,W-268,TOP-24)
panel(4,TOP,476,H-TOP-6)
panel(486,TOP,452,H-TOP-6)
panel(944,TOP,W-950,H-TOP-6)
title(W-22,48,'INPUT/OUTPUT','end')
title(22,H-70,['STRUCTURAL','PROPERTIES'])
title(506,H-70,['TRANSPORT','PROPERTIES'])
title(964,H-70,['THERMODYNAMIC','PROPERTIES'])
for x,flag in ((22,'-g'),(506,'-t'),(964,'-d')):           #the automated workflow runs a whole panel
    cmd(x,H-122,'umd_workflow.py '+flag+' : the whole panel',size=16)

#---------------------------------------------------------------- input / output
md = box(6,6,236,140,None,sw=4.5)
text(124,52,['ab initio /','classical MD'],32,BLACK,'middle',spacing=1.1)
text(124,128,'VASP, QBox, LAMMPS',15,BLACK,'middle')
umd = box(500,140,180,64,'.umd.',24)
line([(242,112),(498,164)])
cmd(274,40,['VaspParser.py','QBoxParser.py','LAMMPSParser2umd.py -f -i -s'])
#repair of old VASP files
line([(522,140),(522,116),(598,116),(598,138)])
cmd(512,112,'repair_umd.py -f',anchor='end',size=15)
#structures built from a snapshot, for new runs
new = box(300,212,150,58,['new','structures'],18)
line([(500,198),(452,230)])
line([(300,241),(124,241),(124,148)])
cmd(286,288,['build_umd.py -f -i -s -a','--supercell --shuffle --insert','--scale --replace --remove'],size=15)
#conversion to other formats
line([(680,172),(858,172)])
cmd(688,163,'umd2out.py -f -r -s -i',size=15.5)
for i,lab in enumerate(['.xyz.','POSCAR','LAMMPS','DeepMD']):
    box(860+i*128,140,118,64,lab,19)
#quality of the run
line([(665,140),(665,52),(798,52)])
cmd(674,40,'check_umd.py -f',size=15.5)
box(800,22,118,58,'.check.',19)
gray(920,51,962,51)
prop(968,45,['QUALITY','OF THE RUN'])

line([(515,204),(515,336+DY),(102,336+DY),(102,370+DY)])
line([(640,204),(640,400+DY)])
line([(540,204),(540,690+DY),(563,690+DY)])
line([(680,196),(1030,370+DY)])
line([(680,188),(1340,280),(1340,592+DY),(1148,592+DY)])
out.append('<g transform="translate(0,%g)">' % DY)
#---------------------------------------------------------------- structural
gofr = box(40,372,124,64,'.gofr.')
cmd(120,328,'gofr_umd.py -f -d -s -i -c -a')
gray(166,404,222,404)
prop(230,398,['INTERATOMIC','BONDING'])
line([(102,436),(102,478)])
cmd(122,463,'analyze_gofr.py -f -v')
mb = stack(30,480,152,62,[['.maximum','bonds.'],'.coordinations.'],19)
gray(184,573,240,573)
prop(248,567,['COORDINATION','NUMBERS'])
line([(102,604),(102,646)])
cmd(122,631,'bonding_umd.py -f -c -l -s -i')
bond = box(30,648,152,62,'.bonding.')
line([(102,710),(102,770)])
cmd(122,735,['speciation_and_angles.py','-f -b -c -a -r -m -s -t -l'])
spec = stack(40,772,124,58,['.step.','.popul.','.stat.'])
gray(166,858,222,858)
prop(230,852,['CHEMICAL','SPECIATION'])
specb = stack(40,958,124,58,['.angles.','.bondslives.'],19)
line([(102,946),(102,956)],head=False)
gray(166,987,222,987)
prop(230,981,'ANGLES')
gray(166,1045,222,1045)
prop(230,1039,['BOND','LIFETIMES'])
#lifetimes and liquid/vapor partitioning from the species
line([(40,859),(18,859),(18,1186),(262,1186)],head=False)
line([(18,1116),(150,1116)],head=False)
cmd(30,1108,'lifetimes.py',size=16)
gray(150,1116,200,1116)
prop(208,1110,['LIFETIMES','OF SPECIES'])
cmd(30,1160,['liquid_and_vapor_speciation.py','stat-concentrate.py'],size=16)
gray(262,1186,312,1186)
prop(320,1180,['LIQUID / VAPOR','PARTITIONING'])

#---------------------------------------------------------------- transport
cmd(652,330,['viscosity_umd.py','-f -i -n -c -v'])
visc = stack(565,402,150,60,['.stress.scf.','.visc.'],19)
line([(715,432),(800,432),(800,526)],head=False)
gray(800,524,800,562)
prop(800,590,'VISCOSITY','middle')
cmd(565,634,['msd_umd.py -f -i -s -m -a','-n -b -e'])
msd = stack(565,660,150,60,['.msd.','.diffusion.'],19)
gray(717,720,770,720)
prop(776,708,['ATOMIC','DIFFUSION'])
text(776,750,'(and tensor)',15,MAGENTA,font=CMDFONT)
cmd(585,810,'analyze_msd.py',size=16)
line([(640,780),(640,795)],head=False)
#species from speciation_and_angles.py
line([(166,830),(470,830),(470,900),(640,900),(640,960)])
cmd(650,946,['msd_clusters_umd.py -f -p'])
msdc = box(565,962,150,60,'.msd.',19)
gray(717,992,770,992)
prop(776,980,['DIFFUSION OF','CHEMICAL','SPECIES'])
cmd(565,1060,'analyze_msd_clusters.py',size=16)
line([(640,1022),(640,1042)],head=False)

#---------------------------------------------------------------- thermodynamic
cmd(1030,330,['averages.py -f -p -a'])
thermo = box(966,372,170,64,'fullthermo.',19)
gray(1138,404,1196,404)
prop(1204,384,['THERMO-','DYNAMICS','(P, T, energy,','heat capacity)'])
cmd(1160,520,['vibrations_umd.py','-f -a -t -w'])
vib = stack(966,562,180,60,['.vels.scf.','.vibr.'],19)
gray(1148,652,1206,652)
prop(1214,646,['VIBRATIONAL','SPECTRUM'])
gray(964,660,884,712)
line([(1056,682),(1056,770)])
cmd(1070,724,['vdos_thermo.py -f -t'])
vdos = stack(966,772,180,60,['.vdos.thermo.','.2pt.'],19)
gray(1148,832,1206,832)
prop(1214,802,['THERMO-','DYNAMICS','(entropy, free energy:','harmonic, 2PT,','quantum corrections)'],size=16.5)
#thermodynamic integration over runs at several lambda
for k in range(3):
    box(990-8*k,960+8*k,120,52,None,sw=2.5)
text(1050-16,992+16,'.umd.',21,BLACK,'middle')
text(1120,990,['runs at','several λ'],15,BLACK)
line([(1034,1028),(1034,1100)])
cmd(1050,1060,['thermo_integration_umd.py','-f -t -k -v'])
ti = box(958,1102,232,60,'thermo_integration.',19)
gray(1192,1132,1240,1132)
prop(1246,1110,['FREE','ENERGY,','ENTROPY'],size=16.5)

out.append('</g>')

svg = ['<?xml version="1.0" encoding="UTF-8"?>',
       '<svg xmlns="http://www.w3.org/2000/svg" width="%dpt" height="%dpt" viewBox="0 0 %d %d">' % (W,H,W,H),
       '<defs><linearGradient id="grayfill" x1="0" y1="0" x2="1" y2="1"><stop offset="0" stop-color="#B2B1B1"/><stop offset="1" stop-color="#7B7979"/></linearGradient></defs>',
       '<rect width="%d" height="%d" fill="white"/>' % (W,H)] + out + ['</svg>']
open(sys.argv[1] if len(sys.argv) > 1 else 'umd-flowchart.svg','w').write('\n'.join(svg)+'\n')
