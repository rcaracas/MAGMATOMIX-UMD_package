###
#AUTHORS: RAZVAN CARACAS
###
#The cases of the validation suite : inputs (files of data/, linked into the case folder ; zip files extracted), commands (scripts of
#src/ with their options, run one after the other in the case folder), probes (functions of the folder and of the outputs of the commands
#returning {name : value} or {name : (value, relative tolerance)}).

import os, re
import probes as P

RTOL = 1e-4             #default relative tolerance of the key numbers : the last printed digit may change (order of sums in threads)
ATOL = 1e-9


class Case:
    def __init__(self,name,inputs,commands,probes=(),full=False,about='',truth=()):
        self.name = name
        self.inputs = inputs
        self.commands = commands
        self.probes = probes
        self.truth = truth                                                  #probes returning (label, value, exact value, tolerance)
        self.full = full
        self.about = about


#----------------------------------------------------------------------------------------------------------------------------- probes

def digest(umd,prefix=''):#key numbers of a UMD file written by a case
    def probe(folder,outs):
        return {prefix+k:v for k,v in P.umd_digest(os.path.join(folder,umd)).items()}
    probe.__name__ = 'digest '+umd
    return probe


STEP = re.compile(r'^   (OK|FAILED|skipped|done before)\b.*?\s(\S+\.py\b.*)$')
HEADER = re.compile(r'^(\S+) :( |$)')

def workflow(index=0,prefix=''):#status of every step of umd_workflow.py (output of the command index), per UMD file ('folder' for the
    #steps over the folder), and, without prefix, the numbers of its summary file
    def probe(folder,outs):
        out = {}
        current = ''
        for line in outs[index].splitlines():
            if line.startswith('Steps :'):#the final count, then the failed steps again with their folder: already read
                break
            head = HEADER.match(line)
            if head:
                current = 'folder' if head.group(1).startswith('/') else head.group(1)
            found = STEP.match(line)
            if found:
                out[prefix+'step '+current+' '+found.group(2).strip()] = found.group(1)
        if prefix:
            return out
        summary = os.path.join(folder,'umd_workflow.summary.dat')
        if os.path.isfile(summary):
            for umd, row in P.table(summary,'file',sep='\t').items():
                for name, value in row.items():
                    if name != 'file':
                        out[umd+' '+name] = value
        return out
    probe.__name__ = 'workflow'
    return probe


def tab(name,first,sep='\t',columns=None,rows=None):#numbers of a table (file name or regular expression) whose header begins with first
    def probe(folder,outs):
        found = P.files(folder,name)
        if not found:
            raise ValueError('no file '+name)
        out = {}
        for f in found:
            for key, row in P.table(os.path.join(folder,f),first,sep).items():
                if rows is None or key in rows:
                    for col, value in row.items():
                        if col != first and (columns is None or col in columns) and isinstance(value,float):
                            out[f+' '+key+' '+col] = value
        return out
    probe.__name__ = 'table '+name
    return probe


def inp(name):#files 'A B value' : bond cutoffs, coordinations
    def probe(folder,outs):
        found = P.files(folder,name)
        if not found:
            raise ValueError('no file '+name)
        return {f+' '+k:v for f in found for k,v in P.pairs(os.path.join(folder,f)).items()}
    probe.__name__ = 'pairs '+name
    return probe


def checks(name):#result of every check of check_umd.py : OK, WARNING, skipped
    def probe(folder,outs):
        out = {}
        for f in P.files(folder,name):
            title = None
            with open(os.path.join(folder,f)) as g:
                for line in g:
                    found = re.match(r'(\d+)\. (.*)',line)
                    if found:
                        title = found.group(1)
                    elif line.strip().startswith('->') and title:
                        out[f+' check '+title] = line.split()[1]
        return out
    probe.__name__ = 'checks '+name
    return probe


def count(name):#number of files whose names match a regular expression
    def probe(folder,outs):
        return {'files '+name:len(P.files(folder,name))}
    probe.__name__ = 'count '+name
    return probe


def stdout(index,label,pattern,tolerance=None):#numbers printed by the command index (a regular expression ; one group per number)
    def probe(folder,outs):
        found = re.search(pattern,outs[index])
        if not found:
            raise ValueError('not found in the output : '+pattern)
        out = {}
        for k, g in enumerate(found.groups()):
            v = P.number(g)
            out[label+('' if len(found.groups()) == 1 else ' %d' % k)] = v if tolerance is None else (v,tolerance)
        return out
    probe.__name__ = 'stdout '+label
    return probe


def digests(name):#the UMD files matching a regular expression
    def probe(folder,outs):
        return {f+' '+k:v for f in P.files(folder,name) for k,v in P.umd_digest(os.path.join(folder,f)).items()}
    probe.__name__ = 'digests '+name
    return probe


#----------------------------------------------------------------------------------------------------------------------------- cases

CASES = []

#Parsers : raw outputs of the simulation codes -> UMD
CP2K = ['cp2k/position.xyz','cp2k/velocity.xyz','cp2k/force.xyz','cp2k/log.out','cp2k/stress.out','cp2k/main_2.out']
CP2K_ARGS = ['-f','position.xyz','-v','velocity.xyz','-r','force.xyz','-e','log.out','-a','stress.out','-c','main_2.out']  #the cell from the MD_INI| lines of
#the main output of a later segment of the same NVT run (steps 75266-76266 ; the trajectory files are steps 9001-10001) : 14.6 A cubic
CASES += [
    Case('parse_vasp_zns',['zns.T1000K.prod.outcar'],[('VaspParser.py',['-f','zns.T1000K.prod.outcar'])],
         [digest('zns.T1000K.prod.outcar.umd.dat')],about='VASP OUTCAR, ZnS 1000 K, on-the-fly ML (10000 steps, 57 DFT)'),
    Case('parse_abinit',['abinit.exMD.output'],[('AbinitParser.py',['-f','abinit.exMD.output'])],
         [digest('abinit.exMD.output.umd.dat')],about='ABINIT 6.4.1 Langevin MD, Mg7Si8O24Fe, non-orthogonal cell'),
    Case('parse_abinit_is',['abinit.exMD.output'],[('AbinitParser.py',['-f','abinit.exMD.output','-i','10','-s','3'])],
         [digest('abinit.exMD.output.umd.dat')],about='ABINIT parser with -i 10 -s 3'),
    Case('parse_lammps_water',['lamps-test1/water.dump','lamps-test1/log.lammps'],
         [('LAMMPSParser2umd.py',['-f','water.dump','-l','log.lammps'])],
         [digest('water.dump.umd.dat')],about='LAMMPS DeePMD water + NaCl, NPT 500 K 10 GPa, no forces (ice VII)'),
    Case('parse_lammps_dump',['lammps-test-dump/test.dump','lammps-test-dump/log.lammps'],
         [('LAMMPSParser2umd.py',['-f','test.dump','-l','log.lammps'])],
         [digest('test.dump.umd.dat')],about='LAMMPS DeePMD water 600 K, velocities, forces, xu'),
    Case('parse_lammps_melt',['LAMMPS-dump-outputs.zip'],
         [('LAMMPSParser2umd.py',['-f','dump.lammpstrj','-l','log.lammps','-a','Pxy_vs_time.dat'])],
         [digest('dump.lammpstrj.umd.dat')],about='LAMMPS MACE melt NPT, dump without id, stress file'),
    Case('parse_lammps_frame',['Lammps_output.zip'],
         [('LAMMPSParser2umd.py',['-f','test.lammpstrj','-a','Pxy_vs_time.dat','-t','1'])],
         [digest('test.lammpstrj.umd.dat')],about='LAMMPS pyrolite melt, one frame without id, velocities or log'),
    Case('parse_cp2k',CP2K,[('CP2KParser.py',CP2K_ARGS)],
         [digest('position.xyz.umd.dat')],about='CP2K NVT, 100 H2O + 2 H2SO4, positions, velocities, forces, energies, stress'),
    Case('parse_cp2k_is',CP2K,[('CP2KParser.py',CP2K_ARGS+['-i','10','-s','3'])],
         [digest('position.xyz.umd.dat')],about='CP2K parser with -i 10 -s 3'),
]

#The automated workflow (check, structure, transport, thermodynamics) on each reference UMD file
for name, umd, full, about in (
        ('workflow_glycine','glycine.umd.dat',False,'glycine crystal, cubic, 100 atoms, 10000 snapshots'),
        ('workflow_pyrolite','pyrolite.a13.T3000K.outcar.umd.dat',False,'pyrolite melt 3000 K'),
        ('workflow_zns','zns.T1000K.prod.outcar.umd.dat',False,'ZnS 1000 K, on-the-fly ML'),
        ('workflow_corundum','al2o3.new.umd.dat',False,'hexagonal corundum + 1 He'),
        ('workflow_andesite','hydrousandesite.a15.0.T3000K.reduced.-A-.outcar.umd.dat',True,'hydrous andesite melt 3000 K, 468 MB, restarts')):
    base = umd[:-8]
    CASES.append(Case(name,[umd],[('umd_workflow.py',['-f',umd])],
                      [workflow(),checks(re.escape(base)+r'\.check\.dat'),inp(re.escape(base)+r'\.maximum_bonds\.inp'),
                       inp(re.escape(base)+r'\.coordinations\.inp'),tab(re.escape(base)+r'\.vacf\.diffusion\.dat','column'),
                       tab(re.escape(base)+r'\.2pt\.dat','element'),tab(r'diffusivities_clusters\.txt','file',columns=['D(m2/s)','alpha'])],
                      full=full,about='umd_workflow.py, '+about))

PYRO = 'pyrolite.a13.T3000K.outcar.umd.dat'
GLY = 'glycine.umd.dat'
COR = 'al2o3.new.umd.dat'

#Parsed files analysed : runs of other codes than VASP through the workflow
CASES += [
    Case('chain_abinit',['abinit.exMD.output'],[('AbinitParser.py',['-f','abinit.exMD.output']),
                                                ('umd_workflow.py',['-f','abinit.exMD.output.umd.dat'])],
         [workflow(1),checks(r'.*\.check\.dat'),inp(r'.*\.coordinations\.inp')],about='ABINIT run parsed, then the workflow'),
    Case('chain_lammps',['lammps-test-dump/test.dump','lammps-test-dump/log.lammps'],
         [('LAMMPSParser2umd.py',['-f','test.dump','-l','log.lammps']),('umd_workflow.py',['-f','test.dump.umd.dat'])],
         [workflow(1),checks(r'.*\.check\.dat'),inp(r'.*\.coordinations\.inp')],about='LAMMPS water parsed, then the workflow'),
    Case('chain_cp2k',CP2K,[('CP2KParser.py',CP2K_ARGS),('umd_workflow.py',['-f','position.xyz.umd.dat'])],
         [workflow(1),checks(r'.*\.check\.dat'),inp(r'.*\.coordinations\.inp')],about='CP2K run parsed, then the workflow'),
]

#Options not used by the workflow
CASES += [
    Case('gofr_select',[PYRO],[('gofr_umd.py',['-f',PYRO,'-c','Si,Mg','-a','O']),
                               ('analyze_gofr.py',['-f','pyrolite.a13.T3000K.outcar.SiMg-O.gofr.dat']),
                               ('bonding_umd.py',['-f',PYRO,'-c','pyrolite.a13.T3000K.outcar.SiMg-O.maximum_bonds.inp','-i','100','-s','2']),
                               ('gofr_umd.py',['-f',PYRO]),
                               ('analyze_gofr.py',['-f','pyrolite.a13.T3000K.outcar.gofr.dat','-c','O','-g','0.05'])],
         [inp(r'.*\.maximum_bonds\.inp'),inp(r'.*\.coordinations\.inp')],
         about='gofr -c/-a, analyze_gofr on it, bonding with -i -s ; all the pairs analyzed with -c O -g'),
    Case('msd_options',[PYRO,GLY],[('msd_umd.py',['-f',PYRO,'-a','tensor','-p','[1,0,0,0,1,0]']),
                                   ('msd_umd.py',['-f',GLY,'-m','atoms','-a','cryst','-i','1000','-s','2','-n','4']),
                                   ('analyze_msd.py',[]),('analyze_msd.py',['-s','100','-e','400'])],
         [tab(r'.*\.diffusion\.dat','column',columns=['D_(m2/s)','D_error_(m2/s)','alpha']),tab(r'diffusivities\.txt','file')],
         about='msd tensor and a plane, per atom along a, b, c with -i -s -n, analyze_msd automatic and with a window'),
    Case('vibrations_options',[COR],[('vibrations_umd.py',['-f',COR,'-a','cryst','-p','none','-g','10','-n','3','-i','200']),
                                     ('vdos_thermo.py',['-f',COR,'-t','300,1500'])],
         [tab(r'.*\.vacf\.diffusion\.dat','column'),tab(r'.*\.vdos\.thermo\.dat','element',columns=['ZPE','E','S','Cv','dE','dS'])],
         about='vibrations along a, b, c of a hexagonal cell, no taper, smoothing ; vdos_thermo at two temperatures'),
    Case('speciation_options',[GLY],
         [('bonding_umd.py',['-f',GLY,'-l','1.6','-i','1000']),
          ('speciation_and_angles.py',['-f',GLY,'-b','glycine.bonding.dat','-c','all','-a','all','-r','0']),
          ('speciation_and_angles.py',['-f',GLY,'-b','glycine.bonding.dat','-c','C,N','-a','H,O','-r','1','-t','1','-l','1','-s','2','-m','5']),
          ('speciation_and_angles.py',['-f',GLY,'-b','glycine.bonding.dat','-c','C','-a','all','-r','2']),
          ('liquid_and_vapor_speciation.py',['-s','glycine.r0.step.dat','-p','glycine.r1.popul.dat','-l','12']),
          ('lifetimes.py',['-f','glycine.r1.popul.dat','-l','12']),
          ('stat-concentrate.py',['-l','12','-c','C,H']),
          ('msd_clusters_umd.py',['-f',GLY,'-p','glycine.r1.popul.dat']),
          ('analyze_msd_clusters.py',[])],
         [tab(r'glycine\.r[012]\.stat\.dat','Cluster',columns=['Percent']),tab(r'glycine\.r1\.(liquid|vapor)\.stat\.dat','Cluster',columns=['Percent']),
          tab(r'diffusivities_clusters\.txt','file',columns=['D(m2/s)','alpha']),count(r'.*\.msd\.dat'),count(r'.*\.pdf')],
         about='bonds of one length from snapshot 1000 ; r0, r1 with -s -m, r2 ; liquid/vapor, lifetimes -l, stat-concentrate, cluster MSDs'),
    Case('averages_options',[PYRO,GLY],[('averages.py',['-f',PYRO+','+GLY]),('averages.py',['-f',PYRO,'-p','StressTensor','-l']),
                                        ('averages.py',['-f',PYRO,'-a','-u','1','-i','100'])],
         [stdout(0,'pyrolite P T E',r'pyrolite\S*\s+#PstdTstdEstd\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)'),
          tab(r'fullthermo\.dat','file')],
         about='averages : default line, a field with -l, all the thermodynamics with units -u 1 and -i'),
    Case('check_options',[GLY,PYRO],[('check_umd.py',['-f',GLY,'-q','C:1','-i','100','-s','5']),('check_umd.py',['-f',PYRO,'-l','5'])],
         [checks(r'.*\.check\.dat')],about='check_umd with charges -q, -i, -s ; overlap level -l'),
    Case('umd2out_formats',[GLY,COR],[('umd2out.py',['-f',GLY,'-s','1000','-r','xyz,poscar,lammps,deepmd,lammpstrj','-u']),
                                      ('umd2out.py',['-f',COR,'-k','10','-r','1,2,3','-o','corundum'])],
         [count(r'.*POSCAR.*'),count(r'.*\.lammps\.data'),count(r'.*\.xyz'),count(r'.*\.lammpstrj')],
         about='umd2out : all the formats, unwrapped, every 1000 ; one snapshot of the hexagonal cell (LAMMPS tilts)'),
    Case('build_structures',[GLY,COR,'molecules.dat','test.poscar.vasp'],
         [('gofr_umd.py',['-f',GLY]),('analyze_gofr.py',['-f','glycine.gofr.dat']),
          ('build_umd.py',['-f',GLY,'--supercell','2,1,1','--name','sc']),
          ('build_umd.py',['-f',GLY,'-i','0','--shuffle','2','--seed','1','--name','sh']),
          ('build_umd.py',['-f',GLY,'-i','0','--scale','1.5','--insert','molecules.dat','-r','1.5','--seed','2','--name','ins']),
          ('build_umd.py',['-f',GLY,'-i','0','--scale','1.05','--density','1.2','--name','sd']),
          ('build_umd.py',['-f',GLY,'-i','0','--vacuum','x:10,y:5','--name','vac']),
          ('build_umd.py',['-f',GLY,'-i','0','--replace','O:S:2','--replace','N:P:0.5','--seed','3','--name','rep']),
          ('build_umd.py',['-f',GLY,'-i','0','--remove','C2H5NO2:2','-b','glycine.maximum_bonds.inp','--seed','4','--name','rem']),
          ('build_umd.py',['-f',GLY,'--extract','0-9,N','-s','100','--name','ext']),
          ('build_umd.py',['-f','test.poscar.vasp','--supercell','2,2,2','-t','1000','--seed','5','--name','pos']),
          ('build_umd.py',['-a','10','--insert','molecules.dat','-t','300','--seed','6','--name','box','-o','umd,poscar,xyz']),
          ('build_umd.py',['-f',COR,'-i','0','--vacuum','y:5','--supercell','1,1,2','-n','2','--seed','7','--name','hex'])],
         [digests(r'.*\.umd\.dat')],
         about='build_umd : every operation, from a UMD file, a POSCAR and an empty box ; hexagonal cell ; extraction'),
    Case('workflow_small',['test.8atoms10stepa.umd.dat','test.5atoms2step.umd.dat','test.glycine-10steps.umd.dat',
                           'test.glycine-1step.umd.dat','test.8atoms1step.umd.dat','short.at1.umd.dat','short.at2.umd.dat',
                           'short.at6.umd.dat','glycine.2steps.umd.dat','test.umd.dat','outcar.umd.dat'],
         [('umd_workflow.py',[]),('umd_clean.py',[])],
         [workflow(),stdout(1,'umd_clean outputs',r'(\d+) output')],
         about='the workflow on small and odd files (1 snapshot, 2 steps, 5 atoms...) ; umd_clean listing'),
    Case('workflow_rerun',[PYRO],[('umd_workflow.py',['-f',PYRO,'-i','500','-s','2','-c','Si,Al','-a','O']),
                                  ('umd_workflow.py',['-f',PYRO,'-i','500','-s','2','-c','Si,Al','-a','O']),
                                  ('umd_workflow.py',['-f',PYRO,'-i','500','-s','2','-g'])],
         [workflow(0,'run 1 '),workflow(1,'run 2 '),workflow(2,'run 3 '),workflow(2)],
         about='the workflow with -i -s -c -a, run again (done before), then -g with other elements (reruns)'),
    Case('repair_corundum',['al2o3.umd.dat'],[('check_umd.py',['-f','al2o3.umd.dat']),('repair_umd.py',['-f','al2o3.umd.dat']),
                                              ('check_umd.py',['-f','al2o3.repaired.umd.dat'])],
         [checks(r'.*\.check\.dat'),digest('al2o3.repaired.umd.dat')],full=True,
         about='the old corundum file (missed crossings, zero velocities) checked, repaired, checked'),
]


#----------------------------------------------------------------------------------------------------------------- synthetic files
#Files written by synthetic.py with exact answers : each result is compared with its exact value (truth probes). Statistical results
#are allowed 8 of their errors when the errors come from 5 blocks (Student's t with 4 degrees of freedom : 99.9 %), 5 with 4 blocks of
#viscosity_umd.py read at 10 correlation times ; exact results (counts, angles, ideal gas) a small tolerance.

import json

def truth_of(folder,kind):
    with open(os.path.join(folder,kind+'.truth.json')) as f:
        return json.load(f)


def t_idealgas(folder,outs):
    t = truth_of(folder,'idealgas')
    path = os.path.join(folder,'idealgas.gofr.dat')
    lines = P.rows(path,'\t')
    names = lines[0]
    data = [[float(x) for x in r] for r in lines[1:]]
    out = []
    for k, name in enumerate(names[1:],1):
        if name.startswith('Int('):
            a, b = name[4:-1].split('-')
            rlast = data[-1][0]+0.5*(data[-1][0]-data[-2][0])                  #outer edge of the last bin
            exact = (t['counts'][b]-(a == b))/t['counts'][b]*t['rho'][b]*4/3*math.pi*rlast**3
            out.append(('idealgas '+name+' at %.2f A / exact' % rlast,data[-1][k]/exact,1.0,0.02))
        else:
            g = [row[k] for row in data if row[0] > 2.0]
            out.append(('idealgas mean g '+name+' beyond 2 A',sum(g)/len(g),1.0,0.02))
    return out


def t_fcc(folder,outs):
    t = truth_of(folder,'fcc')
    out = []
    for cell in ('fcc_cubic','fcc_prim'):
        cn = P.pairs(os.path.join(folder,cell+'.coordinations.inp'))['Cu-Cu']
        out.append((cell+' coordination Cu-Cu',cn,12.0,0.01))
        rmax = P.pairs(os.path.join(folder,cell+'.average_bonds.inp'))['Cu-Cu']
        out.append((cell+' first maximum (A)',rmax,t['nn_distance'],0.02))
        counts = []
        with open(os.path.join(folder,cell+'.bonding.dat')) as f:
            for line in f:
                if re.match(r'\d+\t',line):
                    counts.append(len(line.split())-1)
        out.append((cell+' atoms with 12 bonds (fraction)',sum(c == 12 for c in counts)/len(counts),1.0,0.0))
        stat = P.table(os.path.join(folder,cell+'.r0.stat.dat'),'Cluster',sep='\t')
        out.append((cell+' r0 : one cluster of all the atoms (fraction)',max(r['Percent'] for r in stat.values()),1.0,1e-9))
    return out


def t_water(folder,outs):
    t = truth_of(folder,'water')
    out = []
    for r in ('r1','r0'):
        stat = P.table(os.path.join(folder,'water.'+r+'.stat.dat'),'Cluster',sep='\t')
        out.append(('water '+r+' fraction of O_1H_2',stat.get('O_1H_2',{}).get('Percent',0.0),1.0,1e-9))
    stat = P.table(os.path.join(folder,'water.r1.stat.dat'),'Cluster',sep='\t')
    out.append(('water mean angle H-O-H in the r1 species',stat['O_1H_2']['Mean Angles'] if 'Mean Angles' in stat['O_1H_2'] else
                list(stat['O_1H_2'].values())[-1],t['angle'],1e-3))
    with open(os.path.join(folder,'water.angles.dat')) as f:
        found = re.search(r'Mean of all angles : ([0-9.eE+-]+)',f.read())
    out.append(('water mean of all angles',float(found.group(1)),t['angle'],1e-3))
    return out


def diffusion(folder,name):#{column : (D, error)} of a diffusion.dat or vacf.diffusion.dat file
    return {k:(r['D_(m2/s)'],r['D_error_(m2/s)']) for k,r in P.table(os.path.join(folder,name),'column',sep='\t').items()}


def t_randomwalk(folder,outs):
    t = truth_of(folder,'randomwalk')
    d = diffusion(folder,'randomwalk.diffusion.dat')
    out = []
    exact = {'Ar':t['D']['Ar'],'Kr':t['D']['Kr'],'Xe':t['Xe_isotropic']}
    labels = {'xx':(0,0),'yy':(1,1),'zz':(2,2),'xy':(0,1),'yz':(1,2),'zx':(2,0)}
    for lab,(i,j) in labels.items():
        exact['Xe_'+lab] = t['Xe_tensor'][i][j]
        for el in ('Ar','Kr'):
            exact[el+'_'+lab] = t['D'][el] if i == j else 0.0
    for col, value in exact.items():
        D, err = d[col]
        out.append(('randomwalk D '+col+' (m2/s)',D,value,8*err))
    #principal values and directions of Xe (diffusion.tensor.dat)
    text = open(os.path.join(folder,'randomwalk.diffusion.tensor.dat')).read().split('element Xe')[1]
    for k in range(3):
        found = re.search(r'\nD%d\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)' % (k+1),text)
        D, err = float(found.group(1)), float(found.group(2))
        out.append(('randomwalk Xe principal D%d (m2/s)' % (k+1),D,t['Xe_principal'][k],8*err))
        v = [float(found.group(g)) for g in (3,4,5)]
        cos = abs(sum(a*b for a,b in zip(v,t['Xe_directions'][k])))
        out.append(('randomwalk Xe principal direction %d |cos| with the exact one' % (k+1),cos,1.0,0.1))
    return out


def t_ou(folder,outs):
    t = truth_of(folder,'ou')
    out = []
    for name, label in (('ou.vacf.diffusion.dat','VACF'),('ou.diffusion.dat','MSD')):
        d = diffusion(folder,name)
        for el in ('Ar','Ne'):
            out.append(('ou D '+el+' from the '+label+' (m2/s)',d[el][0],t['D'][el],8*d[el][1]))
    out.append(('ou kinetic temperature (K)',P.find(outs[1],r'Kinetic temperature of the velocities \(3N degrees of freedom\) :\s+(\S+)'),t['T'],10.0))
    for el, n in zip(('Ar','Ne'),(72,72)):
        out.append(('ou integral of the vDOS of '+el,P.find(outs[1],el+r' (\S+) / %d' % n),float(n),0.5))
    return out


def t_oscillators(folder,outs):
    t = truth_of(folder,'oscillators')
    out = []
    lines = P.rows(os.path.join(folder,'oscillators.vibr.dat'),'\t')
    names = lines[0]
    data = [[float(x) for x in r] for r in lines[1:]]
    for el in ('Si','O'):
        k = names.index('VDos_'+el)
        peak = max(data,key=lambda r: r[k])[0]
        out.append(('oscillators vDOS peak of '+el+' (cm-1)',peak,t['nu'][el],3.0))
    rows = {}
    for r in P.rows(os.path.join(folder,'oscillators.vdos.thermo.dat'),'\t'):
        if r[0] == 'T(K)':
            names = r
        else:
            rows[r[1]] = dict(zip(names,r))
    for el in ('Si','O'):
        for q, rel in (('ZPE',0.005),('E',0.005),('F',0.005),('S',0.01),('Cv',0.01)):
            v = float(rows[el][q])
            out.append(('oscillators Einstein '+q+' of '+el,v,t[el][q],rel*abs(t[el][q])))
    return out


def t_stress(folder,outs):
    t = truth_of(folder,'stress')
    out = []
    for k, label in ((1,'all components'),(2,'off-diagonal')):
        eta, err = P.find(outs[k],r'Viscosity =\s+(\S+)\s+\+/-\s+(\S+)\s+Pa\.s'), P.find(outs[k],r'Viscosity =\s+\S+\s+\+/-\s+(\S+)\s+Pa\.s')
        out.append(('stress viscosity, '+label+' (Pa.s)',eta,t['eta'],5*err))
    for name, key in (('Temperature','T'),('Pressure','P')):
        mean = P.find(outs[3],r'Averages of\s+'+name+r'\s+over\s+\d+\s+steps: mean =\s+(\S+)')
        err = P.find(outs[3],r'Averages of\s+'+name+r'\s+over.*?error of the mean =\s+(\S+)')
        out.append(('stress mean '+name,mean,t[key+'_mean'],4*t[key+'_err']))
        out.append(('stress error of the mean '+name+' / exact',err/t[key+'_err'],1.0,0.3))
    return out


def t_ti(folder,outs):
    t = truth_of(folder,'ti')
    out = []
    text = outs[1]
    for line in text.splitlines():
        found = re.match(r'ti_(\S+)\.umd\.dat\s+\S+\s+\S+\s+(\S+)\s+(\S+)',line)
        if found and found.group(1) in t['U']:
            out.append(('ti <U> at lambda '+found.group(1)+' (eV)',float(found.group(2)),t['U'][found.group(1)],5*float(found.group(3))))
    Fig = P.find(text,r'ideal gas free energy at T = \S+ K, V = \S+ A\^3 : (\S+) eV')
    out.append(('ti ideal gas free energy (eV)',Fig,t['F_ideal'],1e-5))
    F, err = P.find(text,r'\nF = (\S+) \+- (\S+) eV per cell',1), P.find(text,r'\nF = (\S+) \+- (\S+) eV per cell',2)
    out.append(('ti free energy F (eV per cell)',F,t['F'],4*err+0.02))
    S, err = P.find(text,r'= (\S+) \+- (\S+) kB/atom',1), P.find(text,r'= (\S+) \+- (\S+) kB/atom',2)
    out.append(('ti entropy S (kB/atom)',S,t['S_kB_per_atom'],4*err+0.02/(KB_T(t))))
    return out


def KB_T(t):#kB T per atom in eV, to turn 0.02 eV per cell into kB/atom
    return 8.617333262e-5*t['T']*t['natom']


def t_faults(folder,outs):
    out = []
    for name, expect in (('faults_clean',{'1':'OK','2':'OK','3':'OK'}),('faults_bad',{'1':'WARNING','3':'WARNING'})):
        status = {}
        title = None
        with open(os.path.join(folder,name+'.check.dat')) as f:
            for line in f:
                found = re.match(r'(\d+)\. ',line)
                if found:
                    title = found.group(1)
                elif line.strip().startswith('->') and title:
                    status[title] = line.split()[1]
        for k, v in expect.items():
            out.append((name+' check '+k,status.get(k,'missing'),v,0))
    return out


import math
CASES += [
    Case('syn_idealgas',[],[('@synthetic.py',['idealgas']),('gofr_umd.py',['-f','idealgas.umd.dat'])],truth=[t_idealgas],
         about='synthetic ideal gas in a triclinic cell deformed at every snapshot : g(r) = 1, running coordination exact'),
    Case('syn_fcc',[],[('@synthetic.py',['fcc'])]+[(s,a) for cell in ('fcc_cubic','fcc_prim') for s,a in
                     (('gofr_umd.py',['-f',cell+'.umd.dat']),('analyze_gofr.py',['-f',cell+'.gofr.dat']),
                      ('bonding_umd.py',['-f',cell+'.umd.dat','-l','3.0856']),
                      ('speciation_and_angles.py',['-f',cell+'.umd.dat','-b',cell+'.bonding.dat','-c','all','-a','all','-r','0']))],
         truth=[t_fcc],about='synthetic fcc Cu in its cubic and its 60 deg primitive cell : 12 neighbours, distances, bonds'),
    Case('syn_water',[],[('@synthetic.py',['water']),('bonding_umd.py',['-f','water.umd.dat','-l','1.2']),
                         ('speciation_and_angles.py',['-f','water.umd.dat','-b','water.bonding.dat','-c','O','-a','H','-r','1','-t','1']),
                         ('speciation_and_angles.py',['-f','water.umd.dat','-b','water.bonding.dat','-c','all','-a','all','-r','0'])],
         truth=[t_water],about='synthetic rigid H2O across the faces of a triclinic cell : species and angle exact'),
    Case('syn_msd',[],[('@synthetic.py',['randomwalk']),('msd_umd.py',['-f','randomwalk.umd.dat','-a','tensor'])],truth=[t_randomwalk],
         about='synthetic random walks, isotropic and anisotropic : D, D tensor, principal values and directions'),
    Case('syn_vacf',[],[('@synthetic.py',['ou']),('vibrations_umd.py',['-f','ou.umd.dat']),('msd_umd.py',['-f','ou.umd.dat'])],
         truth=[t_ou],
         about='synthetic Ornstein-Uhlenbeck velocities : D from the VACF and the MSD, kinetic T, vDOS normalization'),
    Case('syn_vdos',[],[('@synthetic.py',['oscillators']),('vibrations_umd.py',['-f','oscillators.umd.dat']),
                        ('vdos_thermo.py',['-f','oscillators.umd.dat'])],truth=[t_oscillators],
         about='synthetic harmonic oscillators : vDOS peaks, Einstein ZPE, E, F, S, Cv'),
    Case('syn_viscosity',[],[('@synthetic.py',['stress']),('viscosity_umd.py',['-f','stress.umd.dat','-e','250']),
                             ('viscosity_umd.py',['-f','stress.umd.dat','-e','250','-c','offdiag']),
                             ('averages.py',['-f','stress.umd.dat','-l'])],truth=[t_stress],
         about='synthetic stress (Ornstein-Uhlenbeck) : viscosity ; T and P as AR(1) series : averages and their errors'),
    Case('syn_ti',[],[('@synthetic.py',['ti']),('thermo_integration_umd.py',['-t','1500'])],truth=[t_ti],
         about='synthetic harmonic wells at the Gauss-Legendre lambdas : <U>, ideal gas, F, S exact'),
    Case('syn_check',[],[('@synthetic.py',['faults']),('check_umd.py',['-f','faults_clean.umd.dat']),
                         ('check_umd.py',['-f','faults_bad.umd.dat'])],truth=[t_faults],
         about='synthetic random walk, clean and with a duplicated snapshot and a missed crossing : check_umd'),
]
