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
                current = 'folder' if os.path.isabs(head.group(1)) else head.group(1)      #the folder : / or, on Windows, D:\\
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
                       tab(re.escape(base)+r'\.2pt\.dat','element')],
                      full=full,about='umd_workflow.py, '+about))

def comp_numbers(path):#numbers of the summary of composition_umd.py : ratio, error, fraction in voids of each element ; correlations ; voids ; 1/e times
    out = {}
    section = None
    with open(path) as f:
        lines = [l.rstrip('\n') for l in f]
    labels = []
    for line in lines:
        if not line or line.startswith('#'):
            section = None if not line else section
            continue
        cells = line.split('\t')
        if line.startswith('element\tevaporation_(1/ps)\t'):
            section = 'fluxes'
            continue
        if line.startswith('element\tevaporations\t'):
            section = 'exchange'
            continue
        if line.startswith('element\tx\t'):
            section = 'conc'
            continue
        if line.startswith('element\t'):
            section = 'elements'
            continue
        if line.startswith('state\t'):
            section = 'states'
            continue
        if line.startswith('Correlation'):
            section = 'corr'
            continue
        if line.startswith('errors\t'):
            section = 'correrr'
            continue
        if line.startswith('Voids'):
            section = 'voids'
            continue
        if line.startswith('kind\tregions\t'):
            section = 'regions'
            continue
        if line.startswith('series\t'):
            section = 'times'
            continue
        if line.startswith('Times'):
            section = 'timeshead'
            continue
        if section in ('corr','correrr') and cells[0] == '':
            labels = cells[1:]
            continue
        num = lambda x: float(x.lstrip('>')) if x.strip() not in ('','nan') else float('nan')
        if section == 'elements':
            out['ratio '+cells[0]] = num(cells[5])
            out['ratio_err '+cells[0]] = num(cells[6])
            if len(cells) > 7 and cells[7]:
                out['in_voids '+cells[0]] = num(cells[7])
            if len(cells) > 8 and cells[8]:
                out['in_largest_dense '+cells[0]] = num(cells[8])
        elif section == 'conc' and len(cells) >= 10:
            out['conc_ratio '+cells[0]] = num(cells[2])
            out['conc_ratio_err '+cells[0]] = num(cells[3])
            for k in range(4):
                out['conc_quarter%d %s' % (k+1,cells[0])] = num(cells[4+k])
            out['conc_quarter_errs '+cells[0]] = [float(v.lstrip('>')) for v in cells[8].split()]
            out['conc_change '+cells[0]] = num(cells[9])
        elif section == 'states' and len(cells) >= 8:
            out['periods '+cells[0]] = num(cells[1])
            out['cut_periods '+cells[0]] = num(cells[2])
            out['mean_period '+cells[0]] = num(cells[3])
            out['mean_period_err '+cells[0]] = num(cells[4])
            out['time_fraction '+cells[0]] = num(cells[7])
            if len(cells) > 8 and cells[8]:
                out['random_alloy_fraction '+cells[0]] = num(cells[8])
        elif section == 'fluxes' and len(cells) >= 6:
            for k, key in enumerate(('flux_evaporation','flux_condensation','flux_evaporation_area','flux_condensation_area'),1):
                out[key+' '+cells[0]] = num(cells[k])
            out['equilibrium_from '+cells[0]] = num(cells[5]) if cells[5] != 'not reached' else float('inf')
        elif section == 'exchange' and len(cells) >= 13:
            for k, key in enumerate(('evaporations','condensations','k_evaporation','k_evaporation_err','k_condensation','k_condensation_err',
                                     'stay_condensed','stay_condensed_err','stay_gas','stay_gas_err','gas_fraction','net'),1):
                out[key+' '+cells[0]] = num(cells[k])
        elif section in ('corr','correrr'):
            for l, x in zip(labels,cells[1:]):
                out[('corr ' if section == 'corr' else 'corr_err ')+cells[0]+'-'+l] = num(x)
        elif section == 'regions' and len(cells) >= 13:
            for k, key in enumerate(('regions','short_lived','complete','mean_lifetime','mean_lifetime_err','median_lifetime','longest','new','split',
                                     'vanished','merged','mean_volume'),1):
                out['regions_'+key+' '+cells[0]] = num(cells[k])
        elif section == 'times' and len(cells) >= 8:
            for k, key in enumerate(('','_err','_persistent','_persistent_err','_fluctuating','_fluctuating_err','_end'),1):
                out[cells[0]+key] = num(cells[k])
        elif section == 'voids' and len(cells) >= 2:
            out[cells[0]] = num(cells[1])
            if len(cells) > 2:
                out[cells[0]+'_err'] = num(cells[2])
    return out


def comp_summary(name):#numbers of a summary of composition_umd.py, for the cases on real runs
    def probe(folder,outs):
        return {name+' '+k:v for k, v in comp_numbers(os.path.join(folder,name)).items() if isinstance(v,float)}
    probe.__name__ = 'composition '+name
    return probe


def exch_numbers(path):#numbers of the summary of exchange_umd.py : interface, time in each state, paths, exchange, fluxes
    out, section, head = {}, None, []
    num = lambda x: float(x.lstrip('>')) if x.strip() not in ('','nan','not reached') else (float('inf') if x.strip() == 'not reached' else float('nan'))
    for line in open(path):
        line = line.rstrip('\n')
        if not line or line.startswith('#'):
            section = None if not line else section
            continue
        cells = line.split('\t')
        if line.startswith('Interface (means'):
            section = 'interface'
            continue
        if line.startswith('element\t'):
            section, head = 'table', cells
            continue
        if section == 'interface' and len(cells) >= 2:
            out[cells[0]] = num(cells[1])
        elif section == 'table':
            for h, v in zip(head[1:],cells[1:]):
                out[h+' '+cells[0]] = num(v)
    return out


def exch_summary(name):#numbers of a summary of exchange_umd.py, for the cases on real runs
    def probe(folder,outs):
        return {name+' '+k:v for k, v in exch_numbers(os.path.join(folder,name)).items()}
    probe.__name__ = 'exchange '+name
    return probe


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
    Case('rotations_ice',['lamps-test1/water.dump','lamps-test1/log.lammps'],
         [('LAMMPSParser2umd.py',['-f','water.dump','-l','log.lammps']),('gofr_umd.py',['-f','water.dump.umd.dat']),
          ('analyze_gofr.py',['-f','water.dump.gofr.dat']),('bonding_umd.py',['-f','water.dump.umd.dat','-c','water.dump.maximum_bonds.inp']),
          ('rotations_umd.py',['-f','water.dump.umd.dat','-c','O','-a','H'])],
         [tab(r'water\.dump\.OH2\.rot\.dat','quantity')],about='rotations_umd.py on ice VII (LAMMPS DeePMD water + NaCl, 500 K, 10 GPa) : H2O dipoles'),
    Case('composition_h2',['pyrolite+H2.test.umd.dat'],[('composition_umd.py',['-f','pyrolite+H2.test.umd.dat']),
                                                         ('composition_umd.py',['-f','pyrolite+H2.test.umd.dat','-l','7.5','-i','500','-s','2','-c','1:1:1,3:0:2','-a','Si,O,H'])],
         [comp_summary('pyrolite+H2.test.comp.6x6x6/summary.dat'),comp_summary('pyrolite+H2.test.comp.4x4x4/summary.dat')],
         about='composition_umd.py on pyrolite + H2 at 4000 K expanded to 0.2 g/cm3 : melt and gas, fluctuations, voids, partition of H'),
    Case('exchange_h2',['pyrolite+H2.test.umd.dat','hydrousandesite.a15.0.T3000K.reduced.-A-.outcar.umd.dat'],[('gofr_umd.py',['-f','pyrolite+H2.test.umd.dat']),
                                                      ('analyze_gofr.py',['-f','pyrolite+H2.test.gofr.dat']),
                                                      ('bonding_umd.py',['-f','pyrolite+H2.test.umd.dat','-c','pyrolite+H2.test.maximum_bonds.inp']),
                                                      ('exchange_umd.py',['-f','pyrolite+H2.test.umd.dat']),
                                                      ('exchange_umd.py',['-f','hydrousandesite.a15.0.T3000K.reduced.-A-.outcar.umd.dat'])],
         [exch_summary('pyrolite+H2.test.exchange/summary.dat'),exch_summary('hydrousandesite.a15.0.T3000K.reduced.-A-.outcar.exchange/summary.dat')],
         about='exchange_umd.py on pyrolite + H2 (a small melt cluster in its gas ; molecules from its bonds) and on the hydrous andesite (a bulk melt : everything condensed)'),
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
          ('exchange_umd.py',['-S','glycine.r0.step.dat','-P','glycine.r1.popul.dat','-l','12']),
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
import numpy as np

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


def t_rotors(folder,outs):#rotations_umd.py on the synthetic rotors : outputs 5 (rotors), 6 (methane), 7 (sixstate), 8 to 11 (gaps)
    t = truth_of(folder,'rotors')
    out = []
    def col(name,column):
        return P.column(os.path.join(folder,name),column,sep='\t')
    def decay(name,label,pairs_):#time correlations against factor^lag ; 8 errors of 4 blocks, at least 0.01
        time = col(name,'time_(fs)')
        for column, factor in pairs_:
            c, e = col(name,column), col(name,column+'_err')
            for lag in (5,20,50):
                out.append(('%s %s at %g fs' % (label,column,time[lag]),c[lag],factor**lag,max(8*e[lag],0.01)))
    summary = P.table(os.path.join(folder,'rotors.OH2.rot.dat'),'quantity',sep='\t')
    r = t['rotors']
    for key, exact in (('neighbours',r['neighbours']),('cos_first_shell',r['cos_shell']),('P2_first_shell',1.0),('links_per_unit',r['links_per_unit']),
                       ('cos_linked',r['cos_linked'])):
        out.append(('rotors '+key,summary[key]['value'],exact,1e-9))
    dist, pairs_ = col('rotors.OH2.rot.pair.dat','r_(A)'), col('rotors.OH2.rot.pair.dat','pairs_per_snapshot')
    for column, near, far in (('cos',r['cos_3.0'],r['cos_3.6']),('P2',r['P2_3.0'],r['P2_3.6']),('dipolar',r['dipolar_3.0'],r['dipolar_3.6'])):
        values = col('rotors.OH2.rot.pair.dat',column)
        for lo, hi, exact in ((2.8,3.2,near),(3.4,3.8,far)):
            found = [v for x,n,v in zip(dist,pairs_,values) if lo < x < hi and n > 0]
            out.append(('rotors %s of the pairs at %.1f A (worst bin of %d)' % (column,0.5*(lo+hi),len(found)),max(found,key=lambda v: abs(v-exact)) if found else float('nan'),exact,1e-9))
    kirk = col('rotors.OH2.rot.pair.dat','Kirkwood')
    for x0, key in ((3.3,'kirkwood_3.3'),(3.9,'kirkwood_3.9')):
        out.append(('rotors Kirkwood factor within %.1f A' % x0,min(zip(dist,kirk),key=lambda p: abs(p[0]-x0))[1],r[key],1e-9))
    decay('rotors.OH2.rot.time.dat','rotors',[('C1',r['factors']['1']),('C2',r['factors']['2']),('bond_C1',r['factors']['1']),('bond_C2',r['factors']['2'])])
    m = t['methane']
    dist, pairs_ = col('methane.CH4.rot.pair.dat','r_(A)'), col('methane.CH4.rot.pair.dat','pairs_per_snapshot')
    for column, near, far in (('K3',m['K3_4.0'],m['K3_4.6']),('K4',m['K4_4.0'],m['K4_4.6'])):
        values = col('methane.CH4.rot.pair.dat',column)
        for lo, hi, exact in ((3.8,4.2,near),(4.4,4.8,far)):
            found = [v for x,n,v in zip(dist,pairs_,values) if lo < x < hi and n > 0]
            out.append(('methane %s of the pairs at %.1f A (worst bin of %d)' % (column,0.5*(lo+hi),len(found)),max(found,key=lambda v: abs(v-exact)) if found else float('nan'),exact,1e-9))
    found = re.search(r'l=3 ([0-9.]+), l=4 ([0-9.]+) ; orders used : 3, 4',outs[6])
    out.append(('methane mean of P3 between the bonds of a unit',float(found.group(1)) if found else float('nan'),m['self']['3'],1e-3))
    out.append(('methane mean of P4 between the bonds of a unit',float(found.group(2)) if found else float('nan'),m['self']['4'],1e-3))
    decay('methane.CH4.rot.time.dat','methane',[('C3',m['factors']['3']),('C4',m['factors']['4']),('bond_C1',m['factors']['1']),('bond_C2',m['factors']['2'])])
    six = t['sixstate']
    decay('sixstate.OH2.rot.time.dat','sixstate',[('C1',six['factor']),('C2',six['factor'])])
    summary = P.table(os.path.join(folder,'sixstate.OH2.rot.dat'),'quantity',sep='\t')
    out.append(('sixstate entropy of the map + ln(1800 bins)',summary['map_entropy_(kB)']['value']+math.log(1800.0),six['entropy_plus_ln_bins'],0.01))
    out.append(('sixstate cos between first neighbours (independent)',summary['cos_first_shell']['value'],0.0,max(8*summary['cos_first_shell']['error'],0.02)))
    for index, key in ((8,'-g 50'),(9,'-g 5'),(10,'-g 200'),(11,'-g 5 -m 150')):
        out.append(('gaps '+key+' : lives',P.find(outs[index],r'(\d+) lives of'),t['gaps']['lives'][key],0))
        if key in t['gaps']['left_out']:
            out.append(('gaps '+key+' : snapshots of units left out',P.find(outs[index],r': (\d+) snapshots of units left out'),t['gaps']['left_out'][key],0))
    return out


def t_rotors2pt(folder,outs):#rotations_umd.py, entropy (2PT) of the synthetic rigid rotors : the truth is computed from the correlations the run has
    t = truth_of(folder,'rotors2pt')
    out = []
    for kind, unit in (('h2o','OH2'),('ch4','CH4'),('oh','OH')):
        r = t[kind]
        q = P.table(os.path.join(folder,'rotors2pt_%s.%s.rot.dat' % (kind,unit)),'quantity',sep='\t')
        for key, exact, tol in (('symmetry_number',r['symmetry_number'],0),('volume_per_unit_(A3)',r['volume'],1e-5*r['volume']),
                                ('I_1_(amu.A2)',r['moments'][0],1e-6),('I_2_(amu.A2)',r['moments'][1],1e-6),('I_3_(amu.A2)',r['moments'][2],1e-6),
                                ('S_rigid_rotor_(kB)',r['S_rigid_rotor'],1e-5),('T_kin_translation_(K)',r['T_trn'],1e-4),('T_kin_rotation_(K)',r['T_rot'],1e-4),
                                ('fluidicity_translation',r['f_trn'],1e-6),('fluidicity_rotation',r['f_rot'],1e-6),
                                ('S_translation_(kB)',r['S_trn'],2e-3*r['S_trn']),('S_rotation_(kB)',r['S_rot'],2e-3*r['S_rot'])):
            out.append(('%s %s' % (kind,key),q[key]['value'],exact,tol))
        for key, name in (('f_trn','fluidicity_translation'),('f_rot','fluidicity_rotation'),('S_trn','S_translation_(kB)'),('S_rot','S_rotation_(kB)')):
            out.append(('%s %s against exp(-t/tau) : 5 block errors' % (kind,name),q[name]['value'],r['expected'][key],5*q[name]['error']))
    return out


def t_msrd(folder,outs):#msrd_umd.py -c Si -a O -w 120 on the synthetic SiO4 and Si2O7 : every MSRD at the lags of the truth and the static values, computed
    #there by brute force (body frame by Horn's quaternions) ; exact to the rounding of the UMD file
    t = truth_of(folder,'msrd')
    lines = [l for l in open(os.path.join(folder,'msrd.SiO4.msrd.time.dat')) if not l.startswith('#')]
    head = lines[0].rstrip('\n').split('\t')
    data = np.array([[float(x) for x in l.split('\t')] for l in lines[1:]])
    out = []
    for col,values in t['columns'].items():
        for lag,exact in zip(t['lags'],values):
            out.append(('msrd %s at %d fs' % (col,lag),data[lag,head.index(col)],exact,1e-6*abs(exact)+1e-9))
    q = P.table(os.path.join(folder,'msrd.SiO4.msrd.dat'),'quantity',sep='\t')
    for key,exact in t['summary'].items():
        out.append(('msrd '+key,q[key]['value'],exact,1e-6*abs(exact)+1e-9))
    out += [('msrd lives',q['lives']['value'],t['lives'],0),('msrd polyhedra per snapshot',q['polyhedra_per_snapshot']['value'],t['polyhedra_per_snapshot'],1e-12)]
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


def t_composition(folder,outs):#composition_umd.py on the synthetic grids : statistical results within 5 of their errors (autocorrelation time) or
    #8 (4 blocks, 1/e times) ; exact results to the digits of the summary
    t = truth_of(folder,'composition')
    out = []
    for name in ('comp_ideal','comp_segregated'):
        n = comp_numbers(os.path.join(folder,name+'.comp.4x4x4','summary.dat'))
        for el, exact in t[name]['ratio'].items():
            out.append((name+' variance / binomial '+el,n['ratio '+el],exact,max(5*n['ratio_err '+el],1e-4)))
        for pair, exact in t[name]['corr'].items():
            out.append((name+' correlation '+pair,n['corr '+pair],exact,max(5*n['corr_err '+pair],0.01)))
    n = comp_numbers(os.path.join(folder,'comp_ideal.comp.4x4x4','summary.dat'))
    for el in ('Ar','Ne','all'):
        out.append(('comp_ideal 1/e time of the counts of '+el,n['count_'+el],t['comp_ideal']['decay_fs'],8*n['count_'+el+'_err']))
    h, n = t['comp_hole'], comp_numbers(os.path.join(folder,'comp_hole.comp.4x4x4','summary.dat'))
    for key, column in (('void_fraction','void_fraction'),('void_regions','void_regions'),('dense_regions','dense_regions'),
                        ('largest_void','largest_void_(A3)'),('largest_dense','largest_dense_(A3)'),('largest_dense_melt_fraction','largest_dense_melt_fraction'),
                        ('dense_density','dense_density_(g/cm3)'),('void_density','void_density_(g/cm3)'),('melt_density','melt_density_seen_by_its_mass_(g/cm3)')):
        out.append(('comp_hole '+key,n[column],h[key],1e-4*abs(h[key])+1e-9))
    for el, exact in h['in_voids'].items():
        out.append(('comp_hole fraction of '+el+' in the voids',n['in_voids '+el],exact,1e-4))
        out.append(('comp_hole fraction of '+el+' in the largest dense region',n['in_largest_dense '+el],1.0-exact,1e-4))
    return out


def t_alloy(folder,outs):#composition_umd.py -n 3 -e Ni:1.5 on the synthetic Fe-Ni : random alloy and segregation in the second half ;
    #statistical results within 5 of their errors (autocorrelation time), exact ones to the digits of the summary
    t = truth_of(folder,'alloy')
    out = []
    n = comp_numbers(os.path.join(folder,'alloy_random.comp.3x3x3','summary.dat'))
    r = t['alloy_random']
    for el in ('Fe','Ni'):
        out.append(('alloy_random concentration ratio '+el,n['conc_ratio '+el],r['conc_ratio'][el],5*n['conc_ratio_err '+el]))
        out.append(('alloy_random count ratio '+el,n['ratio '+el],r['count_ratio'][el],max(5*n['ratio_err '+el],1e-4)))
        out.append(('alloy_random no change during the run '+el,float(n['conc_change '+el] < 4.0),1.0,0))
    out.append(('alloy_random correlation Fe-Ni',n['corr Fe-Ni'],r['corr_Fe-Ni'],1e-4))
    out.append(('alloy_random fraction enriched in Ni (1.5) expected for a random alloy',n['random_alloy_fraction enriched_Ni_1.5'],r['enriched_Ni_1.5'],1e-4))
    out.append(('alloy_random fraction enriched in Ni (1.5)',n['time_fraction enriched_Ni_1.5'],r['enriched_Ni_1.5'],0.01))
    n = comp_numbers(os.path.join(folder,'alloy_segregated.comp.3x3x3','summary.dat'))
    g = t['alloy_segregated']
    errs = n['conc_quarter_errs Ni']
    for k, exact in ((1,g['first_half']),(2,g['first_half']),(3,g['second_half']),(4,g['second_half'])):
        out.append(('alloy_segregated concentration ratio of Ni, quarter %d' % k,n['conc_quarter%d Ni' % k],exact,5*errs[k-1]))
    out.append(('alloy_segregated change during the run detected',float(n['conc_change Ni'] > 4.0),1.0,0))
    return out


def t_blink(folder,outs):#composition_umd.py -n 4 -g He on the blinking hole : periods exact
    t = truth_of(folder,'blink')
    n = comp_numbers(os.path.join(folder,'blink.comp.4x4x4','summary.dat'))
    out = []
    for state in ('void','dense'):
        out.append(('blink complete '+state+' periods',n['periods '+state],t[state+'_periods'],0))
        out.append(('blink cut '+state+' periods',n['cut_periods '+state],t[state+'_cut'],0))
        out.append(('blink mean '+state+' period (fs)',n['mean_period '+state],t[state+'_mean'],1e-6))
    out.append(('blink void time fraction',n['time_fraction void'],t['void_time_fraction'],1e-4))
    return out


def t_exchange(folder,outs):#exchange_umd.py -g He -t 0 on He jumping between the depths of an Ar slab and of the gas (Markov chain) : rates and stays
    #within 5 of their errors (Poisson ; std / sqrt(stays)) ; the geometry of the slab (interface 8 nm^2, half the volume gas, bulk density)
    t = truth_of(folder,'exchange')
    n = exch_numbers(os.path.join(folder,'exchange.exchange','summary.dat'))
    #the errors : columns named error after each value ; read by position
    rows = [l.rstrip('\n').split('\t') for l in open(os.path.join(folder,'exchange.exchange','summary.dat'))]
    head = next(r for r in rows if r[0] == 'element' and 'evaporations' in r)
    he = next(r for r in rows[rows.index(head):] if r[0] == 'He')
    vals = dict(zip(['evaporations','condensations','k_evaporation','k_evaporation_err','k_condensation','k_condensation_err','stay_condensed',
                     'stay_condensed_err','stay_gas','stay_gas_err','gas_fraction','net'],[float(v) for v in he[1:]]))
    out = [('exchange He '+key,vals[key],t[key],5*vals[key+'_err']) for key in ('k_evaporation','k_condensation','stay_condensed','stay_gas')]
    out.append(('exchange He gas fraction',vals['gas_fraction'],t['gas_fraction'],0.02))
    out.append(('exchange Ar never leaves the slab',n['evaporations Ar'],0.0,0))
    out.append(('exchange interface area (nm2)',n['area_(nm2)'],t['area_nm2'],0.02*t['area_nm2']))
    out.append(('exchange gas volume fraction',n['gas_volume_fraction'],t['gas_volume_fraction'],0.01))
    out.append(('exchange bulk density',n['bulk_density_(atoms/A3)'],t['bulk_density'],0.01*t['bulk_density']))
    out.append(('exchange gas temperature (K)',n['T_gas_(K) He'],t['T_gas'],0.02*t['T_gas']))
    out.append(('exchange partial pressure of He (GPa)',n['P_(GPa) He'],t['P_He_GPa'],5*n['P_err He']+0.01*t['P_He_GPa']))
    out.append(('exchange partition coefficient of He',n['D He'],t['D_He'],5*n['D_err He']))
    out.append(('exchange Hertz-Knudsen flux of He',n['J_HK_(1/ps/nm2) He'],t['J_HK'],0.08*t['J_HK']))
    out.append(('exchange no Ar in the gas',n['N_gas Ar'],0.0,0))
    return out


def t_twophase(folder,outs):#exchange_umd.py -c Fe -t 0 on He jumping between an Fe slab and a Si-O slab of the same density : rates and stays within 5
    #of their errors, D within 5 errors and 2 % (the zones of the rough surfaces), K and D_mass exact from the zone contents, the dissolved atoms in
    #their host ; without -c (-g He) everything condensed
    t = truth_of(folder,'twophase')
    n = exch_numbers(os.path.join(folder,'twophase.exchange.Fe','summary.dat'))
    out = [('twophase He '+k,n[c+' He'],t[k],5*n[e+' He']) for k, c, e in (('k_A_to_B','k_A_to_B_(1/ps)','error'),('stay_A','mean_stay_A_(fs)','error'))]
    rows = [l.rstrip('\n').split('\t') for l in open(os.path.join(folder,'twophase.exchange.Fe','summary.dat'))]
    head = next(r for r in rows if r[0] == 'element' and 'A_to_B' in r)
    he = dict(zip(head[1:],[float(v) for v in next(r for r in rows[rows.index(head):] if r[0] == 'He')[1:]]))
    errs = [float(v) for v in next(r for r in rows[rows.index(head):] if r[0] == 'He')[1:]]
    out.append(('twophase He k_B_to_A',he['k_B_to_A_(1/ps)'],t['k_B_to_A'],5*errs[5]))
    out.append(('twophase He stay_B',he['mean_stay_B_(fs)'],t['stay_B'],5*errs[9]))
    out.append(('twophase He fraction of the time in B',n['B He'],t['B_fraction'],0.02))
    out.append(('twophase Fe dissolved in B',n['B Fe'],t['B_fraction_Fe'],1e-4))
    out.append(('twophase O dissolved in A',n['A O'],t['A_fraction_O'],1e-4))
    out.append(('twophase Fe never changes phase',n['A_to_B Fe']+n['B_to_A Fe'],0.0,0))
    out.append(('twophase O never changes phase',n['A_to_B O']+n['B_to_A O'],0.0,0))
    out.append(('twophase interface area (nm2)',n['area_(nm2)'],t['area_nm2'],0.02*t['area_nm2']))
    out.append(('twophase B volume fraction',n['B_volume_fraction'],t['B_volume_fraction'],0.01))
    out.append(('twophase bulk density of A',n['bulk_density_A_(atoms/A3)'],t['bulk_density'],0.01*t['bulk_density']))
    out.append(('twophase zone of A (A3)',n['A_zone_(A3)'],t['zone'],0.02*t['zone']))
    out.append(('twophase zone of B (A3)',n['B_zone_(A3)'],t['zone'],0.02*t['zone']))
    out.append(('twophase He in A',n['N_A He'],t['N_A_He'],0.01))
    out.append(('twophase He in B',n['N_B He'],t['N_B_He'],0.01))
    out.append(('twophase D of He',n['D He'],t['D_He'],5*n['D_err He']+0.02*t['D_He']))
    for e in ('Fe','O','He'):
        out.append(('twophase K of '+e,n['K '+e],t['K_'+e],1e-3*t['K_'+e]))
        out.append(('twophase D_mass of '+e,n['D_mass '+e],t['D_mass_'+e],1e-3*t['D_mass_'+e]))
    m = exch_numbers(os.path.join(folder,'twophase.exchange','summary.dat'))
    out.append(('twophase without -c : He never leaves the condensed phase',m['condensed He'],1.0,0))
    return out


def t_interface(folder,outs):#exchange_umd.py -g He -t 0 on He scripted through the interface of the slab : passages and times exact
    t = truth_of(folder,'interface')
    n = exch_numbers(os.path.join(folder,'interface.exchange','summary.dat'))
    out = []
    for p in ('G-I-G','G-I-L','L-I-L','L-I-G'):
        out.append(('interface '+p+' passages',n[p+' He'],t[p]['n'],0))
        out.append(('interface '+p+' time (fs)',n[p+'_time_(fs) He'],t[p]['time'],1e-6))
    out.append(('interface sticking',n['sticking He'],t['sticking'],1e-4))
    out.append(('interface escape',n['escape He'],t['escape'],1e-4))
    out.append(('interface area (nm2)',n['area_(nm2)'],t['area_nm2'],0.02*t['area_nm2']))
    out.append(('interface bulk density',n['bulk_density_(atoms/A3)'],t['bulk_density'],0.01*t['bulk_density']))
    return out


def t_persist(folder,outs):#composition_umd.py -n 4 -g none on sub-cells of two fixed weights with atoms drawn anew : persistent fraction within 5 of
    #its error, the decay of the fluctuations within 8 of its error (4 blocks), the correlation at the end close to the persistent fraction
    t = truth_of(folder,'persist')
    n = comp_numbers(os.path.join(folder,'persist.comp.4x4x4','summary.dat'))
    return [('persist persistent fraction of the counts',n['count_Ar_persistent'],t['persistent'],5*n['count_Ar_persistent_err']),
            ('persist decay of the fluctuations (fs)',n['count_Ar_fluctuating'],t['fluctuating_decay'],8*n['count_Ar_fluctuating_err']),
            ('persist correlation at the end',n['count_Ar_end'],t['persistent'],0.05)]


def t_bubbles(folder,outs):#composition_umd.py -n 6 -g none on the scripted bubbles : the regions followed in time, exact
    t = truth_of(folder,'bubbles')
    n = comp_numbers(os.path.join(folder,'bubbles.comp.6x6x6','summary.dat'))
    out = [('bubbles '+key,n['regions_'+key+' void'],t[key],0) for key in ('complete','new','split','vanished','merged','short_lived')]
    out += [('bubbles mean lifetime (fs)',n['regions_mean_lifetime void'],t['mean_lifetime'],1e-6),
            ('bubbles longest lifetime (fs)',n['regions_longest void'],t['longest'],1e-6),
            ('bubbles one dense region',n['regions_regions dense'],t['dense_regions'],0)]
    rows = [l.rstrip('\n').split('\t') for l in open(os.path.join(folder,'bubbles.comp.6x6x6','regions.dat')) if not l.startswith('#')]
    head = rows[0]
    life = sorted(float(r[head.index('lifetime_(fs)')]) for r in rows[1:] if r[1] == 'void' and r[head.index('complete')] == '1')
    out.append(('bubbles lifetimes of the complete bubbles',float(life == sorted(t['lifetimes'])),1.0,0))
    vol = [float(r[head.index('max_volume_(A3)')]) for r in rows[1:] if r[1] == 'void' and r[head.index('birth_(fs)')] == '50.00']
    out.append(('bubbles volume of A (A3)',vol[0] if vol else float('nan'),t['volume_A'],0.01))
    return out


def t_degassing(folder,outs):#exchange_umd.py -n 4 -g He -t 0 -w 150 on He degassing from an Ar slab : the events of every window against their exact
    #expected numbers (within 5 Poisson errors), a net flux in the first window, equilibrium found later
    t = truth_of(folder,'degassing')
    rows = [l.rstrip('\n').split('\t') for l in open(os.path.join(folder,'degassing.exchange','fluxes.dat')) if not l.startswith('#')]
    head, data = rows[0], [[float(v) for v in r] for r in rows[1:]]
    span = t['window_fs']/1000.0
    col = lambda name: head.index(name)
    out = []
    for k, r in enumerate(data):
        for name, exact in (('evaporation',t['evaporations'][k]),('condensation',t['condensations'][k])):
            got = r[col(name+'_(1/ps)_He')]*span
            out.append(('degassing window %d %ss' % (k,name),got,exact,5*np.sqrt(max(exact,1.0))))
    first = data[0]
    out.append(('degassing net flux of the first window positive (> 3 errors)',float(first[col('net_(1/ps)_He')] > 3*first[col('net_err_He')]),1.0,0))
    summary = open(os.path.join(folder,'degassing.exchange','summary.dat')).read()
    found = re.search(r'\nHe\t[^\n]*\t([0-9.]+|not reached)\n',summary.split('equilibrium_from_(fs)')[1])
    out.append(('degassing equilibrium found after the first window',float(found is not None and found.group(1) != 'not reached' and float(found.group(1)) > 0),1.0,0))
    return out


def t_molecules(folder,outs):#exchange_umd.py -g H,He -t 0 on H2 and He jumping between an Ar slab and the gas (bonds from bonding_umd.py) : the molecular
    #events and gas exactly half those of the H atoms, equal for He ; -m clusters : the slab condensed, H2 and He always gas, no exchange
    n = exch_numbers(os.path.join(folder,'molecules.exchange','summary.dat'))
    rows = [l.rstrip('\n').split('\t') for l in open(os.path.join(folder,'molecules.exchange','molecules.dat')) if not l.startswith('#')]
    m = {r[0]:dict(zip(rows[0][1:],[float(v) for v in r[1:]])) for r in rows[1:]}
    out = [('molecules H2 evaporations = H / 2',m['H_2']['evaporations'],n['evaporations H']/2,1e-9),
           ('molecules H2 condensations = H / 2',m['H_2']['condensations'],n['condensations H']/2,1e-9),
           ('molecules He evaporations = He',m['He_1']['evaporations'],n['evaporations He'],1e-9),
           ('molecules H2 in the gas = H / 2',m['H_2']['in_gas'],n['N_gas H']/2,3e-3),          #N_gas is printed with 4 digits
           ('molecules H2 pressure = H / 2',m['H_2']['P_(GPa)'],n['P_(GPa) H']/2,1e-4*n['P_(GPa) H']),
           ('molecules some exchange',float(n['evaporations H'] > 100),1.0,0)]
    #pressure.dat, zone definition : the gas molecules of every snapshot = the H2 and He of the gas (the Xe adsorbed : 1 atom), T the kinetic
    #temperature of their atoms, P = N k T/V - (N_H2 a_H2^1/2 + N_He a_He^1/2)^2/(NA^2 V^2), a of the CRC Handbook (L^2 bar/mol^2 = 0.1 Pa m^6/mol^2)
    truth = json.load(open(os.path.join(folder,'molecules.truth.json')))
    rows = [l.rstrip('\n').split('\t') for l in open(os.path.join(folder,'molecules.exchange','pressure.dat')) if not l.startswith('#')]
    head, data = rows[0], np.array([[float(v) for v in r] for r in rows[1:]])
    col = {h:k for k, h in enumerate(head)}
    counts = np.array(truth['gas_H2_He'],dtype=float)
    T = truth['T_gas_molecules']
    V = data[:,col['V_gas_(A3)_zone']]
    kB, NA, aH2, aHe = 8.617333262e-5, 6.02214076e23, 0.2476, 0.0346
    P = counts.sum(axis=1)*kB*T/V*160.21766208-0.1*(counts[:,0]*np.sqrt(aH2)+counts[:,1]*np.sqrt(aHe))**2/(NA**2*(V*1e-30)**2)*1e-9
    attraction = P/(counts.sum(axis=1)*kB*T/V*160.21766208)
    out += [('pressure : gas molecules of every snapshot = H2 + He of the gas',float(np.abs(data[:,col['molecules_zone']]-counts.sum(axis=1)).max()),0.0,0),
            ('pressure : the Xe adsorbed in every snapshot',float(np.abs(data[:,col['adsorbed_atoms_zone']]-truth['adsorbed_Xe']).max()),0.0,0),
            ('pressure : Van der Waals P of every snapshot (largest relative deviation)',float(np.abs(data[:,col['P_vdW_(GPa)_zone']]/P-1.0).max()),0.0,1e-4),
            ('pressure : the attraction lowers P (mean P_vdW/P_ideal < 1)',float(attraction.mean() < 1.0),1.0,0),
            ('adsorbed : the Xe bonded to the slab never gas',n['gas Xe'],0.0,0),
            ('adsorbed : the Xe in the interface all the time',n['interface Xe'],1.0,0),
            ('adsorbed : the Xe never evaporates',n['evaporations Xe'],0.0,0)]
    c = exch_numbers(os.path.join(folder,'molecules.exchange.clusters','summary.dat'))
    out += [('molecules clusters : Ar condensed',c['condensed Ar'],1.0,0),('molecules clusters : H gas',c['gas H'],1.0,0),
            ('molecules clusters : He gas',c['gas He'],1.0,0),('molecules clusters : no exchange of H',c['evaporations H'],0.0,0)]
    return out


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
    Case('syn_rotations',[],[('@synthetic.py',['rotors'])]+[('bonding_umd.py',['-f',f+'.umd.dat','-l',l]) for f,l in
                             (('rotors','1.2'),('methane','1.3'),('sixstate','1.2'),('gaps','1.2'))]
                            +[('rotations_umd.py',['-f','rotors.umd.dat','-c','O','-a','H','-r','4.0']),
                              ('rotations_umd.py',['-f','methane.umd.dat','-c','C','-a','H','-r','5.0']),
                              ('rotations_umd.py',['-f','sixstate.umd.dat','-c','O','-a','H','-r','5.5'])]
                            +[('rotations_umd.py',['-f','gaps.umd.dat','-c','O','-a','H','-r','7']+o) for o in
                              (['-g','50'],['-g','5'],['-g','200'],['-g','5','-m','150'])],truth=[t_rotors],
         about='synthetic rotors : pair correlations of dipoles and of tetrahedra, time correlations of rotational diffusion and of jumps, map entropy, gaps'),
    Case('syn_rotations2pt',[],[('@synthetic.py',['rotors2pt'])]+[('bonding_umd.py',['-f','rotors2pt_'+m+'.umd.dat','-l','1.2']) for m in ('h2o','ch4','oh')]
                               +[('rotations_umd.py',['-f','rotors2pt_h2o.umd.dat','-c','O','-a','H','-w','300','-t','300','-r','13','-V','bader2pt']),
                                 ('rotations_umd.py',['-f','rotors2pt_ch4.umd.dat','-c','C','-a','H','-w','300','-t','300','-r','13','-V','C:20,H:5']),
                                 ('rotations_umd.py',['-f','rotors2pt_oh.umd.dat','-c','O','-a','H','-w','300','-t','300','-r','13','-V','25'])],
         truth=[t_rotors2pt],about='synthetic rigid rotors (H2O, CH4, OH) with Ornstein-Uhlenbeck translation and rotation : 2PT entropies, fluidicities, symmetry numbers, Bader volumes'),
    Case('syn_msrd',[],[('@synthetic.py',['msrd']),('bonding_umd.py',['-f','msrd.umd.dat','-l','2.1']),
                        ('msrd_umd.py',['-f','msrd.umd.dat','-c','Si','-a','O','-w','120','-p'])],truth=[t_msrd],
         about='synthetic SiO4 and Si2O7 turning, drifting and vibrating in a deforming triclinic cell, with a bridged gap and a broken life : MSRD of the bonds (lab and body frames, bridging and non-bridging) and static values exact ; the plot (-p) written'),
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
    Case('syn_gz',[],[('@synthetic.py',['idealgas_gz']),('gofr_umd.py',['-f','idealgas.umd.dat.gz']),('check_umd.py',['-f','idealgas.umd.dat.gz']),
                      ('composition_umd.py',['-f','idealgas.umd.dat.gz','-n','3'])],truth=[t_idealgas],
         about='the synthetic ideal gas compressed (.umd.dat.gz, decompressed into a temporary file) : g(r) = 1, running coordination exact, outputs named from idealgas.umd.dat'),
    Case('syn_alloy',[],[('@synthetic.py',['alloy']),('composition_umd.py',['-f','alloy_random.umd.dat','-n','3','-e','Ni:1.5,Ni:0.5']),
                         ('composition_umd.py',['-f','alloy_segregated.umd.dat','-n','3','-e','Ni:1.5'])],truth=[t_alloy],
         about='synthetic Fe-Ni on a lattice : random alloy (concentration ratio 1, Fe-Ni correlation -1, enriched fraction hypergeometric) and Ni segregating in the second half (quarters exact, change detected)'),
    Case('syn_exchange',[],[('@synthetic.py',['exchange']),('exchange_umd.py',['-f','exchange.umd.dat','-g','He','-t','0'])],truth=[t_exchange],
         about='synthetic He jumping between the depths of an Ar slab and of the gas (Markov chain) : rates, stays, gas fraction ; interface area, gas volume, bulk density'),
    Case('syn_twophase',[],[('@synthetic.py',['twophase']),('exchange_umd.py',['-f','twophase.umd.dat','-c','Fe','-t','0']),
                            ('exchange_umd.py',['-f','twophase.umd.dat','-g','He','-t','0'])],truth=[t_twophase],
         about='synthetic He jumping between an Fe slab and a Si-O slab of the same density (two condensed phases, -c Fe), O dissolved in the Fe and Fe in the '
               'Si-O : rates, stays, D, K, D_mass exact, the dissolved atoms in their host ; without -c everything condensed'),
    Case('syn_interface',[],[('@synthetic.py',['interface']),('exchange_umd.py',['-f','interface.umd.dat','-g','He','-t','0'])],truth=[t_interface],
         about='synthetic He scripted through the interface of an Ar slab : bounces, condensations, visits, evaporations, times in the interface, sticking and escape exact'),
    Case('syn_persist',[],[('@synthetic.py',['persist']),('composition_umd.py',['-f','persist.umd.dat','-n','4','-g','none'])],truth=[t_persist],
         about='synthetic sub-cells of two fixed weights, atoms drawn anew 10 % per step : persistent fraction of the heterogeneity, decay of the fluctuations 0.9^t'),
    Case('syn_bubbles',[],[('@synthetic.py',['bubbles']),('composition_umd.py',['-f','bubbles.umd.dat','-n','6','-g','none'])],truth=[t_bubbles],
         about='synthetic bubbles of a scripted history (lives, a merge, a split, cut ones) : the regions followed in time exact'),
    Case('syn_degassing',[],[('@synthetic.py',['degassing']),('exchange_umd.py',['-f','degassing.umd.dat','-g','He','-t','0','-w','150'])],
         truth=[t_degassing],about='synthetic He degassing from an Ar slab (Markov chain) : evaporation and condensation fluxes of every window exact, equilibrium found'),
    Case('syn_molecules',[],[('@synthetic.py',['molecules']),('bonding_umd.py',['-f','molecules.umd.dat','-c','molecules.cutoffs.inp']),
                             ('exchange_umd.py',['-f','molecules.umd.dat','-g','H,He,Xe','-t','0']),
                             ('exchange_umd.py',['-f','molecules.umd.dat','-g','H,He,Xe','-t','0','-m','clusters'])],truth=[t_molecules],
         about='synthetic H2 and He jumping between an Ar slab and the gas, a Xe in the gas bonded to the slab, bonds from bonding_umd.py : molecular '
               'events, gas and pressure, adsorbed Xe, modified Van der Waals pressure of every snapshot ; -m clusters'),
    Case('syn_workflow',[],[('@synthetic.py',['exchange']),('umd_workflow.py',['-f','exchange.umd.dat'])],[workflow(1)],
         about='the whole workflow on the synthetic Ar slab with He in the gas (no data needed : run on every system by the tests on GitHub) ; '
               'viscosity FAILED (no stress tensor in the synthetic run)'),
    Case('syn_blink',[],[('@synthetic.py',['blink']),('composition_umd.py',['-f','blink.umd.dat','-n','4','-g','He'])],truth=[t_blink],
         about='synthetic hole opening for 30 fs and closing for 20 fs, ten times : periods of the void and dense states exact'),
    Case('syn_composition',[],[('@synthetic.py',['composition'])]+[('composition_umd.py',['-f',f+'.umd.dat','-n','4']+o) for f, o in
                               (('comp_ideal',[]),('comp_segregated',[]),('comp_hole',['-g','He']))],truth=[t_composition],
         about='synthetic grids : ideal mixing (variance / binomial 1, correlation 0, time correlation 0.9^t), segregated halves (ratio 2), a hole across the periodic faces (voids, densities, partition)'),
    Case('syn_check',[],[('@synthetic.py',['faults']),('check_umd.py',['-f','faults_clean.umd.dat']),
                         ('check_umd.py',['-f','faults_bad.umd.dat'])],truth=[t_faults],
         about='synthetic random walk, clean and with a duplicated snapshot and a missed crossing : check_umd'),
]
