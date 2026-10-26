#!/usr/bin/env python3
###
#AUTHORS: RAZVAN CARACAS
###
#Automated workflow : runs the scripts of the package on UMD files, panel by panel, as in the flowchart (documentation/umd-flowchart.pdf),
#and writes a summary of the results. Each script runs in its own process, in the folder of the UMD file, with its output appended to
#umd_workflow.log of that folder. A step is done when its output files exist, newer than its inputs, written by the same command (kept in
#umd_workflow.done, with the files written) : it is then not run again (-o runs it again). It is skipped when its inputs are missing (a step
#before it failed). When a step runs again, the files it wrote in its former run and not now are removed (named in umd_workflow.done).
#The scripts often end without an error code, so a step succeeds only when it wrote its outputs.

import sys, getopt, os, glob, time, subprocess, json, re
import umd_processes_fast as umdpf
import crystallography as cr

SRC = os.path.dirname(os.path.abspath(__file__))
LOG = 'umd_workflow.log'
DONE = 'umd_workflow.done'
SUMMARY = 'umd_workflow.summary.dat'


class Step:#One script run : command line, input files it needs, output files it writes (names relative to the folder ; pattern :
    #outputs whose names are known only after the run, at least one)
    def __init__(self,script,args,needs,makes,pattern=None):
        self.script = script
        self.args = args
        self.needs = needs
        self.makes = makes
        self.pattern = pattern

    def outputs(self,folder):
        found = [f for f in self.makes if os.path.isfile(os.path.join(folder,f))]
        if self.pattern:
            found += [os.path.basename(f) for f in glob.glob(os.path.join(folder,self.pattern))]
        return found

    def command(self):
        return ' '.join([self.script]+self.args)

    def key(self):#the step in umd_workflow.done
        return self.makes[0] if self.makes else self.pattern

    def label(self,base):
        return ' '.join([self.script]+[x.replace(base,'<base>',1) if base else x for x in self.args])


def base_of(umdfile):
    return umdfile[:-8] if umdfile.endswith('.umd.dat') else os.path.splitext(umdfile)[0]


def read_done(folder):
    try:
        with open(os.path.join(folder,DONE)) as f:
            done = json.load(f)
        return {key : value if isinstance(value,dict) else {'command':value} for key, value in done.items()}
    except (OSError,ValueError):
        return {}


def run(step,folder,base,overwrite,results):#Runs one step ; its status is appended to results
    label = step.label(base)
    short = lambda f: f.replace(base,'<base>',1) if base else f
    mtime = lambda f: os.path.getmtime(os.path.join(folder,f))
    missing = [f for f in step.needs if not os.path.isfile(os.path.join(folder,f))]
    before = step.outputs(folder)
    done = read_done(folder)
    if missing:
        status = 'skipped : no '+', '.join(short(f) for f in missing)
    elif not overwrite and done.get(step.key(),{}).get('command') == step.command() and len(before) >= len(step.makes)+bool(step.pattern) and \
            all(mtime(f) <= min(mtime(g) for g in before) for f in step.needs):
        status = 'done before'                                                  #same command, all the outputs newer than the inputs
    else:
        start = time.time()
        env = dict(os.environ, MPLBACKEND='Agg')                                #no window : the plots are saved by the scripts
        with open(os.path.join(folder,LOG),'a') as log:
            log.write('\n'+'='*100+'\n'+time.strftime('%Y-%m-%d %H:%M:%S')+'   '+step.command()+'\n'+'='*100+'\n')
            log.flush()
            subprocess.run([sys.executable,os.path.join(SRC,step.script)]+step.args,cwd=folder,stdout=log,stderr=subprocess.STDOUT,
                           stdin=subprocess.DEVNULL,env=env)
        written = [f for f in step.outputs(folder) if mtime(f) >= start-1]
        absent = [f for f in step.makes if f not in written]
        if step.pattern and not any(f not in step.makes for f in written):
            absent.append(step.pattern)
        done = read_done(folder)
        stale = [f for f in done.get(step.key(),{}).get('files',[]) if f not in written and os.path.isfile(os.path.join(folder,f)) and mtime(f) < start-1]
        for f in stale:#files written by this step in a former run of the workflow and not now (e.g. species no longer found)
            os.remove(os.path.join(folder,f))
        if absent:
            status = 'FAILED : no '+', '.join(short(f) for f in absent)+' ; see '+LOG
            done.pop(step.key(),None)
        else:
            status = 'OK %7.1f s' % (time.time()-start)
            if step.pattern:
                status += ' : %d files %s' % (len(written)-len(step.makes),short(step.pattern))
            done[step.key()] = {'command':step.command(),'files':written}
        if stale:
            status += (' ; ' if ' : ' in status else ' : ')+'removed %d file(s) of the former run : ' % len(stale)+', '.join(short(f) for f in stale)
        with open(os.path.join(folder,DONE),'w') as f:
            json.dump(done,f,indent=0)
    head, _, detail = status.partition(' : ')
    print('   %-14s %s' % (head,label)+('\n'+' '*18+'('+detail+')' if detail else ''))
    results.append((folder,step.command(),status))


def roles(elements,central,adjacent):#Central and adjacent elements of the speciation : by default the cations around the anions
    if central is None:
        central = [el for el in elements if (cr.Elements2charge(el) or 0) > 0]
    if adjacent is None:
        adjacent = [el for el in elements if (cr.Elements2charge(el) or 0) < 0]
    return central, adjacent


def table(path,first):#Rows of a tab-separated file whose column names are on the line beginning with first : {first column : {name : value}}
    rows = {}
    if not os.path.isfile(path):
        return rows
    names = None
    with open(path) as f:
        for line in f:
            entry = line.rstrip('\n').split('\t')
            if entry[0] == first:
                names = entry
            elif names and len(entry) > 1 and not line.startswith('#'):
                rows[entry[0].strip()] = dict(zip(names,[x.strip() for x in entry]))
    return rows


def summary_row(folder,umd,central,adjacent):#Main results of one UMD file, from the output files of the scripts ; {} if there are none
    base = os.path.join(folder,base_of(umd))
    row = {}
    thermo = table(os.path.join(folder,'fullthermo.dat'),'file').get(umd,{})
    for name, column in (('nsteps','nsteps'),('density(g/cm3)','Density(g/cm3)'),('P(GPa)','Pressure(GPa)'),('P_err','err_P'),
                         ('T(K)','Temperature(K)'),('T_err','err_T'),('E(eV/atom)','InternEnergy(eV/atm)'),('E_err','err_E'),
                         ('Cv(NkB)','Cvm(Nkb)'),('Cv_err','stdev_Cvm')):
        if column in thermo:
            row[name] = thermo[column]
    for el, values in table(base+'.diffusion.dat','column').items():
        if '_' not in el:                                                           #the total MSD of each element, not the axes
            row['D_'+el+'(m2/s)'] = values.get('D_(m2/s)','nan')
            row['D_'+el+'_err'] = values.get('D_error_(m2/s)','nan')
    if os.path.isfile(base+'.visc.dat'):
        with open(base+'.visc.dat') as f:
            for line in f:
                found = re.search(r'viscosity = (\S+) \+/- (\S+) Pa\.s',line)
                if found:
                    row['eta(Pa.s)'], row['eta_err'] = found.group(1), found.group(2)
                    break
    pt = table(base+'.2pt.dat','element').get('all')
    if pt:
        row['S_2PT(kB/atom)'] = pt['S_2PT']
    quantum = table(base+'.vdos.thermo.dat','T(K)')                            #the last line beginning with each temperature : all the atoms
    if quantum:
        allatoms = list(quantum.values())[0]
        for name, column in (('dE_quantum(meV/atom)','dE'),('dS_quantum(kB/atom)','dS'),('dCv_quantum(kB/atom)','dCv')):
            row[name] = allatoms.get(column,'nan')
    if os.path.isfile(base+'.coordinations.inp'):
        with open(base+'.coordinations.inp') as f:
            for line in f:
                entry = line.split()
                if len(entry) == 3 and entry[0] in central and entry[1] in adjacent and entry[0] != entry[1]:
                    row['CN_'+entry[0]+'-'+entry[1]] = entry[2]
    species = table(base+'.r1.stat.dat','Cluster')
    if species:
        main = sorted(species.values(),key=lambda s: -float(s['Percent']))[:5]
        row['main_r1_species(%)'] = ' '.join(s['Cluster']+':'+'%.1f' % (100*float(s['Percent'])) for s in main)
    if row and os.path.isfile(base+'.check.dat'):
        with open(base+'.check.dat') as f:
            row['check_warnings'] = str(sum('WARNING' in line for line in f))
    return row


def write_summary(folder,rows,umds):#One line per UMD file, the columns of all the files (nan where a file has none)
    columns = []
    for row in rows.values():
        columns += [c for c in row if c not in columns]
    groups = ('nsteps','check_warnings','density','P(','P_err','T(','T_err','E(','E_err','Cv','D_','eta','S_2PT','dE','dS','dCv','CN_','main')
    columns.sort(key=lambda c: [c.startswith(g) for g in groups].index(True))  #by quantity, the elements in their order of appearance
    with open(os.path.join(folder,SUMMARY),'w') as f:
        f.write('# umd_workflow.py summary, '+time.strftime('%Y-%m-%d %H:%M')+' : main results of the UMD files of this folder, from the files of the scripts\n')
        f.write('# thermodynamics : fullthermo.dat (averages.py -a, E with the kinetic energy, Cv in N kB) ; D : <base>.diffusion.dat (msd_umd.py) ;\n')
        f.write('# eta : <base>.visc.dat (viscosity_umd.py) ; S_2PT and the quantum corrections (to add to the MD values) : <base>.2pt.dat,\n')
        f.write('# <base>.vdos.thermo.dat (vdos_thermo.py) ; CN : <base>.coordinations.inp (analyze_gofr.py), central - adjacent elements ;\n')
        f.write('# species : <base>.r1.stat.dat (speciation_and_angles.py -r 1), the 5 most frequent, % of the species time ; check_warnings :\n')
        f.write('# number of warnings of check_umd.py (<base>.check.dat) ; errors : _err columns (errors of the mean or of the blocks)\n')
        if os.path.isfile(os.path.join(folder,'thermo_integration.dat')):
            f.write('# free energy and entropy by thermodynamic integration of the runs at several lambda : thermo_integration.dat\n')
        f.write('\t'.join(['file']+columns)+'\n')
        for umd in umds:
            if umd in rows:
                f.write('\t'.join([umd]+[rows[umd].get(c,'nan') for c in columns])+'\n')


def main(argv):
    umdpf.headerumd()
    files = None
    panels = ''
    first = 0
    sampling = 1
    central = None
    adjacent = None
    overwrite = False
    usage = 'umd_workflow.py -f <UMD files> -g -t -d -i <InitialStep> -s <Sampling> -c <central elements> -a <adjacent elements> -o'
    try:
        opts, arg = getopt.getopt(argv,"hf:gtdi:s:c:a:o")
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('umd_workflow.py : automated analysis of UMD files, by panels of the flowchart of the package (documentation/umd-flowchart.pdf)')
            print(usage)
            print('-f : UMD files, separated by commas. Default : every *.umd.dat of the current folder.')
            print('-g : structural properties (geometry) : gofr_umd.py, analyze_gofr.py (bond cutoffs = first minima of the g(r)), bonding_umd.py')
            print('     with these cutoffs, speciation_and_angles.py -r 1 (central atoms and their adjacent neighbours, with the angles -t 1 and')
            print('     the bond lifetimes -l 1) and -r 0 (polymerization of all the atoms : molecules, -c all -a all), lifetimes.py (plot of the')
            print('     r1 species).')
            print('-t : transport properties : msd_umd.py (diffusion coefficients), viscosity_umd.py, and, when the r1 species exist (-g, or a')
            print('     former run), msd_clusters_umd.py and analyze_msd_clusters.py (diffusion of the chemical species).')
            print('-d : thermodynamic properties : averages.py -a (fullthermo.dat, one line per UMD file of the folder), vibrations_umd.py,')
            print('     vdos_thermo.py, and thermo_integration_umd.py on the UMD files of a folder when one of them has lambda_ThermoInt < 1.')
            print('     Without -g, -t, -d : the three panels.')
            print('check_umd.py always runs first ; its warnings are listed at the end. The workflow continues after a warning : read them.')
            print('-i : snapshots skipped at the beginning (equilibration), passed to all the scripts reading the UMD file. Default 0.')
            print('-s : one snapshot out of s is used by check_umd.py (checks of the atoms), gofr_umd.py, bonding_umd.py (hence the speciation)')
            print('     and msd_umd.py. Default 1. viscosity_umd.py and vibrations_umd.py always use every snapshot : their time correlations')
            print('     need the time step of the run (sampling would alias the vibrations and coarsen the Green-Kubo integrals).')
            print('-c, -a : central and adjacent elements of speciation_and_angles.py -r 1, separated by commas. Default : the cations (formal')
            print('     charge > 0, crystallography.Elements2charge) around the anions (charge < 0), e.g. Si,Mg,Al around O.')
            print('-o : runs every step again. Default : a step whose output files are newer than its inputs and were written by the same')
            print('     command (umd_workflow.done) is not run again, so that an interrupted workflow continues where it stopped and a panel')
            print('     can be added later ; a step whose options changed (-i, -s, -c, -a) runs again, and the steps that depend on it. The files')
            print('     a step wrote in a former run and not in the new one (e.g. cluster MSDs of species no longer found) are removed.')
            print('Each script runs with its default options, in the folder of its UMD file ; run it alone for other options. The output of the')
            print('scripts is appended to umd_workflow.log in that folder. No window is opened.')
            print('Output : umd_workflow.summary.dat in each folder, one line per UMD file of the folder with results : number of snapshots,')
            print('warnings of check_umd.py, density, P, T, E, Cv (with errors), D of each element, viscosity, S_2PT, quantum corrections,')
            print('coordination numbers of the central - adjacent pairs, the 5 most frequent r1 species ; nan where a result is missing.')
            sys.exit()
        elif opt == '-f':
            files = [x.strip() for x in arg.split(',') if x.strip()]
        elif opt in ('-g','-t','-d'):
            panels += opt[1]
        elif opt == '-i':
            first = int(arg)
        elif opt == '-s':
            sampling = int(arg)
        elif opt == '-c':
            central = [x.strip() for x in arg.split(',') if x.strip()]
        elif opt == '-a':
            adjacent = [x.strip() for x in arg.split(',') if x.strip()]
        elif opt == '-o':
            overwrite = True
    if files is None:
        files = sorted(glob.glob('*.umd.dat'))
    missing = [f for f in files if not os.path.isfile(f)]
    if missing or not files:
        print('Files not found :',', '.join(missing) if missing else 'no *.umd.dat in the folder')
        sys.exit()
    if sampling < 1:
        print('-s has to be at least 1')
        sys.exit()
    if not panels:
        panels = 'gtd'
    skip = ['-i',str(first)] if first else []
    sample = ['-s',str(sampling)] if sampling > 1 else []
    print('UMD files :',', '.join(files))
    print('panels :',', '.join({'g':'structural','t':'transport','d':'thermodynamic'}[p] for p in 'gtd' if p in panels))

    results = []
    byfolder = {}
    elementroles = {}
    for path in files:
        folder = os.path.dirname(os.path.abspath(path))
        umd = os.path.basename(path)
        byfolder.setdefault(folder,[]).append(umd)
        base = base_of(umd)
        MyCrystal = umdpf.Crystallization(path)[0]
        print('\n'+path+' :',MyCrystal.natom,'atoms,',' '.join(el+str(n) for el,n in zip(MyCrystal.elements,MyCrystal.types)))
        cen, adj = roles(MyCrystal.elements,central,adjacent)
        elementroles[(folder,umd)] = (cen,adj)
        steps = [Step('check_umd.py',['-f',umd]+skip+sample,[umd],[base+'.check.dat'])]
        popul = base+'.r1.popul.dat'

        if 'g' in panels:
            gofr, cutoffs, bonds = base+'.gofr.dat', base+'.maximum_bonds.inp', base+'.bonding.dat'
            steps += [Step('gofr_umd.py',['-f',umd]+skip+sample,[umd],[gofr]),
                      Step('analyze_gofr.py',['-f',gofr],[gofr],[cutoffs,base+'.average_bonds.inp',base+'.coordinations.inp']),
                      Step('bonding_umd.py',['-f',umd,'-c',cutoffs]+skip+sample,[umd,cutoffs],[bonds])]
            if cen and adj:
                steps += [Step('speciation_and_angles.py',['-f',umd,'-b',bonds,'-c',','.join(cen),'-a',','.join(adj),'-r','1','-t','1','-l','1'],
                               [umd,bonds],[popul,base+'.r1.stat.dat',base+'.r1.step.dat',base+'.angles.dat',base+'.bondslives.dat'])]
            else:
                print('   no central or adjacent elements for the r1 speciation (cations',cen,', anions',adj,') : give them with -c and -a')
            steps += [Step('speciation_and_angles.py',['-f',umd,'-b',bonds,'-c','all','-a','all','-r','0'],
                           [umd,bonds],[base+'.r0.popul.dat',base+'.r0.stat.dat',base+'.r0.step.dat'])]
            if cen and adj:
                steps += [Step('lifetimes.py',['-f',popul],[popul],[base+'.r1.lifetimes.pdf'])]

        if 't' in panels:
            steps += [Step('msd_umd.py',['-f',umd]+skip+sample,[umd],[base+'.msd.dat',base+'.msd.blocks.dat',base+'.diffusion.dat']),
                      Step('viscosity_umd.py',['-f',umd]+skip,[umd],[base+'.stress.scf.dat',base+'.visc.dat'])]
            if any(step.makes[0] == popul for step in steps) or os.path.isfile(os.path.join(folder,popul)):
                steps += [Step('msd_clusters_umd.py',['-f',umd,'-p',popul],[umd,popul],[],base+'.r1__*.msd.dat')]

        if 'd' in panels:
            steps += [Step('vibrations_umd.py',['-f',umd]+skip,[umd],[base+'.vels.scf.dat',base+'.vibr.dat',base+'.vibr.err.dat',base+'.vacf.diffusion.dat']),
                      Step('vdos_thermo.py',['-f',umd],[umd,base+'.vibr.err.dat'],[base+'.vdos.thermo.dat',base+'.2pt.dat'])]

        for step in steps:
            run(step,folder,base,overwrite,results)

    for folder, umds in byfolder.items():#Steps over all the UMD files of a folder, and the summary
        print('\n'+folder+' :')
        steps = []
        clusters = sorted(os.path.basename(f) for f in glob.glob(os.path.join(folder,'*.r1__*.msd.dat')))
        if 't' in panels and clusters:
            steps += [Step('analyze_msd_clusters.py',[],clusters,['diffusivities_clusters.txt'])]
        if 'd' in panels:
            steps += [Step('averages.py',['-a','-f',','.join(umds)]+skip,umds,['fullthermo.dat'])]
            if min(umdpf.Crystallization(os.path.join(folder,umd))[0].lambda_ThermoInt for umd in umds) < 1.0:
                steps += [Step('thermo_integration_umd.py',['-f',','.join(umds)]+skip,umds,['thermo_integration.dat'])]
        for step in steps:
            run(step,folder,'',overwrite,results)
        allumds = umds+[os.path.basename(f) for f in sorted(glob.glob(os.path.join(folder,'*.umd.dat'))) if os.path.basename(f) not in umds]
        rows = {}
        for umd in allumds:                                                     #also the UMD files of the folder done before
            if (folder,umd) not in elementroles:
                elementroles[(folder,umd)] = roles(umdpf.Crystallization(os.path.join(folder,umd))[0].elements,central,adjacent)
            row = summary_row(folder,umd,*elementroles[(folder,umd)])
            if row:
                rows[umd] = row
        if rows:
            write_summary(folder,rows,allumds)
            print('   summary of',len(rows),'UMD file(s) in',SUMMARY)

    failed = [r for r in results if r[2].startswith('FAILED')]
    print('\nSteps :',len(results),',',sum(r[2].startswith('OK') for r in results),'run,',
          sum(r[2] == 'done before' for r in results),'done before,',sum(r[2].startswith('skipped') for r in results),'skipped,',len(failed),'failed')
    for folder, command, status in failed:
        print('   FAILED',command,'in',folder)
    for path in files:#Warnings of the quality checks
        report = os.path.join(os.path.dirname(os.path.abspath(path)),base_of(os.path.basename(path))+'.check.dat')
        if os.path.isfile(report):
            with open(report) as f:
                warnings = [line.rstrip() for line in f if 'WARNING' in line]
            print('\ncheck_umd.py,',path,':',len(warnings),'warning(s)' if warnings else 'no warning')
            for line in warnings:
                print('   '+line.strip())


if __name__ == "__main__":
    main(sys.argv[1:])
