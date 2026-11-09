#!/usr/bin/env python3
###
#AUTHORS: RAZVAN CARACAS
###
#Validation suite of the UMD package : runs the scripts of src/ on the reference files of data/ (cases.py), each case in a fresh folder
#of its own with links to its inputs, and compares the results with the reference of the case (reference/<case>.json) :
#   - key numbers (probes.py reads them from the output files : density, P, T, coordination numbers, D, viscosity...) within a tolerance,
#   - a checksum of every output file (dates and the path of the run folder left out),
#   - the time of every command (for information : the 3.0 parallelization).
#A case is CRASH when a command ends with an error code or prints a Python traceback, FAIL when a key number is out of its tolerance or
#missing, CHANGED when only the checksums of some outputs changed (the numbers within tolerance : e.g. the last digit), OK when every
#output is identical. data/ is never written to : it is checked before and after.

import sys, os, getopt, json, time, hashlib, re, subprocess, zipfile, platform, fnmatch, shutil

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SRC = os.path.join(ROOT,'src')
REFERENCE = os.path.join(HERE,'reference')
sys.path.insert(0,HERE)
import cases

VOLATILE = ('*.pdf','*.png','*.log','umd_workflow.done','.DS_Store')        #outputs not compared (figures, logs) ; still listed
DATE = re.compile(rb'\d{4}-\d{2}-\d{2}[ T_]\d{2}:\d{2}(:\d{2})?')


def snapshot_tree(top):#size and modification time of every file below top (to check that data/ is left as it was)
    state = {}
    for folder, dirs, files in os.walk(top):
        for f in files:
            path = os.path.join(folder,f)
            st = os.lstat(path)
            state[os.path.relpath(path,top)] = (st.st_size,st.st_mtime)
    return state


def checksum(path,folder):#sha256 of a file ; for a text file without the dates and the path of the run folder
    with open(path,'rb') as f:
        head = f.read(8192)
    h = hashlib.sha256()
    if b'\0' in head:#binary (numpy arrays of the DeepMD sets)
        with open(path,'rb') as f:
            for block in iter(lambda: f.read(1 << 20),b''):
                h.update(block)
        return h.hexdigest()
    prefix = os.path.abspath(folder).encode()
    with open(path,'rb') as f:
        for line in f:
            h.update(DATE.sub(b'<date>',line.replace(prefix,b'.')))
    return h.hexdigest()


def prepare(case,folder,data):#the case folder : links to the inputs of data/, the members of zip files extracted
    os.makedirs(folder)
    inputs = []
    for item in case.inputs:
        if item.endswith('.zip'):
            with zipfile.ZipFile(os.path.join(data,item)) as z:
                for name in z.namelist():
                    if name.endswith('/') or '/' in name or '\\' in name or name.startswith('.'):
                        continue                                                #flat archives only : nothing outside the case folder
                    with z.open(name) as src, open(os.path.join(folder,name),'wb') as dst:
                        shutil.copyfileobj(src,dst)
                    inputs.append(name)
        else:
            source = os.path.join(data,item)
            if not os.path.exists(source):
                raise FileNotFoundError(source)
            os.symlink(source,os.path.join(folder,os.path.basename(item)))
            inputs.append(os.path.basename(item))
    return inputs


def run_case(case,folder,data):
    result = {'commands':[],'numbers':{},'files':{},'crash':[],'truth':[]}
    try:
        inputs = prepare(case,folder,data)
    except (OSError,zipfile.BadZipFile) as error:
        result['crash'].append('input : '+str(error))
        return result
    env = dict(os.environ,MPLBACKEND='Agg')
    outs = []
    for k, (script, args) in enumerate(case.commands):
        label = ' '.join([script]+args)
        out = os.path.join(folder,'%02d.%s.out' % (k,script.lstrip('@').replace('.py','')))
        start = time.time()
        with open(out,'w') as f:
            path = os.path.join(HERE,script[1:]) if script.startswith('@') else os.path.join(SRC,script)   #@ : a script of validation/
            proc = subprocess.run([sys.executable,path]+args,cwd=folder,stdout=f,stderr=subprocess.STDOUT,
                                  stdin=subprocess.DEVNULL,env=env)
        elapsed = time.time()-start
        with open(out,errors='replace') as f:
            text = f.read()
        outs.append(text)
        result['commands'].append({'command':label,'code':proc.returncode,'time':round(elapsed,2)})
        if proc.returncode != 0:
            result['crash'].append('%s : exit code %d' % (label,proc.returncode))
        if 'Traceback (most recent call last)' in text:
            result['crash'].append('%s : traceback in %s' % (label,os.path.basename(out)))
    for log in sorted(f for f in os.listdir(folder) if f.endswith('.log')):#logs of the workflow (scripts run by it)
        with open(os.path.join(folder,log),errors='replace') as f:
            n = f.read().count('Traceback (most recent call last)')
        if n:
            result['crash'].append('%d traceback(s) in %s' % (n,log))
    for probe in case.probes:
        try:
            for name, value in probe(folder,outs).items():
                result['numbers'][name] = value
        except Exception as error:
            result['numbers']['probe '+getattr(probe,'__name__','?')] = 'ERROR '+type(error).__name__+': '+str(error)
    for probe in case.truth:#comparisons with the exact answers of synthetic files : (label, value, expected, tolerance)
        try:
            for label, value, expected, tol in probe(folder,outs):
                result['numbers'][label] = value
                result['truth'].append((label,value,expected,tol))
        except Exception as error:
            result['truth'].append(('probe '+getattr(probe,'__name__','?'),'ERROR '+type(error).__name__+': '+str(error),'a value',0))
    skip = set(inputs)|{os.path.basename(o) for o in os.listdir(folder) if re.match(r'\d\d\..*\.out$',o)}
    for top, dirs, files in os.walk(folder):
        for f in files:
            rel = os.path.relpath(os.path.join(top,f),folder)
            if rel in skip or os.path.islink(os.path.join(top,f)):
                continue
            if any(fnmatch.fnmatch(f,p) for p in VOLATILE):
                result['files'][rel] = 'not compared'
            else:
                result['files'][rel] = checksum(os.path.join(top,f),folder)
    return result


def value_of(entry):#a number of a probe : value, or (value, relative tolerance)
    return entry if not isinstance(entry,(list,tuple)) else entry[0]


def tolerance_of(entry):
    return cases.RTOL if not isinstance(entry,(list,tuple)) else entry[1]


def compare(result,reference):#status and the lines of the report
    lines = []
    if result['crash']:
        lines += ['   crash : '+c for c in result['crash']]
    wrong = []
    for label, value, expected, tol in result['truth']:
        if isinstance(value,str) or isinstance(expected,str):
            ok = value == expected
        else:
            ok = isinstance(value,(int,float)) and isinstance(expected,(int,float)) and abs(value-expected) <= tol
        if not ok:
            wrong.append(label)
            lines.append('   FAIL  %s = %s, exact %s (tolerance %s)' % (label,value,expected,tol))
    if result['truth']:
        lines.append('   %d of %d results equal to the exact ones within their tolerance' % (len(result['truth'])-len(wrong),len(result['truth'])))
    if reference is None:
        return ('CRASH' if result['crash'] else 'FAIL' if wrong else 'NEW'), lines
    failed = bool(wrong)
    refnumbers = reference.get('numbers',{})
    for name, entry in result['numbers'].items():
        value = value_of(entry)
        if name not in refnumbers:
            lines.append('   new number  %s = %s' % (name,value))
            continue
        ref = value_of(refnumbers[name])
        rtol = tolerance_of(entry)
        if isinstance(value,(int,float)) and isinstance(ref,(int,float)) and not isinstance(value,bool):
            same = value == ref or (value != value and ref != ref) or (value == value and ref == ref and abs(value-ref) <= rtol*max(abs(ref),abs(value))+cases.ATOL)
        else:
            same = value == ref
        if not same:
            failed = True
            dev = ''
            if isinstance(value,(int,float)) and isinstance(ref,(int,float)) and ref:
                dev = '   (%+.2e relative, tolerance %.0e)' % ((value-ref)/abs(ref),rtol)
            lines.append('   FAIL  %s = %s, reference %s%s' % (name,value,ref,dev))
    for name in refnumbers:
        if name not in result['numbers']:
            failed = True
            lines.append('   FAIL  %s missing (reference %s)' % (name,value_of(refnumbers[name])))
    reffiles = reference.get('files',{})
    changed = [f for f in result['files'] if f in reffiles and result['files'][f] != reffiles[f]]
    new = [f for f in result['files'] if f not in reffiles]
    gone = [f for f in reffiles if f not in result['files']]
    for label, group in (('changed',changed),('new',new),('not written',gone)):
        if group:
            lines.append('   %s file(s) : %s' % (label,', '.join(sorted(group))))
    reftimes = {c['command']:c['time'] for c in reference.get('commands',[])}
    for c in result['commands']:
        if c['command'] in reftimes and reftimes[c['command']] > 5 and abs(c['time']/reftimes[c['command']]-1) > 0.3:
            lines.append('   time  %-60s %7.1f s (reference %.1f s)' % (c['command'][:60],c['time'],reftimes[c['command']]))
    if result['crash']:
        return 'CRASH', lines
    if failed or gone:
        return 'FAIL', lines
    if changed or new:
        return 'CHANGED', lines
    return 'OK', lines


def environment():
    import numpy, scipy
    try:
        commit = subprocess.run(['git','-C',ROOT,'describe','--tags','--always','--dirty'],capture_output=True,text=True).stdout.strip()
    except OSError:
        commit = '?'
    return {'commit':commit,'python':platform.python_version(),'numpy':numpy.__version__,'scipy':scipy.__version__,
            'platform':platform.platform(),'cpus':os.cpu_count()}


def main(argv):
    data = os.path.join(ROOT,'data')
    outdir = None
    selected = None
    full = False
    update = False
    usage = 'run_validation.py -d <data folder> -o <run folder> -c <cases> --full --update -l'
    try:
        opts, arg = getopt.getopt(argv,'hd:o:c:l',['full','update'])
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('Validation suite of the UMD package : runs the scripts on the reference files of data/ and compares the results with')
            print('validation/reference/<case>.json (key numbers within their tolerance, checksums of the outputs, times).')
            print(usage)
            print('-d : folder of the reference files. Default : data/ of the package.')
            print('-o : run folder (must not exist), one sub-folder per case. Default : validation/runs/<date-time>/.')
            print('-c : cases to run, names or patterns separated by commas (e.g. workflow_*,parse_abinit). Default : all but the --full ones.')
            print('--full : also the cases on the large files (hydrous andesite, 468 MB ; the old corundum file).')
            print('--update : writes the results of this run as the reference of each case run (not for a case that crashed). Use it when a')
            print('     change of results is understood and wanted ; the reference files are committed with the change.')
            print('-l : lists the cases.')
            print('Status of a case : OK (every output identical), CHANGED (some outputs differ, every key number within its tolerance),')
            print('FAIL (a key number out of tolerance or missing, or an output not written), CRASH (error code or Python traceback),')
            print('NEW (no reference yet).')
            sys.exit()
        elif opt == '-d':
            data = os.path.abspath(arg)
        elif opt == '-o':
            outdir = os.path.abspath(arg)
        elif opt == '-c':
            selected = [x.strip() for x in arg.split(',') if x.strip()]
        elif opt == '--full':
            full = True
        elif opt == '--update':
            update = True
        elif opt == '-l':
            for case in cases.CASES:
                print('%-28s %s%s' % (case.name,'[--full] ' if case.full else '',case.about))
            sys.exit()
    todo = [c for c in cases.CASES if (selected is None and (full or not c.full)) or
            (selected is not None and any(fnmatch.fnmatch(c.name,p) for p in selected))]
    if not todo:
        print('No case matches',selected,'; -l lists them')
        sys.exit(1)
    missing = [lib for lib in ('c_UMDprocess','c_bonds_full','c_clusters','c_gofr')
               if not any(os.path.isfile(os.path.join(SRC,lib+ext)) for ext in ('.dylib','.so','.dll'))]
    if missing:
        print('Compile the C libraries first (README.txt) : missing',', '.join(missing))
        sys.exit(1)
    if outdir is None:
        outdir = os.path.join(HERE,'runs',time.strftime('%Y%m%d-%H%M%S'))
    if os.path.exists(outdir):
        print('The run folder',outdir,'exists : give a new one')
        sys.exit(1)
    if os.path.abspath(outdir).startswith(os.path.abspath(data)+os.sep):
        print('The run folder cannot be inside the data folder')
        sys.exit(1)
    os.makedirs(outdir)
    env = environment()
    print('UMD validation :',len(todo),'case(s), data',data,', run folder',outdir)
    print('   '+', '.join(k+' '+str(v) for k,v in env.items()))
    before = snapshot_tree(data)
    report = []
    statuses = {}
    start = time.time()
    for case in todo:
        t0 = time.time()
        print('%-28s' % case.name,end=' ',flush=True)
        result = run_case(case,os.path.join(outdir,case.name),data)
        path = os.path.join(REFERENCE,case.name+'.json')
        reference = None
        if os.path.isfile(path):
            with open(path) as f:
                reference = json.load(f)
        status, lines = compare(result,reference)
        statuses[case.name] = status
        head = '%-28s %-8s %7.1f s' % (case.name,status,time.time()-t0)
        print(head[29:])
        for line in lines:
            print(line)
        report += [head]+lines
        if update and not result['crash']:
            os.makedirs(REFERENCE,exist_ok=True)
            with open(path,'w') as f:
                json.dump({'about':case.about,'environment':env,'commands':result['commands'],'numbers':result['numbers'],
                           'files':result['files']},f,indent=1,sort_keys=False)
                f.write('\n')
        elif update:
            print('   reference not updated : the case crashed')
    after = snapshot_tree(data)
    touched = sorted(f for f in set(before)|set(after) if before.get(f) != after.get(f))
    tail = ['','%d case(s) in %.0f s : ' % (len(todo),time.time()-start)+', '.join('%d %s' % (list(statuses.values()).count(s),s)
            for s in ('OK','CHANGED','FAIL','CRASH','NEW') if s in statuses.values())]
    if touched:
        tail.append('FAIL : files of the data folder changed during the run : '+', '.join(touched[:20]))
    else:
        tail.append('data folder unchanged')
    if update:
        tail.append('references written in '+REFERENCE)
    for line in tail:
        print(line)
    with open(os.path.join(outdir,'report.txt'),'w') as f:
        f.write('UMD validation '+time.strftime('%Y-%m-%d %H:%M')+'\n'+json.dumps(env)+'\n\n'+'\n'.join(report+tail)+'\n')
    sys.exit(0 if all(s in ('OK','CHANGED','NEW') for s in statuses.values()) and not touched else 1)


if __name__ == '__main__':
    main(sys.argv[1:])
