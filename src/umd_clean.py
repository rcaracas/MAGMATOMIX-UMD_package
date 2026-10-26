#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### Author : Razvan Caracas
### Removes the outputs of the UMD package from a folder, to start the analysis afresh : keeps the UMD files, the OUTCARs, the dumps and
### every file that is not an output of the package (umd_processes_fast.umd_outputs). By default only lists them ; -y moves them to a
### folder umd_cleaned_<date>, from which they can be put back ; --delete removes them for good.


import sys, getopt, os, shutil, time
import umd_processes_fast as umdpf


def remove_empty_folders(folder,top):#the folders left empty by the moves (sets and deepmd folders), up to top excluded
    while folder and os.path.abspath(folder) != os.path.abspath(top) and os.path.isdir(folder) and not os.listdir(folder):
        os.rmdir(folder)
        folder = os.path.dirname(folder)


def main(argv):
    umdpf.headerumd()
    recursive = False
    act = False
    delete = False
    usage = 'umd_clean.py [-r] [-y] [--delete]'
    try:
        opts, arg = getopt.getopt(argv,"hry",["delete"])
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('Removes the outputs of the scripts of the UMD package from the current folder, to start the analysis afresh. Usage:')
            print(usage)
            print('By default it only lists the files it would remove, and changes nothing.')
            print('-y : moves them to a new folder umd_cleaned_<date>_<time> of each folder cleaned ; put them back from there if needed, and')
            print('     delete that folder yourself once you are sure.')
            print('--delete : with -y, deletes them for good instead.')
            print('-r : also the folders below the current one (each one gets its own umd_cleaned_ folder ; those are never entered).')
            print('Removed : the files whose whole name is that of an output of a script for a UMD file <base>.umd.dat of the same folder')
            print('(<base>.gofr.dat, .maximum_bonds.inp, .bonding.dat, .r1.popul.dat, .msd.dat, .vibr.dat, ..., the files umd2out.py writes in')
            print('<base>.deepmd), and the files written once per folder (gofrs_summary.txt, diffusivities.txt, fullthermo.dat, thermo_integration.dat,')
            print('umd_workflow.log/.done/.summary.dat, ...) ; the list is umd_processes_fast.RUN_OUTPUTS and FOLDER_OUTPUTS.')
            print('Kept : all the UMD files (also those of build_umd.py, --extract and repair_umd.py, and the POSCAR and xyz of the builds), the')
            print('OUTCARs, LAMMPS dumps and logs, and every other file (notes, scripts, figures). The outputs of a UMD file no longer in the folder')
            print('are not recognized, hence kept. A cutoff file edited by hand keeps the name of an output (<base>.maximum_bonds.inp) : it is')
            print('listed and moved with the others.')
            sys.exit()
        elif opt == '-r':
            recursive = True
        elif opt == '-y':
            act = True
        elif opt == '--delete':
            delete = True
    if delete and not act:
        print('--delete needs -y : nothing is changed without it')
    stamp = time.strftime('%Y-%m-%d_%H%M%S')
    folders = []
    for dirpath, dirnames, filenames in os.walk(os.curdir):
        dirnames[:] = sorted(d for d in dirnames if not d.startswith('umd_cleaned_'))
        folders.append(dirpath)
        if not recursive:
            break
    total = 0
    for folder in folders:
        if not os.path.isdir(folder):#a set folder of deepmd, emptied and removed by a move above
            continue
        outputs = umdpf.umd_outputs(folder)
        if not outputs:
            continue
        total += len(outputs)
        umds = [n for n in sorted(os.listdir(folder)) if n.endswith('.umd.dat')]
        print('\n'+folder+' : '+str(len(outputs))+' output files of '+str(len(umds))+' UMD files ('+', '.join(umds)+')')
        for rel in outputs:
            print('   '+rel)
        if not act:
            continue
        target = os.path.join(folder,'umd_cleaned_'+stamp)
        for rel in outputs:
            path = os.path.join(folder,rel)
            if path.endswith('.umd.dat'):#never : umd_outputs does not return UMD files
                continue
            if delete:
                os.remove(path)
            else:
                os.makedirs(os.path.dirname(os.path.join(target,rel)),exist_ok=True)
                shutil.move(path,os.path.join(target,rel))
            remove_empty_folders(os.path.dirname(path),folder)
        print('-> '+('deleted' if delete else 'moved to '+target))
    if total == 0:
        print('\nNo output of the UMD package found'+(' here or below' if recursive else ' here (-r : also the folders below)'))
    elif not act:
        print('\n'+str(total)+' files listed, nothing changed. -y moves them to umd_cleaned_<date>_<time> folders, -y --delete deletes them.')


if __name__ == "__main__":
    main(sys.argv[1:])
