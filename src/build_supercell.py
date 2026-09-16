#!/usr/bin/env python3
###
##AUTHORS: RAZVAN CARACAS
###

import sys,getopt,os.path
import numpy as np
import crystallography as cr
import umd_processes_fast as umdpf


def ReadPOSCAR(PoscarFile):
    #reads a VASP POSCAR/CONTCAR file (VASP5 format, with an element-symbols line),
    #such as the ones exported by VESTA
    with open(PoscarFile,'r') as ff:
        lines = ff.readlines()
    comment = lines[0].strip()
    scale = float(lines[1].split()[0])
    if scale < 0.0:
        print('Negative (volume-based) scaling factors in POSCAR files are not supported.')
        sys.exit()
    cellvecs = np.zeros((3,3))
    for ii in range(3):
        cellvecs[ii,:] = [float(xx) for xx in lines[2+ii].split()[0:3]]
    cellvecs = cellvecs * scale
    elemline = lines[5].split()
    if all(xx.lstrip('-').isdigit() for xx in elemline):
        print('The POSCAR file ',PoscarFile,' has no element-symbols line (VASP4 format). Please use a VASP5-style POSCAR, such as the ones exported by VESTA.')
        sys.exit()
    elements = elemline
    counts = [int(xx) for xx in lines[6].split()]
    natom = sum(counts)
    lineno = 7
    if lines[lineno].strip()[0] in ('s','S'):     #optional Selective dynamics line
        lineno += 1
    coordline = lines[lineno].strip()
    lineno += 1
    if coordline[0] in ('d','D'):
        coordtype = 'reduced'
    elif coordline[0] in ('c','C','k','K'):
        coordtype = 'cartesian'
    else:
        print('Could not understand the coordinate-type line in ',PoscarFile,' : ',coordline)
        sys.exit()
    atoms = [cr.Atom() for _ in range(natom)]
    iatom = 0
    for ielem in range(len(elements)):
        for jatom in range(counts[ielem]):
            entry = lines[lineno].split()
            lineno += 1
            atoms[iatom].symbol = elements[ielem]
            coords = [float(xx) for xx in entry[0:3]]
            if coordtype == 'reduced':
                atoms[iatom].xred = coords
                atoms[iatom].xcart = list(np.array(coords) @ cellvecs)
            else:
                atoms[iatom].xcart = coords
                atoms[iatom].xred = list(np.array(coords) @ np.linalg.inv(cellvecs))
            iatom += 1
    print('Read ',natom,' atoms from ',PoscarFile,' with elements ',elements,' and counts ',counts)
    return(comment,cellvecs,elements,counts,natom,atoms)


def BuildSupercell(cellvecs,elements,counts,natom,atoms,nx,ny,nz):
    #replicates the unit cell nx * ny * nz times, working in reduced coordinates
    newcellvecs = np.zeros((3,3))
    newcellvecs[0,:] = cellvecs[0,:] * nx
    newcellvecs[1,:] = cellvecs[1,:] * ny
    newcellvecs[2,:] = cellvecs[2,:] * nz
    newcounts = [cc * nx * ny * nz for cc in counts]
    newnatom = natom * nx * ny * nz
    newatoms = [cr.Atom() for _ in range(newnatom)]
    jatom = 0
    for iatom in range(natom):
        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    newatoms[jatom].symbol = atoms[iatom].symbol
                    newatoms[jatom].xred = [(atoms[iatom].xred[0]+ix)/nx, (atoms[iatom].xred[1]+iy)/ny, (atoms[iatom].xred[2]+iz)/nz]
                    newatoms[jatom].xcart = list(np.array(newatoms[jatom].xred) @ newcellvecs)
                    jatom += 1
    return(newcellvecs,newcounts,newnatom,newatoms)


def WritePOSCAR(OutputName,comment,newcellvecs,elements,newcounts,newnatom,newatoms):
    filename = OutputName + '.vasp'
    with open(filename,'w') as ff:
        ff.write(comment + '\n')
        ff.write('  1.0\n')
        for ii in range(3):
            ff.write('   ' + str(newcellvecs[ii,0]) + '  ' + str(newcellvecs[ii,1]) + '  ' + str(newcellvecs[ii,2]) + '\n')
        ff.write('  ' + '  '.join(elements) + '\n')
        ff.write('  ' + '  '.join(str(cc) for cc in newcounts) + '\n')
        ff.write('Direct\n')
        for iatom in range(newnatom):
            ff.write('  ' + str(newatoms[iatom].xred[0]) + '  ' + str(newatoms[iatom].xred[1]) + '  ' + str(newatoms[iatom].xred[2]) + '\n')
    print('done writing ',filename)


def WriteXYZ(OutputName,comment,newnatom,newatoms):
    filename = OutputName + '.xyz'
    with open(filename,'w') as ff:
        ff.write(str(newnatom) + '\n')
        ff.write(comment + '\n')
        for iatom in range(newnatom):
            ff.write(newatoms[iatom].symbol + '  ' + str(newatoms[iatom].xcart[0]) + '  ' + str(newatoms[iatom].xcart[1]) + '  ' + str(newatoms[iatom].xcart[2]) + '\n')
    print('done writing ',filename)


def main(argv):
    umdpf.headerumd()
    PoscarFile = 'POSCAR'
    OutputName = ''
    nx = 1
    ny = 1
    nz = 1
    try:
        opts, arg = getopt.getopt(argv,"hf:x:y:z:o:",["fPOSCARfile","xExpansion","yExpansion","zExpansion","oOutputName"])
    except getopt.GetoptError:
        print ('build_supercell.py -f <POSCAR_filename> -x <nx> -y <ny> -z <nz> -o <output_name>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('build_supercell.py program to build a supercell from a unit cell read from a VASP POSCAR/CONTCAR file (such as one exported by VESTA)')
            print ('build_supercell.py -f <POSCAR_filename> -x <nx> -y <ny> -z <nz> -o <output_name>')
            print (' -f the input POSCAR file, in reduced or cartesian coordinates. Default: POSCAR')
            print (' -x -y -z the number of repetitions of the unit cell along a, b and c. Default: 1 1 1')
            print (' -o the name used as a prefix of the output .vasp and .xyz files. Default: derived from the input filename')
            sys.exit()
        elif opt in ("-f","--fPOSCARfile"):
            PoscarFile = str(arg)
        elif opt in ("-x","--xExpansion"):
            nx = int(arg)
        elif opt in ("-y","--yExpansion"):
            ny = int(arg)
        elif opt in ("-z","--zExpansion"):
            nz = int(arg)
        elif opt in ("-o","--oOutputName"):
            OutputName = str(arg)
    if not os.path.isfile(PoscarFile):
        print ('the POSCAR file ',PoscarFile,' does not exist')
        sys.exit()
    if nx < 1 or ny < 1 or nz < 1:
        print ('the supercell repetitions -x -y -z must all be >= 1')
        sys.exit()
    if OutputName == '':
        OutputName = os.path.splitext(os.path.basename(PoscarFile))[0] + '_' + str(nx) + 'x' + str(ny) + 'x' + str(nz)

    (comment,cellvecs,elements,counts,natom,atoms) = ReadPOSCAR(PoscarFile)
    (newcellvecs,newcounts,newnatom,newatoms) = BuildSupercell(cellvecs,elements,counts,natom,atoms,nx,ny,nz)
    print('Built a ',nx,'x',ny,'x',nz,' supercell with ',newnatom,' atoms')
    WritePOSCAR(OutputName,comment,newcellvecs,elements,newcounts,newnatom,newatoms)
    WriteXYZ(OutputName,comment,newnatom,newatoms)


if __name__ == "__main__":
    main(sys.argv[1:])
