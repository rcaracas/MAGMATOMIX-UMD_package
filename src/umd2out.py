#!/usr/bin/env python3
###
#AUTHORS: RAZVAN CARACAS
###
#Converts the snapshots of a UMD file to other formats : extended XYZ trajectory, VASP POSCAR, LAMMPS data files and dump trajectory, DeepMD
#training set (deepmd/npy). The UMD file is read line by line (any file size), every -s-th snapshot from -i, or the single snapshot -k.

import sys, getopt, os, re
import numpy as np
import umd_processes_fast as umdpf
from averages import kinetic_included
from vibrations_umd import AMU_A2FS2

FORMATS = {'1':'xyz','2':'poscar','3':'lammps','4':'deepmd','5':'lammpstrj'}
GPA = 1.0/160.21766208                  #1 GPa in eV/A^3
SETSIZE = 5000                          #frames per set.nnn folder of the DeepMD training set (the default of dpdata)
ENERGYFIELDS = ('InternalEnergy','KineticEnergyIons','TotalElectronicEnergy','ElectronicEntropy')


#****Energies and stresses

def kinetic_in_internal(UMDname,natom):#Whether InternalEnergy contains KineticEnergyIons (averages.kinetic_included, from all the snapshots)
    series = {k:[] for k in ENERGYFIELDS}
    for isnap, fields, atoms in umdpf.read_snapshots(UMDname,natom,fields=ENERGYFIELDS,atoms=False):
        for k in ENERGYFIELDS:
            series[k].append(fields[k][0] if fields.get(k) else np.nan)
    IE, K, TEE, ES = (np.array(series[k]) for k in ENERGYFIELDS)
    if np.all(np.isnan(K)) or np.all(np.isnan(IE)):
        return False, 'no KineticEnergyIons'
    return kinetic_included(IE,np.nan_to_num(K),None if np.all(np.isnan(TEE)) else np.nan_to_num(TEE),np.nan_to_num(ES))


def free_energy(fields,included):#Potential energy of the ions minus T S_el (the Mermin free energy, of which the forces are the derivatives), eV
    energy = fields['InternalEnergy'][0]
    if included:
        energy -= fields['KineticEnergyIons'][0]
    if fields.get('ElectronicEntropy'):
        energy -= fields['ElectronicEntropy'][0]
    return energy


def static_stress(fields,atoms,masses,volume):#Stress tensor 3x3 (eV/A^3, positive in compression as the UMD Pressure) without the kinetic part
    #m v v / V, which the StressTensor of the UMD file contains (VASP 'Total+kin.') ; the m v v of the velocities is scaled so that its trace is
    #2 KineticEnergyIons ; isotropic 2/3 KineticEnergyIons when the velocities are 0 or off by more than 2 (first snapshot of the parsers).
    #None when the file has no stress.
    s = fields.get('StressTensor',[])
    if len(s) < 6 or not np.any(s[:6]):
        return None
    xx, yy, zz, xy, yz, zx = s[:6]
    sigma = np.array([[xx,xy,zx],[xy,yy,yz],[zx,yz,zz]])*GPA
    vel = atoms[:,9:12]
    mvv = np.einsum('i,ij,ik->jk',masses,vel,vel)*AMU_A2FS2
    kinetic = fields.get('KineticEnergyIons',[0.0])[0]
    if kinetic > 0.0:
        if 0.5 < np.trace(mvv)/(2.0*kinetic) < 2.0:
            mvv *= 2.0*kinetic/np.trace(mvv)
        else:
            mvv = np.eye(3)*2.0*kinetic/3.0
    return sigma - mvv/volume


#****LAMMPS

def lammps_box(cell):#LAMMPS restricted triclinic box of the same lattice : a along x, b in the xy plane, right-handed (c -> -c otherwise), tilts
    #reduced to at most half a box length (by adding lattice vectors) ; returns the box (rows A, B, C, A) and the rotation (x_lammps = x @ rotation)
    cell = np.array(cell,dtype=float)
    if np.linalg.det(cell) < 0.0:
        cell[2] = -cell[2]
    a, b, c = cell
    lx = np.linalg.norm(a)
    xy = b @ a/lx
    ly = np.sqrt(b @ b - xy**2)
    xz = c @ a/lx
    yz = (b @ c - xy*xz)/ly
    lz = np.sqrt(c @ c - xz**2 - yz**2)
    rotated = np.array([[lx,0.0,0.0],[xy,ly,0.0],[xz,yz,lz]])
    rotation = np.linalg.solve(cell,rotated)
    A, B, C = rotated.copy()
    C -= np.rint(C[1]/ly)*B
    C -= np.rint(C[0]/lx)*A
    B -= np.rint(B[0]/lx)*A
    return np.array([A,B,C]), rotation


def lammps_atoms(atoms,cell):#Box, positions wrapped in the box, image flags (from the unwrapped positions), velocities (A/ps), forces (eV/A)
    box, rotation = lammps_box(cell)
    inverse = np.linalg.inv(box)
    s = atoms[:,3:6] @ rotation @ inverse
    wrapped = np.floor(s)
    image = (wrapped + np.rint((atoms[:,6:9]-atoms[:,3:6]) @ rotation @ inverse)).astype(int)
    return box, (s-wrapped) @ box, image, atoms[:,9:12] @ rotation*1000.0, atoms[:,12:15] @ rotation


def tilts(box):
    return box[1,0], box[2,0], box[2,1]


def write_lammps_data(filename,title,box,elements,masses,types,x,image,vels):#Data file, atom_style atomic, units metal ; type k = elements[k-1]
    xy, xz, yz = tilts(box)
    with open(filename,'w') as ff:
        ff.write(title+'\n\n')
        ff.write('%d atoms\n%d atom types\n\n' % (len(types),len(elements)))
        ff.write('0.0 %.10f xlo xhi\n0.0 %.10f ylo yhi\n0.0 %.10f zlo zhi\n' % (box[0,0],box[1,1],box[2,2]))
        if max(abs(xy),abs(xz),abs(yz)) > 1e-10:
            ff.write('%.10f %.10f %.10f xy xz yz\n' % (xy,xz,yz))
        ff.write('\nMasses\n\n')
        for k,(el,m) in enumerate(zip(elements,masses)):
            ff.write('%d %s  # %s\n' % (k+1,m,el))
        ff.write('\nAtoms  # atomic\n\n')
        for i in range(len(types)):
            ff.write('%d %d %.8f %.8f %.8f %d %d %d\n' % (i+1,types[i]+1,*x[i],*image[i]))
        ff.write('\nVelocities\n\n')
        for i in range(len(types)):
            ff.write('%d %.8e %.8e %.8e\n' % (i+1,*vels[i]))


def lammps_dump_frame(step,time,box,symbols,types,x,image,vels,forces):#One frame of a dump file : id type element x y z ix iy iz vx vy vz fx fy fz
    xy, xz, yz = tilts(box)
    lines = ['ITEM: TIMESTEP',str(step),'ITEM: TIME','%.6f' % time,'ITEM: NUMBER OF ATOMS',str(len(types))]
    if max(abs(xy),abs(xz),abs(yz)) > 1e-10:
        lines.append('ITEM: BOX BOUNDS xy xz yz pp pp pp')
        lines.append('%.10f %.10f %.10f' % (min(0.0,xy,xz,xy+xz),box[0,0]+max(0.0,xy,xz,xy+xz),xy))
        lines.append('%.10f %.10f %.10f' % (min(0.0,yz),box[1,1]+max(0.0,yz),xz))
        lines.append('%.10f %.10f %.10f' % (0.0,box[2,2],yz))
    else:
        lines.append('ITEM: BOX BOUNDS pp pp pp')
        lines += ['0.0 %.10f' % box[k,k] for k in range(3)]
    lines.append('ITEM: ATOMS id type element x y z ix iy iz vx vy vz fx fy fz')
    for i in range(len(types)):
        lines.append('%d %d %s %.8f %.8f %.8f %d %d %d %.8e %.8e %.8e %.6f %.6f %.6f' % (i+1,types[i]+1,symbols[i],*x[i],*image[i],*vels[i],*forces[i]))
    return '\n'.join(lines)+'\n'


#****DeepMD

class DeepMDSet:#deepmd/npy training set : type.raw, type_map.raw, then folders set.000, set.001... of SETSIZE frames : box.npy, coord.npy (A),
    #energy.npy (eV), force.npy (eV/A, if the forces of the set are not all 0), virial.npy (eV, if the file has stresses)
    def __init__(self,folder,elements,types):
        self.folder = folder
        os.makedirs(folder,exist_ok=True)#the files of a former run are overwritten (those not written again removed at the end) ; other files are kept
        np.savetxt(os.path.join(folder,'type.raw'),types,fmt='%d')
        with open(os.path.join(folder,'type_map.raw'),'w') as ff:
            ff.write('\n'.join(elements)+'\n')
        self.nsets = 0
        self.nframes = 0
        self.novirial = 0
        self.noforce = 0
        self.frames = []

    def add(self,box,coord,energy,force,virial):
        self.frames.append((box,coord,energy,force,virial))
        if len(self.frames) == SETSIZE:
            self.flush()

    def flush(self):
        if not self.frames:
            return
        folder = os.path.join(self.folder,'set.%03d' % self.nsets)
        os.makedirs(folder,exist_ok=True)
        box, coord, energy, force, virial = zip(*self.frames)
        np.save(os.path.join(folder,'box.npy'),np.array(box).reshape(len(box),9))
        np.save(os.path.join(folder,'coord.npy'),np.array(coord).reshape(len(box),-1))
        np.save(os.path.join(folder,'energy.npy'),np.array(energy))
        if np.any(np.array(force)):                 #forces all 0 : not in the file (e.g. a LAMMPS dump without fx fy fz), not a training target
            np.save(os.path.join(folder,'force.npy'),np.array(force).reshape(len(box),-1))
        else:
            self.noforce += len(box)
        if all(v is not None for v in virial):
            np.save(os.path.join(folder,'virial.npy'),np.array(virial).reshape(len(box),9))
        else:
            self.novirial += len(box)
        self.nsets += 1
        self.nframes += len(box)
        self.frames = []


#****Main

def output_names(base,formats,ksnapshot):
    tag = '' if ksnapshot < 0 else '.'+str(ksnapshot)
    return {'xyz':base+tag+'.xyz', 'lammpstrj':base+tag+'.lammpstrj', 'deepmd':base+tag+'.deepmd',
            'poscar':base+'._%d_.POSCAR.vasp', 'lammps':base+'._%d_.lammps.data'}


def remove_former(names,formats,ksnapshot):#The outputs of a former run in the formats of this one : the trajectory files, the files per
    folder, name = os.path.split(names['poscar'][:-len('._%d_.POSCAR.vasp')])#snapshot (all of them, another -i/-s writes other snapshots ;
    paths = [names[fmt] for fmt in ('xyz','lammpstrj') if fmt in formats]    #with -k only that one, the others extracted one by one are kept),
    for fmt, suffix in (('poscar',r'\.POSCAR\.vasp'),('lammps',r'\.lammps\.data')):#and in the deepmd folder the files it writes,
        if fmt in formats:                                                          #then the folders left empty
            paths += [names[fmt] % ksnapshot] if ksnapshot >= 0 else umdpf.former_outputs(folder,re.escape(name)+r'\._\d+_'+suffix)
    folders = []
    if 'deepmd' in formats and os.path.isdir(names['deepmd']):
        sets = [os.path.join(names['deepmd'],d) for d in sorted(os.listdir(names['deepmd'])) if re.fullmatch(r'set\.\d{3}',d)]
        paths += [os.path.join(names['deepmd'],f) for f in ('type.raw','type_map.raw')]
        paths += [os.path.join(d,f+'.npy') for d in sets for f in ('box','coord','energy','force','virial')]
        folders = sets+[names['deepmd']]
    umdpf.remove_former_outputs(paths,folders)


def main(argv):
    umdpf.headerumd()
    UMDname = 'output.umd.dat'
    firststep = 0
    sampling = 1
    ksnapshot = -1
    base = None
    formats = ['xyz']
    unwrapped = False
    try:
        opts, arg = getopt.getopt(argv,"hf:i:s:k:o:r:u")
    except getopt.GetoptError:
        print('umd2out.py -f <umdfile> -i <InitialStep> -s <Sampling> -k <Snapshot> -o <OutputName> -r <formats> -u')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('umd2out.py : converts the snapshots of a UMD file to other formats')
            print('umd2out.py -f <umdfile> -i <InitialStep> -s <Sampling> -k <Snapshot> -o <OutputName> -r <formats> -u')
            print('-f : the UMD file (default output.umd.dat)')
            print('-f all : every *.umd.dat file of the current folder, one after the other (-o cannot be given then).')
            print('-i : first snapshot written (from 0 ; default 0)')
            print('-s : every s-th snapshot from -i (default 1)')
            print('-k : only snapshot k (from 0 ; replaces -i and -s)')
            print('-o : name of the outputs (default : the UMD file without .umd.dat)')
            print('-r : output formats, numbers or names separated by commas (default 1) ; the file is read once for all :')
            print('     1 = xyz       : extended XYZ trajectory <name>.xyz : species, positions (A ; wrapped, -u unwrapped) and forces (eV/A) ;')
            print('                     comment line with the cell (Lattice), the time, the temperature and the energy (as for DeepMD). Plain XYZ')
            print('                     readers read the symbol and positions ; ASE, OVITO and the ML potential codes read all.')
            print('     2 = poscar    : one VASP5 POSCAR per snapshot, <name>._<snapshot>_.POSCAR.vasp, direct coordinates, Cartesian velocities (A/fs)')
            print('     3 = lammps    : one LAMMPS data file per snapshot, <name>._<snapshot>_.lammps.data : atom_style atomic, units metal ;')
            print('                     type k = k-th element of the UMD header (listed under Masses) ; positions in the box, image flags from the')
            print('                     unwrapped positions, velocities (A/ps). Read it with  units metal, atom_style atomic, read_data <file>.')
            print('     4 = deepmd    : DeepMD training set (deepmd/npy) in the folder <name>.deepmd : type.raw, type_map.raw (the elements of the')
            print('                     UMD header), set.000, set.001... of '+str(SETSIZE)+' frames : box, coord (A), energy (eV), force (eV/A), virial (eV).')
            print('                     energy = potential energy of the ions - T S_el (InternalEnergy, minus KineticEnergyIons when it contains it,')
            print('                     minus ElectronicEntropy) : the free energy whose derivatives are the forces (VASP TOTEN).')
            print('                     virial = V x (StressTensor - m v v/V) : the stress of the UMD file contains the kinetic part (VASP Total+kin.),')
            print('                     removed with the velocities, scaled to KineticEnergyIons ; not written when the file has no stresses.')
            print('                     On-the-fly ML runs : the UMD file does not tell the DFT steps from the ML ones ; train on the DFT steps')
            print('                     only (e.g. the ML_AB file of VASP).')
            print('     5 = lammpstrj : LAMMPS dump trajectory <name>.lammpstrj : id type element x y z ix iy iz vx vy vz fx fy fz, units metal')
            print('                     (velocities A/ps, time ps), as written by dump custom (with dump_modify time yes) ; for OVITO, VMD, ASE.')
            print('     LAMMPS outputs : restricted triclinic box (a along x, b in the xy plane, tilts at most half a box length) : the')
            print('     structure is rotated, and for a left-handed cell c is replaced by -c ; the lattice and the positions are unchanged.')
            print('-u : unwrapped positions in the XYZ file')
            print('The outputs of a former run in the formats asked are removed first : the trajectory files, all the files per snapshot (with -k')
            print('only that snapshot\'s, so snapshots extracted one by one are kept), and the files umd2out.py writes in the deepmd folder.')
            sys.exit()
        elif opt == '-f':
            UMDname = str(arg)
        elif opt == '-i':
            firststep = int(arg)
        elif opt == '-s':
            sampling = int(arg)
        elif opt == '-k':
            ksnapshot = int(arg)
        elif opt == '-o':
            base = str(arg)
        elif opt == '-r':
            formats = [FORMATS.get(x.strip(),x.strip().lower()) for x in arg.split(',')]
        elif opt == '-u':
            unwrapped = True
    if not os.path.isfile(UMDname):
        umdpf.missing_umd(UMDname)
        sys.exit()
    for fmt in formats:
        if fmt not in FORMATS.values():
            print('Unknown output format',fmt,': use',', '.join(k+' = '+v for k,v in FORMATS.items()))
            sys.exit()
    if sampling < 1 or firststep < 0:
        print('-s must be at least 1 and -i at least 0')
        sys.exit()
    if base is None:
        base = UMDname[:-8] if UMDname.endswith('.umd.dat') else os.path.splitext(UMDname)[0]
    if ksnapshot >= 0:
        firststep, sampling = ksnapshot, 1

    MyCrystal, TimeStep = umdpf.Crystallization(UMDname)
    natom = MyCrystal.natom
    elements = list(MyCrystal.elements)
    types = np.array(MyCrystal.typat[:natom],dtype=int)
    symbols = [elements[t] for t in types]
    masses = np.array([MyCrystal.masses[t] for t in types],dtype=float)
    order = np.argsort(types,kind='stable')                     #POSCAR : atoms grouped by element
    names = output_names(base,formats,ksnapshot)

    included = False
    if 'xyz' in formats or 'deepmd' in formats:
        included, how = kinetic_in_internal(UMDname,natom)
        print('energy : InternalEnergy','- KineticEnergyIons' if included else '','- ElectronicEntropy  (' + how + ')')
    remove_former(names,formats,ksnapshot)
    if 'deepmd' in formats and os.path.exists(names['deepmd']) and not os.path.isdir(names['deepmd']):
        print(names['deepmd'],'is a file, not a folder : remove or rename it, or choose another name with -o')
        sys.exit()
    xyzfile = dumpfile = deepmd = None
    nwritten = 0
    anyforce = False
    isnap = -1
    lefthanded = False
    for isnap, fields, atoms in umdpf.read_snapshots(UMDname,natom,firststep,sampling):
        if atoms is None:
            if ksnapshot >= 0 and isnap > ksnapshot:
                break
            continue
        if nwritten == 0:                                       #outputs opened at the first snapshot selected
            xyzfile = open(names['xyz'],'w') if 'xyz' in formats else None
            dumpfile = open(names['lammpstrj'],'w') if 'lammpstrj' in formats else None
            deepmd = DeepMDSet(names['deepmd'],elements,types) if 'deepmd' in formats else None
        cell = umdpf.cell_of(fields)
        if cell is None:
            cell = np.array(MyCrystal.rprimd,dtype=float)
        lefthanded = lefthanded or np.linalg.det(cell) < 0.0
        anyforce = anyforce or bool(np.any(atoms[:,12:15]))
        time = fields.get('time',[isnap*TimeStep])[0]
        if xyzfile or deepmd:
            energy = free_energy(fields,included) if fields.get('InternalEnergy') else None
            stress = static_stress(fields,atoms,masses,abs(np.linalg.det(cell)))
        if xyzfile:
            info = 'time=%.4f' % time
            if fields.get('Temperature'):
                info += ' temperature=%.2f' % fields['Temperature'][0]
            if energy is not None:
                info += ' energy=%.6f' % energy
            if stress is not None:                              #ASE convention : negative in compression
                info += ' stress="' + ' '.join('%.8f' % x for x in -stress.ravel()) + '"'
            info += ' snapshot=%d source="%s"' % (isnap,os.path.basename(UMDname))
            xyzfile.write(umdpf.xyz_frame(symbols,atoms[:,6:9] if unwrapped else atoms[:,3:6],cell,info,[('forces',atoms[:,12:15])]))
        if 'poscar' in formats:
            umdpf.write_poscar(names['poscar'] % isnap,UMDname+' snapshot '+str(isnap),cell,elements,MyCrystal.types,
                               atoms[order,0:3] % 1.0,atoms[order,9:12])
        if 'lammps' in formats or dumpfile:
            box, x, image, vels, forces = lammps_atoms(atoms,cell)
            if 'lammps' in formats:
                write_lammps_data(names['lammps'] % isnap,'LAMMPS data file from '+UMDname+' snapshot '+str(isnap)+' (umd2out.py) ; types : '
                                  +' '.join(str(k+1)+'='+el for k,el in enumerate(elements)),box,elements,MyCrystal.masses,types,x,image,vels)
            if dumpfile:
                dumpfile.write(lammps_dump_frame(isnap,time/1000.0,box,symbols,types,x,image,vels,forces))
        if deepmd:
            if energy is None:
                print('No InternalEnergy in snapshot',isnap,': no DeepMD set can be written')
                sys.exit()
            deepmd.add(cell,atoms[:,3:6],energy,atoms[:,12:15],None if stress is None else stress*abs(np.linalg.det(cell)))
        nwritten += 1
        if ksnapshot >= 0:
            break
    if nwritten == 0:
        print('No snapshot selected : the file has',isnap+1,'snapshots')
        sys.exit()
    print(nwritten,'snapshots written')
    if not anyforce and (xyzfile or dumpfile or deepmd):
        print('WARNING : all the forces of the file are 0 (not written by the MD code, e.g. a LAMMPS dump without fx fy fz) :'
              +(' no force.npy in the DeepMD set' if deepmd else '')+(' ;' if deepmd and (xyzfile or dumpfile) else '')
              +(' the forces of the xyz/lammpstrj files are 0' if xyzfile or dumpfile else ''))
    if lefthanded and ('lammps' in formats or dumpfile):
        print('left-handed cell : c replaced by -c in the LAMMPS files')
    for fmt in formats:
        if fmt in ('poscar','lammps'):
            print(fmt,':',nwritten,'files',names[fmt].replace('%d','<snapshot>'))
        elif fmt != 'deepmd':
            print(fmt,':',names[fmt])
    if xyzfile:
        xyzfile.close()
    if dumpfile:
        dumpfile.close()
    if deepmd:
        deepmd.flush()
        notes = (['(no virial : the file has no stresses)'] if deepmd.novirial else [])+(['(no force.npy in '+str(deepmd.noforce)+' frames : forces 0)'] if deepmd.noforce else [])
        print(' '.join(['deepmd :',names['deepmd'],':',str(deepmd.nframes),'frames in',str(deepmd.nsets),'sets']+notes))


if __name__ == "__main__":
    if not umdpf.run_on_all(sys.argv[1:],main,single=('-o',)):#-f all : every UMD file of the folder
        main(sys.argv[1:])
