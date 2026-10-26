#!/usr/bin/env python3
###
##AUTHORS: RAZVAN CARACAS
###

import numpy as np

class Atom(object):
    def __init__(self, name = 'XXX', symbol='X', znucl = 1, mass=0.0, magnet=0.0, charge=0.0, xred = None, xcart = None, absxcart = None, vels = None, forces = None, diffcoords = None):
        self.name = name
        self.symbol = symbol
        self.znucl = znucl
        self.mass = mass
        self.magnet = magnet
        self.charge = charge
        if xred is None:
            self.xred = [0.0, 0.0, 0.0]
        else:
            self.xred = xred
        if xcart is None:
            self.xcart = [0.0, 0.0, 0.0]
        else:
            self.xcart = xcart
        if absxcart is None:
            self.absxcart = [0.0, 0.0, 0.0]
        else:
            self.absxcart = absxcart
        if vels is None:
            self.vels = [0.0, 0.0, 0.0]
        else:
            self.vels = vels
        if forces is None:
            self.forces = [0.0, 0.0, 0.0]
        else:
            self.forces = forces
        if diffcoords is None:
            self.diffcoords=[0.0, 0.0, 0.0]
        else:
            self.diffcoords = diffcoords    


class Lattice(object):#Remplacer les vecteurs de vecteurs par des matrices ?
    def __init__(self, name='crystal', acell=None, angles=None,
                 rprim=None, rprimd=None,
                 gprimd=None, stress=None,
                 typat=[],elements=[],masses=[],zelec=[],ntypat=1, natom=1, density=0.0, atoms=[],cellvolume=1.0,enthalpy=0.0,toten=0.0,internalenergy=0.0,electronicentropy = 0.0,kineticenergy=0.0,energywithdrift=0.0,gibbsfreeenergy=0.0,cv=0.0,pressure=0.0,temperature=0.0,magnetization=0.0,elecgap=0.0,noelectrons=0.0,lambda_ThermoInt=1.0):
        #acell/angles/rprim/rprimd/gprimd/stress default to None (rather than a mutable list literal) and are
        #built fresh below, since a mutable default argument is created once and shared by every instance that
        #doesn't override it: in-place edits like MyCrystal.acell[0]=... would otherwise leak into every other
        #Lattice() created afterwards in the same run.
        self.name = name
        self.acell = acell if acell is not None else [0.0 for x in range(3)]
        self.rprim = rprim if rprim is not None else [[0.0 for x in range(3)] for y in range(3)]
        self.rprimd = rprimd if rprimd is not None else [[0.0 for x in range(3)] for y in range(3)]
        self.gprimd = gprimd if gprimd is not None else [[0.0 for x in range(3)] for y in range(3)]
        self.natom = natom                                      #number of atoms
        self.ntypat = ntypat                                    #number of atom types
        self.density = density
        self.typat = [0 for x in range(self.natom)]             #gives the type of each atom
        self.types = [0 for x in range(self.ntypat)]            #gives the no. of atoms of each type
        self.elements = ['X' for x in range(self.ntypat)]       #gives the Chemical Symbol for each atom
        self.masses = ['0.0' for x in range(self.ntypat)]       #gives the atomic mass of each atom
        self.zelec = ['0.0' for x in range(self.ntypat)]        #gives the atomic mass of each atom
        self.atoms = [Atom() for x in range(self.natom)]        #gives the complete Atom() class for each atom
        self.angles = angles if angles is not None else [0.0 for x in range(3)]
        self.cellvolume = cellvolume
        self.enthalpy = enthalpy
        self.internalenergy = internalenergy
        self.electronicentropy = electronicentropy
        self.kineticenergy = kineticenergy
        self.energywithdrift = energywithdrift
        self.gibbsfreeenergy = gibbsfreeenergy
        self.cv = cv
        self.pressure = pressure
        self.elecgap = elecgap
        self.stress = stress if stress is not None else [0.0 for x in range(6)]
        self.toten=toten
        self.temperature=temperature
        self.magnetization=magnetization
        self.elecgap=elecgap
        self.noelectrons=noelectrons
        self.lambda_ThermoInt=lambda_ThermoInt
        self.origin=''                          #the code that produced the run (origin line of the UMD header)
        self.firstsnapshot=0                    #first snapshot of the UMD file used for a bonding file (its first_snapshot line)
    #Coordinates : xcart = xred @ rprimd (rows of rprimd = the vectors a, b, c), xred = xcart @ inverse(rprimd) ; the cell is rprimd as it is
    def xred2xcart(self,red):
        return np.dot(red,self.rprimd).tolist()
    def xcart2xred(self,cart):
        return np.dot(cart,np.linalg.inv(self.rprimd)).tolist()
    def absxcart2xred(self,abscart):#reduced coordinates of an unwrapped position, wrapped into [0,1)
        red = np.dot(abscart,np.linalg.inv(self.rprimd))
        return (red - np.floor(red)).tolist()
    def allxred2xcart(self):
        for atom in self.atoms[:self.natom]:
            atom.xcart = self.xred2xcart(atom.xred)
    def allxcart2xred(self,cart=None):#the cart argument is not used (kept for the former signature)
        for atom in self.atoms[:self.natom]:
            atom.xred = self.xcart2xred(atom.xcart)
    def allabsxcart2xred(self):
        for atom in self.atoms[:self.natom]:
            atom.xred = self.absxcart2xred(atom.absxcart)
    def makerprimd(self):                     #computes the basis metrics as rprimd = acell.rprim
        for ii in range(3):
            for jj in range(3):
                self.rprimd[ii][jj] = self.acell[ii] * self.rprim[ii][jj]
    def makevolume(self):
        cellvolume = np.abs(np.dot(np.cross(self.rprimd[0][:],self.rprimd[1][:]),self.rprimd[2][:]))
        return cellvolume
    def makegprimd(self):#gprimd = inverse(rprimd) transposed : xred[i] = sum_j gprimd[i][j] xcart[j] ; any handedness
        cell = np.array(self.rprimd,dtype=float)
        if abs(np.linalg.det(cell)) < 1.0E-06:
            print('The volume of the cell is 0 : no reduced coordinates')
            return [[0.0 for _ in range(3)] for _ in range(3)]
        return np.linalg.inv(cell).T.tolist()
    def getdensity(self):
        atomicname=''
        atomicnumber=0
        density = 0.0
        self.cellvolume = self.makevolume()
        for itypat in range(self.ntypat):
            (atomicname,atomicsymbol, atomicnumber,self.masses[itypat])=Elements2rest(self.elements[itypat])
            density = density + self.masses[itypat] * self.types[itypat]
        density = density / (6.02214076E+23 * self.cellvolume*1E-24)
        return density


def Elements2rest(identifier):#Name, symbol, atomic number and mass of an element given by any of them ; None if unknown
    entry = element_entry(identifier)
    return None if entry is None else tuple(entry[k] for k in ('atomicname','atomicsymbol','atomicnumber','atomicmass'))


def Elements2charge(identifier):#Formal charge : the most common oxidation state in oxides and ionic compounds (O -2, Si +4, Fe +2, C +4, N -3, S -2 ... ;
    #elements with several states, e.g. Fe 2/3, C, N, S in molecules, need a choice by the user)
    entry = element_entry(identifier)
    return None if entry is None else entry['formalcharge']


def element_entry(identifier):
    PeriodicTable = {}
    PeriodicTable[0]={'atomicname':'Actinium',  'atomicsymbol':'Ac',  'atomicnumber':    89, 'atomicmass':    227.028, 'formalcharge':   3}
    PeriodicTable[1]={'atomicname':'Aluminum',  'atomicsymbol':'Al',  'atomicnumber':    13, 'atomicmass':    26.9815395, 'formalcharge':   3}
    PeriodicTable[2]={'atomicname':'Americium',  'atomicsymbol':'Am',  'atomicnumber':    95, 'atomicmass':    243, 'formalcharge':   3}
    PeriodicTable[3]={'atomicname':'Antimony',  'atomicsymbol':'Sb',  'atomicnumber':    51, 'atomicmass':    121.757, 'formalcharge':   3}
    PeriodicTable[4]={'atomicname':'Argon',  'atomicsymbol':'Ar',  'atomicnumber':    18, 'atomicmass':    39.9481, 'formalcharge':   0}
    PeriodicTable[5]={'atomicname':'Arsenic',  'atomicsymbol':'As',  'atomicnumber':    33, 'atomicmass':    74.921592, 'formalcharge':   5}
    PeriodicTable[6]={'atomicname':'Astatine',  'atomicsymbol':'At',  'atomicnumber':    85, 'atomicmass':    210, 'formalcharge':  -1}
    PeriodicTable[7]={'atomicname':'Barium',  'atomicsymbol':'Ba',  'atomicnumber':    56, 'atomicmass':    137.3277, 'formalcharge':   2}
    PeriodicTable[8]={'atomicname':'Berkelium',  'atomicsymbol':'Bk',  'atomicnumber':    97, 'atomicmass':    247, 'formalcharge':   3}
    PeriodicTable[9]={'atomicname':'Beryllium',  'atomicsymbol':'Be',  'atomicnumber':    4, 'atomicmass':    9.0121823, 'formalcharge':   2}
    PeriodicTable[10]={'atomicname':'Bismuth',  'atomicsymbol':'Bi',  'atomicnumber':    83, 'atomicmass':    208.980373, 'formalcharge':   3}
    PeriodicTable[11]={'atomicname':'Bohrium',  'atomicsymbol':'Bh',  'atomicnumber':    107, 'atomicmass':    262, 'formalcharge':   7}
    PeriodicTable[12]={'atomicname':'Boron',  'atomicsymbol':'B',  'atomicnumber':    5, 'atomicmass':    10.8115, 'formalcharge':   3}
    PeriodicTable[13]={'atomicname':'Bromine',  'atomicsymbol':'Br',  'atomicnumber':    35, 'atomicmass':    79.904, 'formalcharge':  -1}
    PeriodicTable[14]={'atomicname':'Cadmium',  'atomicsymbol':'Cd',  'atomicnumber':    48, 'atomicmass':    112.4118, 'formalcharge':   2}
    PeriodicTable[15]={'atomicname':'Calcium',  'atomicsymbol':'Ca',  'atomicnumber':    20, 'atomicmass':    40.0789, 'formalcharge':   2}
    PeriodicTable[16]={'atomicname':'Californium',  'atomicsymbol':'Cf',  'atomicnumber':    98, 'atomicmass':    251, 'formalcharge':   3}
    PeriodicTable[17]={'atomicname':'Carbon',  'atomicsymbol':'C',  'atomicnumber':    6, 'atomicmass':    12.0111, 'formalcharge':   4}
    PeriodicTable[18]={'atomicname':'Cerium',  'atomicsymbol':'Ce',  'atomicnumber':    58, 'atomicmass':    140.1154, 'formalcharge':   3}
    PeriodicTable[19]={'atomicname':'Cesium',  'atomicsymbol':'Cs',  'atomicnumber':    55, 'atomicmass':    132.905435, 'formalcharge':   1}
    PeriodicTable[20]={'atomicname':'Chlorine',  'atomicsymbol':'Cl',  'atomicnumber':    17, 'atomicmass':    35.45279, 'formalcharge':  -1}
    PeriodicTable[21]={'atomicname':'Chromium',  'atomicsymbol':'Cr',  'atomicnumber':    24, 'atomicmass':    51.99616, 'formalcharge':   3}
    PeriodicTable[22]={'atomicname':'Cobalt',  'atomicsymbol':'Co',  'atomicnumber':    27, 'atomicmass':    58.933201, 'formalcharge':   2}
    PeriodicTable[23]={'atomicname':'Copper',  'atomicsymbol':'Cu',  'atomicnumber':    29, 'atomicmass':    63.5463, 'formalcharge':   2}
    PeriodicTable[24]={'atomicname':'Curium',  'atomicsymbol':'Cm',  'atomicnumber':    96, 'atomicmass':    247, 'formalcharge':   3}
    PeriodicTable[25]={'atomicname':'Dubnium',  'atomicsymbol':'Db',  'atomicnumber':    105, 'atomicmass':    262, 'formalcharge':   5}
    PeriodicTable[26]={'atomicname':'Dysprosium',  'atomicsymbol':'Dy',  'atomicnumber':    66, 'atomicmass':    162.503, 'formalcharge':   3}
    PeriodicTable[27]={'atomicname':'Einsteinium',  'atomicsymbol':'Es',  'atomicnumber':    99, 'atomicmass':    252, 'formalcharge':   3}
    PeriodicTable[28]={'atomicname':'Erbium',  'atomicsymbol':'Er',  'atomicnumber':    68, 'atomicmass':    167.263, 'formalcharge':   3}
    PeriodicTable[29]={'atomicname':'Europium',  'atomicsymbol':'Eu',  'atomicnumber':    63, 'atomicmass':    151.9659, 'formalcharge':   3}
    PeriodicTable[30]={'atomicname':'Fermium',  'atomicsymbol':'Fm',  'atomicnumber':    100, 'atomicmass':    257, 'formalcharge':   3}
    PeriodicTable[31]={'atomicname':'Fluorine',  'atomicsymbol':'F',  'atomicnumber':    9, 'atomicmass':    18.99840329, 'formalcharge':  -1}
    PeriodicTable[32]={'atomicname':'Francium',  'atomicsymbol':'Fr',  'atomicnumber':    87, 'atomicmass':    223, 'formalcharge':   1}
    PeriodicTable[33]={'atomicname':'Gadolinium',  'atomicsymbol':'Gd',  'atomicnumber':    64, 'atomicmass':    157.253, 'formalcharge':   3}
    PeriodicTable[34]={'atomicname':'Gallium',  'atomicsymbol':'Ga',  'atomicnumber':    31, 'atomicmass':    69.7231, 'formalcharge':   3}
    PeriodicTable[35]={'atomicname':'Germanium',  'atomicsymbol':'Ge',  'atomicnumber':    32, 'atomicmass':    72.612, 'formalcharge':   4}
    PeriodicTable[36]={'atomicname':'Gold',  'atomicsymbol':'Au',  'atomicnumber':    79, 'atomicmass':    196.966543, 'formalcharge':   3}
    PeriodicTable[37]={'atomicname':'Hafnium',  'atomicsymbol':'Hf',  'atomicnumber':    72, 'atomicmass':    178.492, 'formalcharge':   4}
    PeriodicTable[38]={'atomicname':'Hassium',  'atomicsymbol':'Hs',  'atomicnumber':    108, 'atomicmass':    265, 'formalcharge':   8}
    PeriodicTable[39]={'atomicname':'Helium',  'atomicsymbol':'He',  'atomicnumber':    2, 'atomicmass':    4.0026022, 'formalcharge':   0}
    PeriodicTable[40]={'atomicname':'Holmium',  'atomicsymbol':'Ho',  'atomicnumber':    67, 'atomicmass':    164.930323, 'formalcharge':   3}
    PeriodicTable[41]={'atomicname':'Hydrogen',  'atomicsymbol':'H',  'atomicnumber':    1, 'atomicmass':    1.007947, 'formalcharge':   1}
    PeriodicTable[42]={'atomicname':'Indium',  'atomicsymbol':'In',  'atomicnumber':    49, 'atomicmass':    114.821, 'formalcharge':   3}
    PeriodicTable[43]={'atomicname':'Iodine',  'atomicsymbol':'I',  'atomicnumber':    53, 'atomicmass':    126.904473, 'formalcharge':  -1}
    PeriodicTable[44]={'atomicname':'Iridium',  'atomicsymbol':'Ir',  'atomicnumber':    77, 'atomicmass':    192.223, 'formalcharge':   4}
    PeriodicTable[45]={'atomicname':'Iron',  'atomicsymbol':'Fe',  'atomicnumber':    26, 'atomicmass':    55.8473, 'formalcharge':   2}
    PeriodicTable[46]={'atomicname':'Krypton',  'atomicsymbol':'Kr',  'atomicnumber':    36, 'atomicmass':    83.801, 'formalcharge':   0}
    PeriodicTable[47]={'atomicname':'Lanthanum',  'atomicsymbol':'La',  'atomicnumber':    57, 'atomicmass':    138.90552, 'formalcharge':   3}
    PeriodicTable[48]={'atomicname':'Lawrencium',  'atomicsymbol':'Lr',  'atomicnumber':    103, 'atomicmass':    262, 'formalcharge':   3}
    PeriodicTable[49]={'atomicname':'Lead',  'atomicsymbol':'Pb',  'atomicnumber':    82, 'atomicmass':    207.21, 'formalcharge':   2}
    PeriodicTable[50]={'atomicname':'Lithium',  'atomicsymbol':'Li',  'atomicnumber':    3, 'atomicmass':    6.9412, 'formalcharge':   1}
    PeriodicTable[51]={'atomicname':'Lutetium',  'atomicsymbol':'Lu',  'atomicnumber':    71, 'atomicmass':    174.9671, 'formalcharge':   3}
    PeriodicTable[52]={'atomicname':'Magnesium',  'atomicsymbol':'Mg',  'atomicnumber':    12, 'atomicmass':    24.30506, 'formalcharge':   2}
    PeriodicTable[53]={'atomicname':'Manganese',  'atomicsymbol':'Mn',  'atomicnumber':    25, 'atomicmass':    54.938051, 'formalcharge':   2}
    PeriodicTable[54]={'atomicname':'Meitnerium',  'atomicsymbol':'Mt',  'atomicnumber':    109, 'atomicmass':    266, 'formalcharge':   3}
    PeriodicTable[55]={'atomicname':'Mendelevium',  'atomicsymbol':'Md',  'atomicnumber':    101, 'atomicmass':    258, 'formalcharge':   3}
    PeriodicTable[56]={'atomicname':'Mercury',  'atomicsymbol':'Hg',  'atomicnumber':    80, 'atomicmass':    200.593, 'formalcharge':   2}
    PeriodicTable[57]={'atomicname':'Molybdenum',  'atomicsymbol':'Mo',  'atomicnumber':    42, 'atomicmass':    95.941, 'formalcharge':   6}
    PeriodicTable[58]={'atomicname':'Neodymium',  'atomicsymbol':'Nd',  'atomicnumber':    60, 'atomicmass':    144.243, 'formalcharge':   3}
    PeriodicTable[59]={'atomicname':'Neon',  'atomicsymbol':'Ne',  'atomicnumber':    10, 'atomicmass':    20.17976, 'formalcharge':   0}
    PeriodicTable[60]={'atomicname':'Neptunium',  'atomicsymbol':'Np',  'atomicnumber':    93, 'atomicmass':    237.048, 'formalcharge':   5}
    PeriodicTable[61]={'atomicname':'Nickel',  'atomicsymbol':'Ni',  'atomicnumber':    28, 'atomicmass':    58.6934, 'formalcharge':   2}
    PeriodicTable[62]={'atomicname':'Niobium',  'atomicsymbol':'Nb',  'atomicnumber':    41, 'atomicmass':    92.906382, 'formalcharge':   5}
    PeriodicTable[63]={'atomicname':'Nitrogen',  'atomicsymbol':'N',  'atomicnumber':    7, 'atomicmass':    14.006747, 'formalcharge':  -3}
    PeriodicTable[64]={'atomicname':'Nobelium',  'atomicsymbol':'No',  'atomicnumber':    102, 'atomicmass':    259, 'formalcharge':   2}
    PeriodicTable[65]={'atomicname':'Osmium',  'atomicsymbol':'Os',  'atomicnumber':    76, 'atomicmass':    190.21, 'formalcharge':   4}
    PeriodicTable[66]={'atomicname':'Oxygen',  'atomicsymbol':'O',  'atomicnumber':    8, 'atomicmass':    15.99943, 'formalcharge':  -2}
    PeriodicTable[67]={'atomicname':'Palladium',  'atomicsymbol':'Pd',  'atomicnumber':    46, 'atomicmass':    106.421, 'formalcharge':   2}
    PeriodicTable[68]={'atomicname':'Phosphorus',  'atomicsymbol':'P',  'atomicnumber':    15, 'atomicmass':    30.9737624, 'formalcharge':   5}
    PeriodicTable[69]={'atomicname':'Platinum',  'atomicsymbol':'Pt',  'atomicnumber':    78, 'atomicmass':    195.083, 'formalcharge':   2}
    PeriodicTable[70]={'atomicname':'Plutonium',  'atomicsymbol':'Pu',  'atomicnumber':    94, 'atomicmass':    244, 'formalcharge':   4}
    PeriodicTable[71]={'atomicname':'Polonium',  'atomicsymbol':'Po',  'atomicnumber':    84, 'atomicmass':    209, 'formalcharge':   4}
    PeriodicTable[72]={'atomicname':'Potassium',  'atomicsymbol':'K',  'atomicnumber':    19, 'atomicmass':    39.09831, 'formalcharge':   1}
    PeriodicTable[73]={'atomicname':'Praseodymium',  'atomicsymbol':'Pr',  'atomicnumber':    59, 'atomicmass':    140.907653, 'formalcharge':   3}
    PeriodicTable[74]={'atomicname':'Promethium',  'atomicsymbol':'Pm',  'atomicnumber':    61, 'atomicmass':    145, 'formalcharge':   3}
    PeriodicTable[75]={'atomicname':'Protactinium',  'atomicsymbol':'Pa',  'atomicnumber':    91, 'atomicmass':    231.0359, 'formalcharge':   5}
    PeriodicTable[76]={'atomicname':'Radium',  'atomicsymbol':'Ra',  'atomicnumber':    88, 'atomicmass':    226.025, 'formalcharge':   2}
    PeriodicTable[77]={'atomicname':'Radon',  'atomicsymbol':'Rn',  'atomicnumber':    86, 'atomicmass':    222, 'formalcharge':   0}
    PeriodicTable[78]={'atomicname':'Rhenium',  'atomicsymbol':'Re',  'atomicnumber':    75, 'atomicmass':    186.2071, 'formalcharge':   7}
    PeriodicTable[79]={'atomicname':'Rhodium',  'atomicsymbol':'Rh',  'atomicnumber':    45, 'atomicmass':    102.905503, 'formalcharge':   3}
    PeriodicTable[80]={'atomicname':'Rubidium',  'atomicsymbol':'Rb',  'atomicnumber':    37, 'atomicmass':    85.46783, 'formalcharge':   1}
    PeriodicTable[81]={'atomicname':'Ruthenium',  'atomicsymbol':'Ru',  'atomicnumber':    44, 'atomicmass':    101.072, 'formalcharge':   4}
    PeriodicTable[82]={'atomicname':'Rutherfordium',  'atomicsymbol':'Rf',  'atomicnumber':    104, 'atomicmass':    261, 'formalcharge':   4}
    PeriodicTable[83]={'atomicname':'Samarium',  'atomicsymbol':'Sm',  'atomicnumber':    62, 'atomicmass':    150.363, 'formalcharge':   3}
    PeriodicTable[84]={'atomicname':'Scandium',  'atomicsymbol':'Sc',  'atomicnumber':    21, 'atomicmass':    44.9559109, 'formalcharge':   3}
    PeriodicTable[85]={'atomicname':'Seaborgium',  'atomicsymbol':'Sg',  'atomicnumber':    106, 'atomicmass':    263, 'formalcharge':   6}
    PeriodicTable[86]={'atomicname':'Selenium',  'atomicsymbol':'Se',  'atomicnumber':    34, 'atomicmass':    78.963, 'formalcharge':  -2}
    PeriodicTable[87]={'atomicname':'Silicon',  'atomicsymbol':'Si',  'atomicnumber':    14, 'atomicmass':    28.08553, 'formalcharge':   4}
    PeriodicTable[88]={'atomicname':'Silver',  'atomicsymbol':'Ag',  'atomicnumber':    47, 'atomicmass':    107.86822, 'formalcharge':   1}
    PeriodicTable[89]={'atomicname':'Sodium',  'atomicsymbol':'Na',  'atomicnumber':    11, 'atomicmass':    22.9897686, 'formalcharge':   1}
    PeriodicTable[90]={'atomicname':'Strontium',  'atomicsymbol':'Sr',  'atomicnumber':    38, 'atomicmass':    87.621, 'formalcharge':   2}
    PeriodicTable[91]={'atomicname':'Sulfur',  'atomicsymbol':'S',  'atomicnumber':    16, 'atomicmass':    32.0666, 'formalcharge':  -2}
    PeriodicTable[92]={'atomicname':'Tantalum',  'atomicsymbol':'Ta',  'atomicnumber':    73, 'atomicmass':    180.94791, 'formalcharge':   5}
    PeriodicTable[93]={'atomicname':'Technetium',  'atomicsymbol':'Tc',  'atomicnumber':    43, 'atomicmass':    98, 'formalcharge':   7}
    PeriodicTable[94]={'atomicname':'Tellurium',  'atomicsymbol':'Te',  'atomicnumber':    52, 'atomicmass':    127.603, 'formalcharge':  -2}
    PeriodicTable[95]={'atomicname':'Terbium',  'atomicsymbol':'Tb',  'atomicnumber':    65, 'atomicmass':    158.925343, 'formalcharge':   3}
    PeriodicTable[96]={'atomicname':'Thallium',  'atomicsymbol':'Tl',  'atomicnumber':    81, 'atomicmass':    204.38332, 'formalcharge':   1}
    PeriodicTable[97]={'atomicname':'Thorium',  'atomicsymbol':'Th',  'atomicnumber':    90, 'atomicmass':    232.03811, 'formalcharge':   4}
    PeriodicTable[98]={'atomicname':'Thulium',  'atomicsymbol':'Tm',  'atomicnumber':    69, 'atomicmass':    168.934213, 'formalcharge':   3}
    PeriodicTable[99]={'atomicname':'Tin',  'atomicsymbol':'Sn',  'atomicnumber':    50, 'atomicmass':    118.7107, 'formalcharge':   4}
    PeriodicTable[100]={'atomicname':'Titanium',  'atomicsymbol':'Ti',  'atomicnumber':    22, 'atomicmass':    47.883, 'formalcharge':   4}
    PeriodicTable[101]={'atomicname':'Tungsten',  'atomicsymbol':'W',  'atomicnumber':    74, 'atomicmass':    183.853, 'formalcharge':   6}
    PeriodicTable[102]={'atomicname':'Uranium',  'atomicsymbol':'U',  'atomicnumber':    92, 'atomicmass':    238.02891, 'formalcharge':   4}
    PeriodicTable[103]={'atomicname':'Vanadium',  'atomicsymbol':'V',  'atomicnumber':    23, 'atomicmass':    50.94151, 'formalcharge':   5}
    PeriodicTable[104]={'atomicname':'Xenon',  'atomicsymbol':'Xe',  'atomicnumber':    54, 'atomicmass':    131.292, 'formalcharge':   0}
    PeriodicTable[105]={'atomicname':'Ytterbium',  'atomicsymbol':'Yb',  'atomicnumber':    70, 'atomicmass':    173.043, 'formalcharge':   3}
    PeriodicTable[106]={'atomicname':'Yttrium',  'atomicsymbol':'Y',  'atomicnumber':    39, 'atomicmass':    88.905852, 'formalcharge':   3}
    PeriodicTable[107]={'atomicname':'Zinc',  'atomicsymbol':'Zn',  'atomicnumber':    30, 'atomicmass':    65.392, 'formalcharge':   2}
    PeriodicTable[108]={'atomicname':'Zirconium',  'atomicsymbol':'Zr',  'atomicnumber':    40, 'atomicmass':    91.2242, 'formalcharge':   4}
    for kk in PeriodicTable.keys():
        entry = PeriodicTable[kk]
        if identifier in (entry['atomicname'],entry['atomicsymbol'],entry['atomicnumber'],entry['atomicmass']):
            return entry



