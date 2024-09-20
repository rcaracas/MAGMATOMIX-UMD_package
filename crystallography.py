#!/usr/bin/env python3
###
##AUTHORS: RAZVAN CARACAS
###

import math
import numpy as np

class Atom(object):
    def __init__(self, name = 'XXX', symbol='X', znucl = 1, mass=0.0, magnet=0.0, charge=0.0, xred = None, xcart = None, absxcart = None, vels = None, forces = None, diffcoords = None):
        self.name = name
        self.symbol = symbol
        self.znucl = znucl
        self.mass = mass
        self.magnet = magnet
        self.charge = charge
        if xred == None:
            self.xred = [0.0, 0.0, 0.0]
        else:
            self.xred = xred
        if xcart == None:
            self.xcart = [0.0, 0.0, 0.0]
        if absxcart == None:
            self.absxcart = [0.0, 0.0, 0.0]
        else:
            self.xcart = xcart
        if vels == None:
            self.vels = [0.0, 0.0, 0.0]
        else:
            self.vels = vels
        if forces == None:
            self.forces = [0.0, 0.0, 0.0]
        else:
            self.forces = forces
        if diffcoords == None:
            self.diffcoords=[0.0, 0.0, 0.0]
        else:
            self.diffcoords = diffcoords    


class Lattice(object):#Remplacer les vecteurs de vecteurs par des matrices ?
    def __init__(self, name='crystal', acell=[0.0 for x in range(3)], angles=[0.0 for x in range(3)],
                 rprim=[[0.0 for x in range(3)] for y in range(3)], rprimd=[[0.0 for x in range(3)] for y in range(3)],
                 gprimd=[[0.0 for x in range(3)] for y in range(3)], stress=[0.0 for x in range(6)],
                 typat=[],elements=[],masses=[],zelec=[],ntypat=1, natom=1, density=0.0, atoms=[],cellvolume=1.0,enthalpy=0.0,internalenergy=0.0,electronicentropy = 0.0,kineticenergy=0.0,energywithdrift=0.0,gibbsfreeenergy=0.0,cv=0.0,pressure=0.0,temperature=0.0,magnetization=0.0,elecgap=0.0,noelectrons=0.0,lambda_ThermoInt=1.0):
        self.name = name
        self.acell = acell
        self.rprim = rprim
        self.rprimd = rprimd
        self.gprimd = gprimd
        self.natom = natom                                      #number of atoms
        self.ntypat = ntypat                                    #number of atom types
        self.density = density
        self.typat = [0 for x in range(self.natom)]             #gives the type of each atom
        self.types = [0 for x in range(self.ntypat)]            #gives the no. of atoms of each type
        self.elements = ['X' for x in range(self.ntypat)]       #gives the Chemical Symbol for each atom
        self.masses = ['0.0' for x in range(self.ntypat)]       #gives the atomic mass of each atom
        self.zelec = ['0.0' for x in range(self.ntypat)]        #gives the atomic mass of each atom
        self.atoms = [Atom() for x in range(self.natom)]        #gives the complete Atom() class for each atom
        self.angles = angles
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
        self.stress=stress
        self.temperature=temperature
        self.magnetization=magnetization
        self.elecgap=elecgap
        self.noelectrons=noelectrons
        self.lambda_ThermoInt=lambda_ThermoInt
    def allxred2xcart(self):
        self.makerprimd()
        for iatom in range(self.natom):
            for ii in range(3):
                self.atoms[iatom].xcart[ii]=0.0
                for jj in range(3):
                    self.atoms[iatom].xcart[ii] = self.atoms[iatom].xcart[ii] + self.rprimd[ii][jj]*self.atoms[iatom].xred[ii]
    def xred2xcart(self,red):
        self.makerprimd()
        cart = [0.0,0.0,0.0]
        for ii in range(3):
            for jj in range(3):
                cart[ii] = cart[ii] + self.rprimd[ii][jj]*red[ii]
        return cart
    def allxcart2xred(self,cart):
        self.makerprimd()
        self.makegprimd()
        for iatom in range(self.natom):
            for ii in range(3):
                self.atoms[iatom].xred[ii]=0.0
                for jj in range(3):
                    self.atoms[iatom].xred[ii] = self.atoms[iatom].xred[ii]  + self.gprimd[ii][jj]*self.atoms[iatom].xcart[ii]
    def xcart2xred(self,cart):
        self.makerprimd()
        self.makegprimd()
        red=[0.0,0.0,0.0]
        for ii in range(3):
            for jj in range(3):
                red[ii] = red[ii] + self.gprimd[ii][jj]*cart[ii]
        return red
    def allabsxcart2xred(self):
        self.makerprimd()
        self.makegprimd()
        for iatom in range(self.natom):
            for ii in range(3):
                self.atoms[iatom].xred[ii]=0.0
                for jj in range(3):
                    self.atoms[iatom].xred[ii] = self.atoms[iatom].xred[ii]  + self.gprimd[ii][jj]*self.atoms[iatom].absxcart[ii]
                while self.atoms[iatom].xred[ii] < 0.0:
                    self.atoms[iatom].xred[ii] += 1.0
                while self.atoms[iatom].xred[ii] >= 1.0:
                    self.atoms[iatom].xred[ii] -= 1.0
    def absxcart2xred(self,abscart):
        self.makerprimd()
        self.makegprimd()
        red=[0.0,0.0,0.0]
        for ii in range(3):
            for jj in range(3):
                red[ii] = red[ii] + self.gprimd[ii][jj]*abscart[ii]
            while red[ii] < 0:
                red[ii] += 1.0
            while red[ii] >= 1.0:
                red[ii] -= 1.0
        return red
    def makerprimd(self):                     #computes the basis metrics as rprimd = acell.rprim
        for ii in range(3):
            for jj in range(3):
                self.rprimd[ii][jj] = self.acell[ii] * self.rprim[ii][jj]
    def makevolume(self):
        cellvolume = np.abs(np.dot(np.cross(self.rprimd[0][:],self.rprimd[1][:]),self.rprimd[2][:]))
        return cellvolume
    def makegprimd(self):
        gprimd = [ [0.0 for _ in range(3)] for _ in range(3)]
        t1 = self.rprimd[1][1] * self.rprimd[2][2] - self.rprimd[2][1] * self.rprimd[1][2]
        t2 = self.rprimd[2][1] * self.rprimd[0][2] - self.rprimd[0][1] * self.rprimd[2][2]
        t3 = self.rprimd[0][1] * self.rprimd[1][2] - self.rprimd[1][1] * self.rprimd[0][2]
        det  = self.rprimd[0][0] * t1 + self.rprimd[1][0] * t2 + self.rprimd[2][0] * t3
        dd = 1.0
        if (det>1.0E-06):
            dd = 1.0/det
        else:
            print ('The norm of the crystallographic space is 0.0 or negative. This is not allowed ')
        gprimd[0][0] = t1 * dd
        gprimd[1][0] = t2 * dd
        gprimd[2][0] = t3 * dd
        gprimd[0][1] = (self.rprimd[2][0]*self.rprimd[1][2]-self.rprimd[1][0]*self.rprimd[2][2]) * dd
        gprimd[1][1] = (self.rprimd[0][0]*self.rprimd[2][2]-self.rprimd[2][0]*self.rprimd[0][2]) * dd
        gprimd[2][1] = (self.rprimd[1][0]*self.rprimd[0][2]-self.rprimd[0][0]*self.rprimd[1][2]) * dd
        gprimd[0][2] = (self.rprimd[1][0]*self.rprimd[2][1]-self.rprimd[2][0]*self.rprimd[1][1]) * dd
        gprimd[1][2] = (self.rprimd[2][0]*self.rprimd[0][1]-self.rprimd[0][0]*self.rprimd[2][1]) * dd
        gprimd[2][2] = (self.rprimd[0][0]*self.rprimd[1][1]-self.rprimd[1][0]*self.rprimd[0][1]) * dd
        return gprimd
    def getdensity(self):
        atomicname=''
        atomicnumber=0
        density = 0.0
        self.cellvolume = self.makevolume()
        for itypat in range(self.ntypat):
            (atomicname,atomicsymbol, atomicnumber,self.masses[itypat])=Elements2rest(self.elements[itypat])
            density = density + self.masses[itypat] * self.types[itypat]
        density = density / (6.022E+23 * self.cellvolume*1E-24)
        return density
    def compangles(self):                   #compute the angles between the lattice vectors
        self.cartlatt()
        index = -1
        for ii in range(3):                 #loop over the first vector, ai   
            for jj in range(ii+1,3):        #loop over the second vector, aj
                cosa = 0.0
                for kk in range(3):         #loop over the terms of the dot product ai.aj
                    cosa = cosa + self.rprimd[ii][kk]*self.rprimd[jj][kk]
                cosa = cosa / (self.acell[ii] * self.acell[jj])
                index = index +1
                self.angles[index] = math.acos(cosa) * 180.0 / math.pi
    def printstruct(self):
        print ('a,b,c =',self.acell)
        print ('rprim =',self.rprim)
        self.cartlatt()
        print ('rprimd =',self.rprimd)
        self.compangles()
        print ('angles =',self.angles)
    def niceprint(self):
        self.cartlatt()
        self.compangles()
        print ('\t',self.acell[0],'\t',self.acell[1],'\t',self.acell[2],'\t',self.angles[0],'\t',self.angles[1],'\t',self.angles[2])
    def printcontcar(self,filetype):
        self.makerprimd()
        filename=self.name+'_.'+filetype
        f=open(filename,'w')
        string = self.name + '\n'
        f.close()
        f=open(filename,'a')
        f.write(string)
        f.write('  1.0\n')
        for ii in range(3):
            string = '      ' + str(self.rprimd[ii][0]) + '  ' + str(self.rprimd[ii][1]) + '  ' + str(self.rprimd[ii][2]) +  '\n'
            f.write(string)
        string = ' '
        for ii in range(self.ntypat):
            string = string + '   ' + str(self.elements[ii])
        string = string+'\n'
        f.write(string)
        string = ' '
        for ii in range(self.ntypat):
            string = string + '   ' + str(self.typat[ii])
        string = string+'\n'
        f.write(string)
        f.write('Direct\n')
        for iatom in range(self.natom):
            string = '  ' + str(self.atoms[iatom].xred[0]) + '  ' + str(self.atoms[iatom].xred[1]) + '  ' + str(self.atoms[iatom].xred[2])+'\n'
            f.write(string)


def TripleLattice(crystal):
    #
    #builds a 3x3x3 supercell centered on the initial cell
    #this allows to take into account all the NN interactions
    #
    bigcrystal=Lattice()
    bigcrystal.name='2x2x2_'+crystal.name
    bigcrystal.acell[0] = crystal.acell[0]*3.0
    bigcrystal.acell[1] = crystal.acell[1]*3.0
    bigcrystal.acell[2] = crystal.acell[2]*3.0
    bigcrystal.rprim = crystal.rprim
    bigcrystal.natom = crystal.natom * 27
    bigcrystal.ntypat = crystal.ntypat * 27
    bigcrystal.typat = [0 for x in range(bigcrystal.ntypat)]
    bigcrystal.elements = ['E' for x in range(bigcrystal.ntypat)]    
    bigcrystal.atoms = [Atom() for iatom in range(bigcrystal.natom)]
    jatom = 0
    box = 0
    for ix in range(-1,2):
        for iy in range(-1,2):
            for iz in range(-1,2):
                bigcrystal.elements[box:box+crystal.ntypat]=crystal.elements
                bigcrystal.typat[box:box+crystal.ntypat]=crystal.typat
                box=box+crystal.ntypat
                for iatom in range(crystal.natom):
                    nx = (crystal.atoms[iatom].xred[0] + float(ix))/3.0
                    ny = (crystal.atoms[iatom].xred[1] + float(iy))/3.0
                    nz = (crystal.atoms[iatom].xred[2] + float(iz))/3.0
                    bigcrystal.atoms[jatom].xred = [nx,ny,nz]
                    bigcrystal.atoms[jatom].xcart = bigcrystal.cartatom(bigcrystal.atoms[jatom].xred)
                    bigcrystal.atoms[jatom].symbol = crystal.atoms[iatom].symbol
                    jatom = jatom + 1
    return bigcrystal


def distance(atom1,atom2):
    dd = 0.0
    for ii in range(3):
        dd = dd + (atom1.xcart[ii]-atom2.xcart[ii])*(atom1.xcart[ii]-atom2.xcart[ii])
    return math.sqrt(dd)


def CheckAtoms(Crystal,cat,an):
    flagcat=0
    flagan =0
    for ii in range(Crystal.ntypat):
        if (Crystal.elements[ii] == cat):
            flagcat = 1
        if (Crystal.elements[ii] == an):
            flagan = 1
    if (flagcat + flagan == 2):
        return 1
    else:
        return 0


def Elements2rest(identifier):
    PeriodicTable = {}
    PeriodicTable[0]={'atomicname':'Actinium',  'atomicsymbol':'Ac',  'atomicnumber':    89, 'atomicmass':    227.028}
    PeriodicTable[1]={'atomicname':'Aluminum',  'atomicsymbol':'Al',  'atomicnumber':    13, 'atomicmass':    26.9815395}
    PeriodicTable[2]={'atomicname':'Americium',  'atomicsymbol':'Am',  'atomicnumber':    95, 'atomicmass':    243}
    PeriodicTable[3]={'atomicname':'Antimony',  'atomicsymbol':'Sb',  'atomicnumber':    51, 'atomicmass':    121.757}
    PeriodicTable[4]={'atomicname':'Argon',  'atomicsymbol':'Ar',  'atomicnumber':    18, 'atomicmass':    39.9481}
    PeriodicTable[5]={'atomicname':'Arsenic',  'atomicsymbol':'As',  'atomicnumber':    33, 'atomicmass':    74.921592}
    PeriodicTable[6]={'atomicname':'Astatine',  'atomicsymbol':'At',  'atomicnumber':    85, 'atomicmass':    210}
    PeriodicTable[7]={'atomicname':'Barium',  'atomicsymbol':'Ba',  'atomicnumber':    56, 'atomicmass':    137.3277}
    PeriodicTable[8]={'atomicname':'Berkelium',  'atomicsymbol':'Bk',  'atomicnumber':    97, 'atomicmass':    247}
    PeriodicTable[9]={'atomicname':'Beryllium',  'atomicsymbol':'Be',  'atomicnumber':    4, 'atomicmass':    9.0121823}
    PeriodicTable[10]={'atomicname':'Bismuth',  'atomicsymbol':'Bi',  'atomicnumber':    83, 'atomicmass':    208.980373}
    PeriodicTable[11]={'atomicname':'Bohrium',  'atomicsymbol':'Bh',  'atomicnumber':    107, 'atomicmass':    262}
    PeriodicTable[12]={'atomicname':'Boron',  'atomicsymbol':'B',  'atomicnumber':    5, 'atomicmass':    10.8115}
    PeriodicTable[13]={'atomicname':'Bromine',  'atomicsymbol':'Br',  'atomicnumber':    35, 'atomicmass':    79.904}
    PeriodicTable[14]={'atomicname':'Cadmium',  'atomicsymbol':'Cd',  'atomicnumber':    48, 'atomicmass':    112.4118}
    PeriodicTable[15]={'atomicname':'Calcium',  'atomicsymbol':'Ca',  'atomicnumber':    20, 'atomicmass':    40.0789}
    PeriodicTable[16]={'atomicname':'Californium',  'atomicsymbol':'Cf',  'atomicnumber':    98, 'atomicmass':    251}
    PeriodicTable[17]={'atomicname':'Carbon',  'atomicsymbol':'C',  'atomicnumber':    6, 'atomicmass':    12.0111}
    PeriodicTable[18]={'atomicname':'Cerium',  'atomicsymbol':'Ce',  'atomicnumber':    58, 'atomicmass':    140.1154}
    PeriodicTable[19]={'atomicname':'Cesium',  'atomicsymbol':'Cs',  'atomicnumber':    55, 'atomicmass':    132.905435}
    PeriodicTable[20]={'atomicname':'Chlorine',  'atomicsymbol':'Cl',  'atomicnumber':    17, 'atomicmass':    35.45279}
    PeriodicTable[21]={'atomicname':'Chromium',  'atomicsymbol':'Cr',  'atomicnumber':    24, 'atomicmass':    51.99616}
    PeriodicTable[22]={'atomicname':'Cobalt',  'atomicsymbol':'Co',  'atomicnumber':    27, 'atomicmass':    58.933201}
    PeriodicTable[23]={'atomicname':'Copper',  'atomicsymbol':'Cu',  'atomicnumber':    29, 'atomicmass':    63.5463}
    PeriodicTable[24]={'atomicname':'Curium',  'atomicsymbol':'Cm',  'atomicnumber':    96, 'atomicmass':    247}
    PeriodicTable[25]={'atomicname':'Dubnium',  'atomicsymbol':'Db',  'atomicnumber':    105, 'atomicmass':    262}
    PeriodicTable[26]={'atomicname':'Dysprosium',  'atomicsymbol':'Dy',  'atomicnumber':    66, 'atomicmass':    162.503}
    PeriodicTable[27]={'atomicname':'Einsteinium',  'atomicsymbol':'Es',  'atomicnumber':    99, 'atomicmass':    252}
    PeriodicTable[28]={'atomicname':'Erbium',  'atomicsymbol':'Er',  'atomicnumber':    68, 'atomicmass':    167.263}
    PeriodicTable[29]={'atomicname':'Europium',  'atomicsymbol':'Eu',  'atomicnumber':    63, 'atomicmass':    151.9659}
    PeriodicTable[30]={'atomicname':'Fermium',  'atomicsymbol':'Fm',  'atomicnumber':    100, 'atomicmass':    257}
    PeriodicTable[31]={'atomicname':'Fluorine',  'atomicsymbol':'F',  'atomicnumber':    9, 'atomicmass':    18.99840329}
    PeriodicTable[32]={'atomicname':'Francium',  'atomicsymbol':'Fr',  'atomicnumber':    87, 'atomicmass':    223}
    PeriodicTable[33]={'atomicname':'Gadolinium',  'atomicsymbol':'Gd',  'atomicnumber':    64, 'atomicmass':    157.253}
    PeriodicTable[34]={'atomicname':'Gallium',  'atomicsymbol':'Ga',  'atomicnumber':    31, 'atomicmass':    69.7231}
    PeriodicTable[35]={'atomicname':'Germanium',  'atomicsymbol':'Ge',  'atomicnumber':    32, 'atomicmass':    72.612}
    PeriodicTable[36]={'atomicname':'Gold',  'atomicsymbol':'Au',  'atomicnumber':    79, 'atomicmass':    196.966543}
    PeriodicTable[37]={'atomicname':'Hafnium',  'atomicsymbol':'Hf',  'atomicnumber':    72, 'atomicmass':    178.492}
    PeriodicTable[38]={'atomicname':'Hassium',  'atomicsymbol':'Hs',  'atomicnumber':    108, 'atomicmass':    265}
    PeriodicTable[39]={'atomicname':'Helium',  'atomicsymbol':'He',  'atomicnumber':    2, 'atomicmass':    4.0026022}
    PeriodicTable[40]={'atomicname':'Holmium',  'atomicsymbol':'Ho',  'atomicnumber':    67, 'atomicmass':    164.930323}
    PeriodicTable[41]={'atomicname':'Hydrogen',  'atomicsymbol':'H',  'atomicnumber':    1, 'atomicmass':    1.007947}
    PeriodicTable[42]={'atomicname':'Indium',  'atomicsymbol':'In',  'atomicnumber':    49, 'atomicmass':    114.821}
    PeriodicTable[43]={'atomicname':'Iodine',  'atomicsymbol':'I',  'atomicnumber':    53, 'atomicmass':    126.904473}
    PeriodicTable[44]={'atomicname':'Iridium',  'atomicsymbol':'Ir',  'atomicnumber':    77, 'atomicmass':    192.223}
    PeriodicTable[45]={'atomicname':'Iron',  'atomicsymbol':'Fe',  'atomicnumber':    26, 'atomicmass':    55.8473}
    PeriodicTable[46]={'atomicname':'Krypton',  'atomicsymbol':'Kr',  'atomicnumber':    36, 'atomicmass':    83.801}
    PeriodicTable[47]={'atomicname':'Lanthanum',  'atomicsymbol':'La',  'atomicnumber':    57, 'atomicmass':    138.90552}
    PeriodicTable[48]={'atomicname':'Lawrencium',  'atomicsymbol':'Lr',  'atomicnumber':    103, 'atomicmass':    262}
    PeriodicTable[49]={'atomicname':'Lead',  'atomicsymbol':'Pb',  'atomicnumber':    82, 'atomicmass':    207.21}
    PeriodicTable[50]={'atomicname':'Lithium',  'atomicsymbol':'Li',  'atomicnumber':    3, 'atomicmass':    6.9412}
    PeriodicTable[51]={'atomicname':'Lutetium',  'atomicsymbol':'Lu',  'atomicnumber':    71, 'atomicmass':    174.9671}
    PeriodicTable[52]={'atomicname':'Magnesium',  'atomicsymbol':'Mg',  'atomicnumber':    12, 'atomicmass':    24.30506}
    PeriodicTable[53]={'atomicname':'Manganese',  'atomicsymbol':'Mn',  'atomicnumber':    25, 'atomicmass':    54.938051}
    PeriodicTable[54]={'atomicname':'Meitnerium',  'atomicsymbol':'Mt',  'atomicnumber':    109, 'atomicmass':    266}
    PeriodicTable[55]={'atomicname':'Mendelevium',  'atomicsymbol':'Md',  'atomicnumber':    101, 'atomicmass':    258}
    PeriodicTable[56]={'atomicname':'Mercury',  'atomicsymbol':'Hg',  'atomicnumber':    80, 'atomicmass':    200.593}
    PeriodicTable[57]={'atomicname':'Molybdenum',  'atomicsymbol':'Mo',  'atomicnumber':    42, 'atomicmass':    95.941}
    PeriodicTable[58]={'atomicname':'Neodymium',  'atomicsymbol':'Nd',  'atomicnumber':    60, 'atomicmass':    144.243}
    PeriodicTable[59]={'atomicname':'Neon',  'atomicsymbol':'Ne',  'atomicnumber':    10, 'atomicmass':    20.17976}
    PeriodicTable[60]={'atomicname':'Neptunium',  'atomicsymbol':'Np',  'atomicnumber':    93, 'atomicmass':    237.048}
    PeriodicTable[61]={'atomicname':'Nickel',  'atomicsymbol':'Ni',  'atomicnumber':    28, 'atomicmass':    58.6934}
    PeriodicTable[62]={'atomicname':'Niobium',  'atomicsymbol':'Nb',  'atomicnumber':    41, 'atomicmass':    92.906382}
    PeriodicTable[63]={'atomicname':'Nitrogen',  'atomicsymbol':'N',  'atomicnumber':    7, 'atomicmass':    14.006747}
    PeriodicTable[64]={'atomicname':'Nobelium',  'atomicsymbol':'No',  'atomicnumber':    102, 'atomicmass':    259}
    PeriodicTable[65]={'atomicname':'Osmium',  'atomicsymbol':'Os',  'atomicnumber':    76, 'atomicmass':    190.21}
    PeriodicTable[66]={'atomicname':'Oxygen',  'atomicsymbol':'O',  'atomicnumber':    8, 'atomicmass':    15.99943}
    PeriodicTable[67]={'atomicname':'Palladium',  'atomicsymbol':'Pd',  'atomicnumber':    46, 'atomicmass':    106.421}
    PeriodicTable[68]={'atomicname':'Phosphorus',  'atomicsymbol':'P',  'atomicnumber':    15, 'atomicmass':    30.9737624}
    PeriodicTable[69]={'atomicname':'Platinum',  'atomicsymbol':'Pt',  'atomicnumber':    78, 'atomicmass':    195.083}
    PeriodicTable[70]={'atomicname':'Plutonium',  'atomicsymbol':'Pu',  'atomicnumber':    94, 'atomicmass':    244}
    PeriodicTable[71]={'atomicname':'Polonium',  'atomicsymbol':'Po',  'atomicnumber':    84, 'atomicmass':    209}
    PeriodicTable[72]={'atomicname':'Potassium',  'atomicsymbol':'K',  'atomicnumber':    19, 'atomicmass':    39.09831}
    PeriodicTable[73]={'atomicname':'Praseodymium',  'atomicsymbol':'Pr',  'atomicnumber':    59, 'atomicmass':    140.907653}
    PeriodicTable[74]={'atomicname':'Promethium',  'atomicsymbol':'Pm',  'atomicnumber':    61, 'atomicmass':    145}
    PeriodicTable[75]={'atomicname':'Protactinium',  'atomicsymbol':'Pa',  'atomicnumber':    91, 'atomicmass':    231.0359}
    PeriodicTable[76]={'atomicname':'Radium',  'atomicsymbol':'Ra',  'atomicnumber':    88, 'atomicmass':    226.025}
    PeriodicTable[77]={'atomicname':'Radon',  'atomicsymbol':'Rn',  'atomicnumber':    86, 'atomicmass':    222}
    PeriodicTable[78]={'atomicname':'Rhenium',  'atomicsymbol':'Re',  'atomicnumber':    75, 'atomicmass':    186.2071}
    PeriodicTable[79]={'atomicname':'Rhodium',  'atomicsymbol':'Rh',  'atomicnumber':    45, 'atomicmass':    102.905503}
    PeriodicTable[80]={'atomicname':'Rubidium',  'atomicsymbol':'Rb',  'atomicnumber':    37, 'atomicmass':    85.46783}
    PeriodicTable[81]={'atomicname':'Ruthenium',  'atomicsymbol':'Ru',  'atomicnumber':    44, 'atomicmass':    101.072}
    PeriodicTable[82]={'atomicname':'Rutherfordium',  'atomicsymbol':'Rf',  'atomicnumber':    104, 'atomicmass':    261}
    PeriodicTable[83]={'atomicname':'Samarium',  'atomicsymbol':'Sm',  'atomicnumber':    62, 'atomicmass':    150.363}
    PeriodicTable[84]={'atomicname':'Scandium',  'atomicsymbol':'Sc',  'atomicnumber':    21, 'atomicmass':    44.9559109}
    PeriodicTable[85]={'atomicname':'Seaborgium',  'atomicsymbol':'Sg',  'atomicnumber':    106, 'atomicmass':    263}
    PeriodicTable[86]={'atomicname':'Selenium',  'atomicsymbol':'Se',  'atomicnumber':    34, 'atomicmass':    78.963}
    PeriodicTable[87]={'atomicname':'Silicon',  'atomicsymbol':'Si',  'atomicnumber':    14, 'atomicmass':    28.08553}
    PeriodicTable[88]={'atomicname':'Silver',  'atomicsymbol':'Ag',  'atomicnumber':    47, 'atomicmass':    107.86822}
    PeriodicTable[89]={'atomicname':'Sodium',  'atomicsymbol':'Na',  'atomicnumber':    11, 'atomicmass':    22.9897686}
    PeriodicTable[90]={'atomicname':'Strontium',  'atomicsymbol':'Sr',  'atomicnumber':    38, 'atomicmass':    87.621}
    PeriodicTable[91]={'atomicname':'Sulfur',  'atomicsymbol':'S',  'atomicnumber':    16, 'atomicmass':    32.0666}
    PeriodicTable[92]={'atomicname':'Tantalum',  'atomicsymbol':'Ta',  'atomicnumber':    73, 'atomicmass':    180.94791}
    PeriodicTable[93]={'atomicname':'Technetium',  'atomicsymbol':'Tc',  'atomicnumber':    43, 'atomicmass':    98}
    PeriodicTable[94]={'atomicname':'Tellurium',  'atomicsymbol':'Te',  'atomicnumber':    52, 'atomicmass':    127.603}
    PeriodicTable[95]={'atomicname':'Terbium',  'atomicsymbol':'Tb',  'atomicnumber':    65, 'atomicmass':    158.925343}
    PeriodicTable[96]={'atomicname':'Thallium',  'atomicsymbol':'Tl',  'atomicnumber':    81, 'atomicmass':    204.38332}
    PeriodicTable[97]={'atomicname':'Thorium',  'atomicsymbol':'Th',  'atomicnumber':    90, 'atomicmass':    232.03811}
    PeriodicTable[98]={'atomicname':'Thulium',  'atomicsymbol':'Tm',  'atomicnumber':    69, 'atomicmass':    168.934213}
    PeriodicTable[99]={'atomicname':'Tin',  'atomicsymbol':'Sn',  'atomicnumber':    50, 'atomicmass':    118.7107}
    PeriodicTable[100]={'atomicname':'Titanium',  'atomicsymbol':'Ti',  'atomicnumber':    22, 'atomicmass':    47.883}
    PeriodicTable[101]={'atomicname':'Tungsten',  'atomicsymbol':'W',  'atomicnumber':    74, 'atomicmass':    183.853}
    PeriodicTable[102]={'atomicname':'Uranium',  'atomicsymbol':'U',  'atomicnumber':    92, 'atomicmass':    238.02891}
    PeriodicTable[103]={'atomicname':'Vanadium',  'atomicsymbol':'V',  'atomicnumber':    23, 'atomicmass':    50.94151}
    PeriodicTable[104]={'atomicname':'Xenon',  'atomicsymbol':'Xe',  'atomicnumber':    54, 'atomicmass':    131.292}
    PeriodicTable[105]={'atomicname':'Ytterbium',  'atomicsymbol':'Yb',  'atomicnumber':    70, 'atomicmass':    173.043}
    PeriodicTable[106]={'atomicname':'Yttrium',  'atomicsymbol':'Y',  'atomicnumber':    39, 'atomicmass':    88.905852}
    PeriodicTable[107]={'atomicname':'Zinc',  'atomicsymbol':'Zn',  'atomicnumber':    30, 'atomicmass':    65.392}
    PeriodicTable[108]={'atomicname':'Zirconium',  'atomicsymbol':'Zr',  'atomicnumber':    40, 'atomicmass':    91.2242}
    for kk in PeriodicTable.keys():
        if PeriodicTable[kk]['atomicname']==identifier:
            #print(PeriodicTable[kk]['atomicname'],PeriodicTable[kk]['atomicnumber'],PeriodicTable[kk]['atomicmass'])
            return(PeriodicTable[kk]['atomicname'],PeriodicTable[kk]['atomicsymbol'],PeriodicTable[kk]['atomicnumber'],PeriodicTable[kk]['atomicmass'])
            break
        if PeriodicTable[kk]['atomicsymbol']==identifier:
            #print(PeriodicTable[kk]['atomicname'],PeriodicTable[kk]['atomicnumber'],PeriodicTable[kk]['atomicmass'])
            return(PeriodicTable[kk]['atomicname'],PeriodicTable[kk]['atomicsymbol'],PeriodicTable[kk]['atomicnumber'],PeriodicTable[kk]['atomicmass'])
            break
        if PeriodicTable[kk]['atomicnumber']==identifier:
            #print(PeriodicTable[kk]['atomicname'],PeriodicTable[kk]['atomicnumber'],PeriodicTable[kk]['atomicmass'])
            return(PeriodicTable[kk]['atomicname'],PeriodicTable[kk]['atomicsymbol'],PeriodicTable[kk]['atomicnumber'],PeriodicTable[kk]['atomicmass'])
            break
        if PeriodicTable[kk]['atomicmass']==identifier:
            #print(PeriodicTable[kk]['atomicname'],PeriodicTable[kk]['atomicnumber'],PeriodicTable[kk]['atomicmass'])
            return(PeriodicTable[kk]['atomicname'],PeriodicTable[kk]['atomicsymbol'],PeriodicTable[kk]['atomicnumber'],PeriodicTable[kk]['atomicmass'])
            break



