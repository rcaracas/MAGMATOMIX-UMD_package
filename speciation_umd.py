#!/usr/bin/env python


###
##AUTHORS: RAZVAN CARACAS
###

import sys,getopt,numpy,os.path,math,itertools
import crystallography as cr
import umd_process as up
from subprocess import call

def analysis_clusters(clusters,MyCrystal,ligands,minlife,Nsteps,FileName,rings):
    #this functions builds the dictionary with all the clusters
    ###labels = {'',[]}
    ###labels = {'Mg2Si4O9',['4137141516171891929394956979899',start,end,duration]}
    ###print ('checking: indexes',indexes)
    population = {}
    #population['clusterindex']={} #'clustername':'', 'formula':'', 'word':'', 'begin':0, 'end':0, 'lifetime':1,'composition':[]}
    #print ('initialization of the population',population,'size of clusters',len(clusters))
    for ii in range(len(clusters)):
        #print ('snapshot no. ',ii)
        #labels.append([[' ' for _ in range(2)] for _ in range(len(clusters[ii]))])
        for jj in range(len(clusters[ii])):
            #print('dealing with cluster',ii,'which is ',clusters[ii][jj])
            word = ''
            molecule = ''
            indexes = [0 for _ in range(len(MyCrystal.elements))]
            #formulalong = ['' for _ in range(len(clusters[ii][jj]))]
            #print (' current ',ii,jj,' cluster  is ', clusters[ii][jj])
            #print ('size of cluster',len(clusters[ii][jj]))
            for kk in range(len(clusters[ii][jj])):
                word = word + str(clusters[ii][jj][kk]) + '-'
                indexes[MyCrystal.typat[clusters[ii][jj][kk]]] +=1
            for kk in range(len(MyCrystal.elements)):
                if indexes[kk]>0:
                    molecule = molecule + str(MyCrystal.elements[kk]) + '_' + str(indexes[kk])
            clustername = molecule + '_' + word
            clusterindex = clustername + '_' + str(ii)
            #print ('     ===>  current cluster to parse is ',clustername,' with index',clusterindex)
            flagalive = 0
                #if len(population)>1:
            for ll in population.keys():
                #print ('comparing currect cluster ',population[ll]['clustername']) #,' with ',population[kk][clustername])
                if clustername == population[ll]['clustername']:
                    #print ('current step',ii,'end of population',population[ll]['end'])
                    if ii - population[ll]['end'] == 1:
                        population[ll]['end'] += 1
                        population[ll]['lifetime'] += 1
                        
                        #print (' cluster already exists since',population[ll]['begin'])
                        flagalive = 1
            if flagalive == 0:
                #print (' cluster doesnt exist')
                population[clusterindex]={'clustername':clustername, 'formula':molecule, 'word':word, 'begin':ii, 'end':ii, 'lifetime':1,'composition':clusters[ii][jj]}
#print ('new cluster added: ',population[clusterindex][molecule],population[clusterindex][begin],population[clusterindex][end])
#print ('new cluster added: ',population[clusterindex])
    FileAll = FileName + '.r' + str(rings) + '.popul.dat'
    #print ('Population will be written in ',FileAll,' file')
    fa = open(FileAll,'a')
    newstring = 'Formula\tBegin\tEnd\tLifetime\t[Composition]\n'
    fa.write(newstring)
    #print('Formula - Begin : End : Lifetime - Composition')
    for kk in population.keys():
        #newstring = population[kk]['formula'] + '\t' + str(population[kk]['begin']*Nsteps) + '\t' + str(population[kk]['end']*Nsteps) + '\t' + str(population[kk]['lifetime']*Nsteps) + '\t' + str(population[kk]['composition']) + '\n'
        if population[kk]['lifetime'] > minlife/float(Nsteps):
            newstring = population[kk]['formula'] + '\t' + str(population[kk]['begin']*Nsteps) + '\t' + str(population[kk]['end']*Nsteps) + '\t' + str(population[kk]['lifetime']*Nsteps) + '\t' + str(population[kk]['composition']) + '\n'
    #print(population[kk]['formula'],population[kk]['begin']*Nsteps,population[kk]['end']*Nsteps,population[kk]['lifetime']*Nsteps,population[kk]['composition'])
            fa.write(newstring)
    
    statclusters = [['',0,0]]
    #print('length of statclusters is',len(statclusters))
    for kk in population.keys():
        #print('treating cluster',population[kk]['formula'])
        ll = 0
        flagnewclust = 0
        while flagnewclust == 0:
            #print('comparing to cluster',statclusters[ll][0])
            if population[kk]['formula'] == statclusters[ll][0]:
                statclusters[ll][1] += population[kk]['lifetime']*Nsteps
                #print('current population is',statclusters[ll][1])
                flagnewclust = 1
            else:
                if ll == len(statclusters)-1:
                    #print('adding new cluster',population[kk]['formula'])
                    statclusters.append([population[kk]['formula'],population[kk]['lifetime']*Nsteps,len(population[kk]['composition'])])
                    flagnewclust = 1
                else:
                    ll +=1

    totalpop = 0
    for ll in range(1,len(statclusters)):
        totalpop += statclusters[ll][1]
    FileStat = FileName + '.r' + str(rings) + '.stat.dat'
#print ('Statistics will be written in ',FileStat,' file')
    fs = open(FileStat,'a')
    newstring = 'Cluster\tTime\tPercent\n'
    fs.write(newstring)
    for ll in range(1,len(statclusters)):
        newstring = statclusters[ll][0] + '\t' + str(statclusters[ll][1]) + '\t' + str(float(statclusters[ll][1])/float(totalpop)) + '\t' + str(statclusters[ll][2]) + '\n'
        fs.write(newstring)


def neighboring(BooleanMap,ligands,iatom,newcluster):
    newcluster.append(ligands[iatom])
    if iatom < len(ligands):
        #        for jatom in range(ligands[iatom+1],len(ligands)):
        for jatom in range(len(ligands)):
            if BooleanMap[ligands[iatom]][ligands[jatom]]==1:
                #print('bond between',ligands[iatom],ligands[jatom])
                BooleanMap[ligands[iatom]][ligands[jatom]]=0
                BooleanMap[ligands[jatom]][ligands[iatom]]=0
                #print ('next atom',ligands[jatom])
                newcluster=neighboring(BooleanMap,ligands,jatom,newcluster)
    return newcluster
def clustering(BooleanMap,ligands):
    #print('     clustering: start')
    #print ('ligands:',ligands)
    icluster = -1
    neighbors = []
    i=[]
    monogas = []
    for iatom in range(len(ligands)):
        #print('sum of bonds for atom ',iatom,ligands[iatom],sum(BooleanMap[ligands[iatom]]))
        #for jatom in range(len(ligands)):
        #    if BooleanMap[ligands[iatom]][ligands[jatom]] == 1:
        #        print ('bonded to atom ',jatom)
        if sum(BooleanMap[ligands[iatom]]) == 0:            #this part has to be treated first, as after clustering a lot of bonds are removed from BooleanMap
                                                            # in the neiboring routine
                                                            # and atoms that are bonded would appear as a monoatomic gas
#            monogas.append(ligands[iatom])
            neighbors.append([ligands[iatom]])
    for iatom in range(len(ligands)):
        if sum(BooleanMap[ligands[iatom]]) > 0:
            #print('current atom',iatom,ligands[iatom],BooleanMap[ligands[iatom]],sum(BooleanMap[ligands[iatom]]))
            newcluster = []
            newcluster1 = neighboring(BooleanMap,ligands,iatom,newcluster)
            #print(' end of cluster. last was: ',newcluster1)
            #print('    same atom after removing bonds',iatom,ligands[iatom],BooleanMap[ligands[iatom]],sum(BooleanMap[ligands[iatom]]))
            newcluster1.sort()
            newcluster = [i[0] for i in itertools.groupby(newcluster1)]
            if len(newcluster)>1:
                neighbors.append(newcluster)
    return neighbors
def clusteringnorings(BooleanMap,centralatoms,outeratoms):
    #print ('ligands:',ligands)
    icluster = -1
    neighbors = []
    for iatom in range(len(centralatoms)):
        newcluster = [centralatoms[iatom]]
        for jatom in range(len(outeratoms)):
            if BooleanMap[centralatoms[iatom]][outeratoms[jatom]]==1:
                newcluster.append(outeratoms[jatom])
        if len(newcluster)>1:
            neighbors.append(newcluster)
    return neighbors


def ComputeBondMapInput(MyCrystal,MySnapshot,ligands,BondTable):
    #Builds the graph of the interatomic bonding
    #print('in ComputeBondMapInput')
    #print('no of atoms = ',MyCrystal.natom)
    #print('unit cell parameters = ',MySnapshot.acell)
    #print('dimension of ligands :',len(ligands))
    BondMap=[[0.0 for _ in range(MyCrystal.natom)] for _ in range(MyCrystal.natom)]
    BooleanMap=[[0 for _ in range(MyCrystal.natom)] for _ in range(MyCrystal.natom)]
    for iatom in range (len(ligands)):
        for jatom in range(iatom+1,len(ligands)):
            #compute distances along x, y and z between the two atoms
            #these lines are taken and adapted from the gofr function
            dx = MySnapshot.atoms[ligands[jatom]].xcart[0] - MySnapshot.atoms[ligands[iatom]].xcart[0]
            dy = MySnapshot.atoms[ligands[jatom]].xcart[1] - MySnapshot.atoms[ligands[iatom]].xcart[1]
            dz = MySnapshot.atoms[ligands[jatom]].xcart[2] - MySnapshot.atoms[ligands[iatom]].xcart[2]
            #if these distances are too large, then we correct them by translating one of the two atoms of one unit cell
            #if we want to work on a non orthogonal cell, then these corrections need to be changed
            #if (dx < -maxlength/2): dx = dx + maxlength  #these two lines are replaced by the one below: if |dx|>acell/2 then int(...) gives +/- 1
            #if (dx >  maxlength/2): dx = dx - maxlength
            #dx = dx-MyCrystal.acell[0]*int(dx/(0.5*MyCrystal.acell[0]))
            #dy = dy-MyCrystal.acell[1]*int(dy/(0.5*MyCrystal.acell[1]))
            #dz = dz-MyCrystal.acell[2]*int(dz/(0.5*MyCrystal.acell[2]))
            valx = min(dx**2, (MyCrystal.acell[0]-dx)**2, (MyCrystal.acell[0]+dx)**2)
            valy = min(dy**2, (MyCrystal.acell[1]-dy)**2, (MyCrystal.acell[1]+dy)**2)
            valz = min(dz**2, (MyCrystal.acell[2]-dz)**2, (MyCrystal.acell[2]+dz)**2)
            distij = valx + valy + valz
            if (distij < BondTable[MyCrystal.typat[ligands[iatom]]][MyCrystal.typat[ligands[jatom]]]): #if the distance between the two atoms is below the cutoff distance between the two atoms
                BondMap[ligands[iatom]][ligands[jatom]] = distij
                BooleanMap[ligands[iatom]][ligands[jatom]] = 1
                BondMap[ligands[jatom]][ligands[iatom]] = distij
                BooleanMap[ligands[jatom]][ligands[iatom]] = 1
        #print (' for atom ',iatom,' bonding scheme is ',BooleanMap[ligands[iatom]],' with total number of sum ',sum(BooleanMap[ligands[iatom]]))
    return (BondMap,BooleanMap)


def read_inputfile(InputFile,MyCrystal,ClusterAtoms):
    BondTable = [[0.0 for _ in range(MyCrystal.ntypat)] for _ in range(MyCrystal.ntypat)]
    nolines = 0
    with open(InputFile) as ff:
        for line in ff:
            line=line.strip()
            entry=line.split()
            if (len(entry)==3):
                nolines +=1
                for ii in range(MyCrystal.ntypat):
                    if MyCrystal.elements[ii]==entry[0]:
                        for jj in range(MyCrystal.ntypat):
                            if MyCrystal.elements[jj]==entry[1]:
                                BondTable[ii][jj]=float(entry[2])*float(entry[2])
                                BondTable[jj][ii]=float(entry[2])*float(entry[2])
    if nolines < len(ClusterAtoms)*len(ClusterAtoms):
        print ('WARNING: missing bonds in the input file. There are ',nolines,'bonds defined instead of ',len(ClusterAtoms)*len(ClusterAtoms))
    if nolines > len(ClusterAtoms)*len(ClusterAtoms):
        print ('WARNING: there are too many bonds in the input file. There are ',nolines,'bonds defined instead of ',len(ClusterAtoms)*len(ClusterAtoms))
            #    bondstoprint = BondTable
            #for ii in range(MyCrystal.ntypat):
            #for jj in range(MyCrystal.ntypat):
            #bondstoprint[ii][jj] = math.sqrt(BondTable[ii][jj])
            #print ('Bond table is: ')
            #for ii in range(MyCrystal.ntypat):
            #print (bondstoprint[ii])
            #return BondTable
    print ('Bond table is:')
    for ii in range(MyCrystal.ntypat):
        print (BondTable[ii])
    return BondTable


def main(argv):
    up.headerumd()
    UMDname='output.umd.dat'
    Nsteps = 1
    InitialStep = 0
    ClusterAtoms = []
    Cations = []
    Anions = []
    maxlength = 3.0
    InputFile = ''
    minlife = 5
    rings = 1
    header = ''
    try:
        opts, arg = getopt.getopt(argv,"hf:s:l:c:a:m:i:r:",["fUMDfile","sSampling_Frequency","lMaxLength","cCations","aAnions","mMinlife","iInputFile","rRings"])
    except getopt.GetoptError:
        print ('speciation.py -f <UMD_filename> -s <Sampling_Frequency> -l <MaxLength> -c <Cations> -a <Anions> -m <MinLife> -i <InputFile> -r <Rings>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('speciation.py program to compute bonding maps and identify speciation')
            print ('speciation.py -f <UMD_filename> -s <Sampling_Frequency> -l <MaxLength> -c <Cations> -a <Anions> -m <MinLife>  -i <InputFile> -r <Rings>')
            print ('  default values: -f output.umd.dat -s 1 -l 3.0 -m 0 -r 1)
            print (' the input file contains the bond lengths for the different atom pairs. \n the values overwrite the option -l')
            print (' rings = 1 default, polymerization, all anions and cations bond to each other; rings = 0 only individual cation-anion groups')
            sys.exit()
        elif opt in ("-f", "--fUMDfile"):
            UMDname = str(arg)
            header = header + 'FILE: -f=' + UMDname
        elif opt in ("-s","--sNsteps"):
            Nsteps = int(arg)
            header = header + ' -s=' + arg
            print('Will sample the MD trajectory every ',Nsteps,' steps')
        elif opt in ("-l","--lMaxLength"):
            maxlength = float(arg)
            print ('Maximum bonding length is',maxlength)
            header = header + ' -l=' + str(maxlength)
        elif opt in ("-m","--mMinlife"):
            minlife = float(arg)
        elif opt in ("-c","--Cations"):
            header = header + arg
            Cations = arg.split(",")
            #print ('Cation list is: ',Cations)
        elif opt in ("-a","--Anions"):
            header = header + arg
            Anions= arg.split(",")
            #print ('Anion list is: ',Anions)
        elif opt in ("-i", "--iInputFile"):
            InputFile = str(arg)
            header = header + ' -i=' + InputFile
            print ('Bonding cutoffs to be read from file ',InputFile)
        elif opt in ("-r","--rRings"):
            rings = int(arg)
            header = header + ' -r=' + arg
            if rings == 0:
                print ('Calculation of non-polymerized coordination polyhedra')
            elif rings == 1:
                print ('Calculation of polymerized coordination')
            else :
                print ('Undefined calculation')


    if not (os.path.isfile(UMDname)):
        print ('the UMD files ',UMDname,' does not exist')            
        sys.exit()

    for ii in range(len(Cations)):
        ClusterAtoms.append(Cations[ii])
    for ii in range(len(Anions)):
        ClusterAtoms.append(Anions[ii])
    if rings == 0:
        print('searching for cations',Cations)
        print('surrounded by anions ',Anions)
    else:
        print('all atoms bonding:',ClusterAtoms)
        
                  
#writes the header of the files containing eventually the clusters
    header = header + '\n'
    FileAll = UMDname + '.r' + str(rings) + '.popul.dat'
    print ('Population will be written in ',FileAll,' file')
    fa = open(FileAll,'w')
    fa.write(header)
    fa.close()
    FileStat = UMDname + '.r' + str(rings) + '.stat.dat'
    print ('Statistics will be written in ',FileStat,' file')
    fa = open(FileStat,'w')
    fa.write(header)
    fa.close()

                  
#reading the xc art coordinates of the atoms from the UMD file. it uses the read_xcart (i.e. only xcart) function from the umd_process library
    MyCrystal = cr.Lattice()
    AllSnapshots = [cr.Lattice]
    (MyCrystal,AllSnapshots,TimeStep)=up.read_xcart(UMDname)
    #print('checks after reading the umd file')
    #print('no of atoms = ',MyCrystal.natom)


#reading the cutoff radii for the bonds
    BondTable = [[maxlength*maxlength for _ in range(MyCrystal.ntypat)] for _ in range(MyCrystal.ntypat)]
    if len(InputFile)>0:
        BondTable = read_inputfile(InputFile,MyCrystal,ClusterAtoms)
#print ('unique bondlength, with square',maxlength)

#defining the ligands, the coordinated and the coordinating atoms
    centralatoms = []
    outeratoms = []
    ligands = []
    for iatom in range(MyCrystal.natom):
        for jatom in range(len(ClusterAtoms)):
            if MyCrystal.elements[MyCrystal.typat[iatom]]==ClusterAtoms[jatom]:
                ligands.append(iatom)           #contains the list with the index of the ligand atoms from the 0 ... natom
    for iatom in range(MyCrystal.natom):
        for jatom in range(len(Cations)):
            if MyCrystal.elements[MyCrystal.typat[iatom]]==Cations[jatom]:
                centralatoms.append(iatom)           #contains the list with the index of the central atoms from the 0 ... natom
        for jatom in range(len(Anions)):
            if MyCrystal.elements[MyCrystal.typat[iatom]]==Anions[jatom]:
                outeratoms.append(iatom)           #contains the list with the index of the coordinating atoms from the 0 ... natom
    print('All ligands are: ',ligands)
    print('Central atoms are :',centralatoms)
    print('Coordinating atoms are :',outeratoms)

#span the entire trajectory and analyze only the Nsteps snapshots
    clusters = []
    for istep in range(0,len(AllSnapshots),Nsteps):
        #print('analyzing new simulation snapshot, step no. ',istep)
#first build the bonding maps
        BondMap=[[0.0 for _ in range(MyCrystal.natom)] for _ in range(MyCrystal.natom)]
        BooleanMap=[[0 for _ in range(MyCrystal.natom)] for _ in range(MyCrystal.natom)]
        (BondMap,BooleanMap) = ComputeBondMapInput(MyCrystal,AllSnapshots[istep],ligands,BondTable)
        if rings == 1:
            clusters.append(clustering(BooleanMap,ligands))
        elif rings == 0:
            clusters.append(clusteringnorings(BooleanMap,centralatoms,outeratoms))
        #print('number of identified clusters ',len(clusters))
        #print(' clusters at this point:\n',clusters)
        else:
            print ('value of rings = ',rings,' is not allowed. Only 0 (polymers) or 1 (coordinating polyhedra) are allowed',)
            sys.exit()
    analysis_clusters(clusters,MyCrystal,ligands,minlife,Nsteps,UMDname,rings)

#        read_umd(UMDname,Nsteps,maxlength,Cations,Anions,ClusterAtoms,minlife,InputFile,rings)
#the read_umd function was replaced with the automatic one plus a bunch of separate funtions
#distance calculations come from gofr



if __name__ == "__main__":
    main(sys.argv[1:])


