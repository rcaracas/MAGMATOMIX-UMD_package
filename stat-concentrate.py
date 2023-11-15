#!/usr/bin/env python3

import sys,getopt,glob,re
import crystallography as cr


def molecularmass(formula,chemistry):
    indexes = [int(s) for s in re.findall('[0-9]+', formula)]
    elements = re.findall('[a-zA-Z]+',formula)
    massM = 0.0
    massA = 0.0
    occurs = 0.0
    totalatomsinmolec = 0.0
    for ielem in range(len(elements)):
        #print('Next element is ',elements[ielem])
        (atn,ats,atno,mass)=cr.Elements2rest(elements[ielem])
        massM = massM + float(mass) * float(indexes[ielem])
        totalatomsinmolec = totalatomsinmolec + float(indexes[ielem])
        if elements[ielem] == chemistry or chemistry=='_':
            massA = massA + float(mass) * float(indexes[ielem])
            occurs = float(indexes[ielem])#How many times any of the elements of chemistry are in the molecule
    return(massM,massA,occurs,totalatomsinmolec)


def concatstatfile(FileName,population,iround,clusterlength,lifetime,chemistry):
    print('\n ++++ reading file ++++',FileName)
    if(chemistry=='_'):
        ch='every element'
    else:
        ch=chemistry
    print('Looking for the chemistry of',ch)
    ff = open(FileName,'r')
    line = ff.readline()
    line = ff.readline()
    totaltime = 0.0                                 #defined as the sum of all lifetimes
    
    totaltimevapor = 0.0                            #defined as the sum of all gas lifetimes
    totalmassvapor = 0.0                            # defined as total mass of the vapor molecules x lifetimes
    totaltimecheminvap = 0.0                        # defined as no of chemistry.in.gas.molecules x lifetime
    totalatomsinvapor = 0.0                         # defined as no of atoms in gas molecules x lifetime
    totaltimemoleculeswithchemistryinvap = 0.0      # defined as sum of lifetimes of gas.molecules.with.chemistry
    totalmasscheminvap = 0.0                        # defined as no of chemistry.in.gas.molecules x chemistry.mass x lifetime

    totaltimeliquid = 0.0                           #defined as the sum of all liquid lifetimes
    totalmassliquid = 0.0                           # defined as total mass of the liquid molecules x lifetimes
    totaltimecheminliq = 0.0                        # defined as no of chemistry.in.liquid.molecules x lifetime
    totalatomsinliquid = 0.0                        # defined as no of atoms in liquid molecules x lifetime
    totaltimemoleculeswithchemistryinliq = 0.0      # defined as sum of lifetimes of liquid.molecules.with.chemistry
    totalmasscheminliq = 0.0                        # defined as no of chemistry.in.liquid.molecules x chemistry.mass x lifetime
    massM = 0.0
    massA = 0.0
    occurs = 0.0
    while True:
        line = ff.readline()
        if not line: break
        line=line.strip()
        entry=line.split(   )
        atomsinmolec=[]
        atomsextract=line.split('[')[1].split(']')[0].split(', ')
        for char in atomsextract:
            if(char.isdigit()):
                atomsinmolec.append(int(char))
        flagcluster = -1        #flag to check the existence of the cluster
        totaltime = totaltime + float(entry[3])
        if int(len(atomsinmolec)) <= clusterlength:          # cluster is in the gas
                                                    # clusters smaller than clusterlength are considered vapor
            totaltimevapor = totaltimevapor + float(entry[3])
            (massM,massA,occurs,totalatomsinmolec)=molecularmass(entry[0],chemistry)
            totalmassvapor = totalmassvapor + massM * float(entry[3])               #defined as total mass of the vapor molecules x lifetimes
            totalatomsinvapor = totalatomsinvapor + totalatomsinmolec * float(entry[3])
            if int(occurs) > 0:
                totaltimemoleculeswithchemistryinvap = totaltimemoleculeswithchemistryinvap + float(entry[3])           # defined as sum of lifetimes of gas.molecules.with.chemistry
                totalmasscheminvap = totalmasscheminvap + massA * float(entry[3])       #defined as no.of chemistry x lifetime x mass
                totaltimecheminvap = totaltimecheminvap + occurs * float(entry[3])      #defined as no.of chemistry x lifetime
            for ii in range(len(population)):
                if population[ii][0] == entry[0]:       #cluster already exists
                    population[ii][iround+2	] = entry[1]
                    flagcluster = 0
            if flagcluster == -1:
                population.append(['' for _ in range(len(population[0]))])
                population[len(population)-1][0] = entry[0]
                population[len(population)-1][iround+2] = entry[1]
                population[len(population)-1][1] = entry[3]
        else:                                   #cluster is in the liquid
            totaltimeliquid = totaltimeliquid + float(entry[3])
            (massM,massA,occurs,totalatomsinmolec)=molecularmass(entry[0],chemistry)
            totalatomsinliquid = totalatomsinliquid + totalatomsinmolec * float(entry[3])
            if entry[0].find(chemistry)>-1:
                    totaltimemoleculeswithchemistryinliq = totaltimemoleculeswithchemistryinliq + float(entry[1])
            totalmassliquid = totalmassliquid + massM * float(entry[3])
            totalmasscheminliq = totalmasscheminliq + massA * float(entry[3])
            totaltimecheminliq = totaltimecheminliq + occurs * float(entry[3])
            

    return (population,totaltime,totaltimevapor,totalmassvapor,totalatomsinvapor,totaltimecheminvap,totaltimemoleculeswithchemistryinvap, totalmasscheminvap,totaltimeliquid,totalmassliquid,totalatomsinliquid,totaltimecheminliq,totaltimemoleculeswithchemistryinliq,totalmasscheminliq)



def relativespecies(FileName,population,header):#not actually used
    ff = open(FileName,'w')
    ff.write(header)
    string = '\n'
    ff.write(string)
    SortedPopul = sorted(population)                 #SortedPopul = population sorted alphabetically over the first column, with the name of the clusters
    element = 'Me'
    oldelem = element
    currstart = 0                #line number at which a new series of central cations start to be listed
    instances = [0.0 for _ in range(len(population[0])-1)]            #stores the number of times clusters starting with a given element occur for all densities
    for ipop in range(len(SortedPopul)):
        if SortedPopul[ipop][0] != 'cluster':                #identifies the first element in each cluster
            cluster = SortedPopul[ipop][0]
            cluster = cluster.split('_')
            ilength = len(cluster[0])           #assumes all the string is only one element
            for ichar in range(1,len(cluster[0])):
                if cluster[0][ichar].isupper():
                    ilength = ichar             #stops at the first capital letter, corresponding to the new element in the cluster
                    break
            element = cluster[0][:ilength]      #reads to the first capital letter or to the end of the cluster
            
            if element == oldelem:      #takes into account all the clusters starting with the same element
                for icol in range(1,len(population[0])):
                    if len(SortedPopul[ipop][icol])>0:               #counts only the entries which exist (as some clusters don't exist at some densities)
                        instances[icol-1] = instances[icol-1] + float(SortedPopul[ipop][icol])
            else:                       #this is a new element
#                if ipop == 0:           #this is the very first time we analyze clusters
                    #                    string = SortedPopul[ipop][0]
                    #                    oldelement = cluster[0]
#                    for icol in range(1,len(population[0])):
#                        if len(SortedPopul[ipop][icol])>0:
                            #print('ipop ',ipop,' column ',icol,' SortedPopul ',SortedPopul[ipop][icol])
                            #instances[icol] = instances[icol] + float(SortedPopul[ipop][icol])
#                            instances[icol-1] = float(SortedPopul[ipop][icol])
#                        print(' for each density instances are:  ',instances)
#                else:
                if ipop > 0:        #th
                    #print('dealing first with the former element',oldelem)
                    #print('   summing up from ',currstart,' till ',ipop-1)
                    #print('   former values of the instances',instances)
                    string = '\n'
                    ff.write(string)
                    for iline in range(currstart,ipop):
                        #print('iline now to be printed is ',iline,' with ',SortedPopul[iline])
                        #print(instances)
                        string = SortedPopul[iline][0]
                        for icol in range(1,len(population[0])):
                            if len(SortedPopul[iline][icol])>0:
                                string = string + '\t' + str( float(SortedPopul[iline][icol]) / instances[icol-1] )
                            else:
                                string = string + '\t'
                        string = string +'\n'
                        ff.write(string)
                        string = SortedPopul[iline][0]
                    for icol in range(2,len(population[0])):
                        if len(SortedPopul[ipop][icol])>0:               #counts only the entries which exist (as some clusters don't exist at some densities)
                            instances[icol-1] = 	float(SortedPopul[ipop][icol])
                        else:
                            instances[icol-1] = 0.0

                    currstart = ipop
                oldelem = element
                for icol in range(2,len(population[0])):
                    if len(SortedPopul[ipop][icol])>0:
                        instances[icol-1] = float(SortedPopul[ipop][icol])

        else:
        #print('hit the last line of the file')
        #print('dealing first with the last element',oldelem)
        #print(' summing up from ',currstart,' till ',ipop-1)
        #print('current values of the instances',instances)
            string = '\n'
            ff.write(string)
            for iline in range(currstart,ipop):
                string = SortedPopul[iline][0]
                for icol in range(1,len(population[0])):
                    if len(SortedPopul[iline][icol])>0:
                        string = string +  '\t' + str( float(SortedPopul[iline][icol]) / instances[icol-1] )
                    else:
                        string = string + '\t'
                string = string +'\n'
                ff.write(string)
    ff.close()





def main(argv):
    logfile = 'cond.spec.'
    lifetime = 0
    clusterlength = 9999
    rings = 0
    chemistry = '_'
    try:
        opts, arg = getopt.getopt(argv,"hf:l:t:r:c:",["flogfile","lLength","tLifeTime","rRings","cChemistry"])
    except getopt.GetoptError:
        print ('stat-concentrate <to concentrate all the stat files into one multi-column file -f logfile -l max_length_gas_molecs> -t lifetime -r polymerization_type -c chemistry')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('stat-concentrate.py program to concentrate all the stat files into one multi-column file  -f logfile -l max_length_gas_molecs> -t lifetime  -r polymerization_type -c Chemistry')
            print ('Defaults are: -f cond.spec. -l 9999 -t 0 -r 0 -c _(all chemistry) ')
            sys.exit()
        elif opt in ("-f","-flogfile"):
            logfile = str(arg)
        elif opt in ("-t","--tLifeTime"):
            lifetime = int(arg)
        elif opt in ("-l","--lMaxLength"):
            clusterlength = float(arg)
        elif opt in ("-r","--rRings"):
            rings = int(arg)
        elif opt in ("-c","--cChemistry"):
            chemistry =str(arg)
    print('\nCondensing all speciation files into ',logfile)
    print('Considering species with min lifetime of ',lifetime)
    print('Maximum vapor species size is ',clusterlength)
    print('L/V partitioning for species containing ',chemistry)
    if rings == 0:
        list_of_files = glob.glob('./*.r0.popul.dat')
    elif rings == 1:
        list_of_files = glob.glob('./*.r1.popul.dat')
    else:
        print('Cluster type should be either 0 or 1, but this was  ',rings,'  This is not allowed')
        sys.exit()
    list_of_files = sorted(list_of_files)
    #print ('files are :',list_of_files)
    if len(list_of_files)>0:
        population = [['' for _ in range(len(list_of_files)+2)]]        #population = matrix containing all the clusters
                            # col 0: cluster name
                            # col 1: cluster length
                            # col 2+: cluster time in each file
                            
                            #the following varialbes contain total over all occurences of vapor/liquid with or without chemistry
                            #they are all defined for each population file
        totaltime = [0.0 for _ in range(len(list_of_files))]                #defined as the sum of all lifetimes
        
        totaltimevapor = [0.0 for _ in range(len(list_of_files))]           # defined as total vapor molecules x lifetimes
        totalmassvapor = [0.0 for _ in range(len(list_of_files))]           # defined as gas.molecules x lifetime x mass
        totalatomsinvapor = [0.0 for _ in range(len(list_of_files))]        # defined as no of atoms in gas molecules x lifetime
        totaltimecheminvap = [0.0 for _ in range(len(list_of_files))]       # defined as no of chemistry.in.gas.molecules x lifetime
        totaltimemoleculeswithchemistryinvap = [0.0 for _ in range(len(list_of_files))]     # defined as sum of lifetimes of gas.molecules.with.chemistry
        totalmasscheminvap = [0.0 for _ in range(len(list_of_files))]       # defined as no of chemistry.in.gas.molecules x chemistry.mass x lifetime

        totaltimeliquid = [0.0 for _ in range(len(list_of_files))]          #defined as the sum of all liquid lifetimes
        totalmassliquid = [0.0 for _ in range(len(list_of_files))]          # defined as liquid molecules x lifetimes x mass
        totalatomsinliquid = [0.0 for _ in range(len(list_of_files))]        # defined as no of atoms in liquid molecules x lifetime
        totaltimecheminliq = [0.0 for _ in range(len(list_of_files))]       # defined as no of chemistry.in.liquid.molecules x lifetime
        totaltimemoleculeswithchemistryinliq = [0.0 for _ in range(len(list_of_files))]       # defined as sum of lifetimes of liquid.molecules.with.chemistry
        totalmasscheminliq = [0.0 for _ in range(len(list_of_files))]       # defined as no of chemistry.in.liquid.molecules x chemistry.mass x lifetime
        
        totalsimulationvaportime = 0.0                          #the entire simulation time of vapor - relevant for several indep confgis at the same P-V-T
        totalsimulationliquidtime = 0.0                         #the entire simulation time of liquid - relevant for several indep confgis at the same P-V-T
        population[0][0] = 'cluster'
        for ifile in range(len(list_of_files)):
            shortfile = list_of_files[ifile]
            population[0][ifile+2] = shortfile.split('stat.dat')[0].split('outcar')[0]
        for ifile in range(len(list_of_files)):
            (population,totaltime[ifile],totaltimevapor[ifile],totalmassvapor[ifile],totalatomsinvapor[ifile],totaltimecheminvap[ifile],totaltimemoleculeswithchemistryinvap[ifile], totalmasscheminvap[ifile],totaltimeliquid[ifile],totalmassliquid[ifile],totalatomsinliquid[ifile],totaltimecheminliq[ifile],totaltimemoleculeswithchemistryinliq[ifile],totalmasscheminliq[ifile]) = concatstatfile(list_of_files[ifile],population,ifile,clusterlength,lifetime,chemistry)
            #print('Total time of chemical species in L+V is ',totaltime[ifile])
            #print('Total time of vapor species is ',totaltimevapor[ifile])
            #print('Total atoms in vapor is ',totalatomsinvapor[ifile])
            #print('Total time of liquid species is ',totaltimeliquid[ifile])
            #print('Total atoms in liquid is ',totalatomsinliquid[ifile])
            totalsimulationvaportime = totalsimulationvaportime + totaltimevapor[ifile]
            totalsimulationliquidtime = totalsimulationliquidtime + totaltimeliquid[ifile]
            if totaltimevapor[ifile] > 0:
                print('VAPOR')
                print('Total time ',chemistry,' atoms in vapor',totaltimecheminvap[ifile])
                print('Mole proportion of ',chemistry,' out of total vapor ',100*totaltimecheminvap[ifile]/totalatomsinvapor[ifile],' %')
                print('Mass proportion of ',chemistry,' out of total vapor mass',100*totalmasscheminvap[ifile]/totalmassvapor[ifile],' %')
                print('Total time molecules with ',chemistry,' is ',totaltimemoleculeswithchemistryinvap[ifile],' which is ',100*totaltimemoleculeswithchemistryinvap[ifile]/totaltimevapor[ifile],' % of total vapor species')
            else:
                print('No ',chemistry,' in the vapor phase')

            if totaltimeliquid[ifile] > 0:
                print('LIQUID')
                print('Total time ',chemistry,' atoms in liquid',totaltimecheminliq[ifile])
                print('Mole proportion of ',chemistry,' out of total liquid ',100*totaltimecheminliq[ifile]/totalatomsinliquid[ifile],' %')
                print('Mass proportion of ',chemistry,' out of total liquid mass',100*totalmasscheminliq[ifile]/totalmassliquid[ifile],' %')
                print('Total time molecules with ',chemistry,' is ',totaltimemoleculeswithchemistryinliq[ifile],' which is ',100*totaltimemoleculeswithchemistryinliq[ifile]/totaltimeliquid[ifile],' % of total liquid species')
            else:
                print('No ',chemistry,' in the liquid phase')

            if totaltimecheminvap[ifile]>0:
                print('L/V PARTITIONING')
                print('L/V mole partitioning of ',chemistry,' is ',totaltimecheminliq[ifile]/totaltimecheminvap[ifile])
            else:
                print("Every",chemistry,"is liquid.")

            
    else:
        print ('There are no *.stat.dat files to concatenate')
        sys.exit()
    
    #Print all sort of stat output files
    outfile = logfile + 'abso.-t' + str(lifetime) + '.-L' + str(clusterlength) + '.-R' + str(rings) + '.dat'
    percfile = logfile + 'perc.-t' + str(lifetime) + '.-L' + str(clusterlength) + '.-R' + str(rings) + '.dat'
    vapfile = logfile + 'vapor.-t' + str(lifetime) + '.-L' + str(clusterlength) + '.-R' + str(rings) + '.-C=' + str(chemistry) + '.dat'
    
        #all ABSOLUTE times for all vapor species
    ff = open(outfile,'w')
    string = 'TotalTimes\tSize\t'
    for icol in range(len(list_of_files)):
        string = string + str(totaltime[icol]) +'\t'
    for icol in range(2,len(population[0])):
        string = string + '%in' + population[0][icol] +'\t'
    string = string + 'CumulativeTime\tCumulativePercent\t'
    string = string + '\n'
    ff.write(string)

    for ipop in range(1,len(population)):
        #print (population[ipop])
        string = ''
        cumultime = 0.0
        for icol in range(len(population[ipop])):
            string = string + str(population[ipop][icol]) +'\t'
        for icol in range(2,len(population[ipop])):
            if len(population[ipop][icol])>0:                        # print percentage of each of the species
                cumultime = cumultime + float(population[ipop][icol])
                string = string + str(100*float(population[ipop][icol])/totaltime[icol-2]) + '\t'
            else:
                string = string + '\t'
        string = string + str(cumultime) + '\t' + str(100*cumultime/totalsimulationvaportime) + '\t'
        string = string +'\n'
        ff.write(string)
    ff.close()
    
        #all VAPORs containing chemistry
    ff = open(vapfile,'w')
    string = 'Molecule\tSize\t'
    for icol in range(len(totaltime)):
        string = string + str(totaltimevapor[icol]) +'\t'
    for icol in range(len(totaltime)):
        string = string + '%inFile' + str(icol) + '\t'
    string = string + 'CumulativeTime\tCumulativePercent\t'
    string = string + '\n'
    ff.write(string)
    for ipop in range(len(population)):
        string = ''
        cumultime = 0.0
        if population[ipop][0].find(chemistry)>-1:
            for icol in range(len(population[ipop])):
                string = string + str(population[ipop][icol]) + '\t'     # print absolute times
            for icol in range(2,len(population[ipop])):
                if len(population[ipop][icol])>0:                        # print percentage of each of the species
                    cumultime = cumultime + float(population[ipop][icol])
                    string = string + str(100*float(population[ipop][icol])/totaltimevapor[icol-2]) + '\t'
                else:
                    string = string + '\t'
            string = string + str(cumultime) + '\t' + str(100*cumultime/totalsimulationvaportime) + '\t'
            string = string +'\n'
        ff.write(string)
    ff.close
    
        #PERCENTAGE times
    header = ''
    for icol in range(len(population[0])):
        header = header + population[0][icol] + '\t'
    header = header + '\n'

if __name__ == "__main__":
    main(sys.argv[1:])
