#!/usr/bin/env python

import sys,getopt,numpy,os.path,math,glob,re
import crystallography as cr
from subprocess import call


def concatstatfile(FileName,population,iround,clusterlength,lifetime):
    print(' ++++ reading file ++++',FileName)
    #print ('in concat round',iround,' with pop of length ',len(population))
    ff = open(FileName,'r')
    line = ff.readline()
    line = ff.readline()
    #    line = ff.readline()
    while True:
        line = ff.readline()
        if not line: break
        line=line.strip()
        entry=line.split()
        flagcluster = -1        #flag to check the existence of the cluster
        #print ('current population is:',population)
        #print ('entries are ',entry)
        if int(entry[3]) <= clusterlength:         #only clusters smaller than clusterlength are considered, to be able to separate gas (i.e. small) from liquid (i..e large polymers) species
            if int(entry[1]) > lifetime:            #only clusters with a total life time larger than lifetime are considered
            #print('reading current line ',entry[0],'--',entry[1],'--',entry[2])
                for ii in range(len(population)):
                #print('previous ',population[ii][0],' compared to new ',entry[0])
                    if population[ii][0] == entry[0]:       #cluster already exists
                        population[ii][iround+1	] = entry[1]
                        flagcluster = 0
            #print('current flagcluster',flagcluster)
                if flagcluster == -1:
                    population.append(['' for _ in range(len(population[0]))])
                    population[len(population)-1][0] = entry[0]
                    population[len(population)-1][iround+1] = entry[1]
                #print('     append new cluster',population[len(population)-1])
    return (population)


def relativespecies(FileName,population,header):
    #for iline in range(len(population)):
    #    print ('line is ',population[iline])
    ff = open(FileName,'w')
    ff.write(header)
    string = '\n'
    ff.write(string)
    sorpop = sorted(population)                 #sorpop = population sorted alphabetically over the first column, with the name of the clusters
    #print('sorpop is ',sorpop,' with length ',len(sorpop),'\n')
    element = 'Me'
    oldelem = element
    currstart = 0                #line number at which a new series of central cations start to be listed
    #print('\n\n')
    #instances = [0.0 for _ in range(len(population[0]))]            #stores the number of times clusters starting with a given element occur for all densities
    instances = [0.0 for _ in range(len(population[0])-1)]            #stores the number of times clusters starting with a given element occur for all densities
    #transprec = [[0.0 for _ in range(len(population[0]))] for _ in range(len(population))]      #stores the percentage of each cluster, normalized to 1 to the first cation
    for ipop in range(len(sorpop)):
        #print()
        #print (' treating current row no. ',ipop,'  :  ',sorpop[ipop])
        if sorpop[ipop][0] != 'cluster':                #identifies the first element in each cluster
            cluster = sorpop[ipop][0]
            #print('current cluster is',cluster)
            cluster = cluster.split('_')
            #print('cluster after _ splitting',cluster,' with length ',len(cluster[0]))
            ilength = len(cluster[0])           #assumes all the string is only one element
            for ichar in range(1,len(cluster[0])):
                #print('\t\t\t current letter no. ichar ',ichar,' is ',cluster[0][ichar])
                if cluster[0][ichar].isupper():
                    ilength = ichar             #stops at the first capital letter, corresponding to the new element in the cluster
                    break
            element = cluster[0][:ilength]      #readds to the first capital letter or to the end of the cluster
            #print('first element of the cluster is ',element)
            
            if element == oldelem:      #takes into account all the clusters starting with the same element
                #print(' same element ',element)
                for icol in range(1,len(population[0])):
                    if len(sorpop[ipop][icol])>0:               #counts only the entries which exist (as some clusters don't exist at some densities)
                        instances[icol-1] = instances[icol-1] + float(sorpop[ipop][icol])
                #print(' for each density instances are:  ',instances)
            else:                       #this is a new element
                #print(' new element ',element)
#                if ipop == 0:           #this is the very first time we analyze clusters
                    #                    string = sorpop[ipop][0]
                    #                    oldelement = cluster[0]
#                    for icol in range(1,len(population[0])):
#                        if len(sorpop[ipop][icol])>0:
                            #print('ipop ',ipop,' column ',icol,' sorpop ',sorpop[ipop][icol])
                            #instances[icol] = instances[icol] + float(sorpop[ipop][icol])
#                            instances[icol-1] = float(sorpop[ipop][icol])
#                        print(' for each density instances are:  ',instances)
#                else:
                if ipop > 0:        #th
                    #print('dealing first with the former element',oldelem)
                    #print('   summing up from ',currstart,' till ',ipop-1)
                    #print('   former values of the instances',instances)
                    string = '\n'
                    ff.write(string)
                    for iline in range(currstart,ipop):
                        #print('iline now to be printed is ',iline,' with ',sorpop[iline])
                        #print(instances)
                        string = sorpop[iline][0]
                        for icol in range(1,len(population[0])):
                            if len(sorpop[iline][icol])>0:
#                                string = string + '\t' + str(sorpop[iline][icol]) + '\t' + str( float(sorpop[iline][icol]) / instances[icol-1] )
                                string = string + '\t' + str( float(sorpop[iline][icol]) / instances[icol-1] )
                            #transprec[icol][iline] = float(sorpop[iline][icol]) / instances[icol]
                            else:
                                string = string + '\t'
                        string = string +'\n'
                        #print('\tline no. ',iline,' oldelement ',oldelem,' percentages: ',string)
                        ff.write(string)
                        string = sorpop[iline][0]
                    for icol in range(1,len(population[0])):
                        if len(sorpop[ipop][icol])>0:               #counts only the entries which exist (as some clusters don't exist at some densities)
                            instances[icol-1] = 	float(sorpop[ipop][icol])
                        else:
                            instances[icol-1] = 0.0

                    currstart = ipop
                oldelem = element
                for icol in range(1,len(population[0])):
                    if len(sorpop[ipop][icol])>0:
                        #print('ipop ',ipop,' column ',icol,' sorpop ',sorpop[ipop][icol])
                        #instances[icol] = instances[icol] + float(sorpop[ipop][icol])
                        instances[icol-1] = float(sorpop[ipop][icol])
                #print(' for each density instances are:  ',instances)

    #print('updated value of the instances ', instances)
        else:
        #print('hit the last line of the file')
        #print('dealing first with the last element',oldelem)
        #print(' summing up from ',currstart,' till ',ipop-1)
        #print('current values of the instances',instances)
            string = '\n'
            ff.write(string)
            for iline in range(currstart,ipop):
            #print('iline now is ',iline,' with ',sorpop[iline])
            #print(instances)
                string = sorpop[iline][0]
                for icol in range(1,len(population[0])):
                    if len(sorpop[iline][icol])>0:
#                        string = string + '\t' + str(sorpop[iline][icol]) +  '\t' + str( float(sorpop[iline][icol]) / instances[icol-1] )
                        string = string +  '\t' + str( float(sorpop[iline][icol]) / instances[icol-1] )
                    else:
                        string = string + '\t'
                string = string +'\n'
            #print('line no. ',iline,' oldelement ',oldelem,' percentages: ',string)
                ff.write(string)
    ff.close()

#transfile = FileName[:-3] + 'transposed.dat'
#    ff.open(transfile)




def main(argv):
    logfile = 'cond.spec.'
    lifetime = 0
    clusterlength = 50
    rings = 0
    try:
        opts, arg = getopt.getopt(argv,"hf:l:t:r:",["flogfile","lLength","tLifeTime","rRings"])
    except getopt.GetoptError:
        print ('stat-concentrate <to concentrate all the stat files into one multi-column file -f logfile -l max_length_gas_molecs> -t lifetime -r polymerization_type')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('stat-concentrate.py program to concentrate all the stat files into one multi-column file  -f logfile -l max_length_gas_molecs> -t lifetime  -r polymerization_type')
            sys.exit()
        elif opt in ("-f","-flogfile"):
            logfile = str(arg)
        elif opt in ("-t","--tLifeTime"):
            lifetime = int(arg)
        elif opt in ("-l","--lMaxLength"):
            clusterlength = float(arg)
        elif opt in ("-r","--rRings"):
            rings = int(arg)
    print('Condensing all speciation files into ',logfile)
    print('Considering species with min lifetime of ',lifetime)
    print('Maximum species size is ',clusterlength)
    if rings == 0:
        list_of_files = glob.glob('./*.r0.stat.dat')
    elif rings == 1:
        list_of_files = glob.glob('./*.r1.stat.dat')
    else:
        print('Cluster type should be either 0 or 1, but this was  ',rings,'  This is not allowed')
        sys.exit()
    list_of_files = sorted(list_of_files)
    #print ('files are :',list_of_files)
    if len(list_of_files)>0:
        population = [['' for _ in range(len(list_of_files)+1)]]                #population = matrix containing all the clusters
                                                                                #  first column is the cluster name, the other columns are (no. of) the files
        population[0][0] = 'cluster'
        for ifile in range(len(list_of_files)):
            shortfile = list_of_files[ifile]
            population[0][ifile+1] = shortfile.split('stat.dat')[0].split('outcar')[0]
        for ifile in range(len(list_of_files)):
            population = concatstatfile(list_of_files[ifile],population,ifile,clusterlength,lifetime)
    else:
        print ('There are no *.stat.dat files to concatenate')
        sys.exit()
    outfile = logfile + 'abso.-T' + str(lifetime) + '.-L' + str(clusterlength) + '.-R' + str(rings) + '.dat'
    percfile = logfile + 'perc.-T' + str(lifetime) + '.-L' + str(clusterlength) + '.-R' + str(rings) + '.dat'
    ff = open(outfile,'w')
    for ipop in range(len(population)):
        #print (population[ipop])
        string = ''
        for icol in range(len(population[ipop])):
            string = string + population[ipop][icol] + '\t'
        string = string + '\n'
        ff.write(string)
    ff.close()
    header = ''
    for icol in range(len(population[0])):
        header = header + population[0][icol] + '\t'
    header = header + '\n'
    #relativespecies(percfile,population[1:len(population)],header)
    relativespecies(percfile,population,header)

if __name__ == "__main__":
    main(sys.argv[1:])
