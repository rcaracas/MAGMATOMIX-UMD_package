#!/usr/bin/env python3
###
##AUTHORS: RAZVAN CARACAS
###

import sys,getopt,os.path,re
import numpy as np
import umd_processes_fast as umdpf
from msd_umd import msd_fft


def msd_trajectory(trajectory):# partialmsd[x] is the mean square displacement between positions i and i+x, averaged over all the origins i
    partialmsd = msd_fft(np.array(trajectory)).sum(axis=-1)
    return partialmsd[:len(trajectory)//2]                                                     #lags up to half the lifetime, as in msd_umd.py


def ana_popul(MyCrystal,AbsXcart,TimeStep,POPname,MinLife,ClusterMaxSize,Nsteps):
    print ('in analysis of the population')
    print ('reading cluster population from ',POPname,' file')
    print ('initial total no of snapshots:',len(AbsXcart))
    masses = np.array([float(MyCrystal.masses[MyCrystal.typat[iatom]]) for iatom in range(MyCrystal.natom)])
    population = {}                                                                             #cluster name -> list of the MSDs of its individual instances
    ff = open(POPname,'r')
    ff.readline()
    ff.readline()
    while True:
        line = ff.readline()
        if not line: break
        line = re.sub(r'[\[,\]]', '',line)
        line = line.strip()
        entry = line.split()
        if len(entry) < 5: continue
        atoms = [int(x) for x in entry[4:]]
        if len(atoms) > ClusterMaxSize: continue
        begin = int(entry[1])
        end = int(entry[2])
        #the Begin/End steps of the popul file are counted in its own sampling ; the lifetime column (fs) gives the conversion to UMD snapshots
        ratio = max(1, int(round(float(entry[3]) / ((end-begin+1)*TimeStep))))
        first = begin*ratio
        last = min(end*ratio, len(AbsXcart)-1)
        if last - first < MinLife: continue
        snapshots = range(first, last+1, Nsteps)
        clusteratoms = AbsXcart[snapshots][:,atoms,:]                                           #(time, atom, xyz)
        trajectory = (clusteratoms*masses[atoms][None,:,None]).sum(axis=1)/masses[atoms].sum()  #center of mass of the current cluster
        partialmsd = msd_trajectory(trajectory)                                                 #calls the MSD on the trajectory
        if len(partialmsd) > 1:                                                                 #there must be at least two steps for diffusion to be added
            population.setdefault(entry[0],[]).append(partialmsd)
    ff.close()
    if population=={}:
        print("No existing cluster in",POPname, "is smaller than the limit size -c and lives longer than -b. Therefore, no output file was produced.")
        sys.exit()
    printfiles(population,POPname,Nsteps*TimeStep)


def printfiles(population,POPname,dt):
    for clustername in population:                                                                #loop over all the cluster types
        filename = POPname[:-9] + '__' + clustername + '.msd.dat'                                    #each one will go into a separate file
        print(filename)
        msds = population[clustername]
        nlags = max(len(msd) for msd in msds)
        ff = open(filename,'w')
        header = 'time_(fs)\t' + '\t'.join('cl'+str(jj) for jj in range(len(msds))) + '\tmean\n'  #one column per individual cluster, then their mean
        ff.write(header)
        for kk in range(nlags):
            values = [msd[kk] for msd in msds if kk < len(msd)]
            row = str(dt*kk) + '\t'
            row = row + '\t'.join(str(msd[kk]) if kk < len(msd) else '' for msd in msds)
            row = row + '\t' + str(sum(values)/len(values)) + '\n'
            ff.write(row)
        ff.close()


def main(argv):
    UMDname = ''
    POPname = ''
    Nsteps = 1
    ClusterMaxSize = 20
    MinLife = 20
    umdpf.headerumd()
    try:
        opts, arg = getopt.getopt(argv,"hf:p:s:b:c:",["fUMDfile","pPOPfile","sSampling_Frequency","bBallistic","cClusterMaxSize"])
    except getopt.GetoptError:
        print ('msd_clusters_umd.py -f <UMD_filename> -p <POPul_filename> -s <Sampling_Frequency> -b <MinLife> -c <ClusterMaxSize>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('Program to compute the mean square displacements of the centers of mass of all the atomic clusters found by speciation_and_angles.py. Usage:')
            print ('msd_clusters_umd.py -f <UMD_filename> -p <POPUL_filename> -s <Sampling_Frequency> -b <MinLife> -c <ClusterMaxSize>')
            print ('UMD_filename = input file with the trajectory, in UMD format.')
            print ('POPul_filename = input file with all the individual atomic clusters and their lifetimes in popul format.')
            print ('Sampling_Frequency = frequency of sampling of the trajectory. Default = 1 (all steps are considered)')
            print ('MinLife = minimum lifetime of a cluster, in steps of the UMD trajectory, for its MSD to be computed. Default = 20.')
            print ('ClusterMaxSize = arbitrary upper size limit for defining chemical species. Default = 20. You can safely assume values up to 10-12.' )
            print ('Output : one <popul base>__<cluster>.msd.dat file per type of cluster, with the time (fs), the MSD (A^2) of each individual cluster, and their mean.')
            sys.exit()
        elif opt in ("-f", "--fUMD_filename"):
            UMDname = str(arg)
            print ('UMDname = ',UMDname)
        elif opt in ("-p", "--pPOP_filename"):
            POPname = str(arg)
            print ('POPname = ',POPname)
        elif opt in ("-s", "--sNsteps"):
            Nsteps = int(arg)
            print ('Density of sampling: every ',Nsteps,' steps')
        elif opt in ("-b", "--bBallistic"):
            MinLife = int(arg)
            print ('Minimum lifetime: ',MinLife,' steps')
        elif opt in ("-c", "--cClusterMaxSize"):
            ClusterMaxSize = int(arg)
            print ('Clusters larger than ',ClusterMaxSize,' atoms are ignored')


#checking the existence of the files before doing any other operation
    if os.path.isfile(UMDname):
        print ('using ',UMDname,' umd file')
    else:
        print ('the UMD file ',UMDname,' does not exist')
        sys.exit()

    if os.path.isfile(POPname):
        print ('using ',POPname,' population file')
    else:
        print ('the population file ',POPname,' does not exist')
        sys.exit()

#read the umd file
    MyCrystal,SnapshotsValues,TimeStep,numsteps = umdpf.read_values(UMDname,"absxcart","line")
    AbsXcart = np.array(SnapshotsValues).reshape(len(SnapshotsValues),MyCrystal.natom,3)
    print()
    ana_popul(MyCrystal,AbsXcart,TimeStep,POPname,MinLife,ClusterMaxSize,Nsteps)



if __name__ == "__main__":
   main(sys.argv[1:])
