###
##AUTHORS: RAZVAN CARACAS
###

import sys,getopt,os.path,re
from . import crystallography as cr
from . import umd_process as umdp


def msd_trajectory(trajectory):
    #print('trajectory size ',len(trajectory))
    partialmsd = [0.0 for _ in range(int(len(trajectory)/2))]
    for ii in range(int(len(trajectory)/2)):
        ref = [0.0 for _ in range(3)]
        for kk in range(3):
            ref = trajectory[ii]
        #print('reference center of mass ',ref)
        for jj in range(int(len(trajectory)/2)):
            partialmsd[jj] = partialmsd[jj] + (trajectory[ii+jj][0]-ref[0])**2 + (trajectory[ii+jj][1]-ref[1])**2 + (trajectory[ii+jj][2]-ref[2])**2
    for ii in range(int(len(trajectory)/2)):
        partialmsd[ii] = partialmsd[ii] / (len(trajectory)/2)
    return(partialmsd)
                        
                        
def ana_popul(MyCrystal,AllSnapshots,POPname,Ballistic,ClusterMaxSize,Nsteps):
    print ('in analysis of the population')
    print ('reading cluster population from ',POPname,' file')
    print ('initial total no of snapshots:',len(AllSnapshots))
    population = {}
    ff = open(POPname,'r')
    ff.readline()
    ff.readline()
    clusterindex = 0
    maxtime = 0
    while True:
        line = ff.readline()
        if not line: break
        line = re.sub('[\[,\]]', '',line)
        line = line.strip()
        entry = line.split()
        #print
        #print ('treating current line ',entry)
        sizeofentry = len(entry)
        if sizeofentry -4 <= ClusterMaxSize :
            #print ('reading line with ', entry[0],' living for ',int(entry[3]),' steps')
            if int(entry[3]) >= Ballistic :
                partialmsd = []
                newmsd = []
                clustername = entry[0]
                timediff = int(entry[3])
                enddiff = int(entry[2])
                begindiff = int(entry[1])
                if timediff > maxtime:
                    maxtime = timediff
                    #print('current maximum time is ',maxtime)
                #print('Cluster ok. Analyzing cluster with name ',clustername, ' living for ',timediff,' total steps')
                trajectory = [ [0.0,0.0,0.0] for _ in range(int(timediff/Nsteps)) ]                 #stores the coordinates of the center of mass for current cluster
                #print ('trajectory will have ',len(trajectory),' instances')
                for itime in range(int(timediff/Nsteps)):
                    #print('current trajectory instance is ',itime,' at actual time ',begindiff + itime*Nsteps)
                    for iatom in range(4,sizeofentry):
                        for kk in range(3):                                                         #computes the coordinates of the center of mass for the current cluster at each timestep
                            trajectory[itime][kk] = trajectory[itime][kk] + AllSnapshots[begindiff + itime*Nsteps].atoms[int(entry[iatom])].xcart[kk]
                #print('trajectory array is',trajectory)
                partialmsd = msd_trajectory(trajectory)                                             #calls the MSD on the trajectory
                #print('corresponding msd is ',partialmsd)
                if len(partialmsd)> 1:                                                              #there must be at least two steps for diffusion to be added
                    flagalive = 0                                                                   #flag to check for previous existence of the cluster
                    for ll in population.keys():                                                    #runs over the ewntire dictionary looking for the current cluster
                        #print ('comparing current cluster ',clustername,' with existent clusters:',population[ll]['clustername']) #,' with ',population[kk][clustername])
                        if clustername == population[ll]['clustername']:                            #cluaster already exists
#                            newmsd = population[ll]['diffusion'] + [partialmsd]                     #appends the msd of the current cluster to the bigger family
                            population[ll]['diffusion'].append(partialmsd)
                            #print('new msd array, after last assignement',newmsd)
                            #population[ll] = {'clustername':clustername,'diffusion':newmsd}         #updates the dictionary
                            flagalive = 1                                                           #resets the flag
                    if flagalive == 0:                                                              #cluster is new
                        population[clusterindex] = {'clustername':clustername,'diffusion':[partialmsd]}       #adds the new cluster and its msd to the dictionary
                        clusterindex = clusterindex + 1                                             #index over the totla number of entries, i.e. clusters, from the dictionary
    ff.close()
    #print(population)
    printfiles(population,POPname,Nsteps,maxtime)


def printfiles(population,POPname,Nsteps,maxtime):
    for ii in range(len(population)):                                                               #loop over all the cluster types
        filename = POPname[:-9] + '__' + population[ii]['clustername'] + '.msd.dat'                 #each one will go into a seprate file
        print(filename)
        ff = open(filename,'w')
        #print('clusters size ',len(population[ii]['diffusion']))
        header = 'time(fs)'
        for jj in range(len(population[ii]['diffusion'])+1):                                           #loop over all the individsual clusters
            header = header + 'cl' + str(jj) + '\t'
        header = header + '\n'
        ff.write(header)
        msd = []
        row = ''
        currmsd = [['' for _ in range(len(population[ii]['diffusion'])+2)] for _ in range(int(maxtime/Nsteps))]
        tlmsd = [[0.0,0] for _ in range(int(maxtime/Nsteps))]
        for jj in range(len(population[ii]['diffusion'])):
            msd = population[ii]['diffusion'][jj]
            for kk in range(len(msd)):
                currmsd[kk][jj] = str(msd[kk])
                tlmsd[kk][0] = tlmsd[kk][0] + msd[kk]
                tlmsd[kk][1] += 1
        for kk in range(int(maxtime/Nsteps)):
            row = str(Nsteps * kk) + '\t'
            for jj in range(len(population[ii]['diffusion'])):
                row = row + str(currmsd[kk][jj]) + '\t'
            if tlmsd[kk][1] > 0 :
                row = row + str(tlmsd[kk][0] / float(tlmsd[kk][1])) + '\t'
            else:
                row = row + '\t'
            row = row + '\n'
            ff.write(row)
    ff.close()




                                 
def main(argv):
    UMDname = ''
    POPname = ''
    RADname = ''
    Nsteps = 1
    ClusterMaxSize = 20
    AtomicRadius = 1.0
    Ballistic = 100
    umdp.headerumd()
    try:
        opts, arg = getopt.getopt(argv,"hf:p:s:b:c:",["fUMDfile","pPOPfile","sSampling_Frequency","bBallistic","cClusterMaxSize"])
    except getopt.GetoptError:
        print ('msd_cluster_umd.py -f <UMD_filename> -p <POPul_filename> -s <Sampling_Frequency> -b <Ballistic> -c <ClusterMaxSize>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('Program to compute the mean square displacements for all the atomic species found by speciation.py. Usage:')
            print ('msd_cluster_umd.py -f <UMD_filename> -p <POPUL_filename> -s <Sampling_Frequency> -b <Ballistic> -c <ClusterMaxSize>')
            print ('UMD_filename = input file with the trajectory, in UMD format.')
            print ('POPul_filename = input file with all the individual atomic clusters and their lifetimes in popul format.')
            print ('Sampling_Frequency = frequency of sampling of the trajectory. Default = 1 (all steps are considered)')
            print ('Ballistic = estimation of the ballistic part of the trajectory. Default is 0. Typical values of 100 are sufficient.')
            print ('ClusterMaxSize = arbitrary upper size limit for defining checmnial species. Default = 5. You can safely assume values up to 10-12.' )
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
        elif opt in ("-s", "--bBallistic"):
            Ballistic = int(arg)
            print ('Ballistic threshold: ',Ballistic,' steps')
        elif opt in ("-c", "--cClusterMaxSize"):
            ClusterMaxSize = int(arg)
            print ('Clusters smaller than ',ClusterMaxSize,' atoms are considered gas')
                
                                 
#checking the existence of the files before doing any other operation
    if os.path.isfile(UMDname):
        print ('using ',UMDname,' umd file')
    else:
        print ('the UMD file ',UMDname,' does not exist')
        sys.exit()
    
    if os.path.isfile(POPname):
        print ('using ',POPname,' population file')
    else:
        print ('the UMD file ',POPname,' does not exist')
        sys.exit()

#read the umd file and allocate structure
    MyCrystal = cr.Lattice()
    AllSnapshots = [cr.Lattice]
    (MyCrystal,AllSnapshots,TimeStep)=umdp.read_absxcart(UMDname)

    ana_popul(MyCrystal,AllSnapshots,POPname,Ballistic,ClusterMaxSize,Nsteps)

 

if __name__ == "__main__":
   main(sys.argv[1:])



