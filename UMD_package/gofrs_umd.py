#!/usr/bin/env python3
###
#AUTHORS: RAZVAN CARACAS, ANAIS KOBSCH
#MODIFICATIONS: LAURENT GILQUIN, 5 September 2020
###

import os
import itertools
import getopt
import sys
import numpy as np
from distutils.util import strtobool
from . import crystallography as cr
from . import umd_process as umdp
from . import c_gofr
from . import gpu_utils


def umd_pdist(X, coeff):
    """
    Interface with C wrapper to compute the umd distance, euclidean distance
    with shift.
    """
    # correct X type if not double
    if (X.dtype != np.double):
        X = np.ascontiguousarray(X, dtype=np.double)
    # initialize output
    res = np.zeros(int(X.shape[0] * (X.shape[0] - 1)/2), dtype=np.double,
                   order='C')
    # call the CPython wrapper
    c_gofr.umd_pdist_wrapper(X, res, coeff)

    return res


    
def print_gofrs(umdfile, MyCrystal, discrete, normalization, maxlength, gofr):
    """
    Normalize gofr, compute integral, create the new file gofrs.dat and write
    data in it."""

    # open output file
    filename = umdfile[:-8] + '.gofr.dat'
    nf = open(filename, 'w')
    # write header
    pairs = itertools.product(MyCrystal.elements, MyCrystal.elements)
    string = ['dist ']
    string.append("\t".join([elem[0]+'-'+elem[1]+'\t Int('+elem[0]+'-'+elem[1]+')' for elem in pairs]))
    string.append('\n')
    nf.write("".join(string))
    nf.close()

    # write data to file
    types = np.array(MyCrystal.types, dtype=int, order='C')
    discrete = np.double(discrete)
    maxlength = np.double(maxlength)
    normalization = np.double(normalization)
    c_gofr.print_gofr_wrapper(gofr, types, MyCrystal.ntypat, maxlength, discrete,
                              normalization, filename)
    return None # empty returns



def computeallgofrs(MyCrystal, atoms_coords, discrete, maxlength, gofr, gpu):
    """
    For every pair of atoms, compute the minimal distance betweens atoms
    and add 1 to the corresponding gofr[r, r + delta_r] slot.
    """
    # check type arguments
    coeffs = np.array(MyCrystal.acell, dtype=np.double, order='C')
    types = np.array(MyCrystal.typat, dtype=np.int, order='C')
    discrete = np.double(discrete)
    maxlength = np.double(maxlength)
    # check if calculation is to be done on gpu
    if gpu is not None:
        # compute the distance matrix through gpu
        atoms_dist = gpu.compute_distance(atoms_coords, coeffs)
        # compute the gofr through C extension
        c_gofr.compute_gofrM_wrapper(atoms_dist.astype(np.double), gofr, types,
                                     maxlength, discrete, MyCrystal.ntypat)
    # or else, compute directly gofr through C extension
    else:
        c_gofr.compute_gofr_wrapper(atoms_coords, gofr, coeffs, types,
                                    maxlength, discrete, MyCrystal.ntypat)
    return None # empty return



def read_umd(umdfile, Nsteps, discrete, InitialStep, **kwargs):
    """
    Read umd file and store data into classes.
    kwargs is meant to pass build flags to the gpu device.
    """

    # convert str to int if possible
    def str_to_int(x):
        
        try:
            return int(x)
        except ValueError:
            return x

    # initialize class storing the data
    MyCrystal = cr.Lattice()
    # initialize header flag and the set of attributes to be reads
    in_head = True
    head_keys = set(["elements", "natom", "ntypat", "typat", "types"])
    
    # open umd file
    with open(umdfile, 'r') as ff:

        # get the head attributes
        while in_head:
            line = ff.readline()
            if not line:
                break
            # pass empty line
            if len(line) > 1:
                # format line
                entry = (line.strip()).split()
                # check if first element is an attribute
                key = entry[0]
                if key in head_keys:
                    # remove the key
                    head_keys.remove(key)
                    values = [str_to_int(val) for val in entry[1:]]
                    if len(values) == 1:
                        setattr(MyCrystal, key, values[0])
                    else:
                        setattr(MyCrystal, key, values)
                    # check if set is empty
                    in_head = len(head_keys) != 0
        
        # read the first acell line
        while key != "acell":
            line = ff.readline()
            if len(line) > 1:
                key = (line.strip()).split()[0]
        entry = (line.strip()).split()
        values = [float(val) for val in entry[1:-1]]
        MyCrystal.acell = values
        
        # initialize gofr structure
        maxlength = min(MyCrystal.acell[:2])
        ndivx = int(maxlength / (2 * discrete))
        print('acell', MyCrystal.acell[0], 'discrete', discrete, 
              'ndivx ', ndivx)
        kvals = int(MyCrystal.ntypat * (MyCrystal.ntypat + 1) / 2)
        gofr = np.zeros((kvals, ndivx+1), dtype=np.intp, order='C')

        # check if gpu must be used and initialize the class
        use_gpu = kwargs["use_gpu"]
        del kwargs["use_gpu"]
        if use_gpu:
            gpu = gpu_utils.gpu((MyCrystal.natom, 3), custom_flags = kwargs)
            print("The following build flags will be used:\n{}".format(gpu.build_flags))
            if gpu.device is None: # catch if gpu memory is not large enough
                gpu = None
        else:
            gpu = None
        

        # get the atoms coordinates and calculate gofr for each time step
        niter = 0
        istep = 0
        while not in_head:
            
            line = ff.readline()
            if not line:
                break
            # pass empty line
            if len(line) > 1:
                entry = (line.strip()).split()       
                if entry[0] == 'atoms:':
                    istep += 1
                    # start computing from InitialStep
                    if istep > InitialStep:
                        # compute gofr only every Nsteps
                        if niter % Nsteps == 0:
                            # instantiate array to store atoms coordinates
                            atoms_coords = np.empty((MyCrystal.natom, 3), order='C',
                                                    dtype=np.double)
                            # read coordinates and fill the structure
                            for iatom in range(MyCrystal.natom):
                                line = ff.readline()
                                entry = (line.strip()).split()
                                for jj in range(3):
                                    atoms_coords[iatom][jj] = np.double(entry[jj+3])
                            # compute the gofr every Nsteps
                            computeallgofrs(MyCrystal, atoms_coords, discrete, maxlength,
                                            gofr, gpu)
                        else:
                            # skip exactly natom lines to speed up the process
                            for _ in range(MyCrystal.natom):
                                next(ff)
                        niter += 1
                    else:
                        # skip exactly natom lines to speed up the process
                          for _ in range(MyCrystal.natom):
                              next(ff)

    # if required, free memory held by gpu
    if gpu is not None:
        gpu.free_memory_pools()

    # number of time steps actually used in the calculation of gofr
    normalization = niter / Nsteps
    print_gofrs(umdfile, MyCrystal, discrete, normalization, maxlength, gofr)

    # empty return
    return None



def main(argv):
    umdp.headerumd()
    umdfile = 'output.umd.dat'
    Nsteps = 1
    discrete = 0.01            #delta_r  = width of bins in histogram
    InitialStep = 0           #in case we want to skip additional timesteps
    gpu_options = {"use_gpu": True}
    try:
        opts, arg = getopt.getopt(argv,"hf:s:d:i:g:", ["fumdfile", "Sampling_Frequency", "ddiscrete", "iInitialStep",
                                                       "ggpu_options"])
    except getopt.GetoptError:
        print ('gofrs_umd.py -f <umdfile>  -s <No.steps>  -d <discretization.interval> -i <InitialStep> -g <gpu_options>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('gofrs_umd.py program to compute pair-distribution fuctions g(r) and print all the results in one file. Usage:')
            print ('gofrs_umd.py -f <umdfile> -s <Sampling_Frequency> -d <discretization.interval> -i <InitialStep> -b <UseGPU>')
            print ('umdfile = name of the umd file. Default = output.umd.dat')
            print ('Sampling_Frequency = frequency of sampling the trajectory. Default = 1')
            print ('discretization.interval = for plotting the g(r). Default = 0.01 Angstroms')
            print ('InitialStep = number of initial steps of the trajectory to be removed. Default = 0')
            print ("""
                   gpu_options = a comma separated list of the following options for the distance calculation by gpu:\n
                       * use_gpu: a boolean specifying whether to use or not a GPU for the calculation.\n
                       * FINT: a string in {"ushort", "int"} specifying the type of integer to use.\n
                       * FREAL: a string in {"double", "float"} specifying the type of floating precision to use.\n
                       * FHALF: a boolean indicating if the half floating point precision should be use
                       for the inputs and ouput. This option offers a compromise between speed and
                       precision.\n
                    By default, only the user_gpu option is used and set to True, the remaining options
                    are automatically deduced.\n
                    """)
            sys.exit()
        elif opt in ("-f", "--fumdfile"):
            umdfile = str(arg)
            #print('umdfile = ',umdfile)
        elif opt in ("-s", "--sNsteps"):
            Nsteps = int(arg)
        elif opt in ("-d", "--ddiscrete"):
            discrete = float(arg)
            #print('the distance delta_r we use to compute the number of atoms around the central one is ',discrete, 'Angstroms')
        elif opt in ("-i", "--iInitialStep"):
            InitialStep = int(arg)
            #print('I will skip  ',initial,' timesteps. ')
        elif opt in ("-g", "--gpu_options"):
            for string in str(arg).split(","):
                key, val = string.split("=")
                if key in ["use_gpu", "FHALF"]:
                    try:
                        gpu_options[key] = bool(strtobool(val))
                    except ValueError:
                        print("In gpu_options, invalid argument {}.".format(key))
                        sys.exit()
                if key == "FINT":
                    if val in ["ushort", "int"]:
                        gpu_options[key] = val
                    else:
                        raise ValueError("In gpu_options, invalid argument 'FINT'.")
                        sys.exit()
                if key == "FREAL":
                    if val in ["double", "float"]:
                        gpu_options[key] = val
                    else:
                        raise ValueError("In gpu_options, invalid argument 'FREAL'.")
                        sys.exit()

    if (os.path.isfile(umdfile)): 
        read_umd(umdfile, Nsteps, discrete, InitialStep, **gpu_options)
    else:
        print ('the umd file ',umdfile,' does not exist')
        sys.exit()


if __name__ == "__main__":
   main(sys.argv[1:])
