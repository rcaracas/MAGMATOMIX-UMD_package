#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### Author : Kevin Jiguet-Covex
### FFT engine and axes decomposition : Razvan Caracas



import sys,getopt,os.path
import umd_processes_fast as umdpf
import time
import ast
import numpy as np


def msd_fft(x):#Mean square displacement of each column of x (time along axis 0), averaged over all the time origins
    #Uses the FFT algorithm of Calandrini et al. (nMoldyn) : MSD(m) = S1(m) - 2 S2(m), in O(N log N) instead of O(N^2)
    #Returns MSD(m) for all the lags m = 0 .. len(x)-1
    n = x.shape[0]
    x = x - x.mean(axis=0)#The MSD is invariant to a constant shift ; centring limits the round-off errors on long unwrapped trajectories
    F = np.fft.rfft(x, n=2*n, axis=0)#Zero-padding to 2n gives the non-cyclic autocorrelation
    S2 = np.fft.irfft(F*F.conj(), n=2*n, axis=0)[:n]#S2(m) = sum_k x(k).x(k+m)
    D = x**2
    cs = np.concatenate([np.zeros((1,)+D.shape[1:]), np.cumsum(D, axis=0)])
    lags = np.arange(n)
    S1 = (cs[n] - cs[lags]) + cs[n-lags]#S1(m) = sum_{k=0}^{n-1-m} x(k+m)^2 + x(k)^2
    counts = (n - lags).reshape((-1,)+(1,)*(x.ndim-1))#Number of time origins for each lag
    msd = (S1 - 2.0*S2)/counts
    msd[0] = 0.0#Exactly zero by definition, instead of a round-off residual
    return msd


def msd_atoms(X,lags,U=None,AxesAreCartesian=False):#MSD of each atom, total and projected on the unit vectors U
    #X : (nsteps,natom,3) absolute (unwrapped) cartesian coordinates
    #Returns the total MSD (nlags,natom) and the projected MSDs (nlags,natom,nAxes), or None if there are no axes
    nsteps,natom,_ = X.shape
    total = np.zeros((len(lags),natom))
    proj = None if U is None else np.zeros((len(lags),natom,len(U)))
    batch = max(1, 5000000//(3*nsteps))#Atoms are processed in batches, to limit the memory used by the FFTs
    for a0 in range(0,natom,batch):
        sys.stdout.write("\rCalculating MSD. Progress : "+str((100*a0)//natom)+"%")
        sys.stdout.flush()
        Xb = X[:,a0:a0+batch,:]
        cart = msd_fft(Xb)[lags]#MSD along x, y and z ; their sum is the total MSD, for any cell shape
        total[:,a0:a0+batch] = cart.sum(axis=-1)
        if U is not None:
            if AxesAreCartesian:
                proj[:,a0:a0+batch,:] = cart
            else:
                proj[:,a0:a0+batch,:] = msd_fft(Xb @ U.T)[lags]#Projection of the trajectory on each axis, then 1D MSD
    sys.stdout.write("\rCalculating MSD. Progress : 100%\n")
    return total, proj


def main(argv):
    vv = 1
    ballistic = 0
    umdpf.headerumd()
    start=time.time()
    umdfile=''
    Mode = "elements"
    AxesMode = None
    Axes = None
    nCores=None
    try:
        opts, arg = getopt.getopt(argv,"hf:z:v:b:m:a:k:",["fumdfile","zHorizontalJump","vVerticalJump","bBallistic","mMode","aAxes","knCores"])
    except getopt.GetoptError:
        print ('msd_umd.py -f <UMD_filename> -v <VerticalJump> -b <Ballistic> -m <Mode> -a <Axes> -k <nCores>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('Program to compute the Mean Square Displacement. Usage: ')
            print ('msd_umd.py -f <UMD_filename> -v <VerticalJump> -b <Ballistic> -m <Mode> -a <Axes> -k <nCores>')
            print ('UMD_filename = input file with the trajectory, in UMD format.')
            print (' The MSD is averaged over all the possible time origins, for time lags up to half the length of the trajectory.')
            print (' It is computed with Fast Fourier Transforms, so no sub-sampling of the time origins is needed.')
            print ('-v : VerticalJump. the MSD is written every VerticalJump time lags. Default is 1.')
            print ('-b : Ballistic. number of steps of the ballistic regime; the MSD is written only for time lags >= Ballistic. Default is 0. Typical values of 100 are sufficient.')
            print ("-m : mode. 'elements' (default) to print the mean MSD of each element, or 'atoms' to print the MSD of each individual atom followed by the mean MSD of each element.")
            print ('-a : Axes. decomposition of the MSD along axes. If absent, only the total MSD is computed. Options are:')
            print ('     cart  : along the cartesian axes x, y, z. The three components add up to the total MSD.')
            print ('     cryst : along the crystallographic axes a, b, c of the initial cell (Auto is accepted as a synonym).')
            print ('             For non-orthogonal cells the three components do not add up to the total MSD.')
            print ('     a list of floats [x1,y1,z1,x2,y2,z2,...] describing any number of axes, e.g. -a [1,0,0,0,1,0,0,1,1]')
            print ('             calculates the MSD along the axes [1,0,0], [0,1,0] and [0,1,1].')
            print ('     Each axis component is the MSD of the displacement projected on the unit vector of that axis,')
            print ('     so that the diffusion coefficient along the axis is D = slope / 2 (D = slope / 6 for the total MSD).')
            print ('-k : nCores. number of cores used to read the UMD file. Default is all.')
            print ('Output : <base>.msd.dat, one time column (fs), then one column per element (or per atom and per element)')
            print (' named El for the total MSD, and El_x, El_a, El_1 ... for the components along the axes. MSD in A^2.')
            sys.exit()
        elif opt in ("-f", "--fumdfile"):
            umdfile = str(arg)
            print ('I will use the ',umdfile,' for input','\n')
        elif opt in ("-z", "--zHorizontalJump"):
            print ('The horizontal jump -z is no longer needed (all the time origins are used) and is ignored.')
        elif opt in ("-v", "--vVerticalJump"):
            vv = int(arg)
            if vv<=0 :
                print('ERROR : vertical jump -v has to be strictly positive.')
                sys.exit()
        elif opt in ("-b", "--bBallistic"):
            ballistic = int(arg)
            if ballistic<0 :
                print('ERROR : ballistic -b has to be positive or zero.')
                sys.exit()
        elif opt in ("-m","--mMode"):
            Mode=str(arg)
            if Mode not in ("elements","atoms"):
                print("ERROR : mode -m has to be 'elements' or 'atoms'.")
                sys.exit()
        elif opt in ("-a","--aAxes"):
            if arg in ("cart","xyz"):
                AxesMode = "cart"
            elif arg in ("cryst","abc","Auto"):
                AxesMode = "cryst"
            else :
                AxesMode = "custom"
                try :
                    Axes = [float(x) for x in ast.literal_eval(arg)]
                except (ValueError,SyntaxError,TypeError):
                    print('ERROR : axes -a have to be cart, cryst, or a list of floats such as [1,0,0,0,1,1].')
                    sys.exit()
                if len(Axes)==0 or len(Axes)%3!=0 :
                    print('ERROR : the list of axes -a must contain 3 components per axis.')
                    sys.exit()
        elif opt in ("-k","--nCores"):
            nCores = int(arg)


    if (os.path.isfile(umdfile)):

        MyCrystal,SnapshotsValues,TimeStep,numsteps = umdpf.read_values(umdfile,"absxcart","line",nCores=nCores)
        X = np.array(SnapshotsValues).reshape(len(SnapshotsValues),MyCrystal.natom,3)#(step, atom, xyz)
        nsteps = len(X)
        print ('\nNumber of atoms of each type is ',MyCrystal.types)

        #Projection axes, as unit vectors (rows of U)
        if AxesMode == "cart":
            U = np.eye(3)
            AxesLabels = ['x','y','z']
        elif AxesMode == "cryst":
            U = np.array(MyCrystal.rprimd,dtype=float)#Lattice vectors a, b, c of the first snapshot
            AxesLabels = ['a','b','c']
        elif AxesMode == "custom":
            U = np.array(Axes).reshape(-1,3)
            AxesLabels = [str(i+1) for i in range(len(U))]
        else :
            U = None
            AxesLabels = []
        if U is not None:
            norms = np.linalg.norm(U,axis=1)
            if np.any(norms==0.0):
                print('ERROR : the axes cannot be null vectors.')
                sys.exit()
            U = U/norms[:,None]
            Axes = U.flatten().tolist()
        nAxes = len(AxesLabels)

        maxlag = nsteps//2
        if ballistic >= maxlag :
            print('ERROR : the ballistic regime -b (',ballistic,' steps) must be shorter than half the trajectory (',maxlag,' steps).')
            sys.exit()
        lags = np.arange(ballistic,maxlag,vv)
        Instants = list(lags*TimeStep)

        total, proj = msd_atoms(X,lags,U,AxesMode=="cart")

        #MSD[species][0] is the total MSD, MSD[species][1+k] its component along axis k
        typat = np.array(MyCrystal.typat)
        MSDel = np.zeros((MyCrystal.ntypat,1+nAxes,len(lags)))
        for itypat in range(MyCrystal.ntypat):
            MSDel[itypat][0] = total[:,typat==itypat].mean(axis=1)
            for kk in range(nAxes):
                MSDel[itypat][1+kk] = proj[:,typat==itypat,kk].mean(axis=1)
        if Mode == "atoms":
            MSDat = np.zeros((MyCrystal.natom,1+nAxes,len(lags)))
            for iatom in range(MyCrystal.natom):
                MSDat[iatom][0] = total[:,iatom]
                for kk in range(nAxes):
                    MSDat[iatom][1+kk] = proj[:,iatom,kk]

        #Column labels and data
        labels = []
        columns = []
        if Mode == "atoms":
            indexat = [0 for _ in range(MyCrystal.ntypat)]
            for iatom in range(MyCrystal.natom):
                name = MyCrystal.elements[MyCrystal.typat[iatom]]+str(indexat[MyCrystal.typat[iatom]])#Atoms are numbered within each element
                indexat[MyCrystal.typat[iatom]] += 1
                labels += [name]+[name+'_'+lab for lab in AxesLabels]
                columns += list(MSDat[iatom])
        for itypat in range(MyCrystal.ntypat):
            labels += [MyCrystal.elements[itypat]]+[MyCrystal.elements[itypat]+'_'+lab for lab in AxesLabels]
            columns += list(MSDel[itypat])

        if umdfile.endswith('.umd.dat'):
            msdfile = umdfile[:-8] + '.msd.dat'
        else :
            msdfile = os.path.splitext(umdfile)[0] + '.msd.dat'
        f = open(msdfile,'w')
        f.write('# MSD from '+umdfile+' : '+str(nsteps)+' snapshots, timestep '+str(TimeStep)+' fs, all time origins, lags from '+str(ballistic)+' to '+str(lags[-1])+' steps every '+str(vv)+' ; mode '+Mode+'\n')
        if U is not None:
            f.write('# axes (unit vectors) : '+' ; '.join(lab+' = '+np.array2string(u,precision=6,separator=',') for lab,u in zip(AxesLabels,U))+'\n')
        f.write('time_(fs)\t'+'\t'.join(labels)+'\n')
        for ii in range(len(lags)):
            f.write(str(Instants[ii])+'\t'+'\t'.join(str(col[ii]) for col in columns)+'\n')
        f.close()
        print ('MSDs for '+Mode+' printed in file ',msdfile)

        end=time.time()
        print("runtime :",end-start)

        if Mode == "atoms":
            return [msdfile,MSDat,Instants,MyCrystal.elements,Axes]
        return [msdfile,MSDel,Instants,MyCrystal.elements,Axes]

    else:
        print ('umd file ',umdfile,'does not exist')
        sys.exit()


if __name__ == "__main__":
    main(sys.argv[1:])
