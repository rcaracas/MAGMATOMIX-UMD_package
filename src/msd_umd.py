#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### Author : Kevin Jiguet-Covex
### FFT engine and axes decomposition : Razvan Caracas



import sys,getopt,os.path
import umd_processes_fast as umdpf
import time
import ast
import numpy as np
from analyze_msd import calculation_diffusivities, block_diffusivity_error, diffusion_exponent, regime_warning, fit_end


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


def msd_atoms(X,lags,U=None,AxesAreCartesian=False,progress=True):#MSD of each atom, total and projected on the unit vectors U
    #X : (nsteps,natom,3) absolute (unwrapped) cartesian coordinates
    #Returns the total MSD (nlags,natom) and the projected MSDs (nlags,natom,nAxes), or None if there are no axes
    nsteps,natom,_ = X.shape
    total = np.zeros((len(lags),natom))
    proj = None if U is None else np.zeros((len(lags),natom,len(U)))
    batch = max(1, 5000000//(3*nsteps))#Atoms are processed in batches, to limit the memory used by the FFTs
    for a0 in range(0,natom,batch):
        if progress:
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
    if progress:
        sys.stdout.write("\rCalculating MSD. Progress : 100%\n")
    return total, proj


def msd_elements(total,proj,typat,ntypat):#Mean over the atoms of each element ; MSD[element][0] is the total MSD, MSD[element][1+k] its component along axis k
    nAxes = 0 if proj is None else proj.shape[2]
    MSD = np.zeros((ntypat,1+nAxes,total.shape[0]))
    for itypat in range(ntypat):
        MSD[itypat][0] = total[:,typat==itypat].mean(axis=1)
        for kk in range(nAxes):
            MSD[itypat][1+kk] = proj[:,typat==itypat,kk].mean(axis=1)
    return MSD


TensorAxes = np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,0,1]])/np.linalg.norm([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,0,1]],axis=1)[:,None]
TensorLabels = ['xx','yy','zz','xy','yz','zx']


def tensor_components(proj):#MSDs projected on x, y, z, (x+y)/sqrt2, (y+z)/sqrt2, (z+x)/sqrt2 (last axis) -> components <dx dx>, <dy dy>, <dz dz>, <dx dy>,
    #<dy dz>, <dz dx> of the MSD tensor : <du^2> = (<di^2> + <dj^2>)/2 + <di dj> along u = (i+j)/sqrt2
    T = proj.copy()
    for k,(i,j) in enumerate(((0,1),(1,2),(2,0))):
        T[...,3+k] = proj[...,3+k]-0.5*(proj[...,i]+proj[...,j])
    return T


def tensor(D6):#3x3 symmetric tensor from its components xx, yy, zz, xy, yz, zx
    return np.array([[D6[0],D6[3],D6[5]],[D6[3],D6[1],D6[4]],[D6[5],D6[4],D6[2]]])


def principal(D6):#Principal values (decreasing) and directions (columns) of a tensor given by its 6 components
    w, v = np.linalg.eigh(tensor(D6))
    order = np.argsort(w)[::-1]
    v = v[:,order]
    v *= np.sign(v[np.argmax(np.abs(v),axis=0),np.arange(3)])          #largest component of each direction positive
    return w[order], v


def plane_directions(arg,step=5):#Unit vectors every step degrees in the plane of two vectors (orthonormalized : true angles)
    P = np.array([float(x) for x in ast.literal_eval(arg)]).reshape(2,3)
    a = P[0]/np.linalg.norm(P[0])
    b = P[1]-(P[1] @ a)*a
    if np.linalg.norm(b) < 1e-8:
        print('ERROR : the two vectors of the plane -p are parallel.')
        sys.exit()
    b /= np.linalg.norm(b)
    thetas = np.arange(0,180,step)
    return thetas, np.cos(np.radians(thetas))[:,None]*a+np.sin(np.radians(thetas))[:,None]*b, a, b


def write_tensor(filename,umdfile,elements,Dfull,Dblocks,tmin,tmax,tend,plane):#Diffusion tensor of each element, principal values and directions, and D in a plane
    #Dfull[element] : the 6 components (m2/s) ; Dblocks[element] : the 6 components fitted on each block (list), or [] ; errors from the spread over the blocks,
    #scaled as block_diffusivity_error : stdev * sqrt(tmax/(nblocks*tend))
    def err(values):#values : one row per block
        return np.std(values,axis=0,ddof=1)*np.sqrt(tmax/(len(values)*tend))
    f = open(filename,'w')
    f.write('# Diffusion tensors from '+umdfile+' : <dr_i dr_j> = 2 D_ij t, fit between '+str(tmin)+' and '+str(tmax)+' fs ; D in m2/s ; '
            +('errors from '+str(len(Dblocks[0]))+' blocks' if len(Dblocks[0]) > 1 else 'no errors (no blocks)')+'\n')
    f.write('# D along a unit vector n : n.D.n ; isotropic D = trace/3 (= D of the total MSD) ; principal directions : cartesian unit vectors\n')
    for el,D6,B in zip(elements,Dfull,Dblocks):
        B = np.array(B)
        blocks = len(B) > 1
        e6 = err(B) if blocks else np.full(6,np.nan)
        w, v = principal(D6)
        ew = err([principal(b)[0] for b in B]) if blocks else np.full(3,np.nan)
        iso = np.trace(tensor(D6))/3.0
        eiso = err([np.trace(tensor(b))/3.0 for b in B]) if blocks else np.nan
        f.write('\nelement '+el+'\n')
        f.write('component\tD_(m2/s)\terror\n')
        for lab,d,e in zip(TensorLabels,D6,e6):
            f.write(lab+'\t'+str(d)+'\t'+str(e)+'\n')
        f.write('isotropic\t'+str(iso)+'\t'+str(eiso)+'\n')
        f.write('principal\tD_(m2/s)\terror\tx\ty\tz\n')
        for k in range(3):
            f.write('D'+str(k+1)+'\t'+str(w[k])+'\t'+str(ew[k])+'\t'+'\t'.join('%.6f' % c for c in v[:,k])+'\n')
        f.write('anisotropy D1/D3\t'+str(w[0]/w[2] if w[2] > 0 else float('nan'))+'\n')
        print('{:<4s} D (m2/s) isotropic {:.4e} +/- {:.1e} ; principal {:.4e} {:.4e} {:.4e} (+/- {:.1e} {:.1e} {:.1e})'.format(el,iso,eiso,*w,*ew))
        for k in range(3):
            print('       D'+str(k+1)+' along ['+', '.join('%.3f' % c for c in v[:,k])+']')
    if plane is not None:
        thetas, dirs, a, b = plane
        f.write('\nD along the directions cos(theta) a + sin(theta) b of the plane a = ['+', '.join('%.6f' % c for c in a)+'], b = ['+', '.join('%.6f' % c for c in b)+'] (orthonormalized)\n')
        f.write('theta_(deg)\t'+'\t'.join(el+'\t'+el+'_error' for el in elements)+'\n')
        for th,n in zip(thetas,dirs):
            row = [str(th)]
            for D6,B in zip(Dfull,Dblocks):
                row.append(str(n @ tensor(D6) @ n))
                row.append(str(err([n @ tensor(bb) @ n for bb in B])) if len(B) > 1 else 'nan')
            f.write('\t'.join(row)+'\n')
    f.close()
    print('Diffusion tensors printed in file ',filename)


def main(argv):
    vv = 1
    ballistic = 0
    FitEnd = None
    nBlocks = 5
    umdpf.headerumd()
    start=time.time()
    umdfile=''
    Mode = "elements"
    AxesMode = None
    Axes = None
    nCores=None
    planearg = None
    first = 0
    sampling = 1
    try:
        opts, arg = getopt.getopt(argv,"hf:i:s:z:v:b:e:n:m:a:k:p:",["fumdfile","iInitialStep","sSampling","zHorizontalJump","vVerticalJump","bBallistic","eFitEnd","nBlocks","mMode","aAxes","knCores","pPlane"])
    except getopt.GetoptError:
        print ('msd_umd.py -f <UMD_filename> -i <InitialStep> -s <Sampling> -v <VerticalJump> -b <Ballistic> -e <FitEnd> -n <nBlocks> -m <Mode> -a <Axes> -k <nCores>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('Program to compute the Mean Square Displacement. Usage: ')
            print ('msd_umd.py -f <UMD_filename> -i <InitialStep> -s <Sampling> -v <VerticalJump> -b <Ballistic> -e <FitEnd> -n <nBlocks> -m <Mode> -a <Axes> -k <nCores>')
            print ('UMD_filename = input file with the trajectory, in UMD format.')
            print ('-i : InitialStep. number of initial snapshots discarded. Default is 0.')
            print ('-s : Sampling. one snapshot out of s is used (the time step becomes s times longer). Default is 1. The steps of -v, -b, -e')
            print ('     count the snapshots used.')
            print (' The MSD is averaged over all the possible time origins, for time lags up to half the length of the trajectory.')
            print (' It is computed with Fast Fourier Transforms, so no sub-sampling of the time origins is needed.')
            print ('-v : VerticalJump. the MSD is written every VerticalJump time lags. Default is 1.')
            print ('-b : Ballistic. number of steps of the ballistic regime; the MSD is written and fitted only for time lags >= Ballistic. Default is 0. Typical values of 100 are sufficient.')
            print ('-e : FitEnd. last time lag, in steps, of the linear fit of the MSD giving the diffusion coefficients.')
            print ('     Default is 0.75 of the last lag of a block (see -n ; with -n 0 or 1, of a block as for 5 blocks), i.e. 7.5% of the trajectory for 5 blocks.')
            print ('     Longer fits, up to lags with few independent time origins, give much noisier D ; increase -e for slow, sub-diffusive systems.')
            print ('-n : nBlocks. the trajectory is also cut into nBlocks consecutive blocks; the spread of the D fitted on each block gives')
            print ('     the error on D, stdev(D_blocks)*sqrt(FitEnd/(nBlocks*end of the block fit)). The block MSDs are written in <base>.msd.blocks.dat, used also by analyze_msd.py.')
            print ('     The fit on the blocks ends at most at 0.75 of their last lag. Default is 5 ; 0 or 1 to switch off.')
            print ("-m : mode. 'elements' (default) to print the mean MSD of each element, or 'atoms' to print the MSD of each individual atom followed by the mean MSD of each element.")
            print ('-a : Axes. decomposition of the MSD along axes. If absent, only the total MSD is computed. Options are:')
            print ('     cart  : along the cartesian axes x, y, z. The three components add up to the total MSD.')
            print ('     cryst : along the crystallographic axes a, b, c of the initial cell (Auto is accepted as a synonym).')
            print ('             For non-orthogonal cells the three components do not add up to the total MSD.')
            print ('     a list of floats [x1,y1,z1,x2,y2,z2,...] describing any number of axes, e.g. -a [1,0,0,0,1,0,0,1,1]')
            print ('             calculates the MSD along the axes [1,0,0], [0,1,0] and [0,1,1].')
            print ('     Each axis component is the MSD of the displacement projected on the unit vector of that axis,')
            print ('     so that the diffusion coefficient along the axis is D = slope / 2 (D = slope / 6 for the total MSD).')
            print ('     tensor : the components xx, yy, zz, xy, yz, zx of the MSD tensor <dr_i dr_j> (from the MSDs along x, y, z and the three')
            print ('             diagonals (i+j)/sqrt2), fitted as 2 D_ij t, and <base>.diffusion.tensor.dat : the diffusion tensor of each element with')
            print ('             its errors, its principal values D1 >= D2 >= D3 and directions, the isotropic D (trace/3) and the anisotropy D1/D3.')
            print ('             D along any unit vector n is n.D.n. Noise makes D1 > D3 even for an isotropic system : compare D1 - D3 with the errors.')
            print ('-p : with -a tensor, D along the directions of a plane, every 5 degrees, given by two vectors [a1,a2,a3,b1,b2,b3] (made orthonormal),')
            print ('     added to <base>.diffusion.tensor.dat (formerly MSDdistrib.py).')
            print ('-k : nCores. number of cores used to read the UMD file. Default is all.')
            print ('Output : <base>.msd.dat, one time column (fs), then one column per element (or per atom and per element)')
            print (' named El for the total MSD, and El_x, El_a, El_1 ... for the components along the axes. MSD in A^2.')
            print (' <base>.diffusion.dat : diffusion coefficient of each element (and along each axis), D = slope/6 for the total MSD')
            print (' and slope/2 along one axis, with its error, the R^2 of the fit, and the exponent alpha of MSD ~ t^alpha over the fit window')
            print (' (1 for normal diffusion, < 1 if still sub-diffusive). For finer fits, plots (-v), or many files at once, use analyze_msd.py.')
            sys.exit()
        elif opt in ("-f", "--fumdfile"):
            umdfile = str(arg)
            print ('I will use the ',umdfile,' for input','\n')
        elif opt in ("-i", "--iInitialStep"):
            first = int(arg)
        elif opt in ("-s", "--sSampling"):
            sampling = int(arg)
            if sampling<=0 :
                print('ERROR : sampling -s has to be strictly positive.')
                sys.exit()
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
        elif opt in ("-e", "--eFitEnd"):
            FitEnd = int(arg)
        elif opt in ("-n", "--nBlocks"):
            nBlocks = int(arg)
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
            elif arg == "tensor":
                AxesMode = "tensor"
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
        elif opt in ("-p","--pPlane"):
            planearg = arg


    if (os.path.isfile(umdfile)):

        MyCrystal,SnapshotsValues,TimeStep,numsteps = umdpf.read_values(umdfile,"absxcart","line",sampling,first,nCores=nCores)
        TimeStep = TimeStep*sampling                                            #time between two snapshots used
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
        elif AxesMode == "tensor":
            U = TensorAxes
            AxesLabels = TensorLabels
        else :
            U = None
            AxesLabels = []
        if planearg is not None and AxesMode != "tensor":
            print('ERROR : the plane -p needs -a tensor.')
            sys.exit()
        plane = plane_directions(planearg) if planearg is not None else None
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

        #Blocks of the trajectory : their MSDs give the errors on D, and their length the default end of the fit
        #(fitting up to long lags, which have few independent time origins, makes D much noisier)
        blocklen = nsteps//(nBlocks if nBlocks > 1 else 5)
        blocklags = np.arange(ballistic,blocklen//2,vv)
        if FitEnd is None:
            if len(blocklags) == 0 :
                print('ERROR : the ballistic regime -b is longer than half a block of ',blocklen,' steps. Use fewer blocks -n, a smaller -b, or set the end of the fit -e.')
                sys.exit()
            FitEnd = int(0.75*blocklags[-1])
            print('The fit of the MSD ends at 0.75 of the last lag of a block of ',blocklen,' steps : ',FitEnd,' steps.')
        elif FitEnd > lags[-1]:
            print('The end of the fit -e is beyond the last lag ; it is set to ',lags[-1],' steps.')
            FitEnd = lags[-1]
        if np.sum((lags >= ballistic) & (lags <= FitEnd)) < 3 :
            print('ERROR : the fit window from -b ',ballistic,' to ',FitEnd,' steps must contain at least 3 time lags. Set a later end of the fit -e, or use fewer blocks -n.')
            sys.exit()

        total, proj = msd_atoms(X,lags,U,AxesMode=="cart")
        if AxesMode == "tensor":
            proj = tensor_components(proj)

        #MSD[species][0] is the total MSD, MSD[species][1+k] its component along axis k
        typat = np.array(MyCrystal.typat)
        MSDel = msd_elements(total,proj,typat,MyCrystal.ntypat)

        #MSD of each block of the trajectory, for the errors on the diffusion coefficients
        BlockMSD = []
        if nBlocks > 1 :
            if len(blocklags) == 0 or np.sum(blocklags <= 0.75*blocklags[-1]) < 3 :
                print('WARNING : the blocks are too short compared to the ballistic regime ; no errors on D. Use fewer blocks -n or a smaller -b.')
            else :
                print('Calculating the MSD of ',nBlocks,' blocks of ',blocklen,' steps')
                for kk in range(nBlocks):
                    tb, pb = msd_atoms(X[kk*blocklen:(kk+1)*blocklen],blocklags,U,AxesMode=="cart",progress=False)
                    if AxesMode == "tensor":
                        pb = tensor_components(pb)
                    BlockMSD.append(msd_elements(tb,pb,typat,MyCrystal.ntypat))
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
        ellabels = []
        for itypat in range(MyCrystal.ntypat):
            ellabels += [MyCrystal.elements[itypat]]+[MyCrystal.elements[itypat]+'_'+lab for lab in AxesLabels]
            columns += list(MSDel[itypat])
        labels += ellabels

        if umdfile.endswith('.umd.dat'):
            base = umdfile[:-8]
        else :
            base = os.path.splitext(umdfile)[0]
        msdfile = base + '.msd.dat'
        f = open(msdfile,'w')
        f.write('# MSD from '+umdfile+' : '+str(nsteps)+' snapshots'+(' from snapshot '+str(first) if first else '')+(' every '+str(sampling) if sampling > 1 else '')+', timestep '+str(TimeStep)+' fs, all time origins, lags from '+str(ballistic)+' to '+str(lags[-1])+' steps every '+str(vv)+' ; mode '+Mode+'\n')
        if U is not None:
            f.write('# axes (unit vectors) : '+' ; '.join(lab+' = '+np.array2string(u,precision=6,separator=',') for lab,u in zip(AxesLabels,U))+'\n')
        f.write('time_(fs)\t'+'\t'.join(labels)+'\n')
        for ii in range(len(lags)):
            f.write(str(Instants[ii])+'\t'+'\t'.join(str(col[ii]) for col in columns)+'\n')
        f.close()
        print ('MSDs for '+Mode+' printed in file ',msdfile)

        if BlockMSD != []:
            blockfile = base + '.msd.blocks.dat'
            f = open(blockfile,'w')
            f.write('# MSD of '+str(nBlocks)+' consecutive blocks of '+str(blocklen)+' steps of '+umdfile+', for the errors on the diffusion coefficients\n')
            f.write('time_(fs)\t'+'\t'.join(lab+'#'+str(kk+1) for lab in ellabels for kk in range(nBlocks))+'\n')
            for ii in range(len(blocklags)):
                f.write(str(blocklags[ii]*TimeStep)+'\t'+'\t'.join(str(BlockMSD[kk][itypat][jj][ii]) for itypat in range(MyCrystal.ntypat) for jj in range(1+nAxes) for kk in range(nBlocks))+'\n')
            f.close()
            print ('MSDs of the blocks printed in file ',blockfile)

        #Diffusion coefficients : linear fit of the MSD between the end of the ballistic regime and FitEnd
        tmin = ballistic*TimeStep
        tmax = FitEnd*TimeStep
        times = np.array(Instants)
        difffile = base + '.diffusion.dat'
        f = open(difffile,'w')
        if BlockMSD != []:
            errorstring = 'error from the spread of D over '+str(nBlocks)+' blocks of '+str(blocklen)+' steps'
        else :
            errorstring = 'error = standard error of the slope of the fit (underestimated : the MSD points are correlated)'
        f.write('# Diffusion coefficients from '+msdfile+' : linear fit of the MSD between '+str(tmin)+' and '+str(tmax)+' fs ; D = slope/6 (total) or slope/2 (along one axis) ; '+errorstring+'\n')
        f.write('column\tD_(m2/s)\tD_error_(m2/s)\tR_squared\talpha\n')
        print ('\nDiffusion coefficients, fit between '+str(tmin)+' and '+str(tmax)+' fs ; '+errorstring)
        print ('{:<10s}{:>14s}{:>14s}{:>10s}{:>8s}'.format('column','D (m2/s)','error','R^2','alpha'))
        icol = 0
        lowalpha = False
        Dtensor = []        #6 components of each element (tensor mode), on the full trajectory and on each block
        Dtensorblocks = []
        blocktimes = blocklags*TimeStep
        tend = min(fit_end(blocktimes),tmax) if BlockMSD != [] else tmax
        for itypat in range(MyCrystal.ntypat):
            for jj in range(1+nAxes):
                lab = ellabels[icol]
                icol += 1
                dimension = 3 if jj == 0 else 1
                D, D_err, R2, slope, intercept = calculation_diffusivities(times,MSDel[itypat][jj],tmin,dimension,tmax)
                if BlockMSD != []:
                    D_err = block_diffusivity_error(blocktimes,[BlockMSD[kk][itypat][jj] for kk in range(nBlocks)],tmin,tmax,dimension)
                offdiagonal = AxesMode == "tensor" and jj > 3          #<dr_i dr_j> with i != j : of any sign, no exponent
                alpha = float('nan') if offdiagonal else diffusion_exponent(times,MSDel[itypat][jj],tmin,tmax)#1 for normal diffusion, < 1 if the window is still sub-diffusive
                f.write(lab+'\t'+str(D)+'\t'+str(D_err)+'\t'+str(R2)+'\t'+str(alpha)+'\n')
                print ('{:<10s}{:>14.4e}{:>14.4e}{:>10.4f}{:>8.2f}'.format(lab,D,D_err,R2,alpha))
                if not offdiagonal and regime_warning(alpha) != '':
                    lowalpha = True
            if AxesMode == "tensor":
                Dtensor.append([calculation_diffusivities(times,MSDel[itypat][1+k],tmin,1,tmax)[0] for k in range(6)])
                Dtensorblocks.append([[calculation_diffusivities(blocktimes,BlockMSD[kk][itypat][1+k],tmin,1,tend)[0] for k in range(6)] for kk in range(nBlocks)] if BlockMSD != [] else [])
        f.close()
        print ('Diffusion coefficients printed in file ',difffile)
        if AxesMode == "tensor":
            print('')
            write_tensor(base+'.diffusion.tensor.dat',umdfile,MyCrystal.elements,Dtensor,Dtensorblocks,tmin,tmax,tend,plane)
        if lowalpha:
            print ('WARNING : alpha is outside 0.9-1.1 for some columns : the fit window is not in the diffusive regime (MSD ~ t^alpha) ;')
            print ('          alpha < 0.9 (sub-diffusive) : consider a larger -e ; alpha > 1.1 (ballistic) : consider a larger -b and -e.')

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
