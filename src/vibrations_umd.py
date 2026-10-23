#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### Authors : Zhi Li, Razvan Caracas, Kevin Jiguet-Covex
### FFT engine, axes decomposition, taper and smoothing, diffusion with block errors : Razvan Caracas
### Velocity autocorrelation functions, vibrational densities of states and diffusion coefficients, by element and along axes


import sys, getopt, os.path
import time
import ast
import numpy as np
from scipy.fft import dct
from scipy.integrate import cumulative_trapezoid
from scipy.ndimage import gaussian_filter1d
import umd_processes_fast as umdpf


CM1 = 33356.40952#1/fs in cm^-1
KB = 8.617333262E-05#Boltzmann constant in eV/K
AMU_A2FS2 = 1.66053906660E-27*1.0E10/1.602176634E-19#1 amu.(A/fs)^2 in eV
A2FS = 1.0E-05#1 A^2/fs in m^2/s


def acf_fft(x,nlags):#Autocorrelation of each column of x (time along axis 0), averaged over all the time origins, for lags 0 .. nlags-1
    n = x.shape[0]
    F = np.fft.rfft(x,n=2*n,axis=0)#Zero-padding to 2n gives the non-cyclic correlation
    acf = np.fft.irfft(F*F.conj(),n=2*n,axis=0)[:nlags]
    return acf/(n-np.arange(nlags)).reshape((-1,)+(1,)*(x.ndim-1))


def vacf_elements(V,nlags,typat,ntypat,U=None,AxesAreCartesian=False):#Mean over the atoms of each element of <v(0).v(t)> and of <v_k(0) v_k(t)> along the unit vectors U
    #V : (nsteps,natom,3) cartesian velocities ; returns VACF[element][0] for the total, VACF[element][1+k] along axis k, (ntypat,1+nAxes,nlags)
    nsteps,natom,_ = V.shape
    nAxes = 0 if U is None else len(U)
    VACF = np.zeros((ntypat,1+nAxes,nlags))
    batch = max(1,5000000//(3*nsteps))#Atoms are processed in batches, to limit the memory used by the FFTs
    for a0 in range(0,natom,batch):
        Vb = V[:,a0:a0+batch,:]
        tb = typat[a0:a0+batch]
        cart = acf_fft(Vb,nlags)#along x, y and z ; their sum is the total, for any cell shape
        if U is not None:
            proj = cart if AxesAreCartesian else acf_fft(Vb @ U.T,nlags)
        for itypat in range(ntypat):
            sel = tb==itypat
            VACF[itypat][0] += cart[:,sel,:].sum(axis=(1,2))
            for kk in range(nAxes):
                VACF[itypat][1+kk] += proj[:,sel,kk].sum(axis=1)
    counts = np.bincount(typat,minlength=ntypat)
    return VACF/counts[:,None,None]


def running_diffusion(VACF,dt):#Running integral of the VACF, in m^2/s : D(t) = 1/3 int <v(0).v(t)> for the total, int <v_k(0) v_k(t)> along one axis
    D = cumulative_trapezoid(VACF,dx=dt,axis=-1,initial=0)*A2FS
    D[:,0] /= 3.0
    return D


def vdos(VACF,dt,npad,taper):#Spectra of the VACFs, 2 int_0^tmax C(t) w(t) cos(2 pi nu t) dt, on the frequencies k/(2 (npad-1) dt)
    nlags = VACF.shape[-1]
    C = VACF*(0.5*(1.0+np.cos(np.pi*np.arange(nlags)/(nlags-1))) if taper == 'hann' else 1.0)#One-sided Hann taper, 1 at t=0 and 0 at tmax
    pad = np.zeros(VACF.shape[:-1]+(npad,))
    pad[...,:nlags] = C
    S = dct(pad,type=1,axis=-1)*dt#DCT-I = trapezoidal integral over -tmax .. tmax of the even function C(t)
    freq = np.arange(npad)/(2.0*(npad-1)*dt)*CM1
    return freq, S


def read_temperatures(umdfile):#Temperature of each snapshot, as written in the UMD file
    temperatures = []
    with open(umdfile,'r') as f:
        for line in f:
            if line.startswith('Temperature'):
                temperatures.append(float(line.split()[1]))
    return np.array(temperatures)


def write_columns(filename,comments,labels,columns):
    f = open(filename,'w')
    for comment in comments:
        f.write('# '+comment+'\n')
    f.write('\t'.join(labels)+'\n')
    for ii in range(len(columns[0])):
        f.write('\t'.join(str(col[ii]) for col in columns)+'\n')
    f.close()


def plot_vibrations(t,vacfnorm,tD,Drun,Dblocks_err,DT,iend,freq,g,elements,title):
    import matplotlib.pyplot as plt
    fig,(ax1,ax2,ax3) = plt.subplots(3,1,figsize=(8,11))
    for itypat,el in enumerate(elements):
        line, = ax1.plot(t,vacfnorm[itypat],lw=0.8,label=el)
        ax2.fill_between(tD,Drun[itypat]-Dblocks_err[itypat],Drun[itypat]+Dblocks_err[itypat],color=line.get_color(),alpha=0.2)
        ax2.plot(tD,Drun[itypat],color=line.get_color(),lw=1.0,label=el+' %.2e' % DT[itypat])
        ax3.plot(freq,g[itypat],color=line.get_color(),lw=0.8,label=el)
    ax3.plot(freq,g.sum(axis=0),'k',lw=1.2,label='total')
    ax1.axhline(0.0,color='k',lw=0.5)
    ax1.set_xlabel('time (fs)')
    ax1.set_ylabel('normalized VACF')
    ax1.set_title(title,fontsize='medium')
    ax1.legend(fontsize='small',ncol=5)
    ax2.axvline(tD[iend],color='r',lw=0.8)
    ax2.set_yscale('symlog',linthresh=1e-10)
    ax2.set_xlabel('time (fs)')
    ax2.set_ylabel('running D (m$^2$/s)')
    ax2.set_title('D read at '+'%.0f' % tD[iend]+' fs (red line) ; shaded : block error',fontsize='medium')
    ax2.legend(fontsize='small',ncol=3)
    gmax = g.sum(axis=0)
    above = np.nonzero(gmax > 1e-3*gmax.max())[0]
    ax3.set_xlim(0,freq[above[-1]]*1.05 if len(above) > 0 else freq[-1])
    ax3.set_xlabel('frequency (cm$^{-1}$)')
    ax3.set_ylabel('vDOS (states/cm$^{-1}$)')
    ax3.legend(fontsize='small',ncol=5)
    fig.tight_layout()
    plt.show()


def main(argv):
    umdpf.headerumd()
    start = time.time()
    umdfile = ''
    temperature = None
    firststep = 0
    window = None
    FitEnd = None
    nBlocks = 4
    taper = 'hann'
    sigma = 0.0
    AxesMode = None
    Axes = None
    plot = False
    nCores = None
    usage = 'vibrations_umd.py -f <UMD_filename> -t <Temperature> -i <InitialStep> -w <Window> -e <FitEnd> -n <nBlocks> -p <Taper> -g <Sigma> -a <Axes> -k <nCores> -v'
    try:
        opts, arg = getopt.getopt(argv,"hf:t:i:w:e:n:p:g:a:k:v",["fumdfile","tTemperature","iInitialStep","wWindow","eFitEnd","nBlocks","pTaper","gSigma","aAxes","knCores","visual"])
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('Program to compute the velocity autocorrelation functions (VACF), the vibrational densities of states (vDOS)')
            print('and the diffusion coefficients of each element, in total and along axes. Usage: ')
            print(usage)
            print('UMD_filename = input file with the trajectory, in UMD format.')
            print(' The VACFs are averaged over all the time origins, computed with Fast Fourier Transforms, and averaged over the atoms of each element.')
            print('-t : Temperature, in K, used to normalize the vDOS. Default is the mean temperature written in the UMD file')
            print('     (with it, the vDOS of all the atoms integrates to the degrees of freedom used by the code for the temperature, 3N or 3N-3).')
            print('-i : InitialStep. number of initial snapshots discarded. Default is 0.')
            print('-w : Window. length of the VACF used for the spectra, in steps. Longer windows give sharper peaks but noisier spectra.')
            print('     Default is half a block (see -n), i.e. 1/8 of the trajectory for 4 blocks ; maximum is half the trajectory.')
            print('-e : FitEnd. time lag, in steps, at which the diffusion coefficients are read on the running integral of the VACF')
            print('     D(t) = 1/3 int_0^t <v(0).v(t)> dt (and int_0^t <v_k(0) v_k(t)> dt along one axis). Default and maximum is half a block.')
            print('-n : nBlocks. the trajectory is cut into nBlocks consecutive blocks ; the spread of their D gives the error on D. Default is 4, minimum 3.')
            print("-p : Taper of the VACF before its Fourier transform, 'hann' (default ; it suppresses the ringing from the truncation of the VACF)")
            print("     or 'none'. It applies to all the spectra.")
            print('-g : Sigma. width, in cm^-1, of a Gaussian smoothing of the spectra. Default is 0 (no smoothing).')
            print('-a : Axes. decomposition of the VACFs, spectra and diffusion coefficients of each element along axes. Options are:')
            print('     cart  : along the cartesian axes x, y, z. The three components add up to the total.')
            print('     cryst : along the crystallographic axes a, b, c of the initial cell (Auto is accepted as a synonym).')
            print('             For non-orthogonal cells the three components do not add up to the total.')
            print('     a list of floats [x1,y1,z1,x2,y2,z2,...] describing any number of axes, e.g. -a [1,0,0,0,1,1].')
            print('-k : nCores. number of cores used to read the UMD file. Default is all.')
            print('-v : plot the VACFs, the running diffusion coefficients and the vDOS of each element.')
            print('Output : <base>.vels.scf.dat : time (fs), VACF of each element (and along each axis), normalized to 1 at t=0, the mass-weighted')
            print('         average of all the atoms, then the running diffusion coefficients D_El (m^2/s).')
            print('         <base>.vibr.dat : frequency (cm^-1), vDOS of each element VDos_El and their sum Total_DOS, in states per cm^-1,')
            print('         mass-weighted and normalized so that the vDOS of an element integrates to 3 N_El (T_kinetic/T).')
            print('         <base>.vibr.axes.dat (with -a) : the same along each axis, VDos_El_x ... and Total_DOS_x ..., integrating to N_El per axis.')
            print('         <base>.vacf.diffusion.dat : diffusion coefficients from the VACF (Green-Kubo), with their block errors.')
            print('         They are an independent check of those fitted on the MSD by msd_umd.py.')
            sys.exit()
        elif opt in ("-f","--fumdfile"):
            umdfile = str(arg)
            print('I will use the ',umdfile,' for input','\n')
        elif opt in ("-t","--tTemperature"):
            temperature = float(arg)
            if temperature <= 0.0:
                print('ERROR : the temperature -t has to be strictly positive.')
                sys.exit()
        elif opt in ("-i","--iInitialStep"):
            firststep = int(arg)
        elif opt in ("-w","--wWindow"):
            window = int(arg)
        elif opt in ("-e","--eFitEnd"):
            FitEnd = int(arg)
        elif opt in ("-n","--nBlocks"):
            nBlocks = int(arg)
            if nBlocks < 3:
                print('ERROR : at least 3 blocks -n are needed for the errors.')
                sys.exit()
        elif opt in ("-p","--pTaper"):
            taper = str(arg)
            if taper not in ('hann','none'):
                print("ERROR : the taper -p has to be 'hann' or 'none'.")
                sys.exit()
        elif opt in ("-g","--gSigma"):
            sigma = float(arg)
        elif opt in ("-a","--aAxes"):
            if arg in ("cart","xyz"):
                AxesMode = "cart"
            elif arg in ("cryst","abc","Auto"):
                AxesMode = "cryst"
            else:
                AxesMode = "custom"
                try:
                    Axes = [float(x) for x in ast.literal_eval(arg)]
                except (ValueError,SyntaxError,TypeError):
                    print('ERROR : axes -a have to be cart, cryst, or a list of floats such as [1,0,0,0,1,1].')
                    sys.exit()
                if len(Axes)==0 or len(Axes)%3!=0:
                    print('ERROR : the list of axes -a must contain 3 components per axis.')
                    sys.exit()
        elif opt in ("-k","--nCores"):
            nCores = int(arg)
        elif opt in ("-v","--visual"):
            plot = True

    if not os.path.isfile(umdfile):
        print('umd file ',umdfile,'does not exist')
        sys.exit()

    MyCrystal,SnapshotsValues,TimeStep,numsteps = umdpf.read_values(umdfile,"velocity","line",nCores=nCores)
    V = np.array(SnapshotsValues).reshape(len(SnapshotsValues),MyCrystal.natom,3)[firststep:]#(step, atom, xyz), in A/fs
    del SnapshotsValues
    nsteps = len(V)
    if firststep > 0:
        print('The first ',firststep,' snapshots are discarded')
    typat = np.array(MyCrystal.typat)
    elements = MyCrystal.elements
    ntypat = MyCrystal.ntypat
    counts = np.bincount(typat,minlength=ntypat)
    masses = np.array(MyCrystal.masses,dtype=float)

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
    else:
        U = None
        AxesLabels = []
    if U is not None:
        norms = np.linalg.norm(U,axis=1)
        if np.any(norms==0.0):
            print('ERROR : the axes cannot be null vectors.')
            sys.exit()
        U = U/norms[:,None]
    nAxes = len(AxesLabels)

    #Blocks give the errors on D up to half their length ; the spectra use the VACF of the whole run up to the window
    blocklen = nsteps//nBlocks
    Dwindow = blocklen//2
    if Dwindow < 10:
        print('ERROR : the trajectory is too short for ',nBlocks,' blocks.')
        sys.exit()
    if FitEnd is None or FitEnd > Dwindow-1:
        if FitEnd is not None:
            print('The end of the integral -e is beyond half a block ; it is set to ',Dwindow-1,' steps. Use fewer blocks -n for a longer integral.')
        FitEnd = Dwindow-1
    if FitEnd < 2:
        print('ERROR : the end of the integral -e must be at least 2 steps.')
        sys.exit()
    if window is None:
        window = Dwindow
    elif window > nsteps//2 or window < 10:
        window = min(max(window,10),nsteps//2)
        print('The window -w has to be between 10 steps and half the trajectory ; it is set to ',window,' steps.')
    nlags = max(window,Dwindow)

    print('Calculating the velocity autocorrelation functions')
    VACF = vacf_elements(V,nlags,typat,ntypat,U,AxesMode=="cart")
    VACFblocks = np.array([vacf_elements(V[kk*blocklen:(kk+1)*blocklen],Dwindow,typat,ntypat,U,AxesMode=="cart") for kk in range(nBlocks)])

    #Temperature : the vDOS is normalized by kB T, and m <v^2> = kB T for each degree of freedom
    Tkin = (counts*masses*VACF[:,0,0]).sum()*AMU_A2FS2/(3.0*MyCrystal.natom*KB)
    Tfile = read_temperatures(umdfile)[firststep:]
    if temperature is None:
        if len(Tfile) > 0 and Tfile.mean() > 0.0:
            temperature = Tfile.mean()
            print('Mean temperature of the UMD file : ','%.1f' % temperature,' K')
        else:
            temperature = Tkin
            print('No temperature in the UMD file : the kinetic temperature of the velocities is used')
    else:
        print('Temperature set by -t : ','%.1f' % temperature,' K')
    print('Kinetic temperature of the velocities (3N degrees of freedom) : ','%.1f' % Tkin,' K')
    print(nsteps,' snapshots of ',TimeStep,' fs ; ',nBlocks,' blocks of ',blocklen,' steps ; spectra from VACFs of ',window,' steps (',window*TimeStep,' fs), taper ',taper)

    #Vibrational densities of states, mass-weighted, in states per cm^-1
    npad = max(4*window,2)#zero-padding interpolates the spectra on a finer frequency grid
    freq, S = vdos(VACF[:,:,:window],TimeStep,npad,taper)
    g = S*(2.0*counts*masses*AMU_A2FS2/(KB*temperature))[:,None,None]/CM1
    if sigma > 0.0:
        g = gaussian_filter1d(g,sigma/(freq[1]-freq[0]),axis=-1,mode='reflect')#the spectra are even in frequency : reflecting at 0 is exact
    integrals = g.sum(axis=-1)*(freq[1]-freq[0]) - 0.5*(g[...,0]+g[...,-1])*(freq[1]-freq[0])
    print('\nIntegral of the vDOS : ','%.1f' % integrals[:,0].sum(),' for ',3*MyCrystal.natom,' degrees of freedom (3N), ',3*MyCrystal.natom-3,' without the centre of mass (3N-3)')
    print('  by element : '+', '.join(el+' '+'%.1f' % integrals[itypat,0]+' / '+str(3*counts[itypat]) for itypat,el in enumerate(elements)))

    #Diffusion coefficients : running integral of the VACF, read at FitEnd, errors from the blocks
    Drun = running_diffusion(VACF[:,:,:Dwindow],TimeStep)
    Dblocks = np.array([running_diffusion(VACFblocks[kk],TimeStep) for kk in range(nBlocks)])
    Derr = Dblocks.std(axis=0,ddof=1)/np.sqrt(nBlocks)
    tD = np.arange(Dwindow)*TimeStep
    ihalf = FitEnd//2

    if umdfile.endswith('.umd.dat'):
        base = umdfile[:-8]
    else:
        base = os.path.splitext(umdfile)[0]
    comment = 'from '+umdfile+' : '+str(nsteps)+' snapshots of '+str(TimeStep)+' fs from snapshot '+str(firststep)+', all time origins'
    axescomment = [] if U is None else ['axes (unit vectors) : '+' ; '.join(lab+' = '+np.array2string(u,precision=6,separator=',') for lab,u in zip(AxesLabels,U))]
    collabels = [[el]+[el+'_'+lab for lab in AxesLabels] for el in elements]

    #VACFs normalized to 1 at t=0, and running diffusion coefficients
    mw = (counts*masses)[:,None,None]*VACF
    average = mw.sum(axis=0)
    labels = ['Vel_AutoCorr_Fct(fs)']
    columns = [np.arange(nlags)*TimeStep]
    for itypat in range(ntypat):
        for jj in range(1+nAxes):
            labels.append(collabels[itypat][jj])
            columns.append(VACF[itypat][jj]/VACF[itypat][jj][0])
    labels += ['Average']+['Average_'+lab for lab in AxesLabels]
    columns += [average[jj]/average[jj][0] for jj in range(1+nAxes)]
    for itypat in range(ntypat):
        for jj in range(1+nAxes):
            labels.append('D_'+collabels[itypat][jj]+'(m2/s)')
            columns.append(np.concatenate([Drun[itypat][jj],np.full(nlags-Dwindow,np.nan)]))
    vaffile = base+'.vels.scf.dat'
    write_columns(vaffile,['VACF '+comment,'normalized to 1 at t=0 ; Average is mass-weighted over all the atoms ; D_El = running integral of the VACF, up to half a block']+axescomment,labels,columns)

    #vDOS : the first file keeps the format read by UMD_GUI.py (one header line, VDos_El columns, Total_DOS last)
    vibrfile = base+'.vibr.dat'
    f = open(vibrfile,'w')
    f.write('Frequency(cm^-1)\t'+'\t'.join('VDos_'+el for el in elements)+'\tTotal_DOS\n')
    total = g[:,0].sum(axis=0)
    for ii in range(npad):
        f.write(str(freq[ii])+'\t'+'\t'.join(str(g[itypat][0][ii]) for itypat in range(ntypat))+'\t'+str(total[ii])+'\n')
    f.close()
    vibraxesfile = None
    if U is not None:
        vibraxesfile = base+'.vibr.axes.dat'
        labels = ['Frequency(cm^-1)']
        columns = [freq]
        for itypat in range(ntypat):
            for kk in range(nAxes):
                labels.append('VDos_'+collabels[itypat][1+kk])
                columns.append(g[itypat][1+kk])
        for kk in range(nAxes):
            labels.append('Total_DOS_'+AxesLabels[kk])
            columns.append(g[:,1+kk].sum(axis=0))
        write_columns(vibraxesfile,['vDOS along axes '+comment,'states per cm^-1, the vDOS of an element integrates to N_El along each axis ; T = '+'%.1f' % temperature+' K, VACF of '+str(window)+' steps, taper '+taper+', smoothing '+str(sigma)+' cm^-1']+axescomment,labels,columns)

    #Diffusion coefficients
    difffile = base+'.vacf.diffusion.dat'
    f = open(difffile,'w')
    f.write('# Diffusion coefficients from the VACF '+comment+'\n')
    f.write('# D = 1/3 int_0^t <v(0).v(t)> dt (total) or int_0^t <v_k(0) v_k(t)> dt (along one axis), read at t = '+str(tD[FitEnd])+' fs ; error from '+str(nBlocks)+' blocks of '+str(blocklen)+' steps\n')
    for line in axescomment:
        f.write('# '+line+'\n')
    f.write('column\tD_(m2/s)\tD_error_(m2/s)\tD_at_half_(m2/s)\n')
    print('\nDiffusion coefficients from the VACF, read at '+str(tD[FitEnd])+' fs ; error from '+str(nBlocks)+' blocks of '+str(blocklen)+' steps')
    print('{:<10s}{:>14s}{:>14s}{:>16s}'.format('column','D (m2/s)','error','D at half'))
    unconverged = []
    for itypat in range(ntypat):
        for jj in range(1+nAxes):
            lab = collabels[itypat][jj]
            D, err, Dhalf = Drun[itypat][jj][FitEnd], Derr[itypat][jj][FitEnd], Drun[itypat][jj][ihalf]
            f.write(lab+'\t'+str(D)+'\t'+str(err)+'\t'+str(Dhalf)+'\n')
            print('{:<10s}{:>14.4e}{:>14.4e}{:>16.4e}'.format(lab,D,err,Dhalf))
            if abs(D-Dhalf) > 2.0*err:
                unconverged.append(lab)
    f.close()
    if unconverged != []:
        print('WARNING : the running D changes by more than 2 errors between half the integral and its end for ',', '.join(unconverged),' :')
        print('          D is not converged ; check the running D in ',vaffile,' (or with -v) and consider a different end -e.')
        print('          In a solid D is zero and the running D oscillates around 0 : the block error does not include this oscillation.')

    print('\nVACFs and running D printed in file ',vaffile)
    print('vDOS printed in file ',vibrfile)
    if vibraxesfile is not None:
        print('vDOS along the axes printed in file ',vibraxesfile)
    print('Diffusion coefficients printed in file ',difffile)
    print('runtime : ','%.1f' % (time.time()-start),' s')

    if plot:
        plot_vibrations(np.arange(nlags)*TimeStep,VACF[:,0]/VACF[:,0,:1],tD,Drun[:,0],Derr[:,0],Drun[:,0,FitEnd],FitEnd,freq,g[:,0],elements,os.path.basename(umdfile))

    return [vibrfile,freq,g,elements,Drun[:,:,FitEnd],Derr[:,:,FitEnd]]


if __name__ == "__main__":
    main(sys.argv[1:])
