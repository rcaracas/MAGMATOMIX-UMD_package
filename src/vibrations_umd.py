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


def taper_weights(nlags,taper):#One-sided Hann taper, 1 at t=0 and 0 at tmax, or no taper
    if taper == 'hann':
        return 0.5*(1.0+np.cos(np.pi*np.arange(nlags)/(nlags-1)))
    return np.ones(nlags)


def vdos(VACF,dt,npad,taper):#Spectra of the VACFs, 2 int_0^tmax C(t) w(t) cos(2 pi nu t) dt, on the frequencies k/(2 (npad-1) dt)
    nlags = VACF.shape[-1]
    C = VACF*taper_weights(nlags,taper)
    pad = np.zeros(VACF.shape[:-1]+(npad,))
    pad[...,:nlags] = C
    S = dct(pad,type=1,axis=-1)*dt#DCT-I = trapezoidal integral over -tmax .. tmax of the even function C(t)
    freq = np.arange(npad)/(2.0*(npad-1)*dt)*CM1
    return freq, S


def vdos_states(VACF,dt,npad,taper,weights,sigma):#vDOS in states per cm^-1 ; weights = 2 N_El m_El/(kB T) of each element, along the axis before the last two
    freq, S = vdos(VACF,dt,npad,taper)
    g = S*weights[:,None,None]/CM1
    if sigma > 0.0:
        g = gaussian_filter1d(g,sigma/(freq[1]-freq[0]),axis=-1,mode='reflect')#the spectra are even in frequency : reflecting at 0 is exact
    return freq, g


def spectral_variance(L,window,taper):#Variance of a spectrum from a run of L steps, up to a constant : the VACF at lag t, averaged over L-t origins, has a variance ~1/(L-t)
    w = taper_weights(window,taper)
    return 1.0/L+(w[1:]**2/(L-np.arange(1,window))).sum()


def vdos_errors(gblocks,g,nsteps,blocklen,window,taper,fwhm):#Error of the vDOS g of the whole run from the scatter of the vDOS of the blocks (along axis 0)
    #The variance of a block spectrum is scaled to the whole run, then averaged over the resolution width of the spectra (fwhm, in frequency points)
    var = gblocks.var(axis=0,ddof=1)*spectral_variance(nsteps,window,taper)/spectral_variance(blocklen,window,taper)
    smooth = lambda x: gaussian_filter1d(x,fwhm/2.3548,axis=-1,mode='reflect')
    if taper == 'hann':#the relative variance is nearly constant across a peak : averaging it keeps the larger error at the top of the peaks
        return np.abs(g)*np.sqrt(smooth(var)/np.maximum(smooth(g**2),1e-300))
    return np.sqrt(smooth(var))#without taper the vDOS crosses zero


def vdos_with_errors(V,gall,window,nBlocks,typat,ntypat,U,AxesAreCartesian,dt,npad,taper,weights,sigma,fwhm,VACFblocks=None):
    #Errors of gall, the vDOS of the elements then of their sum, from blocks at least one window long ; VACFblocks are reused when they fit
    nsteps = len(V)
    nSpec = min(nBlocks,nsteps//window)
    specblock = nsteps//nSpec
    if VACFblocks is not None and len(VACFblocks) == nSpec and VACFblocks.shape[-1] >= window:
        VACFspec = VACFblocks[...,:window]
    else:
        VACFspec = np.array([vacf_elements(V[kk*specblock:(kk+1)*specblock],window,typat,ntypat,U,AxesAreCartesian) for kk in range(nSpec)])
    freq, gblocks = vdos_states(VACFspec,dt,npad,taper,weights,sigma)
    gblocks = np.concatenate([gblocks,gblocks.sum(axis=1,keepdims=True)],axis=1)
    return vdos_errors(gblocks,gall,nsteps,specblock,window,taper,fwhm/(freq[1]-freq[0])), nSpec, specblock


def mean_error(x,with_tau=False):#Statistical error of the mean of the time series x, from its integrated autocorrelation time tau (in steps, also returned with with_tau)
    t = np.arange(len(x))
    r = x-np.polyval(np.polyfit(t,x,1),t)#a drift would inflate the autocorrelation time and hide itself
    if r.std() == 0.0:
        return (0.0,0.5) if with_tau else 0.0
    rho = acf_fft(r,len(r)//2)
    rho = rho/rho[0]
    cut = np.argmax(rho < 0.0) or len(rho)#sum up to the first zero
    tau = 0.5+rho[1:cut].sum()
    err = r.std()*np.sqrt(2.0*tau/len(x))
    return (err,tau) if with_tau else err


def check_zero_velocities(V,firststep):#Snapshots with all the velocities exactly zero : not physical, written by a parser that lost the previous positions
    zero = np.nonzero(~V.any(axis=(1,2)))[0]
    print('\nCheck of the velocities : snapshots with all the velocities exactly zero')
    if len(zero) == 0:
        print('OK : none.')
        return
    print('WARNING : ',len(zero),' snapshots (','%.1f' % (100.0*len(zero)/len(V)),'% of the run) have all their velocities exactly zero, e.g. snapshots (counted from 0, as -i)',
          ', '.join(str(ii+firststep) for ii in zero[:8])+(' ...' if len(zero) > 8 else ''))
    print('          They are written by VaspParser.py versions before October 23, 2026 on the learning steps of on-the-fly ML runs, which')
    print('          also miss the crossings of the cell faces at these steps (wrong unwrapped positions, hence MSD). They lower the VACF,')
    print('          the vDOS and the kinetic temperature, and can look like a drift below. Parse the OUTCAR again with the current VaspParser.py,')
    print('          or repair the UMD file with repair_umd.py.')


def check_temperatures(V,typat,ntypat,masses,elements,nBlocks):#Kinetic temperature of each element by block : drifts, and equipartition between the elements
    counts = np.bincount(typat,minlength=ntypat)
    Tel = np.array([(masses[itypat]*V[:,typat==itypat]**2).sum(axis=(1,2)) for itypat in range(ntypat)])*AMU_A2FS2/(3.0*KB*counts[:,None])
    Tall = (counts[:,None]*Tel).sum(axis=0)/counts.sum()
    blocklen = len(V)//nBlocks
    print('\nCheck of the equilibration : kinetic temperature (K) of each element by block, then over the run')
    print('{:<6s}'.format('')+''.join('{:>10s}'.format('block '+str(kk+1)) for kk in range(nBlocks))+'{:>12s}{:>10s}{:>14s}'.format('run','error','(T-T_all)/err'))
    drifting = []
    hotcold = []
    for itypat,el in enumerate(elements+['all']):
        T = Tel[itypat] if itypat < ntypat else Tall
        Tb = np.array([T[kk*blocklen:(kk+1)*blocklen].mean() for kk in range(nBlocks)])
        err = mean_error(T)
        line = '{:<6s}'.format(el)+''.join('{:>10.1f}'.format(x) for x in Tb)+'{:>12.1f}{:>10.1f}'.format(T.mean(),err)
        if err > 0.0 and np.abs(Tb-T.mean()).max() > 4.0*err*np.sqrt(nBlocks):#error of a block mean = error of the run * sqrt(nBlocks)
            drifting.append(el)
        if itypat < ntypat and counts[itypat] < counts.sum():
            errdiff = mean_error(T-Tall)#the element contributes to T_all : the error of their difference
            z = (T-Tall).mean()/errdiff if errdiff > 0.0 else 0.0
            line += '{:>14.1f}'.format(z)
            if abs(z) > 4.0:
                hotcold.append(el)
        print(line)
    if drifting != []:
        print('WARNING : the kinetic temperature of ',', '.join(drifting),' changes between the blocks by more than 4 errors of a block :')
        print('          the run is not stationary. The vDOS and D are averages over a drift, and their errors do not include it.')
        print('          Discard the beginning of the run with -i, or check the thermostat.')
    if hotcold != []:
        print('WARNING : the kinetic temperature of ',', '.join(hotcold),' differs from that of all the atoms by more than 4 errors :')
        print('          the elements are not in equipartition (not equilibrated, or slow energy exchange as in a nearly harmonic crystal).')
        print('          The vDOS of an element integrates to 3 N_El T_El/T, not to 3 N_El.')
    if drifting == [] and hotcold == []:
        print('OK : no drift between the blocks, and the elements are in equipartition (within 4 errors).')


def check_halves(V,gall,window,nBlocks,typat,ntypat,elements,dt,npad,taper,weights,sigma,fwhm):#vDOS of the two halves of the run : do they differ as much as their errors say ?
    half = len(V)//2
    print('\nCheck of the errors of the vDOS : vDOS of the two halves of the run, with the same window, compared with their errors')
    if half//window < 2:
        print('Skipped : the window is longer than 1/4 of the run, so a half cannot be cut into 2 blocks.')
        return
    g12 = []
    e12 = []
    for part in (V[:half],V[half:2*half]):
        _, g = vdos_states(vacf_elements(part,window,typat,ntypat),dt,npad,taper,weights,sigma)
        g = np.concatenate([g[:,0],g[:,0].sum(axis=0,keepdims=True)])
        gerr = vdos_with_errors(part,g[:,None],window,nBlocks,typat,ntypat,None,False,dt,npad,taper,weights,sigma,fwhm)[0][:,0]
        g12.append(g)
        e12.append(gerr)
    sel = gall[-1][0] > 1.0E-3*gall[-1][0].max()#frequencies where there are vibrations
    ratio = np.sqrt(((g12[0]-g12[1])**2)[:,sel].mean(axis=1)/(e12[0]**2+e12[1]**2)[:,sel].mean(axis=1))
    print('rms(difference)/rms(error), 1 when the errors are right : total '+'%.2f' % ratio[-1]+' ; '+', '.join(el+' '+'%.2f' % ratio[itypat] for itypat,el in enumerate(elements)))
    if ratio[-1] > 1.5:
        print('WARNING : the two halves differ ','%.1f' % ratio[-1],' times more than their errors : the error band of the vDOS is too small by about this factor.')
        print('          Usual causes : modes living as long as the blocks (crystals, sharp molecular modes), or a drift of the run (see above).')
    else:
        print('OK : the two halves agree within their errors.')


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


def plot_vibrations(t,vacfnorm,freq,g,gerr,elements,title):#g, gerr : vDOS and errors of the elements, then of their sum
    import matplotlib.pyplot as plt
    fig,(ax1,ax2) = plt.subplots(2,1,figsize=(8,8))
    for itypat,el in enumerate(elements):
        line, = ax1.plot(t,vacfnorm[itypat],lw=0.8,label=el)
        ax2.fill_between(freq,g[itypat]-gerr[itypat],g[itypat]+gerr[itypat],color=line.get_color(),alpha=0.25,lw=0)
        ax2.plot(freq,g[itypat],color=line.get_color(),lw=0.8,label=el)
    ax2.fill_between(freq,g[-1]-gerr[-1],g[-1]+gerr[-1],color='k',alpha=0.2,lw=0)
    ax2.plot(freq,g[-1],'k',lw=1.2,label='total')
    ax1.axhline(0.0,color='k',lw=0.5)
    ax1.set_xlabel('time (fs)')
    ax1.set_ylabel('normalized VACF')
    ax1.set_title(title,fontsize='medium')
    ax1.legend(fontsize='small',ncol=5)
    gmax = g[-1]
    above = np.nonzero(gmax > 1e-3*gmax.max())[0]
    ax2.set_xlim(0,freq[above[-1]]*1.05 if len(above) > 0 else freq[-1])
    ax2.set_xlabel('frequency (cm$^{-1}$)')
    ax2.set_ylabel('vDOS (states/cm$^{-1}$)')
    ax2.set_title('shaded : one standard error',fontsize='medium')
    ax2.legend(fontsize='small',ncol=5)
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
            print('-f all : every *.umd.dat file of the current folder, one after the other.')
            print('-w : Window. length of the VACF used for the spectra, in steps. Longer windows give sharper peaks but noisier spectra.')
            print('     Default is half a block (see -n), i.e. 1/8 of the trajectory for 4 blocks ; maximum is half the trajectory.')
            print('     The width (FWHM) of the narrowest peaks is ~24000/T cm^-1, T = window in fs (see -p) : to separate two peaks Delta cm^-1')
            print('     apart, T has to be longer than ~24000/Delta fs. Compare the peaks with the error band (-v, <base>.vibr.err.dat) :')
            print('     a peak that does not rise above the band is noise.')
            print('-e : FitEnd. time lag, in steps, at which the diffusion coefficients are read on the running integral of the VACF')
            print('     D(t) = 1/3 int_0^t <v(0).v(t)> dt (and int_0^t <v_k(0) v_k(t)> dt along one axis). Default and maximum is half a block.')
            print('-n : nBlocks. the trajectory is cut into nBlocks consecutive blocks ; the spread of their D gives the error on D. Default is 4, minimum 3.')
            print('     The error of the vDOS comes from the spread of the vDOS of min(nBlocks, trajectory/window) blocks, each at least one window long,')
            print('     scaled to the whole run. It matches the real scatter for melts, but it is too small (by ~2x in a crystal at 1000 K) when the')
            print('     modes live as long as the blocks (crystals, sharp molecular modes), and it does not include drifts of a run that is not')
            print('     equilibrated. With -p none it is unreliable at the peaks. With only 2 blocks (window longer than 1/3 of the run) the band is rough.')
            print('     More blocks give a more precise error, longer blocks a less biased one : the blocks must be much longer than the lifetime')
            print('     of the vibrations. For melts (VACF decaying in a few 100 fs) the default 4 is good. For crystals and molecules with')
            print('     sharp modes fewer blocks help only a little (ZnS at 1000 K : ratio 2.3 with -n 4, 2.1 with -n 3) : the band stays too')
            print('     small, and should be multiplied by the ratio given by the check of the two halves at the end of the run.')
            print("-p : Taper of the VACF before its Fourier transform, 'hann' (default) or 'none'. It applies to all the spectra.")
            print("     none : the VACF is cut abruptly at the end of the window T. Each peak of the spectrum then has ripples on its sides")
            print("            (side lobes up to 22% of the peak, spaced by 1/(2T)), which can look like spurious peaks, and the vDOS can be negative.")
            print("     hann : the VACF is multiplied by 1/2 (1+cos(pi t/T)), going smoothly from 1 at t=0 to 0 at T. The side lobes drop")
            print("            to ~3% and the noisy long-time tail of the VACF gets less weight, but the peaks are ~20% wider.")
            print("     Peak width (FWHM) : ~24000/T cm^-1 with hann, ~20000/T without, T in fs. The normalization of the vDOS is the same.")
            print("     For sharper peaks, use a longer window -w rather than no taper.")
            print('-g : Sigma. width, in cm^-1, of a Gaussian smoothing of the spectra. Default is 0 (no smoothing).')
            print('-a : Axes. decomposition of the VACFs, spectra and diffusion coefficients of each element along axes. Options are:')
            print('     cart  : along the cartesian axes x, y, z. The three components add up to the total.')
            print('     cryst : along the crystallographic axes a, b, c of the initial cell (Auto is accepted as a synonym).')
            print('             For non-orthogonal cells the three components do not add up to the total.')
            print('     a list of floats [x1,y1,z1,x2,y2,z2,...] describing any number of axes, e.g. -a [1,0,0,0,1,1].')
            print('-k : nCores. number of cores used to read the UMD file. Default is all.')
            print('-v : plot the VACFs and the vDOS of each element, with the error band of the vDOS (one standard error).')
            print('Output : <base>.vels.scf.dat : time (fs), VACF of each element (and along each axis), normalized to 1 at t=0, the mass-weighted')
            print('         average of all the atoms, then the running diffusion coefficients D_El (m^2/s).')
            print('         <base>.vibr.dat : frequency (cm^-1), vDOS of each element VDos_El and their sum Total_DOS, in states per cm^-1,')
            print('         mass-weighted and normalized so that the vDOS of an element integrates to 3 N_El (T_kinetic/T).')
            print('         <base>.vibr.axes.dat (with -a) : the same along each axis, VDos_El_x ... and Total_DOS_x ..., integrating to N_El per axis.')
            print('         <base>.vibr.err.dat : frequency, then each vDOS (also along the axes) followed by its error, VDos_El VDos_El_err ...')
            print('         vdos_thermo.py computes the harmonic thermodynamics and the quantum corrections from these vDOS.')
            print('         <base>.vacf.diffusion.dat : diffusion coefficients from the VACF (Green-Kubo), with their block errors.')
            print('         They are an independent check of those fitted on the MSD by msd_umd.py.')
            print('Checks printed in the terminal :')
            print('  velocities : snapshots with all the velocities exactly zero, written by VaspParser.py before October 23, 2026 on the')
            print('     learning steps of on-the-fly ML runs : parse the OUTCAR again, or use repair_umd.py.')
            print('  equilibration : kinetic temperature of each element in each block (-n) and over the run, with its statistical error')
            print('     (from the autocorrelation time, after removing a linear trend). Warnings when the temperature of an element changes')
            print('     between blocks by more than 4 errors (the run is not stationary : use -i), or when an element is hotter or colder')
            print('     than all the atoms by more than 4 errors (no equipartition).')
            print('  errors of the vDOS : the vDOS of the two halves of the run, each with its own errors, are compared. The ratio')
            print('     rms(difference)/rms(error) is ~1 when the errors are right ; a warning above 1.5 means the error band is too small')
            print('     by about this factor. Skipped when the window is longer than 1/4 of the run.')
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
    base = umdfile[:-8] if umdfile.endswith('.umd.dat') else os.path.splitext(umdfile)[0]
    umdpf.remove_former_outputs([base+ext for ext in ('.vels.scf.dat','.vibr.dat','.vibr.axes.dat','.vibr.err.dat','.vacf.diffusion.dat')])#of a former run

    MyCrystal,SnapshotsValues,TimeStep,numsteps = umdpf.read_values(umdfile,"velocity","line",nCores=nCores)
    V = np.array(SnapshotsValues).reshape(len(SnapshotsValues),MyCrystal.natom,3)#(step, atom, xyz), in A/fs
    del SnapshotsValues
    #The parsers have no previous position for the first ionic step : its velocities are zero (current VaspParser) or the positions (older versions)
    v2 = (V[:101]**2).sum(axis=(1,2))
    if firststep == 0 and len(v2) > 1 and not (0.1*np.median(v2[1:]) < v2[0] < 10.0*np.median(v2[1:])):
        print('The velocities of the first snapshot are not those of the trajectory (written by the parser without a previous step) : it is discarded')
        firststep = 1
    V = V[firststep:]
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
    weights = 2.0*counts*masses*AMU_A2FS2/(KB*temperature)
    freq, g = vdos_states(VACF[:,:,:window],TimeStep,npad,taper,weights,sigma)
    integrals = g.sum(axis=-1)*(freq[1]-freq[0]) - 0.5*(g[...,0]+g[...,-1])*(freq[1]-freq[0])
    print('\nIntegral of the vDOS : ','%.1f' % integrals[:,0].sum(),' for ',3*MyCrystal.natom,' degrees of freedom (3N), ',3*MyCrystal.natom-3,' without the centre of mass (3N-3)')
    print('  by element : '+', '.join(el+' '+'%.1f' % integrals[itypat,0]+' / '+str(3*counts[itypat]) for itypat,el in enumerate(elements)))

    #Errors of the vDOS : spread of the vDOS of blocks at least one window long, scaled to the whole run
    fwhm = (24000.0 if taper == 'hann' else 20000.0)/(window*TimeStep)
    gall = np.concatenate([g,g.sum(axis=0,keepdims=True)])#the elements, then their sum
    gerr, nSpec, specblock = vdos_with_errors(V,gall,window,nBlocks,typat,ntypat,U,AxesMode=="cart",TimeStep,npad,taper,weights,sigma,fwhm,VACFblocks)
    print('Resolution of the spectra (FWHM) : ','%.1f' % fwhm,' cm^-1 ; errors of the vDOS from ',nSpec,' blocks of ',specblock,' steps')
    if nSpec < 3:
        print('  with only 2 blocks the error band is rough ; a window -w shorter than 1/3 of the run gives more blocks')

    #Diffusion coefficients : running integral of the VACF, read at FitEnd, errors from the blocks
    Drun = running_diffusion(VACF[:,:,:Dwindow],TimeStep)
    Dblocks = np.array([running_diffusion(VACFblocks[kk],TimeStep) for kk in range(nBlocks)])
    Derr = Dblocks.std(axis=0,ddof=1)/np.sqrt(nBlocks)
    tD = np.arange(Dwindow)*TimeStep
    ihalf = FitEnd//2

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

    #vDOS with their errors
    vibrerrfile = base+'.vibr.err.dat'
    labels = ['Frequency(cm^-1)']
    columns = [freq]
    for itypat in range(ntypat+1):
        for jj in range(1+nAxes):
            lab = ('VDos_'+collabels[itypat][jj]) if itypat < ntypat else ('Total_DOS'+('_'+AxesLabels[jj-1] if jj > 0 else ''))
            labels += [lab,lab+'_err']
            columns += [gall[itypat][jj],gerr[itypat][jj]]
    write_columns(vibrerrfile,['vDOS and their errors '+comment,'states per cm^-1 ; T = '+'%.1f' % temperature+' K, VACF of '+str(window)+' steps, taper '+taper+', smoothing '+str(sigma)+' cm^-1, resolution (FWHM) '+'%.1f' % fwhm+' cm^-1',
                               'error = one standard error, from the spread of the vDOS of '+str(nSpec)+' blocks of '+str(specblock)+' steps scaled to the whole run, averaged over the FWHM']+axescomment,labels,columns)

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
        print('          D is not converged ; check the running D in ',vaffile,' and consider a different end -e.')
        print('          In a solid D is zero and the running D oscillates around 0 : the block error does not include this oscillation.')

    check_zero_velocities(V,firststep)
    check_temperatures(V,typat,ntypat,masses,list(elements),nBlocks)
    check_halves(V,gall,window,nBlocks,typat,ntypat,list(elements),TimeStep,npad,taper,weights,sigma,fwhm)

    print('\nVACFs and running D printed in file ',vaffile)
    print('vDOS printed in file ',vibrfile)
    print('vDOS with their errors printed in file ',vibrerrfile)
    if vibraxesfile is not None:
        print('vDOS along the axes printed in file ',vibraxesfile)
    print('Diffusion coefficients printed in file ',difffile)
    print('runtime : ','%.1f' % (time.time()-start),' s')

    if plot:
        plot_vibrations(np.arange(nlags)*TimeStep,VACF[:,0]/VACF[:,0,:1],freq,gall[:,0],gerr[:,0],elements,os.path.basename(umdfile))

    return [vibrfile,freq,g,elements,Drun[:,:,FitEnd],Derr[:,:,FitEnd],gerr]


if __name__ == "__main__":
    if not umdpf.run_on_all(sys.argv[1:],main):#-f all : every UMD file of the folder
        main(sys.argv[1:])
