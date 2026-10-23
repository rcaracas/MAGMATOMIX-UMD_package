#!/usr/bin/env python3
# -*- coding: utf-8 -*-


### Authors : Razvan Caracas, Adrien Saurety, Emma Stoutenburg
### Green-Kubo with the 5 components of the traceless stress tensor, bias correction and block errors : Razvan Caracas
import sys, getopt, os
import numpy as np
from scipy.integrate import cumulative_trapezoid
import umd_processes_fast as umdpf


COMPONENTS = ['xy','yz','zx','xx','yy','zz']
#Daivis & Evans : eta = V/(10 kB T) int sum_ab <N_ab(0) N_ab(t)> dt, with N the traceless symmetric stress tensor (5 independent components)
#each off-diagonal component appears twice (ab and ba), each diagonal one once ; written as weights on per-component viscosities
#for which the diagonal ones are scaled by 3/4, since <N_aa N_aa> = 4/3 <N_ab N_ab> in an isotropic liquid
WEIGHTS = {'all':np.array([0.2,0.2,0.2,2/15,2/15,2/15]), 'offdiag':np.array([1/3,1/3,1/3,0.0,0.0,0.0])}
DIAGSCALE = np.array([1.0,1.0,1.0,0.75,0.75,0.75])


def gk_prefactor(volume,temperature):#V/(kB T) in Pa.s per GPa^2.fs, for V in A^3 and T in K : 1e-30 * 1e18 * 1e-15 / 1.380649e-23 = 1e-4/1.380649
    return 1.0E-04*volume/(1.380649*temperature)


def shear_stresses(stress):#Traceless stress tensor N = sigma - p I, from the UMD StressTensor (xx yy zz xy yz zx, in GPa)
    p = stress[:,:3].mean(axis=1)
    return np.column_stack([stress[:,3],stress[:,4],stress[:,5],stress[:,0]-p,stress[:,1]-p,stress[:,2]-p])


def autocorrelation(x,nlags):#Autocorrelation of each column of x (time along axis 0), averaged over all the time origins, by FFT
    #x must already be the fluctuations : the mean of the whole trajectory is removed once, also for the blocks, so that they all have the same bias
    n = x.shape[0]
    F = np.fft.rfft(x,n=2*n,axis=0)#Zero-padding to 2n gives the non-cyclic correlation
    acf = np.fft.irfft(F*F.conj(),n=2*n,axis=0)[:nlags]
    return acf/(n-np.arange(nlags))[:,None]


def running_viscosities(acf,dt,pref):#Running Green-Kubo integral of each component, in Pa.s, the diagonal ones scaled to estimate the shear viscosity
    return cumulative_trapezoid(acf,dx=dt,axis=0,initial=0)*pref*DIAGSCALE


def viscosity_estimate(t,eta,etablocks,iend,runtime):#Viscosity at the end of the integral, t[iend], and its error from the blocks
    #removing the mean of a trajectory of length runtime lowers its autocorrelation by about 2 I/runtime at all lags, I being its full integral,
    #so that eta(t) is lowered by 2 eta t/runtime (about 10% at 1 ps for a 16 ps run) ; the self-consistent correction at the end of the integral is
    #eta = eta_raw/(1 - 2 t_end/runtime), and the correction 2 eta t/runtime is the same for all the components and all the times
    factor = 1.0/(1.0-2.0*t[iend]/runtime)
    value = eta[iend]*factor
    correction = 2.0*value*t/runtime
    etaerr = etablocks.std(axis=0,ddof=1)/np.sqrt(len(etablocks))*factor#standard error of the whole-run eta(t), from the spread of the blocks
    #convergence : eta(t) should no longer rise between half the integral and its end
    ihalf = iend//2
    rise = (eta[iend]+correction[iend]) - (eta[ihalf]+correction[ihalf])
    return {'eta':value, 'error':etaerr[iend], 'iend':iend, 'tend':t[iend], 'correction':correction, 'etaerr':etaerr,
            'converged':rise < 2.0*etaerr[iend], 'rise':rise}


def plot_viscosity(t,acf,etacomp,eta,etaerr,res,title):
    import matplotlib.pyplot as plt
    fig,(ax1,ax2) = plt.subplots(2,1,figsize=(8,8),sharex=True)
    for ic,comp in enumerate(COMPONENTS):
        ax1.plot(t,acf[:,ic]*DIAGSCALE[ic],lw=0.8,label=comp+(' (x3/4)' if ic > 2 else ''))
        ax2.plot(t,etacomp[:,ic]*1e3,lw=0.6,alpha=0.6,label=comp)
    ax1.axhline(0.0,color='k',lw=0.5)
    ax1.set_ylabel('stress autocorrelation (GPa$^2$)')
    ax1.legend(fontsize='small',ncol=3)
    ax1.set_title(title,fontsize='medium')
    ax2.fill_between(t,(eta-etaerr)*1e3,(eta+etaerr)*1e3,color='0.8',label='block error')
    ax2.plot(t,eta*1e3,'k',lw=1.5,label='viscosity')
    ax2.axvline(res['tend'],color='r',lw=0.8,label='end of the integral')
    ax2.axhline(res['eta']*1e3,color='r',lw=0.8,ls=':')
    ax2.set_xlabel('time (fs)')
    ax2.set_ylabel('running viscosity (mPa.s)')
    ax2.set_title('viscosity = '+'%.3g' % (res['eta']*1e3)+' +/- '+'%.2g' % (res['error']*1e3)+' mPa.s at '+'%.0f' % res['tend']+' fs',fontsize='medium')
    ax2.legend(fontsize='small',ncol=3)
    fig.tight_layout()
    plt.show()


def main(argv):
    umdpf.headerumd()
    UMDfile = ''
    firststep = 0
    FitEnd = None
    nBlocks = 4
    Components = 'all'
    plot = False
    usage = 'viscosity_umd.py -f <umdfile> -i <InitialSteps> -e <End> -n <nBlocks> -c <Components> -v'
    try:
        opts, args = getopt.getopt(argv,"hf:i:e:n:c:v",["fumdfile","iInitialSteps","eEnd","nBlocks","cComponents","visual"])
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('viscosity_umd.py program to compute the shear viscosity with the Green-Kubo relation, from the stress tensor of a UMD file')
            print(usage)
            print('-f : UMD file')
            print('-i : number of initial snapshots to skip (default 0)')
            print('-e : end of the Green-Kubo integral, in steps, where the viscosity is read (default and maximum: half a block, i.e. 1/8 of the run for 4 blocks)')
            print('-n : number of consecutive blocks of the trajectory, used for the error (default 4, minimum 3) ; fewer blocks allow a longer integral')
            print("-c : 'all' (default) uses the 5 independent components of the traceless stress tensor, for an isotropic liquid (Daivis & Evans) ;")
            print("     'offdiag' uses only the 3 off-diagonal components xy, yz, zx, e.g. for anisotropic systems")
            print('-v : display the stress autocorrelations and the running viscosity, with its error and the end of the integral (default: no display)')
            print('The mean stress is removed, and the running Green-Kubo integral eta(t) is corrected for the bias this causes (+2 eta t / run length).')
            print('The viscosity is eta(t) at the end of the integral, and its error is the standard error of eta(t) from the spread of the blocks.')
            print('By default the integral runs as far as the blocks allow : the stress autocorrelation of melts often has a small but long tail')
            print('that carries much of the viscosity, which automatic windows tend to cut. A warning is printed if eta(t) still rises by more than')
            print('2 errors over the second half of the integral : then use a longer run or fewer blocks. Check the plateau with -v, and use -e')
            print('to stop earlier for fast-relaxing liquids, whose eta(t) is noisier at long times. With few blocks the error may be underestimated')
            print('by 10-30%.')
            print('Output : <base>.visc.dat with the autocorrelation and the running viscosity of each component, the combined running')
            print('viscosity and its error from the blocks, as a function of time.')
            sys.exit()
        elif opt in ("-f","--fumdfile"):
            UMDfile = str(arg)
        elif opt in ("-i","--iInitialSteps"):
            firststep = int(arg)
        elif opt in ("-e","--eEnd"):
            FitEnd = int(arg)
        elif opt in ("-n","--nBlocks"):
            nBlocks = int(arg)
            if nBlocks < 3:
                print('ERROR : at least 3 blocks -n are needed for the error.')
                sys.exit()
        elif opt in ("-c","--cComponents"):
            Components = str(arg)
            if Components not in WEIGHTS:
                print("ERROR : components -c have to be 'all' or 'offdiag'.")
                sys.exit()
        elif opt in ("-v","--visual"):
            plot = True
    if not os.path.isfile(UMDfile):
        print('the umdfile ',UMDfile,' does not exist')
        sys.exit()

    MyCrystal, AllSnapshots, TimeStep = umdpf.read_stresses_4visc(UMDfile)
    if firststep > 0:
        print('The first ',firststep,' snapshots are discarded')
    AllSnapshots = AllSnapshots[firststep:]
    nsteps = len(AllSnapshots)
    stress = np.array([snapshot.stress for snapshot in AllSnapshots])
    temperature = np.mean([snapshot.temperature for snapshot in AllSnapshots])
    volume = np.mean([snapshot.cellvolume for snapshot in AllSnapshots])
    pref = gk_prefactor(volume,temperature)
    weights = WEIGHTS[Components]
    print(nsteps,' snapshots of ',TimeStep,' fs ; mean temperature ','%.1f' % temperature,' K ; mean volume ','%.2f' % volume,' A^3')

    blocklen = nsteps//nBlocks
    window = blocklen//2#the blocks give the error up to half their length
    if FitEnd is None or FitEnd > window-1:
        if FitEnd is not None:
            print('The end of the integral -e is beyond half a block ; it is set to ',window-1,' steps. Use fewer blocks -n for a longer integral.')
        FitEnd = window-1
    if window < 10:
        print('ERROR : the trajectory is too short for ',nBlocks,' blocks.')
        sys.exit()
    t = np.arange(window)*TimeStep

    N = shear_stresses(stress)
    N = N - N.mean(axis=0)#Green-Kubo needs the fluctuations : a residual mean stress would make the integral grow linearly
    acf = autocorrelation(N,window)
    etacomp = running_viscosities(acf,TimeStep,pref)
    etablocks = np.array([running_viscosities(autocorrelation(N[kk*blocklen:(kk+1)*blocklen],window),TimeStep,pref) @ weights for kk in range(nBlocks)])
    res = viscosity_estimate(t,etacomp @ weights,etablocks,FitEnd,nsteps*TimeStep)
    etacomp = etacomp + res['correction'][:,None]#the bias correction, in Pa.s, is the same for all components
    eta = etacomp @ weights
    etaerr = res['etaerr']

    print('Components : ',Components,' ; ',nBlocks,' blocks of ',blocklen,' steps ; time lags up to ',t[-1],' fs')
    print('Running viscosity (mPa.s) : '+', '.join('%.0f' % t[ii]+' fs : '+'%.3g' % (eta[ii]*1e3)+' +/- '+'%.2g' % (etaerr[ii]*1e3) for ii in (FitEnd//4,FitEnd//2,3*FitEnd//4,FitEnd)))
    print('Running viscosity at the end of the integral, by component (mPa.s) : '+', '.join(comp+' '+'%.3g' % (etacomp[FitEnd,ic]*1e3) for ic,comp in enumerate(COMPONENTS) if weights[ic] > 0))
    print('Viscosity = ','%.4g' % res['eta'],' +/- ','%.2g' % res['error'],' Pa.s  (','%.4g' % (res['eta']*1e3),' +/- ','%.2g' % (res['error']*1e3),' mPa.s)')
    print('Viscosity read at the end of the integral, ','%.0f' % res['tend'],' fs')
    if not res['converged']:
        print('WARNING : the running viscosity still rises by ','%.2g' % (res['rise']*1e3),' mPa.s over the second half of the integral (more than 2 errors) :')
        print('          it has not converged, and the value is a lower bound ; use fewer blocks -n (longer integral) or a longer run.')

    if UMDfile.endswith('.umd.dat'):
        viscfile = UMDfile[:-8]+'.visc.dat'
    else:
        viscfile = os.path.splitext(UMDfile)[0]+'.visc.dat'
    with open(viscfile,'w') as f:
        f.write('# Green-Kubo viscosity from '+UMDfile+' : '+str(nsteps)+' snapshots from '+str(firststep)+', T = '+'%.1f' % temperature+' K, V = '+'%.2f' % volume+' A^3, components '+Components+', '+str(nBlocks)+' blocks of '+str(blocklen)+' steps\n')
        f.write('# viscosity = '+str(res['eta'])+' +/- '+str(res['error'])+' Pa.s, at the end of the integral '+str(res['tend'])+' fs'+('' if res['converged'] else ' (NOT converged : lower bound)')+'\n')
        f.write('# acf_ab : autocorrelation of the traceless stress N_ab (GPa^2) ; eta_ab : its running Green-Kubo integral (Pa.s, x3/4 for the diagonal ones, bias-corrected) ; eta : weighted mean, eta_err : error from the blocks\n')
        f.write('time_(fs)\t'+'\t'.join('acf_'+comp for comp in COMPONENTS)+'\t'+'\t'.join('eta_'+comp for comp in COMPONENTS)+'\teta\teta_err\n')
        for ii in range(window):
            f.write(str(t[ii])+'\t'+'\t'.join(str(x) for x in acf[ii])+'\t'+'\t'.join(str(x) for x in etacomp[ii])+'\t'+str(eta[ii])+'\t'+str(etaerr[ii])+'\n')
    print('Autocorrelations and running viscosities printed in file ',viscfile)

    if plot:
        plot_viscosity(t,acf,etacomp,eta,etaerr,res,UMDfile+' : '+Components+' components')
    return res


if __name__ == "__main__":
    main(sys.argv[1:])
