#!/usr/bin/env python3
###
##AUTHORS: RAZVAN CARACAS
###

"""Diffusion coefficients of the clusters (chemical species) from the cluster MSD files written by msd_clusters_umd.py"""
import sys
import getopt
import glob
import os
import numpy as np
from analyze_msd import read_columns, is_cluster_file, calculation_diffusivities, diffusion_exponent, plot_fits, regime_warning


def read_cluster_msd(msdfile):
    """returns the times and the MSD of each cluster instance (nan after the end of its life), from a <popul base>__<cluster>.msd.dat file"""
    nskip, columns = read_columns(msdfile)
    rows = []
    with open(msdfile,'r') as f:
        for iline, line in enumerate(f):
            if iline < nskip:
                continue
            fields = line.rstrip('\n').split('\t')                                  #empty fields once a cluster instance has stopped living
            rows.append([float(x) if x.strip() != '' else np.nan for x in fields])
    data = np.array(rows)
    return data[:,0], data[:,1:-1]                                                  #the last column is the mean, recomputed here


def analyze_cluster_file(msdfile,SkipTemps,EndTemps):
    """fits the MSD averaged over the cluster instances that live through the whole fit window
    the error is the spread of the D fitted on each of these instances, stdev(D_i)/sqrt(n), as they are independent samples"""
    t, msd = read_cluster_msd(msdfile)
    ninst = msd.shape[1]
    lasttime = np.array([t[np.where(~np.isnan(msd[:,ii]))[0][-1]] for ii in range(ninst)])       #last time of the MSD of each instance
    tmax = 0.75*np.median(lasttime) if EndTemps is None else EndTemps                              #by default most instances live through the window
    used = np.where(lasttime >= tmax)[0]
    window = (t >= SkipTemps) & (t <= tmax)
    result = {'ninst':ninst, 'nused':len(used), 'tmin':SkipTemps, 'tmax':tmax}
    if len(used) == 0 or np.sum(window) < 3:
        return result, None
    mean = np.nanmean(msd[:,used],axis=1)
    defined = ~np.isnan(mean)
    D, D_err, R2, slope, intercept = calculation_diffusivities(t[defined],mean[defined],SkipTemps,3,tmax)
    error = 'fit'
    if len(used) > 1:
        Di = [calculation_diffusivities(t[window],msd[window,ii],SkipTemps,3,tmax)[0] for ii in used]
        D_err = np.std(Di,ddof=1)/np.sqrt(len(used))
        error = 'instances'
    alpha = diffusion_exponent(t[defined],mean[defined],SkipTemps,tmax)
    result.update({'D':D, 'D_err':D_err, 'error':error, 'R2':R2, 'slope':slope, 'intercept':intercept, 'alpha':alpha,
                   't':t[defined], 'msd':mean[defined]})
    return result, True


def refit_cluster(fit,tmin,tmax):
    """the fit of a cluster file over another window (the instances used change with it) ; False when none is usable"""
    result, ok = analyze_cluster_file(fit['file'],tmin,tmax)
    if ok is None:
        return False
    fit.update(result)
    fit['title'] = fit['file']+' : '+fit['cluster']+' ('+str(result['nused'])+' instances)'
    return True


def main(argv):
    SkipTemps = 0
    EndTemps = None
    plot = False
    usage = 'analyze_msd_clusters.py -s <SkipTemps>(time in fs, default=0) -e <EndTemps>(time in fs, default=0.75 of the median lifetime) -v'
    try:
        options,arg = getopt.getopt(argv,"hs:e:v",["sSkipTemps","eEndTemps","visual"])
    except getopt.GetoptError:
        print(usage)
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('analyze_msd_clusters.py program to compute the diffusion coefficients of the clusters (chemical species) from the')
            print('<popul base>__<cluster>.msd.dat files of msd_clusters_umd.py, and write them in a diffusivities_clusters.txt file per subfolder')
            print(usage)
            print('The MSD of the center of mass is averaged over the instances of the cluster that live through the whole fit window,')
            print('then fitted linearly between SkipTemps and EndTemps : D = slope/6.')
            print('Default EndTemps : 0.75 of the median lifetime (last time of the MSD) of the instances, so that most of them are used.')
            print('-v : display each averaged MSD with its linear fit and the fit window (grey), one at a time (default: no display). Two clicks')
            print('     on the graph set the fit window of this cluster : the first one its start, the second one its end ; the instances living')
            print('     through it are averaged and fitted again at once. Automatic : back to the window of -s/-e. Previous/Next : the other')
            print('     clusters, one pass ; Next on the last one closes the window, and the files are then written with the windows set.')
            print('')
            print('diffusivities_clusters.txt has one line per cluster file, with:')
            print('     file, cluster, number of instances, number used in the fit, start and end of the fit (fs),')
            print('     D (m^2/s), its error (m^2/s), the origin of the error, R^2, slope, intercept, and alpha')
            print('     - error "instances" : stdev of the D fitted on each instance / sqrt(number of instances used)')
            print('       error "fit" (one instance only) : standard error of the slope, which underestimates the uncertainty')
            print('     - alpha, exponent of MSD ~ t^alpha over the fit window : 1 for normal diffusion, < 1 if the window is still sub-diffusive,')
            print('       > 1 if still ballistic. Clusters living only a few hundred fs rarely reach the diffusive regime : check alpha before using D.')
            sys.exit()
        elif opt in ("-s","--sSkipTemps"):
            SkipTemps = float(arg)
        elif opt in ("-e","--eEndTemps"):
            EndTemps = float(arg)
        elif opt in ('-v','--visual'):
            plot = True
    folders = []                                                                    #(folder, [(file, cluster, fit or the result without fit)])
    for dirpath, dirnames, filenames in os.walk(os.curdir):
        files = [file for file in sorted(glob.glob(dirpath+'/*.msd.dat')) if is_cluster_file(read_columns(file)[1])]
        if files == []:
            continue
        rows = []
        for file in files:
            cluster = os.path.basename(file)[:-len('.msd.dat')].split('__')[-1]
            result, ok = analyze_cluster_file(file,SkipTemps,EndTemps)
            if ok is None:
                print(file,': no instance lives through a fit window of at least 3 points, between ',result['tmin'],' and ',result['tmax'],' fs ; set -e')
                rows.append((file,cluster,result,False))
                continue
            print(file,': ',cluster,', fit between ',result['tmin'],' and ',result['tmax'],' fs on ',result['nused'],' of ',result['ninst'],' instances ; D = ','%.3e' % result['D'],' +/- ','%.2e' % result['D_err'],' m2/s, alpha = ','%.2f' % result['alpha'])
            warning = regime_warning(result['alpha'])
            if warning != '':
                print('    alpha',warning)
            rows.append((file,cluster,dict(result,file=file,cluster=cluster,default=(result['tmin'],result['tmax']),
                                           title=file+' : '+cluster+' ('+str(result['nused'])+' instances)'),True))
        folders.append((dirpath,rows))
    fits = [fit for dirpath, rows in folders for file, cluster, fit, ok in rows if ok]
    if plot and fits != []:
        plot_fits(fits,refit_cluster)
        for fit in fits:
            if (fit['tmin'],fit['tmax']) != fit['default']:
                print(fit['title'],': window set by hand, ','%.0f' % fit['tmin'],' - ','%.0f' % fit['tmax'],' fs : D = ','%.3e' % fit['D'],' +/- ',
                      '%.2e' % fit['D_err'],' m2/s, alpha = ','%.2f' % fit['alpha'])
    for dirpath, rows in folders:                                                   #written after the plots, with the windows set there
        newfilename = dirpath+'/diffusivities_clusters.txt'
        print('The file ',newfilename,' is created')
        with open(newfilename,'w') as f:
            f.write('\t'.join(['file','cluster','n_instances','n_used','fit_start(fs)','fit_end(fs)','D(m2/s)','D_stdev','error','R_squared','slope','intercept','alpha'])+'\n')
            for file, cluster, fit, ok in rows:
                row = [file, cluster, str(fit['ninst']), str(fit['nused']), str(fit['tmin']), str(fit['tmax'])]
                row += [str(fit[key]) for key in ('D','D_err','error','R2','slope','intercept','alpha')] if ok else ['nan']*7
                f.write('\t'.join(row)+'\n')


if __name__ == "__main__":
    main(sys.argv[1:])
