#!/usr/bin/env python3
###
#AUTHORS: ANAIS KOBSCH
###



"""     ********* Importation of the packages and modules used here *********     """
import sys
import getopt
import glob
import os
import numpy as np
from scipy import stats


def read_columns(msdfile):
    """returns the number of lines before the data and the list of column names (without the time column)"""
    with open(msdfile, 'r') as f:
        nskip = 0
        for line in f:
            nskip += 1
            if not line.startswith('#'):
                return nskip, line.split()[1:]
    return nskip, []


def is_cluster_file(columns):
    """msd_clusters_umd.py files have one column per cluster instance, cl0 cl1 ..., then their mean"""
    return len(columns) > 0 and columns[-1] == 'mean' and columns[0] == 'cl0'


def fit_end(t,tmax=None):
    """end of the fit window: tmax if given, 0.75 of the last time otherwise (the long-time MSD is noisy)"""
    return 0.75*t[-1] if tmax is None else tmax


def calculation_diffusivities(t,data,tmin,dimension=3,tmax=None):
    """     ********* Calculation of the self diffusivities (using linear regression between tmin and tmax) *********     """
    window = (t >= tmin) & (t <= fit_end(t,tmax))
    NewTemps = t[window]
    NewData = data[window]
    #least squared linear regression giving the correlation coefficient
    slope, intercept, r_value, p_value, std_err = stats.linregress(NewTemps, NewData) #intercept = y(x=0), std_err = error on the slope
    R_squared = r_value**2
    #MSD = 2 * dimension * D * t ; A^2/fs -> m^2/s is 10^-5
    diffusivity = slope * 10**(-5) / (2.*dimension)
    diff_err = std_err * 10**(-5) / (2.*dimension)
    return diffusivity , diff_err , R_squared , slope , intercept


def diffusion_exponent(t,data,tmin,tmax):
    """exponent alpha of MSD ~ t^alpha over the fit window (slope in log-log scale)
    alpha = 1 for normal diffusion ; alpha < 1 means the window is still in the sub-diffusive (caged) regime, alpha > 1 in the ballistic one"""
    window = (t >= tmin) & (t <= tmax) & (t > 0) & (data > 0)
    if np.sum(window) < 3:
        return float('nan')
    return stats.linregress(np.log(t[window]), np.log(data[window]))[0]


BALLISTIC_ALPHA = 1.2     #the ballistic regime (MSD ~ t^2) ends when the local exponent of the MSD falls below this


def ballistic_end(t,data,tmax=None):
    """end of the ballistic regime : the first time at which the local exponent d ln MSD / d ln t, averaged over a factor 1.3 in time,
    falls below BALLISTIC_ALPHA (2 for free flight, 1 for diffusion, < 1 in the cage of the neighbours) ; None if it does not, or if
    the window up to tmax would then have fewer than 3 points"""
    ok = (t > 0) & (data > 0)
    lt, lm = np.log(t[ok]), np.log(data[ok])
    if len(lt) < 3:
        return None
    a = np.gradient(lm, lt)
    cumul = np.concatenate(([0.0], np.cumsum(a)))
    lo, hi = np.searchsorted(lt, lt-0.13), np.searchsorted(lt, lt+0.13, side='right')
    below = np.where((cumul[hi]-cumul[lo])/(hi-lo) < BALLISTIC_ALPHA)[0]
    if len(below) == 0:
        return None
    tb = float(t[ok][below[0]])
    if tmax is not None and np.sum((t >= tb) & (t <= tmax)) < 3:
        return None
    return tb


def regime_warning(alpha):
    """warning when the fit window is not in the diffusive regime, MSD ~ t^alpha with alpha = 1"""
    if alpha < 0.9:
        return ' < 0.9 : the fit window is still sub-diffusive (caged) ; consider a larger end of the fit -e'
    if alpha > 1.1:
        return ' > 1.1 : the fit window is still ballistic ; consider a larger start -s (or -b in msd_umd.py) and end -e of the fit'
    return ''


def read_blocks(blockfile):
    """returns the times and a dictionary column name -> list of the MSDs of each block, from a <base>.msd.blocks.dat file"""
    nskip, names = read_columns(blockfile)
    alldata = np.loadtxt(blockfile,skiprows=nskip,ndmin=2)
    blocks = {}
    for icol, name in enumerate(names):
        blocks.setdefault(name.split('#')[0],[]).append(alldata[:,icol+1])    #columns are named El#1, El#2 ...
    return alldata[:,0], blocks


def block_diffusivity_error(t,blocks,tmin,tmax,dimension=3):
    """error on the D fitted between tmin and tmax on the full trajectory, from the spread of the D fitted on each block of it
    the block MSDs are shorter, so their fit window ends at most at 0.75 of their last time, tend
    the variance of a fitted D grows with the end of the fit window and decreases with the length of the trajectory, so that
    error = stdev(D_blocks) * sqrt(tmax / (nblocks * tend)) ; it is stdev(D_blocks)/sqrt(nblocks) when both windows are the same"""
    tend = min(fit_end(t), tmax)
    if len(blocks) < 2 or np.sum((t >= tmin) & (t <= tend)) < 3:
        return float('nan')
    D = [calculation_diffusivities(t,data,tmin,dimension,tend)[0] for data in blocks]
    return np.std(D,ddof=1)*np.sqrt(tmax/(len(D)*tend))


def fit_msd(fit,tmin,tmax):
    """fits the MSD of fit (keys t, msd, dimension, blocktimes, blocks) between tmin and tmax : sets tmin, tmax, D, D_err, R2, slope,
    intercept, alpha ; False when the window has fewer than 3 points"""
    t, data = fit['t'], fit['msd']
    if np.sum((t >= tmin) & (t <= tmax)) < 3:
        return False
    D, D_err, R2, slope, intercept = calculation_diffusivities(t,data,tmin,fit['dimension'],tmax)
    note = ''
    if fit['blocks'] is not None:
        D_err = block_diffusivity_error(fit['blocktimes'],fit['blocks'],tmin,tmax,fit['dimension'])
        if np.isnan(D_err):#the error of the linear fit would be much too small : none rather than a wrong one
            note = 'no error : the window must start before the end of the block fits, %.0f fs' % fit_end(fit['blocktimes'])
    fit.update(tmin=tmin, tmax=tmax, D=D, D_err=D_err, R2=R2, slope=slope, intercept=intercept, alpha=diffusion_exponent(t,data,tmin,tmax),
               note=note)
    return True


def plot_fits(fits,refit=None):
    """Plots one MSD at a time with its linear fit and the fit window, one pass with Previous/Next (Next on the last one closes the window).
    Each fit is a dictionary with the keys title, t, msd, tmin, tmax, slope, intercept, D, D_err, R2, alpha, and default (the window of -s/-e).
    With refit(fit,tmin,tmax) (True if fitted), two clicks set the window : the first one its start, the second one its end ; Automatic
    goes back to the default window. The fits are changed in place, the caller writes them after the window is closed."""
    import matplotlib.pyplot as plt                     #imported here, so that msd_umd.py, which imports this module, does not need it
    from matplotlib.widgets import Button
    fig, ax = plt.subplots(figsize=(9,6))
    fig.subplots_adjust(bottom=0.2)
    line, = ax.plot([], [], label='MSD')
    fitline, = ax.plot([], [], 'r--', label='linear fit')
    ax.set_xlabel('time (fs)')
    ax.set_ylabel('MSD (A$^2$)')
    ax.legend(loc='upper left')
    state = {'ifit':0, 'span':None, 'clicks':[], 'marks':[], 'message':''}
    def draw():
        fit = fits[state['ifit']]
        line.set_data(fit['t'], fit['msd'])
        tfit = fit['t'][(fit['t'] >= fit['tmin']) & (fit['t'] <= fit['tmax'])]      #the fit is drawn over its window only, so it does not set the scale
        fitline.set_data(tfit, fit['slope']*tfit + fit['intercept'])
        if state['span'] is not None:
            state['span'].remove()
        state['span'] = ax.axvspan(fit['tmin'], fit['tmax'], color='0.9', zorder=0)    #fit window
        for mark in state['marks']:
            mark.remove()
        state['marks'] = [ax.axvline(x, color='tab:green', linewidth=1) for x in state['clicks']]#start clicked, waiting for the end
        ax.set_title(fit['title'] + '  (' + str(state['ifit']+1) + '/' + str(len(fits)) + ')\nfit %.0f - %.0f fs%s :  D = ' % (fit['tmin'],fit['tmax'],
                     '  (set by hand)' if (fit['tmin'],fit['tmax']) != fit['default'] else '') + '%.3e' % fit['D'] + ' +/- ' + '%.2e' % fit['D_err']
                     + ' m$^2$/s,  R$^2$ = ' + '%.3f' % fit['R2'] + ',  alpha = ' + '%.2f' % fit['alpha']
                     + (('\n'+fit['note']) if fit.get('note') else ''), fontsize='medium')
        if refit is not None:
            todo = 'click the end of the fit window' if state['clicks'] else 'click the start of the fit window, then its end'
        else:
            todo = ''
        if state['ifit'] == len(fits)-1:
            todo += '    (last : Next closes'+(' and writes the files)' if refit is not None else ')')
        fig.suptitle(todo + (('    '+state['message']) if state['message'] else ''), x=0.98, ha='right', fontsize='small')
        ax.relim()
        ax.autoscale(enable=True)      #re-enabled in case the toolbar zoom switched it off
        fig.canvas.draw_idle()
    def step(delta):
        if state['ifit'] + delta >= len(fits):
            plt.close(fig)
            return
        state['ifit'] = max(0, state['ifit'] + delta)
        state['clicks'], state['message'] = [], ''
        draw()
    def onclick(event):
        if refit is None or event.inaxes is not ax or event.button != 1 or (fig.canvas.toolbar is not None and fig.canvas.toolbar.mode != ''):
            return
        state['message'] = ''
        if not state['clicks']:
            state['clicks'] = [event.xdata]
        else:
            tmin, tmax = sorted((state['clicks'][0], event.xdata))
            state['clicks'] = []
            if not refit(fits[state['ifit']], tmin, tmax):
                state['message'] = 'fewer than 3 points in the window : click both ends again'
        draw()
    def automatic(event):
        state['clicks'], state['message'] = [], ''
        if not refit(fits[state['ifit']], *fits[state['ifit']]['default']):
            state['message'] = 'fewer than 3 points in the default window'
        draw()
    buttons = []
    for label, pos, func in ((('Automatic', 0.43, automatic),) if refit is not None else ()) + (('Previous', 0.59, lambda event: step(-1)),
                                                                                             ('Next', 0.75, lambda event: step(1))):
        bt = Button(fig.add_axes([pos, 0.03, 0.15, 0.06]), label)
        bt.on_clicked(func)
        buttons.append(bt)
    fig._umd_buttons = buttons       #keeps the buttons alive, otherwise they are garbage-collected and stop responding
    fig.canvas.mpl_connect('button_press_event', onclick)
    draw()
    plt.show()


def headerfile(atoms, dirpath, headfile):
    """creation of the newfile with correct header"""
    firstline = ['atom','fit','fit','fit']
    secondline = ['file','fit_start(fs)','fit_end(fs)','error']
    for atom in atoms:
        for i in range(0,8):
            firstline.append(atom)
        secondline.extend(['D(m2/s)','D_stdev','R_squared','slope','intercept','alpha','fit_start(fs)','fit_end(fs)'])
    newfilename = dirpath+'/diffusivities.txt'
    print('The file ',newfilename,' is created')
    f = open(newfilename,'w')
    if headfile != '':
        with open(headfile,'r') as hf:
            for line in hf:
                f.write(line)
    f.write("\t".join(x for x in firstline)+ "\n")
    f.write("\t".join(x for x in secondline)+ "\n")
    return f



def main(argv):
    """     ********* Main program *********     """
    SkipTemps = None                                                                #None : the end of the ballistic regime of each column
    EndTemps = None
    headfile = ''
    plot = False
    try:
        options,arg = getopt.getopt(argv,"hs:e:f:v",["sSkipTemps","eEndTemps","file","visual"])
    except getopt.GetoptError:
        print('analyze_msd.py -s <SkipTemps>(time in fs, default : end of the ballistic regime) -e <EndTemps>(time in fs, default see -h) -f <header file (default none)> -v')
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('analyze_msd.py program to compute self diffusivities from the msd.dat files of msd_umd.py and write them in a diffusivities.txt file per subfolder')
            print('analyze_msd.py -s <SkipTemps>(time in fs, default : end of the ballistic regime) -e <EndTemps>(time in fs, default see below) -f <header file (default none)> -v')
            print('The MSD is fitted linearly between SkipTemps (end of the ballistic regime) and EndTemps.')
            print('Default SkipTemps : the end of the ballistic regime of each column, the first time at which the local exponent')
            print('d ln MSD / d ln t (2 in free flight, 1 in diffusion), averaged over a factor 1.3 in time, falls below '+str(BALLISTIC_ALPHA)+'. After it the')
            print('MSD is often still sub-diffusive (alpha < 1 : the atoms rattle in the cage of their neighbours) : check alpha, with -v.')
            print('Default EndTemps : 0.75 of the last time of the blocks in <base>.msd.blocks.dat if it exists, otherwise 0.15 of the last')
            print('time of the MSD (the same for the default 5 blocks of msd_umd.py). Longer fits, up to lags with few independent time')
            print('origins, give much noisier D ; increase EndTemps for slow, sub-diffusive systems.')
            print('-v : display each MSD with its linear fit and the fit window (grey), one at a time (default: no display). Two clicks on the')
            print('     graph set the fit window of this MSD : the first one its start, the second one its end ; D, its error, R^2 and alpha')
            print('     are fitted again at once. Automatic : back to the default window. Previous/Next : the other MSDs, one pass ; Next on')
            print('     the last one closes the window (so does closing it), and the diffusivities.txt files are then written with the')
            print('     windows set. Clicks are ignored while the zoom or pan of the toolbar is active.')
            print('Cluster MSD files from msd_clusters_umd.py are skipped : use analyze_msd_clusters.py for them.')
            print('')
            print('headerfile is a txt file containing the element and number (useful to compute the densities)')
            print('')
            print('This code produces a diffusivities.txt file with:')
            print('     - a column with filename, then the start (\'ballistic\' : the end of the ballistic regime of each column) and end of the fit')
            print('       of -s/-e (fs), and the origin of the error (blocks or fit)')
            print('     - 8 columns per atom type (all the columns found in the files of the folder ; nan if a file does not have it)')
            print('         - Diffusivity (m^2/s), D = slope/6 for the total MSD and slope/2 for the MSD along one axis (columns El_x, El_a, ...)')
            print('         - stdev on diffusivity (m^2/s): if <base>.msd.blocks.dat (written by msd_umd.py) exists, the standard error')
            print('           from the D fitted on each block of the trajectory, stdev(D_blocks)*sqrt(EndTemps/(nblocks*end of the block fit)),')
            print('           i.e. stdev(D_blocks)/sqrt(nblocks) when the block and full fits share the same window ; otherwise the standard')
            print('           error of the slope of the fit, which underestimates the uncertainty as the MSD points are correlated')
            print('         - R^2  coefficient associated to the fit')
            print('         - slope of the fit')
            print('         - intercept (y(x=0))')
            print('         - alpha, exponent of MSD ~ t^alpha over the fit window : 1 for normal diffusion, < 1 if the window is still sub-diffusive,')
            print('           > 1 if still ballistic ; a warning is printed when alpha is outside 0.9-1.1')
            print('         - start and end of the fit window of this column (fs) : those of -s/-e (by default the end of its ballistic regime),')
            print('           or those set by clicks with -v')
            sys.exit()
        elif opt in ("-s","--sSkipTemps"):
            SkipTemps = float(arg)
        elif opt in ("-e","--eEndTemps"):
            EndTemps = float(arg)
        elif opt in ('-f','--file'):
            headfile = str(arg)
        elif opt in ('-v','--visual'):
            plot = True
    folders = []                                                                    #(folder, its columns, its files : (file, window, error, fits))
    for dirpath, dirnames, filenames in os.walk(os.curdir):
        files = []
        atoms = []                                                                  #all the columns found in the files of the folder, in order of appearance
        for file in sorted(glob.glob(dirpath+'/*.msd.dat')):                       #I list every msd files in alphabetic order
            nskip, columns = read_columns(file)
            if is_cluster_file(columns):
                print(file,': cluster MSD, skipped (use analyze_msd_clusters.py)')
                continue
            files.append(file)
            atoms += [atom for atom in columns if atom not in atoms]
        if files == []:
            continue
        rows = []
        for file in files:
            nskip, columns = read_columns(file)
            alldata = np.loadtxt(file,skiprows=nskip,ndmin=2)
            Temps = alldata[:,0]
            blockfile = file[:-len('.msd.dat')]+'.msd.blocks.dat'
            if os.path.isfile(blockfile):
                BlockTemps, blocks = read_blocks(blockfile)
                tmax = fit_end(BlockTemps,EndTemps)                         #by default 0.75 of the last time of a block
                error = 'blocks'
            else:
                BlockTemps, blocks = None, {}
                tmax = 0.15*Temps[-1] if EndTemps is None else EndTemps   #by default as for 5 blocks : 0.75 * (last time / 5)
                error = 'fit'
            print(file,': fit from ',('the end of the ballistic regime of each column' if SkipTemps is None else SkipTemps),' to ',tmax,
                  ' fs ; errors on D from the '+('blocks in '+blockfile if error=='blocks' else 'linear fit'))
            fits = {}
            for atom in columns:
                data = alldata[:,columns.index(atom)+1]
                tmin = SkipTemps
                if tmin is None:#the end of the ballistic regime of the total MSD of the element, also for its axes (El_x ...)
                    total = atom.split('_')[0]
                    tmin = ballistic_end(Temps,alldata[:,columns.index(total)+1] if total in columns else data,tmax)
                    if tmin is None:
                        print('   ',atom,': the ballistic regime (local alpha > ',BALLISTIC_ALPHA,') lasts beyond the fit window : fitted from 0')
                        tmin = 0.0
                fit = {'title':file+' : '+atom, 't':Temps, 'msd':data, 'default':(tmin,tmax),
                       'dimension':1 if '_' in atom else 3,                 #El_x, El_a, El_1 ... are MSDs projected along one axis
                       'blocktimes':BlockTemps, 'blocks':blocks.get(atom)}
                if not fit_msd(fit,tmin,tmax):
                    print('   ',atom,': fewer than 3 points between ',tmin,' and ',tmax,' fs')
                    continue
                warning = regime_warning(fit['alpha'])
                if warning != '':
                    print('   ',atom,': alpha = ','%.2f' % fit['alpha'],warning)
                fits[atom] = fit
            rows.append((file,(SkipTemps,tmax),error,fits))
        folders.append((dirpath,atoms,rows))
    allfits = [fit for dirpath, atoms, rows in folders for file, window, error, fits in rows for fit in fits.values()]
    if plot and allfits != []:
        plot_fits(allfits,fit_msd)
        for fit in allfits:
            if (fit['tmin'],fit['tmax']) != fit['default']:
                print(fit['title'],': window set by hand, ','%.0f' % fit['tmin'],' - ','%.0f' % fit['tmax'],' fs : D = ','%.3e' % fit['D'],' +/- ',
                      '%.2e' % fit['D_err'],' m2/s, R2 = ','%.3f' % fit['R2'],', alpha = ','%.2f' % fit['alpha'],(' ; '+fit['note']) if fit['note'] else '')
    for dirpath, atoms, rows in folders:                                            #written after the plots, with the windows set there
        f = headerfile(atoms, dirpath, headfile)
        for file, window, error, fits in rows:
            results = [file, 'ballistic' if window[0] is None else str(window[0]), str(window[1]), error]
            for atom in atoms:
                if atom not in fits:
                    results.extend(['nan']*8)
                    continue
                fit = fits[atom]
                results.extend([str(fit[key]) for key in ('D','D_err','R2','slope','intercept','alpha','tmin','tmax')])
            f.write("\t".join(x for x in results)+ "\n")                  #we write in the file the result line
        f.close()


#     ********* Execution *********
if __name__ == "__main__":
    main(sys.argv[1:])
