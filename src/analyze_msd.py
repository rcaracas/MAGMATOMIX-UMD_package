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


def plot_fits(fits):
    """Plots one MSD at a time with its linear fit and the fit window, with Previous/Next buttons to move between them.
    Each fit is a dictionary with the keys title, t, msd, tmin, tmax, slope, intercept, D, D_err, R2, alpha"""
    import matplotlib.pyplot as plt                     #imported here, so that msd_umd.py, which imports this module, does not need it
    from matplotlib.widgets import Button
    fig, ax = plt.subplots()
    line, = ax.plot([], [], label='MSD')
    fitline, = ax.plot([], [], 'r--', label='linear fit')
    ax.set_xlabel('time (fs)')
    ax.set_ylabel('MSD (A$^2$)')
    ax.legend(loc='upper left')
    state = {'ifit':0, 'span':None}
    def draw():
        fit = fits[state['ifit']]
        line.set_data(fit['t'], fit['msd'])
        tfit = fit['t'][(fit['t'] >= fit['tmin']) & (fit['t'] <= fit['tmax'])]      #the fit is drawn over its window only, so it does not set the scale
        fitline.set_data(tfit, fit['slope']*tfit + fit['intercept'])
        if state['span'] is not None:
            state['span'].remove()
        state['span'] = ax.axvspan(fit['tmin'], fit['tmax'], color='0.9', zorder=0)    #fit window
        title = fit['title']
        if len(fits) > 1:
            title = title + '  (' + str(state['ifit']+1) + '/' + str(len(fits)) + ')'
        ax.set_title(title + '\nD = ' + '%.3e' % fit['D'] + ' +/- ' + '%.2e' % fit['D_err'] + ' m$^2$/s,  R$^2$ = ' + '%.3f' % fit['R2'] + ',  alpha = ' + '%.2f' % fit['alpha'], fontsize='medium')
        ax.relim()
        ax.autoscale(enable=True)      #re-enabled in case the toolbar zoom switched it off
        fig.canvas.draw_idle()
    draw()
    if len(fits) > 1:
        fig.subplots_adjust(bottom=0.2)
        def step(delta):
            state['ifit'] = (state['ifit'] + delta) % len(fits)
            draw()
        bprev = Button(fig.add_axes([0.59, 0.03, 0.15, 0.06]), 'Previous')
        bnext = Button(fig.add_axes([0.75, 0.03, 0.15, 0.06]), 'Next')
        bprev.on_clicked(lambda event: step(-1))
        bnext.on_clicked(lambda event: step(1))
        fig._umd_buttons = (bprev, bnext)       #keeps the buttons alive, otherwise they are garbage-collected and stop responding
    plt.show()


def headerfile(atoms, dirpath, headfile):
    """creation of the newfile with correct header"""
    firstline = ['atom','fit','fit','fit']
    secondline = ['file','fit_start(fs)','fit_end(fs)','error']
    for atom in atoms:
        for i in range(0,6):
            firstline.append(atom)
        secondline.extend(['D(m2/s)','D_stdev','R_squared','slope','intercept','alpha'])
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
    SkipTemps = 0
    EndTemps = None
    headfile = ''
    plot = False
    try:
        options,arg = getopt.getopt(argv,"hs:e:f:v",["sSkipTemps","eEndTemps","file","visual"])
    except getopt.GetoptError:
        print('analyze_msd.py -s <SkipTemps>(time in fs, default=0) -e <EndTemps>(time in fs, default see -h) -f <header file (default none)> -v')
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('analyze_msd.py program to compute self diffusivities from the msd.dat files of msd_umd.py and write them in a diffusivities.txt file per subfolder')
            print('analyze_msd.py -s <SkipTemps>(time in fs, default=0) -e <EndTemps>(time in fs, default see below) -f <header file (default none)> -v')
            print('The MSD is fitted linearly between SkipTemps (end of the ballistic regime) and EndTemps.')
            print('Default EndTemps : 0.75 of the last time of the blocks in <base>.msd.blocks.dat if it exists, otherwise 0.15 of the last')
            print('time of the MSD (the same for the default 5 blocks of msd_umd.py). Longer fits, up to lags with few independent time')
            print('origins, give much noisier D ; increase EndTemps for slow, sub-diffusive systems.')
            print('-v : display each MSD with its linear fit and the fit window, with Previous/Next buttons (default: no display)')
            print('Cluster MSD files from msd_clusters_umd.py are skipped : use analyze_msd_clusters.py for them.')
            print('')
            print('headerfile is a txt file containing the element and number (useful to compute the densities)')
            print('')
            print('This code produces a diffusivities.txt file with:')
            print('     - a column with filename, then the start and end of the fit (fs) and the origin of the error (blocks or fit)')
            print('     - 6 columns per atom type (all the columns found in the files of the folder ; nan if a file does not have it)')
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
            sys.exit()
        elif opt in ("-s","--sSkipTemps"):
            SkipTemps = float(arg)
        elif opt in ("-e","--eEndTemps"):
            EndTemps = float(arg)
        elif opt in ('-f','--file'):
            headfile = str(arg)
        elif opt in ('-v','--visual'):
            plot = True
    fits = []
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
        if files != []:
            f = headerfile(atoms, dirpath, headfile)
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
                    blocks = {}
                    tmax = 0.15*Temps[-1] if EndTemps is None else EndTemps   #by default as for 5 blocks : 0.75 * (last time / 5)
                    error = 'fit'
                print(file,': fit between ',SkipTemps,' and ',tmax,' fs ; errors on D from the '+('blocks in '+blockfile if error=='blocks' else 'linear fit'))
                results = [file, str(SkipTemps), str(tmax), error]
                for atom in atoms:
                    if atom not in columns:
                        results.extend(['nan']*6)
                        continue
                    data = alldata[:,columns.index(atom)+1]
                    dimension = 1 if '_' in atom else 3  #El_x, El_a, El_1 ... are MSDs projected along one axis
                    D, D_stdev, R, slope, intercept = calculation_diffusivities(Temps,data,SkipTemps,dimension,tmax) #we compute the diffusivities and stdev using linear fit
                    if atom in blocks:
                        D_stdev = block_diffusivity_error(BlockTemps,blocks[atom],SkipTemps,tmax,dimension)
                    alpha = diffusion_exponent(Temps,data,SkipTemps,tmax)
                    results.extend([str(D), str(D_stdev), str(R), str(slope), str(intercept), str(alpha)])
                    warning = regime_warning(alpha)
                    if warning != '':
                        print('   ',atom,': alpha = ','%.2f' % alpha,warning)
                    if plot:
                        fits.append({'title':file+' : '+atom, 't':Temps, 'msd':data, 'tmin':SkipTemps, 'tmax':tmax, 'slope':slope, 'intercept':intercept,
                                     'D':D, 'D_err':D_stdev, 'R2':R, 'alpha':alpha})
                f.write("\t".join(x for x in results)+ "\n")                  #we write in the file the result line
            f.close()
    if plot and fits != []:
        plot_fits(fits)


#     ********* Execution *********
if __name__ == "__main__":
    main(sys.argv[1:])
