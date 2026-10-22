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


def headerfile(firstfile, dirpath, headfile):
    """creation of the newfile with correct header"""
    firstline = ['atom']  #beginning of the first line of the file gofr
    secondline = ['file']
    # creation of the header from the column-header line of the first file (the first line not starting with #)
    nskip, atoms = read_columns(firstfile)  #list of atoms
    for atom in atoms:
        for i in range(0,5):
            firstline.append(atom)
        secondline.extend(['D(m2/s)','D_stdev','R_squared','slope','intercept'])
    newfilename = dirpath+'/diffusivities.txt' 
    print('The file ',newfilename,' is created')
    f = open(newfilename,'w')
    if headfile != '':
        with open(headfile,'r') as hf:
            for line in hf:
                f.write(line)
    f.write("\t".join(x for x in firstline)+ "\n")
    f.write("\t".join(x for x in secondline)+ "\n")
    return f, atoms      #I return the newly created files f along with the list of element couples


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



def main(argv):
    """     ********* Main program *********     """
    SkipTemps = 0
    EndTemps = None
    headfile = ''
    try:
        options,arg = getopt.getopt(argv,"hs:e:f:",["sSkipTemps","eEndTemps","file"])
    except getopt.GetoptError:
        print('analyze_msd.py -s <SkipTemps>(time in fs, default=0) -e <EndTemps>(time in fs, default see -h) -f <header file (default none)>')
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('analyze_msd.py program to compute self diffusivities from msd.dat file and write them in a .txt file per subfolder')
            print('analyze_msd.py -s <SkipTemps>(time in fs, default=0) -e <EndTemps>(time in fs, default see -h) -f <header file (default none)>')
            print('The MSD is fitted linearly between SkipTemps (end of the ballistic regime) and EndTemps.')
            print('Default EndTemps : 0.75 of the last time of the blocks in <base>.msd.blocks.dat if it exists, otherwise 0.15 of the last')
            print('time of the MSD (the same for the default 5 blocks of msd_umd.py). Longer fits, up to lags with few independent time')
            print('origins, give much noisier D ; increase EndTemps for slow, sub-diffusive systems.')
            print('')
            print('headerfile is a txt file containing the element and number (useful to compute the densities)')
            print('')
            print('This code produces a diffusivities.txt file with:')
            print('     - a column with filename')
            print('     - 5 columns per atom type')
            print('         - Diffusivity (m^2/s), D = slope/6 for the total MSD and slope/2 for the MSD along one axis (columns El_x, El_a, ...)')
            print('         - stdev on diffusivity (m^2/s): if <base>.msd.blocks.dat (written by msd_umd.py) exists, the standard error')
            print('           from the D fitted on each block of the trajectory, stdev(D_blocks)*sqrt(EndTemps/(nblocks*end of the block fit)),')
            print('           i.e. stdev(D_blocks)/sqrt(nblocks) when the block and full fits share the same window ; otherwise the standard')
            print('           error of the slope of the fit, which underestimates the uncertainty as the MSD points are correlated')
            print('         - R^2  coefficient associated to the fit')     
            print('         - slope of the fit') 
            print('         - intercept (y(x=0))') 
            sys.exit()
        elif opt in ("-s","--sSkipTemps"):
            SkipTemps = float(arg)
        elif opt in ("-e","--eEndTemps"):
            EndTemps = float(arg)
        elif opt in ('-f','--file'):
            headfile = str(arg)
    for dirpath, dirnames, filenames in os.walk(os.curdir):
        files = sorted(glob.glob(dirpath+'/*.msd.dat')) #I list every gofr files in alphabetic order
        if files != []:
            f, atoms = headerfile(files[0], dirpath, headfile)                          #I create the first newfile for gofr and save the list of element couples 
            for file in files:
                results = [file]
                nskip, columns = read_columns(file)
                alldata = np.loadtxt(file,skiprows=nskip,ndmin=2)
                Temps = alldata[:,0]
                blockfile = file[:-len('.msd.dat')]+'.msd.blocks.dat'
                if os.path.isfile(blockfile):
                    BlockTemps, blocks = read_blocks(blockfile)
                    tmax = fit_end(BlockTemps,EndTemps)                         #by default 0.75 of the last time of a block
                    print(file,': errors on D from the blocks in ',blockfile,' ; fit between ',SkipTemps,' and ',tmax,' fs')
                else:
                    blocks = {}
                    tmax = 0.15*Temps[-1] if EndTemps is None else EndTemps   #by default as for 5 blocks : 0.75 * (last time / 5)
                    print(file,': no blocks file, errors on D from the linear fit ; fit between ',SkipTemps,' and ',tmax,' fs')
                for atom in atoms:
                    data = alldata[:,columns.index(atom)+1]
                    dimension = 1 if '_' in atom else 3  #El_x, El_a, El_1 ... are MSDs projected along one axis
                    D, D_stdev, R, slope, intercept = calculation_diffusivities(Temps,data,SkipTemps,dimension,tmax) #we compute the diffusivities and stdev using linear fit
                    if atom in blocks:
                        D_stdev = block_diffusivity_error(BlockTemps,blocks[atom],SkipTemps,tmax,dimension)
                    results.extend([str(D), str(D_stdev), str(R), str(slope), str(intercept)])
                f.write("\t".join(x for x in results)+ "\n")                  #we write in the file the result line
            f.close()


#     ********* Execution *********
if __name__ == "__main__":
    main(sys.argv[1:])





