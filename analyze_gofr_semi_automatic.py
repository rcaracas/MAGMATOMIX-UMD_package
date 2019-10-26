#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 2017 & Thu Feb 6 2018

@author: Anaïs & Natalia
Langage: Python3

                ****     Analyze the plot of gofr and int gofr       ****
                          and save all the data into a txt file
                          and save the first minimum for all bonds in a .bonds file

This code requires all the gofr.dat files created by the script gofrs_umd.py

WARNING! for efficiency, your data needs to be located in different folders for each T 
and you should launch this script from the folder containing every subfolder temperature'

This code produces a gofrs.txt file per temperature with:
    - a column with filename
    - 4 column per A-B couple of atoms
        - x value of max(gofr)
        - y value associated
        - x value of the first min(gofr)
        - y value of int(gofr) associated to previous x value

And a .bonds file with three colums:
    - atom 1
    - atom 2
    - bondlength

"""

#*********** Importation of the packages and modules used here ************
import glob
import os
import sys
import getopt
import re
import numpy as np
import csv
import matplotlib.pyplot as plt
from matplotlib.widgets import Button



def headerfile(firstfile, dirpath):
    """creation of the newfile with correct header"""
    firstline = ['couple']  #beginning of the first line of the file gofr
    secondline = ['file']
    # creation of the header from the first line of the first file
    with open(firstfile, 'r') as f:
        header = f.readline()
    header = header.strip('dist')
    header = re.sub('(Int\([A-Za-z-]*\))', ' ',header).split()        #I substitute all Int(...)  by a blankspace and split using spaces
    for couple in header:
        for i in range(0,4):
            firstline.append(couple)
        secondline.extend(['xmax','ymax','xmin','coord'])
    newfilename = dirpath+'_gofrs.txt' 
    print('The file ',newfilename,' is created')
    f = open(newfilename,'w')
    f.write("\t".join(x for x in firstline)+ "\n")
    f.write("\t".join(x for x in secondline)+ "\n")
    return f, header      #I return the newly created files f along with the list of element couples

def atoms_columns(firstfile):
    """creation of the newfile for bonds"""
    atom1 = [] #creation of the first array for the file .bonds
    atom2 = [] #creation of the second array for the file .bonds
    with open(firstfile, 'r') as f:
        header = f.readline()
    header = header.strip('dist')
    header = re.sub('(Int\([A-Za-z-]*\))', ' ',header).split()        #I substitute all Int(...)  by a blankspace and split using spaces
    for couple in header:
        atom1.append(couple.split('-')[0])
        atom2.append(couple.split('-')[1])    #at the end of this loop we have created the two first column of our .bonds file
    newfilenamebond = firstfile.split('.gofr.dat')[0] + '.bonds.inp'
    print('The file ', newfilenamebond, 'is also created')
    b = open(newfilenamebond, 'w')
    return b, atom1, atom2      #I return the newly created file b along with the two columns for .bonds file


def interactive_fit(data,file,couples, guess_xmax, guess_xmin, distance, col, bonds, couple):
    """fit max and min of data using interactive plot for one couple of atoms"""
    print(couple, couples.index(couple) + col + 1, couples.index(couple) + col + 2)
    gofr, intgofr = np.loadtxt(file,
                               usecols=(couples.index(couple) + col + 1, couples.index(couple) + col + 2),
                               skiprows=1, unpack=True)
    col = col + 1
    
    orig_size = np.size(distance)  # size of our data
    nonzeros = np.count_nonzero(gofr)
    max_index = np.argmax(gofr)  # find the index of max gofr
    last_index = -1
    cutdistance = distance[max_index:last_index]  # cut off beginning for finding minimum
    cutgofr = gofr[max_index:last_index]  # cut off beginning for finding minimum
    cut_size = np.size(cutdistance)  # size of data
    cut_min_index = np.argmin(cutgofr)  # find the min of gofr in the cut segment
    min_index = orig_size - cut_size + cut_min_index - 1  # min index of the *uncut* gofr
    end_imin = int(min_index + min_index * 0.08)  # 8% to the right
    end_imax = int(max_index + max_index * 0.12)  # 12% to the left

    if nonzeros > 10  and end_imin < (orig_size - 1): #and end_imax < (orig_size - 1):
        plt.close()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.subplots_adjust(bottom=0.2)
        ax.set_title(file.split('.gofr.dat')[0] + '  ' + couple)
        plt.xlabel('Distance ($\AA$)')
        plt.ylabel('g(r)')
        plt.grid()
        ax.plot(distance, gofr, 'o', markersize=2)  # plot the fitted max

        def onclick(event):
            global xmax_click
            xmax_click = event.xdata  # record x value of click
            plt.close()
        fig.canvas.mpl_connect('button_press_event', onclick)
        plt.show()
        max_index = (np.abs(distance - xmax_click)).argmin()  # index of the closest value to clicked value

        # determine start and end cols for fitting the polynomial:
        end_imax = int(max_index + max_index * 0.12)  # 12% to the left (arbitrary)
        start_imax = int(max_index - max_index * 0.08)  # 8% to the right (arbitrary)
        if end_imax > (orig_size - 1):  # for rare cases where the max distance is out of bounds
            end_imax = orig_size - 1  # subtract one because cols start at 0

        # define the regions that need to be fitted for the max area:
        fit_max_dist = distance[int(start_imax):int(end_imax + 1)]  # ranges do not include last element, so add one
        fit_max_gofr = gofr[int(start_imax):int(end_imax + 1)]

        order = 3  # order of the polynomials
        max_fit = np.poly1d(np.polyfit(fit_max_dist, fit_max_gofr, order))
        xaxis_max = np.linspace(distance[int(start_imax)], distance[int(end_imax)],
                                1000)  # create new x values for the fitted max region
        yaxis_max = max_fit(xaxis_max)  # create new y values for fitted max region
        poly_max_i = np.argmax(yaxis_max)  # determine index of max y value

        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.subplots_adjust(bottom=0.2)
        ax.set_title(file.split('.gofr.dat')[0] + '  ' + couple)
        ax.plot(distance, gofr, 'o', xaxis_max, max_fit(xaxis_max), '-', markersize=2)  # plot the fitted max
        plt.xlabel('Distance ($\AA$)')
        plt.ylabel('g(r)')
        plt.grid()

        def onclick(event):
            global xmin_click
            xmin_click = event.xdata  # record x value of click
            plt.close()

        fig.canvas.mpl_connect('button_press_event', onclick)
        plt.show()
        min_index = (np.abs(distance - xmin_click)).argmin()  # index of the closest value to clicked value

        # determine start and end cols for fitting the polynomial for the min region:
        start_imin = int(min_index - min_index * 0.08)  # 8% to the left (arbitrary)
        end_imin = int(min_index + min_index * 0.08)  # 8% to the right (arbitrary)
        if end_imin > (orig_size - 1):  # for rare cases where the min distance is out of bounds
            end_imin = orig_size - 1  # subtract one because cols start at 0
        # define the regions that need to be fitted for min area using the clicked min:
        fit_min_dist = distance[int(start_imin):int(end_imin + 1)]  # ranges do not include last element, so add one
        fit_min_gofr = gofr[int(start_imin):int(end_imin + 1)]
        intgofr = intgofr[int(start_imin):int(end_imin + 1)]

        min_fit = np.poly1d(np.polyfit(fit_min_dist, fit_min_gofr, order))
        intgofr_fit = np.poly1d(np.polyfit(fit_min_dist, intgofr, order))
        xaxis_min = np.linspace(distance[int(start_imin)], distance[int(end_imin)],
                                1000)  # create new x values for the fitted min region
        yaxis_min = min_fit(xaxis_min)  # create new y values for fitted min region
        poly_min_i = np.argmin(yaxis_min)  # determine index of min y value
        yaxis_intgofr = intgofr_fit(xaxis_min)

        xmax_gofr = xaxis_max[int(poly_max_i)]  # x value of the max (x,y) of the max region of gofr
        ymax_gofr = yaxis_max[int(poly_max_i)]  # max y value of the max region of gofr
        xmin_gofr = xaxis_min[int(poly_min_i)]  # x value of the min (x,y) of the min region of gofr
        y_intgofr = yaxis_intgofr[
            int(poly_min_i)]  # y value of the integral of gofr that corresponds to the min x value of gofr

        class Decision(object):
            def good(self, event):
                global value_array
                plt.close(fig)  # close and move to next plot
                print('xmax=', round(xmax_gofr, 2), ' ymax=', round(ymax_gofr, 2), ' xmin=', round(xmin_gofr, 2), ' coord=', round(y_intgofr, 2))
                value_array = (xmax_gofr, ymax_gofr, xmin_gofr, y_intgofr)
                return value_array

            def bad(self, event):
                global value_array
                plt.close(fig)  # close and move to next plot
                print('bad')
                xmax_gofr = 0
                ymax_gofr = 0
                xmin_gofr = 0
                y_intgofr = 0
                value_array = (xmax_gofr, ymax_gofr, xmin_gofr, y_intgofr)
                return value_array

        fig = plt.figure()  # show figure with the newly fitted minimum curve
        ax = fig.add_subplot(111)
        plt.subplots_adjust(bottom=0.2)
        ax.set_title(file.split('.gofr.dat')[0] + '  ' + couple)
        ax.plot(distance, gofr, 'o', xaxis_max, max_fit(xaxis_max), '-', xaxis_min, min_fit(xaxis_min), 'r-',
                markersize=2)
        plt.xlabel('Distance ($\AA$)')
        plt.ylabel('g(r)')

        callback = Decision()
        axbad = plt.axes([0.12, 0.05, 0.1, 0.075])
        axgood = plt.axes([0.8, 0.05, 0.1, 0.075])
        bgood = Button(axgood, 'Good')
        bgood.on_clicked(callback.good)
        bbad = Button(axbad, 'Bad')
        bbad.on_clicked(callback.bad)
        plt.show(fig)

        xmax_gofr = value_array[0]
        ymax_gofr = value_array[1]
        xmin_gofr = value_array[2]
        y_intgofr = value_array[3]

    else: # simply output 0's if the data is no good
        print('skip')
        xmax_gofr = 0
        ymax_gofr = 0
        xmin_gofr = 0
        y_intgofr = 0
    
    data.append(str(xmax_gofr))
    data.append(str(ymax_gofr))
    bonds.append(str(xmin_gofr))  # for each couple of atom we add the corresponding xmin value
    data.append(str(xmin_gofr))
    data.append(str(y_intgofr))
    guess_xmax[couple] = xmax_gofr                  #update the initial guesses for the next acell
    guess_xmin[couple] = xmin_gofr
    
    cutdistance = 0  # re initialization of variables for the next loop
    fit_min_dist = 0
    fit_max_dist = 0
    fit_min_gofr = 0
    fit_max_gofr = 0
    gofr = 0
    intgofr = 0
    return (data, bonds, guess_xmin, guess_xmax, col)



xmin_click = []
xmax_click = []
def analyze_gofrs_interactive(data,file,couples, guess_xmax, guess_xmin, atoms):
    """ extraction of the min, max and coordination number using fits and interactive plot """
    distance = np.loadtxt(file, usecols=(0,), skiprows = 1, unpack=True)
    col = 0
    size = len(data)
    bonds = []
    for couple in couples:
        for atom_couple in atoms:
            if (couple == atom_couple):
                data, bonds, guess_xmin, guess_xmax, col = interactive_fit(data,file,couples, guess_xmax, guess_xmin, distance, col, bonds, couple)
        if len(data) == size: #if data hasn't changed then we put X in the list
            data.extend(['X','X','X','X'])
            bonds.append('0')
            col += 1
        size = len(data) #we update the size of data
    return(data,bonds, guess_xmax, guess_xmin)






def analyze_gofrs_automatic(data,file,couples, guess_xmax, guess_xmin, atoms):
    """ extraction of the min, max and coordination number using fits and previous values as initial guesses """
    distance = np.loadtxt(file, usecols=(0,), skiprows = 1, unpack=True)
    col = 0
    size = len(data)
    bonds = []
    for couple in couples:
        for atom_couple in atoms:
            if (couple == atom_couple):
                if (guess_xmax[couple] == 0 or guess_xmin[couple] == 0):  #if we haven't succeeded to have initial guesses, then we try again with interactive plot
                    data, bonds, guess_xmin, guess_xmax, col = interactive_fit(data,file,couples, guess_xmax, guess_xmin, distance, col, bonds, couple)
                else: #if we have correct initial guesses for this couple of atoms, then we use an automatic fitting process
                    print(couple,couples.index(couple)+col+1,couples.index(couple)+col+2)
                    gofr, intgofr = np.loadtxt(file, usecols = (couples.index(couple)+col+1,couples.index(couple)+col+2), skiprows = 1, unpack = True)
                    col = col + 1
                
                    orig_size = np.size(distance)
                    max_index = (np.abs(distance - guess_xmax[couple])).argmin()  # index of the closest value to guessed value
                    min_index = (np.abs(distance - guess_xmin[couple])).argmin()  # index of the closest value to guessed value
        
                    # determine start and end cols for fitting the polynomial:
                    start_imin = int(min_index - min_index * 0.08) # 8% to the left (arbitrary)
                    end_imin = int(min_index + min_index * 0.08)  # 8% to the right (arbitrary)
                    if end_imin > (orig_size - 1):  # for rare cases where the max distance is out of bounds
                        end_imin = orig_size - 1  # subtract one because cols start at 0
                    end_imax = int(max_index + max_index * 0.12) # 12% to the left (arbitrary)
                    start_imax = int(max_index - max_index * 0.08) # 8% to the right (arbitrary)
                    if end_imax > (orig_size - 1):  # for rare cases where the max distance is out of bounds
                        end_imax = orig_size - 1 # subtract one because cols start at 0
        
                    # define the regions that need to be fitted:
                    # ranges do not include last element, so add one
                    fit_max_dist = distance[int(start_imax):int(end_imax + 1)]
                    fit_min_dist = distance[int(start_imin):int(end_imin + 1)]
                    fit_max_gofr = gofr[int(start_imax):int(end_imax + 1)]
                    fit_min_gofr = gofr[int(start_imin):int(end_imin + 1)]
                    intgofr = intgofr[int(start_imin):int(end_imin + 1)]
        
                    nonzeros = np.count_nonzero(gofr)
                    # fit nth order polynomials to max and min regions
                    # only include gofrs with sufficient data and that doesn't go out of bounds
                    if nonzeros > 10 and end_imax < (orig_size - 1) and end_imin < (orig_size - 1):
                        order = 4
                        max_fit = np.poly1d(np.polyfit(fit_max_dist, fit_max_gofr, order))
                        min_fit = np.poly1d(np.polyfit(fit_min_dist, fit_min_gofr, order))
                        intgofr_fit = np.poly1d(np.polyfit(fit_min_dist, intgofr, order))
                        xaxis_max = np.linspace(distance[int(start_imax)], distance[int(end_imax)], 1000)  # create new x values for the fitted max region
                        xaxis_min = np.linspace(distance[int(start_imin)], distance[int(end_imin)],1000)  # create new x values for the fitted min region
                        yaxis_max = max_fit(xaxis_max)  # create new y values for fitted max region
                        yaxis_min = min_fit(xaxis_min)  # create new y values for fitted min region
                        yaxis_intgofr = intgofr_fit(xaxis_min)
                        poly_max_i = np.argmax(yaxis_max)  # determine index of max y value
                        poly_min_i = np.argmin(yaxis_min)  # determine index of min y value
                        
                        xmax_gofr = xaxis_max[int(poly_max_i)]  # x value of the max (x,y) of the max region of gofr
                        ymax_gofr = yaxis_max[int(poly_max_i)]  # max y value of the max region of gofr
                        xmin_gofr = xaxis_min[int(poly_min_i)]  # x value of the min (x,y) of the min region of gofr
                        y_intgofr = yaxis_intgofr[int(poly_min_i)]  # y value of the integral of gofr that corresponds to the min x value of gofr
                    
                        plt.plot(distance, gofr, 'o', xaxis_max, max_fit(xaxis_max), '-', xaxis_min, min_fit(xaxis_min), 'r-', markersize=2)
                               
                    else: # simply output 0's if the data is no good
                        xmax_gofr = 0
                        ymax_gofr = 0
                        xmin_gofr = 0
                        y_intgofr = 0
        
                    data.append(str(xmax_gofr))
                    data.append(str(ymax_gofr))
                    bonds.append(str(xmin_gofr))                      #for each couple of atom we add the corresponding xmin value
                    data.append(str(xmin_gofr))
                    data.append(str(y_intgofr))
                    guess_xmax[couple] = xmax_gofr                  #update the initial guesses for the next acell
                    guess_xmin[couple] = xmin_gofr
        
                    #re initialization of variables for the next loop
                    fit_min_dist = 0
                    fit_max_dist = 0
                    fit_min_gofr = 0
                    fit_max_gofr = 0
                    gofr = 0
                    intgofr = 0 
        if len(data) == size: #if data hasn't changed then we put X in the list
            data.extend(['X','X','X','X'])
            bonds.append('0')
            col += 1
        size = len(data) #we update the size of data   
    return(data,bonds,guess_xmax, guess_xmin)



def main(argv):
    """     ********* Main program *********     """
    atoms = []
    try:
        options,arg = getopt.getopt(argv,"ha:",["atoms"])
    except getopt.GetoptError:
        print("analyze_gofr.py -a <couples of atoms>(ex: 'Ca-O,Ca-Ca')")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('*******')
            print('analyze_gofr.py program to extract all relevant data from the gofr.dat files created by the script gofrs.py and write them into temperature_gofrs.txt files')
            print("analyze_gofr.py -a <couples of atoms>(ex: 'Ca-O,Ca-Ca')")
            print('WARNING! for efficiency, your data needs to be located in different folders for each T and you should launch this script from the folder containing every subfolder temperature')
            print('WARNING! you have to click FIRST on the maximum, and THEN on the minimum')
            print(' ')
            print('This code produces a gofrs.txt file with:')
            print('    - a column with filename')
            print('    - 4 column per A-B couple of atoms')
            print('        - x value of max(gofr)')
            print('        - y value associated')
            print('        - x value of the first min(gofr)')
            print('        - y value of int(gofr) associated to previous x value')
            print(' ')
            print(' ')
            sys.exit()
        elif opt in ("-a", "--atoms"):
            atoms = arg.split(',')                      #list of atom couples we want to analyze here
    for dirpath, dirnames, filenames in os.walk(os.curdir):
        files = sorted(glob.glob(dirpath+'/*.gofr.dat')) #I list every gofr files in alphabetic order
        if files != []:
            f, couples = headerfile(files[0], dirpath)                          #I create the first newfile for gofr and save the list of element couples 
            if atoms == []:
                atoms = couples[:]
            interactive = 0
            guess_xmax = {}                                                     #dictionnaries for initial guesses (key = couple of atoms, value = xmin or max value)
            guess_xmin = {}
            for file in files:
                data = [file]
                b, atom1, atom2 = atoms_columns(file)                            #I create the newfile for bonds
                if file == files[0]:
                    interactive = 1
                if interactive == 1:
                    data, bonds, guess_xmax, guess_xmin = analyze_gofrs_interactive(data,file,couples, guess_xmax, guess_xmin, atoms)
                    interactive = 0
                else:
                    data, bonds, guess_xmax, guess_xmin = analyze_gofrs_automatic(data,file,couples, guess_xmax, guess_xmin, atoms)       #we compute the min,max etc. from the gofr and int using fit and interactive plot
                f.write("\t".join(x for x in data)+ "\n")                  #we write in the file the result line
                writer = csv.writer(b, delimiter = '\t')  
                writer.writerows(zip(atom1, atom2, bonds))    #we write all our bonds.inp file
                b.close()
                print(b) 
            f.close
    
#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])



