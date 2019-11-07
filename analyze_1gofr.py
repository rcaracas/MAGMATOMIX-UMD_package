#!/usr/bin/env python
# -*- coding: utf-8 -*-

###
#AUTHORS: ANAIS KOBSCH, NATALIA SOLOMATOVA
###

#*********** Importation of the packages and modules used here ************
import sys
import getopt
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button


def headerfile(firstfile):
    """creation of the newfile with correct header"""
    firstline = ['method']  #beginning of the first line of the file gofr
    secondline = ['method']
    # creation of the header from the first line of the first file
    with open(firstfile, 'r') as f:
        header = f.readline()
    header = header.strip('dist')
    header = re.sub('(Int\([A-Za-z-]*\))', ' ',header).split()        #I substitute all Int(...)  by a blankspace and split using spaces
    for couple in header:
        for i in range(0,4):
            firstline.append(couple)
        secondline.extend(['xmax','ymax','xmin','coord','bond'])
    newfilename = 'gofr_' + firstfile.split('/')[-1] +'.txt' 
    print('The file ',newfilename,' is created')
    f = open(newfilename,'w')
    f.write("\t".join(x for x in firstline)+ "\n")
    f.write("\t".join(x for x in secondline)+ "\n")
    return f, header     #I return the newly created files f along with the list of element couples


def average_bond(radius,gofr,xmin_gofr):
    """compute [int(r^3 gofr(r)]/[r^2 gofr(r)] up to the 1st xmin"""
    ii = 0
    r3gofr = 0.0
    r2gofr = 0.0
    
    while radius[ii] < xmin_gofr:
        r3gofr = r3gofr + radius[ii]**3 * gofr[ii]
        r2gofr = r2gofr + radius[ii]**2 * gofr[ii]
        ii += 1
#    print('r3gofr, r2gofr, bond',)
    return r3gofr/r2gofr


def interactive_fit(data, data2, file,couples, distance, col, couple):
    """fit max and min of data using interactive plot for one couple of atoms"""
    print(couple, couples.index(couple) + col + 1, couples.index(couple) + col + 2)
    radius, gofr, intgofr = np.loadtxt(file,
                               usecols=(0, couples.index(couple) + col + 1, couples.index(couple) + col + 2),
                               skiprows=1, unpack=True)
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

    if nonzeros > 10 :
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.subplots_adjust(bottom=0.2)
        ax.set_title(file + ',  ' + couple)
        plt.xlabel('Distance ($\AA$)')
        plt.ylabel('g(r)')
        plt.grid()
        ax.plot(distance, gofr, 'o', markersize=2)  # plot the fitted max

        def onclick(event):
            global xmax_click, ymax_click
            xmax_click = event.xdata  # record x value of click
            ymax_click = event.ydata 
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
        ax.set_title(file + ',  ' + couple)
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
        
        y_intgofr2 = intgofr[min_index]     #the coordination number corresponding to the clicked position
        
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

        bond = average_bond(radius,gofr,xmin_gofr)

        class Decision(object):
            def good(self, event):
                global value_array
                plt.close(fig)  # close and move to next plot
                print('xmax=', round(xmax_gofr, 2), ' ymax=', round(ymax_gofr, 2), ' xmin=', round(xmin_gofr, 2), ' coord=', round(y_intgofr, 2), 'average bond=',round(bond,3))
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
        ax.set_title(file + ',  ' + couple)
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
    
    bond = average_bond(radius,gofr,xmin_gofr)
    
    data.append(str(xmax_gofr))
    data.append(str(ymax_gofr))
    data.append(str(xmin_gofr))
    data.append(str(y_intgofr))
    data.append(str(bond))
    
    data2.append(str(xmax_click))
    data2.append(str(ymax_click))
    data2.append(str(xmin_click))
    data2.append(str(y_intgofr2))
    data2.append(str(bond))
    
    cutdistance = 0  # re initialization of variables for the next loop
    fit_min_dist = 0
    fit_max_dist = 0
    fit_min_gofr = 0
    fit_max_gofr = 0
    gofr = 0
    intgofr = 0
    return (data, data2)



xmin_click = []
xmax_click = []
def analyze_gofrs_interactive(data,data2,file,couples, atoms):
    """ extraction of the min, max and coordination number using fits and interactive plot """
    distance = np.loadtxt(file, usecols=(0,), skiprows = 1, unpack=True)
    col = 0
    size = len(data)
    for couple in couples:
        for atom_couple in atoms:
            if (couple == atom_couple):
                data, data2 = interactive_fit(data,data2,file,couples, distance, col, couple)
        if len(data) == size: #if data hasn't changed then we put X in the list
            data.extend(['X','X','X','X','X'])
            data2.extend(['X','X','X','X','X'])
        size = len(data) #we update the size of data
        col += 1
    return(data, data2)



def main(argv):
    """     ********* Main program *********     """
    try:
        options,arg = getopt.getopt(argv,"hf:a:",["fgofrfile","atoms"])
    except getopt.GetoptError:
        print("analyze-gofr.py -f <gofr_filename> -a <couples of atoms>(ex: 'Ca-O,Ca-Ca')")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('*******')
            print('analyze-gofr.py program to extract all relevant data from 1 gofr.dat file created by the script gofrs.py and write them into .txt (fitted data and clicked values)')
            print("analyze-gofr.py -f <gofr_filename> -a <couples of atoms>(ex: 'Ca-O,Ca-Ca')")
            print('WARNING! you have to click FIRST on the maximum, and THEN on the minimum')
            print('WARNING!!!! If you want precise clicked value, you have to click on the correct position of maximum and minimum')
            print(' ')
            print('This code produces a .txt file with:')
            print('    - a column with the method used (fitted data of raw clicked value)')
            print('    - 4 column per A-B couple of atoms')
            print('        - x value of max(gofr)')
            print('        - y value associated')
            print('        - x value of the first min(gofr)')
            print('        - y value of int(gofr) associated to previous x value')
            print(' ')
            print('*******')
            print(' ')
            sys.exit()
        elif opt in ("-f", "--fgofrfile"):
            file = str(arg)
            #print('gofrfile = ',firstfile)
        elif opt in ("-a", "--atoms"):
            atoms = arg.split(',')                      #list of atom couples we want to analyze here)
    f, couples = headerfile(file)                          #I create the first newfile for gofr and save the list of element couples 
    data = ['fitted']
    data2 = ['clicked']
    data, data2 = analyze_gofrs_interactive(data, data2, file,couples, atoms)
    f.write("\t".join(x for x in data)+ "\n")                  #we write in the file the very last line
    f.write("\t".join(x for x in data2)+ "\n")                  #we write in the file the very last line
    f.close()

        
#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])



