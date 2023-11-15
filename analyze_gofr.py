#!/usr/bin/env python3
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

# *********** Importation of the packages and modules used here ************
import glob
import os
import sys
import getopt
import re
import numpy as np
import csv
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import matplotlib.ticker as ticker
from matplotlib.widgets import Button
import warnings


def headerfile(firstfile):
    """creation of the gofrs.txt file with correct header"""
    firstline = ['file']  # beginning of the first line of the file gofr
    secondline = ['file']
    with open(firstfile, 'r') as f:
        header = f.readline()
    header = header.strip('dist')
    header = re.sub('(Int\([A-Za-z1-9-]*\))', ' ', header).split()
    for couple in header:
        for i in range(0,4):
            firstline.append(couple)
        secondline.extend(['xmax','ymax','xmin','coord'])
    newfilename = 'gofrs_summary.txt'
    print('The file ',newfilename,' is created')
    f = open(newfilename,'w')
    f.write("\t".join(x for x in firstline)+ "\n")
    f.write("\t".join(x for x in secondline)+ "\n")
    return f, header


def atoms_columns_max(firstfile):
    """creation of a new file for maximum bonds"""
    atom1 = []
    atom2 = []
    with open(firstfile, 'r') as f:
        header = f.readline()
    header = header.strip('dist')
    header = re.sub('(Int\([A-Za-z1-9-]*\))', ' ', header).split()
    for couple in header:
        atom1.append(couple.split('-')[0])
        atom2.append(couple.split('-')[1])
    newfilenamebond = firstfile.split('.gofr.dat')[0] + '.maximum_bonds.inp'
    print('The maximum bonds file ', newfilenamebond, 'is created')
    b = open(newfilenamebond, 'w')
    return b, atom1, atom2


def atoms_columns_avg(firstfile):
    """creation of a new file for average bonds"""
    atom1_avg = []
    atom2_avg = []
    with open(firstfile, 'r') as f:
        header = f.readline()
    header = header.strip('dist')
    header = re.sub('(Int\([A-Za-z1-9-]*\))', ' ', header).split()
    for couple in header:
        atom1_avg.append(couple.split('-')[0])
        atom2_avg.append(couple.split('-')[1])
    newfilenamebond_avg = firstfile.split('.gofr.dat')[0] + '.average_bonds.inp'
    print('The average bonds file ', newfilenamebond_avg, 'is created')
    b_avg = open(newfilenamebond_avg, 'w')
    return b_avg, atom1_avg, atom2_avg


def atoms_columns_coord(firstfile):
    """creation of a new file for coordination"""
    atom1_coord = []
    atom2_coord = []
    with open(firstfile, 'r') as f:
        header = f.readline()
    header = header.strip('dist')
    header = re.sub('(Int\([A-Za-z1-9-]*\))', ' ', header).split()
    for couple in header:
        atom1_coord.append(couple.split('-')[0])
        atom2_coord.append(couple.split('-')[1])
    newfilenamebond_coord = firstfile.split('.gofr.dat')[0] + '.coordinations.inp'
    print('The coordinations file ', newfilenamebond_coord, 'is created')
    b_coord = open(newfilenamebond_coord, 'w')
    return b_coord, atom1_coord, atom2_coord

def Approximate_radius(element):
    element_list = ['H','He','Li','Be','B','C','N','O','F','Na','Mg','Al','Si','P','S','Cl','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es']
    radii_list = [2,1,90,59,41,30,132,126,119,116,86,68,54,58,170,167,152,114,89,75,72,94,97,92,89,89,74,87,88,76,87,184,182,166,132,104,86,86,83,75,76,74,100,108,109,94,83,90,207,206,62,181,149,117,115,113,143,111,136,131,108,106,105,104,103,102,101,100,85,86,80,77,77,82,94,151,133,164,133,117,108,76,194,162,126,108,116,116,124,114,140,111,110,109,93]
    if element in element_list:
        element_index = element_list.index(element)
        ionic_radius = radii_list[element_index]
    else:
        ionic_radius = 100
    return(ionic_radius)

def AtomPairsMinMax(atom1,atom2):
    ionic_radius1 = Approximate_radius(atom1)
    ionic_radius2 = Approximate_radius(atom2)
    sum_ionic_radii = ionic_radius1 + ionic_radius2
    if sum_ionic_radii < 140:
        fitting_ranges = [0.05, 0.08, 0.08, 0.12] # left_min, right_min, left_max, right_max
    if 140 <= sum_ionic_radii < 200:
        fitting_ranges = [0.08, 0.1, 0.12, 0.16]  # left_min, right_min, left_max, right_max
    if 200 <= sum_ionic_radii < 230:
        fitting_ranges = [0.08, 0.12, 0.12, 0.18]  # left_min, right_min, left_max, right_max
    if sum_ionic_radii >=230:
        fitting_ranges = [0.12, 0.18, 0.14, 0.25]  # left_min, right_min, left_max, right_max
    return(fitting_ranges)

def interactive_fit(data, file, couples, distance, col, max_bonds, peak_heights, avg_bonds, avg_coord, xmin_clicks, couple, list_atom1, list_atom2):
    """fit max and min of data using interactive plot for one couple of atoms"""
    gofr, intgofr = np.loadtxt(file,
                               usecols=(couples.index(couple) + col + 1, couples.index(couple) + col + 2),
                               skiprows=1, unpack=True)
    col = col + 1
    i = 0
    indices = []
    orig_size = np.size(distance)  # size of our data
    max_index = np.argmax(gofr)  # find the index of max gofr
    cutdistance = distance[max_index:-1]  # cut off beginning for finding minimum
    cutgofr = gofr[max_index:-1]  # cut off beginning for finding minimum
    cut_size = np.size(cutdistance)  # size of data
    cut_min_index = np.argmin(cutgofr)  # find the min of gofr in the cut segment
    min_index = orig_size - cut_size + cut_min_index - 1  # min index of the *uncut* gofr
    atom1 = couple.split('-')[0]
    atom2 = couple.split('-')[1]
    fitting_ranges = AtomPairsMinMax(atom1,atom2)
    left_max = fitting_ranges[0]
    right_max = fitting_ranges[1]
    left_min = fitting_ranges[2]
    right_min = fitting_ranges[3]

    ## Here I check if the gofr actually contains a peak/trough at physically-possible distances:
    # The longest reported bond length is between two atoms is 2.8 Ang (Bi-I)
    # If the max value is larger than 2.8 Angstroms, there is no first coord. shell (no bonding)
    # Thus, we only consider g(r)'s with first coordination shells and record 0's for non-existing bonds
    # If the max value of the gofr is less than 0.01 we don't have a gofr (it is flat), so 0's are recorded as well
    conditions = False
    max_index = np.argmax(gofr)  # determine index of max y value for initial analysis
    max_distance = distance[max_index]
    max_i = np.argmax(gofr)  # find the index of max gofr
    max_d = distance[max_i]
    if max_distance < 2.8 and max(gofr) > 0.01:
        # determine index of min value for initial analysis:
        end_index = (np.abs(distance - 5)).argmin()  # the first minimum should not be larger than 5 Ang.
        temp_cut_distance = distance[max_index:end_index]
        temp_cut_gofr = gofr[max_index:end_index]
        min_i = np.argmin(temp_cut_gofr)  # determine index of min y value
        min_distance = temp_cut_distance[min_i] # this is the estimated location of the first min
        if min_distance < 4:
            conditions = True

    ##  Here I record if the atom pair has already occurred in reverse order (e.g., Si-O vs. O-Si)
    list_atom1.append(atom1)
    list_atom2.append(atom2)
    this_is_a_repeat_pair = False
    for atom in list_atom2:
        # record atom1's that have already occurred as atom2's, omitting same-element pairs (Si-Si)
        if (atom == atom1) and (atom1 != atom2):
            indices.append(i)
        i = i + 1  # assign indices to the atom pairs, so that we can find them for repeat pairs (Si-O -> O-Si)
    for x in indices:
        # in the list of elements that have already occurred, check if the reverse couple has already occurred
        atom1_pair = list_atom1[x]
        if atom1_pair == atom2:
            this_is_a_repeat_pair = True
            matching_index = x

    if conditions == True:  # check if the peak and trough are located at physically-possible distances
        # if the pair has already been processed (e.g., Si-O vs. O-Si):
        # we use the same values for bond distances and max/min positions, but update the coordination
        if this_is_a_repeat_pair is True:
            max_bonds.append(max_bonds[int(matching_index)])
            avg_bonds.append(avg_bonds[int(matching_index)])
            peak_heights.append(peak_heights[int(matching_index)])
            xmin_click_pair = xmin_clicks[matching_index]
            min_index = (np.abs(distance - float(xmin_click_pair))).argmin()  # index of the closest value to clicked value
            # determine start and end cols for fitting the polynomial for the min region:
            start_imin = int(min_index - min_index * left_min)
            end_imin = int(min_index + min_index * right_min)
            if end_imin > (orig_size - 1):  # for rare cases where the min distance is out of bounds
                end_imin = orig_size - 1  # subtract one because cols start at 0
            # define the regions that need to be fitted for min area using the clicked min:
            fit_min_dist = distance[int(start_imin):int(end_imin + 1)]  # ranges do not include last element, so add one
            fit_min_gofr = gofr[int(start_imin):int(end_imin + 1)]
            intgofr = intgofr[int(start_imin):int(end_imin + 1)]
            warnings.filterwarnings('ignore')
            order = 3  # order of the polynomials
            min_fit = np.poly1d(np.polyfit(fit_min_dist, fit_min_gofr, order))
            intgofr_fit = np.poly1d(np.polyfit(fit_min_dist, intgofr, order))
            xaxis_min = np.linspace(distance[int(start_imin)], distance[int(end_imin)], 1000)  # create new x values for the fitted min region
            yaxis_min = min_fit(xaxis_min)  # create new y values for fitted min region
            poly_min_i = np.argmin(yaxis_min)  # determine index of min y value
            yaxis_intgofr = intgofr_fit(xaxis_min)
            y_intgofr = yaxis_intgofr[int(poly_min_i)]  # y value of the integral of gofr that corresponds to the min x value of gofr
            avg_coord.append(str(y_intgofr))
            xmin_clicks.append(xmin_clicks[int(matching_index)])
            data.append(str(avg_bonds[int(matching_index)]))  # equivalent to data.append(str(xmax_gofr))
            data.append(peak_heights[int(matching_index)])   # equivalent to data.append(str(ymax_gofr))
            data.append(str(max_bonds[int(matching_index)])) # equivalent to data.append(str(xmin_gofr))
            data.append(str(y_intgofr)) # equivalent to data.append(str(y_intgofr))

        else:
            plt.close()
            fig = plt.figure()
            ax = fig.add_subplot(111)
            plt.subplots_adjust(bottom=0.2)
            ax.set_title(file.split('.gofr.dat')[0])
            plt.text(0.03, 0.92, couple, fontsize=14, color='black', fontname='Arial', transform=ax.transAxes)
            ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
            plt.xlabel('Distance (Å)', fontsize=11.5, fontname='Arial', labelpad=1)
            plt.ylabel('Pair distribution function, g(r)', fontsize=12, fontname='Arial',labelpad=8)
            plt.ylim(-0.02, max(gofr) * 1.05)
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            ax.tick_params(labelsize=11)
            plt.xticks(fontname="Arial")
            plt.yticks(fontname="Arial")
            fig.subplots_adjust(bottom=0.25)
            fig.subplots_adjust(left=0.15)
            ax.plot(distance[0:-1], gofr[0:-1], 'o', markersize=2.5, color='grey', markerfacecolor='white')  # plot the fitted max

            plt.text(0.9, -0.27, 'Auto', fontsize=11, color='black', fontname='Arial', transform=ax.transAxes)
            ax.plot([0.865, 0.865, 1, 1, 0.865], [-0.195, -0.316, -0.316, -0.195, -0.195], marker='None', linestyle='-', linewidth=0.8,
                    color='black', transform=ax.transAxes, clip_on=False)
            ax.fill_between([0.865, 1], [-0.195, -0.195], [-0.316, -0.316], color='#DBDBDB', transform=ax.transAxes, clip_on=False)

            def onclick(event):
                global xmax_click
                xmax_click = event.xdata  # record x value of click
                if isinstance(xmax_click, float) == False:
                    xmax_click = float(max_distance)
                plt.close()

            fig.canvas.mpl_connect('button_press_event', onclick)
            plt.show()
            max_index = (np.abs(distance - xmax_click)).argmin()  # index of the closest value to clicked value

            # determine start and end cols for fitting the polynomial:
            end_imax = int(max_index + max_index * right_max)
            start_imax = int(max_index - max_index * left_max)
            if end_imax > (orig_size - 1):  # for rare cases where the max distance is out of bounds
                end_imax = orig_size - 1  # subtract one because cols start at 0

            # define the regions that need to be fitted for the max area:
            fit_max_dist = distance[int(start_imax):int(end_imax + 1)]  # ranges do not include last element, so add one
            fit_max_gofr = gofr[int(start_imax):int(end_imax + 1)]

            order = 3  # order of the polynomials
            max_fit = np.poly1d(np.polyfit(fit_max_dist, fit_max_gofr, order))
            xaxis_max = np.linspace(distance[int(start_imax)], distance[int(end_imax)], 1000)  # create new x values
            yaxis_max = max_fit(xaxis_max)  # create new y values for fitted max region
            poly_max_i = np.argmax(yaxis_max)  # determine index of max y value

            xmax_gofr = xaxis_max[int(poly_max_i)]  # x value of the max (x,y) of the max region of gofr
            ymax_gofr = yaxis_max[int(poly_max_i)]  # max y value of the max region of gofr

            fig = plt.figure()
            ax = fig.add_subplot(111)
            plt.subplots_adjust(bottom=0.2)
            ax.set_title(file.split('.gofr.dat')[0])
            plt.text(0.03, 0.92, couple, fontsize=14, color='black', fontname='Arial', transform=ax.transAxes)
            ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
            ax.plot(distance[0:-1], gofr[0:-1], 'o',  markersize=2.5, color='grey', markerfacecolor='white')
            ax.plot(xaxis_max, max_fit(xaxis_max), '-', color='blue', linewidth=2)
            ax.plot(xmax_gofr, ymax_gofr, 'o', markersize=6, color='blue')
            plt.xlabel('Distance (Å)', fontsize=11.5, fontname='Arial', labelpad=1)
            plt.ylabel('Pair distribution function, g(r)', fontsize=12, fontname='Arial',labelpad=8)
            plt.ylim(-0.02, max(gofr) * 1.05)
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            ax.tick_params(labelsize=11)
            plt.xticks(fontname="Arial")
            plt.yticks(fontname="Arial")
            fig.subplots_adjust(bottom=0.25)
            fig.subplots_adjust(left=0.15)

            plt.text(0.9, -0.27, 'Auto', fontsize=11, color='black', fontname='Arial', transform=ax.transAxes)
            ax.plot([0.865, 0.865, 1, 1, 0.865], [-0.195, -0.316, -0.316, -0.195, -0.195], marker='None', linestyle='-', linewidth=0.8,
                    color='black', transform=ax.transAxes, clip_on=False)
            ax.fill_between([0.865, 1], [-0.195, -0.195], [-0.316, -0.316], color='#DBDBDB', transform=ax.transAxes, clip_on=False)

            def onclick(event):
                global xmin_click
                xmin_click = event.xdata  # record x value of click
                if isinstance(xmin_click, float) == False:
                    end_index = (np.abs(distance - 5)).argmin() # the first minimum should not be larger than 5 Ang.
                    temporary_cut_distance = distance[max_index:end_index]
                    temporary_cut_gofr = gofr[max_index:end_index]
                    min_index = np.argmin(temporary_cut_gofr)  # determine index of min y value
                    xmin_click = temporary_cut_distance[min_index]
                plt.close()

            fig.canvas.mpl_connect('button_press_event', onclick)
            plt.show()
            min_index = (np.abs(distance - xmin_click)).argmin()  # index of the closest value to clicked value
            # determine start and end cols for fitting the polynomial for the min region:
            start_imin = int(min_index - min_index * left_min)
            end_imin = int(min_index + min_index * right_min)
            if end_imin > (orig_size - 1):  # for rare cases where the min distance is out of bounds
                end_imin = orig_size - 1  # subtract one because cols start at 0
            # define the regions that need to be fitted for min area using the clicked min:
            fit_min_dist = distance[int(start_imin):int(end_imin + 1)]  # ranges do not include last element, so add one
            fit_min_gofr = gofr[int(start_imin):int(end_imin + 1)]
            cut_intgofr = intgofr[int(start_imin):int(end_imin + 1)]

            min_fit = np.poly1d(np.polyfit(fit_min_dist, fit_min_gofr, order))
            intgofr_fit = np.poly1d(np.polyfit(fit_min_dist, cut_intgofr, order))
            xaxis_min = np.linspace(distance[int(start_imin)], distance[int(end_imin)], 1000)
            yaxis_min = min_fit(xaxis_min)  # create new y values for fitted min region
            poly_min_i = np.argmin(yaxis_min)  # determine index of min y value
            yaxis_intgofr = intgofr_fit(xaxis_min)

            xmin_gofr = xaxis_min[int(poly_min_i)]  # x value of the min (x,y) of the min region of gofr
            ymin_gofr = yaxis_min[int(poly_min_i)]  # x value of the min (x,y) of the min region of gofr
            y_intgofr = yaxis_intgofr[int(poly_min_i)]  # y value of the integral of gofr that corresponds to the min x value of gofr

            class Decision(object):
                def good(self, event):
                    global value_array
                    plt.close(fig)  # close and move to next plot
                    print(atom1, '-', atom2, '  average bond =', "%.2f" % round(xmax_gofr, 2), '  maximum bond =', "%.2f" % round(xmin_gofr, 2))
                    xmin_click_save = xmin_click
                    undo = False
                    manual_click = False
                    value_array = (xmax_gofr, ymax_gofr, xmin_gofr, y_intgofr, xmin_click_save, undo, manual_click)
                    return value_array

                def bad(self, event):
                    global value_array
                    plt.close(fig)  # close and move to next plot
                    print('Zeros will be placed for the bond distances for', couple)
                    xmax_gofr = 0
                    ymax_gofr = 0
                    xmin_gofr = 0
                    y_intgofr = 0
                    xmin_click_save = 0
                    undo = False
                    manual_click = False
                    value_array = (xmax_gofr, ymax_gofr, xmin_gofr, y_intgofr, xmin_click_save, undo, manual_click)
                    return value_array

                def redo(self, event):
                    global value_array
                    plt.close(fig)  # close and move to next plot
                    print('Please redo the fitting procedure for', couple)
                    xmax_gofr = 0
                    ymax_gofr = 0
                    xmin_gofr = 0
                    y_intgofr = 0
                    xmin_click_save = 0
                    undo = True
                    manual_click = False
                    value_array = (xmax_gofr, ymax_gofr, xmin_gofr, y_intgofr, xmin_click_save, undo, manual_click)
                    return value_array

                def manual(self, event):
                    global value_array
                    plt.close(fig)  # close and move to next plot
                    print('Please redo the fitting procedure for', couple)
                    xmax_gofr = 0
                    ymax_gofr = 0
                    xmin_gofr = 0
                    y_intgofr = 0
                    xmin_click_save = 0
                    undo = True
                    manual_click = True
                    value_array = (xmax_gofr, ymax_gofr, xmin_gofr, y_intgofr, xmin_click_save, undo, manual_click)
                    return value_array

            fig = plt.figure()  # show figure with the newly fitted minimum curve
            ax = fig.add_subplot(111)
            plt.subplots_adjust(bottom=0.2)
            ax.set_title(file.split('.gofr.dat')[0])
            plt.text(0.03, 0.92, couple, fontsize=14, color='black', fontname='Arial', transform=ax.transAxes)
            ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
            ax.plot(distance[0:-1], gofr[0:-1], 'o',  markersize=2.5, color='grey', markerfacecolor='white')
            ax.plot(xaxis_max, max_fit(xaxis_max), '-', color='blue', linewidth=2)
            ax.plot(xmax_gofr, ymax_gofr, 'o', markersize=6, color='blue')
            ax.plot(xaxis_min, min_fit(xaxis_min), '-', color='red', linewidth=2)
            ax.plot(xmin_gofr, ymin_gofr, 'o', markersize=6, color='red')

            plt.xlabel('Distance (Å)', fontsize=11.5, fontname='Arial', labelpad=1)
            plt.ylabel('Pair distribution function, g(r)', fontsize=12, fontname='Arial',labelpad=8)
            plt.ylim(-0.02, max(gofr) * 1.05)
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            ax.tick_params(labelsize=11)
            plt.xticks(fontname="Arial")
            plt.yticks(fontname="Arial")
            fig.subplots_adjust(bottom=0.25)
            fig.subplots_adjust(left=0.15)


            class Exit(object):
                def exit(self, event):
                    print('Exiting the fitting procedure.')
                    plt.close(fig)
                    sys.exit()

            callback = Decision()

            click_Exit = Exit()
            bexit = Button(plt.axes([0.15, 0.05, 0.1, 0.078]), 'Quit')
            bexit.on_clicked(click_Exit.exit)

            bmanual = Button(plt.axes([0.27, 0.05, 0.12, 0.078]), 'Redo -\nmanually')
            bmanual.on_clicked(callback.manual)
            bmanual.ax.set_facecolor('#E2E3A0')
            bmanual.label.set_fontsize(9.5)

            bredo = Button(plt.axes([0.41, 0.05, 0.12, 0.078]), 'Redo -\nfitting')
            bredo.on_clicked(callback.redo)
            bredo.ax.set_facecolor('#E2E3A0')
            bredo.label.set_fontsize(9.5)

            bbad = Button(plt.axes([0.55, 0.05, 0.105, 0.078]), 'Discard')
            bbad.on_clicked(callback.bad)
            bbad.ax.set_facecolor('#F1B6A6')

            bgood = Button(plt.axes([0.675, 0.05, 0.105, 0.078]), 'Save') #left, bottom, width, height
            bgood.on_clicked(callback.good)
            bgood.ax.set_facecolor('#98D49A')

            plt.text(0.9, -0.27, 'Auto', fontsize=11, color='grey', fontname='Arial', transform=ax.transAxes)
            ax.plot([0.865, 0.865, 1, 1, 0.865], [-0.195, -0.316, -0.316, -0.195, -0.195], marker='None', linestyle='-', linewidth=0.8,
                    color='grey', transform=ax.transAxes, clip_on=False)
            ax.fill_between([0.865, 1], [-0.195, -0.195], [-0.316, -0.316], color='#DBDBDB', transform=ax.transAxes, clip_on=False)

            plt.show()

            undo = value_array[5]
            manual_click = value_array[6]

            while undo == True:
                if manual_click == False:
                    plt.close()
                    fig = plt.figure()
                    ax = fig.add_subplot(111)
                    plt.subplots_adjust(bottom=0.2)
                    ax.set_title(file.split('.gofr.dat')[0])
                    plt.text(0.03, 0.92, couple, fontsize=14, color='black', fontname='Arial', transform=ax.transAxes)
                    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
                    plt.xlabel('Distance (Å)', fontsize=11.5, fontname='Arial', labelpad=1)
                    plt.ylabel('Pair distribution function, g(r)', fontsize=12, fontname='Arial', labelpad=8)
                    plt.ylim(-0.02, max(gofr) * 1.05)
                    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                    ax.tick_params(labelsize=11)
                    plt.xticks(fontname="Arial")
                    plt.yticks(fontname="Arial")
                    fig.subplots_adjust(bottom=0.25)
                    fig.subplots_adjust(left=0.15)
                    ax.plot(distance[0:-1], gofr[0:-1], 'o', markersize=2.5, color='grey', markerfacecolor='white')  # plot the fitted max

                    plt.text(0.9, -0.27, 'Auto', fontsize=11, color='black', fontname='Arial', transform=ax.transAxes)
                    ax.plot([0.865, 0.865, 1, 1, 0.865], [-0.195, -0.316, -0.316, -0.195, -0.195], marker='None', linestyle='-', linewidth=0.8, color='black', transform=ax.transAxes, clip_on=False)
                    ax.fill_between([0.865, 1], [-0.195, -0.195], [-0.316, -0.316], color='#DBDBDB', transform=ax.transAxes, clip_on=False)

                    def onclick(event):
                        global xmax_click
                        xmax_click = event.xdata  # record x value of click
                        if isinstance(xmax_click, float) == False:
                            xmax_click = float(max_distance)
                        plt.close()

                    fig.canvas.mpl_connect('button_press_event', onclick)
                    plt.show()
                    max_index = (np.abs(distance - xmax_click)).argmin()  # index of the closest value to clicked value

                    # determine start and end cols for fitting the polynomial:
                    end_imax = int(max_index + max_index * right_max)
                    start_imax = int(max_index - max_index * left_max)
                    if end_imax > (orig_size - 1):  # for rare cases where the max distance is out of bounds
                        end_imax = orig_size - 1  # subtract one because cols start at 0

                    # define the regions that need to be fitted for the max area:
                    fit_max_dist = distance[int(start_imax):int(end_imax + 1)]
                    fit_max_gofr = gofr[int(start_imax):int(end_imax + 1)]

                    order = 3  # order of the polynomials
                    max_fit = np.poly1d(np.polyfit(fit_max_dist, fit_max_gofr, order))
                    xaxis_max = np.linspace(distance[int(start_imax)], distance[int(end_imax)], 1000)
                    yaxis_max = max_fit(xaxis_max)  # create new y values for fitted max region
                    poly_max_i = np.argmax(yaxis_max)  # determine index of max y value

                    xmax_gofr = xaxis_max[int(poly_max_i)]  # x value of the max (x,y) of the max region of gofr
                    ymax_gofr = yaxis_max[int(poly_max_i)]  # max y value of the max region of gofr

                    fig = plt.figure()
                    ax = fig.add_subplot(111)
                    plt.subplots_adjust(bottom=0.2)
                    ax.set_title(file.split('.gofr.dat')[0])
                    plt.text(0.03, 0.92, couple, fontsize=14, color='black', fontname='Arial', transform=ax.transAxes)
                    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
                    ax.plot(distance[0:-1], gofr[0:-1], 'o',  markersize=2.5, color='grey', markerfacecolor='white')
                    ax.plot(xaxis_max, max_fit(xaxis_max), '-', color='blue', linewidth=2)
                    ax.plot(xmax_gofr, ymax_gofr, 'o', markersize=6, color='blue')
                    plt.xlabel('Distance (Å)', fontsize=11.5, fontname='Arial', labelpad=1)
                    plt.ylabel('Pair distribution function, g(r)', fontsize=12, fontname='Arial', labelpad=8)
                    plt.ylim(-0.02, max(gofr) * 1.05)
                    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                    ax.tick_params(labelsize=11)
                    plt.xticks(fontname="Arial")
                    plt.yticks(fontname="Arial")
                    fig.subplots_adjust(bottom=0.25)
                    fig.subplots_adjust(left=0.15)

                    plt.text(0.9, -0.27, 'Auto', fontsize=11, color='black', fontname='Arial', transform=ax.transAxes)
                    ax.plot([0.865, 0.865, 1, 1, 0.865], [-0.195, -0.316, -0.316, -0.195, -0.195], marker='None', linestyle='-', linewidth=0.8, color='black', transform=ax.transAxes, clip_on=False)
                    ax.fill_between([0.865, 1], [-0.195, -0.195], [-0.316, -0.316], color='#DBDBDB', transform=ax.transAxes, clip_on=False)

                    def onclick(event):
                        global xmin_click
                        xmin_click = event.xdata  # record x value of click
                        if isinstance(xmin_click, float) == False:
                            end_index = (np.abs(distance - 5)).argmin() # the first minimum should not be larger than 5 Ang.
                            temporary_cut_distance = distance[max_index:end_index]
                            temporary_cut_gofr = gofr[max_index:end_index]
                            min_index = np.argmin(temporary_cut_gofr)  # determine index of min y value
                            xmin_click = temporary_cut_distance[min_index]
                        plt.close()

                    fig.canvas.mpl_connect('button_press_event', onclick)
                    plt.show()
                    min_index = (np.abs(distance - xmin_click)).argmin()  # index of the closest value to clicked value
                    # determine start and end cols for fitting the polynomial for the min region:
                    start_imin = int(min_index - min_index * left_min)
                    end_imin = int(min_index + min_index * right_min)
                    if end_imin > (orig_size - 1):  # for rare cases where the min distance is out of bounds
                        end_imin = orig_size - 1  # subtract one because cols start at 0
                    # define the regions that need to be fitted for min area using the clicked min:
                    fit_min_dist = distance[int(start_imin):int(end_imin + 1)]
                    fit_min_gofr = gofr[int(start_imin):int(end_imin + 1)]
                    cut_intgofr = intgofr[int(start_imin):int(end_imin + 1)]

                    min_fit = np.poly1d(np.polyfit(fit_min_dist, fit_min_gofr, order))
                    intgofr_fit = np.poly1d(np.polyfit(fit_min_dist, cut_intgofr, order))
                    xaxis_min = np.linspace(distance[int(start_imin)], distance[int(end_imin)], 1000)  # create new x values
                    yaxis_min = min_fit(xaxis_min)  # create new y values for fitted min region
                    poly_min_i = np.argmin(yaxis_min)  # determine index of min y value
                    yaxis_intgofr = intgofr_fit(xaxis_min)

                    xmin_gofr = xaxis_min[int(poly_min_i)]  # x value of the min (x,y) of the min region of gofr
                    ymin_gofr = yaxis_min[int(poly_min_i)]  # x value of the min (x,y) of the min region of gofr
                    y_intgofr = yaxis_intgofr[int(poly_min_i)]  # y value of the integral of gofr that corresponds to the min x value of gofr

                    class Decision(object):
                        def good(self, event):
                            global value_array
                            plt.close(fig)  # close and move to next plot
                            print(
                            atom1, '-', atom2, '  average bond =', "%.2f" % round(xmax_gofr, 2), '  maximum bond =',
                            "%.2f" % round(xmin_gofr, 2))
                            xmin_click_save = xmin_click
                            undo = False
                            manual_click = False
                            value_array = (
                            xmax_gofr, ymax_gofr, xmin_gofr, y_intgofr, xmin_click_save, undo, manual_click)
                            return value_array

                        def bad(self, event):
                            global value_array
                            plt.close(fig)  # close and move to next plot
                            print('Zeros will be placed for the bond distances for', couple)
                            xmax_gofr = 0
                            ymax_gofr = 0
                            xmin_gofr = 0
                            y_intgofr = 0
                            xmin_click_save = 0
                            undo = False
                            manual_click = False
                            value_array = (
                            xmax_gofr, ymax_gofr, xmin_gofr, y_intgofr, xmin_click_save, undo, manual_click)
                            return value_array

                        def redo(self, event):
                            global value_array
                            plt.close(fig)  # close and move to next plot
                            print('Please redo the fitting procedure for', couple)
                            xmax_gofr = 0
                            ymax_gofr = 0
                            xmin_gofr = 0
                            y_intgofr = 0
                            xmin_click_save = 0
                            undo = True
                            manual_click = False
                            value_array = (
                            xmax_gofr, ymax_gofr, xmin_gofr, y_intgofr, xmin_click_save, undo, manual_click)
                            return value_array

                        def manual(self, event):
                            global value_array
                            plt.close(fig)  # close and move to next plot
                            print('Please redo the fitting procedure for', couple)
                            xmax_gofr = 0
                            ymax_gofr = 0
                            xmin_gofr = 0
                            y_intgofr = 0
                            xmin_click_save = 0
                            undo = True
                            manual_click = True
                            value_array = (
                            xmax_gofr, ymax_gofr, xmin_gofr, y_intgofr, xmin_click_save, undo, manual_click)
                            return value_array

                    fig = plt.figure()  # show figure with the newly fitted minimum curve
                    ax = fig.add_subplot(111)
                    plt.subplots_adjust(bottom=0.2)
                    ax.set_title(file.split('.gofr.dat')[0])
                    plt.text(0.03, 0.92, couple, fontsize=14, color='black', fontname='Arial', transform=ax.transAxes)
                    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
                    ax.plot(distance[0:-1], gofr[0:-1], 'o',  markersize=2.5, color='grey', markerfacecolor='white')
                    ax.plot(xaxis_max, max_fit(xaxis_max), '-', color='blue', linewidth=2)
                    ax.plot(xmax_gofr, ymax_gofr, 'o', markersize=6, color='blue')
                    ax.plot(xaxis_min, min_fit(xaxis_min), '-', color='red', linewidth=2)
                    ax.plot(xmin_gofr, ymin_gofr, 'o', markersize=6, color='red')

                    plt.xlabel('Distance (Å)', fontsize=11.5, fontname='Arial', labelpad=1)
                    plt.ylabel('Pair distribution function, g(r)', fontsize=12, fontname='Arial', labelpad=8)
                    plt.ylim(-0.02, max(gofr) * 1.05)
                    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                    ax.tick_params(labelsize=11)
                    plt.xticks(fontname="Arial")
                    plt.yticks(fontname="Arial")
                    fig.subplots_adjust(bottom=0.25)
                    fig.subplots_adjust(left=0.15)

                    class Exit(object):
                        def exit(self, event):
                            print('You exited the fitting procedure.')
                            plt.close(fig)
                            sys.exit()

                    callback = Decision()

                    click_Exit = Exit()
                    bexit = Button(plt.axes([0.15, 0.05, 0.1, 0.078]), 'Quit')
                    bexit.on_clicked(click_Exit.exit)

                    bmanual = Button(plt.axes([0.27, 0.05, 0.12, 0.078]), 'Redo -\nmanually')
                    bmanual.on_clicked(callback.manual)
                    bmanual.ax.set_facecolor('#E2E3A0')
                    bmanual.label.set_fontsize(9.5)

                    bredo = Button(plt.axes([0.41, 0.05, 0.12, 0.078]), 'Redo -\nfitting')
                    bredo.on_clicked(callback.redo)
                    bredo.ax.set_facecolor('#E2E3A0')
                    bredo.label.set_fontsize(9.5)

                    bbad = Button(plt.axes([0.55, 0.05, 0.105, 0.078]), 'Discard')
                    bbad.on_clicked(callback.bad)
                    bbad.ax.set_facecolor('#F1B6A6')

                    bgood = Button(plt.axes([0.675, 0.05, 0.105, 0.078]), 'Save')  # left, bottom, width, height
                    bgood.on_clicked(callback.good)
                    bgood.ax.set_facecolor('#98D49A')

                    plt.text(0.9, -0.27, 'Auto', fontsize=11, color='grey', fontname='Arial', transform=ax.transAxes)
                    ax.plot([0.865, 0.865, 1, 1, 0.865], [-0.195, -0.316, -0.316, -0.195, -0.195], marker='None', linestyle='-', linewidth=0.8, color='grey', transform=ax.transAxes, clip_on=False)
                    ax.fill_between([0.865, 1], [-0.195, -0.195], [-0.316, -0.316], color='#DBDBDB', transform=ax.transAxes, clip_on=False)

                    plt.show()

                    undo = value_array[5]
                    manual_click = value_array[6]

                else: ### if manual_click == True
                    plt.close()
                    fig = plt.figure()
                    ax = fig.add_subplot(111)
                    plt.subplots_adjust(bottom=0.2)
                    ax.set_title(file.split('.gofr.dat')[0])
                    plt.text(0.03, 0.92, couple, fontsize=14, color='black', fontname='Arial', transform=ax.transAxes)
                    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
                    plt.xlabel('Distance (Å)', fontsize=11.5, fontname='Arial', labelpad=1)
                    plt.ylabel('Pair distribution function, g(r)', fontsize=12, fontname='Arial', labelpad=8)
                    plt.ylim(-0.02, max(gofr) * 1.05)
                    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                    ax.tick_params(labelsize=11)
                    plt.xticks(fontname="Arial")
                    plt.yticks(fontname="Arial")
                    fig.subplots_adjust(bottom=0.25)
                    fig.subplots_adjust(left=0.15)
                    ax.plot(distance[0:-1], gofr[0:-1], 'o', markersize=2.5, color='grey', markerfacecolor='white')  # plot the fitted max

                    plt.text(0.9, -0.27, 'Auto', fontsize=11, color='black', fontname='Arial', transform=ax.transAxes)
                    ax.plot([0.865, 0.865, 1, 1, 0.865], [-0.195, -0.316, -0.316, -0.195, -0.195], marker='None', linestyle='-', linewidth=0.8, color='black', transform=ax.transAxes, clip_on=False)
                    ax.fill_between([0.865, 1], [-0.195, -0.195], [-0.316, -0.316], color='#DBDBDB', transform=ax.transAxes, clip_on=False)

                    def onclick(event):
                        global xmax_click
                        global ymax_click
                        xmax_click = event.xdata  # record x value of click
                        ymax_click = event.ydata  # record x value of click
                        if isinstance(xmax_click, float) == False:
                            xmax_click = float(max_distance)
                            max_index = (np.abs(distance - xmax_click)).argmin()
                            ymax_click = gofr[max_index]
                        plt.close()

                    fig.canvas.mpl_connect('button_press_event', onclick)
                    plt.show()
                    max_index = (np.abs(distance - xmax_click)).argmin()  # index of the closest value to clicked value

                    xmax_gofr = distance[max_index]  # x value of the max (x,y) of the max region of gofr
                    ymax_gofr = gofr[max_index]  # max y value of the max region of gofr

                    fig = plt.figure()
                    ax = fig.add_subplot(111)
                    plt.subplots_adjust(bottom=0.2)
                    ax.set_title(file.split('.gofr.dat')[0])
                    plt.text(0.03, 0.92, couple, fontsize=14, color='black', fontname='Arial', transform=ax.transAxes)
                    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
                    ax.plot(distance[0:-1], gofr[0:-1], 'o',  markersize=2.5, color='grey', markerfacecolor='white')
                    ax.plot(xmax_click, ymax_click, 'o', markersize=6, color='blue')
                    plt.xlabel('Distance (Å)', fontsize=11.5, fontname='Arial', labelpad=1)
                    plt.ylabel('Pair distribution function, g(r)', fontsize=12, fontname='Arial', labelpad=8)
                    plt.ylim(-0.02, max(gofr) * 1.05)
                    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                    ax.tick_params(labelsize=11)
                    plt.xticks(fontname="Arial")
                    plt.yticks(fontname="Arial")
                    fig.subplots_adjust(bottom=0.25)
                    fig.subplots_adjust(left=0.15)

                    plt.text(0.9, -0.27, 'Auto', fontsize=11, color='black', fontname='Arial', transform=ax.transAxes)
                    ax.plot([0.865, 0.865, 1, 1, 0.865], [-0.195, -0.316, -0.316, -0.195, -0.195], marker='None', linestyle='-', linewidth=0.8, color='black', transform=ax.transAxes, clip_on=False)
                    ax.fill_between([0.865, 1], [-0.195, -0.195], [-0.316, -0.316], color='#DBDBDB', transform=ax.transAxes, clip_on=False)

                    def onclick(event):
                        global xmin_click
                        global ymin_click
                        xmin_click = event.xdata  # record x value of click
                        ymin_click = event.ydata  # record x value of click
                        if isinstance(xmin_click, float) == False:
                            end_index = (np.abs(distance - 5)).argmin() # the first minimum should not be larger than 5 Ang.
                            temporary_cut_distance = distance[max_index:end_index]
                            temporary_cut_gofr = gofr[max_index:end_index]
                            min_index = np.argmin(temporary_cut_gofr)  # determine index of min y value
                            xmin_click = temporary_cut_distance[min_index]
                            ymin_click = temporary_cut_gofr[min_index]
                        plt.close()

                    fig.canvas.mpl_connect('button_press_event', onclick)
                    plt.show()
                    min_index = (np.abs(distance - xmin_click)).argmin()  # index of the closest value to clicked value

                    xmin_gofr = distance[min_index]  # x value of the min (x,y) of the min region of gofr
                    ymin_gofr = gofr[min_index]  # x value of the min (x,y) of the min region of gofr
                    y_intgofr = intgofr[min_index]  # y value of the integral of gofr that corresponds to the min x value of gofr

                    class Decision(object):
                        def good(self, event):
                            global value_array
                            plt.close(fig)  # close and move to next plot
                            print(
                            atom1, '-', atom2, '  average bond =', "%.2f" % round(xmax_gofr, 2), '  maximum bond =',
                            "%.2f" % round(xmin_gofr, 2))
                            xmin_click_save = xmin_click
                            undo = False
                            manual_click = False
                            value_array = (
                            xmax_gofr, ymax_gofr, xmin_gofr, y_intgofr, xmin_click_save, undo, manual_click)
                            return value_array

                        def bad(self, event):
                            global value_array
                            plt.close(fig)  # close and move to next plot
                            print('Zeros will be placed for the bond distances for', couple)
                            xmax_gofr = 0
                            ymax_gofr = 0
                            xmin_gofr = 0
                            y_intgofr = 0
                            xmin_click_save = 0
                            undo = False
                            manual_click = False
                            value_array = (
                            xmax_gofr, ymax_gofr, xmin_gofr, y_intgofr, xmin_click_save, undo, manual_click)
                            return value_array

                        def redo(self, event):
                            global value_array
                            plt.close(fig)  # close and move to next plot
                            print('Please redo the fitting procedure for', couple)
                            xmax_gofr = 0
                            ymax_gofr = 0
                            xmin_gofr = 0
                            y_intgofr = 0
                            xmin_click_save = 0
                            undo = True
                            manual_click = False
                            value_array = (
                            xmax_gofr, ymax_gofr, xmin_gofr, y_intgofr, xmin_click_save, undo, manual_click)
                            return value_array

                        def manual(self, event):
                            global value_array
                            plt.close(fig)  # close and move to next plot
                            print('Please redo the fitting procedure for', couple)
                            xmax_gofr = 0
                            ymax_gofr = 0
                            xmin_gofr = 0
                            y_intgofr = 0
                            xmin_click_save = 0
                            undo = True
                            manual_click = True
                            value_array = (
                            xmax_gofr, ymax_gofr, xmin_gofr, y_intgofr, xmin_click_save, undo, manual_click)
                            return value_array

                    fig = plt.figure()  # show figure with the newly fitted minimum curve
                    ax = fig.add_subplot(111)
                    plt.subplots_adjust(bottom=0.2)
                    ax.set_title(file.split('.gofr.dat')[0])
                    plt.text(0.03, 0.92, couple, fontsize=14, color='black', fontname='Arial', transform=ax.transAxes)
                    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
                    ax.plot(distance[0:-1], gofr[0:-1], 'o',  markersize=2.5, color='grey', markerfacecolor='white')
                    ax.plot(xmax_click, ymax_click, 'o', markersize=6, color='blue')
                    ax.plot(xmin_click, ymin_click, 'o', markersize=6, color='red')
                    plt.xlabel('Distance (Å)', fontsize=11.5, fontname='Arial', labelpad=1)
                    plt.ylabel('Pair distribution function, g(r)', fontsize=12, fontname='Arial', labelpad=8)
                    plt.ylim(-0.02, max(gofr) * 1.05)
                    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                    ax.tick_params(labelsize=11)
                    plt.xticks(fontname="Arial")
                    plt.yticks(fontname="Arial")
                    fig.subplots_adjust(bottom=0.25)
                    fig.subplots_adjust(left=0.15)

                    class Exit(object):
                        def exit(self, event):
                            print('You exited the fitting procedure.')
                            plt.close(fig)
                            sys.exit()

                    callback = Decision()

                    click_Exit = Exit()
                    bexit = Button(plt.axes([0.15, 0.05, 0.1, 0.078]), 'Quit')
                    bexit.on_clicked(click_Exit.exit)

                    bmanual = Button(plt.axes([0.27, 0.05, 0.12, 0.078]), 'Redo -\nmanually')
                    bmanual.on_clicked(callback.manual)
                    bmanual.ax.set_facecolor('#E2E3A0')
                    bmanual.label.set_fontsize(9.5)

                    bredo = Button(plt.axes([0.41, 0.05, 0.12, 0.078]), 'Redo -\nfitting')
                    bredo.on_clicked(callback.redo)
                    bredo.ax.set_facecolor('#E2E3A0')
                    bredo.label.set_fontsize(9.5)

                    bbad = Button(plt.axes([0.55, 0.05, 0.105, 0.078]), 'Discard')
                    bbad.on_clicked(callback.bad)
                    bbad.ax.set_facecolor('#F1B6A6')

                    bgood = Button(plt.axes([0.675, 0.05, 0.105, 0.078]), 'Save')  # left, bottom, width, height
                    bgood.on_clicked(callback.good)
                    bgood.ax.set_facecolor('#98D49A')

                    plt.text(0.9, -0.27, 'Auto', fontsize=11, color='grey', fontname='Arial', transform=ax.transAxes)
                    ax.plot([0.865, 0.865, 1, 1, 0.865], [-0.195, -0.316, -0.316, -0.195, -0.195], marker='None', linestyle='-', linewidth=0.8, color='grey', transform=ax.transAxes, clip_on=False)
                    ax.fill_between([0.865, 1], [-0.195, -0.195], [-0.316, -0.316], color='#DBDBDB', transform=ax.transAxes, clip_on=False)

                    plt.show()

                    undo = value_array[5]
                    manual_click = value_array[6]

            xmax_gofr = value_array[0]
            ymax_gofr = value_array[1]
            xmin_gofr = value_array[2]
            y_intgofr = value_array[3]
            xmin_click_save = value_array[4]
            max_bonds.append(str(xmin_gofr))
            peak_heights.append(str(ymax_gofr))
            avg_bonds.append(str(xmax_gofr))
            avg_coord.append(str(y_intgofr))

            xmin_clicks.append(str(xmin_click_save))
            data.append(str(xmax_gofr))
            data.append(str(ymax_gofr))
            data.append(str(xmin_gofr))
            data.append(str(y_intgofr))

    else:
        xmax_gofr = 0
        ymax_gofr = 0
        xmin_gofr = 0
        y_intgofr = 0
        xmin_click_save = 0
        max_bonds.append(str(xmin_gofr))
        peak_heights.append(str(ymax_gofr))
        avg_bonds.append(str(xmax_gofr))
        avg_coord.append(str(y_intgofr))
        xmin_clicks.append(str(xmin_click_save))
        data.append(str(xmax_gofr))
        data.append(str(ymax_gofr))
        data.append(str(xmin_gofr))
        data.append(str(y_intgofr))

    cutdistance = 0  # re initialization of variables for the next loop
    fit_min_dist = 0
    fit_max_dist = 0
    fit_min_gofr = 0
    fit_max_gofr = 0
    gofr = 0
    intgofr = 0
    indices = []
    return (data, max_bonds, peak_heights, avg_bonds, avg_coord, col, xmin_clicks)


xmin_click = []
xmax_click = []


def analyze_gofrs_interactive(data, file, couples):
    """ extraction of the min, max and coordination number using fits and interactive plot """
    distance = np.loadtxt(file, usecols=(0,), skiprows=1, unpack=True)
    col = 0
    max_bonds = []
    peak_heights = []
    avg_bonds = []
    avg_coord = []
    xmin_clicks = []
    list_atom1 = []
    list_atom2 = []
    for couple in couples:
        data, max_bonds, peak_heights, avg_bonds, avg_coord, col, xmin_clicks = interactive_fit(data, file, couples, distance, col, max_bonds, peak_heights, avg_bonds, avg_coord, xmin_clicks, couple, list_atom1, list_atom2)
    return (data, max_bonds, peak_heights, avg_bonds, avg_coord)


def main(argv):
    """     ********* Main program *********     """
    try:
        options, arg = getopt.getopt(argv, "h")
    except getopt.GetoptError:
        print('analyze_gofr.py')
        sys.exit()
    for opt, arg in options:
        if opt == '-h':
            print('*******')
            print('Analyze g(r) program to extract bond information from the gofr.dat files created by the script gofrs.py')
            print('This program should be run from the folder containing your gofr.dat files')
            print('WARNING! you have to click FIRST on the maximum, and THEN on the minimum')
            print(' ')
            print('This code produces four files: *coordination.inp, *average_bonds.inp, *maximum_bonds.inp and gofrs_summary.txt')
            print('The summary_gofrs.txt file contains the following information for each atom pair: ')
            print('xmax (x-coordinate of the first peak)')
            print('ymax (y-coordinate of the first peak)')
            print('xmin (x-coordinate of the first minimum)')
            print('')
            print('')
            print('*******')
            print(' ')
            sys.exit()
    #print('\n')
    print('********************************************************************************')
    print('Click on the position of the first maximum and then on the first minimum')
    print('OR click on "Auto" to allow the program to automatically find the positions')
    print('********************************************************************************')
    files = sorted(glob.glob('*.gofr.dat'))  # I list every gofr files in alphabetic order
    if files != []:
        f, couples = headerfile(files[0])  # I create the first newfile for gofr and save the list of element couples
        for file in files:
            data = [file]
            b, atom1, atom2 = atoms_columns_max(file)
            b_avg, atom1_avg, atom2_avg = atoms_columns_avg(file)
            b_coord, atom1_coord, atom2_coord = atoms_columns_coord(file)
            data, max_bonds, peak_heights, avg_bonds, avg_coord = analyze_gofrs_interactive(data, file, couples)
            writer = csv.writer(b, delimiter='\t')
            writer.writerows(zip(atom1, atom2, max_bonds))
            b.close()
            writer = csv.writer(b_avg, delimiter='\t')
            writer.writerows(zip(atom1_avg, atom2_avg, avg_bonds))
            b_avg.close()
            writer = csv.writer(b_coord, delimiter='\t')
            writer.writerows(zip(atom1_coord, atom2_coord, avg_coord))
            b_avg.close()
            f.write("\t".join(x for x in data) + "\n")  # we write in the file the result line
        f.close

# ********* Execution *********
if __name__ == "__main__":
    main(sys.argv[1:])



