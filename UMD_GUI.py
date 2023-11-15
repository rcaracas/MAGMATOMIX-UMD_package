#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 12:22:44 2023

@author: Kevin Jiguet-Covex
"""


import sys
import scipy
import numpy as np
import bonding_umd
import speciation_and_angles
import VaspParser2umd
import LAMMPSParser2umd
import vibr_spectrum_umd_fast
import umd_to_lammps
import gofr_umd
import Extract_umd
#import analyze_gofr_forGUI
import msd_umd
import viscosity_new
import averages_forGUI
import umd2poscar
import umd2xyz
import MSDdistrib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt5.QtWidgets import QCheckBox, QApplication, QMainWindow, QLabel, QLineEdit, QTabWidget, QPushButton, QVBoxLayout, QWidget, QFileDialog, QRadioButton, QButtonGroup, QSpacerItem,QHBoxLayout,QGridLayout,QSizePolicy,QFrame
import os
from functools import partial

def isfloat(n):
    try :
        float(n)
        return True
    except ValueError :
        return False

def select_File(LineEdit):
    file_dialog = QFileDialog()
    Line, _ = file_dialog.getOpenFileName()
    LineEdit.setText(Line)

def usefunction(funct,argv,message):#used to call a function outside this script
    try :
        res=funct.main(argv)
        if res!=None :
            return res
        else :
            return True
    except Exception as e :
        message.setText("The error : < "+str(e)+" > has occured. Please check the validity of the file and arguments you provided.")
        return False

def setGraph(self,layout):
    self.display_layout=layout
    self.display_widget.deleteLater()
    self.display_widget=QWidget()
    self.display_widget.setLayout(self.display_layout)
    self.Graph_tab.insertTab(0,self.display_widget,"Visual display")
    self.Graph_tab.setCurrentIndex(0)

def clearLayout(layout):#Clears the layout to display the next graphs
    while layout.count():
        item = layout.takeAt(0)
        if item.widget():
            item.widget().deleteLater()
        elif item.layout():
            while item.layout().count():
                item_int = item.layout().takeAt(0)
                if item_int.widget():
                    item_int.widget().deleteLater()

def create_display_layout(layout,mode="grid",ncol=3,nrow=3):    
    clearLayout(layout)
    if mode=="grid":
        layout = QGridLayout()
        layout.setContentsMargins(0,0,0,0)
        
        for i in range(ncol):            
            for j in range(1,nrow+1):
                layout.setColumnStretch(i,j)
            
    elif mode=="full":
        layout=QVBoxLayout()        
        
    return layout
    
def createMSD_graph(self,MSD,Instants,Elements,axes=None):
    n = len(Elements)
    ncol=int((n+1)/3)
    layout=create_display_layout(self.display_layout,'grid',ncol,int((n-1)/ncol))

    Title = QLabel("x : time ; y : msd (A²)")
    Title.setStyleSheet("font-size: 16px;")
    layout.addWidget(Title,0,0)
    FigList = []
    for i in range(n):
        
        layout_int = QVBoxLayout()
        row=int(i/ncol)+1
        col=i%ncol
        figure = Figure()
        figure.subplots_adjust(left=0.1,right=0.9,top=0.85,bottom=0.15)
        ax = figure.add_subplot(111)

        if axes != None :
            
            ax.plot(Instants[:],MSD[i][0][:],label = "Total MSD")

            axes = [round(axes[i],3) for i in range(len(axes))]

            for n in range(int(len(axes)/3)):
                ax.plot(Instants[:],MSD[i][i][:],label = "Axis "+str(n)+" : "+str(axes[3*n:3*(n+1)]))
        else :
            
            ax.plot(Instants[:],MSD[i][0])
            
        ax.set_title(Elements[i])
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
        ax.legend(fontsize='smaller')
        canvas = FigureCanvas(figure)
        layout_int.addWidget(canvas)
        layout.addLayout(layout_int,row,col)
        FigList.append(figure)

    SaveButton = QPushButton("Save")
    SaveButton.clicked.connect(partial(SaveFunction,FigList,Elements))
    
    layout.addWidget(SaveButton)
    
    setGraph(self,layout)


def createavg_graph(self,data,firststep,TimeStep,unit,parameter):
    layout = create_display_layout(self.display_layout,'full')
    X = [(firststep+i)*TimeStep for i in range(len(data))]
    figure=Figure()
    Title = "Mean "+parameter+" ("+unit+") "+"as a function of time (fs)"
#    Title.setStyleSheet("font-size: 16px;")
#    Title.setFixedSize(500,20)
    ax=figure.add_subplot(111)
    ax.plot(X,data)
    ax.set_title(Title)
    figure.tight_layout()
    canvas=FigureCanvas(figure)
    layout.addWidget(canvas)
    SaveButton = QPushButton("Save")
    SaveButton.clicked.connect(partial(SaveFunction,[figure],[""]))
    layout.addWidget(SaveButton)
    setGraph(self,layout)
            
def creategraph_visc(self,Tau,Visc):#Creates the viscosity graph
    layout = create_display_layout(self.display_layout,"full")
    figure=Figure(figsize=(10,10))
    Title = "Mean viscosity (Pa*s) as a function of lag time (fs)"
    ax=figure.add_subplot(111)
    ax.plot(Tau,Visc)
    ax.set_title(Title)
    figure.tight_layout()
    figure.subplots_adjust(left=0.1,right=0.9,top=0.85,bottom=0.15)
    canvas=FigureCanvas(figure)
    layout.addWidget(canvas)

    SaveButton = QPushButton("Save")
    SaveButton.clicked.connect(partial(SaveFunction,[figure],[""]))    
    layout.addWidget(SaveButton)
    setGraph(self,layout)

def read_visc(viscfile):#Reads the viscosity file
    ff=open(viscfile,"r")
    ff.readline()
    ff.readline()
    Tau=[]
    Visc=[]
    while True :
        line=ff.readline()
        if not line :
            break
        l=line.strip().split()

        Tau.append(float(l[0]))
        Visc.append(float(l[-1]))
    
    return Tau,Visc

def createhist_vib(self,DOS,Freq,Elements,types,limF = "None",wl = "None",OrigDOS = None):#Creates the histograms of vibrational spectrum
    
    #print(DOS)
    #print(len(DOS),len(DOS[0]))
    layout=create_display_layout(self.display_layout,'grid')
    if OrigDOS == None :
        OrigDOS = DOS
    
    limX = len(Freq)
    if limF != "None" :
        limX = 0    
        for i in range(len(Freq)):
            if Freq[i] >= limF :
                limX = i
                break
    
    n = len(Elements)
    natom = sum(types)
    types.append(natom)
    ncol=int((n+3.5)/3)
    #Title = QLabel("x : DOS of vib.spectrum ; y : wavenumber (cm⁻¹)")
    #Title.setStyleSheet("font-size: 16px;")
    #layout.addWidget(Title,0,0)    
    ListFig = []
    for i in range(n):
        DOSnorm = [DOS[i][j]/(3*types[i]) for j in range(len(DOS[i]))]
        OrigDOSnorm = [OrigDOS[i][j]/(3*types[i]) for j in range(len(DOS[i]))]
        layout_int = QVBoxLayout()
        row=int(i/ncol)+1
        col=i%ncol
        figure = Figure(figsize=(3,3))
        subplot = figure.add_subplot(111)
        subplot.plot(Freq[:limX],OrigDOSnorm[:limX],alpha = 0.2)
        subplot.plot(Freq[:limX], DOSnorm[:limX], label=Elements[i], color='red',linewidth=1.5)

        subplot.set_title(Elements[i])
        subplot.set_xlim(left=0)
        subplot.set_yticks([])
        subplot.set_ylabel("arb. units")
        subplot.set_xlabel("wavenumber (cm⁻¹)")
        canvas = FigureCanvas(figure)
        layout_int.addWidget(canvas)
        container = QFrame()
        container.setLayout(layout_int)
        container.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Expanding)
        layout.addWidget(container,row,col)
        ListFig.append(figure)
    
    WindowLength = QLineEdit(str(wl))
    FreqCut = QLineEdit(str(limF))    
    
    message = QLabel("")
    
    
    
    SaveButton = QPushButton("Save the plots")
    SaveButton.clicked.connect(partial(SaveFunction,ListFig,Elements))
    
    layout_buttons = QVBoxLayout()

    FreqCut.setMaximumWidth(200)
    WindowLength.setMaximumWidth(200)
    
    BoxFreq = QHBoxLayout()
    FreqLab = QLabel("Frequency cutoff :")
    FreqLab.setSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed)
    BoxFreq.addWidget(FreqLab)
    BoxFreq.addWidget(FreqCut)
    BoxFreq.setSizeConstraint(QHBoxLayout.SetFixedSize)

    
    BoxFilt = QHBoxLayout()
    FiltLab = QLabel("Window length :")
    FiltLab.setSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed)
    BoxFilt.addWidget(FiltLab)    
    BoxFilt.addWidget(WindowLength)
    BoxFilt.setSizeConstraint(QHBoxLayout.SetFixedSize)

    FilterButton = QPushButton("Display the DOS with these parameters")
    FilterButton.clicked.connect(partial(FilterFunction,self,OrigDOS,Freq,Elements,types[:-1],FreqCut,WindowLength,message))
    
    layout_buttons.addLayout(BoxFreq)
    layout_buttons.addLayout(BoxFilt)    
    layout_buttons.addWidget(FilterButton)
    layout_buttons.addWidget(SaveButton)
    layout_buttons.addWidget(message)
    
    container = QFrame()
    container.setLayout(layout_buttons)
    container.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Expanding)

    
    layout.addWidget(container,n%ncol+1,ncol-1)

    setGraph(self,layout)

def FilterFunction(self,DOS,Freq,Elements,types,Fcut,WindowsLength,message):
   
    DOSfiltered=[]
    flagCalc = True
    WL = WindowsLength.text()
    if WL.isdigit() or WL=="None":
        if WL != "None":
            if int(WL)>3:
                if int(WL)<=len(Freq):
                    for subDOS in DOS:
                        DOSfiltered.append(scipy.signal.savgol_filter(np.array(subDOS),int(WL),3))

                else:
                    message.setText("Error : the Window Length was set to "+str(WL)+", but must not be bigger than "+str(len(Freq))+" (data's length)")
                    WindowsLength.setText("None")
                    flagCalc=False
            else : 
                message.setText("Error : the Window Length was set to "+str(WL)+", but must be at least equal to 3.")
                WindowsLength.setText("None")
                flagCalc = False
                
        elif Fcut.text()!="None" :
            DOSfiltered = DOS
        else:
            flagCalc=False
            
        if flagCalc and isfloat(Fcut.text()):
            createhist_vib(self,DOSfiltered,Freq,Elements,types,float(Fcut.text()),WL,DOS)
        elif flagCalc :
            Fcut.setText("None")
            createhist_vib(self,DOSfiltered,Freq,Elements,types,"None",WL,DOS)
        
    else :
        
        message.setText("Error : the window length must be an integer>3 or 'None'.")        


def displayVar(self, Val,X,title="",xvar="",yvar="",yunit="",timestep=1):
    
    layout=create_display_layout(self.display_layout,'full')
    figure = Figure()
    ax = figure.add_subplot()
    
    ax.plot([timestep*t for t in X],Val)
#    ax.legend()
    ax.set_xlim(left=X[0])
    ax.set_xlabel("time ("+xvar+")")
    ax.set_ylabel(yvar+" ("+yunit+")")
#    ax.set_ylim(bottom=0)
    Title = QLabel(title)
    Title.setStyleSheet("font-size: 16px")
    Title.setFixedSize(500,20)

    canvas = FigureCanvas(figure)
    
    layout.addWidget(Title)
    layout.addWidget(canvas)
    setGraph(self,layout)
        
    
def displaydistrib(self,Thetas,MSDsEl,Elements,Axes):
    A1 = Axes[:3]
    A2 = Axes[3:]
    
    layout=create_display_layout(self.display_layout,"full")
    
    figure = Figure()
    ax = figure.add_subplot(111)
    for el in range(len(Elements)):
        ax.plot(Thetas,MSDsEl[el],label = Elements[el])
        
    ax.legend()
    ax.set_xlabel("Angular fraction")
    ax.set_ylabel("MSD")
    Title = QLabel("Axes : "+str(A1)+", "+str(A2))
    Title.setStyleSheet("font-size: 16px")
    Title.setFixedSize(500,20)
    
    canvas = FigureCanvas(figure)
    
    layout.addWidget(Title)
    layout.addWidget(canvas)
    setGraph(self,layout)    

def SaveFunction(ListFig,ListTitles):

    nameWidget = QWidget()
    saveLayout = QVBoxLayout()
    Name = QLineEdit("")
    saveForRealButton = QPushButton("Save")

    def save():
        name = Name.text()
        for i in range(len(ListFig)) :
            if ListTitles[i]!="":
                ListFig[i].savefig(name+"_"+ListTitles[i], format='png')
            else:
                ListFig[i].savefig(name, format='png')
            
        nameWidget.setVisible(False)
    
    saveForRealButton.clicked.connect(save)
    
    saveLayout.addWidget(QLabel("Enter the base for the name(s) of the file(s) to be saved :"))
    saveLayout.addWidget(Name)
    saveLayout.addWidget(saveForRealButton)
    
    nameWidget.setLayout(saveLayout)
    nameWidget.setVisible(True)


        
def read_vibr(vibfile,ncol):
    ff=open(vibfile,'r')
    DOS=[]
    Freq=[]
    ff.readline()
    while True :
        line = ff.readline().strip().split()
        if not line :
            break
        DOS.append(float(line[ncol]))
        Freq.append(float(line[0]))
    ff.close()
    return DOS,Freq

def read_vibr_allEls(vibfile):#Reads all the elements of the vibr file to create the histograms
    ff=open(vibfile,'r')
    Elements=[]
    Freq=[]
    line = ff.readline().strip().split()
    for label in line[1:-1]:
        el = label.split("_")[-1]
        Elements.append(el)
    Elements.append("Total DOS")
    print(Elements)
    
    DOS=[[] for _ in range(len(Elements))]

    while True :
        line = ff.readline().strip().split()
        if not line :
            break
        for i in range((len(Elements))):
            DOS[i].append(float(line[i+1]))
        Freq.append(float(line[0]))
    ff.close()
    print(len(Freq))
    print(len(DOS[0]))
    return DOS,Freq,Elements



class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Magmatix : a tool for the use of raw VASP and LAMMPS data files")
        self.setMinimumSize(0,0)
        self.resize(600,600)
        
        
        mainWidget = QWidget()        
        mainLayout = QVBoxLayout()
        main_tab = QTabWidget()
        mainLayout.addWidget(main_tab)
        mainWidget.setLayout(mainLayout)

        # Création du widget à onglets
        
        convert_widget=QWidget()
        convert_layout=QVBoxLayout()
        convert_tab=QTabWidget()
        convert_layout.addWidget(convert_tab)
        convert_widget.setLayout(convert_layout)
        
        self.Graph_widget=QWidget()
        self.Graph_layout=QVBoxLayout()
        self.Graph_tab=QTabWidget()
        self.Graph_layout.addWidget(self.Graph_tab)
        self.Graph_widget.setLayout(self.Graph_layout)
                
        operations_widget=QWidget()
        operations_layout=QVBoxLayout()
        operations_tab=QTabWidget()
        operations_layout.addWidget(operations_tab)
        operations_widget.setLayout(operations_layout)
        
        calculations_widget=QWidget()
        calculations_layout=QVBoxLayout()
        calculations_tab=QTabWidget()
        calculations_layout.addWidget(calculations_tab)
        calculations_widget.setLayout(calculations_layout)

        # Création des onglets
        analyze_gofr_widget = QWidget()
        UMD_widget = QWidget()
        Bond_widget = QWidget()
        Speciation_widget = QWidget()
        LAMMPSParser_widget = QWidget()
        Vibration_widget = QWidget()
        viscosity_widget = QWidget()        
        msd_widget = QWidget()
        UMDtoLAMMPS_widget = QWidget()
        gofr_widget = QWidget()
        msd_widget = QWidget()        
        avg_widget = QWidget()
        poscar_widget = QWidget()
        xyz_widget = QWidget()
        msd_dist_widget = QWidget()
        displayVar_widget = QWidget()
        self.display_widget=QWidget()
        


        # Ajout des onglets au widget à onglets
        convert_tab.addTab(UMD_widget, "VASP to UMD")
        convert_tab.addTab(LAMMPSParser_widget, "LAMMPS to UMD")
        convert_tab.addTab(UMDtoLAMMPS_widget, "UMD to LAMMPS")
        convert_tab.addTab(poscar_widget, "UMD to POSCAR")
        convert_tab.addTab(xyz_widget, "UMD to XYZ")
        operations_tab.addTab(gofr_widget,"Gofr file")
#        operations_tab.addTab(analyze_gofr_widget, "Average bonds file!! Not working yet")
        operations_tab.addTab(Bond_widget, "Bonding")
        operations_tab.addTab(Speciation_widget, "Population")
        calculations_tab.addTab(Vibration_widget,"Vibrational spectrum")
        calculations_tab.addTab(viscosity_widget,"Viscosity parameters")
        calculations_tab.addTab(msd_widget,"MSD")
        calculations_tab.addTab(avg_widget, "Average parameters")
        calculations_tab.addTab(displayVar_widget, "Parameters of UMD")
        calculations_tab.addTab(msd_dist_widget, "MSD distribution")
        self.Graph_tab.addTab(self.display_widget,"Visual Display")
        self.Graph_tab.addTab(displayVar_widget,"UMD variables")

        # Création du contenu des onglets
        self.create_Vibration_layout(Vibration_widget)
        self.create_avg_layout(avg_widget)
        self.create_UMD_layout(UMD_widget)
        self.create_LAMMPS_layout(LAMMPSParser_widget)
        self.create_Bond_layout(Bond_widget)
        self.create_Speciation_layout(Speciation_widget)
        self.create_UMD_to_LAMMPS_layout(UMDtoLAMMPS_widget)
        self.create_gofr_layout(gofr_widget)
        self.create_msd_layout(msd_widget)
        self.create_viscosity_layout(viscosity_widget)
        self.create_analyze_gofr_layout(analyze_gofr_widget)
        self.create_poscar(poscar_widget)
        self.create_xyz(xyz_widget)
        self.create_msddistrib(msd_dist_widget)
        self.create_displayVar(displayVar_widget)
        self.display_layout=QVBoxLayout()
        create_display_layout(self.display_layout)
        self.display_widget.setLayout(self.display_layout)

        # Définition du widget à onglets comme widget central de la fenêtre principale
        main_tab.addTab(convert_widget,"Conversions")
        main_tab.addTab(operations_widget,"Operations : file creation")
        main_tab.addTab(calculations_widget,"Calculations")
        main_tab.addTab(self.Graph_widget,"Graphics")
        self.setCentralWidget(main_tab)



    def create_msddistrib(self,tab_widget):
        layout = QVBoxLayout()
        
        ComputeButton = QPushButton("Compute")
        ComputeButton.clicked.connect(self.distrib_msd)
        
        self.MSDdistFile = QLineEdit("")

        SelectButton = QPushButton("Select an umd file")
        SelectButton.clicked.connect(partial(select_File,self.MSDdistFile))
        
        self.MSDdistCoresBox = QCheckBox("Precise the number of cores to be used for the parallelization")        
        self.MSDdistCores = QLineEdit("")
        self.MSDdistCores.setMaximumWidth(50)
        self.MSDdistCores.setVisible(False)
        self.MSDdistCoresBox.stateChanged.connect(self.MSDdistCoresBoxchange)

        self.nVecs = QLineEdit("20")
        
        self.D00 = QLineEdit("1")
        self.D01 = QLineEdit("0")
        self.D02 = QLineEdit("0")
        self.D10 = QLineEdit("0")
        self.D11 = QLineEdit("1")
        self.D12 = QLineEdit("0")

        self.MSDdist_message = QLabel("")        

        hBoxLab = QHBoxLayout()
        hBoxAx0 = QHBoxLayout()
        hBoxAx1 = QHBoxLayout()
        vBoxAxExp = QVBoxLayout()

        
        hBoxLab.addItem(QSpacerItem(50,0))
        hBoxLab.addWidget(QLabel("x"))
        hBoxLab.addWidget(QLabel("y"))
        hBoxLab.addWidget(QLabel("z"))
        hBoxAx0.addWidget(QLabel("Axis 1 :"))
        hBoxAx0.addWidget(self.D00)
        hBoxAx0.addWidget(self.D01)
        hBoxAx0.addWidget(self.D02)
        hBoxAx1.addWidget(QLabel("Axis 2 :"))
        hBoxAx1.addWidget(self.D10)
        hBoxAx1.addWidget(self.D11)
        hBoxAx1.addWidget(self.D12)

        
        vBoxAxExp.addLayout(hBoxLab)
#        vBoxAxExp.addItem(QSpacerItem(0,-75))
        vBoxAxExp.addLayout(hBoxAx0)
        vBoxAxExp.addLayout(hBoxAx1)

        hBoxnVecs = QHBoxLayout()
        hBoxnVecs.addWidget(QLabel("Enter the number of vectors to sample : "))
        hBoxnVecs.addWidget(self.nVecs)

        layout.addWidget(self.MSDdistFile)
        layout.addWidget(SelectButton)
        layout.addItem(QSpacerItem(0,75))
        layout.addLayout(hBoxnVecs)
        layout.addItem(QSpacerItem(0,75))
        layout.addLayout(vBoxAxExp)
        layout.addItem(QSpacerItem(0,75))
        layout.addWidget(self.MSDdistCoresBox)
        layout.addWidget(self.MSDdistCores)
        layout.addWidget(ComputeButton)    
        layout.addWidget(self.MSDdist_message)
        
        tab_widget.setLayout(layout)

    def distrib_msd(self):
        
        if os.path.isfile(self.MSDdistFile.text()) :
            if self.nVecs.text().isnumeric():
                if isfloat(self.D00.text()) and isfloat(self.D01.text()) and isfloat(self.D02.text()) and isfloat(self.D10.text()) and isfloat(self.D11.text()) and isfloat(self.D12.text()) :            
                    if(self.MSDdistCoresBox.isChecked() and self.MSDdistCores.text().isnumeric()) or (not self.MSDCoresBox.isChecked()):
                    
                        Axes = [float(self.D00.text()),float(self.D01.text()),float(self.D02.text()),float(self.D10.text()),float(self.D11.text()),float(self.D12.text())]
                        argv = ["-f",self.MSDdistFile.text(),"-n",self.nVecs.text(),"-a",str(Axes)]
                        if self.MSDdistCoresBox.isChecked() :
                            argv+=["-k",self.MSDdistCores.text()]
                            
                        res = usefunction(MSDdistrib, argv, self.MSDdist_message)
                    
                        if res != False :
                            Theta,MSDEl,Elements,file = res
                            displaydistrib(self,Theta,MSDEl,Elements,Axes)
                            self.MSDdist_message.settext("distribution of MSD successfully calculated and graphically displayed.\nData stored in file : "+file)
                    else :
                        self.MSDdist_message.setText("Error : the number of cores has to be a strictly positive integer.")
                else :
                    self.MSDdist_message.setText("Error : all the values for Axis 1 and 2 must be integer or float. Please check the values you provided.")
            else :
                self.MSDdist_message.setText("Error : please select an integer value for the number of vectors.")
        else :
            self.MSDdist_message.setText("Error : file "+self.MSDdistFile.text()+" is displaced or missing.")

    def MSDdistCoresBoxchange(self,state):
        if state ==2 :
            self.MSDdistCores.setVisible(True)
        else :
            self.MSDdistCores.setVisible(False)
            
    
    def create_poscar(self,tab_widget):
        layout = QVBoxLayout()
        
        self.U2poscar_edit = QLineEdit("")
        self.sampFreqPoscar = QLineEdit("1")
        
        selectUMDButton = QPushButton("Select the UMD file")
        selectUMDButton.clicked.connect(partial(select_File,self.U2poscar_edit))
        submit_button = QPushButton("Convert")
        submit_button.clicked.connect(self.convert_U2poscar)
        
        self.U2poscarMessage = QLabel("")
        
        layout.addItem(QSpacerItem(0,50))
        layout.addWidget(self.U2poscar_edit)
        layout.addWidget(selectUMDButton)
        layout.addItem(QSpacerItem(0,50))
        layout.addWidget(QLabel("Select the sampling frequency :"))
        layout.addWidget(self.sampFreqPoscar)
        layout.addItem(QSpacerItem(0,50))
        layout.addWidget(submit_button)
        layout.addWidget(self.U2poscarMessage)
        
        tab_widget.setLayout(layout)

    def convert_U2poscar(self):
        if os.path.isfile(self.U2poscar_edit.text()):
            if (self.sampFreqPoscar.text()).isnumeric() and self.sampFreqPoscar.text()!="0":
                argv = ["-f",self.U2poscar_edit.text(),"-s",self.sampFreqPoscar.text()]
                res = usefunction(umd2poscar, argv, self.U2poscarMessage)
                if res == True :
                    Name = self.U2poscar_edit.text()[:-7] + 'xyz'
                    self.U2poscarMessage.setText("The file "+Name+" has been successfully created.")
                    
            else :
                self.U2poscarMessage.setText("Error : the sampling frequency must be a strictly positive integer.")
        else :
            self.U2poscarMessage.setText("Error : the file "+self.U2poscar_edit.text()+" is displaced or missing")
            
    def create_xyz(self,tab_widget):
        layout = QVBoxLayout()
        
        self.U2xyz_edit = QLineEdit("")
        self.sampFreqxyz = QLineEdit("1")
        
        selectUMDButton = QPushButton("Select the UMD file")
        selectUMDButton.clicked.connect(partial(select_File,self.U2xyz_edit))
        submit_button = QPushButton("Convert")
        submit_button.clicked.connect(self.convert_U2xyz)
        
        self.U2xyzMessage = QLabel("")
        
        layout.addItem(QSpacerItem(0,50))
        layout.addWidget(self.U2xyz_edit)
        layout.addWidget(selectUMDButton)
        layout.addItem(QSpacerItem(0,50))
        layout.addWidget(QLabel("Select the sampling frequency :"))
        layout.addWidget(self.sampFreqxyz)
        layout.addItem(QSpacerItem(0,50))
        layout.addWidget(submit_button)
        layout.addWidget(self.U2xyzMessage)
        
        tab_widget.setLayout(layout)

    def convert_U2xyz(self):
        if os.path.isfile(self.U2xyz_edit.text()):
            if (self.sampFreqxyz.text()).isnumeric() and self.sampFreqxyz.text()!="0":
                argv = ["-f",self.U2xyz_edit.text(),"-s",self.sampFreqxyz.text()]
                res = usefunction(umd2xyz, argv, self.U2xyzMessage)
                if res == True :
                    Name = self.U2xyz_edit.text()[:-7] + 'xyz'
                    self.U2xyzMessage.setText("The file "+Name+" has been successfully created.")

            else :
                self.U2xyzMessage.setText("Error : the sampling frequency must be a strictly positive integer.")
        else :
            self.U2xyzMessage.setText("Error : the file "+self.U2xyz_edit.text()+" is displaced or missing")

    def create_analyze_gofr_layout(self,tab_widget):
        layout = QVBoxLayout()
        self.gofrFile_edit = QLineEdit("")
        self.anagofr_message = QLabel("")
        SButton = QPushButton("Select the .gofr file")
        SButton.clicked.connect(partial(partial(select_File,self.gofrFile_edit)))
        
        CButton = QPushButton("Compute")
        CButton.clicked.connect(self.analyze_gofr)
        
        layout.addWidget(self.gofrFile_edit)
        layout.addWidget(SButton)
        layout.addItem(QSpacerItem(0,150))
        layout.addWidget(CButton)
        layout.addWidget(self.anagofr_message)

        
        tab_widget.setLayout(layout)

    def analyze_gofr(self):
        if not os.path.isfile(self.gofrFile_edit.text()):
               self.anagofr_message.setText("Error : the UMD file "+self.gofrFile_edit.text()+" is displaced or missing. Pleace check the path or the file name.")
        else :
            argv=["-f",self.gofrFile_edit.text()]
#            testPlot.main()
#            res = analyze_gofr_forGUI.main(argv)
#            if res == True :
#                name = self.gofrFile.text()[:-9]+".average_bonds.dat"
#                print("Average coordination distance successfully calculaterd. The results are stored in the file "+name)

    def create_avg_layout(self,tab_widget):
        layout=QVBoxLayout()
        self.avgUMD_edit = QLineEdit("")
        self.select_parameter = QLineEdit("")
        self.avg_message=QLabel("")
        self.select_init=QLineEdit("0")
        self.select_last=QLineEdit("max")
        Hlayout = QHBoxLayout()
        UMDButton = QPushButton("Select the UMD file")
        ComputeButton = QPushButton("Compute and display")
        UMDButton.clicked.connect(partial(select_File,self.avgUMD_edit))                
        ComputeButton.clicked.connect(self.compute_avg)
        
        Hlayout.addWidget(self.select_init)
        Hlayout.addItem(QSpacerItem(150,0))
        Hlayout.addWidget(self.select_last)
        
        
        layout.addWidget(self.avgUMD_edit)
        layout.addWidget(UMDButton)
        layout.addWidget(QLabel("Enter the parameter of which the average shall be displayed :"))
        layout.addWidget(self.select_parameter)
        layout.addWidget(QLabel("Indexes of the first and last snapshot to be displayed :"))
        layout.addLayout(Hlayout)
        layout.addItem(QSpacerItem(0,150))
        layout.addWidget(ComputeButton)
        layout.addWidget(self.avg_message)
        
        tab_widget.setLayout(layout)
        
    def compute_avg(self):
        argv=["-f",self.avgUMD_edit.text(),"-p",self.select_parameter.text()]
        flag=1
        if self.select_init.text().isnumeric():
            argv+=["-s",self.select_init.text()]
        else :
            self.avg_message.setText("The initial step must be specified and set to a strictly positive integer value.")
            flag=0
            
        if self.select_last.text() != "max":
            if isfloat(self.select_last.text()):
                argv+=["-m",self.select_last.text()]
            else:
                self.select_last.setText("max")
            
        if (not self.select_last.text().isnumeric()) and self.select_last.text()!="max":
            self.avg_message.setText("The last step must be specified and set to a strictly positive integer value (or keyword max).")
            flag=0
        elif flag and self.select_last.text()!="max" and int(self.select_init.text())>int(self.select_last.text()):
            self.avg_message.setText("The initial step must smaller than the last step.")
            flag=0

        if flag :
            if os.path.isfile(self.avgUMD_edit.text()):
                
                res = usefunction(averages_forGUI,argv,self.avg_message)
                if res!=False :
                    data,average,variance,stdev,TimeStep,unit=res
                    if self.select_last.text() == "max": 
                        self.select_last.setText(str(len(data)+int(self.select_init.text())))
                    createavg_graph(self,data[int(self.select_init.text()):int(self.select_last.text())],int(self.select_init.text()),TimeStep,unit,self.select_parameter.text())
                    self.avg_message.setText("Average "+self.select_parameter.text()+" calculated. Evolution graphically displayed. \nAverage "+self.select_parameter.text()+" : "+str(average)+"\nVariance : "+str(variance)+"\nStandard deviation : "+str(stdev))
                    
            else :
               self.avg_message.setText("Error : the UMD file "+self.avgUMD_edit.text()+" is displaced or missing. Pleace check the path or the file name.")

    def create_displayVar(self,tab_widget):
        layout = QVBoxLayout()
        self.umdfile_disp = QLineEdit("")
        selectButton=QPushButton("Select file")
        self.TempButton=QRadioButton("Temperature")
        self.EnButton=QRadioButton("Internal Energy")
        self.PressButton=QRadioButton("Pressure")
        buttongroup=QButtonGroup()
        computeButton = QPushButton("Display")
        self.displaymessage=QLabel("")
        
        selectButton.clicked.connect(partial(select_File,self.umdfile_disp))
        buttongroup.addButton(self.TempButton)
        buttongroup.addButton(self.EnButton)
        buttongroup.addButton(self.PressButton)
        computeButton.clicked.connect(self.display_var)
        
                
        
        layout.addWidget(self.umdfile_disp)
        layout.addWidget(selectButton)
        layout.addWidget(self.TempButton)
        layout.addWidget(self.EnButton)
        layout.addWidget(self.PressButton)
        layout.addWidget(computeButton)
        layout.addWidget(self.displaymessage)
        
        tab_widget.setLayout(layout)
    

    def display_var(self):
        if not os.path.isfile(self.umdfile_disp.text()):
            if self.umdfile_disp.text()=="":
                self.displaymessage.setText("Please select a file.")
            else:
                self.displaymessage.setText("Error : the file "+self.umdfile_disp.text()+" is displaced or missing.")
        else:
            argv=["-f",self.umdfile_disp.text()]
            arg=None
            if self.TempButton.isChecked():
                arg="Temperature"
            elif self.PressButton.isChecked():
                arg="Pressure"
            elif self.EnButton.isChecked():
                arg="InternalEnergy"
            else:
                self.displaymessage.setText("Please select a parameter.")
           
            argv+=["-p",arg]
            if arg!=None:
                res=usefunction(Extract_umd, argv, self.displaymessage)        
                if res !=False:
#                self.create_display_layout()
                    [Var,Times,unit,tunit,timestep]=res
                    displayVar(self,Var,Times,xvar=tunit,yvar=arg,yunit=unit,timestep=timestep)
                    
                    self.displaymessage.setText(arg+" successfully extracted and graphically shown (see the tab 'Visual Display'")

    def create_viscosity_layout(self,tab_widget):
        layout = QVBoxLayout()
        
        self.umdfile_vsc_edit = QLineEdit("")
        selectButton = QPushButton("Select the UMD file")
        computeButton = QPushButton("Compute")
        self.i_vsc = QLineEdit("0")
        self.n_vsc = QLineEdit("2000")
        
        selectButton.clicked.connect(partial(select_File,self.umdfile_vsc_edit))
        computeButton.clicked.connect(self.compute_vsc)
        self.vscmessage=QLabel("")
        
        self.displayVisc=QCheckBox("Graphically Display the mean viscosity after calculation")
        self.displayVisc.setChecked(True)
        
        hboxIS = QHBoxLayout()
        hboxIS.addWidget(QLabel("Enter the index of the initial step :"))
        hboxIS.addWidget(self.i_vsc)
        Precision = QLabel("This value shall not be greater than the total duration of the simulation.")
        Precision.setFixedSize(500,20)
        
        
        hboxLA = QHBoxLayout()
        hboxLA.addWidget(QLabel("Enter the value of lag for autocorrelation in fs (default 2000) :"))
        hboxLA.addWidget(self.n_vsc)
        layout.addWidget(self.umdfile_vsc_edit)
        layout.addWidget(selectButton)
        layout.addItem(QSpacerItem(0,80))
        layout.addLayout(hboxIS)
        layout.addItem(QSpacerItem(0,80))
        layout.addLayout(hboxLA)
        layout.addWidget(Precision)
        layout.addItem(QSpacerItem(0,80))
        layout.addWidget(computeButton)
        layout.addWidget(self.displayVisc)
        layout.addWidget(self.vscmessage)
        
        tab_widget.setLayout(layout)
        
        
    def compute_vsc(self):
        if not os.path.isfile(self.umdfile_vsc_edit.text()):
            self.vscmessage.setText("Error : the UMD file "+self.umdfile_vsc_edit.text()+" is displaced or missing. Pleace check the path or the file name.")
        elif not self.i_vsc.text().isnumeric() or int(self.n_vsc.text())<0:
            self.vscmessage.setText("Error : the index of the initial step must be a positive integer.")
        elif not self.n_vsc.text().isnumeric() or int(self.n_vsc.text())<0:
            self.vscmessage.setText("Error : the vertical step must be a positive integer.")
        else :
            argv = ["-f",self.umdfile_vsc_edit.text(),"-i",self.i_vsc.text(),"-n",self.n_vsc.text()]
            result = usefunction(viscosity_new,argv,self.vscmessage)
            if result==True and self.displayVisc.isChecked():
                Name = self.umdfile_vsc_edit.text().split("/")
                name = Name[-1][:-8]+".visc.dat"
                Tab,Visc=read_visc(name)
                creategraph_visc(self,Tab,Visc)
                self.vscmessage.setText("Viscosity file successfully created under the name "+name+"\nViscosity showed in graphic form (see the tab <Visual Display> above).")
            elif result :
                Name = self.umdfile_vsc_edit.text().split("/")
                name = Name[-1][:-8]+".visc.dat"
                self.vscmessage.setText("Viscosity file successfully created under the name "+name)

        
        
    def create_msd_layout(self,tab_widget):
        layout = QVBoxLayout()
        
        self.umdfile_msd_edit = QLineEdit("")
        selectButton = QPushButton("Select the UMD file")
        computeButton = QPushButton("Compute")
        self.z_msd = QLineEdit("1")
        self.v_msd = QLineEdit("1")
        self.b_msd = QLineEdit("0")
        
        
        self.MSDCoresBox = QCheckBox("Precise the number of cores to be used for the parallelization")        
        self.MSDCores = QLineEdit("")
        self.MSDCores.setMaximumWidth(50)
        self.MSDCores.setVisible(False)
        self.MSDCoresBox.stateChanged.connect(self.MSDCoresBoxchange)
        
        selectButton.clicked.connect(partial(select_File,self.umdfile_msd_edit))
        computeButton.clicked.connect(self.compute_msd)
        self.msdmessage=QLabel("")
        self.atoms = QRadioButton("msd of individual atoms")
        self.elements = QRadioButton("msd of each element")
        buttongroup = QButtonGroup()
        buttongroup.addButton(self.atoms)
        buttongroup.addButton(self.elements)
        self.atoms.setChecked(True)

        self.elements.toggled.connect(self.dispVisual)        

        hboxRB = QHBoxLayout()
        hboxRB.addWidget(self.elements)
        hboxRB.addItem(QSpacerItem(-380,0))
        hboxRB.addWidget(self.atoms)


        hboxLabels = QHBoxLayout()
        hBoxEdits = QHBoxLayout()

        hboxLabels.addWidget(QLabel("Enter the value of the horizontal jump :"))
        hboxLabels.addWidget(QLabel("Enter the value of the vertical jump :"))
        hboxLabels.addWidget(QLabel("Enter the lenght of the ballistic regime :"))

        hBoxEdits.addWidget(self.z_msd)
#        hBoxEdits.addItem(QSpacerItem(150,0))        
        hBoxEdits.addWidget(self.v_msd)
#        hBoxEdits.addItem(QSpacerItem(150,0))                
        hBoxEdits.addWidget(self.b_msd)        

        hBoxAx = QHBoxLayout()
        vBoxAxExp = QVBoxLayout()
        self.ExpWidget = QWidget()
        hBoxLab = QHBoxLayout()
        hBoxAx0 = QHBoxLayout()
        hBoxAx1 = QHBoxLayout()
        hBoxAx2 = QHBoxLayout()
        
        self.Auto = QCheckBox("Auto")
        self.Default = QCheckBox("Default")
        self.Custom = QCheckBox("Custom")
        self.Auto.stateChanged.connect(lambda : self.changeAxes(self.Auto))
        self.Default.stateChanged.connect(lambda : self.changeAxes(self.Default))
        self.Custom.stateChanged.connect(self.CustomBox)
        self.Default.setChecked(True)
        
        self.M00 = QLineEdit('1')
        self.M01 = QLineEdit('0')
        self.M02 = QLineEdit('0')
        self.M10 = QLineEdit('0')
        self.M11 = QLineEdit('1')
        self.M12 = QLineEdit('0')
        self.M20 = QLineEdit('0')
        self.M21 = QLineEdit('0')
        self.M22 = QLineEdit('1')
        
        
        hBoxLab.addWidget(QLabel(""))
        hBoxLab.addWidget(QLabel("x"))
        hBoxLab.addWidget(QLabel("y"))
        hBoxLab.addWidget(QLabel("z"))
        hBoxAx0.addWidget(QLabel("Axis 1 :"))
        hBoxAx0.addWidget(self.M00)
        hBoxAx0.addWidget(self.M01)
        hBoxAx0.addWidget(self.M02)
        hBoxAx1.addWidget(QLabel("Axis 2 :"))
        hBoxAx1.addWidget(self.M10)
        hBoxAx1.addWidget(self.M11)
        hBoxAx1.addWidget(self.M12)
        hBoxAx2.addWidget(QLabel("Axis 3 :"))
        hBoxAx2.addWidget(self.M20)
        hBoxAx2.addWidget(self.M21)
        hBoxAx2.addWidget(self.M22)
        
 
        IntLayout = QHBoxLayout()
        IntLayout.addWidget(QLabel("Enter the explicit expression for the vector axes : "))
        IntLayout.addLayout(vBoxAxExp)
        
        hBoxAx.addWidget(self.Default)
        hBoxAx.addWidget(self.Auto)
        hBoxAx.addWidget(self.Custom)
        vBoxAxExp.addLayout(hBoxLab)
        vBoxAxExp.addLayout(hBoxAx0)
        vBoxAxExp.addLayout(hBoxAx1)
        vBoxAxExp.addLayout(hBoxAx2)
        self.ExpWidget.setLayout(IntLayout)
        hBoxAx.addWidget(self.ExpWidget)
        self.ExpWidget.setVisible(False)

        self.displayMSDbox = QCheckBox("Display the MSD of each element as a graph")
        self.displayMSDbox.setVisible(False)

        layout.addWidget(self.umdfile_msd_edit)
        layout.addWidget(selectButton)
        layout.addItem(QSpacerItem(0,20))
        layout.addLayout(hboxLabels)
        layout.addItem(QSpacerItem(0,-15))
        layout.addLayout(hBoxEdits)
        layout.addItem(QSpacerItem(0,0))
        layout.addWidget(QLabel("Choose the form of the data to be displayed in the output file :"))
        layout.addItem(QSpacerItem(0,-30))
        layout.addLayout(hboxRB)
        layout.addWidget(self.MSDCoresBox)
        layout.addWidget(self.MSDCores)
        layout.addWidget(self.displayMSDbox)
        layout.addItem(QSpacerItem(0,50))
        layout.addWidget(QLabel("Choose the way the msd should be computed :"))
        layout.addLayout(hBoxAx)
        layout.addWidget(computeButton)
        layout.addItem(QSpacerItem(0,-20))
        layout.addWidget(self.msdmessage)
        
        tab_widget.setLayout(layout)
    
    def CustomBox(self):
        state = self.Custom.isChecked()
        if state :
            self.ExpWidget.setVisible(True)
            self.changeAxes(self.Custom)
        else :
            self.ExpWidget.setVisible(False)
    
    def changeAxes(self,checkbox):
        state = checkbox.isChecked()
        if state :
            if checkbox != self.Default :
                self.Default.setChecked(False)
            if checkbox != self.Auto :
                self.Auto.setChecked(False)
            if checkbox != self.Custom :                   
                self.Custom.setChecked(False)
        
    def dispVisual(self):
        state = self.sender().isChecked()
        if state :
            self.displayMSDbox.setVisible(True)
        else :
            self.displayMSDbox.setVisible(False)
            
    def MSDCoresBoxchange(self,state):
        if state ==2 :
            self.MSDCores.setVisible(True)
        else :
            self.MSDCores.setVisible(False)
        
    def compute_msd(self):
        
        if self.MSDCoresBox.isChecked() and not self.MSDCores.text().isnumeric():
            self.msdmessage.setText("Error : the number of cores has to be a strictly positive integer")
        elif not os.path.isfile(self.umdfile_msd_edit.text()):
            self.msdmessage.setText("Error : the UMD file "+self.umdfile_msd_edit.text()+" is displaced or missing. Pleace check the path or the file name.")
        elif not self.z_msd.text().isnumeric() or int(self.z_msd.text())<1:
            self.msdmessage.setText("Error : the horizontal jump must be a strictly positive integer.")
        elif not self.v_msd.text().isnumeric() or int(self.v_msd.text())<1:
            self.msdmessage.setText("Error : the vertical step must be a strictly positive integer.")
        elif not self.b_msd.text().isnumeric() or int(self.b_msd.text())<0:
            self.msdmessage.setText("Error : the length of the ballistic regime must be a positive integer.")
        else :
            Comp = True
            if self.atoms.isChecked():
                argv = ["-f",self.umdfile_msd_edit.text(),"-z",self.z_msd.text(),"-v",self.v_msd.text(),"-b",self.b_msd.text(),"-m","atoms"]
            else :
                argv = ["-f",self.umdfile_msd_edit.text(),"-z",self.z_msd.text(),"-v",self.v_msd.text(),"-b",self.b_msd.text(),"-m","elements"]
            
            if self.Auto.isChecked():
                argv += ["-a","Auto"]
            elif self.Custom.isChecked():
                if isfloat(self.M00.text()) and isfloat(self.M01.text()) and isfloat(self.M02.text()) and isfloat(self.M10.text()) and isfloat(self.M11.text()) and isfloat(self.M12.text()) and isfloat(self.M20.text()) and isfloat(self.M21.text()) and isfloat(self.M22.text()) :
                    axes = "["+str(self.M00.text())+","+str(self.M01.text())+","+str(self.M02.text())+","+str(self.M10.text())+","+str(self.M11.text())+","+str(self.M12.text())+","+str(self.M20.text())+","+str(self.M21.text())+","+str(self.M22.text())+"]"
                    argv+=["-a",axes]
                else :
                    self.msdmessage.setText("Error : each component of the axis vectors must be a float. Please check the values you provided.")
                    Comp = False
                    
            if self.MSDCoresBox.isChecked():
                argv+=["-k",self.MSDCores.text()]
                
            if Comp :
                
                result = usefunction(msd_umd,argv,self.msdmessage)
                if result!=False :
                    [Name,MSD,Instants,Elements,Axes] = result
                    message = "MSD file successfully created under the name : "+Name
                    if self.displayMSDbox.isChecked():
                        self.display_layout=create_display_layout(self.display_layout) 
                        createMSD_graph(self,MSD, Instants,Elements, axes = Axes)
                        message += '\nMSD displayed as graphs in the tab "Graph" -> "Visual display"'
                    self.msdmessage.setText(message)        

    def create_gofr_layout(self,tab_widget):
        layout = QVBoxLayout()
        
        self.umdfileg_edit = QLineEdit("")
        selectButton = QPushButton("Select the UMD file")
        computeButton = QPushButton("Compute")
        self.s_g = QLineEdit("1")
        self.d_g = QLineEdit("0.01")
        self.i_g = QLineEdit("0")


        self.gofrCoresBox = QCheckBox("Precise the number of cores to be used for the parallelization")        
        self.gofrCores = QLineEdit("")
        self.gofrCores.setMaximumWidth(50)
        self.gofrCores.setVisible(False)
        self.gofrCoresBox.stateChanged.connect(self.gofrCoresBoxchange)

        
        selectButton.clicked.connect(partial(select_File,self.umdfileg_edit))
        computeButton.clicked.connect(self.compute_gofr)
        self.gofrmessage=QLabel("")
        
        hboxSF = QHBoxLayout()
        hboxSF.addWidget(QLabel("Enter the sampling frequency :"))
        hboxSF.addWidget(self.s_g)
        
        hboxDI = QHBoxLayout()
        hboxDI.addWidget(QLabel("Enter the discretization interval :"))
        hboxDI.addWidget(self.d_g)
        
        hboxIS = QHBoxLayout()
        hboxIS.addWidget(QLabel("Enter the index of the initial step :"))
        hboxIS.addWidget(self.i_g)
        
        layout.addWidget(self.umdfileg_edit)
        layout.addWidget(selectButton)
        layout.addItem(QSpacerItem(0,50))
        layout.addLayout(hboxSF)
        layout.addItem(QSpacerItem(0,50))
        layout.addLayout(hboxDI)
        layout.addItem(QSpacerItem(0,50))
        layout.addLayout(hboxIS)
        layout.addItem(QSpacerItem(0,50))
        layout.addWidget(self.gofrCoresBox)
        layout.addWidget(self.gofrCores)
        layout.addItem(QSpacerItem(0,50))
        layout.addWidget(computeButton)
#        layout.addItem(QSpacerItem(0,50))
        layout.addWidget(self.gofrmessage)
        
        tab_widget.setLayout(layout)
        
        
    def gofrCoresBoxchange(self,state):
        if state ==2 :
            self.gofrCores.setVisible(True)
        else :
            self.gofrCores.setVisible(False)

    def compute_gofr(self):
        if not os.path.isfile(self.umdfileg_edit.text()):
            self.gofrmessage.setText("Error : the file "+self.umdfileg_edit.text()+" is displaced or missing. Pleace check the path or the file name.")
        elif not self.s_g.text().isnumeric() or int(self.s_g.text())<1:
            self.gofrmessage.setText("Error : the sampling frequency must be a strictly positive integer.")
        elif not isfloat(self.d_g.text()) or float(self.d_g.text())<0:
            self.gofrmessage.setText("Error : the discretization interval must be a positive rational number.")
        elif not self.i_g.text().isnumeric() or int(self.i_g.text())<0:
            self.gofrmessage.setText("Error : the initial step must be a positive integer.")
        else :
            argv = ["-f",self.umdfileg_edit.text(),"-s",self.s_g.text(),"-d",self.d_g.text(),"-i",self.i_g.text(),"-k",self.gofrCores.text()]
            print("argv=",argv)
            result = usefunction(gofr_umd,argv,self.gofrmessage)
            if result==True :
                Name = self.umdfileg_edit.text().split("/")
                name = Name[-1][:-8]+".gofr.dat"
                
                ff=open(name,'r+')
                lines = ff.readlines()
                newFile = [line.replace(",",".") for line in lines]#For some reason the floats in the gofr file when produced from the GUI are written with a comma "," instead of a point "."
                ff.seek(0)
                ff.writelines(newFile)
                ff.close()
                    
                self.gofrmessage.setText("gofr file successfully created under the name "+name)








    def create_Vibration_layout(self,tab_widget):
        self.layoutvib = QVBoxLayout()
        
        self.messageVib=QLabel("")

        self.temperature = QLineEdit("5000")
        self.UMDfileVib_edit = QLineEdit("")
        
        UMDButton = QPushButton("Select the UMD file")
        
        UMDButton.clicked.connect(partial(select_File,self.UMDfileVib_edit))

        hboxT=QHBoxLayout()
        hboxT.addWidget(QLabel("Select the temperature (default 5000 K) :"))
        hboxT.addWidget(self.temperature)
                
        
        ComputeVibButton = QPushButton("Compute vibrational spectrum")
        ComputeVibButton.clicked.connect(self.ComputeVib)

        self.displayVib = QCheckBox("Graphically display the vibrational spectrums after calculation")
        self.displayVib.setChecked(True)


        self.layoutvib.addWidget(self.UMDfileVib_edit)
        self.layoutvib.addWidget(UMDButton)
        self.layoutvib.addItem(QSpacerItem(0,70))
        self.layoutvib.addLayout(hboxT)
        self.layoutvib.addItem(QSpacerItem(0,70))
        self.layoutvib.addWidget(ComputeVibButton)
        self.layoutvib.addWidget(self.displayVib)
        self.layoutvib.addItem(QSpacerItem(0,-10))
        self.layoutvib.addWidget(self.messageVib)
        
        tab_widget.setLayout(self.layoutvib)
        
        
        
    def ComputeVib(self):
        argv = ["-f",self.UMDfileVib_edit.text(),"-t",self.temperature.text()]
        
        if self.UMDfileVib_edit.text()=="":
            self.messageVib.setText("Error : Please select an UMD file.")
        elif not os.path.isfile(self.UMDfileVib_edit.text()):
            self.messageVib.setText("Error : the file "+self.UMDfileVib_edit.text()+" is displaced or missing.") 
        elif not (self.temperature.text().isnumeric() and int(self.temperature.text())>=0):
            self.messageVib.setText("Error : the temperature must be a strictly positive integer.")
        else :
            self.messageVib.setText("Computing...") 
            ff=open(self.UMDfileVib_edit.text(),"r")
            while True :
                line=ff.readline()
                if not line :
                    break
                line=line.strip().split()
                if len(line)>0 and line[0] == "types": 
                    types = [int(line[i]) for i in range (1,len(line))]
                    break
            ff.close()
            result=usefunction(vibr_spectrum_umd_fast,argv,self.messageVib)
            if result and self.displayVib.isChecked() :
                self.messageVib.setText("Vibrational spectrum successfully calculated. Files created under the names "+self.UMDfileVib_edit.text()[:-8]+".vels.scf.dat and "+self.UMDfileVib_edit.text()[:-8]+".vibr.dat\nDOS showed in graphs (see the tab <Visual Display> above)")
                DOS,Freq,Elements = read_vibr_allEls(self.UMDfileVib_edit.text()[:-8]+".vibr.dat")
                createhist_vib(self,DOS, Freq, Elements,types)
            elif result :
                self.messageVib.setText("Vibrational spectrum successfully calculated. Files created under the names "+self.UMDfileVib_edit.text()[:-8]+".vels.scf.dat and "+self.UMDfileVib_edit.text()[:-8]+".vibr.dat")
                



    def create_LAMMPS_layout(self,tab_widget):
        layout = QVBoxLayout()
        
        self.LAMMPSfile = QLineEdit("")
        self.logfile = QLineEdit("")
        self.pressfile = QLineEdit("")
        
        
        selectLAMMPSButton = QPushButton("Select the LAMMPS file")
        selectLAMMPSButton.clicked.connect(partial(select_File,self.LAMMPSfile))
        selectlogButton = QPushButton("Select the log file")
        selectlogButton.clicked.connect(partial(select_File,self.logfile))
        selectpressButton = QPushButton("Select the press file (optional)")
        selectpressButton.clicked.connect(partial(select_File,self.pressfile))
        submit_button = QPushButton("Convert")
        submit_button.clicked.connect(self.convert_LAMMPS)
        
        self.L2Umessage = QLabel("")
        layout.addItem(QSpacerItem(0,50))
        layout.addWidget(self.LAMMPSfile)
        layout.addWidget(selectLAMMPSButton)
        layout.addItem(QSpacerItem(0,50))
        layout.addWidget(self.logfile)
        layout.addWidget(selectlogButton)
        layout.addItem(QSpacerItem(0,50))
        layout.addWidget(self.pressfile)
        layout.addWidget(selectpressButton)
        layout.addItem(QSpacerItem(0,50))
        layout.addWidget(submit_button)
        layout.addWidget(self.L2Umessage)
        
        tab_widget.setLayout(layout)
    
    def convert_LAMMPS(self):
        if not os.path.isfile(self.LAMMPSfile.text()):
            self.L2Umessage.setText("ERROR : LAMMPS file not found. Please check the path or the file name.")        
        elif not os.path.isfile(self.logfile.text()):
            self.L2Umessage.setText("ERROR : log file not found. Please check the path or the file name.")        
        else:            
            if not os.path.isfile(self.pressfile.text()):
                argv=["-f",self.LAMMPSfile.text(),"-l",self.logfile.text()]
            else :                
                argv=["-f",self.LAMMPSfile.text(),"-l",self.logfile.text(),"-a",self.pressfile.text()]
            result = usefunction(LAMMPSParser2umd,argv,self.L2Umessage)
            if result :
                Name = self.LAMMPSfile.text().split("/")
                name = Name[-1]+".umd.dat"
                self.L2Umessage.setText("The UMD file has been successfully created under the name "+name)



    def create_UMD_layout(self,tab_widget):
        layout = QVBoxLayout()
        
        self.VASPfile = QLineEdit("")
        self.istepVasp = QLineEdit("0")
        
        selectVASPButton = QPushButton("Select the VASP file")
        selectVASPButton.clicked.connect(partial(select_File,self.VASPfile))
        submit_button = QPushButton("Convert")
        submit_button.clicked.connect(self.convert_VASP)
        
        self.VASPmessage = QLabel("")
        
        layout.addItem(QSpacerItem(0,50))
        layout.addWidget(self.VASPfile)
        layout.addWidget(selectVASPButton)
        layout.addWidget(QLabel("Initial step :"))
        layout.addWidget(self.istepVasp)
        layout.addItem(QSpacerItem(0,50))
        layout.addWidget(submit_button)
        layout.addWidget(self.VASPmessage)
        
        tab_widget.setLayout(layout)
                
    
    def convert_VASP(self):
        if not os.path.isfile(self.VASPfile.text()):
            self.VASPmessage.setText("ERROR : the VASP file "+self.VASPfile.text()+ " is displaced or missing. Please check the path or filename.")        
        else:
            argv=["-f",self.VASPfile.text(),"-i",self.istepVasp.text()]
            result = usefunction(VaspParser2umd,argv,self.VASPmessage)
            if result :
                Name = self.VASPfile.text().split("/")
                name = Name[-1]+".umd.dat"
                self.VASPmessage.setText("The UMD file has been successfully created under the name "+name)


            
    def create_UMD_to_LAMMPS_layout(self,tab_widget):
        layout = QVBoxLayout()
        
        self.U2Lfile = QLineEdit("")
        self.sampFreqU2Lammps = QLineEdit("1")
        
        selectUMDButton = QPushButton("Select the UMD file")
        selectUMDButton.clicked.connect(partial(select_File,self.U2Lfile))
        submit_button = QPushButton("Convert")
        submit_button.clicked.connect(self.convert_U2L)
        
        self.U2Lmessage = QLabel("")
        
        layout.addItem(QSpacerItem(0,50))
        layout.addWidget(self.U2Lfile)
        layout.addWidget(selectUMDButton)
        layout.addWidget(QLabel("Select the sampling frequency :"))
        layout.addWidget(self.sampFreqU2Lammps)
        layout.addItem(QSpacerItem(0,50))
        layout.addWidget(submit_button)
        layout.addWidget(self.U2Lmessage)
        
        tab_widget.setLayout(layout)
                
    
    def convert_U2L(self):
        if not os.path.isfile(self.U2Lfile.text()):
            self.U2Lmessage.setText("ERROR : UMD file "+self.U2Lfile.text()+" is displaced or missing. Please check the path or the file name.")        
        else:
            argv=["-f",self.U2Lfile.text(),"-s",self.sampFreqU2Lammps.text()]
            result = usefunction(umd_to_lammps,argv,self.U2Lmessage)
            if result :
                Name = self.U2Lfile.text().split("/")
                name = Name[-1][:-8]+".lammps"
                self.U2Lmessage.setText("The LAMMPS file has been successfully created under the name "+name)

    

    def create_Speciation_layout(self, tab_widget):
        layout = QVBoxLayout()
        
        self.bondfile = QLineEdit("")
        self.centralatom = QLineEdit("")
        self.outeratom = QLineEdit("")
        self.umdfileSpec = QLineEdit("")
        self.umdfileSpec.hide()
        self.rings = QCheckBox("All levels of coordination")
        self.angles = QCheckBox("Calculate the angles whithin each polyhedra (needs to have, and sets, a degree of coordination equal to 1)")
        self.rings_edit = QLineEdit("1")
        
        self.rings = QCheckBox("Polymerization (All levels, r = 0)")
        self.angle = False
        
    
        selectbondbutton = QPushButton("Select the bond file")
        selectbondbutton.clicked.connect(partial(select_File,self.bondfile))
        self.specmessage = QLabel("")
        
        submit_button = QPushButton("Compute")
        submit_button.clicked.connect(self.compute_speciation)
        
        self.rings.stateChanged.connect(self.changering)
        self.angles.stateChanged.connect(self.changeangles)
        
        self.rings_edit.textChanged.connect(self.check_boxes)
        
        self.SpecCoresBox = QCheckBox("Precise the number of cores to be used for the parallelization")        
        self.SpecCores = QLineEdit("")
        self.SpecCores.setMaximumWidth(50)
        self.SpecCores.setVisible(False)
        self.SpecCoresBox.stateChanged.connect(self.SpecCoresBoxchange)
        
        self.umdSpecbutton = QPushButton("Select the matching umd file to compute the angles.")
        self.umdSpecbutton.hide()
        self.umdSpecbutton.clicked.connect(partial(select_File,self.umdfileSpec))
        
        hboxCA=QHBoxLayout()
        hboxCA.addWidget(QLabel("Central atoms elements :"))
        hboxCA.addWidget(self.centralatom)
        
        hboxOA=QHBoxLayout()
        hboxOA.addWidget(QLabel("Outer atoms elements :"))
        hboxOA.addWidget(self.outeratom)

        
        layout.addWidget(self.bondfile)
        layout.addWidget(selectbondbutton)
        layout.addItem(QSpacerItem(0,50))
        layout.addLayout(hboxCA)
        layout.addItem(QSpacerItem(0,50))
        layout.addLayout(hboxOA)
        layout.addItem(QSpacerItem(0,50))
        layout.addWidget(QLabel("Enter the degree of coordination :"))
        layout.addWidget(self.rings_edit)
        layout.addWidget(self.rings)
        layout.addWidget(self.angles)
        layout.addWidget(self.SpecCoresBox)
        layout.addWidget(self.SpecCores)
        layout.addItem(QSpacerItem(0,15))
        layout.addWidget(self.umdfileSpec)
        layout.addWidget(self.umdSpecbutton)
        layout.addItem(QSpacerItem(0,30))
        layout.addWidget(submit_button)
        layout.addItem(QSpacerItem(0,30))
        layout.addWidget(self.specmessage)
        
#        layout.setSpacing(2)
        
        tab_widget.setLayout(layout)

    def check_boxes(self):
        if self.rings_edit.text()!="1" and self.angle:
            self.angles.toggle()
        if self.rings_edit.text()!="0" and self.rings_edit.text()!="All levels (0)" and self.rings.isChecked() :
            self.rings.toggle()

    def changering(self,state):
        if state == 2:
            self.ring = 0
            self.rings_edit.setText("All levels (0)")
        else :
            if not self.rings_edit.text().isnumeric():
                self.ring = 1
                self.rings_edit.setText("1")
                
    def changeangles(self,state):
        if state == 2:
            self.ring = 1
            self.angle = True
            self.rings_edit.setText("1")
            self.umdfileSpec.show()
            self.umdSpecbutton.show()
        else :
            self.angle = False
            self.umdfileSpec.hide()
            self.umdSpecbutton.hide()

    def SpecCoresBoxchange(self,state):
        if state ==2 :
            self.SpecCores.setVisible(True)
        else :
            self.SpecCores.setVisible(False)

    
    def compute_speciation(self):
        
        self.specmessage.setText("")
        
        if self.rings_edit.text().isnumeric():
            self.ring = self.rings_edit.text()
        elif self.rings_edit.text() != "All levels (0)":
            self.rings_edit.setText("1")
            self.ring=1
        
        compute = True

        if self.SpecCoresBox.isChecked() and not self.SpecCores.text().isnumeric():
            self.specmessage.setText("Error : the number of cores has to be a strictly positive integer.")
            compute = False
                    
        if self.outeratom.text()=="":
            self.specmessage.setText("Please select at least one element for the outer atoms.")
            compute = False
        elif self.centralatom.text()=="":
            self.specmessage.setText("Please select at least one element for the central atoms.")
            compute = False
            
        if self.bondfile.text()=="":
            self.specmessage.setText("Please select a bond file before launching computation.")            
            compute = False
        elif not os.path.isfile(self.bondfile.text()):
            self.specmessage.setText("ERROR : the file "+self.bondfile.text()+" is displaced or missing")
            compute = False
        

        if compute :
            if self.angle :
                if os.path.isfile(self.umdfileSpec.text()):                    
                    argv=["-f",self.bondfile.text(),"-u",self.umdfileSpec.text(),"-c",self.centralatom.text(),"-a",self.outeratom.text(),"-r", self.ring]
                elif self.umdfileSpec.text()=="" :
                    self.specmessage.setText("If you want to calculate the angles, please select an UMD file before launching computation.")            
                    compute = False
                else :
                    self.specmessage.setText("ERROR : the file "+self.umdfileSpec.text()+" is displaced or missing")
                    compute = False
            else :
                argv=["-f",self.bondfile.text(),"-c",self.centralatom.text(),"-a",self.outeratom.text(),"-r", self.ring]
                if self.SpecCoresBox.isChecked():
                    argv+=["-k",self.SpecCores.text()]

 
        if compute :
            result=usefunction(speciation_and_angles,argv,self.specmessage)
            if result == True :
                self.specmessage.setText("Population file successfully created under the name "+self.bondfile.text().split("/")[-1][:-4]+".r"+str(self.ring)+".popul.dat")
            elif result != False :
                if len(result) ==1 :
                    self.specmessage.setText("ERROR : element "+result[0]+" not present in the provided simulation")
                else :
                    string = ""
                    for el in result[:-2]:
                        string+=el+", "
                    string+=result[-2]+" and "+result[-1]
                    self.specmessage.setText("ERROR : elements "+string+" not present in the provided simulation")



    def create_Bond_layout(self, tab_widget):
        layout = QVBoxLayout()
        # Création des widgets de saisie
        self.umdfileBond = QLineEdit("")
        self.inpfile = QLineEdit("")
                                
        self.bondmessage = QLabel("")
        
        self.l_edit = QLineEdit('')
        self.s_edit = QLineEdit('1')
        
        self.BondCoresBox = QCheckBox("Precise the number of cores to be used for the parallelization")        
        self.BondCores = QLineEdit("")
        self.BondCores.setMaximumWidth(50)
        self.BondCores.setVisible(False)
        self.BondCoresBox.stateChanged.connect(self.BondCoresBoxchange)
        
        self.specificBox = QCheckBox("Restrict bonding to a fraction of the cell")
        
        self.X0 = QLineEdit("0")
        self.X1 = QLineEdit("0")
        self.Y0 = QLineEdit("0")
        self.Y1 = QLineEdit("0")
        self.Z0 = QLineEdit("0")
        self.Z1 = QLineEdit("0")
        
        specificLayout = QHBoxLayout()
        self.specificWidget = QWidget()
        
        XLayout = QVBoxLayout()
        YLayout = QVBoxLayout()
        ZLayout = QVBoxLayout()
        LabLayout = QVBoxLayout()

        specificLayout.addLayout(LabLayout)
        specificLayout.addLayout(XLayout)
        specificLayout.addLayout(YLayout)
        specificLayout.addLayout(ZLayout)

        LabLayout.addWidget(QLabel(""))
        LabLayout.addWidget(QLabel("Lower limit"))
        LabLayout.addWidget(QLabel("Upper limit"))
        XLayout.addWidget(QLabel("x"))
        XLayout.addWidget(self.X0)
        XLayout.addWidget(self.X1)
        YLayout.addWidget(QLabel("y"))
        YLayout.addWidget(self.Y0)
        YLayout.addWidget(self.Y1)
        ZLayout.addWidget(QLabel("z"))
        ZLayout.addWidget(self.Z0)
        ZLayout.addWidget(self.Z1)

        self.specificWidget.setVisible(False)
        
        self.specificWidget.setLayout(specificLayout)

        self.specificBox.stateChanged.connect(self.specBoxchange)

        # Création d'un bouton pour sélectionner un fichier
        select_button = QPushButton("Select the UMD file")
        select_button.clicked.connect(partial(select_File,self.umdfileBond))

        select_buttonI = QPushButton("Select the average bonding input file")
        select_buttonI.clicked.connect(partial(select_File,self.inpfile))

        
        # Création d'un bouton pour soumettre les paramètres
        submit_button = QPushButton("Compute")
        submit_button.clicked.connect(self.compute_Bond)

        hboxBL = QHBoxLayout()
        hboxBL.addWidget(QLabel("Bonding length in Angströms (optional ; will overwrite the values of the input file) :"))
        hboxBL.addWidget(self.l_edit)

        hboxSF = QHBoxLayout()        
        hboxSF.addWidget(QLabel("Sampling frequency :"))
        hboxSF.addWidget(self.s_edit)

        # Création d'un layout vertical pour organiser les widgets
        layout.addWidget(self.umdfileBond)
        layout.addWidget(select_button)
        layout.addItem(QSpacerItem(0,50))
        layout.addWidget(self.inpfile)
        layout.addWidget(select_buttonI)
        layout.addItem(QSpacerItem(0,50))
        layout.addLayout(hboxBL)
        layout.addItem(QSpacerItem(0,50))
        layout.addLayout(hboxSF)
        layout.addWidget(self.BondCoresBox)
        layout.addWidget(self.BondCores)
        layout.addWidget(self.specificBox)
        layout.addWidget(self.specificWidget)
        layout.addItem(QSpacerItem(0,50))
        layout.addWidget(submit_button)
        layout.addWidget(self.bondmessage)

        # Création d'un widget de base et définition du layout comme layout principal
        widget = QWidget()
        widget.setLayout(layout)

        # Définition du widget comme widget central de la fenêtre principale
        self.setCentralWidget(widget)
        
        tab_widget.setLayout(layout)
        
    def compute_Bond(self):
          
        if not os.path.isfile(self.umdfileBond.text()):
            if self.umdfileBond.text=="":
                self.bondmessage.setText("Please select an umd file before launching computation.")
            else :
                self.bondmessage.setText("Error : the file "+self.umdfileBond.text()+" is displaced or missing. Please check the path or file name.") 
        elif not os.path.isfile(self.inpfile.text()) and not(isfloat(self.l_edit.text())):
            if self.inpfile.text()=="":
                self.bondmessage.setText("Please select an input file or set manually the value of the bonding length before launching the computation.")
            else :
                self.bondmessage.setText("Error : the file "+self.inpfile.text()+" is displaced or missing. Please check the path or file name.") 
        elif not(self.s_edit.text().isnumeric()) or int(self.s_edit.text())<1 :        
            self.bondmessage.setText("Error : the value of the sampling frequency has to be a strictly positive integer.")
        elif self.BondCoresBox.isChecked() and not self.BondCores.text().isnumeric():
            self.bondmessage.setText("Error : the number of cores has to be a strictly positive integer.")
        else :
            if not(isfloat(self.l_edit.text())):
                argv=["-f",self.umdfileBond.text(),"-s",self.s_edit.text(),"-i",self.inpfile.text()]
                self.l_edit.setText("")
                self.bondmessage.setText("The .inp file has been used to determine the bonds.")
            else :
                argv=["-f",self.umdfileBond.text(),"-s",self.s_edit.text(),"-l",self.l_edit.text(),"-i",self.inpfile.text()]
                self.bondmessage.setText("The explicit value of "+self.l_edit.text()+" (Angströms) has been used to determine the bonds.")
            
            specified = False
            if self.specificBox.isChecked():
                if isfloat(self.X0.text()) and isfloat(self.X1.text()) and isfloat(self.Y0.text()) and isfloat(self.Y1.text()) and isfloat(self.Z0.text()) and isfloat(self.Z1.text()):
                    specList = "["+self.X0.text()+","+self.Y0.text()+","+self.Z0.text()+","+self.X1.text()+","+self.Y1.text()+","+self.Z1.text()+","+"]"
                    argv+=["-p",specList]
                    specified = True
            if self.BondCoresBox.isChecked():
                argv+=["-k",self.BondCores.text()]


            result = usefunction(bonding_umd,argv,self.bondmessage)                    
                
            if result :
                Text = "Bonds file successfully created under the name "+self.umdfileBond.text().split("/")[-1][:-8]+".bonding.dat"
                if self.specificBox.isChecked() :
                    if specified == False :
                        Text+="\nThe bonds in the whole cell have been computed, since the values provided for the restricting conditions weren't valid."
                    else :
                        Text+="\nThe bonds have been calculated in a restricted space interval."
                self.bondmessage.setText(Text) 
                self.bondfile.setText(self.umdfileBond.text()[:-8]+".bonding.dat")


    def specBoxchange(self,state):
        if state ==2 :
            self.specificWidget.setVisible(True)
        else :
            self.specificWidget.setVisible(False)
    
    def BondCoresBoxchange(self,state):
        if state ==2 :
            self.BondCores.setVisible(True)
        else :
            self.BondCores.setVisible(False)
        
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
