#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 17 07:43:57 2023

@author: Kevin Jiguet-Covex
"""
import sys, getopt, os.path, itertools, math
import numpy as np

def Variables(filename):#Check which variables are present in the file and in which order
    ff=open(filename,'r')
    fline = ff.readline().split()
    dumpstyle='MINI'#If it stays like this, the file is standardized and contains x, y and z only
    indtype=None
    PresList,PosList=None,None
    if len(fline)>0 and fline[0]=='ITEM:':   #Else, it can contain any variable in any order
        dumpstyle='ITEMS'
        VarList=['xs','ys','zs','x','y','z','xa','ya','za','vx','vy','vz','fx','fy','fz','id']
        PresList=[0 for _ in range(16)]
        PosList=[]
        while len(fline)<2 or fline[1]!="ATOMS":
            fline=ff.readline().split()
        for var in fline[2:]:
            if var in VarList :
                PresList[VarList.index(var)]=1
                PosList.append(VarList.index(var))
            elif var=='type':#Information about the element.
                indtype=len(PosList)
                PosList.append(-1)
            else:
                PosList.append(-1)

    ff.close()
    return PresList,PosList,indtype,dumpstyle

def Crystal(filename,logname,indtype):#Extracts the information about the system (# of atoms, elements, etc.)
    natom = 0
    types =  []
    elements = []
    ntypat = 0
    typat = []
    t0 = 0
    style=""#The way information about the parameters (T,P,etc.) is displayed in the log file
    
    fl=open(logname,"r")
    ff=open(filename,"r")

    while True : 
        line=fl.readline()
        if not line : 
            break
        line = line.strip()
        
        if line == "Reading data file ...":
            line=fl.readline().strip().split()
            acellx = float(line[7].strip('('))-float(line[3].strip('('))
            acelly = float(line[8].strip('('))-float(line[4].strip('('))
            acellz = float(line[9].strip(')'))-float(line[5].strip(')'))
        elif line == "# SIMULATION":
            break
        else :
            line = line.split()
            if line !=[]:
                if line[0] == "timestep":
                    if line[1] == "${dt}":
                        line = fl.readline().strip().split()
                        timestepLog = float(line[1])*1000
                    else :
                        timestepLog = float(line[1])*1000
                elif line[0] == "thermo_style":
                    style = line[1]

    isatoms = 0
    fl.close()
    fl=open(logname,'r')

    #Those variables code for the position of each parameter in the log file
    IEindex = None
    Enthindex = None
    Tempindex = None
    Pressindex = None

    while True :
        line=fl.readline()
        if not line : 
            break
        line = line.strip().split()
        if len(line)>0 and line[0]=="group" and isatoms!=2:
            line = fl.readline().strip().split()
            types.append(int(line[0]))
            natom+=types[-1]
            typat+=[ntypat for _ in range(types[-1])]
            ntypat+=1
            isatoms=1
        elif isatoms==1:
            isatoms=2
            
        if len(line)>3 and (line[0]=='dump_modify' or line[0]=='#dump_modify') and line[1]=="d1" and line[2]=="element":
            j=-1
            el=line[j]
            while el!="element":
                elements.append(el)
                j-=1
                el=line[j]
            elements.reverse()
            
        if len(line)>2 and line[0] == "Per" and line[1] == "MPI" and line[2] == "rank":
            if style == "one" or style == "custom":
                Variables = fl.readline().strip().split()
                for i in range(len(Variables)) : 
                    v=Variables[i]
                    if v == "Temp":
                        Tempindex = i
                    elif v == "Press":
                        Pressindex = i
                    elif v == "TotEng":
                        IEindex = i
                    elif v =="Enthalpy":
                        Enthindex = i
            break
    #Determining the timestep between two snapshots in the dump file
    t0=None        
            
    while True :
        l = ff.readline()
        if not l :
            break
        line = l.strip().split()
        if len(line)>1 and line[1] == "TIMESTEP":
            line=ff.readline().strip().split()
            if t0 == None :
                t0=float(line[0])
            else : 
                timestepFile = (float(line[0]) - t0)*timestepLog
                break
        if len(line)>1 and line[1]=="Timestep:":
            if t0 == None :
                t0=float(line[2])
            else : 
                timestepFile = (float(line[2]) - t0)*timestepLog
                break
                
    #If the information about atoms isn't fully contained in the log file, we look in the dump file
    if(natom==0):
        types=[]
        ntypat=0
        typat=[]
        elements=[]
        currentType=-1
        if indtype==None :
            print("error : no way to determine the species")
            sys.exit()
        line=ff.readline().strip().split()
        while len(line)<2 or line[1]!="ATOMS":
            line=ff.readline().strip().split()
        line=ff.readline().strip().split()
        while line[0] != "ITEM:":           
            if currentType != int(line[indtype])-1:
                ntypat+=1
                types.append(1)
                currentType = int(line[indtype])-1 
            else:
                types[-1]+=1
            line=ff.readline().strip().split()
                
            typat.append(currentType)
        icounter=0
        fl.close()
        fl=open(logname,'r')
        line=fl.readline().strip().split()
        while len(line)<3 or (line[0]!="use" or line[1]!="deepmd-kit" or line[2]!="at:"):
            line=fl.readline().strip().split()
            icounter+=1
            if icounter ==1000 :
                break
        j=-1
        while True:
            el=line[j]

            if el =="*":
                break
            else : 
                elements.append(el)
            j-=1
        elements.reverse()
        natom=len(typat)
    

    ff.close()
    fl.close()
    return natom,types,elements,ntypat,typat,[acellx,acelly,acellz],Tempindex,Pressindex,IEindex,Enthindex,timestepLog,timestepFile,style

def main(argv):
    
    filename = ""
    logname = ""
    allname = ""
    press = False
    
    try:
        opts,args = getopt.getopt(argv,"f:l:a:",["fFile","lLog","aAllpress"])
    except getopt.GetoptError:
        print("LAMMPSParser2umd to convert various LAMMPS files into UMD files")
        print("-f : LAMMPS file ; -l : log file ; -a : allpress file (optional)")
        sys.exit()
    for opt,arg in opts : 
        if opt in "-f":
            filename = str(arg)
        elif opt in "-l":
            logname = str(arg)
        elif opt in "-a":
            allname = str(arg)
            press = True
    print("Will use the file <"+filename+"> for positions") 
    print("Will use the file <"+logname+"> for thermodynamic data")
    if press :
        print("Will use the file <"+allname+"> for stress tensor values\n") 
        
    
    PresList,PosList,indtype,dumpstyle=Variables(filename)

    natom,types,elements,ntypat,typat,[acellx,acelly,acellz],Tempindex,Pressindex,IEindex,Enthindex,timestepLog,timestepFile,style = Crystal(filename,logname,indtype)
    if Tempindex == None and style == "one" or style == "custom":
        print("WARNING : Temperature information not present in log file. This can cause problem in some post-processing scripts.")

    fa = open(filename+".umd.dat",'w')    

    fa.write("natom "+str(natom)+"\n")
    fa.write("ntypat "+str(ntypat)+"\n")
    string = ""
    for n in types :
        string+=str(n)+" "
    string+="\n"
    fa.write("types "+string)
    string=""
    for el in elements :
        string+=el+" "
    string+="\n"
    fa.write("elements "+string)
    string=""
    for x in typat :
        string+=" "+str(x)
    string+="\n\n"
    fa.write("typat "+string)
    
    #If a parameter is not present in the log file, we give it the value '0'
    IE = '0'
    Enth = '0'
    Press = '0'
    Temp = '0'
        
    ff=open(filename,'r')
    fl=open(logname,'r')
    
    logline = fl.readline().strip().split()
    while len(logline)<3 or logline[0] != "Per" or logline[1] != "MPI" or logline[2] != "rank":
        logline = fl.readline().strip().split()
    logline = fl.readline().strip().split()
    
    if press :
        fall=open(allname,'r')
    while True :
        line = ff.readline()
        if not line :
            break
        line = line.strip().split()
        if len(line)>0 and (line[0] == "Atoms." or (line[0]=="ITEM:" and line[1]=="TIMESTEP")):
            if dumpstyle=='MINI':
                string="timestep "+str(timestepFile)+" fs\n"
                fa.write(string)
                string = "time "+str(float(line[2])*timestepLog)+" fs\n"
                fa.write(string)
            elif dumpstyle=='ITEMS':
                line = ff.readline().strip().split()
                string="timestep "+str(timestepFile)+" fs\n"
                fa.write(string)
                string = "time "+str(float(line[0])*timestepLog)+" fs\n"
                fa.write(string)
                while(len(line)<2 or line[1]!="ATOMS"):
                    line = ff.readline().strip().split()

            if style == "one" or style == "custom":         #Extraction of the parameters
                logline = fl.readline().strip().split()
                if Tempindex != None :
                    Temp = logline[Tempindex]
                if Pressindex != None :
                    Press = str(float(logline[Pressindex])/10000)
                if Enthindex != None :
                    Enth = logline[Enthindex]
                if IEindex != None :
                    IE = logline[IEindex]

            elif style == "multi":                          #Extraction of the parameters
                logline = fl.readline().strip().split()
                IE = logline[2]
                Temp = logline[8]
                logline = fl.readline().strip().split()
                logline = fl.readline().strip().split()
                logline = fl.readline().strip().split()
                Press = str(float(logline[8])/10000)
                logline = fl.readline().strip().split()

            
            rprimVecs = np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])
            cellVecs = [np.array([1.0,0.0,0.0])*acellx,np.array([0.0,1.0,0.0])*acellx,np.array([0.0,0.0,1.0])*acellx]
                
            rprimVecsstr=["  ".join(str(el) for el in vec) for vec in rprimVecs]
            cellVecsstr=["  ".join(str(el) for el in vec) for vec in cellVecs]
            acell=[acellx,acelly,acellz]
            acellLine = "acell "+str(acellx)+" "+str(acelly)+" "+str(acellz)+" A\n"
            rLines ="rprim_a "+rprimVecsstr[0] +" \nrprim_b "+ rprimVecsstr[1] + " \nrprim_c " +rprimVecsstr[2]+ " \n"
            rdLines ="rprimd_a "+cellVecsstr[0] +" A\nrprimd_b "+ cellVecsstr[1] + " A\nrprimd_c " +cellVecsstr[2]+ " A\n"
            
            string = "Internal Energy "+IE+" eV\nEnthalpy "+Enth+" eV\nTemperature "+Temp+" K\nPressure "+Press+" GPa\n"
            
            if press :
                allLine = fall.readline().strip().split()
            
                while allLine[0]!=line[2]:
                    allLine = fall.readline().strip().split()
                    if not allLine :
                        break
            
                stressline = "StressTensor "+str(float(allLine[1])/10000)+" "+str(float(allLine[2])/10000)+" "+str(float(allLine[3])/10000)+" "+str(float(allLine[4])/10000)+" "+str(float(allLine[5])/10000)+" "+str(float(allLine[6])/10000)+" GPa\n"

                    
            fa.write(string)
            if press :
                fa.write(stressline)
            fa.write(acellLine)
            fa.write(rLines)
            fa.write(rdLines)            
            fa.write("atoms: reduced*3 cartesian*3(A) abs.diff.*3(A) velocity*3(A/fs) force*3(ev/A) \n")
            
            if(dumpstyle=='ITEMS'):
                if PresList[-1]==1 :#If the atoms aren't in order and their index is precised in the variable "id"
                    posId = PosList.index(15)
                    StrList=["" for _ in range(natom)]
                    for at in range(natom):
                        string=""
                        line = ff.readline()
                        line = line.strip().split()
                        CoordsLine=['0.0' for _ in range(15)]
                        ind=0
                        for coord in line :#We parse the line and extract the variables
                            if PosList[ind]!=-1:
                                if ind!=posId :
                                    CoordsLine[PosList[ind]]=coord
                                else:
                                    atom=int(coord)
                            ind+=1
                        #Filling the missing ones with the information we already have
                        if PresList[3]!=0 and PresList[0]==0:
                            for i in range(3):
                                CoordsLine[i]=str(round(float(CoordsLine[3+i])/acell[i],5))


                        if PresList[6]==0:
                            if PresList[3]!=0 :
                                for i in range(3):
                                    CoordsLine[6+i]=CoordsLine[3+i]
                            elif PresList[0]!=0:
                                for i in range(3):
                                    CoordsLine[6+i]=str(round(float(CoordsLine[i])*acell[i],5))
                                    
                        if PresList[0]!=0 and PresList[3]==0:
                            for i in range(3):
                                CoordsLine[3+i]=str(round(float(CoordsLine[i])*acell[i],5))
                                
                        for i in range(15):
                            string+=CoordsLine[i]+" "
                        StrList[atom-1]=string+"\n"
                    
                    for at in range(natom):
                        fa.write(StrList[at])
                                        
                else:
                    for at in range(natom):
                        string=""
                        line = ff.readline()
                        line = line.strip().split()
                        CoordsLine=['0.0' for _ in range(15)]
                        ind=0
                        for coord in line :
                            CoordsLine[PosList[ind]]=coord
                            ind+=1
                        if(PresList[0]==0):
                            for i in range(3):
                                CoordsLine[i]=str(round(float(CoordsLine[3+i])/acell[i],5))
                        if(PresList[6]==0):
                            for i in range(3):
                                CoordsLine[6+i]=CoordsLine[3+i]
                        for i in range(15):
                            string+=CoordsLine[i]+" "
                    
                        fa.write(string+"\n")
    
            elif dumpstyle=='MINI':
                for at in range(natom):
                    line = ff.readline()
                    line = line.strip().split()
                    x,y,z = round(float(line[1]),5),round(float(line[2]),5),round(float(line[3]),5)
                    redx,redy,redz = round(x/acellx,5),round(y/acelly,5),round(z/acellz,5)
                    string = str(redx)+" "+str(redy)+" "+str(redz)+" "+str(round(math.fmod(acellx+x,acellx),5))+" "+str(round(math.fmod(acelly+y,acelly),5))+" "+str(round(math.fmod(acellz+z,acellz),5))+" "+str(x)+" "+str(y)+" "+str(z)+" 0.0 0.0 0.0 0.0 0.0 0.0"
                    fa.write(string+"\n")
                
            fa.write("\n")
    
    ff.close()
    fl.close()
    if press :
        fall.close()
    fa.close()
    print("umd file created under the name <"+filename+".umd.dat>")
            
        

if __name__ =="__main__":
    main(sys.argv[1:])
