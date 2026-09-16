#!/usr/bin/env python3
# coding: utf-8

# In[139]:


#A script that separates the r1 calculated stat and popul files for liquid and vapor species


# In[140]:


import glob
import os
import getopt
import sys
import re
import numpy as np
from itertools import groupby


# In[141]:


####################################Grab user inputs#################
argv=sys.argv[1:]
try:
    opts, args = getopt.getopt(argv, "hs:p:l:")
except:
    print('Error')
for opt, arg in opts:
    if opt=='-h':
        print('-s file.r0.step.dat -p file_r1.popul.dat file -l cutoff length for vapor species (default is 12 atoms)')
    if opt in ['-s']:
        step=str(arg)
    if opt in ['-p']:
        popul=str(arg)
    if opt in ['-l']:
        length=int(arg)
    else:
        length=12


# In[8]:


#Read in data from step.dat file in groups, separated by timestep
def get_groups(seq, group_by):
    data = []
    for line in seq:
        # Here the `startswith()` logic can be replaced with other
        # condition(s) depending on the requirement.
        if line.startswith(group_by):
            if data:
                yield data
                data = []
        data.append(line)

    if data:
        yield data


# In[94]:


#From the r0 step file, get a list of all vapor and liquid atoms at each timestep
vapor_atoms_dict={}
liquid_atoms_dict={}
with open(step, 'r') as f:
    next(f)
    for i, group in enumerate(get_groups(f, 'step'), start=1):
        if i==2:
            ts=float(group[1].split('\t')[1].strip('\n'))
        vapor_atoms=[]
        liquid_atoms=[]
        timestep=float(group[0].split('\t')[1].strip('\n'))
        for j in range(3,len(group)):
            if int(group[j].split('\t')[1])<length:
                for k in group[j].split('\t')[2].strip('\n'+'['+']').split(','):
                    vapor_atoms.append(int(k))
            else:
                for k in group[j].split('\t')[2].strip('\n'+'['+']').split(','):
                    liquid_atoms.append(int(k))
        vapor_atoms_dict[timestep]=vapor_atoms
        liquid_atoms_dict[timestep]=liquid_atoms


# In[85]:


#Make lists for every cluster, the timesteps in which it is alive in the liquid phase, the vapor phase, 
#and all the atoms making up the cluster.
formula_list=[]
liquid_time_list=[]
vapor_time_list=[]
composition_list=[]
with open(popul, 'r') as f:
    next(f)
    next(f)
    for line in f:
        split_line=line.strip('\n').split('\t')
        formula_list.append(split_line[0])
        liquid_time=[]
        vapor_time=[]
        composition_list.append(split_line[4])
        for i in range(int(split_line[1]),int(split_line[2])+1):
            counter=0
            for j in split_line[4].strip('['+']').split(','):
                if int(j) in liquid_atoms_dict[i]:
                    counter=counter+1
            if counter>0:
                liquid_time.append(i)
            else:
                vapor_time.append(i)
        liquid_time_list.append(liquid_time)
        vapor_time_list.append(vapor_time)


# In[106]:


#Write the population file for r1 vapor phases
nf=open(popul[:-9]+'vapor.popul.dat', 'w')
nf.write('Formula\tBegin (step)\tEnd (step)\tLifetime (fs)\tComposition\n')
for i in range(len(formula_list)):
    begin=[]
    end=[]
    if vapor_time_list[i] != []:
        if sorted(vapor_time_list[i]) == list(range(min(vapor_time_list[i]), max(vapor_time_list[i])+1)):
            begin.append(min(vapor_time_list[i]))
            end.append(max(vapor_time_list[i]))
        else:
            out = []
            for _, g in groupby(enumerate(vapor_time_list[i]), lambda x: x[0] - x[1]):
                out.append([v for _, v in g])
            for j in out:
                begin.append(min(j))
                end.append(max(j))
        for j in range(len(begin)):
            lifetime=(end[j]-begin[j]+1)*ts
            nf.write(formula_list[i]+'\t'+str(begin[j])+'\t'+str(end[j])+'\t'+str(lifetime)+'\t'
                     +composition_list[i]+'\n')
nf.close()


# In[136]:


#Write the population file for r1 liquid phases
nf=open(popul[:-9]+'liquid.popul.dat', 'w')
nf.write('Formula\tBegin (step)\tEnd (step)\tLifetime (fs)\tComposition\n')
for i in range(len(formula_list)):
    begin=[]
    end=[]
    if liquid_time_list[i] != []:
        if sorted(liquid_time_list[i]) == list(range(min(liquid_time_list[i]), max(liquid_time_list[i])+1)):
            begin.append(min(liquid_time_list[i]))
            end.append(max(liquid_time_list[i]))
        else:
            out = []
            for _, g in groupby(enumerate(liquid_time_list[i]), lambda x: x[0] - x[1]):
                out.append([v for _, v in g])
            for j in out:
                begin.append(min(j))
                end.append(max(j))
        for j in range(len(begin)):
            lifetime=(end[j]-begin[j]+1)*ts
            nf.write(formula_list[i]+'\t'+str(begin[j])+'\t'+str(end[j])+'\t'+str(lifetime)+'\t'
                     +composition_list[i]+'\n')
nf.close()


# In[137]:


#Write a new stat.dat file for vapor species
species_dict={}
nf=open(popul[:-9]+'vapor.stat.dat', 'w')
nf.write('Cluster\tTime (fs)\tPercent\tNumber of atoms\n')
with open(popul[:-9]+'vapor.popul.dat', 'r') as f:
    next(f)
    time=0
    for line in f:
        split_line=line.strip('\n').split('\t')
        time=time+float(split_line[3])
        if split_line[0] in species_dict:
            species_dict[split_line[0]]=[species_dict[split_line[0]][0]+float(split_line[3]),
            len(split_line[4].split(','))]
        else:
            species_dict[split_line[0]]=[float(split_line[3]), len(split_line[4].split(','))]
for key in species_dict.keys():
    nf.write(key+'\t'+str(species_dict[key][0])+'\t'+str(species_dict[key][0]/time)+'\t'+str(species_dict[key][1])
             +'\n')
nf.close()


# In[138]:


#Write a new stat.dat file for liquid species
species_dict={}
nf=open(popul[:-9]+'liquid.stat.dat', 'w')
nf.write('Cluster\tTime (fs)\tPercent\tNumber of atoms\n')
with open(popul[:-9]+'liquid.popul.dat', 'r') as f:
    next(f)
    time=0
    for line in f:
        split_line=line.strip('\n').split('\t')
        time=time+float(split_line[3])
        if split_line[0] in species_dict:
            species_dict[split_line[0]]=[species_dict[split_line[0]][0]+float(split_line[3]),
            len(split_line[4].split(','))]
        else:
            species_dict[split_line[0]]=[float(split_line[3]), len(split_line[4].split(','))]
for key in species_dict.keys():
    nf.write(key+'\t'+str(species_dict[key][0])+'\t'+str(species_dict[key][0]/time)+'\t'+str(species_dict[key][1])
             +'\n')
nf.close()

