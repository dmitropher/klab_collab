#!/usr/bin/python
from transform import *
import math

atom_list=['O4','Co4']
resname='NPH'
res_conversion={'1':'193','2':'201','1A':'204','2A':'212'}

with open ('test.pdb', 'r') as myfile:
    ligandPdb = myfile.readlines()

with open ('test_full.pdb', 'r') as myfile:
    pdb = myfile.readlines()

atomInfo=[]

for line in ligandPdb:
    if line[12:16].strip(' ') in atom_list:
        coords = [a for a in line[30:56].split(' ') if a != '']
        coords.append(line[12:16].strip(' '))
        coords.append(line[23:27].strip(' ')+line[21:23].strip(' '))
        atomInfo.append(coords)
resPositions=[]
for line in pdb:
    if line[17:21].strip(' ') == resname:
        if line[23:27].strip(' ') not in resPositions:
            resPositions.append(line[23:27].strip(' '))
for i,atom in enumerate(atomInfo):
    for otherAtom in atomInfo[i:]:
        if atom[4] != otherAtom[4]:

            cstLine= 'AtomPair ' + atom[3] + ' ' + res_conversion[atom[4][:-1]+'A'] + atom[4][-1] + ' ' +  otherAtom[3] + ' ' + res_conversion[otherAtom[4][:-1] + 'A'] + otherAtom[4][-1] + ' GAUSSIANFUNC ' +\
            str(
                round(
                    math.sqrt(\
                    math.pow( float(atom[0]) - float(otherAtom[0]), 2) + \
                    math.pow( float(atom[1]) - float(otherAtom[1]), 2) +\
                    math.pow( float(atom[2]) - float(otherAtom[2]), 2)\
                ),2)) +\
            ' 0.05'
            print cstLine
print resPositions
print atomInfo
