# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 18:36:30 2018
pdb2lmp.py
@author: 97588
"""

#Define the filename of your pdb and lmp files
filename_pdb = 'sio2.pdb'
filename_lmp = 'sio2.dat'

#Specify atom types
atom_elements = ['Si4-','O2-']

#Getting the full-text of the pdb file
def load_text_pdb(filename_pdb):
    with open(filename_pdb,'r') as file_pdb:
        text_pdb = [line.rstrip() for line in file_pdb.readlines()]
        return [line.split() for line in text_pdb]
 
#Getting the size of the simulation box
def get_box_size(text_pdb):
    return [float(line) for line in text_pdb[2][1:7]]

#Getting the number of atoms and types
def get_number_atoms(text_pdb):
    atom_count=0
    for line in text_pdb:
        atom_count = atom_count + (1 if line[0] == 'ATOM' else 0) 
    return atom_count

#Specify atom type
def specify_type(atom_type,atom_elements):
    num_elements = len(atom_elements)
    for i in range(num_elements):
        if atom_type == atom_elements[i]:
            return int(i+1)

#Create a type list
def create_typelist(text_pdb, atom_count, atom_elements):
    return [specify_type(line[-1],atom_elements) \
            for line in text_pdb[9:9+atom_count]]

#Create a coordinate list
def create_coordinatelist(text_pdb, atom_count):
    return [line[5:8] for line in text_pdb[9:9+atom_count]]

#Create the lmp file
def create_lmp(filename_lmp,\
               simbox,atom_count,atom_elements,typelist,coordinatelist):
    with open(filename_lmp,'w') as file_lmp:
        file_lmp.write('# {}\n'.format(atom_elements))
        file_lmp.write('{} atoms\n'.format(atom_count))
        file_lmp.write('{} atom types\n'.format(len(atom_elements)))
        file_lmp.write('0.0 {} xlo xhi\n'.format(simbox[0]))
        file_lmp.write('0.0 {} ylo yhi\n'.format(simbox[1]))
        file_lmp.write('0.0 {} zlo zhi\n'.format(simbox[2]))
        file_lmp.write('Atoms\n\n')
        for i in range(atom_count):
            file_lmp.write(str(i+1) + ' ' + str(typelist[i]) + ' 0.0 ' + \
                           coordinatelist[i][0] + ' ' + \
                           coordinatelist[i][1] + ' ' + \
                           coordinatelist[i][2] + '\n')

#Main Script
text_pdb = load_text_pdb(filename_pdb)
simbox = get_box_size(text_pdb)
atom_count = get_number_atoms(text_pdb)
typelist=create_typelist(text_pdb,atom_count,atom_elements)
coordinatelist=create_coordinatelist(text_pdb,atom_count)
create_lmp(filename_lmp,simbox,atom_count,atom_elements,typelist,coordinatelist)






