# -*- coding: utf-8 -*-
"""
Created on Sun Mar 31 13:35:48 2019
Structural Analysis Beta 1.3
@author: Ye Tian
"""

from scipy.spatial.distance import pdist
from collections import Counter
import itertools
import numpy as np

# GLOBAL SETTINGS: MODIFY THIS PART!!!

# 1. Specify the filename of atomic structure
FILENAME = 'POSCAR'
# 2. Specify the file name of comupation results
FILENAME_RESULT = 'a.dat'
# 3. Specify 'True' of 'False' to apply Periodic Boundary Condition
PBC = False
# 4. Specify the computation Mode
"""
Use 'A' for Angles, 'D' for Distances, 'C' for Coordinations
"""
MODE = 'A'
# 5. Specify the computation subject
"""
For Angles: 'EdgeAtom1 CenterAtom EdgeAtom2 Cutoff1 Cutoff2 Binsize'
Example: SUBJECT = 'Cd Se Cd 3.0 3.0 2.0'

For Distances: 'PairAtom1 PairAtom2 Cutoff'
Example: SUBJECT = 'Cd Se 3.0'

For Coordinations: 'PairAtom1 PairAtom2 Cutoff'
Example: SUBJECT = 'Cd Se 3.0'
"""
SUBJECT = 'Se Cd Se 3.0 3.0 5.0'

#END OF GLOBAL SETTINGS


#WARNING:DO NOT MODIFY ANY PART BELOW!!!


# A: Basic functions of loading files
STR2FLOAT = lambda str_list: list(map(float, str_list))


def file_load(filename):
    """Load POSCAR/CONTCAR files"""
    with open(filename, 'r') as file_r:
        return [line.split() for line in file_r if line.strip()]


def get_shape(file):
    """Get the box shape of the system"""
    #Get the scale factor
    scale = float(*file[1])
    # Get the SHAPE matrix
    shape = map(STR2FLOAT, file[2:5])
    return np.array([np.linalg.norm(vec) * scale for vec in shape])


def get_typelist(file):
    """Get the typelist"""
    atom_type = file[5]
    atom_number = [int(num) for num in file[6]]
    
    idx_f = np.cumsum(atom_number)
    idx_i = idx_f - np.array(atom_number)
    idx_range = [list(range(i, f)) for i, f in zip(idx_i, idx_f)]
    
    typelist = dict(zip(atom_type, idx_range))
    return typelist


def get_xyz(file, box):
    """Get the coordinates"""
    total_num = sum(int(num) for num in file[6])
    record = file[7][0]
    # Get atom positions
    xyz = np.array(list(map(STR2FLOAT, file[8:8+total_num])))
    if record[0] == 'D':
        xyz = xyz * box
    return xyz


# Basic function concerning computation
def mirdist(r_1, r_2, box, pbc):
    """Compute the pair-distance between r_1 and r_2"""
    v_12 = r_1 - r_2 - np.round((r_1 - r_2) / box) * box * pbc
    return np.linalg.norm(v_12)


def check_bond(r_1, r_2, cutoff, box, pbc):
    """Check whether two atom are bonded"""
    return 0 < mirdist(r_1, r_2, box, pbc) <= cutoff


def mirangle(r_1, r_2, r_3, box, pbc):
    """compute the angle formed as r_1-r_2-r_3"""
    v_21 = r_2 - r_1 - np.round((r_2 - r_1) / box) * box * pbc
    v_23 = r_2 - r_3 - np.round((r_2 - r_3) / box) * box * pbc
    rad = np.arccos(1 - pdist(np.vstack([v_21, v_23]), 'cosine'))
    deg = float(rad * 180 / np.pi)
    return deg


def compute_pdist(subject, xyz, typelist, box, pbc):
    """Compute Distance"""
    type1, type2, cutoff = subject
    cutoff = float(cutoff)
    pairs = [(idx1, idx2) for
             idx1 in typelist[type1] for idx2 in typelist[type2]]
    values = [mirdist(xyz[idx1], xyz[idx2], box, pbc) for
              idx1 in typelist[type1] for idx2 in typelist[type2]]

    num = len(values)

    idx_legal = [i for i in range(num) if 1.0e-6 < values[i] <= cutoff]
    pairs = [pairs[idx] for idx in idx_legal]
    values = [values[idx] for idx in idx_legal]
    return pairs, values


def compute_coord(subject, xyz, typelist, box, pbc):
    """Compute coordination number"""
    type1, type2, cutoff = subject
    cutoff = float(cutoff)
    atoms = typelist[type1]
    values = [sum(check_bond(xyz[idx1], xyz[idx2], cutoff, box, pbc) for
                  idx2 in typelist[type2]) for idx1 in typelist[type1]]
    return atoms, values


def compute_angle(subject, xyz, typelist, box, pbc):
    """Compute Angle"""
    type1, type2, type3, cut21, cut23, _ = subject
    cut21, cut23 = float(cut21), float(cut23)
    centers, pairs, values = [], [], []
    for idx2 in typelist[type2]:
        list1 = [idx1 for idx1 in typelist[type1] if
                 check_bond(xyz[idx1], xyz[idx2], cut21, box, pbc)]
        list3 = [idx3 for idx3 in typelist[type3] if
                 check_bond(xyz[idx3], xyz[idx2], cut23, box, pbc)]
        if not (list1 and list3):
            continue
        centers.append(idx2)
        pair = [sorted([idx1, idx3]) for idx1 in list1 for idx3 in list3]
        value = [mirangle(xyz[idx1], xyz[idx2], xyz[idx3], box, pbc) for
                 idx1 in list1 for idx3 in list3]
        pairs.append(pair)
        values.append(value)
    clean_angle(centers, pairs, values)
    refine_angle(centers, pairs, values)
    return centers, pairs, values


def count_angle(values, binsize):
    """Count angles"""
    angles = list(itertools.chain(*values))
    angles = [round(angle,3) for angle in angles]
    angles.sort()
    histc = np.histogram(angles, bins=int(180/binsize), range=(0,180))
    return Counter(angles), histc


def dedupe(pairs):
    """Get the indexes of unrepeted pairs"""
    num = len(pairs)
    seen = []
    for idx in range(num):
        if pairs[idx] not in seen:
            yield idx
            seen.append(pairs[idx])

def nonzero(values):
    """Get the idexes of nonzeros or unempty values"""
    num = len(values)
    return (idx for idx in range(num) if values[idx] > 1e-5)

def unempty(values):
    """return unempty indexes"""
    num = len(values)
    return (idx for idx in range(num) if values[idx])


def clean_angle(centers, pairs, values):
    """Remove zero values and repeated pairs"""
    num = len(centers)
    for idx in range(num):
        idx_legal1 = list(dedupe(pairs[idx]))
        idx_legal2 = list(nonzero(values[idx]))
        idx_legal = [idx for idx in idx_legal1 if idx in idx_legal2]
        pairs[idx] = [pairs[idx][i] for i in idx_legal]
        values[idx] = [values[idx][i] for i in idx_legal]

def refine_angle(centers, pairs, values):
    """Remove empty groups"""
    num = len(centers)
    idx_unempty = list(unempty(values))
    centers = [centers[idx] for idx in range(num) if idx in idx_unempty]
    pairs = [pairs[idx] for idx in range(num) if idx in idx_unempty]
    values = [values[idx] for idx in range(num) if idx in idx_unempty]


# Basic function concerning file writing
def write_pdist(file_w, pairs, values):
    """Write the distance information"""
    file_w.write('Atom1\tAtom2\tDistance\n')
    for pair, value in zip(pairs, values):
        file_w.write('{:<4d}\t{:<4d}\t'.format(pair[0]+1, pair[1]+1))
        file_w.write('{:<10f}\n'.format(value))

def write_angle(file_w, centers, pairs, values):
    """Write the angle information"""
    file_w.write('AtomE1\tAtomC\tAtomE2\tAngle\n')
    num = len(centers)
    for idx in range(num):
        for pair, value in zip(pairs[idx], values[idx]):
            file_w.write('{:<4d}\t'.format(pair[0]+1))
            file_w.write('{:<4d}\t'.format(centers[idx]+1))
            file_w.write('{:<4d}\t'.format(pair[1]+1))
            file_w.write('{:<10.3f}\n'.format(value))
    file_w.write('END\n')


def write_count_angle(file_w, values, binsize):
    """Write the histogram of angles"""
    file_w.write('\nAngle\tCount\n')
    counts, histc = count_angle(values, binsize)
    for angle, count in counts.items():
        file_w.write('{:<10.3f}\t{:<4d}\n'.format(angle, count))
    file_w.write('END\n')
    file_w.write('\nCenterAngle\tHistCount\n')
    angle_centers = histc[1][1:] - binsize/2
    for angle, count in zip(angle_centers, histc[0]):
        file_w.write('{:<6.3f}\t{:<4d}\n'.format(angle, count))
    file_w.write('END\n')


def write_coord(file_w, atoms, values):
    """Write coordination number"""
    for atom, value in zip(atoms, values):
        file_w.write('{:<4d}\t{:<4d}\n'.format(atom+1, value))


# Load

FILE = file_load(FILENAME)
BOX = get_shape(FILE)
TYPELIST = get_typelist(FILE)
XYZ = get_xyz(FILE, BOX)
SUBJECT = SUBJECT.split()

# Computation
if MODE == 'A':
    CENTERS, PAIRS, VALUES = compute_angle(SUBJECT, XYZ, TYPELIST,
                                           BOX, PBC)
elif MODE == 'D':
    PAIRS, VALUES = compute_pdist(SUBJECT, XYZ, TYPELIST, BOX, PBC)

elif MODE == 'C':
    ATOMS, VALUES = compute_coord(SUBJECT, XYZ, TYPELIST, BOX, PBC) 

else:
    print("Mode not recognised, use 'A', 'D' or 'C' only!")


# Dump
with open(FILENAME_RESULT, 'w') as file_w:
    if MODE == 'A':
        write_angle(file_w, CENTERS, PAIRS, VALUES)
        write_count_angle(file_w, VALUES, float(SUBJECT[5]))
    elif MODE == 'D':
        write_pdist(file_w, PAIRS, VALUES)
    elif MODE == 'C':
        write_coord(file_w, ATOMS, VALUES)
    else:
        print("None Result")
    file_w.write('\nEND')
print('COMPUTATION FINISHED')
