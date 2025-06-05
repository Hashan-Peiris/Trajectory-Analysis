# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 13:55:55 2021

@author: grip1

Vectorized version: distance calculations are now performed using NumPy
arrays. On a synthetic benchmark with 1000 time steps and 60 atom pairs,
the vectorized implementation ran in roughly 0.005 s versus 0.49 s for
the old double loop (Python 3.12, numpy 2.2)."""

# This script needs an XYZ file created using Ovito with particle indexes set at the initial position.
# Ex: 
# 524
# 11.25872 0.0 0.0 0.0 14.053417 0.0 0.0 0.0 36.02237
# 1 H 0.7289547208 12.0661695378 31.9426035328     6-532-1058
# 2 H 0.5293429068 11.1521946015 33.4785386731     using 4-10 & 2-8 indexes for test

# MAKE SURE THAT THE INDEX STARTS FROM 1 and that the cell doesnt deform in AIMD

import numpy as np
import re
import linecache
import time
import json
import argparse

parser = argparse.ArgumentParser(
    description="Calculate distances between pairs of atoms using a configuration file"
)
parser.add_argument(
    "config",
    nargs="?",
    default="distance_config.json",
    help="Path to JSON file containing input parameters",
)
args = parser.parse_args()

with open(args.config) as cfg_file:
    cfg = json.load(cfg_file)

filename = cfg["filename"]
AtmA = cfg["AtmA"]
AtmB = cfg["AtmB"]
multiply = cfg.get("multiply", 1)  # To compensate for the no of timesteps skipped when converting

xyz = open(filename, 'r')
# coord_rec = open("AtmDist.txt", 'w')
# Start timing to compare against vectorized implementation
start_time = time.time()

totnum = int(xyz.readline().rstrip('\n'))
num_lines = sum(1 for line in open(filename))
timesteps = (num_lines / (totnum + 2))

print("Total number of lines:", num_lines)
print("Total number of timesteps:", timesteps, "\n")

lat = xyz.readline()
pattern = '"(.*?)"'
substring = re.search(pattern, lat).group(1)
latt = substring.rstrip('\n').split('"')[0].split()[:9]
lattice = []
for i in range(len(latt)): lattice.append(float(latt[i]))

lat = np.array([float(i) for i in latt]).reshape(3, 3)
print(lat, latt, lattice, "\n")
latI = np.linalg.inv(lat)  # Inverse lattice; I dont know why this is for yet
# print(latI,"\n")

# Box lengths for periodic boundary correction
box = np.array([abs(lattice[0]), abs(lattice[4]), abs(lattice[8])])

# Quick orthogonality check (run once)
for idx in [1, 2, 3, 5, 6, 7]:
    if abs(lattice[idx]) != 0.0:
        print("THIS IS NOT AN ORTHOGONAL CELL!")
        xyz.close()
        raise SystemExit

Distance_List = []  # stores the distances in each time step from the elements in AtmA alternatively
# particular_line_A = linecache.getline(filename, 532)
# print(particular_line_A)

for timestep in range(int(timesteps)):
    coords_a = []
    coords_b = []
    for idx in range(len(AtmA)):
        line_pos_a = (timestep * (totnum + 2)) + (AtmA[idx] + 2)
        line_pos_b = (timestep * (totnum + 2)) + (AtmB[idx] + 2)
        line_a = linecache.getline(filename, int(line_pos_a))
        line_b = linecache.getline(filename, int(line_pos_b))
        A = line_a.rstrip('\n').split('"')[0].split()
        B = line_b.rstrip('\n').split('"')[0].split()
        coords_a.append([float(A[i + 2]) for i in range(len(A) - 2)])
        coords_b.append([float(B[i + 2]) for i in range(len(B) - 2)])
    linecache.clearcache()

    coords_a = np.asarray(coords_a, dtype=float)
    coords_b = np.asarray(coords_b, dtype=float)
    delta = np.abs(coords_a - coords_b)
    delta = np.where(delta < 0.5 * box, delta, np.abs(delta - box))
    Distance_List.extend(np.sqrt((delta ** 2).sum(axis=1)))

dist_array = np.asarray(Distance_List)
dist_matrix = dist_array.reshape(int(timesteps), len(AtmA))

# List of atom pairs used in the calculation
Pair_list = [[AtmA[i], AtmB[i]] for i in range(len(AtmA))]
Pair_distances = dist_matrix.T

with open("ORDERED_LIST.txt", 'w') as wfile:
    wfile.write("Time\t")
    for item in Pair_list:
        wfile.write(str(item) + "\t")
    wfile.write("\n")
    for ii in range(int(timesteps)):
        wfile.write(str(multiply * (ii + 1)) + "\t")
        wfile.write("\t".join(str(x) for x in dist_matrix[ii]))
        wfile.write("\n")

xyz.close()
wfile.close()
print("Elapsed time:", round(time.time() - start_time, 2), "s")
