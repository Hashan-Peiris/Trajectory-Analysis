# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 13:55:55 2021

@author: grip1
"""

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

filename = 'All.xyz'
AtmA = [57, 58, 59, 60, 49, 50, 51, 52, 57, 58, 59, 60, 45, 46, 47, 48, 37, 38, 39, 40, 45, 46, 47, 48, 33, 34, 35, 36,
        25, 26, 27, 28, 33, 34, 35, 36, 13, 14, 15, 16, 17, 18, 19, 20, 13, 14, 15, 16, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3,
        4]
AtmB = [49, 50, 51, 52, 53, 54, 55, 56, 53, 54, 55, 56, 37, 38, 39, 40, 41, 42, 43, 44, 41, 42, 43, 44, 25, 26, 27, 28,
        29, 30, 31, 32, 29, 30, 31, 32, 17, 18, 19, 20, 21, 22, 23, 24, 21, 22, 23, 24, 5, 6, 7, 8, 9, 10, 11, 12, 9,
        10, 11, 12]
multiply = 1  # To compensate for the no of timesteps skipped when converting

xyz = open(filename, 'r')
# coord_rec = open("AtmDist.txt", 'w')

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

Distance_List = []  # stores the distances in each time step from the elements in AtmA alternatively
# particular_line_A = linecache.getline(filename, 532)
# print(particular_line_A)

for timestep in range(int(timesteps)):
    for index in range(len(AtmA)):
        # separately gets the respective lines for the corresponding elements in the AtmA and AtmB
        line_position_A = ((timestep * (totnum + 2)) + (AtmA[index] + 2))
        particular_line_A = linecache.getline(filename, int(line_position_A))

        line_position_B = ((timestep * (totnum + 2)) + (AtmB[index] + 2))
        particular_line_B = linecache.getline(filename, int(line_position_B))
        if timestep % 1 == 0:
            print("Timestep:", timestep, "Index:", index, "Line_PosA:", line_position_A)

        print(AtmA[index], AtmB[index])
        print("This:", particular_line_A, particular_line_B)
        print("Timestep:", timestep, "Index:", index, "Line_PosB:", line_position_B)

        # stores the [index H x y z] in a list
        A = [];
        B = []
        A = particular_line_A.rstrip('\n').split('"')[0].split()
        B = particular_line_B.rstrip('\n').split('"')[0].split()
        print("A:", A)
        print("B:", B, "\n")
        linecache.clearcache()

        # to get only the x y z values
        A1 = [];
        B1 = []
        for i in range(len(A) - 2): A1.append(float(A[i + 2]))
        for i in range(len(B) - 2): B1.append(float(B[i + 2]))
        print("A1:", A1)
        print("B1:", B1, "\n")

        # CHECKING FOR NON-ORTHOGONAL CELLS
        check = [1, 2, 3, 5, 6, 7]
        for i in check:
            if abs(lattice[i]) == 0.0:
                continue
            else:
                print("THIS IS NOT AN ORTHOGONAL CELL!")
                xyz.close()
                # wfile.close
                break

        # checking for saddling across boundary
        if abs(abs(A1[0]) - abs(B1[0])) < abs(0.5 * lattice[0]):
            X = abs(abs(A1[0]) - abs(B1[0]))

        else:
            X = abs(abs(abs(A1[0]) - abs(B1[0])) - abs(lattice[0]))

        if abs(abs(A1[1]) - abs(B1[1])) < abs(0.5 * lattice[4]):
            Y = abs(abs(A1[1]) - abs(B1[1]))
        else:
            Y = abs(abs(abs(A1[1]) - abs(B1[1])) - abs(lattice[4]))

        if abs(abs(A1[2]) - abs(B1[2])) < abs(0.5 * lattice[8]):
            Z = abs(abs(A1[2]) - abs(B1[2]))
        else:
            Z = abs(abs(abs(A1[2]) - abs(B1[2])) - abs(lattice[8]))

        # print (X,Y,Z,"\n")

        Distance_List.append(np.sqrt(X ** 2 + Y ** 2 + Z ** 2))
        # print("DistList:" , Distance_List, "\n", "--------------------------------", "\n")

        # line=xyz.readline()
        # print(line)

stored = len(Distance_List)
list_items = len(AtmA)
skip = int(stored / list_items)
print(stored, list_items, skip)
counter = 0

# This is for the names of list of pairs
Pair_list = []
for i in range(len(AtmA)):
    Pair_list.append([AtmA[i], AtmB[i]])
# print("Pair_List:", Pair_list)

# This is for the distance of a particular pair in a sublist ex: {[Dist of A1]][Dist of A2]}
Pair_distances = []
for i in range(len(AtmA)):
    Pair_distances.append([AtmA[i], AtmB[i]])
# print("Pair_List:", Pair_distances)

# This replaces the pair values in List_Distances list by the relevant list of distances over all timesteps
for loop in range(len(AtmA)):
    i = []
    for j in range(int(timesteps)):
        print(j)
        i.append(Distance_List[counter])
        counter = counter + list_items
        # print(counter)
        # Pair="A_%i_B_%i" % (AtmA[i], AtmB[i])
    # print("i is:", i)
    Pair_distances[loop] = i
    counter = loop + 1
# print("Pair_List:", Pair_list)
# print("Pair Distances:", Pair_distances)

with open("ORDERED_LIST.txt", 'w') as wfile:
    wfile.write("Time\t")
    for item in Pair_list:
        wfile.write(str(item) + "\t")
    wfile.write("\n")

    for ii in range(int(timesteps)):
        Number_List = Pair_list
        print("WRITING TO FILE")
        print(multiply * (ii + 1))

        wfile.write(str(multiply * (ii + 1)) + "\t")
        for current_line in range(len(Pair_distances)):
            wfile.write(str(Pair_distances[current_line][ii]) + "\t")
        wfile.write("\n")

        # print(Pair_distances[current_line][ii])
        # print("%i\n" %(current_line))

xyz.close()
wfile.close()
