# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 18:06:04 2021

@author: grip1
"""

# This script needs an XYZ file created using Ovito with particle indexes set at the initial position.
# Ex: 
    
# 1048
# Lattice="11.25872 0.0 0.0 0.0 14.053417 0.0 0.0 0.0 72.04474" Properties=id:I:1:species:S:1:pos:R:3
# 1 H 0.4209823429 11.9669238853 32.3066153452
# 2 H 0.7437714215 11.3166162076 33.9709691539


import numpy as np
import re
import linecache
import sys
import time as tm

start = tm.time()
# ------------------------------------------------------------------------------------------
filename='Distances5.xyz'
AtmA=[380,384,389,396,407,409,414,420,423,427,432,438,446,450]
AtmB=["O","O"]
Range=8.0    # Max radius
Limit=0.05   # Bin size
Type=["All","O","C","Li","H","B","F","Ni"]     # To get the RDF relative to all elements or just some of them. Type "All" or one element name
multiply=25  # To compensate for the no of timesteps skipped when converting. Not used here.
# ------------------------------------------------------------------------------------------

xyz=open(filename,'r')
totnum=int(xyz.readline().rstrip('\n'))
num_lines = sum(1 for line in open(filename)) #Total Number of Lines
timesteps = (num_lines/(totnum+2))

discarded=0
counted=0
counter=0

print("Total number of lines:", num_lines)
print("Total number of timesteps:", timesteps)
print("  Range:%f\n  Limit:%f" %(Range,Limit) ,"\n")

lat=xyz.readline()
pattern = '"(.*?)"'
substring = re.search(pattern, lat).group(1)
latt=substring.rstrip('\n').split('"')[0].split()[:9]
lattice=[]
for i in range (len(latt)): lattice.append(float(latt[i]))

lat=np.array([float(i) for i in latt]).reshape(3,3)
#print (lat, latt, lattice, "\n")
latI=np.linalg.inv(lat)  # Inverse lattice; I dont know why this is for yet
#print(latI,"\n")

#CHECKING FOR NON-ORTHOGONAL CELLS
check=[1,2,3,5,6,7]
for i in check:
    if abs(lattice[i]) == 0.0:
        continue
    else:
        print("THIS IS NOT AN ORTHOGONAL CELL!")
        xyz.close 
        #wfile.close
        sys.exit()

# ------------------------------------------------------------------------------------------
# fname="XX"
# def test(fname):
#     if fname == "XXX" : 
#         print ("GOOD")
#     else:
#         #print ("BADDD")
#         return "XXX"
# print (test(fname))
# if test(fname) == "XXX":
#     print("WEELL")
# else:
#     print("NOOOO")
# ------------------------------------------------------------------------------------------    
    
## ---------------------------------------------------------------------------------------------------- ##
## ---------------------------------------------------------------------------------------------------- ##

# A function to get the details in a line from a index number
# Takes in Timestep point, index(no of items in AtmA)
def Read_Line(timestep, indexA):
        #print("STARTING RLA")
        #separately gets the respective lines for the corresponding elements in the AtmA
        line_position_A=((timestep*(totnum+2))+(indexA+2))
        particular_line_A = linecache.getline(filename, int(line_position_A))


        #stores the [index ele_name x y z] in a list
        A=[]
        A=particular_line_A.rstrip('\n').split('"')[0].split()
        #print("A:", A)
        #linecache.clearcache()        
        #sys.exit()
        #print("RETURNING RLA\n")
        return A


# A function to get the details in a line from ALLL index numbers
# Takes in Timestep point, index(specified in the command)
def Read_Line_B(timestep, indexB):
        #separately gets the respective lines for all corresponding atoms
        #print("STARTING RLB")
        line_position_B=((timestep*(totnum+2))+(indexB+2))
        particular_line_B = linecache.getline(filename, int(line_position_B))
        #if timestep % 1 == 0:
            #print ("Timestep:", timestep, "IndexB:", indexB, "Line_PosB:", line_position_B)

        #stores the [index ele_name x y z] in a list
        B=[]
        B=particular_line_B.rstrip('\n').split('"')[0].split()
        #print("B:", B)
        #print("RETURNING RLB\n")
        #linecache.clearcache()        
        return B
        
    
# A function to get the xyz values out of A: ['427', 'O', '10.293610137', '8.117507745', '19.6899517192']
def get_XYZ(A): 
    #print("STARTING G_XYZ")      
    A1=[]
    for i in range (len(A)-2): A1.append(float(A[i+2]))
    #print("A1:", A1)
    #print("RETURNING G_XYZ\n")   
    return A1
    #print("B1:", B1,"\n")
    
def get_Element(A): 
    #print("STARTING Element Check")      
    Element=A[1]
    return str(Element)
    
# A function to get the Distance value from XYZ distances - after get_XYZ    
def get_Distance(A1,B1):
    #print("STARTING GD")
    DistList=[]
    #checking for saddling across boundary & get X Y Z distances  
    if abs(abs(A1[0])-abs(B1[0])) < abs(0.5*lattice[0]):
        X=abs(abs(A1[0])-abs(B1[0]))
    else:
        X=abs(abs(abs(A1[0])-abs(B1[0]))-abs(lattice[0]))
    if abs(abs(A1[1])-abs(B1[1])) < abs(0.5*lattice[4]):
        Y=abs(abs(A1[1])-abs(B1[1]))
    else:
        Y=abs(abs(abs(A1[1])-abs(B1[1]))-abs(lattice[4]))
    if abs(abs(A1[2])-abs(B1[2])) < abs(0.5*lattice[8]):
        Z=abs(abs(A1[2])-abs(B1[2]))
    else:
        Z=abs(abs(abs(A1[2])-abs(B1[2]))-abs(lattice[8]))

    
    DistList.append(np.sqrt(X**2+Y**2+Z**2)) # Appends this to a DistList (not used YET)
    Distance =  np.sqrt(X**2+Y**2+Z**2)
    
    #print (X,Y,Z) # Displacement in X Y Z directions
    #print (Distance)
    #print("RETURNING GD\n")
    return Distance


# # A function to get and assign the number of atoms in each bin for one time step   
# def get_RDF(timestep, indexA, indexB):
#     #print("STARTING gRDF")

#     dist=get_Distance(Line_A,get_XYZ(Read_Line_B(timestep,indexB)))
#     #print(dist)
#     for item in KeysList:
#         #counter+=1
#         if  0.000000  == dist > KeysList[-1]:
#             #print("SAME ATOM IndexA:%s IndexB:%s" %(indexA,indexB))
#             #discarded=discarded+1
#             break
#         elif 0.000000 < dist and item > dist:
#             #print ("WRITING Item:"+str(item))
#             #print (dicSteps[item])
#             dicSteps[item] = int(dicSteps[item]+1)
#             #counted+=1
#             break
#         else:
#             pass
#             #discarded=discarded+1
#             #break
    
#     #print("RETURNING gRDF\n")
    
    
    
    # A function to get and assign the number of atoms in each bin for one time step   
def get_RDF(timestep, indexA, indexB):
    #print("STARTING gRDF")
    
    #Ele1=get_Element(Read_Line(timestep, indexA))
    Ele2=get_Element(Read_Line(timestep, indexB))
    
    if element == Ele2 or element == "All":
        #print("Matching %s" %element)
        print ("Timestep:", timestep, "IndexA:", indexA, "Line_PosA:", indexA, "MATCHING", element)
        dist=get_Distance(Line_A,get_XYZ(Read_Line_B(timestep,indexB)))
        #print(dist)
        for item in KeysList:
            #counter+=1
            if  0.000000  == dist > KeysList[-1]:
                #print("SAME ATOM IndexA:%s IndexB:%s" %(indexA,indexB))
                #discarded=discarded+1
                break
            elif 0.000000 < dist and item > dist:
                #print ("WRITING Item:"+str(item))
                #print (dicSteps[item])
                dicSteps[item] = int(dicSteps[item]+1)
                #counted+=1
                break
            else:
                pass
                #discarded=discarded+1
                #break
    else:
        print ("Timestep:", timestep, "IndexA:", indexA, "Line_PosA:", indexA, "NOT MATCHING", element)
        pass
        
        #print("RETURNING gRDF\n")
## ---------------------------------------------------------------------------------------------------- ##
## ---------------------------------------------------------------------------------------------------- ##

items=int(Range/Limit)
print ("iTEMS:"+str(items))  # Number of items in the KeysList

# Dictionary Comprehension; stores keys with intial 
#dicSteps = {round((x*Limit+Limit),2): 0 for x in range(items)} 
#NewDicSteps = dicSteps

for element in Type:
    rang=totnum  #320
    time=timesteps
    
    for indexA in AtmA:
        dicSteps = {round((x*Limit+Limit),2): 0 for x in range(items)}
        KeysList=[*dicSteps]
        
        with open("%i %s.txt" %(indexA,element),'w') as wfile:
            wfile.write("Bin\tIndex_%i\t\n" %indexA)        
            
            for timestep in range(int(time)):
                Line_A = get_XYZ(Read_Line(timestep, indexA))
    
                for i in range(rang):
                    get_RDF(timestep,indexA,i+1) 
                
            for key , value in dicSteps.items(): 
                #print("WRITING TO FILE"+str(i))
                print (key,value)
                wfile.write('%s\t%s\t\n' % (key, value)) 
                
            
    FileNameList=[]
    for IndexA in AtmA: FileNameList.append("%i %s.txt" %(IndexA,element))
    with open('COMBINED_VALUES_%s.txt' %element, 'w') as writer:
        readers = [open(filename) for filename in FileNameList]
        for lines in zip(*readers):
            print('\t'.join([line.strip() for line in lines]), file=writer)
        
end = tm.time()
TIME=end - start
print("TIME ELAPSED: " + str(TIME/60))        
        
        
    





































