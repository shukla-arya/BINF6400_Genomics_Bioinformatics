#!/usr/bin/env python
# coding: utf-8

# # **PSet4 Trees and Models**
# 
# **Genomics in Bioinformatics**
# <br>
# **Arya Shukla**  
# **3 April 2025**

# **Global Imports**

# In[1]:


import numpy as np
import random


# **Class Assignment**

# In[2]:


distanceMatrix = [[0 ,12,12,13,15,15],
                  [12, 0, 2, 6, 8, 8],
                  [12, 2, 0, 6, 9, 9],
                  [13, 6, 6, 0, 8, 8],
                  [15, 8, 9, 8, 0, 4],
                  [15, 8, 9, 8, 4, 0]]

speciesList = ["M_Spacii", "T_Pain", "G_Unit", "Q_Doba", "R_Mani", "A_Finch"]


# In[3]:


def UPGMA(dM, sp):
    # Open the file in write mode to save the results
    with open("Arya_Shukla_UPGMA.txt", "w") as f:
        while len(dM) > 1:
            leastRow, leastCol = findSmallest(dM)
            cluster_dist = dM[leastRow][leastCol] / 2
            dM = updateMatrix(dM, leastRow, leastCol)
            sp = updateSpecies(sp, leastRow, leastCol, cluster_dist)
            
            # Write the iteration result to the file
            f.write("\nBegin Iteration.\n")
            f.write(str(dM) + "\n")
            f.write(str(sp) + "\n")


# In[4]:


def findSmallest(dM):
    min_val = float('inf')
    row, col = -1, -1
    
    for i in range(len(dM)):
        for j in range(i):  # only check lower triangle of symmetric matrix
            if 0 < dM[i][j] < min_val:
                min_val = dM[i][j]
                row, col = i, j
    
    return row, col


# In[5]:


def updateSpecies(sp, r, c, d):
    sp[r] = f"({sp[r]}:{d},{sp[c]}:{d})"
    del sp[c]
    return sp


# In[6]:


def updateMatrix(dM, row, col):
    newMat = [] # new matrix that is smaller in size than 
    
    for i in range(len(dM)):
        if i != row and i != col:
            new_row = []
            for j in range(len(dM)):
                if j != row and j != col:
                    new_row.append(dM[i][j])
            
            # Compute new distance as average of row and col distances
            new_distance = (dM[i][row] + dM[i][col]) / 2
            new_row.append(new_distance)
            newMat.append(new_row)
    
    # Add last row for the new cluster
    newMat.append([newMat[i][-1] for i in range(len(newMat))] + [0])
    
    return newMat


# In[7]:


UPGMA(distanceMatrix, speciesList)


# **Discussion**

# UPGMA is a greedy clustering algorithm. Once we create a grouping, we do not need to worry about those again. If the updated distances result in 2+ identical values and small, then we can select randomly.
