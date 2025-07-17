#!/usr/bin/env python
# coding: utf-8

# # **PSet5 HMMs and Models**
# 
# **Genomics in Bioinformatics**
# <br>
# **Arya Shukla**  
# **17 April 2025**

# **Global Imports**

# In[1]:


import numpy as np
import random
import matplotlib.pyplot as plt
import math


# **Class Assignment**

# *Nucleotide Composition HMM*

# In[2]:


baseIDx = {"A":0,"C":1, "G":2,"T":3}


# In[3]:


def main():
    spaciiFA = "MSpacii.fa"
    pathogenFA = "pathogen.fa"
    spaciiFA_T = "MSpacii_training.fa"
    pathogenFA_T = "pathogen_training.fa"
    spaciiID2seq = getSeq(spaciiFA)
    pathogenID2seq = getSeq(pathogenFA)
    
    spaciiTrainModel = [[0, 0, 0, 0] for _ in range(4)]
    pathTrainModel = [[0, 0, 0, 0] for _ in range(4)]
    
    spaciiTrainModel = trainModel(spaciiTrainModel, spaciiFA_T)
    pathTrainModel = trainModel(pathTrainModel, pathogenFA_T)
    
    markovScoresSpacii = []
    markovScoresPath = []

    for ID in spaciiID2seq.keys():
        markovScoresSpacii.append(getLogLike(spaciiTrainModel, pathTrainModel, spaciiID2seq[ID]))

    for ID in pathogenID2seq.keys():
        markovScoresPath.append(getLogLike(spaciiTrainModel, pathTrainModel, pathogenID2seq[ID]))
    
    ####----------------------output-------------------------
    plt.hist([markovScoresPath, markovScoresSpacii], bins=20, label=['pathogen', 'spacii'], rwidth=1, density=True)
    plt.legend()
    plt.title("Log-Likelihood Score Distribution")
    plt.xlabel("Log-Likelihood Score")
    plt.ylabel("Density")
    plt.savefig("score_distribution.png")
    plt.show()
    
    scoresOutputText(markovScoresSpacii, markovScoresPath)
    ####----------------------output-------------------------

def scoresOutputText(markovScoresSpacii, markovScoresPath):
    with open("results.tab", "w") as f:
        f.write("SpaciiScores\tpathogenScores\n")
        for i in range(len(markovScoresSpacii)):
            s = markovScoresSpacii[i] if i < len(markovScoresSpacii) else ""
            p = markovScoresPath[i] if i < len(markovScoresPath) else ""
            f.write(f"{s}\t{p}\n")

# Compares how likely a given test sequence is under each model
def getLogLike(model1, model2, seq):
    # Log-likelihood ratio
    loglike = 0.0

    for i in range(1, len(seq)):
        prev = seq[i - 1]
        curr = seq[i]
        if prev in baseIDx and curr in baseIDx:
            prev_idx = baseIDx[prev]
            curr_idx = baseIDx[curr]

            p1 = model1[prev_idx][curr_idx]
            p2 = model2[prev_idx][curr_idx]

            # Avoid log(0) by applying a very small probability
            if p1 == 0:
                p1 = 1e-10
            if p2 == 0:
                p2 = 1e-10

            loglike += math.log(p1 / p2)

    return loglike

def trainModel(model, fasta_file):
    seqs = getSeq(fasta_file)
    counts = [[0 for _ in range(4)] for _ in range(4)]

    for seq in seqs.values():
        for i in range(1, len(seq)):
            prev = seq[i - 1]
            curr = seq[i]
            if prev in baseIDx and curr in baseIDx:
                prev_idx = baseIDx[prev]
                curr_idx = baseIDx[curr]
                counts[prev_idx][curr_idx] += 1

    # Convert counts to probabilities
    for i in range(4):
        total = sum(counts[i])
        if total == 0:
            model[i] = [0.25] * 4  # uniform if no data
        else:
            model[i] = [c / total for c in counts[i]]

    print(model)
    print()
    return model

def getSeq(filename):
    f = open(filename)
    id2seq = {}
    currkey = ""
    for line in f:
        if line.startswith(">"):
            currkey = line.rstrip()[1:]
            id2seq[currkey] = ""
        else:
            id2seq[currkey] += line.rstrip()
    f.close()
    return id2seq

main()


# **Discussion**

# - The goal is to classify sequences as originating from Spacii or the pathogen using their nucleotide patterns.
# - Dinucleotide transitions are point mutation where a purine nucleotide (A or G) changes to another purine, or a pyrimidine nucleotide (C or T) changes to another pyrimidine.
# - If the log-likelihood is positive, the sequence is more likely to be from Spacii. If negative, it's more likely from the pathogen.
# - Each row of the 4x4 matrix sums to 1 because they are conditional probabilities. It represents how one base tends to follow another in an organism’s genome. A higher value means that a dinucleotide is more frequent.
