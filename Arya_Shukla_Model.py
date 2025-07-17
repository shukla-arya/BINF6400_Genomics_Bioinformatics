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


import math
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc


# **Class Assignment**

# *Part One*

# In[2]:


modelCodons = ['TTT', 'TTC', 'TTA', 'TTG', 'CTT', 'CTC', 'CTA',
               'CTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTT', 'GTC',
               'GTA', 'GTG', 'TCT', 'TCC', 'TCA', 'TCG', 'AGT',
               'AGC', 'CCT', 'CCC', 'CCA', 'CCG', 'ACT', 'ACC',
               'ACA', 'ACG', 'GCT', 'GCC', 'GCA', 'GCG', 'TAT',
               'TAC', 'CAT', 'CAC', 'CAA', 'CAG', 'AAT', 'AAC',
               'AAA', 'AAG', 'GAT', 'GAC', 'GAA', 'GAG', 'TGT',
               'TGC', 'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG',
               'GGT', 'GGC', 'GGA', 'GGG', 'TGG', 'TAA', 'TAG',
               'TGA']  ############# This list is the codon labels to the coding and non-coding models


# In[3]:


def scoreModels():
    codingMatrix = getProbs("./codingModel.tab")  # 64x64 matrix of probabilities for coding triplet mutations
    noncodingMatrix = getProbs("./noncodingModel.tab")  # 64x64 matrix of probabilities for non-coding triplet mutations
    
    id2ancestorSeq = getSeq("./Ancestor.fa")  # Reads the ancestor sequences
    id2spaciiSeq = getSeq("./Spacii.fa")  # Reads the M. Spacii sequences
    allID = list(id2ancestorSeq.keys())  # The two sequences above share dictionary indexes

    # Open file for writing results
    with open("Arya_Shukla_Model.txt", "w") as result_file:
        for ID in allID:
            cScore = 0  # Log-likelihood score for the coding model
            nScore = 0  # Log-likelihood score for the non-coding model

            ancestorSeq = id2ancestorSeq[ID]
            spaciiSeq = id2spaciiSeq[ID]

            for i in range(0, len(ancestorSeq), 3):
                ancestorCodon = ancestorSeq[i:i+3]
                spaciiCodon = spaciiSeq[i:i+3]

                if ancestorCodon in modelCodons and spaciiCodon in modelCodons:
                    ancestorIndex = modelCodons.index(ancestorCodon)
                    spaciiIndex = modelCodons.index(spaciiCodon)

                    cProb = codingMatrix[ancestorIndex][spaciiIndex]
                    nProb = noncodingMatrix[ancestorIndex][spaciiIndex]

                    if cProb > 0:
                        cScore += math.log(cProb)
                    if nProb > 0:
                        nScore += math.log(nProb)

            # Save the result to the file
            if cScore > nScore:
                result_file.write(f"{ID} is likely coding {cScore} {nScore}\n")
            else:
                result_file.write(f"{ID} is likely non-coding {cScore} {nScore}\n")

def getProbs(f1):
    f = open(f1)
    pMatrix = []
    for line in f:
        tmp = line.rstrip().split("\t")
        tmp = [float(i) for i in tmp]
        pMatrix.append(tmp)
    return pMatrix

def getSeq(filename):
    f = open(filename)
    id2seq = {}
    currkey = ""
    for line in f:
        if line.find(">") == 0:
            currkey = (line[1:].split("|")[0])
            id2seq[currkey] = ""
        else:
            id2seq[currkey] = id2seq[currkey] + line.rstrip()
    f.close()
    return id2seq

scoreModels()


# *Part Two*

# In[4]:


def scoreModels():
    codingMatrix = getProbs("./codingModel.tab")  # 64x64 matrix of probabilities for coding triplet mutations
    noncodingMatrix = getProbs("./noncodingModel.tab")  # 64x64 matrix of probabilities for non-coding triplet mutations
    
    id2ancestorSeq = getSeq("./Ancestor.fa")  # Reads the ancestor sequences
    id2spaciiSeq = getSeq("./Spacii_2100.fa")  # Reads the M. Spacii sequences
    allID = list(id2ancestorSeq.keys())  # The two sequences above share dictionary indexes

    true_labels = []
    likelihood_ratios = []

    for ID in allID:
        cScore = 0  # Log-likelihood score for the coding model
        nScore = 0  # Log-likelihood score for the non-coding model

        ancestorSeq = id2ancestorSeq[ID]
        spaciiSeq = id2spaciiSeq[ID]

        for i in range(0, len(ancestorSeq), 3):
            ancestorCodon = ancestorSeq[i:i+3]
            spaciiCodon = spaciiSeq[i:i+3]

            if ancestorCodon in modelCodons and spaciiCodon in modelCodons:
                ancestorIndex = modelCodons.index(ancestorCodon)
                spaciiIndex = modelCodons.index(spaciiCodon)

                cProb = codingMatrix[ancestorIndex][spaciiIndex]
                nProb = noncodingMatrix[ancestorIndex][spaciiIndex]

                if cProb > 0:
                    cScore += math.log(cProb)
                if nProb > 0:
                    nScore += math.log(nProb)

        # Compute the likelihood ratio in log-space
        likelihood_ratio = math.exp(cScore - nScore)

        # True labels (_c represents coding)
        if "SimSeq_c" in ID:
            true_labels.append(1)  # coding
        else:
            true_labels.append(0)  # non-coding
        
        likelihood_ratios.append(likelihood_ratio)

    return true_labels, likelihood_ratios


def getProbs(f1):
    with open(f1) as f:
        return [[float(i) for i in line.rstrip().split("\t")] for line in f]


def getSeq(filename):
    id2seq = {}
    with open(filename) as f:
        currkey = ""
        for line in f:
            if line.startswith(">"):
                currkey = line[1:].split("|")[0]
                id2seq[currkey] = ""
            else:
                id2seq[currkey] += line.rstrip()
    return id2seq


# Generate true labels and likelihood ratios
true_labels, likelihood_ratios = scoreModels()

# Calculate ROC curve
fpr, tpr, thresholds = roc_curve(true_labels, likelihood_ratios)
roc_auc = auc(fpr, tpr)

# Plot ROC curve
plt.figure()
plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (area = {roc_auc:.2f})')
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve using Likelihood Ratio')
plt.legend(loc="lower right")
plt.show()

