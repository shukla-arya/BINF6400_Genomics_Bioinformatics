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

# *ORF Generating HMM*

# In[2]:


x = 100 # total emissions per seq
n_sequences = 250 # number of seqs to generate

# States: 0=Start, 1=Coding, 2=Stop
k_0 = np.array([1.0, 0.0, 0.0])  # always start in 'Start'

# Transition matrix (3x3)
A = np.array([
    [0.0, 1.0, 0.0],       # start to coding
    [0.01, 0.94, 0.05],    # coding to start (rare), coding, stop
    [1.0, 0.0, 0.0]        # stop to start (new ORF)
])


# Codon definitions
start_codon = ['ATG']
stop_codons = ['TAA', 'TAG', 'TGA']
all_codons = [a+b+c for a in "ATGC" for b in "ATGC" for c in "ATGC"]
non_special_codons = list(set(all_codons) - set(start_codon) - set(stop_codons))

# Emission probabilities
E_start = {codon: 0.0 for codon in all_codons}
E_start['ATG'] = 1.0  # "ATG" is the only possible start codon

E_coding = {codon: 1/len(non_special_codons) for codon in non_special_codons}
for codon in start_codon + stop_codons:
    E_coding[codon] = 0.0  # Remove start and stop codons from coding

E_stop = {codon: 0.0 for codon in all_codons}
for codon in stop_codons:
    E_stop[codon] = 1/len(stop_codons)  # Only emits stop codons with equal probability

E_k = [E_start, E_coding, E_stop]  # List of emission distributions

# Sampling functions
def sample_emission(state):
    codons = list(E_k[state].keys())
    probs = list(E_k[state].values())
    return np.random.choice(codons, p=probs)

def sample_next_state(current_state):
    return np.random.choice([0, 1, 2], p=A[current_state])

# Storage for sequences and ORF lengths
sequences = []
first_orf_lengths = []

# Generate sequences
for _ in range(n_sequences):
    state = np.random.choice([0, 1, 2], p=k_0)  # Start with a random state
    seq = []  # To store the generated sequence
    orf_len = 0
    in_orf = False
    first_orf_recorded = False

    for _ in range(x):  # Generate x codons (total emissions)
        emission = sample_emission(state)
        seq.append(emission)

        # Start a new ORF
        if state == 0:  # "Start" state
            in_orf = True
            orf_len = 1  # We start an ORF
        elif state == 1 and in_orf:  # "Coding" state
            orf_len += 1  # Increase ORF length
        elif state == 2 and in_orf:  # "Stop" state
            orf_len += 1  # ORF ends here
            if not first_orf_recorded:
                first_orf_lengths.append(orf_len)  # Record the first complete ORF
                first_orf_recorded = True
            in_orf = False  # ORF ended, stop counting
            orf_len = 0

            # After ORF ends, start a new ORF immediately
            state = 0  # Transition back to "Start" state

        state = sample_next_state(state)  # Move to the next state

    # If no ORF was found, append 0
    if not first_orf_recorded:
        first_orf_lengths.append(0)
    sequences.append(''.join(seq))  # Store the generated sequence

# Save sequences to a file
with open("Arya_Shukla_ORFoutput.txt", "w") as f:
    for seq in sequences:
        f.write(seq + "\n")

# Plot ORF length distribution
plt.figure(figsize=(10, 6))
plt.hist(first_orf_lengths, bins=range(0, max(first_orf_lengths)+5, 3),
         edgecolor='black', color='lightpink')
plt.title("Distribution of First Complete ORF Lengths (250 Sequences)")
plt.xlabel("ORF Length (codons)")
plt.ylabel("Frequency")
plt.grid(True)
plt.savefig("Arya_Shukla_ORF_Distribution.pdf")
plt.show()


# **Discussion**

# - 3 States: Start, Coding, Stop
# - The histogram plots the distribution of ORF lengths across sequences. The distribution is extremely positively skewed with a right tail. For a sequence to have a longer ORF, then the start codon ATG must begin earlier in the sequence. 
# - In randomly generated artificial sequences, there is no biological placement of ATG so there are less chances of finding a long ORF downstream. The smaller the sequence length, the greater the chance that the next codon is a stop. Longer ORFs seem to be rarer.
# - If an ORF completes before reaching the designated 100 emissions, then the sequence keeps generating until that amount of codons.
# - The parameters I chose are as follows: the start to coding state will always occur so that has a probability of 1. The coding state likely remains this way majority of the time, so I kept that at 0.94, while the chances of transitioning to a stop is 0.05 and back to a start is 0.01. Additionally, transitioning to a start after a stop will always occur, hence the 1.
