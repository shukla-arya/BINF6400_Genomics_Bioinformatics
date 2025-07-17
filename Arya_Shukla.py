#!/usr/bin/env python
# coding: utf-8

# # **PSet2 ORF Finding and Alignment**
# 
# **Genomics in Bioinformatics**
# <br>
# **Arya Shukla**  
# **20 February 2025**

# **Class Assignment**

# In[1]:


# Motif matrix as given in motifCode.txt
base_idx = {'A': 0, 'T': 1, 'C': 2, 'G': 3}  # letter to number matching for index
motif = [[.5,.5,.5,.5,0,0,0,0,0,2,-99,-99,.5],   # A
         [0,0,0,0,0,0,0,0,0,-99,2,-99,0],        # T
         [0,0,0,0,0,0,0,0,0,-99,-99,-99,0],      # C
         [.5,.5,.5,.5,0,0,0,0,0,.5,-99,2,0]]     # G


# There are 13 positions in the motif matrix, so the scoring is based on the subsequence that matches a motif consistent with the gene start. The values in the matrix are the position weight scores.
# <br>
# <br>
# The start codon is represented at positions 10, 11, and 12 in the motif matrix shown in the file. Out of the three options (ATG, TTG, and GTG), there is a T that appears in all of them. For that reason, we can eliminate positions 1-9, as they contain 0’s for the T nucleotide. The Shine Dalgarno motif occurs in bases A and G in positions 1-4.

# In[2]:


def scanSeq(seq):
    """
    This function scans a given DNA sequence for possible ORFs (open reading frames).
    It returns a list of ORF sequences, their start positions, and their lengths.
    
    The goal here is to determine if there are multiple ORFs in each contig.
    """
    # Initialize lists to store the sequence information
    orfs = []
    start_pos = []
    orf_lengths = []

    # Searching for ORFs that start with ATG or GTG
    for i in range(len(seq) - 2):
        if seq[i:i+3] == 'ATG' or seq[i:i+3] == 'GTG':
            # Found a potential ORF start
            for j in range(i + 3, len(seq) - 2, 3):
                codon = seq[j:j+3]
                if codon == 'TAA' or codon == 'TAG' or codon == 'TGA':  # Stop codons
                    orfs.append(seq[i:j+3])
                    start_pos.append(i)
                    orf_lengths.append(j + 3 - i)
                    break

    return orfs, start_pos, orf_lengths

def scoreMotif(seq, start_idx, window_size=6):
    """
    Scores a 13bp subsequence near the start codon. The start codon (ATG or GTG) is specifically aligned
    with positions 10-12 of the motif matrix. The window also includes the 9 preceding nucleotides and 1 after the start codon.
    """
    max_score = float('-inf')
    
    # Ensures that the start codon is within the motif matrix positions 10-12
    motif_start = start_idx - 9  # considers 9 bases before the start codon
    motif_end = motif_start + 13  # total length of 13bp
    
    if motif_start < 0 or motif_end > len(seq):
        return max_score

    subseq = seq[motif_start:motif_end]
    
    try:
        score = sum(motif[base_idx[base]][i] for i, base in enumerate(subseq))
        max_score = max(max_score, score)  # keep track of best score
    except KeyError:
        pass  # skip sequences with invalid bases
    
    return max_score

def identifyORFs(input_file, output_file, quality_score_cutoff=7.25, min_orf_length=60):
    """
    Identifies ORFs, scores them using the motif matrix, and outputs the results
    in FASTA format. Only ORFs with a length >= min_orf_length are included.
    The start codon is included in the length, but the stop codon is not.
    """
    with open(input_file, 'r') as f:
        contigs = f.read().split('>')[1:]  # split input into individual contig sections
        contigs = [contig.strip() for contig in contigs]

    results = []
    print(f"Number of contigs: {len(contigs)}")

    for i, contig in enumerate(contigs):
        contig_lines = contig.splitlines()
        contig_name = contig_lines[0]
        sequence = ''.join(contig_lines[1:]).upper()

        print(f"\nProcessing Contig {i+1}: {contig_name}, Length: {len(sequence)}")

        orfs, start_positions, orf_lengths = scanSeq(sequence)

        # Filter ORFs based on length, excluding the stop codon
        valid_orfs = []
        valid_start_positions = []
        valid_orf_lengths = []

        for j, orf in enumerate(orfs):
            # Include the start codon (ATG or GTG) in the length but not the stop codon
            # Start codon is methionine
            if orf_lengths[j] >= min_orf_length:
                valid_orfs.append(orf)
                valid_start_positions.append(start_positions[j])
                valid_orf_lengths.append(orf_lengths[j])

        print(f"Found ORFs for Contig {i+1}:")
        print(f"ORFs: {valid_orfs}")
        print(f"Start Positions: {valid_start_positions}")
        print(f"ORF Lengths: {valid_orf_lengths}")

        for j, orf in enumerate(valid_orfs):
            start_pos = valid_start_positions[j]

            # Use the new scanning function
            motif_score = scoreMotif(sequence, start_pos)

            print(f"Best Motif Score for ORF at Position {start_pos}: {motif_score}")

            if motif_score > quality_score_cutoff:
                result = f">contig {i+1}|Length {valid_orf_lengths[j]}|at position {start_pos}\n{orf}"
                results.append(result)

    print(f"\nNumber of ORFs Passed: {len(results)}")

    with open(output_file, 'w') as f:
        f.write("\n\n".join(results))

    print(f"Results saved to {output_file}")


# In[3]:


# Call the main function
identifyORFs('spaceSeq.fa', 'identified_orfs.fa')


# **Discussion**

# In this output, we can see the starting point of each ORF and how long each of them are. While there were ORF regions that exist with less than 60 base pairs (including the start codon but not the stop), they were excluded going forward. There are also ORF regions that are found within others in a nested format. For the purposes of our analysis, we considered all possible options because different start codons may take on various functions and are used for specific isoforms.
# <br>
# <br>
# The next step was to calculate the score where the start codon existed in positions 10-12 for the motif matrix. It is essential that the 13 base pair sequence has T in position 11. This is because its positional score in the motif matrix is not negative. For that reason, we take into consideration the 9 base pairs before the start codon and the 1 base pair after to account for all 13 positions used for calculation. It would NOT take into consideration the next 10 base pairs after the start codon nor would it consider the best group of 13 base pairs to produce the highest score.

# The motif represents the ribosomes ability to interact with the positions of the nucleotides getting fed into the ribosome. If it interacts well, then the tRNA will come in and recognize that spot. For the ribosome to have specific affinity, then the nucleotides must be an ATG or a GTG.

# Based on the criteria, we can deduce that a sequence with GTG will never produce a score of over 7.25 because the G in the motif matrix has a score of 0.5. Even though this is a valid start codon, we will not find any GTG ORFs in our output. Maybe we should reconsider the threshold of 7.25 going forward to determine if it is too high.

# There are also some contigs that result in 0 functional ORFs but some that result in multiple.

# **Simulating Sequences**

# In[4]:


import random
import numpy as np
import matplotlib.pyplot as plt


# In[9]:


# Global variables
l = 150  # sequence length
k = 15000  # number of sequences to generate
bases = ["A", "T", "C", "G"]
stop_codons = {"TAA", "TAG", "TGA"}

# Function to generate a random DNA sequence
def generate_sequence(length):
    return "".join(random.choices(bases, k=length))

# Function to scan a sequence for ORFs - repeated from above code
def scanSeq(seq):
    orf_lengths = []
    for i in range(len(seq) - 2):
        if seq[i:i+3] == 'ATG':
            for j in range(i + 3, len(seq) - 2, 3):
                codon = seq[j:j+3]
                if codon in stop_codons:
                    orf_lengths.append(j + 3 - i)
                    break
    return orf_lengths


# In[10]:


# Generate sequences and analyze ORF length distribution
all_orf_lengths = []
for _ in range(k):
    sequence = generate_sequence(l)
    all_orf_lengths.extend(scanSeq(sequence))

# Plot the distribution of ORF lengths
plt.figure(figsize=(10, 6))
plt.hist(all_orf_lengths, bins=30, color='skyblue', edgecolor='black', alpha=0.7)
plt.xlabel("ORF Length (bp)")
plt.ylabel("Frequency")
plt.title("Distribution of ORF Lengths in Simulated DNA Sequences")
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.show()


# In[11]:


# Generate sequences and count how many have ORFs > 60 bp
count_long_orfs = 0

for _ in range(k):
    sequence = generate_sequence(l)
    orf_lengths = scanSeq(sequence)
    
    # Check if at least one ORF is longer than 60 bp
    if any(length > 60 for length in orf_lengths):
        count_long_orfs += 1

# Compute the fraction
fraction = (count_long_orfs / k) * 100
print(f"Fraction of sequences with ORFs > 60 bp: {fraction:.4f}")

