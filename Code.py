#########Python cod for generating random amino acid sequences

import random

# Define the 20 standard amino acids
amino_acids = "ACDEFGHIKLMNPQRSTVWY"

# Function to generate a random amino acid sequence
def generate_sequence(length):
    return ''.join(random.choice(amino_acids) for _ in range(length))

# Create a FASTA file
with open("SNP.fasta", "w") as fasta_file:
    # Generate 20 sequences
    for i in range(1, 21):
        # Sequence name Prt1, Prt2, ..., Prt20
        fasta_file.write(f">Prt{i}\n")
        # Generate a random sequence of length 200
        sequence = generate_sequence(200)
        fasta_file.write(sequence + "\n")

print("FASTA file generated: amino_acid_sequences.fasta")



#########Python cod for SNP discovery in amino acid sequences
from Bio import AlignIO
import matplotlib.pyplot as plt

# Define the path to your alignment file
alignment_file = "SNP.fasta"
output_file = "snp_plot.pdf"

print(f"Attempting to load file from: {alignment_file}")

# Load the amino acid alignment file
try:
    alignment = AlignIO.read(alignment_file, "fasta")
except FileNotFoundError:
    print(f"Error: File not found at {alignment_file}")
    exit()
except ValueError as e:
    print(f"Error: {e}")
    exit()

# Initialize lists to store SNP positions and corresponding sequences
snp_positions = []
snp_counts = {}  # Dictionary to store SNP counts at each position

# Iterate over the alignment columns to identify SNPs
for i in range(alignment.get_alignment_length()):
    column = alignment[:, i]
    if len(set(column)) > 1:
        snp_positions.append(i)
        snp_counts[i] = len(set(column))

# Prepare data for plotting
positions = list(snp_counts.keys())
counts = list(snp_counts.values())

# Plotting the SNP positions and their frequencies
fig, ax = plt.subplots(figsize=(12, 6))
ax.bar(positions, counts, color='blue', alpha=0.7)

ax.set_title("SNP Locations and Frequencies in Amino Acid Alignment")
ax.set_xlabel("SNP Position")
ax.set_ylabel("Number of SNPs")

# Adjust font size of tick labels
plt.xticks(positions, rotation=90, fontsize=8)
plt.tight_layout()

# Save the plot as a PDF with more than 300 DPI resolution
fig.savefig(output_file, dpi=300)
plt.show()
