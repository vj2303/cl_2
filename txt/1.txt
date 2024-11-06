# DNA Sequence provided
dna_sequence = "ATGCATCGATCGTACGATCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"

dna_sequence = open("datasets\dna_sequence.txt").read().strip()
print("DNA sequence:", dna_sequence)

# Function to calculate GC content
def calculate_gc_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100

# Function to find motifs
def find_motifs(sequence, motif):
    positions = []
    for i in range(len(sequence) - len(motif) + 1):
        if sequence[i:i+len(motif)] == motif:
            positions.append(i)
    return positions

# Function to identify coding regions
def find_coding_regions(sequence):
    coding_regions = []
    i = 0
    while i < len(sequence) - 3:
        # Find start codon
        start = sequence.find("ATG", i)
        if start == -1:
            break  # No more start codons
        # Look for stop codon after start
        for stop_codon in ["TAA", "TAG", "TGA"]:
            stop = sequence.find(stop_codon, start + 3)
            if stop != -1 and (stop - start) % 3 == 0:
                coding_regions.append(sequence[start:stop+3])
                i = stop + 3  # Move to end of this coding region
                break
        else:
            i = start + 3  # Move to next potential start codon
    return coding_regions

# Calculate GC content
gc_content = calculate_gc_content(dna_sequence)

# Find motifs (example motifs TATA and GAATTC)
tata_positions = find_motifs(dna_sequence, "TATAAA")
ecori_positions = find_motifs(dna_sequence, "GAATTC")

# Identify coding regions
coding_regions = find_coding_regions(dna_sequence)

motif = "CGT"
motif_positions = [i for i in range(len(dna_sequence)) if dna_sequence.startswith(motif, i)]
print("Motif", motif, "found at positions:", motif_positions)

start_codon = "ATG"
stop_codons = ["TAA", "TAG", "TGA"]
coding_regions = []

for i in range(len(dna_sequence) - 2):
    if dna_sequence[i:i + 3] == start_codon:
        for j in range(i + 3, len(dna_sequence) - 2, 3):
            if dna_sequence[j:j + 3] in stop_codons:
                coding_regions.append((i, j + 3))
                break

if coding_regions:
    print("Coding regions found:")
    for start, end in coding_regions:
        print("Start:", start, "End:", end, "Sequence:", dna_sequence[start:end])
else:
    print("No coding regions found.")


#optional
print("Summary of DNA Sequence Analysis:")
print("GC Content:", gc_content, "%")
print("Motif positions:", motif_positions)
print("Coding regions:", coding_regions)
