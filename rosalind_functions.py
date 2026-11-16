#DNA 
def count_nucleotides(seq):
    a = seq.count("A")
    c = seq.count("C")
    g = seq.count("G")
    t = seq.count("T")
    return f"{a} {c} {g} {t}"
#RNA
def DNA_transcription(seq):
    RNA = seq.replace('T', 'U')
    return RNA
#REVC
def reverse_comp(seq):
    comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    reversed_seq = seq[::-1]
    return "".join(comp[base] for base in reversed_seq)
#GC
import sys

def read_fasta_text(text):
    sequences = {}
    name = None
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            name = line[1:].split()[0]
            sequences[name] = ""
        else:
            sequences[name] += line.upper()
    return sequences


def gc_content(seq):
    if len(seq) == 0:
        return 0
    return (seq.count("G") + seq.count("C")) / len(seq) * 100


def highest_gc_file(filename):
    with open(filename) as f:
        sequences = read_fasta_text(f.read())

    gc_values = {name: gc_content(seq) for name, seq in sequences.items()}

    best = max(gc_values, key=gc_values.get)
    print(best)
    print(f"{gc_values[best]:.3f}")
#HAMM 
def point_mut(seq1, seq2):
    return sum(a !=b for a, b in zip(seq1, seq2))
#PROT
from itertools import product

# Codon table in the same order as your table, using U
CODON_TABLE = {
    a + b + c: aa
    for (a, b, c), aa in zip(
        product("UCAG", repeat=3),
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    )
}

def translate(rna):
    """Translate an RNA string into a protein string, stopping at a Stop codon."""
    protein = []
    for i in range(0, len(rna) - 2, 3):
        aa = CODON_TABLE[rna[i:i + 3]]
        if aa == "*":
            break
        protein.append(aa)
    return "".join(protein)

def rna_to_protein(filename):
    """Read an RNA string from a file and return its protein string."""
    with open(filename) as f:
        rna = "".join(line.strip() for line in f).upper()
    return translate(rna)

print(rna_to_protein("rna"))
#SUBS
def find_locations(s, t):
    """Return the 1-based starting position of every occurrence of t in s."""
    n = len(t)
    return [i + 1 for i in range(len(s) - n + 1) if s[i:i + n] == t]

def read_two_strings(filename):
    """Read the first two non-blank lines of a file."""
    with open(filename) as f:
        lines = [line.strip() for line in f if line.strip()]
    return lines[0], lines[1]

if __name__ == "__main__":
    s, t = read_two_strings("subs") 
    print(" ".join(str(loc) for loc in find_locations(s, t)))

#PRTM

#REVP

#TRAN

