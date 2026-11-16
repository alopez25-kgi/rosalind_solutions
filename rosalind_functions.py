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

#SUBS

#PRTM

#REVP

#TRAN

